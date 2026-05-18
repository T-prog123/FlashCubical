    using Clock = std::chrono::steady_clock;

    struct VertexKey {
        double value;
        std::size_t index;
        std::uint32_t code;
    };

    using H1Term = std::uint32_t;

    struct H1VertexOrderItem {
        double value;
        std::uint32_t code;
        VertexCode vertex;
    };

    static constexpr std::uint32_t no_rewrite_index =
        std::numeric_limits<std::uint32_t>::max();

    const double* values_;
    const std::uint32_t* value_codes_;
    std::size_t depth_;
    std::size_t rows_;
    std::size_t cols_;
    std::size_t plane_;
    std::size_t n_vertices_;
    std::size_t n_cubes_;
    std::size_t infinity_;
    bool include_apparent_zero_pairs_;
    bool compute_h1_;
    bool profile_h1_;
    std::uint8_t* edge_state_;
    std::uint8_t* square_state_;
    smart_core::UnionFind<VertexCode> h0_parent_;
    smart_core::UnionFind<CubeCode> dual_parent_;
    std::unique_ptr<VertexCode[]> edge_birth_vertex_;
    std::unique_ptr<EdgeCode[]> square_birth_edge_;
    std::unique_ptr<VertexCode[]> square_birth_vertex_;
    std::unique_ptr<std::uint8_t[]> square_second_rank_;
    std::unique_ptr<SquareCode[]> cube_birth_square_;
    std::unique_ptr<VertexCode[]> cube_birth_vertex_;
    std::unique_ptr<std::uint8_t[]> cube_second_rank_;
    const std::uint32_t* zero_lookup_masks_;
    zero_lookup::RunStats* zero_lookup_stats_;

    static double elapsed_ms(Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

    std::size_t vertex_id(std::size_t z, std::size_t y, std::size_t x) const {
        return (z * rows_ + y) * cols_ + x;
    }

    void vertex_coords(std::size_t v,
                       std::size_t& z,
                       std::size_t& y,
                       std::size_t& x) const {
        z = v / plane_;
        const std::size_t rem = v - z * plane_;
        y = rem / cols_;
        x = rem - y * cols_;
    }

    std::size_t cube_id(std::size_t z, std::size_t y, std::size_t x) const {
        return (z * (rows_ - 1) + y) * (cols_ - 1) + x;
    }

    std::size_t cube_source(std::size_t cube) const {
        const std::size_t w = cols_ - 1;
        const std::size_t h = rows_ - 1;
        const std::size_t z = cube / (h * w);
        const std::size_t rem = cube - z * h * w;
        const std::size_t y = rem / w;
        const std::size_t x = rem - y * w;
        return vertex_id(z, y, x);
    }

    VertexKey vertex_key(std::size_t v) const {
        const std::uint32_t code = value_codes_ != nullptr
            ? value_codes_[v]
            : std::numeric_limits<std::uint32_t>::max();
        return {values_[v], v, code};
    }

    static bool key_less(VertexKey a, VertexKey b) {
        if (a.code != std::numeric_limits<std::uint32_t>::max()) {
            if (a.code < b.code) {
                return true;
            }
            if (b.code < a.code) {
                return false;
            }
        } else {
            if (a.value < b.value) {
                return true;
            }
            if (b.value < a.value) {
                return false;
            }
        }
        return a.index < b.index;
    }

    static bool key_greater(VertexKey a, VertexKey b) {
        return key_less(b, a);
    }

    bool vertex_less(std::size_t a, std::size_t b) const {
        if (value_codes_ != nullptr) {
            const std::uint32_t ca = value_codes_[a];
            const std::uint32_t cb = value_codes_[b];
            return ca < cb || (ca == cb && a < b);
        }
        return values_[a] < values_[b] || (values_[a] == values_[b] && a < b);
    }

    bool vertex_greater(std::size_t a, std::size_t b) const {
        return vertex_less(b, a);
    }

    EdgeCode pack_edge(std::size_t source, unsigned dir) const {
        return static_cast<EdgeCode>(3u * source + dir);
    }

    std::size_t edge_source(EdgeCode e) const {
        return static_cast<std::size_t>(e / 3u);
    }

    unsigned edge_dir(EdgeCode e) const {
        return static_cast<unsigned>(e % 3u);
    }

    SquareCode pack_square(std::size_t source, unsigned orient) const {
        return static_cast<SquareCode>(3u * source + orient);
    }

    std::size_t square_source(SquareCode s) const {
        return static_cast<std::size_t>(s / 3u);
    }

    unsigned square_orient(SquareCode s) const {
        return static_cast<unsigned>(s % 3u);
    }

    bool edge_state_has(EdgeCode e, std::uint8_t bit) const {
        return (edge_state_[static_cast<std::size_t>(e)] & bit) != 0;
    }

    void edge_state_mark(EdgeCode e, std::uint8_t bit) const {
        edge_state_[static_cast<std::size_t>(e)] |= bit;
    }

    bool square_state_has(SquareCode s, std::uint8_t bit) const {
        return (square_state_[static_cast<std::size_t>(s)] & bit) != 0;
    }

    void square_state_mark(SquareCode s, std::uint8_t bit) const {
        square_state_[static_cast<std::size_t>(s)] |= bit;
    }

    std::pair<std::size_t, std::size_t> edge_endpoints(EdgeCode e) const {
        const std::size_t s = edge_source(e);
        const unsigned dir = edge_dir(e);
        if (dir == 0) {
            return {s, s + 1};
        }
        if (dir == 1) {
            return {s, s + cols_};
        }
        return {s, s + plane_};
    }

    EdgeCode edge_between(std::size_t a, std::size_t b) const {
        if (a + 1 == b) {
            return pack_edge(a, 0);
        }
        if (b + 1 == a) {
            return pack_edge(b, 0);
        }
        if (a + cols_ == b) {
            return pack_edge(a, 1);
        }
        if (b + cols_ == a) {
            return pack_edge(b, 1);
        }
        if (a + plane_ == b) {
            return pack_edge(a, 2);
        }
        if (b + plane_ == a) {
            return pack_edge(b, 2);
        }
        throw std::logic_error("vertices are not 3D grid neighbors");
    }

    EdgeCode lookup_local_edge_at_vertex(std::size_t v,
                                         std::uint8_t local) const {
        switch (local) {
        case zero_lookup::edge3_x_neg:
            return pack_edge(v - 1, 0);
        case zero_lookup::edge3_x_pos:
            return pack_edge(v, 0);
        case zero_lookup::edge3_y_neg:
            return pack_edge(v - cols_, 1);
        case zero_lookup::edge3_y_pos:
            return pack_edge(v, 1);
        case zero_lookup::edge3_z_neg:
            return pack_edge(v - plane_, 2);
        case zero_lookup::edge3_z_pos:
            return pack_edge(v, 2);
        default:
            throw std::logic_error("invalid local 3D edge id");
        }
    }

    SquareCode lookup_local_square_at_vertex(std::size_t v,
                                             std::uint8_t local) const {
        switch (local) {
        case 0:
            return pack_square(v - cols_ - 1, 0);
        case 1:
            return pack_square(v - cols_, 0);
        case 2:
            return pack_square(v - 1, 0);
        case 3:
            return pack_square(v, 0);
        case 4:
            return pack_square(v - plane_ - 1, 1);
        case 5:
            return pack_square(v - plane_, 1);
        case 6:
            return pack_square(v - 1, 1);
        case 7:
            return pack_square(v, 1);
        case 8:
            return pack_square(v - plane_ - cols_, 2);
        case 9:
            return pack_square(v - plane_, 2);
        case 10:
            return pack_square(v - cols_, 2);
        case 11:
            return pack_square(v, 2);
        default:
            throw std::logic_error("invalid local 3D square id");
        }
    }

    std::size_t lookup_local_cube_at_vertex(std::size_t v,
                                            std::uint8_t local) const {
        const bool xp = (local & 1u) != 0;
        const bool yp = (local & 2u) != 0;
        const bool zp = (local & 4u) != 0;
        const std::size_t source =
            v - (xp ? 0u : 1u) -
            (yp ? 0u : cols_) -
            (zp ? 0u : plane_);
        std::size_t z = 0;
        std::size_t y = 0;
        std::size_t x = 0;
        vertex_coords(source, z, y, x);
        return cube_id(z, y, x);
    }

    void ensure_edge_birth_vertices() {
        if (edge_birth_vertex_ != nullptr) {
            return;
        }
        edge_birth_vertex_.reset(new VertexCode[edge_state_storage_size(
            depth_, rows_, cols_)]);
        for_each_edge([&](EdgeCode e) {
            const auto [u, v] = edge_endpoints(e);
            edge_birth_vertex_[static_cast<std::size_t>(e)] =
                static_cast<VertexCode>(vertex_greater(u, v) ? u : v);
        });
    }

    std::size_t edge_birth_vertex(EdgeCode e) const {
        if (edge_birth_vertex_ != nullptr) {
            return edge_birth_vertex_[static_cast<std::size_t>(e)];
        }
        const auto [u, v] = edge_endpoints(e);
        return vertex_greater(u, v) ? u : v;
    }

    VertexKey edge_birth_key(EdgeCode e) const {
        return vertex_key(edge_birth_vertex(e));
    }

    double edge_birth_value(EdgeCode e) const {
        return values_[edge_birth_vertex(e)];
    }

    bool edge_less(EdgeCode a, EdgeCode b) const {
        const std::size_t va = edge_birth_vertex(a);
        const std::size_t vb = edge_birth_vertex(b);
        if (vertex_less(va, vb)) {
            return true;
        }
        if (vertex_less(vb, va)) {
            return false;
        }
        return a < b;
    }

    void square_boundary_edges(SquareCode s, EdgeCode* edges) const {
        const std::size_t source = square_source(s);
        const unsigned orient = square_orient(s);
        if (orient == 0) {
            edges[0] = pack_edge(source, 0);
            edges[1] = pack_edge(source + cols_, 0);
            edges[2] = pack_edge(source, 1);
            edges[3] = pack_edge(source + 1, 1);
        } else if (orient == 1) {
            edges[0] = pack_edge(source, 0);
            edges[1] = pack_edge(source + plane_, 0);
            edges[2] = pack_edge(source, 2);
            edges[3] = pack_edge(source + 1, 2);
        } else {
            edges[0] = pack_edge(source, 1);
            edges[1] = pack_edge(source + plane_, 1);
            edges[2] = pack_edge(source, 2);
            edges[3] = pack_edge(source + cols_, 2);
        }
    }

    void ensure_square_birth_edges() {
        if (square_birth_edge_ != nullptr) {
            return;
        }
        ensure_edge_birth_vertices();
        const std::size_t square_storage =
            square_state_storage_size(depth_, rows_, cols_);
        square_birth_edge_.reset(new EdgeCode[square_storage]);
        square_birth_vertex_.reset(new VertexCode[square_storage]);
        square_second_rank_.reset(new std::uint8_t[square_storage]);
        std::fill_n(square_second_rank_.get(), square_storage,
                    static_cast<std::uint8_t>(0xffu));
        std::size_t filled = 0;
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                if (local_square == 0xffu) {
                    throw std::logic_error(
                        "lookup square order metadata is incomplete");
                }
                const std::uint8_t local_edge =
                    entry.square_birth_local_edge[local_square];
                if (local_edge == 0xffu) {
                    throw std::logic_error(
                        "lookup square-birth metadata is incomplete");
                }
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                square_birth_edge_[static_cast<std::size_t>(s)] =
                    lookup_local_edge_at_vertex(v, local_edge);
                square_birth_vertex_[static_cast<std::size_t>(s)] =
                    static_cast<VertexCode>(v);
                square_second_rank_[static_cast<std::size_t>(s)] =
                    entry.square_second_rank[local_square];
                ++filled;
            }
        }
        if (filled != primal_square_count()) {
            throw std::logic_error(
                "lookup square-birth metadata did not cover all squares");
        }
    }

    EdgeCode square_birth_edge(SquareCode s) const {
        return square_birth_edge_[static_cast<std::size_t>(s)];
    }

    bool square_less(SquareCode a, SquareCode b) const {
        const std::size_t ai = static_cast<std::size_t>(a);
        const std::size_t bi = static_cast<std::size_t>(b);
        const VertexCode va = square_birth_vertex_[ai];
        const VertexCode vb = square_birth_vertex_[bi];
        if (vertex_less(va, vb)) {
            return true;
        }
        if (vertex_less(vb, va)) {
            return false;
        }
        const std::uint8_t ra = square_second_rank_[ai];
        const std::uint8_t rb = square_second_rank_[bi];
        if (ra != rb) {
            return ra < rb;
        }
        return a < b;
    }

    bool square_greater(SquareCode a, SquareCode b) const {
        return square_less(b, a);
    }

    bool square_before_in_dual_order(SquareCode a, SquareCode b) const {
        if (square_greater(a, b)) {
            return true;
        }
        if (square_greater(b, a)) {
            return false;
        }
        return a > b;
    }

    double square_birth_value(SquareCode s) const {
        return values_[square_birth_vertex_[static_cast<std::size_t>(s)]];
    }

    void ensure_cube_births() {
        if (cube_birth_square_ != nullptr || n_cubes_ == 0) {
            return;
        }
        ensure_square_birth_edges();
        cube_birth_square_.reset(new SquareCode[n_cubes_]);
        cube_birth_vertex_.reset(new VertexCode[n_cubes_]);
        cube_second_rank_.reset(new std::uint8_t[n_cubes_]);
        std::fill_n(cube_second_rank_.get(), n_cubes_,
                    static_cast<std::uint8_t>(0xffu));
        std::size_t filled = 0;
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            for (std::uint8_t i = 0; i < entry.cube_count; ++i) {
                const std::uint8_t local_cube = entry.cube_order_asc[i];
                if (local_cube == 0xffu) {
                    throw std::logic_error(
                        "lookup cube order metadata is incomplete");
                }
                const std::uint8_t local_square =
                    entry.cube_birth_local_square[local_cube];
                if (local_square == 0xffu) {
                    throw std::logic_error(
                        "lookup cube-birth metadata is incomplete");
                }
                const std::size_t cube =
                    lookup_local_cube_at_vertex(v, local_cube);
                cube_birth_square_[cube] =
                    lookup_local_square_at_vertex(v, local_square);
                cube_birth_vertex_[cube] = static_cast<VertexCode>(v);
                cube_second_rank_[cube] =
                    entry.cube_second_rank[local_cube];
                ++filled;
            }
        }
        if (filled != n_cubes_) {
            throw std::logic_error(
                "lookup cube-birth metadata did not cover all cubes");
        }
    }

    bool cube_less(std::size_t a, std::size_t b) const {
        const VertexCode va = cube_birth_vertex_[a];
        const VertexCode vb = cube_birth_vertex_[b];
        if (vertex_less(va, vb)) {
            return true;
        }
        if (vertex_less(vb, va)) {
            return false;
        }
        const std::uint8_t ra = cube_second_rank_[a];
        const std::uint8_t rb = cube_second_rank_[b];
        if (ra != rb) {
            return ra < rb;
        }
        return a < b;
    }

    bool cube_greater(std::size_t a, std::size_t b) const {
        return cube_less(b, a);
    }

    double cube_birth_value(std::size_t cube) const {
        return values_[cube_birth_vertex_[cube]];
    }

    double dual_birth_value(std::size_t dual_vertex) const {
        if (dual_vertex == infinity_) {
            return std::numeric_limits<double>::infinity();
        }
        return cube_birth_value(dual_vertex);
    }

    bool dual_older(std::size_t a, std::size_t b) const {
        if (a == infinity_) {
            return true;
        }
        if (b == infinity_) {
            return false;
        }
        if (cube_greater(a, b)) {
            return true;
        }
        if (cube_greater(b, a)) {
            return false;
        }
        return a > b;
    }

    std::size_t dual_younger_endpoint(std::size_t a, std::size_t b) const {
        return dual_older(a, b) ? b : a;
    }

    std::pair<std::size_t, std::size_t> dual_endpoints(SquareCode s) const {
        const std::size_t source = square_source(s);
        const unsigned orient = square_orient(s);
        std::size_t z = 0;
        std::size_t y = 0;
        std::size_t x = 0;
        vertex_coords(source, z, y, x);

        if (orient == 0) {
            if (z == 0) {
                return {cube_id(0, y, x), infinity_};
            }
            if (z + 1 == depth_) {
                return {cube_id(z - 1, y, x), infinity_};
            }
            return {cube_id(z - 1, y, x), cube_id(z, y, x)};
        }

        if (orient == 1) {
            if (y == 0) {
                return {cube_id(z, 0, x), infinity_};
            }
            if (y + 1 == rows_) {
                return {cube_id(z, y - 1, x), infinity_};
            }
            return {cube_id(z, y - 1, x), cube_id(z, y, x)};
        }

        if (x == 0) {
            return {cube_id(z, y, 0), infinity_};
        }
        if (x + 1 == cols_) {
            return {cube_id(z, y, x - 1), infinity_};
        }
        return {cube_id(z, y, x - 1), cube_id(z, y, x)};
    }

    template <class Fn>
    void for_each_edge(Fn&& fn) const {
        for (std::size_t z = 0; z < depth_; ++z) {
            for (std::size_t y = 0; y < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x + 1 < cols_; ++x) {
                    fn(pack_edge(base + x, 0));
                }
            }
        }
        for (std::size_t z = 0; z < depth_; ++z) {
            for (std::size_t y = 0; y + 1 < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x < cols_; ++x) {
                    fn(pack_edge(base + x, 1));
                }
            }
        }
        for (std::size_t z = 0; z + 1 < depth_; ++z) {
            for (std::size_t y = 0; y < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x < cols_; ++x) {
                    fn(pack_edge(base + x, 2));
                }
            }
        }
    }

    template <class Fn>
    void for_each_square(Fn&& fn) const {
        for (std::size_t z = 0; z < depth_; ++z) {
            for (std::size_t y = 0; y + 1 < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x + 1 < cols_; ++x) {
                    fn(pack_square(base + x, 0));
                }
            }
        }
        for (std::size_t z = 0; z + 1 < depth_; ++z) {
            for (std::size_t y = 0; y < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x + 1 < cols_; ++x) {
                    fn(pack_square(base + x, 1));
                }
            }
        }
        for (std::size_t z = 0; z + 1 < depth_; ++z) {
            for (std::size_t y = 0; y + 1 < rows_; ++y) {
                const std::size_t base = vertex_id(z, y, 0);
                for (std::size_t x = 0; x < cols_; ++x) {
                    fn(pack_square(base + x, 2));
                }
            }
        }
    }

    std::size_t primal_square_count() const {
        return depth_ * (rows_ - 1) * (cols_ - 1) +
               (depth_ - 1) * rows_ * (cols_ - 1) +
               (depth_ - 1) * (rows_ - 1) * cols_;
    }

    int incident_squares_of_edge(EdgeCode e, SquareCode* squares) const {
        const std::size_t source = edge_source(e);
        const unsigned dir = edge_dir(e);
        std::size_t z = 0;
        std::size_t y = 0;
        std::size_t x = 0;
        vertex_coords(source, z, y, x);
        int count = 0;

        if (dir == 0) {
            if (y + 1 < rows_) {
                squares[count++] = pack_square(source, 0);
            }
            if (y > 0) {
                squares[count++] = pack_square(source - cols_, 0);
            }
            if (z + 1 < depth_) {
                squares[count++] = pack_square(source, 1);
            }
            if (z > 0) {
                squares[count++] = pack_square(source - plane_, 1);
            }
        } else if (dir == 1) {
            if (x + 1 < cols_) {
                squares[count++] = pack_square(source, 0);
            }
            if (x > 0) {
                squares[count++] = pack_square(source - 1, 0);
            }
            if (z + 1 < depth_) {
                squares[count++] = pack_square(source, 2);
            }
            if (z > 0) {
                squares[count++] = pack_square(source - plane_, 2);
            }
        } else {
            if (x + 1 < cols_) {
                squares[count++] = pack_square(source, 1);
            }
            if (x > 0) {
                squares[count++] = pack_square(source - 1, 1);
            }
            if (y + 1 < rows_) {
                squares[count++] = pack_square(source, 2);
            }
            if (y > 0) {
                squares[count++] = pack_square(source - cols_, 2);
            }
        }

        return count;
    }

    bool is_edge_square_apparent(EdgeCode e, SquareCode s) const {
        if (square_birth_edge(s) != e) {
            return false;
        }

        SquareCode cofacets[4];
        const int count = incident_squares_of_edge(e, cofacets);
        if (count == 0) {
            return false;
        }

        SquareCode earliest = cofacets[0];
        for (int i = 1; i < count; ++i) {
            if (square_less(cofacets[i], earliest)) {
                earliest = cofacets[i];
            }
        }
        return earliest == s;
    }

    std::size_t mark_primal_h1_apparent_edges() {
        ensure_square_birth_edges();
        std::size_t count = 0;
        for_each_square([&](SquareCode s) {
            const EdgeCode e = square_birth_edge(s);
            if (is_edge_square_apparent(e, s)) {
                edge_state_mark(e, edge_state_h1_apparent | edge_state_h1_positive);
                ++count;
            }
        });
        return count;
    }

    std::size_t mark_dual_h1_apparent_squares() {
        ensure_square_birth_edges();
        std::size_t count = 0;
        for_each_edge([&](EdgeCode e) {
            SquareCode cofacets[4];
            const int cofacet_count = incident_squares_of_edge(e, cofacets);
            if (cofacet_count == 0) {
                return;
            }

            SquareCode earliest = cofacets[0];
            for (int i = 1; i < cofacet_count; ++i) {
                if (square_less(cofacets[i], earliest)) {
                    earliest = cofacets[i];
                }
            }
            if (square_birth_edge(earliest) == e) {
                square_state_mark(earliest,
                                  square_state_dual_h1_apparent |
                                      square_state_dual_h1_positive);
                ++count;
            }
        });
        return count;
    }

    ApparentPrepassStats mark_edge_square_apparent_pairs() {
        ensure_square_birth_edges();
        ApparentPrepassStats stats;
        for_each_edge([&](EdgeCode e) {
            SquareCode cofacets[4];
            const int cofacet_count = incident_squares_of_edge(e, cofacets);
            if (cofacet_count == 0) {
                return;
            }

            SquareCode earliest = cofacets[0];
            for (int i = 1; i < cofacet_count; ++i) {
                if (square_less(cofacets[i], earliest)) {
                    earliest = cofacets[i];
                }
            }
            if (square_birth_edge(earliest) == e) {
                edge_state_mark(e, edge_state_h1_apparent |
                                   edge_state_h1_positive);
                square_state_mark(earliest,
                                  square_state_dual_h1_apparent |
                                      square_state_dual_h1_positive);
                ++stats.primal_h1_apparent_pairs;
                ++stats.dual_h1_apparent_pairs;
                ++stats.primal_h1_positive_edges;
                ++stats.dual_h1_positive_squares;
            }
        });
        return stats;
    }

    static void sort_edges_ascending(EdgeCode* edges, int count) {
        for (int i = 1; i < count; ++i) {
            const EdgeCode x = edges[i];
            int j = i - 1;
            while (j >= 0 && edges[j] > x) {
                edges[j + 1] = edges[j];
                --j;
            }
            edges[j + 1] = x;
        }
    }

    void sort_squares_in_dual_order(SquareCode* squares, int count) const {
        for (int i = 1; i < count; ++i) {
            const SquareCode x = squares[i];
            int j = i - 1;
            while (j >= 0 && square_before_in_dual_order(x, squares[j])) {
                squares[j + 1] = squares[j];
                --j;
            }
            squares[j + 1] = x;
        }
    }
