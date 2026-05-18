#pragma once

#include "algorithms/detail/smart_h0_2d_common.hpp"
#include "zero_lookup/runtime/zero_table.hpp"

namespace smart_h0 {

class H1DualComputer {
public:
    H1DualComputer(const double* values, std::size_t rows, std::size_t cols,
                   std::uint8_t* edge_state = nullptr,
                   const std::uint32_t* value_codes = nullptr)
        : values_(values),
          value_codes_(value_codes),
          rows_(rows),
          cols_(cols),
          n_vertices_(0),
          n_squares_(0),
          infinity_(0),
          edge_state_(edge_state),
          zero_lookup_masks_(nullptr),
          square_birth_edge_(nullptr),
          square_birth_vertex_(nullptr),
          square_second_rank_(nullptr) {
        if (values_ == nullptr) {
            throw std::invalid_argument("image value pointer is null");
        }
        if (rows_ == 0 || cols_ == 0) {
            throw std::invalid_argument("image dimensions must be positive");
        }
        if (rows_ > (std::numeric_limits<std::size_t>::max() / cols_)) {
            throw std::overflow_error("image dimensions overflow size_t");
        }
        n_vertices_ = rows_ * cols_;
        if (n_vertices_ > std::numeric_limits<VertexCode>::max()) {
            throw std::overflow_error("image has too many vertices for packed vertex ids");
        }
        if (rows_ < 2 || cols_ < 2) {
            return;
        }
        n_squares_ = (rows_ - 1) * (cols_ - 1);
        infinity_ = n_squares_;
    }

    void mark_apparent_positive_edges_only() {
        if (n_squares_ == 0 || edge_state_ == nullptr) {
            return;
        }

        ensure_square_births();

        scan_dual_vertices_for_apparent_edges(
            [this](std::size_t x, EdgeCode e) {
                (void)x;
                edge_state_mark(e, edge_state_h1_apparent | edge_state_h1_positive);
            });
    }

    H1StagedRun start_apparent_pass() {
        H1StagedRun staged;

        if (n_squares_ == 0) {
            return staged;
        }

        const auto init_start = Clock::now();
        ensure_square_births();
        parent_.reset(n_squares_ + 1);
        staged.result.pairs.reserve(n_squares_ / 4 + 1);
        staged.result.timing.parent_initialization_ms =
            elapsed_ms(init_start, Clock::now());

        staged.result.dual_vertices = n_squares_ + 1;
        staged.result.primal_edges = edge_count();
        staged.apparent_edge_mask.reset(new std::uint8_t[n_vertices_]());

        const auto scan_start = Clock::now();
        scan_dual_apparent_pairs(staged.apparent_edge_mask.get(),
                                 staged.result);
        staged.apparent_edges_applied = true;
        staged.result.timing.apparent_pair_scan_ms =
            elapsed_ms(scan_start, Clock::now());
        staged.result.timing.apparent_pair_work_ms =
            staged.result.timing.apparent_pair_scan_ms;

        return staged;
    }

    H1StagedRun start_zero_lookup_pass(const std::uint8_t* zero_lookup_masks,
                                       zero_lookup::RunStats* zero_lookup_stats) {
        H1StagedRun staged;
        zero_lookup_masks_ = zero_lookup_masks;

        if (n_squares_ == 0) {
            return staged;
        }

        const auto init_start = Clock::now();
        ensure_square_births();
        parent_.reset(n_squares_ + 1);
        staged.result.pairs.reserve(n_squares_ / 4 + 1);
        staged.result.timing.parent_initialization_ms =
            elapsed_ms(init_start, Clock::now());

        staged.result.dual_vertices = n_squares_ + 1;
        staged.result.primal_edges = edge_count();
        staged.apparent_edge_mask.reset(new std::uint8_t[n_vertices_]());

        const auto scan_start = Clock::now();
        scan_dual_lookup_positive_edges(zero_lookup_masks,
                                        staged.apparent_edge_mask.get(),
                                        staged.result,
                                        zero_lookup_stats);
        staged.apparent_edges_applied = true;
        staged.result.timing.apparent_pair_scan_ms =
            elapsed_ms(scan_start, Clock::now());
        staged.result.timing.apparent_pair_work_ms =
            staged.result.timing.apparent_pair_scan_ms;

        return staged;
    }

    H1Result finish_after_h0(H1StagedRun staged,
                             const std::uint8_t* prefilled_edge_owner_mask = nullptr,
                             std::vector<VertexCode>* prefilled_active_vertices = nullptr) {
        H1Result& result = staged.result;

        if (n_squares_ == 0) {
            return result;
        }

        const auto init_start = Clock::now();
        if (parent_.empty()) {
            parent_.reset(n_squares_ + 1);
        }
        if (result.pairs.capacity() == 0) {
            result.pairs.reserve(n_squares_ / 4 + 1);
        }
        result.timing.parent_initialization_ms +=
            elapsed_ms(init_start, Clock::now());

        const auto fill_start = Clock::now();
        if (!staged.apparent_edges_applied) {
            apply_apparent_edges_from_mask(staged.apparent_edge_mask.get(), result);
        }

        std::unique_ptr<std::uint8_t[]> edge_owner_mask_storage;
        const std::uint8_t* edge_owner_mask = prefilled_edge_owner_mask;
        std::vector<VertexCode> active_vertices_storage;
        std::vector<VertexCode>* active_vertices = prefilled_active_vertices;

        if (edge_owner_mask == nullptr || active_vertices == nullptr) {
            edge_owner_mask_storage.reset(new std::uint8_t[n_vertices_]());
            active_vertices_storage.reserve(n_vertices_ / 4 + 1);
            fill_remaining_dual_edges(staged.apparent_edge_mask.get(),
                                      edge_owner_mask_storage.get(),
                                      active_vertices_storage,
                                      result,
                                      staged.apparent_edges_applied);
            edge_owner_mask = edge_owner_mask_storage.get();
            active_vertices = &active_vertices_storage;
        } else {
            const std::size_t candidate_edges =
                count_mask_edges(*active_vertices, edge_owner_mask);
            result.skipped_known_h0_negative_edges =
                edge_count() - result.apparent_edges - candidate_edges;
        }

        result.timing.apparent_pair_scan_ms +=
            elapsed_ms(fill_start, Clock::now());
        result.timing.apparent_pair_work_ms = result.timing.apparent_pair_scan_ms;
        result.active_edge_owner_vertices = active_vertices->size();

        const bool use_lookup_reverse_order =
            zero_lookup_masks_ != nullptr &&
            edge_owner_mask == prefilled_edge_owner_mask &&
            active_vertices == prefilled_active_vertices;

        const auto sort_start = Clock::now();
        if (use_lookup_reverse_order) {
            validate_lookup_reverse_owner_order(*active_vertices,
                                                edge_owner_mask);
        } else {
            sort_active_vertices_descending(*active_vertices, edge_owner_mask);
        }
        result.timing.sorting_ms = elapsed_ms(sort_start, Clock::now());

        const auto sweep_start = Clock::now();
        if (use_lookup_reverse_order) {
            sweep_dual_edges_lookup_reverse(*active_vertices, edge_owner_mask,
                                            result);
        } else {
            sweep_dual_edges(*active_vertices, edge_owner_mask, result);
        }
        result.timing.union_find_sweep_ms = elapsed_ms(sweep_start, Clock::now());

        const auto final_start = Clock::now();
        result.timing.finalization_ms = elapsed_ms(final_start, Clock::now());
        result.timing.total_ms =
            result.timing.parent_initialization_ms +
            result.timing.apparent_pair_scan_ms +
            result.timing.sorting_ms +
            result.timing.union_find_sweep_ms +
            result.timing.finalization_ms;

        result.memory.parent_bytes = (n_squares_ + 1) * sizeof(VertexCode);
        result.memory.non_apparent_vertex_mask_bytes = n_vertices_ * sizeof(std::uint8_t);
        result.memory.active_vertices_bytes = active_vertices->capacity() * sizeof(VertexCode);
        result.memory.square_birth_bytes =
            n_squares_ * ((square_birth_edge_ != nullptr
                               ? sizeof(EdgeCode)
                               : 0) +
                          sizeof(VertexCode) +
                          (square_second_rank_ != nullptr
                               ? sizeof(std::uint8_t)
                               : 0));
        result.memory.persistence_pairs_bytes =
            result.pairs.capacity() * sizeof(PersistencePair);

        return result;
    }

    H1Result run() {
        return finish_after_h0(start_apparent_pass());
    }

private:
    using Clock = std::chrono::steady_clock;

    const double* values_;
    const std::uint32_t* value_codes_;
    std::size_t rows_;
    std::size_t cols_;
    std::size_t n_vertices_;
    std::size_t n_squares_;
    std::size_t infinity_;
    std::uint8_t* edge_state_;
    const std::uint8_t* zero_lookup_masks_;
    smart_core::UnionFind<VertexCode> parent_;
    struct VertexKey {
        double value;
        std::size_t index;
        std::uint32_t code;
    };

    std::unique_ptr<EdgeCode[]> square_birth_edge_;
    std::unique_ptr<VertexCode[]> square_birth_vertex_;
    std::unique_ptr<std::uint8_t[]> square_second_rank_;

    static constexpr std::uint8_t edge_left_bit = 1u << 0;
    static constexpr std::uint8_t edge_right_bit = 1u << 1;
    static constexpr std::uint8_t edge_up_bit = 1u << 2;
    static constexpr std::uint8_t edge_down_bit = 1u << 3;

    static double elapsed_ms(Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

    void ensure_square_births() {
        if (square_birth_vertex_ != nullptr) {
            return;
        }
        square_birth_vertex_.reset(new VertexCode[n_squares_]);
        if (zero_lookup_masks_ != nullptr) {
            square_second_rank_.reset(new std::uint8_t[n_squares_]);
            fill_square_births_from_lookup();
            return;
        }
        square_birth_edge_.reset(new EdgeCode[n_squares_]);

        for (std::size_t r = 0; r + 1 < rows_; ++r) {
            const std::size_t top = r * cols_;
            const std::size_t bottom = top + cols_;
            const std::size_t square_base = r * (cols_ - 1);
            for (std::size_t c = 0; c + 1 < cols_; ++c) {
                const std::size_t v00 = top + c;
                const std::size_t v01 = v00 + 1;
                const std::size_t v10 = bottom + c;
                const std::size_t v11 = v10 + 1;

                std::size_t youngest = v00;
                if (vertex_greater(v01, youngest)) {
                    youngest = v01;
                }
                if (vertex_greater(v10, youngest)) {
                    youngest = v10;
                }
                if (vertex_greater(v11, youngest)) {
                    youngest = v11;
                }

                EdgeCode a = 0;
                EdgeCode b = 0;
                if (youngest == v00) {
                    a = pack_edge(v00, 0);
                    b = pack_edge(v00, 1);
                } else if (youngest == v01) {
                    a = pack_edge(v00, 0);
                    b = pack_edge(v01, 1);
                } else if (youngest == v10) {
                    a = pack_edge(v10, 0);
                    b = pack_edge(v00, 1);
                } else {
                    a = pack_edge(v10, 0);
                    b = pack_edge(v01, 1);
                }

                const std::size_t square = square_base + c;
                square_birth_edge_[square] = (a > b) ? a : b;
                square_birth_vertex_[square] = static_cast<VertexCode>(youngest);
            }
        }
    }

    void fill_square_births_from_lookup() {
        for (std::size_t r = 0; r + 1 < rows_; ++r) {
            const std::size_t top = r * cols_;
            const std::size_t bottom = top + cols_;
            const std::size_t square_base = r * (cols_ - 1);
            for (std::size_t c = 0; c + 1 < cols_; ++c) {
                const std::size_t v00 = top + c;
                const std::size_t v01 = v00 + 1;
                const std::size_t v10 = bottom + c;
                const std::size_t v11 = v10 + 1;

                std::size_t v = v00;
                std::uint8_t local_square = zero_lookup::square2_se;
                if (vertex_greater(v01, v)) {
                    v = v01;
                    local_square = zero_lookup::square2_sw;
                }
                if (vertex_greater(v10, v)) {
                    v = v10;
                    local_square = zero_lookup::square2_ne;
                }
                if (vertex_greater(v11, v)) {
                    v = v11;
                    local_square = zero_lookup::square2_nw;
                }

                const std::size_t square = square_base + c;
                square_birth_vertex_[square] = static_cast<VertexCode>(v);
                const zero_lookup::Entry2D& entry =
                    zero_lookup::lookup2d(zero_lookup_masks_[v]);
                if (entry.square_second_rank[local_square] == 0xffu) {
                    throw std::logic_error(
                        "lookup square birth rank missing for local square");
                }
                square_second_rank_[square] =
                    entry.square_second_rank[local_square];
            }
        }
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

    static void sort_edges_descending(EdgeCode* edges, int count) {
        for (int i = 1; i < count; ++i) {
            const EdgeCode x = edges[i];
            int j = i - 1;
            while (j >= 0 && edges[j] < x) {
                edges[j + 1] = edges[j];
                --j;
            }
            edges[j + 1] = x;
        }
    }

    void sort_active_vertices_descending(std::vector<VertexCode>& active_vertices,
                                         const std::uint8_t* edge_owner_mask) const {
        if (value_codes_ == nullptr) {
            std::sort(active_vertices.begin(), active_vertices.end(),
                      [this](VertexCode a, VertexCode b) {
                          return key_less(vertex_key(b), vertex_key(a));
                      });
            return;
        }

        std::uint32_t max_code = 0;
        for (VertexCode v : active_vertices) {
            if (value_codes_[v] > max_code) {
                max_code = value_codes_[v];
            }
        }
        if (max_code <= 10000000u && edge_owner_mask != nullptr) {
            std::vector<std::uint32_t> offsets(static_cast<std::size_t>(max_code) + 1u, 0);
            for (VertexCode v : active_vertices) {
                ++offsets[value_codes_[v]];
            }

            std::uint32_t running = 0;
            for (std::uint32_t code = max_code + 1u; code > 0; --code) {
                const std::uint32_t bucket = code - 1u;
                const std::uint32_t bucket_size = offsets[bucket];
                offsets[bucket] = running;
                running += bucket_size;
            }

            std::vector<VertexCode> sorted(active_vertices.size());
            for (std::size_t v = n_vertices_; v > 0; --v) {
                const std::size_t vertex = v - 1;
                if (edge_owner_mask[vertex] != 0) {
                    const std::uint32_t code = value_codes_[vertex];
                    sorted[offsets[code]++] = static_cast<VertexCode>(vertex);
                }
            }
            active_vertices.swap(sorted);
            return;
        }

        std::vector<std::uint64_t> keys;
        keys.reserve(active_vertices.size());
        for (VertexCode v : active_vertices) {
            keys.push_back((static_cast<std::uint64_t>(value_codes_[v]) << 32) |
                           static_cast<std::uint64_t>(v));
        }
        std::sort(keys.begin(), keys.end());
        for (std::size_t i = 0, n = keys.size(); i < n; ++i) {
            active_vertices[i] = static_cast<VertexCode>(keys[n - 1 - i]);
        }
    }

    void validate_lookup_reverse_owner_order(
        const std::vector<VertexCode>& active_vertices,
        const std::uint8_t* edge_owner_mask) const {
        if (n_vertices_ > 4096 || active_vertices.size() < 2) {
            return;
        }
        std::vector<VertexCode> sorted = active_vertices;
        sort_active_vertices_descending(sorted, edge_owner_mask);
        for (std::size_t i = 0, n = active_vertices.size(); i < n; ++i) {
            if (sorted[i] != active_vertices[n - 1 - i]) {
                throw std::logic_error(
                    "lookup H1 active owners are not reverse-sorted");
            }
        }
    }

    std::size_t row(std::size_t v) const {
        return v / cols_;
    }

    std::size_t col(std::size_t v) const {
        return v % cols_;
    }

    EdgeCode pack_edge(std::size_t source, unsigned dir) const {
        return (static_cast<EdgeCode>(source) << 1) | static_cast<EdgeCode>(dir);
    }

    std::size_t edge_source(EdgeCode e) const {
        return static_cast<std::size_t>(e >> 1);
    }

    unsigned edge_dir(EdgeCode e) const {
        return static_cast<unsigned>(e & 1u);
    }

    bool edge_state_has(EdgeCode e, std::uint8_t bit) const {
        return edge_state_ != nullptr &&
               (edge_state_[static_cast<std::size_t>(e)] & bit) != 0;
    }

    void edge_state_mark(EdgeCode e, std::uint8_t bit) const {
        if (edge_state_ != nullptr) {
            edge_state_[static_cast<std::size_t>(e)] |= bit;
        }
    }

    std::pair<std::size_t, std::size_t> endpoints(EdgeCode e) const {
        const std::size_t s = edge_source(e);
        if (edge_dir(e) == 0) {
            return {s, s + 1};
        }
        return {s, s + cols_};
    }

    std::size_t edge_count() const {
        return rows_ * (cols_ - 1) + (rows_ - 1) * cols_;
    }

    double edge_value(EdgeCode e) const {
        return edge_birth_key(e).value;
    }

    VertexKey edge_birth_key(EdgeCode e) const {
        const auto [u, v] = endpoints(e);
        return vertex_greater(u, v) ? vertex_key(u) : vertex_key(v);
    }

    bool edge_before_in_dual_order(EdgeCode a, EdgeCode b) const {
        const VertexKey ka = edge_birth_key(a);
        const VertexKey kb = edge_birth_key(b);
        if (key_greater(ka, kb)) {
            return true;
        }
        if (key_greater(kb, ka)) {
            return false;
        }
        return a > b;
    }

    std::size_t square_id(std::size_t r, std::size_t c) const {
        return r * (cols_ - 1) + c;
    }

    std::size_t square_row(std::size_t square) const {
        return square / (cols_ - 1);
    }

    std::size_t square_col(std::size_t square) const {
        return square % (cols_ - 1);
    }

    std::pair<EdgeCode, VertexCode> compute_square_birth(std::size_t square) const {
        const std::size_t r = square / (cols_ - 1);
        const std::size_t c = square % (cols_ - 1);
        const std::size_t v00 = r * cols_ + c;
        const std::size_t v01 = v00 + 1;
        const std::size_t v10 = v00 + cols_;
        const std::size_t v11 = v10 + 1;

        std::size_t youngest = v00;
        if (vertex_greater(v01, youngest)) {
            youngest = v01;
        }
        if (vertex_greater(v10, youngest)) {
            youngest = v10;
        }
        if (vertex_greater(v11, youngest)) {
            youngest = v11;
        }

        EdgeCode a = 0;
        EdgeCode b = 0;
        if (youngest == v00) {
            a = pack_edge(v00, 0);
            b = pack_edge(v00, 1);
        } else if (youngest == v01) {
            a = pack_edge(v00, 0);
            b = pack_edge(v01, 1);
        } else if (youngest == v10) {
            a = pack_edge(v10, 0);
            b = pack_edge(v00, 1);
        } else {
            a = pack_edge(v10, 0);
            b = pack_edge(v01, 1);
        }

        return {(a > b) ? a : b, static_cast<VertexCode>(youngest)};
    }

    double dual_birth(std::size_t dual_vertex) const {
        if (dual_vertex == infinity_) {
            return std::numeric_limits<double>::infinity();
        }
        return values_[square_birth_vertex_[dual_vertex]];
    }

    bool dual_older(std::size_t a, std::size_t b) const {
        if (a == infinity_) {
            return true;
        }
        if (b == infinity_) {
            return false;
        }
        const VertexCode va = square_birth_vertex_[a];
        const VertexCode vb = square_birth_vertex_[b];
        if (vertex_greater(va, vb)) {
            return true;
        }
        if (vertex_greater(vb, va)) {
            return false;
        }
        if (square_second_rank_ != nullptr) {
            const std::uint8_t ra = square_second_rank_[a];
            const std::uint8_t rb = square_second_rank_[b];
            if (ra > rb) {
                return true;
            }
            if (rb > ra) {
                return false;
            }
            return a > b;
        }
        const EdgeCode ea = square_birth_edge_[a];
        const EdgeCode eb = square_birth_edge_[b];
        if (ea > eb) {
            return true;
        }
        if (eb > ea) {
            return false;
        }
        return a > b;
    }

    std::size_t dual_younger_endpoint(std::size_t a, std::size_t b) const {
        return dual_older(a, b) ? b : a;
    }

    std::size_t find(std::size_t v) {
        return parent_.find(static_cast<VertexCode>(v));
    }

    std::uint8_t owned_edge_bit(EdgeCode e, std::size_t owner) const {
        const std::size_t s = edge_source(e);
        if (edge_dir(e) == 0) {
            return (owner == s) ? edge_right_bit : edge_left_bit;
        }
        return (owner == s) ? edge_down_bit : edge_up_bit;
    }

    static std::uint8_t local_edge_owned_bit(std::uint8_t local) {
        switch (local) {
        case zero_lookup::edge2_n:
            return edge_up_bit;
        case zero_lookup::edge2_e:
            return edge_right_bit;
        case zero_lookup::edge2_s:
            return edge_down_bit;
        case zero_lookup::edge2_w:
            return edge_left_bit;
        default:
            return 0;
        }
    }

    void mark_owned_edge(std::uint8_t* edge_owner_mask,
                         std::vector<VertexCode>& active_vertices,
                         std::size_t owner,
                         EdgeCode e) const {
        const std::uint8_t bit = owned_edge_bit(e, owner);
        if (edge_owner_mask[owner] == 0) {
            active_vertices.push_back(static_cast<VertexCode>(owner));
        }
        edge_owner_mask[owner] |= bit;
    }

    void mark_remaining_dual_edge(std::uint8_t* edge_owner_mask,
                                  std::vector<VertexCode>& active_vertices,
                                  EdgeCode e) const {
        const auto [u, v] = endpoints(e);
        const std::size_t owner =
            vertex_greater(u, v) ? u : v;
        mark_owned_edge(edge_owner_mask, active_vertices, owner, e);
    }

    void scan_primal_edges(std::uint8_t* edge_owner_mask,
                           std::vector<VertexCode>& active_vertices) const {
        auto handle_edge = [&](EdgeCode e) {
            mark_remaining_dual_edge(edge_owner_mask, active_vertices, e);
        };

        for (std::size_t r = 0; r < rows_; ++r) {
            const std::size_t base = r * cols_;
            for (std::size_t c = 0; c + 1 < cols_; ++c) {
                handle_edge(pack_edge(base + c, 0));
            }
        }
        for (std::size_t r = 0; r + 1 < rows_; ++r) {
            const std::size_t base = r * cols_;
            for (std::size_t c = 0; c < cols_; ++c) {
                handle_edge(pack_edge(base + c, 1));
            }
        }
    }

    template <class Fn>
    void for_each_primal_edge(Fn&& fn) const {
        for (std::size_t r = 0; r < rows_; ++r) {
            const std::size_t base = r * cols_;
            for (std::size_t c = 0; c + 1 < cols_; ++c) {
                fn(pack_edge(base + c, 0));
            }
        }
        for (std::size_t r = 0; r + 1 < rows_; ++r) {
            const std::size_t base = r * cols_;
            for (std::size_t c = 0; c < cols_; ++c) {
                fn(pack_edge(base + c, 1));
            }
        }
    }

    template <class Fn>
    void scan_dual_vertices_for_apparent_edges(Fn&& on_apparent) {
        for (std::size_t x = 0; x < n_squares_; ++x) {
            const EdgeCode best = square_birth_edge_[x];
            const auto [a, b] = dual_endpoints(best);
            const std::size_t y = dual_younger_endpoint(a, b);
            if (y == x) {
                on_apparent(x, best);
            }
        }
    }

    bool is_apparent_dual_edge(const std::uint8_t* apparent_edge_mask,
                               EdgeCode e) const {
        return (apparent_edge_mask[edge_source(e)] &
                static_cast<std::uint8_t>(1u << edge_dir(e))) != 0;
    }

    void apply_apparent_edges_from_mask(const std::uint8_t* apparent_edge_mask,
                                        H1Result& result) {
        for (std::size_t source = 0; source < n_vertices_; ++source) {
            const std::uint8_t mask = apparent_edge_mask[source];
            if ((mask & 1u) != 0) {
                apply_apparent_dual_edge(pack_edge(source, 0), result);
            }
            if ((mask & 2u) != 0) {
                apply_apparent_dual_edge(pack_edge(source, 1), result);
            }
        }
    }

    static std::size_t popcount4(std::uint8_t x) {
        x &= 0x0fu;
        return static_cast<std::size_t>((x & 1u) + ((x >> 1) & 1u) +
                                        ((x >> 2) & 1u) + ((x >> 3) & 1u));
    }

    std::size_t count_mask_edges(const std::vector<VertexCode>& active_vertices,
                                 const std::uint8_t* edge_owner_mask) const {
        std::size_t count = 0;
        for (VertexCode v : active_vertices) {
            count += popcount4(edge_owner_mask[v]);
        }
        return count;
    }

    void scan_dual_apparent_pairs(std::uint8_t* apparent_edge_mask,
                                  H1Result& result) {
        auto mark_edge_in_mask = [&](std::uint8_t* mask, EdgeCode e) {
            mask[edge_source(e)] |=
                static_cast<std::uint8_t>(1u << edge_dir(e));
        };

        scan_dual_vertices_for_apparent_edges(
            [&](std::size_t x, EdgeCode best) {
                mark_edge_in_mask(apparent_edge_mask, best);
                edge_state_mark(best, edge_state_h1_apparent | edge_state_h1_positive);
                const auto [a, b] = dual_endpoints(best);
                const std::size_t other = (x == a) ? b : a;
                parent_[x] = static_cast<VertexCode>(other);
                ++result.uf_merge_attempts;
                ++result.uf_successful_merges;
                const double birth = edge_value(best);
                const double death = dual_birth(x);
                if (death > birth && !std::isinf(death)) {
                    result.pairs.push_back({birth, death});
                }
                ++result.apparent_edges;
            });
    }

    EdgeCode local_lookup_edge(std::size_t v, std::uint8_t local) const {
        switch (local) {
        case zero_lookup::edge2_n:
            return pack_edge(v - cols_, 1);
        case zero_lookup::edge2_e:
            return pack_edge(v, 0);
        case zero_lookup::edge2_s:
            return pack_edge(v, 1);
        case zero_lookup::edge2_w:
            return pack_edge(v - 1, 0);
        default:
            throw std::logic_error("invalid local 2D edge id");
        }
    }

    std::size_t local_lookup_square(std::size_t v, std::uint8_t local) const {
        const std::size_t r = row(v);
        const std::size_t c = col(v);
        switch (local) {
        case zero_lookup::square2_ne:
            return square_id(r - 1, c);
        case zero_lookup::square2_se:
            return square_id(r, c);
        case zero_lookup::square2_sw:
            return square_id(r, c - 1);
        case zero_lookup::square2_nw:
            return square_id(r - 1, c - 1);
        default:
            throw std::logic_error("invalid local 2D square id");
        }
    }

    void mark_applied_lookup_edge(std::uint8_t* applied_edge_mask,
                                  EdgeCode e) const {
        applied_edge_mask[edge_source(e)] |=
            static_cast<std::uint8_t>(1u << edge_dir(e));
    }

    void local_lookup_edge_parts(std::size_t v,
                                 std::uint8_t local,
                                 EdgeCode& e,
                                 std::size_t& source,
                                 std::uint8_t& dir) const {
        switch (local) {
        case zero_lookup::edge2_n:
            source = v - cols_;
            dir = 1;
            break;
        case zero_lookup::edge2_e:
            source = v;
            dir = 0;
            break;
        case zero_lookup::edge2_s:
            source = v;
            dir = 1;
            break;
        case zero_lookup::edge2_w:
            source = v - 1;
            dir = 0;
            break;
        default:
            throw std::logic_error("invalid local 2D edge id");
        }
        e = pack_edge(source, dir);
    }

    std::pair<std::size_t, std::size_t> local_dual_endpoints(
        std::size_t r,
        std::size_t c,
        std::uint8_t local) const {
        switch (local) {
        case zero_lookup::edge2_n: {
            const std::size_t er = r - 1;
            if (c == 0) {
                return {square_id(er, 0), infinity_};
            }
            if (c + 1 == cols_) {
                return {square_id(er, c - 1), infinity_};
            }
            return {square_id(er, c - 1), square_id(er, c)};
        }
        case zero_lookup::edge2_e:
            if (r == 0) {
                return {square_id(0, c), infinity_};
            }
            if (r + 1 == rows_) {
                return {square_id(r - 1, c), infinity_};
            }
            return {square_id(r - 1, c), square_id(r, c)};
        case zero_lookup::edge2_s:
            if (c == 0) {
                return {square_id(r, 0), infinity_};
            }
            if (c + 1 == cols_) {
                return {square_id(r, c - 1), infinity_};
            }
            return {square_id(r, c - 1), square_id(r, c)};
        case zero_lookup::edge2_w: {
            const std::size_t ec = c - 1;
            if (r == 0) {
                return {square_id(0, ec), infinity_};
            }
            if (r + 1 == rows_) {
                return {square_id(r - 1, ec), infinity_};
            }
            return {square_id(r - 1, ec), square_id(r, ec)};
        }
        default:
            throw std::logic_error("invalid local 2D edge id");
        }
    }

    void apply_lookup_positive_dual_edge(
        EdgeCode e,
        std::uint8_t* applied_edge_mask,
        H1Result& result,
        zero_lookup::RunStats* zero_lookup_stats) {
        if (applied_edge_mask != nullptr) {
            mark_applied_lookup_edge(applied_edge_mask, e);
        }
        edge_state_mark(e, edge_state_h1_positive);

        const auto [a, b] = dual_endpoints(e);
        std::size_t ra = find(a);
        std::size_t rb = find(b);
        ++result.uf_merge_attempts;
        if (ra == rb) {
            ++result.uf_same_component_attempts;
            return;
        }

        std::size_t older = ra;
        std::size_t younger = rb;
        if (!dual_older(ra, rb)) {
            older = rb;
            younger = ra;
        }

        const double birth = edge_value(e);
        const double death = dual_birth(younger);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
        ++result.apparent_edges;
        if (zero_lookup_stats != nullptr) {
            ++zero_lookup_stats->lookup_h1_positive_edges_applied;
        }
    }

    void apply_lookup_positive_dual_edge(
        std::size_t v,
        std::size_t r,
        std::size_t c,
        std::uint8_t local,
        std::uint8_t* applied_edge_mask,
        H1Result& result,
        zero_lookup::RunStats* zero_lookup_stats) {
        EdgeCode e = 0;
        std::size_t source = 0;
        std::uint8_t dir = 0;
        local_lookup_edge_parts(v, local, e, source, dir);
        if (applied_edge_mask != nullptr) {
            applied_edge_mask[source] |= static_cast<std::uint8_t>(1u << dir);
        }
        edge_state_mark(e, edge_state_h1_positive);

        const auto [a, b] = local_dual_endpoints(r, c, local);
        std::size_t ra = find(a);
        std::size_t rb = find(b);
        ++result.uf_merge_attempts;
        if (ra == rb) {
            ++result.uf_same_component_attempts;
            return;
        }

        std::size_t older = ra;
        std::size_t younger = rb;
        if (!dual_older(ra, rb)) {
            older = rb;
            younger = ra;
        }

        const double birth = values_[v];
        const double death = dual_birth(younger);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
        ++result.apparent_edges;
        if (zero_lookup_stats != nullptr) {
            ++zero_lookup_stats->lookup_h1_positive_edges_applied;
        }
    }

    void scan_dual_lookup_positive_edges(
        const std::uint8_t* zero_lookup_masks,
        std::uint8_t* applied_edge_mask,
        H1Result& result,
        zero_lookup::RunStats* zero_lookup_stats) {
        const auto apply_start = Clock::now();
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const std::size_t r = row(v);
            const std::size_t c = col(v);
            const zero_lookup::Entry2D& entry =
                zero_lookup::lookup2d(zero_lookup_masks[v]);
            if (zero_lookup_stats != nullptr) {
                zero_lookup_stats->lookup_h1_positive_edges +=
                    entry.edge_h1_positive_count;
                zero_lookup_stats->lookup_h1_negative_squares +=
                    entry.square_h1_negative_count;
                zero_lookup_stats->lookup_survivor_edges +=
                    entry.survivor_edge_count;
                zero_lookup_stats->lookup_survivor_squares +=
                    entry.survivor_square_count;
            }
            for (std::uint8_t local : entry.edge_order_desc) {
                if (local == 0xffu ||
                    (entry.edge_h1_positive_mask &
                     static_cast<std::uint8_t>(1u << local)) == 0) {
                    continue;
                }
                apply_lookup_positive_dual_edge(v, r, c, local,
                                                applied_edge_mask, result,
                                                zero_lookup_stats);
            }
            if (zero_lookup_stats != nullptr) {
                zero_lookup_stats->direct_h1_dual_pairs_found +=
                    entry.direct_h1_dual_count;
                zero_lookup_stats->direct_h1_dual_pairs_applied +=
                    entry.direct_h1_dual_count;
            }
        }
        if (zero_lookup_stats != nullptr) {
            zero_lookup_stats->table_apply_ms +=
                elapsed_ms(apply_start, Clock::now());
        }
    }

    void apply_apparent_dual_edge(EdgeCode e, H1Result& result) {
        const auto [a, b] = dual_endpoints(e);
        const std::size_t y = dual_younger_endpoint(a, b);
        const std::size_t other = (y == a) ? b : a;
        parent_[y] = static_cast<VertexCode>(other);
        ++result.uf_merge_attempts;
        ++result.uf_successful_merges;

        const double birth = edge_value(e);
        const double death = dual_birth(y);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
    }

    void fill_remaining_dual_edges(const std::uint8_t* apparent_edge_mask,
                                   std::uint8_t* edge_owner_mask,
                                   std::vector<VertexCode>& active_vertices,
                                   H1Result& result,
                                   bool apparent_edges_already_applied = false) {
        for_each_primal_edge([&](EdgeCode e) {
            if (is_apparent_dual_edge(apparent_edge_mask, e)) {
                if (!apparent_edges_already_applied) {
                    apply_apparent_dual_edge(e, result);
                }
                return;
            }
            if (edge_state_has(e, edge_state_h0_negative)) {
                ++result.skipped_known_h0_negative_edges;
                return;
            }
            mark_remaining_dual_edge(edge_owner_mask, active_vertices, e);
        });
    }

    std::pair<std::size_t, std::size_t> dual_endpoints(EdgeCode e) const {
        const std::size_t s = edge_source(e);
        const std::size_t r = row(s);
        const std::size_t c = col(s);

        if (edge_dir(e) == 0) {
            if (r == 0) {
                return {square_id(0, c), infinity_};
            }
            if (r + 1 == rows_) {
                return {square_id(r - 1, c), infinity_};
            }
            return {square_id(r - 1, c), square_id(r, c)};
        }

        if (c == 0) {
            return {square_id(r, 0), infinity_};
        }
        if (c + 1 == cols_) {
            return {square_id(r, c - 1), infinity_};
        }
        return {square_id(r, c - 1), square_id(r, c)};
    }

    void sweep_dual_edge(EdgeCode e, H1Result& result) {
        if (edge_state_has(e, edge_state_h0_negative)) {
            ++result.skipped_known_h0_negative_edges;
            return;
        }

        const auto [a, b] = dual_endpoints(e);
        std::size_t ra = find(a);
        std::size_t rb = find(b);
        ++result.uf_merge_attempts;
        if (ra == rb) {
            ++result.uf_same_component_attempts;
            return;
        }

        std::size_t older = ra;
        std::size_t younger = rb;
        if (!dual_older(ra, rb)) {
            older = rb;
            younger = ra;
        }

        const double birth = edge_value(e);
        const double death = dual_birth(younger);
        if (death > birth && !std::isinf(death)) {
            edge_state_mark(e, edge_state_h1_positive);
            result.pairs.push_back({birth, death});
        }
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
    }

    void sweep_dual_edge_lookup(std::size_t v,
                                std::size_t r,
                                std::size_t c,
                                std::uint8_t local,
                                H1Result& result) {
        EdgeCode e = 0;
        std::size_t source = 0;
        std::uint8_t dir = 0;
        local_lookup_edge_parts(v, local, e, source, dir);
        (void)source;
        (void)dir;

        const auto [a, b] = local_dual_endpoints(r, c, local);
        std::size_t ra = find(a);
        std::size_t rb = find(b);
        ++result.uf_merge_attempts;
        if (ra == rb) {
            ++result.uf_same_component_attempts;
            return;
        }

        std::size_t older = ra;
        std::size_t younger = rb;
        if (!dual_older(ra, rb)) {
            older = rb;
            younger = ra;
        }

        const double birth = values_[v];
        const double death = dual_birth(younger);
        if (death > birth && !std::isinf(death)) {
            edge_state_mark(e, edge_state_h1_positive);
            result.pairs.push_back({birth, death});
        }
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
    }

    void sweep_dual_edges(const std::vector<VertexCode>& active_vertices,
                          const std::uint8_t* edge_owner_mask,
                          H1Result& result) {
        for (VertexCode owner_code : active_vertices) {
            const std::size_t y = owner_code;
            const std::uint8_t mask = edge_owner_mask[y];
            if (zero_lookup_masks_ != nullptr) {
                const zero_lookup::Entry2D& entry =
                    zero_lookup::lookup2d(zero_lookup_masks_[y]);
                for (std::uint8_t i = 0;
                     i < entry.residual_edge_order_desc_count; ++i) {
                    const std::uint8_t local = entry.residual_edge_order_desc[i];
                    if ((mask & local_edge_owned_bit(local)) != 0) {
                        sweep_dual_edge(local_lookup_edge(y, local), result);
                    }
                }
                continue;
            }

            EdgeCode edges[4];
            int count = 0;

            if ((mask & edge_left_bit) != 0) {
                edges[count++] = pack_edge(y - 1, 0);
            }
            if ((mask & edge_right_bit) != 0) {
                edges[count++] = pack_edge(y, 0);
            }
            if ((mask & edge_up_bit) != 0) {
                edges[count++] = pack_edge(y - cols_, 1);
            }
            if ((mask & edge_down_bit) != 0) {
                edges[count++] = pack_edge(y, 1);
            }

            sort_edges_descending(edges, count);
            for (int i = 0; i < count; ++i) {
                sweep_dual_edge(edges[i], result);
            }
        }
    }

    void sweep_dual_edges_lookup_reverse(
        const std::vector<VertexCode>& active_vertices,
        const std::uint8_t* edge_owner_mask,
        H1Result& result) {
        for (std::size_t i = active_vertices.size(); i > 0; --i) {
            const std::size_t y = active_vertices[i - 1];
            const std::size_t r = row(y);
            const std::size_t c = col(y);
            const std::uint8_t mask = edge_owner_mask[y];
            const zero_lookup::Entry2D& entry =
                zero_lookup::lookup2d(zero_lookup_masks_[y]);
            for (std::uint8_t j = 0;
                 j < entry.residual_edge_order_desc_count; ++j) {
                const std::uint8_t local = entry.residual_edge_order_desc[j];
                if ((mask & local_edge_owned_bit(local)) != 0) {
                    sweep_dual_edge_lookup(y, r, c, local, result);
                }
            }
        }
    }
};

inline H1Result compute_h1_dual(const double* values, std::size_t rows,
                                std::size_t cols,
                                std::uint8_t* edge_state,
                                const std::uint32_t* value_codes = nullptr) {
    return H1DualComputer(values, rows, cols, edge_state, value_codes).run();
}

inline H1Result compute_h1_dual(const double* values, std::size_t rows,
                                std::size_t cols) {
    return compute_h1_dual(values, rows, cols, nullptr);
}

inline void mark_apparent_h1_positive_edges(const double* values,
                                            std::size_t rows,
                                            std::size_t cols,
                                            std::uint8_t* edge_state,
                                            const std::uint32_t* value_codes = nullptr) {
    H1DualComputer(values, rows, cols, edge_state, value_codes)
        .mark_apparent_positive_edges_only();
}

} // namespace smart_h0
