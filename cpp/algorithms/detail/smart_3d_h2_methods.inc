    std::uint8_t h2_owned_square_bit(SquareCode s, std::size_t cube) const {
        const std::size_t source = cube_source(cube);
        if (s == pack_square(source, 2)) {
            return face_x_neg_bit;
        }
        if (s == pack_square(source + 1, 2)) {
            return face_x_pos_bit;
        }
        if (s == pack_square(source, 1)) {
            return face_y_neg_bit;
        }
        if (s == pack_square(source + cols_, 1)) {
            return face_y_pos_bit;
        }
        if (s == pack_square(source, 0)) {
            return face_z_neg_bit;
        }
        if (s == pack_square(source + plane_, 0)) {
            return face_z_pos_bit;
        }
        throw std::logic_error("square is not a face of cube");
    }

    void h2_mark_non_apparent(std::uint8_t* mask,
                              std::vector<CubeCode>& active_cubes,
                              std::size_t cube,
                              SquareCode s,
                              H2Result& result) const {
        const std::uint8_t bit = h2_owned_square_bit(s, cube);
        if ((mask[cube] & bit) != 0) {
            return;
        }
        if (mask[cube] == 0) {
            active_cubes.push_back(static_cast<CubeCode>(cube));
        }
        mask[cube] |= bit;
        ++result.non_apparent_squares;
    }

    void h2_apply_apparent_square(std::size_t cube,
                                  SquareCode s,
                                  std::size_t other,
        H2Result& result) {
        dual_parent_[cube] = static_cast<CubeCode>(other);
        ++result.uf_merge_attempts;
        ++result.uf_successful_merges;
        square_state_mark(s, square_state_h2_apparent |
                             square_state_dual_h0_negative);
        ++result.apparent_squares;

        const double birth = square_birth_value(s);
        const double death = dual_birth_value(cube);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
    }

    void h2_scan_apparent(std::uint8_t* non_apparent_mask,
                          std::vector<CubeCode>& active_cubes,
                          std::vector<H1VertexOrderItem>* lookup_owner_vertices,
                          H2Result& result) {
        h2_scan_lookup(non_apparent_mask, active_cubes,
                       lookup_owner_vertices, result);
    }

    SquareCode h2_local_lookup_square(std::size_t v, std::uint8_t local) const {
        return lookup_local_square_at_vertex(v, local);
    }

    std::size_t h2_local_lookup_cube(std::size_t v, std::uint8_t local) const {
        return lookup_local_cube_at_vertex(v, local);
    }

    static std::uint16_t h2_local_cube_boundary(std::uint8_t local_cube) {
        const bool xp = (local_cube & 1u) != 0;
        const bool yp = (local_cube & 2u) != 0;
        const bool zp = (local_cube & 4u) != 0;
        const std::uint8_t xy =
            static_cast<std::uint8_t>((xp ? 1u : 0u) + 2u * (yp ? 1u : 0u));
        const std::uint8_t xz =
            static_cast<std::uint8_t>(4u + (xp ? 1u : 0u) +
                                      2u * (zp ? 1u : 0u));
        const std::uint8_t yz =
            static_cast<std::uint8_t>(8u + (yp ? 1u : 0u) +
                                      2u * (zp ? 1u : 0u));
        return static_cast<std::uint16_t>((1u << xy) |
                                          (1u << xz) |
                                          (1u << yz));
    }

    static std::uint8_t h2_local_cube_face_bit(std::uint8_t local_cube,
                                               std::uint8_t local_square) {
        const bool xp = (local_cube & 1u) != 0;
        const bool yp = (local_cube & 2u) != 0;
        const bool zp = (local_cube & 4u) != 0;
        if (local_square < 4) {
            return zp ? face_z_neg_bit : face_z_pos_bit;
        }
        if (local_square < 8) {
            return yp ? face_y_neg_bit : face_y_pos_bit;
        }
        return xp ? face_x_neg_bit : face_x_pos_bit;
    }

    void h2_apply_direct_lookup_pair(std::size_t paired_cube,
                                     SquareCode s,
                                     H2Result& result) {
        const auto [a, b] = dual_endpoints(s);
        if (paired_cube != a && paired_cube != b) {
            return;
        }
        const std::size_t other = (paired_cube == a) ? b : a;
        ++result.uf_merge_attempts;
        if (dual_parent_[paired_cube] == paired_cube &&
            dual_older(other, paired_cube)) {
            const double birth = square_birth_value(s);
            const double death = dual_birth_value(paired_cube);
            if (death > birth && !std::isinf(death)) {
                result.pairs.push_back({birth, death});
            }
            square_state_mark(s, square_state_h2_apparent |
                                 square_state_dual_h0_negative);
            dual_parent_[paired_cube] = static_cast<CubeCode>(other);
            ++result.uf_successful_merges;
            ++result.apparent_squares;
            if (zero_lookup_stats_ != nullptr) {
                ++zero_lookup_stats_->direct_h2_dual_pairs_applied;
            }
            return;
        }

        std::size_t ra = dual_parent_.find(static_cast<CubeCode>(a));
        std::size_t rb = dual_parent_.find(static_cast<CubeCode>(b));
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

        const double birth = square_birth_value(s);
        const double death = dual_birth_value(younger);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
        square_state_mark(s, square_state_h2_apparent |
                             square_state_dual_h0_negative);
        dual_parent_[younger] = static_cast<CubeCode>(older);
        ++result.uf_successful_merges;
        ++result.apparent_squares;
        if (zero_lookup_stats_ != nullptr) {
            ++zero_lookup_stats_->direct_h2_dual_pairs_applied;
        }
    }

    void h2_scan_lookup(std::uint8_t* non_apparent_mask,
                        std::vector<CubeCode>& active_cubes,
                        std::vector<H1VertexOrderItem>* lookup_owner_vertices,
                        H2Result& result) {
        const auto apply_start = Clock::now();
        if (lookup_owner_vertices != nullptr) {
            lookup_owner_vertices->clear();
            lookup_owner_vertices->reserve(n_vertices_ / 4 + 1);
        }
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            if (zero_lookup_stats_ != nullptr) {
                zero_lookup_stats_->lookup_h1_negative_squares +=
                    popcount16(entry.square_h1_negative_mask);
                zero_lookup_stats_->lookup_h2_negative_squares +=
                    popcount16(entry.square_h2_negative_mask);
                zero_lookup_stats_->lookup_survivor_squares +=
                    popcount16(entry.survivor_square_mask);
                zero_lookup_stats_->lookup_survivor_cubes +=
                    popcount8(entry.survivor_cube_mask);
            }
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local = entry.square_order_asc[i];
                if ((entry.square_h1_negative_mask &
                     static_cast<std::uint16_t>(1u << local)) == 0) {
                    continue;
                }
                square_state_mark(h2_local_lookup_square(v, local),
                                  square_state_dual_h1_positive);
            }
            for (std::uint8_t i = 0; i < entry.direct_h2_count; ++i) {
                if (zero_lookup_stats_ != nullptr) {
                    ++zero_lookup_stats_->direct_h2_dual_pairs_found;
                }
                const zero_lookup::Direct3DH2Pair pair = entry.direct_h2_pairs[i];
                const SquareCode s = h2_local_lookup_square(v, pair.square);
                const std::size_t cube = h2_local_lookup_cube(v, pair.cube);
                h2_apply_direct_lookup_pair(cube, s, result);
            }
        }

        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            bool has_lookup_work = false;
            for (std::uint8_t i = 0;
                 i < entry.h2_residual_square_asc_count; ++i) {
                const std::uint8_t local_square =
                    entry.h2_residual_square_order_asc[i];
                const SquareCode s = h2_local_lookup_square(v, local_square);
                if (square_state_has(s, square_state_dual_h0_negative)) {
                    continue;
                }
                if (square_state_has(s, square_state_dual_h1_positive)) {
                    ++result.skipped_known_dual_h1_positive_squares;
                    continue;
                }
                const auto [a, b] = dual_endpoints(s);
                const std::size_t owner = (a == infinity_) ? b : a;
                h2_mark_non_apparent(non_apparent_mask, active_cubes,
                                     owner, s, result);
                has_lookup_work = true;
            }
            if (has_lookup_work && lookup_owner_vertices != nullptr) {
                const std::uint32_t code = value_codes_ != nullptr
                    ? value_codes_[v]
                    : std::numeric_limits<std::uint32_t>::max();
                lookup_owner_vertices->push_back({values_[v], code,
                                                  static_cast<VertexCode>(v)});
            }
        }
        if (lookup_owner_vertices != nullptr) {
            result.active_lookup_owner_vertices = lookup_owner_vertices->size();
        }

        if (zero_lookup_stats_ != nullptr) {
            zero_lookup_stats_->table_apply_ms +=
                elapsed_ms(apply_start, Clock::now());
        }
    }

    void h2_sort_lookup_vertex_order(
        std::vector<H1VertexOrderItem>& vertex_order) const {
        std::sort(vertex_order.begin(), vertex_order.end(),
                  [this](const H1VertexOrderItem& a,
                         const H1VertexOrderItem& b) {
                      if (value_codes_ != nullptr) {
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
                      return a.vertex < b.vertex;
                  });
    }

    void h2_sweep_lookup_residual_squares(
        const std::vector<H1VertexOrderItem>& vertex_order,
        H2Result& result) {
        for (auto it = vertex_order.rbegin(); it != vertex_order.rend(); ++it) {
            const std::size_t v = it->vertex;
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            for (std::uint8_t i = entry.h2_residual_square_asc_count;
                 i > 0; --i) {
                const std::uint8_t local_square =
                    entry.h2_residual_square_order_asc[i - 1u];
                const SquareCode s = h2_local_lookup_square(v, local_square);
                if (square_state_has(s, square_state_dual_h0_negative) ||
                    square_state_has(s, square_state_dual_h1_positive)) {
                    continue;
                }
                h2_sweep_square(s, result);
            }
        }
    }

    void h2_sweep_square(SquareCode s, H2Result& result) {
        if (square_state_has(s, square_state_dual_h1_positive)) {
            ++result.skipped_known_dual_h1_positive_squares;
            return;
        }

        const auto [a, b] = dual_endpoints(s);
        std::size_t ra = dual_parent_.find(static_cast<CubeCode>(a));
        std::size_t rb = dual_parent_.find(static_cast<CubeCode>(b));
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

        const double birth = square_birth_value(s);
        const double death = dual_birth_value(younger);
        if (death > birth && !std::isinf(death)) {
            result.pairs.push_back({birth, death});
        }
        square_state_mark(s, square_state_dual_h0_negative);
        dual_parent_[younger] = static_cast<CubeCode>(older);
        ++result.uf_successful_merges;
    }

    H2Result run_h2() {
        H2Result result;
        result.primal_squares = primal_square_count();
        result.dual_vertices = n_cubes_ == 0 ? 0 : n_cubes_ + 1;
        if (n_cubes_ == 0) {
            return result;
        }

        const auto total_start = Clock::now();
        const auto init_start = Clock::now();
        ensure_cube_births();
        dual_parent_.reset(n_cubes_ + 1);
        std::unique_ptr<std::uint8_t[]> non_apparent_mask(
            new std::uint8_t[n_cubes_]());
        std::vector<CubeCode> active_cubes;
        active_cubes.reserve(n_cubes_ / 4 + 1);
        result.pairs.reserve(n_cubes_ / 4 + 1);
        result.timing.parent_initialization_ms = elapsed_ms(init_start, Clock::now());

        std::vector<H1VertexOrderItem> lookup_vertex_order;

        const auto scan_start = Clock::now();
        h2_scan_apparent(non_apparent_mask.get(), active_cubes,
                         &lookup_vertex_order,
                         result);
        result.timing.apparent_pair_scan_ms = elapsed_ms(scan_start, Clock::now());
        result.timing.apparent_pair_work_ms = result.timing.apparent_pair_scan_ms;
        result.active_non_apparent_cubes = active_cubes.size();

        const auto sort_start = Clock::now();
        h2_sort_lookup_vertex_order(lookup_vertex_order);
        result.timing.sorting_ms = elapsed_ms(sort_start, Clock::now());

        const auto sweep_start = Clock::now();
        h2_sweep_lookup_residual_squares(lookup_vertex_order, result);
        result.timing.union_find_sweep_ms = elapsed_ms(sweep_start, Clock::now());
        result.timing.total_ms = elapsed_ms(total_start, Clock::now());

        result.memory.parent_bytes = (n_cubes_ + 1) * sizeof(CubeCode);
        result.memory.non_apparent_vertex_mask_bytes =
            n_cubes_ * sizeof(std::uint8_t);
        result.memory.active_vertices_bytes =
            active_cubes.capacity() * sizeof(CubeCode) +
            lookup_vertex_order.capacity() * sizeof(H1VertexOrderItem);
        result.memory.square_birth_bytes =
            n_cubes_ * sizeof(SquareCode) +
            square_state_storage_size(depth_, rows_, cols_) * sizeof(EdgeCode) +
            edge_state_storage_size(depth_, rows_, cols_) * sizeof(VertexCode);
        result.memory.persistence_pairs_bytes =
            result.pairs.capacity() * sizeof(PersistencePair);
        return result;
    }
