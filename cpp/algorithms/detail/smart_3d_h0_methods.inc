    std::uint8_t h0_owned_edge_bit(EdgeCode e, std::size_t owner) const {
        const std::size_t source = edge_source(e);
        const unsigned dir = edge_dir(e);
        if (dir == 0) {
            return (owner == source) ? edge_x_pos_bit : edge_x_neg_bit;
        }
        if (dir == 1) {
            return (owner == source) ? edge_y_pos_bit : edge_y_neg_bit;
        }
        return (owner == source) ? edge_z_pos_bit : edge_z_neg_bit;
    }

    void h0_mark_non_apparent(std::uint8_t* mask,
                              std::vector<VertexCode>& active_vertices,
                              std::size_t owner,
                              EdgeCode e,
                              H0Result& result) const {
        const std::uint8_t bit = h0_owned_edge_bit(e, owner);
        if (mask[owner] == 0) {
            active_vertices.push_back(static_cast<VertexCode>(owner));
        }
        mask[owner] |= bit;
        ++result.non_apparent_edges;
    }

    void h0_scan_apparent(std::uint8_t* non_apparent_mask,
                          std::vector<VertexCode>& active_vertices,
                          H0Result& result) {
        h0_scan_lookup(non_apparent_mask, active_vertices, result);
    }

    static std::size_t popcount8(std::uint8_t x) {
        std::size_t count = 0;
        while (x != 0) {
            count += x & 1u;
            x = static_cast<std::uint8_t>(x >> 1);
        }
        return count;
    }

    static std::size_t popcount16(std::uint16_t x) {
        std::size_t count = 0;
        while (x != 0) {
            count += x & 1u;
            x = static_cast<std::uint16_t>(x >> 1);
        }
        return count;
    }

    EdgeCode h0_local_lookup_edge(std::size_t v, std::uint8_t local) const {
        return lookup_local_edge_at_vertex(v, local);
    }

    void h0_scan_lookup(std::uint8_t* non_apparent_mask,
                        std::vector<VertexCode>& active_vertices,
                        H0Result& result) {
        const auto apply_start = Clock::now();
        for (std::size_t y = 0; y < n_vertices_; ++y) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[y]);
            const std::uint8_t direct_h0 = entry.direct_h0_edge;
            if (zero_lookup_stats_ != nullptr) {
                if (direct_h0 != 0xffu) {
                    ++zero_lookup_stats_->direct_h0_pairs_found;
                }
                zero_lookup_stats_->lookup_h1_positive_edges +=
                    popcount8(entry.edge_h1_positive_mask);
                zero_lookup_stats_->lookup_survivor_edges +=
                    popcount8(entry.survivor_edge_mask);
            }
            result.skipped_known_h1_positive_edges +=
                popcount8(entry.edge_h1_positive_mask);

            if (direct_h0 != 0xffu) {
                const EdgeCode apparent = h0_local_lookup_edge(y, direct_h0);
                const auto [u, v] = edge_endpoints(apparent);
                h0_parent_[y] = static_cast<VertexCode>((u == y) ? v : u);
                ++result.uf_merge_attempts;
                ++result.uf_successful_merges;
                edge_state_mark(apparent, edge_state_h0_apparent |
                                             edge_state_h0_negative);
                ++result.apparent_edges;
                if (zero_lookup_stats_ != nullptr) {
                    ++zero_lookup_stats_->direct_h0_pairs_applied;
                }
                if (include_apparent_zero_pairs_) {
                    result.pairs.push_back({values_[y], edge_birth_value(apparent)});
                }
            }

            for (std::uint8_t i = 0; i < entry.residual_edge_asc_count; ++i) {
                const EdgeCode e =
                    h0_local_lookup_edge(y, entry.residual_edge_order_asc[i]);
                if (edge_state_has(e, edge_state_h1_positive)) {
                    ++result.skipped_known_h1_positive_edges;
                    continue;
                }
                h0_mark_non_apparent(non_apparent_mask, active_vertices, y,
                                     e, result);
            }
        }
        if (zero_lookup_stats_ != nullptr) {
            zero_lookup_stats_->table_apply_ms +=
                elapsed_ms(apply_start, Clock::now());
        }
    }

    void h0_sweep_edge(EdgeCode e, H0Result& result) {
        if (edge_state_has(e, edge_state_h1_positive)) {
            ++result.skipped_known_h1_positive_edges;
            return;
        }

        const auto [u, v] = edge_endpoints(e);
        std::size_t ru = h0_parent_.find(static_cast<VertexCode>(u));
        std::size_t rv = h0_parent_.find(static_cast<VertexCode>(v));
        ++result.uf_merge_attempts;
        if (ru == rv) {
            ++result.uf_same_component_attempts;
            return;
        }

        std::size_t younger = ru;
        std::size_t older = rv;
        if (vertex_less(ru, rv)) {
            younger = rv;
            older = ru;
        }

        const double birth = values_[younger];
        const double death = edge_birth_value(e);
        if (include_apparent_zero_pairs_ || death > birth) {
            result.pairs.push_back({birth, death});
        }
        edge_state_mark(e, edge_state_h0_negative);
        h0_parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
    }

    void h0_sweep(const std::vector<VertexCode>& active_vertices,
                  const std::uint8_t* non_apparent_mask,
                  H0Result& result) {
        for (VertexCode owner_code : active_vertices) {
            const std::size_t y = owner_code;
            const std::uint8_t mask = non_apparent_mask[y];
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[y]);
            for (std::uint8_t i = 0;
                 i < entry.residual_edge_asc_count; ++i) {
                const std::uint8_t local =
                    entry.residual_edge_order_asc[i];
                const std::uint8_t bit =
                    static_cast<std::uint8_t>(1u << local);
                if ((mask & bit) != 0) {
                    h0_sweep_edge(h0_local_lookup_edge(y, local), result);
                }
            }
        }
    }

    H0Result run_h0() {
        H0Result result;
        const auto total_start = Clock::now();

        const auto init_start = Clock::now();
        h0_parent_.reset(n_vertices_);
        std::unique_ptr<std::uint8_t[]> non_apparent_mask(
            new std::uint8_t[n_vertices_]());
        std::vector<VertexCode> active_vertices;
        active_vertices.reserve(n_vertices_ / 4 + 1);
        result.pairs.reserve(n_vertices_ / 4 + 1);
        result.timing.parent_initialization_ms = elapsed_ms(init_start, Clock::now());

        const auto scan_start = Clock::now();
        h0_scan_apparent(non_apparent_mask.get(), active_vertices, result);
        result.timing.apparent_pair_scan_ms = elapsed_ms(scan_start, Clock::now());
        result.timing.apparent_pair_work_ms = result.timing.apparent_pair_scan_ms;

        const auto sort_start = Clock::now();
        std::sort(active_vertices.begin(), active_vertices.end(),
                  [this](VertexCode a, VertexCode b) {
                      return key_less(vertex_key(a), vertex_key(b));
                  });
        result.timing.sorting_ms = elapsed_ms(sort_start, Clock::now());

        result.active_non_apparent_vertices = active_vertices.size();
        const auto sweep_start = Clock::now();
        h0_sweep(active_vertices, non_apparent_mask.get(), result);
        result.timing.union_find_sweep_ms = elapsed_ms(sweep_start, Clock::now());
        result.timing.total_ms = elapsed_ms(total_start, Clock::now());

        result.memory.parent_bytes = n_vertices_ * sizeof(VertexCode);
        result.memory.non_apparent_vertex_mask_bytes =
            n_vertices_ * sizeof(std::uint8_t);
        result.memory.active_vertices_bytes =
            active_vertices.capacity() * sizeof(VertexCode);
        result.memory.persistence_pairs_bytes =
            result.pairs.capacity() * sizeof(PersistencePair);
        return result;
    }
