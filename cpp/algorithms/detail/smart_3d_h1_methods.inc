    bool is_h1_positive_edge(EdgeCode e) const {
        return !edge_state_has(e, edge_state_h0_negative);
    }

    bool is_h1_negative_square(SquareCode s) const {
        return !square_state_has(s, square_state_dual_h0_negative);
    }

    static bool same_edge(EdgeCode a, EdgeCode b) {
        return a == b;
    }

    H1Term h1_term(EdgeCode e) const {
        return static_cast<H1Term>(e);
    }

    EdgeCode h1_edge(H1Term e) const {
        return static_cast<EdgeCode>(e);
    }

    bool h1_term_less(H1Term a, H1Term b) const {
        return edge_less(h1_edge(a), h1_edge(b));
    }

    void sort_h1_terms_ascending(std::vector<H1Term>& terms) const {
        for (std::size_t i = 1; i < terms.size(); ++i) {
            const H1Term x = terms[i];
            std::size_t j = i;
            while (j > 0 && h1_term_less(x, terms[j - 1])) {
                terms[j] = terms[j - 1];
                --j;
            }
            terms[j] = x;
        }
    }

    static void sort_squares_ascending(SquareCode* squares, int count) {
        for (int i = 1; i < count; ++i) {
            const SquareCode x = squares[i];
            int j = i - 1;
            while (j >= 0 && squares[j] > x) {
                squares[j + 1] = squares[j];
                --j;
            }
            squares[j + 1] = x;
        }
    }

    void xor_with_sorted_terms(std::vector<H1Term>& dst,
                               const H1Term* src,
                               std::size_t src_size,
                               std::vector<H1Term>& scratch) const {
        scratch.clear();
        const std::size_t needed = dst.size() + src_size;
        if (scratch.capacity() < needed) {
            scratch.reserve(needed);
        }

        std::size_t i = 0;
        std::size_t j = 0;
        while (i < dst.size() && j < src_size) {
            if (dst[i] == src[j]) {
                ++i;
                ++j;
            } else if (h1_term_less(dst[i], src[j])) {
                scratch.push_back(dst[i++]);
            } else {
                scratch.push_back(src[j++]);
            }
        }
        while (i < dst.size()) {
            scratch.push_back(dst[i++]);
        }
        while (j < src_size) {
            scratch.push_back(src[j++]);
        }
        dst.swap(scratch);
    }

    void h1_initial_square_column(SquareCode s,
                                  std::vector<H1Term>& column) const {
        column.clear();
        EdgeCode edges[4];
        square_boundary_edges(s, edges);
        for (EdgeCode e : edges) {
            if (is_h1_positive_edge(e)) {
                column.push_back(h1_term(e));
            }
        }
        sort_h1_terms_ascending(column);
    }

    std::size_t h1_initial_apparent_square_rewrite(
        SquareCode s,
        EdgeCode pivot_edge,
        std::vector<H1Term>& column) const {
        column.clear();
        EdgeCode edges[4];
        square_boundary_edges(s, edges);
        bool found_pivot = false;
        std::size_t positive_terms = 0;
        const H1Term pivot = h1_term(pivot_edge);
        for (EdgeCode e : edges) {
            if (same_edge(e, pivot_edge)) {
                found_pivot = true;
                ++positive_terms;
                continue;
            }
            if (!is_h1_positive_edge(e)) {
                continue;
            }
            const H1Term term = h1_term(e);
            if (!h1_term_less(term, pivot)) {
                throw std::logic_error("apparent H1 pivot is not the youngest boundary edge");
            }
            column.push_back(term);
            ++positive_terms;
        }
        if (!found_pivot || !is_h1_positive_edge(pivot_edge)) {
            throw std::logic_error("apparent H1 square does not expose its pivot edge");
        }
        sort_h1_terms_ascending(column);
        return positive_terms;
    }

    // -----------------------------------------------------------------------
    // Combined pivot map encoding (indexed by H1Term / EdgeCode).
    //
    //   bits 31-30 = 00, val < 0x40000000:
    //                    stored rewrite row at rewrite_meta index val
    //   bits 31-30 = 01, pm_in_prog_flag <= val < pm_zero_sq_flag:
    //                    materialisation in progress for sq_index = val & pm_sq_mask
    //   bits 31-30 = 10, pm_zero_sq_flag <= val < pm_no_info:
    //                    unmaterialised zero square at sq_index = val & pm_sq_mask
    //   0xFFFFFFFF (pm_no_info): not a zero pair, no stored row → genuine new pivot
    //
    // The in-progress state encodes dep_sq so that the hot loop can still
    // check dep_sq == current_square for the own-pivot termination case.
    //
    // H1RewriteRow packs offset + length for a single cache-line load instead
    // of two separate SoA arrays.
    // -----------------------------------------------------------------------
    static constexpr std::uint32_t pm_no_info      = 0xFFFFFFFFu;
    static constexpr std::uint32_t pm_zero_sq_flag = 0x80000000u;
    static constexpr std::uint32_t pm_in_prog_flag = 0x40000000u;
    static constexpr std::uint32_t pm_sq_mask      = 0x3FFFFFFFu;

    struct H1RewriteRow {
        std::uint32_t offset;
        std::uint32_t length;
    };

    // Non-demand column reduction — only uses pre-existing rows in pivot_map.
    void h1_reduce_column(
        std::vector<H1Term>& column,
        const std::vector<std::uint32_t>& pivot_map,
        const std::vector<H1RewriteRow>& rewrite_meta,
        const std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& scratch,
        H1Result& result) const {
        while (!column.empty()) {
            const H1Term pivot = column.back();
            const std::uint32_t pm =
                pivot_map[static_cast<std::size_t>(pivot)];
            if (profile_h1_) {
                ++result.profile.reduction_pivot_checks;
            }
            if (pm >= pm_in_prog_flag) {
                break; // no stored row (zero sq, in-progress, or pm_no_info)
            }
            column.pop_back();
            const H1RewriteRow& row = rewrite_meta[pm];
            const std::size_t col_before = column.size();
            if (row.length > 0) {
                xor_with_sorted_terms(column,
                    rewrite_terms.data() + row.offset,
                    row.length, scratch);
            }
            if (profile_h1_) {
                H1Profile& profile = result.profile;
                ++profile.reduction_rewrite_hits;
                profile.reduction_rewrite_terms_read += row.length;
                profile.reduction_column_terms_before_xor += col_before;
                profile.reduction_output_terms_after_xor += column.size();
                profile.reduction_estimated_merge_input_terms +=
                    col_before + row.length;
                profile.max_rewrite_terms_read =
                    std::max(profile.max_rewrite_terms_read,
                             static_cast<std::size_t>(row.length));
                profile.max_column_terms_before_xor =
                    std::max(profile.max_column_terms_before_xor, col_before);
            }
            ++result.reduction_steps;
        }
        result.max_reduced_boundary_size =
            std::max(result.max_reduced_boundary_size, column.size());
    }

    void h1_store_rewrite(
        H1Term pivot,
        const std::vector<H1Term>& column,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        H1Result& result) const {
        if (rewrite_meta.size() >= static_cast<std::size_t>(pm_in_prog_flag)) {
            throw std::overflow_error("too many H1 rewrite rows for 30-bit index");
        }
        if (column.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error("H1 rewrite row is too large");
        }
        if (rewrite_terms.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error("H1 rewrite storage exceeds uint32 offsets");
        }
        const std::uint32_t row_index =
            static_cast<std::uint32_t>(rewrite_meta.size());
        pivot_map[static_cast<std::size_t>(pivot)] = row_index;
        rewrite_meta.push_back({
            static_cast<std::uint32_t>(rewrite_terms.size()),
            static_cast<std::uint32_t>(column.size())
        });
        rewrite_terms.insert(rewrite_terms.end(), column.begin(), column.end());
        if (profile_h1_) {
            H1Profile& profile = result.profile;
            ++profile.rewrite_rows;
            profile.rewrite_terms_stored += column.size();
            profile.max_rewrite_terms_stored =
                std::max(profile.max_rewrite_terms_stored, column.size());
        }
    }

    template <class Fn>
    void for_each_h1_negative_square_in_order(
        const std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result,
        Fn&& fn) const {
        for (const H1VertexOrderItem& item : vertex_order) {
            if (profile_h1_) {
                ++result.profile.stream_vertices;
            }
            const std::size_t v = item.vertex;
            std::size_t z = 0;
            std::size_t y = 0;
            std::size_t x = 0;
            vertex_coords(v, z, y, x);

            EdgeCode edges[6];
            int edge_count = 0;
            auto consider = [&](std::size_t other) {
                if (vertex_less(other, v)) {
                    edges[edge_count++] = edge_between(v, other);
                }
            };

            if (x > 0) {
                consider(v - 1);
            }
            if (x + 1 < cols_) {
                consider(v + 1);
            }
            if (y > 0) {
                consider(v - cols_);
            }
            if (y + 1 < rows_) {
                consider(v + cols_);
            }
            if (z > 0) {
                consider(v - plane_);
            }
            if (z + 1 < depth_) {
                consider(v + plane_);
            }

            sort_edges_ascending(edges, edge_count);
            if (profile_h1_ && edge_count != 0) {
                ++result.profile.stream_vertices_with_lower_edges;
                result.profile.stream_lower_edges +=
                    static_cast<std::size_t>(edge_count);
            }
            for (int i = 0; i < edge_count; ++i) {
                const EdgeCode e = edges[i];
                SquareCode cofacets[4];
                const int cofacet_count = incident_squares_of_edge(e, cofacets);
                if (profile_h1_) {
                    result.profile.stream_incident_square_refs +=
                        static_cast<std::size_t>(cofacet_count);
                }
                SquareCode born_here[4];
                int born_count = 0;
                for (int j = 0; j < cofacet_count; ++j) {
                    const SquareCode s = cofacets[j];
                    if (profile_h1_) {
                        ++result.profile.stream_square_birth_tests;
                    }
                    if (square_birth_edge(s) == e && is_h1_negative_square(s)) {
                        if (profile_h1_) {
                            ++result.profile.stream_square_birth_matches;
                            ++result.profile.stream_h1_negative_squares;
                        }
                        born_here[born_count++] = s;
                    }
                }
                sort_squares_ascending(born_here, born_count);
                if (profile_h1_ && born_count != 0) {
                    ++result.profile.stream_born_square_sort_groups;
                    result.profile.stream_born_square_sort_items +=
                        static_cast<std::size_t>(born_count);
                }
                for (int j = 0; j < born_count; ++j) {
                    fn(born_here[j]);
                }
            }
        }
    }

    template <class Fn>
    void for_each_h1_negative_square_in_lookup_order(
        const std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result,
        Fn&& fn) const {
        for (const H1VertexOrderItem& item : vertex_order) {
            if (profile_h1_) {
                ++result.profile.stream_vertices;
            }
            const std::size_t v = item.vertex;
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            if (profile_h1_ && entry.square_count != 0) {
                ++result.profile.stream_vertices_with_lower_edges;
                result.profile.stream_lower_edges += entry.edge_count;
            }
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                if (!is_h1_negative_square(s)) {
                    continue;
                }
                if (profile_h1_) {
                    ++result.profile.stream_square_birth_matches;
                    ++result.profile.stream_h1_negative_squares;
                }
                fn(s);
            }
        }
    }

    struct H1LookupOwnerContext {
        std::uint8_t h1_positive_local_edges = 0;
        std::uint16_t zero_square_mask = 0;
    };

    H1LookupOwnerContext h1_lookup_owner_context(
        std::size_t owner_vertex,
        const zero_lookup::Entry3D& entry) const {
        H1LookupOwnerContext context;
        for (std::uint8_t local_edge = 0; local_edge < 6; ++local_edge) {
            if ((entry.edge_mask & (1u << local_edge)) == 0) {
                continue;
            }
            const EdgeCode e =
                lookup_local_edge_at_vertex(owner_vertex, local_edge);
            if (is_h1_positive_edge(e)) {
                context.h1_positive_local_edges =
                    static_cast<std::uint8_t>(
                        context.h1_positive_local_edges |
                        (1u << local_edge));
            }
        }
        for (std::uint8_t i = 0; i < entry.zero_h1_count; ++i) {
            const zero_lookup::Direct3DH1Pair pair = entry.zero_h1_pairs[i];
            const EdgeCode e =
                lookup_local_edge_at_vertex(owner_vertex, pair.edge);
            const SquareCode s =
                lookup_local_square_at_vertex(owner_vertex, pair.square);
            if (is_h1_positive_edge(e) && is_h1_negative_square(s)) {
                context.zero_square_mask =
                    static_cast<std::uint16_t>(
                        context.zero_square_mask | (1u << pair.square));
            }
        }
        return context;
    }

    std::uint16_t h1_lookup_owner_zero_square_mask(
        std::size_t owner_vertex,
        const zero_lookup::Entry3D& entry) const {
        std::uint16_t mask = 0;
        for (std::uint8_t i = 0; i < entry.zero_h1_count; ++i) {
            const zero_lookup::Direct3DH1Pair pair = entry.zero_h1_pairs[i];
            const EdgeCode e =
                lookup_local_edge_at_vertex(owner_vertex, pair.edge);
            const SquareCode s =
                lookup_local_square_at_vertex(owner_vertex, pair.square);
            if (is_h1_positive_edge(e) && is_h1_negative_square(s)) {
                mask = static_cast<std::uint16_t>(mask | (1u << pair.square));
            }
        }
        return mask;
    }

    template <class Fn>
    void for_each_h1_lookup_square_in_order(
        const std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result,
        Fn&& fn) const {
        for (const H1VertexOrderItem& item : vertex_order) {
            if (profile_h1_) {
                ++result.profile.stream_vertices;
            }
            const std::size_t v = item.vertex;
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            if (profile_h1_ && entry.square_count != 0) {
                ++result.profile.stream_vertices_with_lower_edges;
                result.profile.stream_lower_edges += entry.edge_count;
            }
            const H1LookupOwnerContext context =
                h1_lookup_owner_context(v, entry);
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                if (!is_h1_negative_square(s)) {
                    continue;
                }
                if (profile_h1_) {
                    ++result.profile.stream_square_birth_matches;
                    ++result.profile.stream_h1_negative_squares;
                }
                fn(v, local_square, s, entry, context);
            }
        }
    }

    template <class Fn>
    void for_each_h1_lookup_nonzero_square_raw_order(
        const std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result,
        Fn&& fn) const {
        for (const H1VertexOrderItem& item : vertex_order) {
            if (profile_h1_) {
                ++result.profile.stream_vertices;
            }
            const std::size_t v = item.vertex;
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            const std::uint16_t zero_square_mask =
                h1_lookup_owner_zero_square_mask(v, entry);
            if (profile_h1_ && entry.square_count != 0) {
                ++result.profile.stream_vertices_with_lower_edges;
                result.profile.stream_lower_edges += entry.edge_count;
            }
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                if (!is_h1_negative_square(s)) {
                    continue;
                }
                if (profile_h1_) {
                    ++result.profile.stream_square_birth_matches;
                    ++result.profile.stream_h1_negative_squares;
                }
                if ((zero_square_mask & (1u << local_square)) != 0) {
                    continue;
                }
                fn(s);
            }
        }
    }

    void h1_collect_lookup_owner_vertices_eager(
        std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result) const {
        vertex_order.clear();
        vertex_order.reserve(n_vertices_ / 4 + 1);
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            bool has_work = false;
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                if (is_h1_negative_square(s)) {
                    has_work = true;
                    break;
                }
            }
            if (!has_work) {
                continue;
            }
            const std::uint32_t code = value_codes_ != nullptr
                ? value_codes_[v]
                : std::numeric_limits<std::uint32_t>::max();
            vertex_order.push_back({values_[v], code,
                                    static_cast<VertexCode>(v)});
        }
        result.active_owner_vertices = vertex_order.size();
        result.h1_zero_squares_skipped_outer = 0;
    }

    void h1_process_negative_square(
        SquareCode s,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& column,
        std::vector<H1Term>& scratch,
        H1Result& result) const {
        ++result.negative_squares;
        const bool apparent =
            square_state_has(s, square_state_dual_h1_apparent);
        if (apparent) {
            ++result.apparent_pairs;
        }

        H1Term pivot = std::numeric_limits<H1Term>::max();
        bool pivot_removed_from_column = false;
        std::size_t initial_positive_terms = 0;
        if (apparent) {
            const EdgeCode pivot_edge = square_birth_edge(s);
            pivot = h1_term(pivot_edge);
            initial_positive_terms =
                h1_initial_apparent_square_rewrite(s, pivot_edge, column);
            pivot_removed_from_column = true;
        } else {
            h1_initial_square_column(s, column);
            initial_positive_terms = column.size();
        }
        if (profile_h1_) {
            H1Profile& profile = result.profile;
            ++profile.initial_columns;
            profile.initial_boundary_edge_tests += 4;
            profile.initial_positive_terms += initial_positive_terms;
            const std::size_t hist_index = std::min<std::size_t>(
                initial_positive_terms,
                profile.initial_column_size_hist.size() - 1);
            ++profile.initial_column_size_hist[hist_index];
        }
        if (apparent) {
            if (profile_h1_) {
                ++result.profile.apparent_fast_path_columns;
            }
        } else {
            if (profile_h1_) {
                ++result.profile.nonapparent_reduced_columns;
            }
            h1_reduce_column(column, pivot_map, rewrite_meta,
                             rewrite_terms, scratch, result);
            if (!column.empty()) {
                pivot = column.back();
            }
        }
        if (pivot == std::numeric_limits<H1Term>::max()) {
            throw std::logic_error("H1-negative square reduced to zero");
        }

        if (pivot_map[static_cast<std::size_t>(pivot)] < pm_in_prog_flag) {
            throw std::logic_error("H1 reduction selected an already paired edge");
        }

        const EdgeCode pivot_edge = h1_edge(pivot);
        ++result.matched_pairs;
        if (apparent && !include_apparent_zero_pairs_) {
            ++result.zero_persistence_pairs;
        } else {
            const double birth = edge_birth_value(pivot_edge);
            const double death = apparent ? birth : square_birth_value(s);
            if (death > birth || include_apparent_zero_pairs_) {
                result.pairs.push_back({birth, death});
            } else {
                ++result.zero_persistence_pairs;
            }
        }

        if (!pivot_removed_from_column) {
            column.pop_back();
        }
        h1_store_rewrite(pivot, column, pivot_map, rewrite_meta,
                         rewrite_terms, result);
    }

    static std::uint8_t h1_local_square_boundary_edges(
        std::uint8_t local_square) {
        if (local_square < 4) {
            const bool xp = (local_square & 1u) != 0;
            const bool yp = (local_square & 2u) != 0;
            return static_cast<std::uint8_t>(
                (1u << (xp ? zero_lookup::edge3_x_pos
                           : zero_lookup::edge3_x_neg)) |
                (1u << (yp ? zero_lookup::edge3_y_pos
                           : zero_lookup::edge3_y_neg)));
        }
        if (local_square < 8) {
            const std::uint8_t t =
                static_cast<std::uint8_t>(local_square - 4u);
            const bool xp = (t & 1u) != 0;
            const bool zp = (t & 2u) != 0;
            return static_cast<std::uint8_t>(
                (1u << (xp ? zero_lookup::edge3_x_pos
                           : zero_lookup::edge3_x_neg)) |
                (1u << (zp ? zero_lookup::edge3_z_pos
                           : zero_lookup::edge3_z_neg)));
        }
        const std::uint8_t t =
            static_cast<std::uint8_t>(local_square - 8u);
        const bool yp = (t & 1u) != 0;
        const bool zp = (t & 2u) != 0;
        return static_cast<std::uint8_t>(
            (1u << (yp ? zero_lookup::edge3_y_pos
                       : zero_lookup::edge3_y_neg)) |
            (1u << (zp ? zero_lookup::edge3_z_pos
                       : zero_lookup::edge3_z_neg)));
    }

    bool h1_lookup_zero_pair_for_square(
        std::size_t owner_vertex,
        const zero_lookup::Entry3D& entry,
        std::uint8_t local_square,
        EdgeCode& pivot_edge) const {
        for (std::uint8_t i = 0; i < entry.zero_h1_count; ++i) {
            const zero_lookup::Direct3DH1Pair pair = entry.zero_h1_pairs[i];
            if (pair.square != local_square) {
                continue;
            }
            const EdgeCode global_pivot =
                lookup_local_edge_at_vertex(owner_vertex, pair.edge);
            const SquareCode global_square =
                lookup_local_square_at_vertex(owner_vertex, pair.square);
            if (!is_h1_positive_edge(global_pivot) ||
                !is_h1_negative_square(global_square)) {
                return false;
            }
            pivot_edge = global_pivot;
            return true;
        }
        return false;
    }

    void h1_toggle_positive_edge(std::vector<H1Term>& column,
                                 EdgeCode e) const {
        if (!is_h1_positive_edge(e)) {
            return;
        }
        const H1Term term = h1_term(e);
        for (auto it = column.begin(); it != column.end(); ++it) {
            if (*it == term) {
                column.erase(it);
                return;
            }
        }
        column.push_back(term);
    }

    void h1_toggle_lookup_nonlocal_boundary(
        std::size_t owner_vertex,
        std::uint8_t local_square,
        std::vector<H1Term>& column) const {
        const std::uint8_t local_boundary =
            h1_local_square_boundary_edges(local_square);
        const SquareCode s =
            lookup_local_square_at_vertex(owner_vertex, local_square);
        EdgeCode boundary[4];
        square_boundary_edges(s, boundary);
        for (EdgeCode e : boundary) {
            bool is_local_boundary_edge = false;
            for (std::uint8_t local_edge = 0; local_edge < 6; ++local_edge) {
                if ((local_boundary & (1u << local_edge)) == 0) {
                    continue;
                }
                if (e == lookup_local_edge_at_vertex(owner_vertex,
                                                     local_edge)) {
                    is_local_boundary_edge = true;
                    break;
                }
            }
            if (!is_local_boundary_edge) {
                h1_toggle_positive_edge(column, e);
            }
        }
    }

    bool h1_lookup_semireduced_template_safe(
        std::uint8_t local_square,
        const zero_lookup::Entry3D& entry,
        const H1LookupOwnerContext& context) const {
        const std::uint16_t substitutions =
            entry.h1_square_substitution_square_mask[local_square];
        if ((substitutions & ~context.zero_square_mask) != 0) {
            return false;
        }
        const std::uint16_t squares =
            static_cast<std::uint16_t>(substitutions | (1u << local_square));
        std::uint8_t required_local_edges = 0;
        for (std::uint8_t sq = 0; sq < 12; ++sq) {
            if ((squares & (1u << sq)) == 0) {
                continue;
            }
            required_local_edges =
                static_cast<std::uint8_t>(
                    required_local_edges |
                    h1_local_square_boundary_edges(sq) & entry.edge_mask &
                    ~entry.edge_h0_negative_mask);
        }
        return (required_local_edges &
                ~context.h1_positive_local_edges) == 0;
    }

    bool h1_initial_lookup_semireduced_column(
        std::size_t owner_vertex,
        std::uint8_t local_square,
        SquareCode s,
        const zero_lookup::Entry3D& entry,
        const H1LookupOwnerContext& context,
        std::vector<H1Term>& column) const {
        if (!h1_lookup_semireduced_template_safe(local_square, entry,
                                                context)) {
            h1_initial_square_column(s, column);
            return false;
        }

        column.clear();
        const std::uint8_t local_edge_mask =
            entry.h1_square_semireduced_boundary[local_square];
        for (std::uint8_t local_edge = 0; local_edge < 6; ++local_edge) {
            if ((local_edge_mask & (1u << local_edge)) == 0) {
                continue;
            }
            const EdgeCode e =
                lookup_local_edge_at_vertex(owner_vertex, local_edge);
            h1_toggle_positive_edge(column, e);
        }

        h1_toggle_lookup_nonlocal_boundary(owner_vertex, local_square, column);
        const std::uint16_t substitutions =
            entry.h1_square_substitution_square_mask[local_square];
        for (std::uint8_t dep_square = 0; dep_square < 12; ++dep_square) {
            if ((substitutions & (1u << dep_square)) != 0) {
                h1_toggle_lookup_nonlocal_boundary(owner_vertex, dep_square,
                                                   column);
            }
        }
        sort_h1_terms_ascending(column);
        return true;
    }

    void h1_record_lookup_initial_profile(
        const std::vector<H1Term>& column,
        H1Result& result) const {
        if (!profile_h1_) {
            return;
        }
        H1Profile& p = result.profile;
        ++p.initial_columns;
        p.initial_boundary_edge_tests += 4;
        p.initial_positive_terms += column.size();
        const std::size_t hist_idx = std::min<std::size_t>(
            column.size(), p.initial_column_size_hist.size() - 1);
        ++p.initial_column_size_hist[hist_idx];
        ++p.nonapparent_reduced_columns;
    }

    void h1_install_lookup_zero_rewrite_eager(
        std::size_t owner_vertex,
        std::uint8_t local_square,
        EdgeCode expected_pivot,
        const zero_lookup::Entry3D& entry,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& column,
        std::vector<H1Term>& scratch,
        H1Result& result) const {
        const auto build_start = profile_h1_ ? Clock::now()
                                             : Clock::time_point{};
        const H1LookupOwnerContext context =
            h1_lookup_owner_context(owner_vertex, entry);
        const bool used_template = h1_initial_lookup_semireduced_column(
            owner_vertex, local_square,
            lookup_local_square_at_vertex(owner_vertex, local_square),
            entry, context, column);
        if (profile_h1_) {
            result.h1_lookup_zero_build_ms +=
                elapsed_ms(build_start, Clock::now());
        }
        if (used_template) {
            ++result.h1_lookup_zero_template_columns;
        } else {
            ++result.h1_lookup_zero_raw_fallback_columns;
        }
        h1_record_lookup_initial_profile(column, result);
        const auto reduce_start = profile_h1_ ? Clock::now()
                                              : Clock::time_point{};
        h1_reduce_column(column, pivot_map, rewrite_meta,
                         rewrite_terms, scratch, result);
        if (profile_h1_) {
            result.h1_lookup_zero_reduce_ms +=
                elapsed_ms(reduce_start, Clock::now());
        }
        if (column.empty()) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "eager lookup H1 zero square reduced to empty column");
        }
        const H1Term pivot = column.back();
        const EdgeCode actual_pivot = h1_edge(pivot);
        if (actual_pivot != expected_pivot) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "eager lookup H1 zero square reduced to unexpected pivot");
        }
        if (pivot_map[static_cast<std::size_t>(pivot)] < pm_in_prog_flag) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "eager lookup H1 zero row selected an already paired edge");
        }
        column.pop_back();
        h1_store_rewrite(pivot, column, pivot_map, rewrite_meta,
                         rewrite_terms, result);
        ++result.h1_zero_rewrites_materialised;
    }

    void h1_process_lookup_square_eager(
        std::size_t owner_vertex,
        std::uint8_t local_square,
        SquareCode s,
        const zero_lookup::Entry3D& entry,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& column,
        std::vector<H1Term>& scratch,
        H1Result& result) const {
        EdgeCode expected_zero_pivot =
            std::numeric_limits<EdgeCode>::max();
        if (h1_lookup_zero_pair_for_square(owner_vertex, entry, local_square,
                                           expected_zero_pivot)) {
            h1_install_lookup_zero_rewrite_eager(
                owner_vertex, local_square, expected_zero_pivot, entry,
                pivot_map, rewrite_meta, rewrite_terms,
                column, scratch, result);
            return;
        }

        ++result.negative_squares;
        const auto build_start = profile_h1_ ? Clock::now()
                                             : Clock::time_point{};
        const H1LookupOwnerContext context =
            h1_lookup_owner_context(owner_vertex, entry);
        const bool used_template = h1_initial_lookup_semireduced_column(
            owner_vertex, local_square, s, entry, context, column);
        if (profile_h1_) {
            result.h1_lookup_nonzero_build_ms +=
                elapsed_ms(build_start, Clock::now());
        }
        if (used_template) {
            ++result.h1_lookup_nonzero_template_columns;
        } else {
            ++result.h1_lookup_nonzero_raw_fallback_columns;
        }
        h1_record_lookup_initial_profile(column, result);
        const auto reduce_start = profile_h1_ ? Clock::now()
                                              : Clock::time_point{};
        h1_reduce_column(column, pivot_map, rewrite_meta,
                         rewrite_terms, scratch, result);
        if (profile_h1_) {
            result.h1_lookup_nonzero_reduce_ms +=
                elapsed_ms(reduce_start, Clock::now());
        }
        if (column.empty()) {
            throw std::logic_error(
                "H1-negative square reduced to zero in eager lookup path: owner " +
                std::to_string(owner_vertex) + " local square " +
                std::to_string(local_square) + " global square " +
                std::to_string(s));
        }
        const H1Term pivot = column.back();
        if (pivot_map[static_cast<std::size_t>(pivot)] < pm_in_prog_flag) {
            throw std::logic_error(
                "H1 eager lookup reduction selected an already paired edge");
        }
        const EdgeCode pivot_edge = h1_edge(pivot);
        ++result.matched_pairs;
        const double birth = edge_birth_value(pivot_edge);
        const double death = square_birth_value(s);
        if (death > birth || include_apparent_zero_pairs_) {
            result.pairs.push_back({birth, death});
        } else {
            ++result.zero_persistence_pairs;
        }
        column.pop_back();
        h1_store_rewrite(pivot, column, pivot_map, rewrite_meta,
                         rewrite_terms, result);
    }

    // -----------------------------------------------------------------------
    // Miscellaneous sentinels kept for non-production code paths.
    // -----------------------------------------------------------------------
    static constexpr bool h1_lookup_use_template_zero_rows = false;
    static constexpr bool h1_lookup_use_template_nonzero_columns = false;

    struct H1MaterialiseBuffers {
        std::deque<std::vector<H1Term>> columns;
        std::deque<std::vector<H1Term>> scratches;
        std::size_t depth = 0;
    };

    // -----------------------------------------------------------------------
    // Build pivot_map from the zero-lookup table.  For each confirmed zero-pair
    // (edge e, square s), encodes pivot_map[e] = pm_zero_sq_flag | s so that
    // the hot-loop recognises it as an unmaterialised zero square in one load.
    // -----------------------------------------------------------------------
    void h1_build_zero_pivot_map(
            std::vector<std::uint32_t>& pivot_map) const {
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            for (std::uint8_t i = 0; i < entry.zero_h1_count; ++i) {
                const zero_lookup::Direct3DH1Pair pair = entry.zero_h1_pairs[i];
                const EdgeCode global_e =
                    lookup_local_edge_at_vertex(v, pair.edge);
                const SquareCode global_s =
                    lookup_local_square_at_vertex(v, pair.square);
                if (!is_h1_positive_edge(global_e) ||
                    !is_h1_negative_square(global_s)) {
                    continue;
                }
                if (global_s > static_cast<SquareCode>(pm_sq_mask)) {
                    throw std::overflow_error(
                        "3D H1 pivot map: square index does not fit in 30 bits");
                }
                pivot_map[static_cast<std::size_t>(global_e)] =
                    pm_zero_sq_flag | static_cast<std::uint32_t>(global_s);
            }
        }
    }

    void h1_collect_lookup_owner_vertices_demand(
        std::vector<H1VertexOrderItem>& vertex_order,
        H1Result& result) const {
        vertex_order.clear();
        vertex_order.reserve(n_vertices_ / 4 + 1);
        std::size_t zero_skipped = 0;
        for (std::size_t v = 0; v < n_vertices_; ++v) {
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[v]);
            const std::uint16_t zero_square_mask =
                h1_lookup_owner_zero_square_mask(v, entry);
            bool has_work = false;
            for (std::uint8_t i = 0; i < entry.square_count; ++i) {
                const std::uint8_t local_square = entry.square_order_asc[i];
                const SquareCode s =
                    lookup_local_square_at_vertex(v, local_square);
                if (!is_h1_negative_square(s)) {
                    continue;
                }
                if ((zero_square_mask & (1u << local_square)) != 0) {
                    ++zero_skipped;
                    continue;
                }
                has_work = true;
            }
            if (!has_work) {
                continue;
            }
            const std::uint32_t code = value_codes_ != nullptr
                ? value_codes_[v]
                : std::numeric_limits<std::uint32_t>::max();
            vertex_order.push_back({values_[v], code,
                                    static_cast<VertexCode>(v)});
        }
        result.active_owner_vertices = vertex_order.size();
        result.h1_zero_squares_skipped_outer = zero_skipped;
    }

    bool h1_lookup_zero_square_owner_local(
        SquareCode square,
        EdgeCode expected_pivot,
        std::size_t& owner_vertex,
        std::uint8_t& local_square) const {
        const std::size_t sidx = static_cast<std::size_t>(square);
        owner_vertex = square_birth_vertex_[sidx];
        const zero_lookup::Entry3D& entry =
            zero_lookup::lookup3d(zero_lookup_masks_[owner_vertex]);
        const H1LookupOwnerContext context =
            h1_lookup_owner_context(owner_vertex, entry);
        for (std::uint8_t i = 0; i < entry.zero_h1_count; ++i) {
            const zero_lookup::Direct3DH1Pair pair = entry.zero_h1_pairs[i];
            const SquareCode global_square =
                lookup_local_square_at_vertex(owner_vertex, pair.square);
            if (global_square != square) {
                continue;
            }
            const EdgeCode global_pivot =
                lookup_local_edge_at_vertex(owner_vertex, pair.edge);
            if (global_pivot != expected_pivot) {
                continue;
            }
            local_square = pair.square;
            return true;
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // Hot-path demand column reduction.
    //
    // Single pivot_map lookup per iteration replaces three separate array
    // reads (rewrite_index + zero_pivot_to_sq + zero_pivot_state).
    // H1RewriteRow AoS loads offset+length in one cache line instead of two.
    // -----------------------------------------------------------------------
    void h1_reduce_column_demand_template(
            std::vector<H1Term>& column,
            std::vector<std::uint32_t>& pivot_map,
            std::vector<H1RewriteRow>& rewrite_meta,
            std::vector<H1Term>& rewrite_terms,
            std::vector<H1Term>& scratch,
            H1MaterialiseBuffers& materialise_buffers,
            SquareCode current_square,
            H1Result& result) const {
        while (!column.empty()) {
            const H1Term pivot_term = column.back();
            const std::uint32_t pm =
                pivot_map[static_cast<std::size_t>(pivot_term)];
            if (profile_h1_) {
                ++result.profile.reduction_pivot_checks;
            }

            if (pm < pm_in_prog_flag) {
                // Path 1: stored row already exists — XOR it in.
                column.pop_back();
                const H1RewriteRow& row = rewrite_meta[pm];
                const std::size_t col_before = column.size();
                if (row.length > 0) {
                    xor_with_sorted_terms(column,
                        rewrite_terms.data() + row.offset,
                        row.length, scratch);
                }
                if (profile_h1_) {
                    H1Profile& p = result.profile;
                    ++p.reduction_rewrite_hits;
                    p.reduction_rewrite_terms_read += row.length;
                    p.reduction_column_terms_before_xor += col_before;
                    p.reduction_output_terms_after_xor += column.size();
                    p.reduction_estimated_merge_input_terms +=
                        col_before + row.length;
                    p.max_rewrite_terms_read =
                        std::max(p.max_rewrite_terms_read,
                                 static_cast<std::size_t>(row.length));
                    p.max_column_terms_before_xor =
                        std::max(p.max_column_terms_before_xor, col_before);
                    const std::size_t c = col_before;
                    std::size_t bi;
                    if      (c == 0)   bi = 0;
                    else if (c == 1)   bi = 1;
                    else if (c == 2)   bi = 2;
                    else if (c == 3)   bi = 3;
                    else if (c <= 5)   bi = 4;
                    else if (c <= 9)   bi = 5;
                    else if (c <= 15)  bi = 6;
                    else if (c <= 31)  bi = 7;
                    else if (c <= 63)  bi = 8;
                    else if (c <= 127) bi = 9;
                    else if (c <= 511) bi = 10;
                    else               bi = 11;
                    ++p.xor_col_size_hist[bi];
                    const std::size_t r = row.length;
                    std::size_t bj;
                    if      (r == 0)  bj = 0;
                    else if (r == 1)  bj = 1;
                    else if (r == 2)  bj = 2;
                    else if (r == 3)  bj = 3;
                    else if (r <= 7)  bj = 4;
                    else if (r <= 15) bj = 5;
                    else if (r <= 63) bj = 6;
                    else              bj = 7;
                    ++p.xor_row_size_hist[bj];
                }
                ++result.reduction_steps;
                continue;
            }

            // pm >= pm_in_prog_flag: in-progress, unmaterialised, or pm_no_info.
            // dep_sq is encoded in the low 30 bits for both in-progress and
            // unmaterialised states, enabling the own-pivot termination check.
            const SquareCode dep_sq =
                static_cast<SquareCode>(pm & pm_sq_mask);

            if (pm < pm_zero_sq_flag) {
                // In-progress: dep_sq == current_square is the own-pivot case.
                if (dep_sq == current_square) {
                    break;
                }
                const EdgeCode pivot_edge = h1_edge(pivot_term);
                ++result.h1_zero_materialisation_failures;
                throw std::logic_error(
                    "template H1 cycle detected at edge " +
                    std::to_string(pivot_edge));
            }

            if (pm == pm_no_info) {
                break; // genuine new pivot
            }

            // Path 2: unmaterialised zero square.
            const EdgeCode pivot_edge = h1_edge(pivot_term);
            if (dep_sq == current_square) {
                break;
            }
            if (!square_less(dep_sq, current_square)) {
                ++result.h1_zero_materialisation_failures;
                throw std::logic_error(
                    "template H1 ordering violation: zero square " +
                    std::to_string(dep_sq) + " is not before current square " +
                    std::to_string(current_square));
            }

            materialise_zero_h1_square_template(
                dep_sq, pivot_edge, pivot_map, rewrite_meta,
                rewrite_terms, materialise_buffers, result);

            const std::uint32_t pm2 =
                pivot_map[static_cast<std::size_t>(pivot_term)];
            if (pm2 >= pm_in_prog_flag) {
                ++result.h1_zero_materialisation_failures;
                throw std::logic_error(
                    "template H1 materialisation did not store a rewrite");
            }
            column.pop_back();
            const H1RewriteRow& row2 = rewrite_meta[pm2];
            const std::size_t col_before2 = column.size();
            if (row2.length > 0) {
                xor_with_sorted_terms(column,
                    rewrite_terms.data() + row2.offset,
                    row2.length, scratch);
            }
            if (profile_h1_) {
                H1Profile& p = result.profile;
                ++p.reduction_rewrite_hits;
                p.reduction_rewrite_terms_read += row2.length;
                p.reduction_column_terms_before_xor += col_before2;
                p.reduction_output_terms_after_xor += column.size();
                p.reduction_estimated_merge_input_terms +=
                    col_before2 + row2.length;
                p.max_rewrite_terms_read =
                    std::max(p.max_rewrite_terms_read,
                             static_cast<std::size_t>(row2.length));
                p.max_column_terms_before_xor =
                    std::max(p.max_column_terms_before_xor, col_before2);
            }
            ++result.reduction_steps;
        }
        result.max_reduced_boundary_size =
            std::max(result.max_reduced_boundary_size, column.size());
    }

    void materialise_zero_h1_square_template(
            SquareCode square,
            EdgeCode expected_pivot,
            std::vector<std::uint32_t>& pivot_map,
            std::vector<H1RewriteRow>& rewrite_meta,
            std::vector<H1Term>& rewrite_terms,
            H1MaterialiseBuffers& materialise_buffers,
            H1Result& result) const {
        const H1Term expected_pivot_term = h1_term(expected_pivot);
        const std::size_t pivot_idx =
            static_cast<std::size_t>(expected_pivot_term);

        const std::uint32_t pm_init = pivot_map[pivot_idx];
        if (pm_init < pm_in_prog_flag) {
            // Already materialised (row stored by a previous call).
            ++result.h1_zero_rewrite_cache_hits;
            return;
        }
        if (pm_init < pm_zero_sq_flag) {
            // pm_in_prog_flag <= pm_init < pm_zero_sq_flag: in-progress.
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "template H1 materialisation called on in-progress row");
        }
        if (pm_init == pm_no_info) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "template H1 materialisation called on non-zero-pair edge");
        }
        // Mark in-progress, preserving dep_sq for cycle detection.
        pivot_map[pivot_idx] =
            pm_in_prog_flag |
            static_cast<std::uint32_t>(square & pm_sq_mask);

        const std::size_t buffer_level = materialise_buffers.depth++;
        struct DepthGuard {
            H1MaterialiseBuffers& buffers;
            ~DepthGuard() { --buffers.depth; }
        } guard{materialise_buffers};
        if (materialise_buffers.columns.size() <= buffer_level) {
            materialise_buffers.columns.emplace_back();
            materialise_buffers.columns.back().reserve(16);
            materialise_buffers.scratches.emplace_back();
            materialise_buffers.scratches.back().reserve(16);
        }
        std::vector<H1Term>& col =
            materialise_buffers.columns[buffer_level];
        std::vector<H1Term>& scratch =
            materialise_buffers.scratches[buffer_level];
        col.clear();
        scratch.clear();

        const auto build_start = profile_h1_ ? Clock::now()
                                             : Clock::time_point{};
        bool used_template = false;
        if constexpr (h1_lookup_use_template_zero_rows) {
            std::size_t owner_vertex = 0;
            std::uint8_t local_square = 0xffu;
            if (!h1_lookup_zero_square_owner_local(square, expected_pivot,
                                                   owner_vertex,
                                                   local_square)) {
                ++result.h1_zero_materialisation_failures;
                throw std::logic_error(
                    "template H1 could not locate zero square in lookup entry");
            }
            const zero_lookup::Entry3D& entry =
                zero_lookup::lookup3d(zero_lookup_masks_[owner_vertex]);
            const H1LookupOwnerContext context =
                h1_lookup_owner_context(owner_vertex, entry);
            used_template = h1_initial_lookup_semireduced_column(
                owner_vertex, local_square, square, entry, context, col);
        } else {
            h1_initial_square_column(square, col);
        }
        if (profile_h1_) {
            result.h1_lookup_zero_build_ms +=
                elapsed_ms(build_start, Clock::now());
        }
        if (used_template) {
            ++result.h1_lookup_zero_template_columns;
        } else {
            ++result.h1_lookup_zero_raw_fallback_columns;
        }
        h1_record_lookup_initial_profile(col, result);

        const auto reduce_start = profile_h1_ ? Clock::now()
                                              : Clock::time_point{};
        h1_reduce_column_demand_template(
            col, pivot_map, rewrite_meta, rewrite_terms, scratch,
            materialise_buffers, square, result);
        if (profile_h1_) {
            result.h1_lookup_zero_reduce_ms +=
                elapsed_ms(reduce_start, Clock::now());
        }

        if (col.empty()) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "template H1 zero square reduced to empty column");
        }
        const H1Term actual_pivot_term = col.back();
        const EdgeCode actual_pivot = h1_edge(actual_pivot_term);
        if (actual_pivot != expected_pivot) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "template H1 zero square reduced to unexpected pivot");
        }
        if (pivot_map[static_cast<std::size_t>(actual_pivot_term)] <
                pm_in_prog_flag) {
            ++result.h1_zero_materialisation_failures;
            throw std::logic_error(
                "template H1 zero row selected an already paired edge");
        }
        col.pop_back();
        h1_store_rewrite(actual_pivot_term, col, pivot_map,
                         rewrite_meta, rewrite_terms, result);
        ++result.h1_zero_rewrites_materialised;
    }

    void h1_process_lookup_square_demand_template(
        std::size_t owner_vertex,
        std::uint8_t local_square,
        SquareCode s,
        const zero_lookup::Entry3D& entry,
        const H1LookupOwnerContext& context,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& column,
        std::vector<H1Term>& scratch,
        H1MaterialiseBuffers& materialise_buffers,
        H1Result& result) const {
        ++result.negative_squares;
        const auto build_start = profile_h1_ ? Clock::now()
                                             : Clock::time_point{};
        bool used_template = false;
        if (h1_lookup_use_template_nonzero_columns) {
            used_template = h1_initial_lookup_semireduced_column(
                owner_vertex, local_square, s, entry, context, column);
        } else {
            h1_initial_square_column(s, column);
        }
        if (profile_h1_) {
            result.h1_lookup_nonzero_build_ms +=
                elapsed_ms(build_start, Clock::now());
        }
        if (used_template) {
            ++result.h1_lookup_nonzero_template_columns;
        } else {
            ++result.h1_lookup_nonzero_raw_fallback_columns;
        }
        h1_record_lookup_initial_profile(column, result);

        const auto reduce_start = profile_h1_ ? Clock::now()
                                              : Clock::time_point{};
        h1_reduce_column_demand_template(
            column, pivot_map, rewrite_meta, rewrite_terms, scratch,
            materialise_buffers, s, result);
        if (profile_h1_) {
            result.h1_lookup_nonzero_reduce_ms +=
                elapsed_ms(reduce_start, Clock::now());
        }

        if (column.empty()) {
            throw std::logic_error(
                "H1-negative square reduced to zero in demand-template path");
        }
        const H1Term pivot = column.back();
        if (pivot_map[static_cast<std::size_t>(pivot)] < pm_in_prog_flag) {
            throw std::logic_error(
                "H1 demand-template reduction selected an already paired edge");
        }
        const EdgeCode pivot_edge = h1_edge(pivot);
        ++result.matched_pairs;
        const double birth = edge_birth_value(pivot_edge);
        const double death = square_birth_value(s);
        if (death > birth || include_apparent_zero_pairs_) {
            result.pairs.push_back({birth, death});
        } else {
            ++result.zero_persistence_pairs;
        }
        column.pop_back();
        h1_store_rewrite(pivot, column, pivot_map, rewrite_meta,
                         rewrite_terms, result);
    }

    void h1_process_lookup_square_raw_demand(
        SquareCode s,
        std::vector<std::uint32_t>& pivot_map,
        std::vector<H1RewriteRow>& rewrite_meta,
        std::vector<H1Term>& rewrite_terms,
        std::vector<H1Term>& column,
        std::vector<H1Term>& scratch,
        H1MaterialiseBuffers& materialise_buffers,
        H1Result& result) const {
        ++result.negative_squares;
        const auto nz_build_start = profile_h1_ ? Clock::now()
                                                 : Clock::time_point{};
        h1_initial_square_column(s, column);
        if (profile_h1_) {
            result.h1_lookup_nonzero_build_ms +=
                elapsed_ms(nz_build_start, Clock::now());
        }
        ++result.h1_lookup_nonzero_raw_fallback_columns;
        h1_record_lookup_initial_profile(column, result);

        const auto nz_reduce_start = profile_h1_ ? Clock::now()
                                                  : Clock::time_point{};
        h1_reduce_column_demand_template(
            column, pivot_map, rewrite_meta, rewrite_terms, scratch,
            materialise_buffers, s, result);
        if (profile_h1_) {
            result.h1_lookup_nonzero_reduce_ms +=
                elapsed_ms(nz_reduce_start, Clock::now());
        }

        if (column.empty()) {
            throw std::logic_error(
                "H1-negative square reduced to zero in raw-demand path");
        }
        const H1Term pivot = column.back();
        if (pivot_map[static_cast<std::size_t>(pivot)] < pm_in_prog_flag) {
            throw std::logic_error(
                "H1 raw-demand reduction selected an already paired edge");
        }
        const EdgeCode pivot_edge = h1_edge(pivot);
        ++result.matched_pairs;
        const double birth = edge_birth_value(pivot_edge);
        const double death = square_birth_value(s);
        if (death > birth || include_apparent_zero_pairs_) {
            result.pairs.push_back({birth, death});
        } else {
            ++result.zero_persistence_pairs;
        }
        column.pop_back();
        h1_store_rewrite(pivot, column, pivot_map, rewrite_meta,
                         rewrite_terms, result);
    }

    H1Result run_h1() {
        H1Result result;
        result.profile.enabled = profile_h1_;
        const auto total_start = Clock::now();

        const auto setup_start = Clock::now();
        ensure_square_birth_edges();
        const std::size_t edge_storage = edge_state_storage_size(
            depth_, rows_, cols_);
        if (edge_storage > std::numeric_limits<H1Term>::max()) {
            throw std::overflow_error(
                "3D H1 matching currently expects edge ids to fit in uint32");
        }

        // Combined pivot map: replaces rewrite_index (12 MB) +
        // zero_pivot_to_sq (12 MB) + zero_pivot_state (3 MB) with a single
        // 12 MB array.  Initialised to pm_no_info; zero-pair entries are
        // populated by h1_build_zero_pivot_map.
        std::vector<std::uint32_t> pivot_map(edge_storage, pm_no_info);
        h1_build_zero_pivot_map(pivot_map);

        result.timing.parent_initialization_ms =
            elapsed_ms(setup_start, Clock::now());

        const std::size_t expected_rewrites = edge_storage / 3u * 2u + 1024u;

        const auto collect_start = Clock::now();
        std::vector<H1VertexOrderItem> vertex_order;
        h1_collect_lookup_owner_vertices_demand(vertex_order, result);
        result.timing.apparent_pair_scan_ms =
            elapsed_ms(collect_start, Clock::now());

        const auto sort_start = Clock::now();
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
        result.timing.sorting_ms = elapsed_ms(sort_start, Clock::now());

        // AoS rewrite metadata: replaces separate rewrite_offsets + rewrite_lengths.
        std::vector<H1RewriteRow> rewrite_meta;
        std::vector<H1Term> rewrite_terms;
        rewrite_meta.reserve(expected_rewrites);
        if (expected_rewrites <=
            std::numeric_limits<std::size_t>::max() / 6u) {
            rewrite_terms.reserve(expected_rewrites * 6u);
        }
        std::vector<H1Term> column;
        std::vector<H1Term> scratch;
        column.reserve(16);
        scratch.reserve(16);
        H1MaterialiseBuffers materialise_buffers;
        result.pairs.reserve(expected_rewrites / 4 + 1);

        const auto reduction_start = Clock::now();
        for_each_h1_lookup_nonzero_square_raw_order(
            vertex_order, result,
            [&](SquareCode s) {
                h1_process_lookup_square_raw_demand(
                    s, pivot_map, rewrite_meta, rewrite_terms,
                    column, scratch, materialise_buffers, result);
            });
        result.timing.union_find_sweep_ms =
            elapsed_ms(reduction_start, Clock::now());
        result.timing.total_ms = elapsed_ms(total_start, Clock::now());
        result.positive_edges = result.matched_pairs;

        result.memory.parent_bytes =
            pivot_map.size() * sizeof(std::uint32_t);
        result.memory.non_apparent_edges_bytes =
            rewrite_meta.capacity() * sizeof(H1RewriteRow) +
            rewrite_terms.capacity() * sizeof(H1Term);
        result.memory.active_vertices_bytes =
            vertex_order.capacity() * sizeof(H1VertexOrderItem) +
            column.capacity() * sizeof(H1Term) +
            scratch.capacity() * sizeof(H1Term);
        result.memory.persistence_pairs_bytes =
            result.pairs.capacity() * sizeof(PersistencePair);
        return result;
    }
