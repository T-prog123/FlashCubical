#pragma once

#include "algorithms/detail/smart_h0_2d_common.hpp"
#include "zero_lookup/runtime/zero_table.hpp"

namespace smart_h0 {

class H0Computer {
public:
    H0Computer(const double* values, std::size_t rows, std::size_t cols,
               bool include_apparent_zero_pairs,
               TimingProfile* timing,
               std::uint8_t* edge_state = nullptr,
               bool skip_known_h1_positive = true,
               const std::uint8_t* h1_positive_edge_mask = nullptr,
               const std::uint32_t* value_codes = nullptr,
               std::uint8_t* h1_candidate_edge_owner_mask = nullptr,
               std::vector<VertexCode>* h1_candidate_active_vertices = nullptr,
               const std::uint8_t* zero_lookup_masks = nullptr,
               zero_lookup::RunStats* zero_lookup_stats = nullptr)
        : values_(values),
          value_codes_(value_codes),
          rows_(rows),
          cols_(cols),
          n_vertices_(0),
          include_apparent_zero_pairs_(include_apparent_zero_pairs),
          timing_(timing),
          edge_state_(edge_state),
          skip_known_h1_positive_(skip_known_h1_positive),
          h1_positive_edge_mask_(h1_positive_edge_mask),
          h1_candidate_edge_owner_mask_(h1_candidate_edge_owner_mask),
          h1_candidate_active_vertices_(h1_candidate_active_vertices),
          zero_lookup_masks_(zero_lookup_masks),
          zero_lookup_stats_(zero_lookup_stats) {
        if (values_ == nullptr) {
            throw std::invalid_argument("image value pointer is null");
        }
        if (rows_ == 0 || cols_ == 0) {
            throw std::invalid_argument("image dimensions must be positive");
        }
        if (rows_ > (std::numeric_limits<std::size_t>::max() / cols_)) {
            throw std::overflow_error("image dimensions overflow size_t");
        }

        const auto start = Clock::now();
        n_vertices_ = rows_ * cols_;
        if (n_vertices_ > std::numeric_limits<VertexCode>::max()) {
            throw std::overflow_error("image has too many vertices for packed vertex ids");
        }
        parent_.reset(n_vertices_);
        if (timing_ != nullptr) {
            timing_->parent_initialization_ms += elapsed_ms(start, Clock::now());
        }
    }

    void mark_apparent_negative_edges_only() {
        if (edge_state_ == nullptr) {
            return;
        }

        for (std::size_t y = 0; y < n_vertices_; ++y) {
            const std::size_t r = row(y);
            const std::size_t c = col(y);

            EdgeCode lower_edges[4];
            int lower_count = 0;

            auto consider = [&](std::size_t z) {
                if (vertex_less(z, y)) {
                    lower_edges[lower_count++] = edge_between(y, z);
                }
            };

            if (c > 0) {
                consider(y - 1);
            }
            if (c + 1 < cols_) {
                consider(y + 1);
            }
            if (r > 0) {
                consider(y - cols_);
            }
            if (r + 1 < rows_) {
                consider(y + cols_);
            }

            if (lower_count == 0) {
                continue;
            }

            int apparent_index = 0;
            for (int i = 1; i < lower_count; ++i) {
                if (lower_edges[i] < lower_edges[apparent_index]) {
                    apparent_index = i;
                }
            }
            edge_state_mark(lower_edges[apparent_index],
                            edge_state_h0_apparent | edge_state_h0_negative);
        }
    }

    Result run() {
        Result result;
        if (timing_ != nullptr) {
            result.timing.parent_initialization_ms = timing_->parent_initialization_ms;
        }

        const auto work_storage_start = Clock::now();
        result.pairs.reserve(include_apparent_zero_pairs_
            ? n_vertices_
            : (n_vertices_ / 4 + 1));

        std::unique_ptr<std::uint8_t[]> non_apparent_vertex_mask(
            new std::uint8_t[n_vertices_]());

        std::vector<VertexCode> active_vertices;
        active_vertices.reserve(n_vertices_ / 4 + 1);
        result.timing.parent_initialization_ms +=
            elapsed_ms(work_storage_start, Clock::now());

        const auto scan_start = Clock::now();
        scan_edges_for_apparent_pairs(non_apparent_vertex_mask.get(),
                                      active_vertices, result);
        result.timing.apparent_pair_scan_ms = elapsed_ms(scan_start, Clock::now());
        result.timing.apparent_pair_work_ms = result.timing.apparent_pair_scan_ms;

        const auto sort_start = Clock::now();
        sort_active_vertices_ascending(active_vertices);
        result.timing.sorting_ms = elapsed_ms(sort_start, Clock::now());

        result.active_non_apparent_vertices = active_vertices.size();
        const auto sweep_start = Clock::now();
        sweep_non_apparent_edges(active_vertices, non_apparent_vertex_mask.get(),
                                 result);
        result.timing.union_find_sweep_ms = elapsed_ms(sweep_start, Clock::now());

        const auto final_start = Clock::now();
        const std::size_t oldest = oldest_vertex();
        const std::size_t surviving_root = find(oldest);
        result.pairs.push_back({values_[surviving_root],
                                std::numeric_limits<double>::infinity()});
        result.timing.finalization_ms = elapsed_ms(final_start, Clock::now());

        result.memory.parent_bytes = n_vertices_ * sizeof(VertexCode);
        result.memory.non_apparent_edges_bytes = 0;
        result.memory.non_apparent_vertex_mask_bytes = n_vertices_ * sizeof(std::uint8_t);
        result.memory.active_vertices_bytes =
            active_vertices.capacity() * sizeof(VertexCode);
        result.memory.persistence_pairs_bytes = result.pairs.capacity() * sizeof(PersistencePair);

        return result;
    }

private:
    using Clock = std::chrono::steady_clock;

    const double* values_;
    const std::uint32_t* value_codes_;
    std::size_t rows_;
    std::size_t cols_;
    std::size_t n_vertices_;
    bool include_apparent_zero_pairs_;
    TimingProfile* timing_;
    std::uint8_t* edge_state_;
    bool skip_known_h1_positive_;
    const std::uint8_t* h1_positive_edge_mask_;
    std::uint8_t* h1_candidate_edge_owner_mask_;
    std::vector<VertexCode>* h1_candidate_active_vertices_;
    const std::uint8_t* zero_lookup_masks_;
    zero_lookup::RunStats* zero_lookup_stats_;
    smart_core::UnionFind<VertexCode> parent_;

    static double elapsed_ms(Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

    struct VertexKey {
        double value;
        std::size_t index;
        std::uint32_t code;
    };

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

    static bool key_equal(VertexKey a, VertexKey b) {
        if (a.code != std::numeric_limits<std::uint32_t>::max()) {
            return a.code == b.code && a.index == b.index;
        }
        return a.value == b.value && a.index == b.index;
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

    static std::size_t popcount4(std::uint8_t x) {
        x &= 0x0fu;
        return static_cast<std::size_t>((x & 1u) + ((x >> 1) & 1u) +
                                        ((x >> 2) & 1u) + ((x >> 3) & 1u));
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

    bool h1_positive_edge_mask_has(EdgeCode e) const {
        return h1_positive_edge_mask_ != nullptr &&
               (h1_positive_edge_mask_[edge_source(e)] &
                static_cast<std::uint8_t>(1u << edge_dir(e))) != 0;
    }

    bool should_skip_known_h1_positive(EdgeCode e) const {
        if (!skip_known_h1_positive_) {
            return false;
        }
        if (h1_positive_edge_mask_ != nullptr) {
            return h1_positive_edge_mask_has(e);
        }
        return edge_state_has(e, edge_state_h1_positive);
    }

    std::pair<std::size_t, std::size_t> endpoints(EdgeCode e) const {
        const std::size_t s = edge_source(e);
        if (edge_dir(e) == 0) {
            return {s, s + 1};
        }
        return {s, s + cols_};
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
        throw std::logic_error("vertices are not grid neighbors");
    }

    VertexKey edge_birth_key(EdgeCode e) const {
        const auto [u, v] = endpoints(e);
        return vertex_greater(u, v) ? vertex_key(u) : vertex_key(v);
    }

    bool edge_less(EdgeCode a, EdgeCode b) const {
        const VertexKey ba = edge_birth_key(a);
        const VertexKey bb = edge_birth_key(b);
        if (key_less(ba, bb)) {
            return true;
        }
        if (key_less(bb, ba)) {
            return false;
        }
        return a < b;
    }

    void sort_active_vertices_ascending(std::vector<VertexCode>& active_vertices) const {
        if (value_codes_ == nullptr) {
            if (active_vertices.size() > 4096) {
                radix_sort_active_vertices_by_value(active_vertices);
                return;
            }
            std::sort(active_vertices.begin(), active_vertices.end(),
                      [this](VertexCode a, VertexCode b) {
                          return key_less(vertex_key(a), vertex_key(b));
                      });
            return;
        }

        std::uint32_t max_code = 0;
        for (VertexCode v : active_vertices) {
            if (value_codes_[v] > max_code) {
                max_code = value_codes_[v];
            }
        }
        if (max_code <= 10000000u) {
            std::vector<std::uint32_t> offsets(static_cast<std::size_t>(max_code) + 1u, 0);
            for (VertexCode v : active_vertices) {
                ++offsets[value_codes_[v]];
            }

            std::uint32_t running = 0;
            for (std::uint32_t& count : offsets) {
                const std::uint32_t bucket_size = count;
                count = running;
                running += bucket_size;
            }

            std::vector<VertexCode> sorted(active_vertices.size());
            for (VertexCode v : active_vertices) {
                sorted[offsets[value_codes_[v]]++] = v;
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
        for (std::size_t i = 0; i < keys.size(); ++i) {
            active_vertices[i] = static_cast<VertexCode>(keys[i]);
        }
    }

    static std::uint64_t sortable_double_key(double x) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &x, sizeof(bits));
        if ((bits & (1ull << 63)) != 0) {
            return ~bits;
        }
        return bits ^ (1ull << 63);
    }

    void radix_sort_active_vertices_by_value(
        std::vector<VertexCode>& active_vertices) const {
        const std::size_t n = active_vertices.size();
        std::vector<VertexCode> tmp_vertices(n);
        std::vector<std::uint64_t> keys(n);
        std::vector<std::uint64_t> tmp_keys(n);

        for (std::size_t i = 0; i < n; ++i) {
            keys[i] = sortable_double_key(values_[active_vertices[i]]);
        }

        std::array<std::uint32_t, 65536> counts{};
        for (unsigned shift = 0; shift < 64; shift += 16) {
            counts.fill(0);
            for (std::uint64_t key : keys) {
                ++counts[static_cast<std::uint16_t>(key >> shift)];
            }

            std::uint32_t running = 0;
            for (std::uint32_t& count : counts) {
                const std::uint32_t bucket_size = count;
                count = running;
                running += bucket_size;
            }

            for (std::size_t i = 0; i < n; ++i) {
                const std::uint16_t bucket =
                    static_cast<std::uint16_t>(keys[i] >> shift);
                const std::uint32_t out = counts[bucket]++;
                tmp_vertices[out] = active_vertices[i];
                tmp_keys[out] = keys[i];
            }

            active_vertices.swap(tmp_vertices);
            keys.swap(tmp_keys);
        }
    }

    bool edge_equal(EdgeCode a, EdgeCode b) const {
        return a == b;
    }

    std::size_t edge_count() const {
        return rows_ * (cols_ - 1) + (rows_ - 1) * cols_;
    }

    std::size_t find(std::size_t v) {
        return parent_.find(static_cast<VertexCode>(v));
    }

    std::size_t oldest_vertex() const {
        std::size_t best = 0;
        for (std::size_t v = 1; v < n_vertices_; ++v) {
            if (vertex_less(v, best)) {
                best = v;
            }
        }
        return best;
    }

    static constexpr std::uint8_t edge_left_bit = 1u << 0;
    static constexpr std::uint8_t edge_right_bit = 1u << 1;
    static constexpr std::uint8_t edge_up_bit = 1u << 2;
    static constexpr std::uint8_t edge_down_bit = 1u << 3;

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

    void local_edge_parts(std::size_t y,
                          std::uint8_t local,
                          EdgeCode& e,
                          std::size_t& other,
                          std::uint8_t& owned_bit) const {
        switch (local) {
        case zero_lookup::edge2_n:
            e = pack_edge(y - cols_, 1);
            other = y - cols_;
            owned_bit = edge_up_bit;
            break;
        case zero_lookup::edge2_e:
            e = pack_edge(y, 0);
            other = y + 1;
            owned_bit = edge_right_bit;
            break;
        case zero_lookup::edge2_s:
            e = pack_edge(y, 1);
            other = y + cols_;
            owned_bit = edge_down_bit;
            break;
        case zero_lookup::edge2_w:
            e = pack_edge(y - 1, 0);
            other = y - 1;
            owned_bit = edge_left_bit;
            break;
        default:
            throw std::logic_error("invalid local 2D edge id");
        }
    }

    void mark_non_apparent_edge(std::uint8_t* non_apparent_vertex_mask,
                                std::vector<VertexCode>& active_vertices,
                                std::size_t owner,
                                EdgeCode e,
                                Result& result) const {
        const std::uint8_t bit = owned_edge_bit(e, owner);
        if (non_apparent_vertex_mask[owner] == 0) {
            active_vertices.push_back(static_cast<VertexCode>(owner));
        }
        non_apparent_vertex_mask[owner] |= bit;
        ++result.non_apparent_edges;
    }

    void mark_h1_candidate_edge(EdgeCode e) const {
        if (h1_candidate_edge_owner_mask_ == nullptr ||
            h1_candidate_active_vertices_ == nullptr) {
            return;
        }

        const auto [u, v] = endpoints(e);
        const std::size_t owner = vertex_greater(u, v) ? u : v;
        const std::uint8_t bit = owned_edge_bit(e, owner);
        if (h1_candidate_edge_owner_mask_[owner] == 0) {
            h1_candidate_active_vertices_->push_back(static_cast<VertexCode>(owner));
        }
        h1_candidate_edge_owner_mask_[owner] |= bit;
    }

    void scan_edges_for_apparent_pairs(std::uint8_t* non_apparent_vertex_mask,
                                       std::vector<VertexCode>& active_vertices,
                                       Result& result) {
        if (zero_lookup_masks_ != nullptr) {
            scan_edges_from_lookup_masks(non_apparent_vertex_mask,
                                         active_vertices, result);
            return;
        }
        for (std::size_t y = 0; y < n_vertices_; ++y) {
            const std::size_t r = row(y);
            const std::size_t c = col(y);

            EdgeCode lower_edges[4];
            std::uint8_t lower_bits[4];
            std::uint8_t lookup_mask = 0;
            int lower_count = 0;

            auto consider = [&](std::uint8_t local_bit,
                                std::size_t z,
                                EdgeCode e) {
                if (vertex_less(z, y)) {
                    lower_bits[lower_count] = local_bit;
                    lower_edges[lower_count++] = e;
                    lookup_mask = static_cast<std::uint8_t>(
                        lookup_mask | (1u << local_bit));
                }
            };

            if (c > 0) {
                consider(zero_lookup::edge2_w, y - 1, pack_edge(y - 1, 0));
            }
            if (c + 1 < cols_) {
                consider(zero_lookup::edge2_e, y + 1, pack_edge(y, 0));
            }
            if (r > 0) {
                consider(zero_lookup::edge2_n, y - cols_, pack_edge(y - cols_, 1));
            }
            if (r + 1 < rows_) {
                consider(zero_lookup::edge2_s, y + cols_, pack_edge(y, 1));
            }

            if (lower_count == 0) {
                continue;
            }

            const std::uint8_t apparent_local =
                zero_lookup::direct_h0_edge2d(lookup_mask);
            int apparent_index = -1;
            for (int i = 0; i < lower_count; ++i) {
                if (lower_bits[i] == apparent_local) {
                    apparent_index = i;
                    break;
                }
            }
            if (apparent_index < 0) {
                apparent_index = 0;
                for (int i = 1; i < lower_count; ++i) {
                    if (lower_edges[i] < lower_edges[apparent_index]) {
                        apparent_index = i;
                    }
                }
            }

            const EdgeCode apparent = lower_edges[apparent_index];
            const auto [u, v] = endpoints(apparent);
            parent_[y] = static_cast<VertexCode>((u == y) ? v : u);
            ++result.uf_merge_attempts;
            ++result.uf_successful_merges;
            edge_state_mark(apparent,
                            edge_state_h0_apparent | edge_state_h0_negative);
            ++result.apparent_edges;
            if (include_apparent_zero_pairs_) {
                result.pairs.push_back({values_[y], edge_birth_key(apparent).value});
            }

            for (int i = 0; i < lower_count; ++i) {
                if (i == apparent_index) {
                    continue;
                }
                if (should_skip_known_h1_positive(lower_edges[i])) {
                    ++result.skipped_known_h1_positive_edges;
                    continue;
                }
                mark_non_apparent_edge(non_apparent_vertex_mask, active_vertices,
                                       y, lower_edges[i], result);
            }
        }
    }

    EdgeCode local_edge(std::size_t y, std::uint8_t local) const {
        switch (local) {
        case zero_lookup::edge2_n:
            return pack_edge(y - cols_, 1);
        case zero_lookup::edge2_e:
            return pack_edge(y, 0);
        case zero_lookup::edge2_s:
            return pack_edge(y, 1);
        case zero_lookup::edge2_w:
            return pack_edge(y - 1, 0);
        default:
            throw std::logic_error("invalid local 2D edge id");
        }
    }

    void scan_edges_from_lookup_masks(std::uint8_t* non_apparent_vertex_mask,
                                      std::vector<VertexCode>& active_vertices,
                                      Result& result) {
        const auto apply_start = Clock::now();
        for (std::size_t y = 0; y < n_vertices_; ++y) {
            const zero_lookup::Entry2D& entry =
                zero_lookup::lookup2d(zero_lookup_masks_[y]);
            const std::uint8_t direct_h0 = entry.direct_h0_edge;
            if (zero_lookup_stats_ != nullptr) {
                if (direct_h0 != 0xffu) {
                    ++zero_lookup_stats_->direct_h0_pairs_found;
                }
            }
            result.skipped_known_h1_positive_edges +=
                entry.edge_h1_positive_count;

            if (direct_h0 != 0xffu) {
                const EdgeCode apparent = local_edge(y, direct_h0);
                const auto [u, v] = endpoints(apparent);
                parent_[y] = static_cast<VertexCode>((u == y) ? v : u);
                ++result.uf_merge_attempts;
                ++result.uf_successful_merges;
                edge_state_mark(apparent,
                                edge_state_h0_apparent | edge_state_h0_negative);
                ++result.apparent_edges;
                if (zero_lookup_stats_ != nullptr) {
                    ++zero_lookup_stats_->direct_h0_pairs_applied;
                }
                if (include_apparent_zero_pairs_) {
                    result.pairs.push_back({values_[y], edge_birth_key(apparent).value});
                }
            }

            std::uint8_t residual = entry.survivor_edge_mask;
            while (residual != 0) {
                std::uint8_t local = 0;
                while ((residual & static_cast<std::uint8_t>(1u << local)) == 0) {
                    ++local;
                }
                residual = static_cast<std::uint8_t>(
                    residual & ~static_cast<std::uint8_t>(1u << local));
                const EdgeCode e = local_edge(y, local);
                if (should_skip_known_h1_positive(e)) {
                    ++result.skipped_known_h1_positive_edges;
                    continue;
                }
                mark_non_apparent_edge(non_apparent_vertex_mask, active_vertices,
                                       y, e, result);
            }
        }
        if (zero_lookup_stats_ != nullptr) {
            zero_lookup_stats_->table_apply_ms +=
                elapsed_ms(apply_start, Clock::now());
        }
    }

    void sweep_edge(EdgeCode e, Result& result) {
        if (should_skip_known_h1_positive(e)) {
            ++result.skipped_known_h1_positive_edges;
            return;
        }

        const auto [u, v] = endpoints(e);
        std::size_t ru = find(u);
        std::size_t rv = find(v);
        ++result.uf_merge_attempts;
        if (ru == rv) {
            ++result.uf_same_component_attempts;
            mark_h1_candidate_edge(e);
            return;
        }

        std::size_t younger = ru;
        std::size_t older = rv;
        if (vertex_less(ru, rv)) {
            younger = rv;
            older = ru;
        }

        const double birth = values_[younger];
        const double death = edge_birth_key(e).value;
        if (include_apparent_zero_pairs_ || death > birth) {
            result.pairs.push_back({birth, death});
        }
        edge_state_mark(e, edge_state_h0_negative);
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
    }

    void mark_h1_candidate_owned(std::size_t owner, std::uint8_t bit) const {
        if (h1_candidate_edge_owner_mask_ == nullptr ||
            h1_candidate_active_vertices_ == nullptr) {
            return;
        }
        if (h1_candidate_edge_owner_mask_[owner] == 0) {
            h1_candidate_active_vertices_->push_back(
                static_cast<VertexCode>(owner));
        }
        h1_candidate_edge_owner_mask_[owner] |= bit;
    }

    void sweep_lookup_edge(std::size_t y,
                           std::uint8_t local,
                           Result& result) {
        EdgeCode e = 0;
        std::size_t z = y;
        std::uint8_t owned_bit = 0;
        local_edge_parts(y, local, e, z, owned_bit);

        const std::size_t ry = find(y);
        const std::size_t rz = find(z);
        ++result.uf_merge_attempts;
        if (ry == rz) {
            ++result.uf_same_component_attempts;
            mark_h1_candidate_owned(y, owned_bit);
            return;
        }

        std::size_t younger = ry;
        std::size_t older = rz;
        if (vertex_less(ry, rz)) {
            younger = rz;
            older = ry;
        }

        const double birth = values_[younger];
        const double death = values_[y];
        if (include_apparent_zero_pairs_ || death > birth) {
            result.pairs.push_back({birth, death});
        }
        edge_state_mark(e, edge_state_h0_negative);
        parent_[younger] = static_cast<VertexCode>(older);
        ++result.uf_successful_merges;
    }

    void sweep_non_apparent_edges(const std::vector<VertexCode>& active_vertices,
                                  const std::uint8_t* non_apparent_vertex_mask,
                                  Result& result) {
        for (VertexCode owner_code : active_vertices) {
            const std::size_t y = owner_code;
            const std::uint8_t mask = non_apparent_vertex_mask[y];
            if (zero_lookup_masks_ != nullptr) {
                const zero_lookup::Entry2D& entry =
                    zero_lookup::lookup2d(zero_lookup_masks_[y]);
                for (std::uint8_t i = 0;
                     i < entry.residual_edge_order_asc_count; ++i) {
                    const std::uint8_t local = entry.residual_edge_order_asc[i];
                    if ((mask & local_edge_owned_bit(local)) != 0) {
                        sweep_lookup_edge(y, local, result);
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

            sort_edges_ascending(edges, count);
            for (int i = 0; i < count; ++i) {
                sweep_edge(edges[i], result);
            }
        }
    }
};

inline Result compute(const double* values, std::size_t rows, std::size_t cols,
                      bool include_apparent_zero_pairs,
                      std::uint8_t* edge_state,
                      bool skip_known_h1_positive = true,
                      const std::uint8_t* h1_positive_edge_mask = nullptr,
                      const std::uint32_t* value_codes = nullptr,
                      std::uint8_t* h1_candidate_edge_owner_mask = nullptr,
                      std::vector<VertexCode>* h1_candidate_active_vertices = nullptr,
                      const std::uint8_t* zero_lookup_masks = nullptr,
                      zero_lookup::RunStats* zero_lookup_stats = nullptr) {
    TimingProfile timing;
    const auto total_start = std::chrono::steady_clock::now();
    H0Computer computer(values, rows, cols, include_apparent_zero_pairs,
                        &timing, edge_state, skip_known_h1_positive,
                        h1_positive_edge_mask, value_codes,
                        h1_candidate_edge_owner_mask,
                        h1_candidate_active_vertices,
                        zero_lookup_masks, zero_lookup_stats);
    Result result = computer.run();
    result.timing.parent_initialization_ms = timing.parent_initialization_ms;
    result.timing.total_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - total_start).count();
    return result;
}

inline Result compute(const double* values, std::size_t rows, std::size_t cols,
                      bool include_apparent_zero_pairs = false) {
    return compute(values, rows, cols, include_apparent_zero_pairs, nullptr);
}

inline void mark_apparent_h0_negative_edges(const double* values,
                                            std::size_t rows,
                                            std::size_t cols,
                                            std::uint8_t* edge_state,
                                            const std::uint32_t* value_codes = nullptr) {
    H0Computer(values, rows, cols, false, nullptr, edge_state, true, nullptr,
               value_codes)
        .mark_apparent_negative_edges_only();
}

} // namespace smart_h0
