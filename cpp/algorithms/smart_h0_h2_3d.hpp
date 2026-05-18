#pragma once

#include "core/states.hpp"
#include "core/types.hpp"
#include "core/union_find.hpp"
#include "zero_lookup/runtime/zero_table.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace smart_h0_h2_3d {

using smart_core::MemoryProfile;
using smart_core::PersistencePair;
using smart_core::TimingProfile;
using smart_core::edge_state_h0_apparent;
using smart_core::edge_state_h0_negative;
using smart_core::edge_state_h1_apparent;
using smart_core::edge_state_h1_positive;
using smart_core::square_state_dual_h0_negative;
using smart_core::square_state_dual_h1_apparent;
using smart_core::square_state_dual_h1_positive;
using smart_core::square_state_h2_apparent;

using VertexCode = std::uint32_t;
using CubeCode = std::uint32_t;
using EdgeCode = std::uint64_t;
using SquareCode = std::uint64_t;

inline std::size_t checked_vertex_count(std::size_t depth,
                                        std::size_t rows,
                                        std::size_t cols) {
    if (depth == 0 || rows == 0 || cols == 0) {
        return 0;
    }
    if (depth > std::numeric_limits<std::size_t>::max() / rows) {
        throw std::overflow_error("3D image dimensions overflow size_t");
    }
    const std::size_t dr = depth * rows;
    if (dr > std::numeric_limits<std::size_t>::max() / cols) {
        throw std::overflow_error("3D image dimensions overflow size_t");
    }
    return dr * cols;
}

inline std::size_t edge_state_storage_size(std::size_t depth,
                                           std::size_t rows,
                                           std::size_t cols) {
    const std::size_t n = checked_vertex_count(depth, rows, cols);
    if (n > std::numeric_limits<std::size_t>::max() / 3u) {
        throw std::overflow_error("3D edge-state table overflows size_t");
    }
    return 3u * n;
}

inline std::size_t square_state_storage_size(std::size_t depth,
                                             std::size_t rows,
                                             std::size_t cols) {
    return edge_state_storage_size(depth, rows, cols);
}

struct H0Result {
    std::vector<PersistencePair> pairs;
    std::size_t apparent_edges = 0;
    std::size_t non_apparent_edges = 0;
    std::size_t active_non_apparent_vertices = 0;
    std::size_t skipped_known_h1_positive_edges = 0;
    std::size_t uf_merge_attempts = 0;
    std::size_t uf_successful_merges = 0;
    std::size_t uf_same_component_attempts = 0;
    TimingProfile timing;
    MemoryProfile memory;
};

struct H2Result {
    std::vector<PersistencePair> pairs;
    std::size_t primal_squares = 0;
    std::size_t dual_vertices = 0;
    std::size_t apparent_squares = 0;
    std::size_t non_apparent_squares = 0;
    std::size_t active_non_apparent_cubes = 0;
    std::size_t active_lookup_owner_vertices = 0;
    std::size_t skipped_known_dual_h1_positive_squares = 0;
    std::size_t uf_merge_attempts = 0;
    std::size_t uf_successful_merges = 0;
    std::size_t uf_same_component_attempts = 0;
    TimingProfile timing;
    MemoryProfile memory;
};

struct H1Profile {
    bool enabled = false;

    std::size_t stream_vertices = 0;
    std::size_t stream_vertices_with_lower_edges = 0;
    std::size_t stream_lower_edges = 0;
    std::size_t stream_incident_square_refs = 0;
    std::size_t stream_square_birth_tests = 0;
    std::size_t stream_square_birth_matches = 0;
    std::size_t stream_h1_negative_squares = 0;
    std::size_t stream_born_square_sort_groups = 0;
    std::size_t stream_born_square_sort_items = 0;

    std::size_t initial_columns = 0;
    std::size_t initial_boundary_edge_tests = 0;
    std::size_t initial_positive_terms = 0;
    std::array<std::size_t, 5> initial_column_size_hist{};

    std::size_t apparent_fast_path_columns = 0;
    std::size_t nonapparent_reduced_columns = 0;

    std::size_t reduction_pivot_checks = 0;
    std::size_t reduction_rewrite_hits = 0;
    std::size_t reduction_rewrite_terms_read = 0;
    std::size_t reduction_column_terms_before_xor = 0;
    std::size_t reduction_output_terms_after_xor = 0;
    std::size_t reduction_estimated_merge_input_terms = 0;
    std::size_t max_rewrite_terms_read = 0;
    std::size_t max_column_terms_before_xor = 0;

    // Column size before XOR at each rewrite hit (buckets: 0,1,2,3,4-5,6-9,
    // 10-15,16-31,32-63,64-127,128-511,512+)
    std::array<std::size_t, 12> xor_col_size_hist{};
    // Rewrite-row size at each XOR (buckets: 0,1,2,3,4-7,8-15,16-63,64+)
    std::array<std::size_t, 8>  xor_row_size_hist{};

    std::size_t rewrite_rows = 0;
    std::size_t rewrite_terms_stored = 0;
    std::size_t max_rewrite_terms_stored = 0;
};

struct H1Result {
    std::vector<PersistencePair> pairs;
    std::size_t positive_edges = 0;
    std::size_t negative_squares = 0;
    std::size_t active_owner_vertices = 0;
    std::size_t apparent_pairs = 0;
    std::size_t matched_pairs = 0;
    std::size_t zero_persistence_pairs = 0;
    std::size_t reduction_steps = 0;
    std::size_t max_reduced_boundary_size = 0;
    // Lookup zero-persistence H1 counters.
    std::size_t h1_zero_squares_skipped_outer = 0;
    std::size_t h1_zero_rewrites_materialised = 0;
    std::size_t h1_zero_rewrite_cache_hits = 0;
    std::size_t h1_zero_materialisation_failures = 0;
    std::size_t h1_lookup_zero_template_columns = 0;
    std::size_t h1_lookup_zero_raw_fallback_columns = 0;
    std::size_t h1_lookup_nonzero_template_columns = 0;
    std::size_t h1_lookup_nonzero_raw_fallback_columns = 0;
    double h1_lookup_zero_build_ms = 0.0;
    double h1_lookup_zero_reduce_ms = 0.0;
    double h1_lookup_nonzero_build_ms = 0.0;
    double h1_lookup_nonzero_reduce_ms = 0.0;
    H1Profile profile;
    TimingProfile timing;
    MemoryProfile memory;
};

struct ApparentPrepassStats {
    std::size_t primal_h1_apparent_pairs = 0;
    std::size_t dual_h1_apparent_pairs = 0;
    std::size_t primal_h1_positive_edges = 0;
    std::size_t dual_h1_positive_squares = 0;
    double primal_h1_ms = 0.0;
    double dual_h1_ms = 0.0;
};

struct FullResult {
    H0Result h0;
    H1Result h1;
    H2Result h2;
    ApparentPrepassStats apparent;
    zero_lookup::RunStats zero_lookup;
    double total_ms = 0.0;
    std::size_t edge_state_bytes = 0;
    std::size_t square_state_bytes = 0;
    bool h1_computed = false;
};

class Computer {
public:
    Computer(const double* values,
             std::size_t depth,
             std::size_t rows,
             std::size_t cols,
             bool include_apparent_zero_pairs,
             const std::uint32_t* value_codes,
             bool compute_h1,
             bool profile_h1,
             const std::uint32_t* zero_lookup_masks)
        : values_(values),
          value_codes_(value_codes),
          depth_(depth),
          rows_(rows),
          cols_(cols),
          plane_(rows * cols),
          n_vertices_(checked_vertex_count(depth, rows, cols)),
          n_cubes_(0),
          infinity_(0),
          include_apparent_zero_pairs_(include_apparent_zero_pairs),
          compute_h1_(compute_h1),
          profile_h1_(profile_h1),
          edge_state_(nullptr),
          square_state_(nullptr),
          edge_birth_vertex_(nullptr),
          square_birth_edge_(nullptr),
          square_birth_vertex_(nullptr),
          square_second_rank_(nullptr),
          cube_birth_square_(nullptr),
          cube_birth_vertex_(nullptr),
          cube_second_rank_(nullptr),
          zero_lookup_masks_(zero_lookup_masks),
          zero_lookup_stats_(nullptr) {
        if (values_ == nullptr) {
            throw std::invalid_argument("3D image value pointer is null");
        }
        if (zero_lookup_masks_ == nullptr) {
            throw std::invalid_argument("3D lookup masks are required");
        }
        if (depth_ == 0 || rows_ == 0 || cols_ == 0) {
            throw std::invalid_argument("3D image dimensions must be positive");
        }
        if (n_vertices_ > std::numeric_limits<VertexCode>::max()) {
            throw std::overflow_error("3D image has too many vertices for packed ids");
        }
        if (depth_ >= 2 && rows_ >= 2 && cols_ >= 2) {
            const std::size_t d = depth_ - 1;
            const std::size_t r = rows_ - 1;
            const std::size_t c = cols_ - 1;
            if (d > std::numeric_limits<std::size_t>::max() / r) {
                throw std::overflow_error("3D cube count overflows size_t");
            }
            const std::size_t dr = d * r;
            if (dr > std::numeric_limits<std::size_t>::max() / c) {
                throw std::overflow_error("3D cube count overflows size_t");
            }
            n_cubes_ = dr * c;
            if (n_cubes_ >= std::numeric_limits<CubeCode>::max()) {
                throw std::overflow_error("3D image has too many cubes for packed ids");
            }
            infinity_ = n_cubes_;
        }
    }

    FullResult run() {
        FullResult result;
        const auto total_start = Clock::now();

        std::vector<std::uint8_t> edge_state(edge_state_storage_size(depth_, rows_, cols_), 0);
        std::vector<std::uint8_t> square_state(square_state_storage_size(depth_, rows_, cols_), 0);
        edge_state_ = edge_state.data();
        square_state_ = square_state.data();
        zero_lookup_stats_ = &result.zero_lookup;

        result.apparent.primal_h1_ms = 0.0;
        result.apparent.dual_h1_ms = 0.0;

        result.h0 = run_h0();
        result.h2 = run_h2();
        if (compute_h1_) {
            result.h1 = run_h1();
            result.h1_computed = true;
        }
        result.edge_state_bytes = edge_state.size() * sizeof(std::uint8_t);
        result.square_state_bytes = square_state.size() * sizeof(std::uint8_t);
        result.total_ms = elapsed_ms(total_start, Clock::now());
        return result;
    }

private:

#include "algorithms/detail/smart_3d_common_methods.inc"
#include "algorithms/detail/smart_3d_h0_methods.inc"
#include "algorithms/detail/smart_3d_h2_methods.inc"
#include "algorithms/detail/smart_3d_h1_methods.inc"

    static constexpr std::uint8_t edge_x_neg_bit = 1u << 0;
    static constexpr std::uint8_t edge_x_pos_bit = 1u << 1;
    static constexpr std::uint8_t edge_y_neg_bit = 1u << 2;
    static constexpr std::uint8_t edge_y_pos_bit = 1u << 3;
    static constexpr std::uint8_t edge_z_neg_bit = 1u << 4;
    static constexpr std::uint8_t edge_z_pos_bit = 1u << 5;

    static constexpr std::uint8_t face_x_neg_bit = 1u << 0;
    static constexpr std::uint8_t face_x_pos_bit = 1u << 1;
    static constexpr std::uint8_t face_y_neg_bit = 1u << 2;
    static constexpr std::uint8_t face_y_pos_bit = 1u << 3;
    static constexpr std::uint8_t face_z_neg_bit = 1u << 4;
    static constexpr std::uint8_t face_z_pos_bit = 1u << 5;
};

inline FullResult compute(const double* values,
                          std::size_t depth,
                          std::size_t rows,
                          std::size_t cols,
                          bool include_apparent_zero_pairs,
                          const std::uint32_t* value_codes,
                          bool compute_h1,
                          bool profile_h1,
                          const std::uint32_t* zero_lookup_masks) {
    return Computer(values, depth, rows, cols, include_apparent_zero_pairs,
                    value_codes, compute_h1, profile_h1, zero_lookup_masks)
        .run();
}

} // namespace smart_h0_h2_3d
