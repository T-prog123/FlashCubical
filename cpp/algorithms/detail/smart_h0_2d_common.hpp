#pragma once

#include "core/states.hpp"
#include "core/types.hpp"
#include "core/union_find.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smart_h0 {

using smart_core::MemoryProfile;
using smart_core::PersistencePair;
using smart_core::TimingProfile;
using smart_core::edge_state_h0_apparent;
using smart_core::edge_state_h0_negative;
using smart_core::edge_state_h1_apparent;
using smart_core::edge_state_h1_positive;

struct Result {
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

struct H1Result {
    std::vector<PersistencePair> pairs;
    std::size_t primal_edges = 0;
    std::size_t dual_vertices = 0;
    std::size_t active_edge_owner_vertices = 0;
    std::size_t apparent_edges = 0;
    std::size_t skipped_known_h0_negative_edges = 0;
    std::size_t uf_merge_attempts = 0;
    std::size_t uf_successful_merges = 0;
    std::size_t uf_same_component_attempts = 0;
    TimingProfile timing;
    MemoryProfile memory;
};

struct H1StagedRun {
    H1Result result;
    std::unique_ptr<std::uint8_t[]> apparent_edge_mask;
    bool apparent_edges_applied = false;
};

using EdgeCode = std::uint64_t;
using VertexCode = std::uint32_t;

inline std::size_t edge_state_storage_size(std::size_t rows, std::size_t cols) {
    if (rows == 0 || cols == 0) {
        return 0;
    }
    if (rows > (std::numeric_limits<std::size_t>::max() / cols)) {
        throw std::overflow_error("image dimensions overflow size_t");
    }
    return 2 * rows * cols;
}

} // namespace smart_h0
