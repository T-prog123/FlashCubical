#pragma once

#include <cstddef>

namespace smart_core {

struct PersistencePair {
    double birth = 0.0;
    double death = 0.0;
};

struct TimingProfile {
    double total_ms = 0.0;
    double parent_initialization_ms = 0.0;
    double apparent_pair_scan_ms = 0.0;
    double apparent_pair_work_ms = 0.0;
    double sorting_ms = 0.0;
    double union_find_sweep_ms = 0.0;
    double finalization_ms = 0.0;
};

struct MemoryProfile {
    std::size_t parent_bytes = 0;
    std::size_t non_apparent_edges_bytes = 0;
    std::size_t non_apparent_vertex_mask_bytes = 0;
    std::size_t active_vertices_bytes = 0;
    std::size_t square_birth_bytes = 0;
    std::size_t persistence_pairs_bytes = 0;

    std::size_t total_bytes() const {
        return parent_bytes + non_apparent_edges_bytes +
               non_apparent_vertex_mask_bytes + active_vertices_bytes +
               square_birth_bytes + persistence_pairs_bytes;
    }
};

} // namespace smart_core
