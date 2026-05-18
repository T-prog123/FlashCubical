#include "compute_2d_all.hpp"

#include "algorithms/smart_h0_2d.hpp"
#include "zero_lookup/runtime/mask_build.hpp"
#include "zero_lookup/runtime/zero_table.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace cubicalp_native {

std::vector<std::array<double, 3>>
compute_2d_all(const double* values, std::size_t rows, std::size_t cols) {
    const std::size_t n_vertices = rows * cols;

    // Shared edge-state array: H1 prepass marks H1-positive edges so H0 can
    // skip them; H0 marks H0-apparent edges so H1 finish can skip them.
    std::vector<std::uint8_t> edge_state(
        smart_h0::edge_state_storage_size(rows, cols), 0);

    // Build per-vertex older-neighbour masks (fast scan over values array).
    zero_lookup::MaskBuild2D masks =
        zero_lookup::build_masks2d(values, rows, cols, nullptr);

    // --- Phase 1: H1 apparent prepass (lookup-driven) ---
    zero_lookup::RunStats stats;
    smart_h0::H1DualComputer h1_computer(
        values, rows, cols, edge_state.data(), nullptr);
    smart_h0::H1StagedRun h1_staged =
        h1_computer.start_zero_lookup_pass(masks.masks.data(), &stats);

    // Buffers produced by H0 and consumed by H1 finish for cross-pruning.
    std::unique_ptr<std::uint8_t[]> h1_candidate_edge_owner_mask(
        new std::uint8_t[n_vertices]());
    std::vector<smart_h0::VertexCode> h1_candidate_active_vertices;
    h1_candidate_active_vertices.reserve(n_vertices / 8 + 1);

    // --- Phase 2: H0 primal union-find ---
    smart_h0::Result h0 = smart_h0::compute(
        values, rows, cols,
        /*include_apparent_zero_pairs=*/false,
        edge_state.data(),
        /*skip_known_h1_positive=*/true,
        h1_staged.apparent_edge_mask.get(),
        /*value_codes=*/nullptr,
        h1_candidate_edge_owner_mask.get(),
        &h1_candidate_active_vertices,
        masks.masks.data(),
        &stats);

    // --- Phase 3: H1 dual union-find (using H0 cross-pruning output) ---
    smart_h0::H1Result h1 = h1_computer.finish_after_h0(
        std::move(h1_staged),
        h1_candidate_edge_owner_mask.get(),
        &h1_candidate_active_vertices);

    std::vector<std::array<double, 3>> pairs;
    pairs.reserve(h0.pairs.size() + h1.pairs.size());

    for (const auto& p : h0.pairs) {
        pairs.push_back({p.birth, p.death, 0.0});
    }
    for (const auto& p : h1.pairs) {
        pairs.push_back({p.birth, p.death, 1.0});
    }

    return pairs;
}

} // namespace cubicalp_native
