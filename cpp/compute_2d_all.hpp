#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace cubicalp_native {

// Run the full 2D cubical persistence pipeline (H0 + H1).
// Returns one entry per pair: {birth, death, dim} where dim is 0 or 1.
// Infinite-death entries (essential H0 class) are included with death = +inf.
// Caller must have called zero_lookup::warm_zero_table2d() exactly once before
// the first call to this function.
std::vector<std::array<double, 3>>
compute_2d_all(const double* values, std::size_t rows, std::size_t cols);

} // namespace cubicalp_native
