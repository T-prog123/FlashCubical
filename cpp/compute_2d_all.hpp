#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace cubicalp_native {

std::vector<std::array<double, 3>>
compute_2d_all(const double* values, std::size_t rows, std::size_t cols);

}
