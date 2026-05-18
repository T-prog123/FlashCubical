#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace zero_lookup {

struct MaskBuild2D {
    std::vector<std::uint8_t> masks;
    std::vector<std::uint32_t> vertex_rank;
    double build_ms = 0.0;
};

struct MaskBuild3D {
    std::vector<std::uint32_t> masks;
    std::vector<std::uint32_t> vertex_rank;
    double build_ms = 0.0;
};

MaskBuild2D build_masks2d(const double* values,
                          std::size_t rows,
                          std::size_t cols,
                          const std::uint32_t* value_codes = nullptr);

MaskBuild3D build_masks3d(const double* values,
                          std::size_t depth,
                          std::size_t rows,
                          std::size_t cols,
                          const std::uint32_t* value_codes = nullptr);

} // namespace zero_lookup
