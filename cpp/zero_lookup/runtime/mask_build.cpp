#include "zero_lookup/runtime/mask_build.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace zero_lookup {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

bool vertex_less(const double* values,
                 const std::uint32_t* value_codes,
                 std::size_t a,
                 std::size_t b) {
    if (value_codes != nullptr) {
        const std::uint32_t ca = value_codes[a];
        const std::uint32_t cb = value_codes[b];
        return ca < cb || (ca == cb && a < b);
    }
    const double va = values[a];
    const double vb = values[b];
    return va < vb || (va == vb && a < b);
}

void set_older_bit2d(std::vector<std::uint8_t>& masks,
                     const double* values,
                     const std::uint32_t* value_codes,
                     std::size_t a,
                     std::size_t b,
                     std::uint8_t bit_b_sees_a,
                     std::uint8_t bit_a_sees_b) {
    if (vertex_less(values, value_codes, a, b)) {
        masks[b] = static_cast<std::uint8_t>(masks[b] | (1u << bit_b_sees_a));
    } else {
        masks[a] = static_cast<std::uint8_t>(masks[a] | (1u << bit_a_sees_b));
    }
}

std::uint32_t bit3d_from_offset(int dx, int dy, int dz) {
    if (dy == 0 && dz == 0) {
        return dx < 0 ? 0u : 1u;
    }
    if (dx == 0 && dz == 0) {
        return dy < 0 ? 2u : 3u;
    }
    if (dx == 0 && dy == 0) {
        return dz < 0 ? 4u : 5u;
    }
    if (dz == 0) {
        const std::uint32_t sx = dx > 0 ? 1u : 0u;
        const std::uint32_t sy = dy > 0 ? 1u : 0u;
        return 6u + sx + 2u * sy;
    }
    if (dy == 0) {
        const std::uint32_t sx = dx > 0 ? 1u : 0u;
        const std::uint32_t sz = dz > 0 ? 1u : 0u;
        return 10u + sx + 2u * sz;
    }
    if (dx == 0) {
        const std::uint32_t sy = dy > 0 ? 1u : 0u;
        const std::uint32_t sz = dz > 0 ? 1u : 0u;
        return 14u + sy + 2u * sz;
    }
    const std::uint32_t sx = dx > 0 ? 1u : 0u;
    const std::uint32_t sy = dy > 0 ? 1u : 0u;
    const std::uint32_t sz = dz > 0 ? 1u : 0u;
    return 18u + sx + 2u * sy + 4u * sz;
}

void set_older_bit3d(std::vector<std::uint32_t>& masks,
                     const double* values,
                     const std::uint32_t* value_codes,
                     std::size_t a,
                     std::size_t b,
                     std::uint32_t bit_b_sees_a,
                     std::uint32_t bit_a_sees_b) {
    if (vertex_less(values, value_codes, a, b)) {
        masks[b] |= (1u << bit_b_sees_a);
    } else {
        masks[a] |= (1u << bit_a_sees_b);
    }
}

void propagate_offset3d(std::vector<std::uint32_t>& masks,
                        const double* values,
                        const std::uint32_t* value_codes,
                        std::size_t depth,
                        std::size_t rows,
                        std::size_t cols,
                        int dx,
                        int dy,
                        int dz) {
    const std::uint32_t bit_b_sees_a = bit3d_from_offset(-dx, -dy, -dz);
    const std::uint32_t bit_a_sees_b = bit3d_from_offset(dx, dy, dz);

    const std::size_t z_begin = dz < 0 ? static_cast<std::size_t>(-dz) : 0;
    const std::size_t y_begin = dy < 0 ? static_cast<std::size_t>(-dy) : 0;
    const std::size_t x_begin = dx < 0 ? static_cast<std::size_t>(-dx) : 0;
    const std::size_t z_end =
        depth - (dz > 0 ? static_cast<std::size_t>(dz) : 0);
    const std::size_t y_end =
        rows - (dy > 0 ? static_cast<std::size_t>(dy) : 0);
    const std::size_t x_end =
        cols - (dx > 0 ? static_cast<std::size_t>(dx) : 0);

    const std::ptrdiff_t delta =
        static_cast<std::ptrdiff_t>(dz) *
            static_cast<std::ptrdiff_t>(rows * cols) +
        static_cast<std::ptrdiff_t>(dy) * static_cast<std::ptrdiff_t>(cols) +
        static_cast<std::ptrdiff_t>(dx);

    for (std::size_t z = z_begin; z < z_end; ++z) {
        for (std::size_t y = y_begin; y < y_end; ++y) {
            std::size_t a = (z * rows + y) * cols + x_begin;
            for (std::size_t x = x_begin; x < x_end; ++x, ++a) {
                const std::size_t b =
                    static_cast<std::size_t>(
                        static_cast<std::ptrdiff_t>(a) + delta);
                set_older_bit3d(masks, values, value_codes, a, b,
                                bit_b_sees_a, bit_a_sees_b);
            }
        }
    }
}

}

MaskBuild2D build_masks2d(const double* values,
                          std::size_t rows,
                          std::size_t cols,
                          const std::uint32_t* value_codes) {
    const auto start = Clock::now();
    const std::size_t n_vertices = rows * cols;
    MaskBuild2D out;
    out.masks.assign(n_vertices, 0);

    for (std::size_t r = 0; r < rows; ++r) {
        const std::size_t base = r * cols;
        for (std::size_t c = 0; c < cols; ++c) {
            const std::size_t v = base + c;
            if (c + 1 < cols) {
                set_older_bit2d(out.masks, values, value_codes, v, v + 1, 3, 1);
            }
            if (r + 1 < rows) {
                set_older_bit2d(out.masks, values, value_codes, v, v + cols, 0, 2);
                if (c + 1 < cols) {
                    set_older_bit2d(out.masks, values, value_codes, v, v + cols + 1,
                                    7, 5);
                }
                if (c > 0) {
                    set_older_bit2d(out.masks, values, value_codes, v, v + cols - 1,
                                    4, 6);
                }
            }
        }
    }
    out.build_ms = elapsed_ms(start, Clock::now());
    return out;
}

MaskBuild3D build_masks3d(const double* values,
                          std::size_t depth,
                          std::size_t rows,
                          std::size_t cols,
                          const std::uint32_t* value_codes) {
    const auto start = Clock::now();
    const std::size_t plane = rows * cols;
    const std::size_t n_vertices = depth * plane;
    MaskBuild3D out;
    out.masks.assign(n_vertices, 0);

    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       1, 0, 0);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       -1, 1, 0);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       0, 1, 0);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       1, 1, 0);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       -1, -1, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       0, -1, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       1, -1, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       -1, 0, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       0, 0, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       1, 0, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       -1, 1, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       0, 1, 1);
    propagate_offset3d(out.masks, values, value_codes, depth, rows, cols,
                       1, 1, 1);

    out.build_ms = elapsed_ms(start, Clock::now());
    return out;
}

}
