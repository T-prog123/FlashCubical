#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "compute_2d_all.hpp"
#include "algorithms/smart_h0_h2_3d.hpp"
#include "zero_lookup/runtime/mask_build.hpp"
#include "zero_lookup/runtime/zero_table.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <mutex>
#include <string>
#include <vector>

namespace py = pybind11;

static std::once_flag g_init_2d;
static std::once_flag g_init_3d;

static py::array_t<double> _compute(
    py::array_t<double, py::array::c_style | py::array::forcecast> x,
    bool h1)
{
    const auto buf = x.request();
    const int ndim = static_cast<int>(buf.ndim);

    if (ndim != 2 && ndim != 3) {
        throw py::value_error(
            "cubicalp: input must be 2D or 3D, got " +
            std::to_string(ndim) + "D");
    }
    if (ndim == 2 && !h1) {
        throw py::value_error(
            "cubicalp: h1=False is only valid for 3D inputs; "
            "2D inputs always compute H0 and H1");
    }

    const double* values = static_cast<const double*>(buf.ptr);
    std::vector<std::array<double, 3>> pairs;

    if (ndim == 2) {
        std::call_once(g_init_2d, []() { zero_lookup::warm_zero_table2d(); });

        const std::size_t rows = static_cast<std::size_t>(buf.shape[0]);
        const std::size_t cols = static_cast<std::size_t>(buf.shape[1]);
        pairs = cubicalp_native::compute_2d_all(values, rows, cols);

    } else {
        std::call_once(g_init_3d, []() { zero_lookup::warm_zero_table3d(); });

        const std::size_t depth = static_cast<std::size_t>(buf.shape[0]);
        const std::size_t rows  = static_cast<std::size_t>(buf.shape[1]);
        const std::size_t cols  = static_cast<std::size_t>(buf.shape[2]);

        zero_lookup::MaskBuild3D masks =
            zero_lookup::build_masks3d(values, depth, rows, cols, nullptr);

        smart_h0_h2_3d::FullResult result = smart_h0_h2_3d::compute(
            values, depth, rows, cols,
            false,
            nullptr,
            h1,
            false,
            masks.masks.data());

        pairs.reserve(
            result.h0.pairs.size() +
            result.h1.pairs.size() +
            result.h2.pairs.size());

        for (const auto& p : result.h0.pairs) {
            pairs.push_back({p.birth, p.death, 0.0});
        }
        {
            const std::size_t n_verts = depth * rows * cols;
            const double birth = *std::min_element(values, values + n_verts);
            pairs.push_back({birth,
                             std::numeric_limits<double>::infinity(),
                             0.0});
        }
        if (result.h1_computed) {
            for (const auto& p : result.h1.pairs) {
                pairs.push_back({p.birth, p.death, 1.0});
            }
        }
        for (const auto& p : result.h2.pairs) {
            pairs.push_back({p.birth, p.death, 2.0});
        }
    }

    const std::size_t n = pairs.size();
    py::array_t<double> out({n, std::size_t{3}});
    auto r = out.mutable_unchecked<2>();
    for (std::size_t i = 0; i < n; ++i) {
        r(i, 0) = pairs[i][0];
        r(i, 1) = pairs[i][1];
        r(i, 2) = pairs[i][2];
    }
    return out;
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Cubical persistent homology – C++ backend";
    m.def("compute", &_compute,
          py::arg("x"), py::arg("h1"),
          "Compute cubical persistence pairs.\n\n"
          "Parameters\n"
          "----------\n"
          "x : ndarray, float64, C-contiguous, 2D or 3D\n"
          "h1 : bool – compute H1 (for 3D; always True for 2D)\n\n"
          "Returns\n"
          "-------\n"
          "ndarray of shape (n, 3), columns: birth, death, dim\n"
          "Infinite death values (+inf) represent essential classes.");
}
