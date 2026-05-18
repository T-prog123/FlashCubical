import numpy as np
import pytest

import flash_cubical as fc
from flash_cubical import Persistence


RNG = np.random.default_rng(42)


def make_2d(rows=8, cols=8):
    return RNG.random((rows, cols))


def make_3d(depth=4, rows=4, cols=4):
    return RNG.random((depth, rows, cols))


def test_1d_raises():
    with pytest.raises(ValueError, match="2D or 3D"):
        fc.compute(np.ones(10))


def test_4d_raises():
    with pytest.raises(ValueError, match="2D or 3D"):
        fc.compute(np.ones((2, 2, 2, 2)))


def test_2d_h1_false_raises():
    with pytest.raises(ValueError, match="h1=False"):
        fc.compute(make_2d(), h1=False)


def test_3d_h1_false_runs():
    ph = fc.compute(make_3d(), h1=False)
    assert ph.values.shape[1] == 3
    dims = set(ph.values[:, 2].astype(int).tolist())
    assert 1 not in dims


def test_negative_min_persistence_raises():
    with pytest.raises(ValueError, match="non-negative"):
        fc.compute(make_2d(), min_persistence=-1.0)


def test_nonzero_min_persistence_filters():
    ph_all = fc.compute(make_2d())
    ph_filtered = fc.compute(make_2d(), min_persistence=0.1)
    finite_all = np.sum(np.isfinite(ph_all.values[:, 1]))
    finite_filtered = np.sum(np.isfinite(ph_filtered.values[:, 1]))
    assert finite_filtered <= finite_all
    assert np.sum(np.isinf(ph_filtered.values[:, 1])) == 1


def test_2d_runs():
    ph = fc.compute(make_2d())
    assert isinstance(ph, Persistence)


def test_2d_values_shape():
    ph = fc.compute(make_2d())
    assert ph.values.ndim == 2
    assert ph.values.shape[1] == 3


def test_2d_values_dtype():
    ph = fc.compute(make_2d())
    assert ph.values.dtype == np.float64


def test_2d_dims_are_0_and_1():
    ph = fc.compute(make_2d(16, 16))
    dims = set(ph.values[:, 2].astype(int).tolist())
    assert dims <= {0, 1}


def test_2d_has_essential_h0():
    ph = fc.compute(make_2d())
    h0_pairs = ph.values[ph.values[:, 2] == 0]
    essential = h0_pairs[np.isinf(h0_pairs[:, 1])]
    assert len(essential) == 1


def test_3d_runs():
    ph = fc.compute(make_3d())
    assert isinstance(ph, Persistence)


def test_3d_values_shape():
    ph = fc.compute(make_3d())
    assert ph.values.ndim == 2
    assert ph.values.shape[1] == 3


def test_3d_has_essential_h0():
    ph = fc.compute(make_3d())
    h0_pairs = ph.values[ph.values[:, 2] == 0]
    essential = h0_pairs[np.isinf(h0_pairs[:, 1])]
    assert len(essential) == 1


def test_3d_dims_are_subset_of_0_1_2():
    ph = fc.compute(make_3d())
    dims = set(ph.values[:, 2].astype(int).tolist())
    assert dims <= {0, 1, 2}


def test_values_attribute():
    ph = fc.compute(make_2d())
    assert hasattr(ph, "values")
    assert isinstance(ph.values, np.ndarray)


def test_np_asarray_equals_values():
    ph = fc.compute(make_2d())
    np.testing.assert_array_equal(np.asarray(ph), ph.values)


def test_plot_returns_axes():
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    ph = fc.compute(make_2d())
    ax = ph.plot()
    assert hasattr(ax, "scatter")
    plt.close("all")


def test_float32_input_accepted():
    x = RNG.random((6, 6)).astype(np.float32)
    ph = fc.compute(x)
    assert ph.values.dtype == np.float64


def test_fortran_order_accepted():
    x = np.asfortranarray(RNG.random((6, 6)))
    ph = fc.compute(x)
    assert ph.values.shape[1] == 3


def test_list_input_accepted():
    x = [[float(i + j) / 10 for j in range(4)] for i in range(4)]
    ph = fc.compute(x)
    assert isinstance(ph, Persistence)
