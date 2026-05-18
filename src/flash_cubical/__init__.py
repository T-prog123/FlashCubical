"""flash_cubical – fast cubical persistent homology for 2D and 3D scalar grids."""

import warnings
from importlib import resources

import numpy as np

from ._core import compute as _compute_native
from ._core import initialize_lookup_tables as _initialize_lookup_tables
from .persistence import Persistence

__all__ = ["compute", "Persistence"]


def _load_packaged_lookup_tables() -> None:
    data = resources.files(__package__).joinpath("data")
    path2d = data.joinpath("zero_lookup_2d.bin")
    path3d = data.joinpath("zero_lookup_3d.bin")
    with resources.as_file(path2d) as p2d, resources.as_file(path3d) as p3d:
        _initialize_lookup_tables(str(p2d), str(p3d))


_load_packaged_lookup_tables()


def compute(x, *, min_persistence: float = 0.0, h1: bool = True) -> Persistence:
    """Compute cubical persistent homology of a 2D or 3D scalar array.

    Parameters
    ----------
    x : array-like
        Input scalar field.  Converted to a C-contiguous float64 ndarray.
        Must be 2D (rows × cols) or 3D (depth × rows × cols).
    min_persistence : float, optional
        Minimum persistence threshold.  Currently not implemented by the
        backend; any value other than 0.0 emits a RuntimeWarning and is
        treated as 0.0.
    h1 : bool, optional
        Whether to compute H1 (loops / tunnels).
        - For 3D inputs: defaults to True; pass False to skip H1.
        - For 2D inputs: must be True (raising ValueError otherwise), because
          the 2D pipeline always computes H0 and H1 together.

    Returns
    -------
    Persistence
        Object with a ``.values`` attribute (float64 array of shape (n, 3),
        columns: birth, death, homology dimension) and a ``.plot()`` method.

    Raises
    ------
    ValueError
        If ``x`` is not 2D or 3D, if ``min_persistence < 0``, or if
        ``h1=False`` is requested for a 2D input.
    RuntimeWarning
        If ``min_persistence != 0.0`` (backend does not yet support this).
    """
    x = np.ascontiguousarray(np.asarray(x, dtype=np.float64))

    if x.ndim < 2 or x.ndim > 3:
        raise ValueError(
            f"Input must be 2D or 3D, got {x.ndim}D array."
        )

    if min_persistence < 0.0:
        raise ValueError(
            f"min_persistence must be non-negative, got {min_persistence}."
        )

    if min_persistence != 0.0:
        warnings.warn(
            "min_persistence != 0.0 is not yet supported by the backend; "
            "running as if min_persistence=0.0.",
            RuntimeWarning,
            stacklevel=2,
        )

    if x.ndim == 2 and not h1:
        raise ValueError(
            "h1=False is only valid for 3D inputs.  The 2D backend always "
            "computes H0 and H1 together."
        )

    arr = _compute_native(x, h1)
    return Persistence(arr)
