import numpy as np

from ._core import compute as _compute_native
from .persistence import Persistence

__all__ = ["compute", "Persistence"]


def compute(x, *, min_persistence: float = 0.0, h1: bool = True) -> Persistence:
    x = np.ascontiguousarray(np.asarray(x, dtype=np.float64))

    if x.ndim < 2 or x.ndim > 3:
        raise ValueError(
            f"Input must be 2D or 3D, got {x.ndim}D array."
        )

    if min_persistence < 0.0:
        raise ValueError(
            f"min_persistence must be non-negative, got {min_persistence}."
        )

    if x.ndim == 2 and not h1:
        raise ValueError(
            "h1=False is only valid for 3D inputs.  The 2D backend always "
            "computes H0 and H1 together."
        )

    arr = _compute_native(x, h1)
    if min_persistence > 0.0:
        finite = np.isfinite(arr[:, 1])
        keep = finite & ((arr[:, 1] - arr[:, 0]) >= min_persistence)
        arr = arr[keep | ~finite]
    return Persistence(arr)
