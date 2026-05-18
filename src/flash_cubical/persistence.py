import numpy as np


class Persistence:
    """Cubical persistence diagram.

    Attributes
    ----------
    values : ndarray of shape (n, 3), dtype float64
        Columns: birth, death, homology dimension.
        Essential classes (connected components) have death = +inf.
    """

    def __init__(self, values: np.ndarray) -> None:
        arr = np.asarray(values, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(
                f"Persistence values must have shape (n, 3), got {arr.shape}"
            )
        self._values = arr

    @property
    def values(self) -> np.ndarray:
        """float64 array of shape (n, 3): birth, death, dimension."""
        return self._values

    def __array__(self, dtype=None):
        if dtype is not None:
            return self._values.astype(dtype)
        return self._values

    def __repr__(self) -> str:
        n = len(self._values)
        dims = sorted(int(d) for d in np.unique(self._values[:, 2])) if n else []
        dim_str = ", ".join(f"H{d}" for d in dims)
        return f"Persistence({n} pairs [{dim_str}])"

    def plot(self, ax=None):
        """Plot the persistence diagram (birth vs death, coloured by dimension).

        Infinite-death pairs (essential classes) are omitted from the scatter.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; created if not provided.

        Returns
        -------
        ax : matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        _COLORS = {0: "tab:blue", 1: "tab:orange", 2: "tab:green"}

        if ax is None:
            _, ax = plt.subplots()

        v = self._values
        if len(v) == 0:
            ax.set_xlabel("Birth")
            ax.set_ylabel("Death")
            ax.set_title("Persistence Diagram")
            return ax

        finite_mask = np.isfinite(v[:, 1])
        finite = v[finite_mask]

        dims = sorted(int(d) for d in np.unique(v[:, 2]))
        for d in dims:
            sel = finite[:, 2] == d
            ax.scatter(
                finite[sel, 0],
                finite[sel, 1],
                color=_COLORS.get(d, f"C{d}"),
                label=f"H{d}",
                s=12,
                alpha=0.7,
                linewidths=0,
            )

        if len(finite) > 0:
            lo = min(finite[:, 0].min(), finite[:, 1].min())
            hi = max(finite[:, 0].max(), finite[:, 1].max())
            ax.plot([lo, hi], [lo, hi], "k--", linewidth=0.5, alpha=0.4)

        ax.set_xlabel("Birth")
        ax.set_ylabel("Death")
        ax.set_title("Persistence Diagram")
        ax.legend()
        return ax
