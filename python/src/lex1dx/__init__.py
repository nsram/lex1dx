"""lex1dx: fast 1-D k-means with a 0.01% relative-excess-SSE tolerance.

    import numpy as np, lex1dx
    r = lex1dx.cluster(x, 64)
    s = lex1dx.cluster(x, 64, mode="sig5")

``cluster`` sorts the data, optionally groups to 5 significant digits, runs
tuned batch Collapse to ``M = default_M(k) = max(4096, 32 * k)`` (or an explicit
``M``), a weighted wilber-interp DP, and Refine. Errors (empty or non-finite
input, k < 1, or a bad mode) raise ``ValueError``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import _core

__all__ = ["Result", "cluster", "default_M"]
__version__ = "0.1.0"


def default_M(k: int) -> int:
    """The default Collapse target for k bins: max(4096, 32 * k)."""
    if int(k) < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    return int(_core.default_M(int(k)))


@dataclass(frozen=True)
class Result:
    """A fitted partition of the sorted data. Bin i holds sorted values bounds[i]:bounds[i+1]."""

    bounds: np.ndarray   # k + 1 cut indices into the sorted data (int64)
    sizes: np.ndarray    # values per bin (int64)
    centers: np.ndarray  # bin means of the raw values
    upper: np.ndarray    # largest raw value in each bin
    sse: float           # sum of squared deviations from the bin means (raw values)
    D: int               # distinct values after sort and -0 canonicalization
    Dprime: int          # atoms after the grouping actually used (equals D for tuned)
    sig: int = 0         # significant digits used for grouping; 0 = full precision
    sweeps: int
    moves: int
    presorted: bool      # the input was already sorted, so the sort was skipped

    @property
    def k(self) -> int:
        """Bins used: min(requested k, atoms)."""
        return int(self.sizes.shape[0])

    def assign(self, x) -> np.ndarray:
        """Bin index for each value: the first bin whose upper value is >= x (last bin if none)."""
        idx = np.searchsorted(self.upper, np.asarray(x, dtype=np.float64), side="left")
        return np.minimum(idx, self.k - 1)


def _vec(a, name: str) -> np.ndarray:
    v = np.ascontiguousarray(a, dtype=np.float64)
    if v.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {v.shape}")
    return v


def cluster(x, k: int, *, mode: str = "tuned", M: Optional[int] = None,
            tolerance: float = 1e-11) -> Result:
    """Cluster ``x`` into ``k`` contiguous bins of the sorted values.

    ``mode`` is ``"tuned"`` (full-precision atoms) or ``"sig5"`` (group to 5
    significant digits first). ``M=None`` uses ``default_M(k)``.
    """
    if mode not in ("tuned", "sig5"):
        raise ValueError("mode must be 'tuned' or 'sig5'")
    if M is not None and int(M) < 0:
        raise ValueError(f"M must be >= 0, got {M}")
    m = 0 if M is None else int(M)
    return Result(**_core.cluster(_vec(x, "x"), int(k), mode, m, float(tolerance)))
