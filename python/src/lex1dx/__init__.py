"""lex1dx: fast 1-D k-means with a 0.01% relative-excess-SSE tolerance.

    import numpy as np, lex1dx
    r = lex1dx.cluster(x, 64)
    s = lex1dx.cluster(x, 64, sig=5)
    e = lex1dx.cluster(x, 64, mode="exact", solver="ckmeans")

``cluster`` sorts the data, optionally groups to ``sig`` significant digits
(0 is full precision), runs tuned batch Collapse to
``M = default_M(k) = max(4096, 32 * k)`` (or an explicit ``M``), a DP on those
groups (wilber-interp by default, or the Fisher / ckmeans DP), and Refine.
``mode="exact"`` skips Collapse. ``mode="sig5"`` is an alias for ``sig=5``.
Errors (empty or non-finite input, k < 1, a bad mode, solver, sig, or weights)
raise ``ValueError``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import _core
from . import datasets

__all__ = ["Result", "cluster", "default_M", "datasets", "histogram", "bin_table"]
__version__ = "0.1.2"


def default_M(k: int) -> int:
    """The default Collapse target for k bins: max(4096, 32 * k)."""
    if int(k) < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    return int(_core.default_M(int(k)))


@dataclass(frozen=True)
class Result:
    """A fitted partition of the sorted data. Bin i holds sorted values bounds[i]:bounds[i+1]."""

    bounds: np.ndarray   # k + 1 cut indices into the sorted data (int64)
    sizes: np.ndarray    # values per bin (int64 element counts)
    centers: np.ndarray  # bin means of the raw values (weighted if weights were given)
    upper: np.ndarray    # largest raw value in each bin
    sse: float           # sum of squared deviations from the bin means (weighted if weights were given)
    D: int               # distinct values after sort and -0 canonicalization
    Dprime: int          # atoms after the grouping actually used (equals D when sig is 0)
    sig: int             # significant digits used for grouping; 0 = full precision
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


def _sig_arg(sig) -> int:
    if sig is None:
        return -1
    if isinstance(sig, bool) or not isinstance(sig, (int, np.integer)):
        raise ValueError(f"sig must be an integer in 0..9, got {sig!r}. Use sig=0 for full precision")
    s = int(sig)
    if s < 0 or s > 9:
        raise ValueError(f"sig must be an integer in 0..9, got {s}. Use sig=0 for full precision")
    return s


def cluster(x, k: int, *, mode: str = "tuned", sig: Optional[int] = None,
            M: Optional[int] = None, tolerance: float = 1e-11,
            weights=None, solver: str = "wilber") -> Result:
    """Cluster ``x`` into ``k`` contiguous bins of the sorted values.

    ``mode`` is ``"tuned"`` (the default pipeline), ``"sig5"`` (an alias for
    ``sig=5``), or ``"exact"`` (no Collapse; the chosen solver runs on the
    atoms). ``sig`` is an integer 0..9. 0 means full precision. 1..9 groups
    to that many significant digits. If that grouping yields fewer than
    ``k`` groups, precision is raised (up to 15, then raw atoms) so ``k``
    is kept when the raw data has at least ``k`` distinct values.
    ``Result.sig`` is the sig actually used (0 if the fallback reached
    raw atoms). ``None`` follows ``mode``. ``M=None``
    uses ``default_M(k)``. ``solver`` is ``"wilber"`` (default) or
    ``"ckmeans"``. Both solvers are exact and give the same optimum.
    Wilber is faster as k grows, because wilber-interp is O(k N) in
    style while ckmeans is a divide-and-conquer O(k N log N) DP.
    ckmeans is the original LEX exact DP and serves as a
    cross-check. ``weights`` is an optional non-negative finite vector of
    the same length as ``x``.
    """
    if mode not in ("tuned", "sig5", "exact"):
        raise ValueError("mode must be 'tuned', 'sig5', or 'exact'")
    if solver not in ("wilber", "ckmeans"):
        raise ValueError("solver must be 'wilber' or 'ckmeans'")
    if M is not None and int(M) < 0:
        raise ValueError(f"M must be >= 0, got {M}")
    xv = _vec(x, "x")
    wv = None
    if weights is not None:
        wv = _vec(weights, "weights")
        if wv.shape[0] != xv.shape[0]:
            raise ValueError("weights size must match input")
    m = 0 if M is None else int(M)
    return Result(**_core.cluster(xv, int(k), mode, m, float(tolerance),
                                  _sig_arg(sig), wv, solver))

from .histogram import histogram, bin_table
