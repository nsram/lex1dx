"""Quantized bin table and histogram. matplotlib is imported only when drawing."""
from __future__ import annotations

from typing import Optional

import numpy as np

_THEMES = {
    "default": {"bg": "#005fff", "fg": "#ffff00", "bar": "#ffff00"},
    "light": {"bg": "#ffffff", "fg": "#000000", "bar": "#555555"},
}
_MAX_BLOCKS = 60
_LABEL_WIDTH = 26  # matches the 2024 Scala audio layout (indent + lo - hi)


def _sig5(v: float) -> str:
    return f"{v:.5g}"


def _fmt_pct(pct: float) -> str:
    s = f"{pct:.5f}".rstrip("0").rstrip(".")
    return f"{s}%"


def _fmt_n(n) -> str:
    if isinstance(n, (int, np.integer)) or float(n) == int(n):
        return str(int(n))
    return f"{n:g}"


def bin_table(result, x, weights=None):
    """Rows ``(lo, hi, n, pct, mean)`` for each bin of ``result`` on ``x``.

    ``lo`` and ``hi`` are the minimum and maximum raw values in the bin.
    ``n`` is the count, or the weight mass when ``weights`` is given.
    ``pct`` is 100 * n / N (or 100 * mass / total mass).
    ``mean`` is the (weighted) mean of the values in the bin.
    """
    xv = np.ascontiguousarray(x, dtype=np.float64)
    if xv.ndim != 1:
        raise ValueError(f"x must be one-dimensional, got shape {xv.shape}")
    wv = None
    if weights is not None:
        wv = np.ascontiguousarray(weights, dtype=np.float64)
        if wv.ndim != 1:
            raise ValueError(f"weights must be one-dimensional, got shape {wv.shape}")
        if wv.shape[0] != xv.shape[0]:
            raise ValueError("weights size must match the input length")

    order = np.argsort(xv, kind="mergesort")
    sx = xv[order]
    sw = None if wv is None else wv[order]
    bounds = np.asarray(result.bounds)
    if int(bounds[0]) != 0 or int(bounds[-1]) != xv.shape[0]:
        raise ValueError(
            f"result.bounds must cover x: got {bounds.tolist()} for length {xv.shape[0]}"
        )
    total = float(xv.shape[0] if sw is None else sw.sum())
    rows = []
    for i in range(result.k):
        a = int(bounds[i])
        b = int(bounds[i + 1])
        lo = float(sx[a])
        hi = float(sx[b - 1])
        if sw is None:
            n: int | float = int(b - a)
            mean = float(sx[a:b].mean()) if b > a else 0.0
        else:
            mass = float(sw[a:b].sum())
            n = mass
            mean = float(np.dot(sw[a:b], sx[a:b]) / mass) if mass else 0.0
        pct = 100.0 * (n / total) if total else 0.0
        rows.append((lo, hi, n, pct, mean))
    return rows


def _r2_mse(result, x, weights=None):
    xv = np.ascontiguousarray(x, dtype=np.float64)
    if weights is None:
        mu = float(xv.mean())
        sst = float(np.sum((xv - mu) ** 2))
        denom = float(xv.shape[0])
    else:
        wv = np.ascontiguousarray(weights, dtype=np.float64)
        wsum = float(wv.sum())
        mu = float(np.dot(wv, xv) / wsum) if wsum else 0.0
        sst = float(np.sum(wv * (xv - mu) ** 2))
        denom = wsum
    sse = float(result.sse)
    r2 = 1.0 if sst == 0 else 1.0 - sse / sst
    mse = sse / denom if denom else 0.0
    return r2, mse


def _lines(result, x, weights=None, title: str = ""):
    rows = bin_table(result, x, weights=weights)
    r2, mse = _r2_mse(result, x, weights=weights)
    k = result.k
    lines = []
    if title:
        lines.append(str(title))
    lines.append(f"Squared correlation r² of {k}-bin model with data = {r2:.8f}")
    lines.append(f"MSE of {k}-bin model with data = {mse:.8f}")
    lines.append("")
    n_max = max((row[2] for row in rows), default=0)
    for lo, hi, n, pct, mean in rows:
        label = f"{_sig5(lo)} - {_sig5(hi)}"
        left = ("    " + label).ljust(_LABEL_WIDTH)
        if n_max > 0 and n > 0:
            blocks = max(1, int(round(_MAX_BLOCKS * float(n) / float(n_max))))
        else:
            blocks = 0
        annot = f"[{_fmt_n(n)}, {_fmt_pct(pct)}, {_sig5(mean)}]"
        lines.append(left + ("█" * blocks) + (" " if blocks else "") + annot)
        lines.append("")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def histogram(
    result,
    x,
    file: Optional[str] = "lex1dx_histogram.png",
    ax=None,
    weights=None,
    title: str = "",
    theme: str = "default",
):
    """Draw a horizontal bin table for ``result`` on ``x``.

    One row per bin, top to bottom from the lowest bin. Each row has the
    bin's min and max (5 significant digits), a bar whose length is
    proportional to n (or weight mass), and ``[n, pct%, mean]``.
    The header shows ``title``, then r² = 1 - SSE/SST and MSE = SSE/N.

    ``theme="default"`` is a blue background with yellow bars and text,
    after the 2024 Scala audio figures. ``theme="light"`` is a white
    background with dark grey bars and black text.

    By default the figure is written to ``file`` (``lex1dx_histogram.png``)
    and that path is returned. Pass ``file=None`` to skip writing and
    return the Axes. If ``ax`` is given, the plot is drawn on it.

    matplotlib is imported only inside this function. Install it with
    ``pip install lex1dx[plot]``.
    """
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from matplotlib.font_manager import FontProperties
    except ImportError as e:
        raise ImportError(
            "histogram() requires matplotlib. Install it with: pip install lex1dx[plot]"
        ) from e

    if theme not in _THEMES:
        raise ValueError(f"theme must be 'default' or 'light', got {theme!r}")
    colors = _THEMES[theme]
    lines = _lines(result, x, weights=weights, title=title)
    nlines = max(len(lines), 1)
    width_ch = max((len(L) for L in lines), default=40) + 4
    # DejaVu Sans Mono is ~0.6 em wide; 11pt, 12in canvas fits the 60-block rows.
    fig_w = max(10.0, min(14.0, 0.095 * width_ch))
    fig_h = 0.40 + 0.22 * nlines

    if ax is None:
        fig = Figure(figsize=(fig_w, fig_h), facecolor=colors["bg"])
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        owned = True
    else:
        fig = ax.figure
        owned = False

    ax.set_facecolor(colors["bg"])
    ax.set_xlim(0, width_ch)
    ax.set_ylim(nlines + 0.3, -0.3)
    ax.set_axis_off()
    for s in ax.spines.values():
        s.set_visible(False)
    font = FontProperties(family="DejaVu Sans Mono", size=10)
    for i, line in enumerate(lines):
        # Paint the bar segment in the bar color (same as text in default theme).
        ax.text(
            0,
            i,
            line,
            fontproperties=font,
            color=colors["fg"],
            ha="left",
            va="top",
            clip_on=False,
        )
    if owned:
        fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    else:
        fig.tight_layout()

    if file is not None:
        fig.savefig(file, dpi=120, facecolor=colors["bg"])
        return file
    return ax
