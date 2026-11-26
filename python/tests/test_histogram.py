from pathlib import Path

import numpy as np
import pytest

import lex1dx

matplotlib = pytest.importorskip("matplotlib")

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def test_bin_table_matches_sizes_and_centers():
    x = lex1dx.datasets.load("athletics20")
    r = lex1dx.cluster(x, 5)
    rows = lex1dx.bin_table(r, x)
    assert len(rows) == r.k
    n = np.array([row[2] for row in rows])
    means = np.array([row[4] for row in rows])
    assert n.sum() == len(x)
    np.testing.assert_array_equal(n, r.sizes)
    np.testing.assert_allclose(means, r.centers, rtol=0, atol=1e-12)
    pct = np.array([row[3] for row in rows])
    np.testing.assert_allclose(pct.sum(), 100.0, atol=1e-9)
    sx = np.sort(x)
    for i, (lo, hi, _n, _pct, _mean) in enumerate(rows):
        a, b = int(r.bounds[i]), int(r.bounds[i + 1])
        assert lo == sx[a]
        assert hi == sx[b - 1]


def test_histogram_writes_png(tmp_path):
    x = lex1dx.datasets.load("athletics20")
    r = lex1dx.cluster(x, 5)
    path = tmp_path / "lex1dx_histogram.png"
    out = lex1dx.histogram(r, x, file=str(path), title="athletics20")
    assert Path(out) == path or Path(out).exists()
    data = path.read_bytes()
    assert data[:8] == PNG_MAGIC
    assert path.stat().st_size > 0


def test_histogram_file_none_returns_axes():
    from matplotlib.axes import Axes

    x = lex1dx.datasets.load("athletics20")
    r = lex1dx.cluster(x, 5)
    ax = lex1dx.histogram(r, x, file=None, theme="light")
    assert isinstance(ax, Axes)
