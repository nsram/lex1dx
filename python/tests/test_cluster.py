import struct
from pathlib import Path

import numpy as np
import pytest

import lex1dx


def bits(x: float) -> str:
    return format(struct.unpack("<Q", struct.pack("<d", x))[0], "016x")


def from_bits(h: str) -> float:
    return struct.unpack("<d", struct.pack("<Q", int(h, 16)))[0]


def builtin():
    # Same 800 points as cpp/tests/test_agree.cpp, stored as C++ IEEE bits
    # so the comparison is not sensitive to libm differences.
    path = Path(__file__).resolve().parent / "data" / "builtin_x_bits.txt"
    return np.array([from_bits(h) for h in path.read_text().split()], dtype=np.float64)


def load_golden(path: Path):
    cases = []
    cur = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("case="):
            if cur:
                cases.append(cur)
            meta = dict(tok.split("=", 1) for tok in line.split())
            cur = {"mode": meta["mode"], "k": int(meta["k"]), "D": int(meta["D"]),
                   "Dprime": int(meta["Dprime"])}
        elif line.startswith("bounds="):
            cur["bounds"] = np.array([int(v) for v in line.split("=", 1)[1].split(",") if v],
                                     dtype=np.int64)
        elif line.startswith("sse_bits="):
            cur["sse"] = from_bits(line.split("=", 1)[1])
        elif line.startswith("centers_bits="):
            cur["centers"] = np.array([from_bits(h) for h in line.split("=", 1)[1].split(",") if h],
                                      dtype=np.float64)
    if cur:
        cases.append(cur)
    return cases


GOLDEN = Path(__file__).resolve().parent / "data" / "builtin_cpp.txt"


def test_bit_identical_to_cpp_builtin():
    x = builtin()
    cases = load_golden(GOLDEN)
    assert len(cases) == 8
    for c in cases:
        r = lex1dx.cluster(x, c["k"], mode=c["mode"])
        np.testing.assert_array_equal(r.bounds, c["bounds"])
        assert r.sse == c["sse"] and bits(r.sse) == bits(c["sse"])
        assert r.centers.shape == c["centers"].shape
        for a, b in zip(r.centers, c["centers"]):
            assert a == b and bits(a) == bits(b)
        assert r.D == c["D"]
        assert r.Dprime == c["Dprime"]


def test_empty_rejected():
    with pytest.raises(ValueError):
        lex1dx.cluster(np.array([], dtype=np.float64), 2)


def test_nan_rejected():
    with pytest.raises(ValueError):
        lex1dx.cluster(np.array([1.0, np.nan, 2.0]), 2)


def test_inf_rejected():
    with pytest.raises(ValueError):
        lex1dx.cluster(np.array([1.0, np.inf, 2.0]), 2)


def test_default_M():
    assert lex1dx.default_M(1) == 4096
    assert lex1dx.default_M(64) == 4096
    assert lex1dx.default_M(256) == 8192
    with pytest.raises(ValueError):
        lex1dx.default_M(0)


def test_assign_and_properties():
    x = np.array([10.0, 1.0, 2.0, 11.0, 1.5, 10.5])
    r = lex1dx.cluster(x, 2)
    assert r.k == 2
    assert not r.presorted
    labels = r.assign(x)
    assert labels.shape == x.shape
    assert list(labels) == [1, 0, 0, 1, 0, 1]
    r_sorted = lex1dx.cluster(np.sort(x), 2)
    assert r_sorted.presorted
    np.testing.assert_array_equal(r.bounds, r_sorted.bounds)


def test_k_clamped_to_distinct():
    x = np.array([1.0, 1.0, 2.0, 2.0, 3.0])
    r = lex1dx.cluster(x, 10)
    assert r.k == 3
    assert r.D == 3
    assert r.sse == 0.0


def test_explicit_M_and_sig5_fallback():
    # v0.1.0 linspace case: 5-sig grouping would collapse the 20 points, so
    # the fallback raises precision (or uses raw atoms) and keeps k=5.
    x = np.linspace(1.0, 1.0 + 1e-9, 20)
    r = lex1dx.cluster(x, 5, mode="sig5", M=16)
    assert r.k == 5
    assert r.D == 20
    # Tight pairs share a 5-sig key. The fallback retries at 6 digits,
    # keeps k=4, and Result.sig reports 6, the sig actually used.
    tight = np.concatenate([
        np.zeros(50),
        1.000012 + 1e-10 * (np.arange(50) - 25),
        1.000017 + 1e-10 * (np.arange(50) - 25),
        np.full(50, 2.0),
    ])
    r5 = lex1dx.cluster(tight, 4, sig=5, M=16)
    assert r5.k == 4 and r5.sig == 6 and r5.Dprime == 4
    r_alias = lex1dx.cluster(tight, 4, mode="sig5")
    assert r_alias.sig == 6 and r_alias.k == 4
    np.testing.assert_array_equal(r_alias.bounds, r5.bounds)


def test_invalid_arguments():
    x = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        lex1dx.cluster(x, 0)
    with pytest.raises(ValueError, match="mode must be"):
        lex1dx.cluster(x, 2, mode="invalid")
    with pytest.raises(ValueError):
        lex1dx.cluster(x, 2, M=-1)
    with pytest.raises(ValueError):
        lex1dx.cluster(np.zeros((2, 2)), 2)
    with pytest.raises(ValueError, match="0..9"):
        lex1dx.cluster(x, 2, sig=10)
    with pytest.raises(ValueError, match="0..9"):
        lex1dx.cluster(x, 2, sig=5.5)
    with pytest.raises(ValueError, match="solver must be"):
        lex1dx.cluster(x, 2, solver="other")


def test_weights_match_expanded_counts_and_validation():
    vals = np.array([1.0, 2.0, 3.0])
    wts = np.array([1.0, 3.0, 1.0])
    expanded = np.array([1.0, 2.0, 2.0, 2.0, 3.0])
    rw = lex1dx.cluster(vals, 2, weights=wts)
    re = lex1dx.cluster(expanded, 2)
    assert rw.sse == re.sse
    np.testing.assert_array_equal(rw.centers, re.centers)
    mass = [float(wts[lo:hi].sum()) for lo, hi in zip(rw.bounds[:-1], rw.bounds[1:])]
    np.testing.assert_array_equal(mass, re.sizes)
    with pytest.raises(ValueError, match="weights must be finite and >= 0"):
        lex1dx.cluster(vals, 2, weights=[1.0, -1.0, 1.0])
    with pytest.raises(ValueError, match="weights must be finite and >= 0"):
        lex1dx.cluster(vals, 2, weights=[1.0, np.nan, 1.0])
    with pytest.raises(ValueError, match="size must match"):
        lex1dx.cluster(vals, 2, weights=[1.0, 1.0])
    with pytest.raises(ValueError, match="positive finite total"):
        lex1dx.cluster(vals, 2, weights=[0.0, 0.0, 0.0])


def test_exact_solvers_agree_on_sse():
    rng = np.random.default_rng(1)
    x = rng.lognormal(size=2000)
    rw = lex1dx.cluster(x, 8, mode="exact", solver="wilber")
    rc = lex1dx.cluster(x, 8, mode="exact", solver="ckmeans")
    denom = max(rw.sse, rc.sse)
    rel = 0.0 if denom == 0 else abs(rw.sse - rc.sse) / denom
    assert rel < 1e-12


def test_sig_monotone_dprime_on_synthetic():
    rng = np.random.default_rng(1)
    x = rng.lognormal(size=10_000)
    dprimes = []
    prev = 0
    for d in range(1, 10):
        r = lex1dx.cluster(x, 8, sig=d)
        assert r.sig == d
        dprimes.append(r.Dprime)
        assert r.Dprime >= prev
        prev = r.Dprime
    r0 = lex1dx.cluster(x, 8, sig=0)
    rt = lex1dx.cluster(x, 8)
    assert r0.sig == 0 and r0.Dprime == rt.Dprime
    np.testing.assert_array_equal(r0.bounds, rt.bounds)
    assert r0.sse == rt.sse


def test_athletics20_k5_solvers_same_partition():
    x = lex1dx.datasets.load("athletics20")
    for mode in ("exact", "tuned"):
        rw = lex1dx.cluster(x, 5, mode=mode, solver="wilber")
        rc = lex1dx.cluster(x, 5, mode=mode, solver="ckmeans")
        assert list(rw.bounds) == list(rc.bounds), (mode, list(rw.bounds), list(rc.bounds))
        denom = max(rw.sse, rc.sse)
        rel = 0.0 if denom == 0 else abs(rw.sse - rc.sse) / denom
        assert rel < 1e-12, (mode, rw.sse, rc.sse, rel)

