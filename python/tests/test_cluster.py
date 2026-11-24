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
