from importlib.resources import files

import numpy as np
import pytest

import lex1dx


def test_names_and_unknown():
    assert lex1dx.datasets.names() == (
        "athletics20",
        "worldathletics-mens-100m-2379",
    )
    with pytest.raises(ValueError, match="unknown dataset"):
        lex1dx.datasets.load("missing")


def test_athletics20():
    x = lex1dx.datasets.load("athletics20")
    assert x.dtype == np.float64
    assert x.shape == (20,)
    assert np.all(np.isfinite(x))
    assert np.all(np.diff(x) >= 0)
    assert x[0] == 9.58
    assert x[-1] == 9.84
    assert len(np.unique(x)) == 12


def test_worldathletics_2379():
    x = lex1dx.datasets.load("worldathletics-mens-100m-2379")
    assert x.shape == (2379,)
    assert np.all(np.isfinite(x))
    assert np.all(np.diff(x) >= 0)
    assert x[0] == 9.58
    assert x[19] == 9.84
    assert np.array_equal(x[:20], lex1dx.datasets.load("athletics20"))
    assert len(np.unique(x)) == 58


def test_package_data_files_are_shipped():
    root = files("lex1dx") / "data"
    assert (root / "athletics20.csv").is_file()
    assert (root / "worldathletics-mens-100m-2379.csv").is_file()
    assert (root / "README.md").is_file()
