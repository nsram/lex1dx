"""Bundled example datasets.

    import lex1dx
    x = lex1dx.datasets.load("athletics20")

Source: World Athletics 100m all-time lists, retrieved January 2024.

These files are sample inputs for testing and demonstrating the
algorithms in this package.

Sample Data Attribution: The sample 100m timing datasets (top 20 performances and 2379 all time bests) contain historical athletic performance records compiled from publicly available statistics published by World Athletics (https://worldathletics.org/). These datasets are provided strictly as static benchmark/example inputs for testing and demonstrating the algorithms in this package. This project is not affiliated with, endorsed by, or sponsored by World Athletics.
"""
from __future__ import annotations

from importlib.resources import files

import numpy as np

_DATASETS = {
    "athletics20": "athletics20.csv",
    "worldathletics-mens-100m-2379": "worldathletics-mens-100m-2379.csv",
}


def names() -> tuple[str, ...]:
    """The dataset names accepted by ``load``."""
    return tuple(_DATASETS)


def load(name: str) -> np.ndarray:
    """Return a 1-D float64 vector for a bundled dataset.

    ``name`` is ``"athletics20"`` or ``"worldathletics-mens-100m-2379"``.
    """
    if name not in _DATASETS:
        known = ", ".join(_DATASETS)
        raise ValueError(f"unknown dataset {name!r}; known names: {known}")
    text = (files(__package__) / "data" / _DATASETS[name]).read_text(encoding="utf-8")
    rows = []
    for i, line in enumerate(text.splitlines()):
        line = line.strip()
        if not line:
            continue
        if i == 0 and line.lower() == "time":
            continue
        rows.append(float(line))
    if not rows:
        raise ValueError(f"dataset {name!r} is empty")
    return np.asarray(rows, dtype=np.float64)
