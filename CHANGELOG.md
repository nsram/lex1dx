# Changelog

Dates are the GitHub tag dates, in America/New_York.

## 0.1.4 — 2026-11-27

DOI: https://doi.org/10.5281/zenodo.24607940

When the default Wilber solver's λ search fails (an empty interval, typically on collinear SSE*(k), for example a regular grid or `cluster([0, 1, 2, 3], 3)`), lex1dx falls back to the Fisher / Ckmeans DP instead of throwing. If that fallback's back-pointer table would exceed 2 GiB, the table is not allocated and a clear error is raised. The result now includes `solver_used` (`"wilber"` or `"ckmeans"`), which records which DP ran (when the DP is skipped, it mirrors the requested solver).

The R `tolerance` argument is now checked to lie in (0, 1), matching the C++ core. The maintainer ORCID is recorded. Example output headers that still said "from lex1dx 0.1.1" now say 0.1.4. The main README now uses absolute URL links.

## 0.1.3 — 2026-11-26

DOI: https://doi.org/10.5281/zenodo.24572680

Search keywords and k-means descriptions; no code changes.

## 0.1.2 — 2026-11-26

DOI: https://doi.org/10.5281/zenodo.24571810

`histogram()` in Python no longer writes a file by default; R and Python now match.

## 0.1.1 — 2026-11-26

DOI: https://doi.org/10.5281/zenodo.24571093

First release on PyPI (`pip install lex1dx`; add `[plot]` for histograms).

- `histogram()`: bin chart in the style of the original audio output (blue background, one yellow bar per bin, with n, percent and mean per row; r² and MSE in the header). `theme="light"` for print.
- `bin_table()`: per-bin lo, hi, n, percent and mean, in Python and R.
- Sample histogram PNGs in `docs/img` and the R vignette.
- README rewritten: C++ core with Python and R bindings and a Scala port, solver options, weights, histogram usage, citation.
- R histogram draws bars with `rect()`, so it works on any PNG device and locale (fixes macOS R CMD check).
- Builds on older toolchains: removed `std::to_chars(double)`; marked overriding methods in `kmeans.h`.
- Scala agreement test tolerates 1-ulp differences and CRLF dumps from Windows.
- Wheels for Python 3.10–3.13 on Linux (manylinux_2_28), macOS 10.15+ (Intel and Apple Silicon) and Windows, published by trusted publishing.
- Cites Sriram, N. (2026). LEX: Fast scalar quantization by grouping and local exchanges. https://doi.org/10.5281/zenodo.24543220
- Author name N. Sriram in R and PyPI metadata.

## 0.1.0 — 2026-11-25

DOI: https://doi.org/10.5281/zenodo.24544282

First public release of lex1dx, a fast least-squares scalar quantizer (1D k-means) using grouping and local exchanges.

- C++ core with Python and R bindings, plus the table-free exact reference solver used for benchmarking.
- Full-precision and precision-grouped (sig5) modes; precision-grouped mode keeps peak memory near 0.1 GiB on ~8M values at k = 512.
- Tests for C++, Python and R; CI on all three.
- Fixes the Python `Result` dataclass field order that broke import in the previous commit.

Companion paper and replication package: https://github.com/nsram/lex1dx-paper
