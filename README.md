# lex1dx

lex1dx is a 1-D k-means library (least-squares partition of a sorted array into k contiguous bins).

It continues [lex1d](https://github.com/nsram/lex1d) ([Zenodo](https://zenodo.org/records/24451573)).

lex1d stays the paper-faithful exact-quality reference.

lex1dx is a new repository because the objectives are different: it targets performance against a 0.01% relative-excess-SSE (REE) tolerance, so its defaults and algorithms are not the same as lex1d.

## Defaults

The default method is the tuned pipeline.

It sorts the data (and canonicalizes `-0` to `+0`).

Optional 5-significant-digit grouping (`Mode::Sig5`) can be turned on; the default `Mode::Tuned` keeps full-precision atoms (exact ties only).

Collapse uses a 0.75-quantile batch rule and finishes with a heap once at most 16384 merges remain.

The Collapse target is `M = max(4096, 32k)`.

The exact dynamic program on those M groups is replaced by a weighted wilber-interp DP, derived from Grønlund's kmeans1d (commit `585386d`) with a weighted interval-cost patch.

Refine is unchanged in role: it locally improves the cuts.

k larger than the number of atoms is clamped to that number.

## Build and use (C++)

From the repository root, with CMake 3.21+ and a C++17 compiler:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Link the `lex1dx` library and include `<lex1dx/lex1dx.hpp>`:

```cpp
#include <lex1dx/lex1dx.hpp>
lex1dx::Result r = lex1dx::cluster(values, 64);                 // tuned, default M
lex1dx::Options o;
o.mode = lex1dx::Mode::Sig5;
lex1dx::Result s = lex1dx::cluster(values, 64, o);              // 5-significant-digit grouping
```

`Result` holds `bounds`, `sizes`, `centers`, `upper`, `sse`, `D`, `Dprime`, and `assign(x)`.

The `lex1dx_dump` tool prints bounds, centers, and SSE for a little-endian float64 file:

```
./build/lex1dx_dump tuned data.f64 64
```

## Build and use (Scala)

Scala 3.3 and sbt 1.10, from the `scala/` directory:

```
sbt test
```

```scala
import lex1dx.{Lex, Options, Mode}
val r = Lex.cluster(values, 64)
val s = Lex.cluster(values, 64, Options(mode = Mode.Sig5))
```

## Build and use (Python)

From the `python/` directory, with a C++17 compiler and CMake 3.21+:

```
pip install .
pip install pytest
pytest
```

```python
import numpy as np
import lex1dx
r = lex1dx.cluster(x, 64)                         # tuned, default M
s = lex1dx.cluster(x, 64, mode="sig5")            # 5-significant-digit grouping
t = lex1dx.cluster(x, 64, M=8192)                 # override Collapse target
```

`M=None` (the default) uses `default_M(k) = max(4096, 32 * k)`.

## Build and use (R)

From the `r/` directory, with a C++17 compiler. `configure` copies the C++ core into `src/`. Run `sh tools/sync_core.sh` before `R CMD build` as well:

```
sh tools/sync_core.sh
R CMD build .
R CMD check --no-manual --no-vignettes lex1dx_0.1.0.tar.gz
```

```r
library(lex1dx)
r <- cluster(x, 64L)                              # tuned, default M
s <- cluster(x, 64L, mode = "sig5")
t <- cluster(x, 64L, M = 8192L)
```

`M = NULL` (the default) uses `default_M(k) = max(4096, 32 * k)`.

## Results (v7 Mac pack)

On the v7 scaling pack, a post-sort comparison against Grønlund wilber-interp on the raw points gave speedups up to 156× (sig5) and 37× (tuned), with REE at most 0.0056% on all 270/270 cells.

Those times exclude the sort, because both methods need the same sorted input.

Those figures come from the `bot_v7_scaling_mac_pack` Mac run (270 cells, all ok).

## Citation

lex1d is at https://github.com/nsram/lex1d and is archived at https://zenodo.org/records/24451573.

## Attribution

lex1dx is MIT licensed (see `LICENSE`).

The wilber DP is derived from kmeans1d by Allan Grønlund, commit `585386d4678000691706acdb954911be8fea1c5c`, also MIT; see `NOTICE` and `third_party/kmeans1d/`.
