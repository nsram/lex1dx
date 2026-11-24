#!/usr/bin/env bash
# Same commands as .github/workflows/ci.yml, for local validation.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export LEX1DX_AGREE_DIR="${LEX1DX_AGREE_DIR:-/tmp/lex1dx-agree}"
mkdir -p "$LEX1DX_AGREE_DIR"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
./build/test_agree
if ! command -v sbt >/dev/null; then
  echo "sbt not found (install sbt 1.10+ and JDK 17, and put sbt on PATH)" >&2
  exit 1
fi
(cd scala && sbt test)

# Python job: pip install . + pytest, from python/
if ! command -v python >/dev/null && ! command -v python3 >/dev/null; then
  echo "python not found" >&2
  exit 1
fi
PYTHON="${PYTHON:-python3}"
(cd python && "$PYTHON" -m pip install . && "$PYTHON" -m pip install pytest && "$PYTHON" -m pytest)

# R job: sync, R CMD build, R CMD check
if ! command -v R >/dev/null; then
  echo "R not found (install R and the Rcpp, testthat packages)" >&2
  exit 1
fi
(cd r && sh tools/sync_core.sh && R CMD build . && R CMD check --no-manual --no-vignettes lex1dx_0.1.0.tar.gz)
