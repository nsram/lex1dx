#!/usr/bin/env bash
# Stage a self-contained CRAN tarball (with C++ core synced into src/ and
# dev-only configure scripts removed so CRAN does not see no-op configure files).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

(cd "$ROOT" && sh tools/sync_core.sh)
cp -R "$ROOT" "$STAGE/lex1dx"
rm -f "$STAGE/lex1dx/configure" "$STAGE/lex1dx/configure.win"
rm -rf "$STAGE/lex1dx/tools" "$STAGE/lex1dx/src"/*.o "$STAGE/lex1dx/src"/*.so "$STAGE/lex1dx/src"/*.dll

(cd "$STAGE" && R CMD build lex1dx && R CMD check --as-cran lex1dx_*.tar.gz)
cp "$STAGE"/lex1dx_*.tar.gz "$ROOT/"
echo "Built CRAN tarball in $ROOT/"
