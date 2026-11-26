#!/bin/sh
# Copies the C++ core from the repository into src/, so the package is self-contained.
# configure and configure.win run this before every install from a checkout. Run it yourself, from
# the r/ directory, before `R CMD build`. In a source tarball there is no parent tree and it does nothing.
set -e
root=..
[ -f "$root/cpp/include/lex1dx/lex1dx.hpp" ] || exit 0
cp "$root"/cpp/src/*.cpp "$root"/cpp/src/*.hpp src/
mkdir -p src/lex src/lex1dx
cp "$root/cpp/include/lex/lex.hpp" src/lex/
cp "$root/cpp/include/lex1dx/lex1dx.hpp" src/lex1dx/
# Flatten kmeans1d sources and headers next to the rest of the core.
cp "$root"/third_party/kmeans1d/*.cpp src/
cp "$root"/third_party/kmeans1d/*.h src/
cp "$root"/third_party/kmeans1d/*.hpp src/

# R CMD check warns on .hpp names directly in src/ ("unlikely file names for src
# files"). Public C++ headers stay as .hpp. Only the copied src-root headers are
# renamed, and includes in the copied tree are rewritten to match.
for f in detail interval_sum siggroup tuned_collapse; do
  if [ -f "src/${f}.hpp" ]; then
    mv "src/${f}.hpp" "src/${f}.h"
  fi
done
for file in src/*.cpp src/*.h src/*.hpp; do
  [ -f "$file" ] || continue
  tmp="$file.tmp$$"
  sed -e 's/"detail.hpp"/"detail.h"/g' \
      -e 's/"interval_sum.hpp"/"interval_sum.h"/g' \
      -e 's/"siggroup.hpp"/"siggroup.h"/g' \
      -e 's/"tuned_collapse.hpp"/"tuned_collapse.h"/g' \
      "$file" > "$tmp" && mv "$tmp" "$file"
done
