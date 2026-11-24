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
