#!/usr/bin/env bash
# Builds the sequential and OpenMP-parallel binaries from the SAME source file.
#
# The source (src/ConstructorRedes2D_C4_Genetico_Final.c) contains #pragma omp
# directives. Compiled WITHOUT -fopenmp, those pragmas are silently ignored by
# the compiler and the program runs single-threaded (the "sequential" build
# used as the baseline in Fig. 14 of the paper). Compiled WITH -fopenmp, the
# same pragmas are honored and the program runs multi-threaded (the "parallel"
# build). This gives a single, byte-identical source for both configurations,
# which is what the paper itself describes ("the parallel code was designed
# to maintain the same behavior ... as the sequential version").
#
# Usage:
#   scripts/build.sh
#
# Produces:
#   bin/sequential/ga_seq
#   bin/parallel/ga_omp

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src/ConstructorRedes2D_C4_Genetico_Final.c"
CC="${CC:-gcc}"

mkdir -p "$REPO_ROOT/bin/sequential" "$REPO_ROOT/bin/parallel"

echo "Using compiler: $("$CC" --version | head -n 1)"

echo "Building sequential binary (no -fopenmp)..."
"$CC" -O2 -std=c99 -Wall "$SRC" -o "$REPO_ROOT/bin/sequential/ga_seq" -lm

echo "Building parallel (OpenMP) binary (-fopenmp)..."
"$CC" -O2 -std=c99 -Wall -fopenmp "$SRC" -o "$REPO_ROOT/bin/parallel/ga_omp" -lm

echo "Done:"
echo "  $REPO_ROOT/bin/sequential/ga_seq"
echo "  $REPO_ROOT/bin/parallel/ga_omp"
