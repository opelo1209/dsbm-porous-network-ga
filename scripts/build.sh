#!/usr/bin/env bash
# Builds the sequential and OpenMP-parallel GA binaries, plus the pure Monte
# Carlo baseline, from their source files.
#
# The GA source (src/genetico/ConstructorRedes2D_C4_Genetico_Final.c) contains #pragma
# omp directives. Compiled WITHOUT -fopenmp, those pragmas are silently
# ignored by the compiler and the program runs single-threaded (the
# "sequential" build used as the baseline in Fig. 14 of the paper). Compiled
# WITH -fopenmp, the same pragmas are honored and the program runs
# multi-threaded (the "parallel" build). This gives a single, byte-identical
# source for both configurations, which is what the paper itself describes
# ("the parallel code was designed to maintain the same behavior ... as the
# sequential version").
#
# The Monte Carlo baseline (src/montecarlo/ConstructorRedes2D_C4_MonteCarlo.c, Sec. 2.2's
# classical baseline) is inherently sequential - each exchange depends on the
# outcome of the previous one - so it has no OpenMP variant.
#
# Usage:
#   scripts/build.sh
#
# Produces:
#   bin/sequential/ga_seq
#   bin/parallel/ga_omp
#   bin/montecarlo/mc_puro

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_GA="$REPO_ROOT/src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"
SRC_MC="$REPO_ROOT/src/montecarlo/ConstructorRedes2D_C4_MonteCarlo.c"
CC="${CC:-gcc}"

mkdir -p "$REPO_ROOT/bin/sequential" "$REPO_ROOT/bin/parallel" "$REPO_ROOT/bin/montecarlo"

echo "Using compiler: $("$CC" --version | head -n 1)"

echo "Building sequential GA binary (no -fopenmp)..."
"$CC" -O3 -std=c99 -Wall "$SRC_GA" -o "$REPO_ROOT/bin/sequential/ga_seq" -lm

echo "Building parallel (OpenMP) GA binary (-fopenmp)..."
"$CC" -O3 -std=c99 -Wall -fopenmp "$SRC_GA" -o "$REPO_ROOT/bin/parallel/ga_omp" -lm

echo "Building pure Monte Carlo baseline binary..."
"$CC" -O3 -std=c99 -Wall "$SRC_MC" -o "$REPO_ROOT/bin/montecarlo/mc_puro" -lm

echo "Done:"
echo "  $REPO_ROOT/bin/sequential/ga_seq"
echo "  $REPO_ROOT/bin/parallel/ga_omp"
echo "  $REPO_ROOT/bin/montecarlo/mc_puro"
