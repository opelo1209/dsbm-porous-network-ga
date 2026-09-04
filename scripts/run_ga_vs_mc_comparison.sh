#!/usr/bin/env bash
# Compares the pure Monte Carlo baseline (src/ConstructorRedes2D_C4_MonteCarlo.c,
# Sec. 2.2's classical approach) against the constraint-preserving genetic
# algorithm (src/ConstructorRedes2D_C4_Genetico_Final.c, Sec. 3), on the same
# lattice sizes and the same (mediaS, mediaE, desviacion) parameters - i.e.
# the same statistical inventory and the same Omega (see the Omega derivation
# in README.md).
#
# This directly reproduces, quantitatively, the qualitative claim in Sec. 1.2/
# 2.2 of the paper that blind Monte Carlo search needs a very large number of
# random exchanges compared to the population-based GA.
#
# CAVEAT: unlike the GA, the Monte Carlo baseline has no relaxation fallback
# (by design - see the header comment in ConstructorRedes2D_C4_MonteCarlo.c),
# so for large L it may hit its --max-attempts cap without reaching zero
# error. When that happens this script records Status=MAX_ATTEMPTS and the
# residual error, rather than a convergence time, which is itself a
# meaningful (and expected) data point for the comparison.
#
# Usage:
#   scripts/run_ga_vs_mc_comparison.sh
#
# Environment overrides:
#   L_VALUES         space-separated list    (default: "20 50 100")
#   POP_SIZE         GA population size      (default: 64)
#   THREADS          GA thread count          (default: 32)
#   MEDIA_S / MEDIA_E / DESVIACION            (default: 410 / 400 / 50)
#   MC_MAX_ATTEMPTS  MC max attempts (0=unlimited, dangerous) (default: 20000000)
#   TIMEOUT_SECONDS  per-run wall-clock cap   (default: 3600)
#   CC               compiler                 (default: gcc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_GA="$REPO_ROOT/src/ConstructorRedes2D_C4_Genetico_Final.c"
SRC_MC="$REPO_ROOT/src/ConstructorRedes2D_C4_MonteCarlo.c"
CC="${CC:-gcc}"

L_VALUES="${L_VALUES:-20 50 100}"
POP_SIZE="${POP_SIZE:-64}"
THREADS="${THREADS:-32}"
MEDIA_S="${MEDIA_S:-410}"
MEDIA_E="${MEDIA_E:-400}"
DESVIACION="${DESVIACION:-50}"
MC_MAX_ATTEMPTS="${MC_MAX_ATTEMPTS:-20000000}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-3600}"

OUT_DIR="$REPO_ROOT/results/ga_vs_mc"
mkdir -p "$OUT_DIR"

echo "Building pop=$POP_SIZE GA parallel binary..."
PATCHED_SRC="$OUT_DIR/patched_pop${POP_SIZE}.c"
sed -E "s/int numCromosomas = [0-9]+;/int numCromosomas = ${POP_SIZE};/" "$SRC_GA" > "$PATCHED_SRC"
BIN_GA="$OUT_DIR/ga_omp_pop${POP_SIZE}"
"$CC" -O2 -std=c99 -fopenmp "$PATCHED_SRC" -o "$BIN_GA" -lm

echo "Building Monte Carlo baseline binary..."
BIN_MC="$OUT_DIR/mc_puro"
"$CC" -O2 -std=c99 "$SRC_MC" -o "$BIN_MC" -lm

CSV="$OUT_DIR/ga_vs_mc_raw.csv"
echo "Method,L,Attempts,ExecutionTimeSeconds,FinalError,Status" > "$CSV"

for L in $L_VALUES; do
    echo "=== L=$L ==="

    # --- Monte Carlo ---
    LOG_MC="$OUT_DIR/log_mc_L${L}.txt"
    if timeout "$TIMEOUT_SECONDS" bash -c \
        "printf '%d\n%f\n%f\n%f\n%ld\n' '$L' '$MEDIA_S' '$MEDIA_E' '$DESVIACION' '$MC_MAX_ATTEMPTS' | '$BIN_MC'" \
        > "$LOG_MC" 2>&1; then
        STATUS_MC="OK"
    else
        STATUS_MC="TIMEOUT_OR_ERROR"
    fi
    ATTEMPTS_MC=$(grep -oE 'Intentos totales: [0-9]+' "$LOG_MC" | grep -oE '[0-9]+$' || true)
    TIME_MC=$(grep -oE 'Tiempo total de ejecución: [0-9.]+' "$LOG_MC" | grep -oE '[0-9.]+$' || true)
    if grep -q "Maximo de intentos alcanzado" "$LOG_MC"; then
        STATUS_MC="MAX_ATTEMPTS"
    fi
    ERR_MC=$(grep -oE 'Optimo\[[0-9]+\] = [0-9]+' "$LOG_MC" | tail -n 1 | grep -oE '[0-9]+$' || true)
    echo "MonteCarlo,${L},${ATTEMPTS_MC:-NA},${TIME_MC:-NA},${ERR_MC:-NA},${STATUS_MC}" >> "$CSV"
    echo "  [MonteCarlo] status=$STATUS_MC attempts=${ATTEMPTS_MC:-NA} time=${TIME_MC:-NA}s finalError=${ERR_MC:-NA}"

    # --- Genetic Algorithm ---
    LOG_GA="$OUT_DIR/log_ga_L${L}.txt"
    if OMP_NUM_THREADS="$THREADS" timeout "$TIMEOUT_SECONDS" bash -c \
        "printf '%d\n%f\n%f\n%f\n' '$L' '$MEDIA_S' '$MEDIA_E' '$DESVIACION' | '$BIN_GA' dummy dummy" \
        > "$LOG_GA" 2>&1; then
        STATUS_GA="OK"
    else
        STATUS_GA="TIMEOUT_OR_ERROR"
    fi
    GEN_GA=$(grep -oE 'Optimo\[[0-9]+\] = 0$' "$LOG_GA" | tail -n 1 | grep -oE '[0-9]+' | head -n 1 || true)
    TIME_GA=$(grep -oE 'Tiempo total de ejecución: [0-9.]+' "$LOG_GA" | grep -oE '[0-9.]+$' || true)
    echo "GeneticAlgorithm,${L},${GEN_GA:-NA},${TIME_GA:-NA},0,${STATUS_GA}" >> "$CSV"
    echo "  [GeneticAlgorithm] status=$STATUS_GA generations=${GEN_GA:-NA} time=${TIME_GA:-NA}s"
done

echo
echo "Raw data written to $CSV"
echo "Plot with: python analysis/plot_ga_vs_mc.py $CSV"
