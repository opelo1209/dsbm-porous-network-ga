#!/usr/bin/env bash
# Reproduces the population-size experiment behind Table 1 (Sec. 4.2 of the
# paper): population sizes {4, 8, 16, 32, 64, 128} on a 50x50 lattice with
# mu_S = 410, mu_B = 400, sigma = 50 (i.e. Omega = 0.9, see the Omega/parameter
# derivation in the README).
#
# IMPORTANT METHODOLOGICAL DIFFERENCE FROM THE PAPER (read before trusting
# numbers from this script):
#   The paper capped every Table 1 run at a FIXED budget of 13000 iterations
#   and reported the (generally non-zero) minimum error still present at that
#   point. The current public source has no such iteration cap: main()'s loop
#   is `while (mejorFitness != 0)` and only stops once the lattice is fully
#   free of Construction Principle violations (helped along by the relaxation
#   fallback, relajaSitios, once the error rate drops below ~0.001%). So this
#   script cannot reproduce a "minimum error at 13000 iterations" column - it
#   will always converge to 0 - and instead reports GENERATIONS-TO-CONVERGENCE
#   and EXECUTION-TIME-TO-CONVERGENCE per population size, which is the
#   closest faithful comparison obtainable from the shipped interface. See
#   README.md "Reproducibility notes" for the full explanation.
#
# `numCromosomas` (population size) is a compile-time constant in main(), not
# a runtime argument (the .c source in src/ is intentionally left unmodified -
# see README). This script therefore patches a TEMPORARY copy of the source
# with sed for each population size and compiles that copy; src/ itself is
# never touched.
#
# Usage:
#   scripts/run_population_sweep.sh
#
# Environment overrides:
#   POP_SIZES        space-separated list (default: "4 8 16 32 64 128")
#   L                lattice size            (default: 50)
#   MEDIA_S          mean site radius        (default: 410)
#   MEDIA_E          mean bond radius        (default: 400)
#   DESVIACION       half-width / std. dev.  (default: 50)
#   TIMEOUT_SECONDS  per-run wall-clock cap  (default: 3600)
#   CC               compiler                (default: gcc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"
CC="${CC:-gcc}"

POP_SIZES="${POP_SIZES:-4 8 16 32 64 128}"
L="${L:-50}"
MEDIA_S="${MEDIA_S:-410}"
MEDIA_E="${MEDIA_E:-400}"
DESVIACION="${DESVIACION:-50}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-3600}"

OUT_DIR="$REPO_ROOT/results/table1"
mkdir -p "$OUT_DIR"
CSV="$OUT_DIR/table1_convergence.csv"
echo "PopulationSize,GenerationsToConvergence,ExecutionTimeSeconds,Status" > "$CSV"

for P in $POP_SIZES; do
    echo "=== Population size = $P ==="

    PATCHED_SRC="$OUT_DIR/patched_pop${P}.c"
    sed -E "s/int numCromosomas = [0-9]+;/int numCromosomas = ${P};/" "$SRC" > "$PATCHED_SRC"

    BIN="$OUT_DIR/ga_pop${P}"
    "$CC" -O2 -std=c99 -fopenmp "$PATCHED_SRC" -o "$BIN" -lm

    LOG="$OUT_DIR/log_pop${P}.txt"
    # stdin order matches the program's scanf sequence: L, mediaS, mediaE, desviacion.
    # Two dummy positional args are passed because the program reads argv[2]
    # unconditionally (see README "Known limitations").
    if timeout "$TIMEOUT_SECONDS" bash -c \
        "printf '%d\n%f\n%f\n%f\n' '$L' '$MEDIA_S' '$MEDIA_E' '$DESVIACION' | '$BIN' dummy dummy" \
        > "$LOG" 2>&1; then
        STATUS="OK"
    else
        STATUS="TIMEOUT_OR_ERROR"
    fi

    GEN=$(grep -oE 'Optimo\[[0-9]+\] = 0$' "$LOG" | tail -n 1 | grep -oE '[0-9]+' | head -n 1 || true)
    TIME=$(grep -oE 'Tiempo total de ejecución: [0-9.]+' "$LOG" | grep -oE '[0-9.]+$' || true)

    echo "${P},${GEN:-NA},${TIME:-NA},${STATUS}" >> "$CSV"
    echo "  status=$STATUS generations=${GEN:-NA} time=${TIME:-NA}s (log: $LOG)"
done

echo
echo "Summary written to $CSV"
echo "Plot it with: python analysis/plot_table1.py $CSV"
