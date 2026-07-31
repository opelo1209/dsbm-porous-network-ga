#!/usr/bin/env bash
# Reproduces the thread-count / lattice-size sweep behind Fig. 12 (execution
# time) and, with the caveat below, Fig. 13 (minimum error) of the paper
# (Sec. 4.3): population = 64, Omega = 0.9 (mu_S=410, mu_B=400, sigma=50),
# thread counts {4,8,16,32,64,128}, lattice sizes L in {50,100,...,500}.
#
# IMPORTANT METHODOLOGICAL DIFFERENCE FROM THE PAPER:
#   Just like the population sweep (see run_population_sweep.sh), the shipped
#   main() loop runs until the lattice is fully free of violations (no fixed
#   13000-iteration budget), so every run's FINAL error converges to 0
#   regardless of thread count or L. This script can therefore faithfully
#   reproduce Fig. 12 (execution time vs. L per thread count) but NOT Fig. 13
#   as originally plotted (a non-zero "minimum error at a fixed budget"
#   curve) - the "MinError" column recorded here is always 0 and is kept only
#   for completeness/sanity-checking, not as a Fig. 13 substitute.
#
# Population size (64) is compile-time in main(); this script patches a
# TEMPORARY copy of the source with sed, exactly like run_population_sweep.sh.
# src/ is never modified. Thread count IS a runtime setting here, via the
# standard OpenMP OMP_NUM_THREADS environment variable (no source change
# needed for that part).
#
# Usage:
#   scripts/run_thread_sweep.sh
#
# Environment overrides:
#   THREAD_COUNTS    space-separated list   (default: "4 8 16 32 64 128")
#   L_VALUES         space-separated list   (default: "50 100 150" - the full
#                                             paper range is "50 100 150 200
#                                             250 300 350 400 450 500", which
#                                             is very heavy: set L_VALUES to
#                                             that explicitly to reproduce
#                                             Fig. 12/13 in full)
#   POP_SIZE         population size        (default: 64, per Sec. 4.3)
#   MEDIA_S          mean site radius        (default: 410)
#   MEDIA_E          mean bond radius        (default: 400)
#   DESVIACION       half-width / std. dev.  (default: 50)
#   TIMEOUT_SECONDS  per-run wall-clock cap  (default: 21600, i.e. 6h - large
#                                             L with few threads is genuinely
#                                             slow, see the paper's own
#                                             ~28,000s sequential figure)
#   CC               compiler                (default: gcc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src/ConstructorRedes2D_C4_Genetico_Final.c"
CC="${CC:-gcc}"

THREAD_COUNTS="${THREAD_COUNTS:-4 8 16 32 64 128}"
L_VALUES="${L_VALUES:-50 100 150}"
POP_SIZE="${POP_SIZE:-64}"
MEDIA_S="${MEDIA_S:-410}"
MEDIA_E="${MEDIA_E:-400}"
DESVIACION="${DESVIACION:-50}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-21600}"

OUT_DIR="$REPO_ROOT/results/fig12_13"
mkdir -p "$OUT_DIR"

echo "Building the pop=$POP_SIZE parallel binary used for the whole sweep..."
PATCHED_SRC="$OUT_DIR/patched_pop${POP_SIZE}.c"
sed -E "s/int numCromosomas = [0-9]+;/int numCromosomas = ${POP_SIZE};/" "$SRC" > "$PATCHED_SRC"
BIN="$OUT_DIR/ga_omp_pop${POP_SIZE}"
"$CC" -O2 -std=c99 -fopenmp "$PATCHED_SRC" -o "$BIN" -lm

CSV="$OUT_DIR/fig12_13_raw.csv"
echo "Threads,L,ExecutionTimeSeconds,FinalMinError,Status" > "$CSV"

for T in $THREAD_COUNTS; do
    for L in $L_VALUES; do
        echo "=== threads=$T L=$L ==="
        LOG="$OUT_DIR/log_t${T}_L${L}.txt"
        if OMP_NUM_THREADS="$T" timeout "$TIMEOUT_SECONDS" bash -c \
            "printf '%d\n%f\n%f\n%f\n' '$L' '$MEDIA_S' '$MEDIA_E' '$DESVIACION' | '$BIN' dummy dummy" \
            > "$LOG" 2>&1; then
            STATUS="OK"
        else
            STATUS="TIMEOUT_OR_ERROR"
        fi

        TIME=$(grep -oE 'Tiempo total de ejecución: [0-9.]+' "$LOG" | grep -oE '[0-9.]+$' || true)
        ERR=$(grep -oE 'Optimo\[[0-9]+\] = [0-9]+' "$LOG" | tail -n 1 | grep -oE '[0-9]+$' || true)

        echo "${T},${L},${TIME:-NA},${ERR:-NA},${STATUS}" >> "$CSV"
        echo "  status=$STATUS time=${TIME:-NA}s finalError=${ERR:-NA} (log: $LOG)"
    done
done

echo
echo "Raw data written to $CSV"
echo "Plot execution time (Fig. 12) with: python analysis/plot_execution_time.py $CSV"
echo "Plot min error       (Fig. 13) with: python analysis/plot_min_error.py $CSV"
