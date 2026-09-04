#!/usr/bin/env bash
# Reproduces the sequential-vs-parallel execution time comparison behind
# Fig. 14 (Sec. 4.4): population = 64, Omega = 0.9 (mu_S=410, mu_B=400,
# sigma=50), 32 OpenMP threads for the parallel build, across lattice sizes L.
#
# Requires scripts/build.sh to have been run first (uses bin/sequential/ga_seq
# and bin/parallel/ga_omp). Population size (64) is compile-time in main();
# this script builds its own pop=64 patched copies of both binaries the same
# way run_population_sweep.sh / run_thread_sweep.sh do - src/ is never
# modified.
#
# Same convergence-loop caveat as the other sweep scripts applies to any error
# metric (always converges to 0); only ExecutionTimeSeconds is a faithful
# reproduction target here, which is exactly what Fig. 14 plots.
#
# Usage:
#   scripts/run_seq_vs_parallel.sh
#
# Environment overrides:
#   L_VALUES         space-separated list   (default: "50 100 150" - use the
#                                             full "50 100 150 200 250 300 350
#                                             400 450 500" to reproduce Fig. 14
#                                             in full; expect multi-hour runs
#                                             for the sequential build at large L)
#   POP_SIZE         population size        (default: 64)
#   THREADS_PARALLEL thread count            (default: 32, per Sec. 4.3/4.4)
#   MEDIA_S / MEDIA_E / DESVIACION           (default: 410 / 400 / 50)
#   TIMEOUT_SECONDS  per-run wall-clock cap  (default: 21600)
#   CC               compiler                (default: gcc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"
CC="${CC:-gcc}"

L_VALUES="${L_VALUES:-50 100 150}"
POP_SIZE="${POP_SIZE:-64}"
THREADS_PARALLEL="${THREADS_PARALLEL:-32}"
MEDIA_S="${MEDIA_S:-410}"
MEDIA_E="${MEDIA_E:-400}"
DESVIACION="${DESVIACION:-50}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-21600}"

OUT_DIR="$REPO_ROOT/results/fig14"
mkdir -p "$OUT_DIR"

PATCHED_SRC="$OUT_DIR/patched_pop${POP_SIZE}.c"
sed -E "s/int numCromosomas = [0-9]+;/int numCromosomas = ${POP_SIZE};/" "$SRC" > "$PATCHED_SRC"

BIN_SEQ="$OUT_DIR/ga_seq_pop${POP_SIZE}"
BIN_OMP="$OUT_DIR/ga_omp_pop${POP_SIZE}"
echo "Building sequential (no -fopenmp) pop=$POP_SIZE binary..."
"$CC" -O2 -std=c99 "$PATCHED_SRC" -o "$BIN_SEQ" -lm
echo "Building parallel (-fopenmp) pop=$POP_SIZE binary..."
"$CC" -O2 -std=c99 -fopenmp "$PATCHED_SRC" -o "$BIN_OMP" -lm

CSV="$OUT_DIR/fig14_raw.csv"
echo "Mode,Threads,L,ExecutionTimeSeconds,Status" > "$CSV"

run_one () {
    local MODE="$1" BIN="$2" THREADS="$3" L="$4"
    local LOG="$OUT_DIR/log_${MODE}_L${L}.txt"
    if OMP_NUM_THREADS="$THREADS" timeout "$TIMEOUT_SECONDS" bash -c \
        "printf '%d\n%f\n%f\n%f\n' '$L' '$MEDIA_S' '$MEDIA_E' '$DESVIACION' | '$BIN' dummy dummy" \
        > "$LOG" 2>&1; then
        STATUS="OK"
    else
        STATUS="TIMEOUT_OR_ERROR"
    fi
    local TIME
    TIME=$(grep -oE 'Tiempo total de ejecución: [0-9.]+' "$LOG" | grep -oE '[0-9.]+$' || true)
    echo "${MODE},${THREADS},${L},${TIME:-NA},${STATUS}" >> "$CSV"
    echo "  [$MODE] L=$L status=$STATUS time=${TIME:-NA}s (log: $LOG)"
}

for L in $L_VALUES; do
    echo "=== L=$L ==="
    run_one "sequential" "$BIN_SEQ" 1 "$L"
    run_one "parallel"   "$BIN_OMP" "$THREADS_PARALLEL" "$L"
done

echo
echo "Raw data written to $CSV"
echo "Plot with: python analysis/plot_seq_vs_parallel.py $CSV"
