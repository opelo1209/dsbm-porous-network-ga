#!/usr/bin/env bash
# Reproduces the three overlap snapshots behind Fig. 15 (Omega=0.3), Fig. 16
# (Omega=0.6) and Fig. 17 (Omega=0.9), whose output CSVs are also the input
# to the spatial correlation function C(r) in Fig. 18.
#
# WHERE THE (mediaS, mediaE, desviacion) TRIPLES BELOW COME FROM:
#   The current public interface takes mu_S, mu_B and a shared half-width
#   sigma directly (uniform sampling in [mu - sigma, mu + sigma]); it does
#   NOT take Omega as an input and does not print/compute it. Sec. 4.2 of the
#   paper states mu_S=410, mu_B=400, sigma=50 for Omega=0.9. Treating both
#   distributions as same-width uniform windows, the overlap fraction of one
#   window covered by the other is:
#       Omega = 1 - |mu_S - mu_B| / (2 * sigma)
#   which reproduces 0.9 exactly for (410, 400, 50) and was cross-checked
#   against two independently logged historical runs of an earlier version of
#   this codebase that printed Omega directly: (440, 400, 50) -> Omega=0.6,
#   and (440, 370, 50) -> Omega=0.3. All three triples below are therefore
#   grounded in the paper and in real prior program output, not guessed.
#
# Population size (64, Sec. 4.2's optimum) is patched into a temporary copy
# of the source exactly as in the other sweep scripts; src/ is never modified.
#
# The program always writes its CSV to a fixed relative filename,
# "red_colores.csv" (Sec. 4.5 in the paper / exportarRedConColores in the
# code); this script works around that by running each configuration in its
# own working directory.
#
# Usage:
#   scripts/run_overlap_snapshots.sh
#
# Environment overrides:
#   L                lattice size            (default: 200; the paper's Figs.
#                                             15-17 use L=500, which is heavy -
#                                             set L=500 to match the paper exactly)
#   POP_SIZE         population size         (default: 64)
#   TIMEOUT_SECONDS  per-run wall-clock cap  (default: 21600)
#   CC               compiler                (default: gcc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"
CC="${CC:-gcc}"

L="${L:-200}"
POP_SIZE="${POP_SIZE:-64}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-21600}"

OUT_DIR="$REPO_ROOT/results/fig15_16_17_18"
mkdir -p "$OUT_DIR"

PATCHED_SRC="$OUT_DIR/patched_pop${POP_SIZE}.c"
sed -E "s/int numCromosomas = [0-9]+;/int numCromosomas = ${POP_SIZE};/" "$SRC" > "$PATCHED_SRC"
BIN="$OUT_DIR/ga_omp_pop${POP_SIZE}"
"$CC" -O2 -std=c99 -fopenmp "$PATCHED_SRC" -o "$BIN" -lm

# name:mediaS:mediaE:desviacion  (desviacion=50 for all three, see header comment)
CONFIGS=(
    "omega_0.3:440:370:50"
    "omega_0.6:440:400:50"
    "omega_0.9:410:400:50"
)

for CFG in "${CONFIGS[@]}"; do
    IFS=':' read -r NAME MS ME SD <<< "$CFG"
    RUN_DIR="$OUT_DIR/$NAME"
    mkdir -p "$RUN_DIR"
    echo "=== $NAME (mediaS=$MS mediaE=$ME desviacion=$SD L=$L) ==="

    LOG="$RUN_DIR/log.txt"
    if (cd "$RUN_DIR" && timeout "$TIMEOUT_SECONDS" bash -c \
        "printf '%d\n%f\n%f\n%f\n' '$L' '$MS' '$ME' '$SD' | '$BIN' dummy dummy" \
        > "$LOG" 2>&1); then
        echo "  OK -> $RUN_DIR/red_colores.csv"
    else
        echo "  TIMEOUT_OR_ERROR (see $LOG)"
    fi
done

echo
echo "Snapshots written under $OUT_DIR/<omega_*>/red_colores.csv"
echo "Plot site distributions (Fig. 15-17) with:"
echo "  python analysis/plot_site_distribution.py $OUT_DIR/omega_0.3/red_colores.csv --title 'Omega = 0.3'"
echo "  python analysis/plot_site_distribution.py $OUT_DIR/omega_0.6/red_colores.csv --title 'Omega = 0.6'"
echo "  python analysis/plot_site_distribution.py $OUT_DIR/omega_0.9/red_colores.csv --title 'Omega = 0.9'"
echo "Plot the correlation function (Fig. 18) with:"
echo "  python analysis/plot_correlation.py $OUT_DIR/omega_0.3/red_colores.csv $OUT_DIR/omega_0.6/red_colores.csv $OUT_DIR/omega_0.9/red_colores.csv"
