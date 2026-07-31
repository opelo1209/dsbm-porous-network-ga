#!/usr/bin/env bash
# Captures the exact toolchain/hardware environment a run was produced with,
# so results can be attributed to a specific compiler version and machine.
# This is what should be committed alongside any results/ artifact you share,
# since the compiler version and -fopenmp support are not otherwise recorded
# anywhere in the source.
#
# Usage:
#   scripts/env_info.sh > results/environment.txt

set -euo pipefail

CC="${CC:-gcc}"

echo "=== Date ==="
date -u

echo
echo "=== Compiler ==="
"$CC" --version || echo "$CC not found"
"$CC" -dumpmachine || true
echo "OpenMP support check (expands to _OPENMP date macro if enabled):"
echo | "$CC" -fopenmp -dM -E - 2>/dev/null | grep -i openmp || echo "  _OPENMP macro not found"

echo
echo "=== OS ==="
uname -a || true

echo
echo "=== CPU ==="
if command -v lscpu >/dev/null 2>&1; then
    lscpu
elif command -v wmic >/dev/null 2>&1; then
    wmic cpu get Name,NumberOfCores,NumberOfLogicalProcessors
else
    echo "nproc: $(nproc 2>/dev/null || echo unknown)"
fi

echo
echo "=== OMP_NUM_THREADS (current shell) ==="
echo "${OMP_NUM_THREADS:-<unset, OpenMP will default to all logical cores>}"
