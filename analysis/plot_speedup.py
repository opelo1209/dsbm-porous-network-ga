#!/usr/bin/env python3
"""Computes and plots the speedup factor S_N(L) = T_1(L) / T_N(L) for each
thread count N, from the raw CSV produced by scripts/run_thread_sweep.sh
(columns: Threads, L, ExecutionTimeSeconds, FinalMinError, Status).

This directly answers the reviewer request to compare thread configurations
via relative speedup rather than absolute execution time alone, which
depends on the specific hardware used.

Usage:
    python plot_speedup.py results/fig12_13/fig12_13_raw.csv --out fig12_speedup.png

Requires that the sweep include a single-thread (N=1) or otherwise a
designated baseline column; if the raw CSV's smallest thread count is not 1
(e.g., the default sweep starts at N=4), pass --baseline-threads to specify
which column to divide by (must be one of the thread counts present in the
CSV). An ideal-speedup reference line (S_N = N) is included for comparison.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_path", help="fig12_13_raw.csv produced by run_thread_sweep.sh")
    parser.add_argument("--out", default="fig12_speedup.png", help="Output image path")
    parser.add_argument("--baseline-threads", type=int, default=None,
                         help="Thread count to use as the T_1 baseline (default: smallest thread count present)")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    df = df[df["Status"] == "OK"].copy()
    df["ExecutionTimeSeconds"] = pd.to_numeric(df["ExecutionTimeSeconds"], errors="coerce")
    df = df.dropna(subset=["ExecutionTimeSeconds"])

    baseline_n = args.baseline_threads or df["Threads"].min()
    baseline = df[df["Threads"] == baseline_n].set_index("L")["ExecutionTimeSeconds"]

    plt.figure(figsize=(7, 6))
    for n in sorted(df["Threads"].unique()):
        sub = df[df["Threads"] == n].set_index("L")["ExecutionTimeSeconds"]
        common_L = sorted(set(sub.index) & set(baseline.index))
        speedup = [baseline[L] / sub[L] for L in common_L]
        plt.plot(common_L, speedup, "o-", label=f"N = {n} threads")

    # Ideal linear speedup reference for the largest thread count tested.
    all_L = sorted(df["L"].unique())
    max_n = df["Threads"].max()
    plt.plot(all_L, [max_n / baseline_n] * len(all_L), "k--", alpha=0.4,
              label=f"Ideal S={max_n // baseline_n}x (N={max_n})")

    plt.xlabel("Lattice size L")
    plt.ylabel(f"Speedup $S_N = T_{{{baseline_n}}}/T_N$")
    plt.title("Speedup factor across lattice sizes and thread configurations")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
