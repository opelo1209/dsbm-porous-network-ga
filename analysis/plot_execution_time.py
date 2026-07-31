#!/usr/bin/env python3
"""Reproduces Fig. 12 (execution time vs. lattice size, one line per thread
count) from the raw CSV produced by scripts/run_thread_sweep.sh.

Usage:
    python plot_execution_time.py results/fig12_13/fig12_13_raw.csv
"""
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="fig12_13_raw.csv produced by run_thread_sweep.sh")
    parser.add_argument("--out", default="fig12_execution_time.png", help="Output image path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    bad = df[df["Status"] != "OK"]
    if not bad.empty:
        print(f"Warning: {len(bad)} run(s) did not finish OK and are excluded:", file=sys.stderr)
        print(bad, file=sys.stderr)
    df = df[df["Status"] == "OK"].sort_values(["Threads", "L"])

    fig, ax = plt.subplots(figsize=(8, 6))
    for threads, group in df.groupby("Threads"):
        ax.plot(group["L"], group["ExecutionTimeSeconds"], "o-", label=f"{threads} Threads")

    ax.set_xlabel("Size of Porous Networks (L)")
    ax.set_ylabel("Execution time (sec)")
    ax.set_title("Execution time across lattice sizes and thread configurations")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
