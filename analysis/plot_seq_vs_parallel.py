#!/usr/bin/env python3
"""Reproduces Fig. 14 (sequential vs. parallel execution time vs. lattice
size) from the raw CSV produced by scripts/run_seq_vs_parallel.sh.

Usage:
    python plot_seq_vs_parallel.py results/fig14/fig14_raw.csv
"""
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="fig14_raw.csv produced by run_seq_vs_parallel.sh")
    parser.add_argument("--out", default="fig14_seq_vs_parallel.png", help="Output image path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    bad = df[df["Status"] != "OK"]
    if not bad.empty:
        print(f"Warning: {len(bad)} run(s) did not finish OK and are excluded:", file=sys.stderr)
        print(bad, file=sys.stderr)
    df = df[df["Status"] == "OK"].sort_values(["Mode", "L"])

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {"sequential": "tab:orange", "parallel": "black"}
    for mode, group in df.groupby("Mode"):
        ax.plot(group["L"], group["ExecutionTimeSeconds"], "o-",
                 label=("Sequential GA" if mode == "sequential" else "Parallel GA (N Threads)"),
                 color=colors.get(mode))

    ax.set_xlabel("Size of Porous Networks (L)")
    ax.set_ylabel("Execution time (sec)")
    ax.set_title("Sequential vs. parallel GA execution time")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
