#!/usr/bin/env python3
"""Attempts to reproduce Fig. 13 (minimum error vs. lattice size, one line per
thread count) from the raw CSV produced by scripts/run_thread_sweep.sh.

READ THIS FIRST: the shipped main() loop runs until the lattice is fully free
of Construction Principle violations (no fixed 13000-iteration budget), so
FinalMinError in the CSV is expected to be 0 for every row - unlike Fig. 13 in
the paper, which reports the non-zero error remaining after a fixed budget.
This script will plot whatever is in the CSV and loudly warn if everything is
zero (the expected outcome with the current public code), rather than pretend
it reproduced the original figure. See scripts/run_thread_sweep.sh and
README.md "Reproducibility notes" for the full explanation.

Usage:
    python plot_min_error.py results/fig12_13/fig12_13_raw.csv
"""
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="fig12_13_raw.csv produced by run_thread_sweep.sh")
    parser.add_argument("--out", default="fig13_min_error.png", help="Output image path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    df = df[df["Status"] == "OK"].sort_values(["Threads", "L"])

    if (df["FinalMinError"].fillna(0) == 0).all():
        print(
            "NOTE: every FinalMinError is 0, as expected from the current "
            "unbounded convergence loop (no 13000-iteration cap). This plot "
            "will be flat at zero and is NOT a reproduction of the paper's "
            "Fig. 13. See this script's module docstring.",
            file=sys.stderr,
        )

    fig, ax = plt.subplots(figsize=(8, 6))
    for threads, group in df.groupby("Threads"):
        ax.plot(group["L"], group["FinalMinError"], "o-", label=f"{threads} Threads")

    ax.set_xlabel("Size of Porous Networks (L)")
    ax.set_ylabel("Minimum error (Type 1 violations at convergence)")
    ax.set_title("Minimum error across lattice sizes and thread configurations")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
