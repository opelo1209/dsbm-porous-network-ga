#!/usr/bin/env python3
"""Plots the population-size sweep produced by scripts/run_population_sweep.sh
(Table 1 analogue, Sec. 4.2). See that script's header for why this reports
generations/time TO CONVERGENCE instead of "minimum error at a fixed 13000
iteration budget": the shipped main() loop has no iteration cap, so every run
converges to error = 0 regardless of population size.

Usage:
    python plot_table1.py results/table1/table1_convergence.csv
"""
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="table1_convergence.csv produced by run_population_sweep.sh")
    parser.add_argument("--out", default="table1_convergence.png", help="Output image path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    bad = df[df["Status"] != "OK"]
    if not bad.empty:
        print(f"Warning: {len(bad)} run(s) did not finish OK and will show gaps:", file=sys.stderr)
        print(bad, file=sys.stderr)

    df = df.sort_values("PopulationSize")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1.plot(df["PopulationSize"], df["GenerationsToConvergence"], "o-", color="tab:blue")
    ax1.set_xlabel("Population size")
    ax1.set_ylabel("Generations to convergence")
    ax1.set_xscale("log", base=2)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.set_title("Convergence speed vs. population size")

    ax2.plot(df["PopulationSize"], df["ExecutionTimeSeconds"], "o-", color="tab:red")
    ax2.set_xlabel("Population size")
    ax2.set_ylabel("Execution time (s)")
    ax2.set_xscale("log", base=2)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.set_title("Execution time vs. population size")

    fig.suptitle("Population-size sweep (Table 1 analogue, L=50, Omega=0.9)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
