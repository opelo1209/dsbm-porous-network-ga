#!/usr/bin/env python3
"""Plots the pure Monte Carlo baseline vs. the genetic algorithm comparison
produced by scripts/run_ga_vs_mc_comparison.sh - a quantitative version of
the qualitative claim in Sec. 1.2/2.2 of the paper that blind Monte Carlo
search needs far more random exchanges than the population-based GA.

Runs with Status != OK (timed out, or MC hit its --max-attempts cap without
reaching zero error) are shown as open/hollow markers rather than dropped,
since "Monte Carlo did not converge in the budget" is itself part of the
comparison's point.

Usage:
    python plot_ga_vs_mc.py results/ga_vs_mc/ga_vs_mc_raw.csv
"""
import argparse

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="ga_vs_mc_raw.csv produced by run_ga_vs_mc_comparison.sh")
    parser.add_argument("--out", default="ga_vs_mc.png", help="Output image path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path).sort_values(["Method", "L"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"MonteCarlo": "tab:red", "GeneticAlgorithm": "tab:blue"}

    for method, group in df.groupby("Method"):
        ok = group[group["Status"] == "OK"]
        not_ok = group[group["Status"] != "OK"]

        ax1.plot(ok["L"], ok["ExecutionTimeSeconds"], "o-", label=method, color=colors.get(method))
        if not not_ok.empty:
            ax1.scatter(not_ok["L"], not_ok["ExecutionTimeSeconds"], facecolors="none",
                        edgecolors=colors.get(method), marker="o", s=80,
                        label=f"{method} (did not fully converge)")

        ax2.plot(ok["L"], ok["Attempts"], "o-", label=method, color=colors.get(method))
        if not not_ok.empty:
            ax2.scatter(not_ok["L"], not_ok["Attempts"], facecolors="none",
                        edgecolors=colors.get(method), marker="o", s=80)

    ax1.set_xlabel("Lattice size (L)")
    ax1.set_ylabel("Execution time (s)")
    ax1.set_yscale("log")
    ax1.set_title("Execution time to convergence")
    ax1.legend(fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.5)

    ax2.set_xlabel("Lattice size (L)")
    ax2.set_ylabel("Attempts / generations to convergence")
    ax2.set_yscale("log")
    ax2.set_title("Search effort to convergence")
    ax2.grid(True, linestyle="--", alpha=0.5)

    fig.suptitle("Pure Monte Carlo baseline vs. constraint-preserving GA")
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
