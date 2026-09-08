#!/usr/bin/env python3
"""Plots the GA's fitness (residual Type-1 violation count) as a function of
generation number, from the raw stdout of a single run of the GA binary.

The program prints a line "Optimo[<gen>] = <value>" every time the
population-best fitness improves (Sec. 3 elitist replacement); this is a
step function (best-so-far), not a per-generation sample, since fitness is
piecewise-constant between improvements under elitism.

Usage:
    python plot_convergence.py run_L50.log --out fig_convergence.png
"""
import argparse
import re

import matplotlib.pyplot as plt


def parse_trace(log_path):
    pattern = re.compile(r"Optimo\[(\d+)\]\s*=\s*(\d+)")
    gens, vals = [], []
    with open(log_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                gens.append(int(m.group(1)))
                vals.append(int(m.group(2)))
    return gens, vals


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("log_path", help="Raw stdout log of a single GA run")
    parser.add_argument("--out", default="fig_convergence.png")
    parser.add_argument("--title", default="GA Convergence: Best Fitness vs. Generation")
    args = parser.parse_args()

    gens, vals = parse_trace(args.log_path)
    if not gens:
        raise SystemExit(f"No 'Optimo[gen] = value' lines found in {args.log_path}")

    plt.figure(figsize=(7, 5))
    plt.step(gens, vals, where="post", linewidth=1.5)
    plt.xlabel("Generation")
    plt.ylabel("Best fitness (residual Type-1 violations)")
    plt.title(args.title)
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Parsed {len(gens)} improvement events, generations {gens[0]}-{gens[-1]}, "
          f"fitness {vals[0]}->{vals[-1]}")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
