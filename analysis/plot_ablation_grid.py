#!/usr/bin/env python3
"""Combines several ablation (unrestricted crossover+mutation) site-distribution
snapshots at different lattice sizes into a single multi-panel figure, each
panel annotated with its same-category-neighbor clustering fraction, to show
visually and quantitatively that the result is statistically indistinguishable
from spatial randomness (~1/3) regardless of lattice size - written for the
Ablation Study (Sec. 4.7) rebuttal figure.

Usage:
    python plot_ablation_grid.py \
        --entry 20  red_colores_ablacion_L20.csv \
        --entry 50  red_colores_ablacion_L50.csv \
        --entry 200 red_colores_ablacion_L200.csv \
        --entry 500 red_colores_ablacion_L500.csv \
        --out ablation_grid.png
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import gaussian_filter


def load_grid(csv_path, L):
    df = pd.read_csv(csv_path)
    cmap = {"yellow": 0, "blue": 1, "red": 2}
    grid = np.zeros((L, L), dtype=int)
    for _, row in df.iterrows():
        grid[int(row["y"]), int(row["x"])] = cmap[row["color"]]
    return grid


def clustering_fraction(grid, L):
    same = 0
    for i in range(L):
        for j in range(L):
            if grid[i, j] == grid[i, (j + 1) % L]:
                same += 1
            if grid[i, j] == grid[(i + 1) % L, j]:
                same += 1
    return same / (2 * L * L)


def smoothed_labels(grid, L, sigma=1.0):
    class_matrices = np.zeros((3, L, L))
    for c in range(3):
        class_matrices[c] = (grid == c).astype(float)
    smoothed = np.array([gaussian_filter(class_matrices[c], sigma=sigma) for c in range(3)])
    return np.argmax(smoothed, axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--entry", nargs=2, action="append", required=True, metavar=("L", "CSV"),
                         help="Lattice size followed by its red_colores.csv path")
    parser.add_argument("--out", default="ablation_grid.png")
    parser.add_argument("--sigma", type=float, default=1.0, help="Cosmetic smoothing sigma (kept small so small lattices stay legible)")
    args = parser.parse_args()

    entries = [(int(L), path) for L, path in args.entry]
    n = len(entries)
    ncols = 2
    nrows = (n + 1) // 2

    cmap = ListedColormap(["yellow", "blue", "red"])
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = np.array(axes).reshape(-1)

    for ax, (L, path) in zip(axes, entries):
        grid = load_grid(path, L)
        frac = clustering_fraction(grid, L)
        labels = smoothed_labels(grid, L, sigma=args.sigma)
        ax.pcolormesh(labels, cmap=cmap, shading="auto")
        ax.set_title(f"L = {L}  (same-category neighbors: {frac*100:.1f}%)", fontsize=11)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_aspect("equal")

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle("Unrestricted crossover+mutation ablation: site distribution vs. lattice size\n"
                  "(random-mixing baseline = 33.3%)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.out, dpi=300)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
