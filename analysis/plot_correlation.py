#!/usr/bin/env python3
"""Reproduces Fig. 18 (spatial correlation function C(r) on log-log axes for
several Omega values) from one or more red_colores.csv files, implementing
Eq. 4 of the paper directly:

    C(r) = < (R_Si - R_S_bar)(R_Sj - R_S_bar) > / sigma^2

averaged over all site pairs (i, j) separated by lattice distance r, using
the site radii only (not the small/medium/large color labels). The lattice is
periodic (Sec. 3.1), so the correlation is computed as a circular
autocorrelation via FFT and radially averaged over integer-rounded Euclidean
distance bins under periodic (wrap-around) distance.

This replaces the project's earlier ad hoc correlacion_articulo.py, which
read the CSV with header=None as if it were a raw numeric matrix - it does
not match the program's actual "x,y,r_Sitio,color" CSV format and was not
reused here.

Usage:
    python plot_correlation.py results/fig15_16_17_18/omega_0.3/red_colores.csv \\
                                results/fig15_16_17_18/omega_0.6/red_colores.csv \\
                                results/fig15_16_17_18/omega_0.9/red_colores.csv
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compute_correlation(csv_path, max_r=15):
    df = pd.read_csv(csv_path)
    L = int(max(df["x"].max(), df["y"].max()) + 1)

    grid = np.zeros((L, L))
    xs = df["x"].to_numpy(dtype=int)
    ys = df["y"].to_numpy(dtype=int)
    rs = df["r_Sitio"].to_numpy(dtype=float)
    grid[ys, xs] = rs

    mean = grid.mean()
    std = grid.std()
    norm = (grid - mean) / std

    # Circular (periodic) autocorrelation via FFT; ac[0,0] before normalizing
    # equals the r=0 variance term, used here as C(0) = 1.
    spectrum = np.fft.fft2(norm)
    ac = np.fft.ifft2(spectrum * np.conj(spectrum)).real / (L * L)
    ac /= ac[0, 0]

    # Periodic (wrap-around) distance along each axis, then radial bins.
    idx = np.arange(L)
    axis_dist = np.minimum(idx, L - idx)
    dx, dy = np.meshgrid(axis_dist, axis_dist, indexing="ij")
    dist = np.sqrt(dx.astype(float) ** 2 + dy.astype(float) ** 2)

    max_r = min(max_r, L // 2 - 1)
    r_bins = np.arange(1, max_r + 1)
    c_r = np.array([ac[(dist >= r - 0.5) & (dist < r + 0.5)].mean() for r in r_bins])
    return r_bins, c_r


def label_for(csv_path):
    match = re.search(r"omega_([0-9.]+)", csv_path.replace(os.sep, "/"))
    if match:
        return f"Empirical C(r) (Overlap {match.group(1)})"
    return f"Empirical C(r) ({os.path.basename(csv_path)})"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_paths", nargs="+", help="One or more red_colores.csv files")
    parser.add_argument("--out", default="fig18_spatial_correlation.png", help="Output image path")
    parser.add_argument("--max-r", type=int, default=15, help="Maximum radial distance to plot")
    args = parser.parse_args()

    plt.figure(figsize=(7, 6))
    for csv_path in args.csv_paths:
        r_bins, c_r = compute_correlation(csv_path, max_r=args.max_r)
        plt.loglog(r_bins, c_r, "o-", label=label_for(csv_path))

    plt.xlabel("Radial Distance r")
    plt.ylabel("Spatial Correlation C(r)")
    plt.title("Empirical Spatial Correlation for Different Overlap Values (Log-Log Scale)")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
