#!/usr/bin/env python3
"""Reproduces the small/medium/large site-size snapshots (Fig. 15 Omega=0.3,
Fig. 16 Omega=0.6, Fig. 17 Omega=0.9) from a red_colores.csv file produced by
the program (columns: x,y,r_Sitio,color). Adapted from the project's original
exploratory script (Visualizacion_Zonas_Coloreadas.py) into a reusable,
parameterized CLI tool.

Usage:
    python plot_site_distribution.py red_colores.csv --title "Omega = 0.9"
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import gaussian_filter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="red_colores.csv (columns: x,y,r_Sitio,color)")
    parser.add_argument("--out", default=None, help="Output image path (default: <csv_path>.png)")
    parser.add_argument("--title", default="Spatial Distribution of Site Types in the Network")
    parser.add_argument("--sigma", type=float, default=2.0,
                         help="Gaussian smoothing sigma used to pick a dominant class per pixel (cosmetic only)")
    args = parser.parse_args()
    out_path = args.out or (args.csv_path.rsplit(".", 1)[0] + "_distribution.png")

    df = pd.read_csv(args.csv_path)
    color_dict = {"yellow": 0, "blue": 1, "red": 2}
    df["color_val"] = df["color"].map(color_dict)
    L = int(max(df["x"].max(), df["y"].max()) + 1)

    class_matrices = np.zeros((3, L, L))
    for _, row in df.iterrows():
        x, y = int(row["x"]), int(row["y"])
        class_matrices[int(row["color_val"]), y, x] = 1  # note: y first, matches the array's row axis

    smoothed_classes = np.array([gaussian_filter(class_matrices[i], sigma=args.sigma) for i in range(3)])
    smoothed_labels = np.argmax(smoothed_classes, axis=0)

    cmap = ListedColormap(["yellow", "blue", "red"])
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(smoothed_labels, cmap=cmap, shading="auto")
    cbar = plt.colorbar(ticks=[0.5, 1.5, 2.5])
    cbar.ax.set_yticklabels(["Small", "Medium", "Large"])
    cbar.set_label("Site Type")

    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
