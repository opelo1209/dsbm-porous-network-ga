#!/usr/bin/env python3
"""Averages C(r) over multiple independent realizations per Omega value,
reusing compute_correlation() from plot_correlation.py, and reports the
mean +/- SD across realizations at each radial distance r. Written to
answer Reviewer 2's request for how many independent realizations are
averaged in the spatial correlation analysis (Sec. 4.6).

Usage:
    python plot_correlation_averaged.py \
        --omega 0.3 results/correlation_avg/omega_0.3/red_colores_rep1.csv \
                     results/correlation_avg/omega_0.3/red_colores_rep2.csv \
                     results/correlation_avg/omega_0.3/red_colores_rep3.csv \
        --omega 0.6 results/correlation_avg/omega_0.6/red_colores_rep1.csv \
                     results/correlation_avg/omega_0.6/red_colores_rep2.csv \
                     results/correlation_avg/omega_0.6/red_colores_rep3.csv \
        --omega 0.9 results/correlation_avg/omega_0.9/red_colores_rep1.csv \
                     results/correlation_avg/omega_0.9/red_colores_rep2.csv \
                     results/correlation_avg/omega_0.9/red_colores_rep3.csv
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt

from plot_correlation import compute_correlation


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--omega", nargs="+", action="append", required=True,
                         help="Omega value followed by one or more red_colores.csv realizations, e.g. --omega 0.3 rep1.csv rep2.csv rep3.csv")
    parser.add_argument("--out", default="fig18_spatial_correlation_averaged.png")
    parser.add_argument("--max-r", type=int, default=15)
    args = parser.parse_args()

    plt.figure(figsize=(7, 6))
    for group in args.omega:
        omega_val, csv_paths = group[0], group[1:]
        curves = []
        for csv_path in csv_paths:
            r_bins, c_r = compute_correlation(csv_path, max_r=args.max_r)
            curves.append(c_r)
        curves = np.array(curves)  # shape (n_realizations, n_bins)
        mean_c = curves.mean(axis=0)
        std_c = curves.std(axis=0)

        print(f"=== Omega = {omega_val} ({len(csv_paths)} realizations) ===")
        for r, m, s in zip(r_bins, mean_c, std_c):
            print(f"  r={r:2d}  C(r) mean={m: .5f}  std={s: .5f}")

        # Only plot bins where the mean is positive (log axis requirement).
        mask = mean_c > 0
        plt.errorbar(r_bins[mask], mean_c[mask], yerr=std_c[mask], fmt="o-",
                     label=f"Omega = {omega_val} (n={len(csv_paths)})", capsize=3)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radial Distance r")
    plt.ylabel("Spatial Correlation C(r)")
    plt.title("Averaged Spatial Correlation for Different Overlap Values (Log-Log Scale)")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
