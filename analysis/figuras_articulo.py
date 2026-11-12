#!/usr/bin/env python3
"""Figures of the revised manuscript (PDF), generated from the experiment folders in results/.

Usage:
  python3 analysis/figuras_articulo.py OUTPUT_DIR [--resultados results]

Generates (when the required data exist):
  fig_population.pdf  final violation rate per seed and variant (L=200, Omega 0.9 and 0.95; evidence 02)
  fig_convergence.pdf best violation count vs. evaluations (L=200, Omega=0.95; evidence 02, 04, 05)
  fig_violations.pdf  maps of the sites that still violate the Construction Principle (seed 1)
  fig_methods_maps.pdf raw site-category maps for the same runs (seed 1)
  fig_omega.pdf       raw maps of the final method for each Omega (L=100) and C(r) averaged over seeds
  tabla_omega.csv     per-Omega summary of the final method (used by the manuscript table)
Category colours follow the convention of the manuscript (small = yellow, medium = blue, large = red).
"""
import argparse
import csv
import glob
import os

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# Categorical palette validated with the dataviz checks (light surface); low-contrast
# hues always carry a second encoding (marker or line style) and a legend.
C = {"ga64": "#2a78d6", "ga16": "#6250d6", "ga4": "#eb6834", "ea": "#1baf7a", "eaf": "#eda100"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
CAT = ListedColormap(["#f2c80f", "#2156c6", "#d0312d"])

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.4, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def resumen(carpeta):
    ruta = os.path.join(carpeta, "resumen.csv")
    return list(csv.DictReader(open(ruta))) if os.path.exists(ruta) else []


def traza(ruta):
    f = list(csv.DictReader(open(ruta)))
    return np.array([float(x["evaluaciones"]) for x in f]), np.array([float(x["mejor_error"]) for x in f])


def red(ruta):
    """Returns (site radii L x L, errors L x L) from a *_completa.csv file."""
    f = list(csv.DictReader(open(ruta)))
    L = max(int(x["x"]) for x in f) + 1
    r, e = np.zeros((L, L)), np.zeros((L, L), dtype=int)
    for x in f:
        i, j = int(x["x"]), int(x["y"])
        r[i, j], e[i, j] = float(x["r_Sitio"]), int(x["errT1"])
    return r, e


def categorias(r):
    mn, mx = r.min(), r.max()
    x1, x2 = mn + (mx - mn) / 3, mx - (mx - mn) / 3
    return np.where(r <= x1, 0, np.where(r <= x2, 1, 2))


def misma_cat(k):
    return 0.5 * (np.mean(k == np.roll(k, 1, 0)) + np.mean(k == np.roll(k, 1, 1)))


def correlacion_radial(r):
    """C(r) of the site radii on the periodic lattice (FFT autocorrelation, integer distance bins)."""
    L = r.shape[0]
    z = r - r.mean()
    ac = np.fft.ifft2(np.abs(np.fft.fft2(z)) ** 2).real / (z.size * z.var())
    d = np.minimum(np.arange(L), L - np.arange(L))
    dist = np.sqrt(d[:, None] ** 2 + d[None, :] ** 2)
    rb = np.rint(dist).astype(int)
    rmax = L // 2
    return np.arange(1, rmax + 1), np.array([ac[rb == k].mean() for k in range(1, rmax + 1)])


def fig_population(res, salida):
    filas = resumen(os.path.join(res, "experimento_regimen_grande"))
    if not filas:
        return
    orden = [("busquedaLocal_P1", "(1+1)-EA", C["ea"], "o"), ("GA4_rango_rep2", "GA $\\mu$=4", C["ga4"], "s"),
             ("GA_P16", "GA $\\mu$=16", C["ga16"], "^"), ("GA_P64", "GA $\\mu$=64", C["ga64"], "D")]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4), sharey=True)
    for ax, om in zip(axes, ["0.900", "0.950"]):
        for k, (var, etq, col, mk) in enumerate(orden):
            v = [float(f["error_pct"]) for f in filas if f["variante"] == var and f["omega"] == om]
            x = k + np.linspace(-0.12, 0.12, len(v))
            ax.scatter(x, v, s=18, color=col, marker=mk, edgecolor="white", linewidth=0.6, zorder=3)
            ax.hlines(np.median(v), k - 0.25, k + 0.25, color=INK, linewidth=1.0, zorder=4)
        ax.set_yscale("log")
        ax.set_xticks(range(len(orden)), [o[1] for o in orden])
        ax.set_title(f"$\\Omega$ = {float(om):g}, L = 200")
        ax.grid(True, axis="y", color=GRID, linewidth=0.5, which="both")
    axes[0].set_ylabel("Residual violations (% of incidences)")
    fig.savefig(os.path.join(salida, "fig_population.pdf"))
    plt.close(fig)


def metodos_L200(res):
    E, F, V = (os.path.join(res, d) for d in ("experimento_regimen_grande", "piloto_ga_focal", "piloto_v3_L200"))
    return [("GA $\\mu$=64", E, "GA_P64", C["ga64"], "-"),
            ("GA $\\mu$=4, conflict-directed", F, "GA4_pmx_focal", C["ga4"], "--"),
            ("(1+1)-EA", E, "busquedaLocal_P1", C["ea"], "-."),
            ("(1+1)-EA, conflict-directed", V, "focal0.75_P1", C["eaf"], "-")]


def fig_convergence(res, salida):
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for etq, d, var, col, ls in metodos_L200(res):
        rutas = sorted(glob.glob(os.path.join(d, "trazas", f"w0.95_{var}_s[123].csv")))
        if not rutas:
            continue
        tr = [traza(p) for p in rutas]
        rej = np.logspace(5, np.log10(min(e[-1] for e, _ in tr)), 300)
        cur = np.array([y[np.clip(np.searchsorted(e, rej, side="right") - 1, 0, None)] for e, y in tr])
        ax.plot(rej, np.median(cur, 0), color=col, linestyle=ls, label=etq)
        ax.fill_between(rej, cur.min(0), cur.max(0), color=col, alpha=0.18, linewidth=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Site evaluations (cost)")
    ax.set_ylabel("Best violation count")
    ticks = [3000, 5000, 10000, 20000, 40000]
    ax.set_yticks(ticks, [f"{t:,}" for t in ticks])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid(True, color=GRID, linewidth=0.5, which="major")
    ax.legend(frameon=False, loc="lower left")
    fig.savefig(os.path.join(salida, "fig_convergence.pdf"))
    plt.close(fig)


def fig_maps(res, salida):
    redes = []
    for etq, d, var, _, _ in metodos_L200(res):
        p = os.path.join(d, "redes", f"w0.95_{var}_s1_completa.csv")
        if os.path.exists(p):
            redes.append((etq, red(p)))
    if not redes:
        return
    fig, axes = plt.subplots(1, len(redes), figsize=(1.75 * len(redes), 2.0))
    for ax, (etq, (r, e)) in zip(axes, redes):
        ax.imshow(e > 0, cmap=ListedColormap(["white", INK]), interpolation="nearest")
        ax.set_title(f"{etq}\n{int((e > 0).sum())} sites ({100 * (e > 0).mean():.1f}%)", fontsize=6.5)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(True); s.set_color(MUTED); s.set_linewidth(0.5)
    fig.savefig(os.path.join(salida, "fig_violations.pdf"))
    plt.close(fig)
    fig, axes = plt.subplots(1, len(redes), figsize=(1.75 * len(redes), 2.0))
    for ax, (etq, (r, e)) in zip(axes, redes):
        k = categorias(r)
        ax.imshow(k, cmap=CAT, vmin=0, vmax=2, interpolation="nearest")
        ax.set_title(f"{etq}\nsame category {100 * misma_cat(k):.1f}%", fontsize=6.5)
        ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(os.path.join(salida, "fig_methods_maps.pdf"), dpi=300)
    plt.close(fig)


def fig_omega(res, salida):
    d = os.path.join(res, "omega_final")
    filas = resumen(d)
    if not filas:
        return
    omegas = sorted({f["omega"] for f in filas}, key=float)
    var = "focal0.75_P1"
    fig = plt.figure(figsize=(7.0, 2.0), layout="constrained")
    gs = fig.add_gridspec(1, len(omegas) + 1, width_ratios=[1] * len(omegas) + [1.9])
    estilos = ["-", "--", "-.", ":"]
    axc = fig.add_subplot(gs[0, -1])
    tabla = []
    for k, om in enumerate(omegas):
        cs, mismas = [], []
        for p in sorted(glob.glob(os.path.join(d, "redes", f"w{float(om):g}_{var}_s*_completa.csv"))):
            r, e = red(p)
            rr, c = correlacion_radial(r)
            cs.append(c)
            mismas.append(misma_cat(categorias(r)))
            if p.endswith("_s1_completa.csv"):
                ax = fig.add_subplot(gs[0, k])
                ax.imshow(categorias(r), cmap=CAT, vmin=0, vmax=2, interpolation="nearest")
                ax.set_title(f"$\\Omega$ = {float(om):g}", fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
        if cs:
            cs = np.array(cs)
            m = cs.mean(0)
            sel = rr <= 25
            axc.plot(rr[sel], m[sel], color=[C["ea"], C["ga64"], C["ga4"], C["eaf"]][k % 4],
                     linestyle=estilos[k % 4], marker="o", markersize=2.5, label=f"$\\Omega$ = {float(om):g}")
            fo = [f for f in filas if f["omega"] == om and f["variante"] == var]
            tabla.append({"omega": om, "n": len(fo),
                          "error_pct_median": np.median([float(f["error_pct"]) for f in fo]),
                          "llegaron_a_cero": sum(int(f["error"]) == 0 for f in fo),
                          "misma_cat_mean": 100 * np.mean([float(f["misma_cat"]) for f in fo]),
                          "misma_cat_sd": 100 * np.std([float(f["misma_cat"]) for f in fo], ddof=1),
                          "C1_mean": np.mean([float(f["C1"]) for f in fo]),
                          "C1_sd": np.std([float(f["C1"]) for f in fo], ddof=1),
                          "C5_mean": m[rr == 5][0], "C10_mean": m[rr == 10][0]})
    axc.axhline(0, color=MUTED, linewidth=0.5)
    axc.set_xlabel("Distance $r$ (lattice units)")
    axc.set_ylabel("C(r)")
    axc.grid(True, color=GRID, linewidth=0.5)
    axc.legend(frameon=False)
    fig.savefig(os.path.join(salida, "fig_omega.pdf"), dpi=300)
    plt.close(fig)
    with open(os.path.join(salida, "tabla_omega.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tabla[0].keys()))
        w.writeheader(); w.writerows(tabla)
    for t in tabla:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in t.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("salida")
    ap.add_argument("--resultados", default="results")
    a = ap.parse_args()
    os.makedirs(a.salida, exist_ok=True)
    fig_population(a.resultados, a.salida)
    fig_convergence(a.resultados, a.salida)
    fig_maps(a.resultados, a.salida)
    fig_omega(a.resultados, a.salida)
    print("Figuras en", a.salida, ":", sorted(os.listdir(a.salida)))


if __name__ == "__main__":
    main()
