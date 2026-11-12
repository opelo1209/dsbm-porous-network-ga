#!/usr/bin/env python3
"""Dibuja la distribucion de sitios (pequeno/mediano/grande) SIN suavizado.

Cada celda del PNG es exactamente un sitio de la red (interpolacion 'nearest'),
asi que lo que se ve es lo que hay: no se aplica gaussian_filter ni argmax como
en plot_site_distribution.py. En el titulo de cada panel se reportan las
metricas crudas (% de vecinos de misma categoria y C(1)).

Uso:
  # un PNG por archivo
  python3 analysis/plot_red_cruda.py omega_0.3_colores.csv omega_0.9_colores.csv

  # ademas, un solo PNG con todos los paneles lado a lado
  python3 analysis/plot_red_cruda.py a.csv b.csv c.csv --panel comparacion.png

  # titulos personalizados (uno por CSV)
  python3 analysis/plot_red_cruda.py a.csv b.csv --titulos "Omega=0.3" "Omega=0.9"
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

COLORES = {"yellow": 0, "blue": 1, "red": 2}
CMAP = ListedColormap(["#f2d600", "#1f4fd8", "#d62020"])
LEYENDA = [Patch(color="#f2d600", label="Pequeno"),
           Patch(color="#1f4fd8", label="Mediano"),
           Patch(color="#d62020", label="Grande")]


def leer(ruta):
    filas = list(csv.DictReader(open(ruta, newline="")))
    L = max(max(int(f["x"]), int(f["y"])) for f in filas) + 1
    cat = np.zeros((L, L), int)
    rad = np.zeros((L, L))
    for f in filas:
        i, j = int(f["x"]), int(f["y"])
        cat[i, j] = COLORES[f["color"]]
        rad[i, j] = float(f["r_Sitio"])
    return cat, rad


def metricas(cat, rad):
    misma = ((cat == np.roll(cat, 1, 0)).mean() + (cat == np.roll(cat, 1, 1)).mean()) / 2
    z = (rad - rad.mean()) / rad.std()
    c1 = ((z * np.roll(z, 1, 0)).mean() + (z * np.roll(z, 1, 1)).mean()) / 2
    return misma, c1


def dibujar(ax, cat, rad, titulo):
    misma, c1 = metricas(cat, rad)
    ax.imshow(cat, cmap=CMAP, vmin=0, vmax=2, interpolation="nearest", origin="lower")
    ax.set_title(f"{titulo}\nL={cat.shape[0]} | misma categoria={100 * misma:.1f}% | C(1)={c1:+.2f}", fontsize=10)
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    return misma, c1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", help="Archivos *_colores.csv (x,y,r_Sitio,color)")
    ap.add_argument("--titulos", nargs="+", help="Un titulo por CSV (por defecto, el nombre del archivo)")
    ap.add_argument("--panel", help="Ruta de un PNG adicional con todos los paneles juntos")
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    titulos = a.titulos or [os.path.basename(p).replace("_colores.csv", "") for p in a.csv]
    if len(titulos) != len(a.csv):
        ap.error("--titulos debe tener un titulo por CSV")

    datos = [leer(p) for p in a.csv]

    for ruta, (cat, rad), titulo in zip(a.csv, datos, titulos):
        fig, ax = plt.subplots(figsize=(6, 6.4))
        dibujar(ax, cat, rad, titulo)
        ax.legend(handles=LEYENDA, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False)
        fig.tight_layout()
        salida = ruta.rsplit(".", 1)[0] + "_cruda.png"
        fig.savefig(salida, dpi=a.dpi)
        plt.close(fig)
        print(f"Guardado {salida}")

    if a.panel:
        n = len(datos)
        fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.8), squeeze=False)
        for ax, (cat, rad), titulo in zip(axes[0], datos, titulos):
            dibujar(ax, cat, rad, titulo)
        fig.legend(handles=LEYENDA, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.suptitle("Distribucion de sitios SIN suavizado (cada pixel = un sitio)", fontsize=12)
        fig.tight_layout(rect=[0, 0.07, 1, 0.95])
        fig.savefig(a.panel, dpi=a.dpi)
        plt.close(fig)
        print(f"Guardado {a.panel}")


if __name__ == "__main__":
    main()
