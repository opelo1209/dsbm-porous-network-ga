#!/usr/bin/env python3
"""Verificacion de la organizacion espacial SIN suavizado (hallazgo A2).

Dos modos:

  python3 analysis/verificar_red.py --aleatoria 200
      Genera una red de categorias 100% aleatoria y muestra cuanto "agrupamiento"
      aparece solo por aplicar el mismo suavizado que plot_site_distribution.py
      (gaussian_filter + argmax).

  python3 analysis/verificar_red.py red_colores.csv [otro.csv ...]
      Para cada CSV producido por el GA/MC (x,y,r_Sitio,color) reporta, sobre los
      datos CRUDOS: % de vecinos de misma categoria, el valor esperado si el
      arreglo fuera aleatorio (sum p_k^2), y C(r) para r = 1..5; y, para comparar,
      el % de misma categoria tras el suavizado usado en las figuras.

Solo requiere numpy y scipy.
"""
import argparse
import csv
import numpy as np
from scipy.ndimage import gaussian_filter

COLORES = {"yellow": 0, "blue": 1, "red": 2}


def misma_categoria(c):
    return ((c == np.roll(c, 1, 0)).mean() + (c == np.roll(c, 1, 1)).mean()) / 2


def suavizar(c, sigma):
    capas = [gaussian_filter((c == k).astype(float), sigma=sigma) for k in range(3)]
    return np.argmax(capas, axis=0)


def correlacion(r_sitio, rmax=5):
    z = (r_sitio - r_sitio.mean()) / r_sitio.std()
    ac = np.fft.ifft2(np.abs(np.fft.fft2(z)) ** 2).real / z.size
    L = z.shape[0]
    d = np.minimum(np.arange(L), L - np.arange(L))
    dist = np.hypot(*np.meshgrid(d, d, indexing="ij"))
    return [ac[(dist >= r - 0.5) & (dist < r + 0.5)].mean() for r in range(1, rmax + 1)]


def leer_csv(ruta):
    filas = list(csv.DictReader(open(ruta, newline="")))
    L = max(max(int(f["x"]), int(f["y"])) for f in filas) + 1
    cat = np.zeros((L, L), int)
    rad = np.zeros((L, L))
    for f in filas:
        i, j = int(f["x"]), int(f["y"])
        cat[i, j] = COLORES[f["color"]]
        rad[i, j] = float(f["r_Sitio"])
    return cat, rad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="*")
    ap.add_argument("--aleatoria", type=int, metavar="L")
    ap.add_argument("--semilla", type=int, default=1)
    a = ap.parse_args()

    if a.aleatoria:
        rng = np.random.default_rng(a.semilla)
        c = rng.integers(0, 3, (a.aleatoria, a.aleatoria))
        print(f"Red ALEATORIA L={a.aleatoria}")
        print(f"  crudo:                 {100 * misma_categoria(c):5.1f}% vecinos de misma categoria")
        for s in (1.0, 2.0):
            print(f"  suavizado sigma={s}:    {100 * misma_categoria(suavizar(c, s)):5.1f}%"
                  "   <- solo por el filtro de la grafica")

    for ruta in a.csv:
        cat, rad = leer_csv(ruta)
        p = np.bincount(cat.ravel(), minlength=3) / cat.size
        print(f"\n{ruta}  (L={cat.shape[0]})")
        print(f"  proporciones S/M/L: {p[0]:.3f} / {p[1]:.3f} / {p[2]:.3f}")
        print(f"  CRUDO   misma categoria: {100 * misma_categoria(cat):5.1f}%   "
              f"(esperado si fuera aleatorio: {100 * (p ** 2).sum():5.1f}%)")
        print(f"  SUAVIZADO sigma=2:       {100 * misma_categoria(suavizar(cat, 2.0)):5.1f}%   (como en las figuras)")
        print("  C(r) crudo, r=1..5:      " + "  ".join(f"{v:+.3f}" for v in correlacion(rad)))


if __name__ == "__main__":
    main()
