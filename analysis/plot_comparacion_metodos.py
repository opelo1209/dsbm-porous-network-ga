#!/usr/bin/env python3
"""Figuras de comparacion de metodos que reunen corridas de varias carpetas de results/.

Genera en la carpeta de salida:
  convergencia.png   mediana (banda = min-max) del mejor error contra evaluaciones, por metodo
  violaciones.png    mapa de los sitios que aun violan el Principio de Construccion (semilla dada)

Uso:
  python3 analysis/plot_comparacion_metodos.py SALIDA --semilla 1 \
      --metodo "Etiqueta" results/carpeta/trazas/w0.95_variante_s{s}.csv results/carpeta/redes/w0.95_variante_s{s}_completa.csv \
      [--metodo ...]
  {s} se sustituye por cada semilla de --semillas (por defecto 1 2 3).
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def leer_traza(ruta):
    f = list(csv.DictReader(open(ruta)))
    return np.array([float(x["evaluaciones"]) for x in f]), np.array([float(x["mejor_error"]) for x in f])


def leer_errores(ruta):
    f = list(csv.DictReader(open(ruta)))
    L = max(int(x["x"]) for x in f) + 1
    m = np.zeros((L, L))
    for x in f:
        m[int(x["x"]), int(x["y"])] = int(x["errT1"])
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("salida")
    ap.add_argument("--metodo", nargs=3, action="append", metavar=("ETIQUETA", "TRAZA", "RED"), required=True)
    ap.add_argument("--semillas", nargs="+", default=["1", "2", "3"])
    ap.add_argument("--semilla", default="1", help="semilla de los mapas de violaciones")
    ap.add_argument("--titulo", default="")
    a = ap.parse_args()
    os.makedirs(a.salida, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5.2))
    for etq, traza, _ in a.metodo:
        rutas = [traza.format(s=s) for s in a.semillas if os.path.exists(traza.format(s=s))]
        if not rutas:
            continue
        trazas = [leer_traza(r) for r in rutas]
        emax = min(e[-1] for e, _ in trazas)
        rej = np.logspace(5, np.log10(emax), 300)
        curv = np.array([np.where(np.searchsorted(e, rej, side="right") - 1 >= 0,
                                  y[np.clip(np.searchsorted(e, rej, side="right") - 1, 0, None)], y[0])
                         for e, y in trazas])
        ax.plot(rej, np.median(curv, 0), label=f"{etq} (n={len(rutas)}, final {np.median(curv[:, -1]):.0f})")
        ax.fill_between(rej, curv.min(0), curv.max(0), alpha=0.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Evaluaciones de sitio (costo, igual para todos)")
    ax.set_ylabel("Violaciones del PC (mejor hasta el momento)")
    ax.set_title(a.titulo or "Convergencia con igual presupuesto")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.salida, "convergencia.png"), dpi=200)
    plt.close(fig)

    redes = [(etq, red.format(s=a.semilla)) for etq, _, red in a.metodo if os.path.exists(red.format(s=a.semilla))]
    fig, axes = plt.subplots(1, len(redes), figsize=(4.2 * len(redes), 4.6), squeeze=False)
    for ax, (etq, ruta) in zip(axes[0], redes):
        m = leer_errores(ruta)
        n = int((m > 0).sum())
        ax.imshow(m > 0, cmap="Greys", interpolation="nearest")
        ax.set_title(f"{etq}\n{int(m.sum())} violaciones en {n} sitios ({100 * n / m.size:.1f}%)", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"Sitios que aún violan el Principio de Construcción (negro), semilla {a.semilla}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(a.salida, "violaciones.png"), dpi=200)
    plt.close(fig)
    print("Guardado", os.path.join(a.salida, "convergencia.png"), "y", os.path.join(a.salida, "violaciones.png"))


if __name__ == "__main__":
    main()
