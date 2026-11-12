#!/usr/bin/env python3
"""Analiza el EXPERIMENTO 01 (aporte del cruce) producido por scripts/experimento_cruce.sh.

Genera dentro de la carpeta del experimento:
  tabla_resumen.md        tabla por (Omega, variante): error final, % misma categoria,
                          C(1), corridas que llegaron a cero, evaluaciones y tiempo
  pruebas_estadisticas.md Mann-Whitney U (dos colas) de la referencia (GA_P64, o la variable de entorno REFERENCIA) contra cada variante
  convergencia_omega_X.png  mediana (y rango intercuartil) del mejor error vs. evaluaciones
  cajas_error.png         distribucion del error final por variante

Uso: python3 analysis/resumen_cruce.py results/experimento_cruce
"""
import csv
import glob
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

REF = os.environ.get("REFERENCIA", "GA_P64")
ORDEN = ["GA_P64", "GA_P16", "sinCruce_P64", "busquedaLocal_P1"]
ETIQUETA = {
    "GA_P64": "GA P=64 (cruce+mutacion)",
    "GA_P16": "GA P=16 (cruce+mutacion)",
    "sinCruce_P64": "P=64 sin cruce",
    "busquedaLocal_P1": "Busqueda local P=1",
}


def leer_resumen(carpeta):
    filas = list(csv.DictReader(open(os.path.join(carpeta, "resumen.csv"), newline="")))
    datos = defaultdict(list)
    for f in filas:
        datos[(f["omega"], f["variante"])].append(f)
        if f["variante"] not in ORDEN:          # variantes nuevas: al final, en orden de aparicion
            ORDEN.append(f["variante"])
    return datos


def fmt(v, d=1):
    return f"{np.mean(v):.{d}f} ± {np.std(v, ddof=1) if len(v) > 1 else 0:.{d}f}"


def tabla(datos, carpeta):
    lineas = ["| Ω | Variante | n | Error final (media ± DE) | Mediana | Error % incidencias | Llegaron a 0 | "
              "% misma categoría | C(1) | Evaluaciones usadas | Tiempo (s) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for om in sorted({k[0] for k in datos}):
        for var in ORDEN:
            f = datos.get((om, var))
            if not f:
                continue
            err = [int(x["error"]) for x in f]
            pct = [float(x["error_pct"]) for x in f]
            mis = [100 * float(x["misma_cat"]) for x in f]
            c1 = [float(x["C1"]) for x in f]
            ev = [float(x["evaluaciones"]) for x in f]
            t = [float(x["tiempo_s"]) for x in f]
            ok = all(x["inventario_ok"] == "1" and x["fitness_ok"] == "1" for x in f)
            lineas.append(f"| {om} | {ETIQUETA.get(var, var)} | {len(f)} | {fmt(err)} | {np.median(err):.0f} | {fmt(pct, 3)} | "
                          f"{sum(e == 0 for e in err)}/{len(f)} | {fmt(mis)} | {fmt(c1, 3)} | {np.mean(ev):.3g} | {fmt(t)} |"
                          + ("" if ok else " ⚠ verificación falló"))
    texto = "\n".join(lineas) + "\n"
    open(os.path.join(carpeta, "tabla_resumen.md"), "w", encoding="utf-8").write(texto)
    return texto


def pruebas(datos, carpeta):
    lineas = [f"| Ω | Comparación | Mediana {REF} | Mediana otra | U | p (dos colas) | Conclusión (α=0.05) |",
              "|---|---|---|---|---|---|---|"]
    for om in sorted({k[0] for k in datos}):
        base = [int(x["error"]) for x in datos.get((om, REF), [])]
        if not base:
            continue
        for var in [v for v in ORDEN if v != REF]:
            otra = [int(x["error"]) for x in datos.get((om, var), [])]
            if not otra:
                continue
            if len(set(base + otra)) == 1:
                u, p = float("nan"), 1.0
            else:
                u, p = mannwhitneyu(base, otra, alternative="two-sided")
            if p >= 0.05:
                concl = "sin diferencia significativa"
            elif np.median(base) < np.median(otra):
                concl = f"{REF} mejor"
            else:
                concl = f"{var} mejor"
            lineas.append(f"| {om} | {REF} vs {var} | {np.median(base):.0f} | {np.median(otra):.0f} | {u:.1f} | {p:.4g} | {concl} |")
    texto = "\n".join(lineas) + "\n"
    open(os.path.join(carpeta, "pruebas_estadisticas.md"), "w", encoding="utf-8").write(texto)
    return texto


def curvas(carpeta, datos):
    rutas = []
    for om in sorted({k[0] for k in datos}):
        fig, ax = plt.subplots(figsize=(7.5, 5))
        for var in ORDEN:
            trazas = sorted(glob.glob(os.path.join(carpeta, "trazas", f"w{float(om):g}_{var}_s*.csv")))
            if not trazas:
                continue
            emax = max(float(list(csv.DictReader(open(t)))[-1]["evaluaciones"]) for t in trazas)
            rejilla = np.logspace(np.log10(1e5), np.log10(emax), 200)
            curv = []
            for t in trazas:
                f = list(csv.DictReader(open(t)))
                e = np.array([float(x["evaluaciones"]) for x in f])
                y = np.array([float(x["mejor_error"]) for x in f])
                idx = np.searchsorted(e, rejilla, side="right") - 1
                curv.append(np.where(idx >= 0, y[np.clip(idx, 0, None)], y[0]))
            curv = np.array(curv)
            med = np.median(curv, 0)
            q1, q3 = np.percentile(curv, 25, 0), np.percentile(curv, 75, 0)
            ax.plot(rejilla, med + 1, label=f"{ETIQUETA.get(var, var)} (n={len(trazas)})")
            ax.fill_between(rejilla, q1 + 1, q3 + 1, alpha=0.2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Evaluaciones de sitio (costo)")
        ax.set_ylabel("Mejor error + 1 (mediana, banda = RIC)")
        ax.set_title(f"Convergencia con igual presupuesto, Ω = {om}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        ruta = os.path.join(carpeta, f"convergencia_omega_{om}.png")
        fig.savefig(ruta, dpi=200)
        plt.close(fig)
        rutas.append(ruta)
    return rutas


def cajas(carpeta, datos):
    omegas = sorted({k[0] for k in datos})
    fig, axes = plt.subplots(1, len(omegas), figsize=(6 * len(omegas), 4.8), squeeze=False)
    for ax, om in zip(axes[0], omegas):
        vals, etq = [], []
        for var in ORDEN:
            f = datos.get((om, var))
            if f:
                vals.append([int(x["error"]) for x in f])
                etq.append(ETIQUETA.get(var, var).replace(" (", "\n("))
        ax.boxplot(vals, labels=etq)
        ax.set_title(f"Error final, Ω = {om}")
        ax.set_ylabel("Violaciones residuales")
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    ruta = os.path.join(carpeta, "cajas_error.png")
    fig.savefig(ruta, dpi=200)
    plt.close(fig)
    return ruta


def main():
    carpeta = sys.argv[1] if len(sys.argv) > 1 else "results/experimento_cruce"
    datos = leer_resumen(carpeta)
    print(tabla(datos, carpeta))
    print(pruebas(datos, carpeta))
    for r in curvas(carpeta, datos):
        print("Guardado", r)
    print("Guardado", cajas(carpeta, datos))


if __name__ == "__main__":
    main()
