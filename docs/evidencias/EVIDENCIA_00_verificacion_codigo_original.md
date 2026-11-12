# EVIDENCIA 00 — Verificación del código original (antes de las correcciones)

- **Fecha:** 2026-11-11
- **Código verificado:** `src/genetico/ConstructorRedes2D_C4_Genetico_Final.c` y `src/ablacion/ConstructorRedes2D_C4_Ablacion.c` (commit `3a50809`)
- **Herramientas:**
  - `tests/verificar_operadores.c`: arnés que incluye el `.c` original sin modificarlo y llama directamente a sus operadores.
  - `analysis/verificar_red.py`: métricas sobre datos crudos.
  - `scripts/verificar.sh`: corre todo.
- **Reproducir (Ubuntu/WSL):** `bash scripts/verificar.sh`. Resultados en `results/verificacion/`.
- **Plataformas:** Windows 11 (MinGW-w64 GCC 16.1.0, UCRT) y WSL2 Ubuntu 20.04 (gcc 9.4.0).

Configuración del arnés: L = 50, P = 16, 300 generaciones de cruce + mutación, μS = 410, μB = 400, a = 50 (Ω = 0.9), semilla 12345.

## A1 — La mutación con restricción de categoría congela el mapa de categorías de sitios

La mutación solo intercambia sitios de la **misma** categoría, y el cruce solo mueve enlaces. Por eso la
categoría de cada posición queda fija desde la inicialización aleatoria.

| Medida (Ubuntu) | GA original | Ablación (sin restricción) |
|---|---|---|
| Cromosomas finales cuyo mapa de categorías es idéntico a uno inicial | **16 de 16** | 0 de 16 |
| % de vecinos de la misma categoría, inicial → final | 33.2% → 33.7% | 33.2% → 33.6% |

CSV del GA original guardados en `results/correlation_avg/` (Ω = 0.3, 0.6, 0.9; 3 réplicas cada uno):
todos dan 32–35% de vecinos de la misma categoría y C(1) ≈ 0, es decir, sin correlación para cualquier Ω.

## A2 — Los dominios de las figuras eran un artefacto del suavizado

`plot_site_distribution.py` aplica `gaussian_filter` + `argmax` antes de graficar.

| Red | % de vecinos de la misma categoría |
|---|---|
| Aleatoria L=200, datos crudos | 33.5% |
| La misma, suavizada con σ=1 | 68.5% |
| La misma, suavizada con σ=2 | 83.2% |
| `correlation_avg/omega_0.9/red_colores_rep1.csv`, cruda | 34.9% (C(1) = +0.03) |
| La misma, suavizada con σ=2 | 82.0% |

## A4 — El cruce no conserva el inventario de enlaces de cada cromosoma

El cruce intercambia enlaces de la misma categoría pero de **distinto valor** entre dos cromosomas, y la
selección elitista favorece al que recibió enlaces más chicos.

| Medida (Ubuntu, 300 gen) | GA original | Ablación |
|---|---|---|
| Cromosomas con el mismo multiconjunto de enlaces que M_base | **0 de 16** | **0 de 16** |
| Media de enlaces: M_base → población final | 399.93 → 395.70 | 399.93 → **360.56** |
| Cromosomas con el mismo multiconjunto de sitios que M_base | 16 de 16 | 16 de 16 |

La ablación de la Tabla 3 del artículo "converge más rápido" en parte porque sus enlaces se encogen ~10%,
lo que cambia el Ω efectivo. La Tabla 3 solo reportaba la media de sitios.

## B1 / B2 — Generador aleatorio

| Medida | Windows (MinGW) | Linux (glibc) |
|---|---|---|
| `RAND_MAX` | 32 767 | 2 147 483 647 |
| `rand_r(&semilla)` | ignora la semilla (llama a `rand()`) | usa la semilla |
| Dos mutaciones desde el mismo estado coinciden en | 2 003 / 2 500 posiciones | **2 500 / 2 500** (mutación idéntica en cada generación) |
| Sitios que regresan a su valor tras 2 mutaciones seguidas | n/a | 204 de 245 |

Consecuencias:
- En Linux la mutación repite los mismos intercambios en cada generación y deshace la mayoría.
- En Windows, con L ≥ 182, el Fisher–Yates no es uniforme, y el cruce solo alcanza las primeras 32 768
  posiciones de cada categoría (en L = 500, solo las filas 0–196).
- El mismo código da resultados distintos según la plataforma.

## A3 — Discrepancias entre el artículo y los logs del repositorio

| Caso | Artículo | Logs del repositorio |
|---|---|---|
| GA, Ω = 0.9, L = 50 | 6 ± 3 errores | ~2 443 errores (`results/ga_vs_mc/log_ga_L50.txt`); Fig. 8 ≈ 2 650 |
| GA, L = 20 | 0 errores | 362 y 401 sin converger |
| MC, L = 20 | 3 errores, 2.10 s | 0 errores, 6.95 s |
| MC, L = 50 | 178 errores, 2.03 s | 136 errores, 12.71 s |
| C(r) en r = 1 | ≈ 0.85 para todo Ω (Fig. 14) | ≈ 0 en los CSV del repositorio |

## Otros defectos del código original (revisión estática)

- **B3:** `mejorCromosoma` vale 0 salvo en las generaciones con mejora. La inyección de diversidad protege la
  ranura 0 en lugar del mejor, y la relajación se aplica al cromosoma 0.
- **B4:** la mutación también altera al mejor de la población (no hay elitismo real).
- **B5:** `relajaSitios` modifica radios (rompe el inventario) sin actualizar la categoría.
- **B6:** se lee `argv[2]` aunque se ejecute sin argumentos (comportamiento indefinido).
