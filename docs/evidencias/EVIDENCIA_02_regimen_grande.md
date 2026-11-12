# EVIDENCIA 02: régimen grande (L = 200, Ω = 0.9 y 0.95)

**Fecha:** 2026-11-12. **Estado:** completo (5 semillas × 4 variantes × 2 Ω = 40 corridas).

## Pregunta

La evidencia 01 (L = 50) mostró que la población y el cruce no aportan. ¿Se mantiene en una red 16 veces más grande
y con un traslape más alto (Ω = 0.95), donde la búsqueda local podría tener más dificultades y la diversidad de una
población podría servir? Este experimento decide entre los escenarios A y B del plan de reformulación.

## Diseño

| Parámetro | Valor |
|---|---|
| Código | `src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c`, sha256 `ac9d6cad…89f63610` (commit base `3a50809` con cambios sin commit). Las versiones posteriores (recocido, focal) dan resultados idénticos bit a bit con estas opciones |
| Red | L = 200 (N = 40 000 sitios, 80 000 enlaces) |
| Distribuciones | Uniformes, semiancho 50. Ω = 0.9: μS = 410, μB = 400. Ω = 0.95: μS = 405, μB = 400 |
| Presupuesto | 2 × 10¹⁰ evaluaciones de sitio por corrida (5 × 10⁵ por sitio, similar al de la evidencia 01) |
| Semillas | 1–5 |
| Hilos | 1 por corrida, 8 corridas simultáneas |
| Hardware | Intel i7-7700HQ, WSL2 Ubuntu 20.04, gcc 9.4.0 `-O3` |

Variantes:

| Nombre | P | Operadores |
|---|---|---|
| GA_P64 | 64 | PMX por bloques + mutación guiada |
| GA_P16 | 16 | PMX por bloques + mutación guiada |
| GA4_rango_rep2 | 4 | Cruce por rangos + reparación de frontera (2 barridos) + mutación guiada. Es la mejor variante con cruce del piloto (evidencia 03) |
| busquedaLocal_P1 | 1 | Solo mutación guiada |

Para reproducirlo:

```bash
L=200 PRESUPUESTO=2e10 SEMILLAS="1 2 3 4 5" \
OMEGAS="0.95:405:400:50 0.9:410:400:50" \
VARIANTES="GA_P64:64 GA_P16:16 GA4_rango_rep2:4:--cruce=rango,--reparar=2 busquedaLocal_P1:1:--sin-cruce" \
OUT=results/experimento_regimen_grande TRABAJOS=8 bash scripts/experimento_cruce.sh
REFERENCIA=busquedaLocal_P1 python3 analysis/resumen_cruce.py results/experimento_regimen_grande
```

Datos crudos: `results/experimento_regimen_grande/` (resumen.csv, logs/, trazas/, convergencia_omega_0.900.png,
convergencia_omega_0.950.png, cajas_error.png).

## Resultados

Errores tipo 1 al agotar el presupuesto (de 160 000 incidencias sitio–enlace). Ninguna corrida llegó a 0.

| Ω | Variante | Mediana | Media ± DE | Error % | % misma categoría | C(1) | Tasa de aceptación de hijos |
|---|---|---|---|---|---|---|---|
| 0.90 | **Búsqueda local P=1** | **35** | 49.0 ± 27.9 | **0.031** | 85.9 | 0.948 | — |
| 0.90 | GA4_rango_rep2 | 1 688 | 1 702 ± 156 | 1.064 | 80.6 | 0.832 | ~0.2% |
| 0.90 | GA_P16 | 4 060 | 4 074 ± 137 | 2.546 | 74.1 | 0.689 | — |
| 0.90 | GA_P64 | 7 633 | 7 647 ± 164 | 4.779 | 65.2 | 0.519 | 0.02–0.03% |
| 0.95 | **Búsqueda local P=1** | **4 328** | 4 342 ± 174 | **2.714** | 86.2 | 0.840 | — |
| 0.95 | GA4_rango_rep2 | 6 490 | 6 479 ± 192 | 4.049 | 81.3 | 0.739 | 0.16–0.20% |
| 0.95 | GA_P16 | 9 133 | 9 158 ± 170 | 5.724 | 74.9 | 0.625 | 0.006–0.008% |
| 0.95 | GA_P64 | 13 646 | 13 644 ± 206 | 8.527 | 65.2 | 0.468 | — |

Pruebas de Mann-Whitney U, dos colas, n = 5 contra 5, referencia búsqueda local P=1:

| Ω | Comparación | Mediana P=1 | Mediana otra | U | p | Conclusión |
|---|---|---|---|---|---|---|
| 0.90 | P=1 vs GA_P64 | 35 | 7 633 | 0 | 0.0079 | P=1 mejor |
| 0.90 | P=1 vs GA_P16 | 35 | 4 060 | 0 | 0.0079 | P=1 mejor |
| 0.90 | P=1 vs GA4_rango_rep2 | 35 | 1 688 | 0 | 0.0079 | P=1 mejor |
| 0.95 | P=1 vs GA_P64 | 4 328 | 13 646 | 0 | 0.0079 | P=1 mejor |
| 0.95 | P=1 vs GA_P16 | 4 328 | 9 133 | 0 | 0.0079 | P=1 mejor |
| 0.95 | P=1 vs GA4_rango_rep2 | 4 328 | 6 490 | 0 | 0.012 | P=1 mejor |

U = 0 en todas: **las 5 corridas de la búsqueda local son mejores que las 5 de cada GA**. p = 0.0079 es el mínimo
posible con 5 contra 5. Con Bonferroni para 6 comparaciones (α = 0.0083), cinco de las seis quedan justo bajo el
umbral. La sexta (p = 0.012, por un empate) no, aunque también tiene U = 0.

## Interpretación

1. **La búsqueda local gana en todo el rango probado, incluido el régimen difícil.** No hay ningún régimen en el que
   la población convenga. En Ω = 0.9 la diferencia es de dos órdenes de magnitud (0.03% contra 1–5%).
2. **El orden P=1 < P=4 < P=16 < P=64 es perfectamente monótono** en los dos valores de Ω, igual que en L = 50
   (evidencias 01 y 03). Repartir un presupuesto fijo entre más individuos solo diluye el esfuerzo.
3. **El cruce es todavía menos útil en redes grandes.** Con PMX sobrevive 1 hijo de cada 10 000 en L = 200 (0.006–0.03%),
   contra 0.5% en L = 50. El cruce por rangos con reparación sube a ~0.2%, pero sigue siendo marginal. Un bloque
   trasplantado crea violaciones en su frontera, y en una red casi factible eso casi nunca compensa.
4. **En Ω = 0.95 la búsqueda local no está atascada: converge lentamente.** Su mejor error baja de 5 972 (25% del
   presupuesto) a 5 101 (50%), 4 649 (75%) y 4 367 (100%) en la semilla 1, y sigue bajando al agotar el presupuesto.
   El residuo de 2.7% lo limita el presupuesto. Este hallazgo motivó la mutación focalizada de la evidencia 04.
5. **La correlación acompaña al error:** la variante con menos violaciones tiene más vecinos de la misma categoría y
   mayor C(1). Con Ω = 0.9 la búsqueda local llega a C(1) = 0.95.

## Implicaciones para el artículo

- **Decide el escenario B** del plan de reformulación (`docs/PLAN_REFORMULACION_ARTICULO.md`): la población de redes
  completas no aporta en ningún régimen probado (L = 50 y 200, Ω = 0.6, 0.9 y 0.95). Lo que funciona es la
  representación que conserva el inventario y el movimiento guiado.
- La tabla y la curva de convergencia de Ω = 0.95 son la figura central de la sección "Comparación con igual
  presupuesto".
- La caída de la tasa de aceptación con L (0.5% → 0.01%) es un argumento cuantitativo de *por qué* el cruce no sirve
  en este problema.
- La evidencia 04 explora alternativas que podrían cambiar la conclusión del artículo: el recocido (no ayuda), las
  poblaciones por tesela (v3, en prueba) y la mutación focalizada (prometedora).

## Limitaciones

- n = 5 por celda. Es suficiente para estas diferencias (U = 0 en todas), pero en el artículo conviene n ≥ 10 en
  las configuraciones finales.
- Un solo presupuesto. No se probó si con presupuestos mucho mayores el GA alcanza a la búsqueda local; las curvas
  de convergencia no lo sugieren.
- Solo distribuciones uniformes con el mismo semiancho.
- El tiempo de pared no es comparable entre variantes: 8 corridas simultáneas en 4 núcleos físicos (8 hilos lógicos)
  compiten por memoria. La comparación válida es por evaluaciones.
- 30 de las 40 corridas se completaron en una reanudación con el binario recompilado. Las opciones por defecto dan
  resultados idénticos bit a bit, verificado.
