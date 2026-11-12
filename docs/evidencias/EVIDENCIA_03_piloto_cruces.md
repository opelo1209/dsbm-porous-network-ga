# EVIDENCIA 03: piloto de variantes de cruce (opción 2 del plan)

**Fecha:** 2026-11-11. **Estado:** piloto completo (n = 3 por variante). Es exploratorio: con n = 3 la prueba de
Mann-Whitney no puede dar p < 0.1, así que aquí no se reportan pruebas. Sirve para elegir qué variantes entran a un
experimento formal.

## Pregunta

¿Se puede rediseñar el cruce del GA v2 para que aporte frente a la búsqueda local con el mismo presupuesto? Se
probaron dos cambios:

1. **Cruce por rangos** (`--cruce=rango`). Dentro del bloque, el hijo conserva sus propios valores pero los ordena
   según el orden relativo del donante. Transmite el *patrón espacial* (qué posiciones tienen radios grandes o
   chicos) sin mover inventario entre cromosomas. El PMX por bloques, en cambio, transmite valores concretos y
   desplaza genes fuera del bloque.
2. **Reparación de frontera** (`--reparar=R`). Después del cruce se hacen R barridos de búsqueda local guiada sobre
   un anillo de ancho 1 dentro y fuera del bloque y sobre las posiciones desplazadas por el PMX. Su costo se cobra al
   presupuesto.

También se probó una población más chica (P = 4), porque la evidencia 01 mostró que P = 16 es mejor que P = 64.

## Diseño

| Parámetro | Valor |
|---|---|
| Código | `src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c`, sha256 `ac9d6cad…89f63610`, commit base `3a50809` con cambios sin commit |
| Red | L = 50 (N = 2 500 sitios, 5 000 enlaces) |
| Distribuciones | Uniformes, μS = 410, μB = 400, semiancho 50 (Ω = 0.9) |
| Presupuesto | 1 × 10⁹ evaluaciones de sitio por corrida (igual para todas) |
| Semillas | 1, 2, 3 |
| Hilos | 1 por corrida, 3 corridas simultáneas |
| Hardware | Intel i7-7700HQ, WSL2 Ubuntu 20.04, gcc 9.4.0 `-O3` |

Variantes:

| Nombre | P | Cruce | Reparación |
|---|---|---|---|
| GA16_pmx | 16 | PMX por bloques | 0 |
| GA16_pmx_rep2 | 16 | PMX por bloques | 2 barridos |
| GA16_rango | 16 | rangos | 0 |
| GA16_rango_rep2 | 16 | rangos | 2 barridos |
| GA4_rango_rep2 | 4 | rangos | 2 barridos |
| busquedaLocal_P1 | 1 | sin cruce | — |

Para reproducirlo:

```bash
L=50 PRESUPUESTO=1e9 SEMILLAS="1 2 3" OMEGAS="0.9:410:400:50" \
VARIANTES="GA16_pmx:16 GA16_pmx_rep2:16:--reparar=2 GA16_rango:16:--cruce=rango \
GA16_rango_rep2:16:--cruce=rango,--reparar=2 GA4_rango_rep2:4:--cruce=rango,--reparar=2 \
busquedaLocal_P1:1:--sin-cruce" OUT=results/piloto_cruce TRABAJOS=3 bash scripts/experimento_cruce.sh
python3 analysis/resumen_cruce.py results/piloto_cruce
```

Datos crudos: `results/piloto_cruce/` (resumen.csv, logs/, trazas/, convergencia_omega_0.900.png, cajas_error.png).

## Resultados

Errores tipo 1 al agotar el presupuesto (de 10 000 incidencias sitio–enlace):

| Variante | Mediana | Media ± DE | Error % | Llegaron a 0 | % misma categoría | C(1) | Tasa de aceptación de hijos |
|---|---|---|---|---|---|---|---|
| **busquedaLocal_P1** | **0** | 13.3 ± 23.1 | 0.13 | **2/3** | 83.4 | 0.925 | — |
| GA4_rango_rep2 | 92 | 126.7 ± 68.0 | 1.27 | 0/3 | 77.5 | 0.793 | 1.59–1.64% |
| GA16_rango | 193 | 225.7 ± 64.5 | 2.26 | 0/3 | 73.6 | 0.701 | 0.72–0.83% |
| GA16_pmx | 197 | 229.0 ± 59.8 | 2.29 | 0/3 | 73.7 | 0.704 | 0.55–0.57% |
| GA16_rango_rep2 | 263 | 286.3 ± 67.1 | 2.86 | 0/3 | 71.5 | 0.650 | 1.08–1.10% |
| GA16_pmx_rep2 | 423 | 458.0 ± 72.2 | 4.58 | 0/3 | 64.3 | 0.529 | 0.43–0.51% |

La búsqueda local que llegó a 0 se detuvo antes de agotar el presupuesto (6.2 × 10⁸ evaluaciones en promedio).

La **tasa de aceptación** es la fracción de hijos que sobreviven a la selección elitista 2 de 4 (padres + hijos).
Se lee de la línea `CRUCE,...` de cada log.

## Interpretación

1. **Ninguna variante con cruce se acerca a la búsqueda local.** La mejor (P = 4, rangos, reparación) deja una
   mediana de 92 errores, contra 0 de P = 1.
2. **Cruce por rangos ≈ PMX** (193 contra 197). Conservar el patrón espacial en lugar de los valores no cambia el
   resultado. La aceptación sube un poco (0.8% contra 0.56%), pero sigue siendo marginal.
3. **La reparación de frontera empeora con P = 16** (rangos: 193 → 263; PMX: 197 → 423). El costo de reparar se
   come el presupuesto que, sin cruce, se gastaría en mutación guiada sobre toda la red.
4. **Reducir la población es lo que más ayuda** (P = 16 → 4: 263 → 92). La tendencia es monótona: P = 64 > 16 > 4 >
   1 (evidencia 01 y este piloto). Mientras menos presupuesto se reparte entre individuos, mejor. Esto apunta a que la
   población **diluye** el esfuerzo de búsqueda local en vez de aportar diversidad útil.
5. **Por qué el cruce no ayuda:** más del 98% de los hijos se descartan. Un hijo hereda un bloque de un padre que,
   en su contexto, tenía pocas violaciones. Al pegarlo en otro contexto se crean violaciones en la frontera más
   rápido de lo que se arreglan. La función objetivo es local (cada sitio depende de sus 4 enlaces) y no hay
   "bloques constructivos" separables que el cruce pueda recombinar.
6. **La correlación final acompaña al error:** las variantes con menos errores tienen más vecinos de la misma
   categoría y mayor C(1). Es decir, la correlación es consecuencia de satisfacer el Principio de Construcción, no
   de un artefacto del operador.

## Implicaciones para el artículo

- Respalda el **escenario B** del plan de reformulación (`docs/PLAN_REFORMULACION_ARTICULO.md`): la contribución es
  la representación que conserva el inventario y el movimiento guiado, no la población.
- La tasa de aceptación < 2% es el dato concreto para explicar *por qué* (sección 5.3 propuesta).
- El cruce por rangos y la reparación pueden describirse como "intentos razonables de adaptar el cruce al problema"
  que no cambian la conclusión. Esto fortalece el estudio frente a la objeción "usaron un mal cruce".

## Limitaciones

- n = 3: piloto exploratorio, sin pruebas estadísticas. Las diferencias grandes (P1 contra cualquier GA; P = 4
  contra P = 16 con reparación) son consistentes en las tres semillas, pero hay que confirmarlas con n ≥ 10 si se
  van a reportar.
- Solo L = 50, Ω = 0.9. El régimen grande (L = 200, Ω = 0.95) está en la evidencia 02.
- Un solo ancho de reparación (2 barridos) y un solo rango de tamaños de bloque (L/8 a L/2).
- Las 3 corridas que faltaban (fallo de `xargs`, ver evidencia 01) se completaron con un binario recompilado del
  mismo fuente; con las opciones por defecto el resultado es idéntico bit a bit.
