# EVIDENCIA 01 — ¿Aporta algo el cruce? GA v2 frente a variantes sin cruce

- **Fecha de ejecución:** 2026-11-11 (22:06, hora local)
- **Código:** `src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c`
  - SHA-256: `4946defb0e4bd3dfccc7b7ee12e07157a119fee28be95d53fc878bfb3bdb747e`
  - Base git: `3a50809` + cambios sin commit (v2 todavía no está en el repositorio)
- **Datos crudos:** `results/experimento_cruce/` (`resumen.csv`, `trazas/`, `logs/`, `redes/`, `entorno.txt`)
- **Reproducir:**
  ```bash
  bash scripts/experimento_cruce.sh                       # ~15 min con 8 núcleos
  python3 analysis/resumen_cruce.py results/experimento_cruce
  ```

## 1. Pregunta

Con la mutación guiada de v2 (que acepta un intercambio solo si no aumenta las violaciones), cada cromosoma
ya realiza una búsqueda local por sí mismo. ¿La población y el cruce PMX por bloques mejoran el resultado
respecto a la misma búsqueda local **con el mismo esfuerzo de cómputo**?

## 2. Diseño

| Elemento | Valor |
|---|---|
| Red | L = 50 (2 500 sitios, 10 000 incidencias sitio–enlace), periódica, coordinación 4 |
| Distribuciones | Uniformes de semiancho a = 50. Ω = 0.6: μS=440, μB=400. Ω = 0.9: μS=410, μB=400 |
| Semillas | 1–10 por cada (Ω, variante): 80 corridas en total |
| Presupuesto | **4×10⁹ evaluaciones de error de sitio** por corrida (idéntico para todas las variantes) |
| Paro | Error = 0 o presupuesto agotado (sin límite de generaciones) |
| Hilos | 1 por corrida (8 corridas simultáneas), así que el tiempo también es comparable |
| Mutación | Guiada, pm = 0.15, sin restricción de categoría (sitios y enlaces) |
| Hardware | Intel Core i7-7700HQ, WSL2 Ubuntu 20.04, gcc 9.4.0 `-O3` |

**Medida de costo.** Una "evaluación" es una llamada a `errSitio()`, que revisa los 4 enlaces de un sitio.
Es la operación básica de todos los operadores:
- Evaluar un cromosoma completo cuesta N evaluaciones.
- Un intercambio de sitios en la mutación cuesta 4.
- Un intercambio de enlaces cuesta 2 × (sitios afectados).

Es una medida independiente del hardware y del número de hilos.

**Variantes:**

| Variante | Población | Cruce | Comentario |
|---|---|---|---|
| `GA_P64` | 64 | PMX por bloques + selección elitista 2 de 4 | método propuesto (configuración del artículo) |
| `GA_P16` | 16 | igual | población más pequeña |
| `sinCruce_P64` | 64 | ninguno | 64 búsquedas locales independientes (aísla el efecto del cruce) |
| `busquedaLocal_P1` | 1 | ninguno | una sola búsqueda local. Equivale al Monte Carlo "acepta si no empeora", con los mismos movimientos y el mismo generador aleatorio |

Todas las variantes conservan exactamente el inventario de sitios y enlaces. Se verificó en las 80 corridas
(`inventario_ok = 1`, y el fitness recalculado coincide con el registrado).

## 3. Resultados

### 3.1 Tabla resumen (media ± DE sobre 10 semillas)

| Ω | Variante | Error final | Mediana | Error % incidencias | Llegaron a 0 | % misma categoría | C(1) | Evaluaciones usadas | Tiempo (s) |
|---|---|---|---|---|---|---|---|---|---|
| 0.6 | GA P=64 | 0.0 ± 0.0 | 0 | 0.000 | 10/10 | 43.7 ± 1.0 | 0.279 ± 0.017 | 9.2×10⁷ | 4.3 ± 0.5 |
| 0.6 | GA P=16 | 0.0 ± 0.0 | 0 | 0.000 | 10/10 | 43.8 ± 0.6 | 0.285 ± 0.014 | 2.5×10⁷ | 1.1 ± 0.1 |
| 0.6 | P=64 sin cruce | 0.0 ± 0.0 | 0 | 0.000 | 10/10 | 43.3 ± 1.0 | 0.274 ± 0.016 | 6.9×10⁷ | 3.2 ± 0.4 |
| 0.6 | Búsqueda local P=1 | 0.0 ± 0.0 | 0 | 0.000 | 10/10 | 44.2 ± 0.8 | 0.293 ± 0.018 | **1.5×10⁶** | **0.2 ± 0.1** |
| 0.9 | GA P=64 | 250.2 ± 48.3 | 230 | 2.50 ± 0.48 | 0/10 | 73.4 ± 0.8 | 0.686 ± 0.029 | 4×10⁹ | 156.5 ± 21.8 |
| 0.9 | GA P=16 | 77.0 ± 42.1 | 64 | 0.77 ± 0.42 | 0/10 | 80.7 ± 0.9 | 0.854 ± 0.033 | 4×10⁹ | 152.2 ± 20.9 |
| 0.9 | P=64 sin cruce | 211.3 ± 44.4 | 196 | 2.11 ± 0.44 | 0/10 | 74.9 ± 0.7 | 0.719 ± 0.023 | 4×10⁹ | 196.5 ± 15.5 |
| 0.9 | Búsqueda local P=1 | **0.0 ± 0.0** | **0** | **0.000** | **10/10** | 84.6 ± 1.5 | 0.943 ± 0.009 | **1.4×10⁹** | **57.8 ± 41.1** |

"Evaluaciones usadas" es la media. Las corridas que llegan a 0 se detienen antes de agotar el presupuesto.

### 3.2 Pruebas estadísticas (Mann-Whitney U, dos colas, n = 10 contra 10)

| Ω | Comparación | Mediana GA_P64 | Mediana otra | U | p | Conclusión (α = 0.05) |
|---|---|---|---|---|---|---|
| 0.9 | GA_P64 vs GA_P16 | 230 | 64 | 100.0 | 1.8×10⁻⁴ | GA_P16 mejor |
| 0.9 | GA_P64 vs sinCruce_P64 | 230 | 196 | 75.5 | 0.059 | sin diferencia significativa (el cruce no ayuda) |
| 0.9 | GA_P64 vs búsquedaLocal_P1 | 230 | 0 | 100.0 | 6.4×10⁻⁵ | búsqueda local P=1 mejor |

En Ω = 0.6 todas las variantes llegan a 0, así que el error final no las distingue. La diferencia está
en el **costo para llegar a 0**: la búsqueda local necesita ~60 veces menos evaluaciones que el GA P=64
(1.5×10⁶ frente a 9.2×10⁷).

### 3.3 Figuras

- `results/experimento_cruce/convergencia_omega_0.900.png`: mediana del mejor error contra las evaluaciones (banda = rango intercuartil).
- `results/experimento_cruce/convergencia_omega_0.600.png`: lo mismo para Ω = 0.6.
- `results/experimento_cruce/cajas_error.png`: distribución del error final.

En Ω = 0.9 la curva de la búsqueda local P=1 está por debajo de las tres variantes con población **durante
todo el presupuesto**, y llega a 0 alrededor de 10⁹ evaluaciones. Las curvas de GA_P64 y sinCruce_P64
prácticamente se superponen.

## 4. Interpretación

1. **El cruce no aporta.** Con P = 64, quitar el cruce no cambia significativamente el resultado
   (p = 0.059), y la mediana incluso es algo mejor sin cruce (196 frente a 230).
2. **La población diluye el presupuesto.** El rendimiento empeora conforme crece P (P=1 > P=16 > P=64):
   con el mismo costo total, cada cromosoma recibe 1/P del esfuerzo de búsqueda local, y el cruce no
   compensa esa pérdida.
3. **La mejor estrategia observada es la búsqueda local simple** con los movimientos de v2
   (intercambios internos que conservan el inventario, aceptados si no empeoran). Resuelve Ω = 0.9 sin
   violaciones en 10/10 semillas, usando ~35% del presupuesto.
4. **La organización espacial depende del método y del error residual**, no solo de Ω. En Ω = 0.9 el %
   de vecinos de la misma categoría va de 73% (GA P=64, con 2.5% de violaciones) a 85% (P=1, sin
   violaciones). Esto refuerza la necesidad de un criterio de paro y equilibrado común (punto 2 del plan
   de mejoras) antes de reportar correlaciones.

**Hipótesis (no verificadas aún) de por qué el cruce no ayuda:**
- Al insertar un bloque de otro padre se crean violaciones en la frontera del bloque. El hijo casi nunca
  supera a los padres y la selección 2 de 4 lo descarta.
- La inyección de diversidad mete permutaciones aleatorias (~4 000 errores) que no sobreviven.

Se pueden comprobar midiendo la tasa de hijos que entran a la población.

## 5. Implicaciones para el artículo

- **Tal como está, el artículo no puede afirmar que el GA (población + cruce) sea la razón del buen
  desempeño.** Un revisor que haga esta misma comparación obtendrá este resultado.
- **Opciones:**
  1. **Reformular la contribución** como un operador de búsqueda que conserva el inventario (intercambios
     internos guiados por el Principio de Construcción), con el GA como uno de los marcos evaluados.
     Reportar honestamente que la búsqueda local simple es la más eficiente en estas condiciones.
  2. **Rediseñar el cruce** para que sí aporte (reparación local de la frontera del bloque, cruce solo
     entre individuos ya optimizados, o un esquema memético con pocos individuos) y repetir este
     experimento como prueba.
  3. **Buscar el régimen donde la población ayude:** L mayores, Ω > 0.9, o distribuciones con menor
     margen, donde la búsqueda local pueda estancarse en mínimos locales.
- Esta comparación (mismo presupuesto, varias semillas, prueba estadística) debería quedar en el artículo
  como **tabla de ablación**, sustituyendo a la Tabla 3 actual.

## 6. Limitaciones

- Una sola L (50) y dos Ω. Falta L = 100–500 y Ω cercanos a 1.
- Un solo presupuesto (4×10⁹). La conclusión podría cambiar con presupuestos mucho mayores si la
  búsqueda local se estancara, aunque aquí resolvió todo antes del 35% del presupuesto.
- Parámetros del GA sin ajustar (pm = 0.15, bloque L/8–L/2, estancamiento = 50).
- Hardware distinto al del artículo (i7-7700HQ frente a i9-13900HX). No afecta la comparación en
  evaluaciones, solo los tiempos absolutos.
