# EVIDENCIA 07: experimento formal del método final (n = 10, L = 50 y L = 200)

**Fecha:** 2026-11-12. **Estado:** completo (2 tamaños × 4 Ω × 5 variantes × 10 semillas = 400 corridas; 16 de
ellas reutilizadas de evidencias anteriores, ver `reutilizadas.txt`). Sustituye los valores piloto (n = 2–5) de las
evidencias 04 y 05 en el manuscrito revisado (`Paper_TCJ_EA_v2`, sección 6.4).

## Pregunta

Con 10 semillas por celda:
1. ¿La mutación dirigida por conflictos (CD) mejora al (1+1)-EA con mutación uniforme (U) de forma significativa?
2. ¿Qué tan sensible es al parámetro f?
3. ¿El GA μ = 4 con la misma mutación CD alcanza al (1+1)-EA?

## Diseño

| Parámetro | Valor |
|---|---|
| Código | v2, sha256 `f36f445e…` (commit `2a9bc5a`, ver nota del README) |
| Redes | L = 50 (presupuesto 4 × 10⁹) y L = 200 (presupuesto 2 × 10¹⁰); la corrida se detiene antes si llega a f = 0 |
| Configuraciones | Ω = 0.3 (440/370), 0.6 (440/400), 0.9 (410/400), 0.95 (405/400); semiancho 50 |
| Variantes | `busquedaLocal_P1` = (1+1)-EA U (`--sin-cruce`); `focal{0.5,0.75,0.9}_P1` = (1+1)-EA CD (`--sin-cruce --focal=f`); `GA4_pmx_focal` = GA μ = 4, PMX, CD f = 0.75 (`--focal=0.75`) |
| Semillas | 1–10 |
| Hardware | i7-7700HQ, WSL2, gcc 9.4.0 `-O3`; 1 hilo por corrida, 8 simultáneas |
| Duración | 09:17 – 15:00 (L = 200); L = 50 antes |

Para reproducirlo: `bash scripts/experimento_formal.sh`. Las pruebas adicionales (evaluaciones hasta f = 0 y
comparaciones entre valores de f) se calcularon con Mann–Whitney de dos colas y, cuando las semillas son pareadas,
también con Wilcoxon de rangos con signo.

Las 400 corridas pasaron la verificación de inventario y la de fitness.

## Resultados

Datos completos: `results/formal_L50/` y `results/formal_L200/` (`tabla_resumen.md`, `pruebas_estadisticas.md`,
`resumen.csv`, logs y trazas). Mediana de violaciones finales y, cuando todas llegan a 0, mediana de evaluaciones
hasta f = 0.

### L = 200

| Ω | U | CD f = 0.5 | CD f = 0.75 | CD f = 0.9 | GA μ=4 CD |
|---|---|---|---|---|---|
| 0.3 | 0 (10/10), 4.8e6 ev | 0, 4.0e5 | 0, 3.2e5 | 0, 2.7e5 | 0, 1.6e6 |
| 0.6 | 0 (10/10), 3.1e7 ev | 0, 4.1e6 | 0, 3.2e6 | 0, 2.9e6 | 0, 1.7e7 |
| 0.9 | **38 (0/10)**, rango 28–136 | **0 (10/10)**, 5.5e9 | **0 (10/10)**, 6.5e9 | **0 (10/10)**, 8.9e9 | 922 (0/10) |
| 0.95 | 4 422 (2.77%) | 3 452 (2.18%) | **3 142 (2.01%)** | **3 124 (1.97%)** | 5 354 (3.37%) |

Pruebas contra U (Mann–Whitney, n = 10 contra 10):
- Ω = 0.9: CD gana con las tres f (U = 100, p = 6.3 × 10⁻⁵); el GA pierde (U = 0, p = 1.8 × 10⁻⁴).
- Ω = 0.95: CD gana con las tres f (U = 100, p = 1.8 × 10⁻⁴); el GA pierde (U = 0, p = 1.8 × 10⁻⁴). Reducción del
  residuo: −22% (f = 0.5), −29% (f = 0.75), −29% (f = 0.9). Las 30 corridas CD terminan por debajo del mejor U
  (4 104).
- Ω = 0.3 y 0.6: todas llegan a 0; CD llega de 8 a 18 veces más rápido (p = 1.8 × 10⁻⁴ en todas).

Comparación entre valores de f (semillas pareadas):

| Ω | Comparación | Medianas | Mann–Whitney | Wilcoxon | Mejor |
|---|---|---|---|---|---|
| 0.9 (evaluaciones) | 0.5 vs 0.75 | 5.5e9 vs 6.5e9 | p = 0.009 | p = 0.002 | f = 0.5 en 10/10 |
| 0.9 (evaluaciones) | 0.75 vs 0.9 | 6.5e9 vs 8.9e9 | p = 3.3 × 10⁻⁴ | p = 0.002 | f = 0.75 en 10/10 |
| 0.95 (violaciones) | 0.75 vs 0.5 | 3 142 vs 3 452 | p = 0.011 | p = 0.002 | f = 0.75 en 10/10 |
| 0.95 (violaciones) | 0.75 vs 0.9 | 3 142 vs 3 124 | p = 0.57 | p = 0.014 | f = 0.9 en 8/10 (≈ 1.7%) |

### L = 50

| Ω | U | CD f = 0.5 | CD f = 0.75 | CD f = 0.9 | GA μ=4 CD |
|---|---|---|---|---|---|
| 0.3 | 0, 2.1e5 ev | 0, 2.3e4 | 0, 1.8e4 | 0, 1.6e4 | 0, 9.9e4 |
| 0.6 | 0, 1.6e6 ev | 0, 2.5e5 | 0, 1.9e5 | 0, 1.8e5 | 0, 1.0e6 |
| 0.9 | 0, 1.0e9 ev | 0, 2.5e8 | 0, 3.0e8 | 0, 4.0e8 | 0, 1.4e9 |
| 0.95 | 182 (1.82%) | 113 | **100** | 106 | 168 |

- Ω = 0.95: CD mejora a U (p = 0.045, 0.038 y 0.013 para f = 0.5, 0.75 y 0.9; −38% a −45%). El GA no se distingue
  de U (p = 0.50). Entre valores de f no hay diferencias (p ≥ 0.06).
- Ω = 0.9: CD llega a 0 de 2.6 a 4.2 veces más rápido (p ≤ 0.026); el GA no se distingue de U (p = 0.34).

### Organización espacial (% de vecinos de la misma categoría, U contra CD f = 0.75)

| L | Ω = 0.3 | Ω = 0.6 | Ω = 0.9 | Ω = 0.95 |
|---|---|---|---|---|
| 50 | 34.1 / 34.3 (p = 0.82) | 44.1 / 40.6 (p = 1.8e-4) | 84.4 / 80.8 (p = 0.001) | 87.3 / 88.3 |
| 200 | 34.4 / 34.0 (p = 3e-4) | 44.6 / 41.6 (p = 1.8e-4) | 86.0 / 82.0 (p = 1.8e-4) | 86.2 / 87.0 |

## Interpretación

1. **La mutación dirigida por conflictos es la mejora decisiva, y es significativa con n = 10.** En L = 200 con
   Ω = 0.9 convierte un residuo de 28–136 violaciones en redes **sin ninguna violación en 30/30 corridas**,
   con una mediana de entre el 28% y el 44% del presupuesto. En Ω = 0.95 reduce el residuo un 29%.
2. **El método es robusto en f ∈ [0.5, 0.9]**: todas las f mejoran significativamente a U en todos los casos. Hay
   un compromiso: una f más baja llega antes a la factibilidad en Ω = 0.9, y una f más alta deja menos residuo en
   Ω = 0.95. **f = 0.75 es el valor recomendado**: nunca es claramente la peor, empata con f = 0.9 en Ω = 0.95
   (diferencia ≈ 1.7%) y llega a f = 0 un 27% antes que f = 0.9 en Ω = 0.9.
3. **La población no ayuda, ni siquiera con la mutación CD.** En L = 200 el GA μ = 4 CD es significativamente peor
   que el (1+1)-EA *uniforme* (Ω = 0.9: 922 contra 38; Ω = 0.95: +21%). En L = 50 empata con U.
4. **La correlación depende de la trayectoria** (confirma la evidencia 06): cuando ambos llegan a f = 0, CD se
   detiene antes y con menos organización espacial (Ω = 0.6 y 0.9, p ≤ 0.001). En Ω = 0.95, donde ninguno llega a
   0, CD termina con algo más de correlación.

## Implicaciones para el artículo

- Reemplazar los valores piloto de la sección 6.4 y de la Tabla de la mutación focal por estos.
- Resumen y contribuciones: "elimina todas las violaciones en L = 200, Ω = 0.9 (30/30 corridas)" y "reduce el
  residuo un 29% en Ω = 0.95"; f = 0.75 como valor por defecto, con la sensibilidad a f.
- Eliminar la nota `\pend{}` del experimento formal.

## Limitaciones

- Un solo par de tamaños (50 y 200) y una sola familia de distribuciones (uniformes de semiancho 50).
- En Ω = 0.95 ninguna variante llega a 0 con el presupuesto usado; el residuo sigue bajando al final (la búsqueda
  está limitada por presupuesto, evidencia 04).
- Las correlaciones de las corridas que llegan a f = 0 corresponden a la primera configuración factible, sin
  equilibrar.
