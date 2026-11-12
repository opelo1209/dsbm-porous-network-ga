# EVIDENCIA 04: recocido simulado y poblaciones por tesela (v3)

**Fecha:** 2026-11-12. **Estado:** completo (calibración en L = 50 y piloto en L = 200, n = 2–3; exploratorio).

## Motivación

Las evidencias 01–03 muestran que la población de redes completas no aporta: la búsqueda local P=1 gana en todos los
casos. En el régimen difícil (L=200, Ω=0.95) la búsqueda local termina con ~2.7% de violaciones. Al diseñar este
experimento se supuso que se estancaba; la calibración mostró después que solo converge lentamente (ver abajo). Las
dos hipótesis iniciales, cada una con su remedio:

1. **Mínimo local.** La búsqueda local solo acepta movimientos que no empeoran. → **Recocido simulado.**
2. **La población de réplicas desperdicia la localidad de la función objetivo.** Cada sitio depende solo de sus 4
   enlaces, así que la red es casi separable por regiones. Una población de redes completas reparte el presupuesto
   entre réplicas, y su cruce pega bloques en contextos ajenos (>98% de hijos rechazados). → **Poblaciones por
   tesela** (propuesta del autor, 2026-11-12): una sola configuración, partida en teselas, cada tesela con su propia
   población.

## Implementación

### Recocido simulado en v2 (`--recocido=T0,Tf`)

- Criterio de Metropolis en la mutación guiada: un intercambio que empeora en δ se acepta con probabilidad
  exp(−δ/T).
- La temperatura baja geométricamente de T0 a Tf según la fracción de presupuesto usada:
  T = T0·(Tf/T0)^(evals/presupuesto).
- El mejor histórico se guarda aparte, así que el resultado reportado es el mejor visitado, revisado al final de
  cada generación.
- **Regresión:** sin `--recocido` no se consume ningún número aleatorio extra. Se verificó salida idéntica bit a bit
  a la versión anterior en tres configuraciones (GA PMX P=16, búsqueda local P=1, GA rango+reparación P=4; L=40,
  semilla 7, presupuesto 3e7).

### v3: poblaciones por tesela (`src/genetico_v3/ConstructorRedes2D_C4_Teselas_v3.c`)

Una sola configuración L×L que evoluciona por **épocas**:

1. **Cortar** en T×T teselas de lado S, con un corrimiento de rejilla (oi, oj). Cada tesela recibe K versiones,
   copia de su pedazo. Las versiones de una tesela comparten exactamente el mismo inventario, guardado como
   permutación de ids locales.
2. **Evolucionar** durante `--epoca` generaciones. En cada generación y para cada tesela:
   - se evalúan sus versiones en el contexto de las mejores versiones instaladas de las vecinas (coevolución
     cooperativa);
   - cruce PMX de un sub-bloque entre parejas de versiones, con selección elitista 2 de 4;
   - mutación guiada dentro de la tesela (con recocido opcional);
   - se instala la mejor versión.
3. **Migrar** con intercambios guiados entre un sitio o enlace y otro de una tesela vecina (`--migracion=f`
   barridos).
4. **Desplazar** la rejilla S/2 en diagonal, para que las fronteras no queden fijas, y volver a 1.

Detalles que garantizan la corrección:

- **Costo de una tesela:** errores de sus sitios más los de la columna oeste y la fila norte vecinas, que son los
  únicos otros sitios que dependen de sus enlaces. El delta global de cambiar una tesela es exacto. En cada época se
  compara el error acumulado por deltas con un recálculo completo. Si difieren se marca `fitness_ok=0`; no ha
  ocurrido en ninguna prueba.
- **Paralelismo determinista:** las teselas se procesan en 4 fases, con coloreo 2×2 por paridad. Dos teselas del
  mismo color no comparten sitios, enlaces ni incidencias. Verificado: con 1 y 4 hilos se obtiene el mismo resultado
  en 5 configuraciones.
- **Comparabilidad:** con la misma semilla, M_base y la configuración inicial son idénticas a las de v2. El costo se
  mide igual, en llamadas a errSitio().
- **Tasa de aceptación:** se reportan dos. La de v2 (hijo que entra a la población) y una estricta (hijo que entra
  *y* es mejor que el peor padre). Al inicio de cada época las versiones son idénticas, y los clones inflan la
  primera tasa.

Con K=1 el método es búsqueda local restringida a cada tesela, más la migración y el desplazamiento de rejilla. Esto
aísla el efecto de la población por tesela.

### Infraestructura

- `scripts/experimento_cruce.sh`: una variante usa v3 si incluye `--motor=v3`. Se compilan los dos binarios y
  `entorno.txt` registra los sha256 de ambos fuentes.
- `analysis/resumen_cruce.py`: la referencia de las pruebas de Mann-Whitney se elige con
  `REFERENCIA=<variante>`. Por defecto sigue siendo GA_P64.

## Calibración (L=50, Ω=0.95, presupuesto 1e9, semillas 1–2)

Datos: `results/calibracion_v3/`. Con n = 2 no hay significancia posible (p mínimo = 0.33); sirve para elegir
parámetros.

| Variante | Errores (s1 / s2) | Mediana | % misma cat. | C(1) | Tasa estricta de hijos |
|---|---|---|---|---|---|
| **Búsqueda local P=1** | 191 / 202 | **196** | 84.2 | 0.827 | — |
| Recocido T0=0.3, Tf=0.01 | 204 / 225 | 214 | 85.0 | 0.860 | — |
| v3 K=1 (S=25, 2×2 teselas) | 205 / 231 | 218 | 83.0 | 0.808 | — |
| Recocido T0=1, Tf=0.01 | 228 / 249 | 238 | 84.0 | 0.821 | — |
| v3 K=4 (S=25) | 318 / 327 | 322 | 77.6 | 0.685 | 0.28% |

### Hallazgo: la búsqueda local no está atascada, converge lentamente

Mejor error de la búsqueda local P=1 en función de la fracción del presupuesto usada:

| Fracción | L=50 (calibración, s1) | L=200 (evidencia 02, s1) |
|---|---|---|
| 0.10 | 411 | 7 547 |
| 0.25 | 312 | 5 972 |
| 0.50 | 248 | 5 101 |
| 0.75 | 208 | 4 649 |
| 1.00 | 191 | 4 367 |

Sigue bajando al agotar el presupuesto. **Corrige lo que se dijo antes ("la búsqueda local se estanca en 2.7%"):**
el residuo en Ω=0.95 lo limita el presupuesto, no un mínimo local. Esto explica que el recocido no ayude: gasta la
primera parte del presupuesto explorando (1 641 errores al 5%, contra 507 de la búsqueda local), y al final apenas
recupera el terreno perdido. El cuello de botella es la **eficiencia de cada propuesta**. Cerca de la factibilidad,
casi todos los intercambios al azar involucran sitios sin violaciones y se rechazan.

### Interpretación preliminar

- El recocido, con estas temperaturas, no mejora. Una temperatura más baja lo acercaría a la búsqueda local, pero
  no lo haría mejor.
- v3 K=1 es comparable a la búsqueda local, ligeramente peor. Restringir los intercambios a la tesela no estorba
  mucho.
- v3 K=4 repite el patrón de la población de redes completas: repartir el presupuesto entre K versiones cuesta más
  de lo que aporta la diversidad. Solo 0.28% de los hijos mejora al peor padre.
- Con L=50 solo caben 2×2 teselas de 25, así que la prueba real de v3 es L=200 (4×4 teselas de 50).

## Mutación focalizada (`--focal=f`, v2)

Responde al hallazgo anterior: si el cuello de botella es la eficiencia de cada propuesta, conviene proponer donde
hay violaciones.

- Cada configuración mantiene el conjunto de sitios con al menos una violación (lista e índice inverso). Se
  actualiza en O(1) con los errores que el intercambio aceptado ya calculó, así que no hay evaluaciones extra.
  Reconstruirlo tras un cruce o una copia cuesta N evaluaciones, y se cobran.
- Mismo número de intentos por generación que la mutación normal (pm·N de sitios y pm·2N de enlaces). Con
  probabilidad f, el primer elemento sale del conjunto de violados: el sitio mismo, o uno de sus 4 enlaces al azar.
  La pareja sigue siendo uniforme en toda la red, y la aceptación es la misma (guiada).
- **Regresión:** con f = 0 el resultado es idéntico bit a bit en 4 configuraciones (GA PMX, búsqueda local, GA
  rango+reparación, recocido). Con f > 0 es determinista entre 1 y 4 hilos.

Calibración (L=50, Ω=0.95, 1e9, semillas 1–2; misma tabla que arriba):

| Variante | Errores (s1 / s2) | Mediana | % misma cat. | C(1) | Tiempo (s) |
|---|---|---|---|---|---|
| Búsqueda local P=1 (f = 0) | 191 / 202 | 196 | 84.2 | 0.827 | 73.7 |
| **Focal f = 0.5** | 108 / 132 | **120** | 86.0 | 0.881 | 41.6 |
| Focal f = 0.9 | 125 / 133 | 129 | 84.7 | 0.841 | 40.8 |
| Focal f = 1.0 | 189 / 213 | 201 | 79.0 | 0.763 | 39.9 |

- Con f = 0.5 el residuo baja ~39% con el mismo presupuesto, y la corrida tarda ~44% menos. Las propuestas
  focalizadas tocan menos memoria al azar.
- f = 1 deja de funcionar: hacen falta también intercambios entre sitios sin violaciones para "hacer espacio" (mover
  radios grandes hacia donde se necesitan). El óptimo está en una mezcla, alrededor de f ≈ 0.5–0.9.
- La correlación final sube con el mejor resultado (C(1) = 0.881 con f = 0.5), igual que en las evidencias
  anteriores: la correlación acompaña a la satisfacción del Principio de Construcción.

## Piloto (L=200, Ω=0.95, presupuesto 2e10, semillas 1–3)

Datos: `results/piloto_v3_L200/`. El recocido se descartó tras la calibración. La búsqueda local P=1 son las corridas
de la evidencia 02 (mismas semillas, M_base y presupuesto; código idéntico bit a bit con las opciones por defecto).
v3 usa S = 50 (4×4 teselas), épocas de 20 generaciones y migración 0.15.

| Variante | Mediana | Media ± DE | Error % | % misma cat. | C(1) | vs P=1 |
|---|---|---|---|---|---|---|
| Búsqueda local P=1 | 4 328 | 4 266 ± 142 | 2.666 | 86.1 | 0.842 | — |
| **Focal f = 0.75** | **3 103** | 3 026 ± 158 | **1.891** | 87.1 | 0.864 | **−28%**; 3/3 semillas mejores (U = 9) |
| Focal f = 0.5 | 3 373 | 3 301 ± 162 | 2.063 | 87.2 | 0.860 | −22%; 3/3 semillas mejores (U = 9) |
| v3 K=1 | 4 607 | 4 510 ± 174 | 2.819 | 85.4 | 0.825 | +6% |
| v3 K=2 | 5 264 | 5 240 ± 147 | 3.275 | 83.1 | 0.774 | +22%; tasa estricta de hijos 0.1% |

Con n = 3 contra 3 el p mínimo posible es 0.1, así que no hay significancia formal. Pero en la focal **las 3
corridas son mejores que las 3 de la búsqueda local** (U máximo), y la diferencia es de 5–8 desviaciones estándar.

Convergencia (mejor error, semilla 1, por fracción del presupuesto):

| Fracción | 0.10 | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|---|
| Búsqueda local P=1 | 7 547 | 5 972 | 5 101 | 4 649 | 4 367 |
| Focal f = 0.5 | 6 652 | 5 069 | 4 136 | 3 671 | 3 414 |
| v3 K=1 | 8 171 | 6 398 | 5 403 | 4 931 | 4 614 |

La focal va por delante durante toda la corrida: con el 50% del presupuesto ya está por debajo del resultado final
de la búsqueda local.

## Conclusiones de la evidencia 04

1. **La mutación focalizada es la mejora más clara encontrada hasta ahora:** −39% en L = 50 y −28% en L = 200 con el
   mismo presupuesto, sin costo adicional por evaluación. El mejor f en L = 200 es 0.75 (en L = 50 fue 0.5). Hay que
   barrer f en el experimento formal.
2. **El recocido no ayuda** con este presupuesto: el problema no son los mínimos locales, sino la eficiencia de las
   propuestas.
3. **Las poblaciones por tesela (v3) no superan a la búsqueda local.** K = 1 queda 6% peor: restringir los
   intercambios a la tesela cuesta un poco. K = 2 queda 22% peor: la población por tesela repite el problema de la
   dilución, y solo 0.1% de los hijos mejora al peor padre. La conclusión del escenario B se mantiene: en este
   problema la población no aporta, ni de redes completas ni de teselas.
4. **Lo que v3 sí aporta es una vía de paralelización determinista de una sola red.** El coloreo 2×2 permite procesar
   las teselas en paralelo con resultados independientes del número de hilos, a un costo de calidad de ~6% con K = 1.
   Para redes muy grandes (L ≥ 500), una v3 con K = 1 y mutación focalizada podría ser la forma de usar varios núcleos
   en una sola red. Falta medir su escalamiento. La referencia es la búsqueda local P=1 de la evidencia 02: mismos logs,
mismas semillas, M_base y presupuesto, y código idéntico bit a bit con las opciones por defecto.
