# Plan de reformulación del artículo (opción 1)

**Estado:** borrador de trabajo, actualizado el 2026-11-12. **El experimento 02 decide el escenario B**: la búsqueda
local gana en las 6 comparaciones de L=200 (evidencia 02). Hay dos líneas abiertas que pueden enriquecerlo (evidencia
04):
- **Mutación focalizada:** en L=50 reduce el residuo ~39%. Si se confirma en L=200, pasa a ser el método principal.
- **Poblaciones por tesela (v3):** si en L=200 superan a la búsqueda local, el artículo vuelve hacia el escenario A
  con un esquema poblacional nuevo. Si no, entran como un intento más que no cambia la conclusión.
**Revista objetivo:** *The Computer Journal* (OUP/BCS), Sección C (Computational Intelligence, Machine Learning
and Data Analytics), con la Sección A (algoritmos) como alternativa.

## 1. Por qué hay que reformular

Las evidencias 00 y 01 (`docs/evidencias/`) muestran que la versión actual no se sostiene:

| Afirmación del manuscrito actual | Lo que muestran las evidencias |
|---|---|
| Los operadores con restricción de categoría producen dominios homotácticos que crecen con Ω (Fig. 11–14, Sec. 4.7) | El mapa de categorías nunca cambia (A1). Los dominios de las figuras salían del suavizado (A2). |
| Los operadores preservan exactamente F_B(R_B) | El cruce encoge los enlaces: −1% en el GA y −10% en la ablación (A4) |
| 6 ± 3 errores en L=50, Ω=0.9; Tablas 1–3 | Los logs del repositorio no reproducen esos números (A3) |
| El GA con población y cruce es la razón del buen desempeño | Con igual presupuesto, el cruce no aporta (p = 0.059) y la búsqueda local P=1 es mejor (evidencia 01). Solo 0.5–1.6% de los hijos sobreviven a la selección (piloto) |

Mantener el texto actual expondría el artículo a rechazo, y a un problema de integridad si un revisor ejecuta el
repositorio, que es público y está citado. **Hay que retirar las afirmaciones no sostenidas y rehacer todos los
resultados con v2.**

## 2. Qué sí es una contribución sólida (en ambos escenarios)

1. **Representación como permutación del inventario.** Cada configuración es una permutación de los radios de
   M_base. Todos los operadores son intercambios internos, así que F_S y F_B se conservan *exactamente por
   construcción* (se verifica en cada corrida). Esto corrige un defecto real de los operadores "por categoría".
2. **Movimientos guiados por el Principio de Construcción con evaluación local O(1)** de sitios y enlaces, sin
   relajación de radios.
3. **Emergencia verificada de la correlación dependiente de Ω**, medida sobre datos crudos (% de vecinos de la
   misma categoría, C(r)), sin suavizado: Ω=0.3 → 35%, Ω=0.6 → 44%, Ω=0.9 → 77–86%.
4. **Protocolo de comparación justo:** presupuesto medido en evaluaciones de sitio (independiente del hardware),
   varias semillas, prueba de Mann-Whitney, reproducibilidad bit a bit entre plataformas y número de hilos.
5. **Una ablación rigurosa** que responde con datos si la población y el cruce aportan, y en qué régimen.

## 3. Escenario A — el GA gana en el régimen difícil

**Condición:** en L=200, Ω=0.95, alguna variante con población (P=4 o P=16) deja significativamente menos
violaciones que la búsqueda local P=1 con el mismo presupuesto.

- **Tesis:** "un algoritmo memético que conserva el inventario". La búsqueda local basta en el régimen fácil,
  pero la población evita estancamientos cuando el traslape es alto.
- **Título tentativo:** *An Inventory-Preserving Memetic Algorithm for Constructing Highly Correlated
  Two-Dimensional Dual Site-Bond Porous Networks*
- **Resultados clave:** curva de "régimen". Ventaja del GA frente a P=1 en función de Ω y L, mostrando dónde
  empieza a convenir la población.

## 4. Escenario B — la búsqueda local gana en todo el rango probado

**Condición:** la búsqueda local P=1 es igual o mejor en todos los casos.

- **Tesis:** "la representación y el movimiento importan más que la población". Se propone una búsqueda local
  que conserva el inventario y está guiada por el PC. Se muestra, con un estudio controlado, que agregar
  población y cruce genético **no** mejora (y explicamos por qué: tasa de aceptación de hijos < 2%). Es un
  resultado metodológico valioso para la comunidad de metaheurísticas aplicadas a materiales.
- **Título tentativo:** *Inventory-Preserving Local Search versus Genetic Algorithms for Constructing Correlated
  Dual Site-Bond Porous Networks: A Controlled Comparison*
- **Paralelismo:** en lugar de OpenMP dentro del GA (eficiencia de 19%), réplicas independientes en paralelo
  (escalado casi perfecto) y, si se requiere, paralelización por dominios de una sola red grande.
- **Riesgo:** el artículo deja de ser "un GA". Es honesto, pero hay que redactar la motivación con cuidado: se
  partió de un GA, se diagnosticaron sus defectos y se llegó a un método mejor.

## 5. Estructura propuesta (común, con variantes A/B)

1. **Introducción:** la construcción DSBM como problema combinatorio con restricción de inventario. Brecha:
   los métodos previos alteran radios (relajación) o no conservan la distribución.
2. **DSBM y Principio de Construcción.** Definición de Ω (fórmula para uniformes del mismo ancho; generalizar si
   se usan otras distribuciones).
3. **Método**
   - 3.1 Representación como permutación del inventario.
   - 3.2 Función objetivo por incidencia y evaluación local.
   - 3.3 Movimientos (sitio–sitio, enlace–enlace), aceptación guiada.
   - 3.4 Marco poblacional: cruce PMX/rango y reparación (A: como componente; B: como variante evaluada).
   - 3.5 Generador aleatorio y reproducibilidad.
4. **Diseño experimental:** presupuesto en evaluaciones, semillas, métricas crudas, pruebas estadísticas,
   hardware.
5. **Resultados**
   - 5.1 Validación: inventario exacto, emergencia de correlación vs Ω (figura cruda de 3 paneles + C(r)
     promediado sobre semillas).
   - 5.2 Comparación con igual presupuesto: evidencias 01 y 02, con curvas de convergencia.
   - 5.3 Por qué el cruce aporta poco o mucho: tasa de aceptación de hijos (piloto/evidencia 03).
   - 5.4 Escalamiento en L y tiempos.
   - 5.5 (Recomendado) Equilibrado posterior a la factibilidad y comparación con la teoría del DSBM (punto 2 del
     plan de mejoras).
6. **Discusión:** limitaciones (2D, distribuciones uniformes, Ω efectivo), extensión a 3D.
7. **Conclusiones.**

## 6. Qué se elimina o reemplaza del manuscrito actual

| Elemento actual | Acción |
|---|---|
| Fig. 8 (convergencia ~2 650 errores) | Reemplazar por las curvas de convergencia con igual presupuesto |
| Fig. 9–10 y speedup 6.22 (OpenMP de v1) | Rehacer con v2 o sustituir por escalamiento de réplicas independientes |
| Fig. 11–13 (suavizadas) | Reemplazar por distribuciones crudas (`analysis/plot_red_cruda.py`) |
| Fig. 14 C(r) | Rehacer sobre datos crudos, promediado sobre semillas |
| Fig. 15 y Tabla 3 (ablación con deriva de enlaces) | Reemplazar por la ablación de la evidencia 01/02 |
| Tabla 1 (tamaño de población) | Rehacer con igual presupuesto (es parte de la ablación) |
| Tabla 2 (MC) | La búsqueda local P=1 es el MC con los mismos movimientos. Opcionalmente mantener el MC original como referencia histórica |
| Relajación de radios (Sec. 3.4) | Eliminar del método. v2 no modifica radios |
| "σ = 50" | Corregir: es el semiancho de una uniforme (σ = 28.9) |
| Citas "[?]" sin resolver | Resolver todas antes de enviar |

## 7. Pendientes antes de redactar

- [x] Experimento 02 completo (L=200, Ω=0.9/0.95) → escenario B (evidencia 02).
- [x] Evidencia 03: piloto de cruces (tasa de aceptación) con n completa.
- [x] Piloto L=200 de mutación focalizada y v3 (evidencia 04) → método final: búsqueda local que conserva el
      inventario, guiada por el PC y **focalizada** (f ≈ 0.5–0.75). v3 queda como variante evaluada y como vía de
      paralelización.
- [ ] Experimento formal con el método final (n ≥ 10, Ω = 0.3/0.6/0.9/0.95, L = 50/200).
- [ ] Equilibrado posterior y criterio de paro común para reportar correlaciones.
- [ ] C(r) promediado sobre ≥ 10 semillas para Ω = 0.3/0.6/0.9 (datos crudos).
- [ ] Escalamiento: L = 50, 100, 200, 500 con el método final.
- [ ] Subir v2, scripts y evidencias al repositorio (rama nueva) y actualizar el README.
- [ ] Respaldar `results/` (no se versiona).
