# Evidencias para la revisión del artículo

Cada archivo documenta un experimento o una verificación con lo necesario para citarlo en el artículo:
- versión del código (hash),
- parámetros, hardware y comandos para reproducirlo,
- resultados numéricos, pruebas estadísticas y figuras,
- interpretación y limitaciones.

| # | Archivo | Tema | Conclusión principal |
|---|---|---|---|
| 00 | [EVIDENCIA_00_verificacion_codigo_original.md](EVIDENCIA_00_verificacion_codigo_original.md) | Verificación del GA original | Los operadores no crean correlación (A1), los dominios de las figuras venían del suavizado (A2), el cruce encoge los enlaces (A4), la mutación es determinista en Linux (B1) |
| 01 | [EVIDENCIA_01_aporte_del_cruce.md](EVIDENCIA_01_aporte_del_cruce.md) | GA v2 frente a variantes sin cruce, con igual presupuesto | El cruce no aporta (p = 0.059). La búsqueda local P=1 es la más eficiente: llega a 0 errores en Ω = 0.9 en 10/10 semillas, mientras que el GA P=64 deja ~2.5% |
| 02 | [EVIDENCIA_02_regimen_grande.md](EVIDENCIA_02_regimen_grande.md) | Régimen grande: L = 200, Ω = 0.9 / 0.95 | La búsqueda local gana en las 6 comparaciones (U = 0, p ≤ 0.012). Orden monótono P=1 < 4 < 16 < 64. Ω = 0.9: 0.03% contra 1–5%. Ω = 0.95: 2.7% contra 4–8.5%. Decide el escenario B |
| 03 | [EVIDENCIA_03_piloto_cruces.md](EVIDENCIA_03_piloto_cruces.md) | Piloto de variantes de cruce (PMX/rango, reparación, P) | Ningún rediseño del cruce alcanza a la búsqueda local. Rango ≈ PMX, la reparación empeora y menos población es mejor. Solo 0.4–1.6% de los hijos sobreviven |
| 04 | [EVIDENCIA_04_recocido_y_teselas.md](EVIDENCIA_04_recocido_y_teselas.md) | Recocido simulado, poblaciones por tesela (v3) y mutación focalizada | **La mutación focalizada reduce el residuo 39% (L=50) y 28% (L=200, f=0.75).** El recocido no ayuda. v3 queda 6% (K=1) y 22% (K=2) peor que la búsqueda local. La búsqueda local en Ω=0.95 no está atascada, converge lentamente |
| 05 | [EVIDENCIA_05_ga_con_mutacion_focal.md](EVIDENCIA_05_ga_con_mutacion_focal.md) | El GA con la misma mutación focalizada | La focal mejora al GA (−15% con P=4, −10% con P=16), menos que a P=1 (−29%). El mejor GA focal (5 232) queda peor incluso que la búsqueda local sin focal (4 328). El orden P=1 < 4 < 16 se mantiene |
| 06 | [EVIDENCIA_06_omega_metodo_final.md](EVIDENCIA_06_omega_metodo_final.md) | Dependencia de Ω con el método final (L=100) | La correlación crece con Ω (34% de misma categoría con Ω=0.3, 84% con Ω=0.95). La mutación focal llega a 0 en Ω=0.9 en 4/5 semillas (0/5 sin ella). Con Ω bajo llega a 0 de 8 a 14 veces más rápido, pero con menos correlación: la correlación depende de la trayectoria |

Plan de reformulación del artículo: [../PLAN_REFORMULACION_ARTICULO.md](../PLAN_REFORMULACION_ARTICULO.md)

**Nota sobre los sha256 del código.** Los `sha256 fuente` de cada `entorno.txt` se calcularon sobre la copia de
trabajo en Windows. La de v2 tenía fin de línea CRLF, y git guarda los archivos con LF, así que en un clon el sha de
v2 no coincide directamente. Para comprobarlo:
`git show <commit>:src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c | sed 's/$/\r/' | sha256sum`.
El sha de la versión final de v2 (`f36f445e…`) se reproduce así con el commit `2a9bc5a`. El de v3 (`8df40210…`)
coincide directamente. Las evidencias 01–04 se corrieron con versiones intermedias de v2 que no se versionaron por
separado. Cada cambio posterior se verificó con una prueba de regresión: con las opciones de esas evidencias, la
versión final da resultados idénticos bit a bit. Para reproducirlas, basta la versión final.

Los datos de cada evidencia están en `results/`. Se versionan los resúmenes, las tablas, las pruebas estadísticas,
las figuras, los logs por corrida y las trazas de convergencia. Las redes finales (`results/*/redes/`, ~175 MB) no
se versionan; hay que respaldarlas aparte (por ejemplo, en Zenodo junto con el artículo).
