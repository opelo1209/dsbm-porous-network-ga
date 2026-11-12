# EVIDENCIA 05: el GA con mutación focalizada

**Fecha:** 2026-11-12. **Estado:** piloto completo (n = 3, exploratorio).

## Pregunta

La evidencia 04 mostró que la mutación focalizada mejora a la búsqueda local P=1 en 28% (L=200, Ω=0.95). Pero en la
evidencia 02 el GA se había comparado con la mutación *sin* focalizar. Para que la comparación sea justa: **si el GA
recibe la misma mutación focalizada, ¿alcanza o supera a la búsqueda local focalizada?** (Pregunta del autor,
2026-11-12.)

## Diseño

| Parámetro | Valor |
|---|---|
| Código | `src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c`, sha256 `f36f445e…90fc63fb`. Con las opciones por defecto es idéntico bit a bit a las versiones de las evidencias 01–03 |
| Red | L = 200, Ω = 0.95 (μS = 405, μB = 400, semiancho 50) |
| Presupuesto | 2 × 10¹⁰ evaluaciones de sitio |
| Semillas | 1, 2, 3 (mismas M_base y configuración inicial que las evidencias 02 y 04) |
| Mutación focal | f = 0.75, el mejor valor en L=200 (evidencia 04) |
| Hilos | 1 por corrida, 8 simultáneas; i7-7700HQ, WSL2, gcc 9.4.0 `-O3` |

Variantes nuevas: GA P=4 con PMX + focal; GA P=4 con rangos + reparación (2 barridos) + focal; GA P=16 con PMX +
focal. Las referencias se copiaron de `results/piloto_v3_L200`, con las mismas semillas: búsqueda local P=1 con y sin
focal. Las variantes del GA sin focal provienen de la evidencia 02.

Para reproducirlo:

```bash
L=200 PRESUPUESTO=2e10 SEMILLAS="1 2 3" OMEGAS="0.95:405:400:50" \
VARIANTES="busquedaLocal_P1:1:--sin-cruce focal0.75_P1:1:--sin-cruce,--focal=0.75 \
GA4_pmx_focal:4:--focal=0.75 GA4_rango_rep2_focal:4:--cruce=rango,--reparar=2,--focal=0.75 \
GA16_pmx_focal:16:--focal=0.75" OUT=results/piloto_ga_focal TRABAJOS=8 bash scripts/experimento_cruce.sh
REFERENCIA=focal0.75_P1 python3 analysis/resumen_cruce.py results/piloto_ga_focal
```

## Resultados

Errores al agotar el presupuesto (de 160 000 incidencias), semillas 1 / 2 / 3:

| Variante | Sin focal | Con focal (f = 0.75) | Mejora por la focal | Error % con focal | C(1) con focal |
|---|---|---|---|---|---|
| **Búsqueda local P=1** | 4 367 / 4 328 / 4 104 | **3 131 / 3 103 / 2 845** | **−28 a −31%** | **1.89** | 0.864 |
| GA P=4, PMX | — | 5 345 / 5 232 / 5 003 | — | 3.25 | 0.773 |
| GA P=4, rangos + reparación | 6 490 / 6 453 / 6 210 | 5 555 / 5 486 / 5 179 | −14 a −17% | 3.38 | 0.761 |
| GA P=16, PMX | 9 177 / 9 132 / 8 936 | 8 323 / 8 250 / 7 984 | −9 a −11% | 5.12 | 0.656 |

Tasa de aceptación de hijos con focal: PMX P=16, 0.01%; PMX P=4, 0.11–0.13%; rangos + reparación P=4, 0.14–0.17%.

En cada semilla, la búsqueda local focal es mejor que todas las variantes del GA focal (U = 0 en las tres
comparaciones; con n = 3 el p mínimo es 0.1).

## Interpretación

1. **La focalización también mejora al GA, pero menos:** 14–17% con P=4 y 9–11% con P=16, contra 28–31% con P=1. La
   focalización rinde más cerca de la factibilidad, cuando quedan pocas violaciones. Con P individuos, cada uno recibe
   1/P del presupuesto y pasa menos tiempo en esa fase final, donde la focal es más útil.
2. **El GA focal no alcanza ni a la búsqueda local *sin* focal.** El mejor GA focal (P=4 PMX, 5 232) queda 21% peor
   que la búsqueda local normal (4 328) y 69% peor que la focal (3 103).
3. **El orden por tamaño de población no cambia:** P=1 < P=4 < P=16, con o sin focal. La causa de fondo no la ataca
   la mutación. Es la dilución del presupuesto entre individuos y un cruce cuyos hijos casi nunca sobreviven (≤ 0.17%).
4. **El cruce PMX y el de rangos con reparación quedan prácticamente empatados con P=4** (5 232 contra 5 486).

## Implicaciones para el artículo

- Cierra la objeción "comparan un GA débil contra una búsqueda local mejorada": con exactamente los mismos operadores
  de mutación, la población sigue perdiendo.
- Encuadre recomendado dentro de la computación evolutiva: el método final es un **(1+1)-EA** (un individuo,
  mutación y reemplazo elitista) con mutación que conserva el inventario, guiada por el PC y focalizada en las
  violaciones. El estudio controlado muestra que aumentar μ (P) y agregar recombinación empeora el desempeño para este
  problema, y explica por qué: dilución del presupuesto y tasa de aceptación de hijos < 0.2% en L=200.
- La tabla "sin focal / con focal" por tamaño de población es una figura natural de la sección de ablación: muestra
  que las dos mejoras (focalizar y reducir la población) se suman.

## Limitaciones

- n = 3. Las diferencias son grandes y consistentes semilla por semilla, pero en el experimento formal conviene
  n ≥ 10.
- Un solo valor de f (0.75), elegido para P=1. El GA podría preferir otro valor, aunque la brecha es demasiado
  grande para que eso la cierre.
- Solo L = 200, Ω = 0.95.
