# EVIDENCIA 06: dependencia de Ω con el método final (L = 100)

**Fecha:** 2026-11-12. **Estado:** completo (5 semillas × 2 variantes × 4 Ω = 40 corridas). Es la base de la Fig. 9 y
la Tabla 7 del manuscrito revisado (`Paper_TCJ_EA_v2`).

## Pregunta

¿Cómo cambian la organización espacial y el residuo con Ω usando el método final, el (1+1)-EA con mutación dirigida
por conflictos? ¿Y cómo se compara con el (1+1)-EA con mutación uniforme?

## Diseño

| Parámetro | Valor |
|---|---|
| Código | v2, mismo fuente que la evidencia 05 (ver `results/omega_final/entorno.txt`) |
| Red | L = 100 |
| Configuraciones | Ω = 0.3 (440/370), 0.6 (440/400), 0.9 (410/400), 0.95 (405/400); semiancho 50 |
| Variantes | `focal0.75_P1` (`--sin-cruce --focal=0.75`) y `busquedaLocal_P1` (`--sin-cruce`) |
| Presupuesto | 2 × 10⁹ evaluaciones; la corrida se detiene antes si llega a f = 0 |
| Semillas | 1–5 |
| Hardware | i7-7700HQ, WSL2, gcc 9.4.0 `-O3`; 1 hilo por corrida, 8 simultáneas |

Para reproducirlo:

```bash
L=100 PRESUPUESTO=2e9 SEMILLAS="1 2 3 4 5" \
OMEGAS="0.3:440:370:50 0.6:440:400:50 0.9:410:400:50 0.95:405:400:50" \
VARIANTES="focal0.75_P1:1:--sin-cruce,--focal=0.75 busquedaLocal_P1:1:--sin-cruce" \
OUT=results/omega_final TRABAJOS=8 bash scripts/experimento_cruce.sh
REFERENCIA=focal0.75_P1 python3 analysis/resumen_cruce.py results/omega_final
python3 analysis/figuras_articulo.py <carpeta_de_figuras>     # fig_omega.pdf y tabla_omega.csv
```

## Resultados

U = mutación uniforme; CD = dirigida por conflictos (f = 0.75). Evaluaciones: mediana hasta f = 0 o hasta agotar el
presupuesto.

| Ω | Mutación | Mediana f (% de 4N) | f = 0 | Evaluaciones | % misma categoría | C(1) |
|---|---|---|---|---|---|---|
| 0.3 | U | 0 (0) | 5/5 | 1.1e6 | 34.6 | 0.031 |
| 0.3 | CD | 0 (0) | 5/5 | 8.0e4 | 34.1 | 0.017 |
| 0.6 | U | 0 (0) | 5/5 | 6.3e6 | 44.5 | 0.298 |
| 0.6 | CD | 0 (0) | 5/5 | 8.0e5 | 41.4 | 0.218 |
| 0.9 | U | 112 (0.28) | 0/5 | 2.0e9 | 83.9 | 0.904 |
| 0.9 | CD | 0 (0) | 4/5 | 1.2e9 | 81.3 | 0.897 |
| 0.95 | U | 1 256 (3.14) | 0/5 | 2.0e9 | 83.8 | 0.788 |
| 0.95 | CD | 927 (2.32) | 0/5 | 2.0e9 | 84.4 | 0.812 |

C(r) promediado sobre las 5 semillas (CD): C(5) ≈ 0.36 (Ω = 0.9) y 0.38 (Ω = 0.95). A partir de r ≈ 15 es
prácticamente 0. Con Ω ≤ 0.6, C(r) ≈ 0 desde r = 2.

Pruebas de Mann-Whitney, CD contra U:
- error final: p = 0.085 (Ω = 0.9) y p = 0.095 (Ω = 0.95);
- evaluaciones hasta f = 0: p = 0.008 (Ω = 0.3 y 0.6);
- % de misma categoría: p = 0.02 (Ω = 0.3) y p = 0.008 (Ω = 0.6).

Las 40 corridas pasaron la verificación de inventario y la de fitness.

## Interpretación

1. **La correlación crece con Ω:** con Ω = 0.3 la red es indistinguible de una aleatoria. Con Ω = 0.6 hay
   correlación débil de corto alcance. Con Ω ≥ 0.9 se forman dominios extensos.
2. **La mutación dirigida por conflictos también mejora en L = 100**, pero con n = 5 no es significativa:
   - Ω = 0.9: 4/5 corridas llegan a 0, contra 0/5 con mutación uniforme;
   - Ω = 0.95: −26% de residuo;
   - la única corrida CD que no llegó a 0 en Ω = 0.9 terminó con 185 violaciones.
3. **Hallazgo nuevo: la correlación depende de la trayectoria.** Con Ω bajo, las dos versiones llegan a f = 0. La
   CD lo hace de 8 a 14 veces más rápido y, por eso, se detiene con *menos* organización espacial (Ω = 0.6: 41.4%
   contra 44.5%; p = 0.008). La correlación medida en la primera configuración factible no es una propiedad solo de
   Ω, sino también del algoritmo y del momento en que se detiene.

## Implicaciones para el artículo

- Es la figura y la tabla de dependencia con Ω del manuscrito (sección 6.6).
- El hallazgo 3 entra en la Discusión. Para aplicaciones físicas hace falta una etapa común posterior a la
  factibilidad: intercambios restringidos a Δf = 0 que muestreen las configuraciones factibles. Solo después de
  ella son comparables las correlaciones entre métodos. Queda como trabajo futuro, y coincide con el punto
  "equilibrado posterior" del plan.

## Limitaciones

- n = 5 y un solo valor de f.
- L = 100 con presupuesto 2 × 10⁹, menor por sitio que en L = 200.
- Las correlaciones de Ω ≤ 0.6 corresponden a la primera configuración factible, que no está equilibrada (ver
  hallazgo 3).
