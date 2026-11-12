#!/usr/bin/env bash
# Prueba de la version 2 del GA (correcciones A1 y A4).
#
#   bash scripts/verificar_v2.sh                  # L=50, 3000 gen, unos 2-3 minutos
#   L=100 G=13000 bash scripts/verificar_v2.sh     # configuracion completa (mas lento)
#
# Comprueba:
#   1) Reproducibilidad: misma semilla con 1 y con N hilos -> archivos identicos.
#   2) Conservacion exacta del inventario de sitios y enlaces (A4).
#   3) Organizacion espacial cruda vs Omega = 0.3 / 0.6 / 0.9 (A1).
#
# Requisitos en Ubuntu: sudo apt install build-essential python3-numpy python3-scipy
set -euo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$RAIZ/results/verificacion_v2"
mkdir -p "$OUT"
CC="${CC:-gcc}"
L="${L:-50}"; P="${P:-64}"; G="${G:-3000}"; SEMILLA="${SEMILLA:-2026}"
HILOS="${HILOS:-$(nproc 2>/dev/null || echo 4)}"
PY="$(command -v python3 || command -v python)"

echo "Compilando GA v2..."
"$CC" -O3 -std=gnu99 -Wall -fopenmp "$RAIZ/src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c" -o "$OUT/ga_v2" -lm

echo
echo "################ 1) REPRODUCIBILIDAD (L=30, 300 gen, 1 vs $HILOS hilos) ################"
( cd "$OUT" && OMP_NUM_THREADS=1      ./ga_v2 30 410 400 50 16 300 7 rep1 > /dev/null
              OMP_NUM_THREADS="$HILOS" ./ga_v2 30 410 400 50 16 300 7 repN > /dev/null )
if cmp -s "$OUT/rep1_completa.csv" "$OUT/repN_completa.csv"; then
    echo "  OK: resultados identicos con 1 y $HILOS hilos (misma semilla)."
else
    echo "  FALLA: los resultados dependen del numero de hilos."
fi

echo
echo "################ 2-3) OMEGA = 0.3 / 0.6 / 0.9  (L=$L, P=$P, G=$G, semilla=$SEMILLA) ################"
for cfg in "0.3 440 370" "0.6 440 400" "0.9 410 400"; do
    set -- $cfg
    echo "--- Omega = $1 ---"
    ( cd "$OUT" && OMP_NUM_THREADS="$HILOS" ./ga_v2 "$L" "$2" "$3" 50 "$P" "$G" "$SEMILLA" "omega_$1" \
        | grep -vE "^Optimo|^Estancamiento" | tee "log_omega_$1.txt" \
        | grep -E "residual|Solucion|Inventario|incidencia|misma categoria|Tiempo" )
done

echo
echo "################ Analisis crudo (sin suavizado) ################"
"$PY" "$RAIZ/analysis/verificar_red.py" "$OUT/omega_0.3_colores.csv" "$OUT/omega_0.6_colores.csv" "$OUT/omega_0.9_colores.csv"
echo
echo "################ Imagenes PNG de la distribucion (sin suavizado) ################"
"$PY" "$RAIZ/analysis/plot_red_cruda.py" \
    "$OUT/omega_0.3_colores.csv" "$OUT/omega_0.6_colores.csv" "$OUT/omega_0.9_colores.csv" \
    --titulos "Omega = 0.3" "Omega = 0.6" "Omega = 0.9" \
    --panel "$OUT/comparacion_omega.png"

echo
echo "Esperado si A1 quedo corregido: '% misma categoria' y C(1) CRECEN con Omega (version original: ~33% y ~0 siempre)."
echo "Resultados en $OUT"
