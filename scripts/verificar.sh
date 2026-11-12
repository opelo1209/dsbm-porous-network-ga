#!/usr/bin/env bash
# Prueba rapida de los hallazgos A1, A2, A4, B1 y B2 (y opcionalmente A3).
#
#   bash scripts/verificar.sh            # ~1 minuto: arnes de operadores + demo de suavizado
#   CORRIDA_GA=1 bash scripts/verificar.sh   # ademas corre el GA real (L=50, Omega=0.9, 13000 gen)
#
# Requisitos en Ubuntu:  sudo apt install build-essential python3-numpy python3-scipy
set -euo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$RAIZ/results/verificacion"
mkdir -p "$OUT"
CC="${CC:-gcc}"
L="${L:-50}"; P="${P:-16}"; G="${G:-300}"; SEMILLA="${SEMILLA:-12345}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

compilar () {   # $1 = fuente del GA relativa a tests/ (portable Linux/MinGW), $2 = binario
    "$CC" -O2 -std=c99 -fopenmp -w -DGA_SRC="\"$1\"" "$RAIZ/tests/verificar_operadores.c" -o "$2" -lm
}

echo "################ 1) GA PROPUESTO (operadores con restriccion de categoria) ################"
compilar "../src/genetico/ConstructorRedes2D_C4_Genetico_Final.c" "$OUT/verificar_ga"
"$OUT/verificar_ga" "$L" "$P" "$G" 410 400 50 "$SEMILLA" | tee "$OUT/verificar_ga.txt"

echo
echo "################ 2) VARIANTE DE ABLACION (sin restriccion de categoria) ################"
compilar "../src/ablacion/ConstructorRedes2D_C4_Ablacion.c" "$OUT/verificar_ablacion"
"$OUT/verificar_ablacion" "$L" "$P" "$G" 410 400 50 "$SEMILLA" | tee "$OUT/verificar_ablacion.txt"

echo
echo "################ 3) ARTEFACTO DEL SUAVIZADO EN LAS FIGURAS (A2) ################"
python3 "$RAIZ/analysis/verificar_red.py" --aleatoria 200 | tee "$OUT/verificar_suavizado.txt"

if [[ "${CORRIDA_GA:-0}" == "1" ]]; then
    echo
    echo "################ 4) CORRIDA REAL DEL GA: L=50, Omega=0.9, P=64, 13000 gen (A3) ################"
    "$CC" -O3 -std=c99 -fopenmp "$RAIZ/src/genetico/ConstructorRedes2D_C4_Genetico_Final.c" -o "$OUT/ga_omp" -lm
    ( cd "$OUT" && printf '50\n410\n400\n50\n' | ./ga_omp x x > corrida_ga.txt )
    grep -E "Presupuesto|Solución|Tiempo total" "$OUT/corrida_ga.txt" || true
    echo "  (el articulo reporta 6 +/- 3 errores residuales para esta configuracion)"
    python3 "$RAIZ/analysis/verificar_red.py" "$OUT/red_colores.csv"
fi

echo
echo "Resultados guardados en $OUT"
