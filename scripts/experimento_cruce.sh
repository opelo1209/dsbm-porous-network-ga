#!/usr/bin/env bash
# Experimentos de comparacion de variantes del GA v2 con IGUAL presupuesto.
#
# Todas las variantes reciben exactamente el mismo presupuesto de costo
# (evaluaciones de errSitio, independiente del hardware) y corren con 1 hilo
# cada una, de modo que el tiempo tambien es comparable. Se lanzan varias
# corridas en paralelo como procesos independientes.
#
# EXPERIMENTO 01 (valores por defecto):
#   bash scripts/experimento_cruce.sh
# EXPERIMENTO 02 (regimen grande, ver docs/evidencias):
#   L=200 PRESUPUESTO=2e10 SEMILLAS="1 2 3 4 5" \
#   OMEGAS="0.9:410:400:50 0.95:405:400:50" \
#   VARIANTES="GA_P64:64 GA_P16:16 busquedaLocal_P1:1:--sin-cruce" \
#   OUT=results/experimento_regimen_grande bash scripts/experimento_cruce.sh
#
# Formatos:
#   OMEGAS    = lista de omega:mediaS:mediaE:semiancho
#   VARIANTES = lista de nombre:P[:opciones extra separadas por comas]
#   BINARIO   = ejecutable alternativo (por defecto se compila GA v2)
#
# Salida en $OUT: resumen.csv, trazas/, logs/, redes/, entorno.txt
# Analisis:       python3 analysis/resumen_cruce.py $OUT
set -euo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${OUT:-$RAIZ/results/experimento_cruce}"
case "$OUT" in /*) ;; *) OUT="$RAIZ/$OUT" ;; esac
SRC="${SRC:-$RAIZ/src/genetico_v2/ConstructorRedes2D_C4_Genetico_v2.c}"
CC="${CC:-gcc}"
L="${L:-50}"
PRESUPUESTO="${PRESUPUESTO:-4e9}"
SEMILLAS="${SEMILLAS:-1 2 3 4 5 6 7 8 9 10}"
OMEGAS="${OMEGAS:-0.6:440:400:50 0.9:410:400:50}"
VARIANTES="${VARIANTES:-GA_P64:64 GA_P16:16 sinCruce_P64:64:--sin-cruce busquedaLocal_P1:1:--sin-cruce}"
TRABAJOS="${TRABAJOS:-$(( $(nproc 2>/dev/null || echo 4) / 2 ))}"
[ "$TRABAJOS" -lt 1 ] && TRABAJOS=1

mkdir -p "$OUT/trazas" "$OUT/logs" "$OUT/redes"
BIN="$OUT/ga_v2"
# v3 (poblaciones por tesela): se usa en las variantes que incluyen la opcion --motor=v3.
SRC3="${SRC3:-$RAIZ/src/genetico_v3/ConstructorRedes2D_C4_Teselas_v3.c}"
BIN3="$OUT/ga_v3"
# Si se reanuda un experimento, se guarda el entorno anterior antes de sobrescribirlo.
[ -f "$OUT/entorno.txt" ] && cat "$OUT/entorno.txt" >> "$OUT/entorno_anterior.txt"
"$CC" -O3 -std=gnu99 -Wall -fopenmp ${CFLAGS_EXTRA:-} "$SRC" -o "$BIN" -lm 2>&1 || exit 1
if [ -f "$SRC3" ]; then "$CC" -O3 -std=gnu99 -Wall -fopenmp ${CFLAGS_EXTRA:-} "$SRC3" -o "$BIN3" -lm 2>&1 || exit 1; fi

{
    echo "fecha: $(date -Iseconds)"
    echo "compilador: $("$CC" --version | head -n 1)"
    echo "flags: -O3 -std=gnu99 -fopenmp ${CFLAGS_EXTRA:-} ; OMP_NUM_THREADS=1 por corrida ; $TRABAJOS corridas simultaneas"
    echo "sistema: $(uname -a)"
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 | sed 's/^ //' || echo desconocido)"
    echo "fuente: $SRC"
    echo "sha256 fuente: $(sha256sum "$SRC" | cut -d' ' -f1)"
    [ -f "$SRC3" ] && echo "fuente v3: $SRC3" && echo "sha256 fuente v3: $(sha256sum "$SRC3" | cut -d' ' -f1)"
    echo "git commit: $(git -C "$RAIZ" rev-parse HEAD 2>/dev/null || echo NA) (cambios sin commit: $(git -C "$RAIZ" status --porcelain 2>/dev/null | wc -l))"
    echo "L=$L presupuesto=$PRESUPUESTO semillas=[$SEMILLAS]"
    echo "omegas (omega:mediaS:mediaE:semiancho): $OMEGAS"
    echo "variantes (nombre:P:opciones): $VARIANTES"
} > "$OUT/entorno.txt"
cat "$OUT/entorno.txt"

# Una linea por corrida: id P mS mE a semilla variante opciones
TAREAS="$OUT/tareas.txt"
: > "$TAREAS"
for o in $OMEGAS; do
    IFS=: read -r om mS mE a <<< "$o"
    for v in $VARIANTES; do
        IFS=: read -r nombre P extra <<< "$v"
        for s in $SEMILLAS; do
            # Sin espacio final: xargs -L trata un blanco al final de linea como
            # continuacion de la siguiente y fusionaria tareas.
            echo "w${om}_${nombre}_s${s} $P $mS $mE $a $s ${nombre}${extra:+ $extra}" >> "$TAREAS"
        done
    done
done
echo "Corridas: $(wc -l < "$TAREAS")  (simultaneas: $TRABAJOS)"

correr () {
    local id="$1" P="$2" mS="$3" mE="$4" a="$5" s="$6" nombre="$7" extra="${8:-}"
    # Reanudable: si la corrida ya termino (tiene linea RESUMEN), no se repite.
    if [ -f "$OUT/logs/$id.txt" ] && grep -q -E '^RESUMEN,[^,]+,[0-9]' "$OUT/logs/$id.txt"; then
        echo "  ya hecho $id"; return 0
    fi
    local bin="$BIN"
    case ",$extra," in *,--motor=v3,*) bin="$BIN3" ;; esac
    extra="$(echo "$extra" | tr ',' '\n' | grep -v -x -e '--motor=v3' | tr '\n' ' ')"
    # shellcheck disable=SC2046,SC2086
    OMP_NUM_THREADS=1 "$bin" "$L" "$mS" "$mE" "$a" "$P" 0 "$s" "$OUT/redes/$id" \
        $extra \
        --presupuesto="$PRESUPUESTO" --traza="$OUT/trazas/$id.csv" --variante="$nombre" \
        > "$OUT/logs/$id.txt" 2>&1
    echo "  listo $id: error%,misma_cat,C1 = $(grep -E '^RESUMEN,[^,]+,[0-9]' "$OUT/logs/$id.txt" | cut -d, -f12,13,14)"
}
export -f correr
export BIN BIN3 L PRESUPUESTO OUT

xargs -P "$TRABAJOS" -L 1 bash -c 'correr "$@"' _ < "$TAREAS"

{
    echo "variante,L,omega,semilla,P,cruce,generaciones,evaluaciones,tiempo_s,error,error_pct,misma_cat,C1,inventario_ok,fitness_ok"
    grep -h -E '^RESUMEN,[^,]+,[0-9]' "$OUT"/logs/*.txt | sed 's/^RESUMEN,//'
} > "$OUT/resumen.csv"
echo
echo "Resumen: $OUT/resumen.csv ($(($(wc -l < "$OUT/resumen.csv") - 1)) corridas)"
echo "Analisis: python3 analysis/resumen_cruce.py $OUT"
