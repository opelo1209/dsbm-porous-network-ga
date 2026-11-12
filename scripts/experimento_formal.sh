#!/usr/bin/env bash
# EXPERIMENTO FORMAL (evidencia 07): confirmacion con n = 10 semillas del metodo final.
#
#   L = 50  (presupuesto 4e9)  y  L = 200 (presupuesto 2e10)
#   Omega = 0.3, 0.6, 0.9, 0.95 ; semillas 1..10
#   Variantes: (1+1)-EA uniforme, (1+1)-EA dirigido por conflictos con f = 0.5/0.75/0.9,
#              GA mu=4 PMX con mutacion dirigida (f = 0.75)
#
# Reutiliza corridas previas con la MISMA semilla, presupuesto y opciones (el codigo es
# identico bit a bit con esas opciones; ver docs/evidencias/README.md):
#   L=50 : experimento_cruce (busquedaLocal_P1, Omega 0.6 y 0.9)
#   L=200: experimento_regimen_grande (busquedaLocal_P1, Omega 0.9 y 0.95),
#          piloto_v3_L200 (focal0.75_P1, Omega 0.95), piloto_ga_focal (GA4_pmx_focal, Omega 0.95)
# Es reanudable: las corridas con linea RESUMEN no se repiten.
#
# Uso:  bash scripts/experimento_formal.sh        (unas 9-10 h con 8 corridas simultaneas)
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RAIZ"
R=results
OMEGAS="0.3:440:370:50 0.6:440:400:50 0.9:410:400:50 0.95:405:400:50"
VARIANTES="busquedaLocal_P1:1:--sin-cruce focal0.5_P1:1:--sin-cruce,--focal=0.5 focal0.75_P1:1:--sin-cruce,--focal=0.75 focal0.9_P1:1:--sin-cruce,--focal=0.9 GA4_pmx_focal:4:--focal=0.75"
SEMILLAS="1 2 3 4 5 6 7 8 9 10"

reutilizar () {   # origen variante omega semillas destino
    local o="$1" v="$2" om="$3" ss="$4" d="$5"
    mkdir -p "$d/logs" "$d/trazas"
    for s in $ss; do
        local id="w${om}_${v}_s${s}"
        if [ -f "$o/logs/$id.txt" ] && grep -q -E '^RESUMEN,[^,]+,[0-9]' "$o/logs/$id.txt" && [ ! -f "$d/logs/$id.txt" ]; then
            cp "$o/logs/$id.txt" "$d/logs/"; cp "$o/trazas/$id.csv" "$d/trazas/" 2>/dev/null
            echo "$id <- $o" >> "$d/reutilizadas.txt"
        fi
    done
}

D50=$R/formal_L50
reutilizar $R/experimento_cruce busquedaLocal_P1 0.6 "$SEMILLAS" $D50
reutilizar $R/experimento_cruce busquedaLocal_P1 0.9 "$SEMILLAS" $D50
L=50 PRESUPUESTO=4e9 SEMILLAS="$SEMILLAS" OMEGAS="$OMEGAS" VARIANTES="$VARIANTES" \
  OUT=$D50 TRABAJOS=8 bash scripts/experimento_cruce.sh > $R/formal_L50_consola.txt 2>&1
REFERENCIA=busquedaLocal_P1 python3 analysis/resumen_cruce.py $D50 >> $R/formal_L50_consola.txt 2>&1

D200=$R/formal_L200
reutilizar $R/experimento_regimen_grande busquedaLocal_P1 0.9 "1 2 3 4 5" $D200
reutilizar $R/experimento_regimen_grande busquedaLocal_P1 0.95 "1 2 3 4 5" $D200
reutilizar $R/piloto_v3_L200 focal0.75_P1 0.95 "1 2 3" $D200
reutilizar $R/piloto_ga_focal GA4_pmx_focal 0.95 "1 2 3" $D200
L=200 PRESUPUESTO=2e10 SEMILLAS="$SEMILLAS" OMEGAS="$OMEGAS" VARIANTES="$VARIANTES" \
  OUT=$D200 TRABAJOS=8 bash scripts/experimento_cruce.sh > $R/formal_L200_consola.txt 2>&1
REFERENCIA=busquedaLocal_P1 python3 analysis/resumen_cruce.py $D200 >> $R/formal_L200_consola.txt 2>&1
echo "Experimento formal terminado: $(date)" >> $R/formal_L200_consola.txt
