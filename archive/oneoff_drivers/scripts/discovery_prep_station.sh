#!/bin/bash
# Discovery prep for a wave station: fetch coords -> PNSN-driven candidate detection -> picker score
# (baseline) -> ADAPTIVE-P threshold -> data/<sta>_cand_filtered.parquet. CPU-only.
# Usage: scripts/discovery_prep_station.sh <STA>
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]'); sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) $*" >> "$STATUS"; echo "$(date +%H:%M) $*"; }
CAND=data/${sta}_pnsn_candidates_100km.parquet

if [ ! -f data/${sta}_cand_filtered.parquet ]; then
  # station coords (IRIS) -> bbox +-0.9 lat / +-1.35 lon
  coords=$($PY -c "
from obspy.clients.fdsn import Client
inv=Client('IRIS').get_stations(network='PB', station='$STA', level='station')
st=inv[0][0]; print(f'{st.latitude:.4f} {st.longitude:.4f}')" 2>/dev/null)
  [ -n "$coords" ] || { log "$STA prep ABORT: no coords"; exit 1; }
  set -- $coords; LA=$1; LO=$2
  la0=$($PY -c "print(round($LA-0.9,2))"); la1=$($PY -c "print(round($LA+0.9,2))")
  lo0=$($PY -c "print(round($LO-1.35,2))"); lo1=$($PY -c "print(round($LO+1.35,2))")

  if [ ! -f "$CAND" ]; then
    log "$STA candidate detect start (bbox $la0 $la1 $lo0 $lo1)"
    PYTHONPATH=src $PY scripts/discover_nllb_pnsn_driven.py --station "$STA" --wfdir data/waveforms \
       --bbox $la0 $la1 $lo0 $lo1 --pnsn catalogs/pnsn_tremor_cascadia_full.csv --fs 100 \
       --candidates-only --candidates-out "$CAND" --workers 12 > logs/cand_${sta}.log 2>&1 \
       || { log "$STA prep ABORT: candidate detect failed"; exit 2; }
  fi

  cp -f lfe_features/models/tremor_picker_b011.joblib lfe_features/models/tremor_picker_${sta}.joblib
  if [ ! -f data/${sta}_cand_baseline.parquet ]; then
    log "$STA picker score start"
    PYTHONPATH=src $PY lfe_features/score_candidates.py --net PB --sta "$STA" --cand "$CAND" \
       --y0 2010 --y1 2026 --thr 0.5 --workers 12 > logs/score_${sta}.log 2>&1 \
       || { log "$STA prep ABORT: picker score failed"; exit 3; }
  fi

  $PY scripts/adaptive_cand_threshold.py "$STA" 15000 > logs/adaptive_${sta}.log 2>&1 \
     || { log "$STA prep ABORT: adaptive threshold failed"; exit 4; }
  log "$STA discovery prep done: $(tail -1 logs/adaptive_${sta}.log)"
fi
