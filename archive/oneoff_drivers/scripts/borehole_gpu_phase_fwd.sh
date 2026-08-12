#!/bin/bash
# FORWARD-ONLY GPU phase (retired-reverse method): cluster -> appB -> select -> resample -> FORWARD densify
# -> verify. No reversed densify (causality certifies reliability at the end; a sampled reverse fake-rate
# control comes later, separately). Halves per-station GPU vs the concurrent fwd+rev variant.
# Usage: scripts/borehole_gpu_phase_fwd.sh <STA>
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]'); sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; mkdir -p logs
log(){ echo "$(date +%H:%M) $*" >> "$STATUS"; echo "$(date +%H:%M) $*"; }
gpu_env(){
  export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
  export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
}
anon_gb(){ awk '/^anon /{print int($2/1073741824)}' /sys/fs/cgroup/memory.stat; }
wait_gpu_free(){
  while true; do
    busy=$(ps -eo pid,cmd | grep -E 'densify_gnw_gpu|discover_gpu' | grep -v grep)
    mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
    [ -z "$busy" ] && [ "${mem:-99999}" -lt 2000 ] && break
    sleep 120
  done
}
wait_mem_ok(){ while [ "$(anon_gb)" -gt 130 ]; do echo "$(date +%H:%M) $STA hold: anon $(anon_gb)GB>130" >> "$STATUS"; sleep 120; done; }

DISC=data/${sta}_disc_p70_2010_2026_m3; CAND=data/${sta}_cand_filtered.parquet
PICK=data/family_picker_${sta}_p70_2010_2026_m3.csv; SEL=${DISC}_sel300.summary.csv
NPZ40=${DISC}_40hz.npz; FLOG=logs/densify_${sta}p90f40.log

[ -f "$CAND" ] || { log "$STA ABORT: missing $CAND"; exit 2; }

if [ ! -f "${DISC}.summary.csv" ]; then
  wait_gpu_free; gpu_env; log "$STA cluster start"
  $PY scripts/discover_gpu.py --station "$STA" --candidates "$CAND" --out "${DISC}.npz" \
     --fs 100 --cc-threshold 0.80 --min-family-members 3 --min-years 1 --workers 12 \
     > logs/cluster_${sta}.log 2>&1 || { log "$STA ABORT: cluster failed"; exit 3; }
  nf=$(wc -l < "${DISC}.summary.csv"); log "$STA cluster done: $((nf-1)) families"
fi
if [ ! -f "$PICK" ]; then
  log "$STA appB start"
  PYTHONPATH=src $PY lfe_features/score_family_stacks_picker.py --net PB --sta "$STA" \
     --members "${DISC}.members.parquet" --summary "${DISC}.summary.csv" --out "$PICK" --workers 12 \
     > logs/appb_${sta}.log 2>&1 || { log "$STA ABORT: appB failed"; exit 4; }
  nl=$($PY -c "import pandas as pd;print((pd.read_csv('$PICK').pred=='LFE').sum())" 2>/dev/null); log "$STA appB done: ${nl:-?} LFE-label families"
fi
# early-exit (Merlin 2026-07-13): <60 LFE-label families cannot reach ~20 causality-certified (~1/3 survival) -> skip densify
nlfe=$($PY -c "import pandas as pd;print(int((pd.read_csv('$PICK').pred=='LFE').sum()))" 2>/dev/null)
if [ "${nlfe:-0}" -lt 60 ]; then
  log "$STA EARLY-EXIT: ${nlfe:-0} LFE-label families (<60) — too weak to reach ~20 certified; skipping densify"
  exit 8
fi
if [ ! -f "$SEL" ]; then
  $PY scripts/select_families_coverage.py "$DISC" "$PICK" 300 > logs/select_${sta}.log 2>&1 \
     || { log "$STA ABORT: select failed"; exit 5; }
  ns=$(wc -l < "$SEL"); log "$STA select done: $((ns-1)) families"
fi
if [ ! -f "$NPZ40" ]; then
  $PY scripts/make_40hz_npz.py "$DISC" > logs/resample_${sta}.log 2>&1 \
     || { log "$STA ABORT: resample failed"; exit 6; }
  log "$STA resample 40hz done"
fi

if ! grep -q 'ALL YEARS DONE' "$FLOG" 2>/dev/null; then
  wait_mem_ok; wait_gpu_free; gpu_env
  log "$STA FORWARD densify start (forward-only; anon $(anon_gb)GB)"
  $PY scripts/densify_gnw_gpu.py --templates-npz "$NPZ40" --summary-csv "$SEL" --min-snr 0 \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 \
     --top-n 0 --despike-mad 8 --workers 8 --start-year 2010 --end-year 2026 \
     --out-dir data --out-prefix "mf_${sta}p90f40_" > "$FLOG" 2>&1
  grep -q 'ALL YEARS DONE' "$FLOG" || { log "$STA ABORT: forward densify incomplete"; exit 7; }
  tc=$(grep -oE '[0-9]+ templates,' "$FLOG" | head -1 | grep -oE '[0-9]+'); want=$(($(wc -l < "$SEL")-1))
  mf=$(ls data/mf_${sta}p90f40_20*.csv 2>/dev/null | head -1)
  hrs=$($PY -c "import pandas as pd; d=pd.read_csv('$mf',usecols=['time'],nrows=400000); h=pd.to_datetime(d.time).dt.hour; print(int(h.min()),int(h.max()))" 2>/dev/null)
  if [ "${tc:-0}" = "$want" ] && [ "$hrs" = "0 23" ]; then log "$STA fwd densify VERIFY PASS: $tc templates (=sel), mf hours $hrs"
  else log "$STA fwd densify VERIFY WARN: templates=$tc want=$want hours=[$hrs]"; fi
fi
log "$STA GPU PHASE done"
