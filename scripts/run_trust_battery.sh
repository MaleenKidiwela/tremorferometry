#!/bin/bash
# Family Trust Battery — one command per station.
# Usage: run_trust_battery.sh <NET> <STA> [REV_YEARS="2015 2016"] [KEEP_WF=0]
# Chain: (waveforms must exist or be downloading) -> reversed-template npz -> reversed densify (GPU,
# sample years) -> Tier-1 v2 (stack-vs-random + day/night vs reversed-fake calibration) -> verdict merge
# into data/family_trust_master.csv -> waveform cleanup (unless KEEP_WF=1) -> note stub.
set -u
NET=$1; STA=$2; REV_YEARS=${3:-"2015 2016"}; KEEP_WF=${4:-0}
S=$(echo $STA | tr A-Z a-z)
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/trust_${S}.log; : > $L
say(){ echo "=== $(date +%H:%M:%S) [$STA] $* ===" | tee -a $L; }

say "0. wait for waveforms (dir must appear, then download must finish)"
for i in $(seq 1 60); do [ -d data/waveforms/$NET.$STA ] && break; sleep 60; done
[ -d data/waveforms/$NET.$STA ] || { say "NO WAVEFORMS after 60 min — abort"; exit 1; }
until ! pgrep -f "[d]ownload_station.*$STA" >/dev/null 2>&1; do sleep 60; done
say "waveforms ready: $(du -sh data/waveforms/$NET.$STA | cut -f1)"

say "1. reversed-template npz"
[ -f data/${S}_reversed_templates.npz ] || $PY - <<EOF >> $L 2>&1
import numpy as np
z=np.load('data/${S}_pnsn_families_100km.npz',allow_pickle=True)
np.savez('data/${S}_reversed_templates.npz',**{k:np.asarray(z[k],float)[::-1].copy() for k in z.files})
print('reversed',len(z.files))
EOF

set -- $REV_YEARS; Y0=$1; Y1=${2:-$1}
say "2. reversed-template densify $Y0-$Y1 (GPU; the guaranteed-fake calibration)"
if [ $(ls data/mf_${S}rev_[12]*.csv 2>/dev/null | wc -l) -eq 0 ]; then
  $PY scripts/densify_launcher.py --templates-npz data/${S}_reversed_templates.npz \
    --summary-csv data/${S}_coverage_selection.summary.csv --min-snr 0 --network $NET --station $STA \
    --out-prefix mf_${S}rev_ --workers 20 --top-n 100 --max-raw-det 3000000 --despike-mad 8 \
    --start-year $Y0 --end-year $Y1 >> $L 2>&1
fi

say "3. Tier-1 v2 battery"
$PY scripts/family_trust_tier1.py --net $NET --sta $STA --rev-mf "data/mf_${S}rev_*.csv" \
  --nsamp 250 --nrand 2500 --nboot 200 --workers 16 >> $L 2>&1
grep -E "REAL fams|FAKES|verdicts" $L | tail -3

say "4. merge into master trust table"
$PY - <<EOF >> $L 2>&1
import pandas as pd, os
t=pd.read_csv('data/family_trust_tier1_${STA}.csv'); t=t[t.kind=='F'].copy(); t['station']='${STA}'
mst='data/family_trust_master.csv'
if os.path.exists(mst):
    m=pd.read_csv(mst); m=m[m.station!='${STA}']; t=pd.concat([m,t],ignore_index=True)
t.to_csv(mst,index=False)
print('master:',len(t),'rows,',t.station.nunique(),'stations')
EOF

if [ "$KEEP_WF" != "1" ] && [ ! -f data/.keep_wf_$STA ]; then
  say "5. waveform cleanup"
  rm -rf data/waveforms/$NET.$STA
else
  say "5. waveform cleanup SKIPPED (keep flag set for $STA)"
fi
say "TRUST BATTERY DONE — verdicts in data/family_trust_tier1_${STA}.csv (+ master)"
