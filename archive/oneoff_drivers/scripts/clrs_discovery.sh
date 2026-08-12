#!/bin/bash
# CLRS picker-discovery validation — 40 Hz THROUGHOUT (broadband picker trained at 40 Hz; native 100 Hz HHZ
# would give wrong-Nyquist features). Stops after scoring so the P(LFE) distribution can be validated before
# clustering. Usage: nohup bash scripts/clrs_discovery.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
log(){ echo "$(date +%H:%M) [clrs] $*" | tee -a logs/clrs_discovery.log; }
log "stage1: PNSN-driven detection @ fs40 (redo — was mistakenly fs100)"
PYTHONPATH=src $PY scripts/discover_nllb_pnsn_driven.py --station CLRS --wfdir data/waveforms \
  --bbox 47.92 49.72 -125.48 -122.78 --pnsn catalogs/pnsn_tremor_cascadia_full.csv --fs 40 \
  --candidates-only --candidates-out data/clrs_pnsn_candidates_100km.parquet --workers 16 \
  > logs/cand_clrs.log 2>&1 || { log "detection FAILED"; exit 1; }
log "stage1 done: $($PY -c "import pandas as pd;print(len(pd.read_parquet('data/clrs_pnsn_candidates_100km.parquet')))") candidates"
cp -f lfe_features/models/picker_broadband_pgc.joblib lfe_features/models/tremor_picker_clrs.joblib
log "stage2: broadband picker score @ target-fs 40 (thr 0.5)"
PYTHONPATH=src $PY lfe_features/score_candidates.py --net CN --sta CLRS \
  --cand data/clrs_pnsn_candidates_100km.parquet --y0 2014 --y1 2026 --thr 0.5 --target-fs 40 --workers 16 \
  > logs/score_clrs.log 2>&1 || { log "score FAILED"; exit 2; }
log "stage2 done: $(grep -E 'scored|BASELINE' logs/score_clrs.log | tail -1)"
log "CLRS SCORING COMPLETE @40Hz — validate P(LFE) vs B926, then cluster"
