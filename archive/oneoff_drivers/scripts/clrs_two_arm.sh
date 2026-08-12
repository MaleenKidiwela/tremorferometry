#!/bin/bash
# Merlin's decisive test: does the picker beat plain PNSN-windowing? Cluster the picker arm (top-30k by P(LFE))
# AND a control arm (top-30k by SNR, picker-blind) through discover_gpu; compare families. Sequential (one GPU).
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
log(){ echo "$(date +%H:%M) [clrs-2arm] $*" | tee -a logs/clrs_two_arm.log; }
run_arm(){
  local NAME=$1 CAND=$2 OUT=$3
  log "$NAME arm: discover_gpu clustering @40Hz"
  PYTHONPATH=src $PY scripts/discover_gpu.py --station CLRS --wfdir data/waveforms \
    --candidates "$CAND" --out "$OUT" --fs 40 --workers 24 > logs/clrs_gpu_${NAME}.log 2>&1 \
    || { log "$NAME arm FAILED (see logs/clrs_gpu_${NAME}.log)"; return 1; }
  local nf=$(($(wc -l < "${OUT%.npz}.summary.csv" 2>/dev/null)-1))
  log "$NAME arm DONE: $nf families (>=3 members, >=3 years) -> ${OUT%.npz}.summary.csv"
}
run_arm picker  data/clrs_cand_filtered.parquet    data/clrs_picker_families.npz
run_arm control data/clrs_control_filtered.parquet data/clrs_control_families.npz
log "CLRS TWO-ARM CLUSTERING COMPLETE"
