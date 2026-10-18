#!/bin/bash
# S12c scale-up (PLAN §5): 2003-2026 preprocessing for PGC, VGZ, SNB (skips finished station-years), then visibility
# and the superposed-epoch analysis on all ETS. Idempotent: stages with DONE markers are skipped.
cd /home/jovyan/tremorferometry
WS=svi_science/results/s12_lvl/acf_pgc
log() { echo "[$(TZ=America/Los_Angeles date '+%F %H:%M PT')] $*"; }
while [ ! -f $WS/PILOT_ANA_DONE ]; do sleep 60; done   # keep total CPU <= 8
if [ ! -f $WS/FULL_PRE_DONE ]; then
  for S in CN.PGC CN.VGZ CN.SNB; do
    Y1=2026; [ $S = CN.SNB ] && Y1=2021
    log "pre $S"; python $WS/scripts/acf01_daily.py --stations $S --years 2003-$Y1 --workers 8 || { log "FAIL pre $S"; PREFAIL=1; }
  done
  # guard added 2026-10-14: no marker if a run failed or any station-year came out empty (0 days)
  if [ -z "$PREFAIL" ] && ! find $WS/data -name "*.npz" -size -10k | grep -q .; then date > $WS/FULL_PRE_DONE
  else log "pre NOT complete -> stop"; exit 1; fi
fi
if [ ! -f $WS/FULL_VIS_DONE ]; then
  for M in tremorfree none; do
    log "vis $M"; python $WS/scripts/acf02_visibility.py --stations CN.PGC,CN.VGZ,CN.SNB --years 2003-2026 --mask $M --tag full
  done
  date > $WS/FULL_VIS_DONE
fi
run() { K=$WS/EPOCH_full_$1_$2_$3_DONE; [ -f $K ] && return
  Y1=2026; [ $1 = CN.SNB ] && Y1=2021
  log "epoch $1 $2 $3"
  python $WS/scripts/acf03_epoch.py --station $1 --years 2003-$Y1 --band $2 --mask $3 --tag full --ndraw 2000 \
    > $WS/logs/epoch_full_$1_$2_$3.log 2>&1 && date > $K; }
for S in CN.PGC CN.VGZ CN.SNB; do     # 4 epoch runs at a time (4 CPUs, ~4 GB each)
  for B in b05_2 b1_4; do run $S $B tremorfree & run $S $B none & done
  wait
done
python $WS/scripts/acf04_compare.py --tag full > $WS/logs/compare_full.log 2>&1
log "queue finished"; date > $WS/FULL_QUEUE_DONE
