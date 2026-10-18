#!/bin/bash
# S12c PLAN amendment B: echo amplitude / similarity stage, after the full run (waits for FULL_QUEUE_DONE).
# 3 stations x 2 bands x 2 masks, 4 at a time; markers AMP_full_*_DONE; all done = AMP_QUEUE_DONE (after the
# auto-report is appended to RESULTS.md). Idempotent.
cd /home/jovyan/tremorferometry
WS=svi_science/results/s12_lvl/acf_pgc
log() { echo "[$(TZ=America/Los_Angeles date '+%F %H:%M PT')] $*"; }
while [ ! -f $WS/FULL_QUEUE_DONE ]; do sleep 60; done
run() { K=$WS/AMP_full_$1_$2_$3_DONE; [ -f $K ] && return
  Y1=2026; [ $1 = CN.SNB ] && Y1=2021
  log "amp $1 $2 $3"
  python $WS/scripts/acf05_amp.py --station $1 --years 2003-$Y1 --band $2 --mask $3 --tag full --ndraw 2000 \
    > $WS/logs/amp_full_$1_$2_$3.log 2>&1 && date > $K; }
for S in CN.PGC CN.VGZ CN.SNB; do
  for B in b05_2 b1_4; do run $S $B tremorfree & run $S $B none & done
  wait
done
if [ ! -f $WS/AMP_QUEUE_DONE ]; then python $WS/scripts/acf06_report.py > $WS/logs/report_full.log 2>&1; fi
log "amp queue finished"; date > $WS/AMP_QUEUE_DONE
