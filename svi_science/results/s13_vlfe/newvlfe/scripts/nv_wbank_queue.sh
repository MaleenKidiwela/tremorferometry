#!/bin/bash
# Western-group candidate templates, relay acceptance on 2018 (report only). Marker newvlfe/STACKW_DONE; FAIL -> logs/queue.log.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; N=$D/newvlfe; cd $D
export OMP_NUM_THREADS=1 PROD_YEAR=2018 PROD_RELAY=1 PILOT_SFX=_v2 SCAN_YEAR=2018 SCAN_RUN=real_v2 SCAN_BP=$D/data/scan_2018/bp_2018v2 \
  SCAN_GRID=$D/data/scan_2018/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 CAT_NW=2
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
for m in mf loo battery verdict; do
  [ -f $N/data/wbank/STEP_${m}_DONE ] && continue
  log "wbank start $m"
  python $N/scripts/nv_wbank.py $m >> $N/logs/nv_wbank.out 2>&1 || { log "FAIL wbank $m"; exit 1; }
  touch $N/data/wbank/STEP_${m}_DONE; log "wbank end $m"
done
TZ=America/Los_Angeles date > $N/STACKW_DONE; log "wbank done"
