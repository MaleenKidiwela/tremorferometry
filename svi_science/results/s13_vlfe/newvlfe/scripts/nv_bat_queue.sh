#!/bin/bash
# newvlfe stage 2: battery for real + 7 control catalogues per year. Waits for each control run's DONE.
# 1 worker while the GPU scan queue runs, else 2 (rule: <= 2 CPU workers). Marker newvlfe/BATTERY_DONE.
N=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/newvlfe; cd $N
export OMP_NUM_THREADS=1
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
for Y in 2013 2018; do
  for R in real polB shA3 shA7 shA11; do
    [ -f data/bat_${Y}_${R}.parquet ] && continue
    if [ $R != real ]; then while [ ! -f data/runs/$Y/$R/DONE ]; do
      grep -q "FAIL" logs/queue.log && { log "battery: scan FAIL seen, stop"; exit 1; }; sleep 30; done; fi
    W=2; [ -f $N/SCANCTRL_DONE ] || W=1
    log "battery start $Y $R (workers $W)"
    SCAN_YEAR=$Y python scripts/nv_battery.py $Y $R --workers $W >> logs/nv_battery_$Y.out 2>&1 || { log "FAIL battery $Y $R"; exit 1; }
    log "battery end $Y $R"
  done
done
TZ=America/Los_Angeles date > $N/BATTERY_DONE; log "battery done"
