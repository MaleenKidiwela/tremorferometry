#!/bin/bash
# NV3 batteries (1 worker while the NV3 GPU queue runs, else 2). Marker NV3_BAT_DONE.
N=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/newvlfe; cd $N
export OMP_NUM_THREADS=1
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
for Y in 2016 2019 2011 2013; do
  for R in realW polBW shA3W shA7W shA11W; do
    [ -f data/bat_${Y}_${R}.parquet ] && continue
    RD=$R; [ $R = realW ] && RD=polBW
    while [ ! -f data/runs/$Y/$RD/DONE ]; do grep -q "FAIL NV3" logs/queue.log && { log "NV3 battery: upstream FAIL, stop"; exit 1; }; sleep 60; done
    W=2; [ -f $N/NV3_SCAN_DONE ] || W=1
    log "NV3 battery start $Y $R (workers $W)"
    SCAN_YEAR=$Y python scripts/nv_battery.py $Y $R --workers $W >> logs/nv_battery_$Y.out 2>&1 || { log "FAIL NV3 battery $Y $R"; exit 1; }
    log "NV3 battery end $Y $R"
  done
done
TZ=America/Los_Angeles date > $N/NV3_BAT_DONE; log "NV3 batteries done"
