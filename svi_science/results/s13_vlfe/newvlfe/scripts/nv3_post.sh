#!/bin/bash
# NV3 post: after NV3_BAT_DONE -> 2018 batteries with the NV3 columns (realW/polBW/shA3/shA7/shA11; 2 workers) ->
# nv3_summary.py -> nv3_labels.py (blind figures + key). Marker NV3_DONE; FAIL lines in logs/queue.log.
N=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/newvlfe; cd $N
export OMP_NUM_THREADS=1
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
while [ ! -f NV3_BAT_DONE ]; do grep -q "FAIL NV3" logs/queue.log && { log "NV3 post: upstream FAIL, stop"; exit 1; }; sleep 60; done
for R in realW polBW shA3 shA7 shA11; do
  [ -f data/bat_2018_${R}.parquet ] && continue
  log "NV3 battery start 2018 $R"
  SCAN_YEAR=2018 python scripts/nv_battery.py 2018 $R --workers 2 >> logs/nv_battery_2018.out 2>&1 || { log "FAIL NV3 battery 2018 $R"; exit 1; }
  log "NV3 battery end 2018 $R"
done
log "NV3 summary start"; python scripts/nv3_summary.py > logs/nv3_summary.out 2>&1 || { log "FAIL NV3 summary"; exit 1; }
log "NV3 labels start"; python scripts/nv3_labels.py > logs/nv3_labels.out 2>&1 || { log "FAIL NV3 labels"; exit 1; }
TZ=America/Los_Angeles date > NV3_DONE; log "NV3 done"
