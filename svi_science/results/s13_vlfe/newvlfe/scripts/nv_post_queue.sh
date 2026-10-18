#!/bin/bash
# newvlfe stage 3-4: after BATTERY_DONE -> nv_summary.py (separation, verdict, raw figures) -> nv_survivors.py
# (new-candidate columns + waveform figures). Marker newvlfe/POST_DONE; FAIL lines in logs/queue.log.
N=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/newvlfe; cd $N
export OMP_NUM_THREADS=1
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
while [ ! -f $N/BATTERY_DONE ]; do grep -q "FAIL" logs/queue.log && { log "post: FAIL seen upstream, stop"; exit 1; }; sleep 60; done
log "post: summary start"
python scripts/nv_summary.py > logs/nv_summary.out 2>&1 || { log "FAIL summary"; exit 1; }
log "post: survivors start"
python scripts/nv_survivors.py > logs/nv_survivors.out 2>&1 || { log "FAIL survivors"; exit 1; }
TZ=America/Los_Angeles date > $N/POST_DONE; log "post done"
