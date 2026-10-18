#!/bin/bash
# Amendment Q3/Q4 queue: after Q1_DONE -> scan year 2013 -> pilot 2013 -> scan year 2008 -> pilot 2008.
# Markers: catalogue/PILOT2013_DONE, PILOT2008_DONE. FAIL lines -> catalogue/logs/q_main.log. Resume: re-run.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') pilots $*" >> $L; }
until [ -f Q1_DONE ]; do pgrep -f q1_queue.sh > /dev/null || { log "FAIL q1 queue not running and no Q1_DONE"; exit 1; }; sleep 60; done
for Y in 2013 2008; do
  [ -f PILOT${Y}_DONE ] && continue
  bash ../scripts/scan_year_queue.sh $Y || { log "FAIL scan year $Y"; exit 1; }
  mkdir -p data/pilot$Y
  for m in cont bank mf detect battery inj_list inj summary; do
    [ -f data/pilot$Y/STEP_${m}_DONE ] && continue
    log "start $Y $m"
    PILOT=$Y SCAN_YEAR=$Y python scripts/cat_pilot.py $m >> logs/pilot$Y.out 2>&1 || { log "FAIL pilot $Y $m (see logs/pilot$Y.out)"; exit 1; }
    touch data/pilot$Y/STEP_${m}_DONE; log "end $Y $m"
  done
  TZ=America/Los_Angeles date > PILOT${Y}_DONE; log "pilot $Y done"
done
