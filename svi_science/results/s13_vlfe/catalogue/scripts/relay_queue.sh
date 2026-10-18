#!/bin/bash
# Amendment T2 relay banks 2016-2018 (after the 2010-2018 production queue). Markers RELAY_<Y>_DONE.
# FAIL lines -> logs/q_main.log. Resume: re-run (STEP markers in data/prod<Y>/relay/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') relay $*" >> $L; }
until [ -f PROD_2010_DONE ] && [ -f PROD_2016_DONE ] && [ -f PROD_2017_DONE ] && [ -f PROD_2018_DONE ]; do
  pgrep -f "scripts/prod_queue.sh" > /dev/null || { log "FAIL prod_queue not running and production markers missing"; exit 1; }; sleep 120; done
for Y in 2018 2017 2016; do
  [ -f RELAY_${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; RD=data/prod$Y/relay; mkdir -p $RD
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  st(){ local m=$1; shift; [ -f $RD/STEP_${m}_DONE ] && return 0; log "start $Y $m"
    "$@" >> logs/relay$Y.out 2>&1 || { log "FAIL $Y $m (see logs/relay$Y.out)"; exit 1; }; touch $RD/STEP_${m}_DONE; log "end $Y $m"; }
  st bank python scripts/cat_relay.py bank
  if [ ! -f tables/prod${Y}_relay_bank.csv ]; then log "$Y: no relay bank (see tables/prod${Y}_relay_verdict.json)"; TZ=America/Los_Angeles date > RELAY_${Y}_DONE; continue; fi
  for m in mf loo battery verdict; do st $m python scripts/cat_relay.py $m; done
  TZ=America/Los_Angeles date > RELAY_${Y}_DONE; log "relay $Y done"
done
