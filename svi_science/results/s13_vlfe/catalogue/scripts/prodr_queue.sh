#!/bin/bash
# Amendment U: production 2016-2018 with bank = 2011 bank + 2018 relay bank. Markers PRODR_<Y>_DONE.
# FAIL lines -> logs/q_main.log. Resume: re-run (STEP markers in data/prodr<Y>/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') prodr $*" >> $L; }
if [ ! -f data/inj/STEP_relay_sources_DONE ]; then
  python scripts/cat_inj.py relay_sources >> logs/prodr.out 2>&1 || { log "FAIL relay_sources"; exit 1; }
  touch data/inj/STEP_relay_sources_DONE; log "relay sources done"
fi
for Y in 2018 2017 2016; do
  [ -f PRODR_${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; PD=data/prodr$Y; mkdir -p $PD
  export SCAN_YEAR=$Y PROD_YEAR=$Y PROD_RELAY=1 PILOT_SFX=_v2 SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  st(){ local m=$1; shift; [ -f $PD/STEP_${m}_DONE ] && return 0; log "start $Y $m"
    "$@" >> logs/prodr$Y.out 2>&1 || { log "FAIL $Y $m (see logs/prodr$Y.out)"; exit 1; }; touch $PD/STEP_${m}_DONE; log "end $Y $m"; }
  for m in bank mf detect battery inj_list inj summary; do st $m python scripts/cat_pilot.py $m; done
  st R python scripts/cat_r.py
  st psum python scripts/cat_prod.py summary $Y
  TZ=America/Los_Angeles date > PRODR_${Y}_DONE; log "year $Y done"
done
