#!/bin/bash
# Amendment W2b: diagnostics for all production years, exclusion table, re-runs 2013 (check), 2016, 2017 (2011 bank and relay bank).
# Markers W_<NAME><Y>_DONE, W_DONE. FAIL lines -> logs/q_main.log. Resume: re-run.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') W $*" >> $L; }
for y in 2011 2012 2014 2015 2018; do [ -f tables/diag_noise_$y.csv ] || python scripts/cat_diag.py $y > /dev/null 2>&1 || { log "FAIL diag $y"; exit 1; }; done
python scripts/qc_exclusions.py > logs/qc_exclusions.out 2>&1 || { log "FAIL qc_exclusions"; exit 1; }
log "exclusions written"
for spec in "2013 0" "2016 1" "2017 1" "2016 0" "2017 0"; do
  set -- $spec; Y=$1; R=$2
  NM=prod; [ $R = 1 ] && NM=prodr
  [ -f W_${NM}q${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; PD=data/${NM}q$Y; mkdir -p $PD
  export SCAN_YEAR=$Y PROD_YEAR=$Y PROD_QC=1 PILOT_SFX=_v2 SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  if [ $R = 1 ]; then export PROD_RELAY=1; else unset PROD_RELAY; fi
  st(){ local m=$1; shift; [ -f $PD/STEP_${m}_DONE ] && return 0; log "start ${NM}q$Y $m"
    "$@" >> logs/${NM}q$Y.out 2>&1 || { log "FAIL ${NM}q$Y $m (see logs/${NM}q$Y.out)"; exit 1; }; touch $PD/STEP_${m}_DONE; log "end ${NM}q$Y $m"; }
  for m in bank mf detect battery inj_list inj summary; do st $m python scripts/cat_pilot.py $m; done
  st R python scripts/cat_r.py
  st psum python scripts/cat_prod.py summary $Y
  TZ=America/Los_Angeles date > W_${NM}q${Y}_DONE; log "${NM}q$Y done"
done
TZ=America/Los_Angeles date > W_DONE
