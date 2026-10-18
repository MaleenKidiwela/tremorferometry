#!/bin/bash
# Amendment S production 2010-2018 (scan v2, 2011 bank, weighted, n = 800). Waits for SCANV2_2018_DONE (CPU cap),
# runs carriers_prep.sh, then years in order. Per-year marker PROD_<Y>_DONE; FAIL lines -> logs/q_main.log.
# Resume: re-run (stages with markers in data/prod<Y>/ are skipped).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') prod $*" >> $L; }
until [ -f SCANV2_2018_DONE ]; do pgrep -f "v2_year.sh 2018" > /dev/null || { log "FAIL v2_year 2018 not running and no SCANV2_2018_DONE"; exit 1; }; sleep 60; done
[ -f CARRIERS_DONE ] || bash scripts/carriers_prep.sh || { log "FAIL carriers_prep"; exit 1; }
for Y in 2011 2012 2013 2018 2014 2015 2016 2017 2010; do
  [ -f PROD_${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  if [ ! -f $SD/runs/real/DONE ]; then bash $D/scripts/scan_year_queue.sh $Y || { log "FAIL scan year $Y"; exit 1; }; fi
  PD=data/prod$Y; mkdir -p $PD
  st(){ local m=$1; shift; [ -f $PD/STEP_${m}_DONE ] && return 0; log "start $Y $m"
    "$@" >> logs/prod$Y.out 2>&1 || { log "FAIL $Y $m (see logs/prod$Y.out)"; exit 1; }; touch $PD/STEP_${m}_DONE; log "end $Y $m"; }
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2
  unset SCAN_REFBP; { [ $Y = 2011 ] || [ $Y = 2012 ]; } && export SCAN_REFBP=$SD/bp2_ref_v2
  MISS=$(python scripts/cat_prod.py missing $Y)
  if [ -n "$MISS" ] && [ ! -f $PD/STEP_prepmiss_DONE ]; then
    log "start $Y prep missing carriers $MISS"
    (cd $D && S13_DISP=$D/data/disp1sps_v3 python scripts/scan_prep_ref.py $MISS $Y-01-01 $Y-12-31 --workers 8 > /dev/null 2>&1) || { log "FAIL $Y prepmiss"; exit 1; }
    touch $PD/STEP_prepmiss_DONE; [ -f $PD/STEP_stations_DONE ] && mv $PD/STEP_stations_DONE $PD/STEP_stations_DONE.before_prepmiss
    log "end $Y prepmiss"
  fi
  st stations python scripts/cat_prod.py stations $Y
  if [ ! -f $SD/runs/real_v2/DONE ] || [ $Y = 2011 ] || [ $Y = 2012 ]; then
    st v2setup python scripts/cat_scanv2.py setup $Y
    st v2detect bash -c "cd $D && SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
  fi
  export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  for m in cont bank mf detect battery inj_list inj summary; do st $m python scripts/cat_pilot.py $m; done
  st R python scripts/cat_r.py
  st psum python scripts/cat_prod.py summary $Y
  unset SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP
  TZ=America/Los_Angeles date > PROD_${Y}_DONE; log "year $Y done"
done
log "production 2010-2018 done"
