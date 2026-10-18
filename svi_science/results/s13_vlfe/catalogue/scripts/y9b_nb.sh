#!/bin/bash
# Amendment Y9b: similarity-matched noise banks 2014-2017 (Y9 step-4 method), waits for Y10B_DONE; marker Y9B_DONE.
# Per-stage markers data/y9nb/<run><Y>/STEP_<m>_DONE; Y9_STEP4_DONE at the end. FAIL -> logs/q_main.log. Resume: re-run.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=2
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y9nbb $*" >> $L; }
for Y in 2014 2015 2016 2017; do
  R=prod; { [ $Y = 2016 ] || [ $Y = 2017 ]; } && R=prodq
  ND=data/y9nb/$R$Y; mkdir -p $ND
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  unset PROD_QC SCAN_REFBP; [ $R = prodq ] && export PROD_QC=1
  { [ $Y = 2011 ] || [ $Y = 2012 ]; } && export SCAN_REFBP=$SD/bp2_ref_v2
  export PROD_YEAR=$Y PILOT_SFX=_v2 SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  [ -f tables/y9nb_${Y}_banks.csv ] && [ -f $ND/banks.npz ] && touch $ND/STEP_build_DONE
  for m in build mf detect battery; do
    [ -f $ND/STEP_${m}_DONE ] && continue
    log "start $Y $m"
    python scripts/cat_y9nb.py $m >> logs/y9nb.out 2>&1 || { log "FAIL $Y $m (logs/y9nb.out)"; exit 1; }
    touch $ND/STEP_${m}_DONE; log "end $Y $m"
  done
done
unset PROD_YEAR PROD_QC SCAN_YEAR SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP
python scripts/cat_y9nb.py summary > logs/y9nb_summary.out 2>&1 || { log "FAIL summary"; exit 1; }
TZ=America/Los_Angeles date > Y9B_DONE; log "done"
