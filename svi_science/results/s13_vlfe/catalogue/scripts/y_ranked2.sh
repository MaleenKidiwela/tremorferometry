#!/bin/bash
# Amendment Y8b ranked catalogue (fix of Y8): full battery columns for every agreement row (real, noise templates,
# scrambled draws A1/A2 0-1) under each year's scan-v2 env, monotone score, curves, sanity checks.
# Marker Y_RANKED2_DONE; per-year markers data/ranked/R2BAT_<Y>_DONE; FAIL lines -> logs/q_main.log. Resume: re-run.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=2
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Yranked2 $*" >> $L; }
for Y in 2018 2013 2016 2017 2010 2014 2015 2011 2012; do
  [ -f data/ranked/R2BAT_${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  unset SCAN_REFBP; { [ $Y = 2011 ] || [ $Y = 2012 ]; } && export SCAN_REFBP=$SD/bp2_ref_v2
  export SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  log "start battery $Y"
  python scripts/cat_ranked2.py battery $Y >> logs/y_ranked2.out 2>&1 || { log "FAIL battery $Y (see logs/y_ranked2.out)"; exit 1; }
  touch data/ranked/R2BAT_${Y}_DONE; log "end battery $Y"
done
unset SCAN_YEAR SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP
log "start score"
python scripts/cat_ranked2.py score > logs/y_ranked2_score.out 2>&1 || { log "FAIL score (see logs/y_ranked2_score.out)"; exit 1; }
TZ=America/Los_Angeles date > Y_RANKED2_DONE; log "done"
