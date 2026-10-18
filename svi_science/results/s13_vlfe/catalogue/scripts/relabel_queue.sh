#!/bin/bash
# Amendment U fix: lfe_lists lower bound (was 2011-05-31) -> re-score LFE labels for 2010 and 2011 (summary, R, prod summary).
# Waits for PRODR_2016_DONE (CPU). Marker RELABEL_DONE. FAIL lines -> logs/q_main.log.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') relabel $*" >> $L; }
until [ -f PRODR_2016_DONE ]; do pgrep -f prodr_queue.sh > /dev/null || { log "FAIL prodr_queue not running and no PRODR_2016_DONE"; exit 1; }; sleep 60; done
for Y in 2010 2011; do
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  unset SCAN_REFBP; [ $Y = 2011 ] && export SCAN_REFBP=$SD/bp2_ref_v2
  for f in tables/prod${Y}_summary_v2.json tables/prod${Y}_episodes_v2.csv tables/prod${Y}_R_verdict_v2.json; do [ -f $f ] && cp $f ${f%.*}_preRelabel.${f##*.}; done
  python scripts/cat_pilot.py summary >> logs/relabel.out 2>&1 || { log "FAIL $Y summary"; exit 1; }
  python scripts/cat_r.py >> logs/relabel.out 2>&1 || { log "FAIL $Y R"; exit 1; }
  python scripts/cat_prod.py summary $Y >> logs/relabel.out 2>&1 || { log "FAIL $Y psum"; exit 1; }
  log "$Y relabelled"
done
TZ=America/Los_Angeles date > RELABEL_DONE
