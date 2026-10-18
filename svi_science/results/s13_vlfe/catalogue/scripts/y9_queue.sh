#!/bin/bash
# Amendment Y9 steps 5a (emergence columns), 2-3 (ladder + tables), 5b (level-5 figures). Markers Y9_STEP2_DONE, Y9_STEP3_DONE,
# Y9_STEP5_DONE (figures; eye labels added afterwards). FAIL -> logs/q_main.log. Resume: re-run (done years skipped).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y9 $*" >> $L; }
log "start emergence"
for Y in 2010 2011 2012 2013 2014 2015 2016 2017 2018; do [ -f data/ranked/y9_emerge_$Y.parquet ] || echo $Y; done | \
  xargs -P 2 -I{} sh -c 'python scripts/cat_y9.py emerge {} >> logs/y9.out 2>&1 || echo "$(TZ=America/Los_Angeles date "+%F %H:%M:%S PT") Y9 FAIL emerge {}" >> logs/q_main.log'
for Y in 2010 2011 2012 2013 2014 2015 2016 2017 2018; do [ -f data/ranked/y9_emerge_$Y.parquet ] || { log "FAIL emerge $Y missing"; exit 1; }; done
log "end emergence; start ladder"
python scripts/cat_y9.py ladder > logs/y9_ladder.out 2>&1 || { log "FAIL ladder (logs/y9_ladder.out)"; exit 1; }
TZ=America/Los_Angeles date > Y9_STEP2_DONE; TZ=America/Los_Angeles date > Y9_STEP3_DONE; log "ladder done"
python scripts/cat_y9.py figs > logs/y9_figs.out 2>&1 || { log "FAIL figs (logs/y9_figs.out)"; exit 1; }
TZ=America/Los_Angeles date > Y9_STEP5_DONE; log "figs done"
