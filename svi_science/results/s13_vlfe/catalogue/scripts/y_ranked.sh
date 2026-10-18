#!/bin/bash
# Amendment Y8 ranked catalogue product. Waits for Y_CHECKS_DONE (2018 W2b run needed). Marker Y_RANKED_DONE.
# FAIL lines -> logs/q_main.log. Resume: re-run (per-year outputs are rewritten).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Yranked $*" >> $L; }
until [ -f Y_CHECKS_DONE ]; do pgrep -f "scripts/y_checks.sh" > /dev/null || { log "FAIL y_checks not running and no Y_CHECKS_DONE"; exit 1; }; sleep 60; done
[ -f data/ranked/score_model.json ] || python scripts/cat_ranked.py fit > logs/y_ranked.out 2>&1 || { log "FAIL fit"; exit 1; }
for Y in 2010 2011 2012 2013 2014 2015 2016 2017 2018; do
  python scripts/cat_ranked.py year $Y >> logs/y_ranked.out 2>&1 || { log "FAIL year $Y"; exit 1; }; log "year $Y done"
done
TZ=America/Los_Angeles date > Y_RANKED_DONE; log "done"
