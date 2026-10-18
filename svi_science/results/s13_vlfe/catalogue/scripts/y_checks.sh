#!/bin/bash
# Amendment Y3/Y4 checks after Y_2018W2B_DONE. Marker Y_CHECKS_DONE. FAIL lines -> logs/q_main.log.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Ychecks $*" >> $L; }
until [ -f Y_2018W2B_DONE ]; do pgrep -f "scripts/y_queue.sh" > /dev/null || { log "FAIL y_queue not running and no Y_2018W2B_DONE"; exit 1; }; sleep 60; done
OMP_NUM_THREADS=1 python scripts/cat_ychecks.py > logs/y_checks.out 2>&1 || { log "FAIL cat_ychecks (see logs/y_checks.out)"; exit 1; }
TZ=America/Los_Angeles date > Y_CHECKS_DONE; log "done"
