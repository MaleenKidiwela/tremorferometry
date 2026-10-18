#!/bin/bash
# Master <= 2-worker queue (after Y10B_DONE): Y11 inj38 -> Y12 pilot 2021 -> Y11 ladder 2019-21 -> Y9b banks 2014-17.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
L=logs/q_main.log; log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') master $*" >> $L; }
until [ -f Y10B_DONE ]; do pgrep -f "scripts/y10b_queue.sh" > /dev/null || { log "FAIL y10b not running and no Y10B_DONE"; exit 1; }; sleep 60; done
bash scripts/y11_queue.sh --inj-only || { log "FAIL y11 inj"; exit 1; }
bash scripts/y12_pilot.sh || { log "FAIL y12 pilot"; exit 1; }
bash scripts/y11_queue.sh || { log "FAIL y11 ladder/y9b"; exit 1; }
log "master done"
