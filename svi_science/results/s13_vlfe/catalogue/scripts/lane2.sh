#!/bin/bash
# 2-worker lane after Y9B_DONE: Y12b 2022 re-run, then Y13.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
L=logs/q_main.log; log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') lane2 $*" >> $L; }
until [ -f Y9B_DONE ]; do pgrep -f "scripts/yq_master.sh" > /dev/null || { [ -f Y9B_DONE ] || { log "FAIL master queue not running and no Y9B_DONE"; exit 1; }; }; sleep 120; done
[ -f Y12B_2022_DONE ] || bash scripts/y12b_2022.sh || { log "FAIL y12b 2022"; exit 1; }
bash scripts/y13_queue.sh || { log "FAIL y13"; exit 1; }
