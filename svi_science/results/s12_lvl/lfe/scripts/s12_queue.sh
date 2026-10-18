#!/bin/bash
# S12-LVL LFE-side queue (PLAN.md): reflection scan -> epoch SVI -> pool SVI -> epoch WA-NW -> pool WA-NW.
# Idempotent: stages with a DONE marker are skipped. Log ../logs/queue.log
cd /home/jovyan/tremorferometry/svi_science/results/s12_lvl/lfe/scripts
O=..; L=$O/logs/queue.log
echo $$ > $O/queue.pid
run() { local mk=$1; shift; if [ -f $O/$mk ]; then echo "[$(date)] skip $mk" >> $L; return 0; fi
        echo "[$(date)] start $mk: $*" >> $L; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?
        echo "[$(date)] end $mk rc=$rc" >> $L; return $rc; }
run REFLECT_DONE python s12_reflect.py --nscr 1000
S11_REGION=SVI run EPOCH_SVI_DONE env S11_REGION=SVI python s12_epoch.py --workers 12 --nscr 200
run POOL_SVI_DONE env S11_REGION=SVI python s12_epoch_pool.py
run EPOCH_WA_NW_DONE env S11_REGION=WA_NW python s12_epoch.py --workers 12 --nscr 200
run POOL_WA_NW_DONE env S11_REGION=WA_NW python s12_epoch_pool.py
touch $O/QUEUE_DONE; echo "[$(date)] QUEUE_DONE" >> $L
