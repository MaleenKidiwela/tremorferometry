#!/bin/bash
# S12-LVL LFE side, coordinator tasks A-D (PLAN.md §6): novel measurement per region -> analysis. Idempotent (DONE markers).
cd /home/jovyan/tremorferometry/svi_science/results/s12_lvl/lfe/scripts
O=..; L=$O/logs/queue.log
echo $$ > $O/queue2.pid
run() { local mk=$1; shift; if [ -f $O/$mk ]; then echo "[$(date)] skip $mk" >> $L; return 0; fi
        echo "[$(date)] start $mk: $*" >> $L; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?
        echo "[$(date)] end $mk rc=$rc" >> $L; return $rc; }
run REFLECT_INJECT_DONE python s12_reflect_inject.py
run NOVEL_SVI_DONE env S11_REGION=SVI python s12_novel_measure.py --workers 12
run NOVEL_WA_NW_DONE env S11_REGION=WA_NW python s12_novel_measure.py --workers 12
run NOVELAN_SVI_DONE env S11_REGION=SVI python s12_novel.py
run NOVELAN_WA_NW_DONE env S11_REGION=WA_NW python s12_novel.py
touch $O/QUEUE2_DONE; echo "[$(date)] QUEUE2_DONE" >> $L
