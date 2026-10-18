#!/bin/bash
# S12-LVL LFE side, PLAN §7: WA-NW recharge confirmation (+ SVI by era) -> 1-3 Hz reflection restack -> figures.
# Idempotent (DONE markers in ../). Log ../logs/queue.log
cd /home/jovyan/tremorferometry/svi_science/results/s12_lvl/lfe/scripts
O=..; L=$O/logs/queue.log
echo $$ > $O/queue3.pid
run() { local mk=$1; shift; if [ -f $O/$mk ]; then echo "[$(date)] skip $mk" >> $L; return 0; fi
        echo "[$(date)] start $mk: $*" >> $L; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?
        echo "[$(date)] end $mk rc=$rc" >> $L; return $rc; }
run BMEAS_WA_NW_DONE env S11_REGION=WA_NW python s12_b_measure.py --workers 12
run R2_WA_NW_DONE env S11_REGION=WA_NW python s12_recharge2.py
run BMEAS_SVI_DONE env S11_REGION=SVI python s12_b_measure.py --workers 12
run R2_SVI_DONE env S11_REGION=SVI python s12_recharge2.py
run LOWBAND_BUILD_DONE python s12_lowband.py build --workers 7 --nmax 100
run REFLECT_LB_DONE python s12_lowband.py scan --nscr 1000
run REFLECT_INJECT_LB_DONE env S12_TAG=_lb python s12_reflect_inject.py
python s12_figures.py recharge lowband >> $O/logs/queue_stage.out 2>&1
touch $O/QUEUE3_DONE; echo "[$(date)] QUEUE3_DONE" >> $L
