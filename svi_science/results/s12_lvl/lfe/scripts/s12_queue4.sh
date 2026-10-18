#!/bin/bash
# S12-LVL LFE side, PLAN §7: rerun of the 1-3 Hz reflection scan only (crashed 10-13 on NaN S-P pairs; fixed 10-14 in
# s12_reflect.resid_logenv_P). Build (LOWBAND_BUILD_DONE) is reused. Idempotent (DONE markers in ../). Log ../logs/queue.log
cd /home/jovyan/tremorferometry/svi_science/results/s12_lvl/lfe/scripts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
O=..; L=$O/logs/queue.log
echo $$ > $O/queue4.pid
run() { local mk=$1; shift; if [ -f $O/$mk ]; then echo "[$(date)] skip $mk" >> $L; return 0; fi
        echo "[$(date)] start $mk: $*" >> $L; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?
        echo "[$(date)] end $mk rc=$rc" >> $L; return $rc; }
run REFLECT_LB_DONE python s12_lowband.py scan --nscr 1000
run REFLECT_INJECT_LB_DONE env S12_TAG=_lb python s12_reflect_inject.py
python s12_figures.py lowband >> $O/logs/queue_stage.out 2>&1
touch $O/QUEUE4_DONE; echo "[$(date)] QUEUE4_DONE" >> $L
