#!/bin/bash
# S11 round 6 (svi_science/results/s11_depth/PLAN.md §10): task-1 forms + cause analysis per region -> adoption ->
# task-2 size regressions per region -> figures. Every stage writes a DONE marker in round6/ and is skipped when the
# marker exists (safe to restart after a reboot). <= 12 workers (serial here), 1 BLAS thread each.
# Usage: nohup bash s11_r6_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round6
mkdir -p $B
L=$B/queue.log
echo $$ > $B/queue.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
say "queue start (pid $$)"
run(){ m=$1; shift; if [ -f $B/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $B/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $B/$m ]; }
for R in SVI WA_NW; do run MODEL_${R}_DONE env S11_REGION=$R python s11_r6_model.py --stage forms; done
run ADOPT_DONE python s11_r6_model.py --stage adopt
for R in SVI WA_NW; do run SIZE_${R}_DONE env S11_REGION=$R python s11_r6_size.py; done
run FIGURES_DONE python s11_r6_figures.py
say "QUEUE_DONE"; touch $B/QUEUE_DONE
