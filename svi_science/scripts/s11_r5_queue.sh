#!/bin/bash
# S11 round 5 (svi_science/results/s11_depth/PLAN.md §9 + Amendment 1): task-1 forms per region -> selection -> smooth Sc +
# Pcal + Q3 band per region -> Sc adoption -> Q1 per region (task 3 + Q1 band) -> 2019-11 per instrument (SVI) -> figures.
# Every stage writes a DONE marker in round5/ and is skipped when the marker exists (safe to restart after a reboot).
# <= 12 workers, 1 BLAS thread each. Usage: nohup bash s11_r5_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round5
mkdir -p $B
L=$B/queue.log
echo $$ > $B/queue.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
say "queue start (pid $$)"
run(){ m=$1; shift; if [ -f $B/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $B/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $B/$m ]; }
for R in SVI WA_NW; do S11_REGION=$R run MODEL_FORMS_${R}_DONE env S11_REGION=$R python s11_r5_model.py --stage forms; done
run SELECT_DONE python s11_r5_model.py --stage select
for R in SVI WA_NW; do run MODEL_SC_${R}_DONE env S11_REGION=$R python s11_r5_model.py --stage sc; done
run ADOPT_DONE python s11_r5_model.py --stage adopt
for R in SVI WA_NW; do run Q1_${R}_DONE env S11_REGION=$R python s11_r5_q1.py --nperm 500; done
run PGC_DONE env S11_REGION=SVI python s11_r5_pgc.py --workers 12
run FIGURES_DONE python s11_r5_figures.py
say "QUEUE_DONE"; touch $B/QUEUE_DONE
