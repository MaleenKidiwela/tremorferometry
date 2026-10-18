#!/bin/bash
# S11 round 6, §10 Amendment 1 (post-result for task 1): forms E4lu / E4yu per region -> adoption -> task-2 column -> figures redo.
# DONE markers in round6/ (skipped when present). Usage: nohup bash s11_r6a_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round6
L=$B/queue.log
echo $$ > $B/queue_amend.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
say "amendment queue start (pid $$)"
run(){ m=$1; shift; if [ -f $B/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $B/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $B/$m ]; }
for R in SVI WA_NW; do run AMEND_${R}_DONE env S11_REGION=$R python s11_r6_model.py --stage amend; done
run ADOPT_AMEND_DONE python s11_r6_model.py --stage adopt_amend
for R in SVI WA_NW; do run SIZE_AMEND_${R}_DONE env S11_REGION=$R python s11_r6_size.py --amend; done
say "AMEND_QUEUE_DONE"; touch $B/AMEND_QUEUE_DONE
