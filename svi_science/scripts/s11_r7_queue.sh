#!/bin/bash
# S11 round 7 (svi_science/results/s11_depth/PLAN.md §11): per region front timing -> zones -> joint -> hint; then the hint
# verdict and figures. Every stage writes a DONE marker in round7/ and is skipped when the marker exists (safe to restart
# after a reboot). <= 12 workers, 1 BLAS thread each.
# Usage: nohup bash s11_r7_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round7
mkdir -p $B
L=$B/queue.log
echo $$ > $B/queue.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
say "queue start (pid $$)"
run(){ m=$1; shift; if [ -f $B/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $B/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $B/$m ]; }
for R in SVI WA_NW; do
  run FRONT_${R}_near_DONE env S11_REGION=$R python s11_r7_front.py --part near --workers 12
  run FRONT_${R}_far_DONE env S11_REGION=$R python s11_r7_front.py --part far --workers 12
  run FRONT_${R}_DONE env S11_REGION=$R python s11_r7_front.py --part analyse
  run ZONES_${R}_DONE env S11_REGION=$R python s11_r7_zones.py
  run JOINT_${R}_DONE env S11_REGION=$R python s11_r7_joint.py
  run HINT_${R}_DONE env S11_REGION=$R python s11_r7_hint.py
done
run HINT_VERDICT_DONE python s11_r7_hint.py --stage verdict
run FIGURES_DONE python s11_r7_figures.py
say "QUEUE_DONE"; touch $B/QUEUE_DONE
