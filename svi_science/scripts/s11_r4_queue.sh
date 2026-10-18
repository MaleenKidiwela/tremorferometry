#!/bin/bash
# S11 round 4 (svi_science/results/s11_depth/PLAN.md §8): per region (SVI, WA_NW) S1 measure (near [+ext] + far) -> S2 far
# kernels -> S3 lapse-selection recovery -> S4 combined model + injection suite; then S5 2019-11 deep dive (SVI) -> S6 figures.
# Every stage writes a DONE marker in round4/ and is skipped when the marker exists (safe to restart after a reboot).
# <= 12 workers, 1 BLAS thread per worker (lapse stage: 4 threads, runs alone). Usage: nohup bash s11_r4_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round4
mkdir -p $B
L=$B/queue.log
echo $$ > $B/queue.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
while pgrep -f r4_s12.sh > /dev/null; do sleep 20; done        # first-launch helper for S1/S2 (same scripts) still running
say "queue start (pid $$)"
run(){ m=$1; shift; if [ -f $B/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $B/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $B/$m ]; }
for R in SVI WA_NW; do
  export S11_REGION=$R
  run MEASURE_${R}_DONE python s11_r4_measure.py --set near --workers 12
  if [ "$R" = "SVI" ]; then run MEASURE_${R}_ext_DONE python s11_r4_measure.py --set ext --workers 12; fi
  run MEASURE_${R}_far_DONE python s11_r4_measure.py --set far --workers 12
  if ! ls $B/FARKERNEL_${R}_* > /dev/null 2>&1; then run FARKERNEL_${R}_DONE python s11_r4_farkernel.py --workers 12; fi
  run LAPSE_${R}_DONE env OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 python s11_r4_lapse.py --ndraw 300 --nperm 500
  run COMBINED_${R}_DONE python s11_r4_combined.py --workers 12
done
export S11_REGION=SVI
run DIVE_DONE python s11_r4_dive.py --workers 12
unset S11_REGION
run FIGURES_DONE python s11_r4_figures.py
say "QUEUE_DONE"; touch $B/QUEUE_DONE
