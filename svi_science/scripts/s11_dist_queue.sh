#!/bin/bash
# S11 distant-station lever (svi_science/results/s11_depth/distant/PLAN.md): per region (SVI, WA_NW)
# inventory -> [DIST_STOP: region ends] -> windows -> measure -> kernels -> [KERNEL_STOP: region ends] -> invert -> gates;
# then figures. Stages skip if their DONE marker exists (safe to restart after a reboot). <= 10 workers.
# Usage: nohup bash s11_dist_queue.sh > /dev/null 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/scripts
B=/home/jovyan/tremorferometry/svi_science/results/s11_depth/distant
mkdir -p $B
L=$B/queue.log
echo $$ > $B/queue.pid
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
say "queue start (pid $$)"
for R in SVI WA_NW; do
  SUB=$(echo $R | tr 'A-Z' 'a-z'); O=$B/$SUB; mkdir -p $O/logs
  export S11_REGION=$R
  run(){ m=$1; shift; if [ -f $O/$m ]; then say "$R skip $m"; return 0; fi; say "$R start $*"; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?; say "$R end $* rc=$rc"; [ -f $O/$m ]; }
  run INVENTORY_DONE python s11_dist_inventory.py --workers 10 || { say "$R inventory failed -> next region"; continue; }
  if [ -f $O/DIST_STOP ]; then say "$R DIST_STOP (expected usable cells < 100) -> region ends"; continue; fi
  run WINDOWS_DONE python s11_dist_windows.py --workers 10 || { say "$R windows failed -> next region"; continue; }
  run MEASURE_DONE python s11_dist_measure.py --workers 10 || { say "$R measure failed -> next region"; continue; }
  if [ ! -f $O/KERNEL_DONE ] && [ ! -f $O/KERNEL_STOP ]; then run KERNEL_DONE python s11_dist_kernel.py --workers 10 --npart 100000; fi
  if [ -f $O/KERNEL_STOP ]; then say "$R KERNEL_STOP -> region ends"; continue; fi
  run INVERT_DONE python s11_dist_invert.py --workers 10 || { say "$R invert failed -> next region"; continue; }
  run GATES_DONE python s11_dist_gates.py || say "$R gates failed"
done
unset S11_REGION
if [ ! -f $B/FIGURES_DONE ]; then say "start figures"; python s11_dist_figures.py >> $B/figures.out 2>&1; say "end figures rc=$?"; fi
say "QUEUE_DONE"; touch $B/QUEUE_DONE
