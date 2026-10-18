#!/bin/bash
# S11 Revision 2 queue: added-family windows -> measure -> kernels -> Q3 analysis. Stages skip if DONE exists.
cd /home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/s11_depth
L=$O/logs/queue_q3.log
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
run(){ m=$1; shift; if [ -f $O/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $O/logs/queue_q3_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $O/$m ] || { say "missing $m -> halt"; exit 1; }; }
say "queue start (pid $$)"
run EXT_WINDOWS_DONE python s11_ext_windows.py --workers 12
run MEASURE_ext_DONE python s11_measure.py --workers 12 --tag _ext --src $O/tables/ext
run KERNEL_ext_DONE python s11_kernel_mc.py --workers 12 --npart 100000 --tag _ext --ext
run Q3_DONE python s11_q3.py --workers 12
say "QUEUE_Q3_DONE"; touch $O/QUEUE_Q3_DONE
