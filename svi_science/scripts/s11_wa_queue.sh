#!/bin/bash
# S11 step 7 (PLAN §5) region queue: inventory -> windows -> measure -> kernels -> v2 -> major list -> q3 -> WA extras.
# Usage: S11_REGION=WA_NW nohup bash s11_wa_queue.sh > /dev/null 2>&1 &     Stages skip if their DONE marker exists.
cd /home/jovyan/tremorferometry/svi_science/scripts
export S11_REGION=${S11_REGION:-WA_NW}
SUB=$(echo $S11_REGION | tr 'A-Z' 'a-z')
O=/home/jovyan/tremorferometry/svi_science/results/s11_depth/$SUB
mkdir -p $O/logs
L=$O/logs/queue.log
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
run(){ m=$1; shift; if [ -f $O/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $O/$m ] || { say "missing $m -> halt"; exit 1; }; }
echo $$ > $O/logs/queue.pid
say "queue start (pid $$, region $S11_REGION)"
run WA_INVENTORY_DONE python s11_wa_inventory.py
run WA_WINDOWS_DONE python s11_wa_windows.py --workers 14
run MEASURE_DONE python s11_measure.py --workers 14
run KERNEL_DONE python s11_kernel_mc.py --workers 14 --npart 100000
if [ -f $O/KERNEL_STOP ]; then say "KERNEL_STOP -> halt"; exit 1; fi
run V2_DONE python s11_v2.py --workers 14
run MAJLIST_DONE python s11_wa_majlist.py
run Q3_DONE python s11_q3.py --workers 8
if [ -f s11_wa_compare.py ]; then run WA_COMPARE_DONE python s11_wa_compare.py --workers 8; fi
say "QUEUE_DONE"; touch $O/QUEUE_DONE
