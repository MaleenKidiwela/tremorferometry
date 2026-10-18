#!/bin/bash
# WA-SW column (PLAN §5): starts after the WA-NW queue ends (QUEUE_DONE or a halt), runs inventory -> windows -> measure ->
# kernels only (Q1/Q3 decided on cell counts afterwards). Markers in s11_depth/wa_sw/.
cd /home/jovyan/tremorferometry/svi_science/scripts
NW=/home/jovyan/tremorferometry/svi_science/results/s11_depth/wa_nw
O=/home/jovyan/tremorferometry/svi_science/results/s11_depth/wa_sw
mkdir -p $O/logs; L=$O/logs/queue.log
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
echo $$ > $O/logs/chain.pid
say "chain waiting for WA-NW queue (pid $(cat $NW/logs/queue.pid))"
while kill -0 $(cat $NW/logs/queue.pid) 2>/dev/null; do sleep 60; done
export S11_REGION=WA_SW
run(){ m=$1; shift; if [ -f $O/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $O/$m ] || { say "missing $m -> halt"; exit 1; }; }
run WA_INVENTORY_DONE python s11_wa_inventory.py
run WA_WINDOWS_DONE python s11_wa_windows.py --workers 14
run MEASURE_DONE python s11_measure.py --workers 14
run KERNEL_DONE python s11_kernel_mc.py --workers 14 --npart 100000
say "SW_CHAIN_DONE"; touch $O/SW_CHAIN_DONE
