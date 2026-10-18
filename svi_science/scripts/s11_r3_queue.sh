#!/bin/bash
# S11 round 3 queue (PLAN.md §6): kernels -> REDPy measurement -> task 1 (kernel-robust deep term) -> task 2 (Q3 with
# REDPy estimators) -> task 3 (along-strike zones), SVI and WA-NW. Idempotent: every stage skips if its DONE marker exists.
# Usage: nohup bash s11_r3_queue.sh > /dev/null 2>&1 &      Log: s11_depth/round3/queue.log
cd /home/jovyan/tremorferometry/svi_science/scripts
L=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round3
mkdir -p $L/svi/logs $L/wa_nw/logs
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L/queue.log; }
echo $$ > $L/queue.pid
say "queue start (pid $$)"
run(){ reg=$1; mk=$2; shift 2; sub=svi; [ $reg = WA_NW ] && sub=wa_nw
  if [ -f $L/$sub/$mk ]; then say "skip $reg $mk"; return 0; fi
  say "start $reg $*"; S11_REGION=$reg "$@" >> $L/$sub/logs/stage.out 2>&1; rc=$?; say "end $reg $* rc=$rc"
  [ -f $L/$sub/$mk ] || { say "missing $sub/$mk -> halt"; exit 1; }; }
# 0) kernels: wait for a running kernel chain, then (re)run it (it skips finished parts)
while pgrep -f s11_r3_kernel_chain.sh >/dev/null; do sleep 60; done
bash s11_r3_kernel_chain.sh
for f in svi/KERNEL_DONE svi/KERNEL_ext_DONE wa_nw/KERNEL_DONE; do [ -f $L/$f ] || { say "missing $f (KERNEL_STOP?) -> halt"; exit 1; }; done
# 1) REDPy measurement (cc weights, ccw store, graph + bootstrap)
run SVI MEASURE_DONE python s11_r3_measure.py --workers 12
run SVI MEASURE_ext_DONE python s11_r3_measure.py --workers 12 --ext
run WA_NW MEASURE_DONE python s11_r3_measure.py --workers 12
# 2) task 1
run SVI KERR_DONE python s11_r3_kerr.py --workers 12
run WA_NW KERR_DONE python s11_r3_kerr.py --workers 12
# 3) task 2
run SVI R3Q3_DONE python s11_r3_q3.py --workers 12
run WA_NW R3Q3_DONE python s11_r3_q3.py --workers 12
# 4) task 3
run SVI ZONES_DONE python s11_r3_zones.py --workers 12
run WA_NW ZONES_DONE python s11_r3_zones.py --workers 12
say "QUEUE_DONE"; touch $L/QUEUE_DONE
