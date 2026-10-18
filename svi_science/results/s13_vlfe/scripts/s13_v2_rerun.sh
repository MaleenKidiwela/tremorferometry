#!/bin/bash
# S13 Amendment 5 addendum (2026-10-14 ~7:35 pm PT): re-run the two wa2_ stacks started before the instrument-veto patch
# of s13_stack2_v2.py (wa2_between, wa2_2012-02; killed 7:30 pm PT). Waits for s13_v2_queue.sh to end, then runs them
# (2 parallel), then writes V2_QUEUE_DONE if every v2 stack exists. Run: nohup bash scripts/s13_v2_rerun.sh > logs/v2_rerun.out 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2
LOG=logs/v2_queue.log
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
QP=$(cat logs/v2_queue.pid)
while kill -0 $QP 2>/dev/null; do sleep 60; done
WA_STA=$(grep '^WA_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
echo "$(ts) rerun (veto patch) start" >> $LOG
for b in between 2012-02; do
  tag=wa2_$b
  [ -f STACK_${tag}_DONE ] && continue
  ( python scripts/s13_stack2_v2.py data/rows_v2/wa_$b.parquet $tag $WA_STA --kmax 5 > logs/${tag}.out 2>&1 && echo "$(ts) $tag done (rerun)" >> $LOG || echo "$(ts) $tag FAILED (rerun)" >> $LOG ) &
done
wait
miss=0; for f in data/rows_v2/wa_*.parquet; do b=$(basename $f .parquet); [ -f STACK_wa2_${b#wa_}_DONE ] || miss=$((miss+1)); done
for t in pbin_svi_2004-07 pbin_svi_2005-09 pbin_svi_2012-08 pbin_svi_2013-09 svibe2; do [ -f STACK_${t}_DONE ] || miss=$((miss+1)); done
[ $miss -eq 0 ] && touch V2_QUEUE_DONE || echo "$(ts) NOT done after rerun: $miss stacks missing" >> $LOG
echo "$(ts) rerun end" >> $LOG
