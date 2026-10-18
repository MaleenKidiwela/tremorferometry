#!/bin/bash
# S13 Amendment 4 v2 guarded SVI re-runs (2026-10-14 ~6:25 pm PT; WA withdrawn pending origin fix). <=4 parallel.
# 1 non-fill SVI row files with a v2 drop (tables/g2_svi_refit_list.txt) -> tag refit_*, first-pass station list.
# 2 fill-affected SVI row files (tables/g2_svi_fill_list.txt) -> tag scalef_*, fill station list (adds CN.GOBB), filled data.
# Run: nohup bash scripts/s13_refit_queue.sh > logs/refit_queue.out 2>&1 &   Marker REFIT_QUEUE_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2
LOG=logs/refit_queue.log
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
SVI_STA=$(sed -n 9p scripts/s13_scale_queue.sh | cut -d= -f2)
SVI_FILL=$(grep '^SVI_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
echo "$(ts) refit queue (v2) start" >> $LOG
n=0
for job in $(sed 's/^/refit:/' tables/g2_svi_refit_list.txt) $(sed 's/^/scalef:/' tables/g2_svi_fill_list.txt); do
  pre=${job%%:*}; b=${job#*:}; tag=${pre}_$b
  [ -f STACK_${tag}_DONE ] && continue
  STA=$([ $pre = refit ] && echo $SVI_STA || echo $SVI_FILL)
  ( python scripts/s13_stack2.py data/rows/$b.parquet $tag $STA > logs/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG ) &
  n=$((n+1)); if [ $n -ge 4 ]; then wait -n; n=$((n-1)); fi
done
wait
miss=0
for b in $(cat tables/g2_svi_refit_list.txt); do [ -f STACK_refit_${b}_DONE ] || miss=$((miss+1)); done
for b in $(cat tables/g2_svi_fill_list.txt); do [ -f STACK_scalef_${b}_DONE ] || miss=$((miss+1)); done
[ $miss -eq 0 ] && touch REFIT_QUEUE_DONE || echo "$(ts) NOT done: $miss stacks missing" >> $LOG
echo "$(ts) refit queue (v2) end" >> $LOG
