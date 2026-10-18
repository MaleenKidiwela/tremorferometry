#!/bin/bash
# Q3 step 0 re-stacks (PLAN Amendment 4). Same stacker/stations/data as the S13 table run of each episode
# (2012-08 = refit_, first-pass SVI_STA; 2018-06, 2024-08 = scalef_, fill list). <= 3 in parallel. Marker Q3_RESTACK_DONE.
S13=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
S15=/home/jovyan/tremorferometry/svi_science/results/s15_vlfe_science
cd $S13
export OMP_NUM_THREADS=1
export S13_DISP=$S13/data/disp1sps_v2
LOG=$S15/tables/q3_restack.log
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
SVI_STA=$(sed -n 9p scripts/s13_scale_queue.sh | cut -d= -f2)
SVI_FILL=$(grep '^SVI_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
echo "$(ts) start" >> $LOG
for kind in ex in; do
  for e in 2012-08 2018-06 2024-08; do
    tag=s15${kind}_svi_$e
    [ -f STACK_${tag}_DONE ] && continue
    STA=$([ $e = 2012-08 ] && echo $SVI_STA || echo $SVI_FILL)
    ( python scripts/s13_stack2.py $S15/data/rows/$tag.parquet $tag $STA > $S15/tables/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG ) &
    while [ $(jobs -rp | wc -l) -ge 3 ]; do sleep 10; done
  done
done
wait
echo "$(ts) end" >> $LOG
touch $S15/Q3_RESTACK_DONE
