#!/bin/bash
# S13 Amendment 6 queue (2026-10-14 PT). Idempotent (DONE markers). Run:
#   nohup bash scripts/s13_v3_queue.sh > logs/v3_queue.out 2>&1 &      log: logs/v3_queue.log, PID logs/v3_queue.pid
# WA step 4: wa3_* (built majors in data/rows_v3/wa_2*.parquet + wa3_between), kmax 5, grids '' (+-8 s) + wide_ + fix2_
# SVI step 5: nmt_svi_* (net_mad terciles within P >= 0.9), kmax 1, 200 draws, neighbour pool = the episode's pooled rows
# <= 6 stacks in parallel, OMP_NUM_THREADS=1. Done when V3_QUEUE_DONE exists; then the queue runs
# python scripts/s13_v3_summary.py (logs/v3_summary.out) if that script exists.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2
LOG=logs/v3_queue.log
echo $$ > logs/v3_queue.pid
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
SVI_STA=$(grep '^SVI_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
WA_STA=$(grep '^WA_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
echo "$(ts) v3 queue start" >> $LOG
run(){  # tag rows stations extra-args...
  local tag=$1 rows=$2 sta=$3; shift 3
  [ -f STACK_${tag}_DONE ] && return 0
  python scripts/s13_stack2_v2.py $rows $tag $sta "$@" > logs/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG
}
JOBS=()
JOBS+=("wa3_between|data/rows_v3/wa_between.parquet|WA|--kmax 5 --grids ,wide_,fix2_")
JOBS+=("nmt_svi_2004-07|data/rows_v3/svi_nmt_2004-07.parquet|SVI|--kmax 1 --nb-pool data/rows/svi_2004-07.parquet")
JOBS+=("nmt_svi_2005-09|data/rows_v3/svi_nmt_2005-09.parquet|SVI|--kmax 1 --nb-pool data/rows/svi_2005-09.parquet")
for f in $(ls data/rows_v3/wa_2*.parquet); do
  b=$(basename $f .parquet)
  JOBS+=("wa3_${b#wa_}|$f|WA|--kmax 5 --grids ,wide_,fix2_")
done
JOBS+=("nmt_svi_2012-08|data/rows_v3/svi_nmt_2012-08.parquet|SVI|--kmax 1 --nb-pool data/rows/svi_2012-08.parquet")
JOBS+=("nmt_svi_2013-09|data/rows_v3/svi_nmt_2013-09.parquet|SVI|--kmax 1 --nb-pool data/rows/svi_2013-09.parquet")
for j in "${JOBS[@]}"; do
  IFS='|' read tag rows reg extra <<< "$j"
  [ -f STACK_${tag}_DONE ] && continue
  if [ $reg = WA ]; then sta=$WA_STA; else sta=$SVI_STA; fi
  run $tag $rows $sta $extra &
  sleep 2
  while [ $(jobs -rp | wc -l) -ge 6 ]; do sleep 15; done
done
wait
miss=0
for j in "${JOBS[@]}"; do IFS='|' read tag rows reg extra <<< "$j"; [ -f STACK_${tag}_DONE ] || miss=$((miss+1)); done
if [ $miss -eq 0 ]; then
  touch V3_QUEUE_DONE; echo "$(ts) all ${#JOBS[@]} stacks done -> V3_QUEUE_DONE" >> $LOG
  [ -f scripts/s13_v3_summary.py ] && python scripts/s13_v3_summary.py > logs/v3_summary.out 2>&1 && echo "$(ts) summary ok" >> $LOG || echo "$(ts) summary not run or failed (logs/v3_summary.out)" >> $LOG
else
  echo "$(ts) NOT done: $miss of ${#JOBS[@]} stacks missing (re-run the same nohup command)" >> $LOG
fi
echo "$(ts) v3 queue end" >> $LOG
