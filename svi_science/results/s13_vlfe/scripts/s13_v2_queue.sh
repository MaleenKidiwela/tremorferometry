#!/bin/bash
# S13 Amendment 5 queue (2026-10-14 PT). Idempotent (DONE markers). Run:
#   nohup bash scripts/s13_v2_queue.sh > logs/v2_queue.out 2>&1 &      log: logs/v2_queue.log, PID logs/v2_queue.pid
# A: SVI step-4 P-bin/thinning (pbin_*) and step-5 between variants (svibe2), <=2 parallel, start at once (no prep needed)
# B: prep the extra WA days (+-6 d for the +-5-day controls) -> data/disp1sps_v2 (skips existing) -> PREP_wa2_DONE
# C: WA stacks wa2_* (54 episodes + between), kmax 5; <=4 of mine in parallel while the fill/refit queues run, <=6 after
# Done when V2_QUEUE_DONE exists. Then: python scripts/s13_v2_summary.py
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2
LOG=logs/v2_queue.log
echo $$ > logs/v2_queue.pid
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
SVI_STA=$(grep '^SVI_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
WA_STA=$(grep '^WA_STA=' scripts/s13_fill_queue.sh | cut -d= -f2)
echo "$(ts) v2 queue start" >> $LOG
run(){  # tag rows stations extra-args...
  local tag=$1 rows=$2 sta=$3; shift 3
  [ -f STACK_${tag}_DONE ] && return 0
  python scripts/s13_stack2_v2.py $rows $tag $sta "$@" > logs/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG
}
limit(){ if pgrep -f "s13_fill_queue.sh|s13_refit_queue.sh" > /dev/null; then echo 4; else echo 6; fi; }
# ---- A (background, <=2 parallel)
(
  n=0
  for ep in 2004-07 2005-09 2012-08 2013-09; do
    run pbin_svi_$ep data/rows_v2/svi_pbin_$ep.parquet $SVI_STA --kmax 1 --nb-pool data/rows/svi_$ep.parquet &
    n=$((n+1)); if [ $n -ge 2 ]; then wait -n; n=$((n-1)); fi
  done
  run svibe2 data/rows_v2/svi_be_v2.parquet $SVI_STA --kmax 1 &
  wait
) &
PA=$!
# ---- B
if [ ! -f PREP_wa2_DONE ]; then
  echo "$(ts) prep wa2 start" >> $LOG
  python scripts/s13_prep.py $WA_STA x x --days-file tables/wa2_prep_days.csv --workers 12 > logs/prep_wa2.out 2>&1 && ! grep -q "ModuleNotFoundError\|Traceback" logs/prep_wa2.out && touch PREP_wa2_DONE
  echo "$(ts) prep wa2 end ($(grep -c ' ok ' logs/prep_wa2.out) ok)" >> $LOG
fi
[ -f PREP_wa2_DONE ] || { echo "$(ts) prep wa2 FAILED, stop" >> $LOG; wait $PA; exit 1; }
# ---- C (between first: the longest)
for f in data/rows_v2/wa_between.parquet $(ls data/rows_v2/wa_2*.parquet); do
  b=$(basename $f .parquet); tag=wa2_${b#wa_}
  [ -f STACK_${tag}_DONE ] && continue
  run $tag $f $WA_STA --kmax 5 &
  sleep 2
  while [ $(jobs -rp | grep -v "^$PA\$" | wc -l) -ge $(limit) ]; do sleep 20; done
done
wait
miss=0; for f in data/rows_v2/wa_*.parquet; do b=$(basename $f .parquet); [ -f STACK_wa2_${b#wa_}_DONE ] || miss=$((miss+1)); done
for t in pbin_svi_2004-07 pbin_svi_2005-09 pbin_svi_2012-08 pbin_svi_2013-09 svibe2; do [ -f STACK_${t}_DONE ] || miss=$((miss+1)); done
[ $miss -eq 0 ] && touch V2_QUEUE_DONE || echo "$(ts) NOT done: $miss stacks missing" >> $LOG
echo "$(ts) v2 queue end" >> $LOG
