#!/bin/bash
# S13 Amendment 7 queue (2026-10-15 PT): between-ETS stacks at candidate times. Idempotent (markers). <= 3 CPU workers.
#   nohup bash scripts/s13_be3_queue.sh > logs/be3_queue.out 2>&1 &     log logs/be3_queue.log, PID logs/be3_queue.pid
# Steps: PNSN envelopes + origins -> rows -> prep (v2 recipe, data/disp1sps_v2, skips existing) -> 4 stacks (<= 3 parallel)
# Done when BE3_QUEUE_DONE exists. Then: python scripts/s13_be3_summary.py
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
cd $D
export OMP_NUM_THREADS=1 S13_DISP=$D/data/disp1sps_v2
L=$D/logs/be3_queue.log
echo $$ > logs/be3_queue.pid
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
STA=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,UW.LRIV
W=${BE3_WORKERS:-3}
log "queue start (workers $W)"
[ -s tables/be3_pnsn_origins.csv ] || { python scripts/s13_be3_env.py --workers $W > logs/be3_env.out 2>&1; log "env rc $?"; }
[ -s data/rows_be3/be3_pnsn.parquet ] || { python scripts/s13_units_be3.py rows > logs/be3_rows.out 2>&1; log "rows rc $?"; }
if [ ! -f PREP_be3_DONE ]; then
  log "prep start"
  python scripts/s13_prep.py $STA x x --days-file tables/be3_prep_days.csv --workers $W > logs/prep_be3.out 2>&1 && \
    ! grep -q "Traceback" logs/prep_be3.out && touch PREP_be3_DONE
  log "prep end ($(grep -c ' ok ' logs/prep_be3.out) ok, $(grep -c 'exists' logs/prep_be3.out) existing)"
fi
[ -f PREP_be3_DONE ] || { log "prep FAILED, stop"; exit 1; }
run(){ # tag
  [ -f STACK_$1_DONE ] && return 0
  python scripts/s13_stack2_v2.py data/rows_be3/$1.parquet $1 $STA --kmax 1 > logs/$1.out 2>&1 && log "$1 done" || log "$1 FAILED"
}
n=0
for t in be3_c010 be3_c030 be3_hi be3_pnsn; do
  run $t &
  n=$((n+1)); if [ $n -ge $W ]; then wait -n; n=$((n-1)); fi
done
wait
miss=0; for t in be3_c010 be3_c030 be3_hi be3_pnsn; do [ -f STACK_${t}_DONE ] || miss=$((miss+1)); done
[ $miss -eq 0 ] && touch BE3_QUEUE_DONE || log "NOT done: $miss stacks missing"
log "queue end"
