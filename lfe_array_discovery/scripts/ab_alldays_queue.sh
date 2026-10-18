#!/usr/bin/env bash
# =====================================================================================
# ALL-DAYS detection pass -- follow-on queue to ab_royer_queue.sh (never edits or signals it).
#
# Launch (foreground Bash call):
#   cd /home/jovyan/tremorferometry/lfe_array_discovery && \
#     nohup setsid -f bash scripts/ab_alldays_queue.sh >> logs/ab_alldays_queue.log 2>&1 < /dev/null
#
# What it waits for (cheap `sleep 600` loop): each job's FINAL GROWTH output from the chain,
#   real  : results/ab_royer/<SUB>/grow/summary.json          (written once stage 5 of <SUB> ends)
#   nullK : results/ab_royer/<SUB>/null<K>/grow/summary.json  (stage 7, realisation K)
# plus the all-days prefetch (results/ab_royer/alldays/PREFETCH_DONE).  Jobs run in the order
# real NW, SW, NVI, then null1 NW/SW/NVI, null2 ..., null3 ...; each as soon as it is ready.  A
# job whose growth ended without final templates is recorded as skipped.  After each job the report
# is rewritten (ab_alldays_report.py).  The queue ends when all 12 jobs are done/skipped, or when
# AB_ROYER_DONE exists and nothing else can become ready.
#
# Markers (results/):  AB_ALLDAYS_PID, AB_ALLDAYS_HEARTBEAT, AB_ALLDAYS_DONE, AB_ALLDAYS_FAILED
# Stop it with `kill -TERM $(cat results/AB_ALLDAYS_PID)` -- never pkill -f.
# =====================================================================================
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
RES=$PWD/results
OUT=$RES/ab_royer
AD=$OUT/alldays
LOCKS=$PWD/locks
mkdir -p "$AD/jobs" "$LOCKS" logs
PY=$(command -v python)
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export AB_NULL=0

exec 8>"$LOCKS/ab_alldays.lock"
if ! flock -n 8; then echo "[alldays] another instance holds the lock; exiting"; exit 0; fi
echo $$ > "$RES/AB_ALLDAYS_PID"
rm -f "$RES/AB_ALLDAYS_FAILED"
STAGE=init

# heartbeat: a background subshell that does NOT keep the lock fd (closed first, so its `sleep`
# children never hold the flock) and is killed, with its sleep, in the EXIT trap
( exec 8>&-; trap 'kill $SP 2>/dev/null; exit 0' TERM
  while :; do date -u +%FT%TZ > "$RES/AB_ALLDAYS_HEARTBEAT"; sleep 300 & SP=$!; wait $SP; done ) &
HB=$!

finish() {
  rc=$?
  kill -TERM "$HB" 2>/dev/null
  for c in $(pgrep -P "$HB" 2>/dev/null); do kill "$c" 2>/dev/null; done
  if [ ! -f "$RES/AB_ALLDAYS_DONE" ]; then
    { echo "stage=$STAGE"; echo "rc=$rc"; date -u +%FT%TZ; } > "$RES/AB_ALLDAYS_FAILED"
    echo "[alldays] FAILED in $STAGE (rc=$rc)"
  fi
}
trap finish EXIT
trap 'STAGE="$STAGE(signalled)"; exit 143' TERM INT HUP

log() { echo "[alldays] $(date -u +%FT%TZ) $*"; }

JOBS="real:NW real:SW real:NVI null1:NW null1:SW null1:NVI null2:NW null2:SW null2:NVI null3:NW null3:SW null3:NVI"
growdir() { if [ "$1" = real ]; then echo "$OUT/$2/grow"; else echo "$OUT/$2/$1/grow"; fi; }

# ---- stage P: prefetch + 40 sps cache of every pool station-day in the array lifetime
#      (<= 8 FDSN connections, <= 12 cache processes; the chain shares the 32-CPU pod)
if [ ! -f "$AD/PREFETCH_DONE" ]; then
  STAGE=prefetch
  log "start prefetch"
  "$PY" scripts/ab_alldays_fetch.py all --workers 8 --cache-workers 12 >> logs/ab_alldays_fetch.log 2>&1 \
    || { log "prefetch failed"; exit 1; }
  log "done prefetch"
fi

log "waiting for the chain's growth summaries"
WAITED=0
while :; do
  STAGE=wait
  pending=0; ran=0
  for J in $JOBS; do
    K=${J%%:*}; S=${J##*:}
    [ -e "$AD/jobs/${K}_${S}.done" ] || [ -e "$AD/jobs/${K}_${S}.skipped" ] && continue
    pending=$((pending+1))
    GD=$(growdir "$K" "$S")
    [ -f "$GD/summary.json" ] || continue
    [ -f "$AD/PREFETCH_DONE" ] || continue
    if ! ls "$GD"/tpl_it*.npz >/dev/null 2>&1; then
      log "skip $K $S: growth has no templates"; date -u +%FT%TZ > "$AD/jobs/${K}_${S}.skipped"
      continue
    fi
    STAGE="alldays_${K}_${S}"
    log "start $STAGE"
    "$PY" scripts/ab_alldays.py run --sub "$S" --kind "$K" >> "logs/ab_alldays_${K}_${S}.log" 2>&1 \
      || { log "$STAGE failed"; exit 1; }
    date -u +%FT%TZ > "$AD/jobs/${K}_${S}.done"
    STAGE="report_after_${K}_${S}"
    "$PY" scripts/ab_alldays_report.py >> logs/ab_alldays_report.log 2>&1 || { log "report failed"; exit 1; }
    log "done alldays_${K}_${S}"
    ran=1; pending=$((pending-1))
    break          # re-scan from the top: real jobs keep priority over nulls
  done
  if [ "$pending" -eq 0 ]; then break; fi
  if [ "$ran" -eq 1 ]; then continue; fi
  if [ -f "$RES/AB_ROYER_DONE" ]; then
    # the chain is finished: anything still without a growth summary will never get one
    for J in $JOBS; do
      K=${J%%:*}; S=${J##*:}
      [ -e "$AD/jobs/${K}_${S}.done" ] || [ -e "$AD/jobs/${K}_${S}.skipped" ] && continue
      GD=$(growdir "$K" "$S")
      if [ ! -f "$GD/summary.json" ]; then
        log "skip $K $S: chain DONE without its growth summary"
        date -u +%FT%TZ > "$AD/jobs/${K}_${S}.skipped"
      fi
    done
  fi
  WAITED=$((WAITED+600))
  if [ "$WAITED" -ge $((45*86400)) ]; then STAGE="wait_timeout_45d"; exit 1; fi
  sleep 600 8>&-
done

STAGE=final_report
"$PY" scripts/ab_alldays_report.py >> logs/ab_alldays_report.log 2>&1 || exit 1
date -u +%FT%TZ > "$RES/AB_ALLDAYS_DONE"
log "ALL-DAYS PASS DONE"
exit 0
