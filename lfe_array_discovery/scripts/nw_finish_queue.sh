#!/usr/bin/env bash
# =====================================================================================
# NW finish queue (2026-09-26): (1) converge growth from it8 (ab_grow_cont.py, Bostock 2012 rule),
# (2) all-station final templates on the literal NW 25-station pool (ab_allsta.py), (3) ab_hypo.py
# on the all-station common-time stacks, (4) comparison + figure.
#
# Launch:  cd /home/jovyan/tremorferometry/lfe_array_discovery && \
#   nohup setsid -f bash scripts/nw_finish_queue.sh >> logs/nw_finish_queue.log 2>&1 < /dev/null
# Markers: results/NW_FINISH_DONE | NW_FINISH_FAILED | NW_FINISH_HEARTBEAT, pid in NW_FINISH_PID.
# Stop:    kill -TERM $(cat results/NW_FINISH_PID)   (never pkill -f)
# Stages are guarded by their outputs, so a relaunch resumes.
# =====================================================================================
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
RES=$PWD/results
NW=$RES/ab_royer/NW
GROW=$NW/grow
AS=$NW/hypo_allsta
mkdir -p logs locks "$AS/null317"
PY=$(command -v python)
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}

exec 8>locks/nw_finish.lock
if ! flock -n 8; then echo "[nwf] lock held; exiting"; exit 0; fi
echo $$ > "$RES/NW_FINISH_PID"
rm -f "$RES/NW_FINISH_FAILED"

STAGE=init
HB=""
finish() {
  rc=$?
  [ -n "$HB" ] && kill -TERM "$HB" 2>/dev/null
  if [ "$rc" -ne 0 ] || [ ! -f "$RES/NW_FINISH_DONE" ]; then
    { echo "stage=$STAGE"; echo "rc=$rc"; date -u +%FT%TZ; } > "$RES/NW_FINISH_FAILED"
    echo "[nwf] FAILED in $STAGE rc=$rc"
  fi
}
trap finish EXIT
trap 'STAGE="$STAGE(signalled)"; exit 143' TERM INT HUP
# heartbeat subshell: does NOT inherit the lock fd
( while :; do date -u +%FT%TZ > "$RES/NW_FINISH_HEARTBEAT"; sleep 120; done ) 8>&- &
HB=$!

log() { echo "[nwf] $(date -u +%FT%TZ) $*"; }

STAGE=grow_cont
if [ ! -f "$GROW/FINAL_ITER" ]; then
  log "growth continuation from it8"
  "$PY" scripts/ab_grow_cont.py --sub NW --from-iter 8 --max-more 6 --workers 20 || exit 1
fi
LAST=$(cat "$GROW/FINAL_ITER")
log "final growth iteration it$LAST"

STAGE=allsta_stack
if [ ! -f "$AS/stacks_loc.npz" ] || [ ! -f "$AS/allsta_stack_stats.csv" ]; then
  "$PY" scripts/ab_allsta.py --sub NW --det "$GROW/det_it${LAST}.parquet" \
    --tpl-sta "$GROW/tpl_it${LAST}_stations.csv" --out "$AS" --passes 0 --workers 20 || exit 1
fi
# per-detection x-corr refinement REJECTED in the smoke test (a 317-s time-shifted null gets the
# same SNR gain: 2.23 -> 7.99 vs real 2.24 -> 7.72; lags uniform in +-0.5 s) -> --passes 0, i.e.
# the common-time stacks (moveout embedded) are the templates.  Station-level NULL stacks (every
# detection moved +317 s) as a real-vs-null column per (family, station):
STAGE=allsta_null
if [ ! -f "$AS/null317/allsta_stack_stats.csv" ]; then
  "$PY" scripts/ab_allsta.py --sub NW --det "$GROW/det_it${LAST}.parquet" \
    --tpl-sta "$GROW/tpl_it${LAST}_stations.csv" --out "$AS/null317" --passes 0 --workers 20 \
    --null-shift 317 || exit 1
fi

STAGE=hypo_allsta
if [ ! -f "$AS/loc_summary.csv" ]; then
  "$PY" scripts/ab_hypo.py --sub NW --det "$GROW/det_it${LAST}.parquet" --out "$AS" --boot 60 || exit 1
fi
STAGE=report
"$PY" scripts/ab_allsta_report.py --sub NW --main allsta \
  --runs "it8_hypo=$NW/hypo,allsta=$AS" \
  --fig "$AS/NW_catalogue_allsta.png" > "$AS/allsta_report.txt" 2>&1 || exit 1
cat "$AS/allsta_report.txt"
date -u +%FT%TZ > "$RES/NW_FINISH_DONE"
STAGE=done
log "DONE"
