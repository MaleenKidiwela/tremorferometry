#!/usr/bin/env bash
# 2026-09-29: full-lag (±15 s) Royer chain for NVI then SW, exactly as scripts/nw_lag15_queue.sh (NW_lag15, 65 real):
#   ab_cluster.py --lag-max 15 --lag-mode zeropad on the *_mute autocorrelation detections (auto_window reused)
#   -> ab_grow_fast.py --merge most, 14 iters, tol 0.01, prep 6, restack 4, runaway guard (>5 M detections)
#   -> ab_hypo --boot 60 -> ab_restack28 (workers 4) -> review sheets -> ab_nomerge_score.py (auto_score.csv, n_sites_3km).
# Mute masks reused from results/ab_royer/{NVI,SW}_mute/mute (AB_MUTE_DIR). Outputs results/ab_royer/{NVI,SW}_lag15/.
# First waits for the Arm B scrambled control (results/ARMB_NULL_DONE|FAILED; proceeds if no marker and no
# ARMB_NULL_PID after 2 h, or if the ARMB_NULL_PID process is gone). Markers results/LAG15_{NVI,SW}_{DONE,FAILED,RUNAWAY},
# results/LAG15_HEARTBEAT, pid results/LAG15_QUEUE_PID, lock locks/lag15_nvi_sw.lock, log logs/lag15_nvi_sw_queue.log.
# Stages are guarded: a relaunch skips finished regions/stages.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/lag15_nvi_sw.lock; flock -n 9 || exit 0
R=results
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/LAG15_HEARTBEAT; sleep 120; done ) &
HB=$!
GU=
trap 'kill $HB $GU 2>/dev/null' EXIT
echo "$$" > $R/LAG15_QUEUE_PID
echo "[queue] start $(date -u +%FT%TZ)"

# ---- wait for the Arm B scrambled control ----
T0=$(date +%s)
while true; do
  [ -f $R/ARMB_NULL_DONE ] && { echo "[queue] ARMB_NULL_DONE seen"; break; }
  [ -f $R/ARMB_NULL_FAILED ] && { echo "[queue] ARMB_NULL_FAILED seen"; break; }
  if [ -f $R/ARMB_NULL_PID ]; then
    P=$(cat $R/ARMB_NULL_PID 2>/dev/null)
    if [ -n "$P" ] && ! kill -0 "$P" 2>/dev/null; then
      sleep 60; [ -f $R/ARMB_NULL_DONE ] || [ -f $R/ARMB_NULL_FAILED ] || { echo "[queue] ARMB pid $P gone, no marker; proceeding"; break; }
      continue
    fi
  elif [ $(( $(date +%s) - T0 )) -gt 7200 ]; then
    echo "[queue] no ARMB marker/pid after 2 h; proceeding"; break
  fi
  sleep 300
done
echo "[queue] Arm B wait over $(date -u +%FT%TZ)"

run_region(){
  local S=$1; local M=results/ab_royer/${S}_mute; local B=results/ab_royer/${S}_lag15; local K=$R/LAG15_${S}
  [ -f ${K}_DONE ] && { echo "[$S] already done"; return 0; }
  rm -f ${K}_FAILED ${K}_RUNAWAY
  mkdir -p $B
  export AB_MUTE_DIR=$PWD/$M/mute
  echo "[$S] start $(date -u +%FT%TZ) AB_MUTE_DIR=$AB_MUTE_DIR"
  [ -f "$M/mute/params.json" ] || { echo "stage=mask_missing $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  [ -f "$M/auto_window/_params.json" ] || { echo "stage=auto_missing $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  [ -f "$B/clusters/summary.json" ] || \
    python scripts/ab_cluster.py --sub $S --auto "$M/auto_window" --out "$B/clusters" \
      --lag-max 15 --lag-mode zeropad \
      --headline-from results/ab_royer/NW/clusters/summary.json 9>&- || { echo "stage=cluster $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  echo "[$S] cluster done $(date -u +%FT%TZ)"
  if [ ! -f "$B/grow/GROW_DONE" ]; then
    python scripts/ab_grow_fast.py --sub $S --tpl "$B/clusters/tpl_it0.npz" --out "$B/grow" \
      --merge most --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 4 9>&- &
    local GP=$!; echo "$GP" > ${K}_GROW_PID
    ( exec 9>&-; while kill -0 $GP 2>/dev/null; do
        if [ -f "$B/grow/iters.csv" ] && python -c "import pandas as p,sys; d=p.read_csv('$B/grow/iters.csv'); sys.exit(0 if d.n_detections.max()>5e6 else 1)"; then
          echo "runaway $(date -u +%FT%TZ)" > ${K}_RUNAWAY; kill $(pgrep -P $GP) $GP 2>/dev/null; break; fi
        sleep 60; done ) &
    GU=$!
    wait $GP; local RC=$?
    kill $GU 2>/dev/null; GU=
    [ -f ${K}_RUNAWAY ] && { echo "stage=runaway $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
    [ $RC -eq 0 ] || { echo "stage=grow $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
    date -u +%FT%TZ > "$B/grow/GROW_DONE"
  fi
  local LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
  echo "[$S] grow done, last iteration $LAST $(date -u +%FT%TZ)"
  [ -f "$B/hypo/loc_headline.csv" ] || \
    python scripts/ab_hypo.py --sub $S --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 9>&- || true
  [ -f "$B/hypo/loc_headline.csv" ] || { echo "stage=hypo $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  [ -f "$B/hypo_28/stacks_loc.npz" ] || \
    python scripts/ab_restack28.py --sub $S --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/stacks_loc.npz" --workers 4 9>&- \
    || { echo "stage=restack28 $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  [ -f "$B/hypo_28/review/SHEETS_DONE" ] || { python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
    --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" 9>&- && date -u +%FT%TZ > "$B/hypo_28/review/SHEETS_DONE"; } \
    || { echo "stage=sheets $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  python scripts/ab_nomerge_score.py "$B/hypo_28/stacks_loc.npz" "$B/hypo_28/review/auto_score.csv" $S 9>&- \
    || { echo "stage=score $(date -u +%FT%TZ)" > ${K}_FAILED; return 1; }
  date -u +%FT%TZ > ${K}_DONE
  echo "[$S] done $(date -u +%FT%TZ)"
}

run_region NVI || echo "[queue] NVI failed; continuing to SW"
run_region SW  || echo "[queue] SW failed"
echo "[queue] end $(date -u +%FT%TZ)"
