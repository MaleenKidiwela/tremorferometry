#!/usr/bin/env bash
# 2026-09-28 (user-approved Merlin step 3): NW clustering with the FULL lag window, then the burst-muted chain with the
# ORIGINAL per-iteration merge (as NW_mute, which gave 66 real families; not fair, not no-merge).
#   ab_cluster.py --lag-max 15 --lag-mode zeropad (literal linear xcorr of the 15 s super traces; see ab_cluster.py header)
#   -> ab_grow_fast.py (--merge most = original; 14 iterations, tol 0.01, prep 6, restack 4 for the 200 GB pod)
#   -> ab_hypo (default GSC headline, --boot 60, as NW_mute) -> ab_restack28 -> review sheets
#   -> ab_nomerge_score.py site-count score (orig grouping + stations < 3 km apart).
# Same autocorrelation detections + mute masks as NW_mute. Outputs results/ab_royer/NW_lag15/;
# markers results/NW_LAG15_{DONE,FAILED,HEARTBEAT,RUNAWAY}; pids NW_LAG15_PID / NW_LAG15_GROW_PID;
# lock locks/nw_lag15.lock (held only by this shell); log logs/nw_lag15_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nw_lag15.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/NW_lag15; M=results/ab_royer/NW_mute
mkdir -p $B
export AB_MUTE_DIR=$PWD/$M/mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/NW_LAG15_HEARTBEAT; sleep 120; done ) &
HB=$!
GU=
trap 'kill $HB $GU 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NW_LAG15_FAILED; exit 1; }
rm -f $R/NW_LAG15_FAILED $R/NW_LAG15_DONE $R/NW_LAG15_RUNAWAY
echo "$$" > $R/NW_LAG15_PID
echo "[queue] start $(date -u +%FT%TZ) AB_MUTE_DIR=$AB_MUTE_DIR"
[ -f "$M/mute/params.json" ] || fail mask_missing
[ -f "$B/clusters/summary.json" ] || \
  python scripts/ab_cluster.py --sub NW --auto "$M/auto_window" --out "$B/clusters" \
    --lag-max 15 --lag-mode zeropad \
    --headline-from results/ab_royer/NW/clusters/summary.json 9>&- || fail cluster
echo "[queue] cluster done $(date -u +%FT%TZ)"
python scripts/ab_grow_fast.py --sub NW --tpl "$B/clusters/tpl_it0.npz" --out "$B/grow" \
  --merge most --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 4 9>&- &
GP=$!; echo "$GP" > $R/NW_LAG15_GROW_PID
# runaway guard: POST-merge detection total (iters.csv n_detections) > 5 M -> stop growth by exact pid
( exec 9>&-; while kill -0 $GP 2>/dev/null; do
    if [ -f "$B/grow/iters.csv" ] && python -c "import pandas as p,sys; d=p.read_csv('$B/grow/iters.csv'); sys.exit(0 if d.n_detections.max()>5e6 else 1)"; then
      echo "runaway $(date -u +%FT%TZ)" > $R/NW_LAG15_RUNAWAY; kill $(pgrep -P $GP) $GP 2>/dev/null; break; fi
    sleep 60; done ) &
GU=$!
wait $GP || { [ -f $R/NW_LAG15_RUNAWAY ] && fail runaway; fail grow; }
[ -f $R/NW_LAG15_RUNAWAY ] && fail runaway
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
echo "[queue] grow done, last iteration $LAST $(date -u +%FT%TZ)"
python scripts/ab_hypo.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 9>&- || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/stacks_loc.npz" --workers 4 9>&- || fail restack28
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" 9>&- || fail sheets
python scripts/ab_nomerge_score.py "$B/hypo_28/stacks_loc.npz" "$B/hypo_28/review/auto_score.csv" NW 9>&- || fail score
date -u +%FT%TZ > $R/NW_LAG15_DONE
echo "[queue] done $(date -u +%FT%TZ)"
