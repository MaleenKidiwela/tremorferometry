#!/usr/bin/env bash
# 2026-09-28 (user-approved): NW growth with the FAIR per-iteration merge (ab_grow_fast.py --merge fair: +-15 s modal lag
# with +-60 s shift control, excess > 10 %; representative = higher 2-8 Hz site-count score, tie more detections; star;
# never merge across the 3-site line) from NW_mute it0 cluster stacks, 14 iterations, tol 0.01; then ab_hypo,
# ab_restack28, review sheets, site-count score -> NW_fair/hypo_28/review/auto_score.csv.
# Outputs results/ab_royer/NW_fair/; markers results/NW_FAIR_{DONE,FAILED,HEARTBEAT,RUNAWAY}; pid results/NW_FAIR_PID;
# grow pid results/NW_FAIR_GROW_PID; lock locks/nw_fair.lock (held only by this shell); log logs/nw_fair_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nw_fair.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/NW_fair
mkdir -p $B
export AB_MUTE_DIR=$PWD/results/ab_royer/NW_mute/mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/NW_FAIR_HEARTBEAT; sleep 120; done ) &
HB=$!
GU=
trap 'kill $HB $GU 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NW_FAIR_FAILED; exit 1; }
rm -f $R/NW_FAIR_FAILED $R/NW_FAIR_DONE $R/NW_FAIR_RUNAWAY
echo "$$" > $R/NW_FAIR_PID
echo "[queue] start $(date -u +%FT%TZ) AB_MUTE_DIR=$AB_MUTE_DIR"
python scripts/ab_grow_fast.py --sub NW --tpl results/ab_royer/NW_mute/clusters/tpl_it0.npz --out "$B/grow" \
  --merge fair --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 9>&- &
GP=$!; echo "$GP" > $R/NW_FAIR_GROW_PID
# runaway guard: stop the growth (exact pid) if any iteration's detection total AFTER the merge exceeds 5 M (raw pre-merge counts include the duplicates the merge removes)
( exec 9>&-; while kill -0 $GP 2>/dev/null; do
    if [ -f "$B/grow/iters.csv" ] && python -c "import pandas as p,sys; d=p.read_csv('$B/grow/iters.csv'); sys.exit(0 if d.n_detections.max()>5e6 else 1)"; then
      echo "runaway $(date -u +%FT%TZ)" > $R/NW_FAIR_RUNAWAY; kill $(pgrep -P $GP) $GP 2>/dev/null; break; fi
    sleep 60; done ) &
GU=$!
wait $GP || { [ -f $R/NW_FAIR_RUNAWAY ] && fail runaway; fail grow; }
[ -f $R/NW_FAIR_RUNAWAY ] && fail runaway
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
echo "[queue] grow done, last iteration $LAST $(date -u +%FT%TZ)"
python scripts/ab_hypo.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/stacks_loc.npz" || fail restack28
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" || fail sheets
python scripts/ab_nomerge_score.py "$B/hypo_28/stacks_loc.npz" "$B/hypo_28/review/auto_score.csv" NW || fail score
date -u +%FT%TZ > $R/NW_FAIR_DONE
echo "[queue] done $(date -u +%FT%TZ)"
