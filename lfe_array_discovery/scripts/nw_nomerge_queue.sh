#!/usr/bin/env bash
# Merlin 09-28 step 2 (user-approved): full NW no-merge growth from NW_mute it0 cluster stacks, 14 iterations,
# then 2-8 Hz restack + site-count score of all families -> ONE final dedupe (ab_final_dedupe.py) ->
# ab_hypo + review sheets + scores on the representatives. Outputs results/ab_royer/NW_nomerge/full/.
# Markers results/NW_NOMERGE_{DONE,FAILED,HEARTBEAT}; pid results/NW_NOMERGE_PID; log logs/nw_nomerge_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nw_nomerge.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/NW_nomerge/full
mkdir -p $B
export AB_MUTE_DIR=$PWD/results/ab_royer/NW_mute/mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/NW_NOMERGE_HEARTBEAT; sleep 120; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NW_NOMERGE_FAILED; exit 1; }
rm -f $R/NW_NOMERGE_FAILED $R/NW_NOMERGE_DONE $R/NW_NOMERGE_RUNAWAY
echo "$$" > $R/NW_NOMERGE_PID
echo "[queue] start $(date -u +%FT%TZ) AB_MUTE_DIR=$AB_MUTE_DIR"
python scripts/ab_grow_fast.py --sub NW --tpl results/ab_royer/NW_mute/clusters/tpl_it0.npz --out "$B/grow" \
  --no-merge --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 9>&- &
GP=$!; echo "$GP" > $R/NW_NOMERGE_GROW_PID
# runaway guard: stop the growth (exact pid) if any iteration's detection total exceeds 5 M
( exec 9>&-; while kill -0 $GP 2>/dev/null; do
    if [ -f "$B/grow/iters.csv" ] && python -c "import pandas as p,sys; d=p.read_csv('$B/grow/iters.csv'); sys.exit(0 if d.n_det_raw.max()>5e6 else 1)"; then
      echo "runaway $(date -u +%FT%TZ)" > $R/NW_NOMERGE_RUNAWAY; kill $GP; break; fi
    sleep 60; done ) &
GU=$!
trap 'kill $HB $GU 2>/dev/null' EXIT
wait $GP || fail grow
[ -f $R/NW_NOMERGE_RUNAWAY ] && fail runaway
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
python scripts/ab_restack28.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/stacks28_all.npz" || fail restack28
python scripts/ab_nomerge_score.py "$B/stacks28_all.npz" "$B/score_all.csv" NW || fail score
python scripts/ab_final_dedupe.py --det "$B/grow/det_it${LAST}.parquet" --score "$B/score_all.csv" \
  --score-col n_sites_3km --out "$B/dedupe" || fail dedupe
python scripts/ab_hypo.py --sub NW --det "$B/dedupe/det_dedup.parquet" --out "$B/hypo" --boot 60 || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub NW --det "$B/dedupe/det_dedup.parquet" --out "$B/hypo_28/stacks_loc.npz" || fail restack28_dedup
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/dedupe/det_dedup.parquet" --out "$B/hypo_28/review" || fail sheets
python scripts/ab_nomerge_score.py "$B/hypo_28/stacks_loc.npz" "$B/hypo_28/score_dedup.csv" NW || fail score_dedup
date -u +%FT%TZ > $R/NW_NOMERGE_DONE
