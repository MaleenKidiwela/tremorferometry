#!/usr/bin/env bash
# NVI with BURST MUTING + line notching (2026-09-27, coordinator task: NVI after SW). Same chain and settings as
# scripts/nw_mute_queue.sh (K=15, qmax=2, qlong=1.3, pad 3 s). Masks built for the NVI 27-station pool + scan days first.
# Outputs results/ab_royer/NVI_mute/ only. Markers results/NVI_MUTE_{DONE,FAILED,HEARTBEAT}; log logs/nvi_mute_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nvi_mute.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/NVI_mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/NVI_MUTE_HEARTBEAT; sleep 120; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NVI_MUTE_FAILED; exit 1; }
rm -f $R/NVI_MUTE_FAILED $R/NVI_MUTE_DONE
echo "$$" > $R/NVI_MUTE_PID
echo "[queue] start $(date -u +%FT%TZ)"
# wait for scripts/ab_refetch_nvi.py (if still running), then the 40 sps cache for the NVI station-days re-fetched by scripts/ab_refetch_nvi.py (existing files are skipped)
while pgrep -f "scripts/ab_refetch_nvi.py" >/dev/null; do sleep 30; done
[ -f "$B/CACHE_REFETCH_DONE" ] || { mkdir -p "$B" && python scripts/ab_cache40.py --sub NVI --workers 16 \
  && date -u +%FT%TZ > "$B/CACHE_REFETCH_DONE"; } || fail cache40
# autocorrelation stations: sec.3.1 "different combinations of 7 of the 10" = per hour the first 7 of the
# listed 10 with 3-C data (ab_auto.py default; diagnosis in results/ab_royer/AB_ROYER_LOG.md, NVI 7-of-10)
for M in energy lines mask; do
  [ -f "$B/mute/${M}_DONE" ] || { python scripts/ab_mute.py $M --sub NVI --out "$B/mute" --K 15 --qmax 2 --qlong 1.3 --pad 3 --workers 16 \
    && date -u +%FT%TZ > "$B/mute/${M}_DONE"; } || fail mute_$M
done
[ -f "$B/mute/params.json" ] || fail mask_missing
export AB_MUTE_DIR=$PWD/$B/mute
[ -f "$B/auto_window/_params.json" ] || \
  python scripts/ab_auto.py --sub NVI --cols window --out "$B/auto_window" || fail auto
[ -f "$B/clusters/summary.json" ] || \
  python scripts/ab_cluster.py --sub NVI --auto "$B/auto_window" --out "$B/clusters" \
    --headline-from results/ab_royer/NW/clusters/summary.json || fail cluster
python scripts/ab_grow_fast.py --sub NVI --tpl "$B/clusters/tpl_it0.npz" --out "$B/grow" \
  --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 || fail grow
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
python scripts/ab_hypo.py --sub NVI --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub NVI --det "$B/grow/det_it${LAST}.parquet" \
  --out "$B/hypo_28/stacks_loc.npz" || fail restack28
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" || fail sheets
date -u +%FT%TZ > $R/NVI_MUTE_DONE
