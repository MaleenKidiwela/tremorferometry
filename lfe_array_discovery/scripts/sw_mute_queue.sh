#!/usr/bin/env bash
# SW with BURST MUTING + line notching (2026-09-27, user-approved "ok SW then"). Same chain and settings as
# scripts/nw_mute_queue.sh (K=15, qmax=2, qlong=1.3, pad 3 s). Masks built for the SW pool + scan days first.
# Outputs results/ab_royer/SW_mute/ only. Markers results/SW_MUTE_{DONE,FAILED,HEARTBEAT}; log logs/sw_mute_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/sw_mute.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/SW_mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/SW_MUTE_HEARTBEAT; sleep 120; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/SW_MUTE_FAILED; exit 1; }
rm -f $R/SW_MUTE_FAILED $R/SW_MUTE_DONE
echo "$$" > $R/SW_MUTE_PID
echo "[queue] start $(date -u +%FT%TZ)"
for M in energy lines mask; do
  [ -f "$B/mute/${M}_DONE" ] || { python scripts/ab_mute.py $M --sub SW --out "$B/mute" --K 15 --qmax 2 --qlong 1.3 --pad 3 --workers 16 \
    && date -u +%FT%TZ > "$B/mute/${M}_DONE"; } || fail mute_$M
done
[ -f "$B/mute/params.json" ] || fail mask_missing
export AB_MUTE_DIR=$PWD/$B/mute
[ -f "$B/auto_window/_params.json" ] || \
  python scripts/ab_auto.py --sub SW --cols window --out "$B/auto_window" || fail auto
[ -f "$B/clusters/summary.json" ] || \
  python scripts/ab_cluster.py --sub SW --auto "$B/auto_window" --out "$B/clusters" \
    --headline-from results/ab_royer/NW/clusters/summary.json || fail cluster
python scripts/ab_grow_fast.py --sub SW --tpl "$B/clusters/tpl_it0.npz" --out "$B/grow" \
  --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 || fail grow
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
python scripts/ab_hypo.py --sub SW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub SW --det "$B/grow/det_it${LAST}.parquet" \
  --out "$B/hypo_28/stacks_loc.npz" || fail restack28
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" || fail sheets
date -u +%FT%TZ > $R/SW_MUTE_DONE
