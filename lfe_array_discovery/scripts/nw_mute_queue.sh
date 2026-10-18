#!/usr/bin/env bash
# NW rerun with BURST MUTING + line notching (2026-09-26, user-approved). Mute masks from scripts/ab_mute.py
# (K=15, qmax=2, qlong=1.3, pad 3 s; results/ab_royer/NW_mute/mute/), consumed by ab_common.load_day via
# AB_MUTE_DIR (muted samples = gap samples). Chain from the autocorrelation stage, same settings as run 1;
# full 25-station pool (W060/FACB back in). Outputs results/ab_royer/NW_mute/ only.
# Markers results/NW_MUTE_{DONE,FAILED,HEARTBEAT}; log logs/nw_mute_queue.log.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nw_mute.lock; flock -n 9 || exit 0
R=results; B=results/ab_royer/NW_mute
export AB_MUTE_DIR=$PWD/$B/mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
( exec 9>&-; while true; do date -u +%FT%TZ > $R/NW_MUTE_HEARTBEAT; sleep 120; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NW_MUTE_FAILED; exit 1; }
rm -f $R/NW_MUTE_FAILED $R/NW_MUTE_DONE
echo "$$" > $R/NW_MUTE_PID
[ -f "$B/mute/params.json" ] || fail mask_missing
echo "[queue] start $(date -u +%FT%TZ) AB_MUTE_DIR=$AB_MUTE_DIR"
[ -f "$B/auto_window/_params.json" ] || \
  python scripts/ab_auto.py --sub NW --cols window --out "$B/auto_window" || fail auto
[ -f "$B/clusters/summary.json" ] || \
  python scripts/ab_cluster.py --sub NW --auto "$B/auto_window" --out "$B/clusters" \
    --headline-from results/ab_royer/NW/clusters/summary.json || fail cluster
python scripts/ab_grow_fast.py --sub NW --tpl "$B/clusters/tpl_it0.npz" --out "$B/grow" \
  --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 || fail grow
LAST=$(ls $B/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
python scripts/ab_hypo.py --sub NW --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo" --boot 60 || true
[ -f "$B/hypo/loc_headline.csv" ] || fail hypo
python scripts/ab_restack28.py --sub NW --det "$B/grow/det_it${LAST}.parquet" \
  --out "$B/hypo_28/stacks_loc.npz" || fail restack28
python scripts/ab_review_sheets.py --stacks "$B/hypo_28/stacks_loc.npz" \
  --det "$B/grow/det_it${LAST}.parquet" --out "$B/hypo_28/review" || fail sheets
date -u +%FT%TZ > $R/NW_MUTE_DONE
