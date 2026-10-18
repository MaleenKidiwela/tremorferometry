#!/bin/bash
# Thomas & Bostock (2015) chain AFTER the cross-station detections exist (results/xs/*_det.parquet).
# Launch-and-exit: nohup bash scripts/xs_chain_queue.sh <TAG> [extra xs_grow args] &
# Markers: results/<TAG>/CHAIN_DONE or results/<TAG>/CHAIN_FAILED (+ the failing step name).
# GPU: xs_grow / xs_scan46 use ab_grow.matched_filter on cuda -- do NOT start while another GPU
# run owns the card unless the coordinator says so.
set -u
cd "$(dirname "$0")/.."
TAG=${1:-chain_v1}; shift || true
R=results/$TAG; mkdir -p $R logs
export OMP_NUM_THREADS=4
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
step() { local name=$1; shift; echo "[$(date +%T)] START $name" >> $R/queue.log
  if ! "$@" >> logs/${TAG}_$name.log 2>&1; then echo "$name" > $R/CHAIN_FAILED; echo "[$(date +%T)] FAIL $name" >> $R/queue.log; exit 1; fi
  echo "[$(date +%T)] DONE $name" >> $R/queue.log; }
echo $$ > $R/QUEUE_PID
step grow    python scripts/xs_grow.py --out $R/grow "$@"
step scan46  python scripts/xs_scan46.py --grow $R/grow --out $R/scan46
step locate  python scripts/xs_locate.py --grow $R/grow --scan $R/scan46 --out $R/hypo
step review  python scripts/xs_review.py --mode stacks --stacks $R/hypo/stacks_loc.npz --picks $R/hypo/picks.csv --out $R/review
date +%FT%T > $R/CHAIN_DONE
