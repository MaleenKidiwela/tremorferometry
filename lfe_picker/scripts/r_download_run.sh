#!/usr/bin/env bash
# P3D download runner.  Resumable: re-running never re-fetches a day that is on disk.
#   bash scripts/r_download_run.sh P0            # SSE windows first
#   bash scripts/r_download_run.sh P1,P2 6       # the fill, 6 workers
# Logs -> logs/p3d_download_<prios>_<stamp>.log ; progress -> results/p3d/download_progress.json
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRIOS="${1:-P0}"
WORKERS="${2:-8}"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
TAG="$(echo "$PRIOS" | tr ',' '-')"
LOG="$ROOT/logs/p3d_download_${TAG}_${STAMP}.log"
mkdir -p "$ROOT/logs" "$ROOT/results/p3d" "$ROOT/traces"
cd "$ROOT"
echo "[run] prios=$PRIOS workers=$WORKERS log=$LOG"
nohup python -u scripts/r_download_worker.py --workers "$WORKERS" --priorities "$PRIOS" \
      >> "$LOG" 2>&1 &
echo $! > "$ROOT/logs/p3d_download_${TAG}.pid"
echo "[run] pid $(cat "$ROOT/logs/p3d_download_${TAG}.pid")"
