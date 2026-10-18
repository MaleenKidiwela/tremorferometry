#!/usr/bin/env bash
# Track B runner - launch the resumable download worker in the background.
#   bash scripts/download_run.sh              # everything, P0 -> P3
#   bash scripts/download_run.sh P0,P1        # only those priorities
# Safe to re-run: a second instance refuses to start while one is alive; a restart after a
# crash resumes from results/download_manifest.csv and never re-fetches a present day.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRIOS="${1:-P0,P1,P2,P2b,P3}"
WORKERS="${2:-8}"
PY=/opt/conda/bin/python
mkdir -p "$HERE/logs" "$HERE/results" "$HERE/traces" "$HERE/figures"
PIDF="$HERE/logs/download_run.pid"
if [ -f "$PIDF" ] && kill -0 "$(cat "$PIDF")" 2>/dev/null; then
  echo "already running as PID $(cat "$PIDF") - tail $HERE/logs/download_run.log"; exit 1
fi
cd "$HERE"
nohup "$PY" scripts/download_worker.py --workers "$WORKERS" --priorities "$PRIOS" \
  >> "$HERE/logs/download_run.log" 2>&1 &
echo $! > "$PIDF"
echo "started PID $(cat "$PIDF")  prios=$PRIOS workers=$WORKERS"
echo "  log:      tail -f $HERE/logs/download_run.log"
echo "  progress: cat $HERE/results/download_progress.json"
echo "  stop:     kill \$(cat $PIDF)   # resumable"
