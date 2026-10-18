#!/usr/bin/env bash
# P4-0 (Puget / Olympic site-noise + registration-proxy inventory) - ONE queue, start to finish.
#
#   nohup bash scripts/u_p4_0_queue.sh > logs/p4_probe_queue.log 2>&1 &
#
# Steps (each resumable; re-running never re-fetches a day that is on disk):
#   1 u_select.py     candidates, IRIS StationXML (response), the 4+4 probe days, the queue
#   2 u_download.py   the queue at <= 8 concurrent FDSN requests -> traces/
#   3 u_cache.py      frozen 20 Hz / 1-8 Hz chain (rational variant for 50/125 Hz) -> cache20hz/
#   4 u_stats.py      per station-day velocity PSD + incoherent 6-s energy
#   5 u_inventory.py  noise_inventory / shortlist_v0 / geometry / figure
#
# Finishes by writing results/p4/P4_0_DONE (success) or results/p4/P4_0_FAILED (traceback).
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
P4="$ROOT/results/p4"
LOG="$ROOT/logs/p4_probe_queue.log"
DONE="$P4/P4_0_DONE"
FAILED="$P4/P4_0_FAILED"
mkdir -p "$P4" "$ROOT/logs" "$ROOT/figures" "$ROOT/traces" "$ROOT/cache20hz"
rm -f "$DONE" "$FAILED"
T0=$(date -u +%s)

step () {                       # step <name> <command...>
  local name="$1"; shift
  echo "=== [$(date -u +%FT%TZ)] STEP $name: $* ==="
  if ! "$@"; then
    local rc=$?
    {
      echo "P4-0 FAILED in step $name (exit $rc) at $(date -u +%FT%TZ)"
      echo "command: $*"
      echo "--- last 120 lines of $LOG ---"
      tail -120 "$LOG" 2>/dev/null
    } > "$FAILED"
    echo "!!! FAILED in $name (exit $rc) -> $FAILED"
    exit "$rc"
  fi
  echo "=== [$(date -u +%FT%TZ)] STEP $name OK ==="
}

step select    python -u scripts/u_select.py
step download  python -u scripts/u_download.py --workers 8
step cache     python -u scripts/u_cache.py --workers 24
step stats     python -u scripts/u_stats.py --workers 24
step inventory python -u scripts/u_inventory.py

{
  echo "P4-0 DONE at $(date -u +%FT%TZ)  (elapsed $(( ($(date -u +%s) - T0) / 60 )) min)"
  echo "traces:    $(du -sh "$ROOT/traces" 2>/dev/null | cut -f1)"
  echo "cache20hz: $(du -sh "$ROOT/cache20hz" 2>/dev/null | cut -f1)"
  for f in probe_days.csv candidates.csv station_probe_days.csv noise_inventory.csv \
           shortlist_v0.csv station_rank_all.csv shortlist_cell_geometry.csv \
           shortlist_geometry.json download_manifest.csv station_day_stats.csv; do
    [ -f "$P4/$f" ] && echo "  $f  $(wc -l < "$P4/$f") lines"
  done
  echo "figure: $ROOT/figures/u_noise_inventory.png"
  echo "--- inventory summary ---"
  sed -n '1,40p' "$P4/u_inventory.log" 2>/dev/null
} > "$DONE"
cat "$DONE"
