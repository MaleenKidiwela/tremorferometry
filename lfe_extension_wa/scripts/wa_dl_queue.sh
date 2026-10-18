#!/usr/bin/env bash
# WA production downloads, tiers A + B (results/WA_PROD_PLAN.md §4; user-approved 2026-10-02). Tiers C/D NOT approved.
#   nohup setsid -f bash scripts/wa_dl_queue.sh > logs/wa_dl_queue.out 2>&1 < /dev/null
#   results/wa_dl/PID        queue pid           results/wa_dl/HEARTBEAT  "<utc> <stage> <progress>" every 60 s
#   results/wa_dl/DONE / results/wa_dl/FAILED (stage + log tail); results/wa_dl/STOPPED_STORAGE if the 0.8 TB cap hit
#   resumable: re-launch skips every station-day already in results/wa_dl/manifest.csv as data/nodata/linked
#   stop: kill $(cat results/wa_dl/PID); then kill the python child (ps --ppid <pid>)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
O=$S/results/wa_dl; L=$S/logs/wa_dl.log
mkdir -p "$O" "$S/logs"
exec 9>"$O/.lock"
if ! flock -n 9; then echo "another wa_dl_queue holds the lock"; exit 1; fi
echo $$ > "$O/PID"
rm -f "$O/DONE" "$O/FAILED"
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2
( while true; do echo "$(date -u +%FT%TZ) $(cat $O/.stage 2>/dev/null) $(cat $O/PROGRESS 2>/dev/null)" > "$O/HEARTBEAT"; sleep 60; done ) 9>&- &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT

run() {  # run <name> <cmd...>
  local name=$1; shift
  echo "$name" > "$O/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" >> "$L"
  "$@" >> "$S/logs/wa_dl_$name.out" 2>&1   # python writes its own log lines to $L
  local rc=$?
  if [ $rc -ne 0 ]; then
    { echo "stage $name failed rc=$rc $(date -u +%FT%TZ)"; cat "$O/STOPPED_STORAGE" 2>/dev/null; tail -40 "$S/logs/wa_dl_$name.out"; } > "$O/FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $name rc=$rc" >> "$L"
    exit 1
  fi
  echo "$(date -u +%FT%TZ) [ok] $name" >> "$L"
}

df -h /home/jovyan | tail -1 >> "$L"
run fetch   python scripts/wa_dl.py run --workers 4 --rate 2 --cap-tb 1.0
run retry1  python scripts/wa_dl.py run --workers 4 --rate 2 --cap-tb 1.0 --retry-errors-only
run retry2  python scripts/wa_dl.py run --workers 2 --rate 1 --cap-tb 1.0 --retry-errors-only
run summary python scripts/wa_dl.py summary
df -h /home/jovyan | tail -1 >> "$L"
echo "done" > "$O/.stage"
date -u +%FT%TZ > "$O/DONE"
echo "$(date -u +%FT%TZ) [DONE]" >> "$L"
