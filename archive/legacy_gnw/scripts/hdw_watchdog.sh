#!/bin/bash
# Standalone anon-memory watchdog for an already-running HDW discovery. Kills it if real (anon)
# memory exceeds 130 GB (cap 187 GB), ignoring reclaimable page cache. Exits when discovery ends.
L=logs/hdw_watchdog.log
while pgrep -f "[d]iscover_nllb_pnsn" >/dev/null; do
  a=$(awk '/^anon /{print $2}' /sys/fs/cgroup/memory.stat)
  echo "$(date '+%H:%M:%S') anon $((a/1000000000))G" >> "$L"
  if [ "$a" -gt 130000000000 ]; then
    echo "$(date '+%H:%M:%S') WATCHDOG KILL anon $((a/1000000000))G" >> "$L"
    pkill -9 -f "[m]ultiprocessing.spawn"; pkill -9 -f "[d]iscover_nllb_pnsn"
    break
  fi
  sleep 20
done
echo "$(date '+%H:%M:%S') watchdog exit (discovery ended)" >> "$L"
