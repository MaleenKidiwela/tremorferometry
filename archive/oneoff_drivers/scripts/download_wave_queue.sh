#!/bin/bash
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
log(){ echo "$(date +%H:%M) [dlq] $*" | tee -a logs/rollout_status.log; }
MAXCONC=3
for S in B012 B024 B027 B031 B035 B943; do
  # throttle to MAXCONC concurrent downloads
  while [ "$(ps -eo cmd | grep -c '[d]ownload_borehole_3comp')" -ge "$MAXCONC" ]; do sleep 60; done
  [ "$(ls data/waveforms/PB.$S/*/*.mseed 2>/dev/null | wc -l)" -gt 6000 ] && { log "$S already downloaded, skip"; continue; }
  log "$S download START"
  nohup $PY scripts/download_borehole_3comp.py --net PB --sta $S --start 2007-01-01 --end 2026-12-31 --workers 6 >> logs/download_${S,,}.log 2>&1 &
  sleep 5
done
log "all wave downloads launched/queued (B012 B024 B027 B031 B035 B943)"
