#!/bin/bash
# Download-gated wave processor (CAUSALITY rollout, forward-only). For each station: wait for its
# download to finish -> discovery prep (candidate detect + adaptive-P) -> forward-only GPU phase
# (cluster/appB/select/resample/forward densify) -> CPU tail (fwd stacks + Z dv/v + causality finalize + <20 flag).
# GPU serialized via wait_gpu_free inside the GPU phase; the CPU tail runs in background so the next
# station's densify overlaps it. Usage: scripts/process_wave.sh [STA ...]  (default: full wave)
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log
WAVE="${*:-B003 B009 B010 B012 B024 B027 B031 B035 B943}"
echo "$(date +%H:%M) WAVE start (causality forward-only): $WAVE" >> "$STATUS"

dl_running(){ ps -eo args | grep -E "download_borehole.*--sta ${1}( |\$)" | grep -qv grep; }

for STA in $WAVE; do
  # wait until the download process for this station is gone
  waited=0
  until [ -d data/waveforms/PB.$STA ]; do sleep 300; waited=$((waited+5)); [ $waited -gt 1440 ] && break; done
  while dl_running "$STA"; do sleep 300; done
  yrs=$(ls -d data/waveforms/PB.$STA/20?? 2>/dev/null | wc -l)
  if [ ! -d data/waveforms/PB.$STA ] || [ "$yrs" -lt 12 ]; then
    echo "$(date +%H:%M) $STA WAVE skip: download absent/incomplete ($yrs yr dirs)" >> "$STATUS"; continue
  fi
  echo "$(date +%H:%M) === $STA WAVE begin (download complete, $yrs yr dirs) ===" >> "$STATUS"
  bash scripts/discovery_prep_station.sh "$STA" || { echo "$(date +%H:%M) $STA prep failed -> skip" >> "$STATUS"; continue; }
  bash scripts/borehole_gpu_phase_fwd.sh "$STA" || { echo "$(date +%H:%M) $STA gpu_phase failed -> skip" >> "$STATUS"; continue; }
  nohup bash scripts/borehole_cpu_tail.sh "$STA" > logs/tail_${STA}.log 2>&1 &
  echo "$(date +%H:%M) $STA -> CPU tail launched (background); next station's densify may overlap" >> "$STATUS"
done
wait
echo "$(date +%H:%M) WAVE complete" >> "$STATUS"
