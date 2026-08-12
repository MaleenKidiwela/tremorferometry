#!/bin/bash
# Rolling-buffer wave (user policy 2026-07-11): downloads are HELD; process each ON-DISK station fully,
# then DELETE its traces before moving on. Adds the BATTERY GATE after dv/v. Serial (no cross-station
# overlap) so peak memory stays low (OOM-safe) and traces free immediately. Fully guarded/idempotent:
# waits for any in-flight same-station work (crash orphans), skips completed steps.
# Per station: [wait-idle] -> discovery_prep -> gpu_phase_fwd -> cpu_tail(stacks+dvv+finalize)
#              -> battery_gate -> DELETE TRACES.  Usage: scripts/rolling_wave.sh <STA...>
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) [roll] $*" | tee -a "$STATUS"; }
QUEUE="${*:?usage: rolling_wave.sh STA [STA...]}"
log "ROLLING WAVE start (downloads HELD; battery+delete tail): $QUEUE"

# wait while any pipeline process for station $1 is still running (handles crash-orphaned in-flight work)
wait_station_idle(){
  local S=$1 sl; sl=$(echo "$S" | tr '[:upper:]' '[:lower:]'); local waited=0
  while ps -eo cmd | grep -vE 'grep|rolling_wave' \
        | grep -qiE "[ /=]$S( |/|\$)|PB[.]$S|${sl}p90f40|${sl}_(cand|disc|pnsn)"; do
    [ $waited -eq 0 ] && log "$S waiting for in-flight work to finish"
    sleep 60; waited=$((waited+1))
  done
}

for STA in $QUEUE; do
  sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]'); TAG=${STA}p90f40
  log "=== $STA begin ==="
  wait_station_idle "$STA"

  bash scripts/discovery_prep_station.sh "$STA"  || { log "$STA prep FAILED -> skip";      continue; }
  bash scripts/borehole_gpu_phase_fwd.sh "$STA"  || { log "$STA gpu_phase FAILED -> skip";  continue; }
  bash scripts/borehole_cpu_tail.sh "$STA"       || { log "$STA cpu_tail FAILED -> skip";   continue; }
  # battery gate SKIPPED (user 2026-07-11): battery can't separate real/ringing on cap-off forward-only
  # data (ringing coda_sigma >= certified; no threshold separates) -> causality-only is the product.

  # DELETE TRACES (rolling buffer) — only after stacks + causality finalize verified present
  if [ -f data/long_window_daily_${TAG}_Z.npz ] && [ -f data/${sta}_3comp_summary.json ]; then
    sz=$(du -sh data/waveforms/PB.$STA 2>/dev/null | cut -f1)
    rm -rf data/waveforms/PB.$STA && log "$STA traces DELETED ($sz freed; stacks+finalize verified)"
  else
    log "$STA traces KEPT (safety): stacks or finalize missing"
  fi
  log "=== $STA COMPLETE ==="
done
log "ROLLING WAVE complete: $QUEUE"
