#!/bin/bash
# Autonomous "continue through 29" driver (user 2026-07-12). Downloads lifted from hold. Brings the
# borehole fleet to 29 completed dv/v stations. Two parts:
#   FEEDER (bg): downloads queue stations in order, throttled to keep <=4 unprocessed trace-sets buffered.
#   PROCESSOR (fg, serial): per station wait-for-download -> prep -> gpu_phase -> cpu_tail (battery-free)
#                           -> delete traces. Gap audit + map refresh happen via map_refresh_daemon.
# Stops when 29 stations have a *_3comp_summary.json. Resumable (dl markers + guarded steps + skip-if-done).
# Usage: nohup bash scripts/continue_29.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) [29] $*" | tee -a "$STATUS"; }
# B031/B943 (finish partials) then prepped boreholes ranked by GOLD count; extra tail for failure headroom.
# Endgame redirect (Merlin 2026-07-13): refuse GOLD-1/2 padding (B935/B204/B028: duplicates/no new interface
# coverage); finish B941+B932 and RECOVER the high-quality B045(GOLD52)+B030(GOLD23) instead. Floor now 0.2.
QUEUE="B017 B935"
TARGET=99
ndone(){ ls data/*_3comp_summary.json 2>/dev/null | wc -l; }
log "CONTINUE-THROUGH-29 start ($(ndone) done). queue: $QUEUE"

# count queue stations with traces on disk but not yet finalized (buffer depth)
buffered(){ local n=0 S sl; for S in $QUEUE; do sl=${S,,}; [ -d data/waveforms/PB.$S ] && [ ! -f data/${sl}_3comp_summary.json ] && n=$((n+1)); done; echo $n; }

feeder(){
  for S in $QUEUE; do
    [ "$(ndone)" -ge "$TARGET" ] && { log "feeder: target reached, stop downloading"; break; }
    sl=${S,,}; [ -f data/${sl}_3comp_summary.json ] && continue      # already done
    [ -f data/.dl_done_$S ] && continue                              # already downloaded
    while [ "$(buffered)" -ge 4 ]; do sleep 300; done                # throttle: <=4 buffered
    log "download $S start"
    $PY scripts/download_borehole_3comp.py --net PB --sta $S --start 2007-01-01 --end 2026-12-31 --workers 6 \
       >> logs/download_${sl}.log 2>&1
    nd=$(ls data/waveforms/PB.$S/*/*.mseed 2>/dev/null | wc -l)
    touch data/.dl_done_$S
    [ "${nd:-0}" -gt 3000 ] && log "download $S done ($nd days)" || log "download $S WEAK (${nd:-0} days) — proceeding"
  done
  log "feeder: all dispatched"
}
feeder & FEEDER=$!

wait_station_idle(){ local S=$1 sl=${1,,}; while ps -eo cmd | grep -vE 'grep|continue_29|rolling_wave' \
  | grep -qiE "[ /=]$S( |/|\$)|PB[.]$S|${sl}p90f40|${sl}_(cand|disc|pnsn)"; do sleep 60; done; }

for STA in $QUEUE; do
  sta=${STA,,}; TAG=${STA}p90f40
  [ "$(ndone)" -ge "$TARGET" ] && { log "reached $TARGET stations -> STOP"; break; }
  [ -f data/${sta}_3comp_summary.json ] && { log "$STA already done -> skip"; continue; }
  log "=== $STA begin (await download) ==="
  until [ -f data/.dl_done_$STA ]; do sleep 120; done
  yrs=$(ls -d data/waveforms/PB.$STA/20?? 2>/dev/null | wc -l)
  [ "$yrs" -lt 8 ] && { log "$STA SKIP: too few year-dirs ($yrs) — no usable data"; continue; }
  wait_station_idle "$STA"
  bash scripts/discovery_prep_station.sh "$STA" || { log "$STA prep FAILED -> skip"; continue; }
  bash scripts/borehole_gpu_phase_fwd.sh "$STA" || { log "$STA gpu_phase FAILED -> skip"; continue; }
  bash scripts/borehole_cpu_tail.sh "$STA"      || { log "$STA cpu_tail FAILED -> skip"; continue; }
  if [ -f data/.keep_traces_$STA ]; then
    log "$STA traces KEPT (.keep_traces marker — needed for catalog-template augmentation)"
  elif [ -f data/long_window_daily_${TAG}_Z.npz ] && [ -f data/${sta}_3comp_summary.json ]; then
    sz=$(du -sh data/waveforms/PB.$STA 2>/dev/null | cut -f1); rm -rf data/waveforms/PB.$STA
    log "$STA traces DELETED ($sz freed)"
  fi
  zc=$($PY -c "import json;print(json.load(open('data/${sta}_3comp_summary.json'))['z_certified'])" 2>/dev/null)
  log "=== $STA COMPLETE (Z ${zc} fam; $(ndone) stations done) ==="
done
kill $FEEDER 2>/dev/null
log "CONTINUE-THROUGH-29 finished: $(ndone) stations done"
