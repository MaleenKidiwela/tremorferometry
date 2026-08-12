#!/bin/bash
# Vigilant supervisor for the continue-through-29 run (2026-07-13). Keeps the driver + guards alive across
# crashes (all resumable) and flags stalls. Loops every 5 min until 29 stations are done.
#   - restarts continue_29.sh if it died and <29 done (resumes via markers + guarded steps)
#   - restarts oom_guard / map_refresh_daemon if they died
#   - flags a stall: a single densify/build_long_window running >4 h (likely hung)
#   - logs a heartbeat with count + active station
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) [sup] $*" | tee -a "$STATUS"; }
ndone(){ ls data/*_3comp_summary.json 2>/dev/null | wc -l; }
alive(){ ps -eo cmd | grep -q "[${1:0:1}]${1:1}"; }
log "supervisor start ($(ndone)/29 done)"
while [ "$(ndone)" -lt 29 ]; do
  # keep guards alive
  alive oom_guard.sh          || { nohup bash scripts/oom_guard.sh >> logs/oom_watch.log 2>&1 & log "restarted oom_guard"; }
  alive map_refresh_daemon.sh || { nohup bash scripts/map_refresh_daemon.sh >> logs/map_refresh_daemon.log 2>&1 & log "restarted map_daemon"; }
  # keep the driver alive (resumable)
  if ! alive continue_29.sh; then
    nohup bash scripts/continue_29.sh >> logs/continue_29.log 2>&1 &
    log "RESTARTED continue_29 (was dead; resuming, $(ndone)/29 done)"
  fi
  # stall detector: densify/stack running > 4 h
  ps -eo etimes,cmd | grep -E '[d]ensify_gnw|[b]uild_long_window' | while read -r et cmd; do
    if [ "${et:-0}" -gt 14400 ]; then
      st=$(echo "$cmd" | grep -oE 'station [A-Z0-9]+' | head -1)
      log "STALL? $st running $((et/3600))h (>4h) — check for hang"
    fi
  done
  act=$(ps -eo cmd | grep -oE 'station [A-Z0-9]+' | grep -v grep | head -1)
  log "heartbeat: $(ndone)/29 | active: ${act:-idle}"
  sleep 300
done
log "SUPERVISOR DONE — 29/29 stations complete"
