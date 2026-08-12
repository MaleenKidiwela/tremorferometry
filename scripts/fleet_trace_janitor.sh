#!/bin/bash
# Trace-janitor for the broadband fleet (rolling buffer). While fleet_broadband_all is alive, every 5 min free
# the RAW TRACES of every COMPLETED station (has "STATION DONE" in its log) EXCEPT the currently-processing one.
# Deletes ONLY data/waveforms/<NET>.<STA>/. NEVER touches logs, dv/v csvs, long_window stacks, cert/summary.
# Idempotent; safe to relaunch after a crash. Usage: nohup bash scripts/fleet_trace_janitor.sh &
cd /home/jovyan/tremorferometry
while ps -eo cmd | grep -q "[f]leet_broadband_all"; do
  cur=$(ps -eo cmd | grep -oE "fleet_station.sh [A-Z]+ [A-Z0-9]+" | awk "{print \$3}")   # currently processing (keep)
  for dpath in data/waveforms/*/; do
    [ -d "$dpath" ] || continue
    ns=$(basename "$dpath"); sta=${ns#*.}; s=$(echo "$sta" | tr A-Z a-z)
    [ "$sta" = "$cur" ] && continue
    if grep -q "STATION DONE: $sta" logs/fleet_${s}.log 2>/dev/null; then
      sz=$(du -sm "$dpath" 2>/dev/null | cut -f1)
      rm -rf "$dpath" && echo "$(date +%H:%M) janitor freed $ns (${sz}MB)" >> logs/trace_cleanup.log
    fi
  done
  sleep 300
done
echo "$(date +%H:%M) janitor exit (batch done)" >> logs/trace_cleanup.log
