#!/bin/bash
# Auto-refresh the borehole dv/v map whenever a station finalizes. DECOUPLED from rolling_wave (safe to
# run alongside it and across waves / manual runs) — it watches the set of per-station causality summaries
# and rebuilds figures/borehole_dvv_map.html on any change. Usage: nohup bash scripts/map_refresh_daemon.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
sig(){ ls -la --time-style=+%s data/*_3comp_summary.json 2>/dev/null | awk '{print $6,$NF}' | md5sum; }
last=""
echo "$(date +%H:%M:%S) [map] daemon start (watching *_3comp_summary.json)"
while true; do
  cur=$(sig)
  if [ "$cur" != "$last" ]; then
    n=$(ls data/daily_dvv_*p90f40_Z_2to4.csv 2>/dev/null | wc -l)
    echo "$(date +%H:%M:%S) [map] change detected -> rebuilding ($n stations)"
    if $PY scripts/build_borehole_dvv_map.py >> logs/map_refresh.log 2>&1; then
      echo "$(date +%H:%M:%S) [map] rebuilt figures/borehole_dvv_map.html ($n stations)"
    else
      echo "$(date +%H:%M:%S) [map] REBUILD FAILED (see logs/map_refresh.log)"
    fi
    # gap audit on every finalize so a truncated/gapped year never slips past unnoticed
    $PY scripts/gap_audit.py 2>/dev/null | grep '^FLAG' | while read -r ln; do
      echo "$(date +%H:%M:%S) [gap] $ln" | tee -a logs/rollout_status.log
    done
    last="$cur"
  fi
  sleep 120
done
