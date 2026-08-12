#!/bin/bash
# Full broadband tier (BH/HH stations) in region-stratified order, sequential (GPU serialized by wait_gpu_free
# in fleet_station.sh). Rolling buffer: free each station's traces after it finishes. Usage: nohup bash scripts/fleet_broadband_all.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
log(){ echo "$(date +%H:%M) [bb-fleet] $*" | tee -a logs/fleet_broadband.log; }
mapfile -t STNS < <($PY -c "
import pandas as pd
c=pd.read_csv('data/broadband_fleet_order.csv')
bb=c[c.chans.str.contains('BHZ|HHZ') & ~c.sta.isin(['SHB','CLRS','PGC'])].sort_values('order')
[print(r.net, r.sta) for _,r in bb.iterrows()]")
log "BROADBAND FLEET START: ${#STNS[@]} stations, region-stratified order"
for ns in "${STNS[@]}"; do
  set -- $ns; net=$1; sta=$2; s=$(echo $sta|tr A-Z a-z)
  grep -q "STATION DONE: $sta" logs/fleet_${s}.log 2>/dev/null && { log "$sta already done, skip"; continue; }
  log "=== $sta start ==="
  bash scripts/fleet_station.sh $net $sta 2009-01-01 auto > logs/fleet_${s}.out 2>&1
  gate=$(grep -E 'GATE:' logs/fleet_${s}.log 2>/dev/null | tail -1)
  log "$sta RESULT: ${gate#*GATE: }"
  # rolling buffer: free traces after EVERY completed station (INCLUDE or FLAG). Re-run candidates (band/
  # provider bugs) re-download from scratch anyway, so nothing worth keeping.
  if grep -q "STATION DONE: $sta" logs/fleet_${s}.log 2>/dev/null; then rm -rf data/waveforms/$net.$sta 2>/dev/null; log "$sta traces freed"; fi
done
log "BROADBAND FLEET COMPLETE"
