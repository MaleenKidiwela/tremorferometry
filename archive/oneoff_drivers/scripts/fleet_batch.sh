#!/bin/bash
# Batch fleet run — sequential (GPU serialized), region-stratified top of data/broadband_fleet_order.csv.
# Each station: fleet_station.sh with auto band-pick. Reports a scoreboard. Usage: nohup bash scripts/fleet_batch.sh &
cd /home/jovyan/tremorferometry
STATIONS="UW.SMW UW.HEBO UW.BBO NC.LCSB CN.TOFB UW.FL2 UW.MPO UW.RNO"
log(){ echo "$(date +%H:%M) [fleet-batch] $*" | tee -a logs/fleet_batch.log; }
log "BATCH START: $STATIONS"
for ns in $STATIONS; do
  net=${ns%.*}; sta=${ns#*.}; s=$(echo $sta | tr A-Z a-z)
  log "=== $sta start ==="
  bash scripts/fleet_station.sh $net $sta 2009-01-01 auto > logs/fleet_${s}.out 2>&1
  gate=$(grep 'GATE:' logs/fleet_${s}.log 2>/dev/null | tail -1)
  log "$sta RESULT: ${gate#*] }"
done
log "BATCH COMPLETE — scoreboard:"; grep 'RESULT:' logs/fleet_batch.log
