#!/bin/bash
# Download-ahead: pre-download broadband stations in fleet order, staying <=K ahead of processing, so the
# network works during each station's GPU densify instead of sitting idle. The processing driver's own
# download step then becomes a fast no-op. Usage: nohup bash scripts/fleet_download_ahead.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
K=3   # stay at most 3 stations ahead of processing
log(){ echo "$(date +%H:%M) [dl-ahead] $*" | tee -a logs/fleet_dl_ahead.log; }
mapfile -t STNS < <($PY -c "
import pandas as pd
c=pd.read_csv('data/broadband_fleet_order.csv')
bb=c[c.chans.str.contains('BHZ|HHZ') & ~c.sta.isin(['SHB','CLRS','PGC'])].sort_values('order')
[print(r.net, r.sta) for _,r in bb.iterrows()]")
log "start: ${#STNS[@]} stations, staying $K ahead"
i=0
for ns in "${STNS[@]}"; do
  set -- $ns; net=$1; sta=$2; s=$(echo $sta|tr A-Z a-z)
  while :; do
    j=$(grep -c 'RESULT:' logs/fleet_broadband.log 2>/dev/null || echo 0)
    [ $((i - j)) -lt $K ] && break
    sleep 60
  done
  CHAN=$($PY scripts/pick_band.py $net $sta 2>/dev/null); [ -n "$CHAN" ] || CHAN="BH?,HH?,EH?,SH?"
  log "pre-download $sta (chan $CHAN)  [ahead=$((i-j))]"
  $PY scripts/download_broadband.py --net $net --sta $sta --start 2009-01-01 --end 2026-08-01 --workers 5 --chan "$CHAN" >> logs/dl_${s}.log 2>&1
  i=$((i+1))
done
log "download-ahead COMPLETE"
