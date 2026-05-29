#!/bin/bash
# HDW discovery stage-2 (clustering) with MEMORY-SAFE settings + a watchdog that kills it before
# it can OOM the pod. Reuses cached candidates (stage 1 already done). The all-pairs CC scales as
# N^2*shifts per bin, so lower --max-bin-candidates (700) + fewer --workers (10) keeps RAM at a
# few GB. Watchdog hard-kills if cgroup memory ever exceeds 150 GB (cap is 187 GB).
cd /home/jovyan/tremorferometry
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=logs/hdw_discovery2.log; : > "$L"
echo "$(date '+%H:%M:%S') launching stage-2 (max-bin 700, workers 10, watchdog @150G)" >> "$L"

python scripts/discover_nllb_pnsn_driven.py --station HDW \
  --bbox 46.5 48.7 -124.3 -122.4 \
  --candidates-out data/hdw_pnsn_candidates.parquet --use-cached-candidates \
  --out data/hdw_pnsn_families.npz \
  --max-bin-candidates 700 --workers 10 >> "$L" 2>&1 &
DPID=$!

# memory watchdog — watches ANON (real, OOM-risk) memory, NOT memory.current (which includes
# reclaimable page cache from reading day files and would false-trip).
while kill -0 "$DPID" 2>/dev/null; do
  mem=$(awk '/^anon /{print $2}' /sys/fs/cgroup/memory.stat)
  g=$((mem/1000000000))
  echo "$(date '+%H:%M:%S') watchdog: anon ${g}G" >> "$L"
  if [ "$mem" -gt 130000000000 ]; then
    echo "$(date '+%H:%M:%S') WATCHDOG TRIPPED at anon ${g}G — killing discovery to prevent OOM" >> "$L"
    pkill -9 -f "[m]ultiprocessing.spawn"; kill -9 "$DPID" 2>/dev/null
    echo "WATCHDOG_KILLED" >> "$L"
    exit 1
  fi
  sleep 20
done
echo "$(date '+%H:%M:%S') stage-2 finished (watchdog saw no breach)" >> "$L"
