#!/bin/bash
# Targeted 2017 gap-fill for B011 (user-authorized 2026-07-11). The 2017 gap is a DOWNLOAD truncation
# (only days 1-59 were ever fetched; traces since deleted). Re-download full 2017 -> re-densify 2017
# (reliable families, forward-only) -> rebuild 2017 stack -> merge into main npz -> recompute dv/v ->
# delete 2017 traces -> rebuild map. GPU serializes behind the running wave via wait_gpu_free.
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
PY=/home/jovyan/envs/tremorferometry/bin/python
S=B011; s=b011; TAG=B011p90f40
log(){ echo "$(date +%H:%M) [2017fix-$S] $*" | tee -a logs/rollout_status.log; }
wait_gpu_free(){ while ps -eo cmd | grep -E 'densify_gnw_gpu|discover_gpu' | grep -qv grep; do sleep 60; done; }

# 1. Download full 2017 (resumable; skips existing days)
log "download 2017 start"
$PY scripts/download_borehole_3comp.py --net PB --sta $S --start 2017-01-01 --end 2017-12-31 --workers 6 \
   >> logs/download_${s}_2017.log 2>&1
yrs=$(ls data/waveforms/PB.$S/2017/*.mseed 2>/dev/null | wc -l)
log "download 2017 done ($yrs day-files)"
[ "$yrs" -lt 200 ] && { log "ABORT: 2017 download incomplete ($yrs days)"; exit 1; }

# 2. Remove the stale 59-day mf so densify re-runs the full year
rm -f data/mf_${s}p90f40_2017.csv; log "removed stale 59-day mf_${s}p90f40_2017.csv"

# 3. Reliable-families summary (fwd-vs-rev ratio>1.5 = the frozen B011 cert)
RSEL=data/${s}_reliable2017.summary.csv
$PY - data/${s}_disc_p70_2010_2026_m3_sel300.summary.csv data/${s}_fwd_vs_rev_coda.csv $RSEL <<'PY'
import sys, pandas as pd
sel, fvr, out = sys.argv[1:4]
s = pd.read_csv(sel); cert = set(pd.read_csv(fvr).query("ratio>1.5").fam)
s[s.family_id.isin(cert)].to_csv(out, index=False); print(f"reliable: {int(s.family_id.isin(cert).sum())} fams")
PY

# 4. Re-densify 2017 (reliable, forward-only)
wait_gpu_free; log "densify 2017 fwd (reliable) start"
$PY scripts/densify_gnw_gpu.py --templates-npz data/${s}_disc_p70_2010_2026_m3_40hz.npz --summary-csv $RSEL \
   --min-snr 0 --network PB --station $S --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 --top-n 0 \
   --despike-mad 8 --workers 8 --start-year 2017 --end-year 2017 --out-dir data --out-prefix mf_${s}p90f40_ \
   >> logs/densify_${s}p90f40.log 2>&1
nd=$(wc -l < data/mf_${s}p90f40_2017.csv 2>/dev/null); log "densify 2017 done (${nd:-0} rows)"

# 5. Build 2017-only stack -> merge into main npz (Z/H1/H2)
$PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${s}p90f40_2017.csv" --network PB --station $S \
   --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
   --out-prefix data/long_window_daily_${TAG}_2017only >> logs/stack_${s}_2017.log 2>&1
for C in Z H1 H2; do
  $PY scripts/merge_year_into_npz.py data/long_window_daily_${TAG}_${C}.npz \
     data/long_window_daily_${TAG}_2017only_${C}.npz 2017 >> logs/stack_${s}_2017.log 2>&1
done
log "2017 stack merged into main npz"

# 6. Recompute dv/v (Z, H2)
for C in Z H2; do
  $PY scripts/dvv_roll30cal.py --station $S --npz data/long_window_daily_${TAG}_${C}.npz \
     --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_${C}_2to4.csv >> logs/dvv_${s}_2017.log 2>&1
done
log "dv/v recomputed with 2017"

# 7. Delete 2017 traces (rolling buffer) + temp stacks
rm -rf data/waveforms/PB.$S/2017 data/long_window_daily_${TAG}_2017only_*.npz; log "2017 traces + temp stacks deleted"

# 8. Rebuild map (touch summary to also nudge the auto-refresh daemon)
touch data/${s}_3comp_summary.json 2>/dev/null
$PY scripts/build_borehole_dvv_map.py >> logs/map_refresh.log 2>&1 && log "map rebuilt (B011 2017-filled)"

# 9. Verify the gap is gone
$PY - <<PY
import pandas as pd
dv = pd.read_csv("data/daily_dvv_${TAG}_Z_2to4.csv"); dv['date'] = pd.to_datetime(dv.date)
yr = dv.groupby(dv.date.dt.year).date.agg(lambda x: x.dt.normalize().nunique())
print(f"[verify] B011 daily_dvv days/yr 2016={yr.get(2016)} 2017={yr.get(2017)} 2018={yr.get(2018)}")
PY
log "COMPLETE"
