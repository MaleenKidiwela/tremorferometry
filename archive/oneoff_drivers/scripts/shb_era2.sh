#!/bin/bash
# SHB era-2 (HHZ 2020-2026): densify the SAME 331 templates on the HHZ record, PER-ERA reference (PGC-style),
# causality. Combined with era-1 (BHZ 2010-2019, 56 certified) = SHB full two-segment record.
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
log(){ echo "$(date +%H:%M) [shb-e2] $*" | tee -a logs/shb_era2.log; }
log "download SHB HHZ 2020-2026"
$PY scripts/download_broadband.py --net CN --sta SHB --start 2020-01-01 --end 2026-08-01 --workers 6 --chan "HH?" \
  > logs/dl_shbhhz.log 2>&1 || { log "download FAILED"; exit 1; }
log "download HHZ: $(grep -oE 'DONE.*' logs/dl_shbhhz.log | tail -1)"
log "densify 331 templates on HHZ era 2020-2026 (100->40)"
$PY scripts/densify_gnw_gpu.py --templates-npz data/shb_disc.npz --summary-csv data/shb_disc.summary.csv \
  --min-snr 0 --network CN --station SHB --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 --top-n 0 \
  --despike-mad 8 --workers 8 --start-year 2020 --end-year 2026 --out-dir data --out-prefix mf_shbhhz_ \
  > logs/densify_shbhhz.log 2>&1 || { log "densify FAILED"; exit 2; }
log "stacks + dv/v (PER-ERA ref) + causality"
$PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_shbhhz_*.csv" --network CN --station SHB --fs 40 \
  --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
  --out-prefix data/long_window_daily_SHBhhz > logs/stack_shbhhz.log 2>&1 || { log "stacks FAILED"; exit 3; }
$PY scripts/dvv_roll30cal.py --station SHB --npz data/long_window_daily_SHBhhz_Z.npz --window 2 4 --origin-anchor \
  --workers 10 --out data/daily_dvv_SHBhhz_Z_2to4.csv > logs/dvv_shbhhz.log 2>&1 || { log "dvv FAILED"; exit 4; }
$PY scripts/finalize_causality.py SHBhhz SHBhhz > logs/finalize_shbhhz.log 2>&1 || { log "finalize FAILED"; exit 5; }
zc=$($PY -c "import json;print(json.load(open('data/shbhhz_3comp_summary.json'))['z_certified'])" 2>/dev/null)
log "SHB ERA-2 DONE: ${zc} certified (HHZ 2020-2026) | era-1 = 56 (BHZ 2010-2019) -> SHB full record"
