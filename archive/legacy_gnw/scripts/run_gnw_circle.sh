#!/bin/bash
# Near-epicenter ("circle") test: densify the 22 SNR>=12 families within 40 km of the 2001
# Nisqually epicenter, stack with resample_poly (no-deconv; within-era test), per-era dv/v.
cd /home/jovyan/tremorferometry
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=logs/gnw_circle.log; mkdir -p logs figures; : > "$L"
say(){ echo "$(date '+%H:%M:%S') $*" >> "$L"; }

say "=== DENSIFY circle families (GPU) ==="
python scripts/densify_gnw_gpu.py --summary-csv data/gnw_circle.summary.csv --min-snr 12 \
  --out-dir data --out-prefix mf_gnwcircle_ --workers 12 >> "$L" 2>&1
grep -q "ALL YEARS DONE" "$L" || { say "ABORT densify"; exit 1; }

say "=== CONCAT ==="
python - >> "$L" 2>&1 <<'PY'
import glob,os,pandas as pd
fs=sorted(glob.glob("data/mf_gnwcircle_[12]*.csv")); out="data/mf_gnwcircle_all.csv"
if os.path.exists(out): os.remove(out)
first=True; tot=0
for f in fs:
    d=pd.read_csv(f); tot+=len(d); d.to_csv(out,mode="w" if first else "a",header=first,index=False); first=False
print(f"concat {len(fs)} files -> {out}: {tot:,} rows")
PY
[ -s data/mf_gnwcircle_all.csv ] || { say "ABORT concat"; exit 1; }

say "=== STACK (resample_poly, no-deconv) ==="
python scripts/build_long_window_resp.py --mf-csv data/mf_gnwcircle_all.csv --no-deconv \
  --network UW --station GNW --out data/long_window_daily_GNWcircle.npz --workers 22 >> "$L" 2>&1
[ -f data/long_window_daily_GNWcircle.npz ] || { say "ABORT stack"; exit 1; }

say "=== DV/V per-era ==="
python scripts/dvv_coda_perera.py --station "GNW near-epicenter circle" \
  --npz data/long_window_daily_GNWcircle.npz --window 1.0 4.0 --era-bounds 2010-09-10,2019-05-07 \
  --workers 22 --out-csv data/daily_dvv_GNWcircle_perera.csv \
  --out-fig figures/smoke_dvv_GNWcircle_perera.png >> "$L" 2>&1
say "=== CIRCLE PIPELINE DONE ==="
