#!/bin/bash
# Full GNW pipeline: GPU densify (all years, resumable) -> concat -> stack -> dv/v plot.
# Survives python-level crashes via the densify retry loop (densify skips done years).
cd /home/jovyan/tremorferometry
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p logs figures
L=logs/gnw_gpu_pipeline.log
say(){ echo "$(date '+%H:%M:%S') $*" >> "$L"; }

say "=================== PIPELINE START ==================="

# 1) DENSIFY all years (GPU), retry-resume on crash
ok=0
for att in 1 2 3 4 5 6 7 8; do
  say "--- densify attempt $att ---"
  python scripts/densify_gnw_gpu.py --workers 24 >> "$L" 2>&1
  if grep -q "ALL YEARS DONE" "$L"; then ok=1; say "densify COMPLETE"; break; fi
  say "densify exited early (attempt $att); resuming after 5s"
  sleep 5
done
[ $ok -eq 1 ] || { say "ABORT: densify incomplete after retries"; exit 1; }

# 2) CONCAT per-year csvs -> one file for the stacker
say "--- concat year csvs ---"
python - >> "$L" 2>&1 <<'PY'
import glob, os, pandas as pd
fs=sorted(glob.glob("data/mf_gnw_[12]*.csv"))
out="data/mf_gnw_all.csv"
if os.path.exists(out): os.remove(out)
first=True; tot=0
for f in fs:
    df=pd.read_csv(f); tot+=len(df)
    df.to_csv(out, mode="w" if first else "a", header=first, index=False); first=False
    print(f"  + {f}: {len(df):,}")
print(f"concat {len(fs)} files -> {out}: {tot:,} rows")
PY
[ -s data/mf_gnw_all.csv ] || { say "ABORT: concat produced no file"; exit 1; }

# 3) STACK -> long-window daily npz
say "--- stack ---"
python scripts/build_long_window_daily_all51.py \
  --mf-csv data/mf_gnw_all.csv --network UW --station GNW \
  --out data/long_window_daily_GNW.npz --workers 30 >> "$L" 2>&1
[ -f data/long_window_daily_GNW.npz ] || { say "ABORT: stack produced no npz"; exit 1; }

# 4) dv/v coda 1-3 s plot (NLLB-style)
say "--- dv/v ---"
python scripts/dvv_coda_51.py \
  --npz data/long_window_daily_GNW.npz --window 1.0 3.0 --station GNW \
  --out-csv data/daily_dvv_GNW_coda_1to3.csv \
  --out-fig figures/smoke_dvv_GNW_coda_1to3.png >> "$L" 2>&1

say "=================== PIPELINE COMPLETE ==================="
ls -la figures/smoke_dvv_GNW_coda_1to3.png data/daily_dvv_GNW_coda_1to3.csv >> "$L" 2>&1
