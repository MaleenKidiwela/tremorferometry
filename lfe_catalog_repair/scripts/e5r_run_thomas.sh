#!/bin/bash
# E5 step 5 STEP F -- Thomas metrics on OUR 2012 catalog with the UNCHANGED e5_thomas_harness.py, plus a Bostock random-decimation
# reference at our measured overall recall (same estimators, same nulls).  Magnitude proxy is uncalibrated -> --mc none only.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=8
set -e
python scripts/e5_thomas_harness.py --catalog results/e5_catalog_2012_sVI.parquet --tag ours2012_L8 --mc none > logs/e5r_thomas_L8.log 2>&1 &
python scripts/e5_thomas_harness.py --catalog results/e5_catalog_2012_sVI_L12.parquet --tag ours2012_L12 --mc none > logs/e5r_thomas_L12.log 2>&1 &
python - <<'PY'
import json, numpy as np, pandas as pd
RES="results"; s=json.load(open(f"{RES}/e5r_score_e5r_snr.json"))
frac=float(s["L8_T8"]["bars"]["recall_overall"]["cat"])
c=pd.read_csv(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
t=pd.to_datetime(dict(year=c.year,month=c.month,day=c.day,hour=c.hour,minute=c["min"],second=c.sec.astype(int)))+pd.to_timedelta(c.sec%1,unit="s")
c=c[(t>="2012-08-15")&(t<"2012-10-16")]
rng=np.random.default_rng(5); keep=rng.random(len(c))<frac
c[keep].to_csv(f"{RES}/e5r_bostock2012_random_decimated.csv", index=False)
print(f"[decimation] Bostock 2012 randomly decimated to our catalog-level recall {frac:.3f}: {int(keep.sum())} of {len(c)} events")
PY
python scripts/e5_thomas_harness.py --catalog results/e5r_bostock2012_random_decimated.csv --tag bostock2012_dec --mc none > logs/e5r_thomas_dec.log 2>&1 &
wait
echo "[$(date -u +%FT%TZ)] STEP F DONE"
