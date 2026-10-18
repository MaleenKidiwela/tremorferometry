#!/usr/bin/env python
"""SH_V2 prep: burst mute masks, spectral lines and guarded channel refs for every station-day the
Shelly detector / stacker will read (union of build + test days).  Writes results/sh_v2/prep/<date>.npz
and results/sh_v2/prep_summary.csv.  See sh_v2_common.py.

  python scripts/sh_v2_prep.py --workers 12
"""
import os, sys, argparse, json, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v2_common import *
from concurrent.futures import ProcessPoolExecutor

_G = {}
def _job(d):
    warnings.filterwarnings("ignore")
    try:
        return prep_day(d, _G["stas"], _G["C"], _G["nref"])
    except Exception as e:
        log(f"  prep {d}: FAILED {type(e).__name__}: {e}"); raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--templates", default=f"{RES}/aa_shelly/templates6_build")
    ap.add_argument("--days", default=f"{RES}/aa_shelly/stage2_days_build.json,{RES}/aa_shelly/stage2_days_test.json")
    ap.add_argument("--t0", default="2018-01-01"); ap.add_argument("--t1", default="2024-12-01")
    ap.add_argument("--len", type=float, default=6.0); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(PREP, exist_ok=True)
    M = pd.read_csv(a.templates + ".csv"); M = M[M.source == "Shelly2025"]
    _G["stas"] = sorted(M.netsta.unique()); _G["C"] = coords(); _G["nref"] = int(round(a.len * FS))
    days = set()
    for p in a.days.split(","):
        for k, v in json.load(open(p)).items():
            if k.split("|")[0] == "Shelly2025": days |= {pd.Timestamp(x) for x in v}
    days = sorted(d for d in days if pd.Timestamp(a.t0) <= d < pd.Timestamp(a.t1))
    if a.limit: days = days[:a.limit]
    log(f"[prep] {len(days)} days, {len(_G['stas'])} stations, workers={a.workers}")
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(_job, [str(d.date()) for d in days])):
            rows += r
            if i % 50 == 0: log(f"  {i + 1}/{len(days)}")
    old = f"{V2}/prep_summary.csv"
    df = pd.DataFrame(rows)
    if os.path.exists(old) and len(df):
        df = pd.concat([pd.read_csv(old), df]).drop_duplicates(["day", "netsta"], keep="last")
    elif os.path.exists(old):
        df = pd.read_csv(old)
    df.to_csv(old, index=False)
    if len(df):
        log(f"[prep] station-days {len(df)}; muted frac mean {df.muted_frac.mean():.4f}; "
            f"station-days with >=1 dead channel {int((df.dead_ch > 0).sum())}; dead stations {int(df.station_dead.sum())}; "
            f"lines {int(df.n_lines.sum())}")
if __name__ == "__main__": main()
