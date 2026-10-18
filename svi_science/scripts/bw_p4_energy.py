#!/usr/bin/env python
"""P4 quiet flag, step 1: CN4 hourly tremor-band energy (lfe_array_discovery/scripts/w_energy_hours.py statistic):
per station-hour the median over the hour of the 3-C summed squared amplitude of the 20 Hz 1-8 Hz cache (>= 1000 finite
samples), UTC hours of the date.  Caches: backward first, then the forward chain (same order as P5).  CPU only.

  python bw_p4_energy.py --years 2010,2011 [--workers 8]   -> p4/tables/energy/cn4_hourly_<Y>.parquet
"""
import os, sys, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

LFE = "/home/jovyan/tremorferometry/lfe_extension"
CACHES = [f"{LFE}/cache20hz_backward", f"{LFE}/cache20hz_r1", f"{LFE}/cache20hz",
          "/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz", "/home/jovyan/tremorferometry/lfe_picker/cache20hz"]
CN4 = {"PGC": "CN.PGC", "VGZ": "CN.VGZ", "NLLB": "CN.NLLB", "SNB": "CN.SNB"}
OUT = "/home/jovyan/tremorferometry/svi_science/results/backward/p4/tables/energy"
FS = 20.0


def one(args):
    sta, date = args
    d = pd.Timestamp(date); y, doy = d.year, d.dayofyear
    p = next((q for q in (f"{c}/{CN4[sta]}/{y}/{doy:03d}.npz" for c in CACHES) if os.path.exists(q)), None)
    if p is None:
        return []
    try:
        z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"]); fs = float(z["fs"]) if "fs" in z.files else FS
    except Exception:
        return [dict(date=date, sta=sta, hour=-1, e=np.nan, n=0, cache=p.split("/")[-4])]
    e = (x ** 2).sum(0)
    day0 = d.value / 1e9
    out = []
    for h in range(24):
        i0 = int(round((day0 + 3600 * h - t0) * fs)); i1 = int(round((day0 + 3600 * (h + 1) - t0) * fs))
        i0 = max(i0, 0); i1 = min(i1, len(e))
        v = e[i0:i1] if i1 > i0 else np.zeros(0)
        v = v[np.isfinite(v) & (v > 0)]                     # exact zeros = gap fill, not data
        out.append(dict(date=date, sta=sta, hour=h, e=float(np.median(v)) if v.size >= 1000 else np.nan, n=int(v.size),
                        cache=p.split("/")[-4]))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--years", required=True); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for y in [int(v) for v in a.years.split(",") if v]:
        f = f"{OUT}/cn4_hourly_{y}.parquet"
        if os.path.exists(f) and not a.force:
            print(f"[energy] {y} exists, skip", flush=True); continue
        days = [str(d.date()) for d in pd.date_range(f"{y}-01-01", f"{y}-12-31")]
        jobs = [(s, d) for d in days for s in CN4]
        rows = []
        with ProcessPoolExecutor(a.workers) as ex:
            for r in ex.map(one, jobs, chunksize=8):
                rows += r
        E = pd.DataFrame(rows)
        E.to_parquet(f + ".tmp", index=False); os.replace(f + ".tmp", f)
        ok = E[E.hour >= 0].groupby("sta").apply(lambda g: g.date[g.e.notna()].nunique())
        print(f"[energy] {y}: {len(E)} rows; days with any hour per station {ok.to_dict()}", flush=True)


if __name__ == "__main__":
    main()
