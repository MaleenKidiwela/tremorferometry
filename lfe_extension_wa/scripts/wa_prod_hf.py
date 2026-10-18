#!/usr/bin/env python
"""WA PRODUCTION: HF ratio column for EVERY candidate (data, scrambled controls, time-reversed control).

hfr = E(10-19 Hz) / E(2-8 Hz) over the 12-s window at the family-day's HF station (wa_prod_mf: clear template first,
then w; native >= 40 Hz), from the NATIVE-rate raw trace, 3 components summed, Butterworth-4 zero-phase designed at the
native rate (step1_features.hf_station_day, unchanged); hfr_rel = log10(hfr) - median log10 of 300 random windows of the
same station-day.  Scrambled-control rows are read at the HF station's ROLLED position (wa_prod_mf hf_pos).
Status per (date, HF station) -> OUT/status/hf_<run>_<shard>.csv; guarded (> 3 % error / missing fails).

  python scripts/wa_prod_hf.py --run prod --workers 8 [--shard 0 --nshard 1] [--wait-mf]
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.patch_wa()
import step1_features as S1


def day_job(args):
    run, ds = args
    dp = f"{P.OUT}/mf/{run}/det_{ds}.parquet"
    if "k" not in pq.read_schema(dp).names:                 # empty day (no observable family): wa_prod_mf writes only fid, t
        P.ensure(f"{P.OUT}/hf/{run}")                       # (seen in the fixed-10 run, 2026-10-07; production never had one)
        pd.DataFrame(dict(k=pd.Series(dtype="int64"), hfr=pd.Series(dtype="float64"),
                          hfr_rel=pd.Series(dtype="float64"))).to_parquet(f"{P.OUT}/hf/{run}/hf_{ds}.parquet", index=False)
        return ds, [(ds, "", "empty", "no detections file content")], 0
    D = pd.read_parquet(dp, columns=["k", "hf_sta", "hf_pos"])
    J = json.load(open(f"{P.OUT}/mf/{run}/run_{ds}.json"))
    nodes = J.get("nodes", {})
    out = pd.DataFrame(dict(k=D.k.to_numpy(), hfr=np.nan, hfr_rel=np.nan))
    st = []
    for ns, gi in D.groupby("hf_sta").groups.items():
        if not ns:
            continue
        node = nodes.get(ns)
        if node is None:
            st.append((ds, ns, "error", "HF station not in the day's nodes")); continue
        if P.raw_path(ns, ds) is None:
            st.append((ds, ns, "missing", "no raw for the HF station")); continue
        try:
            hr, ref = S1.hf_station_day(ns, ds, node, D.loc[gi, "hf_pos"].to_numpy(float))
        except Exception as e:
            st.append((ds, ns, "error", f"{type(e).__name__}: {str(e)[:80]}")); continue
        out.loc[gi, "hfr"] = hr
        with np.errstate(divide="ignore", invalid="ignore"):
            out.loc[gi, "hfr_rel"] = np.log10(hr) - ref
        st.append((ds, ns, "ok" if np.isfinite(hr).any() else "nodata", ""))
    P.ensure(f"{P.OUT}/hf/{run}")
    out.to_parquet(f"{P.OUT}/hf/{run}/hf_{ds}.parquet.tmp", index=False)
    os.replace(f"{P.OUT}/hf/{run}/hf_{ds}.parquet.tmp", f"{P.OUT}/hf/{run}/hf_{ds}.parquet")
    return ds, st, len(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--wait-mf", default="", help="poll until this marker exists, processing day files as they appear")
    a = ap.parse_args()
    P.ensure(f"{P.OUT}/hf/{a.run}")
    sp_ = f"{P.OUT}/status/hf_{a.run}.csv"
    status = [tuple(x) for x in pd.read_csv(sp_, dtype=str).fillna("").itertuples(index=False)] if os.path.exists(sp_) else []
    t0 = time.time(); ndone = 0
    with ProcessPoolExecutor(a.workers) as ex:
        while True:
            fin = (not a.wait_mf) or os.path.exists(a.wait_mf)
            have = sorted(os.path.basename(p)[4:14] for p in glob.glob(f"{P.OUT}/mf/{a.run}/det_*.parquet"))
            have = [d for d in have if os.path.exists(f"{P.OUT}/mf/{a.run}/run_{d}.json")]
            todo = [d for d in have if not os.path.exists(f"{P.OUT}/hf/{a.run}/hf_{d}.parquet")]
            if todo:
                fut = [ex.submit(day_job, (a.run, d)) for d in todo]
                for f in as_completed(fut):
                    ds, st, nr = f.result()
                    status += st; ndone += 1
                    if ndone % 50 == 0:
                        P.write_status("hf", a.run, status)
                        P.log(f"[hf {a.run}] {ndone} days {time.time() - t0:.0f}s")
            if fin and not todo:
                break
            if not todo:
                time.sleep(120)
    P.write_status("hf", a.run, status)
    P.log(f"[hf {a.run}] done {ndone} days {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
