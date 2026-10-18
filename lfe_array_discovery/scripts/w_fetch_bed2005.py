#!/usr/bin/env python
"""Step 1c bed -- download + frozen-chain cache of the 2005-09 POLARIS bed.

Stations, window and quiet days come from results/step1c_bed.json (w_setup1c.py), which is frozen
before launch.  Download `r_download_day.fetch_day` verbatim, <= 4 concurrent (agent D keeps its
8), both providers tried in the probe's winning order.  Cache `o_cache_month.cache_one` verbatim.
"""
from __future__ import annotations
import argparse, json, os, sys
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import r_download_day as RD                                                 # noqa: E402
import o_cache_month as OC                                                  # noqa: E402
from w_fetch_bed2012b import PROVIDERS, fetch_any                           # noqa: E402
import w_fetch_bed2012b as FB                                               # noqa: E402

TRACES = f"{W.STUDY}/traces/bed2005_polaris"
FB.TRACES = TRACES                      # fetch_any writes through the module global
CACHE = W.STUDY_CACHE


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cache-workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    cfg = json.load(open(W.STEP1C_BED))
    os.makedirs(TRACES, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)
    D = [date.fromisoformat(cfg["day0"]) + timedelta(days=k)
         for k in range((date.fromisoformat(cfg["day1"])
                         - date.fromisoformat(cfg["day0"])).days + 1)]
    D += [date.fromisoformat(x) for x in cfg["quiet_days"]]
    R = cfg["stations_download"]
    jobs = []
    for r in R:
        first = r.get("prov") or "IRIS"
        provs = [first] + [p for p in PROVIDERS if p != first]
        for d in D:
            jobs.append((r["net"], r["sta"], r.get("band", "HH"), d, provs))
    if a.limit:
        jobs = jobs[: a.limit]
    W.log(f"[fetch1c] {len(R)} stations x {len(D)} days = {len(jobs)} station-days -> {TRACES}")
    recs = []
    with ThreadPoolExecutor(max_workers=min(a.workers, 4)) as ex:
        futs = [ex.submit(fetch_any, n, s, b, d, pv) for (n, s, b, d, pv) in jobs]
        for k, f in enumerate(futs):
            recs.append(f.result())
            if (k + 1) % 100 == 0:
                W.log(f"[fetch1c] {k+1}/{len(jobs)} "
                      f"({sum(1 for x in recs if x['status']=='data')} data)")
    M = pd.DataFrame(recs)
    M.to_csv(f"{W.RES}/step1c_download_manifest.csv", index=False)
    W.log("[fetch1c] " + str(M.status.value_counts().to_dict()))
    W.log("[fetch1c] per station: " + str(M[M.status == "data"].groupby("sta").size().to_dict()))
    cjobs = []
    for r in M[M.status == "data"].itertuples():
        d = date.fromisoformat(r.date)
        cjobs.append((f"{r.net}.{r.sta}", d.year, int(d.timetuple().tm_yday),
                      RD.day_path(TRACES, r.net, r.sta, d), CACHE, False,
                      float(r.coverage_frac) if r.coverage_frac != "" else 0.0))
    W.log(f"[cache1c] {len(cjobs)} station-days -> {CACHE}")
    out = []
    with ThreadPoolExecutor(max_workers=a.cache_workers) as ex:
        for k, rr in enumerate(ex.map(OC.cache_one, cjobs)):
            out.append(rr)
            if (k + 1) % 200 == 0:
                W.log(f"[cache1c] {k+1}/{len(cjobs)}")
    C = pd.DataFrame(out) if out else pd.DataFrame(columns=["netsta", "doy", "err"])
    if "err" not in C.columns:
        C["err"] = ""
    C.to_csv(f"{W.RES}/step1c_cache_manifest.csv", index=False)
    bad = C[C.err.astype(str) != ""]
    W.log(f"[cache1c] ok {len(C)-len(bad)} / {len(C)}; errors {len(bad)}")
    if len(bad):
        W.log("[cache1c] errors: " + str(bad.head(8)[["netsta", "doy", "err"]].to_dict("records")))
    W.log("[cache1c] per station: "
          + str(C[C.err.astype(str) == ""].groupby("netsta").size().to_dict()))


if __name__ == "__main__":
    main()
