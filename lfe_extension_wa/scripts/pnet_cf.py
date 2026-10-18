#!/usr/bin/env python
"""PNET S characteristic function for every (station, day) of the WA permanent network.

Exactly lfe_picker's p_scan.cf_day (60-s windows, step 20 s, central third kept, per-window 3-C peak
normalisation, bf16 autocast) on the WA cache rotated to ZNE (pnet_common.to_zne).
  python scripts/pnet_cf.py --model <pt> --tag w05 --windows V2008,V2010,A2011
  -> results/pnet/cf_<tag>/<NET.STA>/<date>.npy  (float16, NDAY samples) + cf_<tag>/index.csv
"""
from __future__ import annotations
import os, sys, time, argparse
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC
sys.path.insert(0, PC.PICKER)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--windows", default="V2008,V2010,A2011")
    ap.add_argument("--days", default="")
    ap.add_argument("--stations", default="")
    a = ap.parse_args()
    import torch
    import p_scan
    model, meta = p_scan.load_model(a.model)
    S = W.load_stations()
    stas = a.stations.split(",") if a.stations else list(S.netsta)
    days = [d for w in a.windows.split(",") for d, _ in W.window_days(w)]
    if a.days:
        days = [d for d in days if d in set(a.days.split(","))]
    jobs = [(ns, d) for d in days for ns in stas if os.path.exists(W.cache_path(ns, d))
            and not os.path.exists(PC.cf_path(a.tag, ns, d))]
    W.log(PC.LOG, f"[cf:{a.tag}] model {a.model} ({meta.get('arch')}, win {meta.get('win_s')} s); {len(jobs)} station-days")
    t0 = time.time()
    rows = []
    ex = ThreadPoolExecutor(4)
    futs = [ex.submit(PC.load_zne, ns, d) for ns, d in jobs[:8]]
    for k, (ns, d) in enumerate(jobs):
        if k + 8 < len(jobs):
            futs.append(ex.submit(PC.load_zne, *jobs[k + 8]))
        z = futs[k].result()
        futs[k] = None
        if z is None:
            continue
        cf = p_scan.cf_day(model, meta, z["x"])
        p = PC.cf_path(a.tag, ns, d)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        np.save(p + ".tmp.npy", cf.astype(np.float16))
        os.replace(p + ".tmp.npy", p)
        rows.append(dict(netsta=ns, date=d, node=z["node"], n_horiz=z["n_horiz"], frac_valid=float(z["m"].mean()),
                         cf_mean=float(cf.mean()), cf_p99=float(np.quantile(cf[::20], 0.99)), cf_frac03=float((cf > 0.3).mean())))
        if (k + 1) % 200 == 0:
            W.log(PC.LOG, f"[cf:{a.tag}] {k+1}/{len(jobs)} {time.time()-t0:.0f}s")
    ex.shutdown()
    ip = f"{PC.OUT}/cf_{a.tag}/index.csv"
    I = pd.DataFrame(rows)
    if os.path.exists(ip) and len(I):
        I = pd.concat([pd.read_csv(ip), I]).drop_duplicates(["netsta", "date"], keep="last")
    elif os.path.exists(ip):
        I = pd.read_csv(ip)
    os.makedirs(os.path.dirname(ip), exist_ok=True)
    I.to_csv(ip, index=False)
    W.log(PC.LOG, f"[cf:{a.tag}] done {len(rows)} station-days {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
