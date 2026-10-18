#!/usr/bin/env python
"""Run a trained P1b model over every station-day of a year and store the declustered peaks.

  python scripts/p1b_scan.py --tag pnet_w05 --year 2012
  python scripts/p1b_scan.py --tag pnet_w05 --year 2011      # ONE shot, referee-gated

Output: results/p1/scan_<tag>_<year>.npz with, per station-day, the declustered peak sample indices
and values (top 6000), plus the day's valid-sample count.  Scoring is a separate step so a scan is
never repeated for a new metric.
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_scan as S

OUT = f"{P.RES}/p1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--model", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--sep-s", type=float, default=2.0)
    ap.add_argument("--save-cf", default="", help="optional directory for full CF arrays (float16)")
    ap.add_argument("--doys", default="")
    ap.add_argument("--stations", default="")
    a = ap.parse_args()
    assert a.year != 2013, "2013 is the sealed referee year — forbidden here"
    mp = a.model or f"{OUT}/models/{a.tag}.pt"
    model, meta = S.load_model(mp)
    st = P.DayStore()
    di = st.di[st.di.year == a.year]
    if a.doys:
        di = di[di.doy.isin([int(x) for x in a.doys.split(",")])]
    if a.stations:
        di = di[di.sta.isin(a.stations.split(","))]
    di = di.reset_index(drop=True)
    print(f"[scan] {a.tag} {a.year}: {len(di)} station-days, arch={meta['arch']} "
          f"win={meta['win_s']}s", flush=True)
    if a.save_cf:
        os.makedirs(a.save_cf, exist_ok=True)

    rows, store = [], {}
    t0 = time.time()
    for k, r in enumerate(di.itertuples()):
        i = st.idx(r.sta, r.year, r.doy)
        n = int(st.nvalid[i])
        x = np.asarray(st.get_i(i)[:, :n], np.float32)
        cf = S.cf_day(model, meta, x)
        idx, val = S.decluster(cf, sep_s=a.sep_s)
        key = f"{r.sta}_{r.doy}"
        store[f"i_{key}"] = idx.astype(np.int32)
        store[f"v_{key}"] = val.astype(np.float32)
        rows.append(dict(sta=r.sta, year=r.year, doy=r.doy, t0=float(st.t0[i]), n=n,
                         npeak=len(idx), cf_med=float(np.median(cf)), cf_p99=float(np.quantile(cf, 0.99))))
        if a.save_cf:
            np.save(f"{a.save_cf}/{r.sta}_{r.year}_{r.doy}.npy", cf.astype(np.float16))
        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{len(di)}  {time.time()-t0:.0f}s", flush=True)
    meta_df = pd.DataFrame(rows)
    out = a.out or f"{OUT}/scan_{a.tag}_{a.year}.npz"
    np.savez_compressed(out, **store)
    meta_df.to_csv(out.replace(".npz", "_days.csv"), index=False)
    print(f"[done] {len(di)} station-days in {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
