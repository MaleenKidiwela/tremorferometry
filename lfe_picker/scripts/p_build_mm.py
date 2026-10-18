#!/usr/bin/env python
"""Build the memory-mapped ZNE day store used by every P1/P2 script.

For every cached station-day of the S15 roster (2008-2013), read the frozen 20 Hz / 1-8 Hz cache,
ROTATE the two horizontals to N/E with the StationXML azimuths (boreholes: EH1/EH2 are NOT at 0/90),
scale by the day's robust amplitude so float16 has headroom, and write into one
(Ndays, 3, 1728000) float16 memmap.  ~39 GB; the pod has 1.5 TB of RAM, so it stays in page cache.

Data loading, not GPU, is the bottleneck (Merlin F11) -- this pays for it once.

  python scripts/p_build_mm.py --workers 24
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P


def one(job):
    irow, sta, year, doy, path, out = job
    try:
        z = np.load(path)
        x = z["x"].astype(np.float64)
        order = [str(c) for c in z["order"]]
        t0 = float(z["t0"])
        if len(order) != 3 or x.shape[0] != 3:
            return dict(i=irow, ok=0, err=f"shape {x.shape} order {order}")
        # cache order is (Z, N|1, E|2) resolved by o_chanorder -> rotate slots 1,2 to true N,E
        Rm = P.rot_matrix(sta)
        h = Rm @ x[1:3]
        y = np.vstack([x[0:1], h])
        n = min(y.shape[1], P.NDAY)
        sc = np.median(np.abs(y[:, :n])) * 1.4826
        if not np.isfinite(sc) or sc <= 0:
            sc = 1.0
        buf = np.zeros((3, P.NDAY), np.float16)
        buf[:, :n] = np.clip(y[:, :n] / sc, -3e4, 3e4).astype(np.float16)
        mm = np.memmap(out, dtype=np.float16, mode="r+", shape=(3, P.NDAY),
                       offset=irow * 3 * P.NDAY * 2)
        mm[:] = buf
        mm.flush(); del mm
        frac = float((np.abs(y[:, :n]).sum(0) > 0).mean())
        return dict(i=irow, ok=1, t0=t0, n=int(n), scale=float(sc), frac_nonzero=frac,
                    order=",".join(order), err="")
    except Exception as e:
        return dict(i=irow, ok=0, err=f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--years", default="2008,2009,2010,2011,2012,2013")
    ap.add_argument("--tag", default="", help="suffix for a SECOND store (e.g. _ssl)")
    ap.add_argument("--stride", type=int, default=1, help="keep every Nth day per station")
    a = ap.parse_args()
    years = [int(y) for y in a.years.split(",")]
    os.makedirs(P.MM, exist_ok=True)

    idx = P.cache_index(years=years)
    keys = sorted(idx)
    if a.stride > 1:
        keep, seen = [], {}
        for k in keys:
            c = seen.get(k[0], 0)
            seen[k[0]] = c + 1
            if c % a.stride == 0:
                keep.append(k)
        keys = keep
    di = pd.DataFrame([dict(sta=s, netsta=P.NETSTA[s], year=y, doy=d, path=idx[(s, y, d)])
                       for (s, y, d) in keys])
    out = f"{P.MM}/days_zne_f16{a.tag}.dat"
    nbytes = len(di) * 3 * P.NDAY * 2
    print(f"[mm] {len(di)} station-days -> {nbytes/1e9:.1f} GB at {out}", flush=True)
    with open(out, "wb") as f:                       # sparse allocate
        f.truncate(nbytes)

    jobs = [(i, r.sta, r.year, r.doy, r.path, out) for i, r in enumerate(di.itertuples())]
    t0 = time.time()
    res = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, r in enumerate(ex.map(one, jobs, chunksize=4)):
            res[r["i"]] = r
            if (k + 1) % 250 == 0:
                print(f"  {k+1}/{len(jobs)}  {time.time()-t0:.0f}s", flush=True)
    for c, d in [("ok", 0), ("t0", np.nan), ("n", 0), ("scale", np.nan),
                 ("frac_nonzero", np.nan), ("order", ""), ("err", "")]:
        di[c] = [res[i].get(c, d) for i in range(len(di))]
    bad = di[di.ok == 0]
    if len(bad):
        print(f"[warn] {len(bad)} failed:\n{bad[['sta','year','doy','err']].head(20).to_string()}")
    di = di[di.ok == 1].drop(columns=["ok"]).reset_index(drop=True)
    # rows are positional in the memmap -> keep the ORIGINAL row id
    ok_rows = [i for i in range(len(res)) if res[i]["ok"] == 1]
    di["row"] = ok_rows
    di.to_parquet(f"{P.MM}/day_index{a.tag}.parquet", index=False)
    print(f"[done] {len(di)} good station-days in {time.time()-t0:.0f}s")
    print(di.groupby(["year"]).size())


if __name__ == "__main__":
    main()
