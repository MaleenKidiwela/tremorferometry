#!/usr/bin/env python
"""P3D: the picker's memory-mapped ZNE day store for the P3D roster and years.

`p_build_mm.py` is `p_common.cache_index`-driven and that index SILENTLY DROPS every station that is
not in S15 (Merlin flaw 3).  This is the same builder with `s_common.cache_index` (no S15 filter)
and the P3D azimuth table, written to its OWN tag so `results/mm/day_index.parquet` and
`days_zne_f16.dat` -- the frozen P1/P2 stores -- are never touched.

  python scripts/s_build_mm.py --workers 20 --tag _p3d
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

NDAY = S.NDAY


def one(job):
    irow, sta, year, doy, path, out = job
    try:
        z = np.load(path)
        x = z["x"].astype(np.float64)
        order = [str(c) for c in z["order"]]
        t0 = float(z["t0"])
        nc = x.shape[0]
        if nc == 3:
            x = np.vstack([x[0:1], S.rot_matrix(sta) @ x[1:3]])
        elif nc != 1:
            return dict(i=irow, ok=0, err=f"shape {x.shape} order {order}")
        n = min(x.shape[1], NDAY)
        sc = np.median(np.abs(x[:, :n])) * 1.4826
        if not np.isfinite(sc) or sc <= 0:
            sc = 1.0
        buf = np.zeros((3, NDAY), np.float16)
        buf[:nc, :n] = np.clip(x[:, :n] / sc, -3e4, 3e4).astype(np.float16)
        mm = np.memmap(out, dtype=np.float16, mode="r+", shape=(3, NDAY),
                       offset=irow * 3 * NDAY * 2)
        mm[:] = buf
        mm.flush()
        del mm
        return dict(i=irow, ok=1, t0=t0, n=int(n), scale=float(sc), n_comp=nc,
                    frac_nonzero=float((np.abs(x[:, :n]).sum(0) > 0).mean()),
                    order=",".join(order), err="")
    except Exception as e:
        return dict(i=irow, ok=0, err=f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--years", default=",".join(str(y) for y in range(2019, 2027)))
    ap.add_argument("--tag", default="_p3d")
    ap.add_argument("--stations", default="")
    a = ap.parse_args()
    years = [int(y) for y in a.years.split(",")]
    stas = a.stations.split(",") if a.stations else None
    os.makedirs(S.MM, exist_ok=True)
    idx = S.cache_index(years=years, stas=stas)
    keys = sorted(idx)
    di = pd.DataFrame([dict(sta=s, netsta=S.NETSTA[s], year=y, doy=d, path=idx[(s, y, d)])
                       for (s, y, d) in keys])
    out = f"{S.MM}/days_zne_f16{a.tag}.dat"
    nbytes = len(di) * 3 * NDAY * 2
    print(f"[mm{a.tag}] {len(di)} station-days -> {nbytes/1e9:.1f} GB at {out}", flush=True)
    with open(out, "wb") as f:
        f.truncate(nbytes)
    jobs = [(i, r.sta, r.year, r.doy, r.path, out) for i, r in enumerate(di.itertuples())]
    t0, res = time.time(), {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, r in enumerate(ex.map(one, jobs, chunksize=4)):
            res[r["i"]] = r
            if (k + 1) % 500 == 0:
                print(f"  {k+1}/{len(jobs)}  {time.time()-t0:.0f}s", flush=True)
    for c, d in [("ok", 0), ("t0", np.nan), ("n", 0), ("scale", np.nan), ("n_comp", 0),
                 ("frac_nonzero", np.nan), ("order", ""), ("err", "")]:
        di[c] = [res[i].get(c, d) for i in range(len(di))]
    bad = di[di.ok == 0]
    if len(bad):
        print(f"[warn] {len(bad)} failed\n{bad[['sta','year','doy','err']].head(10).to_string()}")
    ok_rows = [i for i in range(len(res)) if res[i]["ok"] == 1]
    di = di[di.ok == 1].drop(columns=["ok"]).reset_index(drop=True)
    di["row"] = ok_rows
    di.to_parquet(f"{S.MM}/day_index{a.tag}.parquet", index=False)
    print(f"[done] {len(di)} station-days in {time.time()-t0:.0f}s")
    print(di.groupby(["sta"]).size().to_string())


if __name__ == "__main__":
    main()
