#!/usr/bin/env python
"""Track D, D1(b) -- daily noise (median 8-s window energy per component) for the REFEREE years,
computed from the 20 Hz caches with EXACTLY the estimator Track C uses in production
(`c_cache_month.med_win_energy`, 8 s at 20 Hz, float64 prefix sums).  Output has the same schema as
`results/trackC/daily_noise_<YYYY-MM>.parquet` so the day term of the completeness model is formed the
same way in the referee years and in the production record."""
from __future__ import annotations
import argparse, os
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

CACHES = ["/home/jovyan/tremorferometry/lfe_extension/cache20hz",
          "/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz"]
FS, WIN_S = 20.0, 8.0
NETSTA = dict(PGC="CN.PGC", VGZ="CN.VGZ", NLLB="CN.NLLB", SNB="CN.SNB", LRIV="UW.LRIV", B001="PB.B001",
              B003="PB.B003", B004="PB.B004", B005="PB.B005", B006="PB.B006", B007="PB.B007",
              B009="PB.B009", B010="PB.B010", B011="PB.B011", B926="PB.B926")


def med_win_energy(X, n):
    c = np.concatenate([np.zeros((X.shape[0], 1)), np.cumsum(X.astype(np.float64) ** 2, axis=1)], 1)
    return np.median(c[:, n:] - c[:, :-n], axis=1)


def one(job):
    netsta, year, doy = job
    p = next((f"{c}/{netsta}/{year}/{doy:03d}.npz" for c in CACHES
              if os.path.exists(f"{c}/{netsta}/{year}/{doy:03d}.npz")), None)
    if p is None:
        return dict(netsta=netsta, year=year, doy=doy, status="missing", e8_Z=np.nan, e8_N=np.nan,
                    e8_E=np.nan, frac_nonzero=np.nan, order="")
    z = np.load(p); X = z["x"].astype(np.float64)
    me = med_win_energy(X, int(round(WIN_S * FS)))
    return dict(netsta=netsta, year=year, doy=doy, status="cached", e8_Z=me[0], e8_N=me[1], e8_E=me[2],
                frac_nonzero=float((np.abs(X).sum(0) > 0).mean()), order=",".join(list(z["order"])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/home/jovyan/tremorferometry/lfe_extension/results/trackD/daily_noise_ref.parquet")
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    jobs = []
    for s in NETSTA.values():
        jobs += [(s, 2012, d) for d in range(228, 280)]
        jobs += [(s, 2010, d) for d in range(222, 257)]
    with ProcessPoolExecutor(a.workers) as ex:
        rows = list(ex.map(one, jobs, chunksize=8))
    D = pd.DataFrame(rows)
    D["date"] = [str((pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=int(d) - 1)).date())
                 for y, d in zip(D.year, D.doy)]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    D.to_parquet(a.out, index=False)
    print(D.groupby(["year", "status"]).size().to_string())
    print(f"[write] {a.out}  ({len(D)} station-days, {int((D.status=='cached').sum())} cached)")


if __name__ == "__main__":
    main()
