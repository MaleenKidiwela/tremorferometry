#!/usr/bin/env python
"""P3D DISCOVERY ARM B -- seeds: the top 1 % of tremor-hour windows by INCOHERENT S energy.

Merlin's Arm B starts from LOUD windows, not from linked pairs.  The seed statistic is the MEDIAN
per-site z of the 3-C energy in the S window (placed at the grid-centre moveout, the same placement
the pair statistic uses), z being relative to that station-day's own median window energy.  The
MEDIAN over sites -- not the sum -- is deliberate and pre-registered: Step 3 showed the site SUM is
dominated by the PGC cluster, so a sum would let one loud site nominate every seed.

  python scripts/s_discover_seeds.py --net B --top-frac 0.01
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_discover_pairs as DP                                              # noqa: E402

OUT = S.P3D
WIN_S = 12.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="B")
    ap.add_argument("--top-frac", type=float, default=0.01)
    ap.add_argument("--out", default=f"{OUT}/disc_B_seeds.parquet")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    ST = pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta")
    stas = DP.NETS[a.net]
    site_of = {s: ST.site[s] for s in stas}
    sites = sorted(set(site_of.values()))
    C, idx, days, cat_end = DP.candidates(a.net, DP.EPISODES)
    cst, ctr, cdep = DP.station_constants(stas, ST)
    n = int(round(WIN_S * S.FS))
    Z = np.full((len(C), len(stas)), np.nan, np.float32)
    dc = S.DayCache(maxn=4, years=list(range(2019, 2027)), stas=stas)
    for si, s in enumerate(stas):
        for (y, doy), g in C.groupby(["year", "doy"]):
            v = dc.get(s, int(y), int(doy))
            if v is None:
                continue
            t0d, X = v
            e = (np.asarray(X, np.float64) ** 2).sum(0)
            cs = np.concatenate([[0.0], np.cumsum(e)])
            i0 = np.rint((g.t.to_numpy() + cst[s] - t0d) * S.FS).astype(np.int64)
            ok = (i0 >= 0) & (i0 + n < len(e))
            w = np.full(len(g), np.nan)
            w[ok] = cs[i0[ok] + n] - cs[i0[ok]]
            # day background: the median over a 1-s-stride sliding window of the same length
            st = np.arange(0, len(e) - n, int(S.FS))
            bg = np.median(cs[st + n] - cs[st])
            mad = 1.4826 * np.median(np.abs((cs[st + n] - cs[st]) - bg))
            Z[g.index.to_numpy(), si] = (w - bg) / max(mad, 1e-9)
        print(f"  [seed] {s}: {np.isfinite(Z[:, si]).sum()}/{len(C)}  {time.time()-t0:.0f}s",
              flush=True)
    # site-median (a site's stations are averaged first, then the median over sites)
    per_site = np.full((len(C), len(sites)), np.nan, np.float32)
    for qi, q in enumerate(sites):
        cols = [i for i, s in enumerate(stas) if site_of[s] == q]
        per_site[:, qi] = np.nanmean(Z[:, cols], axis=1)
    med = np.nanmedian(per_site, axis=1)
    C = C.assign(seed_z=med)
    ok = C[np.isfinite(C.seed_z)]
    k = max(1, int(round(a.top_frac * len(ok))))
    Sd = ok.nlargest(k, "seed_z").sort_values("t").reset_index(drop=True)
    Sd[["t", "episode", "seed_z"]].to_parquet(a.out, index=False)
    json.dump(dict(net=a.net, n_windows=int(len(C)), n_scored=int(len(ok)), top_frac=a.top_frac,
                   n_seeds=int(len(Sd)), seed_z_min=float(Sd.seed_z.min()),
                   seed_z_max=float(Sd.seed_z.max()), sites=sites,
                   statistic="median over sites of the per-site 3-C S-window energy z",
                   runtime_s=time.time() - t0),
              open(f"{OUT}/disc_B_seeds.json", "w"), indent=1, default=str)
    print(f"[seed] {len(Sd)} seeds (top {a.top_frac:.1%} of {len(ok)}), "
          f"seed_z {Sd.seed_z.min():.2f}..{Sd.seed_z.max():.2f} -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
