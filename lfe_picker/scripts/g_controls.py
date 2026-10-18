#!/usr/bin/env python
"""Gate-1 POSITIVE CONTROL table: at each withheld station, does a KNOWN Bostock family at matched N
and comparable distance pass the Gate-1 test?  If not, a proposal that fails there is UNTESTABLE, not
FAIL (Merlin: "a known Bostock family at matched N must pass, else 'untestable'").

Members are Bostock's own 2012 ETS event times (truth, independent of any detector of ours).
Writes results/gates/g_gate1_controls.csv (one row per station x family x N).
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
import g_gates as GG                                             # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--withheld", default="VGZ,SNB,B926")
    ap.add_argument("--n-families", type=int, default=26)
    ap.add_argument("--N", default="20,50,100")
    ap.add_argument("--n-time-null", type=int, default=200)
    ap.add_argument("--n-moveout-null", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--out", default=f"{G.GATES}/g_gate1_controls.csv")
    a = ap.parse_args()
    cfg = GG.GateConfig(withheld=tuple(a.withheld.split(",")),
                        g1_n_time_null=a.n_time_null, g1_n_moveout_null=a.n_moveout_null)
    ctx = GG.GateContext(cfg)
    rng = np.random.default_rng(a.seed)
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv").set_index("bkey")
    gt = pd.read_parquet(f"{G.REP}/results/e5r_gt_bostock2012.parquet")
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    n = gt.groupby("doy").size()
    gt = gt[gt.doy.isin(n[n >= 100].index)]
    Ns = [int(x) for x in a.N.split(",")]
    rows = []
    for sta in cfg.withheld:
        la, lo = float(ctx.S.loc[sta, "lat"]), float(ctx.S.loc[sta, "lon"])
        cnt = gt.groupby("bkey").size()
        elig = [b for b in cnt.index if cnt[b] >= max(Ns)]
        d = {b: float(G.epi_km(F.loc[b, "lat"], F.loc[b, "lon"], la, lo)) for b in elig}
        # spread the control families over the distance range seen at this station
        order = sorted(elig, key=lambda b: d[b])
        pick = [order[i] for i in np.unique(np.linspace(0, len(order) - 1, a.n_families).astype(int))]
        for b in pick:
            t_all = gt[gt.bkey == b].t_origin.to_numpy(float)
            for N in Ns:
                t = np.sort(rng.choice(t_all, size=N, replace=False))
                p = GG.Proposal(pid=f"ctl_{sta}_f{b}_N{N}", times=t, lat=float(F.loc[b, "lat"]),
                                lon=float(F.loc[b, "lon"]), dep=float(F.loc[b, "dep"]))
                t0 = time.time()
                r = GG._gate1_one_station(ctx, p, sta, rng)
                r.update(bkey=int(b), N=int(N), runtime_s=time.time() - t0)
                rows.append(r)
                print(f"  {sta} f{b} d={r['dist_km']:.0f} km N={N}: {r['verdict']} "
                      f"excess {r.get('excess_db', np.nan):.2f} dB ({r['runtime_s']:.1f} s)", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    print(R.groupby(["sta", "N"]).verdict.value_counts().to_string())
    print("->", a.out)


if __name__ == "__main__":
    main()
