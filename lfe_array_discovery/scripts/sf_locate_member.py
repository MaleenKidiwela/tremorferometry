#!/usr/bin/env python
"""Locate templates that are picked on MEMBER stations (Arm B) instead of on the array SITES.

`n_locate.py` reads its station table from the bed adapter, which on the array bed is the 22-site
table -- so it crashes on Arm B, whose correlation set is "3 stations per array plus the singles"
(41 member stations).  That is why agent O's Arm-B locate stage died with `KeyError: 'B005'` and
why 739 Chestler-chain families were never located.  This module is `n_locate` with the MEMBER
table substituted and nothing else changed: the same wide grid (46.5-48.8 N / -125.2..-121.8 /
15-60 km), the same picking convention (S at `off + 2 s`), the same >= 5-station bar, the same
4-MAD stack-peak bar, the same demeaned-residual rms fit, the same 50-realisation shuffled-station
control, the same intraslab S station terms and the same Slab2 reference.

  python scripts/sf_locate_member.py --tpl ... --det ... --out ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import k_locate as KL                                                       # noqa: E402
import n_locate as NL                                                       # noqa: E402
import sf_common as SF                                                      # noqa: E402


def member_corrections(stas):
    p = f"{SF.RES}/x_station_corrections.csv"
    corr, flagged = {}, []
    C = pd.read_csv(p).set_index("sta") if os.path.exists(p) else pd.DataFrame()
    for s in stas:
        if len(C) and s in C.index and bool(C.loc[s, "usable"]):
            v = float(C.loc[s, "corr_S"])
            corr[s] = v if np.isfinite(v) else 0.0
            if not np.isfinite(v):
                flagged.append(s)
        else:
            corr[s] = 0.0
            flagged.append(s)
    return corr, flagged


def tt_table(stas, lats, lons, deps, corr, ST):
    import s_common as S
    import p_vel1d as V
    LA, LO, DP = np.meshgrid(lats, lons, deps, indexing="ij")
    LA, LO, DP = LA.ravel(), LO.ravel(), DP.ravel()
    dgrid = np.arange(0.0, 400.1, 0.5)
    tabs = {float(dp): np.array([V.tt(float(x), float(dp), "S") for x in dgrid])
            for dp in np.unique(DP)}
    out = np.zeros((LA.size, len(stas)), np.float32)
    for k, s in enumerate(stas):
        d = S.epi_km(LA, LO, float(ST.lat[s]), float(ST.lon[s]))
        for dp, tab in tabs.items():
            m = DP == dp
            out[m, k] = np.interp(d[m], dgrid, tab) + float(corr.get(s, 0.0))
    return np.column_stack([LA, LO, DP]), out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--det", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    import o_common as O
    import g_common as G
    ST = O.member_table()

    z = np.load(a.tpl, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    P = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        P.setdefault(int(g), {})[s] = (float(v[0]) + KL.TPL_PRE, float(v[1]))
    if a.limit:
        keep = sorted(P)[:a.limit]
        P = {g: P[g] for g in keep}
    D = pd.read_parquet(a.det)
    nd = D.groupby("grp").size()
    ndays = D.groupby("grp").date.nunique()
    stas = sorted({s for g in P for s in P[g] if s in ST.index})
    missing = sorted({s for g in P for s in P[g]} - set(stas))
    if missing:
        SF.log(f"WARNING: {len(missing)} template stations absent from the member table: "
               f"{missing[:10]}")
    corr, flagged = member_corrections(stas)
    SF.log(f"{len(stas)} member stations; intraslab S terms for {len(stas) - len(flagged)}, "
           f"0.0-FLAGGED for {len(flagged)}")
    lats, lons, deps = NL.grid()
    SF.log(f"{len(P)} templates; grid {lats.size}x{lons.size}x{deps.size} = "
           f"{lats.size * lons.size * deps.size} nodes")
    NODE, TTg = tt_table(stas, lats, lons, deps, corr, ST)
    si = {s: i for i, s in enumerate(stas)}
    rng = np.random.default_rng(3)
    rows = []
    for n, g in enumerate(sorted(P)):
        ss = [s for s, (o, h) in P[g].items() if h >= KL.MIN_PEAK and s in si]
        if len(ss) < KL.MIN_STA:
            rows.append(dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)),
                             n_days=int(ndays.get(g, 0)), located=False))
            continue
        obs = np.array([P[g][s][0] for s in ss])
        idx = np.array([si[s] for s in ss])
        i, r, _ = KL.fit(obs, TTg, idx)
        rs = np.array([KL.fit(obs, TTg, idx[rng.permutation(len(idx))])[1]
                       for _ in range(KL.N_SHUF)])
        la, lo, dp = NODE[i]
        sz = float(np.asarray(G.slab2_depth(np.array([la]), np.array([lo]))).ravel()[0])
        rows.append(dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)),
                         n_days=int(ndays.get(g, 0)), located=True,
                         lat=float(la), lon=float(lo), dep=float(dp), rms=r,
                         rms_shuf_med=float(np.median(rs)),
                         rms_shuf_p05=float(np.quantile(rs, .05)),
                         ctrl_pass=bool(r < np.quantile(rs, .05)),
                         slab2=sz, dz=(float(dp) - sz) if np.isfinite(sz) else np.nan,
                         edge=NL.is_edge(la, lo, dp, lats, lons, deps),
                         n_corr=int(sum(1 for s in ss if s not in flagged))))
        if (n + 1) % 100 == 0:
            SF.log(f"  located {n + 1}/{len(P)}")
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    Lo = R[R.located == True]                                               # noqa: E712
    ed = Lo.edge.astype(bool) if len(Lo) else pd.Series(dtype=bool)
    SF.log(f"{len(R)} templates; located {len(Lo)}; interior {int((~ed).sum())} "
           f"(on a grid face {int(ed.sum())}); control pass "
           f"{int(Lo.ctrl_pass.sum()) if len(Lo) else 0}; "
           f"depth 15-58 km {int(((Lo.dep >= 15) & (Lo.dep <= 58)).sum()) if len(Lo) else 0}")


if __name__ == "__main__":
    main()
