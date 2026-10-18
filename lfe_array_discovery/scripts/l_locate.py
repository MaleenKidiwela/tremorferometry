#!/usr/bin/env python
"""Locate every grown WA template from its own stack.  Prefix `l_`.

Agent K's validated locator (`k_locate.py`) with three WA-bed substitutions and NOTHING else:
  * the search grid is the Olympic sub-region, not southern Vancouver Island;
  * the station table is the WA one (`l_common.station_table`);
  * the S station terms are the INTRASLAB-EARTHQUAKE corrections of
    `results/x_station_corrections.csv` (x_corrections.py); a station without a usable term gets
    0.0 and is counted in `n_corr_flagged`, never silently imputed.
The picking convention (S at `off + TPL_PRE` after the node time), the minimum-station bar (5),
the minimum stack-peak bar (4 MAD), the demeaned-residual rms fit and the 50-realisation
SHUFFLED-STATION control are imported from `k_locate` unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l_common as L                                                        # noqa: E402
import k_locate as KL                                                       # noqa: E402

LAT = (47.50, 48.80, 0.03)
LON = (-124.80, -122.20, 0.04)
DEP = (20.0, 60.0, 2.0)


def tt_table(stas, lats, lons, deps, corr):
    """(n_node, n_sta) predicted S times -- k_locate.tt_table's computation, with the per-depth
    distance table built once instead of once per station (96 stations here, 15 there)."""
    import s_common as S
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
    import p_vel1d as V
    ST = L.station_table()
    LA, LO, DP = np.meshgrid(lats, lons, deps, indexing="ij")
    LA, LO, DP = LA.ravel(), LO.ravel(), DP.ravel()
    dgrid = np.arange(0.0, 300.1, 0.5)
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
    a = ap.parse_args()
    L.activate()
    import g_common as G
    z = np.load(a.tpl, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    P = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        P.setdefault(int(g), {})[s] = (float(v[0]) + KL.TPL_PRE, float(v[1]))
    D = pd.read_parquet(a.det)
    nd = D.groupby("grp").size()
    ndays = D.groupby("grp").date.nunique()
    stas = sorted({s for g in P for s in P[g]})
    corr, flagged = L.station_corrections()
    L.log(f"{len(stas)} template stations; intraslab S corrections used for "
          f"{len([s for s in stas if s not in flagged])}, 0.0-FLAGGED for "
          f"{len([s for s in stas if s in flagged])}")
    lats = np.arange(*LAT[:2], LAT[2])
    lons = np.arange(*LON[:2], LON[2])
    deps = np.arange(*DEP[:2], DEP[2])
    L.log(f"grid {lats.size}x{lons.size}x{deps.size} = {lats.size*lons.size*deps.size} nodes")
    NODE, TTg = tt_table(stas, lats, lons, deps, corr)
    si = {s: i for i, s in enumerate(stas)}
    rng = np.random.default_rng(3)
    rows = []
    for g in sorted(P):
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
                         n_corr=int(sum(1 for s in ss if s not in flagged))))
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    Lo = R[R.located == True]                                               # noqa: E712
    L.log(f"{len(R)} templates; located {len(Lo)}; control pass {int(Lo.ctrl_pass.sum())}; "
          f"|dz|<=8 km {int((Lo.dz.abs() <= 8).sum())}; "
          f"both + >=3 days {int(((Lo.ctrl_pass) & (Lo.dz.abs() <= 8) & (Lo.n_days >= 3)).sum())}")


if __name__ == "__main__":
    main()
