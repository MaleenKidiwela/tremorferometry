#!/usr/bin/env python
"""Locate every grown template from its own stack, and classify the unmatched ones.

Picks: each template is cut TPL_PRE (2 s) before the maximum of its station's stack envelope, so
the S arrival of the template at station s sits at `off_s + TPL_PRE` after the group's node time.
Those are relative S times; the origin time is a free parameter.

Grid search over lat/lon/depth on the Crosson 1-D S travel times (p_vel1d) with the a5 station
corrections where one exists, minimising the rms of the demeaned residual; a coarse pass then a
fine pass.  Control: the same picks with the STATION LABELS SHUFFLED, 50 realisations - a real
moveout must fit far better than a shuffled one.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ks_common as K                                                        # noqa: E402

TPL_PRE = 2.0
MIN_STA = 5
MIN_PEAK = 4.0
N_SHUF = 50


def tt_table(stas, lats, lons, deps):
    """(n_node, n_sta) predicted S times on the grid."""
    import s_common as S
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
    import p_vel1d as V
    ST = K.station_table()
    try:
        CR = S.station_corr(stas)
    except Exception:
        CR = {s: 0.0 for s in stas}
    LA, LO, DP = np.meshgrid(lats, lons, deps, indexing="ij")
    LA, LO, DP = LA.ravel(), LO.ravel(), DP.ravel()
    dgrid = np.arange(0.0, 200.1, 0.5)
    out = np.zeros((LA.size, len(stas)), np.float32)
    for k, s in enumerate(stas):
        d = S.epi_km(LA, LO, float(ST.lat[s]), float(ST.lon[s]))
        for j, dp in enumerate(np.unique(DP)):
            m = DP == dp
            tab = np.array([V.tt(float(x), float(dp), "S") for x in dgrid])
            out[m, k] = np.interp(d[m], dgrid, tab) + float(CR.get(s, 0.0) or 0.0)
    return np.column_stack([LA, LO, DP]), out


def fit(obs, TTg, idx):
    """rms of the demeaned residual at every grid node; returns (best node, rms, rms_all)."""
    R = obs[None, :] - TTg[:, idx]
    R = R - R.mean(1, keepdims=True)
    r = np.sqrt((R ** 2).mean(1))
    i = int(np.argmin(r))
    return i, float(r[i]), r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--det", required=True)
    ap.add_argument("--score", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    z = np.load(a.tpl, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    P = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        P.setdefault(int(g), {})[s] = (float(v[0]) + TPL_PRE, float(v[1]))
    D = pd.read_parquet(a.det)
    nd = D.groupby("grp").size()
    ndays = D.groupby("grp").date.nunique()
    stas = sorted({s for g in P for s in P[g]})
    lats = np.arange(47.9, 49.31, 0.03)
    lons = np.arange(-124.9, -122.79, 0.04)
    deps = np.arange(18.0, 56.1, 2.0)
    K.log(f"grid {lats.size}x{lons.size}x{deps.size} = {lats.size*lons.size*deps.size} nodes, "
          f"{len(stas)} stations")
    NODE, TTg = tt_table(stas, lats, lons, deps)
    si = {s: i for i, s in enumerate(stas)}
    sl = None
    try:
        sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
        import p_common as P_
        sl = P_.slab2_interp()
    except Exception as e:
        K.log(f"slab2 unavailable: {e}")
    rng = np.random.default_rng(3)
    rows = []
    for g in sorted(P):
        ss = [s for s, (o, h) in P[g].items() if h >= MIN_PEAK and s in si]
        if len(ss) < MIN_STA:
            rows.append(dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)),
                             n_days=int(ndays.get(g, 0)), located=False))
            continue
        obs = np.array([P[g][s][0] for s in ss])
        idx = np.array([si[s] for s in ss])
        i, r, _ = fit(obs, TTg, idx)
        rs = []
        for _ in range(N_SHUF):
            pi = rng.permutation(len(idx))
            rs.append(fit(obs, TTg, idx[pi])[1])
        rs = np.array(rs)
        la, lo, dp = NODE[i]
        sz = np.nan
        if sl is not None:
            try:
                sz = float(np.asarray(sl(np.array([la]), np.array([lo]))).ravel()[0])
            except Exception:
                pass
        rows.append(dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)),
                         n_days=int(ndays.get(g, 0)), located=True,
                         lat=float(la), lon=float(lo), dep=float(dp), rms=r,
                         rms_shuf_med=float(np.median(rs)), rms_shuf_p05=float(np.quantile(rs, .05)),
                         ctrl_pass=bool(r < np.quantile(rs, .05)),
                         slab2=sz, dz=(float(dp) - sz) if np.isfinite(sz) else np.nan))
    R = pd.DataFrame(rows)
    if a.score and os.path.exists(a.score):
        S = pd.read_csv(a.score).set_index("grp")
        R["bkey"] = [int(S.bkey.get(g, -1)) if S.frac.get(g, 0) >= 0.30 else -1 for g in R.grp]
        R["frac"] = [float(S.frac.get(g, np.nan)) for g in R.grp]
    R.to_csv(a.out, index=False)
    L = R[R.located == True]                                                # noqa: E712
    K.log(f"{len(R)} templates; located {len(L)}; control pass {int(L.ctrl_pass.sum())}; "
          f"|dz|<=8 km {int((L.dz.abs() <= 8).sum())}")
    if "bkey" in R.columns:
        un = L[(L.bkey < 0)]
        cand = un[(un.ctrl_pass) & (un.dz.abs() <= 8) & (un.n_days >= 3)]
        K.log(f"unmatched-to-kept located templates {len(un)}; "
              f"CANDIDATE (control pass, |dz|<=8 km, >=3 days): {len(cand)}")


if __name__ == "__main__":
    main()
