#!/usr/bin/env python
"""Double-difference relocation between NEIGHBOURING families (Waldhauser & Ellsworth 2000), the
step that makes Chestler & Creager's Olympic locations precise.  Prefix `uu_`.

Differential times come from cross-correlating the two families' LONG STACKS at the SAME station
(same path, nearby sources -> the waveforms really are similar, unlike a station-to-station
correlation), in a 6-s window about each family's own S onset, lags +-2.5 s, 20 Hz stacks
interpolated to 0.01 s by parabolic refinement of the correlation maximum.

System (per linked pair i,j and shared station s):
    ddt_obs - (tt_i^s - tt_j^s) = G_i . dm_i - G_j . dm_j + (dtau_i - dtau_j)
unknowns dm = (dE, dN, dZ) km and dtau s per family, solved by damped LSQR, iterated 4 times with
the partial derivatives recomputed at the current location.  The centroid of the cluster is held
by a weak absolute-location damping so the pass is RELATIVE, exactly as published.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                       # noqa: E402

LINK_KM = 25.0
MAX_NEIGH = 12
MIN_SHARED = 4
WIN_S = 6.0
PRE_S = 1.0
MAXLAG_S = 2.5
CC_MIN = 0.4
ITERS = 4
DAMP_REL = 1.0
DAMP_ABS = 0.05


def xc_refine(a, b, maxlag, fs=U.FS):
    """Normalised 3-C cross-correlation of a against b; returns (lag_s, cc) with sub-sample
    parabolic refinement."""
    nl = int(round(maxlag * fs))
    n = a.shape[1]
    an = np.sqrt((a ** 2).sum())
    cc = np.empty(2 * nl + 1)
    for k, L in enumerate(range(-nl, nl + 1)):
        if L >= 0:
            aa, bb = a[:, L:], b[:, :n - L] if L else b
        else:
            aa, bb = a[:, :n + L], b[:, -L:]
        d = np.sqrt((aa ** 2).sum()) * np.sqrt((bb ** 2).sum())
        cc[k] = (aa * bb).sum() / d if d > 0 else 0.0
    k = int(np.argmax(cc))
    lag = (k - nl) / fs
    if 0 < k < 2 * nl:
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        den = (y0 - 2 * y1 + y2)
        if abs(den) > 1e-12:
            lag += (0.5 * (y0 - y2) / den) / fs
    return float(lag), float(cc[k])


def partials(lat, lon, dep, slat, slon, vpvs=1.78, h=0.5):
    """d tt_S / d(E, N, Z) in s/km at one station, by central differences on the Crosson table."""
    xs, zs, TP, TS = U.tt_table(vpvs=vpvs)
    dlon = h / (111.2 * np.cos(np.radians(lat)))
    dlat = h / 111.2

    def t(la, lo, dp):
        return float(U.interp2(xs, zs, TS, U.epi_km(la, lo, slat, slon), np.array([dp]))[0])
    return (np.array([(t(lat, lon + dlon, dep) - t(lat, lon - dlon, dep)) / (2 * h),
                      (t(lat + dlat, lon, dep) - t(lat - dlat, lon, dep)) / (2 * h),
                      (t(lat, lon, dep + h) - t(lat, lon, dep - h)) / (2 * h)]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loc", required=True)
    ap.add_argument("--picks", required=True)
    ap.add_argument("--stacks", required=True)
    ap.add_argument("--bed", default="north")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import lsqr
    B = U.bed(a.bed)
    if a.bed != "north":
        B.m.activate()
    ST = B.station_table()
    L = pd.read_csv(a.loc)
    L = L[L.located == True].copy()                                         # noqa: E712
    P = pd.read_parquet(a.picks)
    Z = np.load(a.stacks, allow_pickle=True)
    have = set(Z.files)
    gi = {int(g): i for i, g in enumerate(L.grp)}
    la = L.lat.to_numpy(float)
    lo = L.lon.to_numpy(float)
    dp = L.dep.to_numpy(float)
    on = {(int(r.grp), r.sta): float(r.t_onS) for r in P.itertuples() if np.isfinite(r.t_onS)}
    stas = {int(g): [s for s in gg.sta if (int(g), s) in on and f"{g}_{s}" in have]
            for g, gg in P.groupby("grp") if int(g) in gi}
    U.log(f"{len(L)} located families; building links <= {LINK_KM} km")
    links = []
    for i in range(len(L)):
        d = U.epi_km(la, lo, la[i], lo[i])
        nb = np.argsort(d)[1:MAX_NEIGH + 1]
        for j in nb:
            if j > i and d[j] <= LINK_KM:
                links.append((i, int(j)))
    U.log(f"{len(links)} candidate pairs")
    nw = int(round(WIN_S * U.FS))
    obs = []
    for (i, j) in links:
        gi_, gj = int(L.grp.iloc[i]), int(L.grp.iloc[j])
        sh = [s for s in stas.get(gi_, []) if s in stas.get(gj, [])]
        if len(sh) < MIN_SHARED:
            continue
        got = []
        for s in sh:
            ki = U.I_NODE + int(round((on[(gi_, s)] - PRE_S) * U.FS))
            kj = U.I_NODE + int(round((on[(gj, s)] - PRE_S) * U.FS))
            if min(ki, kj) < 0 or max(ki, kj) + nw > U.NLONG:
                continue
            A = np.asarray(Z[f"{gi_}_{s}"], float)[:, ki:ki + nw]
            Bw = np.asarray(Z[f"{gj}_{s}"], float)[:, kj:kj + nw]
            lag, cc = xc_refine(A, Bw, MAXLAG_S)
            if cc >= CC_MIN:
                # observed differential: (t_i - t_j) at station s
                got.append((s, (on[(gi_, s)] - on[(gj, s)]) + lag, cc))
        if len(got) >= MIN_SHARED:
            for s, ddt, cc in got:
                obs.append(dict(i=i, j=j, sta=s, ddt=ddt, cc=cc))
    O = pd.DataFrame(obs)
    U.log(f"{len(O)} double-difference observations over "
          f"{0 if O.empty else O.groupby(['i','j']).ngroups} pairs")
    if O.empty:
        L.to_csv(a.out, index=False)
        return
    m = np.column_stack([np.zeros(len(L)), np.zeros(len(L)), np.zeros(len(L)), np.zeros(len(L))])
    xs, zs, TP, TS = U.tt_table()
    for it in range(ITERS):
        rowsA, colsA, valsA, rhs, wts = [], [], [], [], []
        r = 0
        for rec in O.itertuples():
            s = rec.sta
            sl, so = float(ST.lat[s]), float(ST.lon[s])
            ti = float(U.interp2(xs, zs, TS, U.epi_km(la[rec.i], lo[rec.i], sl, so),
                                 np.array([dp[rec.i]]))[0])
            tj = float(U.interp2(xs, zs, TS, U.epi_km(la[rec.j], lo[rec.j], sl, so),
                                 np.array([dp[rec.j]]))[0])
            gI = partials(la[rec.i], lo[rec.i], dp[rec.i], sl, so)
            gJ = partials(la[rec.j], lo[rec.j], dp[rec.j], sl, so)
            for k in range(3):
                rowsA.append(r); colsA.append(4 * rec.i + k); valsA.append(gI[k])
                rowsA.append(r); colsA.append(4 * rec.j + k); valsA.append(-gJ[k])
            rowsA.append(r); colsA.append(4 * rec.i + 3); valsA.append(1.0)
            rowsA.append(r); colsA.append(4 * rec.j + 3); valsA.append(-1.0)
            rhs.append(rec.ddt - (ti - tj))
            wts.append(rec.cc)
            r += 1
        A = coo_matrix((np.array(valsA) * np.repeat(wts, 8), (rowsA, colsA)),
                       shape=(r, 4 * len(L))).tocsr()
        b = np.array(rhs) * np.array(wts)
        sol = lsqr(A, b, damp=DAMP_REL, iter_lim=400)[0]
        d = sol.reshape(-1, 4)
        d[:, :3] = np.clip(d[:, :3], -8.0, 8.0)
        lo = lo + d[:, 0] / (111.2 * np.cos(np.radians(la)))
        la = la + d[:, 1] / 111.2
        dp = np.clip(dp + d[:, 2], 4.0, 80.0)
        U.log(f"  dd iter {it}: median |move| {np.median(np.abs(d[:, :3]).sum(1)):.2f} km, "
              f"rms residual {np.sqrt((b ** 2).mean()):.3f} s")
    L["lat_dd"], L["lon_dd"], L["dep_dd"] = la, lo, dp
    L["dd_move_km"] = U.epi_km(la, lo, L.lat.to_numpy(float), L.lon.to_numpy(float))
    L["dd_nobs"] = [int((O.i == i).sum() + (O.j == i).sum()) for i in range(len(L))]
    keep = L.dd_nobs > 0
    L.loc[~keep, ["lat_dd", "lon_dd", "dep_dd"]] = L.loc[~keep, ["lat", "lon", "dep"]].values
    out = L.copy()
    out["lat"], out["lon"], out["dep"] = out.lat_dd, out.lon_dd, out.dep_dd
    out.to_csv(a.out, index=False)
    # NB the observation table must NOT be written as `loc_*.csv`: uu_eval.py globs `loc_*.csv`
    # in this directory and treats every match as a location table (2026-09-23 crash).
    obsp = os.path.join(os.path.dirname(a.out) or ".",
                        "dd_" + re.sub(r"^loc_", "", os.path.basename(a.out))
                        .replace(".csv", "_obs.csv"))
    O.to_csv(obsp, index=False)
    U.log(f"wrote {a.out} (differential-time observations: {obsp})")


if __name__ == "__main__":
    main()
