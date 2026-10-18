#!/usr/bin/env python
"""S11 tomo adaptation (PLAN.md Amendment 2): stack the stacks of families within R km, aligned on S, re-pick.

  python s11_tomo_cluster.py --set B --radius 3
writes tables/picks_<set>cl<R>S1.csv and ..S2.csv (one 'family' per cluster), tables/clusters_<set>_<R>.csv,
and the cluster stacks npz (tables/clstacks_<set>_<R>_S{1,2}.npz) for figures.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
os.chdir(f"{ROOT}/lfe_array_discovery")
import a4f_hypo as H  # noqa: E402

T = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
FS = 40.0


def shift(x, d):
    """shift array (3, n) by d seconds (positive = later) with zero fill, integer samples."""
    k = int(round(d * FS))
    y = np.zeros_like(x)
    if k > 0:
        y[:, k:] = x[:, :-k]
    elif k < 0:
        y[:, :k] = x[:, -k:]
    else:
        y = x.copy()
    return y


def twoway(S):
    """robust two-way fit tS(f,s) = tau_s + delta_f on a long table (fam, sta, tS); delta median-centred."""
    d = pd.Series(0.0, index=S.fam.unique())
    for _ in range(30):
        tau = (S.tS - S.fam.map(d)).groupby(S.sta).median()
        dn = (S.tS - S.sta.map(tau)).groupby(S.fam).median()
        dn = dn.reindex(d.index).fillna(0.0)
        dn -= dn.median()
        if np.allclose(dn, d, atol=1e-3):
            d = dn
            break
        d = dn
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--radius", type=float, default=3.0)
    a = ap.parse_args()
    P = pd.read_csv(f"{T}/picks_{a.set}.csv", low_memory=False)
    P["fam"] = P.fam.astype(str)
    F = P.groupby("fam").agg(lat=("lat", "first"), lon=("lon", "first"), depth_km=("depth_km", "first"),
                             n_detections=("n_detections", "first"), era_start=("era_start", "first"),
                             era_end=("era_end", "first"), src=("src", "first"), grp=("grp", "first"))
    F = F[F.lat.notna()]
    xy = np.c_[(F.lon - F.lon.mean()) * 111.19 * np.cos(np.radians(F.lat.mean())), (F.lat - F.lat.mean()) * 111.19]
    F["cl"] = fcluster(linkage(xy, "single"), a.radius, "distance")
    F.to_csv(f"{T}/clusters_{a.set}_{a.radius:g}.csv")
    print(f"[{a.set} R{a.radius:g}] {len(F)} families -> {F.cl.nunique()} clusters; members: "
          f"{F.cl.value_counts().value_counts().sort_index().to_dict()}", flush=True)
    stk = {}
    for src in F.src.unique():
        z = np.load(f"{ROOT}/{src}")
        want = set(int(g) for g in F[F.src == src].grp)
        for k in z.files:
            g, s = k.split("_", 1)
            if int(g) in want:
                stk[(src, int(g), s)] = z[k].astype(np.float64)
    out = {"S1": [], "S2": []}
    stacks = {"S1": {}, "S2": {}}
    conf = P[(P.rS >= 1.2) & P.tS.notna()]
    for cl, M in F.groupby("cl"):
        S = conf[conf.fam.isin(M.index)][["fam", "sta", "tS"]]
        d = twoway(S) if len(S) else pd.Series(0.0, index=M.index)
        d = d.reindex(M.index).fillna(0.0)
        tau = (S.tS - S.fam.map(d)).groupby(S.sta).median() if len(S) else pd.Series(dtype=float)
        stas = sorted({k[2] for k in stk if (k[0], k[1]) in set(zip(M.src, M.grp.astype(int)))})
        for mode in ("S1", "S2"):
            ST = {}
            for s in stas:
                acc, nm = 0, 0
                for f, r in M.iterrows():
                    key = (r.src, int(r.grp), s)
                    if key not in stk:
                        continue
                    x = stk[key]
                    m = np.abs(x).max()
                    if not np.isfinite(x).all() or m <= 0:
                        continue
                    sh = -d[f]
                    if mode == "S2":
                        c = S[(S.fam == f) & (S.sta == s)]
                        if len(c) and s in tau.index:
                            sh = tau[s] - float(c.tS.iloc[0])
                    acc = acc + shift(x / m, sh)
                    nm += 1
                if nm:
                    ST[s] = acc / nm
            stacks[mode][cl] = ST
            Q = H.make_picks({int(cl): ST})
            if not len(Q):
                continue
            Q["n_members"] = len(M)
            Q["members"] = ";".join(M.index)
            Q["lat"], Q["lon"], Q["depth_km"] = M.lat.mean(), M.lon.mean(), M.depth_km.median()
            Q["n_detections"] = M.n_detections.sum()
            Q["era_start"], Q["era_end"] = M.era_start.min(), M.era_end.max()
            Q["src"] = f"cluster R{a.radius:g} {mode}"
            Q["fam"] = f"cl{cl}"
            Q["delta_sd"] = float(d.std()) if len(d) > 1 else 0.0
            out[mode].append(Q)
    for mode in ("S1", "S2"):
        Q = pd.concat(out[mode], ignore_index=True)
        tag = f"{a.set}cl{a.radius:g}{mode}"
        Q.to_csv(f"{T}/picks_{tag}.csv", index=False)
        np.savez_compressed(f"{T}/clstacks_{a.set}_{a.radius:g}_{mode}.npz",
                            **{f"{c}_{s}": v.astype(np.float32) for c in stacks[mode] for s, v in stacks[mode][c].items()})
        ok = (Q.rP >= 1.2) & (Q.rS >= 1.2)
        print(f"[{tag}] {Q.fam.nunique()} clusters, rows {len(Q)}, rP&rS>=1.2 rows {ok.sum()} "
              f"(clusters with >=6: {(Q[ok].groupby('fam').size() >= 6).sum()})", flush=True)


if __name__ == "__main__":
    main()
