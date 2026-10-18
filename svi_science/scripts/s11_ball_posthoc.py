#!/usr/bin/env python
"""s11_ball_posthoc -- POST-HOC diagnostics for the S11 ballistic test (labelled post hoc in RESULTS.md):
(1) scatter of demeaned dS in bins of the bootstrap sigma (co-located vs far) -- where does the heavy tail come from;
(2) sensitivity: inject an episode-common term of SD s (ms) into the real demeaned values, rerun the decomposition +
    1,000-draw permutation; detection rate at p < 0.01 -> the smallest common SD this data set would catch;
(3) largest episode-common terms of the co-located pairs (station timing, not path).
"""
from __future__ import annotations
import numpy as np, pandas as pd
import s11_ball_analysis as A

TB = A.TB
rng = np.random.default_rng(20261013)
pr = pd.read_parquet(f"{TB}/lags_pair.parquet"); pr = pr[pr.kind == "data"].copy()
pr.loc[pr.railS, "dS"] = np.nan
st = pd.read_parquet(f"{TB}/lags_station.parquet"); st = st[st.kind == "data"].copy()
st["SmP"] = st.lagS - st.lagP; st.loc[st.railS | st.railP, "SmP"] = np.nan
PK = ["fam", "pair", "epA", "epB"]; SK = ["fam", "netsta", "epoch"]
geo = pd.read_csv(f"{TB}/series_geom.csv"); spmap = geo.groupby(["fam", "netsta"]).SmP_pick.median()

# (1) scatter by sigma bin
bins = [0, 2, 4, 8, 16, 1e9]; rows = []
items = {"co-located (PGC ref)": A.COLOC["CN.PGC"], "co-located (B011 ref)": A.COLOC["PB.B011"]}
items.update({f"far {r}-{f}": [f"{r}-{f}"] for r in ("CN.PGC", "PB.B011") for f in A.FAR})
for lab, pl in items.items():
    d = A.demean(pr[pr.pair.isin(pl)], "dS", PK)
    d["b"] = pd.cut(d.sig_dS, bins)
    for b, g in d.groupby("b", observed=True):
        rows.append(dict(item=lab, sig_bin=str(b), n=len(g), sd_ms=np.sqrt(np.mean(g.r ** 2 * g.K / (g.K - 1))),
                         robust_ms=1.4826 * np.median(np.abs(g.r)), sig_med=g.sig_dS.median()))
B = pd.DataFrame(rows); B.to_csv(f"{TB}/posthoc_scatter_by_sigma.csv", index=False)
print(B.round(2).to_string())


# (2) injection sensitivity
def detect(d, s, ntrial=40, nperm=1000):
    e_n = d.groupby("episode").size(); keep = e_n[e_n >= 3].index
    q = d[d.episode.isin(keep)].reset_index(drop=True)
    g = q.episode.astype("category").cat.codes.to_numpy(); G = g.max() + 1
    sid = q.sid.to_numpy(); order = np.argsort(sid, kind="stable"); sid_s = sid[order]
    hits = 0
    for t in range(ntrial):
        a = rng.normal(0, s, G)
        y0 = q.r.to_numpy() + a[g]
        y = y0 - pd.Series(y0).groupby(sid).transform("mean").to_numpy()      # re-demean per series
        F = A.anova(y, g)[3]; c = 0
        for i in range(nperm):
            perm = order[np.argsort(sid_s + rng.random(len(sid_s)), kind="stable")]
            yp = np.empty_like(y); yp[order] = y[perm]
            c += A.anova(yp, g)[3] >= F
        hits += (1 + c) / (nperm + 1) < 0.01
    return hits / ntrial


S_GRID = [2, 4, 6, 8, 12]
rows = []
for r_ in ("CN.PGC", "PB.B011"):
    for f in A.FAR:
        pp = f"{r_}-{f}"; d = A.demean(pr[pr.pair == pp], "dS", PK)
        TSX = np.nanmean([spmap.get((fa, f), np.nan) * 1.75 / 0.75 for fa in d.fam.unique()])
        for s in S_GRID:
            rows.append(dict(item=f"dS {pp}", s_ms=s, detect_rate=detect(d, s), TS_X_s=TSX,
                             dbb_Xleg_pct=s / 1e3 / TSX * 100)); print(rows[-1], flush=True)
for ns in ["CN.PGC", "PB.B011", "PB.B926", "CN.SNB", "CN.VGZ"]:
    d = A.demean(st[st.netsta == ns], "SmP", SK)
    SP = np.nanmean([spmap.get((fa, ns), np.nan) for fa in d.fam.unique()])
    for s in S_GRID + [16, 24]:
        rows.append(dict(item=f"S-P {ns}", s_ms=s, detect_rate=detect(d, s), SP_s=SP,
                         dbb_constVpVs_pct=s / 1e3 / SP * 100, dbb_VsOnly_pct=s / 1e3 / (SP * 1.75 / 0.75) * 100))
        print(rows[-1], flush=True)
I = pd.DataFrame(rows); I.to_csv(f"{TB}/posthoc_sensitivity.csv", index=False)

# (3) co-located episode-common terms
E = pd.read_csv(f"{TB}/episode_common.csv")
E = E[(E.kind == "data") & (E.subset == "all") & E["item"].str.startswith("floor dS")]
E.reindex(E.common.abs().sort_values(ascending=False).index).head(20).to_csv(f"{TB}/posthoc_floor_top_episodes.csv", index=False)
print(E.reindex(E.common.abs().sort_values(ascending=False).index).head(20).round(2).to_string())
