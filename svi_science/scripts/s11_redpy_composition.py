#!/usr/bin/env python
"""s11_redpy_composition -- post-hoc diagnostic (logged as such in RESULTS): what drives the within-episode chunk-to-chunk
dv/v that co-located PGC/B009/B010/B011 share on identical events? Candidates per chunk: mean S realign lag (site-averaged
over the 4 instruments = source-position proxy), share of the high-cc sub-family (SNR), mean detector P, mean single-event
cc with the grand stack. Correlated with the chunk graph values (within-episode residuals) from coloc_nodes_fam<F>.parquet.

    python s11_redpy_composition.py --fam 10
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R
import s10_measure as M
import s10_common as S

STA = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]


def lags(args):
    fam, ns = args
    rows = pd.read_parquet(f"{R.PIL}/tables/rows.parquet")
    rs = pd.read_parquet(f"{R.PIL}/tables/row_station.parquet"); rs = rs[(rs.netsta == ns) & (rs.epoch == "E1")][["row_id"]]
    ser = pd.read_csv(f"{R.PIL}/tables/series.csv"); tS = float(ser[(ser.fam == fam) & (ser.netsta == ns) & (ser.epoch == "E1")].tS_pick.iloc[0])
    z = np.load(f"{R.PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
    d = rows[(rows.fam == fam) & (rows.series == "data") & (rows.P >= 0.5)].merge(rs, on="row_id"); d = d[d.row_id.isin(pos)]
    W = z["W"][[pos[r] for r in d.row_id]].astype(np.float64)
    ok = np.isfinite(W[:, 0, 0]); d, W = d[ok].reset_index(drop=True), W[ok]
    G = R.wstack(W, d.P.to_numpy())
    _, lag, cc = M.realign(W, G, tS)
    F = R.feat(M.realign(W, G, tS)[0].astype(np.float64), tS); Gs = F.sum(0)
    ccg = (F @ Gs - (F * F).sum(1)) / np.linalg.norm(Gs[None] - F, axis=1)
    return ns, pd.DataFrame(dict(row_id=d.row_id, lag=lag, realign_cc=cc, ccg=ccg))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fam", type=int, default=10); a = ap.parse_args()
    from multiprocessing import Pool
    with Pool(4) as p:
        res = dict(p.map(lags, [(a.fam, ns) for ns in STA]))
    shared = set.intersection(*[set(v.row_id) for v in res.values()])
    rows = pd.read_parquet(f"{R.PIL}/tables/rows.parquet").set_index("row_id")
    L = pd.read_parquet(f"{R.TB}/subfam_labels_fam{a.fam}.parquet"); lab = L[(L.netsta == "CN.PGC") & (L.epoch == "E1")].set_index("row_id")
    g0 = lab.loc[lab.index.isin(shared)].reset_index()[["row_id", "t", "date", "cluster", "P", "labA"]]
    g0 = g0[g0.date < "2015-01-12"].sort_values("t").reset_index(drop=True)
    for ns in STA:
        q = res[ns].set_index("row_id"); g0[f"lag_{ns}"] = q.loc[g0.row_id, "lag"].to_numpy(); g0[f"ccg_{ns}"] = q.loc[g0.row_id, "ccg"].to_numpy()
    g0["lag_site"] = g0[[f"lag_{ns}" for ns in STA]].mean(1)
    g0["ccg_site"] = g0[[f"ccg_{ns}" for ns in STA]].mean(1)
    V = pd.read_parquet(f"{R.TB}/coloc_nodes_fam{a.fam}.parquet")
    out = []
    for kind, N in (("N20", 20), ("N50", 50)):
        nodes = R.make_nodes(g0, "chunk", N)
        feats = pd.DataFrame(dict(node=[l for l, _ in nodes], lag_site=[g0.lag_site.iloc[ix].mean() for _, ix in nodes],
                                  sub1_share=[(g0.labA.iloc[ix] == 1).mean() for _, ix in nodes], P_mean=[g0.P.iloc[ix].mean() for _, ix in nodes],
                                  ccg_site=[g0.ccg_site.iloc[ix].mean() for _, ix in nodes],
                                  cluster=[g0.cluster.iloc[ix].mode().iloc[0] for _, ix in nodes]))
        for c in ["lag_site", "sub1_share", "P_mean", "ccg_site"]:
            feats[c + "_res"] = feats[c] - feats.groupby("cluster")[c].transform("mean")
        q = V[(V.kind == kind) & (V.subset == "all")]
        for ns in STA:
            x = q[q.netsta == ns].set_index("node")[["v", "cluster"]].join(feats.set_index("node"), rsuffix="_f").dropna(subset=["v"])
            x["v_res"] = x.v - x.groupby("cluster").v.transform("mean")
            for c in ["lag_site", "sub1_share", "P_mean", "ccg_site"]:
                for lvl, yv, xv in (("within-episode", "v_res", c + "_res"), ("raw", "v", c)):
                    r = float(np.corrcoef(x[yv], x[xv])[0, 1]); slope = float(np.polyfit(x[xv], x[yv], 1)[0])
                    out.append(dict(kind=kind, netsta=ns, feature=c, level=lvl, n=len(x), rho=r, slope_pct_per_unit=-100 * slope))
    O = pd.DataFrame(out); O.to_csv(f"{R.TB}/composition_drivers_fam{a.fam}.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(O.round(3).to_string())
        print(g0[["lag_site"] + [f"lag_{ns}" for ns in STA]].corr().round(3))


if __name__ == "__main__":
    main()
