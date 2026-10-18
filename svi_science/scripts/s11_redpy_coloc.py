#!/usr/bin/env python
"""s11_redpy_coloc -- PLAN §5 T6: co-located PGC site check. Chunk nodes (N=20/50) built on the events shared by
CN.PGC E1 and PB.B009/B010/B011 (2005-11-02 .. 2015-01-12), graph-solved at each station on its own windows; Pearson rho
between stations for identical event sets, permutation p. Also with PGC-E1 sub-family labels (each sub-family own zero).

    python s11_redpy_coloc.py --fam 10
"""
from __future__ import annotations
import os, sys, argparse, itertools
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R

STA = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]


def load(args):
    fam, ns = args
    s = R.load_series(fam, ns)["E1"]
    return ns, s["g"][["row_id", "t", "date", "cluster", "P"]], s["W"], s["tS"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fam", type=int, default=10); a = ap.parse_args()
    from multiprocessing import Pool
    with Pool(4) as p:
        dat = {ns: (g, W, tS) for ns, g, W, tS in p.map(load, [(a.fam, ns) for ns in STA])}
    shared = set(dat[STA[0]][0].row_id)
    for ns in STA[1:]:
        shared &= set(dat[ns][0].row_id)
    g0 = dat["CN.PGC"][0]; g0 = g0[g0.row_id.isin(shared) & (g0.date < "2015-01-12")].sort_values("t").reset_index(drop=True)
    L = pd.read_parquet(f"{R.TB}/subfam_labels_fam{a.fam}.parquet")
    lab = L[(L.netsta == "CN.PGC") & (L.epoch == "E1")].set_index("row_id")["sub"]
    g0["sub"] = lab.loc[g0.row_id].to_numpy()
    rng = np.random.default_rng(11)
    vals, res = [], []
    for kind, N in (("N20", 20), ("N50", 50), ("episode", None)):
        sets = [("all", g0)] + [(f"sub{k}", g0[g0["sub"] == k]) for k in sorted(g0["sub"].unique()) if k > 0 and (g0["sub"] == k).sum() >= 40]
        for sname, gs in sets:
            gs = gs.reset_index(drop=True)
            nodes0 = R.make_nodes(gs, "chunk" if N else "episode", N)
            if len(nodes0) < 4:
                continue
            for ns in STA:
                g, W, tS = dat[ns]; pos = dict(zip(g.row_id, range(len(g))))
                nodes = [(l, np.array([pos[r] for r in gs.row_id.iloc[ix]])) for l, ix in nodes0]
                gr = R.graph_run(W, g.P.to_numpy(), nodes, tS)
                for i, (l, ix) in enumerate(nodes0):
                    vals.append(dict(kind=kind, subset=sname, netsta=ns, node=l, cluster=gs.cluster.iloc[ix].mode().iloc[0],
                                     n=len(ix), v=gr["v"][i], sv=gr["sv"][i], chi2=gr["chi2"]))
    V = pd.DataFrame(vals); V.to_parquet(f"{R.TB}/coloc_nodes_fam{a.fam}.parquet", index=False)
    for (kind, sname), q in V.groupby(["kind", "subset"]):
        P = q.pivot_table(index="node", columns="netsta", values="v")
        for s1, s2 in itertools.combinations(STA, 2):
            x = P[[s1, s2]].dropna()
            if len(x) < 4:
                continue
            r0 = float(np.corrcoef(x[s1], x[s2])[0, 1])
            perm = np.array([np.corrcoef(rng.permutation(x[s1].to_numpy()), x[s2])[0, 1] for _ in range(2000)])
            res.append(dict(kind=kind, subset=sname, sta1=s1, sta2=s2, n_nodes=len(x), rho=r0,
                            p=float((np.sum(perm >= r0) + 1) / 2001), sd1=float(x[s1].std()), sd2=float(x[s2].std()),
                            sd_diff=float((x[s1] - x[s2]).std())))
    RES = pd.DataFrame(res); RES.to_csv(f"{R.TB}/coloc_rho_fam{a.fam}.csv", index=False)
    print(RES.round(3).to_string())


if __name__ == "__main__":
    main()
