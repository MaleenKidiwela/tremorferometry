#!/usr/bin/env python
"""s11_redpy_dvv -- PLAN §2-5 per (family, station, epoch): selective-stacking weights, nodes (day / N=20 / N=50 chunks /
episode), all-pairs graph vs reference-stack method, event bootstrap, waveform-level step injection, scrambled-control nodes,
time-shuffled data. Writes tables/nodes_fam<F>.parquet, tables/halfsplit_fam<F>.csv, tables/ctrl_fam<F>.parquet,
tables/graph_meta_fam<F>.csv.

    python s11_redpy_dvv.py --fam 10 [--workers 6] [--nboot 50]
"""
from __future__ import annotations
import os, sys, argparse, zlib, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R

STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
KINDS = {"day": ("day", None), "N20": ("chunk", 20), "N50": ("chunk", 50), "episode": ("episode", None)}
INJ = (0.001, 0.003)
NBOOT = 50


def nodes_for(g, kind):
    k, N = KINDS[kind]
    return R.make_nodes(g, k, N)


def run_cfg(W, w, g, nodes, tS, rng, meta, boot=0, ref_boot=0):
    """graph + reference on one node set -> DataFrame of nodes"""
    if len(nodes) < 3:
        return None
    gr = R.graph_run(W, w, nodes, tS)
    rf = R.ref_run(W, w, nodes, tS, rng, ref_boot)
    n = len(nodes)
    mpc = np.full(n, np.nan)
    for i in range(n):
        m = (gr["pi"] == i) | (gr["pj"] == i); mpc[i] = np.nanmedian(gr["cc"][m]) if m.any() else np.nan
    sbg = np.full(n, np.nan)
    if boot:
        V = []
        for _ in range(boot):
            nb = [(lab, rng.choice(ix, len(ix))) for lab, ix in nodes]
            V.append(R.graph_run(W, w, nb, tS)["v"])
        sbg = np.nanstd(np.array(V), 0, ddof=1)
    T = pd.DataFrame(dict(node=[l for l, _ in nodes], n=[len(ix) for _, ix in nodes],
                          t_med=[float(np.median(g.t.to_numpy()[ix])) for _, ix in nodes],
                          cluster=[g.cluster.iloc[ix].mode().iloc[0] for _, ix in nodes],
                          v=gr["v"], sv=gr["sv"], sboot=sbg, pair_cc_med=mpc,
                          eps_ref=rf["eps"], eps_ref_raw=rf["eps_raw"], cc_ref=rf["cc"], sb_ref=rf["sb"]))
    T["chi2"], T["scale"], T["npair"] = gr["chi2"], gr["scale"], gr["npair"]
    T["nrail"] = int((np.abs(gr["d"]) >= R.RAIL).sum()); T["pair_cc_all_med"] = float(np.median(gr["cc"]))
    T["trip_rms"] = R.triplet_misclosure(n, gr["pi"], gr["pj"], gr["d"], rng)
    for k, v in meta.items():
        T[k] = v
    return T


def boundary(g):
    tmed = np.median(g.t); cl = g.cluster.iloc[int(np.argmin(np.abs(g.t - tmed)))]
    b = g.date[g.cluster == cl].min()
    if b == g.date.min():                                   # median episode is the first one -> next episode start
        later = g.date[g.cluster != cl]; later = later[later > b]
        b = later.min() if len(later) else b
    return b


def run(args):
    fam, ns, ep, nboot, labcol = args
    t0 = time.time()
    D = R.load_series(fam, ns)
    if ep not in D:
        return None
    s = D[ep]; g, W, tS = s["g"], s["W"], s["tS"]
    L = pd.read_parquet(f"{R.TB}/subfam_labels_fam{fam}.parquet")
    L = L[(L.netsta == ns) & (L.epoch == ep)].set_index("row_id")
    g["sub"] = L.loc[g.row_id, labcol].to_numpy(); g["cc_own"] = L.loc[g.row_id, "cc_own"].to_numpy(); path = L.path.iloc[0] if labcol == "sub" else "B"
    rng = np.random.default_rng(zlib.crc32(f"{fam}|{ns}|{ep}".encode()))
    wts = {"P": g.P.to_numpy(), "cc": np.clip(g.cc_own.to_numpy(), 0, None), "cc2": np.clip(g.cc_own.to_numpy(), 0, None) ** 2}
    for k in ("cc", "cc2"):
        wts[k] = np.where(wts[k] > 0, wts[k], 1e-6)
    base = dict(fam=fam, netsta=ns, epoch=ep, path=path)
    out = []
    subs = [k for k in sorted(g["sub"].unique()) if k > 0 and (g["sub"] == k).sum() >= 20]
    if path == "A" and (g["sub"] == 0).sum() >= 20:
        subs = subs + [0]                                                     # unclustered pool, kept as its own set
    if path == "B" and len(subs) < 2:
        subs = []
    bnd = boundary(g); post = (g.date >= bnd).to_numpy()
    # ---------- natural + injected, all events
    for tag, a in [("nat", 0.0)] + [(f"inj{a*100:.1f}", a) for a in INJ]:
        Wx = W.copy()
        if a:
            Wx[post] = R.inject(W[post], tS, a)
        for kind in KINDS:
            nodes = nodes_for(g, kind)
            for wn in (("P", "cc", "cc2") if kind in ("N50", "episode") or tag == "nat" else ("P",)):
                T = run_cfg(Wx, wts[wn], g, nodes, tS, rng, dict(base, subset="all", kind=kind, weight=wn, tag=tag),
                            boot=nboot if (tag == "nat" and wn == "P") else 0, ref_boot=100 if (tag == "nat" and wn == "P") else 0)
                if T is not None:
                    out.append(T)
        # ---------- sub-family node sets (P weights)
        for k in subs:
            gk = g[g["sub"] == k]
            for kind in KINDS:
                nodes = [(l, gk.index.to_numpy()[ix]) for l, ix in nodes_for(gk.reset_index(drop=True), kind)]
                T = run_cfg(Wx, wts["P"], g, nodes, tS, rng, dict(base, subset=f"sub{k}", kind=kind, weight="P", tag=tag),
                            boot=nboot if (tag == "nat" and kind in ("N50", "episode")) else 0,
                            ref_boot=100 if tag == "nat" else 0)
                if T is not None:
                    out.append(T)
    # ---------- time-shuffled data (5 shuffles; events permuted across the time-ordered slots)
    for r_ in range(5):
        pm = rng.permutation(len(g)); Wp = W[pm]; wp = wts["P"][pm]
        for kind in ("N20", "N50"):
            T = run_cfg(Wp, wp, g, nodes_for(g, kind), tS, rng, dict(base, subset="all", kind=kind, weight="P", tag=f"shuf{r_}"))
            if T is not None:
                out.append(T)
    NODES = pd.concat(out, ignore_index=True); NODES["boundary"] = bnd
    # ---------- scrambled-control nodes vs the data series reference
    ctrl = []
    gc, Wc = s["gc"], s["Wc"]
    Gall = R.wstack(W, wts["P"]); Rf = R.Ref3(Gall, tS)
    if len(gc) >= 3:
        gc = gc.sort_values("t").reset_index(drop=True)
        for kind in ("day", "N20"):
            nodes = nodes_for(gc, kind)
            if not nodes:
                continue
            St = R.node_stacks(Wc, np.ones(len(gc)), nodes)
            e, c = Rf.measure(St)
            sb = []
            for (_, ix) in nodes:
                B = np.array([R.wstack(Wc[b], np.ones(len(b))) for b in (rng.choice(ix, len(ix)) for _ in range(nboot))])
                sb.append(float(np.std(Rf.measure(B)[0], ddof=1)))
            gr = R.graph_run(Wc, np.ones(len(gc)), nodes, tS) if len(nodes) >= 3 else None
            for i, (lab, ix) in enumerate(nodes):
                ctrl.append(dict(base, kind=kind, node=lab, n=len(ix), cluster=gc.cluster.iloc[ix].mode().iloc[0], eps=e[i], cc=c[i],
                                 sboot=sb[i], ctrl_pair_cc_med=float(np.median(gr["cc"])) if gr else np.nan,
                                 v_graph=gr["v"][i] if gr else np.nan, sv_graph=gr["sv"][i] if gr else np.nan))
    CT = pd.DataFrame(ctrl)
    # ---------- half-split precision by weight (busiest episode)
    hs = []
    ce = g.cluster.value_counts(); cl = ce.index[0]; ixe = g.index[g.cluster == cl].to_numpy(); rest = np.setdiff1d(np.arange(len(g)), ixe)
    for wn, w in wts.items():
        Rr = R.Ref3(R.wstack(W[rest], w[rest]), tS)
        for nn in (20, 50, 100):
            if 2 * nn > len(ixe):
                continue
            diffs, rails = [], []
            for _ in range(50):
                pp = rng.permutation(ixe); A_, B_ = pp[:nn], pp[nn:2 * nn]
                e = Rr.measure(np.array([R.wstack(W[A_], w[A_]), R.wstack(W[B_], w[B_])]))[0]
                diffs.append(e[0] - e[1]); rails.append(bool(np.any(np.abs(e) >= R.RAIL)))
            diffs, rails = np.array(diffs), np.array(rails)
            hs.append(dict(base, episode=cl, weight=wn, n=nn, sigma=float(np.std(diffs, ddof=1) / np.sqrt(2)), rail_frac=float(rails.mean()),
                           sigma_unrailed=float(np.std(diffs[~rails], ddof=1) / np.sqrt(2)) if (~rails).sum() > 5 else np.nan))
    print(f"{fam} {ns} {ep}: {len(g)} events, {len(NODES)} node rows, {time.time()-t0:.0f} s", flush=True)
    return NODES, CT, pd.DataFrame(hs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fam", type=int, default=10); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--nboot", type=int, default=NBOOT); ap.add_argument("--only", default="")
    ap.add_argument("--labcol", default="sub", help="sub (pre-registered path) or labB (path-B column)")
    a = ap.parse_args(); R.ensure()
    L = pd.read_parquet(f"{R.TB}/subfam_labels_fam{a.fam}.parquet")
    jobs = [(a.fam, ns, ep, a.nboot, a.labcol) for (ns, ep) in L[["netsta", "epoch"]].drop_duplicates().itertuples(index=False)]
    if a.only:
        jobs = [j for j in jobs if f"{j[1]}_{j[2]}" in a.only.split(",")]
    jobs.sort(key=lambda j: -(L[(L.netsta == j[1]) & (L.epoch == j[2])].shape[0]))
    from multiprocessing import Pool
    with Pool(min(a.workers, 8)) as p:
        res = [r for r in p.imap_unordered(run, jobs) if r is not None]
    sfx = f"_fam{a.fam}" + (f"_{a.only.replace(',', '+')}" if a.only else "") + ("" if a.labcol == "sub" else f"_{a.labcol}")
    pd.concat([r[0] for r in res]).to_parquet(f"{R.TB}/nodes{sfx}.parquet", index=False)
    pd.concat([r[1] for r in res]).to_parquet(f"{R.TB}/ctrl{sfx}.parquet", index=False)
    pd.concat([r[2] for r in res]).to_csv(f"{R.TB}/halfsplit{sfx}.csv", index=False)


if __name__ == "__main__":
    main()
