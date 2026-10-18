#!/usr/bin/env python
"""s11_redpy_summary -- PLAN §5 tests T1-T6 and §6 figures from the s11_redpy_dvv / _coloc / _cluster tables.

    python s11_redpy_summary.py --fam 10
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

KEY = ["netsta", "epoch", "subset", "kind", "weight"]
INJ = {"inj0.1": 0.001, "inj0.3": 0.003}


def wmean(x, s):
    ok = np.isfinite(x)
    if not ok.any():
        return np.nan, np.nan
    w = np.where(np.isfinite(s[ok]) & (s[ok] > 0), 1 / np.maximum(s[ok], 1e-5) ** 2, np.nan)
    if not np.isfinite(w).all():
        return float(np.mean(x[ok])), np.nan
    return float(np.sum(w * x[ok]) / w.sum()), float(1 / np.sqrt(w.sum()))


def add_weights(N):
    """per node: graph value v / sv, reference eps / sb_ref; injected and shuffled runs take the natural run's sigma"""
    nat = N[N.tag == "nat"].set_index(KEY + ["node"])[["sv", "sb_ref", "sboot"]].rename(columns={"sv": "sv_nat", "sb_ref": "sb_nat", "sboot": "sg_nat"})
    N = N.join(nat, on=KEY + ["node"])
    N["sv_nat"] = N.sv_nat.fillna(N.sv); N["sb_nat"] = N.sb_nat.fillna(N.sv_nat)
    return N


def episode_means(q, est):
    col, sig = {"graph": ("v", "sv_nat"), "graph_boot": ("v", "sg_nat"), "reference": ("eps_ref", "sb_nat")}[est]
    out = []
    for cl, e in q.groupby("cluster"):
        m, s = wmean(e[col].to_numpy(), e[sig].to_numpy())
        out.append(dict(cluster=cl, t=float(e.t_med.median()), m=m, s=s, nn=len(e), cc=float(e.pair_cc_med.median() if est != "reference" else e.cc_ref.median())))
    E = pd.DataFrame(out)
    if len(E):
        mm, _ = wmean(E.m.to_numpy(), E.s.to_numpy()); E["m"] -= mm                    # series own zero
    return E


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fam", type=int, default=10); ap.add_argument("--sfx", default=""); a = ap.parse_args(); F = a.fam; X = a.sfx
    N = pd.read_parquet(f"{R.TB}/nodes_fam{F}{a.sfx}.parquet"); N = add_weights(N)
    CT = pd.read_parquet(f"{R.TB}/ctrl_fam{F}{a.sfx}.parquet")
    HS = pd.read_csv(f"{R.TB}/halfsplit_fam{F}{a.sfx}.csv")
    rep = []
    # ---------------- T1 injection
    t1 = []
    for k, q in N[N.tag.isin(["nat", *INJ])].groupby(KEY):
        bnd = pd.Timestamp(q.boundary.iloc[0]).timestamp()
        nat = q[q.tag == "nat"].set_index("node")
        if len(nat) < 3:
            continue
        for tag, amp in INJ.items():
            inj = q[q.tag == tag].set_index("node")
            if len(inj) != len(nat):
                continue
            for est, col, sig in (("graph", "v", "sv_nat"), ("reference", "eps_ref", "sb_nat")):
                post = nat.t_med >= bnd; pre = ~post
                if post.sum() == 0 or pre.sum() == 0:
                    continue
                d = (inj[col] - nat[col]).to_numpy(); s = nat[sig].to_numpy()
                rec = (wmean(d[post.to_numpy()], s[post.to_numpy()])[0] - wmean(d[pre.to_numpy()], s[pre.to_numpy()])[0]) / amp
                cl_pre = nat.cluster[pre]; lastcl = cl_pre[nat.t_med[pre].idxmax()]
                lp = (pre & (nat.cluster == lastcl)).to_numpy(); ep_ = (pre & (nat.cluster != lastcl)).to_numpy()
                leak = (wmean(d[lp], s[lp])[0] - wmean(d[ep_], s[ep_])[0]) / amp if lp.any() and ep_.any() else np.nan
                t1.append(dict(zip(KEY, k), amp=amp, est=est, recovery=rec, leakage=leak, n_nodes=len(nat)))
    T1 = pd.DataFrame(t1, columns=KEY + ["amp", "est", "recovery", "leakage", "n_nodes"]); T1["subset_type"] = np.where(T1.subset == "all", "all", "sub-family")
    T1.to_csv(f"{R.TB}/T1_injection_fam{F}{X}.csv", index=False)
    g1 = T1.groupby(["subset_type", "kind", "weight", "est", "amp"]).agg(rec_med=("recovery", "median"), rec_min=("recovery", "min"),
                                                                      rec_max=("recovery", "max"), leak_med=("leakage", "median"),
                                                                      absleak_med=("leakage", lambda x: np.nanmedian(np.abs(x))), n_series=("recovery", "size"),
                                                                      n_pass=("recovery", lambda x: int(((x > 0.9) & (x < 1.1)).sum()))).reset_index()
    g1["pass"] = (g1.rec_med.between(0.9, 1.1)) & (g1.leak_med.abs() < 0.2)
    g1.to_csv(f"{R.TB}/T1_summary_fam{F}{X}.csv", index=False)
    # ---------------- T2 misclosure vs bootstrap
    q = N[(N.tag == "nat") & N.sboot.notna() & (N.sboot > 0)].copy(); q["ratio"] = q.sv / q.sboot
    T2 = q.groupby(["netsta", "epoch", "subset", "kind", "weight"]).agg(ratio_med=("ratio", "median"), sv_med=("sv", "median"), sboot_med=("sboot", "median"),
                                                                       sbref_med=("sb_ref", "median"), chi2=("chi2", "first"), scale=("scale", "first"),
                                                                       trip_rms=("trip_rms", "first"), n=("ratio", "size")).reset_index()
    T2["pass"] = T2.ratio_med.between(0.8, 1.25); T2.to_csv(f"{R.TB}/T2_misclosure_fam{F}{X}.csv", index=False)
    # ---------------- T3 scrambled-control nodes
    t3 = []
    datacc = N[(N.tag == "nat") & (N.subset == "all") & (N.kind == "N20") & (N.weight == "P")].groupby(["netsta", "epoch"]).pair_cc_all_med.first()
    for (ns, ep, kind), c in CT.groupby(["netsta", "epoch", "kind"]):
        rail = c.eps.abs() >= R.RAIL; nr = c[~rail]
        z = (nr.eps.abs() / nr.sboot)
        t3.append(dict(netsta=ns, epoch=ep, kind=kind, n_nodes=len(c), rail_frac=float(rail.mean()), n_nonrail=len(nr),
                       frac_within2=float((z < 2).mean()) if len(nr) else np.nan, cc_med=float(c.cc.median()),
                       ctrl_pair_cc_med=float(c.ctrl_pair_cc_med.iloc[0]), data_pair_cc_med=float(datacc.get((ns, ep), np.nan))))
    T3 = pd.DataFrame(t3); T3["cc_gap"] = T3.data_pair_cc_med - T3.ctrl_pair_cc_med
    T3["pass"] = (T3.frac_within2 >= 0.9) & (T3.cc_gap >= 0.2); T3.to_csv(f"{R.TB}/T3_control_fam{F}{X}.csv", index=False)
    # ---------------- T3b time-shuffled data
    t3b = []
    for (ns, ep, kind, tag), q in N[N.tag.str.startswith("shuf")].groupby(["netsta", "epoch", "kind", "tag"]):
        nb = N[(N.tag == "nat") & (N.netsta == ns) & (N.epoch == ep) & (N.kind == kind) & (N.subset == "all") & (N.weight == "P")].set_index("node").sboot
        E = episode_means(q.assign(sv_nat=q.sv), "graph"); E = E[np.isfinite(E.s)]
        Eb = episode_means(q.assign(sg_nat=q.node.map(nb)), "graph_boot"); Eb = Eb[np.isfinite(Eb.s)]
        if len(E) >= 3:
            t3b.append(dict(netsta=ns, epoch=ep, kind=kind, tag=tag, n_ep=len(E), chi2_red=float(np.sum((E.m / E.s) ** 2) / (len(E) - 1)),
                            chi2_red_boot=float(np.sum((Eb.m / Eb.s) ** 2) / (len(Eb) - 1)) if len(Eb) >= 3 else np.nan, ep_sd=float(E.m.std())))
    T3b = pd.DataFrame(t3b, columns=["netsta", "epoch", "kind", "tag", "n_ep", "chi2_red", "chi2_red_boot", "ep_sd"]); T3b.to_csv(f"{R.TB}/T3b_shuffled_fam{F}{X}.csv", index=False)
    T3bs = T3b.groupby(["netsta", "epoch", "kind"]).agg(chi2_med=("chi2_red", "median"), chi2_boot_med=("chi2_red_boot", "median"), ep_sd_shuf=("ep_sd", "median"), n_ep=("n_ep", "first")).reset_index()
    T3bs["pass"] = T3bs.chi2_med.between(0.5, 1.5)
    # ---------------- T4 scatter (natural) and T5 sub-family resolved
    t4, EPS_ALL = [], []
    nat = N[N.tag == "nat"]
    for k, q in nat.groupby(KEY):
        for est in ("graph", "graph_boot", "reference"):
            if est == "graph_boot" and not q.sg_nat.notna().all():
                continue
            E = episode_means(q, est); E[KEY] = list(k); E["est"] = est; EPS_ALL.append(E)
    EA = pd.concat(EPS_ALL, ignore_index=True)
    # sub-family resolved: combine each subset's episode means per episode (each subset its own zero)
    sub = EA[EA.subset != "all"]
    comb = []
    for (ns, ep, kind, w, est, cl), q in sub.groupby(["netsta", "epoch", "kind", "weight", "est", "cluster"]):
        m, s = wmean(q.m.to_numpy(), q.s.to_numpy())
        comb.append(dict(netsta=ns, epoch=ep, kind=kind, weight=w, est=est, cluster=cl, t=q.t.median(), m=m, s=s, nn=q.nn.sum(), cc=q.cc.median(), subset="subres"))
    EA = pd.concat([EA, pd.DataFrame(comb)], ignore_index=True)
    EA["station"] = EA.netsta
    EA.to_parquet(f"{R.TB}/episode_means_fam{F}{X}.parquet", index=False)
    # pooled per station (PGC epochs together; each epoch its own zero already); common episode set = 'episode' kind episodes
    for (st, subset, kind, w, est), q in EA[EA.subset.isin(["all", "subres"])].groupby(["station", "subset", "kind", "weight", "est"]):
        ref_eps = set(EA[(EA.station == st) & (EA.subset == "all") & (EA.kind == "episode") & (EA.weight == "P") & (EA.est == "graph")].cluster)
        q = q[q.cluster.isin(ref_eps) & np.isfinite(q.m)]
        if len(q) < 3:
            continue
        sd = float(q.m.std()); ms2 = float(np.nanmean(q.s ** 2))
        t4.append(dict(station=st, subset=subset, kind=kind, weight=w, est=est, n_ep=len(q), ep_sd=sd, mean_sig=float(np.sqrt(ms2)),
                       excess=float(np.sqrt(max(sd ** 2 - ms2, 0))), chi2_red=float(np.sum((q.m / q.s) ** 2) / (len(q) - 1))))
    T4 = pd.DataFrame(t4)
    rec = g1[(g1.amp == 0.003)].set_index(["subset_type", "kind", "weight", "est"]).rec_med
    T4["recovery"] = [rec.get(("all" if r.subset == "all" else "sub-family", r.kind, r.weight, "reference" if r.est == "reference" else "graph"), np.nan)
                      for r in T4.itertuples()]
    T4["ep_sd_calibrated"] = T4.ep_sd / T4.recovery
    # within-episode day-to-day SD (day nodes)
    wd = []
    for (ns, ep, subset, w), q in nat[(nat.kind == "day")].groupby(["netsta", "epoch", "subset", "weight"]):
        for est, col in (("graph", "v"), ("reference", "eps_ref")):
            sds = q.groupby("cluster")[col].agg(lambda x: x.std() if x.notna().sum() >= 3 else np.nan).dropna()
            wd.append(dict(station=ns, epoch=ep, subset=subset, weight=w, est=est, within_sd_med=float(sds.median()) if len(sds) else np.nan,
                           n_ep=len(sds), day_sd=float(q[col].std()), n_days=int(q[col].notna().sum())))
    WD = pd.DataFrame(wd); WD.to_csv(f"{R.TB}/T4_within_episode_fam{F}{X}.csv", index=False)
    T4.to_csv(f"{R.TB}/T4_scatter_fam{F}{X}.csv", index=False)
    # T5: excess ratio subres vs all, same kind/weight/est
    T5 = T4[T4.subset == "subres"].merge(T4[T4.subset == "all"], on=["station", "kind", "weight", "est"], suffixes=("_sub", "_all"))
    T5["excess_ratio"] = T5.excess_sub / T5.excess_all.replace(0, np.nan); T5["sd_ratio"] = T5.ep_sd_sub / T5.ep_sd_all
    T5.to_csv(f"{R.TB}/T5_subfamily_fam{F}{X}.csv", index=False)
    # composition: share of the high-cc sub-family (path A sub 1) per episode, episode-locked? (V and permutation p)
    Lb = pd.read_parquet(f"{R.TB}/subfam_labels_fam{F}.parquet"); comp = []
    rng = np.random.default_rng(3)
    for (ns, ep), q in Lb.groupby(["netsta", "epoch"]):
        ce = q.cluster.value_counts(); q = q[q.cluster.isin(ce.index[ce >= 20])]
        if q.cluster.nunique() < 2:
            continue
        x = (q.labA == 1).astype(int).to_numpy(); c = q.cluster.to_numpy()
        def V(x):
            t = pd.crosstab(x, c).to_numpy().astype(float); e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
            return float(np.sqrt(((t - e) ** 2 / e).sum() / t.sum()))
        v0 = V(x); perm = np.array([V(rng.permutation(x)) for _ in range(2000)])
        fr = q.groupby("cluster").labA.apply(lambda s: (s == 1).mean())
        comp.append(dict(netsta=ns, epoch=ep, n_ep=len(fr), V_sub1_vs_rest=v0, p=float((np.sum(perm >= v0) + 1) / 2001),
                         share_min=float(fr.min()), share_max=float(fr.max()), share_med=float(fr.median())))
    COMP = pd.DataFrame(comp); COMP.to_csv(f"{R.TB}/composition_sub1_share_fam{F}{X}.csv", index=False)
    print("== composition (sub1 share per episode)\n", COMP.round(3).to_string())
    HSp = HS.groupby(["netsta", "epoch", "weight", "n"]).sigma.first().unstack("weight").reset_index(); HSp.to_csv(f"{R.TB}/halfsplit_pivot_fam{F}{X}.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_rows", 400):
        print("== T1 (all, P)\n", g1[(g1.weight == "P")].round(3).to_string())
        print("== T1 (weights)\n", g1[(g1.weight != "P")].round(3).to_string())
        print("== T2\n", T2[T2.subset == "all"].round(5).to_string())
        print("== T2 sub\n", T2[T2.subset != "all"].groupby(["kind"]).ratio_med.describe().round(3).to_string())
        print("== T3\n", T3.round(3).to_string()); print("== T3b\n", T3bs.round(3).to_string())
        print("== T4\n", T4[T4.weight == "P"].assign(ep_sd=lambda d: d.ep_sd * 100, excess=lambda d: d.excess * 100, mean_sig=lambda d: d.mean_sig * 100).round(3).to_string())
        print("== T4 weights\n", T4[(T4.weight != "P")].assign(ep_sd=lambda d: d.ep_sd * 100, excess=lambda d: d.excess * 100).round(3).to_string())
        print("== within\n", WD[WD.subset == "all"].round(5).to_string())
        print("== T5\n", T5[["station", "kind", "weight", "est", "n_ep_sub", "ep_sd_sub", "ep_sd_all", "excess_sub", "excess_all", "excess_ratio"]].round(4).to_string())
        print("== halfsplit\n", HSp.round(4).to_string())
    figures(F, N, EA, T1, X)


def figures(F, N, EA, T1, X):
    sts = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
    # --- Fig 1: graph episode series (N50 nodes) vs reference, nodes N20 as small dots coloured by cc
    fig, axs = plt.subplots(4, 1, figsize=(14, 13), sharex=True)
    for ax, st in zip(axs, sts):
        q = N[(N.tag == "nat") & (N.netsta == st) & (N.subset == "all") & (N.weight == "P") & (N.kind == "N20")]
        for ep_, qq in q.groupby("epoch"):
            ax.scatter(pd.to_datetime(qq.t_med, unit="s"), -100 * qq.v, c=qq.pair_cc_med, cmap="viridis", vmin=0.2, vmax=0.8, s=8, lw=0, alpha=0.6)
        for est, ls, col in (("graph", "-", "k"), ("reference", "--", "tab:red")):
            E = EA[(EA.station == st) & (EA.subset == "all") & (EA.kind == "N50") & (EA.weight == "P") & (EA.est == est)].sort_values("t")
            for ep_, e in E.groupby("epoch"):
                ax.plot(pd.to_datetime(e.t, unit="s"), -100 * e.m, ls, color=col, lw=1.1, label=f"{est} (N=50 nodes → episode)" if ep_ == E.epoch.iloc[0] else None)
                if est == "graph":
                    sc = ax.scatter(pd.to_datetime(e.t, unit="s"), -100 * e.m, c=e.cc, cmap="viridis", vmin=0.2, vmax=0.8, s=40, edgecolor="k", lw=0.5, zorder=4)
        ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.set_ylim(-0.6, 0.6); ax.set_ylabel("dv/v (%)")
        g = EA[(EA.station == st) & (EA.subset == "all") & (EA.kind == "N50") & (EA.weight == "P")]
        sdg, sdr = (g[g.est == e_].m.std() * 100 for e_ in ("graph", "reference"))
        ax.set_title(f"Family {F} {st}: small dots = N=20 graph nodes, big dots + black line = graph episode means, red dashed = reference-stack method"
                     f"  (episode SD graph {sdg:.3f} %, reference {sdr:.3f} %)", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
    fig.colorbar(sc, ax=axs, shrink=0.4, label="median pair cc (graph)")
    fig.savefig(f"{R.FIG}/fam{F}_graph_vs_reference{X}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
    # --- Fig 2: sub-family series (episode kind, graph), one row per sub-family set, columns stations
    subs = sorted([s for s in N.subset.unique() if s != "all"], key=lambda s: (s == "sub0", s))
    have = set(EA[(EA.kind == "N50") & (EA.est == "graph") & (EA.weight == "P")].groupby("subset").filter(lambda d: len(d) >= 3).subset)
    rows_ = [s_ for s_ in ["all"] + subs + ["subres"] if s_ in have]
    fig, axs = plt.subplots(len(rows_), 4, figsize=(18, 2.0 * len(rows_)), sharex=True, sharey=True)
    for j, st in enumerate(sts):
        for i, ss in enumerate(rows_):
            ax = axs[i, j]; E = EA[(EA.station == st) & (EA.subset == ss) & (EA.kind == "N50") & (EA.weight == "P") & (EA.est == "graph")].sort_values("t")
            for ep_, e in E.groupby("epoch"):
                ax.plot(pd.to_datetime(e.t, unit="s"), -100 * e.m, "-", color="k", lw=0.9)
                sc = ax.scatter(pd.to_datetime(e.t, unit="s"), -100 * e.m, c=e.cc, cmap="viridis", vmin=0.2, vmax=0.8, s=22, edgecolor="k", lw=0.4, zorder=3)
            ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.set_ylim(-0.8, 0.8)
            lab = {"all": "all events", "subres": "sub-family resolved", "sub0": "unclustered pool"}.get(ss, ss)
            ax.set_title(f"{st} — {lab} (episode SD {E.m.std() * 100:.3f} %)" if len(E) else f"{st} — {lab}", fontsize=8)
            if j == 0:
                ax.set_ylabel("dv/v (%)")
    fig.colorbar(sc, ax=axs, shrink=0.3, label="median pair cc")
    fig.suptitle(f"Family {F}: graph dv/v per episode (N=50 nodes), by sub-family ({'path-B labels (column)' if X else 'pre-registered labels'}; each set its own zero)", fontsize=11)
    fig.savefig(f"{R.FIG}/fam{F}_subfamily_series{X}.png", dpi=105, bbox_inches="tight"); plt.close(fig)
    # --- Fig 3: injection recovery
    q = T1[(T1.weight == "P") & (T1.subset_type == "all")]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
    kinds = ["day", "N20", "N50", "episode"]
    for ax, col in zip(axs, ["recovery", "leakage"]):
        for k_, est in enumerate(["graph", "reference"]):
            for a_, amp in enumerate([0.001, 0.003]):
                for i, kd in enumerate(kinds):
                    v = q[(q.est == est) & (q.amp == amp) & (q.kind == kd)][col]
                    x = i + (k_ * 2 + a_ - 1.5) * 0.15
                    ax.plot(np.full(len(v), x), v, "o", ms=4, color=["k", "tab:red"][k_], alpha=0.4 + 0.5 * a_,
                            label=f"{est}, +{amp*100:.1f} %" if i == 0 else None)
        ax.set_xticks(range(4)); ax.set_xticklabels(kinds); ax.set_title(f"T1 {col} per series (all events, P weights)", fontsize=9)
        if col == "recovery":
            ax.axhspan(0.9, 1.1, color="g", alpha=0.1)
        else:
            ax.axhspan(-0.2, 0.2, color="g", alpha=0.1)
        ax.legend(fontsize=7)
    fig.savefig(f"{R.FIG}/fam{F}_injection{X}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
