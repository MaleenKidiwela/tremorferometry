#!/usr/bin/env python
"""s11_ball3_figures -- S11 ballistic round 3 figures (round3/PLAN.md). User rules: points coloured by cc, no error bars,
line graphs; dates in the figures are calendar dates (PT not needed at year scale).

    python s11_ball3_figures.py [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L
import s11_ball2_longterm as LT

FIG = f"{L.ROOT3}/figures"; TBL = f"{L.ROOT3}/tables"
CMAP = plt.get_cmap("viridis"); NORM = Normalize(0.3, 1.0)


def short(s):
    return str(s).replace("CN.", "").replace("PB.", "").replace("UW.", "").replace("TA.", "").replace("US.", "").replace("C8.", "")


def fig_pairs(tag):
    p = f"{L.dirs('WA_NW')[1]}/pairs_decade{tag}.csv"
    if not os.path.exists(p):
        return
    P = pd.read_csv(p); P = P[P.ref.isin(["G1BH", "G2"]) & (P.n_series >= 5) & P.leg_bound2s_pct.notna()]
    order = P.groupby("X").TS_X.mean().sort_values().index.tolist()
    xi = {x: i for i, x in enumerate(order)}
    fig, axs = plt.subplots(2, 1, figsize=(max(10, 0.28 * len(order)), 8), sharex=True)
    for ref, mk in (("G1BH", "o"), ("G2", "s")):
        q = P[P.ref == ref]
        for ax, col in zip(axs, ("decade_ms", "leg_bound2s_pct")):
            sc = ax.scatter([xi[x] for x in q.X], q[col], c=q.cc_med, cmap=CMAP, norm=NORM, marker=mk, s=28, lw=0,
                            label=f"reference {ref}")
    axs[0].axhline(0, color="0.5", lw=0.6); axs[0].set_ylabel("decade contrast of S difference\n(≥2015 − ≤2014, ms)")
    axs[1].axhline(0.03, color="r", lw=0.8, ls="--", label="rule 1 (0.03 %)")
    axs[1].axhline(0.014, color="k", lw=0.8, ls=":", label="SVI B926 leg (0.014 %)")
    axs[1].set_yscale("log"); axs[1].set_ylabel("2σ bound on δVs/Vs\non the station's leg (%)")
    axs[1].set_xticks(range(len(order))); axs[1].set_xticklabels([short(x) for x in order], rotation=90, fontsize=7)
    axs[1].legend(fontsize=7, loc="upper left"); axs[0].legend(fontsize=7)
    fig.colorbar(sc, ax=axs, shrink=0.5, label="median S cc")
    fig.suptitle("WA-NW: decade change of the S-time difference between each station and the borehole sites "
                 "(stations sorted by S travel time)", fontsize=10)
    fig.savefig(f"{FIG}/fig_wa_pairs_decade{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


def series_panel(ax, d, title):
    d = d.copy(); d["yr"] = d.t + 0.5
    for _, g in d.groupby("sid"):
        g = g.sort_values("yr"); ax.plot(g.yr, g.r, color="0.85", lw=0.5, zorder=1)
    sc = ax.scatter(d.yr, d.r, c=d.cc, cmap=CMAP, norm=NORM, s=7, lw=0, zorder=2)
    cy = d.groupby("yr").apply(lambda x: np.average(x.r, weights=1 / x.sig ** 2)).sort_index()
    ax.plot(cy.index, cy.values, color="k", lw=1.6, zorder=3)
    lim = np.nanpercentile(np.abs(d.r), 97) * 1.2 if len(d) else 5
    ax.set_ylim(-lim, lim); ax.axhline(0, color="0.5", lw=0.5); ax.axvline(2015, color="r", lw=0.6, ls=":")
    ax.set_title(title, fontsize=8, loc="left"); ax.set_ylabel("ms", fontsize=8)
    return sc


def fig_wa_years(tag):
    reg = "WA_NW"; TB = L.dirs(reg)[1]
    p = f"{TB}/pair_corrected{tag}.parquet"
    if not os.path.exists(p):
        return
    Dc = pd.read_parquet(p); Dc = Dc[Dc.bin == "YR"]; Lg = L.oriented(Dc)
    st = L.load(reg, "lags_station", tag); s = st[(st.kind == "data") & (st.bin == "YR")].copy()
    s["SmP"] = np.where(s.railS | s.railP3, np.nan, s.lagS - s.lagP3)
    P = pd.read_csv(f"{TB}/pairs_decade{tag}.csv")
    best = P[(P.ref == "G1BH") & (P.n_series >= 10)].sort_values("leg_bound2s_pct").X.head(3).tolist()
    R = pd.read_csv(f"{TB}/longterm{tag}.csv")
    panels = [("dS", "G1BH", "PB.B013"), ("dS", "G1BH", "PB.B943")] + [("dS", "G1BH", x) for x in best if x not in ("PB.B013", "PB.B943")] + \
             [("S-P3", None, "PB.B005"), ("S-P3", None, "PB.B013"), ("S-P3", None, "UW.LRIV")]
    fig, axs = plt.subplots(len(panels), 1, figsize=(12, 2.0 * len(panels)), sharex=True)
    sc = None
    for ax, (obs, ref, X) in zip(axs, panels):
        if obs == "dS":
            F = L.ref_frame(Lg, L.REFSETS[reg][ref], "rc"); F = F[F.X == X]
            d = LT.series_frame(F, "v", "sig", ["fam", "X", "epX", "bin"], "ccX"); name = f"dS {ref}-{X} corrected"
        else:
            d = LT.series_frame(s[s.netsta == X], "SmP", "sigSmP3", ["fam", "netsta", "epoch", "bin"], "ccP3"); name = f"S-P3 {X}"
        if not len(d):
            continue
        row = R[(R["item"] == name) & (R.bin == "YR") & (R.kind == "data")]
        txt = f"  decade {row.decade_ms.iloc[0]:+.1f} ± {row.decade_sig_ms.iloc[0]:.1f} ms (p {row.decade_perm_p.iloc[0]:.2f})" if len(row) and np.isfinite(row.decade_ms.iloc[0]) else ""
        lab = f"S({short(X)}) − S(B005/6/7 mean)" if obs == "dS" else f"S − P at {short(X)}"
        sc = series_panel(ax, d, lab + txt)
    if sc is not None:
        fig.colorbar(sc, ax=axs, shrink=0.25, label="cc")
    fig.suptitle("WA-NW year stacks: grey = one family (demeaned), black = family-common; red dotted = 2015 split", fontsize=10)
    fig.savefig(f"{FIG}/fig_wa_longterm_years{tag}.png", dpi=105, bbox_inches="tight"); plt.close(fig)


def fig_phase(tag):
    p = f"{TBL}/phase_curves{tag}.csv"
    if not os.path.exists(p):
        return
    Cv = pd.read_csv(p); K = pd.read_csv(f"{TBL}/phase_contrasts{tag}.csv")
    dists = np.load(f"{TBL}/phase_perm_dists{tag}.npz") if os.path.exists(f"{TBL}/phase_perm_dists{tag}.npz") else {}
    regs = [r for r in ("SVI", "WA_NW") if (Cv.region == r).any()]
    fig, axs = plt.subplots(len(regs), 4, figsize=(17, 3.8 * len(regs)), squeeze=False)
    for i, reg in enumerate(regs):
        C = Cv[Cv.region == reg]
        ax = axs[i, 0]
        for kind, ls in (("data", "-"), ("ctrl", "--")):
            q = C[(C.obs == "S-P3") & (C.scheme == "PH") & (C.unit == "constVpVs") & (C.kind == kind)].set_index("lab").reindex(["BTW", "ARR", "AFT"])
            ax.plot(range(3), q.value, ls, marker="o", color="k" if kind == "data" else "0.6",
                    label="data" if kind == "data" else "scrambled control")
        ax.set_xticks(range(3)); ax.set_xticklabels(["between ETS\n(bursts)", "front arrival\n(≤1 d)", "behind front\n(>1 d)"], fontsize=8)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_ylabel("δ(S−P)/(S−P) → dv/v at constant Vp/Vs (%)"); ax.legend(fontsize=7)
        ax.set_title(f"{reg}: region-wide S−P through the ETS cycle", fontsize=9, loc="left")
        ax = axs[i, 1]
        for scheme, labs, col in (("PH5", ["BPRE", "BMID", "BPOST", "EARLY", "ARR0", "AFT1", "AFT3"], "C0"), ("PHE", ["BTW", "E1", "E2", "E3"], "C1")):
            q = C[(C.obs == "S-P3") & (C.scheme == scheme) & (C.unit == "constVpVs") & (C.kind == "data")].set_index("lab").reindex(labs)
            ax.plot(range(len(labs)), q.value, marker="o", color=col, label=scheme)
            for j, l in enumerate(labs):
                ax.annotate(l, (j, q.value.iloc[j]) if np.isfinite(q.value.iloc[j]) else (j, 0), fontsize=6, color=col)
        ax.axhline(0, color="0.5", lw=0.5); ax.legend(fontsize=7); ax.set_title("finer phases (columns)", fontsize=9, loc="left")
        ax = axs[i, 2]
        for (obs, ref), g in C[(C.obs == "dS") & (C.scheme == "PH") & (C.kind == "data")].groupby(["obs", "ref"]):
            q = g.set_index("lab").reindex(["BTW", "ARR", "AFT"]); ax.plot(range(3), q.value, marker="o", label=f"S(X) − S({ref}), all X")
        ax.set_xticks(range(3)); ax.set_xticklabels(["between", "arrival", "behind"], fontsize=8)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_ylabel("ms"); ax.legend(fontsize=7)
        ax.set_title("inter-station S difference", fontsize=9, loc="left")
        ax = axs[i, 3]
        for c, col in (("C1", "C0"), ("C2", "C3")):
            key = f"{reg}|S-P3|-|constVpVs|{c}|date"
            if key in dists:
                x = dists[key]; ax.hist(x, bins=30, histtype="step", color=col, label=f"{c} date-permuted")
                ob = K[(K.region == reg) & (K.obs == "S-P3") & (K.X == "all") & (K.unit == "constVpVs") & (K.kind == "data")][c]
                if len(ob):
                    ax.axvline(ob.iloc[0], color=col, lw=1.5)
        ax.set_xlabel("pooled contrast (%)"); ax.legend(fontsize=7)
        ax.set_title("episode-date permutation (C1 = arrival − between, C2 = behind − arrival)", fontsize=8, loc="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_phase_composite{tag}.png", dpi=105); plt.close(fig)


def fig_common(tag):
    p = f"{TBL}/common_year_curves{tag}.csv"
    if not os.path.exists(p):
        return
    Cv = pd.read_csv(p); Mj = pd.read_csv(f"{TBL}/common_majors{tag}.csv")
    Sp = pd.read_csv(f"{TBL}/common_majors_split{tag}.csv") if os.path.exists(f"{TBL}/common_majors_split{tag}.csv") else pd.DataFrame()
    regs = [r for r in ("SVI", "WA_NW") if (Cv.region == r).any()]
    fig, axs = plt.subplots(len(regs), 3, figsize=(18, 3.8 * len(regs)), squeeze=False)
    for i, reg in enumerate(regs):
        C = Cv[Cv.region == reg].copy(); C["yr"] = C.lab.astype(str).str[:4].astype(int) + 0.5
        ax = axs[i, 0]; q = C[(C.obs == "S-P3") & (C.unit == "constVpVs")].sort_values("yr")
        ax.plot(q.yr, q.value, color="k", lw=1); sc = ax.scatter(q.yr, q.value, c=q.cc, cmap=CMAP, norm=NORM, s=25, zorder=3)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_ylabel("dv/v at constant Vp/Vs (%)")
        ax.set_title(f"{reg}: region-wide S−P per year (all stations, all families)", fontsize=9, loc="left")
        ax = axs[i, 1]
        for obs, col in (("absS", "C0"), ("absP3", "C3"), ("S-P3", "k")):
            q = C[(C.obs == obs) & (C.unit == "ms")].sort_values("yr")
            ax.plot(q.yr, q.value, color=col, lw=1, label={"absS": "S lag vs catalogue time", "absP3": "P lag vs catalogue time", "S-P3": "S − P"}[obs])
            ax.scatter(q.yr, q.value, c=q.cc, cmap=CMAP, norm=NORM, s=14, zorder=3)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_ylabel("ms"); ax.legend(fontsize=7)
        ax.set_title("catalogue timing (moves S and P alike) vs S − P", fontsize=9, loc="left")
        ax = axs[i, 2]; m = Mj[(Mj.region == reg) & (Mj.kind == "data")].sort_values("date"); m["d"] = pd.to_datetime(m.date)
        ax.plot(m.d, m.dvv_constVpVs_pct, color="0.4", lw=0.8, label="all stations")
        ax.scatter(m.d, m.dvv_constVpVs_pct, c=m.cc, cmap=CMAP, norm=NORM, s=18, zorder=3)
        if len(Sp):
            for half, col, lab in (("A", "C0", "borehole/co-located set"), ("B", "C1", "other stations")):
                h = Sp[(Sp.region == reg) & (Sp.kind == "data") & (Sp.half == half) & (Sp.split == "station")].merge(m[["episode", "d"]], on="episode").sort_values("d")
                ax.plot(h.d, h.value, color=col, lw=0.9, label=lab)
        ax.axhline(0, color="0.5", lw=0.5); ax.legend(fontsize=7); ax.set_ylabel("dv/v at constant Vp/Vs (%)")
        ax.set_title("per major ETS (demeaned over majors)", fontsize=9, loc="left")
    fig.colorbar(sc, ax=axs, shrink=0.4, label="median P cc")
    fig.savefig(f"{FIG}/fig_common_years{tag}.png", dpi=105, bbox_inches="tight"); plt.close(fig)


def fig_jointp(tag):
    p = f"{TBL}/jointp_quality{tag}.csv"
    if not os.path.exists(p):
        return
    R = pd.read_csv(p); R = R[R.kind == "data"]
    fig, ax = plt.subplots(figsize=(9, 4.5)); bins = ["EP", "YR", "Y2", "PH", "PH5", "PHE"]
    for (reg, grp), g in R.groupby(["region", "group"]):
        g = g.set_index("bin").reindex(bins)
        l, = ax.plot(range(len(bins)), g.joint_sigP3, marker="o", label=f"{reg} {grp}: site-joint P")
        ax.plot(range(len(bins)), g.best_member_sigP3_same, marker="s", ls="--", color=l.get_color(), label=f"{reg} {grp}: best single instrument")
    ax.set_xticks(range(len(bins))); ax.set_xticklabels(bins); ax.set_ylabel("median bootstrap σ of the P lag (ms)")
    ax.legend(fontsize=7); ax.set_title("Joint multi-instrument P alignment vs the best single instrument", fontsize=9)
    fig.savefig(f"{FIG}/fig_jointp{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    for f in (fig_pairs, fig_wa_years, fig_phase, fig_common, fig_jointp):
        try:
            f(a.tag); print("ok", f.__name__)
        except Exception as e:
            import traceback; print("FIG ERROR", f.__name__, traceback.format_exc()[-600:])


if __name__ == "__main__":
    main()
