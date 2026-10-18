#!/usr/bin/env python
"""s14_raw -- S14 step A: the existing LFE dv/v and its stretching CC side by side (raw), ETS majors shaded, and the raw
list of 'CC low while dv/v flat' cells. PLAN §4 step A, §5 C7 (raw form here; REG form in s14_tests.py).
    python s14_raw.py [--pilot]
"""
import os, sys, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import s14_common as C

CMAP, VMIN, VMAX = "viridis", 0.3, 1.0
MK = {"CN.PGC": "o", "PB.B009": "s", "PB.B010": "^", "PB.B011": "D"}
LC = {"CN.PGC": "#1f4e79", "PB.B009": "#c55a11", "PB.B010": "#548235", "PB.B011": "#7030a0"}


def shade(ax, U, t0=None, t1=None):
    for _, r in U[U.cls == "major"].iterrows():
        ax.axvspan(r.start, r.end + pd.Timedelta(days=1), color="#f4b183", alpha=0.35, lw=0, zorder=0)


def rob(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def fig_fam10(U):
    T = pd.read_csv(f"{C.RES}/s10_dvv/pilot/tables/fam10_PGC_3c25_daily.csv", parse_dates=["date"])
    Ep = pd.read_csv(f"{C.RES}/s10_dvv/pilot/tables/fam10_PGC_episode_dvv.csv", parse_dates=["date"])
    fig, axs = plt.subplots(4, 1, figsize=(13, 10), sharex=True)
    for k, (pre, lab) in enumerate([("raw", "raw day stacks"), ("svd", "SVD-Wiener day stacks")]):
        a, b = axs[2 * k], axs[2 * k + 1]
        shade(a, U); shade(b, U)
        sc = a.scatter(T.date, T[f"{pre}_3c"], c=T[f"{pre}_cc3c26"], cmap=CMAP, vmin=VMIN, vmax=VMAX, s=14, lw=0)
        a.set_ylabel("dv/v (%)"); a.set_title(f"family 10 at PGC, {lab}: dv/v (3C 2–5 s, as stored) and its stretching CC", fontsize=10)
        b.plot(T.date, T[f"{pre}_cc3c26"], "-", color="0.6", lw=0.6)
        b.scatter(T.date, T[f"{pre}_cc3c26"], c=T[f"{pre}_cc3c26"], cmap=CMAP, vmin=VMIN, vmax=VMAX, s=14, lw=0)
        b.set_ylabel("stretch CC")
        if k == 0:
            a.scatter(Ep.date, Ep.dvv, c=Ep.cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=60, marker="*", edgecolor="k", lw=0.5, zorder=4,
                      label="episode stack (S10 pilot)")
            a.legend(loc="upper left", fontsize=8)
            b.scatter(Ep.date, Ep.cc, c=Ep.cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=60, marker="*", edgecolor="k", lw=0.5, zorder=4)
    fig.colorbar(sc, ax=axs, shrink=0.5, label="stretch CC")
    fig.savefig(f"{C.FIG}/figA1_fam10_PGC_dvv_cc_raw.png", dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_pilot(E, U):
    fams = C.PILOT_FAMS
    fig, axs = plt.subplots(2 * len(fams), 1, figsize=(14, 4.2 * len(fams)), sharex=True)
    for i, f in enumerate(fams):
        a, b = axs[2 * i], axs[2 * i + 1]
        shade(a, U); shade(b, U)
        for ns in C.PILOT_STAS:
            q = E[(E.fam == f) & (E.netsta == ns)].sort_values("date")
            if not len(q):
                continue
            sc = a.scatter(q.date, q.dvv25, c=q.cc3c_25, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=10 + 0.15 * q.n.clip(upper=300),
                           marker=MK[ns], lw=0.3, edgecolor="0.3", label=ns)
            for ep, qq in q.groupby("epoch"):
                b.plot(qq.date, qq.cc3c_25, "-", color=LC[ns], lw=0.7, alpha=0.8)
            b.scatter(q.date, q.cc3c_25, c=q.cc3c_25, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=12, marker=MK[ns], lw=0)
        a.set_ylabel("dv/v (%)"); b.set_ylabel("stretch CC"); a.set_ylim(-0.9, 0.9); b.set_ylim(-0.1, 1.0)
        a.set_title(f"family {f}: per-episode dv/v (3C 2–5 s, leave-episode-out reference) and its stretching CC; "
                    f"marker = instrument, size = n repeats; shaded = ETS majors", fontsize=9)
        if i == 0:
            a.legend(loc="lower left", fontsize=7, ncol=4)
    fig.colorbar(sc, ax=axs, shrink=0.4, label="stretch CC")
    fig.savefig(f"{C.FIG}/figA2_pilot_episode_dvv_cc_raw.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def pooled(E, U, region):
    Q = E[np.isfinite(E.dvv25) | np.isfinite(E.cc3c_25)].copy()
    Q["dv_a"] = Q.dvv25 - Q.groupby("series").dvv25.transform("median")
    Q["cc_a"] = Q.cc3c_25 - Q.groupby("series").cc3c_25.transform("median")
    Q["c2_a"] = Q.c2_cc3c_25 - Q.groupby("series").c2_cc3c_25.transform("median")
    Q["n_a"] = Q.lnn - Q.groupby("series").lnn.transform("median")
    P = Q.groupby("episode").agg(date=("date", "median"), cls=("cls", "first"), dv=("dv_a", "median"), cc=("cc_a", "median"),
                                 c2=("c2_a", "median"), lnn=("n_a", "median"), ccabs=("cc3c_25", "median"), ncell=("cc_a", "size")).reset_index()
    P = P[P.ncell >= 10].sort_values("date")
    P.to_csv(f"{C.TB}/raw_pooled_{region}.csv", index=False)
    return P


def fig_pooled(Ps, Us):
    fig, axs = plt.subplots(2 * len(Ps), 1, figsize=(14, 4 * len(Ps)), sharex=True)
    for i, (region, P) in enumerate(Ps.items()):
        a, b = axs[2 * i], axs[2 * i + 1]
        shade(a, Us[region]); shade(b, Us[region])
        a.plot(P.date, P.dv, "-", color="0.7", lw=0.6)
        sc = a.scatter(P.date, P.dv, c=P.ccabs, cmap=CMAP, vmin=VMIN, vmax=0.8, s=8 + 0.1 * P.ncell, lw=0.3, edgecolor="0.3")
        b.plot(P.date, P.cc, "-", color="0.4", lw=0.6, label="data CC (median over series, series-demeaned)")
        b.scatter(P.date, P.cc, c=P.ccabs, cmap=CMAP, vmin=VMIN, vmax=0.8, s=8 + 0.1 * P.ncell, lw=0.3, edgecolor="0.3")
        b.plot(P.date, P.c2, "-", color="#c00000", lw=0.6, alpha=0.6, label="c2 noise-window CC (same)")
        b2 = b.twinx(); b2.plot(P.date, P.lnn, ":", color="#2e75b6", lw=0.8); b2.set_ylabel("ln n (demeaned)", color="#2e75b6", fontsize=8)
        a.set_ylabel("dv/v (%)"); b.set_ylabel("Δ stretch CC"); b.legend(loc="lower left", fontsize=7)
        a.set_title(f"{region}: per-episode median over all family-station series (series-demeaned); size = number of cells; "
                    f"colour = median CC", fontsize=9)
    fig.colorbar(sc, ax=axs, shrink=0.4, label="median stretch CC")
    fig.savefig(f"{C.FIG}/figA3_pooled_dvv_cc_raw.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def raw_list(E, region):
    """cells with series-demeaned CC below -3 robust sigma while |dv/v - series median| < its bootstrap sigma (raw form)"""
    Q = E[np.isfinite(E.cc3c_25)].copy()
    g = Q.groupby("series")
    Q["cc_a"] = Q.cc3c_25 - g.cc3c_25.transform("median")
    Q["cc_s"] = g.cc_a.transform(rob)
    Q["dv_a"] = Q.dvv25 - g.dvv25.transform("median")
    Q["flat"] = Q.dv_a.abs() < Q.sig25
    Q["low"] = Q.cc_a < -3 * Q.cc_s
    Q["n_ser"] = g.cc3c_25.transform("size")
    Q = Q[Q.n_ser >= 8]
    L = Q[Q.low]
    L = L.assign(rel_n=L.lnn - Q.groupby("series").lnn.transform("median").loc[L.index],
                 rel_r=L.lnr - Q.groupby("series").lnr.transform("median").loc[L.index])
    cols = ["fam", "netsta", "epoch", "episode", "cls", "date", "n", "cc3c_25", "cc_a", "cc_s", "dvv25", "dv_a", "sig25", "flat",
            "rel_n", "rel_r", "dcc_shift", "shift2p_25", "glitch"]
    L[cols].sort_values(["episode", "netsta"]).to_csv(f"{C.TB}/rawlist_cells_{region}.csv", index=False)
    # episode clustering: number of families flagging the same episode at a station
    agg = L[L.flat].groupby(["netsta", "episode"]).agg(nfam=("fam", "nunique"), date=("date", "median"), cls=("cls", "first"),
                                                       rel_n=("rel_n", "median"), rel_r=("rel_r", "median"),
                                                       dcc_shift=("dcc_shift", "median"), glitch=("glitch", "any")).reset_index()
    tot = Q.groupby(["netsta", "episode"]).fam.nunique().rename("nfam_measured")
    agg = agg.join(tot, on=["netsta", "episode"]).sort_values(["nfam", "date"], ascending=[False, True])
    agg.to_csv(f"{C.TB}/rawlist_episodes_{region}.csv", index=False)
    s = dict(region=region, cells=len(Q), low=int(Q.low.sum()), low_flat=int((Q.low & Q.flat).sum()),
             frac_low=float(Q.low.mean()), frac_low_flat=float((Q.low & Q.flat).mean()),
             frac_flat_all=float(Q.flat.mean()),
             low_flat_relN_med=float(L[L.flat].rel_n.median()), low_flat_relR_med=float(L[L.flat].rel_r.median()),
             ep_ge3=int((agg.nfam >= 3).sum()))
    return s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pilot", action="store_true"); a = ap.parse_args()
    os.makedirs(C.FIG, exist_ok=True); os.makedirs(C.TB, exist_ok=True)
    Us = {r: C.units(r) for r in ("SVI", "WA_NW")}
    fig_fam10(Us["SVI"])
    Es = {r: C.load_episode(r) for r in ("SVI", "WA_NW")}
    fig_pilot(Es["SVI"][Es["SVI"].fam.isin(C.PILOT_FAMS) & Es["SVI"].netsta.isin(C.PILOT_STAS)], Us["SVI"])
    Ps = {r: pooled(Es[r], Us[r], r) for r in Es}
    fig_pooled(Ps, Us)
    S = [raw_list(Es["SVI"][Es["SVI"].fam.isin(C.PILOT_FAMS) & Es["SVI"].netsta.isin(C.PILOT_STAS)], "SVI_pilot")]
    if not a.pilot:
        S += [raw_list(Es[r], r) for r in Es]
    pd.DataFrame(S).to_csv(f"{C.TB}/rawlist_summary.csv", index=False)
    print(pd.DataFrame(S).to_string())
    open(f"{C.OUT}/RAW_DONE", "w").write(C.pt_now() + "\n")


if __name__ == "__main__":
    main()
