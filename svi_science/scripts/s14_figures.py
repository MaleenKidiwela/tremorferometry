#!/usr/bin/env python
"""s14_figures -- S14 figures: CC vs stack size (the confound), contrasts vs scrambled band, C3 day profile, co-located
residuals. No error bars (scrambled 95 % bands as shading).
    python s14_figures.py TAG [TAG ...]       (TAG = SVI_pilot, SVI, WA_NW)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import s14_common as C

COMPS = ["C1_major_minus_burst", "C2_post_minus_pre", "C3_late_minus_early", "C3_pre_minus_early"]
WINS = ["S", "2035", "25", "3560", "36", "Z25", "H25"]
WLAB = {"S": "direct S", "2035": "2–3.5 s", "25": "2–5 s", "3560": "3.5–6 s", "36": "3–6 s", "Z25": "Z 2–5", "H25": "H 2–5"}


def fig_confound(tag):
    region = tag.replace("_pilot", "")
    E = C.load_episode(region, pilot=tag.endswith("pilot"))
    E = E[np.isfinite(E.z_25) & np.isfinite(E.lnn) & np.isfinite(E.lnr)]
    g = E.groupby("series")
    x = E.lnn - g.lnn.transform("mean"); y = E.z_25 - g.z_25.transform("mean"); r = E.lnr - g.lnr.transform("mean")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    sc = ax[0].scatter(x, y, c=E.cc3c_25, cmap="viridis", vmin=0.2, vmax=1, s=6, alpha=0.6, lw=0)
    ax[0].set_xlabel("ln n repeats (series-demeaned)"); ax[0].set_ylabel("Fisher z of stretch CC (series-demeaned)")
    ax[1].scatter(r, y, c=E.cc3c_25, cmap="viridis", vmin=0.2, vmax=1, s=6, alpha=0.6, lw=0)
    ax[1].set_xlabel("ln (pre-P RMS / coda RMS) (series-demeaned)")
    for a, xx in zip(ax, (x, r)):
        b = np.polyfit(xx, y, 1); xs = np.linspace(np.nanpercentile(xx, 1), np.nanpercentile(xx, 99), 10)
        a.plot(xs, np.polyval(b, xs), "k-", lw=1.2); a.text(0.03, 0.93, f"slope {b[0]:+.3f}, r = {np.corrcoef(xx, y)[0, 1]:.2f}",
                                                             transform=a.transAxes, fontsize=9)
    fig.colorbar(sc, ax=ax, label="stretch CC"); fig.suptitle(f"{tag}: episode stretch CC vs stack size and noise (2–5 s)", fontsize=10)
    fig.savefig(f"{C.FIG}/figB1_confound_{tag}.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def fig_contrasts(tags):
    fig, axs = plt.subplots(len(COMPS), 1, figsize=(10, 2.6 * len(COMPS)), sharex=True)
    mk = {"SVI_pilot": "o", "SVI": "s", "WA_NW": "^"}; cl = {"SVI_pilot": "#7f7f7f", "SVI": "#1f4e79", "WA_NW": "#c55a11"}
    for i, comp in enumerate(COMPS):
        a = axs[i]
        for k, tag in enumerate(tags):
            f = f"{C.TB}/tests_{tag}.csv"
            if not os.path.exists(f):
                continue
            T = pd.read_csv(f); T = T[(T.comp == comp) & (T.get("kind", pd.Series("data", index=T.index)).fillna("data") == "data")]
            xs = [WINS.index(w) + 0.25 * (k - 1) for w in T.win.astype(str)]
            for x, (_, r) in zip(xs, T.iterrows()):
                a.fill_between([x - 0.1, x + 0.1], r.perm_med - 1.96 * r.perm_sd, r.perm_med + 1.96 * r.perm_sd, color=cl[tag], alpha=0.18, lw=0)
            a.scatter(xs, T.gamma, marker=mk[tag], color=cl[tag], s=30, label=tag, zorder=3)
        a.axhline(0, color="k", lw=0.6); a.set_ylabel("Δ Fisher z"); a.set_title(comp + "  (shading = scrambled-control 95 % band)", fontsize=9)
        if i == 0:
            a.legend(fontsize=8, ncol=3)
    axs[-1].set_xticks(range(len(WINS))); axs[-1].set_xticklabels([WLAB[w] for w in WINS])
    fig.savefig(f"{C.FIG}/figB2_contrasts_{'_'.join(tags)}.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def fig_profile(tags):
    fig, ax = plt.subplots(figsize=(7, 4))
    for tag in tags:
        f = f"{C.TB}/c3_profile_{tag}.csv"
        if os.path.exists(f):
            P = pd.read_csv(f)
            ax.plot(P.lb.astype(str), P.res, "-o", label=f"{tag} (n days {int(P.n.sum())})")
    ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("days since the slip front reached the family")
    ax.set_ylabel("mean REG residual, Fisher z (2–5 s, leave-day-out)"); ax.legend(fontsize=8)
    ax.set_title("Day stretch CC through an episode after removing stack size and noise", fontsize=9)
    fig.savefig(f"{C.FIG}/figB3_day_profile_{'_'.join(tags)}.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def fig_coloc(tag):
    region = tag.replace("_pilot", "")
    if region != "SVI":
        return
    sys.path.insert(0, HERE)
    import s14_tests as Tm
    E = C.load_episode(region, pilot=tag.endswith("pilot"))
    des = Tm.Design(E, "z_25", Tm.COVS, "series")
    Q = des.D.assign(res=des.yr)
    P = Q[Q.netsta.isin(C.COLOC["SVI"][0])].pivot_table(index=["fam", "episode"], columns="netsta", values="res")
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for a, b in zip(axs, ["PB.B009", "PB.B010", "PB.B011"]):
        ok = P["CN.PGC"].notna() & P[b].notna()
        a.scatter(P["CN.PGC"][ok], P[b][ok], s=6, alpha=0.5, color="#1f4e79", lw=0)
        a.set_xlabel("PGC residual z"); a.set_ylabel(f"{b} residual z")
        a.set_title(f"r = {np.corrcoef(P['CN.PGC'][ok], P[b][ok])[0, 1]:.2f}, n = {ok.sum()}", fontsize=9)
        a.axhline(0, color="k", lw=0.4); a.axvline(0, color="k", lw=0.4)
    fig.suptitle(f"{tag}: same family-episode, CC residual after stack size and noise (2–5 s)", fontsize=10)
    fig.savefig(f"{C.FIG}/figB4_coloc_{tag}.png", dpi=120, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    tags = sys.argv[1:] or ["SVI_pilot"]
    for t in tags:
        fig_confound(t); fig_coloc(t)
    fig_contrasts(tags); fig_profile(tags)
