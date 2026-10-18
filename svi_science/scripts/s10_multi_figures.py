#!/usr/bin/env python
"""s10_multi_figures -- S10 multi figures (multi/PLAN.md §6). User plot rules: dv/v points coloured by CC with the
reference, marker = station, NO error bars, averages CC^2-weighted; uncertainties only in text/tables."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
import s10_multi_analysis as A

STA = M.CORE5
MK = {"CN.PGC": "s", "PB.B011": "o", "PB.B009": "^", "PB.B010": "v", "PB.B926": "D", "CN.VGZ": "P", "CN.SNB": "X"}
CMAP, VMIN, VMAX = "viridis", 0.4, 1.0
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
                     "grid.linewidth": 0.5, "font.size": 9})


def sta_legend(stas=STA, extra=()):
    return [Line2D([], [], marker=MK[s], ls="", color="0.6", markeredgecolor="k", markeredgewidth=0.4, label=s.split(".")[1])
            for s in stas] + list(extra)


def episode_panels(D, U):
    maj = U[U.cls == "major"].set_index("episode")
    X = D[D.netsta.isin(STA) & (D.cls == "major") & D.v24.notna()].copy()
    X["x"] = (pd.to_datetime(X.date) - pd.to_datetime(X.episode.map(maj.start))).dt.days
    eps = [e for e in maj.sort_values("start").index if e in set(X.episode)]
    half = int(np.ceil(len(eps) / 2))
    for part, sub in (("a", eps[:half]), ("b", eps[half:])):
        nc = 6; nr = int(np.ceil(len(sub) / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc, 2.6 * nr + 0.9), sharey=True, squeeze=False, layout="constrained")
        for k, e in enumerate(sub):
            ax = axs.flat[k]; g = X[X.episode == e]
            for j, ns in enumerate(STA):
                h = g[g.netsta == ns]
                if len(h):
                    sc = ax.scatter(h.x + (j - 2) * 0.1, h.v24, c=h.cc24, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=7, marker=MK[ns],
                                    edgecolor="none", alpha=0.75)
            m = np.average(g.v24, weights=g.w24)
            ax.axhline(0, color="0.5", lw=0.5); ax.axhline(m, color="tab:red", lw=1.2)
            ax.set_title(f"{e} ({maj.start[e]}): {g.fam.nunique()} fam, mean {m:+.3f} %", fontsize=7.5)
            ax.set_ylim(-1.05, 1.05); ax.tick_params(labelsize=7)
        for ax in axs.flat[len(sub):]:
            ax.axis("off")
        for ax in axs[:, 0]:
            ax.set_ylabel("dv/v (%)", fontsize=8)
        for ax in axs[-1, :]:
            ax.set_xlabel("days since episode start", fontsize=8)
        fig.legend(handles=sta_legend(extra=[Line2D([], [], color="tab:red", label="CC²-weighted mean (all families, stations, days)")]),
                   loc="outside lower center", ncol=7, fontsize=8)
        cb = fig.colorbar(sc, ax=axs, shrink=0.5, pad=0.01, aspect=30); cb.set_label("CC with reference (2–4 s)")
        fig.suptitle("Daily dv/v of all families (2–4 s after S, raw day stacks), one panel per major episode — colour = CC, "
                     "marker = station", fontsize=10)
        fig.savefig(f"{M.FIG}/fig_episode_panels_{part}.png", dpi=110); plt.close(fig)


def heatmap(U, tag="dm"):
    F = pd.read_csv(f"{M.TB}/fe_{'core5_d24_all' if tag == 'dm' else 'es_core5_d24_all'}.csv").merge(U[["episode", "start", "cls"]], on="episode")
    order = U.sort_values("start")
    majc = [e for e in order[order.cls == "major"].episode if e in set(F.episode)]
    burc = [e for e in order[order.cls == "burst"].episode if (F.episode == e).sum() >= 3]
    V = F.pivot(index="fam", columns="episode", values="value")
    fo = V[majc].notna().sum(1).sort_values(ascending=False).index
    lim = float(np.nanquantile(np.abs(V.to_numpy()), 0.95))
    fig, axs = plt.subplots(1, 2, figsize=(18, 12), gridspec_kw=dict(width_ratios=[len(majc), max(len(burc), 1)]), layout="constrained")
    for ax, cols, ttl in ((axs[0], majc, "major episodes"), (axs[1], burc, "between-episode bursts (≥ 3 families)")):
        im = ax.imshow(V.reindex(index=fo, columns=cols).to_numpy(), aspect="auto", cmap="RdBu", vmin=-lim, vmax=lim, interpolation="nearest")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.lstrip("B")[:7] for c in cols], rotation=90, fontsize=6)
        ax.set_yticks(range(len(fo))); ax.set_yticklabels(fo, fontsize=5.5); ax.set_title(ttl); ax.grid(False)
        ax.set_facecolor("0.93")
    axs[0].set_ylabel("family (sorted by number of major episodes)")
    cb = fig.colorbar(im, ax=axs, shrink=0.5, pad=0.01); cb.set_label("family-episode dv/v (%), CC²-weighted mean; blue = faster")
    unit = "mean of daily values" if tag == "dm" else "per-episode stacks"
    fig.suptitle(f"Family × episode dv/v (2–4 s, core 5 stations, {unit}); grey = no value", fontsize=11)
    fig.savefig(f"{M.FIG}/fig_heatmap_{tag}.png", dpi=110); plt.close(fig)


def pooled_series(T, tag):
    T = T.copy(); T["t"] = pd.to_datetime(T.start)
    fig, axs = plt.subplots(3, 1, figsize=(15, 11), sharex=True, layout="constrained")
    ax = axs[0]
    for cls, mk in (("major", "o"), ("burst", "^")):
        h = T[T.cls == cls]
        sc = ax.scatter(h.t, h.core5_d24_value, c=h.core5_d24_mean_cc, cmap=CMAP, vmin=VMIN, vmax=VMAX,
                        s=8 + 4 * h.core5_d24_n_fam, marker=mk, edgecolor="k", linewidth=0.4, alpha=0.9)
    ax.set_title("Pooled episode dv/v, all families and core stations (CC²-weighted); marker size = number of families; "
                 "o major episode, △ between-episode burst", fontsize=9)
    ax = axs[1]
    for k, mk, lab in (("pgc_d24", "s", "PGC only"), ("bh_d24", "o", "boreholes only (B009/B010/B011/B926)")):
        h = T[T.cls == "major"]
        ax.scatter(h.t, h[f"{k}_value"], c=h[f"{k}_mean_cc"], cmap=CMAP, vmin=VMIN, vmax=VMAX, s=8 + 4 * h[f"{k}_n_fam"].fillna(0),
                   marker=mk, edgecolor="k", linewidth=0.4, label=lab)
    ax.legend(loc="upper left", fontsize=8); ax.set_title("Major episodes: PGC-only vs borehole-only pooled series", fontsize=9)
    ax = axs[2]
    h = T[T.cls == "major"]
    ax.scatter(h.t, h.core5_c24_value, c=h.core5_c24_mean_cc, cmap=CMAP, vmin=0, vmax=VMAX, s=8 + 4 * h.core5_c24_n_fam.fillna(0),
               marker="o", edgecolor="k", linewidth=0.4)
    ax.set_title("Scrambled control (noise windows through the same path), major episodes; colour = its CC (scale 0–1)", fontsize=9)
    for ax in axs:
        ax.axhline(0, color="0.5", lw=0.6); ax.set_ylabel("dv/v (%)")
        for d, lab in (("2015-01-12", "PGC E2"), ("2017-08-10", "PGC E3"), ("2005-11-02", "boreholes start")):
            ax.axvline(pd.Timestamp(d), color="0.6", lw=0.6, ls=":")
    lim = np.nanquantile(np.abs(T.core5_d24_value), 0.98) * 1.15
    for ax in axs[:2]:
        ax.set_ylim(-lim, lim)
    cb = fig.colorbar(sc, ax=axs[:2], shrink=0.6, pad=0.01); cb.set_label("mean CC with reference (2–4 s)")
    unit = "CC²-weighted mean of daily values" if tag == "dm" else "CC²-weighted mean of per-episode stacks"
    fig.suptitle(f"SVI 2003–2026: dv/v per slow-slip episode pooled over LFE families (2–4 s after S; {unit})", fontsize=11)
    fig.savefig(f"{M.FIG}/fig_pooled_series_{tag}.png", dpi=110); plt.close(fig)


def corr_matrix(key, tag):
    P = pd.read_csv(f"{M.TB}/pairs_{key}.csv")
    fams = sorted(set(P.fam_i) | set(P.fam_j)); ix = {f: i for i, f in enumerate(fams)}
    R = np.full((len(fams), len(fams)), np.nan)
    for r in P[P.eligible].itertuples():
        R[ix[r.fam_i], ix[r.fam_j]] = R[ix[r.fam_j], ix[r.fam_i]] = r.r
    o = np.argsort(-np.nanmean(R, 1))
    fig, ax = plt.subplots(figsize=(10, 9), layout="constrained")
    im = ax.imshow(R[np.ix_(o, o)], cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest"); ax.grid(False); ax.set_facecolor("0.93")
    ax.set_xticks(range(len(fams))); ax.set_xticklabels(np.array(fams)[o], rotation=90, fontsize=5)
    ax.set_yticks(range(len(fams))); ax.set_yticklabels(np.array(fams)[o], fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.7).set_label("Pearson r of episode dv/v (pairs sharing ≥ 8 major episodes)")
    unit = "means of daily values" if tag == "dm" else "per-episode stacks"
    ax.set_title(f"Cross-family agreement, core 5 stations, 2–4 s, {unit} ({int(P.eligible.sum())} pairs; grey < 8 shared); "
                 "families sorted by mean r", fontsize=9)
    fig.savefig(f"{M.FIG}/fig_corr_matrix_{tag}.png", dpi=110); plt.close(fig)


def perm_hists():
    AG = pd.read_csv(f"{M.TB}/agreement.csv").set_index("set")
    sets = [("core5_d24", "day means: core 5, 2–4 s (headline)"), ("pgc_d24", "day means: PGC only"), ("bh_d24", "day means: boreholes only"),
            ("core5_d35", "day means: core 5, 3–5 s"), ("es_core5_d24", "episode stacks: core 5, 2–4 s"), ("es_pgc_d24", "episode stacks: PGC only"),
            ("es_bh_d24", "episode stacks: boreholes only"), ("es_core5_d35", "episode stacks: core 5, 3–5 s"),
            ("es_core5_c24", "episode stacks: scrambled control"), ("dm_core5_d24_deseason", "day means, annual cycle removed"),
            ("es_core5_d24_deseason", "episode stacks, annual cycle removed"), ("es_bh_d24_deseason", "episode stacks, boreholes, annual removed")]
    fig, axs = plt.subplots(3, 4, figsize=(18, 10), layout="constrained")
    for ax, (k, lab) in zip(axs.flat, sets):
        p = f"{M.TB}/perm_{k}.npy"
        if not os.path.exists(p) or k not in AG.index:
            ax.axis("off"); continue
        Ms = np.load(p)
        ax.hist(Ms, bins=60, color="0.6")
        ax.axvline(AG.loc[k, "perm_p99"], color="0.2", ls="--", lw=1, label="99th percentile of shuffle")
        ax.axvline(AG.loc[k, "M"], color="tab:red", lw=2, label=f"observed median r = {AG.loc[k, 'M']:.3f}")
        ax.set_title(f"{lab}\n{int(AG.loc[k, 'n_pairs'])} pairs, p = {AG.loc[k, 'perm_p']:.4f}", fontsize=8.5)
        ax.set_xlabel("median pairwise r"); ax.legend(fontsize=7)
    fig.suptitle("Episode labels shuffled within each family, 10,000 draws", fontsize=11)
    fig.savefig(f"{M.FIG}/fig_permutation.png", dpi=110); plt.close(fig)


def within(D):
    FRt = pd.read_csv(f"{S.ROOT}/svi_science/results/s5/step3_front_arrivals.csv")
    FRt["front"] = pd.to_datetime(FRt.t_front).dt.normalize()
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharey=True, layout="constrained")
    for ax, (keys, ttl) in zip(axs, ((("pgc_d24", "bh_d24"), "data (2–4 s)"), (("core5_c24",), "scrambled control (2–4 s)"))):
        for k in keys:
            var, grp = k.split("_")[1], k.split("_")[0]
            X = A.within_episode(D, var, grp, FRt)
            v, w, s, c = A.VARS[var]
            for ns, h in X.groupby("netsta"):
                b = h.assign(wy=h[w] * h.y).groupby("x")
                prof = b.wy.sum() / b[w].sum(); cc = b[c].mean(); n = b.size()
                sc = ax.scatter(prof.index + 0.08 * STA.index(ns) - 0.16, prof, c=cc, cmap=CMAP, vmin=0,
                                vmax=VMAX, marker=MK[ns], s=10 + 2 * np.sqrt(n), edgecolor="k", linewidth=0.3)
        ax.axhline(0, color="0.5", lw=0.6); ax.axvline(0, color="0.5", lw=0.6, ls=":")
        ax.set_xlabel("days since the front reached the family (S5 fronts 2015–24; else episode start)")
        ax.set_title(ttl, fontsize=9)
    axs[0].set_ylabel("dv/v minus its family-station-episode mean (%)")
    axs[0].legend(handles=sta_legend(), fontsize=7, loc="lower right")
    fig.colorbar(sc, ax=axs, shrink=0.7).set_label("mean CC with reference")
    fig.suptitle("Within-episode change: CC²-weighted mean per day, per station, all families pooled (major episodes)", fontsize=10)
    fig.savefig(f"{M.FIG}/fig_within_episode.png", dpi=110); plt.close(fig)


def lapse(T, tag):
    pre = "" if tag == "dm" else "es_"
    F24 = pd.read_csv(f"{M.TB}/fe_{pre}core5_d24.csv"); F35 = pd.read_csv(f"{M.TB}/fe_{pre}core5_d35.csv")
    F = F24.merge(F35, on=["fam", "episode"], suffixes=("", "_35"))
    h = T[T.cls == "major"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.5), layout="constrained")
    axs[0].scatter(F.value, F.value_35, c=F.mean_cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=8, alpha=0.6, edgecolor="none")
    axs[0].set_title("family-episode values (core 5)", fontsize=9)
    sc = axs[1].scatter(h.core5_d24_value, h.core5_d35_value, c=h.core5_d24_mean_cc, cmap=CMAP, vmin=VMIN, vmax=VMAX,
                        s=8 + 4 * h.core5_d24_n_fam, edgecolor="k", linewidth=0.4)
    axs[1].set_title("pooled major-episode values; size = number of families", fontsize=9)
    for ax in axs:
        lim = np.nanquantile(np.abs(np.r_[F.value, F.value_35]), 0.98) if ax is axs[0] else np.nanmax(np.abs(np.r_[h.core5_d24_value, h.core5_d35_value])) * 1.1
        ax.plot([-lim, lim], [-lim, lim], color="0.4", lw=0.8, ls="--"); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel("dv/v 2–4 s (%)"); ax.set_ylabel("dv/v 3–5 s (%)"); ax.set_aspect("equal")
    fig.colorbar(sc, ax=axs, shrink=0.7).set_label("mean CC with reference (2–4 s)")
    fig.suptitle("Lapse check: does the later coda (3–5 s) show the same episode-to-episode change?", fontsize=10)
    fig.savefig(f"{M.FIG}/fig_lapse_{tag}.png", dpi=110); plt.close(fig)


def seasonal(Ts):
    R = json.load(open(f"{M.TB}/result.json"))
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5), sharey=True, layout="constrained")
    for ax, (tag, T) in zip(axs, Ts.items()):
        h = T[(T.cls == "major") & (T.core5_d24_n_fam >= 3)].copy()
        h["doy"] = pd.to_datetime(h.start).dt.dayofyear
        for k, mk, lab in (("core5_d24", "o", "all core stations"), ("pgc_d24", "s", "PGC only"), ("bh_d24", "^", "boreholes only")):
            sc = ax.scatter(h.doy, h[f"{k}_value"], c=h[f"{k}_mean_cc"], cmap=CMAP, vmin=VMIN, vmax=VMAX, marker=mk,
                            s=10 + (3 * h[f"{k}_n_fam"].fillna(0) if k == "core5_d24" else 25), edgecolor="k", linewidth=0.4,
                            label=lab, alpha=0.9 if k == "core5_d24" else 0.6)
        b0, b1 = R[f"annual_model_{tag}"]; dd = np.arange(1, 366); ph = 2 * np.pi * (dd + 10) / 365.25
        ax.plot(dd, b0 * np.cos(ph) + b1 * np.sin(ph), color="tab:red", lw=1.5, label="annual sinusoid fit (all core stations)")
        B = R[f"battery_{tag}"]
        ax.set_title(f"{'means of daily values' if tag == 'dm' else 'per-episode stacks'}: amplitude {B['core5_d24_annual_amp']:.3f} %, "
                     f"peak day {B['core5_d24_annual_peak_doy']:.0f}, permutation p = {B['core5_d24_annual_perm_p']:.4f}", fontsize=9)
        ax.axhline(0, color="0.5", lw=0.6); ax.set_xlabel("day of year of episode start")
        ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335])
        ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    axs[0].set_ylabel("pooled episode dv/v (%)"); axs[0].legend(fontsize=7, loc="lower left")
    fig.colorbar(sc, ax=axs, shrink=0.7).set_label("mean CC with reference (2–4 s)")
    fig.suptitle("Major-episode dv/v against time of year (episodes with ≥ 3 families; 2003–2026)", fontsize=10)
    fig.savefig(f"{M.FIG}/fig_seasonal.png", dpi=110); plt.close(fig)


def main():
    D, U = A.load()
    Ts = {t: pd.read_csv(f"{M.TB}/pooled_series_{t}.csv", index_col=0) for t in ("dm", "es")}
    episode_panels(D, U)
    for t, T in Ts.items():
        heatmap(U, t); pooled_series(T, t); lapse(T, t)
    corr_matrix("core5_d24", "dm"); corr_matrix("es_core5_d24", "es"); perm_hists(); within(D); seasonal(Ts)
    open(f"{M.OUT}/FIGURES_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
