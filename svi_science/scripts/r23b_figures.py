#!/usr/bin/env python
"""r23b_figures (2026-10-11) -- figures for R23b (findings 1 and 6 on the 23-year SVI record).
Outputs svi_science/results/r23b/figures/r23b_{f1_neighbours,f1b_why,f6_size,f6b_ets_sizes}.png"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
TB = f"{ROOT}/svi_science/results/r23b/tables"
FIG = f"{ROOT}/svi_science/results/r23b/figures"
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"     # validated default categorical order (dataviz palette)
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "lines.linewidth": 2})
BLK = {"C_eq|2003-13": ("2003–13 (10+3 majors, new)", C1), "C_eq|2014-26": ("2014–26 (39 majors)", C2),
       "C_eq|2003-26": ("2003–26 together (52)", C3)}


def f1_neighbours():
    B = pd.read_csv(f"{TB}/f1/neighbour_corr_bins.csv")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    cen = {"(0, 3]": 1.7, "(3, 6]": 4.2, "(6, 10]": 7.7, "(10, 20]": 14.1, "(20, 40]": 28.3, "(40, 80]": 56.6, "(80, 400]": 120}
    for k, (lab, col) in BLK.items():
        g = B[B.block == k]
        x = g.b.map(cen)
        ax[0].errorbar(x, g["mean"], yerr=g.se, color=col, marker="o", ms=5, capsize=2, label=lab)
    g = B[B.block == "once|S15f_2010-13"]
    ax[0].plot(g.b.map(cen), g["mean"], color=C4, ls="--", lw=1.5, marker="s", ms=4, label="2010–13 counted once, one network")
    ax[0].axhline(0, color=MUTED, lw=0.8); ax[0].set_xscale("log")
    ax[0].set_xlabel("family separation (km, bin centre)"); ax[0].set_ylabel("mean residual correlation across ETS (± s.e.)")
    ax[0].set_title("Neighbouring families change together, in both decades", fontsize=10.5, color=INK)
    ax[0].legend(frameon=False, fontsize=8.5)
    for k, (lab, col) in BLK.items():
        P = pd.read_csv(f"{TB}/f1/pairs_{k.replace('|', '_')}.csv")
        q = pd.qcut(P.jaccard, 6, duplicates="drop")
        s = P.groupby(q, observed=True).agg(j=("jaccard", "mean"), c=("corr", "mean"), sd=("corr", "std"), n=("corr", "size"))
        ax[1].errorbar(s.j, s.c, yerr=s.sd / np.sqrt(s.n), color=col, marker="o", ms=5, capsize=2, label=lab)
    ax[1].axhline(0, color=MUTED, lw=0.8)
    ax[1].set_xlabel("shared template stations (Jaccard)"); ax[1].set_ylabel("mean residual correlation (± s.e.)")
    ax[1].set_title("…and against shared stations", fontsize=10.5, color=INK)
    fig.tight_layout(); fig.savefig(f"{FIG}/r23b_f1_neighbours.png", dpi=150); plt.close(fig)


def f1b_why():
    R = json.load(open(f"{TB}/f1/f1_results.json"))
    preds = [("in_reach", "inside tremor extent"), ("km_outside", "km outside tremor extent"), ("u_along_path", "position along front path"),
             ("dir_mode", "same-direction pattern"), ("T1", "days since last active"), ("T2", "previous big-ETS residual"),
             ("T3", "cumulative past residual"), ("detectability", "station availability"), ("kinematics_set", "ALL front terms"),
             ("budget_set", "ALL budget terms"), ("all", "everything")]
    des = [("A_ceq_2003-13", "2003–13 (10 big ETS, new)", C1), ("A_ceq_2014-26", "2014–26 (13)", C2), ("A_ceq_2003-26", "2003–26 together (23)", C3)]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw=dict(width_ratios=[1.05, 1]))
    y = np.arange(len(preds))[::-1]; h = 0.26
    for i, (k, lab, col) in enumerate(des):
        pr = {x["predictor"]: x["r2_within"] for x in R[k]["partial_r2"]}
        ax[0].barh(y + (1 - i) * h, [pr.get(p, np.nan) for p, _ in preds], height=h - 0.03, color=col, label=lab)
    ax[0].set_yticks(y); ax[0].set_yticklabels([l for _, l in preds]); ax[0].set_xlabel("share of within-episode variance explained (R²)")
    ax[0].set_title("What explains which patches are strong in a big ETS", fontsize=10.5, color=INK)
    ax[0].legend(frameon=False, fontsize=8.5, loc="center right"); ax[0].grid(axis="y", visible=False)
    ev = R["A_ceq_2003-26"]["episode_level_values"]
    inr = json.load(open(f"{TB}/f1/posthoc_persistence.json"))["2003-26_in_reach_frac"]
    for e, U, w in zip(ev["eps"], ev["U"], ev["trem_width"]):
        col = C1 if e < "2014" else C2
        half = inr.get(e, 1) < 0.9
        ax[1].scatter(w, U, s=46 if half else 34, color=col, marker="D" if half else "o", edgecolor=SURF, linewidth=1.2, zorder=3)
        ax[1].annotate(e, (w, U), xytext=(3, 3), textcoords="offset points", fontsize=7, color=MUTED)
    ax[1].set_yscale("log"); ax[1].set_xlabel("along-strike extent of tremor in the ETS (km, 5–95 %)")
    ax[1].set_ylabel("how uneven the ETS is (χ²/df of patch counts)")
    rho = R["A_ceq_2003-26"]["episode_level"]["trem_width"]["rho"]
    ax[1].set_title(f"2003–13 ETS all covered the whole box (Spearman 2003–26 {rho:+.2f})", fontsize=10.5, color=INK)
    from matplotlib.lines import Line2D
    ax[1].legend(handles=[Line2D([], [], marker="o", ls="", color=C1, label="2003–13"), Line2D([], [], marker="o", ls="", color=C2, label="2014–26"),
                          Line2D([], [], marker="D", ls="", color=MUTED, label="half-box ETS (in-reach < 0.9)")], frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/r23b_f1b_why.png", dpi=150); plt.close(fig)


SEGSPAN = [("era detector", "2003-02-01", "2006-06-30", "#f0efec"), ("frozen B2/B3", "2006-07-01", "2009-08-26", "#e4e3df"),
           ("frozen S15", "2009-08-27", "2013-12-31", "#f0efec"), ("forward", "2014-01-01", "2026-12-31", SURF)]


def shade(ax):
    for lab, a, b, c in SEGSPAN:
        a_, b_ = pd.Timestamp(a).year + pd.Timestamp(a).dayofyear / 365.25, pd.Timestamp(b).year + pd.Timestamp(b).dayofyear / 365.25
        ax.axvspan(a_, b_, color=c, zorder=0, lw=0)


def f6_size():
    Y = pd.read_csv(f"{TB}/f6/yearly_group_medians.csv")
    T = json.load(open(f"{TB}/f6/f6_size.json"))
    fig, axs = plt.subplots(3, 4, figsize=(14, 8.2), sharex=True, sharey=True)
    for g, ax in zip(range(1, 13), axs.ravel()):
        y = Y[(Y.group == g) & (Y.n >= 10)]
        shade(ax)
        ax.plot(y.year + 0.5, y.dM, color=C1, marker="o", ms=3.5, lw=1.6, label="ΔM above shared completeness floor")
        ax.plot(y.year + 0.5, y.dN, color=C2, marker="o", ms=3.0, lw=1.2, label="noise ΔN")
        ax.axhline(0, color=MUTED, lw=0.7)
        ax.set_title(f"group {g} (n={int(y.n.sum())})", fontsize=9, color=INK)
    axs[0, 0].legend(frameon=False, fontsize=7.5, loc="lower left")
    fig.supylabel("yearly median log10 amplitude (family-relative, fixed CN stations)")
    fig.supxlabel("year (shading: era detector 2003–06 · frozen B2/B3 2006–09 · frozen S15 2009–13 · forward 2014–)")
    a = T["T6a"]
    fig.suptitle(f"LFE size 2003–2026 per family group, ETS only, stations PGC/SNB/VGZ/NLLB — trend {a['mean_slope_per_decade']:+.3f} "
                 f"log units/decade (95 % {a['ci95'][0]:+.3f}…{a['ci95'][1]:+.3f}, p {a['p_perm']:.2f})", fontsize=10.5, color=INK)
    fig.tight_layout(); fig.savefig(f"{FIG}/r23b_f6_size.png", dpi=150); plt.close(fig)


def f6b_ets_sizes():
    E = pd.read_csv(f"{TB}/f6/episode_sizes.csv")
    R = json.load(open(f"{TB}/f6/f6_rates.json")); Rp = json.load(open(f"{TB}/f6/f6_rates_posthoc.json"))
    B = E[E.big].copy(); B["yr"] = pd.to_datetime(B.start).dt.year + pd.to_datetime(B.start).dt.dayofyear / 365.25
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.0), gridspec_kw=dict(width_ratios=[1.3, 1.0, 1.0]))
    shade(ax[0])
    ax[0].plot(B.yr, B.S, color=C1, marker="o", ms=6, lw=0, label="our Σ rate (Bostock-equivalent LFEs)")
    pre = B[B.n_bostock > 0]
    ax[0].plot(pre.yr, pre.n_bostock, color=C2, marker="s", ms=5, lw=0, label="Bostock's own catalogue")
    ax[0].plot(B.yr, B.S_kc, color=C4, marker="^", ms=4, lw=0, label="Σ rate_kc (chained completeness)")
    for r in B.itertuples():
        ax[0].annotate(r.episode, (r.yr, r.S), xytext=(2, 4), textcoords="offset points", fontsize=6.5, color=MUTED, rotation=60)
    ax[0].set_yscale("log"); ax[0].set_ylabel("LFEs in the ETS (115 families)")
    ax[0].set_title("Big-ETS sizes 2003–2026", fontsize=10.5, color=INK); ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    p = B[B.start >= "2010-01-06"]
    cols = [C1 if s == "S15f" else (C2 if f >= 0.9 else C3) for s, f in zip(p.seg, p.frac_S15)]
    ax[1].scatter(p.yr, p.S_per_pnsn_expday, color=cols, s=40, edgecolor=SURF, linewidth=1.2, zorder=3)
    for r in p.itertuples():
        ax[1].annotate(r.episode, (r.yr, r.S_per_pnsn_expday), xytext=(3, 2), textcoords="offset points", fontsize=6.5, color=MUTED)
    ax[1].set_ylabel("Σ rate per PNSN tremor family-day"); ax[1].set_title("Per tremor (PNSN 2010–26)", fontsize=10.5, color=INK)
    from matplotlib.lines import Line2D
    ax[1].legend(handles=[Line2D([], [], marker="o", ls="", color=C1, label="2010–13 frozen S15"),
                          Line2D([], [], marker="o", ls="", color=C2, label="2014–18, same model + stations"),
                          Line2D([], [], marker="o", ls="", color=C3, label="2019–26, S14/S13 stations")], frameon=False, fontsize=8, loc="upper right")
    rows = [("2014–26 / 2003–13\n(rate)", R["DEC_rate"]), ("2014–26 / 2003–13\n(rate_kc)", R["DEC_rate_kc"]),
            ("2014–18 / 2010–13\nsame model + stations", R["L4L_rate"]), ("same, without 2014-11\nand 2017-03 (post hoc)", Rp["L4L_rate_clean"])]
    yy = np.arange(len(rows))[::-1]
    for y_, (lab, v) in zip(yy, rows):
        ax[2].plot(v["ci95_ratio"], [y_, y_], color=C1, lw=2.5); ax[2].plot(v["ratio"], y_, "o", color=C1, ms=7)
        ax[2].annotate(f"{v['ratio']:.2f}  (p {v['p_perm']:.3f})", (v["ratio"], y_), xytext=(0, 7), textcoords="offset points", fontsize=8, ha="center", color=INK)
    ax[2].axvline(1, color=MUTED, lw=0.9); ax[2].set_xscale("log"); ax[2].set_yticks(yy); ax[2].set_yticklabels([r[0] for r in rows], fontsize=8)
    ax[2].set_xlabel("per-family ETS size ratio\n(median, 95 % CI; p = episode-label permutation)"); ax[2].set_title("Did family ETS sizes change?", fontsize=10.5, color=INK)
    ax[2].grid(axis="y", visible=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/r23b_f6b_ets_sizes.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    f1_neighbours(); f1b_why(); f6b_ets_sizes(); f6_size()
