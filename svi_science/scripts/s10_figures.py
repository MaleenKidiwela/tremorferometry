#!/usr/bin/env python
"""s10_figures -- S10 pilot figures (PLAN §8): daily eps vs time with ETS shading + metadata lines (gate f),
within-episode eps vs days since the front, precision vs number of repeats, POLARIS vs PGC same-day scatter."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

BLUE, ORANGE, AQUA, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
D = pd.read_parquet(f"{S.TB}/daily.parquet")
EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")


def shade(ax):
    for r in EPI.itertuples():
        a, b = pd.Timestamp(r.start), pd.Timestamp(r.end) + pd.Timedelta(days=1)
        ax.axvspan(a, b, color="#d9d8d3" if r.big else "#eeede9", lw=0, zorder=0)


def meta(ax, ns):
    for eid, a, b in S.EPOCHS.get(ns, []):
        for d in (a, b):
            if "2002" < d < "2027":
                ax.axvline(pd.Timestamp(d), color=INK, lw=0.9, ls="--", zorder=1)
    for d, lab in S.CHANGES:
        ax.axvline(pd.Timestamp(d), color=GRAY, lw=0.6, ls=":", zorder=1)


def daily_fig(fam):
    fig, axs = plt.subplots(len(S.CORE), 1, figsize=(13, 14), sharex=True)
    for ax, ns in zip(axs, S.CORE):
        shade(ax); meta(ax, ns)
        g = D[(D.fam == fam) & (D.netsta == ns)]
        t = pd.to_datetime(g.date)
        ax.errorbar(t, 100 * g.eps24, yerr=100 * g.sig24, fmt="o", ms=3, color=BLUE, ecolor=BLUE, elinewidth=0.6,
                    alpha=0.85, label="day stack ε (2–4 s) ± bootstrap σ", zorder=3)
        if "c2_eps24" in g:
            ax.plot(t, 100 * g.c2_eps24, "o", ms=3.5, mfc="none", mec=GRAY, mew=0.8, label="scrambled control (frame time)", zorder=2)
        rl = g.rail24 == True
        ax.plot(t[rl], np.zeros(rl.sum()), "x", ms=3, color=ORANGE, label="no peak within ±1 % (not a measurement)", zorder=2)
        ax.axhline(0, color=MUTED, lw=0.5)
        ax.set_ylim(-1.15, 1.15); ax.set_ylabel(f"{ns}\nε = −dv/v (%)")
        ax.text(0.005, 0.9, f"days {len(g)}, measured {g.eps24.notna().sum()}, median σ {100 * g.sig24.median():.2f} %",
                transform=ax.transAxes, fontsize=8, color=MUTED)
    axs[0].legend(loc="upper right", fontsize=7, ncol=3, frameon=False)
    axs[0].set_title(f"Family {fam}: daily active-day stacks, Z, each response epoch zero-mean (dashed = epoch boundary, "
                     f"dotted = instrument/network change, grey bands = ETS episodes, darker = big)", fontsize=9, loc="left")
    axs[-1].set_xlim(pd.Timestamp("2003-01-01"), pd.Timestamp("2026-10-01"))
    fig.tight_layout(); fig.savefig(f"{S.FIG}/s10_daily_fam{fam}.png", dpi=130); plt.close(fig)


def within_fig():
    P = pd.read_parquet(f"{S.TB}/within_episode_points.parquet")
    P = P[P.netsta.isin(S.CORE)]
    bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 6.5, 10.5, 40]
    labs = ["0", "1", "2", "3", "4–6", "7–10", "11+"]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, fam in zip(axs, S.FAMS):
        for est, col, lab, off in (("d_", BLUE, "data (day stacks)", -0.12), ("c2_", GRAY, "scrambled control", 0.12)):
            g = P[(P.fam == fam) & (P.est == est)].copy()
            g["b"] = pd.cut(g.dsf, bins, labels=False)
            ax.plot(g.b + off + np.random.default_rng(0).uniform(-0.08, 0.08, len(g)), 100 * g.de, ".", ms=2,
                    color=col, alpha=0.25)
            m, e, x = [], [], []
            for b, h in g.groupby("b"):
                w = 1 / h.sig ** 2
                m.append(100 * np.sum(w * h.de) / w.sum()); e.append(100 / np.sqrt(w.sum())); x.append(b + off)
            ax.errorbar(x, m, yerr=e, fmt="o-", ms=5, color=col, lw=1.5, label=lab, capsize=2)
        ax.axhline(0, color=MUTED, lw=0.5)
        ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs)
        ax.set_xlabel("days since the front reached the family"); ax.set_title(f"family {fam}", loc="left")
        ax.set_ylim(-1.0, 1.0)
    axs[0].set_ylabel("ε − episode mean (%)  [ε = −dv/v]"); axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Within-episode change at the core stations (2–4 s): weighted bin means ± 1σ", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{S.FIG}/s10_within_episode.png", dpi=130); plt.close(fig)


def precision_fig():
    P = pd.read_csv(f"{S.TB}/precision.csv")
    P["kind"] = np.select([P.netsta.str.startswith("PB."), P.netsta.isin(["CN.PGC", "CN.VGZ", "CN.SNB"])],
                          ["PB boreholes", "CN PGC/VGZ/SNB"], "POLARIS + LZB")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, win in zip(axs, (24, 35)):
        for kind, col in (("PB boreholes", BLUE), ("CN PGC/VGZ/SNB", ORANGE), ("POLARIS + LZB", AQUA)):
            g = P[(P.win == win) & (P.kind == kind)]
            med = g.groupby("n").sigma.median(); rf = g.groupby("n").rail_frac.mean(); cnt = g.groupby("n").size()
            ax.plot(med.index, 100 * med, "o-", color=col, lw=2, ms=6, label=kind)
            for n_, v in med.items():
                if rf[n_] >= 0.05:
                    ax.plot(n_, 100 * v, "o", ms=11, mfc="none", mec=col, mew=1)
                ax.annotate(f"{cnt[n_]}", (n_, 100 * v), textcoords="offset points", xytext=(4, 4), fontsize=7, color=MUTED)
        for lv in (0.1, 0.05):
            ax.axhline(lv, color=MUTED, lw=0.8, ls="--"); ax.text(10.5, lv * 1.05, f"{lv} %", fontsize=7, color=MUTED)
        ax.axhline(0.577, color=GRAY, lw=0.6, ls=":")
        ax.text(10.5, 0.6, "grid-limited (uniform in ±1 % gives 0.58 %)", fontsize=7, color=GRAY)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("repeats per stack (n, equal weights)")
        ax.set_title({24: "2–4 s after S", 35: "3–5 s after S"}[win], loc="left")
        ax.set_ylim(0.03, 1.2)
    axs[0].set_ylabel("σ of ε per stack (%) — disjoint half-split"); axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Precision vs number of repeats in the busiest days (median over series-days; numbers = series-days; "
                 "ringed = ≥ 5 % of draws at the ±1 % edge)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{S.FIG}/s10_precision_vs_n.png", dpi=130); plt.close(fig)


def polaris_fig():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.3), sharex=True, sharey=True)
    for ax, fam in zip(axs, S.FAMS):
        Df = D[(D.fam == fam) & (D.date < "2007-01-01")]
        ref = Df[(Df.netsta == "CN.PGC")][["date", "eps24", "sig24"]]
        for k, ns in enumerate(S.POLARIS):
            m = Df[Df.netsta == ns].merge(ref, on="date", suffixes=("", "_pgc"))
            m = m[np.isfinite(m.eps24) & np.isfinite(m.eps24_pgc)]
            if not len(m):
                continue
            ax.errorbar(100 * m.eps24_pgc, 100 * m.eps24, xerr=100 * m.sig24_pgc, yerr=100 * m.sig24, fmt="o", ms=3,
                        color=BLUE, ecolor="#c9d9ee", elinewidth=0.6, alpha=0.8, label="POLARIS / LZB day" if k == 0 else None)
        for ns, col in (("CN.VGZ", ORANGE), ("CN.SNB", AQUA)):
            m = Df[Df.netsta == ns].merge(ref, on="date", suffixes=("", "_pgc"))
            m = m[np.isfinite(m.eps24) & np.isfinite(m.eps24_pgc)]
            ax.plot(100 * m.eps24_pgc, 100 * m.eps24, "s", ms=4, color=col, label=f"{ns} day")
        ax.plot([-1, 1], [-1, 1], color=MUTED, lw=0.8); ax.axhline(0, color=GRAY, lw=0.4); ax.axvline(0, color=GRAY, lw=0.4)
        ax.set_title(f"family {fam}, 2003–2006", loc="left"); ax.set_xlabel("PGC ε same day (%)")
    axs[0].set_ylabel("other station ε (%)"); axs[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Same-day agreement with PGC (2–4 s): a shared change would fall on the 1:1 line", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{S.FIG}/s10_polaris_vs_pgc.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    for f in S.FAMS:
        daily_fig(f)
    within_fig(); precision_fig(); polaris_fig()
    print("figures written")
