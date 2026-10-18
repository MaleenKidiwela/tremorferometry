"""Q6 figures: front extents along the margin over time; stop probability vs time since last activity; Strait tests."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q6_front_arrest"
FIG = f"{ROOT}/svi_science/figures"
C_NW, C_SE, C_NONE = "#2a78d6", "#eb6834", "#8a8984"
INK, INK2, GRID, STRAIT_BG = "#0b0b0b", "#52514e", "#e4e3df", "#f3e7c9"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 9.5, "figure.facecolor": "white"})
T0 = pd.Timestamp("2010-01-01")
PLACES = [("Victoria", 9.5), ("Port Angeles", -16.7), ("Lake Cowichan", 75), ("Port Alberni", 146), ("Hood Canal S", -110),
          ("Olympia", -147), ("Columbia R.", -243)]


def fronts():
    E = pd.read_csv(f"{OUT}/episodes_pnsn.csv"); B = E[E.big].copy()
    B["date"] = T0 + pd.to_timedelta(B.t_first, unit="D")
    fig, axs = plt.subplots(2, 1, figsize=(12, 10), constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 1.15]))
    for ax, (lo, hi), ttl in ((axs[0], (-970, 330), "Whole margin (PNSN tremor, 2010-01 → 2026-05): every big episode"),
                              (axs[1], (-260, 260), "Vancouver Island + Washington (zoom; Strait of Juan de Fuca shaded)")):
        ax.axhspan(-15, 15, color=STRAIT_BG, zorder=0, lw=0)
        sub = B[(B.stop_nw_s >= lo) & (B.stop_se_s <= hi)]
        sz = np.log10(sub.n.values)
        lw = 0.8 + 2.2 * (sz - 2.3) / (4.2 - 2.3)
        for (_, r), w in zip(sub.iterrows(), lw):
            ax.plot([r.date, r.date], [r.stop_se_s, r.stop_nw_s], color=C_NONE, lw=w, solid_capstyle="butt", zorder=1)
        ax.scatter(sub.date, sub.start_s, s=14, marker="o", color=INK, zorder=3, lw=0, label="start (first 10 % of epicentres)")
        ax.scatter(sub.date, sub.stop_nw_s, s=16, marker="v", color=C_NW, zorder=2, lw=0, label="stop, NW end")
        ax.scatter(sub.date, sub.stop_se_s, s=16, marker="^", color=C_SE, zorder=2, lw=0, label="stop, SE end")
        ax.set_ylim(lo, hi); ax.set_ylabel("along-margin km (0 = mid-Strait; + toward NW)")
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_title(ttl, loc="left")
        for nm, s in PLACES:
            if lo < s < hi:
                ax.annotate(nm, (1.0, s), xycoords=("axes fraction", "data"), xytext=(4, 0), textcoords="offset points",
                            va="center", fontsize=7.5, color=INK2, annotation_clip=False)
    axs[0].annotate("N. California", (1.0, -900), xycoords=("axes fraction", "data"), xytext=(4, 0), textcoords="offset points", fontsize=7.5, color=INK2)
    axs[0].annotate("Oregon", (1.0, -550), xycoords=("axes fraction", "data"), xytext=(4, 0), textcoords="offset points", fontsize=7.5, color=INK2)
    axs[0].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False, ncol=3)
    axs[1].set_xlabel("episode start (PT)")
    fig.suptitle("Q6 — where each big tremor episode started and where its fronts stopped (grey bar = front run; thicker = more tremor)",
                 x=0.01, ha="left", fontsize=10)
    fig.savefig(f"{FIG}/q6_fronts_margin.png", dpi=130); plt.close(fig)


def arrest():
    t = json.load(open(f"{OUT}/tests.json"))
    q = t["q61_stopfrac_quintiles"]; tt = t["q61"]
    W = pd.read_csv(f"{OUT}/walk_rows_pnsn.csv")
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    ax = axs[0]
    lt = np.log10(W.T_any[W.T_any > 0]); ed = [lt.min()] + q["edges_log10T"] + [lt.max()]
    mids = [10 ** ((ed[i] + ed[i + 1]) / 2) for i in range(5)]
    ax.fill_between(mids, 100 * np.array(q["null_q05"]), 100 * np.array(q["null_q95"]), color=GRID, lw=0,
                    label="scrambled control (episode order shuffled), 5–95 %")
    ax.plot(mids, 100 * np.array(q["null_mean"]), color=C_NONE, lw=1.5, ls="--")
    ax.plot(mids, 100 * np.array(q["real"]), color=C_NW, lw=2, marker="o", ms=7, label="real")
    for x, y in zip(mids, q["real"]):
        ax.annotate(f"{100 * y:.1f} %", (x, 100 * y), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8, color=INK)
    ax.set_xscale("log"); ax.set_xlabel("days since that 5-km stretch last had tremor (quintile medians)")
    ax.set_ylabel("chance the front stops at the next 5-km bin (%)")
    ax.set_title("Fronts stop more often where the fault slipped recently", loc="left")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax = axs[1]
    names = ["walk_Tany_full", "walk_Tbig_full", "walk_Tany_north", "walk_Tany7_full", "walk_Tany_E12", "walk_Tany_E3",
             "enter_Tany_full", "enter_Tbig_full", "enter_Tany_north"]
    labels = ["stop hazard, any episode (PRIMARY)", "stop hazard, since last big ETS", "stop hazard, Vancouver Is.+WA only",
              "stop hazard, ignore last 7 days", "stop hazard, 2010–2014", "stop hazard, 2015–2026",
              "enter vs not, any episode", "enter vs not, last big ETS", "enter vs not, VI+WA only"]
    y = np.arange(len(names))[::-1]
    for yi, n in zip(y, names):
        r = tt[n]
        ax.plot([r["null_q05"], r["null_q95"]], [yi, yi], color=C_NONE, lw=5, alpha=0.5, solid_capstyle="butt")
        ax.plot(r["beta"], yi, "o", color=C_NW, ms=7)
        ax.annotate(f"p {r['p_pos']:.4f}" if r["p_pos"] >= 0.001 else "p < 0.001", (r["beta"], yi), xytext=(8, -3), textcoords="offset points", fontsize=7.5, color=INK2)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("effect of log(time since last tremor) on 'front continues'\n(logit per 1 SD; grey = scrambled control 5–95 %)")
    ax.set_title("Same sign in every version", loc="left"); ax.grid(color=GRID, lw=0.6, axis="x"); ax.set_axisbelow(True)
    fig.savefig(f"{FIG}/q6_arrest_vs_tsince.png", dpi=130); plt.close(fig)


def strait():
    t = json.load(open(f"{OUT}/tests.json"))["q62"]
    E = pd.read_csv(f"{OUT}/episodes_pnsn.csv"); B = E[E.big]
    ep = pd.read_parquet(f"{OUT}/epicentres.parquet", columns=["s_km"])
    s = ep.s_km.values; s = s[(s >= -200) & (s <= 200)]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    bins = np.arange(-200, 201, 15)
    for ax, vals, n, ttl, key in ((axs[0], B.start_s[(B.start_s >= -200) & (B.start_s <= 200)], None, "Starts of big episodes", "starts_big_north_hw15"),
                                  (axs[1], np.r_[B.stop_se_s[(B.start_s.between(-200, 200)) & B.stop_se_s.between(-200, 200)],
                                                 B.stop_nw_s[(B.start_s.between(-200, 200)) & B.stop_nw_s.between(-200, 200)]], None,
                                   "Front stops of big episodes", "stops_big_north_hw15")):
        ax.axvspan(-15, 15, color=STRAIT_BG, lw=0, zorder=0)
        h, _ = np.histogram(vals, bins)
        ax.bar(bins[:-1] + 7.5, h, width=13, color=C_NW, zorder=2, label="real")
        dens, _ = np.histogram(s, bins); dens = dens / dens.sum() * len(vals)
        ax.step(bins, np.r_[dens, dens[-1]], where="post", color=C_SE, lw=1.6, zorder=3, label="matched to where tremor is")
        ax.axhline(len(vals) / (len(bins) - 1), color=C_NONE, ls="--", lw=1.2, label="uniform", zorder=3)
        r = t[key]
        ax.set_title(f"{ttl}: {r['k_strait']} in Strait zone vs {r['uniform_expect']:.0f} uniform, {r['matched_expect']:.0f} matched",
                     loc="left", fontsize=8.8)
        ax.set_xlabel("along-margin km (0 = mid-Strait)"); ax.set_ylabel("count per 15 km")
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.legend(fontsize=7.5, frameon=False)
    CZ = pd.read_csv(f"{OUT}/crossing_zones.csv")
    ax = axs[2]
    ax.axvspan(-15, 15, color=STRAIT_BG, lw=0, zorder=0)
    m = (CZ.z0 + 15)
    ax.plot(CZ.z0 + 15, 100 * CZ.rate, color=C_NW, lw=1.8, marker="o", ms=3)
    c = t["crossing"]
    ax.axhline(100 * c["tiles_rate"], color=C_NONE, ls="--", lw=1.2)
    ax.annotate(f"Strait: {c['strait_crossed']}/{c['strait_arrived']} crossed ({100 * c['strait_rate']:.0f} %)\n"
                f"other 30-km stretches: {100 * c['tiles_rate']:.0f} % (Fisher p {c['fisher_p']:.2f})", (0.02, 0.04), xycoords="axes fraction", fontsize=8)
    ax.set_xlabel("centre of 30-km stretch (km)"); ax.set_ylabel("fronts that arrived and crossed (%)")
    ax.set_title("Do fronts that reach a stretch get across it?", loc="left"); ax.set_ylim(0, 100)
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    fig.suptitle("Q6 — Strait of Juan de Fuca (shaded, ±15 km): PNSN tremor 2010–2026, Vancouver Island + Washington (−200…+200 km)",
                 x=0.01, ha="left", fontsize=10)
    fig.savefig(f"{FIG}/q6_strait.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    fronts(); arrest(); strait()
