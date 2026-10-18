"""Q2 figures: band width per episode 2003-2026 by region; duration vs speed; depth thirds; thinning."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
R = f"{ROOT}/svi_science/results/q2_band"
FIG = f"{ROOT}/svi_science/figures"
EP = pd.read_csv(f"{R}/episodes.csv")
U = pd.read_csv(f"{R}/units.csv")
COL = {"SVI-B 2003-13": "#1f77b4", "SVI-O 2014-26": "#17becf", "NW-E 2007-14": "#d62728", "NW-O 2015-26": "#ff7f0e"}
LAB = {"SVI-B 2003-13": "SVI, Bostock catalogue 2003–13", "SVI-O 2014-26": "SVI, ours 2014–26",
       "NW-E 2007-14": "Washington-NW, ours 2007–14", "NW-O 2015-26": "Washington-NW, ours 2015–26"}

# ---------------- fig 1: band width per episode over time
fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
for j, reg in enumerate(["SVI", "NW"]):
    for i, (col, ttl) in enumerate([("W_dur", "band width from duration: |front speed| × median source D80 (km)"),
                                    ("W_snap", "direct: along-strike extent of activity in 12-h windows (km)")]):
        ax = axs[i, j]
        for era, g in EP[EP.region == reg].groupby("era"):
            sf = g[g.single_front] if col == "W_dur" else g
            lo = sf[col] * np.exp(-sf.sig_logW) if col == "W_dur" else None
            hi = sf[col] * np.exp(sf.sig_logW) if col == "W_dur" else None
            if col == "W_dur":
                ax.errorbar(sf.year, sf[col], yerr=[sf[col] - lo, hi - sf[col]], fmt="o", color=COL[era], label=LAB[era], capsize=2)
            else:
                ax.plot(sf.year, sf[col], "o", color=COL[era], label=LAB[era])
            ax.plot(sf.year, sf[f"{col}_ctrl_med"], "x", color=COL[era], alpha=0.5, ms=5)
            if era.startswith("SVI-B"):
                ax.plot(sf.year, sf[f"{col}_thinM"], "s", mfc="none", color=COL[era], ms=6, label="Bostock thinned to our count (largest M)")
            if col == "W_dur":
                ns = g[~g.single_front]
                ax.plot(ns.year, ns[col], "o", mfc="none", color="grey", ms=5)
        ax.set_title(f"{reg}: {ttl}", fontsize=9, loc="left")
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        if i == 1:
            ax.set_xlabel("year of episode")
axs[0, 0].legend(fontsize=7, loc="upper left")
axs[0, 1].legend(fontsize=7, loc="upper left")
fig.text(0.5, 0.005, "x = scrambled control (source positions shuffled among events within the episode, front speed kept); "
         "grey open = episodes without a single front (not tested)", ha="center", fontsize=8)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(f"{FIG}/q2_band_width_by_year.png", dpi=140)

# ---------------- fig 2: duration vs speed
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
sf = EP[EP.single_front]
vv = np.linspace(1.5, 14, 100)
for era, g in sf.groupby("era"):
    axs[0].plot(g.v.abs(), g.medD80, "o", color=COL[era], label=LAB[era])
for w in (10, 15, 20):
    axs[0].plot(vv, w / vv, "k:", lw=0.8); axs[0].text(vv[-1], w / vv[-1], f"{w} km", fontsize=7)
axs[0].set_xlabel("|front speed| (km/day)"); axs[0].set_ylabel("median source active duration D80 (days)")
axs[0].set_title("duration vs speed; dotted = constant band width", fontsize=9, loc="left"); axs[0].legend(fontsize=7)
axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].grid(alpha=0.3)
# depth
for k, era in enumerate(COL):
    g = sf[sf.era == era]
    m = [g[f"W_dur_{c}"].median() for c in ("shallow", "mid", "deep")]
    axs[1].plot([0, 1, 2], m, "o-", color=COL[era], label=LAB[era])
axs[1].set_xticks([0, 1, 2]); axs[1].set_xticklabels(["up-dip third", "middle", "down-dip third"])
axs[1].set_ylabel("median band width (km)"); axs[1].set_title("band width by down-dip position", fontsize=9, loc="left"); axs[1].grid(alpha=0.3)
# unit-level D80 vs dd (within episode, centred)
for era, g in U.dropna(subset=["dd"]).groupby("era"):
    g = g.copy(); g["c"] = g.D80_n30 / g.groupby("ep").D80_n30.transform("median")
    gm = g.groupby(pd.cut(g.dd, 8), observed=True).c.median()
    axs[2].plot([iv.mid for iv in gm.index], gm.values, "o-", color=COL[era], label=LAB[era])
axs[2].axhline(1, color="k", lw=0.6)
axs[2].set_xlabel("down-dip coordinate of source (km, frame of each region)"); axs[2].set_ylabel("source D80 / episode median")
axs[2].set_title("source duration vs down-dip position (within episode)", fontsize=9, loc="left"); axs[2].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(f"{FIG}/q2_duration_vs_speed.png", dpi=140)

# ---------------- fig 3: thinning (duration width uses the 15-event floor: random thinning leaves no source with >= 30)
T = pd.read_csv(f"{R}/thinning.csv"); TP = pd.read_csv(f"{R}/thinning_post.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
b = EP[EP.era == "SVI-B 2003-13"].set_index("ep")
specs = [("W_dur15", TP[TP.draw >= 0].rename(columns={}), TP[TP.draw == -1].set_index("ep")),
         ("W_snap", T[T.kind == "random"], T[T.kind == "top_mag"].set_index("ep"))]
for (c, rnd, tm), ax in zip(specs, axs):
    for ep, g in rnd.groupby("ep"):
        ax.plot([b.loc[ep, c]] * len(g), g[c], ".", color="0.6", ms=3)
    ax.plot(b[c], tm.loc[b.index, c], "s", color="#1f77b4", label="thinned: largest-M events kept")
    ax.plot(b[c], rnd.groupby("ep")[c].median().loc[b.index], "o", color="k", label="thinned at random (median of 50)")
    for ep in b.index[~b.single_front]:
        ax.annotate(f"{ep} (no single front)", (b.loc[ep, c], tm.loc[ep, c]), fontsize=6)
    lim = [0, max(b[c].max(), np.nanmax(rnd[c])) * 1.05]
    ax.plot(lim, lim, "k--", lw=0.7); ax.set_xlabel(f"{c}, full Bostock catalogue (km)"); ax.set_ylabel(f"{c}, thinned to our count (km)")
    ax.set_title(f"Bostock 2003–13 thinned to 639 LFEs per episode ({'15-event floor' if c == 'W_dur15' else 'direct'})", fontsize=9, loc="left")
    ax.grid(alpha=0.3); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/q2_thinning.png", dpi=140)
print("figures written")
