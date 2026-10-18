#!/usr/bin/env python
"""Figure for §6 Amendment 4: envelope width vs distance (LFE stacks, earthquakes), distance-limited κ, synthetic bias.
  python s11_eq_swidth_fig.py  -> tomo/eq/figures/fig_swidth_kappa.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EQ = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/tomo/eq"
L = pd.read_csv(f"{EQ}/tables/swidth2_lfe_rows.csv"); Q = pd.read_csv(f"{EQ}/tables/swidth2_eq_rows.csv")
K = pd.read_csv(f"{EQ}/tables/swidth2_dist_kappa.csv"); S = pd.read_csv(f"{EQ}/tables/synthw_kappa.csv")
COL = {"WA_NW": "#1f77b4", "WA_SW": "#d62728", "SVI": "#2ca02c"}
EQK = {"WA_NW": 1.775, "WA_SW": 1.735, "SVI": 1.755}
bins = np.arange(0, 101, 10); mid = bins[:-1] + 5


def binned(d, col):
    g = d.groupby(pd.cut(d.epi_km, bins), observed=False)[col]
    m, n = g.median().to_numpy(copy=True), g.size().to_numpy()
    m[n < 8] = np.nan
    return m


fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
for reg, c in COL.items():
    l = L[(L.reg == reg) & (L.rS >= 1.2)]; q = Q[(Q.reg == reg) & (Q.rS >= 1.2)]
    ax[0, 0].plot(mid, binned(l, "wS"), "-o", color=c, label=f"{reg} LFE stacks")
    ax[0, 0].plot(mid, binned(q, "wS"), "--s", color=c, mfc="none", label=f"{reg} earthquakes")
    l2 = L[(L.reg == reg) & (L.rS >= 1.2) & (L.rP >= 1.2)].copy(); l2["d"] = l2.wS - l2.wP
    q2 = Q[(Q.reg == reg) & (Q.rS >= 1.2) & (Q.rPv3 >= 1.2)].copy(); q2["d"] = q2.wS - q2.wP
    ax[0, 1].plot(mid, binned(l2, "d"), "-o", color=c, label=f"{reg} LFE stacks")
    ax[0, 1].plot(mid, binned(q2, "d"), "--s", color=c, mfc="none", label=f"{reg} earthquakes")
ax[0, 0].set(xlabel="epicentral distance (km)", ylabel="S envelope width, FWHM (s), 1–8 Hz",
             title="S envelope width vs distance (bin medians, ≥ 8 rows)")
ax[0, 1].set(xlabel="epicentral distance (km)", ylabel="S width − P width (s)", title="S minus P width vs distance")
ax[0, 0].legend(fontsize=7, ncol=2); ax[0, 1].legend(fontsize=7, ncol=2)
lab = {40: "≤ 40", 50: "≤ 50", 1000: "all"}
for i, (reg, c) in enumerate(COL.items()):
    k = K[(K.reg == reg) & (K.dist == "epi") & (K.nmin == np.where(K.dmax == 1000, 4, 3))]
    x = np.arange(len(k)) + 0.2 * (i - 1)
    ax[1, 0].errorbar(x, k.kappa, yerr=k.sig, fmt="o", color=c, label=f"{reg} LFE κ")
    for xx, r in zip(x, k.itertuples()):
        if np.isfinite(r.scr_kappa_med):
            ax[1, 0].plot([xx + 0.07] * 2, [r.scr_kappa_p16, r.scr_kappa_p84], color=c, alpha=0.4, lw=5)
    ax[1, 0].axhline(EQK[reg], color=c, ls=":", lw=1)
ax[1, 0].set_xticks(range(3)); ax[1, 0].set_xticklabels(["≤ 40 km", "≤ 50 km", "all stations"])
ax[1, 0].set(ylabel="Wadati κ (clipped tP-on-tS)", xlabel="epicentral distance limit (family ≥ 3 rows; all: ≥ 4)",
             title="LFE κ by distance limit")
ax[1, 0].text(0.02, 0.02, "shaded: scrambled control 16–84 %\ndotted: near-interface earthquake κ", transform=ax[1, 0].transAxes, fontsize=7)
ax[1, 0].legend(fontsize=7)
s = S[(S.draw > 0) & (S.kappa_true == 1.75) & (S.picker == "V3")].groupby(["reg", "variant"]).bias.mean().unstack(0)
order = ["W0", "W1", "W2", "W3", "W2s"]
s = s.reindex(order)
xx = np.arange(len(order))
for j, reg in enumerate(["WA_NW", "WA_SW"]):
    ax[1, 1].bar(xx + 0.35 * (j - 0.5), s[reg], 0.35, color=COL[reg], label=reg)
ax[1, 1].axhline(0.17, color="k", ls="--", lw=1); ax[1, 1].text(0, 0.172, "observed WA excess +0.17", fontsize=7, va="bottom")
ax[1, 1].axhline(0.10, color="gray", ls=":", lw=1)
ax[1, 1].set_xticks(xx); ax[1, 1].set_xticklabels(["W0 const", "W1 trend", "W2 own", "W3 3×", "W2s scr."])
ax[1, 1].set(ylabel="recovered κ − true κ (true 1.75)", title="Synthetic κ bias, measured broadening")
ax[1, 1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{EQ}/figures/fig_swidth_kappa.png", dpi=130)
print("ok")
