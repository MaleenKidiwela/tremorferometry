#!/usr/bin/env python
"""Track D §3b-2 — the summary figure for the purity-matched front estimator: front speed per SSE,
both tiers, with bootstrap CIs, the position-shuffle significance, and the activation yield."""
from __future__ import annotations
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
FIG = f"{EXT}/figures"
BOSTOCK = -4.694908


def main():
    os.makedirs(FIG, exist_ok=True)
    R = pd.read_csv(f"{TD}/thomas_front_v2.csv").sort_values(["t0", "tier"])
    order = list(dict.fromkeys(R.sse_id))
    xi = {s: i for i, s in enumerate(order)}
    fig, ax = plt.subplots(3, 1, figsize=(13.5, 10.5), sharex=True,
                           gridspec_kw=dict(height_ratios=[2.2, 2.2, 1.2]))
    for k, (tier, ttl) in enumerate((("cat_core", "cat_core"), ("cat_ext", "cat_ext"))):
        g = R[R.tier == tier].copy()
        g["xi"] = g.sse_id.map(xi)
        f = g[g.speed_kmday.notna()]
        sig = f[f.plausible & (f.null_xshuf_p <= 0.05)]
        wk = f[~(f.plausible & (f.null_xshuf_p <= 0.05))]
        ax[k].axhline(0, color="k", lw=0.9)
        ax[k].axhline(BOSTOCK, color="tab:blue", ls="--", lw=1.1,
                      label=f"Bostock 2012 main front ({BOSTOCK:.2f} km/day)")
        ax[k].axhspan(BOSTOCK * 1.2, BOSTOCK * 0.8, color="tab:blue", alpha=0.10, label="±20 %")
        for d, col, lab, ms in ((wk, "0.65", "fitted, not significant", 5),
                                (sig, "tab:red", "plausible + position-shuffle p ≤ 0.05", 7)):
            if not len(d):
                continue
            ax[k].errorbar(d.xi, d.speed_kmday,
                           yerr=[np.clip(d.speed_kmday - d.ci_lo, 0, None),
                                 np.clip(d.ci_hi - d.speed_kmday, 0, None)],
                           fmt="o", ms=ms, color=col, ecolor=col, elinewidth=1.1, capsize=2.5,
                           lw=0, label=lab)
        miss = g[g.speed_kmday.isna()]
        if len(miss):
            ax[k].plot(miss.xi, np.zeros(len(miss)), "kx", ms=7, label="no fit")
        ax[k].set_ylim(-14, 14)
        ax[k].set_ylabel("front speed (km/day)\n− = SE (southward)")
        ax[k].set_title(f"({'ab'[k]}) purity-matched main front, {ttl} — "
                        f"{len(f)}/{len(g)} windows fitted, {len(sig)} significant", fontsize=11)
        ax[k].legend(fontsize=8, ncol=2, loc="upper left")
        ax[k].grid(alpha=0.25)
    g = R[R.tier == "cat_core"].copy(); g["xi"] = g.sse_id.map(xi)
    ge = R[R.tier == "cat_ext"].copy(); ge["xi"] = ge.sse_id.map(xi)
    ax[2].bar(g.xi - 0.2, g.n_fam_activated.fillna(0), width=0.4, color="tab:purple", label="cat_core")
    ax[2].bar(ge.xi + 0.2, ge.n_fam_activated.fillna(0), width=0.4, color="tab:orange", label="cat_ext")
    ax[2].axhline(5, color="k", ls=":", lw=1, label="MIN_FAM = 5")
    ax[2].set_ylabel("families activated\nabove their own λ")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.25)
    ax[2].set_xticks(range(len(order)))
    ax[2].set_xticklabels([s + ("*" if bool(R[(R.sse_id == s)].major.iloc[0]) else "")
                           for s in order], rotation=70, fontsize=7)
    ax[2].set_xlabel("SSE window (* = major, ≥ 10 active tremor days)")
    fig.suptitle("Track D §3b-2 — purity-matched main-front estimator across every SSE, 2014–2026",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{FIG}/d3b2_front_v2_summary.png", dpi=135)
    plt.close(fig)
    print(f"[write] {FIG}/d3b2_front_v2_summary.png")


if __name__ == "__main__":
    main()
