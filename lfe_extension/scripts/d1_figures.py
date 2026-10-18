#!/usr/bin/env python
"""Track D, D1 -- the completeness-model figures: (a) the out-of-sample validation (predicted vs
measured own recall by magnitude bin on the held-out 2010 year), (b) per-family recall by epoch,
(c) the Mc / day-term distributions over the production record, (d) the false-alarm floor by season."""
from __future__ import annotations
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
FIG = f"{EXT}/figures"


def main():
    os.makedirs(FIG, exist_ok=True)
    fig = plt.figure(figsize=(13.5, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.28)

    # (a) validation
    ax = fig.add_subplot(gs[0, 0])
    V = pd.read_csv(f"{TD}/mc_model_validation.csv")
    Vh = V[~V["bin"].str.contains("in-sample")]
    for tier, c, m in (("T_core", "tab:blue", "o"), ("T_ext", "tab:red", "s")):
        g = Vh[Vh.tier == tier]
        ax.plot(g.measured, g.predicted, m, color=c, ms=7, label=tier)
        for lab, xx, yy in zip(g["bin"], g.measured, g.predicted):
            ax.annotate(lab, (xx, yy), fontsize=6, xytext=(3, 3), textcoords="offset points")
    lim = [0, max(0.35, Vh[["measured", "predicted"]].to_numpy().max() * 1.15)]
    ax.plot(lim, lim, "k-", lw=0.8)
    ax.fill_between(lim, [l - 0.05 for l in lim], [l + 0.05 for l in lim], color="k", alpha=0.08,
                    label="bar ±0.05")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("measured own recall, held-out 2010"); ax.set_ylabel("predicted")
    ax.set_title("(a) completeness model, out of sample", fontsize=10); ax.legend(fontsize=7)
    ax.grid(alpha=0.25)

    # (b) per-family recall by epoch
    ax = fig.add_subplot(gs[0, 1])
    C = pd.read_csv(f"{TD}/completeness_family_epoch.csv")
    eps = [e for e in ("E1", "E2", "E3") if e in set(C.epoch_id)]
    for i, tier in enumerate(("T_core", "T_ext")):
        d = [C[(C.epoch_id == e) & (C.tier == tier)].recall.dropna() for e in eps]
        bp = ax.boxplot(d, positions=np.arange(len(eps)) + (i - 0.5) * 0.3, widths=0.26,
                        patch_artist=True, showfliers=False)
        for b in bp["boxes"]:
            b.set_facecolor(["tab:blue", "tab:red"][i]); b.set_alpha(0.55)
    ax.set_xticks(range(len(eps))); ax.set_xticklabels(eps)
    ax.set_ylabel("per-family own recall (2012 test half)")
    ax.set_title("(b) station-set ablation by epoch\n(blue T_core, red T_ext)", fontsize=10)
    ax.grid(alpha=0.25)

    # (c) recall vs number of template stations
    ax = fig.add_subplot(gs[0, 2])
    g = C[(C.epoch_id == "E1") & (C.tier == "T_core")]
    ax.scatter(g.n_template_sta, g.recall, s=16, c="tab:blue", alpha=0.6)
    ax.set_xlabel("stations with a usable 2012 template")
    ax.set_ylabel("own recall, E1, T_core")
    ax.set_title("(c) family-dependence of completeness", fontsize=10); ax.grid(alpha=0.25)

    # (d) Mc distribution
    ax = fig.add_subplot(gs[1, 0])
    p = f"{TD}/mc_family_day.parquet"
    if os.path.exists(p):
        M = pd.read_parquet(p)
        for tier, c in (("Mc_T_core", "tab:blue"), ("Mc_T_ext", "tab:red")):
            v = M[tier].replace([np.inf, -np.inf], np.nan).dropna()
            v = v[(v > 0) & (v < 6)]
            ax.hist(v, bins=60, histtype="step", color=c, lw=1.6, label=tier)
        ax.set_xlabel("Mc(family, day)"); ax.set_ylabel("family-days")
        ax.set_title("(d) per-family-day completeness magnitude", fontsize=10)
        ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # (e) day terms
    ax = fig.add_subplot(gs[1, 1])
    if os.path.exists(p):
        ax.hist(M.u.dropna(), bins=60, histtype="step", color="k", lw=1.5, label="u (MAD scale)")
        ax.hist(M.v.dropna(), bins=60, histtype="step", color="tab:green", lw=1.5, label="v (noise)")
        ax.set_xlabel("log10 day-term anomaly"); ax.set_ylabel("family-days")
        ax.set_title("(e) day terms of the sensitivity model", fontsize=10)
        ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # (f) floor by season
    ax = fig.add_subplot(gs[1, 2])
    S = pd.read_csv(f"{TD}/floor_epoch_season_summary.csv")
    order = ["DJF", "MAM", "JJA", "SON"]
    for tier, c in (("T_core", "tab:blue"), ("T_ext", "tab:red")):
        for eid, ls in zip(sorted(S.epoch_id.unique()), ["-", "--", ":"]):
            g = S[(S.tier == tier) & (S.epoch_id == eid)].set_index("season").reindex(order)
            ax.plot(order, g.lambda_pooled, ls, marker="o", color=c, label=f"{tier} {eid}")
    ax.set_ylabel("λ (false alarms / template / quiet day)")
    ax.set_title("(f) false-alarm floor by season", fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.25)

    fig.suptitle("Track D — completeness model (Merlin Step 6)", fontsize=13)
    fig.savefig(f"{FIG}/d1_completeness_model.png", dpi=135, bbox_inches="tight")
    print(f"[write] {FIG}/d1_completeness_model.png")


if __name__ == "__main__":
    main()
