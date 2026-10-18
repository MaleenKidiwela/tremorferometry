#!/usr/bin/env python
"""P3 figure: the pooled (2010 + 2012) Step-1a repeat and the site-subset ablation."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
FIG = f"{P.PICK}/figures"
THRN = ["q_0.001", "q_0.0001", "q_1e-5", "q_1e-6"]
COLS = {"q_0.001": "#9ecae1", "q_0.0001": "#4292c6", "q_1e-5": "#2171b5", "q_1e-6": "#08306b"}


def main():
    os.makedirs(FIG, exist_ok=True)
    R = pd.read_csv(f"{P3}/pool_isolated_hidden.csv")
    S = json.load(open(f"{P3}/pool_summary.json"))
    AB = pd.read_csv(f"{P3}/pool_site_ablation.csv")
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))

    # (a) per isolated hidden family, pooled vs 2012-only, at the 1e-6 and 1e-5 quantiles
    x = np.arange(len(R))
    w = 0.2
    ax[0].bar(x - 1.5 * w, R["q_0.0001"], w, color=COLS["q_0.0001"], label="pooled, null 1e-4")
    ax[0].bar(x - 0.5 * w, R["q_1e-5"], w, color=COLS["q_1e-5"], label="pooled, null 1e-5")
    ax[0].bar(x + 0.5 * w, R["q_1e-6"], w, color=COLS["q_1e-6"], label="pooled, null 1e-6")
    ax[0].bar(x + 1.5 * w, R["q_1e-6_2012only"], w, color="#bdbdbd",
              label="2012 only, null 1e-6")
    ax[0].axhline(30, color="#d62728", ls="--", lw=1.2)
    ax[0].text(0.02, 33, "Merlin's bar: 30 pairs", color="#d62728", fontsize=8)
    ax[0].set_xticks(x); ax[0].set_xticklabels(R.bkey.astype(str))
    ax[0].set_yscale("symlog", linthresh=1)
    ax[0].set_xlabel("hidden family (isolated stratum)")
    ax[0].set_ylabel("same-family pairs above threshold")
    ax[0].set_title(f"pooled 2010 + 2012 — {S['isolated_at_bar']}/10 reach the bar\n"
                    f"({S['pooled_events']:,} events, {S['pooled_pairs']/1e6:.0f}M pairs)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=.25, lw=.5, axis="y")

    # (b) site-subset ablation: same-family rate above each threshold
    lbl = {"protocol, 7 sites": "7 discovery sites\n(protocol)",
           "4 nearest discovery sites": "4 nearest of those 7\n(oracle site choice)",
           "oracle 4 incl. withheld B926": "PGCq+B926+VGZ+nearest\nfar borehole (oracle+withheld)"}
    xs = np.arange(len(AB))
    for i, k in enumerate(THRN):
        r = AB[f"iso_{k}"] / AB["iso_same_pairs"].clip(lower=1)
        ax[1].bar(xs + (i - 1.5) * 0.2, r, 0.2, color=COLS[k], label=k.replace("q_", "null "))
    ax[1].set_xticks(xs); ax[1].set_xticklabels([lbl[v] for v in AB.row], fontsize=8)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("fraction of isolated-hidden same-family pairs above threshold")
    ax[1].set_title("do the far sites dilute the pair statistic?")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25, lw=.5, axis="y")

    # (c) the site-count extrapolation
    g = S["site_scaling"]
    ns = sorted(int(k) for k in g)
    need = [float(g[str(n)]["pairs_for_30"]) for n in ns]
    ev = [np.sqrt(2 * v) for v in need]
    ax[2].plot(ns, ev, "o-", color="#08306b")
    ax[2].axhline(561, color="#d62728", ls="--", lw=1.2)
    ax[2].text(ns[0], 620, "biggest 2012 family: 561 events", color="#d62728", fontsize=8)
    ax[2].axhline(S["pooled_events"] / 130, color="#2ca02c", ls=":", lw=1.2)
    ax[2].set_yscale("log")
    ax[2].set_xlabel("number of independent sites in the network sum")
    ax[2].set_ylabel("events of ONE family needed for 30 linkable pairs")
    ax[2].set_title("$\\sqrt{N_{sites}}$ scaling of the network sum\n"
                    f"(same-family shift {S['delta_sigma_same_family']:.2f}$\\sigma$ at 7 sites)")
    ax[2].grid(alpha=.25, lw=.5)
    fig.tight_layout()
    fig.savefig(f"{FIG}/q_pairlink_pooled.png", dpi=140)
    print(f"-> {FIG}/q_pairlink_pooled.png")


if __name__ == "__main__":
    main()
