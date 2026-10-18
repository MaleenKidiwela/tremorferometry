#!/usr/bin/env python
"""P3 figures: the Step-1a pair-linkability kill test."""
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
C = dict(null="#6b7280", true="#1f77b4", same="#d62728", near="#ff7f0e", far="#2ca02c")


def centres(lo, hi, n):
    e = np.linspace(lo, hi, n + 1)
    return 0.5 * (e[1:] + e[:-1])


def fig_dist():
    zn = np.load(f"{P3}/pairlink_null_hist.npz")
    zt = np.load(f"{P3}/pairlink_true_hist.npz")
    zc = np.load(f"{P3}/pairlink_true_hist_class.npz")
    jn = json.load(open(f"{P3}/pairlink_null.json"))
    x = centres(float(zn["lo"]), float(zn["hi"]), len(zn["hist"]))
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for h, k, lab in [(zn["hist"], "null", "NULL (stations circularly shifted)"),
                      (zt["hist"], "true", "TRUE (all pairs)")]:
        y = h / h.sum()
        ax[0].step(x, np.where(y > 0, y, np.nan), where="mid", color=C[k], lw=1.4, label=lab)
    for k in ["same", "near", "far"]:
        h = zc[k]
        lab = dict(same="same family", near="different family < 5 km",
                   far="different family > 5 km")[k]
        y = h / h.sum()
        ax[1].step(x, np.where(y > 0, y, np.nan), where="mid", color=C[k], lw=1.4,
                   label=f"{lab}  (n={int(h.sum()):,})")
    for a in ax:
        a.axvline(jn["q_1e-6"], color="k", ls="--", lw=1.0)
        a.text(jn["q_1e-6"], 3e-1, "  null 1e-6", rotation=90, va="top", fontsize=8)
        a.set_yscale("log"); a.set_xlim(0, 32); a.set_ylim(1e-9, 1)
        a.set_xlabel("network pair cc  (sum of 7 site MAD-z, max over $\\pm$3 s lag)")
        a.set_ylabel("fraction of pairs"); a.legend(fontsize=8, frameon=False)
        a.grid(alpha=.25, lw=.5)
    ax[0].set_title("Step 1a — true vs null, 2012 ETS days, 3.8e8 Bostock-event pairs")
    ax[1].set_title("true run, split by the pair's family relation")
    fig.tight_layout()
    fig.savefig(f"{FIG}/q_pairlink_distributions.png", dpi=140)
    print(f"-> {FIG}/q_pairlink_distributions.png")


def fig_family():
    F = pd.read_csv(f"{P3}/pairlink_true_by_family_thresholds.csv")
    ho = json.load(open(f"{P.RES}/holdouts.json"))
    st = {h["bkey"]: h["stratum"] for h in ho["hidden_families"]}
    F["stratum"] = F.bkey.map(st).fillna("kept")
    order = ["isolated_ge5km", "mid_2.5_5km", "close_lt2.5km", "kept"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    cols = {"q_0.001": "#9ecae1", "q_0.0001": "#4292c6", "q_1e-5": "#2171b5", "q_1e-6": "#08306b"}
    Fi = F[F.stratum == "isolated_ge5km"].sort_values("bkey")
    w = 0.2
    for i, (c, col) in enumerate(cols.items()):
        ax[0].bar(np.arange(len(Fi)) + (i - 1.5) * w, Fi[c], width=w, color=col,
                  label=f"null {c.replace('q_','')} quantile")
    ax[0].axhline(30, color="#d62728", ls="--", lw=1.2)
    ax[0].text(0.02, 30 * 1.15, "Merlin's bar: 30 pairs", color="#d62728", fontsize=8)
    ax[0].set_xticks(np.arange(len(Fi))); ax[0].set_xticklabels(Fi.bkey.astype(str))
    ax[0].set_yscale("symlog", linthresh=1)
    ax[0].set_xlabel("hidden family (isolated stratum, NN $\\geq$ 5 km)")
    ax[0].set_ylabel("same-family pairs above threshold")
    ax[0].set_title("Step 1a readout — 0/10 reach the bar at the 1e-6 threshold")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=.25, lw=.5, axis="y")

    # enrichment of same-family over far-field cross-family, by threshold
    j = json.load(open(f"{P3}/pairlink_true.json"))
    ks = list(cols)
    n_same = F.n_same_pairs.sum()
    n_far = j["far_tot"]
    e = [(j["same_hit_by_thr"][k] / n_same) / max(1e-15, j["far_hit_by_thr"][k] / n_far)
         for k in ks]
    ax[1].bar(range(len(ks)), e, color=[cols[k] for k in ks])
    for i, (k, v) in enumerate(zip(ks, e)):
        ax[1].text(i, v * 1.03, f"{v:.0f}$\\times$", ha="center", fontsize=9)
    ax[1].axhline(1, color="k", lw=1)
    ax[1].set_xticks(range(len(ks)))
    ax[1].set_xticklabels([k.replace("q_", "") for k in ks])
    ax[1].set_xlabel("threshold (null quantile)")
    ax[1].set_ylabel("same-family / far-cross-family rate above threshold")
    ax[1].set_title("the linkage signal is REAL but tiny\n"
                    "(0.6 % of same-family pairs reach even the 1e-3 quantile)")
    ax[1].grid(alpha=.25, lw=.5, axis="y")
    fig.tight_layout()
    fig.savefig(f"{FIG}/q_pairlink_by_family.png", dpi=140)
    print(f"-> {FIG}/q_pairlink_by_family.png")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_dist()
    fig_family()
