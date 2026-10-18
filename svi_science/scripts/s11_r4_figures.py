#!/usr/bin/env python
"""s11_r4_figures -- S11 round 4 figures (PLAN.md §8). User rules: dv/v points coloured by weight or CC, no error bars on
dv/v (sigma in the tables), line graphs.
  fig_r4_dD.png      DeltaD per episode, headline combined model (NOM), SVI and WA-NW; majors outlined; composite line
  fig_r4_suite.png   injection rules per generating kernel (T1 deep-pulse recovery, T2 shallow leak, T5 annual -> composite)
  fig_r4_lapse.png   lapse-selection confusion matrices (N1, N2) per region
  fig_r4_2019_11.png the 2019-11 dive: carriers by station, subsets, sub-episode split, realignment lag / residual shift
    python s11_r4_figures.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
B = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/round4"
FIG = f"{B}/figures"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                     "grid.linewidth": 0.5, "font.size": 9})
REG = {"svi": "SVI", "wa_nw": "WA-NW"}
V6 = ["ell10", "ell20", "ell40", "ell80", "two5_40", "two3_20", "AVG6"]


def fig_dD():
    fig, ax = plt.subplots(2, 1, figsize=(15, 8.5))
    for a, (sub, name) in zip(ax, REG.items()):
        p = f"{B}/{sub}/tables/r4_dD_H.csv"
        if not os.path.exists(p):
            a.set_title(f"{name}: not run"); continue
        T = pd.read_csv(p); T = T[T.label == "H_ell40"].copy(); T["mid"] = pd.to_datetime(T["mid"]); T = T.sort_values("mid")
        w = np.log10(1 / T.sig ** 2)
        a.plot(T.mid, T.dD, "-", color="#9a9a96", lw=0.8, zorder=1)
        sc = a.scatter(T.mid, T.dD, c=w, cmap="viridis", s=16, lw=0, zorder=2)
        M = T[T.major]
        a.scatter(M.mid, M.dD, c=np.log10(1 / M.sig ** 2), cmap="viridis", vmin=w.min(), vmax=w.max(), s=70, edgecolor="k", lw=1, zorder=3)
        for r in M.itertuples():
            a.annotate(r.episode, (r.mid, r.dD), xytext=(0, 7), textcoords="offset points", fontsize=6.5, ha="center")
        a.axhline(0, color="grey", lw=0.6, ls=":")
        a.axhline(M.dD.mean(), color="#eb6834", lw=1.0, ls="--", label=f"composite of {len(M)} majors {M.dD.mean():+.3f} %")
        a.legend(fontsize=8, loc="upper left")
        a.set(ylabel="ΔD (%)", title=f"{name}: deep-weighted ΔD per episode, combined model (near + far B2, Sc, receiver, M/B, sum-to-zero); "
                                         f"σ per major median {M.sig.median():.3f} % (tables), not drawn")
        cb = fig.colorbar(sc, ax=a, pad=0.01); cb.set_label("log10 weight 1/σ²")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r4_dD.png", dpi=130); plt.close(fig)


def fig_suite():
    fig, ax = plt.subplots(2, 3, figsize=(15, 7.5))
    for i, (sub, name) in enumerate(REG.items()):
        p = f"{B}/{sub}/tables/r4_rules.csv"
        if not os.path.exists(p):
            continue
        T = pd.read_csv(p); T = T[T.label == "H"]
        for j, (col, lo, hi, lab) in enumerate((("T1_Pd_rec", 0.8, 1.1, "T1 deep-pulse recovery"), ("T2_Ps_leak", 0, 0.35, "T2 common shallow-pulse leak"),
                                                 ("T5_B_comp", 0, 0.35, "T5 shallow annual → composite"))):
            a = ax[i, j]
            a.axhspan(lo, hi, color="#1baf7a", alpha=0.12, lw=0)
            for k, (setup, g) in enumerate(T.groupby("setup")):
                g = g.set_index("gen").reindex(V6)
                a.plot(range(len(V6)), g[col], "-o", ms=4, lw=1.2, color=["#2a78d6", "#eb6834"][k % 2], label=f"inverted with {setup}")
            a.set_xticks(range(len(V6))); a.set_xticklabels(V6, rotation=30, fontsize=7)
            a.set(title=f"{name}: {lab}", xlabel="generating kernel")
            if j == 0:
                a.legend(fontsize=7)
    fig.suptitle("Waveform-level injections through the combined model (green band = pass rule)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r4_suite.png", dpi=130); plt.close(fig)


def fig_lapse():
    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    for i, (sub, name) in enumerate(REG.items()):
        p = f"{B}/{sub}/tables/lapse_confusion.csv"
        if not os.path.exists(p):
            continue
        CF = pd.read_csv(p); vj = json.load(open(f"{B}/{sub}/tables/lapse_verdict.json"))
        for j, noise in enumerate(("N1", "N2")):
            M = CF[(CF.noise == noise) & CF.sel.isin(V6[:6])].pivot(index="true", columns="sel", values="p").reindex(index=V6[:6], columns=V6[:6]).fillna(0)
            a = ax[i, j]; im = a.imshow(M.to_numpy(), cmap="Blues", vmin=0, vmax=1)
            for (r, c), x in np.ndenumerate(M.to_numpy()):
                a.text(c, r, f"{x:.2f}", ha="center", va="center", fontsize=7, color="k" if x < 0.6 else "w")
            a.set_xticks(range(6)); a.set_xticklabels(V6[:6], rotation=30, fontsize=7); a.set_yticks(range(6)); a.set_yticklabels(V6[:6], fontsize=7)
            a.set(xlabel="kernel chosen by the lapse fit", ylabel="kernel that made the data",
                  title=f"{name}, noise {noise}: data choose {vj['vstar']}; adopted {vj.get('adopted')}")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r4_lapse.png", dpi=130); plt.close(fig)


def fig_dive():
    T = f"{B}/svi/tables"
    if not os.path.exists(f"{T}/dive_verdict.json"):
        return
    V = json.load(open(f"{T}/dive_verdict.json"))
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    CA = pd.read_csv(f"{T}/dive_carriers.csv"); st = CA[CA.by == "netsta"].sort_values("contrib")
    a = ax[0, 0]; a.barh(st.key, st.contrib, color=np.where(st.contrib > 0, "#2a78d6", "#eb6834"))
    a.set(xlabel="contribution to ΔD(2019-11) (%)", title=f"carriers by station (ΔD {V['H']['dD']:+.3f} %)")
    SB = pd.read_csv(f"{T}/dive_subsets.csv"); SB = SB[np.isfinite(SB.dD)]
    a = ax[0, 1]; z = SB.dD / SB.sig
    sc = a.scatter(SB.dD, range(len(SB)), c=1 / SB.sig ** 2, cmap="viridis", s=40, zorder=3)
    a.plot(SB.dD, range(len(SB)), "-", color="#bbbbbb", lw=0.6, zorder=1)
    a.set_yticks(range(len(SB))); a.set_yticklabels([f"{t} (z {x:+.1f})" for t, x in zip(SB.item, z)], fontsize=7)
    a.axvline(0, color="grey", lw=0.6, ls=":"); a.set(xlabel="ΔD(2019-11) (%)", title="subsets, controls, alternative models (σ in table)")
    fig.colorbar(sc, ax=a, pad=0.01).set_label("weight 1/σ²")
    a = ax[1, 0]
    if os.path.exists(f"{T}/dive_realign_lag_all.csv"):
        L_ = pd.read_csv(f"{T}/dive_realign_lag_all.csv")
        U = __import__("s11_common").episodes()
        for ns, g in L_.groupby("netsta"):
            g = g[g.episode.isin(U.index)].assign(mid=lambda q: q.episode.map(U.mid)).sort_values("mid")
            a.plot(g.mid, g.lag_rel, "-", lw=0.9, label=ns)
        a.axvline(U.mid.get("2019-11"), color="k", lw=0.8, ls="--")
        a.legend(fontsize=6, ncol=4); a.set(ylabel="median realignment lag − series median (ms)", title="per-event realignment lag by episode (2019-11 dashed)")
    a = ax[1, 1]
    if os.path.exists(f"{T}/dive_split.csv"):
        SP = pd.read_csv(f"{T}/dive_split.csv").set_index("part").reindex(["early", "late", "all"]).dropna(how="all")
        a.plot(range(len(SP)), SP.dvv, "-o", color="#2a78d6")
        a.set_xticks(range(len(SP))); a.set_xticklabels(["UTC ≤ 11-10 (burst phase)", "UTC ≥ 11-11", "whole episode"][:len(SP)], fontsize=8)
        a.axhline(0, color="grey", lw=0.6, ls=":")
        a.set(ylabel="dv/v minus neighbouring episodes (%)", title="2019-11 split in time (near cells, kernel-free; σ in table)")
    fig.suptitle(f"SVI 2019-11: verdict = {V['verdict']['call']}", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r4_2019_11.png", dpi=130); plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    for f in (fig_dD, fig_suite, fig_lapse, fig_dive):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            import traceback
            print("FIGURE ERROR", f.__name__, e, traceback.format_exc()[-600:])
    open(f"{B}/FIGURES_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
