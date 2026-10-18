#!/usr/bin/env python
"""s11_dist_figures -- figures of the S11 distant-station lever (distant/PLAN.md). User plot rules: dv/v points coloured by
CC (or weight), NO error bars on dv/v, line graphs. Reads distant/<region>/tables; writes distant/figures/*.png.
    python s11_dist_figures.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/distant"
FIG = f"{B}/figures"
NEAR, FAR, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
REG = {"svi": "Southern Vancouver Island", "wa_nw": "Washington NW"}
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "font.size": 9, "axes.grid": True, "grid.color": "#e6e5e1",
                     "grid.linewidth": 0.6})


def fig_inventory():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, r in zip(axs, REG):
        p = f"{B}/{r}/tables/far_pairs.csv"
        if not os.path.exists(p):
            ax.set_visible(False); continue
        F = pd.read_csv(p); thr = json.load(open(f"{B}/{r}/tables/inventory_summary.json"))["ctrl_snr_p99"]
        sc = ax.scatter(F.hypo_km, F.pick_snr, c=F.causal3c.clip(0.8, 3), cmap="viridis", s=10, linewidths=0)
        ax.scatter(F.hypo_km[F.usable], F.pick_snr[F.usable], s=26, facecolors="none", edgecolors=FAR, linewidths=1.0, label="usable pair")
        ax.axhline(thr, color=MUTED, lw=1, ls="--"); ax.text(148, thr * 0.82, "scrambled-control 99th pct", ha="right", color=MUTED, fontsize=8)
        ax.set_yscale("log"); ax.set_xlabel("hypocentral distance (km)"); ax.set_title(f"{REG[r]}: {int(F.usable.sum())} usable of {len(F)} far pairs", fontsize=10)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
    axs[0].set_ylabel("S pick SNR of the sample stack")
    cb = fig.colorbar(sc, ax=axs, shrink=0.8); cb.set_label("causality (coda 2-5 s / pre-P RMS)")
    fig.savefig(f"{FIG}/fig_dist_inventory.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_annual():
    """kernel-free annual cycle: series-demeaned dv/v vs day of year, far vs near, cells coloured by cc; weighted monthly lines"""
    regs = [r for r in REG if os.path.exists(f"{B}/{r}/tables/inv_cells.csv")]
    if not regs:
        return
    fig, axs = plt.subplots(len(regs), 2, figsize=(11, 3.6 * len(regs)), squeeze=False, sharey=True)
    for i, r in enumerate(regs):
        Cc = pd.read_csv(f"{B}/{r}/tables/inv_cells.csv")
        w = 1 / Cc.sig ** 2; key = Cc.fam.astype(str) + "|" + Cc.netsta + "|" + Cc.epoch
        Cc["dm"] = Cc.d - (Cc.d * w).groupby(key).transform("sum") / w.groupby(key).transform("sum")
        Cc["doy"] = pd.to_datetime(Cc["mid"]).dt.dayofyear
        for j, (far, col, lab) in enumerate(((False, NEAR, "near (<= 80 km)"), (True, FAR, "far (80-150 km)"))):
            ax = axs[i, j]; q = Cc[(Cc.far == far) & Cc.doy.notna()]
            ax.scatter(q.doy, q.dm, c=q.cc, cmap="viridis", vmin=0.35, vmax=1, s=5, linewidths=0, alpha=0.6)
            mb = (q.doy // 30.5).clip(0, 11); ww = 1 / q.sig ** 2
            mm = (q.dm * ww).groupby(mb).sum() / ww.groupby(mb).sum()
            ax.plot(mm.index * 30.5 + 15, mm.values, color=col, lw=2, marker="o", ms=4, label=f"{lab}: weighted monthly mean")
            ax.axhline(0, color=MUTED, lw=0.8); ax.set_ylim(-0.6, 0.6); ax.set_xlim(0, 366)
            ax.set_title(f"{REG[r]} — {lab}, {len(q)} cells", fontsize=10); ax.legend(frameon=False, fontsize=8, loc="lower left")
            ax.set_xlabel("day of year (episode)")
        axs[i, 0].set_ylabel("dv/v, series-demeaned (%)")
    fig.savefig(f"{FIG}/fig_dist_annual.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_q3():
    regs = [r for r in REG if os.path.exists(f"{B}/{r}/tables/inv_q3.csv")]
    if not regs:
        return
    fig, axs = plt.subplots(len(regs), 1, figsize=(11, 3.4 * len(regs)), squeeze=False)
    for i, r in enumerate(regs):
        T = pd.read_csv(f"{B}/{r}/tables/inv_q3.csv"); T = T[T.major15]; ax = axs[i, 0]
        for (cfg, setup), col, lab in ((("N1", "NOM"), NEAR, "near only (N1)"), (("F1", "NOM"), FAR, "near + far (F1)"),
                                       (("F1", "NOM_Sc"), "#1baf7a", "near + far, regional shallow free (F1, T+Sc)")):
            q = T[(T.config == cfg) & (T.setup == setup)].sort_values("mid")
            if not len(q):
                continue
            x = pd.to_datetime(q["mid"])
            ax.plot(x, q.dD, color=col, lw=1.5, alpha=0.9, label=lab)
            ax.scatter(x, q.dD, c=1 / q.sig ** 2, cmap="viridis", s=22, zorder=3, linewidths=0)
        ax.axhline(0, color=MUTED, lw=0.8); ax.set_ylabel("ΔD per major (%)")
        ax.set_title(f"{REG[r]}: episode-locked deep-weighted change vs neighbours (points coloured by weight 1/σ²)", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=3)
    fig.savefig(f"{FIG}/fig_dist_q3.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_lever():
    regs = [r for r in REG if os.path.exists(f"{B}/{r}/tables/inv_summary.json")]
    if not regs:
        return
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
    mets = (("shallow-pulse leak (median |·|)", lambda s, c: s[c]["shallow_pulse_median_abs"]),
            ("σ ΔD per major, median (%)", lambda s, c: s[c]["q3_sig_major_median"]),
            ("B shallow-annual → Q1 leak (two5_40 kernel)", lambda s, c: s[c]["B_Q1_leak"]["two5_40"]))
    cfgs = ["N0", "N1", "F1"]
    for ax, (t, f) in zip(axs, mets):
        for r, col in zip(regs, ("#1baf7a", "#4a3aa7")):
            s = json.load(open(f"{B}/{r}/tables/inv_summary.json"))
            ys = [f(s, c) if c in s else np.nan for c in cfgs]
            ax.plot(cfgs, ys, color=col, lw=2, marker="o", ms=6, label=REG[r])
            sc = s.get("F1", {}).get("Sc_model") if t.startswith("shallow") or t.startswith("σ") else None
            if sc:
                v = sc["shallow_pulse_median_abs"] if t.startswith("shallow") else sc["q3_sig_major_median"]
                ax.plot(["F1"], [v], marker="D", ms=7, color=col, markerfacecolor="white", linestyle="none")
        ax.set_title(t, fontsize=9); ax.set_ylim(bottom=0)
    axs[0].legend(frameon=False, fontsize=8); fig.text(0.5, -0.04, "N0 near, existing model; N1 near + mid-crust/below-45 km terms; F1 near + far (headline); "
                                                         "open diamond = F1 with a free regional shallow term per episode", ha="center", color=MUTED, fontsize=8)
    fig.savefig(f"{FIG}/fig_dist_lever.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_b2():
    """per-major dD with the round-3 B2 graph estimator for near cells: near only, near + far, near + far with Sc"""
    regs = [r for r in REG if os.path.exists(f"{B}/{r}/tables/b2_q3.csv")]
    if not regs:
        return
    fig, axs = plt.subplots(len(regs), 1, figsize=(11, 3.4 * len(regs)), squeeze=False)
    for i, r in enumerate(regs):
        T = pd.read_csv(f"{B}/{r}/tables/b2_q3.csv"); T = T[T.major15]; ax = axs[i, 0]
        for (cfg, sc), col, lab in ((("N1_B2", False), NEAR, "near only, graph estimator"), (("F1_B2", False), FAR, "near (graph) + far"),
                                    (("F1_B2", True), "#1baf7a", "near (graph) + far, regional shallow free"),
                                    (("FO_B0", False), "#4a3aa7", "far only")):
            q = T[(T.config == cfg) & (T.Sc == sc)].sort_values("mid")
            if not len(q):
                continue
            x = pd.to_datetime(q["mid"])
            ax.plot(x, q.dD, color=col, lw=1.4, alpha=0.9, label=lab)
            ax.scatter(x, q.dD, c=1 / q.sig ** 2, cmap="viridis", s=20, zorder=3, linewidths=0)
        ax.axhline(0, color=MUTED, lw=0.8); ax.set_ylabel("ΔD per major (%)")
        ax.set_title(f"{REG[r]}: per-major deep-weighted change, graph estimator for near cells (points coloured by weight 1/σ²)", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=4)
    fig.savefig(f"{FIG}/fig_dist_b2.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_station():
    p = f"{B}/station_near_far_annual.csv"
    if not os.path.exists(p):
        return
    A = pd.read_csv(p); A = A[A.shared]
    fig, ax = plt.subplots(figsize=(11, 3.8))
    st = sorted(A.netsta.unique()); xs = {s: i for i, s in enumerate(st)}
    for far, col, lab, dx in ((False, NEAR, "near sources (≤ 80 km)", -0.12), (True, FAR, "far sources (80–150 km)", 0.12)):
        q = A[A.far == far]
        ax.scatter([xs[s] + dx for s in q.netsta], q.amp, s=np.clip(4 * np.sqrt(q.n), 10, 140), color=col, label=lab, zorder=3)
    ax.set_xticks(range(len(st))); ax.set_xticklabels([s.split(".")[1] for s in st])
    ax.set_ylabel("annual amplitude, kernel-free (%)"); ax.set_ylim(bottom=0)
    ax.set_title("Same station, near vs far families: free-phase annual amplitude (marker size ∝ √cells; σ in station_near_far_annual.csv)", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(f"{FIG}/fig_dist_station.png", dpi=150, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    for f in (fig_inventory, fig_annual, fig_q3, fig_lever, fig_b2, fig_station):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            import traceback; print("FAILED", f.__name__, e, traceback.format_exc()[-500:])
    open(f"{B}/FIGURES_DONE", "w").write("done\n")
