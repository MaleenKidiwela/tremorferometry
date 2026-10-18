#!/usr/bin/env python
"""S11 tomo/eq figures: (1) κ per zone, earthquakes (crust / intraslab / layer) and LFEs; (2) κ of the layer
between 25 km and the LFE sources along strike, with the receiver-function range.

  python s11_eq_figures.py [--kappa eq_kappa]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = "/home/jovyan/tremorferometry"
TT = f"{ROOT}/svi_science/results/s11_depth/tomo"
EQ = f"{TT}/eq"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_eq_diff as D  # noqa: E402

COL = {"lfe": "#2a78d6", "crust": "#eb6834", "slab": "#1baf7a", "layer": "#eda100", "diff": "#4a3aa7"}
INK, MUTED = "#0b0b0b", "#52514e"
REGS = ["WA_SW", "WA_NW", "SVI"]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#e4e3df", lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8)


def pt(ax, x, k, s, color, ok, marker="o", label=None):
    ax.errorbar(x, k, yerr=s if np.isfinite(s) else None, fmt=marker, ms=6, lw=0.8, color=color,
                mfc=color if ok else "white", mec=color, capsize=0, label=label, zorder=3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", default="eq_kappa")
    ap.add_argument("--prefix", default="", help="diff table prefix and figure prefix (tests)")
    a = ap.parse_args()
    Z = pd.read_csv(f"{EQ}/tables/{a.kappa}_zones.csv")
    Z = Z[(Z.est == "clip_inv") & (Z.variant == "d100")]
    os.makedirs(f"{EQ}/figures", exist_ok=True)
    # ---- figure 1: kappa per zone
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, reg in zip(axs, REGS):
        off = {"lfe": -4, "crust": -1.3, "slab": 1.3, "layer": 4}
        for kind in ("crust", "slab", "layer"):
            q = Z[(Z.kind == kind) & (Z.reg == reg) & (Z.zone != "all") & Z.kappa.notna()]
            for _, r in q.iterrows():
                pt(ax, int(r.zone) + 15 + off[kind], r.kappa, r.sig, COL[kind], bool(r["pass"]))
            r = Z[(Z.kind == kind) & (Z.reg == reg) & (Z.zone == "all")]
            if len(r) and np.isfinite(r.iloc[0].kappa):
                ax.axhline(r.iloc[0].kappa, color=COL[kind], lw=1.0, ls="--", alpha=0.8)
        geo = D.lfe_geometry(reg)
        for z in sorted(geo.zone.unique()):
            for est, k, s, ok in D.lfe_zone(reg, int(z)):
                if np.isfinite(k):
                    pt(ax, int(z) + 15 + off["lfe"] + (0 if est == "clip" else -2), k, s, COL["lfe"], ok,
                       marker="s" if est == "clip" else "D")
        lab, k, s = D.LFE[reg][0]
        ax.axhspan(k - s, k + s, color=COL["lfe"], alpha=0.15, lw=0)
        ax.set_title(reg.replace("_", " "), fontsize=10, color=INK, loc="left")
        ax.set_xlabel("along-strike zone centre (km, az 330° from 47.0N 123.0W)", fontsize=8, color=MUTED)
        style(ax)
    axs[0].set_ylabel("Vp/Vs κ (Wadati, tP on tS)", fontsize=9, color=INK)
    axs[0].set_ylim(1.4, 2.3)
    hs = [plt.Line2D([], [], marker=m, ls="", color=c, label=l) for m, c, l in
          (("s", COL["lfe"], "LFE families (clipped; ◆ profile)"), ("o", COL["crust"], "crustal earthquakes ≤ 25 km"),
           ("o", COL["slab"], "intraslab earthquakes"), ("o", COL["layer"], "earthquakes 25–45 km"))]
    hs.append(plt.Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, label="open = fails a gate"))
    fig.legend(handles=hs, loc="upper center", ncol=5, fontsize=8, frameon=False)
    fig.text(0.01, 0.01, "Dashed lines: whole-region earthquake κ by class; blue band: whole-region LFE κ (headline). "
             "Bars: max(event, station) bootstrap σ.", fontsize=7, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 0.92))
    fig.savefig(f"{EQ}/figures/{a.prefix}fig_eq_kappa_zones.png", dpi=150)
    plt.close(fig)
    # ---- figure 2: layer kappa along strike
    L = pd.read_csv(f"{EQ}/tables/{a.prefix}diff_layer.csv")
    L = L[L.eq_est == "clip_inv"]
    fig, ax = plt.subplots(figsize=(10, 4.4))
    ax.axhspan(2.35, 2.8, color="#e4e3df", lw=0)
    ax.set_xlim(-70, 250)
    ax.text(-65, 2.75, "receiver-function κ of the thin low-velocity layer (Audet/Bostock 2.35–2.8)", fontsize=7.5,
            color=MUTED, va="center")
    xr = {}
    for reg in REGS:
        q = L[L.reg == reg]
        h = q[(q.scope == "all") & q.lfe_est.str.endswith("(headline)")]
        geo = D.lfe_geometry(reg)
        xmid = float(np.median(geo.along)); xr[reg] = xmid
        if len(h):
            r = h.iloc[0]
            ok = bool(r.crust_pass) and bool(r.lfe_pass)
            ax.errorbar(xmid, r.kappa_layer, yerr=[[r.kappa_layer - r.kl_p16], [r.kl_p84 - r.kappa_layer]], fmt="o",
                        ms=9, color=COL["diff"], mfc=COL["diff"] if ok else "white", lw=1.2, zorder=4)
            ax.annotate(f"{reg.replace('_', ' ')}\n{r.kappa_layer:.2f}", (xmid, r.kappa_layer), xytext=(8, 6),
                        textcoords="offset points", fontsize=8, color=INK)
            ax.plot([xmid - 20, xmid + 20], [r.kappa_bg_cvm] * 2, color=MUTED, lw=1.5)
        for _, r in q[(q.scope != "all") & q.lfe_est.str.contains("clip")].iterrows():
            ok = bool(r.crust_pass) and bool(r.lfe_pass)
            ax.errorbar(int(r.scope) + 15, r.kappa_layer, yerr=[[max(r.kappa_layer - r.kl_p16, 0)],
                                                                  [max(r.kl_p84 - r.kappa_layer, 0)]],
                        fmt="s", ms=5, color=COL["diff"], mfc=COL["diff"] if ok else "white", lw=0.7, alpha=0.8)
    ax.set_ylim(1.4, 3.0)
    ax.set_xlabel("along strike (km, az 330° from 47.0N 123.0W)   SW WA → NW WA → S Vancouver Island", fontsize=8,
                  color=MUTED)
    ax.set_ylabel("κ of the layer 25 km → LFE depth", fontsize=9, color=INK)
    ax.set_title("Differential Wadati: LFE κ minus crustal-earthquake κ, ray-weighted in CVM v1.7\n"
                 "circles: whole region (16–84 %); squares: 30-km zones; grey ticks: CVM v1.7 κ of the same layer; "
                 "open = an input fails a gate", fontsize=8.5, color=INK, loc="left")
    style(ax)
    fig.tight_layout()
    fig.savefig(f"{EQ}/figures/{a.prefix}fig_eq_layer_kappa.png", dpi=150)
    plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main()
