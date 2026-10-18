#!/usr/bin/env python
"""S11 tomo figures (tomo/figures/fig_tomo_*.png). Colour = quality (bootstrap spread / rms); no error bars."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo"
CAT = {"BV3": "#2a78d6", "CV3": "#eb6834", "S1V3": "#1baf7a", "B": "#4a3aa7", "A": "#e87ba4", "AV3": "#e87ba4",
       "C": "#eda100", "EallV3": "#e34948", "E60V3": "#008300"}
NAME = {"BV3": "WA NW (V3 P)", "CV3": "WA SW (V3 P)", "S1V3": "SVI 2012 (V3 P)", "B": "WA NW (frozen P)",
        "C": "WA SW (frozen P)", "EallV3": "SVI 2003-06 POLARIS era, all stations",
        "E60V3": "SVI 2003-06 era, <= 60 km", "A": "WA NW run-1 eye-real (frozen P)", "AV3": "WA NW run-1 eye-real (V3 P)"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6})


def fig_family():
    sets = [("A", "spec set A, OLS"), ("B", "catalogue NW, OLS"), ("C", "catalogue SW, OLS"),
            ("BV3_V1inv", "NW, V3 P, clipped, tP-on-tS"), ("CV3_V1inv", "SW, V3 P, clipped, tP-on-tS"),
            ("S1V3_V1inv", "SVI 2012, V3 P, clipped, tP-on-tS")]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True)
    for ax, (s, lab) in zip(axs.flat, sets):
        K = pd.read_csv(f"{T}/tables/kappa_{s}.csv")
        K = K[K.measured]
        ax.axvspan(0.0, 0.03, color="#cde2fb", alpha=0.6, lw=0)
        ax.axhspan(1.6, 2.2, color="#f0efec", lw=0)
        ax.scatter(K.scr1_sig.clip(0.002, 5), K.scr1_kappa.clip(-1, 4), s=10, c="#a8a7a0", marker="x", lw=0.8,
                   label="scrambled control (1 draw)")
        sc = ax.scatter(K.sig_k.clip(0.002, 5), K.kappa.clip(-1, 4), s=22, c=K.rms.clip(0, 3), cmap="Blues_r",
                        vmin=0, vmax=3, edgecolor="#52514e", lw=0.4, label="families")
        ax.set_xscale("log"); ax.set_xlim(0.005, 5); ax.set_ylim(-1, 4)
        ax.set_title(f"{lab}\n{len(K)} measured, {int(((K.sig_k <= 0.03) & (K.rms <= 0.25)).sum())} pass G1",
                     fontsize=8.5)
    for ax in axs[1]:
        ax.set_xlabel("bootstrap σκ (log)")
    for ax in axs[:, 0]:
        ax.set_ylabel("per-family κ")
    axs[0, 0].legend(loc="lower left", fontsize=7, frameon=False)
    fig.colorbar(sc, ax=axs, label="fit rms (s)", shrink=0.6)
    fig.suptitle("Per-family Wadati κ: none reaches the G1 box (blue band σκ ≤ 0.03 s; grey band = plausible κ)",
                 fontsize=10)
    fig.savefig(f"{T}/figures/fig_tomo_family.png", dpi=110, bbox_inches="tight")


def fig_profiles():
    z = np.load(f"{T}/tables/profile_V3wa_curves.npz")
    z2 = np.load(f"{T}/tables/profile_other_curves.npz")
    zc = np.load(f"{T}/tables/profile_CV3_curves.npz")
    ze = np.load(f"{T}/tables/profile_E_curves.npz")
    kg = z["kgrid"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for key, src in (("BV3|all", z), ("CV3|all", zc), ("B|all", z2), ("AV3|eye_lfe", z2), ("S1V3|all", z2), ("EallV3|all", ze),
                     ("E60V3|all", ze)):
        s = src[key]; sc = src[key + "|scr"]
        st = key.split("|")[0]
        ax.plot(kg, s / np.percentile(sc, 99), color=CAT[st], lw=2, label=NAME[st])
    ax.axhline(1, color="#52514e", lw=1, ls="--")
    ax.text(1.31, 1.03, "scrambled-control 99th percentile of the peak", fontsize=8, color="#52514e")
    ax.set_xlabel("κ = Vp/Vs"); ax.set_ylabel("inlier count / control p99")
    ax.set_title("Inlier profile per region (|residual| ≤ 0.25 s; families ≥ 4 stations)")
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(f"{T}/figures/fig_tomo_profiles.png", dpi=110, bbox_inches="tight")


def fig_zones():
    Z = pd.concat([pd.read_csv(f"{T}/tables/profile_V3wa_zones.csv"),
                   pd.read_csv(f"{T}/tables/profile_other_zones.csv")])
    Z = Z[~((Z.set == "CV3"))]
    Z = pd.concat([Z, pd.read_csv(f"{T}/tables/profile_CV3_zones.csv"), pd.read_csv(f"{T}/tables/profile_E_zones.csv")])
    Z = Z[Z.kappa.notna()]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, kind, xc, xl in ((axs[0], "along", "along_med", "along-strike position (km, azimuth 330° from 47N 123W)"),
                             (axs[1], "depth", "depth_med", "median catalogue depth of the zone's families (km)")):
        q = Z[Z.zone.str.startswith(kind)]
        for st, g in q.groupby("set"):
            if st not in ("BV3", "CV3", "S1V3", "B", "EallV3"):
                continue
            g = g.sort_values(xc)
            ax.plot(g[xc], g.kappa, color=CAT[st], lw=1.2, alpha=0.7)
            sc = ax.scatter(g[xc], g.kappa, s=30 + g.n_fam, c=g.sig.clip(0, 0.3), cmap="Blues_r", vmin=0, vmax=0.3,
                            edgecolor=CAT[st], lw=1.6, zorder=3)
            ax.annotate(NAME[st], (g[xc].iloc[-1], g.kappa.iloc[-1]), fontsize=7, color="#52514e",
                        xytext=(4, 0), textcoords="offset points")
        allz = Z[Z.zone == "all"]
        ax.set_xlabel(xl)
    axs[0].set_ylabel("zone κ (inlier-profile peak)")
    for st in ("BV3", "CV3", "S1V3", "EallV3"):
        r = Z[(Z.set == st) & (Z.zone == "all")]
        if len(r):
            axs[0].axhline(r.kappa.iloc[0], color=CAT[st], lw=0.8, ls=":")
    fig.colorbar(sc, ax=axs, label="bootstrap σκ (darker = tighter)", shrink=0.8)
    fig.suptitle("Zone κ along strike and with depth (marker edge = region; dotted = whole-region value)", fontsize=10)
    fig.savefig(f"{T}/figures/fig_tomo_zones.png", dpi=110, bbox_inches="tight")


def fig_distance():
    sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
    cwd = os.getcwd(); os.chdir(f"{ROOT}/lfe_array_discovery")
    import a4f_common as AC
    os.chdir(cwd)
    I = pd.read_csv(f"{T}/tables/profile_V3wa_inliers.csv", low_memory=False)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, st in zip(axs, ("BV3", "CV3")):
        P = pd.read_csv(f"{T}/tables/picks_{st}.csv", low_memory=False)
        P = P[(P.rP >= 1.2) & (P.rS >= 1.2) & P.tS.notna()].copy()
        for D in (P,):
            D["slat"] = D.sta.map(lambda s: AC.COORD[s][0]); D["slon"] = D.sta.map(lambda s: AC.COORD[s][1])
            D["hyp"] = np.hypot(np.hypot((D.slon - D.lon) * 111.19 * np.cos(np.radians(47.3)),
                                         (D.slat - D.lat) * 111.19), D.depth_km.clip(20, 50))
        key = set(zip(I[I.set == st].fam.astype(str), I[I.set == st].sta))
        P["inl"] = [(str(f), s) in key for f, s in zip(P.fam, P.sta)]
        ax.scatter(P.hyp[~P.inl], (P.tS - P.tP)[~P.inl], s=7, c="#c3c2b7", label="outside the 0.25-s band")
        ax.scatter(P.hyp[P.inl], (P.tS - P.tP)[P.inl], s=9, c=CAT[st], label="inlier at the peak κ")
        h = np.linspace(20, 120, 20)
        for k, ls in ((1.73, "--"), (1.9, "-")):
            ax.plot(h, h * (k - 1) / 6.8, color="#0b0b0b", ls=ls, lw=1, label=f"κ {k}, Vp 6.8 km/s (guide only)")
        ax.set_title(NAME[st]); ax.set_xlabel("hypocentral distance from catalogue location (km)")
        ax.set_xlim(15, 130)
    axs[0].set_ylabel("tS − tP (s)"); axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("S−P vs distance (distance is only for display; κ is estimated without locations)", fontsize=10)
    fig.savefig(f"{T}/figures/fig_tomo_sp_distance.png", dpi=110, bbox_inches="tight")


if __name__ == "__main__":
    for f in (fig_family, fig_profiles, fig_zones, fig_distance):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            print("FAILED", f.__name__, repr(e))
