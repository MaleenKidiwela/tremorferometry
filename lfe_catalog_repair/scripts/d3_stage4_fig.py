#!/usr/bin/env python
"""D3 stage 4 figure: map + depth section + uncertainty ellipses + Lin/Bostock comparison
+ example network stacks with the automatic P/S picks (the manual-style sanity check)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.patches import Ellipse              # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from d3_stage4_stacks import STA_LL, FS, NPRE       # noqa: E402

R = Path("lfe_catalog_repair/results")
LIN_FAM = "catalogs/lin_families_ets_2010_vi.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fam", default=str(R / "d3_stage4_families.csv"))
    ap.add_argument("--picks", default=str(R / "d3_stage4_picks.csv"))
    ap.add_argument("--events", default=str(R / "d3_stage4_events.csv"))
    ap.add_argument("--stacks", default=str(R / "d3_stage4_stacks"))
    ap.add_argument("--out", default="lfe_catalog_repair/figures/d3_stage4_locations.png")
    a = ap.parse_args()
    def _read(p_):
        try:
            d_ = pd.read_csv(p_)
        except Exception:
            return pd.DataFrame()
        return d_
    F = _read(a.fam)
    ok = (F.dropna(subset=["lat"]).copy() if "lat" in F.columns else pd.DataFrame())
    lin = pd.read_csv(LIN_FAM)
    P = _read(a.picks)
    E = _read(a.events)
    Z = np.load(f"{a.stacks}.npz", allow_pickle=True) if Path(f"{a.stacks}.npz").exists() else None
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(17, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.33, wspace=0.28, height_ratios=[1.25, 1, 1])

    # ---- A map ----
    ax = fig.add_subplot(gs[0, 0:2])
    ax.scatter(lin.lon, lin.lat, s=14, marker="s", facecolor="none", edgecolor="#888",
               label=f"Lin/Bostock families (n={len(lin)})", lw=0.7)
    if len(ok):
        ax.scatter(ok.cell_lon, ok.cell_lat, s=10, c="#c8531b", marker="x",
                   label="discovery cell (template id)")
        sc = ax.scatter(ok.lon, ok.lat, s=44, c=ok.depth, cmap="viridis", vmin=25, vmax=50,
                        edgecolor="k", lw=0.4, zorder=4, label="D3 family hypocentre")
        plt.colorbar(sc, ax=ax, label="depth (km)", pad=0.01)
        for _, r in ok.iterrows():
            if np.isfinite(r.sigma_h) and r.sigma_h > 0:
                dx = r.sigma_h / (111.2 * np.cos(np.radians(r.lat)))
                dy = r.sigma_h / 111.2
                ax.add_patch(Ellipse((r.lon, r.lat), 2 * dx, 2 * dy, fill=False,
                                     ec="k", lw=0.5, alpha=0.5, zorder=3))
            ax.plot([r.cell_lon, r.lon], [r.cell_lat, r.lat], "-", color="0.6", lw=0.4, zorder=2)
    for s, (la, lo) in STA_LL.items():
        ax.plot(lo, la, "v", ms=8, color="#1b6ca8", zorder=5)
        ax.annotate(s, (lo, la), fontsize=7, xytext=(2, 4), textcoords="offset points")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("A  family hypocentres (1-sigma bootstrap ellipses) vs Lin/Bostock families "
                 "and the discovery cells", fontsize=10, loc="left")
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(alpha=0.3)

    # ---- B depth vs Slab2 ----
    ax = fig.add_subplot(gs[0, 2])
    if len(ok):
        ax.errorbar(ok.slab_depth, ok.depth, yerr=ok.sigma_z, fmt="o", ms=4, lw=0.6,
                    color="#1b7f4b", label="with Slab2 prior")
        ax.plot(ok.slab_depth, ok.depth_noprior, "x", ms=5, color="#c0392b", label="no prior")
        lim = [min(ok.slab_depth.min(), ok.depth.min()) - 4,
               max(ok.slab_depth.max(), ok.depth.max()) + 4]
        ax.plot(lim, lim, "k--", lw=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("Slab2 interface depth (km)")
        ax.set_ylabel("located depth (km)")
        ax.legend(fontsize=7)
    ax.set_title("B  depth vs the Slab2 interface", fontsize=10, loc="left")
    ax.grid(alpha=0.3)

    # ---- C depth section (E-W) ----
    ax = fig.add_subplot(gs[1, 0])
    if len(ok):
        ax.scatter(ok.lon, ok.depth, s=30, c="#1b7f4b", edgecolor="k", lw=0.3, label="located")
        ax.errorbar(ok.lon, ok.depth, yerr=ok.sigma_z, fmt="none", ecolor="0.5", lw=0.6)
        s = ok.sort_values("lon")
        ax.plot(s.lon, s.slab_depth, "-", color="#c0392b", lw=1.4, label="Slab2 at the same point")
        ax.invert_yaxis()
        ax.set_xlabel("longitude")
        ax.set_ylabel("depth (km)")
        ax.legend(fontsize=7)
    ax.set_title("C  E-W depth section", fontsize=10, loc="left")
    ax.grid(alpha=0.3)

    # ---- D offset distributions ----
    ax = fig.add_subplot(gs[1, 1])
    if len(ok):
        for col, c, lab in [("dist_to_Lin_km", "#1b6ca8", "to nearest Lin family"),
                            ("dist_to_cell_km", "#c8531b", "to the discovery cell"),
                            ("sigma_h", "0.4", "bootstrap sigma_h")]:
            v = ok[col].dropna().values
            if not len(v):
                continue
            x = np.sort(v)
            ax.plot(x, np.arange(len(x)) / len(x), color=c, lw=1.7,
                    label=f"{lab} (med {np.median(v):.1f} km)")
        ax.axvline(3, color="k", ls=":", lw=1)
        ax.text(3.1, 0.05, "3 km target", fontsize=7, rotation=90)
        ax.set_xlabel("km")
        ax.set_ylabel("cumulative fraction of families")
        ax.legend(fontsize=7)
    ax.set_title("D  offsets and uncertainties", fontsize=10, loc="left")
    ax.grid(alpha=0.3)

    # ---- E event-level scatter ----
    ax = fig.add_subplot(gs[1, 2])
    if len(E) and "de_km" in E.columns:
        g = E.groupby("family")[["de_km", "dn_km", "dz_km"]].std()
        for col, c in [("de_km", "#1b6ca8"), ("dn_km", "#1b7f4b"), ("dz_km", "#c0392b")]:
            v = g[col].dropna().values
            if not len(v):
                continue
            x = np.sort(v)
            ax.plot(x, np.arange(len(x)) / len(x), color=c, lw=1.6,
                    label=f"{col[1].upper()} (med {np.median(v):.2f} km)")
        ax.set_xlabel("within-family scatter of relative event locations (km)")
        ax.set_ylabel("cumulative fraction of families")
        ax.legend(fontsize=7)
    else:
        ax.text(.5, .5, "no event-level relocations", ha="center")
    ax.set_title("E  event-level (double-difference) scatter", fontsize=10, loc="left")
    ax.grid(alpha=0.3)

    # ---- F-H example network stacks with picks ----
    if Z is not None and len(P) and "family" in P.columns and len(ok):
        best = (ok.sort_values("n_stations", ascending=False).family.tolist()[:3]
                if len(ok) else [])
        for j, fam in enumerate(best):
            ax = fig.add_subplot(gs[2, j])
            g = P[P.family == fam].sort_values("dist_km")
            for i, (_, r) in enumerate(g.iterrows()):
                kz = f"{fam}|{r.station}|Z"
                k1, k2 = f"{fam}|{r.station}|H1", f"{fam}|{r.station}|H2"
                if kz not in Z.files:
                    continue
                t = (np.arange(len(Z[kz])) - NPRE) / FS
                z = np.asarray(Z[kz], float)
                h = (np.sqrt(np.asarray(Z[k1], float) ** 2 + np.asarray(Z[k2], float) ** 2)
                     if k1 in Z.files else np.abs(z))
                sc = 0.42 / (np.abs(z).max() + 1e-12)
                sh = 0.42 / (np.abs(h).max() + 1e-12)
                ax.plot(t, i + z * sc, color="0.35", lw=0.6)
                ax.plot(t, i + h * sh, color="#1b6ca8", lw=0.6, alpha=0.8)
                if np.isfinite(r.i_P):
                    ax.plot([(r.i_P - NPRE) / FS] * 2, [i - .45, i + .45], "-",
                            color="#c0392b", lw=1.3)
                if np.isfinite(r.i_S):
                    ax.plot([(r.i_S - NPRE) / FS] * 2, [i - .45, i + .45], "-",
                            color="#1b7f4b", lw=1.3)
                ax.annotate(f"{r.station} {r.dist_km:.0f}km", (t[0] + 0.3, i + 0.18), fontsize=6)
            ax.axvline(0, color="k", ls=":", lw=0.8)
            ax.set_xlabel("s relative to the predicted S")
            ax.set_yticks([])
            ax.set_title(f"{'FGH'[j]}  {fam}  (grey Z, blue |H|; red P, green S)",
                         fontsize=8.5, loc="left")
    fig.suptitle("D3 stage 4 — LFE family locations from network-consistent template stacks "
                 "(Crosson 1976 1-D model, Vp/Vs 1.78, Slab2 weak prior)", fontsize=12)
    fig.savefig(a.out, dpi=125, bbox_inches="tight")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
