#!/usr/bin/env python
"""Figures for P1b and P2.

  python scripts/p_figs.py p1b --tags pnet_w05,lin_w05,pnet_w10 --year 2012
  python scripts/p_figs.py p2map --tag p2core_2011 --doys 240,243
  python scripts/p_figs.py p2diag --tag p2core_2011
"""
from __future__ import annotations
import os, sys, json, argparse, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

FIG = f"{P.PICK}/figures"
NEAR = ["PGC", "B009", "B010", "B011", "B926", "NLLB", "SNB", "VGZ"]
FAR = ["B001", "B003", "B004", "B005", "B006", "B007", "LRIV"]


def p1b(a):
    tags = a.tags.split(",")
    fig, axs = plt.subplots(1, len(tags), figsize=(4.6 * len(tags), 4.4), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, t in zip(axs, tags):
        f = f"{P.RES}/p1/curve_{t}_{a.year}_tol1.csv"
        if not os.path.exists(f):
            ax.set_title(f"{t}\n(missing)")
            continue
        c = pd.read_csv(f)
        for s in sorted(c.sta.unique()):
            g = c[c.sta == s].sort_values("rate")
            col = "tab:blue" if s in NEAR else "tab:red"
            ls = "--" if s in P.HELD_OUT_STATIONS else "-"
            ax.plot(g.rate, g.recall, ls, color=col, lw=1.2, alpha=.8,
                    label=s if s in ("PGC", "B003") else None)
        ax.axvline(100, color="k", lw=.8, ls=":")
        ax.axhline(0.30, color="g", lw=.8, ls=":")
        ax.axhline(0.119, color="purple", lw=.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("declustered peaks / station-day")
        ax.set_title(t)
        ax.grid(alpha=.3)
    axs[0].set_ylabel(r"single-station S recall, $\pm$1 s")
    axs[0].legend(fontsize=8, title="blue = near sites\nred = far boreholes\ndashed = held out")
    fig.suptitle(f"P1b single-station S recall vs peaks/station-day ({a.year} development)\n"
                 "green = Merlin's 0.30 bar   purple = template cc at 100 peaks/TEMPLATE-day (0.119)"
                 "   dotted vertical = 100 peaks/station-day", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{FIG}/p1b_recall_curves_{a.year}.png", dpi=140)
    print(f"-> {FIG}/p1b_recall_curves_{a.year}.png")


def p2map(a):
    for d in a.doys.split(","):
        f = f"{P.RES}/p2/maps_{a.tag}_{int(d)}.npz"
        if not os.path.exists(f):
            print("missing", f)
            continue
        z = np.load(f)
        Z, sh = z["z"].astype(np.float32), tuple(z["shape"])
        lat, lon = z["lat"], z["lon"]
        det = pd.read_parquet(f"{P.RES}/p2/det_{a.tag}.parquet")
        det = det[det.doy == int(d)]
        peak = Z.max(1).reshape(sh)
        occ = (Z > 6).sum(1).reshape(sh).astype(float)
        FAMs = P.families()
        fig, axs = plt.subplots(1, 2, figsize=(11, 5))
        for ax, M, ttl in [(axs[0], peak, "peak Z over the day"),
                           (axs[1], np.log10(occ + 1), r"log$_{10}$(1 + # 0.25 s slots with Z > 6)")]:
            im = ax.pcolormesh(lon.reshape(sh), lat.reshape(sh), M, shading="auto", cmap="magma")
            ax.plot(FAMs.lon, FAMs.lat, "c+", ms=4, mew=.7, label="Bostock families")
            ST = P.station_table()
            ax.plot(ST.lon, ST.lat, "w^", ms=5, mec="k", label="stations")
            plt.colorbar(im, ax=ax)
            ax.set_title(ttl)
            ax.set_xlabel("lon")
        axs[0].set_ylabel("lat")
        axs[0].legend(fontsize=7, loc="lower left")
        fig.suptitle(f"P2 back-projection, {a.tag} doy {d} ({len(det)} stored peaks)")
        fig.tight_layout()
        fig.savefig(f"{FIG}/p_p2_map_{a.tag}_{int(d)}.png", dpi=140)
        print(f"-> {FIG}/p_p2_map_{a.tag}_{int(d)}.png")


def p2diag(a):
    det = pd.read_parquet(f"{P.RES}/p2/det_scored_{a.tag}.parquet")
    FAMs = P.families()
    fa = FAMs[["lat", "lon"]].to_numpy()
    dd = np.stack([P.epi_km(det.lat.values, det.lon.values, fa[i, 0], fa[i, 1])
                   for i in range(len(fa))], 1).min(1)
    cal = pd.read_csv(f"{P.RES}/p2/calib_{a.tag}.csv")
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    axs[0].hist(dd, bins=np.arange(0, 40, 1), color="tab:blue")
    axs[0].axvline(5, color="r", ls="--")
    axs[0].set_xlabel("km from the nearest Bostock family")
    axs[0].set_ylabel("detections")
    axs[0].set_title(f"centroid collapse: {np.mean(dd>5):.1%} > 5 km from every known family")
    ST = P.station_table()
    axs[1].hexbin(det.lon, det.lat, gridsize=45, cmap="magma", bins="log")
    axs[1].plot(FAMs.lon, FAMs.lat, "c+", ms=4, mew=.7)
    axs[1].plot(ST.lon, ST.lat, "w^", ms=5, mec="k")
    axs[1].set_title("where detections back-project")
    axs[1].set_xlabel("lon"); axs[1].set_ylabel("lat")
    axs[2].plot([0, 1], [0, 1], "k--", lw=.8)
    axs[2].plot(cal.p_mean, cal.bostock, "o-", label="Bostock match rate")
    if "withheld" in cal and np.isfinite(cal.withheld).any():
        axs[2].plot(cal.p_mean, cal.withheld, "s-", label="withheld-station excess rate")
    axs[2].set_xlabel("p_real (isotonic, fitted on the dev year)")
    axs[2].set_ylabel("observed")
    axs[2].legend(fontsize=8)
    axs[2].set_title("reliability")
    fig.suptitle(f"P2 diagnostics — {a.tag}")
    fig.tight_layout()
    fig.savefig(f"{FIG}/p_p2_diag_{a.tag}.png", dpi=140)
    print(f"-> {FIG}/p_p2_diag_{a.tag}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["p1b", "p2map", "p2diag"])
    ap.add_argument("--tags", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--doys", default="")
    a = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    {"p1b": p1b, "p2map": p2map, "p2diag": p2diag}[a.what](a)
