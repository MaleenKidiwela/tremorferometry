#!/usr/bin/env python
"""P3D DISCOVERY figures: figures/s_discovery_*.png"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

FIG, OUT = S.FIG, S.P3D
C_OLD, C_NEW, C_PQ, C_HID = "#4c6ef5", "#f76707", "#12b886", "#e03131"
ARMS = [("B_true", "ARM A  pairwise"), ("Bseed_true", "ARM B  seed-and-grow"),
        ("C_true", "ARM C  CF-triggered")]


def _load(tag):
    p = f"{OUT}/disc_{tag}_clusters.csv"
    if not os.path.exists(p):
        return None, None, None
    try:
        C = pd.read_csv(p)
    except Exception:
        return None, None, None
    M = (pd.read_parquet(f"{OUT}/disc_{tag}_members.parquet")
         if os.path.exists(f"{OUT}/disc_{tag}_members.parquet") else pd.DataFrame())
    try:
        Sc = pd.read_csv(f"{OUT}/disc_{tag}_scored.csv")
    except Exception:
        Sc = pd.DataFrame()
    return C, M, Sc


def fig_spacetime():
    """Space-time: cluster members in time, against the hidden families' T1 events."""
    T1 = S.hp_events(0.8, episodes=["S18", "S21", "S22"])
    hid, iso, _ = S.holdouts()
    rows = [(t, C, M, Sc) for t, lab in ARMS for C, M, Sc in [_load(t)] if C is not None]
    if not rows:
        print("no arm products yet")
        return
    fig, ax = plt.subplots(len(rows), 1, figsize=(13, 3.1 * len(rows)), sharex=True, squeeze=False)
    for k, (tag, C, M, Sc) in enumerate(rows):
        a = ax[k][0]
        lab = dict(ARMS)[tag]
        th = pd.to_datetime(T1[T1.bkey.isin(iso)].t, unit="s")
        a.plot(th, np.zeros(len(th)) - 0.6, "|", color=C_HID, ms=6, alpha=.5,
               label="T1: isolated hidden families")
        if len(M):
            keep = set(C[C.pass_excess].cid) if "pass_excess" in C else set()
            mm = M[M.cid.isin(keep)]
            for i, (cid, g) in enumerate(mm.groupby("cid")):
                a.plot(pd.to_datetime(g.t, unit="s"), np.full(len(g), i % 25 + 1), ".",
                       ms=2.5, color=C_NEW, alpha=.7)
            a.set_ylabel("cluster (stacked)")
        a.set_title(f"{lab} — {0 if C is None else int(C.pass_excess.sum()) if 'pass_excess' in C else 0} "
                    f"clusters past the S-excess clause", loc="left", fontsize=10)
        a.grid(alpha=.25)
        a.legend(fontsize=7, loc="upper right")
    ax[-1][0].set_xlabel("time (pooled S18 + S21 + S22)")
    fig.tight_layout()
    fig.savefig(f"{FIG}/s_discovery_spacetime.png", dpi=150)
    print(f"-> {FIG}/s_discovery_spacetime.png")


def fig_map():
    """Map: Bostock families (kept vs hidden), the network, and any located cluster."""
    F = S.families()
    hid, iso, _ = S.holdouts()
    ST = pd.read_csv(f"{OUT}/stations_p3d_v2.csv").set_index("sta")
    SI = S.site_table(ST)
    fig, ax = plt.subplots(figsize=(8.2, 7.4))
    k = F[~F.index.isin(hid)]
    ax.scatter(k.lon, k.lat, s=16, c="0.72", label="kept Bostock families", zorder=2)
    h = F[F.index.isin(hid) & ~F.index.isin(iso)]
    ax.scatter(h.lon, h.lat, s=34, facecolor="none", edgecolor=C_HID, lw=1.1,
               label="hidden (other strata)", zorder=3)
    i = F[F.index.isin(iso)]
    ax.scatter(i.lon, i.lat, s=64, c=C_HID, marker="o", ec="k", lw=.5,
               label="hidden, ISOLATED stratum", zorder=4)
    for b in i.index:
        ax.annotate(str(b), (i.lon[b], i.lat[b]), fontsize=7, xytext=(4, 3),
                    textcoords="offset points", color=C_HID)
    for r in SI.itertuples():
        pq = r.site in ("LHLYB", "JRDRB", "ALBHB", "DAOB")
        c = C_PQ if pq else (C_NEW if r.role == "new" else C_OLD)
        m = "^" if not r.withheld else "v"
        ax.plot(r.lon, r.lat, m, ms=11, color=c, mec="k", mew=.6, zorder=5)
        ax.annotate(r.site, (r.lon, r.lat), fontsize=7, xytext=(3, -9),
                    textcoords="offset points")
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ax.set_title("P3D discovery: the V2 network and the truth\n"
                 "blue = frozen sites, orange = P3D new, green = PQ; v = withheld (never summed)")
    ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(f"{FIG}/s_discovery_map.png", dpi=150)
    print(f"-> {FIG}/s_discovery_map.png")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_map()
    fig_spacetime()
