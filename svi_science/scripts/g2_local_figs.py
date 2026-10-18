#!/usr/bin/env python
"""g2_local_figs (2026-10-11) -- Study G2 summary figures: catalogue map + classes, Q1 superposed epoch (hours and
days), Q2 rate ratios with control ranges.  Light, single process.

    python svi_science/scripts/g2_local_figs.py  -> results/g2_local/figures/g2_{map,superposed,q2_rates}.png
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import g2_local_tests as LT

ROOT = "/home/jovyan/tremorferometry"
G2 = f"{ROOT}/svi_science/results/g2_local"
T = f"{G2}/tables"; FIG = f"{G2}/figures"
NS = 10 ** 9; H = 3600 * NS
COL = {"crustal": "#4a7fb5", "intraslab": "#c0504d", "near_interface": "#e5a33d", "overriding_deep": "#8a8a8a",
       "offshore_noslab": "#bbbbbb", "unknown_depth": "#dddddd"}


def fig_map():
    A = pd.read_parquet(f"{T}/eq_local.parquet")
    A = A[A.in_analysis & A.gk_main & (A.mag >= 1.5)]
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    W = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv"); W = W[W.ab & W.copy_of.isna()]
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw=dict(width_ratios=[1.25, 1]))
    for c in ("crustal", "overriding_deep", "near_interface", "intraslab"):
        q = A[A.cls == c]
        ax[0].scatter(q.lon, q.lat, s=2 + 2 * (q.mag - 1.5) ** 2, c=COL[c], alpha=0.5, lw=0, label=f"{c} ({len(q)})")
    ax[0].scatter(F.lon, F.lat, s=10, marker="^", c="k", label="SVI families")
    ax[0].scatter(W.lon, W.lat, s=10, marker="v", c="g", label="WA-NW families")
    ax[0].set_xlim(-126.8, -120.4); ax[0].set_ylim(45.9, 50.5); ax[0].set_aspect(1 / np.cos(np.radians(48)))
    ax[0].legend(fontsize=7, loc="lower left", markerscale=3); ax[0].set_title("GK mainshocks M ≥ 1.5 within 150 km of the families")
    q = A[A.WA_dist_km.le(150) | A.SVI_dist_km.le(150)]
    for c in ("crustal", "overriding_deep", "near_interface", "intraslab"):
        z = q[q.cls == c]
        ax[1].scatter(z.SVI_downdip_km, z.depth, s=3, c=COL[c], alpha=0.4, lw=0)
    ax[1].scatter(F.downdip_km, F.dep, s=12, marker="^", c="k")
    ax[1].invert_yaxis(); ax[1].set_xlabel("down-dip distance from SVI family centre (km, Slab2 az)")
    ax[1].set_ylabel("depth (km)"); ax[1].set_title("SVI profile: classes vs depth"); ax[1].set_xlim(-150, 150)
    fig.tight_layout(); fig.savefig(f"{FIG}/g2_map.png", dpi=110); plt.close(fig)


def fig_superposed():
    fig, ax = plt.subplots(2, 2, figsize=(12, 7))
    for j, reg in enumerate(("SVI", "WA")):
        TT = pd.read_parquet(f"{T}/q1_trials_{reg}.parquet")
        K, W, days, dP, ep = LT.load_lfe(reg)
        for i, (lab, sel) in enumerate((("all classes", TT.status.eq("between")),
                                         ("intraslab + near-interface", TT.status.eq("between") & TT.cls.isin(["intraslab", "near_interface"])))):
            t0 = TT.t0[sel].to_numpy()
            edges = np.arange(-24, 25, 1)
            yd = np.array([LT.wcount(W["data"], t0 + a * H, t0 + (a + 1) * H).sum() for a in edges[:-1]])
            yc = np.array([LT.wcount(W["c1"], t0 + a * H, t0 + (a + 1) * H).sum() for a in edges[:-1]])
            ax[i, j].bar(edges[:-1] + 0.5, yd - yc, width=0.9, color="C0", label="Σ P data − scrambled")
            ax[i, j].axhline((yd - yc)[:24].mean(), color="k", ls="--", lw=0.8, label="mean of the 24 h before")
            ax[i, j].axvline(0, color="C1")
            ax[i, j].set_title(f"{reg}: {lab}, between ETS (n = {sel.sum()})", fontsize=9)
            ax[i, j].set_xlabel("hours after local earthquake origin"); ax[i, j].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/g2_superposed.png", dpi=110); plt.close(fig)


def fig_q2():
    Q = pd.read_csv(f"{T}/q2_rates.csv")
    tests = ["ets_all", "ets_major", "ets_all_after10", "tremor_days", "lfe_active_days", "bursts"]
    sets = ["all", "crustal", "intraslab", "near_interface", "overriding_deep", "beneath_intraslab"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for j, reg in enumerate(("SVI", "WA")):
        for i, t in enumerate(tests):
            for k, s in enumerate(sets):
                r = Q[(Q.region == reg) & (Q.test == t) & (Q.set == s)]
                if not len(r):
                    continue
                r = r.iloc[0]; x = i + (k - 2.5) * 0.13
                ax[j].plot([x, x], [r.get("RR_lb99_exact", np.nan), r.get("RR_ub99_exact", np.nan)], color=COL.get(s, "k"), lw=1)
                ax[j].scatter(x, r.RR, s=18, color=COL.get(s, "k"), zorder=3, label=s if i == 0 else None)
        ax[j].axhline(1, color="0.5", ls="--"); ax[j].set_yscale("log"); ax[j].set_ylim(0.1, 10)
        ax[j].set_xticks(range(len(tests))); ax[j].set_xticklabels(tests, rotation=30, fontsize=8)
        ax[j].set_title(f"{reg}: local seismicity rate ratio (exact 99 % range)", fontsize=9)
    ax[0].set_ylabel("rate in state / rate outside"); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/g2_q2_rates.png", dpi=110); plt.close(fig)


def fig_cases():
    """Two M >= 3.5 deep earthquakes near the families that fell inside between-ETS LFE bursts (post-hoc chance in
    RESULTS): hourly LFE counts (data, scrambled), PNSN tremor <= 100 km, local earthquakes, +-3 days."""
    cases = [("SVI", "2006-01-15 12:29:00", "SVI 2006-01-15 04:29 PST M3.6 near-interface, 7 km from a family"),
             ("WA", "2018-11-19 11:09:00", "WA-NW 2018-11-19 03:09 PST M4.1 intraslab (39 km deep), 10 km from a family")]
    A = pd.read_parquet(f"{T}/eq_local.parquet"); A = A[A.in_analysis]
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["tns"] = pd.to_datetime(tr.time).astype("datetime64[ns]").astype("int64")
    fig, ax = plt.subplots(2, 1, figsize=(12, 7))
    for a, (reg, t, lab) in zip(ax, cases):
        K, W, days, dP, ep = LT.load_lfe(reg)
        t0 = pd.Timestamp(t).value
        e = np.arange(-72, 73, 1)
        yd = [np.diff(W["data"][1][np.searchsorted(W["data"][0], [t0 + h * H, t0 + (h + 1) * H])])[0] for h in e[:-1]]
        yc = [np.diff(W["c1"][1][np.searchsorted(W["c1"][0], [t0 + h * H, t0 + (h + 1) * H])])[0] for h in e[:-1]]
        a.bar(e[:-1] + 0.5, yd, width=1, color="C0", label="LFEs, Σ P per hour (data)")
        a.bar(e[:-1] + 0.5, yc, width=1, color="C3", label="scrambled control")
        la, lo = (48.4837, -123.7485) if reg == "SVI" else (47.8432, -123.0282)
        dk = 111.195 * np.hypot(tr.lat - la, (tr.lon - lo) * np.cos(np.radians(la)))
        q = tr.tns[(dk <= 100) & (abs(tr.tns - t0) < 72 * H)]
        if len(q):
            a.plot((q - t0) / H, np.full(len(q), -1), "|", color="0.4", label=f"PNSN tremor ≤ 100 km ({len(q)})")
        x = A[(abs(A.t.astype("datetime64[ns]").astype("int64") - t0) < 72 * H) & (A[f"{reg}_dist_km"] <= 60) & (A.mag >= 1.5)]
        a2 = a.twinx()
        a2.scatter((x.t.astype("datetime64[ns]").astype("int64") - t0) / H, x.mag, color="C1", s=20, label="local EQs ≤ 60 km (M)")
        a2.set_ylabel("M"); a2.set_ylim(1, 5)
        a.axvline(0, color="C1", lw=1)
        a.set_title(lab, fontsize=9); a.set_ylabel("LFEs / h"); a.legend(fontsize=7, loc="upper left")
        a.set_xlabel("hours from the earthquake (PT labels in title)")
    fig.tight_layout(); fig.savefig(f"{FIG}/g2_cases_bursts.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    import sys as _s
    if len(_s.argv) > 1:
        for f in _s.argv[1:]:
            globals()[f]()
    else:
        fig_map(); fig_superposed(); fig_q2(); fig_cases()
