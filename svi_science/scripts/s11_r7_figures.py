#!/usr/bin/env python
"""s11_r7_figures -- S11 round 7 figures (PLAN.md §11.5): line graphs, points coloured by weight, no error bars on dv/v."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
BASE = "/home/jovyan/tremorferometry/svi_science/results/s11_depth"
R7 = f"{BASE}/round7"; FIG = f"{R7}/figures"; os.makedirs(FIG, exist_ok=True)
REG = [("svi", "SVI"), ("wa_nw", "WA-NW")]
BINS = ["le-1", "0", "1", "2", "3", "4-6", "ge7"]
BLAB = ["≤−1", "0", "1", "2", "3", "4–6", "≥7"]


def tb(sub, f):
    p = f"{R7}/{sub}/tables/{f}"
    return pd.read_csv(p) if os.path.exists(p) else None


def wscatter(ax, x, y, w, cmap="viridis", s=36, label=None):
    w = np.asarray(w, float); c = np.log10(np.clip(w, 1e-12, None))
    return ax.scatter(x, y, c=c, cmap=cmap, s=s, zorder=3, edgecolor="k", linewidth=0.3, label=label)


def fig_front():
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    for i, (sub, name) in enumerate(REG):
        cv = tb(sub, "r7_front_curves.csv"); sc = tb(sub, "r7_front_scr_curves.csv")
        if cv is None:
            continue
        x = np.arange(len(BINS))
        for j, (col, ttl) in enumerate((("shape", "within-episode shape"), ("absolute", "minus between-ETS baseline"))):
            ax = axs[i, j]
            if j == 0 and sc is not None:
                for d, g in sc.groupby("draw"):
                    g = g.set_index("bin").reindex(BINS)
                    ax.plot(x, g["shape"], color="0.75", lw=0.7, zorder=1)
            h = cv[cv.label == "headline"].set_index("bin").reindex(BINS)
            ax.plot(x, h[col], color="C3", lw=1.6, zorder=2, label="data (headline fronts)")
            sca = wscatter(ax, x, h[col], h.weight)
            c = cv[cv.label == "ctrl_timeofday"].set_index("bin").reindex(BINS)
            ax.plot(x, c[col], color="C0", lw=1.0, ls="--", zorder=2, label="time-of-day scrambled control")
            ax.axhline(0, color="k", lw=0.5)
            ax.set_xticks(x); ax.set_xticklabels(BLAB); ax.set_xlabel("days since the front reached the family")
            ax.set_ylabel("dv/v (%)"); ax.set_title(f"{name}: {ttl}" + (" (grey: front-scrambled draws)" if j == 0 else ""), fontsize=9)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, loc="best")
        fig.colorbar(sca, ax=axs[i, :], label="log10 weight", shrink=0.8)
    fig.savefig(f"{FIG}/fig_r7_front.png", dpi=150); plt.close(fig)


def fig_zones():
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    for i, (sub, name) in enumerate(REG):
        mids = pd.read_csv(f"{BASE}/round6/{sub}/tables/r6_sizes.csv").set_index("episode").mid_PT
        Umid = pd.to_datetime(mids)
        pm = tb(sub, "r7_zones_permajor.csv"); sl = tb(sub, "r7_zones_slip_permajor.csv")
        if pm is None:
            continue
        pm = pm[pm.setup == "two3_20"].copy()
        pm["t"] = pm.major.map(Umid)
        pm = pm.sort_values("t")
        ax = axs[i, 0]
        zs = [c for c in pm.columns if c.startswith("z") and "-" not in c and not c.endswith("_sig")]
        for k, z in enumerate(zs):
            ax.plot(pm.t, pm[z], lw=1.0, color=f"C{k}", label=f"zone {z} (SE→NW)", zorder=2)
            sca = wscatter(ax, pm.t, pm[z], 1 / pm[f"{z}_sig"] ** 2, s=26)
        ax.axhline(0, color="k", lw=0.5); ax.set_ylabel("ΔD per zone (dv/v %)"); ax.set_title(f"{name}: zone ΔD per major (E4, two3_20)", fontsize=9)
        ax.legend(fontsize=7)
        ax = axs[i, 1]
        if sl is not None and len(sl):
            s2 = sl[sl.setup == "two3_20"].copy(); s2["t"] = s2.major.map(Umid)
            for k, (lab, g) in enumerate(s2.groupby("labels")):
                g = g.sort_values("t")
                ax.plot(g.t, g.contrast, lw=1.0, color=f"C{k + 3}", label=f"labels {lab}", zorder=2)
                sca = wscatter(ax, g.t, g.contrast, 1 / g.sig ** 2, s=30)
        ax.axhline(0, color="k", lw=0.5); ax.set_ylabel("slipped − not slipped (dv/v %)"); ax.set_title(f"{name}: slipped minus not-slipped zones", fontsize=9)
        ax.legend(fontsize=7)
        fig.colorbar(sca, ax=axs[i, :], label="log10 weight (1/σ²)", shrink=0.8)
    fig.savefig(f"{FIG}/fig_r7_zones.png", dpi=150); plt.close(fig)


def fig_joint():
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), constrained_layout=True)
    for i, (sub, name) in enumerate(REG):
        su = tb(sub, "r7_joint_summary.csv"); de = tb(sub, "r7_joint_decade.csv")
        if su is None:
            continue
        ax = axs[2 * i]
        for k, f in enumerate(["coda", "ballistic", "joint", "joint_dS", "joint_Vs"]):
            t = su[su.fit == f].set_index("group").reindex(["s", "m", "d"])
            ax.plot([0, 1, 2], t.sig_major_med, marker="o", color=f"C{k}", label=f)
        ax.set_yscale("log"); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["shallow 0–5", "mid 5–25", "deep 25–45"])
        ax.set_ylabel("σ per major (%)"); ax.set_title(f"{name}: per-episode σ", fontsize=9); ax.legend(fontsize=7)
        ax = axs[2 * i + 1]
        if de is not None:
            for k, f in enumerate(["coda", "ballistic", "joint"]):
                t = de[de.fit == f].set_index("group").reindex(["s", "m", "d"])
                ax.plot([0, 1, 2], t.sig, marker="o", color=f"C{k}", label=f)
        ax.set_yscale("log"); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["shallow", "mid", "deep"])
        ax.set_ylabel("σ decade step (%)"); ax.set_title(f"{name}: decade (≥2015 − ≤2014) σ", fontsize=9); ax.legend(fontsize=7)
    fig.savefig(f"{FIG}/fig_r7_joint.png", dpi=150); plt.close(fig)


def fig_hint():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    labs = ["full", "fam_odd", "fam_even", "sta_A", "sta_B"]
    for j, meas in enumerate(["Du", "T"]):
        ax = axs[j]; xt = []; k = 0
        for sub, name in REG:
            t = tb(sub, "r7_hint_contrasts.csv")
            if t is None:
                continue
            t = t[t.measure == meas].set_index("label")
            xs = []
            for lab in labs:
                if lab not in t.index:
                    continue
                r = t.loc[lab]
                sca = wscatter(ax, [k], [r.contrast], [1 / r.sig ** 2], s=50)
                ax.plot([k - 0.3, k + 0.3], [r.scr_c95, r.scr_c95], color="0.6", lw=0.8)
                ax.plot([k - 0.3, k + 0.3], [-r.scr_c95, -r.scr_c95], color="0.6", lw=0.8)
                xt.append((k, f"{name}\n{lab}")); xs.append((k, r.contrast)); k += 1
            if xs:
                ax.plot([a for a, _ in xs], [b for _, b in xs], color="C3", lw=0.8, zorder=2)
            k += 1
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks([a for a, _ in xt]); ax.set_xticklabels([b for _, b in xt], fontsize=6, rotation=45)
        ax.set_ylabel("top third − rest, ΔD (dv/v %)"); ax.set_title(f"top third by {'duration' if meas == 'Du' else 'tremor'} (grey: scrambled 95 %)", fontsize=9)
    fig.colorbar(sca, ax=axs, label="log10 weight (1/σ²)", shrink=0.8)
    fig.savefig(f"{FIG}/fig_r7_hint.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    for f in (fig_front, fig_zones, fig_joint, fig_hint):
        try:
            f()
        except Exception as ex:
            import traceback; print("FIGURE ERROR", f.__name__, ex, traceback.format_exc()[-600:])
    open(f"{R7}/FIGURES_DONE", "w").write(pd.Timestamp.now(tz="America/Los_Angeles").strftime("%Y-%m-%d %-I:%M %p PT") + "\n")
