#!/usr/bin/env python
"""q4_figures (2026-10-07) -- recharge curves by era (WA-NW) with the SVI overlay.
    python svi_science/scripts/q4_figures.py
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q4_recharge_wa"
FIG = f"{ROOT}/svi_science/figures/q4_recharge_curves.png"
COL = {"E1": "#2a78d6", "E2": "#eb6834", "E3": "#1baf7a", "pooled": "#0b0b0b", "SVI": "#e87ba4", "SVIq4": "#4a3aa7"}
INK2 = "#52514e"

C = pd.concat([pd.read_csv(f"{OUT}/recharge_curves_wa.csv"), pd.read_csv(f"{OUT}/recharge_curves_posthoc_daycount.csv")])
O = pd.concat([pd.read_csv(f"{OUT}/recharge_fits_wa.csv"), pd.read_csv(f"{OUT}/recharge_fits_posthoc_daycount.csv")]).set_index("name")
P3 = pd.read_csv(f"{ROOT}/svi_science/results/p3_bostock/recharge_curves.csv")
svi = P3[P3.set == "LFE deep groups 1/4/10/12 | ref big list | D1"].copy()
svi["norm_rate"] = svi.rate / (svi.events.sum() / svi.exposure.sum())
EDGES = [0.5, 7.5, 14.5, 30.5, 60.5, 90.5, 150.5, 240.5, 400.5]
MID = np.sqrt(np.array(EDGES[:-1]) * np.array(EDGES[1:]))


def lab(name, short):
    o = O.loc[name]
    p = o.get("p_D30", np.nan)
    ps = f", p {p:.3f}" if np.isfinite(p) else ""
    return f"{short}: D30 {o.D30:+.2f}{ps}, β {o.beta:+.2f}"


def panel(ax, prefix, title, svi_on):
    for e in ("E1", "E2", "E3"):
        nm = f"{prefix} | {e}"
        g = C[C.set == nm]
        if not len(g):
            continue
        ax.plot(MID, g.norm_rate, "-o", color=COL[e], lw=2, ms=6, label=lab(nm, e))
    pn = [s for s in C.set.unique() if s.startswith(prefix + " | pooled") or s.startswith(prefix + " | 2009-08")]
    if pn:
        g = C[C.set == pn[0]]
        ax.plot(MID, g.norm_rate, "-s", color=COL["pooled"], lw=2.5, ms=6, label=lab(pn[0], "all eras" if "pooled-n" in pn[0] else "2009–26"))
    if svi_on:
        ax.plot(MID, svi.norm_rate, "--D", color=COL["SVI"], lw=2, ms=6, label="SVI P3 (big clock, as coded)")
        nm = "SVI deep groups | clock big | D1 | drop all majors (Q4 definition)"
        g = C[C.set == nm]
        g = g.assign(norm_rate=g.rate / (g.events.sum() / g.exposure.sum()))
        ax.plot(MID, g.norm_rate, ":D", color=COL["SVIq4"], lw=2, ms=6, label=lab(nm, "SVI, all majors dropped"))
    ax.axhline(1, color=INK2, lw=0.8, alpha=0.6)
    ax.axvline(30.5, color=INK2, lw=0.8, ls=":", alpha=0.6)
    ax.set_xscale("log"); ax.set_xlim(0.8, 420)
    ax.set_xticks([1, 7, 14, 30, 60, 90, 150, 240, 400]); ax.set_xticklabels(["1", "7", "14", "30", "60", "90", "150", "240", "400"])
    ax.set_xlabel("days since the reference ETS ended")
    ax.set_title(title, fontsize=10, loc="left")
    ax.grid(alpha=0.2); ax.legend(fontsize=7, loc="upper left", frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


fig, axs = plt.subplots(2, 2, figsize=(13, 9.5), sharey=False)
panel(axs[0, 0], "LFE deep | clock big | D1", "a) WA deep LFEs, big-ETS clock (headline), with SVI", True)
panel(axs[0, 1], "LFE deep | clock all majors | D1", "b) WA deep LFEs, all-majors clock", False)
panel(axs[1, 0], "PNSN deep tremor | clock big | D1", "c) PNSN deep tremor near WA deep families, big-ETS clock", False)
panel(axs[1, 1], "POSTHOC active days | clock big | D1", "d) not pre-registered: days with >= 1 deep LFE (a burst counts once)", False)
g = C[C.set == "POSTHOC SVI active days | clock big | D1 | drop all majors"]
g = g.assign(norm_rate=g.rate / (g.events.sum() / g.exposure.sum()))
axs[1, 1].plot(MID, g.norm_rate, ":D", color=COL["SVIq4"], lw=2, ms=6, label=lab("POSTHOC SVI active days | clock big | D1 | drop all majors", "SVI, all majors dropped"))
axs[1, 1].legend(fontsize=7, loc="upper left", frameon=False)
for ax in axs[:, 0]:
    ax.set_ylabel("rate ÷ that curve's mean rate")
for ax in axs.flat:
    ax.set_ylim(0, None)
fig.suptitle("Q4: deep activity outside major ETS vs time since the last big ETS — Washington-NW 2007–2026", fontsize=12)
fig.tight_layout()
fig.savefig(FIG, dpi=140)
print("wrote", FIG)
