"""Figures for study D (fixed-station cross-era test). Reads svi_science/results/d_fixed_station/*.csv."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = "/home/jovyan/tremorferometry/svi_science/results/d_fixed_station"
FIG = "/home/jovyan/tremorferometry/svi_science/figures"
COL = {"prod": "#eb6834", "core9": "#2a78d6", "all10": "#1baf7a"}
LAB = {"prod": "production network (12/37/19 stations)", "core9": "fixed 9 boreholes", "all10": "fixed 9 boreholes + GNW"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})

R = pd.read_csv(f"{OUT}/era_rates.csv")
R = R[R.variant == "rev075"]
series = [("prod", "all"), ("core9", "all"), ("core9", "strict"), ("all10", "all"), ("all10", "strict")]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
panels = [("m1_per_famday", "net repeats per family-tremor-day\n(all tremor days; E1 incl. 2007–09 WECC windows)"),
          ("m1_per_famday_post2009", "net repeats per family-tremor-day\n(days from 2009-08-06 only)"),
          ("m2_per_epicentre", "net repeats per PNSN epicentre ≤ 30 km\n(from 2009-08-06)")]
for ax, (m, title) in zip(axs, panels):
    for j, (net, ds) in enumerate(series):
        x = R[(R.measure == m) & (R.network == net) & (R.days == ds)].set_index("era").loc[["E1", "E2", "E3"]]
        xs = np.arange(3) + (j - 2) * 0.13
        filled = ds == "strict" or net == "prod"
        ax.errorbar(xs, x.rate, yerr=[x.rate - x.lo95, x.hi95 - x.rate], fmt="o" if net != "prod" else "s",
                    color=COL[net], mfc=COL[net] if filled else "white", mec=COL[net], ms=7, lw=1.5, capsize=0,
                    label=f"{LAB[net]}" + (", strict days" if ds == "strict" else (", all days" if net != "prod" else "")))
        if net == "prod" or (net == "core9" and ds == "strict"):
            ax.plot(xs, x.rate, color=COL[net], lw=1, alpha=0.6)
    ax.set_xticks(range(3)); ax.set_xticklabels(["E1 2007–11", "E2 2012–14", "E3 2015–26"])
    ax.set_title(title, fontsize=9, color=INK); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_ylim(bottom=0)
axs[0].legend(fontsize=7, frameon=False, loc="lower left")
fig.suptitle("WA-NW LFE activity per era, P_own_WAVE ≥ 0.75 minus time-reversed control (95 % cluster bootstrap). "
             "Filled = network held fixed on every day (strict) or production", fontsize=9, color=INK)
fig.tight_layout()
fig.savefig(f"{FIG}/d_era_rates.png", dpi=150)

S = pd.read_csv(f"{OUT}/episode_series.csv")
S = S[S.major]
fig, axs = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
for ax, m, yl in ((axs[0], "m1", "net repeats per family-tremor-day"), (axs[1], "m2", "net repeats per PNSN epicentre ≤ 30 km")):
    for net, ds in (("prod", "all"), ("core9", "strict")):
        x = S[(S.network == net) & (S.days == ds) & (S.unit_days >= 100)]
        t = pd.to_datetime(x.start)
        ax.scatter(t, x[m], s=np.clip(x.unit_days / 25, 10, 90), color=COL[net], alpha=0.75 if net == "core9" else 0.55,
                   edgecolor="white", lw=0.5, label=f"{LAB[net]}" + (" (strict days)" if ds == "strict" else ""), zorder=3)
    for b in ("2012-01-01", "2015-01-01"):
        ax.axvline(pd.Timestamp(b), color=MUTED, lw=0.8, ls="--")
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_ylabel(yl); ax.grid(axis="y", color=GRID, lw=0.8)
axs[0].legend(fontsize=8, frameon=False, loc="upper right")
axs[0].set_title("Per MAJOR ETS episode (≥ 100 family-tremor-days; marker size = family-tremor-days). "
                 "Dashed = era boundaries; production network changes there, the fixed network does not", fontsize=9)
fig.tight_layout()
fig.savefig(f"{FIG}/d_episode_series.png", dpi=150)
print("ok")
