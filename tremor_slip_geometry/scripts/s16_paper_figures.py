#!/usr/bin/env python
"""S16 PART B — the paper figure set (Figs 1-4), Tables 1-2, from the SAVED results only.

Everything plotted here is read from tremor_slip_geometry/results/*.csv written by S1-S16.
Nothing is re-derived except (i) group medians / IQRs over per-event exponents that the earlier
steps saved individually (marked "aggregated here" in the captions and in the `source` columns of
Table 2), and (ii) the daily (R7, C, E) tracks of Fig 2a, which no earlier step saved and which
`scripts/s16_daily_tracks.py` regenerated with the S16 PART A code (validated: every recomputed
per-event exponent reproduces the stored S9/S11/S16 value to 8e-17).

INPUTS (all inside tremor_slip_geometry/results/)
  Fig 1  s7_occupancy_triangles.csv, s1_events.csv, s1_controls.csv, s1_trends.csv,
         s2_yearly.csv, s2_trends.csv
  Fig 2  s16_daily_tracks.csv, s9_rate_scaling.csv, s11_exponents.csv, s11_energy_events.csv,
         s11b_exponents.csv, s15_cand_exponents.csv, s16_fb_aggregation.csv, s16_fb_events.csv
  Fig 3  s4_fronts_events.csv, s4_fronts_summary.csv, s9b_phase10.csv, s9b_downdip_stack.csv,
         s9b_sensitivity.csv, s9_depth_migration_stack.csv
  Fig 4  s12_representability_triangles.csv, s12_misfit_by_latitude.csv, s7_occupancy_summary.csv,
         s12b_data_vs_model.csv, s6_interets_conf.csv, s3_profiles.csv
  Tables s1_events.csv, s15_events_check.csv, s9_rate_scaling.csv, s11_energy_events.csv,
         s16_fb_events.csv, s4_fronts_events.csv  (Table 1);  the exponent CSVs above (Table 2)

OUTPUTS (new files only)
  figures/s16_paper_fig1.png .. s16_paper_fig4.png   300 dpi, colourblind-safe (Okabe-Ito)
  results/s16_table1.csv                             the 27 confirmed events
  results/s16_table2.csv                             every exponent with CI, null and source

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s16_paper_figures.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path("/home/jovyan/tremorferometry/tremor_slip_geometry")
RES = HERE / "results"
FIG = HERE / "figures"

# ---- Okabe-Ito colourblind-safe palette -------------------------------------------------
BLUE = "#0072B2"      # north / thrust / observed
VERM = "#D55E00"      # central
GREEN = "#009E73"     # second measure
ORANGE = "#E69F00"    # third measure
SKY = "#56B4E9"       # fourth measure
PURPLE = "#CC79A7"    # F&B reference
GREY = "#999999"      # nulls / controls
DARK = "#333333"
SEGC = {"north": BLUE, "central": VERM, "south": GREY}

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.linewidth": 0.8, "axes.grid": True, "grid.alpha": 0.25,
    "grid.linewidth": 0.5, "lines.linewidth": 1.3, "figure.dpi": 110,
    "savefig.dpi": 300, "savefig.bbox": "tight", "axes.titlelocation": "left",
    "axes.titleweight": "bold", "legend.frameon": False,
})
LAT_ASPECT = 1.0 / np.cos(np.deg2rad(47.0))


def rd(name):
    return pd.read_csv(RES / name)


def panel(ax, letter, title):
    ax.set_title(f"{letter}  {title}", loc="left")


def decyear(s):
    t = pd.to_datetime(pd.Series(s)).values.astype("datetime64[ns]")
    y = t.astype("datetime64[Y]").astype(int) + 1970
    y0 = t.astype("datetime64[Y]").astype("datetime64[ns]")
    y1 = (t.astype("datetime64[Y]") + 1).astype("datetime64[ns]")
    return y + (t - y0) / (y1 - y0)


def theil_line(ax, x, y, slope, color, lw=1.4, ls="-", label=None):
    xm, ym = np.nanmedian(x), np.nanmedian(y)
    xx = np.array([np.nanmin(x), np.nanmax(x)])
    ax.plot(xx, ym + slope * (xx - xm), color=color, lw=lw, ls=ls, label=label, zorder=4)


def strip(ax, xpos, v, color, width=0.16, rng=None, ms=14):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    rng = rng or np.random.default_rng(7)
    ax.scatter(xpos + rng.uniform(-width, width, v.size), v, s=ms, color=color,
               alpha=0.75, edgecolor="none", zorder=3)
    ax.plot([xpos - 0.28, xpos + 0.28], [np.median(v)] * 2, color=DARK, lw=2.0, zorder=5)
    q1, q3 = np.percentile(v, [25, 75])
    ax.plot([xpos, xpos], [q1, q3], color=DARK, lw=0.9, zorder=4)
    return float(np.median(v)), v.size


# =====================================================================================
# FIGURE 1 — co-registration and depth occupancy (S1, S2, S7)
# =====================================================================================
tri7 = rd("s7_occupancy_triangles.csv")
ev = rd("s1_events.csv")
ctl = rd("s1_controls.csv")
trd = rd("s1_trends.csv")
yr = rd("s2_yearly.csv")
ytr = rd("s2_trends.csv")

cmap_c = np.log10(tri7["N_det_tot"].values + 1.0)
dy = decyear(ev["start"])
rQ3 = trd[(trd.set == "confirmed") & (trd.seg == "combined")
          & (trd.quantity == "Q3_depth_km")].iloc[0]

fig = plt.figure(figsize=(7.4, 8.8), layout="constrained")
gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.25], height_ratios=[1.15, 0.78, 0.78])

# (a) map
ax = fig.add_subplot(gs[:2, 0])
sc = ax.scatter(tri7["lon"], tri7["lat"], c=cmap_c, s=12, marker="s", cmap="viridis",
                vmin=0, vmax=np.nanpercentile(cmap_c, 99), linewidths=0)
ss = 16 + 240 * ev["potency_km3"] / ev["potency_km3"].max()
ax.scatter(ev["slip_lon"], ev["slip_lat"], s=ss, facecolor="none", edgecolor="#FFFFFF",
           linewidths=1.5, zorder=5)
ax.scatter(ev["slip_lon"], ev["slip_lat"], s=ss, facecolor="none", edgecolor=VERM,
           linewidths=0.9, zorder=6)
for lat in (44.0, 46.0):
    ax.axhline(lat, color=DARK, lw=0.8, ls="--", alpha=0.55, zorder=4)
bb = dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=0.85)
ax.text(-127.9, 49.3, "north", color=DARK, fontsize=7.5, ha="left", bbox=bb)
ax.text(-127.9, 45.1, "central", color=DARK, fontsize=7.5, ha="left", bbox=bb)
ax.text(-127.9, 41.3, "south (excluded)", color=DARK, fontsize=7.5, ha="left", bbox=bb)
ax.set_aspect(LAT_ASPECT)
ax.set_xlabel("longitude (deg E)")
ax.set_ylabel("latitude (deg N)")
ax.set_xticks([-127, -125, -123])
ax.grid(alpha=0.15)
cb = fig.colorbar(sc, ax=ax, location="bottom", shrink=0.92, pad=0.02, aspect=26)
cb.set_label("log$_{10}$ PNSN detections per triangle, 2010-2026", fontsize=7)
cb.ax.tick_params(labelsize=7)
panel(ax, "a", "tremor and the 27 events")

# (d) completeness control
ax = fig.add_subplot(gs[2, 0])
for seg, g in yr[~yr.partial].groupby("seg"):
    if seg == "south":
        continue
    ax.plot(g["year"], g["completeness"], marker="o", ms=3, color=SEGC[seg], label=seg)
ax.axvline(2019, color=GREY, lw=0.8, ls="--")
ax.set_xlabel("year")
ax.set_ylabel("completeness proxy")
ax.set_xticks([2010, 2014, 2018, 2022, 2025])
ax.legend(loc="upper left", ncol=2)
panel(ax, "d", "detection growth (control)")

# (b) thrust vs mirror
ax = fig.add_subplot(gs[0, 1])
m1, n1 = strip(ax, 0, ev["Q1_downdip_km"], BLUE, rng=np.random.default_rng(1))
m2, n2 = strip(ax, 1, ctl["Q1_downdip_km"], GREY, rng=np.random.default_rng(2))
m3, n3 = strip(ax, 2.5, ev["Q3_depth_km"], BLUE, rng=np.random.default_rng(3))
m4, n4 = strip(ax, 3.5, ctl["Q3_depth_km"], GREY, rng=np.random.default_rng(4))
ax.axhline(0, color=DARK, lw=0.8, ls=":")
ax.set_xticks([0, 1, 2.5, 3.5])
ax.set_xticklabels([f"thrust\nn={n1}", f"back-slip\nn={n2}", f"thrust\nn={n3}",
                    f"back-slip\nn={n4}"], fontsize=7.5)
ax.set_ylabel("tremor $-$ slip centroid (km)")
ax.set_ylim(-95, 78)
ax.text(0.22, 0.97, "across-strike Q1", transform=ax.transAxes, ha="center", va="top",
        fontsize=8, color=DARK)
ax.text(0.83, 0.97, "depth Q3", transform=ax.transAxes, ha="center", va="top",
        fontsize=8, color=DARK)
for x, m in zip([0, 1, 2.5, 3.5], [m1, m2, m3, m4]):
    ax.annotate(f"{m:.1f}", (x, m), textcoords="offset points", xytext=(11, 3),
                fontsize=7.5, color=DARK)
panel(ax, "b", "co-registration vs back-slip control")

# (c) Q3 vs time
ax = fig.add_subplot(gs[1, 1])
for seg, g in ev.groupby("seg"):
    ax.scatter(decyear(g["start"]), g["Q3_depth_km"], s=26, color=SEGC[seg], label=seg, zorder=3)
theil_line(ax, dy, ev["Q3_depth_km"].values, rQ3["slope"], DARK, ls="--")
ax.axhline(0, color=GREY, lw=0.8, ls=":")
ax.set_xlabel("event start (year)")
ax.set_ylabel("Q3 depth offset (km)")
ax.set_ylim(-4.8, 7.6)
ax.text(0.02, 0.03, f"Theil-Sen {rQ3['slope']:+.2f} km/yr [{rQ3['lo']:+.2f}, {rQ3['hi']:+.2f}], "
                    f"p$_{{perm}}$ {rQ3['p_perm']:.2f}; 80 % MDD 0.2 km/yr",
        transform=ax.transAxes, fontsize=7, color=DARK)
ax.legend(loc="upper left", ncol=2)
panel(ax, "c", "no drift of the depth offset")

# (e) depth occupancy
ax = fig.add_subplot(gs[2, 1])
gn = yr[(yr.seg == "north") & (~yr.partial)]
gc = yr[(yr.seg == "central") & (~yr.partial)]
ax.fill_between(gn["year"], gn["p10"], gn["p90"], color=BLUE, alpha=0.18,
                label="north tremor p10-p90")
ax.plot(gn["year"], gn["p50"], marker="o", ms=3, color=BLUE, label="north tremor median")
ax.plot(gn["year"], gn["slip_depth"], marker="s", ms=3, color=GREEN, ls="--",
        label="north slip-weighted")
ax.plot(gc["year"], gc["p50"], marker="^", ms=3, color=VERM, label="central tremor median")
ax.invert_yaxis()
ax.set_xlabel("year")
ax.set_ylabel("Slab2 depth (km)")
pn = ytr[(ytr.seg == "north") & (ytr.quantity == "p50")].iloc[0]
ax.text(0.02, 0.03, f"north median trend {pn['theilsen_slope_km_per_yr']:+.3f} km/yr "
                    f"[{pn['lo90']:+.3f}, {pn['hi90']:+.3f}]", transform=ax.transAxes,
        fontsize=7, color=DARK)
ax.legend(loc="upper left", ncol=2, fontsize=6.3)
panel(ax, "e", "depth occupancy per year")

fig.savefig(FIG / "s16_paper_fig1.png")
plt.close(fig)
print("wrote figures/s16_paper_fig1.png")


# =====================================================================================
# FIGURE 2 — tremor as a sub-linear function of slip rate (S9, S11, S11b, S15B, S16)
# =====================================================================================
trk = rd("s16_daily_tracks.csv")
s9 = rd("s9_rate_scaling.csv")
s9e = s9[s9.id > 0]
s11x = rd("s11_exponents.csv")
s11e = rd("s11_energy_events.csv")
s11b = rd("s11b_exponents.csv")
s15c = rd("s15_cand_exponents.csv")
s16a = rd("s16_fb_aggregation.csv")
s16e = rd("s16_fb_events.csv")

FB = 1.0 / 3.17
FB_LO, FB_HI = 1.0 / 3.46, 1.0 / 2.88


def s16row(meas, direction="forward_measure_on_rate", level="per_event_theilsen"):
    r = s16a[(s16a.measure == meas) & (s16a.direction == direction) & (s16a.level == level)]
    return r.iloc[0]


fig = plt.figure(figsize=(7.4, 8.2), layout="constrained")
gs = fig.add_gridspec(3, 2)

# (a) example daily tracks
ax = fig.add_subplot(gs[0, 0])
for pid, col in [(1863, BLUE), (1601, GREEN), (1361, VERM)]:
    g = trk[(trk.id == pid) & trk.in_fit]
    x = np.log10(g["R7"].values)
    y = np.log10(g["C"].values + 1.0)
    b = float(s9e.loc[s9e.id == pid, "exponent"].iloc[0])
    st = str(g["ev_start"].iloc[0])[:7]
    sg = str(g["seg"].iloc[0])[0].upper()
    ax.scatter(x, y, s=12, color=col, alpha=0.8, edgecolor="none")
    theil_line(ax, x, y, b, col, label=f"{pid} {st} ({sg})  b={b:.2f}")
ax.set_xlabel("log$_{10}$ potency rate R7 (mm km$^2$ d$^{-1}$)")
ax.set_ylabel("log$_{10}$ (tremor count + 1)")
ax.legend(loc="lower right", fontsize=6.5)
panel(ax, "a", "three example events (daily)")

# (b) per-event count exponents
ax = fig.add_subplot(gs[0, 1])
g = s9e.sort_values("exponent").reset_index(drop=True)
rC = s11x[(s11x.measure == "C") & (s11x.slip_var == "R7") & (s11x.subset == "all")].iloc[0]
ax.axvspan(rC["null_mean"] - rC["null_sd"], rC["null_mean"] + rC["null_sd"], color=GREY,
           alpha=0.25, zorder=1, label="shifted null $\\pm$1 sd")
ax.axvspan(rC["b_boot_lo90"], rC["b_boot_hi90"], color=DARK, alpha=0.12, zorder=1)
ax.axvline(rC["b_median"], color=DARK, lw=1.4, label=f"median {rC['b_median']:.2f}")
ax.axvline(1.0, color=PURPLE, lw=1.0, ls="--", label="b = 1 (proportional)")
ax.axvline(0.0, color=GREY, lw=0.8, ls=":")
ax.scatter(g["exponent"], np.arange(len(g)), s=20, c=[SEGC[s] for s in g["seg"]], zorder=3)
ax.set_yticks([])
ax.set_ylim(-2.5, len(g) + 1)
ax.set_xlabel("per-event exponent b (count vs R7)")
ax.set_ylabel("27 events, sorted")
ax.text(0.03, 0.03, f"{int(rC['n_below_1'])}/27 below 1; {int(rC['n_z_positive'])}/"
                    f"{int(rC['n_events_z'])} above null (p {rC['sign_p']:.0e})",
        transform=ax.transAxes, fontsize=7, color=DARK)
ax.legend(loc="upper left", fontsize=6.5)
panel(ax, "b", "count exponent, event by event")

# (c) measure ladder
ax = fig.add_subplot(gs[1, 0])
lad = [
    ("count, 27 ev (S11)", rC["b_median"], rC["b_boot_lo90"], rC["b_boot_hi90"], BLUE),
    ("count, 18 ev (S16)", *[s16row("C")[k] for k in ("b_median", "b_boot_lo90", "b_boot_hi90")],
     BLUE),
    ("summed energy (S11)", *[s16row("E_sum_epsS11")[k] for k in
                              ("b_median", "b_boot_lo90", "b_boot_hi90")], GREEN),
    ("summed amplitude", *[s16row("Asum_lin")[k] for k in
                           ("b_median", "b_boot_lo90", "b_boot_hi90")], ORANGE),
    ("median energy / det.", *[s16row("E_med")[k] for k in
                               ("b_median", "b_boot_lo90", "b_boot_hi90")], SKY),
    ("F&B amplitude A1", *[s16row("A1_sqrt_Emed")[k] for k in
                           ("b_median", "b_boot_lo90", "b_boot_hi90")], PURPLE),
    ("median magnitude", *[s16row("M_med")[k] for k in
                           ("b_median", "b_boot_lo90", "b_boot_hi90")], GREY),
]
ypos = np.arange(len(lad))[::-1]
for yp, (lab, bm, lo, hi, col) in zip(ypos, lad):
    ax.errorbar(bm, yp, xerr=[[bm - lo], [hi - bm]], fmt="o", ms=5, color=col, capsize=2.5,
                lw=1.2)
ax.axvline(1.0, color=DARK, lw=1.0, ls="--")
ax.axvline(0.0, color=GREY, lw=0.8, ls=":")
ax.axvspan(FB_LO, FB_HI, color=PURPLE, alpha=0.15)
ax.axvline(FB, color=PURPLE, lw=1.0, ls="-.")
ax.annotate("F&B 1/3.17", (FB, -0.45), fontsize=6.5, color=PURPLE, ha="center", va="center")
ax.set_yticks(ypos)
ax.set_yticklabels([l[0] for l in lad], fontsize=7.5)
ax.set_ylim(-0.9, len(lad) - 0.3)
ax.set_xlabel("exponent b (measure vs R7), median, 90 % CI")
ax.set_xlim(-0.05, 1.15)
panel(ax, "c", "one relation, four aggregations")

# (d) model-free check
ax = fig.add_subplot(gs[1, 1])
prox = [("vbICA R7", "R7_vbICA"), ("ICA displ.\n31 d", "Gica31_proj"),
        ("raw GNSS\n31 d", "Graw31_proj"), ("raw GNSS\n1 d", "Graw_proj")]
w = 0.32
for i, (meas, mk, col) in enumerate([("C", "o", BLUE), ("E", "s", GREEN)]):
    for j, (lab, p) in enumerate(prox):
        rr = s11b[(s11b.measure == meas) & (s11b.proxy == p) & (s11b.subset == "all")].iloc[0]
        ax.errorbar(j + (i - 0.5) * w, rr["b_median"],
                    yerr=[[rr["b_median"] - rr["b_boot_lo90"]],
                          [rr["b_boot_hi90"] - rr["b_median"]]],
                    fmt=mk, ms=5, color=col, capsize=2.5, lw=1.2)
ax.axhline(1.0, color=DARK, lw=1.0, ls="--")
ax.axhline(0.0, color=GREY, lw=0.8, ls=":")
ax.set_xticks(range(len(prox)))
ax.set_xticklabels([p[0] for p in prox], fontsize=7)
ax.set_ylabel("exponent b, median and 90 % CI")
ax.set_ylim(-0.05, 1.32)
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="", label="tremor count"),
                   Line2D([], [], color=GREEN, marker="s", ls="", label="tremor energy")],
          loc="upper left", fontsize=7)
panel(ax, "d", "model-free regressors agree (S11b)")

# (e) constancy
ax = fig.add_subplot(gs[2, 0])
groups = [("north", s9e[s9e.seg == "north"]["exponent"], BLUE),
          ("central", s9e[s9e.seg == "central"]["exponent"], VERM),
          ("2010-14", s9e[s9e.epoch == "2010-2014"]["exponent"], DARK),
          ("2015-19", s9e[s9e.epoch == "2015-2019"]["exponent"], DARK),
          ("2020-26", s9e[s9e.epoch == "2020-2026"]["exponent"], DARK)]
rng = np.random.default_rng(11)
for i, (lab, v, col) in enumerate(groups):
    strip(ax, i, v, col, rng=rng, ms=16)
ax.axhline(rC["b_median"], color=DARK, lw=1.0, ls="--", label="all 27 events")
ax.axhline(1.0, color=PURPLE, lw=0.9, ls=":")
ax.set_xticks(range(len(groups)))
ax.set_xticklabels([f"{g[0]}\nn={np.isfinite(g[1]).sum()}" for g in groups], fontsize=7)
ax.set_ylabel("per-event exponent b (count)")
ax.legend(loc="lower right", fontsize=7)
panel(ax, "e", "no segment or epoch dependence")

# (f) F&B aggregation ladder
ax = fig.add_subplot(gs[2, 1])
xs, ys, los, his = [], [], [], []
for i, (lab, lv) in enumerate([("1 d", "bin1d_theilsen"), ("3 d", "bin3d_theilsen"),
                               ("7 d", "bin7d_theilsen")]):
    rr = s16row("A1_sqrt_Emed", "forward_measure_on_rate", lv)
    xs.append(i)
    ys.append(rr["b_median"])
    los.append(rr["b_median"] - rr["b_boot_lo90"])
    his.append(rr["b_boot_hi90"] - rr["b_median"])
ce = s16row("A1_sqrt_Emed", "cross_event_amp_on_rate", "n18_theilsen")
ax.errorbar(xs, ys, yerr=[los, his], fmt="o-", ms=5, color=PURPLE, capsize=2.5, lw=1.3,
            label="forward, 90 % CI")
ax.plot([3], [ce["b_median"]], marker="D", ms=6, color=PURPLE, ls="none",
        label="cross-event (n = 18)")
rma = s16row("A1_sqrt_Emed", "rma_amp_on_rate", "per_event_rma")
ax.errorbar([4], [rma["b_median"]],
            yerr=[[rma["b_median"] - rma["b_boot_lo90"]], [rma["b_boot_hi90"] - rma["b_median"]]],
            fmt="^", ms=6, color=GREEN, capsize=2.5, lw=1.2, label="RMA, daily")
rv = s16row("A1_sqrt_Emed", "reverse_rate_on_amplitude_FB", "per_event_theilsen")
ax.plot([5], [1.0 / rv["b_median"]], marker="v", ms=6, color=ORANGE, ls="none",
        label="1/beta (F&B direction), daily")
ax.axhspan(FB_LO, FB_HI, color=PURPLE, alpha=0.15)
ax.axhline(FB, color=PURPLE, lw=1.0, ls="-.")
ax.annotate("F&B 1/3.17", (5.35, FB), fontsize=6.5, color=PURPLE, ha="right", va="bottom")
ax.set_xticks(range(6))
ax.set_xticklabels(["1 d", "3 d", "7 d", "event", "RMA", "1/$\\beta$"], fontsize=7)
ax.set_xlim(-0.5, 5.5)
ax.set_xlabel("aggregation / estimator")
ax.set_ylabel("amplitude exponent (A1 vs R7)")
ax.set_ylim(-0.05, 1.2)
ax.legend(loc="upper left", fontsize=6.2)
panel(ax, "f", "matching Frank & Brodsky")

fig.savefig(FIG / "s16_paper_fig2.png")
plt.close(fig)
print("wrote figures/s16_paper_fig2.png")


# =====================================================================================
# FIGURE 3 — within-event geometry (S4, S9, S9b)
# =====================================================================================
fr = rd("s4_fronts_events.csv")
frs = rd("s4_fronts_summary.csv")
ph = rd("s9b_phase10.csv")
dd = rd("s9b_downdip_stack.csv")
sn = rd("s9b_sensitivity.csv")
dms = rd("s9_depth_migration_stack.csv")


def frac(q, group="combined"):
    r = frs[(frs.group == group) & (frs.quantity == q)]
    return float(r["median"].iloc[0])


fig = plt.figure(figsize=(7.4, 6.4), layout="constrained")
gs = fig.add_gridspec(2, 3)

# (a) lead/lag
ax = fig.add_subplot(gs[0, 0])
ax.hist(fr["L_days"], bins=np.arange(-25.5, 26.5, 3), color=BLUE, alpha=0.9, edgecolor="w",
        linewidth=0.4)
ax.axvspan(-5, 5, color=GREEN, alpha=0.15)
ax.axvline(0, color=DARK, lw=0.9, ls=":")
ax.set_ylim(0, 10.5)
ax.set_xlabel("slip $-$ tremor lag L (d)")
ax.set_ylabel("events")
ax.text(0.03, 0.97, f"|L| $\\leq$ 5 d: {frac('frac_absL_le5_observed')*100:.0f} % obs\n"
                    f"vs {frac('frac_absL_le5_null')*100:.1f} % null\n"
                    f"binomial p {frac('p_binom_absL_le5_obs_gt_null'):.0e}",
        transform=ax.transAxes, va="top", fontsize=7, color=DARK)
panel(ax, "a", "day-scale synchrony")

# (b) speeds
ax = fig.add_subplot(gs[0, 1])
for seg, g in fr.groupby("seg"):
    ax.scatter(g["v_t_km_d"], g["v_s_km_d"], s=22, color=SEGC[seg], label=seg, zorder=3)
lim = 11
ax.plot([-lim, lim], [-lim, lim], color=GREY, lw=0.9, ls="--")
ax.axhline(0, color=GREY, lw=0.7, ls=":")
ax.axvline(0, color=GREY, lw=0.7, ls=":")
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xlabel("tremor front $v_t$ (km/d)")
ax.set_ylabel("slip centroid $v_s$ (km/d)")
ax.text(0.03, 0.97, f"direction agrees {frac('dir_agree_frac')*100:.0f} %\n"
                    f"$v_s/v_t$ {frac('v_ratio'):.2f} (mixture)",
        transform=ax.transAxes, va="top", fontsize=7, color=DARK)
ax.legend(loc="lower right", fontsize=7)
panel(ax, "b", "same direction")

# (c) depth vs phase (10 bins)
ax = fig.add_subplot(gs[0, 2])
g = ph[ph.series == "depth"]
xc = 0.5 * (g["phase_lo"] + g["phase_hi"])
ax.fill_between(xc, g["null_lo90"], g["null_hi90"], color=GREY, alpha=0.3, label="null 90 %")
ax.errorbar(xc, g["obs_median"], yerr=[g["obs_median"] - g["boot_lo90"],
                                       g["boot_hi90"] - g["obs_median"]],
            fmt="o-", ms=4, color=BLUE, capsize=2, lw=1.2, label="observed")
ax.axhline(0, color=DARK, lw=0.8, ls=":")
ax.set_xlabel("phase within event")
ax.set_ylabel("depth anomaly (km)")
ax.legend(loc="upper left", fontsize=7)
panel(ax, "c", "the tail is deeper")

# (d) stacked migration
ax = fig.add_subplot(gs[1, 0])
ok = dd["n_ev_slip_start"] >= 10
ax.plot(dd["lag_days"][ok], dd["depth_resid_from_slip_start"][ok], color=BLUE, lw=1.5,
        label="depth")
ax.plot(dd["lag_days"][ok], dd["ddist_resid_from_slip_start"][ok], color=ORANGE, ls="--",
        lw=1.0, alpha=0.9, label="down-dip distance")
ax.axhline(0, color=DARK, lw=0.8, ls=":")
ax.axvline(0, color=GREY, lw=0.8, ls="--")
ax.set_xlabel("days from slip start")
ax.set_ylabel("residual (km)")
ax.text(0.03, 0.03, "$\\geq$ 10 events per lag", transform=ax.transAxes, fontsize=6.5,
        color=DARK)
ax.legend(loc="upper left", fontsize=7)
panel(ax, "d", "stacked migration")

# (e) radius / anchor control
ax = fig.add_subplot(gs[1, 1])
for i, (anch, col, mk) in enumerate([("tremor", BLUE, "o"), ("catalog", GREEN, "s")]):
    g = sn[sn.anchor == anch]
    ax.errorbar(g["radius_km"] + (i - 0.5) * 3.0, g["dz_depth_km"], yerr=g["null_sd"],
                fmt=mk + "-", ms=5, color=col, capsize=2.5, lw=1.2, label=f"anchor: {anch}")
    for _, rr in g.iterrows():
        ax.annotate(f"{rr['p_null']:.2f}", (rr["radius_km"] + (i - 0.5) * 3.0,
                                            rr["dz_depth_km"]),
                    textcoords="offset points", xytext=(0, 8 if i == 0 else -14),
                    fontsize=6.2, color=col, ha="center")
ax.axhline(0, color=DARK, lw=0.8, ls=":")
ax.set_xlabel("tremor-selection radius (km)")
ax.set_ylabel("tail $-$ main depth (km)")
ax.text(0.03, 0.03, "numbers = null p", transform=ax.transAxes, fontsize=6.5, color=DARK)
ax.legend(loc="upper right", fontsize=7)
panel(ax, "e", "radius / anchor control")

# (f) five-bin version
ax = fig.add_subplot(gs[1, 2])
g = dms[dms.profile == "tremor"]
xc = 0.5 * (g["phase_lo"] + g["phase_hi"])
ax.fill_between(xc, g["null_lo90"], g["null_hi90"], color=GREY, alpha=0.3, label="null 90 %")
ax.errorbar(xc, g["obs_median"], yerr=[g["obs_median"] - g["boot_lo90"],
                                       g["boot_hi90"] - g["obs_median"]],
            fmt="s-", ms=4, color=VERM, capsize=2, lw=1.2, label="observed")
ax.axhline(0, color=DARK, lw=0.8, ls=":")
ax.set_xlabel("phase within event")
ax.set_ylabel("tremor depth anomaly (km)")
ax.legend(loc="upper left", fontsize=7)
panel(ax, "f", "S9 five-bin version")

fig.savefig(FIG / "s16_paper_fig3.png")
plt.close(fig)
print("wrote figures/s16_paper_fig3.png")


# =====================================================================================
# FIGURE 4 — model checks and central Oregon (S3, S6, S7, S12, S12b)
# =====================================================================================
rep = rd("s12_representability_triangles.csv")
mis = rd("s12_misfit_by_latitude.csv")
occ = rd("s7_occupancy_summary.csv")
dvm = rd("s12b_data_vs_model.csv")
ie = rd("s6_interets_conf.csv")
pro = rd("s3_profiles.csv")

fig = plt.figure(figsize=(7.4, 7.4), layout="constrained")
gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 0.85])

# (a) representability map
ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(rep["lon"], rep["lat"], c=rep["r30"], s=12, marker="s", cmap="cividis",
                vmin=0, vmax=np.nanpercentile(rep["r30"], 98), linewidths=0)
cl = rep[rep["s7_cluster"].notna()]
ax.scatter(cl["lon"], cl["lat"], s=24, marker="s", facecolor="none", edgecolor=VERM,
           linewidths=0.9, label="central cluster")
ax.set_aspect(LAT_ASPECT)
ax.set_xlabel("longitude (deg E)")
ax.set_ylabel("latitude (deg N)")
ax.set_xticks([-127, -125, -123])
ax.grid(alpha=0.15)
cb = fig.colorbar(sc, ax=ax, location="bottom", shrink=0.92, pad=0.02, aspect=26)
cb.set_label("$r_{30}$: 30-km patch recovered by the 13 patterns", fontsize=6.8)
cb.ax.tick_params(labelsize=7)
ax.legend(loc="lower left", fontsize=6.5)
panel(ax, "a", "pattern representability")

# (b) coverage vs latitude
ax = fig.add_subplot(gs[0, 1])
ax.plot(mis["median_r30_deep_tri"], mis["lat_band"] + 0.5, "o-", ms=4, color=BLUE,
        label="$r_{30}$ (deep tri.)")
ax.plot(mis["median_frac_misfit"], mis["lat_band"] + 0.5, "s--", ms=4, color=GREEN,
        label="GNSS frac. misfit")
ax2 = ax.twiny()
ax2.plot(mis["median_nearest_station_km_deep_tri"], mis["lat_band"] + 0.5, "^:", ms=4,
         color=ORANGE, label="nearest station")
ax2.set_xlabel("nearest GNSS station (km)", color=ORANGE, fontsize=7.5)
ax2.tick_params(axis="x", colors=ORANGE, labelsize=7)
ax.axhspan(44.3, 45.4, color=VERM, alpha=0.12)
ax.set_xlabel("fraction")
ax.set_ylabel("latitude (deg N)")
ax.set_ylim(39.5, 50)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=6.3)
panel(ax, "b", "not coverage-starved")

# (c) rank correlation per segment
ax = fig.add_subplot(gs[0, 2])
rk = occ[occ.kind == "rank_diagnostic"].set_index("seg")
ep = ["spearman_S_D_E1", "spearman_S_D_E2", "spearman_S_D_E3"]
w = 0.26
for i, seg in enumerate(["north", "central", "south"]):
    ax.bar(np.arange(3) + (i - 1) * w, [float(rk.loc[seg, e]) for e in ep], width=w,
           color=SEGC[seg], label=seg)
ax.set_xticks(range(3))
ax.set_xticklabels(["2010-14", "2015-19", "2020-26"], fontsize=7.5)
ax.set_ylabel("Spearman(slip, tremor-days)")
ax.set_ylim(0, 0.95)
ax.legend(loc="upper right", fontsize=7)
panel(ax, "c", "weakest co-location: centre")

# (d) data vs model
ax = fig.add_subplot(gs[1, 0])
g = dvm[dvm.estimator == "JOINT"].set_index("target")
targets = ["CLUSTER", "NORTH_ST"]
w = 0.26
for i, (key, col, lab) in enumerate([("raw_GNSS", GREY, "raw GNSS"),
                                     ("ICA_13comp", SKY, "13-comp. ICA"),
                                     ("MODEL_slip", GREEN, "model slip")]):
    ax.bar(np.arange(2) + (i - 1) * w, [float(g.loc[t, key]) for t in targets], width=w,
           color=col, label=lab)
for j, t in enumerate(targets):
    ax.annotate(f"model/data\n{float(g.loc[t, 'model_over_data']):.2f}", (j, 1.0),
                textcoords="offset points", xytext=(0, 3), ha="center", fontsize=6.5,
                color=DARK)
ax.set_xticks(range(2))
ax.set_xticklabels(["central\ncluster", "north\ncontrol"], fontsize=7.5)
ax.set_ylabel("transient amplitude (mm)")
ax.set_ylim(0, 33)
ax.legend(loc="upper left", fontsize=6.5)
panel(ax, "d", "model vs data amplitude")

# (e) inter-ETS confirmation
ax = fig.add_subplot(gs[1, 1])
sub = ie[(ie.group == "north/central") & (ie.burst_kind == "slip_burst")].reset_index(drop=True)
xs = np.arange(len(sub))
ax.bar(xs - 0.18, sub["obs_conf_rate"].fillna(0), width=0.34, color=BLUE, label="observed")
ax.bar(xs + 0.18, sub["null_mean"].fillna(0), width=0.34, color=GREY, label="chance")
for i, rr in sub.iterrows():
    ax.annotate(f"n={int(rr['n_bursts'])}", (i, 0.02), ha="center", fontsize=6.8, color=DARK)
ax.set_xticks(xs)
ax.set_xticklabels(sub["definition"], fontsize=8)
ax.set_xlabel("inter-ETS definition")
ax.set_ylabel("tremor-confirmed fraction")
ax.set_ylim(0, 1.0)
ax.text(0.5, 0.72, "A0 leaks real ETS days;\nA1/A2 (strict) leave 0-1\nbursts north of 44N",
        transform=ax.transAxes, fontsize=6.5, color=DARK, ha="center")
ax.legend(loc="upper right", fontsize=6.8)
panel(ax, "e", "inter-ETS slip: unconfirmed")

# (f) slip per tremor along strike
ax = fig.add_subplot(gs[1, 2])
for lab, key, col in [("2010-14", "E1", DARK), ("2015-19", "E2", SKY), ("2020-26", "E3", GREEN)]:
    g = pro[pro.epoch == key].dropna(subset=["median_eff_norm"])
    ax.plot(g["median_eff_norm"], g["lat_bin_centre"], "o-", ms=3, color=col, label=lab)
ax.set_xscale("log")
ax.axhspan(44.3, 45.4, color=VERM, alpha=0.12)
ax.set_xlabel("slip per tremor day (norm.)")
ax.set_ylabel("latitude (deg N)")
ax.set_ylim(39.5, 50)
ax.legend(loc="lower right", fontsize=6.8)
panel(ax, "f", "same pattern every epoch")

fig.savefig(FIG / "s16_paper_fig4.png")
plt.close(fig)
print("wrote figures/s16_paper_fig4.png")


# =====================================================================================
# TABLE 1 — the 27 confirmed events
# =====================================================================================
chk = rd("s15_events_check.csv")[["id", "PASS", "overlap", "jaccard", "peak_mm_per30d"]]
t1 = ev[["id", "start", "end", "seg", "Mw", "dur_d", "n_tri_F", "potency_km3", "N",
         "Q1_downdip_km", "Q3_depth_km", "log10Q5", "z1", "z3"]].copy()
t1 = t1.merge(chk, on="id", how="left")
t1 = t1.merge(s9e[["id", "n_days_fit", "exponent", "null_exp_p"]].rename(
    columns={"exponent": "b_count", "null_exp_p": "b_count_null_p"}), on="id", how="left")
t1 = t1.merge(s11e[["id", "b_E_R7"]].rename(columns={"b_E_R7": "b_energy_sum"}), on="id",
              how="left")
t1 = t1.merge(s16e[["id", "b_E_med", "b_A1_sqrt_Emed", "beta_A1_sqrt_Emed"]].rename(
    columns={"b_E_med": "b_energy_median", "b_A1_sqrt_Emed": "b_amplitude_FB",
             "beta_A1_sqrt_Emed": "beta_FB_direction"}), on="id", how="left")
t1 = t1.merge(fr[["id", "L_days", "v_t_km_d", "v_s_km_d", "dir_agree", "corr_at_L"]], on="id",
              how="left")
t1 = t1.sort_values("start").reset_index(drop=True)
t1.insert(0, "no", np.arange(1, len(t1) + 1))
# round to paper precision (full precision stays in the source CSVs of S1/S9/S11/S15/S16)
rnd = {"Mw": 2, "potency_km3": 3, "Q1_downdip_km": 1, "Q3_depth_km": 1, "log10Q5": 2, "z1": 2,
       "z3": 2, "overlap": 3, "jaccard": 3, "peak_mm_per30d": 1, "b_count": 3,
       "b_count_null_p": 3, "b_energy_sum": 3, "b_energy_median": 3, "b_amplitude_FB": 3,
       "beta_FB_direction": 2, "v_t_km_d": 2, "v_s_km_d": 2, "corr_at_L": 2}
t1 = t1.round(rnd)
t1.to_csv(RES / "s16_table1.csv", index=False)
print(f"wrote results/s16_table1.csv ({len(t1)} rows, {t1.shape[1]} cols)")


# =====================================================================================
# TABLE 2 — every exponent with CI, null and provenance
# =====================================================================================
rows = []


def add(block, measure, aggregation, regressor, src_step, src_csv, r, note="",
        cols=None):
    cols = cols or {}
    g = lambda k, d=np.nan: float(r[cols.get(k, k)]) if cols.get(k, k) in r and pd.notna(
        r[cols.get(k, k)]) else d
    rows.append(dict(
        block=block, measure=measure, aggregation=aggregation, regressor=regressor,
        n_events=g("n_events"), b_median=g("b_median"), b_q25=g("b_q25"), b_q75=g("b_q75"),
        b_ci90_lo=g("b_boot_lo90"), b_ci90_hi=g("b_boot_hi90"), n_below_1=g("n_below_1"),
        pooled_ols=g("pooled_ols"), pooled_theilsen=g("pooled_theilsen"),
        null_mean=g("null_mean"), null_sd=g("null_sd"), z_median=g("z_median"),
        n_z_positive=g("n_z_positive"), n_events_z=g("n_events_z"), sign_p=g("sign_p"),
        n_events_null_p05=g("n_events_p05"), source_step=src_step, source_csv=src_csv,
        note=note))


# --- S9/S11 model-based
for meas, lab, note in [("C", "tremor count", "R9 headline; 27 events"),
                        ("U", "unique tremor locations", "network-growth control"),
                        ("CperU", "count per unique location", ""),
                        ("C_ns8", "count, >= 8 stations", "detection-threshold control"),
                        ("E", "summed tremor energy", "2017+; R9 headline"),
                        ("EperDet", "energy per detection", "= E_mean in S16"),
                        ("Dur", "summed duration", ""),
                        ("Mmax", "daily max magnitude", "")]:
    r = s11x[(s11x.measure == meas) & (s11x.slip_var == "R7") & (s11x.subset == "all")]
    if len(r):
        add("model-based (S11)", lab, "daily sum over the 60-km neighbourhood", "R7 (vbICA)",
            "S11", "s11_exponents.csv", r.iloc[0], note)
for sv, lab in [("R1", "R1 (unsmoothed rate)"), ("R3", "R3"), ("Rcum", "cumulative slip")]:
    r = s11x[(s11x.measure == "C") & (s11x.slip_var == sv) & (s11x.subset == "all")]
    if len(r):
        add("bandwidth ladder (S11)", "tremor count", "daily", lab, "S11",
            "s11_exponents.csv", r.iloc[0], "regressor bandwidth, errors-in-variables ladder")

# --- S11b model-free
for meas, lab in [("C", "tremor count"), ("E", "summed tremor energy")]:
    for p, plab, note in [("R7_vbICA", "R7 (vbICA, reference)", "reproduces S11 exactly"),
                          ("Gica31_proj", "ICA displacement, 31 d", "matched bandwidth"),
                          ("Graw31_proj", "raw GNSS displacement, 31 d",
                           "MODEL-FREE, matched bandwidth"),
                          ("Graw_proj", "raw GNSS displacement, 1 d", "MODEL-FREE, unsmoothed")]:
        r = s11b[(s11b.measure == meas) & (s11b.proxy == p) & (s11b.subset == "all")]
        if len(r):
            add("model-free (S11b)", lab, "daily", plab, "S11b", "s11b_exponents.csv",
                r.iloc[0], note)

# --- S15B candidate measures
s15s = s15c[s15c.kind == "summary"]
for meas, lab, note in [
        ("C_cand", "picker-ranked candidate count", "NOT independent of PNSN (daily rho 0.97)"),
        ("W_cand", "candidate count, p_lfe-weighted", "rank-weighted, not probability-weighted"),
        ("C_q75", "candidates, top-quartile p_lfe", ""),
        ("C_s10", "candidates, 10 nearest stations", "network-growth control"),
        ("M_cand", "multiplicity C_cand / C_pnsn", "the genuinely new content of S15B"),
        ("C_pnsn", "PNSN count (reference)", "same as S11 C")]:
    r = s15s[s15s.measure == meas]
    if len(r):
        add("station-side (S15B)", lab, "daily", "R7 (vbICA)", "S15B",
            "s15_cand_exponents.csv", r.iloc[0], note,
            cols={"null_mean": "null_pooled_mean", "null_sd": "null_pooled_sd",
                  "z_median": "null_pooled_z_median", "n_z_positive": "null_pooled_n_positive",
                  "n_events_z": "null_pooled_n_events", "sign_p": "null_pooled_sign_p",
                  "n_events_p05": "null_pooled_n_p05"})

# --- S16 amplitude aggregations
for meas, lab, note in [
        ("C", "tremor count (18 events)", "cross-check: reproduces S11 per event"),
        ("E_sum", "summed energy, eps = 0", ""),
        ("E_sum_epsS11", "summed energy, S11 eps", "reproduces S11 b_E exactly"),
        ("E_med", "MEDIAN energy per detection", "the Frank & Brodsky aggregation"),
        ("E_mean", "mean energy per detection", "= S11 EperDet"),
        ("E_p10", "10th percentile energy", ""),
        ("E_p90", "90th percentile energy", ""),
        ("E_med_ns8", "median energy, >= 8 stations", "threshold control"),
        ("M_med", "median magnitude", "units: magnitude per decade of R"),
        ("A1_sqrt_Emed", "F&B amplitude A1 = sqrt(median E)", "PRIMARY like-for-like proxy"),
        ("A2_mag_0p5", "A2 = 10^(0.5 M_med)", "0.561x A1 BY CONSTRUCTION"),
        ("A3_mag_EMfit", "A3 = 10^(0.891 M_med)", "A1 via the fitted E-M slope"),
        ("A4_mag_1p0", "A4 = 10^(1.0 M_med)", ""),
        ("Asum_sqrt_Esum", "sqrt(daily summed energy)", "= b(E_sum)/2 identically, not an amplitude"),
        ("Asum_lin", "daily summed sqrt(energy)", "literal summed amplitude = count x amplitude")]:
    add("F&B aggregation (S16)", lab, "daily, per-detection", "R7 (vbICA)", "S16",
        "s16_fb_aggregation.csv", s16row(meas), note)
for lv, lab in [("bin1d_theilsen", "1-d bins"), ("bin3d_theilsen", "3-d bins"),
                ("bin7d_theilsen", "7-d bins")]:
    add("F&B bin ladder (S16)", "F&B amplitude A1", lab, "R7 (vbICA)", "S16",
        "s16_fb_aggregation.csv", s16row("A1_sqrt_Emed", "forward_measure_on_rate", lv),
        "forward slope rises as the bin lengthens (attenuation)")
add("F&B bin ladder (S16)", "F&B amplitude A1", "one point per event (n = 18)", "R7 (vbICA)",
    "S16", "s16_fb_aggregation.csv",
    s16row("A1_sqrt_Emed", "cross_event_amp_on_rate", "n18_theilsen"),
    "closest to F&B's 8-SSE granularity; no bootstrap CI at this level")
add("F&B direction (S16)", "beta: rate on A1 (F&B's direction)", "daily", "A1", "S16",
    "s16_fb_aggregation.csv", s16row("A1_sqrt_Emed", "reverse_rate_on_amplitude_FB",
                                     "per_event_theilsen"),
    "compare F&B beta = 3.17 +- 0.29; all 18 events below 2.88")
add("F&B direction (S16)", "RMA slope, A1 on rate", "daily, direction-symmetric", "R7 (vbICA)",
    "S16", "s16_fb_aggregation.csv", s16row("A1_sqrt_Emed", "rma_amp_on_rate", "per_event_rma"),
    "beta_RMA = 2.54 brackets F&B 3.17")

# --- group medians over saved per-event exponents (aggregated here)
def grp(name, sub, col, block, note):
    v = np.asarray(sub[col], float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return
    rows.append(dict(block=block, measure=name, aggregation="median of saved per-event exponents",
                     regressor="R7 (vbICA)", n_events=v.size, b_median=float(np.median(v)),
                     b_q25=float(np.percentile(v, 25)), b_q75=float(np.percentile(v, 75)),
                     b_ci90_lo=np.nan, b_ci90_hi=np.nan, n_below_1=int((v < 1).sum()),
                     pooled_ols=np.nan, pooled_theilsen=np.nan, null_mean=np.nan, null_sd=np.nan,
                     z_median=np.nan, n_z_positive=np.nan, n_events_z=np.nan, sign_p=np.nan,
                     n_events_null_p05=np.nan, source_step="S9/S11/S16 (aggregated here)",
                     source_csv="s9_rate_scaling.csv / s11_energy_events.csv / s16_fb_events.csv",
                     note=note))


for seg in ["north", "central"]:
    grp(f"tremor count, {seg}", s9e[s9e.seg == seg], "exponent", "stationarity (S9)",
        "segment split of the 27-event count exponent")
    grp(f"F&B amplitude A1, {seg}", s16e[s16e.seg == seg], "b_A1_sqrt_Emed",
        "stationarity (S16)", "segment split of the 18-event amplitude exponent")
for epo in ["2010-2014", "2015-2019", "2020-2026"]:
    grp(f"tremor count, {epo}", s9e[s9e.epoch == epo], "exponent", "stationarity (S9)",
        "epoch split; Kruskal p 0.94 (NOTES S11)")
for epo in ["2015-2019", "2020-2026"]:
    grp(f"summed energy, {epo}", s11e[s11e.epoch == epo], "b_E_R7", "stationarity (S11)",
        "epoch split of the energy exponent (2017+ only)")
    grp(f"F&B amplitude A1, {epo}", s16e[s16e.epoch == epo], "b_A1_sqrt_Emed",
        "stationarity (S16)", "epoch split of the amplitude exponent")

T2 = pd.DataFrame(rows)
num = [c for c in T2.columns if T2[c].dtype.kind == "f"]
T2[num] = T2[num].round(4)
T2.to_csv(RES / "s16_table2.csv", index=False)
print(f"wrote results/s16_table2.csv ({len(T2)} rows, {T2.shape[1]} cols)")
