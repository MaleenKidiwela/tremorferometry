#!/usr/bin/env python
"""bw_firstlook_plot -- figure + tables for the backward first look (reads first_look/daily_counts.csv).
Context only: Bostock ETS window DATES (truth_episodes.csv); no event-level matching. Confidence model P not applied."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
OUT = f"{BW}/first_look"
D = pd.read_csv(f"{OUT}/daily_counts.csv"); D["date"] = pd.to_datetime(D.date)
E = pd.read_csv(f"{BW}/../backward_scope/truth_episodes.csv")
ets = np.zeros(len(D), bool); dist = np.full(len(D), 10 ** 6)
for r in E.itertuples():
    a, b = pd.Timestamp(r.start), pd.Timestamp(r.end)
    ets |= ((D.date >= a) & (D.date <= b)).to_numpy()
    dd = np.where(D.date < a, (a - D.date).dt.days, np.where(D.date > b, (D.date - b).dt.days, 0))
    dist = np.minimum(dist, dd)
D["ets"] = ets; D["far30"] = dist > 30; D["yr"] = D.date.dt.year
BOUND = [("2003-08-01", "POLARIS/CN array on (era)"), ("2005-07-01", "B1"), ("2005-11-02", "B2"),
         ("2006-07-01", "POLARIS off"), ("2007-11-13", "B3"), ("2009-08-27", "S15")]

# detector variants: (label, mask, data col, control col)
fz = D.detector == "frozen"; er = D.detector == "era"
VAR = [("frozen @ T_core(set)", fz, "n_data_Tcore", "n_c1_Tcore"),
       ("era @ 8.75 (provisional)", er, "n_data_8.75", "n_c1_8.75"),
       ("era @ 9.0 (provisional)", er, "n_data_9", "n_c1_9")]

# ---- table 1: ETS days vs days > 30 d from any window, per year x detector
rows = []
for lab, m, cd, cc in VAR:
    X = D[m]
    for y, g in X.groupby("yr"):
        for grp, k in (("ETS window", g.ets), ("> 30 d from window", g.far30)):
            h = g[k]
            if not len(h):
                continue
            ex = h[cd] - h[cc]
            rows.append(dict(detector=lab, year=y, days=grp, n_days=len(h), median_excess=ex.median(),
                             median_data=h[cd].median(), median_control=h[cc].median(),
                             sum_data=int(h[cd].sum()), sum_control=int(h[cc].sum()),
                             ratio_data_over_control=h[cd].sum() / max(h[cc].sum(), 1),
                             data_per_famday=h[cd].sum() / max(h.fd_data.sum(), 1),
                             control_per_famday=h[cc].sum() / max(h.fd_c1.sum(), 1)))
T1 = pd.DataFrame(rows); T1.to_csv(f"{OUT}/table_ets_vs_between.csv", index=False)

# ---- table 2: steps at station-set changes (frozen; between-episode days, 60 d each side)
rows = []
for b, lab in BOUND:
    t = pd.Timestamp(b)
    for detlab, m, cd, cc in VAR[:2]:
        X = D[m & D.far30]
        for side, k in (("before", (X.date < t) & (X.date >= t - pd.Timedelta(days=60))),
                        ("after", (X.date >= t) & (X.date < t + pd.Timedelta(days=60)))):
            h = X[k]
            if len(h) < 10:
                continue
            rows.append(dict(boundary=b, change=lab, detector=detlab, side=side, n_days=len(h),
                             T=h["T_core"].iloc[0] if detlab.startswith("frozen") else 8.75,
                             control_per_famday=h[cc].sum() / max(h.fd_c1.sum(), 1),
                             data_per_famday=h[cd].sum() / max(h.fd_data.sum(), 1),
                             median_excess=(h[cd] - h[cc]).median(), obs_fam_per_day=h.fd_data.median()))
T2 = pd.DataFrame(rows); T2.to_csv(f"{OUT}/table_steps.csv", index=False)

# ---- figure
C = {"frozen": "#2a78d6", "era": "#eb6834"}
fig, ax = plt.subplots(3, 1, figsize=(14, 10.5), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 2.2, 2]))
for a in ax:
    for r in E.itertuples():
        a.axvspan(pd.Timestamp(r.start), pd.Timestamp(r.end) + pd.Timedelta(days=1), color="0.85", lw=0, zorder=0)
    for bd, lab in BOUND:
        a.axvline(pd.Timestamp(bd), color="0.35", ls=":", lw=1, zorder=1)
    a.grid(axis="y", color="0.92", lw=0.6); a.set_axisbelow(True)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
for i, lab, m, cd, cc, col, ls, lw in [(0, *VAR[0], C["frozen"], "-", 0.7), (1, *VAR[1], C["era"], "-", 0.7),
                                         (1, *VAR[2], "#8a3410", "-", 0.5)]:
    X = D[m].sort_values("date").set_index("date")
    sr = (X[cd] - X[cc]).asfreq("D")   # gaps stay gaps
    ax[i].plot(sr.index, sr.values, color=col, ls=ls, lw=lw, label=lab)
for i, t in ((0, "(a) frozen detector (S15 sites), T_core per station set: CN4 10.75, B1 10.5, B2/B3 9.0, S15 8.75"),
             (1, "(b) era detector (S15ERA, with KELB; POLARIS era only, runs to 2006-06): T = 8.75 and 9.0, PROVISIONAL (no P4-3 value)")):
    ax[i].set_yscale("symlog", linthresh=5); ax[i].axhline(0, color="0.2", lw=0.8)
    ax[i].set_ylabel("data minus scrambled\ncontrol, per day")
    ax[i].set_title(t, fontsize=10, loc="left")
ax[1].legend(loc="upper right", frameon=False, fontsize=9)
for bd, lab in BOUND:
    ax[0].text(pd.Timestamp(bd) + pd.Timedelta(days=10), -18, lab, rotation=90, va="bottom", fontsize=7.5, color="0.3")
fig.suptitle("SVI backward extension, first look (provisional): candidates above threshold, data minus scrambled control.\n"
             "Confidence model P not yet applied. Grey = Bostock ETS window dates; dotted = station-set changes.",
             fontsize=10.5, x=0.01, ha="left")
for (lab, m, cd, cc), col in zip(VAR[:2], [C["frozen"], C["era"]]):
    X = D[m & D.far30].sort_values("date").set_index("date").asfreq("D")
    num_d = X[cd].rolling(30, min_periods=10).sum(); num_c = X[cc].rolling(30, min_periods=10).sum()
    den_d = X.fd_data.rolling(30, min_periods=10).sum(); den_c = X.fd_c1.rolling(30, min_periods=10).sum()
    nm = lab.split(" @")[0]
    ax[2].plot(X.index, num_d / den_d, color=col, lw=1.4, label=f"{nm} data")
    ax[2].plot(X.index, num_c / den_c, color=col, lw=1.0, ls=(0, (2, 2)), label=f"{nm} scrambled control")
ax[2].set_yscale("log"); ax[2].set_ylabel("candidates per observable\nfamily-day (30-d mean)")
ax[2].set_title("(c) between episodes only (days > 30 d from any ETS window); era detector at 8.75", fontsize=10, loc="left")
ax[2].legend(loc="lower right", frameon=False, fontsize=9, ncol=2)
ax[2].set_xlabel("date (UTC day; the PT day starts 7-8 h later)")
fig.tight_layout(rect=(0, 0, 1, 0.955)); fig.savefig(f"{OUT}/fig_daily_excess_2003_2013.png", dpi=150)
pd.set_option("display.width", 220)
print(T1.round(3).to_string()); print(T2.round(4).to_string())
