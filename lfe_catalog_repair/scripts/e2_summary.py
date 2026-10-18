#!/usr/bin/env python
"""E2 summary tables (2026-09-10, resumed run) -> results/e2_summary.txt (and stdout).
Every number quoted in NOTES.md / RESULTS.md for E2 comes from this file, so it can be regenerated.
"""
from __future__ import annotations

import glob
import io
import os

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
buf = io.StringIO()


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    buf.write(s + "\n")


# ---- 1. picks
p("=" * 100)
p("1. PICKS ON DISK (e2_picks_*.parquet; the parquet holds the UNION over the threshold sweep,")
p("   each pick tagged with the lowest sweep threshold that fires it -> any threshold is exact)")
rows = []
for f in sorted(glob.glob(f"{RES}/e2_picks_*.parquet")):
    d = pd.read_parquet(f, columns=["set", "model", "sta", "phase", "prob", "thr_first", "day"])
    n = dict(file=os.path.basename(f), set=d["set"].iloc[0], model=d.model.iloc[0],
             stations=d.sta.nunique(), station_days=int(d.groupby(["sta", "day"]).ngroups),
             picks_all=len(d))
    for thr in (0.1, 0.3, 0.5):
        n[f"S_thr{thr}"] = int(((d.phase == "S") & (d.thr_first <= thr) & (d.prob >= thr)).sum())
    n["S_per_sta_day_thr0.3"] = round(n["S_thr0.3"] / max(n["station_days"], 1), 1)
    rows.append(n)
pk = pd.DataFrame(rows)
p(pk.to_string(index=False))
pk.to_csv(f"{RES}/e2_picks_inventory.csv", index=False)

# ---- 2. throughput
p("")
p("=" * 100)
p("2. GPU THROUGHPUT (e2_throughput.csv, NVIDIA L40S, 4 concurrent queues x 2 torch threads)")
tp = pd.read_csv(f"{RES}/e2_throughput.csv")
g = tp.groupby(["set", "model"]).agg(runs=("station", "count"), station_days=("n_days", "sum"),
                                     skipped=("n_skipped", "sum"), picks=("n_picks", "sum"),
                                     seconds=("seconds", "sum"))
g["sta_days_per_min_per_queue"] = (g.station_days / g.seconds * 60).round(1)
p(g.to_string())
tot_d, tot_s = tp.n_days.sum(), tp.seconds.sum()
rate = tot_d / tot_s * 60
p(f"TOTAL {tot_d} station-days in {tot_s/3600:.2f} queue-hours -> {rate:.1f} station-days/min "
  f"per queue (median per run {tp.station_days_per_min.median():.1f}, "
  f"range {tp.station_days_per_min.min():.1f}-{tp.station_days_per_min.max():.1f})")
p(f"  = {365/rate:.2f} GPU-min per station-year per queue per model; "
  f"{4*rate*60/365:.1f} station-years/hour on 4 queues")
if os.path.exists(f"{RES}/f1_cost_estimate.csv"):
    ce = pd.read_csv(f"{RES}/f1_cost_estimate.csv")
    for scope in ce.scope.unique():
        t = ce[(ce.scope == scope) & (ce.year == "TOTAL")]
        if not len(t):
            continue
        sy = float(t.station_years_effective.iloc[0])
        qh = sy * 365 / rate / 60
        p(f"  F1 {scope}: {sy:,.0f} effective station-years -> {qh:,.0f} GPU queue-hours/model "
          f"= {qh/4:,.0f} h wall-clock on 4 queues ({qh/4/24:.1f} days); both weight sets "
          f"= {qh/2:,.0f} h wall-clock")

# ---- 3. PGC single-station validation
p("")
p("=" * 100)
p("3. SINGLE-STATION VALIDATION, CN.PGC 2010-05-01..2010-11-01 (184 days)")
p("   3a. fixed +-3 s window, next to the time-shifted null (e2_validation_pgc2010.csv)")
v = pd.read_csv(f"{RES}/e2_validation_pgc2010.csv")
s = v[(v.phase == "S") & (v.tol_s == 3.0)
      & (v["gt"].isin(["lin_hq", "lin_near40", "bostock_stack", "bostock_near40"]))]
p(s.pivot_table(index=["gt", "threshold"], columns="model",
                values=["picks_per_day", "recall", "recall_null", "precision", "precision_null",
                        "offset_s"]).round(4).to_string())
p("")
p("   3b. chance-free PEAK EXCESS (+-1 s around the measured peak; e2_resid_excess.csv)")
ex = pd.read_csv(f"{RES}/e2_resid_excess.csv")
e1 = ex[(ex.half_width_s == 1.0) & (ex["gt"].isin(["lin_hq", "bostock_near40", "lin_1d"]))]
p(e1.pivot_table(index=["gt", "threshold"], columns="model",
                 values=["peak_lag_s", "n_picks", "excess", "snr", "recall_excess",
                         "precision_excess"]).round(4).to_string())
p("")
o = pd.read_csv(f"{RES}/e2_pgc2010_offsets.csv")
p("   3c. catalogue time-base offsets (peak of the pick-minus-predicted histogram):")
p(o.round(3).to_string(index=False))

# ---- 4. association
p("")
p("=" * 100)
p("4. ASSOCIATION (PyOcto, Crosson 1-D; e2_assoc_scores.csv, tolerance 2 s on ORIGIN times)")
sc = pd.read_csv(f"{RES}/e2_assoc_scores.csv")
sc["min_sta"] = sc.tag.str.extract(r"_n(\d)$").astype(int)
for gk in ("lin_hq", "bostock"):
    t = sc[(sc.tol_s == 2.0) & (sc["gt"] == gk) & (sc["set"] == "vi2010")]
    if not len(t):
        continue
    p(f"   4a. vi2010 vs {gk}")
    p(t.pivot_table(index=["threshold", "min_sta"], columns="model",
                    values=["events_per_day", "recall", "recall_null", "precision",
                            "precision_null", "recall_excess_2s", "precision_excess_2s"]
                    ).round(4).to_string())
    p("")
p("   4b. event rate vs PNSN tremor (from the per-config logs, recomputed here)")
rows = []
for f in sorted(glob.glob(f"{RES}/e2_assoc_daily_*.csv")):
    tag = os.path.basename(f)[len("e2_assoc_daily_"):-4]
    d = pd.read_csv(f, index_col=0, parse_dates=True).fillna(0.0).sort_index()
    roll = d.tremor.rolling(30, min_periods=5, center=True).sum()
    hi, lo = np.nanpercentile(roll, 90), np.nanpercentile(roll, 25)
    from scipy.stats import spearmanr
    rho, pv = spearmanr(d.events, d.tremor)
    ets, qui = d.events[roll >= hi].mean(), d.events[roll <= lo].mean()
    rows.append(dict(tag=tag, n_days=len(d), events_per_day=d.events.mean(),
                     ets=ets, quiet=qui, ets_quiet=ets / max(qui, 1e-9), rho=rho, p=pv))
r = pd.DataFrame(rows)
r.to_csv(f"{RES}/e2_assoc_rates.csv", index=False)
p(r.round(3).to_string(index=False))

# ---- 5. locations
p("")
p("=" * 100)
p("5. LOCATIONS (NonLinLoc, Crosson 1-D, S 0.4 s / P 1.5 s; e2_loc_compare_*.csv)")
for f in sorted(glob.glob(f"{RES}/e2_loc_compare_*.csv")):
    lc = pd.read_csv(f)
    tag = os.path.basename(f)[len("e2_loc_compare_"):-4]
    p(f"   {tag}: {len(lc)} located events")
    for c, lab in (("dh_nearest_family", "epicentre to nearest Bostock family"),
                   ("dh_lin", "to time-matched Lin (<=2 s)"),
                   ("dh_bostock", "to time-matched Bostock (<=2 s)")):
        vv = lc[c].dropna()
        if len(vv) >= 5:
            p(f"     {lab:38s} n={len(vv):5d} median {np.median(vv):6.1f} km "
              f"IQR {np.percentile(vv,25):5.1f}-{np.percentile(vv,75):6.1f}")
    for c, lab in (("dz_lin", "depth - Lin depth"), ("dz_bostock", "depth - Bostock depth")):
        vv = lc[c].dropna()
        if len(vv) >= 5:
            p(f"     {lab:38s} n={len(vv):5d} median {np.median(vv):+6.1f} km "
              f"IQR {np.percentile(vv,25):+5.1f}..{np.percentile(vv,75):+6.1f}")
    p(f"     depth median {lc.depth.median():.1f} km, frac 25-45 km "
      f"{float(lc.depth.between(25,45).mean()):.3f}, frac <5 km {float((lc.depth<5).mean()):.3f}")
    p(f"     formal errors: ERH median {lc.err_h_km.median():.1f} km, ERZ median "
      f"{lc.err_z_km.median():.1f} km, rms median {lc.rms.median():.2f} s")
    if "n_stations" in lc:
        p("     by station count: " + "; ".join(
            f"N={int(n)} (n={len(gg)}) dh_fam {np.nanmedian(gg.dh_nearest_family):.1f} km, "
            f"ERH {gg.err_h_km.median():.0f} km"
            for n, gg in lc.groupby("n_stations") if len(gg) >= 20))

# ---- 6. off ground truth
p("")
p("=" * 100)
p("6. OFF-GROUND-TRUTH STATIONS 2017-2018 (e2_offgt_summary.csv, e2_offgt_degen.csv)")
if os.path.exists(f"{RES}/e2_offgt_summary.csv"):
    su = pd.read_csv(f"{RES}/e2_offgt_summary.csv")
    p(su.round(3).to_string(index=False))
if os.path.exists(f"{RES}/e2_offgt_coincidence.csv"):
    p("")
    p("   WISH x B943 same-source coincidence (S-pick cross-correlation):")
    p(pd.read_csv(f"{RES}/e2_offgt_coincidence.csv").round(3).to_string(index=False))
if os.path.exists(f"{RES}/e2_offgt_degen.csv"):
    p("")
    p("   degenerate-window QC (sampled days; e2_offgt_degen.csv):")
    p(pd.read_csv(f"{RES}/e2_offgt_degen.csv").round(4).to_string(index=False))

# ---- 7. CPU/GPU profile and concurrency scaling
p("")
p("=" * 100)
p("7. IS IT GPU-BOUND?  (e2_profile.csv, e2_scaling.csv)")
if os.path.exists(f"{RES}/e2_profile.csv"):
    pr = pd.read_csv(f"{RES}/e2_profile.csv")
    med = pr.median(numeric_only=True)
    p("   per-station-day wall time of ONE queue, last profiled station:")
    for k in ("read", "pre", "gpu_and_window", "post", "total"):
        p(f"     {k:16s} {med[k]:7.2f} s  {100*med[k]/med['total']:5.1f} %")
    if "annotate_total_cpu" in pr:
        p(f"     annotate on CPU-only: {pr.annotate_total_cpu.median():.2f} s vs "
          f"{med['annotate_total']:.2f} s on the L40S")
if os.path.exists(f"{RES}/e2_scaling.csv"):
    sc = pd.read_csv(f"{RES}/e2_scaling.csv")
    p("")
    p(sc.to_string(index=False))
    for sta, g in sc.groupby("station"):
        g = g.sort_values("workers")
        p(f"   {sta}: per-queue rate {g.per_queue_sta_days_per_min.iloc[0]:.1f} -> "
          f"{g.per_queue_sta_days_per_min.iloc[-1]:.1f} station-days/min going "
          f"{int(g.workers.iloc[0])} -> {int(g.workers.iloc[-1])} queues; steady-state aggregate "
          f"{g.steady_state_aggregate_per_min.iloc[0]:.0f} -> "
          f"{g.steady_state_aggregate_per_min.iloc[-1]:.0f}; GPU util mean <= "
          f"{g.gpu_util_mean.max():.0f} % (max {g.gpu_util_max.max():.0f} %), "
          f"GPU memory {g.gpu_mem_max_MiB.max():.0f} MiB")
    if os.path.exists(f"{RES}/f1_cost_estimate.csv"):
        ce = pd.read_csv(f"{RES}/f1_cost_estimate.csv")
        b = sc[(sc.station == "PB.B943")]
        if len(b):
            for scope in ce.scope.unique():
                t = ce[(ce.scope == scope) & (ce.year == "TOTAL")]
                if not len(t):
                    continue
                sd = float(t.station_years_effective.iloc[0]) * 365
                p(f"   fleet {scope} ({sd/365:,.0f} station-years = {sd:,.0f} station-days), "
                  "100-Hz-station rates, per weight set:")
                for _, r in b.sort_values("workers").iterrows():
                    p(f"     {int(r.workers):2d} queues @ "
                      f"{r.steady_state_aggregate_per_min:6.1f} station-days/min -> "
                      f"{sd/r.steady_state_aggregate_per_min/60:6.1f} h wall-clock "
                      f"({sd/r.steady_state_aggregate_per_min/60/24:4.1f} d)")

with open(f"{RES}/e2_summary.txt", "w") as f:
    f.write(buf.getvalue())
print(f"\n[write] {RES}/e2_summary.txt")
