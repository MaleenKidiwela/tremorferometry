#!/usr/bin/env python
"""D3 stage 3 — validation of the re-detected catalogue, NEW (MAD) rule vs OLD (fixed cc>=0.80) rule.

Consumes the stage-2 outputs (`<prefix>_detections.parquet`, `<prefix>_daily.parquet`) and produces
  results/d3_detections_<STA>.parquet   final catalogue (one row per detection, per-k flags)
  results/d3_thresholds_<STA>.csv       per template x 30-day window thresholds (k = 8/9/10)
  results/d3_validation_<STA>.csv       per template: ETS/quiet, Spearman vs tremor, counts, old vs new
  results/d3_station_summary.csv        one row per station x rule
Network coincidence and the Lin ground-truth test are separate sub-commands.

ETS/quiet definition is D1's, verbatim: PNSN tremor within 40 km of the template's SOURCE CELL,
trailing 30-day count; ETS day = count > 90th pct of that cell's series, quiet day = < 25th pct.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PNSN = "catalogs/pnsn_tremor_cascadia_full.csv"
KS = (8.0, 9.0, 10.0)
OLD_THR = 0.80
OLD_GAP_S = 6.0
FS = 40.0


def cell_of(fam):
    a = fam.split("__")[0].split("_")
    return float(a[0]), float(a[1])


def load_tremor():
    t = pd.read_csv(PNSN, usecols=["time", "lat", "lon"], parse_dates=["time"])
    t["day"] = t.time.dt.floor("D")
    return t


def cell_tremor_daily(tremor, lat, lon, radius_km=40.0):
    d = np.hypot((tremor.lat.values - lat) * 111.2,
                 (tremor.lon.values - lon) * 111.2 * np.cos(np.radians(lat)))
    s = tremor.loc[d <= radius_km].groupby("day").size()
    return s


def ets_quiet_mask(daily_counts, all_days):
    """D1 rule: trailing 30-d tremor count; ETS = >p90, quiet = <p25."""
    c = daily_counts.reindex(all_days, fill_value=0).astype(float)
    roll = c.rolling(30, min_periods=1).sum()
    p90, p25 = np.percentile(roll.values, 90), np.percentile(roll.values, 25)
    return (roll.values > p90), (roll.values < p25), roll.values


def windows_30d(days, d0=None):
    d0 = days.min() if d0 is None else d0
    return ((pd.to_datetime(days) - d0).days // 30 if hasattr(days, "days")
            else ((pd.DatetimeIndex(days) - d0).days // 30)).astype(int)


def build_thresholds(daily, ks=KS):
    """Per template x 30-day calendar window: med/mad from the daily robust stats -> thr per k."""
    d = daily.copy()
    d["day"] = pd.to_datetime(d["day"])
    d = d[d.qc_drop == 0].copy()
    d["template"] = d["template"].astype(str)
    d0 = d["day"].min()
    d["win"] = windows_30d(d["day"].values, d0)
    g = d.groupby(["template", "win"], observed=True).agg(
        med=("med", "median"), mad=("mad", "median"), n_days=("day", "nunique"),
        n_comp=("n_comp", "median"), day0=("day", "min"), day1=("day", "max"))
    for k in ks:
        g[f"thr{int(k)}"] = g["med"] + k * g["mad"]
        g[f"attainable{int(k)}"] = (g[f"thr{int(k)}"] <= 1.0).astype(int)
    g = g.reset_index()
    g["template"] = g["template"].astype(str)
    g.attrs["day0"] = str(d0)
    return g


def apply_rules(det, thr, day0, ks=KS):
    """Flag each detection with the new rule at each k, and with the old fixed rule."""
    det = det.copy()
    det["dt"] = pd.to_datetime(det["t"], unit="s")
    det["day"] = det["dt"].dt.floor("D")
    det["win"] = ((det["day"] - pd.Timestamp(day0)).dt.days // 30).astype(int)
    key = thr.set_index(["template", "win"])
    idx = pd.MultiIndex.from_arrays([det["template"].astype(str), det["win"]])
    for k in ks:
        col = f"thr{int(k)}"
        det[col] = key[col].reindex(idx).values
        det[f"new{int(k)}"] = (det["cc_stack"].values >= det[col].values)
    det["old"] = det["cc_stack"].values >= OLD_THR
    return det


def decluster_old(det):
    """The OLD rule used a 6 s gap; re-apply it on the 0.80 subset (monotone-safe superset in hand)."""
    sub = det[det["old"]].sort_values(["template", "t"])
    keep = np.zeros(len(sub), bool)
    last = {}
    tv = sub["t"].values
    tm = sub["template"].astype(str).values
    for i in range(len(sub)):
        t0 = last.get(tm[i], -1e18)
        if tv[i] - t0 >= OLD_GAP_S:
            keep[i] = True
            last[tm[i]] = tv[i]
    return sub.index.values[keep]


def per_template_stats(det, daily, tremor, ks=KS):
    """Per-template ETS/quiet, tremor correlation and counts for the old and new rules.

    Written with groupby dictionaries: the naive per-template boolean scan over a ~1e6-row
    detection table x ~3e3 templates is O(1e9) and does not finish.
    """
    daily = daily.copy()
    daily["day"] = pd.to_datetime(daily["day"])
    daily["fam"] = daily["template"].astype(str)
    det = det.copy()
    det["fam"] = det["template"].astype(str)
    all_days = pd.DatetimeIndex(sorted(daily["day"].unique()))
    day_index = {d: i for i, d in enumerate(all_days)}
    old_keep = set(decluster_old(det).tolist())
    det["oldkeep"] = det.index.isin(old_keep)
    dgroups = {k: v for k, v in det.groupby("fam", sort=False)}
    lgroups = {k: v for k, v in daily.groupby("fam", sort=False)}
    cells = {}
    rows = []
    for fam in sorted(lgroups):
        lat, lon = cell_of(fam)
        ck = (round(lat, 3), round(lon, 3))
        if ck not in cells:
            cells[ck] = ets_quiet_mask(cell_tremor_daily(tremor, lat, lon), all_days)
        ets_m, qui_m, roll = cells[ck]
        dd = lgroups[fam].set_index("day").reindex(all_days)
        live = (dd["qc_drop"].fillna(1).values == 0)
        n_ets_d = int((ets_m & live).sum())
        n_qui_d = int((qui_m & live).sum())
        r = dict(template=fam, lat=lat, lon=lon, n_ets_days=n_ets_d, n_quiet_days=n_qui_d,
                 n_days=int(live.sum()),
                 mad_med=float(np.nanmedian(dd["mad"].values)),
                 med_med=float(np.nanmedian(dd["med"].values)),
                 n_comp=float(np.nanmedian(dd["n_comp"].values)))
        sub = dgroups.get(fam)
        for label in ["old"] + [f"new{int(k)}" for k in ks]:
            if sub is None or not len(sub):
                s = sub
                cnt = np.zeros(len(all_days))
                nsel = 0
                snrm = ccm = np.nan
            else:
                m = sub["oldkeep"].values if label == "old" else sub[label].values
                s = sub[m]
                nsel = len(s)
                cnt = np.zeros(len(all_days))
                if nsel:
                    vc = s["day"].value_counts()
                    for d, n in vc.items():
                        j = day_index.get(d)
                        if j is not None:
                            cnt[j] = n
                snrm = float(np.nanmedian(s["snr"].values)) if nsel else np.nan
                ccm = float(np.nanmedian(s["cc_stack"].values)) if nsel else np.nan
            rate_ets = cnt[ets_m & live].mean() if n_ets_d else np.nan
            rate_qui = cnt[qui_m & live].mean() if n_qui_d else np.nan
            ratio = (rate_ets / rate_qui) if (rate_qui and rate_qui > 0) else np.nan
            ok = live & np.isfinite(roll)
            rho, pv = (np.nan, np.nan)
            if ok.sum() > 30 and np.nanstd(roll[ok]) > 0 and np.nanstd(cnt[ok]) > 0:
                rho, pv = spearmanr(cnt[ok], roll[ok])
            r.update({f"n_{label}": int(nsel),
                      f"rate_ets_{label}": rate_ets, f"rate_quiet_{label}": rate_qui,
                      f"ratio_{label}": ratio, f"rho_{label}": rho, f"p_{label}": pv,
                      f"snr_med_{label}": snrm, f"cc_med_{label}": ccm})
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True)
    ap.add_argument("--prefix", required=True, help="stage-2 out prefix (without _detections.parquet)")
    ap.add_argument("--outdir", default="lfe_catalog_repair/results")
    a = ap.parse_args()
    det = pd.read_parquet(f"{a.prefix}_detections.parquet")
    daily = pd.read_parquet(f"{a.prefix}_daily.parquet")
    print(f"[{a.sta}] {len(det):,} stored detections, {len(daily):,} template-days, "
          f"{daily.template.nunique()} templates", flush=True)
    # COMPONENT-COUNT GUARD. On a day where a horizontal is missing, cc_stack collapses onto the
    # single-channel cc, whose MAD is ~2x larger, so the extraction floor falls back to the old 0.80
    # while the 30-day window threshold (dominated by the 3-C days) stays near 0.8-0.9 -> such days
    # would dump a burst of spurious "new" detections. Template-days below the run's modal component
    # count are therefore QC-dropped, and their detections removed, from BOTH numerator and
    # denominator.
    nmode = int(daily.loc[daily.qc_drop == 0, "n_comp"].mode().iloc[0])
    short = (daily["n_comp"] < nmode) & (daily["qc_drop"] == 0)
    print(f"[{a.sta}] modal n_comp = {nmode}; QC-dropping {int(short.sum()):,} "
          f"({100*short.mean():.2f} %) template-days with fewer components", flush=True)
    daily.loc[short, "qc_drop"] = 1
    # DEGENERATE-WINDOW GUARD. A day with a zero-filled / constant segment yields detections whose
    # 2 s window has literally zero amplitude (peak |Z| = 0 -> SNR = 0). At YACT one such day held
    # 50 % of all k=9 detections and dragged the median SNR to ~1e-8. Drop them.
    if "snr" in det:
        n0 = len(det)
        det = det[det["snr"] > 1e-3]
        if n0 - len(det):
            print(f"[{a.sta}] dropped {n0-len(det):,} zero-amplitude detections "
                  f"({100*(n0-len(det))/max(n0,1):.2f} %, degenerate/zero-filled windows)", flush=True)
    if "n_comp" in det:
        n0 = len(det)
        det = det[det["n_comp"] >= nmode]
        print(f"[{a.sta}] dropped {n0-len(det):,} detections on component-short days "
              f"({100*(n0-len(det))/max(n0,1):.2f} %)", flush=True)
    thr = build_thresholds(daily)
    # TRUNCATION DIAGNOSTIC: stage 2 extracted everything above min(med_day + 8*MAD_day, 0.80).
    # A 30-day-window threshold BELOW that floor would mean some window-rule detections were never
    # stored. Measure how often that happens instead of assuming it does not.
    dq = daily[daily.qc_drop == 0].copy()
    dq["template"] = dq["template"].astype(str)
    dq["day"] = pd.to_datetime(dq["day"])
    dq["win"] = windows_30d(dq["day"].values, pd.Timestamp(thr.attrs["day0"]))
    dq = dq.merge(thr[["template", "win", "thr8"]], on=["template", "win"], how="left")
    floor_day = np.minimum(dq["med"] + 8.0 * dq["mad"], OLD_THR)
    trunc = float((dq["thr8"] < floor_day - 1e-6).mean())
    print(f"[{a.sta}] truncation check: {100*trunc:.3f} % of template-days have "
          f"thr8(30-d window) BELOW the stage-2 extraction floor "
          f"(median gap {float(np.nanmedian(dq['thr8'] - floor_day)):+.4f})", flush=True)
    det = apply_rules(det, thr, thr.attrs["day0"])
    tremor = load_tremor()
    stats = per_template_stats(det, daily, tremor)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    thr.to_csv(out / f"d3_thresholds_{a.sta}.csv", index=False)
    keep_cols = [c for c in det.columns if c not in ("samp",)]
    det[keep_cols].to_parquet(out / f"d3_detections_{a.sta}.parquet", index=False)
    stats.to_csv(out / f"d3_validation_{a.sta}.csv", index=False)
    # station summary
    rows = []
    for label in ["old"] + [f"new{int(k)}" for k in KS]:
        resp = stats[f"ratio_{label}"]
        rows.append(dict(
            station=a.sta, rule=label,
            n_det=int(stats[f"n_{label}"].sum()),
            det_per_tmpl_day=float(stats[f"n_{label}"].sum() / max(stats["n_days"].sum(), 1)),
            n_tmpl_with_det=int((stats[f"n_{label}"] > 0).sum()),
            ratio_median=float(resp.median()),
            frac_ratio_ge2=float((resp >= 2).mean()),
            frac_ratio_ge5=float((resp >= 5).mean()),
            rho_median=float(stats[f"rho_{label}"].median()),
            frac_rho_pos=float((stats[f"rho_{label}"] > 0).mean()),
            frac_rho_sig=float((stats[f"p_{label}"] < 0.05).mean()),
            snr_median=float(stats[f"snr_med_{label}"].median()),
            cc_median=float(stats[f"cc_med_{label}"].median()),
        ))
        if label != "old":
            k = int(label[3:])
            rows[-1]["frac_thr_attainable"] = float(thr[f"attainable{k}"].mean())
            rows[-1]["thr_median"] = float(thr[f"thr{k}"].median())
        rows[-1]["frac_trunc_template_days"] = trunc
    summ = pd.DataFrame(rows)
    f = out / "d3_station_summary.csv"
    if f.exists():
        prev = pd.read_csv(f)
        prev = prev[prev.station != a.sta]
        summ = pd.concat([prev, summ], ignore_index=True)
    summ.to_csv(f, index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
