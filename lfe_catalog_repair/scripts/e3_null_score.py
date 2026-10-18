#!/usr/bin/env python
"""E3 step 1 — score the chance-association null against the observed PyOcto catalogue.

For every (set, model, threshold) and every station cut N>=3/4/5 it reports the OBSERVED event rate
and the NULL event rate (mean +- sd over the 5 shifted-pick seeds), separately on ETS days and on
quiet days, plus the null/observed ratio that the audit's branch rule is read from.

ETS / quiet definition -- Merlin F10 fix.  E2 built the 30-d rolling tremor sum INSIDE the 34-d
association window with min_periods=5, so the window's first and last two weeks had incomplete
sums and "quiet" degenerated to the window edges.  Here the daily tremor count within
`--radius-km` (default 60) of the NETWORK CENTROID is built over the FULL PNSN record
(2010-01-06..2026-05-27), the 30-d centred rolling sum is taken on that full series (so days inside
the pilot window see their true 30-d neighbourhood, including data outside the window), and the
p90 / p25 thresholds are percentiles of the FULL-RECORD rolling series.  Days in the pilot window
are then labelled against those fixed thresholds.

Usage:
  python lfe_catalog_repair/scripts/e3_null_score.py                     # all sets, writes csv
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
TREMOR = "catalogs/pnsn_tremor_cascadia_full.csv"
R_KM_DEG = 111.195
# the station sets, from e2_pick.SETS (kept literal so this script needs no torch import)
SET_STATIONS = {
    "vi2010": ["PGC", "SNB", "VGZ", "LZB", "PFB", "YOUB", "NLLB", "B011", "B926", "LRIV"],
    "vi2017": ["PGC", "SYMB", "SNB", "VGZ", "CLRS", "PTRF", "B011", "LRIV"],
}
NCUTS = (3, 4, 5)


def centroid(setname):
    st = pd.read_csv(f"{RES}/e2_stations.csv").set_index("sta")
    s = st.loc[SET_STATIONS[setname]]
    return float(s.lat.mean()), float(s.lon.mean())


def ets_quiet_labels(lat0, lon0, radius_km=60.0, win=30, p_hi=90, p_lo=25):
    """FULL-RECORD 30-d rolling tremor count within radius_km of (lat0, lon0) -> per-day label."""
    tr = pd.read_csv(TREMOR, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time).astype("datetime64[ns]")
    d = np.hypot((tr.lon.values - lon0) * R_KM_DEG * np.cos(np.radians(lat0)),
                 (tr.lat.values - lat0) * R_KM_DEG)
    tr = tr[d <= radius_km]
    cnt = tr.groupby(tr.time.dt.normalize()).size()
    full = pd.date_range(cnt.index.min(), cnt.index.max(), freq="D")
    cnt = cnt.reindex(full, fill_value=0)
    roll = cnt.rolling(win, min_periods=win, center=True).sum()   # full 30-d windows only
    hi, lo = np.nanpercentile(roll.dropna(), p_hi), np.nanpercentile(roll.dropna(), p_lo)
    lab = pd.Series(np.where(roll >= hi, "ets", np.where(roll <= lo, "quiet", "mid")),
                    index=full)
    lab[roll.isna()] = "na"
    return lab, cnt, roll, float(hi), float(lo)


def rates(evfile, lab, ncut, t0=None, t1=None):
    ev = pd.read_csv(evfile, usecols=["time", "n_stations"])
    ev = ev[ev["n_stations"] >= ncut]
    if len(ev) == 0:
        return dict(n_events=0, per_day=0.0, ets_per_day=np.nan, quiet_per_day=np.nan,
                    mid_per_day=np.nan, n_days=0, n_ets_days=0, n_quiet_days=0, n_mid_days=0)
    day = pd.to_datetime(ev.time, unit="s").dt.normalize()
    c = day.value_counts()
    lo = t0 if t0 is not None else day.min()
    hi = t1 if t1 is not None else day.max()
    days = pd.date_range(lo, hi, freq="D")
    c = c.reindex(days, fill_value=0)
    L = lab.reindex(days).fillna("na")
    out = dict(n_events=int(c.sum()), per_day=float(c.mean()),
               n_days=len(days), n_ets_days=int((L == "ets").sum()),
               n_quiet_days=int((L == "quiet").sum()), n_mid_days=int((L == "mid").sum()))
    out["ets_per_day"] = float(c[L == "ets"].mean()) if out["n_ets_days"] else np.nan
    out["quiet_per_day"] = float(c[L == "quiet"].mean()) if out["n_quiet_days"] else np.nan
    out["mid_per_day"] = float(c[L == "mid"].mean()) if out["n_mid_days"] else np.nan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius-km", type=float, default=60.0)
    ap.add_argument("--prefix", default="e3_nullevents",
                    help="basename prefix of the null event files to score "
                         "(e3_nullevents = the per-station-constant null; e3_nullhr = the "
                         "per-station-per-hour sensitivity variant)")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    # observed configurations to score = the (set, model, thr) for which nulls exist
    out = args.out or f"{RES}/{args.prefix.replace('e3_', 'e3_')}_rates.csv"
    nulls = sorted(glob.glob(f"{RES}/{args.prefix}_*.csv"))
    keys = {}
    for p in nulls:
        m = re.match(args.prefix + r"_(vi\d+)_(\w+?)_t([\d.]+)_s(\d+)\.csv$",
                     os.path.basename(p))
        keys.setdefault((m.group(1), m.group(2), m.group(3)), []).append((int(m.group(4)), p))
    labs = {}
    rows = []
    for (setname, model, thr), seeds in sorted(keys.items()):
        if setname not in labs:
            labs[setname] = ets_quiet_labels(*centroid(setname), radius_km=args.radius_km)
        lab, cnt, roll, hi, lo = labs[setname]
        obsf = f"{RES}/e2_events_{setname}_{model}_t{thr}.csv"
        # the observed window sets the day range for both observed and null (the null shifts move
        # a handful of picks past the window edge; using one common range keeps rates comparable)
        oev = pd.read_csv(obsf, usecols=["time"])
        d0 = pd.to_datetime(oev.time.min(), unit="s").normalize()
        d1 = pd.to_datetime(oev.time.max(), unit="s").normalize()
        for ncut in NCUTS:
            o = rates(obsf, lab, ncut, d0, d1)
            ns = [rates(p, lab, ncut, d0, d1) for _, p in sorted(seeds)]
            r = dict(set=setname, model=model, threshold=float(thr), n_min=ncut,
                     n_seeds=len(ns), tremor_radius_km=args.radius_km,
                     roll_p90=hi, roll_p25=lo,
                     n_days=o["n_days"], n_ets_days=o["n_ets_days"],
                     n_quiet_days=o["n_quiet_days"], n_mid_days=o["n_mid_days"],
                     obs_events=o["n_events"], obs_per_day=o["per_day"],
                     obs_ets_per_day=o["ets_per_day"], obs_quiet_per_day=o["quiet_per_day"])
            r["obs_mid_per_day"] = o["mid_per_day"]
            for k, ok in (("all", "per_day"), ("ets", "ets_per_day"),
                          ("quiet", "quiet_per_day"), ("mid", "mid_per_day")):
                v = np.array([x[ok] for x in ns], float)
                r[f"null_{k}_per_day"] = float(np.nanmean(v))
                r[f"null_{k}_sd"] = float(np.nanstd(v, ddof=1)) if len(v) > 1 else np.nan
                r[f"ratio_{k}"] = float(np.nanmean(v) / o[ok]) if o[ok] and np.isfinite(o[ok]) \
                    else np.nan
            r["obs_ets_quiet"] = (o["ets_per_day"] / o["quiet_per_day"]
                                  if o["quiet_per_day"] else np.nan)
            r["null_ets_quiet"] = (r["null_ets_per_day"] / r["null_quiet_per_day"]
                                   if r["null_quiet_per_day"] else np.nan)
            rows.append(r)
    df = pd.DataFrame(rows)
    df["null_kind"] = args.prefix
    df.to_csv(out, index=False)
    cols = ["set", "model", "threshold", "n_min", "n_ets_days", "n_quiet_days", "n_mid_days",
            "obs_per_day", "null_all_per_day", "null_all_sd", "ratio_all",
            "obs_ets_per_day", "null_ets_per_day", "null_ets_sd", "ratio_ets",
            "obs_quiet_per_day", "null_quiet_per_day", "null_quiet_sd", "ratio_quiet",
            "obs_mid_per_day", "null_mid_per_day", "ratio_mid",
            "obs_ets_quiet", "null_ets_quiet"]
    print(df[cols].round(3).to_string(index=False))
    print(f"\n[write] {out}")
    for s, (lab, cnt, roll, hi, lo) in labs.items():
        print(f"[{s}] full-record 30-d rolling tremor within {args.radius_km:.0f} km of the "
              f"network centroid: p90={hi:.0f}, p25={lo:.0f} events/30 d "
              f"(record {cnt.index.min().date()}..{cnt.index.max().date()}, "
              f"{len(cnt)} d; labels ets/quiet/mid/na = "
              f"{lab.value_counts().reindex(['ets','quiet','mid','na']).fillna(0).astype(int).to_dict()})")


if __name__ == "__main__":
    main()
