#!/usr/bin/env python
"""Royer & Bostock (2014) sec.3.3 "scanning through all available tremor dates" -- the scan-day list.

The paper does not list the dates.  sec.2.1: NW recorded "major tremor episodes in 2007, 2008, 2009,
2010 and 2011", SW "only data from the 2007 and 2008 tremor episodes".  Our rule (documented in
results/ab_royer/POOL_RECONSTRUCTION.md):

  1. every sec.3.1 network-autocorrelation day of the subarray is a scan day;
  2. where a machine-readable tremor LOCATION catalogue exists (PNSN/Wech, from 2009-08-06), a UTC
     day is a tremor day if it holds >= MIN_EPI PNSN epicentres (2.5-min windows, so >= 1 hour of
     tremor) within RADIUS_KM of any station of the subarray's 25-station pool.  NW: 2009-08-06 ..
     2011-12-31.
  3. before 2009-08-06 no location catalogue is public, so the scan days are the PNSN published ETS
     windows (Wech's WECC, https://assets.pnsn.org/oldets/compare_maps/Density.html):
     winter 2007 Jan 12 - Feb 2, spring 2008 Apr 30 - Jun 2, spring 2009 Apr 28 - Jun 3 (NW only).
     These windows are for the whole N-Washington / S-Vancouver-Island ETS, so they cannot be
     restricted to the subarray region (over-inclusive: extra days only cost compute).

  python scripts/ab_tremor_days.py  ->  results/ab_royer/tremor/scan_days_{NW,SW}.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

RADIUS_KM = 25.0
MIN_EPI = 24
ETS_WINDOWS = {"NW": [("2007-01-12", "2007-02-02"), ("2008-04-30", "2008-06-02"),
                      ("2009-04-28", "2009-06-03")],
               "SW": [("2007-01-12", "2007-02-02"), ("2008-04-30", "2008-06-02")],
               # NVI, sec.2.2: "Two tremor episodes in 2006 and 2007 were available for analysis,
               # each episode lasting under 1 week" -> the scan days ARE the sec.3.1 days
               "NVI": []}
CAT_SPAN = {"NW": ("2009-08-06", "2011-12-31"), "SW": None, "NVI": None}
OUT = f"{A.RES}/tremor"


def catalogue():
    a = pd.read_csv(f"{OUT}/pnsn_tremor_2009.csv", skipinitialspace=True)
    a.columns = [c.strip() for c in a.columns]
    a = a.rename(columns={"starttime": "time"})[["time", "lat", "lon"]]
    b = pd.read_csv(f"{A.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv",
                    usecols=["time", "lat", "lon"])
    T = pd.concat([a, b], ignore_index=True)
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    return T.dropna().drop_duplicates()


def main():
    os.makedirs(OUT, exist_ok=True)
    T = catalogue()
    report = []
    for sub in ("NW", "SW", "NVI"):
        rows = {d: dict(date=d, autocorr=True, n_epi=np.nan) for d in A.days(sub)}
        # before 2009-08-06: OUR Wech & Creager (2008)-style detector (ab_wecc.py), calibrated on
        # the PNSN catalogue era, when its day list exists; otherwise the PNSN ETS windows
        wp = f"{OUT}/wecc_days_{sub}.csv"
        W = pd.read_csv(wp) if os.path.exists(wp) and sub in ("NW", "SW") else None
        ets = {d for a, b in ETS_WINDOWS[sub] for d in pd.date_range(a, b).strftime("%Y-%m-%d")}
        if W is not None:
            cut = CAT_SPAN[sub][0] if CAT_SPAN[sub] else "2099-01-01"
            W = W[W.date < cut]
            for r in W.itertuples():
                if r.tremor_day:
                    rows.setdefault(r.date, dict(date=r.date, autocorr=False, n_epi=np.nan))
                    rows[r.date]["wecc_tremor_day"] = True
                    rows[r.date]["wecc_windows"] = int(r.n_tremor_windows)
            for d in rows:
                rows[d]["ets_window"] = d in ets
            report.append(dict(subarray=sub, rule="WECC days before 2009-08-06",
                               n_days=int(W.tremor_day.sum())))
        else:
            for d in sorted(ets):
                rows.setdefault(d, dict(date=d, autocorr=False, n_epi=np.nan))
                rows[d]["ets_window"] = True
        if CAT_SPAN[sub]:
            a, b = CAT_SPAN[sub]
            C = T[(T.time >= a) & (T.time < pd.Timestamp(b) + pd.Timedelta("1D"))]
            la = np.array([A.COORD[s][0] for s in A.POOL[sub]])
            lo = np.array([A.COORD[s][1] for s in A.POOL[sub]])
            dk = A.epi_km(C.lat.to_numpy()[:, None], C.lon.to_numpy()[:, None], la, lo).min(1)
            n = C[dk <= RADIUS_KM].groupby(C.time[dk <= RADIUS_KM].dt.strftime("%Y-%m-%d")).size()
            for thr in (1, 12, 24, 48):
                report.append(dict(subarray=sub, rule=f">={thr} epicentres within {RADIUS_KM} km",
                                   n_days=int((n >= thr).sum())))
            for d, k in n.items():
                if k >= MIN_EPI:
                    rows.setdefault(d, dict(date=d, autocorr=False, n_epi=k))
                if d in rows:
                    rows[d]["n_epi"] = int(k)
        D = pd.DataFrame(sorted(rows.values(), key=lambda r: r["date"]))
        D["year"] = D.date.str[:4]
        D.to_csv(f"{OUT}/scan_days_{sub}.csv", index=False)
        A.log(f"[tremor days {sub}] {len(D)} scan days; by year "
              + D.groupby("year").size().to_dict().__repr__())
    pd.DataFrame(report).to_csv(f"{OUT}/rule_sensitivity.csv", index=False)
    A.log("\n" + pd.DataFrame(report).to_string(index=False))


if __name__ == "__main__":
    main()
