#!/usr/bin/env python
"""WA preparation (2) -- the HOUR-MATCHED, EARTHQUAKE-SCREENED quiet set.  Prefix `b_`.

Merlin's Step-1 audit defect F2: `w_common.quiet_hours` with `QUIET_HOURS_PER_DAY = 2` evaluates
`linspace(0, 23, 2)` = clock hours {0, 23} ONLY, while the candidate (tremor) hours span all 24
clock hours -- so daytime transients are never sampled by the null -- and no earthquake screen was
applied, which let one regional event at 2005-07-09 23:48 set T = 12.50 single-handed.
`results/freeze_wa_v1.json` inherits both defects (quiet_hours_per_day 2).

This script builds the replacement for Washington, for freeze_wa_v2 to adopt verbatim:

  * 40 hours on the 20 quiet days of the window,
  * clock hours CYCLING 0-23 so every clock hour appears at least once, with the remaining slots
    allocated by largest remainder to match the CANDIDATE-hour clock histogram,
  * hours containing a catalogued earthquake excluded: ComCat M >= 2.5 within 300 km of the
    sub-region centre (the hour +- 10 min), or M >= 5.5 anywhere (the hour and the two preceding
    hours, which covers teleseismic body- and surface-wave arrivals).  A catalog query is not a
    waveform download,
  * hours where more than half the bed's station-hours are vetoed by the health census
    (`b_wa_health.py`) skipped, when that table exists.

2010 (N-Olympic): quiet days and candidate hours from the PNSN rule, i.e. the same definitions
`freeze_wa_v1` used.  2008 (Puget/S-WA): there is no tremor catalogue before 2010-01-06, so the
candidate hours are the ARRAY-ENERGY hot hours (`results/hour_energy_2008_raw.csv` with the
frozen p95 of `results/probe_days_2008.json`) and the quiet days are the literature-quiet days of
`results/hour_energy_flag_2008.csv`.

Usage
  python scripts/b_wa_quiet.py --year 2010
  python scripts/b_wa_quiet.py --year 2008
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
CACHE = f"{STUDY}/cache20hz"
N_HOURS = 40
EQ_LOCAL_MAG, EQ_LOCAL_KM = 2.5, 300.0
EQ_GLOBAL_MAG = 5.5
EQ_PAD_MIN = 10.0
EQ_GLOBAL_LEAD_H = 2
VETO_FRAC_MAX = 0.5


def comcat(t0, t1, minmag, lat=None, lon=None, km=None):
    p = dict(format="csv", starttime=str(t0), endtime=str(t1), minmagnitude=minmag,
             orderby="time-asc")
    if lat is not None:
        p.update(latitude=lat, longitude=lon, maxradiuskm=km)
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urllib.parse.urlencode(p)
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                raw = r.read().decode()
            break
        except Exception as e:                                              # noqa: BLE001
            print(f"[eq] try {k+1}: {type(e).__name__}: {str(e)[:120]}")
            time.sleep(5 * (k + 1))
    else:
        raise RuntimeError("ComCat unreachable")
    if not raw.strip() or "\n" not in raw.strip():
        return pd.DataFrame(columns=["time", "latitude", "longitude", "mag"])
    E = pd.read_csv(pd.io.common.StringIO(raw))
    E["time"] = pd.to_datetime(E.time, utc=True, format="mixed")
    return E


def hours_2010():
    """(quiet days, candidate-hour frame, sub-region centre, bed stations) for the WA 2010 bed."""
    F = json.load(open(f"{RES}/freeze_wa_v1.json"))
    import x_common as X
    X.activate()
    import w_common as W
    H = W.candidate_hours()
    stas = [m for v in F["sites"]["site_members"].values() for m in v]
    return (list(F["window"]["quiet_days"]), H,
            (F["sub_region"]["centre"]["lat"], F["sub_region"]["centre"]["lon"]), stas)


def hours_2008():
    """Same for the Puget/S-WA 2008 bed, from the ARRAY-ENERGY flag (no tremor catalogue)."""
    P = json.load(open(f"{RES}/probe_days_2008.json"))
    Fl = pd.read_csv(f"{RES}/hour_energy_flag_2008.csv")
    R = pd.read_csv(f"{RES}/hour_energy_2008_raw.csv")
    R = R[R.status == "ok"]
    hcols = [f"h{h:02d}" for h in range(24)]
    quiet_days = list(Fl[Fl.is_quiet].date)
    # per station: normalise by its own median quiet-day hourly value; network median per hour
    Q = R[R.date.isin(quiet_days)]
    med = Q.groupby("netsta")[hcols].median().median(axis=1)
    rows = []
    for r in R.itertuples():
        m = float(med.get(r.netsta, np.nan))
        if not np.isfinite(m) or m <= 0:
            continue
        for h, c in enumerate(hcols):
            rows.append(dict(netsta=r.netsta, date=r.date, hour=h,
                             rel=float(getattr(r, c)) / m))
    D = pd.DataFrame(rows)
    N = D.groupby(["date", "hour"], as_index=False).rel.median()
    p95 = float(P["p95"])
    tremor_days = set(Fl[Fl.tremor_flag].date)
    H = N[(N.rel > p95) & (N.date.isin(tremor_days))][["date", "hour"]].copy()
    H["n_epi"] = np.nan
    stas = sorted({p.split("/")[-3].split(".")[-1] for p in glob.glob(f"{CACHE}/*/2008/*.npz")})
    lat = np.nan
    try:
        import d_common as D2
        ST = pd.read_csv(f"{RES}/sites_array_era_final.csv", keep_default_na=False)
        sub = ST[(ST.lat >= 46.2) & (ST.lat <= 47.3)]
        lat, lon = float(sub.lat.median()), float(sub.lon.median())
    except Exception:                                                       # noqa: BLE001
        lat, lon = 46.75, -122.9
    return quiet_days, H.reset_index(drop=True), (lat, lon), stas


def veto_fraction(year, stas):
    """{(date, hour): fraction of the bed's cached station-hours vetoed}, or {} if no census."""
    p = f"{RES}/wa_station_hour_veto_{year}.csv"
    if not os.path.exists(p):
        return {}
    V = pd.read_csv(p)
    keep = {s.split(".")[-1] for s in stas}
    V = V[V.netsta.str.split(".").str[-1].isin(keep)]
    if not len(V):
        return {}
    G = V.groupby(["date", "hour"]).veto.mean()
    return {(d, int(h)): float(v) for (d, h), v in G.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True, choices=[2008, 2010])
    ap.add_argument("--n-hours", type=int, default=N_HOURS)
    a = ap.parse_args()
    quiet_days, H, centre, stas = hours_2010() if a.year == 2010 else hours_2008()
    quiet_days = sorted(quiet_days)
    print(f"[bed] {a.year}: {len(quiet_days)} quiet days, {len(H)} candidate hours, "
          f"centre {centre}, {len(stas)} bed stations")

    # ---------------------------------------------------------------- candidate histogram -----
    hist = np.zeros(24)
    for h, k in H.groupby("hour").size().items():
        hist[int(h)] = k
    p = hist / max(hist.sum(), 1)

    # ---------------------------------------------------------------- earthquake screen -------
    t0 = pd.Timestamp(quiet_days[0]) - pd.Timedelta(hours=3)
    t1 = pd.Timestamp(quiet_days[-1]) + pd.Timedelta(days=1)
    EL = comcat(t0.isoformat(), t1.isoformat(), EQ_LOCAL_MAG, centre[0], centre[1], EQ_LOCAL_KM)
    EG = comcat(t0.isoformat(), t1.isoformat(), EQ_GLOBAL_MAG)
    print(f"[eq] {len(EL)} local M>={EQ_LOCAL_MAG} within {EQ_LOCAL_KM:.0f} km, "
          f"{len(EG)} global M>={EQ_GLOBAL_MAG}")
    excl = {}
    for r in EL.itertuples():
        t = r.time
        for d in (-EQ_PAD_MIN, 0.0, EQ_PAD_MIN):
            tt = t + pd.Timedelta(minutes=d)
            excl.setdefault((str(tt.date()), int(tt.hour)), []).append(
                f"local M{r.mag:.1f} {t.strftime('%H:%M')}")
    for r in EG.itertuples():
        t = r.time
        for k in range(0, EQ_GLOBAL_LEAD_H + 1):
            tt = t + pd.Timedelta(hours=k)
            excl.setdefault((str(tt.date()), int(tt.hour)), []).append(
                f"global M{r.mag:.1f} {t.strftime('%m-%d %H:%M')}")
    vf = veto_fraction(a.year, stas)
    cached = {}
    for s in stas:
        for pth in glob.glob(f"{CACHE}/*.{s}/{a.year}/*.npz"):
            doy = int(os.path.basename(pth).split(".")[0])
            d = str((pd.Timestamp(year=a.year, month=1, day=1)
                     + pd.Timedelta(days=doy - 1)).date())
            cached[d] = cached.get(d, 0) + 1

    # ---------------------------------------------------------------- quota per clock hour ----
    n = a.n_hours
    quota = np.ones(24, int)
    rem = n - 24
    if rem > 0:
        want = p * rem
        base = np.floor(want).astype(int)
        left = rem - base.sum()
        order = np.argsort(-(want - base))
        base[order[:left]] += 1
        quota += base
    elif rem < 0:
        raise SystemExit("n_hours must be >= 24 for the cycling rule")

    # ---------------------------------------------------------------- pick day x hour ---------
    rows, skipped = [], []
    used_per_day = {d: 0 for d in quiet_days}
    nd = len(quiet_days)
    for h in range(24):
        for k in range(int(quota[h])):
            # rotate the day order so the hours spread evenly over the quiet days
            cand = sorted(quiet_days, key=lambda d: (used_per_day[d],
                                                     (quiet_days.index(d) - (h + k)) % nd))
            placed = False
            for d in cand:
                key = (d, h)
                if any(r["date"] == d and r["hour"] == h for r in rows):
                    continue
                if cached.get(d, 0) == 0:
                    skipped.append(dict(date=d, hour=h, why="no cached station-day"))
                    continue
                if key in excl:
                    skipped.append(dict(date=d, hour=h, why="earthquake: " + ";".join(excl[key])))
                    continue
                f = vf.get(key, 0.0)
                if f > VETO_FRAC_MAX:
                    skipped.append(dict(date=d, hour=h, why=f"health veto frac {f:.2f}"))
                    continue
                rows.append(dict(date=d, hour=h, clock_quota=int(quota[h]),
                                 candidate_frac=float(p[h]), veto_frac=f,
                                 n_cached_stations=int(cached.get(d, 0))))
                used_per_day[d] += 1
                placed = True
                break
            if not placed:
                skipped.append(dict(date="", hour=h, why="no eligible quiet day for this hour"))
    Q = pd.DataFrame(rows).sort_values(["date", "hour"]).reset_index(drop=True)
    out = f"{RES}/wa_quiet_hours_{a.year}.csv"
    Q.to_csv(out, index=False)
    pd.DataFrame(skipped).to_csv(f"{RES}/wa_quiet_hours_{a.year}_skipped.csv", index=False)
    meta = dict(year=a.year, n_hours=int(len(Q)), n_target=n, n_quiet_days=len(quiet_days),
                n_days_used=int(Q.date.nunique()) if len(Q) else 0,
                clock_hours_present=sorted(set(int(x) for x in Q.hour)) if len(Q) else [],
                candidate_hours=int(len(H)),
                candidate_hour_histogram={int(k): int(v) for k, v in
                                          zip(range(24), hist.astype(int))},
                quota={int(k): int(v) for k, v in zip(range(24), quota)},
                centre=dict(lat=float(centre[0]), lon=float(centre[1])),
                eq_rule=f"ComCat M>={EQ_LOCAL_MAG} within {EQ_LOCAL_KM:.0f} km of the centre "
                        f"(hour +- {EQ_PAD_MIN:.0f} min) or M>={EQ_GLOBAL_MAG} anywhere (the hour "
                        f"and the {EQ_GLOBAL_LEAD_H} following hours)",
                n_eq_local=int(len(EL)), n_eq_global=int(len(EG)),
                n_excluded_hour_slots=int(len(excl)),
                health_census=os.path.exists(f"{RES}/wa_station_hour_veto_{a.year}.csv"),
                veto_frac_max=VETO_FRAC_MAX,
                supersedes="freeze_wa_v1 quiet_hours_per_day=2 -> clock hours {0,23} only (F2)")
    json.dump(meta, open(f"{RES}/wa_quiet_hours_{a.year}.json", "w"), indent=1, default=str)
    print(f"[quiet] {len(Q)} hours on {meta['n_days_used']} days -> {out}")
    print("        clock hours:", meta["clock_hours_present"])
    print("        per-day counts:", Q.groupby("date").size().to_dict() if len(Q) else {})
    print(f"        excluded slots (earthquake): {len(excl)}; skipped attempts: {len(skipped)}")


if __name__ == "__main__":
    main()
