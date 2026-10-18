#!/usr/bin/env python
"""STEP 1c' STAGE 3 -- the HOUR-MATCHED, EARTHQUAKE-SCREENED QUIET SET (Merlin F2).

Step 1c drew its 40 quiet hours as `linspace(0, 23, 2)` = clock hours {0, 23} on 20 literature-quiet
days, while the 299 candidate hours span all 24 clock hours (5-18 per hour): the null never sampled
daytime transients.  Worse, one hour (2005-07-09 23h) contained the teleseismic arrivals of an
M 5.8 near Japan and set T = 12.50 on its own.

This stage:
  * allocates 40 quiet hours over the 20 quiet days so that EVERY clock hour appears at least once
    and the remaining slots follow the candidate-hour histogram (largest remainder);
  * screens every hour against an earthquake CATALOG (ComCat + NRCan; a catalog query, never a
    waveform download): a LOCAL event M >= 2.5 within 300 km with its origin inside the hour, or a
    GLOBAL event M >= 5.5 whose predicted P (iasp91) or 3.5 km/s surface-wave arrival falls in the
    hour;
  * requires >= MIN_SITES sites to survive the station-hour health veto in the hour;
  * assigns each clock hour to the day with the fewest slots used, so no day carries the same hour
    twice.

Writes results/step1c_prime/quiet_hours.csv (the list the engine reads), _quiet_screen.csv (every
day x hour with its verdict) and _quiet_events.csv (the catalog rows that vetoed something).
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

LOCAL_MAG, LOCAL_KM = 2.5, 300.0
GLOBAL_MAG = 5.5
PRE_PAD_S = 300.0            # a local event just before the hour still rings inside it
SURF_V = 3.5                 # km/s surface-wave proxy
COMCAT = "https://earthquake.usgs.gov/fdsnws/event/1/query"
NRCAN = "https://earthquakescanada.nrcan.gc.ca/fdsnws/event/1/query"


def _get(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def comcat(t0, t1, **kw):
    q = "&".join(f"{k}={v}" for k, v in kw.items())
    txt = _get(f"{COMCAT}?format=geojson&starttime={t0}&endtime={t1}&{q}")
    d = json.loads(txt)
    rows = []
    for f in d.get("features", []):
        p, g = f["properties"], f["geometry"]["coordinates"]
        rows.append(dict(src="comcat", t=float(p["time"]) / 1000.0, mag=p.get("mag"),
                         lat=float(g[1]), lon=float(g[0]), dep=float(g[2] or 0.0),
                         place=str(p.get("place"))[:80]))
    return pd.DataFrame(rows)


def nrcan(t0, t1, minmag):
    txt = _get(f"{NRCAN}?format=text&starttime={t0}&endtime={t1}&minmagnitude={minmag}")
    rows = []
    for ln in txt.splitlines():
        if not ln or ln.startswith("#"):
            continue
        p = ln.split("|")
        if len(p) < 7:
            continue
        try:
            rows.append(dict(src="nrcan", t=float(pd.Timestamp(p[1]).timestamp()),
                             mag=float(p[6]), lat=float(p[2]), lon=float(p[3]),
                             dep=float(p[4] or 0.0), place=p[-1][:80]))
        except Exception:                                                   # noqa: BLE001
            continue
    return pd.DataFrame(rows)


def arrivals(E, lat0, lon0):
    """Add epicentral distance and the predicted first-arrival / surface-wave times."""
    from obspy.taup import TauPyModel
    m = TauPyModel("iasp91")
    d_km, t_p, t_s = [], [], []
    for r in E.itertuples():
        dk = float(S.epi_km(lat0, lon0, r.lat, r.lon))
        d_km.append(dk)
        deg = dk / 111.195
        tp = np.nan
        try:
            a = m.get_travel_times(source_depth_in_km=max(float(r.dep), 0.0),
                                   distance_in_degree=deg, phase_list=["p", "P", "Pn", "Pdiff"])
            if a:
                tp = float(a[0].time)
        except Exception:                                                   # noqa: BLE001
            pass
        t_p.append(tp)
        t_s.append(dk / SURF_V)
    E = E.copy()
    E["dist_km"] = d_km
    E["t_p_s"] = t_p
    E["t_surf_s"] = t_s
    return E


def screen_hour(E, t_start, t_end):
    """(vetoed, [reasons], [event indices]) for one hour."""
    why, idx = [], []
    for i, r in enumerate(E.itertuples()):
        loc = (r.dist_km <= LOCAL_KM) and (r.mag is not None) and (r.mag >= LOCAL_MAG)
        glo = (r.mag is not None) and (r.mag >= GLOBAL_MAG)
        if loc and (t_start - PRE_PAD_S) <= r.t <= t_end:
            why.append(f"local M{r.mag:.1f} {r.dist_km:.0f} km")
            idx.append(i)
            continue
        if glo:
            for lab, dt in (("P", r.t_p_s), ("surf", r.t_surf_s)):
                if np.isfinite(dt) and (t_start - 60.0) <= (r.t + dt) <= t_end:
                    why.append(f"teleseismic M{r.mag:.1f} {lab} arrival ({r.dist_km:.0f} km)")
                    idx.append(i)
                    break
    return bool(why), why, idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/step1c_prime/quiet_hours.csv")
    ap.add_argument("--n-hours", type=int, default=W.QUIET_MIN_HOURS)
    ap.add_argument("--veto", default=f"{W.RES}/step1c_prime/station_hour_veto.csv")
    ap.add_argument("--candidate-hours", default=f"{W.RES}/step1c_candidate_hours.csv")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    mem, of, order, det = W.sites()
    days = sorted(W.STEP1C_QUIET_DAYS)
    lat0, lon0 = W.GRID_LAT0, W.GRID_LON0

    # ---------------------------------------------------------------- catalogs
    t0 = (pd.Timestamp(days[0]) - pd.Timedelta(days=1)).strftime("%Y-%m-%d")
    t1 = (pd.Timestamp(days[-1]) + pd.Timedelta(days=2)).strftime("%Y-%m-%d")
    parts, srcs = [], []
    try:
        parts.append(comcat(t0, t1, minmagnitude=LOCAL_MAG, latitude=lat0, longitude=lon0,
                            maxradiuskm=int(LOCAL_KM)))
        srcs.append("comcat_local")
    except Exception as e:                                                  # noqa: BLE001
        W.log(f"[quiet] ComCat local query FAILED: {e}")
    try:
        parts.append(comcat(t0, t1, minmagnitude=GLOBAL_MAG))
        srcs.append("comcat_global")
    except Exception as e:                                                  # noqa: BLE001
        W.log(f"[quiet] ComCat global query FAILED: {e}")
    try:
        N = nrcan(t0, t1, LOCAL_MAG)
        if len(N):
            parts.append(N)
        srcs.append("nrcan")
    except Exception as e:                                                  # noqa: BLE001
        W.log(f"[quiet] NRCan query FAILED: {e}")
    if not parts:
        raise SystemExit("no earthquake catalog could be queried -- the quiet set cannot be "
                         "screened; stop and report")
    E = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["t", "lat", "lon"])
    E = arrivals(E, lat0, lon0)
    E = E.sort_values("t").reset_index(drop=True)
    E.to_csv(a.out.replace(".csv", "_events.csv"), index=False)
    W.log(f"[quiet] {len(E)} catalog events from {srcs} "
          f"({int((E.dist_km <= LOCAL_KM).sum())} within {LOCAL_KM:.0f} km)")

    # ---------------------------------------------------------------- health availability
    V = set()
    if os.path.exists(a.veto):
        D = pd.read_csv(a.veto)
        D = D[D.veto.astype(bool)]
        V = {(str(r.sta), str(r.date), int(r.hour)) for r in D.itertuples()}
    idx = W.cache_idx()

    def n_sites_ok(d, h):
        ts = pd.Timestamp(d)
        ss = [s for s in det if (s, int(ts.year), int(ts.dayofyear)) in idx
              and (s, d, h) not in V]
        return len({of[s] for s in ss})

    # ---------------------------------------------------------------- the target histogram
    CH = pd.read_csv(a.candidate_hours)
    hist = CH.groupby("hour").size().reindex(range(24), fill_value=0).to_numpy(float)
    n = int(a.n_hours)
    base = np.ones(24)
    rem = n - 24
    share = hist / hist.sum() * rem if hist.sum() > 0 else np.full(24, rem / 24.0)
    add = np.floor(share).astype(int)
    left = rem - int(add.sum())
    for k in np.argsort(-(share - add))[:max(left, 0)]:
        add[k] += 1
    target = (base + add).astype(int)
    W.log(f"[quiet] target clock-hour counts {dict(enumerate(target.tolist()))}")

    # ---------------------------------------------------------------- screen every day x hour
    rows = []
    for d in days:
        ts = pd.Timestamp(d)
        for h in range(24):
            t_s = float(ts.timestamp()) + h * 3600.0
            veto_eq, why, ei = screen_hour(E, t_s, t_s + 3600.0)
            ns = n_sites_ok(d, h)
            rows.append(dict(date=d, hour=h, n_sites_healthy=ns, eq_veto=veto_eq,
                             eq_reason=";".join(why), usable=bool((not veto_eq)
                                                                  and ns >= W.MIN_SITES)))
    SC = pd.DataFrame(rows)
    SC.to_csv(a.out.replace(".csv", "_screen.csv"), index=False)

    # ---------------------------------------------------------------- allocate
    want = [h for h in range(24) for _ in range(int(target[h]))]
    used = collections.Counter()
    taken = set()
    chosen = []
    unfilled = []
    ok = {(r.date, int(r.hour)) for r in SC[SC.usable].itertuples()}
    for h in want:
        cand = [d for d in days if (d, h) in ok and (d, h) not in taken]
        if not cand:
            unfilled.append(h)
            continue
        d = min(cand, key=lambda x: (used[x], x))
        used[d] += 1
        taken.add((d, h))
        chosen.append(dict(date=d, hour=h))
    Q = pd.DataFrame(chosen).sort_values(["date", "hour"]).reset_index(drop=True)
    Q.to_csv(a.out, index=False)

    got = Q.groupby("hour").size().reindex(range(24), fill_value=0)
    out = dict(
        n_quiet_days=len(days), n_quiet_hours=int(len(Q)), requested=n,
        clock_hours_covered=int((got > 0).sum()),
        target_counts={int(k): int(v) for k, v in enumerate(target)},
        achieved_counts={int(k): int(v) for k, v in got.items()},
        unfilled_clock_hours=unfilled,
        hours_per_day={k: int(v) for k, v in Q.groupby("date").size().items()},
        n_hours_eq_vetoed=int(SC.eq_veto.sum()),
        n_hours_health_short=int((SC.n_sites_healthy < W.MIN_SITES).sum()),
        eq_vetoed_hours=[f"{r.date} {int(r.hour):02d}h: {r.eq_reason}"
                         for r in SC[SC.eq_veto].itertuples()],
        catalog_sources=srcs, n_catalog_events=int(len(E)),
        rule=dict(local_mag=LOCAL_MAG, local_km=LOCAL_KM, global_mag=GLOBAL_MAG,
                  pre_pad_s=PRE_PAD_S, surface_v_km_s=SURF_V,
                  hour_match="every clock hour >= 1; remainder by largest remainder on the "
                             "candidate-hour histogram"),
        verdict=("ACCEPT" if (len(Q) >= 24 and int((got > 0).sum()) == 24) else "BELOW_ACCEPT"),
        reason=("" if (len(Q) >= 24 and int((got > 0).sum()) == 24)
                else f"only {int((got > 0).sum())}/24 clock hours could be filled"))
    json.dump(out, open(a.out.replace(".csv", ".json"), "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
