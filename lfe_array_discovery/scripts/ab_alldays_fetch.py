#!/usr/bin/env python
"""ALL-DAYS detection pass -- prefetch + 40 sps cache of every pool station-day in each subarray's
array lifetime (no tremor-day selection).  Prefix `ab_alldays_`.

The literal Royer & Bostock (2014) chain (ab_royer_queue.sh) builds its templates on tremor days
only.  The all-days pass (ab_alldays.py) then scans EVERY day on which a final template's stations
recorded, so the station-days of the pool outside the scan days have to exist in the same cache.

This script only READS the chain's code paths: it imports `ab_fetch.fetch_one` and
`ab_cache40.cache_one` unchanged (same FDSN providers, same one-band pick, same 1-8 Hz / 40 sps /
validity-mask cache format, same tmp+rename writes, so a station-day the chain caches at the same
time is harmless).  Nothing in the chain's files is modified.

Array lifetimes (documented rule; the matched filter itself needs >= 4 template stations per day):
  NW   2006-07-21 .. 2012-05-31   CAFE XU lifetime of the NW pool (BS11/W020/W040/PL11/W060)
  SW   2006-07-09 .. 2012-05-31   CAFE XU lifetime; the SW lines end 2008-12-31, after which only
                                  the permanent PB/UW/US and the XU PL11/W060 of the pool remain
  NVI  2005-06-01 .. 2007-08-31   POLARIS PO NVI deployment (monthly probe: VI10/VI05 data
                                  2005-06..2006-01 and 2006-05..2007-07; a later 2009-10.. segment
                                  is a different deployment and is not included)
Per station the window is intersected with its FDSN station epoch (results/fdsn_stations_40-50N_
2003-2014.txt; NVI/CN/C8 have no usable epochs there) and with a MONTHLY PROBE (10 min on the 15th,
12:00 UTC): a month is fetched when its probe or a neighbouring month's probe returned data.

  python scripts/ab_alldays_fetch.py plan                  # probe + station-day lists
  python scripts/ab_alldays_fetch.py run [--workers 8 --cache-workers 12] [--sub NW,SW,NVI]
  python scripts/ab_alldays_fetch.py run --only-days 2010-07-01,2010-09-20 --sub NW   # smoke
-> results/ab_royer/alldays/probe.csv, station_days_<SUB>.csv
-> results/ab_royer/alldays/fetch_manifest.csv, cache40_alldays.csv
-> results/ab_royer/alldays/PREFETCH_DONE
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

OUT = f"{A.RES}/alldays"
os.makedirs(OUT, exist_ok=True)
LIFETIME = {"NW": ("2006-07-21", "2012-05-31"), "SW": ("2006-07-09", "2012-05-31"),
            "NVI": ("2005-06-01", "2007-08-31")}
FMAN = f"{OUT}/fetch_manifest.csv"
CMAN = f"{OUT}/cache40_alldays.csv"


def epochs():
    p = f"{A.STUDY}/results/fdsn_stations_40-50N_2003-2014.txt"
    F = pd.read_csv(p, sep="|")
    F.columns = [c.strip("# ").strip() for c in F.columns]
    F["Station"] = F.Station.str.strip()
    F["Network"] = F.Network.str.strip()
    return F


def station_window(sta, sub, F):
    a, b = LIFETIME[sub]
    r = F[(F.Station == sta) & (F.Network.isin([A.NET[sta]] + [e[0] for e in A.EPOCHS.get(sta, [])]))]
    if len(r):
        s = str(r.StartTime.min())[:10]
        ends = r.EndTime.astype(str)
        e = "2099-12-31" if ends.str.contains("nan|None").any() or ends.str.strip().eq("").any() \
            else str(ends.max())[:10]
        a, b = max(a, s), min(b, e)
    return a, b


def probe_one(job):
    sta, month = job
    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    d = f"{month}-15"
    t = UTCDateTime(f"{d}T12:00:00")
    import ab_fetch
    for prov in ab_fetch.PROVIDERS:
        try:
            st = Client(prov, timeout=60).get_waveforms(A.net_for(sta, d), sta, "*", ab_fetch.CHAN,
                                                        t, t + 600)
            if len(st):
                return dict(sta=sta, month=month, data=True, provider=prov)
        except Exception:
            pass
    return dict(sta=sta, month=month, data=False, provider="")


def plan(subs, workers):
    F = epochs()
    pp = f"{OUT}/probe.csv"
    P = pd.read_csv(pp) if os.path.exists(pp) else pd.DataFrame(columns=["sta", "month", "data"])
    jobs = set()
    win = {}
    for sub in subs:
        for s in A.POOL[sub]:
            a, b = station_window(s, sub, F)
            win[(sub, s)] = (a, b)
            if a > b:
                continue
            for m in pd.period_range(a[:7], b[:7], freq="M").strftime("%Y-%m"):
                jobs.add((s, m))
    done = set(zip(P.sta, P.month))
    todo = sorted(jobs - done)
    A.log(f"[probe] {len(jobs)} station-months, {len(todo)} to probe")
    rows = []
    with ThreadPoolExecutor(min(workers, 8)) as ex:
        for i, r in enumerate(ex.map(probe_one, todo), 1):
            rows.append(r)
            if i % 200 == 0:
                A.log(f"  probe {i}/{len(todo)}")
    if rows:
        P = pd.concat([P, pd.DataFrame(rows)], ignore_index=True)
        P.to_csv(pp, index=False)
    P["data"] = P.data.astype(str).str.lower().eq("true")
    have = {(r.sta, r.month) for r in P.itertuples() if r.data}
    for sub in subs:
        rows = []
        for s in A.POOL[sub]:
            a, b = win[(sub, s)]
            if a > b:
                continue
            for d in pd.date_range(a, b):
                m = d.to_period("M")
                ok = any((s, str(m + k)) in have for k in (-1, 0, 1))
                if ok:
                    rows.append(dict(sta=s, date=d.strftime("%Y-%m-%d")))
        SD = pd.DataFrame(rows)
        SD.to_csv(f"{OUT}/station_days_{sub}.csv", index=False)
        A.log(f"[plan {sub}] {len(SD)} station-days, {SD.date.nunique() if len(SD) else 0} days, "
              f"{SD.sta.nunique() if len(SD) else 0} stations")


def known_nodata():
    out = set()
    for p in [f"{A.RES}/fetch_manifest.csv", f"{A.RES}/fetch_manifest_NVI.csv",
              f"{A.RES}/tremor/wecc_fetch_manifest.csv", f"{A.RES}/ext/fetch_manifest.csv", FMAN]:
        if os.path.exists(p):
            M = pd.read_csv(p, usecols=["sta", "date", "status"])
            out |= set(zip(M[M.status == "nodata"].sta, M[M.status == "nodata"].date))
    return out


def _append(path, rows):
    if not rows:
        return
    D = pd.DataFrame(rows)
    D.to_csv(path, mode="a", header=not os.path.exists(path), index=False)


def run(subs, workers, cache_workers, only_days=None):
    import ab_cache40
    import ab_fetch
    warnings.filterwarnings("ignore")
    jobs = []
    for sub in subs:
        p = f"{OUT}/station_days_{sub}.csv"
        if not os.path.exists(p):
            raise SystemExit(f"{p} missing: run `plan` first")
        SD = pd.read_csv(p)
        if only_days:
            SD = SD[SD.date.isin(only_days)]
        jobs += list(zip(SD.sta, SD.date))
    jobs = sorted(set(jobs), key=lambda x: (x[1], x[0]))
    nod = known_nodata()
    need_fetch = [(s, d) for s, d in jobs
                  if not os.path.exists(A.raw_path(s, d)) and (s, d) not in nod]
    need_cache = [(s, d) for s, d in jobs
                  if os.path.exists(A.raw_path(s, d)) and not os.path.exists(A.cache_path(s, d))]
    A.log(f"[alldays fetch] {len(jobs)} station-days; {len(need_fetch)} to fetch, "
          f"{len(need_cache)} raw on disk to cache, "
          f"{sum(os.path.exists(A.cache_path(s, d)) for s, d in jobs)} already cached, "
          f"{sum((s, d) in nod for s, d in jobs)} known nodata")
    ctx = mp.get_context("fork")
    pool = ctx.Pool(min(cache_workers, 12), maxtasksperchild=20)
    pend = [pool.apply_async(ab_cache40.cache_one, ((s, d, False),)) for s, d in need_cache]
    frows, crows, t0 = [], [], time.time()

    def drain(block=False):
        nonlocal pend
        keep = []
        for r in pend:
            if block or r.ready():
                crows.append(r.get())
            else:
                keep.append(r)
        pend = keep
        if len(crows) >= 200 or (block and crows):
            _append(CMAN, crows)
            crows.clear()

    with ThreadPoolExecutor(min(workers, 8)) as ex:
        futs = {ex.submit(ab_fetch.fetch_one, (s, d, A.raw_path(s, d))): (s, d)
                for s, d in need_fetch}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            frows.append(r)
            if r["status"] == "data":
                s, d = futs[f]
                pend.append(pool.apply_async(ab_cache40.cache_one, ((s, d, False),)))
            if i % 100 == 0 or i == len(futs):
                _append(FMAN, frows)
                frows.clear()
                drain()
                A.log(f"  fetched {i}/{len(futs)}  cache pending {len(pend)}  "
                      f"{time.time() - t0:.0f}s")
    _append(FMAN, frows)
    while pend:
        drain()
        if pend:
            time.sleep(30)
            A.log(f"  cache pending {len(pend)}")
    drain(block=True)
    pool.close()
    pool.join()
    if os.path.exists(CMAN):
        C = pd.read_csv(CMAN)
        A.log("[alldays cache] " + C.groupby("status").size().to_dict().__repr__())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "run", "all"])
    ap.add_argument("--sub", default="NW,SW,NVI")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cache-workers", type=int, default=12)
    ap.add_argument("--only-days", default="")
    a = ap.parse_args()
    subs = a.sub.split(",")
    if a.cmd in ("plan", "all"):
        plan(subs, a.workers)
    if a.cmd in ("run", "all"):
        run(subs, a.workers, a.cache_workers,
            only_days=a.only_days.split(",") if a.only_days else None)
        if a.cmd == "all" and not a.only_days:
            open(f"{OUT}/PREFETCH_DONE", "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                                   time.gmtime()) + "\n")


if __name__ == "__main__":
    main()
