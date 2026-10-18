#!/usr/bin/env python
"""Tranche 6 - intraslab earthquake windows for the station corrections.

(a) USGS ComCat (https://earthquake.usgs.gov/fdsnws/event/1/query): M >= 2.0, depth 25-50 km,
    inside the box, 2006-06-01..2012-06-30 -> results/intraslab_eq_2006_2012.csv
(b) +- 3 min around every origin, for every 3-C station alive on that day (all networks of
    tranche 1) -> traces_eq/<NET.STA>/<EVID>.mseed, manifest rows tranche=T6 in
    results/download_manifest.csv.

  python scripts/d_eq.py --catalog-only
  python scripts/d_eq.py
"""
from __future__ import annotations
import os
import sys
import csv
import time
import argparse
import threading
import warnings
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "lfe_picker", "scripts")))
import d_common as D                                                   # noqa: E402
import d_queues as Q                                                   # noqa: E402
from r_download_day import get_client                                  # noqa: E402

CAT = os.path.join(D.RES, "intraslab_eq_2006_2012.csv")
PAD_S = 180.0
MMIN, DEP = 2.0, (25.0, 50.0)
COLS = ["tranche", "net", "sta", "day", "status", "bytes", "chain", "linked", "loc", "band",
        "fs", "three_c", "n_channels", "coverage_frac", "n_frag", "outroot", "group",
        "updated_utc", "note"]
_lock = threading.Lock()
_n = {"done": 0, "data": 0, "nodata": 0, "error": 0, "bytes": 0}
_t0 = time.time()


def fetch_catalog(force=False):
    if os.path.exists(CAT) and not force:
        E = pd.read_csv(CAT)
        D.log(f"[eq] catalog cached: {len(E)} events -> {CAT}")
        return E
    import urllib.request
    import urllib.parse
    t0, t1 = D.WINDOWS["t6_intraslab"]
    p = dict(format="csv", starttime=t0, endtime=t1, minmagnitude=MMIN,
             mindepth=DEP[0], maxdepth=DEP[1], minlatitude=D.LAT0, maxlatitude=D.LAT1,
             minlongitude=D.LON0, maxlongitude=D.LON1, orderby="time-asc")
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urllib.parse.urlencode(p)
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                raw = r.read().decode()
            break
        except Exception as e:
            D.log(f"[eq] ComCat try {k+1}: {type(e).__name__}: {str(e)[:120]}")
            time.sleep(5 * (k + 1))
    else:
        raise RuntimeError("ComCat unreachable")
    E = pd.read_csv(pd.io.common.StringIO(raw))
    E = E.rename(columns={"latitude": "lat", "longitude": "lon", "depth": "depth_km"})
    E["evid"] = E.id.astype(str)
    E["origin_utc"] = pd.to_datetime(E.time, utc=True, format="mixed").astype("datetime64[ns, UTC]")
    E = E[["evid", "origin_utc", "lat", "lon", "depth_km", "mag", "magType", "place"]]
    E.to_csv(CAT, index=False)
    D.log(f"[eq] {len(E)} events M>={MMIN} {DEP[0]}-{DEP[1]} km in box -> {CAT}")
    return E


def one(job):
    from obspy import UTCDateTime
    net, sta, loc, band, fs, chain, three_c, evid, t0 = job
    dst = os.path.join(D.TRACES_EQ, f"{net}.{sta}", f"{evid}.mseed")
    rec = dict(tranche="T6", net=net, sta=sta, day=evid, status="", bytes=0, chain=chain,
               linked="", loc=loc, band=band, fs=fs, three_c=three_c, n_channels="",
               coverage_frac="", n_frag="", outroot=D.TRACES_EQ, group="intraslab",
               updated_utc="", note="")
    if os.path.exists(dst):
        rec.update(status="data", bytes=os.path.getsize(dst), note="resumed")
    elif os.path.exists(dst + ".nd"):
        rec.update(status="nodata", note="resumed")
    else:
        t = UTCDateTime(t0)
        last = ""
        for k in range(3):
            try:
                st = get_client("IRIS", 120).get_waveforms(
                    net, sta, loc if loc else "*", band + "?", t - PAD_S, t + PAD_S)
                if loc:
                    st = st.select(location=loc)
                st = st.select(channel=band + "?")
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                if len(st) == 0:
                    open(dst + ".nd", "w").close()
                    rec.update(status="nodata")
                else:
                    st.write(dst + ".tmp", format="MSEED")
                    os.replace(dst + ".tmp", dst)
                    rec.update(status="data", bytes=os.path.getsize(dst),
                               n_channels=len({tr.stats.channel for tr in st}))
                break
            except Exception as e:
                last = f"{type(e).__name__}:{str(e)[:100]}"
                if "No data" in str(e) or "204" in str(e) or "FDSNNoData" in last:
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    open(dst + ".nd", "w").close()
                    rec.update(status="nodata")
                    break
                time.sleep(3 * (k + 1) ** 2)
        else:
            rec.update(status="error", note=last)
    rec["updated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with _lock:
        _n["done"] += 1
        _n[rec["status"]] = _n.get(rec["status"], 0) + 1
        _n["bytes"] += int(rec["bytes"] or 0)
        with open(D.MANIFEST, "a", newline="") as f:
            csv.DictWriter(f, COLS, extrasaction="ignore").writerow(rec)
        if _n["done"] % 500 == 0:
            D.log(f"  [eq] {_n['done']} done data={_n.get('data',0)} "
                  f"nodata={_n.get('nodata',0)} err={_n.get('error',0)} "
                  f"{_n['bytes']/2**30:.2f} GiB {(time.time()-_t0)/60:.1f} min")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog-only", action="store_true")
    ap.add_argument("--workers", type=int, default=D.MAX_CONCURRENT_FDSN)
    ap.add_argument("--force-catalog", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    E = fetch_catalog(a.force_catalog)
    if a.catalog_only or not len(E):
        return
    F = D.station_flags()
    C = pd.read_csv(os.path.join(D.RES, "channels_all.csv"), keep_default_na=False,
                    dtype={"loc": str})
    C["loc"] = C["loc"].astype(str).replace("nan", "")
    C["fs"] = pd.to_numeric(C.fs, errors="coerce").fillna(0.0)
    stations = [(r.net, r.sta) for r in F[F.three_c].itertuples()]
    E["origin_utc"] = pd.to_datetime(E.origin_utc, utc=True)
    E["day"] = E.origin_utc.dt.strftime("%Y-%m-%d")
    dlist = sorted(set(E.day))
    jobs = []
    for net, sta in stations:
        r = Q.resolve_days(C, net, sta, dlist)
        if not r:
            continue
        for e in E.itertuples():
            vs = r.get(e.day)
            if not vs or not vs[0][4]:
                continue
            v = vs[0]
            jobs.append((net, sta, v[0], v[1], v[2], v[3], True, e.evid,
                         e.origin_utc.tz_convert(None).isoformat()))
    done = set()
    if os.path.exists(D.MANIFEST):
        M = pd.read_csv(D.MANIFEST, keep_default_na=False)
        M = M[(M.tranche == "T6") & (M.status.isin(["data", "nodata"]))]
        done = set(zip(M.net, M.sta, M.day.astype(str)))
    else:
        with open(D.MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
    todo = [j for j in jobs if (j[0], j[1], j[7]) not in done]
    D.log(f"[eq] {len(E)} events x 3-C stations alive = {len(jobs)} windows, {len(todo)} to fetch")
    with ThreadPoolExecutor(max_workers=min(a.workers, D.MAX_CONCURRENT_FDSN)) as ex:
        list(ex.map(one, todo))
    D.log(f"[eq] finished: {_n} in {(time.time()-_t0)/60:.1f} min")


if __name__ == "__main__":
    main()
