#!/usr/bin/env python
"""Download CN.PGC broadband 3-comp (BH?) daily mseed for a date range. Uses the working EarthScope FDSN
URL (the 'IRIS' short name currently 404s at service.earthscope.org). Resumable (skips existing days).
Usage: python scripts/download_pgc_broadband.py --start 2005-01-01 --end 2009-01-01 [--workers 8]"""
import os, argparse
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

import argparse as _ap0
_A = None                                # set in main(); workers read module globals after spawn re-import
NET, STA = "CN", "PGC"
CHAN = "BH?,HH?"                          # fetch whichever band exists that day (BH 40Hz / HH 100Hz)


def _root():
    return f"data/waveforms/{NET}.{STA}"


def _client():
    import time
    from obspy.clients.fdsn import Client
    # the flappy EarthScope migration: short name 'IRIS' + explicit URLs rotate in/out; retry the discovery
    for prov in ["IRIS", "https://service.iris.edu", "EARTHSCOPE"]:
        for _ in range(3):
            try:
                return Client(prov, timeout=60)
            except Exception:
                time.sleep(3)
    raise RuntimeError("no FDSN service reachable")


def fetch_day(day):
    import time
    from obspy import UTCDateTime
    p = f"{ROOT}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"
    if os.path.exists(p):
        return "skip"
    t0 = UTCDateTime(day)
    for attempt in range(3):
        try:
            c = _client()
            st = c.get_waveforms("CN", "PGC", "*", CHAN, t0, t0 + 86400)
            if not st:
                return "nodata"
            os.makedirs(os.path.dirname(p), exist_ok=True)
            st.write(p + ".tmp", format="MSEED"); os.replace(p + ".tmp", p)
            return "ok"
        except Exception as e:
            if "No data" in str(e) or "204" in str(e):
                return "nodata"
            time.sleep(5)
    return "fail"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    d0 = datetime.fromisoformat(a.start); d1 = datetime.fromisoformat(a.end)
    days = []; d = d0
    while d < d1:
        days.append(d); d += timedelta(days=1)
    print(f"[pgc-dl] {len(days)} days {a.start}..{a.end}", flush=True)
    counts = {}; done = 0
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(fetch_day, days, chunksize=4):
            counts[r] = counts.get(r, 0) + 1; done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(days)} {counts}", flush=True)
    print(f"[pgc-dl] DONE {a.start}..{a.end}: {counts}")


if __name__ == "__main__":
    main()
