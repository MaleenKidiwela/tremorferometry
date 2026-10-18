#!/usr/bin/env python
"""D3 stage 1 — 3-component continuous download into lfe_catalog_repair/traces/.

Same logic as scripts/download_broadband.py (retry-wrapped FDSN, resumable, `.nd` no-data markers,
spawn-safe) but with a CONFIGURABLE output root so nothing is written outside this study folder.

Usage:
  python lfe_catalog_repair/scripts/d3_download.py --net UW --sta WISH --start 2017-01-01 --end 2019-01-01 \
      --chan 'BH?' --outroot lfe_catalog_repair/traces --workers 6
Layout: <outroot>/<NET>.<STA>/<YYYY>/<JJJ>.mseed   (identical to data/waveforms, so the project's
        readers work unchanged with --wfdir <outroot>)
"""
import os
import argparse
import multiprocessing as mp
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor


def _client(net):
    import time
    from obspy.clients.fdsn import Client
    provs = (["NCEDC", "https://service.ncedc.org"] if net in ("BK", "NC")
             else ["IRIS", "https://service.iris.edu", "EARTHSCOPE"])
    for prov in provs:
        for _ in range(3):
            try:
                return Client(prov, timeout=90)
            except Exception:
                time.sleep(3)
    raise RuntimeError("no FDSN service reachable")


def fetch_day(item):
    import time
    from obspy import UTCDateTime
    net, sta, chan, y, m, d, outroot = item
    day = datetime(y, m, d)
    p = f"{outroot}/{net}.{sta}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"
    if os.path.exists(p) or os.path.exists(p + ".nd"):
        return "skip"
    t0 = UTCDateTime(day)
    for _ in range(4):
        try:
            st = _client(net).get_waveforms(net, sta, "*", chan, t0, t0 + 86400)
            if not st:
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p + ".nd", "w").close()
                return "nodata"
            os.makedirs(os.path.dirname(p), exist_ok=True)
            st.write(p + ".tmp", format="MSEED")
            os.replace(p + ".tmp", p)
            return "ok"
        except Exception as e:
            if "No data" in str(e) or "204" in str(e):
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p + ".nd", "w").close()
                return "nodata"
            time.sleep(5)
    return "fail"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--chan", required=True, help="e.g. 'BH?' / 'HH?' / 'EH?' (ONE band -> unambiguous 3-C)")
    ap.add_argument("--outroot", default="lfe_catalog_repair/traces")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    outroot = os.path.abspath(a.outroot)
    d0 = datetime.fromisoformat(a.start)
    d1 = datetime.fromisoformat(a.end)
    items = []
    d = d0
    while d < d1:
        items.append((a.net, a.sta, a.chan, d.year, d.month, d.day, outroot))
        d += timedelta(days=1)
    print(f"[d3-dl {a.net}.{a.sta}] {len(items)} days {a.start}..{a.end} chan={a.chan} -> {outroot}", flush=True)
    counts, done = {}, 0
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(fetch_day, items, chunksize=2):
            counts[r] = counts.get(r, 0) + 1
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(items)} {counts}", flush=True)
    print(f"[d3-dl {a.net}.{a.sta}] DONE: {counts}", flush=True)


if __name__ == "__main__":
    main()
