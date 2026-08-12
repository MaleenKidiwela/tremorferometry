#!/usr/bin/env python
"""Download a broadband/land station's daily 3-comp mseed (BH?/HH?, whatever band exists per day) for a
date range, from the EarthScope FDSN archive (retry-wrapped for the flappy migration). Decimation to 40 Hz
happens downstream in the pipeline. Resumable (skips existing days). Spawn-safe (net/sta travel in the item).
Usage: python scripts/download_broadband.py --net CN --sta CLRS --start 2007-01-01 --end 2026-08-01 [--workers 6]"""
import os, argparse
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# broadband (BH 40Hz / HH 100Hz) AND short-period land stations (EH?/SH?, often EHZ-only). The coda dv/v is
# Z-only, so a single vertical (EHZ) is a valid fleet receiver; fetch whatever vertical-bearing band exists.
CHAN = "BH?,HH?,EH?,SH?"


def _client(net):
    import time
    from obspy.clients.fdsn import Client
    # BK (Berkeley) + NC (USGS N. California) live at NCEDC, NOT IRIS/EarthScope; everyone else at IRIS.
    provs = (["NCEDC", "https://service.ncedc.org"] if net in ("BK", "NC")
             else ["IRIS", "https://service.iris.edu", "EARTHSCOPE"])
    for prov in provs:
        for _ in range(3):
            try:
                return Client(prov, timeout=60)
            except Exception:
                time.sleep(3)
    raise RuntimeError("no FDSN service reachable")


def fetch_day(item):
    import time
    from obspy import UTCDateTime
    net, sta, chan, y, m, d = item
    day = datetime(y, m, d)
    p = f"data/waveforms/{net}.{sta}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"
    if os.path.exists(p) or os.path.exists(p + ".nd"):   # .nd = cached no-data day -> skip on re-download
        return "skip"
    t0 = UTCDateTime(day)
    for _ in range(3):
        try:
            st = _client(net).get_waveforms(net, sta, "*", chan, t0, t0 + 86400)
            if not st:
                os.makedirs(os.path.dirname(p), exist_ok=True); open(p + ".nd", "w").close()
                return "nodata"
            os.makedirs(os.path.dirname(p), exist_ok=True)
            st.write(p + ".tmp", format="MSEED"); os.replace(p + ".tmp", p)
            return "ok"
        except Exception as e:
            if "No data" in str(e) or "204" in str(e):
                os.makedirs(os.path.dirname(p), exist_ok=True); open(p + ".nd", "w").close()
                return "nodata"
            time.sleep(5)
    return "fail"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--chan", default=CHAN, help="channel selector; restrict to ONE band for multi-band "
                    "stations (e.g. BH? for SHB) so the downstream *Z pick is unambiguous")
    a = ap.parse_args()
    d0 = datetime.fromisoformat(a.start); d1 = datetime.fromisoformat(a.end)
    items = []; d = d0
    while d < d1:
        items.append((a.net, a.sta, a.chan, d.year, d.month, d.day)); d += timedelta(days=1)
    print(f"[dl {a.net}.{a.sta}] {len(items)} days {a.start}..{a.end}", flush=True)
    counts = {}; done = 0
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(fetch_day, items, chunksize=4):
            counts[r] = counts.get(r, 0) + 1; done += 1
            if done % 300 == 0:
                print(f"  {done}/{len(items)} {counts}", flush=True)
    print(f"[dl {a.net}.{a.sta}] DONE: {counts}")


if __name__ == "__main__":
    main()
