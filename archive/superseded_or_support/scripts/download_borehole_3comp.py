#!/usr/bin/env python
"""Download a PB-style borehole's 3-component EH1/EH2/EHZ daily mseed into
data/waveforms/{NET}.{STA}/{year}/{jday}.mseed (same layout as existing data).
Resumable: skips days that already have all 3 channels. Keeps ALL EH channels.

Usage: python scripts/download_borehole_3comp.py --net PB --sta B926 \
         --start 2018-01-01 --end 2026-05-31 --workers 10
"""
import argparse, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client

WANT = {"EH1", "EH2", "EHZ"}


def one(day, client, net, sta, root):
    p = root / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if p.exists():
        try:
            if WANT.issubset({tr.stats.channel for tr in read(str(p))}):
                return "skip"
        except Exception:
            pass
    t0 = UTCDateTime(day); t1 = t0 + 86400
    for attempt in (1, 2):
        try:
            st = client.get_waveforms(net, sta, "*", "EH?", t0, t1); break
        except Exception:
            if attempt == 2:
                return "nodata"
            time.sleep(2)
    st = st.select(channel="EH?")
    if len(st) == 0:
        return "nodata"
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        st.write(str(p), format="MSEED")
    except Exception:
        return "fail"
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    root = Path(f"data/waveforms/{a.net}.{a.sta}")
    start = datetime.fromisoformat(a.start); end = datetime.fromisoformat(a.end)
    days = []; d = start
    while d <= end:
        days.append(d); d += timedelta(days=1)
    print(f"[dl] {a.net}.{a.sta} 3-comp {len(days)} days {start.date()}..{end.date()}, "
          f"{a.workers} workers", flush=True)
    client = Client("IRIS", timeout=60)
    counts = {"ok": 0, "skip": 0, "nodata": 0, "fail": 0}; done = 0; t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one, dd, client, a.net, a.sta, root): dd for dd in days}
        for f in as_completed(futs):
            counts[f.result()] += 1; done += 1
            if done % 200 == 0:
                el = time.time() - t0
                print(f"[dl] {done}/{len(days)} ({done/len(days)*100:.0f}%) "
                      f"ok={counts['ok']} skip={counts['skip']} nodata={counts['nodata']} "
                      f"fail={counts['fail']} | {done/el:.1f} day/s | {el/60:.0f} min", flush=True)
    print(f"[dl] DONE: {counts}", flush=True)


if __name__ == "__main__":
    main()
