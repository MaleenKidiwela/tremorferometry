#!/usr/bin/env python
"""Download B011 (PB borehole) 3-component EH1/EH2/EHZ daily mseed for 2018-2026,
into data/waveforms/PB.B011/{year}/{jday}.mseed (same layout as 2007-2017).
Resumable: skips days that already have all 3 channels. Keeps ALL EH channels
(unlike download_station.py which picks one vertical)."""
import sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client

NET, STA = "PB", "B011"
ROOT = Path(f"data/waveforms/{NET}.{STA}")
WANT = {"EH1", "EH2", "EHZ"}


def one(day, client):
    p = ROOT / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if p.exists():
        try:
            have = {tr.stats.channel for tr in read(str(p))}
            if WANT.issubset(have):
                return "skip"
        except Exception:
            pass
    t0 = UTCDateTime(day); t1 = t0 + 86400
    for attempt in (1, 2):
        try:
            st = client.get_waveforms(NET, STA, "*", "EH?", t0, t1)
            break
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
    start = datetime.fromisoformat("2018-01-01")
    end = datetime.fromisoformat("2026-05-31")
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    days = []
    d = start
    while d <= end:
        days.append(d); d += timedelta(days=1)
    print(f"[dl] {NET}.{STA} 3-comp {len(days)} days {start.date()}..{end.date()}, "
          f"{workers} workers", flush=True)
    client = Client("IRIS", timeout=60)
    counts = {"ok": 0, "skip": 0, "nodata": 0, "fail": 0}
    done = 0; t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, dd, client): dd for dd in days}
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
