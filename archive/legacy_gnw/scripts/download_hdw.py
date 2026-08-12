"""Download UW.HDW vertical-component daily mseed (1995->now) from IRIS FDSN into the
data/waveforms/UW.HDW/<year>/<jday>.mseed layout the pipeline expects.

Robust + resumable: skips days already on disk, retries once, picks one preferred vertical
channel per day (EHZ>HHZ>BHZ>SHZ), logs progress. Parallel via threads (FDSN is I/O-bound).
"""
from __future__ import annotations
import os, sys, time, warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
warnings.filterwarnings("ignore")
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException

NET, STA = "UW", "HDW"
OUT = Path("data/waveforms/UW.HDW")
PREF = ["EHZ", "HHZ", "BHZ", "SHZ"]
START = datetime(1995, 1, 1)
END = datetime.now(timezone.utc).replace(tzinfo=None)  # today
WORKERS = 8

def day_path(d): return OUT / f"{d.year}" / f"{d.timetuple().tm_yday:03d}.mseed"

def fetch(d, client):
    p = day_path(d)
    if p.exists() and p.stat().st_size > 0:
        return "skip"
    t0 = UTCDateTime(d); t1 = t0 + 86400
    for attempt in (1, 2):
        try:
            st = client.get_waveforms(NET, STA, "*", "?HZ", t0, t1)
            break
        except FDSNNoDataException:
            return "nodata"
        except Exception:
            if attempt == 2:
                return "fail"
            time.sleep(2)
    codes = {tr.stats.channel for tr in st}
    pick = next((c for c in PREF if c in codes), None)
    if pick is None:
        return "nodata"
    st = st.select(channel=pick)
    if len(st) == 0:
        return "nodata"
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        st.write(str(p), format="MSEED")
    except Exception:
        return "fail"
    return "ok"

def main():
    days = []
    d = START
    while d <= END:
        days.append(d); d += timedelta(days=1)
    print(f"[hdw-dl] {len(days)} days {START.date()}..{END.date()}, {WORKERS} workers", flush=True)
    client = Client("IRIS", timeout=60)
    counts = {"ok": 0, "skip": 0, "nodata": 0, "fail": 0}
    done = 0; t0 = time.time(); fails = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch, dd, client): dd for dd in days}
        for f in as_completed(futs):
            r = f.result(); counts[r] += 1; done += 1
            if r == "fail":
                fails.append(futs[f])
            if done % 500 == 0:
                el = time.time() - t0
                print(f"[hdw-dl] {done}/{len(days)} ({done/len(days)*100:.0f}%) "
                      f"ok={counts['ok']} skip={counts['skip']} nodata={counts['nodata']} "
                      f"fail={counts['fail']} | {done/el:.1f} day/s | {el/60:.0f} min", flush=True)
    print(f"[hdw-dl] DONE: {counts}", flush=True)
    if fails:
        print(f"[hdw-dl] {len(fails)} failed days (re-run to retry): {[str(x.date()) for x in fails[:10]]}", flush=True)

if __name__ == "__main__":
    main()
