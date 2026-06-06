"""Download one station's vertical-component daily mseed from FDSN into the
data/waveforms/<NET>.<STA>/<year>/<jday>.mseed layout the pipeline expects.

Generic + resumable: skips days already on disk, retries once, picks one preferred
vertical channel per day (EHZ>HHZ>BHZ>SHZ), threaded. Reusable for any station.

  python scripts/download_station.py --network UW --station CPW \
      --start 1996-01-01 --end 2022-08-01 --workers 8
"""
from __future__ import annotations
import argparse
import time
import warnings
from datetime import datetime, timedelta
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

warnings.filterwarnings("ignore")
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException

PREF = ["EHZ", "HHZ", "BHZ", "SHZ"]


def fetch(d, client, net, sta, out):
    p = out / f"{d.year}" / f"{d.timetuple().tm_yday:03d}.mseed"
    if p.exists() and p.stat().st_size > 0:
        return "skip"
    t0 = UTCDateTime(d); t1 = t0 + 86400
    for attempt in (1, 2):
        try:
            st = client.get_waveforms(net, sta, "*", "?HZ", t0, t1)
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--network", required=True)
    ap.add_argument("--station", required=True)
    ap.add_argument("--start", required=True, help="YYYY-MM-DD")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--client", default="IRIS")
    args = ap.parse_args()

    out = Path(f"data/waveforms/{args.network}.{args.station}")
    start = datetime.fromisoformat(args.start)
    end = datetime.fromisoformat(args.end)
    days = []
    d = start
    while d <= end:
        days.append(d); d += timedelta(days=1)
    print(f"[dl] {args.network}.{args.station} {len(days)} days "
          f"{start.date()}..{end.date()}, {args.workers} workers", flush=True)
    client = Client(args.client, timeout=60)
    counts = {"ok": 0, "skip": 0, "nodata": 0, "fail": 0}
    done = 0; t0 = time.time(); fails = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch, dd, client, args.network, args.station, out): dd
                for dd in days}
        for f in as_completed(futs):
            r = f.result(); counts[r] += 1; done += 1
            if r == "fail":
                fails.append(futs[f])
            if done % 500 == 0:
                el = time.time() - t0
                print(f"[dl] {done}/{len(days)} ({done/len(days)*100:.0f}%) "
                      f"ok={counts['ok']} skip={counts['skip']} nodata={counts['nodata']} "
                      f"fail={counts['fail']} | {done/el:.1f} day/s | {el/60:.0f} min", flush=True)
    print(f"[dl] DONE: {counts}", flush=True)
    if fails:
        print(f"[dl] {len(fails)} failed days (re-run to retry)", flush=True)


if __name__ == "__main__":
    main()
