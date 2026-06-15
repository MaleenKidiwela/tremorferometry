"""Append EH1/EH2 horizontals to B927's existing EHZ-only daily mseed files.

B927 (WA borehole) currently has EHZ on disk for 2008-2026. We add the two
horizontals (for polarization features) WITHOUT re-downloading EHZ: read the
existing day file, fetch only the missing channels, write the merged stream back.
"""
import sys, time
from pathlib import Path
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed

from obspy import read, UTCDateTime
from obspy.clients.fdsn import Client

cl = Client("IRIS")
ROOT = Path("data/waveforms/PB.B927")
START, END = datetime(2008, 1, 1), datetime(2026, 6, 7)


def one(day):
    out = ROOT / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if not out.exists():
        return "nozday"                      # only complete days we already have Z for
    try:
        existing = read(str(out))
    except Exception:
        return "readerr"
    have = {tr.stats.channel for tr in existing}
    need = [c for c in ["EH1", "EH2"] if c not in have]
    if not need:
        return "skip"
    t0 = UTCDateTime(day)
    try:
        st = cl.get_waveforms("PB", "B927", "*", ",".join(need), t0, t0 + 86400)
    except Exception:
        return "nodata"
    if len(st) == 0:
        return "nodata"
    (existing + st).write(str(out), format="MSEED")
    return "done"


def main():
    days = []
    d = START
    while d < END:
        days.append(d); d += timedelta(days=1)
    t0 = time.time()
    counts = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(one, dd): dd for dd in days}
        for i, f in enumerate(as_completed(futs)):
            r = f.result()
            counts[r] = counts.get(r, 0) + 1
            if (i + 1) % 500 == 0:
                print(f"  {i+1}/{len(days)}  {counts}  ({time.time()-t0:.0f}s)", flush=True)
    print(f"B927 horizontals: {counts}  total {time.time()-t0:.0f}s")
    print("B927 HORIZONTALS DONE", flush=True)


if __name__ == "__main__":
    main()
