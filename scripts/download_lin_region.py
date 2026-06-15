"""Download 3-component waveforms for the Lin-region validation set.

Targets (in Lin et al. 2023 southern-VI box, co-located with the LFE catalog):
  CN.PGC  (Pat Bay broadband, BH? 40 Hz)        2005-2017
  PB.B011 (Pat Bay borehole,  EH? 100 Hz)       2007-2017  -- already graded (38 GOLD)
  PB.B926 (VI borehole,       EH? 100 Hz)       2007-2017

Window = Lin catalog overlap (catalog ends 2017-02-21). 3-component on purpose:
keeps the multi-component / network-moveout option open and gives richer LFE
features at PGC. Writes data/waveforms/{net}.{sta}/{year}/{jday}.mseed.
"""
import sys, time, logging
from pathlib import Path
from datetime import datetime

sys.path.insert(0, "src")
from tremorferometry.waveforms import fetch_many

logging.basicConfig(level=logging.ERROR, format="%(asctime)s %(message)s")
root = Path("data/waveforms")
END = datetime(2017, 3, 1)
JOBS = [
    ("CN", "PGC",  datetime(2005, 1, 1), "BH?"),
    ("PB", "B011", datetime(2007, 1, 1), "EH?"),
    ("PB", "B926", datetime(2007, 1, 1), "EH?"),
]
for net, sta, start, chan in JOBS:
    t0 = time.time()
    print(f"=== {net}.{sta} {start.date()}..{END.date()} chan={chan} ===", flush=True)
    paths = fetch_many(root, net, [sta], start, END,
                       provider="IRIS", max_workers=16, channel=chan)
    print(f"    {net}.{sta}: {len(paths)} day-files in {time.time()-t0:.0f}s", flush=True)
print("ALL DOWNLOADS DONE", flush=True)
