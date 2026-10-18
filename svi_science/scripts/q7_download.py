#!/usr/bin/env python
"""q7_download (2026-10-07 PT) -- public EarthScope FDSN download of PBO borehole strainmeter gauges.

Level-2 processed strain (data.earthscope.org/archive/strain/bsm/processed/) redirects to a login page, so per the
user's rule (no accounts / credentials) we use the PUBLIC FDSN dataselect service instead:
  https://service.earthscope.org/fdsnws/dataselect/1/query?net=PB&sta=<STA>&cha=RS*,RDO&start=<Y>-01-01&end=<Y+1>-01-01
  RS1-RS4 = the four gauge readings at 10-min sampling (raw counts), RDO = barometric pressure (30-min).
Output: data/bsm/raw/<STA>/<STA>.<YEAR>.mseed ; a manifest data/bsm/MANIFEST.csv (url, bytes, http code).
"""
import os, subprocess, csv, time
ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/data/bsm/raw"
STAS = ["B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B012", "B013", "B014", "B017", "B018",
        "B019", "B020", "B201", "B926", "B927", "B928", "B941", "B943"]
YEARS = range(2009, 2027)
import sys
if len(sys.argv) > 1:
    STAS = sys.argv[1].split(",")
man = f"{ROOT}/data/bsm/MANIFEST.csv"
new = not os.path.exists(man)
with open(man, "a", newline="") as fh:
    w = csv.writer(fh)
    if new:
        w.writerow(["station", "year", "url", "http", "bytes", "fetched_utc"])
    for s in STAS:
        os.makedirs(f"{OUT}/{s}", exist_ok=True)
        for y in YEARS:
            f = f"{OUT}/{s}/{s}.{y}.mseed"
            if os.path.exists(f) and os.path.getsize(f) > 0:
                continue
            end = f"{y+1}-01-01T00:00:00" if y < 2026 else "2026-10-07T00:00:00"
            url = (f"https://service.earthscope.org/fdsnws/dataselect/1/query?net=PB&sta={s}&cha=RS*,RDO"
                   f"&start={y}-01-01T00:00:00&end={end}&nodata=404")
            r = subprocess.run(["curl", "-sS", "--max-time", "600", "-o", f, "-w", "%{http_code}", url],
                               capture_output=True, text=True)
            code = r.stdout.strip()
            size = os.path.getsize(f) if os.path.exists(f) else 0
            if code != "200" and os.path.exists(f):
                os.rename(f, f + f".http{code}")
            w.writerow([s, y, url, code, size, time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())]); fh.flush()
            print(s, y, code, size, flush=True)
