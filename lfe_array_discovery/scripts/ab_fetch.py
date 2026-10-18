#!/usr/bin/env python
"""STAGE 1 of the Royer queue: download the station-days of Royer & Bostock's 14 named stations
that are not on disk as raw miniSEED.  Prefix `ab_`.

Why raw and not the 20 Hz cache: the paper resamples to 40 sps (sec.3) and the existing cache is
20 Hz, so every station-day has to be re-cached from the native rate.  638 of the 833 needed
station-days still have their raw miniSEED under `traces/`; the rest were freed and are re-fetched
here.

  python scripts/ab_fetch.py --plan          # print the list and the volume, download nothing
  python scripts/ab_fetch.py --workers 8     # fetch (<= 8 concurrent FDSN connections)

-> traces/<NET.STA>/<YYYY>/<DOY>.mseed
-> results/ab_royer/fetch_manifest.csv   one row per station-day: status, bytes, fs, channels, err
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

# 3-C band preference.  PBO boreholes are EH?; CAFE (XU) and PNSN short period are BH?/EH?/HH?.
# The channel list is resolved per station-day from the FDSN answer, not assumed.
CHAN = "BH?,HH?,EH?,SH?"
PROVIDERS = ["IRIS", "EARTHSCOPE", "https://earthquakescanada.nrcan.gc.ca"]   # CN data: NRCan
MANIFEST = f"{A.RES}/fetch_manifest.csv"


def wanted(ext=False):
    """[(sta, date, subarrays)] for every station-day the chain reads, with what is on disk.
    ext=True: the EXTENSION station-days (results/ab_royer/ext/station_days_<SUB>.csv)."""
    need = {}
    for sub in ("NW", "SW", "NVI"):
        if ext:
            SD = pd.read_csv(f"{A.EXT_DIR}/station_days_{sub}.csv")
            pairs = list(zip(SD.sta, SD.date))
        else:
            pairs = A.station_days(sub)
        for s, d in pairs:
            need.setdefault((s, d), set()).add(sub)
    rows = []
    for (s, d), subs in sorted(need.items()):
        p = A.raw_path(s, d)
        have = os.path.exists(p)
        rows.append(dict(sta=s, netsta=A.netsta(s), date=d, subarray="+".join(sorted(subs)),
                         raw_path=p, on_disk=have,
                         bytes_on_disk=(os.path.getsize(p) if have else 0)))
    return pd.DataFrame(rows)


def fetch_one(job):
    sta, date, out = job
    t0 = time.time()
    try:
        from obspy import UTCDateTime
        from obspy.clients.fdsn import Client
        t = UTCDateTime(date)
        last = None
        for prov in PROVIDERS:
            try:
                c = Client(prov, timeout=180)
                st = c.get_waveforms(A.net_for(sta, date), sta, "*", CHAN, t, t + 86400,
                                     attach_response=False)
                if len(st) == 0:
                    last = "empty stream"
                    continue
                st, _band = A.pick_band(st)
                os.makedirs(os.path.dirname(out), exist_ok=True)
                tmp = out + f".tmp{os.getpid()}"
                st.write(tmp, format="MSEED")
                os.replace(tmp, out)
                chans = sorted({tr.stats.channel for tr in st})
                fs = sorted({float(tr.stats.sampling_rate) for tr in st})
                cov = sum(float(tr.stats.npts) / float(tr.stats.sampling_rate)
                          for tr in st) / (86400.0 * max(len(chans), 1))
                return dict(sta=sta, date=date, status="data", provider=prov,
                            bytes=os.path.getsize(out), chans=",".join(chans),
                            fs=",".join(f"{v:g}" for v in fs), coverage_frac=round(cov, 4),
                            secs=round(time.time() - t0, 1), err="")
            except Exception as e:                       # try the next provider, then report
                last = f"{type(e).__name__}: {e}"
        return dict(sta=sta, date=date, status="nodata", provider="", bytes=0, chans="", fs="",
                    coverage_frac=0.0, secs=round(time.time() - t0, 1), err=str(last)[:300])
    except Exception as e:
        return dict(sta=sta, date=date, status="error", provider="", bytes=0, chans="", fs="",
                    coverage_frac=0.0, secs=round(time.time() - t0, 1),
                    err=f"{type(e).__name__}: {e}"[:300])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8, help="<= 8 concurrent FDSN connections")
    ap.add_argument("--plan", action="store_true", help="report the gap and exit")
    ap.add_argument("--only", default="", help="restrict to these station codes")
    ap.add_argument("--ext", action="store_true", help="the EXTENSION station-days")
    a = ap.parse_args()
    global MANIFEST
    if a.ext:
        MANIFEST = f"{A.EXT_DIR}/fetch_manifest.csv"
    W = wanted(a.ext)
    if a.only:
        W = W[W.sta.isin(a.only.split(","))]
    miss = W[~W.on_disk].copy()
    have = W[W.on_disk]
    mb = have.bytes_on_disk.sum() / 1e6
    per_day = (mb / max(len(have), 1))
    A.log(f"Royer station-days: {len(W)} needed, {len(have)} raw on disk "
          f"({mb/1e3:.1f} GB, mean {per_day:.0f} MB/day), {len(miss)} to fetch "
          f"(~{len(miss)*per_day/1e3:.1f} GB at that mean)")
    if len(miss):
        A.log("to fetch, by station:\n" +
              miss.groupby(["netsta", "subarray"]).size().to_string())
    W.to_csv(f"{A.EXT_DIR if a.ext else A.RES}/fetch_plan.csv", index=False)
    if a.plan or len(miss) == 0:
        return
    jobs = [(r.sta, r.date, r.raw_path) for r in miss.itertuples()]
    rows = []
    with ThreadPoolExecutor(max_workers=min(a.workers, 8)) as ex:
        for i, r in enumerate(ex.map(fetch_one, jobs), 1):
            rows.append(r)
            if i % 10 == 0 or i == len(jobs):
                got = sum(x["bytes"] for x in rows) / 1e6
                A.log(f"  {i}/{len(jobs)}  {got:.0f} MB")
    M = pd.DataFrame(rows)
    if os.path.exists(MANIFEST):
        M = pd.concat([pd.read_csv(MANIFEST), M], ignore_index=True) \
              .drop_duplicates(["sta", "date"], keep="last")
    M.to_csv(MANIFEST, index=False)
    A.log(M.groupby("status").agg(n=("sta", "size"), MB=("bytes", lambda v: round(v.sum()/1e6)))
           .to_string())
    for r in M[M.status != "data"].head(20).itertuples():
        A.log(f"  MISS {r.sta} {r.date}: {r.err}")


if __name__ == "__main__":
    main()
