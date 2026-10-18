#!/usr/bin/env python
"""Re-fetch the NVI station-days that ab_fetch.py recorded as `nodata` (2026-09-27 NVI diagnosis).

ab_fetch.py asks for 24 h in one request.  NRCan answers 404 to some day-long requests that it serves
in shorter pieces (CN.MAYB 2006-09-07 / 2007-06-15 hour 5 returned EHZ data), and EarthScope has
PO.VI56 2006-09-07 that the manifest calls nodata.  Here every nodata station-day is asked for in
6 x 4 h pieces, from every provider, merged, one 3-C band kept (A.pick_band), and written to the
same raw path.  Station-days that stay empty get a manifest row (no file).

  python scripts/ab_refetch_nvi.py [--partial]   (--partial: 1 h pieces, 6 tries, for days < 99 % covered)
-> traces/<NET.STA>/<YYYY>/<DOY>.mseed ; results/ab_royer/refetch_NVI.csv
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

PROVIDERS = ["EARTHSCOPE", "https://www.earthquakescanada.nrcan.gc.ca",
             "https://earthquakescanada.nrcan.gc.ca"]
CHAN = "BH?,HH?,EH?,SH?"


PIECE_S = 14400
TRIES = 3


def one(job):
    sta, date = job
    from obspy import Stream, UTCDateTime
    from obspy.clients.fdsn import Client
    out = A.raw_path(sta, date)
    t = UTCDateTime(date)
    net = A.net_for(sta, date)
    t0 = time.time()
    best, prov_used, errs = Stream(), "", []
    for prov in PROVIDERS:
        st = Stream()
        try:
            c = Client(prov, timeout=180)
        except Exception as e:
            errs.append(f"{prov}: {e}"[:80]); continue
        for k in range(86400 // PIECE_S):
            for attempt in range(TRIES):
                try:
                    st += c.get_waveforms(net, sta, "*", CHAN, t + k * PIECE_S, t + (k + 1) * PIECE_S,
                                          attach_response=False)
                    break
                except Exception as e:
                    msg = str(e)
                    if "204" in msg or "No data" in msg:
                        break
                    time.sleep(2)
                    if attempt == TRIES - 1:
                        errs.append(f"{prov} piece{k}: {msg[:60]}")
        if len(st) and sum(tr.stats.npts for tr in st) > sum(tr.stats.npts for tr in best):
            best, prov_used = st, prov
    if os.path.exists(out):                    # --partial pass: union with what is on disk
        from obspy import read
        best += read(out)
        prov_used = prov_used + "+existing"
    if len(best) == 0:
        return dict(sta=sta, date=date, status="nodata", provider="", chans="", coverage_frac=0.0,
                    secs=round(time.time() - t0, 1), err="; ".join(errs)[:300])
    best.merge(method=1)
    best = best.split()
    best, _ = A.pick_band(best)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + f".tmp{os.getpid()}"
    best.write(tmp, format="MSEED")
    os.replace(tmp, out)
    chans = sorted({tr.stats.channel for tr in best})
    cov = sum(tr.stats.npts / tr.stats.sampling_rate for tr in best) / (86400.0 * max(len(chans), 1))
    return dict(sta=sta, date=date, status="data", provider=prov_used, chans=",".join(chans),
                coverage_frac=round(cov, 4), secs=round(time.time() - t0, 1), err="")


def main():
    global PIECE_S, TRIES
    M = pd.read_csv(f"{A.RES}/fetch_manifest_NVI.csv")
    miss = M[M.status != "data"]
    jobs = [(r.sta, r.date) for r in miss.itertuples() if not os.path.exists(A.raw_path(r.sta, r.date))]
    tag = ""
    if "--partial" in sys.argv:                # second pass: 1 h pieces, 6 tries, days < 99 % covered
        PIECE_S, TRIES, tag = 3600, 6, "_partial"
        R0 = pd.read_csv(f"{A.RES}/refetch_NVI.csv")
        jobs = [(r.sta, r.date) for r in R0.itertuples() if r.status == "data" and r.coverage_frac < 0.99]
    A.log(f"[refetch NVI] {len(jobs)} nodata station-days, 6 x 4 h pieces x {len(PROVIDERS)} providers")
    with ThreadPoolExecutor(max_workers=6) as ex:
        rows = list(ex.map(one, jobs))
    R = pd.DataFrame(rows)
    R.to_csv(f"{A.RES}/refetch_NVI{tag}.csv", index=False)
    A.log(R.groupby(["sta", "status"]).size().to_string())


if __name__ == "__main__":
    main()
