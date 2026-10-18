#!/usr/bin/env python
"""Step 1b bed -- download + frozen-chain cache of the CLOSE 2012 stations.

Adds to the Step-1a bed the stations that actually exist in the 2012-08-15..10-05 window within
~60 km epicentral of the Bostock family field:

    C8.SOKB   epi_med 24.8 km, 67/130 families <= 25 km   (station #1 of Bostock's own 21)
    CN.LZB    epi_med 28.0 km, 54 families                (#4 of Bostock's 21)
    C8.PA01..PA05, PA12   co-located with PGC (< 0.6 km)  -> 6 more members for the PGCq SITE
    C8.BPCB   epi_med 53.0 km, 19 families
    CN.PFB    epi_med 55.8 km,  9 families
    CN.YOUB   epi_med 59.8 km,  6 families   (WITHHELD -- the dense Gate-1 confirmer)

Download: `r_download_day.fetch_day` verbatim (retry-wrapped FDSN, `.nd` no-data markers,
resumable, NATIVE rate), <= 4 concurrent so agent D's 8-worker queue is not starved.
Cache: `o_cache_month.cache_one` verbatim -- merge(method=1) -> ZNE order -> Butterworth-4
zero-phase 1-8 Hz at the native rate -> integer decimation to 20 Hz (HH 100 Hz / 5, BH 40 Hz / 2:
both integer, so NO rational twin is needed; any station that is not an integer multiple raises
and is reported).

  python scripts/w_fetch_bed2012b.py --workers 4
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import r_download_day as RD                                                 # noqa: E402
import o_cache_month as OC                                                  # noqa: E402

TRACES = f"{W.STUDY}/traces/bed2012_dense"
STEP1B_CSV = W.STEP1B_STATIONS_CSV
CACHE = W.STUDY_CACHE
PROVIDERS = {"IRIS": "IRIS", "NRCAN": "https://earthquakescanada.nrcan.gc.ca"}


def roster():
    """(net, sta, band, [providers in order]) -- the provider column comes from the delivery probe
    (results/step1b_stations.csv); both providers are still TRIED, because metadata is not
    delivery and delivery is day-dependent (C8 delivers partially from NRCan only).

    Reads its OWN csv, not `W.extra_stations()`: that helper is the UNION of the Step-1b and
    Step-1c rosters (it has to be, so the station table knows every bed's coordinates), and using
    it here made the Step-1b fetch try all 27 stations x 52 days of 2012 -- including stations
    that only exist in 2005."""
    E = pd.read_csv(STEP1B_CSV)
    out = []
    for r in E.itertuples():
        first = getattr(r, "prov", "IRIS")
        provs = [first] + [p for p in PROVIDERS if p != first]
        out.append((r.net, r.sta, getattr(r, "band", "HH"), provs))
    return out


def fetch_any(net, sta, band, d, provs):
    """fetch_day against each provider in turn; the first that returns data wins."""
    rec = None
    for pn in provs:
        rec = RD.fetch_day(net, sta, band, d, "", TRACES, PROVIDERS[pn])
        if rec.get("status") == "data":
            rec["note"] = (str(rec.get("note", "")) + f";prov={pn}").strip(";")
            return rec
        # a `.nd` marker from a provider that simply does not serve this network would block the
        # other provider on the next pass -- remove it unless every provider has been tried
        nd = RD.day_path(TRACES, net, sta, d) + ".nd"
        if pn != provs[-1] and os.path.exists(nd):
            os.remove(nd)
    if rec is not None:
        rec["note"] = (str(rec.get("note", "")) + ";all_providers_nodata").strip(";")
    return rec


def days():
    d0 = date.fromisoformat(W.STEP1_DAY0)
    d1 = date.fromisoformat(W.STEP1_DAY1)
    return [d0 + timedelta(days=k) for k in range((d1 - d0).days + 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cache-workers", type=int, default=12)
    ap.add_argument("--prov", default="IRIS")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(TRACES, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)
    R = roster()
    D = days()
    jobs = [(n, s, b, d, pv) for (n, s, b, pv) in R for d in D]
    if a.limit:
        jobs = jobs[: a.limit]
    W.log(f"[fetch] {len(R)} stations x {len(D)} days = {len(jobs)} station-days -> {TRACES}")
    recs = []
    with ThreadPoolExecutor(max_workers=min(a.workers, 4)) as ex:
        futs = [ex.submit(fetch_any, n, s, b, d, pv) for (n, s, b, d, pv) in jobs]
        for k, f in enumerate(futs):
            recs.append(f.result())
            if (k + 1) % 50 == 0:
                ok = sum(1 for r in recs if r["status"] == "data")
                W.log(f"[fetch] {k+1}/{len(jobs)} ({ok} data)")
    M = pd.DataFrame(recs)
    M.to_csv(f"{W.RES}/step1b_download_manifest.csv", index=False)
    W.log("[fetch] " + M.status.value_counts().to_dict().__str__())
    per = M[M.status == "data"].groupby("sta").size()
    W.log(f"[fetch] station-days with data: {per.to_dict()}")

    cjobs = []
    for r in M[M.status == "data"].itertuples():
        netsta = f"{r.net}.{r.sta}"
        d = date.fromisoformat(r.date)
        src = RD.day_path(TRACES, r.net, r.sta, d)
        cjobs.append((netsta, d.year, int(d.timetuple().tm_yday), src, CACHE, False,
                      float(r.coverage_frac) if r.coverage_frac != "" else 0.0))
    W.log(f"[cache] {len(cjobs)} station-days -> {CACHE}")
    out = []
    with ThreadPoolExecutor(max_workers=a.cache_workers) as ex:
        for k, rr in enumerate(ex.map(OC.cache_one, cjobs)):
            out.append(rr)
            if (k + 1) % 100 == 0:
                W.log(f"[cache] {k+1}/{len(cjobs)}")
    C = pd.DataFrame(out) if out else pd.DataFrame(columns=["netsta", "doy", "err"])
    if "err" not in C.columns:
        C["err"] = ""
    C.to_csv(f"{W.RES}/step1b_cache_manifest.csv", index=False)
    bad = C[C.err.astype(str) != ""]
    W.log(f"[cache] ok {len(C) - len(bad)} / {len(C)}; errors {len(bad)}")
    if len(bad):
        W.log("[cache] error sample: " + str(bad.head(5)[["netsta", "doy", "err"]].to_dict("records")))
    per = C[C.err.astype(str) == ""].groupby("netsta").size()
    W.log(f"[cache] cached station-days: {per.to_dict()}")


if __name__ == "__main__":
    main()
