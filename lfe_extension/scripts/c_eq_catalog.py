#!/usr/bin/env python
"""Track C — the regional earthquake catalog used for the EQ flag (fetched ONCE).

Reuses ../data/eq_catalog_pgc.csv (PNSN regional, 2009-2017, already on disk) and tops it up from
USGS FDSN for the whole 2009-2026 window over the station box, so the flag has the same meaning in
every epoch.  Writes results/trackC/eq_catalog.csv (ot, lat, lon, depth, mag, src).
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

OLD = "/home/jovyan/tremorferometry/data/eq_catalog_pgc.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t1", default="2009-01-01")
    ap.add_argument("--t2", default="2026-12-31")
    ap.add_argument("--minlat", type=float, default=46.5)
    ap.add_argument("--maxlat", type=float, default=50.5)
    ap.add_argument("--minlon", type=float, default=-126.5)
    ap.add_argument("--maxlon", type=float, default=-121.0)
    ap.add_argument("--minmag", type=float, default=-1.0)
    ap.add_argument("--out", default=f"{C.TRC}/eq_catalog.csv")
    a = ap.parse_args()
    os.makedirs(C.TRC, exist_ok=True)
    rows = []
    if os.path.exists(OLD):
        d = pd.read_csv(OLD)
        d["src"] = "eq_catalog_pgc"
        rows.append(d)
        print(f"[eq] {len(d)} events from {OLD}")
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    got = []
    for yr in range(int(a.t1[:4]), int(a.t2[:4]) + 1):
        for cl in ("USGS", "IRIS"):
            try:
                cat = Client(cl, timeout=180).get_events(
                    starttime=UTCDateTime(f"{yr}-01-01"), endtime=UTCDateTime(f"{yr+1}-01-01"),
                    minlatitude=a.minlat, maxlatitude=a.maxlat,
                    minlongitude=a.minlon, maxlongitude=a.maxlon, minmagnitude=a.minmag)
                for ev in cat:
                    o = ev.preferred_origin() or (ev.origins[0] if ev.origins else None)
                    if o is None:
                        continue
                    m = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
                    got.append(dict(ot=str(o.time), lat=float(o.latitude), lon=float(o.longitude),
                                    depth=float(o.depth or 0) / 1000.0,
                                    mag=float(m.mag) if m else np.nan, src=cl))
                print(f"[eq] {yr} {cl}: {len(cat)} events", flush=True)
                break
            except Exception as e:
                print(f"[eq] {yr} {cl} FAILED: {type(e).__name__}: {e}", flush=True)
    if got:
        rows.append(pd.DataFrame(got))
    D = pd.concat(rows, ignore_index=True)
    D["ot"] = pd.to_datetime(D.ot, errors="coerce", utc=True).dt.tz_localize(None)
    D = D.dropna(subset=["ot"]).sort_values("ot")
    # de-duplicate: same origin within 3 s and 20 km
    t = D.ot.astype("datetime64[ns]").astype("int64") / 1e9
    keep = np.ones(len(D), bool)
    lat = D.lat.to_numpy(); lon = D.lon.to_numpy(); tv = t.to_numpy()
    j = 0
    for i in range(len(D)):
        if not keep[i]:
            continue
        j = i + 1
        while j < len(D) and tv[j] - tv[i] <= 3.0:
            if keep[j] and C.epi_km(lat[j], lon[j], lat[i], lon[i]) < 20:
                keep[j] = False
            j += 1
    D = D[keep]
    D.to_csv(a.out, index=False)
    print(f"[eq] {len(D)} unique events {D.ot.min()} .. {D.ot.max()} -> {a.out}")


if __name__ == "__main__":
    main()
