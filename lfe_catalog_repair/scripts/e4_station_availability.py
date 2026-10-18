#!/usr/bin/env python
"""E4 step 6 — METADATA ONLY: which stations were recording 3-C at >= 40 Hz, and when?

No waveforms are requested anywhere in this script (level="channel" station metadata only).

Three question sets:
  A. the 21 SeisBench template stations (the stack metadata labels them all "PO"; POLARIS stations
     can appear under PO or CN, so every network is queried and the answer recorded) for
     2010-08-01..2010-09-30 and 2012-08-15..2012-10-15;
  B. the current outer ring + PB boreholes for 2019-01-01..2026-09-01;
  C. data volume of the two candidate downloads at the NATIVE sample rate.

Writes results/e4_station_availability.csv and results/e4_station_availability_summary.md
"""
from __future__ import annotations
import sys, time
import numpy as np, pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"

TEMPLATE21 = ["GLBC", "GOWB", "JRBC", "KELB", "KHVB", "LCBC", "LZB", "MGCB", "PFB", "PGC",
              "SHVB", "SILB", "SNB", "SOKB", "SSIB", "TSJB", "TWBB", "TWGB", "TWKB", "VGZ", "YOUB"]
CURRENT = [("CN", "PGC"), ("CN", "LZB"), ("CN", "VGZ"), ("CN", "SNB"), ("CN", "YOUB"),
           ("CN", "PFB"), ("C8", "SOKB"), ("CN", "NLLB"), ("PB", "B011"), ("PB", "B926"),
           ("UW", "LRIV"), ("UW", "OCP")]
WINDOWS = {"2010ETS": ("2010-08-01", "2010-09-30"),
           "2012ETS": ("2012-08-15", "2012-10-15"),
           "post2019": ("2019-01-01", "2026-09-01")}
NETS = "CN,PO,C8,PB,UW,UO,NV,BK,TA,US"


def query(client, name, net, sta, t0, t1):
    try:
        inv = client.get_stations(network=net, station=sta, channel="?H?,?N?",
                                  starttime=UTCDateTime(t0), endtime=UTCDateTime(t1),
                                  level="channel")
    except Exception as e:
        return [], f"{type(e).__name__}: {str(e)[:90]}"
    rows = []
    for nw in inv:
        for st in nw:
            for ch in st:
                rows.append(dict(source=name, net=nw.code, sta=st.code, loc=ch.location_code,
                                 chan=ch.code, fs=float(ch.sample_rate or 0),
                                 start=str(ch.start_date)[:10],
                                 end=str(ch.end_date)[:10] if ch.end_date else "open",
                                 lat=float(st.latitude), lon=float(st.longitude)))
    return rows, ""


def main():
    clients = {}
    for name in ("EARTHSCOPE", "IRIS", "NRCAN"):
        try:
            clients[name] = Client(name, timeout=60)
            print(f"[fdsn] {name} ok", flush=True)
        except Exception as e:
            print(f"[fdsn] {name} FAILED: {e}", flush=True)
    rows, errs = [], []
    tasks = []
    for w, (t0, t1) in WINDOWS.items():
        if w == "post2019":
            for net, sta in CURRENT:
                tasks.append((w, net, sta, t0, t1))
        else:
            for sta in TEMPLATE21:
                tasks.append((w, NETS, sta, t0, t1))
            for net, sta in CURRENT:
                tasks.append((w, net, sta, t0, t1))
    for w, net, sta, t0, t1 in tasks:
        got = False
        for name, cl in clients.items():
            r, e = query(cl, name, net, sta, t0, t1)
            if r:
                for d in r:
                    d["window"] = w
                rows.extend(r)
                got = True
            elif e:
                errs.append((name, w, net, sta, e))
        if not got:
            rows.append(dict(source="none", net=net, sta=sta, loc="", chan="", fs=np.nan,
                             start="", end="", lat=np.nan, lon=np.nan, window=w))
    D = pd.DataFrame(rows).drop_duplicates(subset=["window", "net", "sta", "loc", "chan",
                                                   "start", "end"])
    D.to_csv(f"{RES}/e4_station_availability.csv", index=False)
    print(f"[write] {RES}/e4_station_availability.csv  ({len(D)} channel epochs)")
    for e in errs[:15]:
        print("  [err]", e)
    return D


if __name__ == "__main__":
    main()
