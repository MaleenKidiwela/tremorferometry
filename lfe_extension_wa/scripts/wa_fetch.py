#!/usr/bin/env python
"""Stage 1: download the permanent station-days that are not already on disk.

Reuses every raw day already under lfe_array_discovery/traces, lfe_extension/traces,
lfe_catalog_repair/traces (results/setup/station_days.csv `raw`).  Band = the band of the station's
metadata node(s) live that day (data/extension_stations/*_bridge_nodes_strict.csv).  Provider IRIS
(EarthScope); CN falls back to NRCan.  <= 8 concurrent connections.

  -> traces/<NET.STA>/<YYYY>/<DOY>.mseed  (.nd = verified empty)
  -> results/fetch_manifest.csv
"""
from __future__ import annotations
import os, sys, time, threading, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

NRCAN = "https://earthquakescanada.nrcan.gc.ca"
_CL = {}
LOG = f"{W.STUDY}/logs/fetch.log"


def client(prov):
    from obspy.clients.fdsn import Client
    k = (prov, threading.get_ident())
    if k not in _CL:
        _CL[k] = Client(prov, timeout=180)
    return _CL[k]


def fetch(ns, date, bands):
    from obspy import UTCDateTime
    net, sta = ns.split(".")
    ts = pd.Timestamp(date)
    p = f"{W.TRACES}/{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
    if os.path.exists(p):
        return dict(netsta=ns, date=date, status="data", bytes=os.path.getsize(p), note="resumed")
    if os.path.exists(p + ".nd"):
        return dict(netsta=ns, date=date, status="nodata", bytes=0, note="resumed")
    t0 = UTCDateTime(date)
    provs = ["EARTHSCOPE"] + ([NRCAN] if net in ("CN", "C8") else [])
    last = ""
    for prov in provs:
        for k in range(3):
            try:
                st = client(prov).get_waveforms(net, sta, "*", ",".join(b + "?" for b in bands), t0, t0 + 86400)
                if len(st) == 0:
                    break
                os.makedirs(os.path.dirname(p), exist_ok=True)
                st.write(p + ".tmp", format="MSEED")
                os.replace(p + ".tmp", p)
                return dict(netsta=ns, date=date, status="data", bytes=os.path.getsize(p),
                            note=f"{prov}:" + ",".join(sorted({tr.id for tr in st})))
            except Exception as e:
                last = f"{type(e).__name__}:{str(e)[:100]}"
                if "No data" in str(e) or "204" in str(e) or "NoData" in last or (prov != "EARTHSCOPE" and "404" in str(e)):
                    break
                time.sleep(5 * (k + 1))
    if last and not ("No data" in last or "204" in last or "NoData" in last or "404" in last):   # NRCan: 404 = no data
        return dict(netsta=ns, date=date, status="error", bytes=0, note=last)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p + ".nd", "w").close()
    return dict(netsta=ns, date=date, status="nodata", bytes=0, note=last)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--windows", default="")
    ap.add_argument("--days", default="", help="comma list of dates (smoke)")
    a = ap.parse_args()
    SD = pd.read_csv(f"{W.SETUP}/station_days.csv", keep_default_na=False)
    if a.windows:
        SD = SD[SD.window.isin(a.windows.split(","))]
    if a.days:
        SD = SD[SD.date.isin(a.days.split(","))]
    SD["raw"] = [W.raw_path(ns, d) or "" for ns, d in zip(SD.netsta, SD.date)]
    todo = SD[SD.meta_live.astype(str).isin(["True", "true", "1"]) & (SD.raw == "")]
    todo = todo.drop_duplicates(["netsta", "date"])
    if a.limit:
        todo = todo.head(a.limit)
    W.log(LOG, f"[fetch] {len(todo)} station-days to fetch, workers {a.workers}")
    rows = []
    t0 = time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        fut = {}
        for r in todo.itertuples():
            bands = sorted({n.split(".")[3] for n in r.meta_nodes.split(";") if n})
            fut[ex.submit(fetch, r.netsta, r.date, bands)] = r
        for i, f in enumerate(as_completed(fut)):
            rows.append(f.result())
            if (i + 1) % 100 == 0:
                gb = sum(x["bytes"] for x in rows) / 1e9
                W.log(LOG, f"[fetch] {i+1}/{len(todo)} {gb:.1f} GB {time.time()-t0:.0f}s")
    M = pd.DataFrame(rows)
    old = f"{W.RES}/fetch_manifest.csv"
    if os.path.exists(old) and len(M):
        M = pd.concat([pd.read_csv(old), M]).drop_duplicates(["netsta", "date"], keep="last")
    if len(M):
        M.to_csv(old, index=False)
        W.log(LOG, f"[fetch] done: {M.status.value_counts().to_dict()} {M.bytes.sum()/1e9:.1f} GB")
    else:
        W.log(LOG, "[fetch] nothing to do")


if __name__ == "__main__":
    main()
