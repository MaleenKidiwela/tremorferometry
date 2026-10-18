#!/usr/bin/env python
"""rf01_download: P-window waveforms (P-150 s .. P+150 s) for every arrival in data/rf_teleseis/arrivals_*.csv.

CN from NRCAN/CNDC FDSN, others from EarthScope FDSN; dataselect bulk requests, 40 events per request, <= 3 stations in parallel (one request
in flight per worker). Cached + resumable: one file per event (raw/<NET>.<STA>/<evid>.mseed), a DONE marker per
chunk (raw/<NET>.<STA>/_chunks/cNNNN.done, holds n traces returned), a station DONE marker
(raw/<NET>.<STA>.DONE) and raw/DOWNLOAD_DONE when all stations are finished.
Usage: python -P rf01_download.py [--workers 3] [--stations PGC,SNB]
"""
import argparse, os, sys, time, traceback
from multiprocessing import Pool
import pandas as pd
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client

ROOT = "/home/jovyan/tremorferometry"
D = f"{ROOT}/data/rf_teleseis"
RAW = f"{D}/raw"
LOGD = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin/logs"
PRE, POST, NCH = 150.0, 150.0, 40


def log(fh, m):
    fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n"); fh.flush()


def do_station(key):
    net, sta = key
    sd = f"{RAW}/{net}.{sta}"; cd = f"{sd}/_chunks"
    os.makedirs(cd, exist_ok=True)
    fh = open(f"{LOGD}/download_{net}.{sta}.log", "a")
    if os.path.exists(f"{RAW}/{net}.{sta}.DONE"):
        return key, "already"
    # CN data: Canadian National Data Centre FDSN (IRIS copy of CN lacks most pre-2010 data); others: EarthScope
    cl = Client(base_url="https://earthquakescanada.nrcan.gc.ca", timeout=300) if net == "CN" else \
        Client("EARTHSCOPE", timeout=300)
    xml = f"{RAW}/{net}.{sta}.xml"
    if not os.path.exists(xml):
        for k in range(4):
            try:
                inv = cl.get_stations(network=net, station=sta, channel="BH?,HH?", level="response")
                inv.write(xml + ".tmp", format="STATIONXML"); os.replace(xml + ".tmp", xml); break
            except Exception as e:
                log(fh, f"stationxml retry {k}: {e}"); time.sleep(20 * (k + 1))
    A = pd.read_csv(f"{D}/arrivals_{net}.{sta}.csv", keep_default_na=False)
    A = A.sort_values("tP").reset_index(drop=True)
    nch = (len(A) + NCH - 1) // NCH
    log(fh, f"start {net}.{sta}: {len(A)} events, {nch} chunks")
    ntot = 0
    for c in range(nch):
        mk = f"{cd}/c{c:04d}.done"
        if os.path.exists(mk):
            continue
        S = A.iloc[c * NCH:(c + 1) * NCH]
        bulk = []
        for _, r in S.iterrows():
            loc, cha = r["chan"].split("|")
            if loc.replace(".", "").isdigit():   # 2026-10-14: CSV round-trip turned loc "00" into "0.0" (US.NLWA) -> server 422
                loc = f"{int(float(loc)):02d}"
            tp = UTCDateTime(str(r["tP"]))
            bulk.append((net, sta, loc if loc else "--", cha, tp - PRE, tp + POST))
        st = None
        for k in range(5):
            try:
                st = cl.get_waveforms_bulk(bulk)
                break
            except Exception as e:
                msg = str(e)
                if "No data" in msg or "204" in msg or "HTTP code: 404" in msg:  # NRCAN answers 404 for no data
                    st = Stream(); break
                log(fh, f"chunk {c} retry {k}: {msg[:200]}"); time.sleep(30 * (k + 1))
        if st is None:
            log(fh, f"chunk {c} FAILED after retries; left without marker (rerun resumes)")
            continue
        n = 0
        for (_, r), b in zip(S.iterrows(), bulk):
            sub = st.slice(b[4], b[5]).copy()
            if len(sub) == 0:
                continue
            fn = f"{sd}/{r['evid']}.mseed"
            try:
                sub.write(fn + ".tmp", format="MSEED"); os.replace(fn + ".tmp", fn); n += 1
            except Exception as e:
                log(fh, f"write fail {r['evid']}: {e}")
        open(mk, "w").write(f"{len(st)} traces, {n} events\n")
        ntot += n
        if c % 10 == 0:
            log(fh, f"chunk {c}/{nch}: {n} events written (session total {ntot})")
    left = [c for c in range(nch) if not os.path.exists(f"{cd}/c{c:04d}.done")]
    if not left:
        open(f"{RAW}/{net}.{sta}.DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
        log(fh, "DONE")
        return key, "done"
    log(fh, f"incomplete: {len(left)} chunks left")
    return key, f"incomplete {len(left)}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--stations", default="")
    a = ap.parse_args()
    S = pd.read_csv(f"{D}/stations.csv")
    keys = [(r.net, r.sta) for r in S.itertuples()]
    if a.stations:
        want = a.stations.split(","); keys = [k for k in keys if k[1] in want]
    os.makedirs(LOGD, exist_ok=True)
    with Pool(min(a.workers, 3)) as P:
        res = list(P.imap_unordered(do_station, keys))
    print(res, flush=True)
    if all(r[1] in ("done", "already") for r in res) and not a.stations:
        open(f"{RAW}/DOWNLOAD_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
