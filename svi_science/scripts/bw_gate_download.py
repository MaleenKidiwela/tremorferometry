#!/usr/bin/env python
"""bw_gate_download -- G0 / G2 downloader for the SVI backward extension (2003-2013).

Layout = the forward convention (lfe_extension/scripts/download_day.py): <root>/<NET.STA>/<YYYY>/<DOY>.mseed, native
rate, no filtering / response removal / gap merging.  Root = lfe_extension/traces_backward (new folder for 2003-2013).

CN from NRCan/CNDC FDSN (https://www.earthquakescanada.nrcan.gc.ca/fdsnws).  NRCan answers the SAME request at random
with 200 / 204 / 404 (SCOPE R1: ~1/3 first-pass success), so a 204 or 404 is NEVER taken as "no data": every channel
is requested ON ITS OWN (explicit code, e.g. BHZ; the band wildcard failed more often) up to --tries times with
exponential backoff + jitter (1 s * 1.3^k, capped at 10 s; connect timeout 15 s, read 120 s).  No `.nd` marker is written; the ledger says what happened.
PB / UW from EarthScope (204 there is stable = no data; errors retried 5x).

Ledger (append): svi_science/results/backward_gates/tables/download_ledger.csv, one row per station-day:
  net, sta, band, date, centre, status (data / partial / nodata / error), bytes, n_channels, fs, coverage_frac
  (min over channels), n_frag, chans_ok, requests, secs, chan_log.

  python bw_gate_download.py --stations PGC,VGZ,NLLB,SNB --start 2008-05-01 --end 2008-06-04 [--lanes-nrcan 4]
"""
from __future__ import annotations
import os, sys, io, time, random, argparse, threading
from datetime import date, timedelta
import numpy as np
import pandas as pd
import requests
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
sys.path.insert(0, G.LFE_SCRIPTS)
from download_day import verify          # forward header check (n_channels, fs, coverage_frac, n_frag, note)

BASE = {"NRCAN": "https://www.earthquakescanada.nrcan.gc.ca/fdsnws", "EARTHSCOPE": "https://service.earthscope.org/fdsnws"}
COMP = {"CN": "ZNE", "UW": "ZNE", "PB": "Z12", "C8": "ZNE", "PO": "ZNE"}
LOCK = threading.Lock()
REPAIR = [False]


def band_on(sta, d):
    for net, s, band, fs, a, b in G.BW_EPOCHS:
        if s == sta and a <= d < b:
            return net, band
    F = pd.read_csv(f"{G.LFE}/results/station_epochs.csv")
    F = F[(F.sta == sta) & (F.start <= d) & (d < F.end)]
    return (F.net.iloc[0], F.band.iloc[0]) if len(F) else (None, None)


def es_first(net, d):
    return net in ("C8", "PO") or (net == "CN" and d >= "2010-01-01")


def cached_elsewhere(net, sta, d, dirs):
    x = pd.Timestamp(d)
    return any(os.path.exists(f"{c}/{net}.{sta}/{x.year}/{x.dayofyear:03d}.npz") for c in dirs)


def get(centre, net, sta, cha, d, timeout=(15, 120)):
    t0 = f"{d}T00:00:00"; t1 = f"{(date.fromisoformat(d) + timedelta(days=1)).isoformat()}T00:00:00"
    url = f"{BASE[centre]}/dataselect/1/query?net={net}&sta={sta}&loc=*&cha={cha}&starttime={t0}&endtime={t1}"
    try:
        r = requests.get(url, timeout=timeout)
    except Exception as e:
        return "exc", None
    if r.status_code != 200 or not r.content:
        return r.status_code, None
    try:
        from obspy import read
        return 200, read(io.BytesIO(r.content))
    except Exception:
        return "parse", None


def fetch(net, sta, band, d, root, tries):
    from obspy import Stream
    yr = d[:4]; doy = date.fromisoformat(d).timetuple().tm_yday
    p = f"{root}/{net}.{sta}/{yr}/{doy:03d}.mseed"
    rec = dict(net=net, sta=sta, band=band, date=d, centre="", status="", bytes=0, n_channels=0, fs="", coverage_frac="",
               n_frag="", chans_ok="", requests=0, secs=0.0, chan_log="")
    t_start = time.time()
    if os.path.exists(p) and REPAIR[0] and net in ("CN", "PO", "C8"):
        from obspy import read as _read
        have = _read(p, headonly=True)
        miss = [band + c for c in COMP[net] if not len(have.select(channel=band + c))]
        if miss:                                  # repair pass: re-request ONLY the missing channels, merge, rewrite
            st = _read(p); logs = []; nreq = 0
            for cha in miss:
                got = False
                for k in range(tries):
                    code, s = get("NRCAN", net, sta, cha, d); nreq += 1
                    if code == 200 and s is not None and len(s.select(channel=cha)):
                        st += s.select(channel=cha); got = True; logs.append(f"{cha}:ok@{k + 1}"); break
                    time.sleep(min(1.0 * 1.3 ** k, 10.0) * (0.75 + 0.5 * random.random()))
                if not got:
                    logs.append(f"{cha}:fail@{tries}({code})")
            st.write(p + ".tmp", format="MSEED"); os.replace(p + ".tmp", p)
            n, fs, cov, nf, note = verify(p, band)
            rec.update(centre="NRCAN", status="data" if n >= 3 else "partial", bytes=os.path.getsize(p), n_channels=n,
                       fs=f"{fs:g}", coverage_frac=f"{cov:.4f}", n_frag=nf, requests=nreq, secs=round(time.time() - t_start, 1),
                       chan_log="repair " + " ".join(logs))
            return rec
    if REPAIR[0] and not os.path.exists(p):         # repair passes only complete days on disk: a day with no file failed all
        rec.update(status="gap_skip", chan_log="repair: no file after first pass (3 x tries failed = a real gap)")   # channels
        return rec                                   # tries x 3 in the first pass (P(data exists) ~ 1e-5): do not re-ask
    if os.path.exists(p):
        n, fs, cov, nf, note = verify(p, band)
        rec.update(centre="disk", status="data" if n >= 3 else "partial", bytes=os.path.getsize(p), n_channels=n,
                   fs=f"{fs:g}", coverage_frac=f"{cov:.4f}", n_frag=nf, chan_log="resumed")
        return rec
    centre = "NRCAN" if net in ("CN", "C8", "PO") else "EARTHSCOPE"
    st = Stream(); logs = []; nreq = 0; ok = []
    if es_first(net, d):                          # POLARIS / CN >= 2010: EarthScope once (stable 204), else NRCan per channel
        code, s = get("EARTHSCOPE", net, sta, band + "?", d); nreq += 1
        if code == 200 and s is not None and len({tr.stats.channel for tr in s.select(channel=band + "?")}) >= 3:
            st = s.select(channel=band + "?"); centre = "EARTHSCOPE"; ok = sorted({tr.stats.channel for tr in st})
            locs = sorted({tr.stats.location for tr in st})
            if len(locs) > 1: st = st.select(location=locs[0])
        logs.append(f"ES:{code}")
    rec["centre"] = centre
    if centre == "NRCAN":
        for c in COMP[net]:
            if c and len(st.select(channel=band + c)):
                continue
            cha = band + c; got = False
            for k in range(tries):
                code, s = get(centre, net, sta, cha, d); nreq += 1
                if code == 200 and s is not None and len(s.select(channel=cha)):
                    st += s.select(channel=cha); got = True; logs.append(f"{cha}:ok@{k + 1}"); ok.append(cha); break
                time.sleep(min(1.0 * 1.3 ** k, 10.0) * (0.75 + 0.5 * random.random()))
            if not got:
                logs.append(f"{cha}:fail@{tries}({code})")
    else:
        code = None
        for k in range(5):
            code, s = get(centre, net, sta, band + "?", d); nreq += 1
            if code == 200 and s is not None:
                st = s.select(channel=band + "?"); break
            if code == 204:
                break
            time.sleep(3 * (k + 1))
        locs = sorted({tr.stats.location for tr in st})
        if len(locs) > 1:
            st = st.select(location=locs[0])
        ok = sorted({tr.stats.channel for tr in st})
        logs.append(f"{band}?:{code}:{','.join(ok)}")
    rec.update(requests=nreq, chans_ok=",".join(ok), chan_log=" ".join(logs))
    if len(st) == 0:
        rec.update(status="nodata" if centre == "EARTHSCOPE" else "error", secs=round(time.time() - t_start, 1))
        return rec
    os.makedirs(os.path.dirname(p), exist_ok=True)
    st.write(p + ".tmp", format="MSEED")
    n, fs, cov, nf, note = verify(p + ".tmp", band)
    os.replace(p + ".tmp", p)
    rec.update(status="data" if n >= 3 else "partial", bytes=os.path.getsize(p), n_channels=n, fs=f"{fs:g}",
               coverage_frac=f"{cov:.4f}", n_frag=nf, secs=round(time.time() - t_start, 1))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True); ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--root", default=G.TRACES_BW); ap.add_argument("--tries", type=int, default=20)
    ap.add_argument("--lanes-nrcan", type=int, default=8); ap.add_argument("--lanes-es", type=int, default=4)
    ap.add_argument("--skip-cached-in", default="", help="comma list of forward cache dirs: skip station-days cached there")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--tag", default=""); ap.add_argument("--repair", action="store_true",
                    help="second pass: for CN days on disk missing a channel, re-request only that channel")
    a = ap.parse_args()
    REPAIR[0] = a.repair
    days = [str(x.date()) for x in pd.date_range(a.start, a.end, freq="D")]
    jobs = []
    for sta in a.stations.split(","):
        for d in days:
            net, band = band_on(sta, d)
            if net is None:
                continue
            jobs.append((net, sta, band, d))
    skipd = [x for x in a.skip_cached_in.split(",") if x]
    if skipd:
        n0 = len(jobs); jobs = [j for j in jobs if not cached_elsewhere(j[0], j[1], j[3], skipd)]
        G.log(f"[dl {a.tag}] {n0 - len(jobs)} station-days already cached in forward caches: skipped")
    if a.dry_run:
        J = pd.DataFrame(jobs, columns=["net", "sta", "band", "date"])
        print(J.assign(nrcan=[not es_first(n_, d_) and n_ in ("CN", "PO", "C8") for n_, d_ in zip(J.net, J.date)]).groupby(["net", "nrcan"]).size().to_string())
        return
    os.makedirs(os.path.dirname(G.LEDGER), exist_ok=True)
    G.log(f"[dl {a.tag}] {len(jobs)} station-days {a.start}..{a.end} -> {a.root}")
    from concurrent.futures import ThreadPoolExecutor
    t0 = time.time(); n = [0]

    def run(j):
        rec = fetch(*j, a.root, a.tries); rec["tag"] = a.tag
        with LOCK:
            hdr = not os.path.exists(G.LEDGER)
            pd.DataFrame([rec]).to_csv(G.LEDGER, mode="a", header=hdr, index=False)
            n[0] += 1
            if n[0] % 25 == 0 or rec["status"] not in ("data",):
                G.log(f"  {n[0]}/{len(jobs)} {rec['net']}.{rec['sta']} {rec['date']} {rec['status']} req {rec['requests']} "
                      f"{rec['chan_log'][:80]} {time.time() - t0:.0f}s")
        return rec
    cn = [j for j in jobs if j[0] in ("CN", "PO") and not (j[0] == "CN" and es_first(j[0], j[3]))]
    es = [j for j in jobs if j not in cn]                      # EarthScope pool: PB, UW, C8 and CN >= 2010 (ES first)
    with ThreadPoolExecutor(a.lanes_nrcan) as e1, ThreadPoolExecutor(a.lanes_es) as e2:
        f1 = [e1.submit(run, j) for j in cn]; f2 = [e2.submit(run, j) for j in es]
        res = [f.result() for f in f1 + f2]
    R = pd.DataFrame(res)
    G.log(f"[dl {a.tag}] done {time.time() - t0:.0f}s: {R.groupby(['net', 'status']).size().to_dict()}; "
          f"{R.bytes.astype(float).sum() / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
