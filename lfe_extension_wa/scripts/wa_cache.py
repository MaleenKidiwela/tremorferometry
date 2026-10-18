#!/usr/bin/env python
"""Stage 2: 20 Hz cache on the FROZEN chain (lfe_extension e4_cache20.py / c_cache_month.py):
merge(method=1, fill_value=0) with the CN per-fragment fs snap -> band HH > BH > EH -> component
order Z / N|1 / E|2 resolved, NOT rotated (missing slots are zero and flagged not live) -> 1-8 Hz
Butterworth-4 zero-phase designed at the NATIVE rate -> resample_poly to 20 Hz (50 Hz: up 2 / down 5)
-> placed on a midnight-aligned 1,728,000-sample grid.  `m` = sample inside recorded data, dilated
5 s at gap edges (the band-pass transient) -- the Royer-chain validity mask.

  -> cache20hz/<NET.STA>/<YYYY>/<DOY>.npz : x (3,N) f32, m (N,) bool, live (3,) bool, order, t0,
     fs_native, band, node
  -> results/cache_manifest.csv
"""
from __future__ import annotations
import os, sys, argparse, time
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

LOG = f"{W.STUDY}/logs/cache.log"
BANDS = ["HH", "BH", "EH"]
_N = None


def nodes():
    global _N
    if _N is None:
        _N = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]]
                       ).drop_duplicates("node")
        extra = os.environ.get("WA_EXTRA_NODES", "")      # step 4b: epochs of the added stations
        if extra:
            X = pd.read_csv(extra)
            _N = pd.concat([_N[~(_N.net + "." + _N.sta).isin(set(X.net + "." + X.sta))], X]).drop_duplicates("node")
    return _N


def node_for(ns, band, fs, date, zlive=True):
    net, sta = ns.split(".")
    N = nodes()
    n = N[(N.net == net) & (N.sta == sta) & (N.band == band) & (np.isclose(N.fs, fs)) & (N.start <= date) & (N.end > date)]
    if zlive and n.comps.astype(str).str.contains("Z").any():
        n = n[n.comps.astype(str).str.contains("Z")]
    if len(n):
        return n.sort_values("start").node.iloc[-1]
    return f"{ns}..{band}.{fs:g}@data"


def slot(ch):
    c = ch[-1].upper()
    return 0 if c == "Z" else 1 if c in "N1" else 2 if c in "E2" else -1


def cache_one(job):
    ns, date, raw = job
    out = W.cache_path(ns, date)
    if os.path.exists(out):
        z = np.load(out)
        return dict(netsta=ns, date=date, status="cached", node=str(z["node"]), band=str(z["band"]),
                    fs_native=float(z["fs_native"]), live="".join("1" if v else "0" for v in z["live"]),
                    frac_valid=float(z["m"].mean()))
    from obspy import read, UTCDateTime
    try:
        st = read(raw)
    except Exception as e:
        return dict(netsta=ns, date=date, status=f"unreadable:{type(e).__name__}")
    D0 = UTCDateTime(date)
    for tr in st:
        fs = tr.stats.sampling_rate
        if abs(fs - round(fs)) / fs < 1e-4:
            tr.stats.sampling_rate = float(round(fs))
    chans = {tr.stats.channel for tr in st}
    band = next((b for b in BANDS if any(c.startswith(b) and slot(c) >= 0 for c in chans)), None)
    if band is None:
        return dict(netsta=ns, date=date, status="no_band:" + ",".join(sorted(chans)))
    st = st.select(channel=band + "?")
    locs = sorted({tr.stats.location for tr in st})
    loc = "" if "" in locs else locs[0]
    st = st.select(location=loc)
    fss = sorted({tr.stats.sampling_rate for tr in st})
    fs0 = max(fss, key=lambda f: sum(tr.stats.npts for tr in st if tr.stats.sampling_rate == f))
    st = st.select(sampling_rate=fs0)
    st.trim(D0, D0 + 86400 - 1.0 / fs0)
    X = np.zeros((3, W.NDAY), np.float32)
    live = np.zeros(3, bool)
    order = ["", "", ""]
    mk = np.zeros((3, W.NDAY), bool)
    fr = Fraction(W.FS / fs0).limit_denominator(100)
    for ch in sorted({tr.stats.channel for tr in st}):
        k = slot(ch)
        if k < 0 or live[k]:
            continue
        s1 = st.select(channel=ch).copy()
        s1.merge(method=1, fill_value=None)
        tr = s1[0]
        data = tr.data
        if np.ma.isMaskedArray(data):
            gap = np.ma.getmaskarray(data).copy(); data = data.filled(0)
        else:
            gap = np.zeros(len(data), bool)
        data = data.astype(np.float64)
        data = data - (np.median(data[~gap]) if (~gap).any() else 0.0)
        data[gap] = 0.0
        if len(data) < fs0 * 600:
            continue
        y = W.bp(data, fs0, 1.0, 8.0)
        if fr != 1:
            from scipy.signal import resample_poly
            y = resample_poly(y, fr.numerator, fr.denominator)
            g = gap[::int(round(1 / float(fr)))] if float(1 / fr).is_integer() else \
                np.interp(np.arange(len(y)) / W.FS, np.arange(len(gap)) / fs0, gap.astype(float)) > 0
        else:
            g = gap
        i0 = int(round((tr.stats.starttime - D0) * W.FS))
        a, b = max(i0, 0), min(i0 + len(y), W.NDAY)
        if b <= a:
            continue
        X[k, a:b] = y[a - i0:b - i0]
        mk[k, a:b] = ~g[a - i0:b - i0][: b - a]
        live[k] = True
        order[k] = ch
    if not live.any():
        return dict(netsta=ns, date=date, status="empty")
    m = mk[live].all(0)
    # dilate the invalid region by 5 s (band-pass transient at gap edges)
    bad = ~m
    if bad.any():
        n5 = int(5 * W.FS)
        c = np.concatenate([[0], np.cumsum(bad)])
        idx = np.arange(W.NDAY)
        lo = np.clip(idx - n5, 0, W.NDAY); hi = np.clip(idx + n5 + 1, 0, W.NDAY)
        m = (c[hi] - c[lo]) == 0
    node = node_for(ns, band, fs0, date, bool(live[0]))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out + ".tmp.npz", x=X, m=m, live=live, order=np.array(order), t0=np.float64(D0.timestamp),
             fs_native=np.float64(fs0), band=band, node=node, fs=np.float64(W.FS))
    os.replace(out + ".tmp.npz", out)
    return dict(netsta=ns, date=date, status="ok", node=node, band=band, fs_native=fs0,
                live="".join("1" if v else "0" for v in live), frac_valid=float(m.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--windows", default="")
    ap.add_argument("--days", default="", help="comma list of dates (smoke)")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    SD = pd.read_csv(f"{W.SETUP}/station_days.csv", keep_default_na=False)
    if a.windows:
        SD = SD[SD.window.isin(a.windows.split(","))]
    if a.days:
        SD = SD[SD.date.isin(a.days.split(","))]
    SD = SD.drop_duplicates(["netsta", "date"])
    jobs = []
    for r in SD.itertuples():
        p = W.raw_path(r.netsta, r.date)
        if p:
            jobs.append((r.netsta, r.date, p))
    if a.limit:
        jobs = jobs[: a.limit]
    W.log(LOG, f"[cache] {len(jobs)} station-days with raw data, workers {a.workers}")
    rows, t0 = [], time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(cache_one, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            try:
                rows.append(f.result())
            except Exception as e:
                rows.append(dict(status=f"crash:{type(e).__name__}:{str(e)[:80]}"))
            if (i + 1) % 200 == 0:
                W.log(LOG, f"[cache] {i+1}/{len(jobs)} {time.time()-t0:.0f}s")
    M = pd.DataFrame(rows)
    p = f"{W.RES}/cache_manifest.csv"
    if os.path.exists(p) and len(M):
        M = pd.concat([pd.read_csv(p), M]).drop_duplicates(["netsta", "date"], keep="last")
    M.to_csv(p, index=False)
    W.log(LOG, f"[cache] done {M.status.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
