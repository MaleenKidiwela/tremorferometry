#!/usr/bin/env python
"""Stage 5: HF-ratio column for every data detection >= T_ext (the SVI HF veto, WA equivalent).

SVI: E(10-19 Hz)/E(2-8 Hz) at CN.PGC from the NATIVE-rate trace over the scored 12-s window,
Butterworth-4 zero-phase, summed over 3 components; cut = 99th percentile of build-half detections
matched to an own-family catalogue event, per instrument epoch (a4_hf_veto_cuts.csv).
WA equivalent (no single station sees all 329 families): the HF station of a family is its
highest-weight clear template station live that day (all our stations are >= 40 Hz native, so the
10-19 Hz band is below Nyquist); the cut is calibrated per HF station NODE by the same 99th-percentile
rule (wa_score.py).  Filters are designed with fs = the native rate.

  python scripts/wa_hf.py --mf val --tpl build
  -> results/mf_<mf>/hf.parquet  (fid, t, hf_sta, hf_node, hf_fs, hfr)
"""
from __future__ import annotations
import os, sys, argparse, glob, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

LOG = f"{W.STUDY}/logs/hf.log"


def ratios(job):
    ns, date, node, band, rows = job
    from obspy import read, UTCDateTime
    p = W.raw_path(ns, date)
    out = []
    if p is None:
        return out
    try:
        st = read(p).select(channel=band + "?")
    except Exception:
        return out
    locs = sorted({tr.stats.location for tr in st})
    if not locs:
        return out
    st = st.select(location="" if "" in locs else locs[0])
    for tr in st:
        fs = tr.stats.sampling_rate
        if abs(fs - round(fs)) / fs < 1e-4:
            tr.stats.sampling_rate = float(round(fs))
    D0 = UTCDateTime(date)
    E28, E1019, fs0 = None, None, None
    for ch in sorted({tr.stats.channel for tr in st}):
        s1 = st.select(channel=ch).copy()
        s1.merge(method=1, fill_value=0)
        tr = s1[0]
        fs0 = tr.stats.sampling_rate
        if fs0 < 40:
            return out
        x = tr.data.astype(np.float64); x -= np.median(x)
        i0 = int(round((tr.stats.starttime - D0) * fs0))
        N = int(86400 * fs0)
        y = np.zeros(N)
        a, b = max(i0, 0), min(i0 + len(x), N)
        if b <= a:
            continue
        y[a:b] = x[a - i0:b - i0]
        e1 = W.bp(y, fs0, 2.0, 8.0) ** 2
        e2 = W.bp(y, fs0, W.HF_BAND[0], W.HF_BAND[1]) ** 2
        E28 = e1 if E28 is None else E28 + e1
        E1019 = e2 if E1019 is None else E1019 + e2
    if E28 is None:
        return out
    c1 = np.concatenate([[0], np.cumsum(E28)]); c2 = np.concatenate([[0], np.cumsum(E1019)])
    nw = int(W.TPL_S * fs0)
    D0t = pd.Timestamp(date).timestamp()
    for fid, t, off in rows:
        j = int(round((t + off - D0t) * fs0))
        if j < 0 or j + nw > len(E28):
            continue
        s28 = c1[j + nw] - c1[j]; s1019 = c2[j + nw] - c2[j]
        out.append((fid, t, ns, node, fs0, float(s1019 / s28) if s28 > 0 else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mf", required=True)
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--cuts", default="", help="hf_cuts.csv: prefer HF stations whose node has a calibrated cut")
    ap.add_argument("--all-templates", action="store_true", help="Merlin T1: HF station may be any >= 30-event "
                    "template station (clear first, then w)")
    a = ap.parse_args()
    od = f"{W.RES}/mf_{a.mf}"
    DET = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/det_*.parquet"))], ignore_index=True)
    DET = DET[(DET.series == "data") & (DET.net_mad >= W.T_EXT)]
    FD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/famday_*.parquet"))], ignore_index=True)
    Q = pd.concat([pd.read_csv(f"{W.RES}/templates/tpl_qc_{t}.csv") for t in a.tpl.split(",")])
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    if a.all_templates:        # Merlin T1: HF station from every >= 30-event template, clear ones preferred
        Q = Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna() & Q.w.notna()].copy()
    else:
        Q = Q[Q.clear].copy()
    Q["has_cut"] = False
    if a.cuts and os.path.exists(a.cuts):
        C = pd.read_csv(a.cuts)
        Q["has_cut"] = Q.node.isin(set(C[C.cut.notna()].hf_node))
    Q = Q.sort_values(["has_cut", "clear", "w"], ascending=[False, False, False])
    import json
    NODES = {}
    for p in glob.glob(f"{od}/run_*.json"):
        r = json.load(open(p)); NODES[r["date"]] = r["nodes"]
    jobs = {}
    FST = {(int(f), d): v for f, d, v in zip(FD.fid, FD.date, FD.stations)}
    QG = dict(tuple(Q.groupby("fid")))
    for (f, ds), g in DET.groupby(["fid", "date"]):
        v = FST.get((int(f), ds), "")
        stas = set(v.split(",")) if v else set()
        qf = QG.get(f, Q.iloc[:0])
        cand = qf[qf.netsta.isin(stas)]
        cand = cand[[NODES.get(ds, {}).get(ns) == nd for ns, nd in zip(cand.netsta, cand.node)]]
        if not len(cand):
            continue
        r = cand.iloc[0]
        band = r.node.split(".")[3]
        jobs.setdefault((r.netsta, ds, r.node, band), []).extend((int(f), float(t), float(r.offset_s)) for t in g.t)
    W.log(LOG, f"[hf:{a.mf}] {len(DET)} detections >= T_ext, {len(jobs)} station-day jobs")
    rows, t0 = [], time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(ratios, (k[0], k[1], k[2], k[3], v)) for k, v in jobs.items()]
        for i, f in enumerate(as_completed(fut)):
            rows += f.result()
            if (i + 1) % 200 == 0:
                W.log(LOG, f"[hf:{a.mf}] {i+1}/{len(fut)} {time.time()-t0:.0f}s")
    H = pd.DataFrame(rows, columns=["fid", "t", "hf_sta", "hf_node", "hf_fs", "hfr"])
    H.to_parquet(f"{od}/hf.parquet", index=False)
    W.log(LOG, f"[hf:{a.mf}] wrote {len(H)} ratios; median {H.hfr.median():.3f}")


if __name__ == "__main__":
    main()
