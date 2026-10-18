#!/usr/bin/env python
"""s11_wa_windows -- S11 step 7 (PLAN §5): windows for the WA inventory with the UNCHANGED S10 reader/processing
(s10_windows.task: 2-8 Hz Butterworth-4 zero-phase, 40 Hz, despike 8 x MAD, RMS-normalised, 3C, -14..+10 s about S,
read span S-15..S+30 s). Only the raw-file location changes: s10_windows.Raw is subclassed to read the WA raw roots
(wa_prod_common.prod_roots(): wa_common.TRACE_ROOTS + traces_dl2 + the production extra roots; first readable file with a
Z channel wins). Data rows and scrambled-control rows (c1 rows read at frame time t + S offset = "c2") are separate
s10_windows.task calls (same despike reference as SVI) on a per-worker shared day cache, so each day file is read once.
Writes C.MWIN/<NET.STA>.npz and <NET.STA>__c2.npz, C.MTB/window_sizes.csv.
    S11_REGION=WA_NW python s11_wa_windows.py --workers 14
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import s10_common as S
import s10_windows as SW
import s11_common as C
import wa_prod_common as WP

ROOTS = WP.prod_roots()
_CACHE = {}


class WARaw(SW.Raw):
    """s10_windows.Raw with the WA raw roots; day cache shared by all calls in one worker process."""
    def __init__(self, netsta):
        self.ns = netsta
        self.cache = _CACHE.setdefault(netsta, {})

    def day(self, d):
        if d in self.cache:
            return self.cache[d]
        import obspy
        y, doy = d.year, d.timetuple().tm_yday
        out = None
        for base in ROOTS:
            p = f"{base}/{self.ns}/{y}/{doy:03d}.mseed"
            if not (os.path.exists(p) and os.path.getsize(p) > 0):
                continue
            try:
                st = obspy.read(p)
            except Exception:
                continue
            bands = {}
            for tr in st:
                c = tr.stats.channel
                if c[-1] in SW.CORD:
                    bands.setdefault(c[:2], set()).add(SW.CORD[c[-1]])
            full = [b for b, cs in bands.items() if 0 in cs]
            if not full:
                continue
            band = sorted(full, key=lambda b: (-len(bands[b]), {"HH": 0, "EH": 1, "BH": 2}.get(b, 3)))[0]
            st = st.select(channel=f"{band}?")
            comps = {}
            for c in set(tr.stats.channel for tr in st):
                s2 = st.select(channel=c).copy()
                try:
                    s2.merge(method=0, fill_value=0)
                except Exception:
                    continue
                tr = s2[0]
                comps[SW.CORD[c[-1]]] = (tr.data.astype(np.float64), float(tr.stats.sampling_rate), float(tr.stats.starttime.timestamp))
            out = comps
            break
        self.cache[d] = out
        while len(self.cache) > 4:
            self.cache.pop(next(iter(self.cache)))
        return out


SW.Raw = WARaw          # s10_windows.task resolves Raw at call time


def task2(arg):
    ns, d, Rd, Rc = arg
    out = []
    for R in (Rd, Rc):
        out.append(SW.task((ns, d, R)) if len(R) else (ns, d, np.zeros(0, np.int64), None, None, None, "empty"))
    _CACHE.pop(ns, None)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=14); ap.add_argument("--stations", default="")
    a = ap.parse_args()
    os.makedirs(C.MWIN, exist_ok=True)
    log = C.Log(f"{C.LOGD}/wa_windows.log")
    RS = pd.read_parquet(f"{C.MTB}/row_station.parquet")
    X = pd.read_parquet(f"{C.MTB}/rows.parquet", columns=["row_id", "det_id", "series"])
    RS = RS.merge(X, on="row_id")
    import multiprocessing as mp
    stas = a.stations.split(",") if a.stations else C.STATIONS
    sizes = []
    t00 = time.time()
    for ns in stas:
        outs = {"data": f"{C.MWIN}/{ns}.npz", "c2": f"{C.MWIN}/{ns}__c2.npz"}
        if all(os.path.exists(p) for p in outs.values()):
            log(f"{ns}: exists"); continue
        R = RS[RS.netsta == ns]
        tasks = []
        for d, g in R.groupby("date"):
            tasks.append((ns, d, g[g.series == "data"].reset_index(drop=True), g[g.series == "c1"].assign(det_id=-1).reset_index(drop=True)))
        tasks.sort(key=lambda t: -(len(t[2]) + len(t[3])))
        acc = {"data": [], "c2": []}; nerr = 0; t0 = time.time()
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=200) as p:
            for res in p.imap_unordered(task2, tasks, chunksize=1):
                for tag, (_, d, rid, W, src, rr, err) in zip(("data", "c2"), res):
                    if err == "empty":
                        continue
                    if err:
                        nerr += 1; continue
                    acc[tag].append((rid, W, src, rr))
        for tag, L in acc.items():
            if not L:
                continue
            rid = np.concatenate([x[0] for x in L]); W = np.concatenate([x[1] for x in L])
            src = np.concatenate([x[2] for x in L]); rr = np.concatenate([x[3] for x in L])
            np.savez(outs[tag] + ".tmp.npz", row_id=rid, W=W, src=src, rms_raw=rr, tau=S.TAU); os.replace(outs[tag] + ".tmp.npz", outs[tag])
            okz = np.isfinite(W[:, 0, 0]); okh = np.isfinite(W[:, 1, 0]) | np.isfinite(W[:, 2, 0])
            sizes.append(dict(netsta=ns, series=tag, n_windows=len(rid), n_Z_ok=int(okz.sum()), n_H_ok=int(okh.sum()),
                              station_days=len(tasks), bytes=os.path.getsize(outs[tag]), errors=nerr, secs=round(time.time() - t0, 1)))
            log(f"{ns} {tag}: {len(rid)} windows, Z ok {okz.sum()}, H ok {okh.sum()}, {len(tasks)} station-days, errors {nerr}, "
                f"{os.path.getsize(outs[tag]) / 1e6:.0f} MB, {time.time() - t0:.0f}s")
        del acc
        pd.DataFrame(sizes).to_csv(f"{C.MTB}/window_sizes{'_' + a.stations.replace(',', '_') if a.stations else ''}.csv", index=False)
    log(f"windows done {time.time() - t00:.0f}s")
    if not a.stations:
        open(f"{C.OUT}/WA_WINDOWS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
