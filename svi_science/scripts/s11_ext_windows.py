#!/usr/bin/env python
"""s11_ext_windows -- S11 Revision 2: windows for the 37 added SVI families with the unchanged S10 reader
(s10_windows.task: cut store where the row is in it, raw day files otherwise; 2-8 Hz, 40 Hz, despike, RMS), same two
passes as s10_multi_windows (data rows; c1 rows read at frame time t + offset = c2). Writes tables/ext/win/<NET.STA>[__c2].npz.
    python s11_ext_windows.py --workers 12
"""
from __future__ import annotations
import os, sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_windows as SW
import s11_common as C

EXT = f"{C.TB}/ext"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    os.makedirs(f"{EXT}/win", exist_ok=True)
    log = C.Log(f"{C.LOGD}/ext_windows.log")
    RS = pd.read_parquet(f"{EXT}/row_station.parquet")
    X = pd.read_parquet(f"{EXT}/rows.parquet", columns=["row_id", "det_id", "series"])
    RS = RS.merge(X, on="row_id")
    import multiprocessing as mp
    for ns in C.STATIONS:
        for tag, series, sfx in (("data", "data", ""), ("c2", "c1", "__c2")):
            out = f"{EXT}/win/{ns}{sfx}.npz"
            if os.path.exists(out):
                log(f"{ns} {tag}: exists"); continue
            R = RS[(RS.netsta == ns) & (RS.series == series)].copy()
            if tag == "c2":
                R["det_id"] = -1
            tasks = [(ns, d, g.reset_index(drop=True)) for d, g in R.groupby("date")]
            tasks.sort(key=lambda t: -len(t[2]))
            acc = []; nerr = 0; t0 = time.time()
            with mp.get_context("fork").Pool(a.workers) as p:
                for k, (_, d, rid, W, src, rr, err) in enumerate(p.imap_unordered(SW.task, tasks, chunksize=1)):
                    if err:
                        nerr += 1; continue
                    acc.append((rid, W, src, rr))
            if not acc:
                continue
            rid = np.concatenate([x[0] for x in acc]); W = np.concatenate([x[1] for x in acc])
            src = np.concatenate([x[2] for x in acc]); rr = np.concatenate([x[3] for x in acc])
            np.savez(out + ".tmp.npz", row_id=rid, W=W, src=src, rms_raw=rr, tau=S.TAU); os.replace(out + ".tmp.npz", out)
            log(f"{ns} {tag}: {len(rid)} windows, Z ok {np.isfinite(W[:, 0, 0]).sum()}, cut {(src == 1).sum()}, raw "
                f"{(src == 2).sum()}, station-day errors {nerr}, {time.time() - t0:.0f}s")
    open(f"{C.OUT}/EXT_WINDOWS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
