#!/usr/bin/env python
"""s10_multi_windows -- S10 multi step 2 (multi/PLAN.md §0-1): windows for the new families with the pilot reader
(s10_windows.task, unchanged: cut store where the row is in it, raw day files otherwise; 2-8 Hz, 40 Hz, despike, RMS).
Two passes per station, as in the pilot: data rows (P >= 0.5), and scrambled-control c1 rows read at frame time
t + offset (c2). Pilot families 10/87/17 are reused from the pilot windows by row key. One station at a time.

    python s10_multi_windows.py --workers 12 [--stations NET.STA,...]
"""
from __future__ import annotations
import os, sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
import s10_windows as SW


def run_pass(ns, RS, workers, log, tag):
    tasks = [(ns, d, g.reset_index(drop=True)) for d, g in RS.groupby("date")]
    tasks.sort(key=lambda t: -len(t[2]))
    log(f"{ns} {tag}: {len(tasks)} station-days, {len(RS)} windows")
    acc = []; nerr = 0; t0 = time.time()
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers) as p:
        for k, (_, d, rid, W, src, rr, err) in enumerate(p.imap_unordered(SW.task, tasks, chunksize=1)):
            if err:
                nerr += 1; log(f"ERROR {ns} {d}: {err}"); continue
            acc.append((rid, W, src, rr))
            if k % 500 == 0:
                log(f"  {ns} {tag} {k + 1}/{len(tasks)} {time.time() - t0:.0f}s")
    if not acc:
        return None
    return (np.concatenate([x[0] for x in acc]), np.concatenate([x[1] for x in acc]),
            np.concatenate([x[2] for x in acc]), np.concatenate([x[3] for x in acc]), nerr)


def reuse(ns, X, sfx, series):
    """pilot windows of families 10/87/17 (pilot <ns>.npz data rows or <ns>__c1f.npz) -> new row ids"""
    p = f"{M.PILOT}/win/{ns}{sfx}.npz"
    if not os.path.exists(p):
        return None
    z = np.load(p)
    mp_ = X[(X.pilot_row_id >= 0) & (X.series == series)].set_index("pilot_row_id").row_id
    keep = np.isin(z["row_id"], mp_.index.to_numpy())
    rid = mp_.reindex(z["row_id"][keep]).to_numpy()
    return rid, z["W"][keep], z["src"][keep], z["rms_raw"][keep]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stations", default="")
    a = ap.parse_args()
    M.ensure()
    log = S.Log(f"{M.LOGD}/windows.log")
    RS = pd.read_parquet(f"{M.TB}/row_station.parquet")
    X = pd.read_parquet(f"{M.TB}/rows.parquet", columns=["row_id", "det_id", "series", "fam", "pilot_row_id"])
    RS = RS.merge(X[["row_id", "det_id", "series", "fam"]], on="row_id")
    stas = a.stations.split(",") if a.stations else M.STATIONS
    for ns in stas:
        for tag, series, sfx_out, sfx_pilot in (("data", "data", "", ""), ("c2", "c1", "__c2", "__c1f")):
            out = f"{M.WIN}/{ns}{sfx_out}.npz"
            if os.path.exists(out):
                log(f"{ns} {tag}: exists, skipped"); continue
            R = RS[(RS.netsta == ns) & (RS.series == series) & ~RS.fam.isin(M.PILOT_FAMS)].drop(columns="fam").copy()
            if tag == "c2":
                R["det_id"] = -1                     # pilot Rev. 1: c2 windows always from the raw day files
            res = run_pass(ns, R, a.workers, log, tag) if len(R) else None
            parts = []
            if res is not None:
                parts.append(res[:4])
            ru = reuse(ns, X, sfx_pilot, series)
            if ru is not None:
                parts.append(ru)
            if not parts:
                continue
            rid = np.concatenate([p[0] for p in parts]); W = np.concatenate([p[1] for p in parts])
            src = np.concatenate([p[2] for p in parts]); rr = np.concatenate([p[3] for p in parts])
            np.savez(out + ".tmp.npz", row_id=rid, W=W, src=src, rms_raw=rr, tau=S.TAU)
            os.replace(out + ".tmp.npz", out)
            log(f"{ns} {tag}: saved {len(rid)} windows (new {0 if res is None else len(res[0])}, reused "
                f"{0 if ru is None else len(ru[0])}), Z ok {np.isfinite(W[:, 0, 0]).sum()}, cut {(src == 1).sum()}, raw {(src == 2).sum()}"
                + ("" if res is None else f", station-day errors {res[4]}"))
            del W, parts, res
    open(f"{M.OUT}/WINDOWS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
