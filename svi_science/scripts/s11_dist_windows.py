#!/usr/bin/env python
"""s11_dist_windows -- S11 distant-station lever, step 2a (distant/PLAN.md §4): windows of every data and c1 row of the
usable far pairs, at the corrected S (predicted + the pair's sample-stack pick shift), on the wide grid -28..+10 s, with
the unchanged S10 processing (s10_windows.task; data and c1 rows as separate calls = same despike reference as the near
runs). Writes distant/<region>/tables/{far_rows.parquet, far_offsets.csv, win/<NET.STA>[__c2].npz, window_sizes.csv}.
    S11_REGION=SVI python s11_dist_windows.py --workers 10
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
import s11_dist_common as D
D.set_grid(D.FAR_GRID)
import s10_windows as SW
import s11_dist_inventory as INV          # DistRaw reader (sets SW.Raw); its grid switch is undone below
D.set_grid(D.FAR_GRID)


def task2(arg):
    ns, d, Rd, Rc = arg
    out = []
    for R in (Rd, Rc):
        out.append(SW.task((ns, d, R)) if len(R) else (ns, d, np.zeros(0, np.int64), None, None, None, "empty"))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10); ap.add_argument("--stations", default="")
    a = ap.parse_args()
    D.ensure(); os.makedirs(f"{D.DTB}/win", exist_ok=True)
    log = C.Log(f"{D.DLOG}/windows.log"); t00 = time.time()
    FP = pd.read_csv(f"{D.DTB}/far_pairs.csv")
    U = FP[FP.usable].copy()
    U["offset"] = U.off_pred + U.pick_shift
    U[["fam", "netsta", "offset", "sp", "off_pred", "pick_shift", "epi_km", "hypo_km", "pick_snr", "causal3c", "cc20"]].to_csv(
        f"{D.DTB}/far_offsets.csv", index=False)
    AV = pd.read_parquet(f"{D.DTB}/avail.parquet")
    INV._RTS = D.roots(); INV._AV = {(r.netsta, r.date): r.path for r in AV.itertuples()}
    NT = D.near_tables(); X = NT["X"]
    dead = D.dead_days()
    offk = U.set_index(["fam", "netsta"]).offset
    rows = []
    for ns in sorted(U.netsta.unique()):
        fams = set(U.fam[U.netsta == ns])
        R = X[X.fam.isin(fams)]
        R = R[[(ns, d) in INV._AV and (ns, d) not in dead for d in R.date]]
        ep = D.epoch_of(D.epochs_table(ns), R.date)
        R = R.assign(netsta=ns, epoch=ep, tS=R.t.to_numpy() + np.array([offk[(f, ns)] for f in R.fam]))
        rows.append(R[R.epoch != ""][["row_id", "netsta", "epoch", "tS", "date", "fam", "series", "P"]])
    RS = pd.concat(rows, ignore_index=True)
    RS.to_parquet(f"{D.DTB}/far_rows.parquet", index=False)
    inv = RS.groupby(["netsta", "series"]).size().unstack(fill_value=0)
    log(f"usable pairs {len(U)} at {U.netsta.nunique()} stations; row-stations {len(RS)}; station-days "
        f"{RS.groupby('netsta').date.nunique().sum()}\n{inv.to_string()}")
    import multiprocessing as mp
    stas = a.stations.split(",") if a.stations else sorted(RS.netsta.unique())
    sizes = []
    for ns in stas:
        outs = {"data": f"{D.DTB}/win/{ns}.npz", "c2": f"{D.DTB}/win/{ns}__c2.npz"}
        if all(os.path.exists(p) for p in outs.values()):
            log(f"{ns}: exists"); continue
        R = RS[RS.netsta == ns].assign(det_id=-1)
        tasks = [(ns, d, g[g.series == "data"].reset_index(drop=True), g[g.series == "c1"].reset_index(drop=True))
                 for d, g in R.groupby("date")]
        tasks.sort(key=lambda t: -(len(t[2]) + len(t[3])))
        acc = {"data": [], "c2": []}; nerr = 0; t0 = time.time()
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=100) as p:
            for res in p.imap_unordered(task2, tasks, chunksize=1):
                for tag, (_, d, rid, W, src, rr, err) in zip(("data", "c2"), res):
                    if err == "empty":
                        continue
                    if err:
                        nerr += 1; continue
                    acc[tag].append((rid, W, src, rr))
        for tag, L in acc.items():
            if not L:
                np.savez(outs[tag], row_id=np.zeros(0, np.int64), W=np.zeros((0, 3, S.NW), np.float32), src=np.zeros(0, np.int8),
                         rms_raw=np.zeros((0, 3), np.float32), tau=S.TAU)
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
        pd.DataFrame(sizes).to_csv(f"{D.DTB}/window_sizes{'_' + a.stations.replace(',', '_') if a.stations else ''}.csv", index=False)
    log(f"windows done {time.time() - t00:.0f}s")
    if not a.stations:
        open(f"{D.DOUT}/WINDOWS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
