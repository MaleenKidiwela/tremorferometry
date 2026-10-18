#!/usr/bin/env python
"""s10_tide_prep -- prepare and cache every (family, station, epoch) series for the tide / sub-daily follow-up
(s10_measure's grand stack, S pick and realignment, unchanged). Checks tS against the pilot's series.csv.
    python s10_tide_prep.py --workers 6
"""
import os, sys, argparse, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_tide_common as T

X = RS = None


def _job(ns):
    return ns, T.prepare_station(ns, X, RS)


def main():
    global X, RS
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    T.ensure(); log = S.Log(f"{T.LOGD}/prep.log")
    X = pd.read_parquet(f"{T.PILOT_TB}/rows.parquet")
    RS = pd.read_parquet(f"{T.PILOT_TB}/row_station.parquet")
    stas = sorted(RS.netsta.unique())
    import multiprocessing as mp
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as p:
        for ns, fl in p.imap_unordered(_job, stas):
            log(f"{ns}: {len(fl)} series {time.time() - t0:.0f}s")
    ser = pd.read_csv(f"{T.PILOT_TB}/series.csv")
    rows = []
    for f in sorted(os.listdir(T.CACHE)):
        z = np.load(f"{T.CACHE}/{f}", allow_pickle=False)
        s = ser[(ser.fam == int(z["fam"])) & (ser.netsta == str(z["netsta"])) & (ser.epoch == str(z["epoch"]))]
        rows.append(dict(fam=int(z["fam"]), netsta=str(z["netsta"]), epoch=str(z["epoch"]), n=len(z["row_id"]), n2=len(z["row_id2"]),
                         n_days=len(np.unique(z["date"])), tS=float(z["tS"]),
                         tS_pilot=float(s.tS_pick.iloc[0]) if len(s) else np.nan))
    D = pd.DataFrame(rows); D["dtS"] = D.tS - D.tS_pilot
    D.to_csv(f"{T.TB}/series_prep.csv", index=False)
    log(f"prep done: {len(D)} series, max |tS - pilot| = {np.nanmax(np.abs(D.dtS)):.4f} s")
    open(f"{T.OUT}/PREP_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
