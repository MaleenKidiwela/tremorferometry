#!/usr/bin/env python
"""Equivalence + speed test: ab_grow_fast.matched_filter vs ab_grow.matched_filter.

  nice -n 10 python scripts/ab_grow_fast_equiv.py [--days D1,D2,D3] [--every 9] [--speed-days 8]
Part 1: ~300 templates (every --every-th of NW tpl_it0) on 3 NW scan days: detection lists
(grp, t) identical? max |network function| difference over all samples (NaN pattern must match).
Part 2: all templates on --speed-days consecutive days, wall time per day for both."""
import argparse, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import ab_common as A, ab_grow as G, ab_grow_fast as F

POOL = "STW,B007,OPC,B04A,LRIV,B001,BLN,BS11,W020,B943,B013,W030,W040,DOSE,C04A,HDW,GNW,SQM,PL11,SMW,N050,N060,N080"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", default=f"{A.RES}/NW/clusters/tpl_it0.npz")
    ap.add_argument("--days", default="2008-05-08,2008-05-12,2008-05-16")
    ap.add_argument("--every", type=int, default=9)
    ap.add_argument("--speed-days", type=int, default=8)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    pool = POOL.split(",")
    TPLall, _ = G.cut_templates(G.load_stacks(a.tpl))
    keys = sorted(TPLall)[::a.every]
    TPL = {g: TPLall[g] for g in keys}
    days = a.days.split(",")
    A.log(f"[equiv] {len(TPL)} templates x {days}")
    d1, d2 = {}, {}
    t = time.time(); r1 = G.matched_filter(A.DayCache(pool), TPL, days, _dump=d1); tr = time.time() - t
    t = time.time(); r2 = F.matched_filter(A.DayCache(pool), TPL, days, _dump=d2, prep_workers=a.workers); tf = time.time() - t
    k1 = set(zip(r1.grp, np.round(r1.t * A.FS).astype(np.int64)))
    k2 = set(zip(r2.grp, np.round(r2.t * A.FS).astype(np.int64)))
    m = r1.merge(r2, on=["grp", "t", "date"], suffixes=("_r", "_f"))
    dmax, nanmis = 0.0, 0
    for k in d1:
        x, y = d1[k], d2[k]
        nanmis += int((np.isnan(x) != np.isnan(y)).sum())
        ok = ~np.isnan(x) & ~np.isnan(y)
        if ok.any():
            dmax = max(dmax, float(np.abs(x[ok] - y[ok]).max()))
    res = dict(n_tpl=len(TPL), days=days, n_det_ref=len(r1), n_det_fast=len(r2),
               det_identical=bool(k1 == k2), only_ref=len(k1 - k2), only_fast=len(k2 - k1),
               fn_maxabs=dmax, nan_mismatch=nanmis,
               cc_maxabs=float((m.cc_r - m.cc_f).abs().max()) if len(m) else np.nan,
               nmad_maxabs=float((m.nmad_r - m.nmad_f).abs().max()) if len(m) else np.nan,
               nsta_equal=bool((m.n_sta_r == m.n_sta_f).all()), nch_equal=bool((m.n_ch_r == m.n_ch_f).all()),
               thr_maxabs=float((m.thr_r - m.thr_f).abs().max()) if len(m) else np.nan,
               t_ref_s=tr, t_fast_s=tf)
    print("EQUIV", res, flush=True)
    del d1, d2
    if a.speed_days:
        sd = [d for d in A.scan_days("NW") if d >= "2008-05-01"][:a.speed_days]
        t = time.time(); s1 = G.matched_filter(A.DayCache(pool), TPLall, sd); tr = time.time() - t
        prof = []
        t = time.time(); s2 = F.matched_filter(A.DayCache(pool), TPLall, sd, prep_workers=a.workers, _prof=prof); tf = time.time() - t
        k1 = set(zip(s1.grp, np.round(s1.t * A.FS).astype(np.int64)))
        k2 = set(zip(s2.grp, np.round(s2.t * A.FS).astype(np.int64)))
        print("SPEED", dict(n_tpl=len(TPLall), days=len(sd), ref_s_per_day=tr / len(sd),
                            fast_s_per_day=tf / len(sd), speedup=tr / tf, det_identical=bool(k1 == k2),
                            n_det=len(s1), n_det_fast=len(s2)), flush=True)
        print(pd.DataFrame(prof).to_string(), flush=True)



def restack_check(n_dates=4, workers=8):
    """ab_grow.restack vs ab_grow_fast.restack (parallel) on n_dates of NW grow_qc det_it0."""
    pool = POOL.split(",")
    det = pd.read_parquet(f"{A.RES}/NW/grow_qc/det_it0.parquet")
    ds = sorted(det.date.unique())[::max(1, len(det.date.unique()) // n_dates)][:n_dates]
    det = det[det.date.isin(ds)]
    t = time.time(); S1 = G.restack(A.DayCache(pool), pool, det); t1 = time.time() - t
    t = time.time(); S2 = F.restack(A.DayCache(pool), pool, det, workers=workers); t2 = time.time() - t
    k1 = {(g, s) for g in S1 for s in S1[g]}; k2 = {(g, s) for g in S2 for s in S2[g]}
    d = max(float(np.abs(S1[g][s] - S2[g][s]).max()) for g, s in k1 & k2)
    rel = max(float(np.abs(S1[g][s] - S2[g][s]).max() / max(np.abs(S1[g][s]).max(), 1e-30)) for g, s in k1 & k2)
    TA, _ = G.cut_templates(S1); TB, _ = G.cut_templates(S2)
    same_sel = all(set(TA[g]) == set(TB.get(g, {})) for g in TA) and set(TA) == set(TB)
    print("RESTACK", dict(dates=ds, n_det=len(det), keys_equal=k1 == k2, n_keys=len(k1), maxabs=d,
                          maxrel=rel, same_template_stations=same_sel, t_ref_s=t1, t_fast_s=t2), flush=True)


if __name__ == "__main__":
    if "--restack" in sys.argv:
        restack_check()
    else:
        main()
