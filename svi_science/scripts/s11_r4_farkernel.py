#!/usr/bin/env python
"""s11_r4_farkernel -- S11 round 4, stage S2 (PLAN.md §8.0): far-pair kernels for the round-3 variants ell10 and two3_20,
plus an ell40 re-run (reproduction check against distant/<region>/tables/kernel_coef.parquet). Unchanged s11_kernel_mc physics
(VARIANTS replaced by these three; same job order and seeds as s11_dist_kernel -> common random numbers), lapse weights from the
far P store grand stacks (as s11_dist_kernel). 3C, every window. Stop rules: |sum a - 1| <= 0.02, a_d < 0 for <= 10 % of pairs
per variant, ell40 re-run within 1e-6 of the stored far ell40.
    S11_REGION=... python s11_r4_farkernel.py [--workers 12] [--npart 100000] [--pairs N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_dist_common as D
import s11_kernel_mc as K
import s11_dist_kernel as DK
import s11_r4_common as R4

K.VARIANTS = {"ell40": (40.0, 40.0, 0), "ell10": (10.0, 10.0, 0), "two3_20": (3.0, 20.0, 0)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--npart", type=int, default=100000); ap.add_argument("--pairs", type=int, default=0)
    a = ap.parse_args()
    lg = R4.log("farkernel"); t0 = time.time()
    NT = D.near_tables(); F = NT["F"]
    ser = pd.read_csv(f"{D.DTB}/series.csv"); ser = ser[(ser.n_data > 0) & (ser.causal25 > 1.5)]
    pairs = sorted(set(zip(ser.fam.astype(int), ser.netsta)))
    XY = D.station_coords(); dep = D.sensor_depths(sorted({p[1] for p in pairs}), lg)
    lat0, lon0 = float(F.lat.mean()), float(F.lon.mean())
    pj = lambda la, lo: np.array([(lo - lon0) * 111.19 * np.cos(np.radians(lat0)), (la - lat0) * 111.19])
    Fx = {int(r.fam): np.r_[pj(r.lat, r.lon), r.depth_km] for r in F.itertuples()}
    Sx = {ns: np.r_[pj(XY.lat[ns], XY.lon[ns]), dep[ns] / 1000.0] for ns in {p[1] for p in pairs}}
    jobs = [(2000 + i, f, ns, Fx[f], Sx[ns], a.npart) for i, (f, ns) in enumerate(pairs)]     # = s11_dist_kernel job order/seeds
    if a.pairs:
        jobs = jobs[:a.pairs]
    lg(f"region {C.REGION}: {len(jobs)} far pairs x {list(K.VARIANTS)}, {a.npart} particles, {a.workers} workers")
    import multiprocessing as mp
    res = {}
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (i, fam, ns, out) in enumerate(p.imap_unordered(DK._pair, jobs)):
            res[(fam, ns)] = out
            if (k + 1) % 25 == 0 or k + 1 == len(jobs):
                lg(f"  {k + 1}/{len(jobs)} pairs {time.time() - t0:.0f}s")
    W = DK.lapse_weights("")
    rows = []
    for key, Eb in W.items():
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in res:
            continue
        for v in K.VARIANTS:
            k = res[(fam, ns)][v]
            for wn, win in C.WINS.items():
                wb = K.window_weights(Eb, win, "3c")
                if wb.sum() <= 0:
                    continue
                av = K.coef(k, wb)
                rows.append(dict(series=key, fam=fam, netsta=ns, epoch=ep, variant=v, win=wn, comp="3c", unit_sum=float(np.nansum(av)),
                                 a_s=float(C.group_coef(av, "s")), a_m=float(C.group_coef(av, "m")), a_d=float(C.group_coef(av, "d")),
                                 a_b=float(C.group_coef(av, "b"))))
    KC = pd.DataFrame(rows)
    tag = "_test" if a.pairs else ""
    KC.to_parquet(f"{R4.TB4}/kernel_coef_far_r4{tag}.parquet", index=False)
    old = pd.read_parquet(f"{D.DTB}/kernel_coef.parquet"); old = old[(old.variant == "ell40") & (old.comp == "3c")][["series", "win", "a_s", "a_d"]]
    m = KC[KC.variant == "ell40"].merge(old, on=["series", "win"], suffixes=("", "_old"))
    rep = float(np.nanmax(np.abs(np.r_[m.a_s - m.a_s_old, m.a_d - m.a_d_old]))) if len(m) else np.nan
    ut = float((KC.unit_sum - 1).abs().max())
    prim = KC[KC.win == C.HEAD].groupby(["variant", "fam", "netsta"]).a_d.first().reset_index()
    neg = prim.groupby("variant").a_d.apply(lambda x: float((x < 0).mean())).to_dict()
    med = KC[KC.win == C.HEAD].groupby("variant")[["a_s", "a_m", "a_d", "a_b"]].median().round(3).to_dict("index")
    lg(f"unit test max |sum-1| {ut:.2e}; a_d<0 fraction {neg}; ell40 reproduction max |diff| {rep:.2e} over {len(m)} rows; medians 2-5 s {med}")
    stop = (ut > 0.02) or any(v > 0.10 for v in neg.values()) or not (rep <= 1e-6)
    if not a.pairs:
        mk = f"{R4.R4}/FARKERNEL_{C.REGION}_{'STOP' if stop else 'DONE'}"
        open(mk, "w").write(f"{R4.pt_now()}\nunit {ut:.3e} neg {neg} repro {rep:.3e}\n")
        lg(f"wrote {mk} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
