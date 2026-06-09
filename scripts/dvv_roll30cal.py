#!/usr/bin/env python
"""Per-family dv/v from a TRUE 30-CALENDAR-DAY trailing rolling stack (not count-based).
For each date that has a daily stack, stack ONLY the daily stacks whose date is within the 30
calendar days ENDING on that date (require >= MIN_STACKS in the window, else the date gets no point).
Then SVD-Wiener filter the rolling-stack matrix and stretch each on the coda window vs the filtered
all-time reference.  This drops the gap-spanning "stitched" values the count-based 30-stack window made.
Usage: dvv_roll30cal.py --station PGC --npz data/long_window_daily_PGC_despike.npz --window 2.0 4.0 --out ...
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402

_C = {}
def _init(t, fs, w1, w2, days, minst, mindays, neps):
    _C.update(t=t, fs=fs, tmin=w1 - t[0], tmax=w2 - t[0], days=int(days), minst=int(minst),
              mindays=int(mindays), neps=int(neps), prenoise=(t < -1.0))

def _family(arg):
    patch, S, dstr = arg
    if len(S) < _C["mindays"]:
        return []
    order = np.argsort(dstr); S = S[order].astype(np.float64); dstr = np.asarray(dstr)[order]
    Dn = pd.to_datetime(dstr).values.astype('datetime64[D]')
    n = len(S); idx = np.arange(n)
    cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
    los = np.searchsorted(Dn, Dn - np.timedelta64(_C["days"], 'D'), side='right')  # first idx within 30 cal days
    counts = idx + 1 - los
    keep = counts >= _C["minst"]
    if int(keep.sum()) < 8:
        return []
    R = (cs[idx + 1] - cs[los]) / counts[:, None]          # mean of the in-window daily stacks
    R = R[keep]; Dr = dstr[keep]; ck = counts[keep]
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    ne, ns = Rn.shape; sig = Rn[:, _C["prenoise"]].std()
    try:
        U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
        nmf = sig * (np.sqrt(ne) + np.sqrt(ns)); w = sv ** 2 / (sv ** 2 + nmf ** 2); Rf = (U * (sv * w)) @ Vt
    except np.linalg.LinAlgError:
        Rf = Rn
    ref = Rf.mean(0)
    rows = []
    for i in range(ne):
        try:
            r = stretch_dvv(ref, Rf[i], fs=_C["fs"], t_min=_C["tmin"], t_max=_C["tmax"], eps_max=0.02, n_eps=_C["neps"])
            rows.append((patch, Dr[i], float(r.dvv), float(r.cc_max), int(ck[i])))
        except Exception:
            pass
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--station", required=True); p.add_argument("--npz", required=True)
    p.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    p.add_argument("--days", type=int, default=30); p.add_argument("--min-stacks", type=int, default=5)
    p.add_argument("--mindays", type=int, default=60); p.add_argument("--n-eps", type=int, default=201)
    p.add_argument("--workers", type=int, default=28); p.add_argument("--out", required=True)
    args = p.parse_args()
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    d = np.load(args.npz, allow_pickle=True)
    t = d["t"]; fs = float(d["fs"]); stacks = d["stacks"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).strftime("%Y-%m-%d").values
    uniq = pd.unique(pat)
    print(f"[{args.station}] {len(stacks)} daily stacks, {len(uniq)} families; window {args.window[0]}-{args.window[1]}s, "
          f"{args.days}-cal-day trailing, min {args.min_stacks} stacks", flush=True)
    tasks = [(fam, stacks[pat == fam], dates[pat == fam]) for fam in uniq]
    ctx = mp.get_context("spawn"); out = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(t, fs, args.window[0], args.window[1], args.days, args.min_stacks,
                                       args.mindays, args.n_eps)) as ex:
        futs = [ex.submit(_family, tk) for tk in tasks]
        for f in as_completed(futs):
            out.extend(f.result())
    df = pd.DataFrame(out, columns=["patch", "date", "dvv", "cc_max", "n_stack"])
    df.to_csv(args.out, index=False)
    if len(df):
        print(f"[{args.station}] DONE {len(df):,} rows, {df.patch.nunique()} families, mean cc {df.cc_max.mean():.3f}, "
              f"median window {df.n_stack.median():.0f} stacks -> {args.out}", flush=True)
    else:
        print(f"[{args.station}] EMPTY -> {args.out}", flush=True)

if __name__ == "__main__":
    main()
