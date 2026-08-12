#!/usr/bin/env python
"""Re-measure per-family dv/v in the 2-4 s coda window using the method we converged on
for B927 c660 (off the diluting direct-S, inside the coherent coda):

  per family:  (1) 30-day rolling WAVEFORM stack (step 1 day; count-based, cumsum)
               (2) SVD-Wiener filter the rolling-stack matrix (random-matrix noise floor)
               (3) stretch each filtered rolling stack vs the filtered all-time reference on 2-4 s

No medianing of dv/v values anywhere -- the stacking IS the smoothing.  Parallel across families.
Usage:  dvv_2to4_roll30svd.py --station B927 --npz data/long_window_daily_B927.npz --out data/daily_dvv_B927_2to4_roll.csv
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402

_C = {}
def _init(t, fs, w1, w2, W, mindays, neps):
    _C.update(t=t, fs=fs, tmin=w1 - t[0], tmax=w2 - t[0], W=W, mindays=mindays, neps=neps,
              prenoise=(t < -1.0))

def _family(arg):
    patch, S, dstr = arg
    t = _C["t"]; W = _C["W"]
    if len(S) < _C["mindays"]:
        return []
    order = np.argsort(dstr); S = S[order].astype(np.float64); dstr = np.asarray(dstr)[order]
    # (1) 30-stack rolling waveform stack via cumsum  -> rows end at index W-1..n-1
    cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
    R = (cs[W:] - cs[:-W]) / W
    Dr = dstr[W - 1:]
    if len(R) < 8:
        return []
    # (2) SVD-Wiener filter (row-normalised; random-matrix noise floor)
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    ne, ns = Rn.shape
    sig = Rn[:, _C["prenoise"]].std()
    try:
        U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
        nmf = sig * (np.sqrt(ne) + np.sqrt(ns))
        w = sv ** 2 / (sv ** 2 + nmf ** 2)
        Rf = (U * (sv * w)) @ Vt
    except np.linalg.LinAlgError:
        Rf = Rn
    ref = Rf.mean(0)
    # (3) stretch each filtered rolling stack on 2-4 s
    rows = []
    for i in range(ne):
        try:
            r = stretch_dvv(ref, Rf[i], fs=_C["fs"], t_min=_C["tmin"], t_max=_C["tmax"],
                            eps_max=0.02, n_eps=_C["neps"])
            rows.append((patch, Dr[i], float(r.dvv), float(r.cc_max)))
        except Exception:
            pass
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--station", required=True)
    p.add_argument("--npz", required=True)
    p.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    p.add_argument("--roll", type=int, default=30)
    p.add_argument("--mindays", type=int, default=60)
    p.add_argument("--n-eps", type=int, default=201)
    p.add_argument("--workers", type=int, default=28)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    d = np.load(args.npz, allow_pickle=True)
    t = d["t"]; fs = float(d["fs"]); stacks = d["stacks"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).strftime("%Y-%m-%d").values
    uniq = pd.unique(pat)
    print(f"[{args.station}] {len(stacks)} daily stacks, {len(uniq)} families; "
          f"window {args.window[0]}-{args.window[1]}s, roll {args.roll}", flush=True)
    tasks = [(fam, stacks[pat == fam], dates[pat == fam]) for fam in uniq]
    ctx = mp.get_context("spawn")
    out = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(t, fs, args.window[0], args.window[1], args.roll,
                                       args.mindays, args.n_eps)) as ex:
        futs = [ex.submit(_family, tk) for tk in tasks]
        done = 0
        for f in as_completed(futs):
            out.extend(f.result()); done += 1
            if done % 20 == 0:
                print(f"  [{args.station}] {done}/{len(tasks)} families, {len(out):,} rows", flush=True)
    df = pd.DataFrame(out, columns=["patch", "date", "dvv", "cc_max"])
    df.to_csv(args.out, index=False)
    if len(df):
        nf = df.patch.nunique()
        print(f"[{args.station}] DONE {len(df):,} rows, {nf} families, mean cc {df.cc_max.mean():.3f}, "
              f"dv/v std {df.dvv.std()*100:.3f}%  -> {args.out}", flush=True)
    else:
        print(f"[{args.station}] DONE but EMPTY -> {args.out}", flush=True)

if __name__ == "__main__":
    main()
