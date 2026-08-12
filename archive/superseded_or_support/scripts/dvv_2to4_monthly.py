#!/usr/bin/env python
"""Per-family 2-4 s dv/v from MONTHLY WAVEFORM STACKS (no median of dv/v anywhere).
For each family: stack (mean WAVEFORM) of each calendar month's daily stacks -> SVD-Wiener the
monthly-stack matrix -> stretch each monthly stack vs the filtered all-time reference on 2-4 s.
One dv/v per (family, month), measured from stacked waveforms. Parallel across families.
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
def _init(t, fs, w1, w2, mindays, neps):
    _C.update(t=t, fs=fs, tmin=w1 - t[0], tmax=w2 - t[0], mindays=mindays, neps=neps, prenoise=(t < -1.0))

def _family(arg):
    patch, S, dstr = arg
    if len(S) < _C["mindays"]:
        return []
    ym = pd.to_datetime(np.asarray(dstr)).to_period('M').astype(str)
    df = pd.DataFrame(S.astype(np.float64)); df['ym'] = ym
    M = df.groupby('ym').mean()                       # MEAN WAVEFORM per month = monthly stack (NOT a dv/v median)
    nday = df.groupby('ym').size()
    months = M.index.values; R = M.values
    if len(R) < 8:
        return []
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    ne, ns = Rn.shape; sig = Rn[:, _C["prenoise"]].std()
    try:
        U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
        nmf = sig * (np.sqrt(ne) + np.sqrt(ns)); w = sv ** 2 / (sv ** 2 + nmf ** 2)
        Rf = (U * (sv * w)) @ Vt
    except np.linalg.LinAlgError:
        Rf = Rn
    ref = Rf.mean(0)
    rows = []
    for i in range(ne):
        try:
            r = stretch_dvv(ref, Rf[i], fs=_C["fs"], t_min=_C["tmin"], t_max=_C["tmax"], eps_max=0.02, n_eps=_C["neps"])
            rows.append((patch, str(months[i]) + '-15', float(r.dvv), float(r.cc_max), int(nday.values[i])))
        except Exception:
            pass
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--station", required=True); p.add_argument("--npz", required=True)
    p.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    p.add_argument("--mindays", type=int, default=60); p.add_argument("--n-eps", type=int, default=201)
    p.add_argument("--workers", type=int, default=28); p.add_argument("--out", required=True)
    args = p.parse_args()
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    d = np.load(args.npz, allow_pickle=True)
    t = d["t"]; fs = float(d["fs"]); stacks = d["stacks"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).strftime("%Y-%m-%d").values
    uniq = pd.unique(pat)
    print(f"[{args.station}] {len(stacks)} daily stacks, {len(uniq)} families; monthly-stack {args.window[0]}-{args.window[1]}s", flush=True)
    tasks = [(fam, stacks[pat == fam], dates[pat == fam]) for fam in uniq]
    ctx = mp.get_context("spawn"); out = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(t, fs, args.window[0], args.window[1], args.mindays, args.n_eps)) as ex:
        futs = [ex.submit(_family, tk) for tk in tasks]
        for f in as_completed(futs):
            out.extend(f.result())
    df = pd.DataFrame(out, columns=["patch", "date", "dvv", "cc_max", "n_day"])
    df.to_csv(args.out, index=False)
    print(f"[{args.station}] DONE {len(df):,} monthly rows, {df.patch.nunique() if len(df) else 0} families, "
          f"mean cc {df.cc_max.mean():.3f} -> {args.out}" if len(df) else f"[{args.station}] EMPTY", flush=True)

if __name__ == "__main__":
    main()
