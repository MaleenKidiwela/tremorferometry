#!/usr/bin/env python
"""s12_reflect_inject -- detection threshold of the LVL-base reflection scan (task 2; PLAN §2 calibration column).
A copy of each pair's own direct-S pulse (horizontals, -0.2..+0.4 s about S, cosine-tapered) scaled by A (relative to the
direct S amplitude) is added at the nominal predicted SrS lag (H 3.5 km, Vs 3.0 km/s), with polarity +1 or -1 (both runs),
to the grand stacks; the residual log-envelope and the nominal-model moveout peak are recomputed and compared with the
scrambled control (lags permuted within station, 500 draws) of the injected data. Output tables/refl_inject.csv.
Physical scale: SH reflection coefficient of a sharp LVL base (Vs 3.0 -> 3.9 km/s, rho 2.9 -> 3.0) ~0.15 at normal
incidence; ~0.13 x spreading at the observed lags.
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s12_reflect as RF

TB = RF.TB


def main():
    tag = os.environ.get("S12_TAG", "")
    if tag:
        RF.TAG = tag; RF.SMN, RF.CEN, RF.FL0, RF.FL1 = 13, 0.15, 0.45, 0.9
    z = np.load(f"{TB}/refl_traces{tag}.npz")
    T = z["TREL"]; X0 = z["X"].astype(np.float64); st = z["netsta"]; w = np.sqrt(z["n"])
    HG, VG = RF.HG, RF.VG
    lag = z["lag_SrS_b"][:, HG.index(3.5), VG.index(3.0)]
    win = (T >= -0.2) & (T <= 0.4) if not tag else (T >= -0.4) & (T <= 0.8)
    nw = win.sum(); tap = np.hanning(nw)
    groups = {s: np.flatnonzero(st == s) for s in np.unique(st)}
    rng = np.random.default_rng(7)
    rows = []
    for pol in (1, -1):
        for A in (0.0, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25):
            RH = []
            for i in range(len(X0)):
                X = X0[i].copy()
                k0 = int(round((lag[i] - (0.4 if tag else 0.2) - T[0]) * RF.FS))
                for c in (1, 2):
                    if np.isfinite(X[c]).all() and k0 + nw < len(T):
                        X[c, k0:k0 + nw] += pol * A * tap * X0[i, c][win]
                EZ, EH = RF.envelopes(X)
                RH.append(RF.resid_logenv(EH))
            RH = np.array(RH)
            M = RF.pair_metric(RH)
            ix = np.clip(np.round((lag - T[0]) * RF.FS).astype(int), 0, len(T) - 1)
            def val(perm):
                v = M[np.arange(len(M)), ix[perm]]; ok = np.isfinite(v)
                return float(np.sum(w[ok] * v[ok]) / np.sum(w[ok]))
            pk = val(np.arange(len(M)))
            sc = []
            for k in range(500):
                perm = np.arange(len(M))
                for s_, idx in groups.items():
                    perm[idx] = rng.permutation(idx)
                sc.append(val(perm))
            sc = np.array(sc)
            r = dict(polarity=pol, A=A, peak=pk, scr_q95=float(np.quantile(sc, .95)), scr_q99=float(np.quantile(sc, .99)),
                     p=float((np.sum(sc >= pk) + 1) / 501))
            rows.append(r); print(r, flush=True)
    pd.DataFrame(rows).to_csv(f"{TB}/refl_inject{tag}.csv", index=False)
    open(f"{RF.OUT}/REFLECT_INJECT{tag.upper()}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
