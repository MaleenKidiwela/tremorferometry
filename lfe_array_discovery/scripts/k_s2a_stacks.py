#!/usr/bin/env python
"""Step 2a: corrected-time STACKS at every 2005 station, and the EMPIRICAL S offset.

The frozen chain cuts at `t_origin + fam_tt` (Crosson 1-D + a5 corrections).  STEP1CP_RESULTS a.
measured that this prediction is off by a median 2.8 s (p90 6.0 s).  With the Step-1 corrected
times the stacks are coherent, so the true S offset can simply be READ OFF the stack: that is what
every later cut in this study uses instead of fam_tt.

Writes results/repro/s2a_stacks.npz (per family/station 3-C stack, 30 s long, starting at
t_corr + fam_tt - 10 s) and s2a_offsets.csv (envelope-peak offset, MAD height, coherence).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

FS = K.FS
PRE = 10.0
LEN = 30.0


def envelope(x, smooth_s=0.4, fs=FS):
    from scipy.signal import hilbert
    e = np.sqrt((np.abs(hilbert(x, axis=-1)) ** 2).sum(0))
    k = max(int(round(smooth_s * fs)), 1)
    return np.convolve(e, np.ones(k) / k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--times", default=f"{K.OUT}/s1_times.csv")
    ap.add_argument("--stas", default=",".join(K.DET_SITES))
    ap.add_argument("--out", default=f"{K.OUT}/s2a")
    ap.add_argument("--loud-only", action="store_true")
    ap.add_argument("--max-ev", type=int, default=1200)
    a = ap.parse_args()
    stas = a.stas.split(",")
    C = pd.read_csv(a.times)
    if a.loud_only:
        C = C[C.loud]
    C["t_corr"] = C.t_origin + C.dt
    TT = K.fam_tt(stas)
    fams = sorted(C.bkey.unique())
    n = int(round(LEN * FS))
    dc = K.DayCache(stas)
    acc = {(b, s): [np.zeros((3, n)), np.zeros((3, n)), 0] for b in fams for s in stas}
    for date, g in C.groupby("date"):
        dc.day(date)
        for r in g.itertuples():
            b = int(r.bkey)
            for s in stas:
                c = dc.cut(s, r.t_corr + TT[(b, s)] - PRE, LEN, date)
                if c is None or (c ** 2).sum() <= 0:
                    continue
                c = K.norm3(c)
                A = acc[(b, s)]
                A[0] += c
                A[1] += c ** 2
                A[2] += 1
        K.log(f"[stack] {date} {len(g)}")
    out, arrs = [], {}
    for (b, s), (Sm, Sq, N) in acc.items():
        if N < 20:
            continue
        st = Sm / N
        env = envelope(st)
        # search the peak within +- 8 s of the prediction (index PRE*FS)
        i0 = int(PRE * FS)
        w = int(8 * FS)
        seg = env[i0 - w:i0 + w]
        k = int(np.argmax(seg)) + i0 - w
        bg = np.concatenate([env[:i0 - w], env[i0 + w:]])
        med, mad = float(np.median(bg)), float(1.4826 * np.median(np.abs(bg - np.median(bg))))
        coh = float(N * (st ** 2).mean() / max((Sq / N).mean(), 1e-300))
        arrs[f"{b}_{s}"] = st.astype(np.float32)
        out.append(dict(bkey=b, sta=s, n=N, dist_km=K.epi_km(b, s),
                        off_s=(k - i0) / FS, peak_mad=(env[k] - med) / mad if mad > 0 else np.nan,
                        coh=coh, tt_pred=TT[(b, s)]))
    O = pd.DataFrame(out)
    O.to_csv(f"{a.out}_offsets.csv", index=False)
    np.savez_compressed(f"{a.out}_stacks.npz", **arrs)
    K.log(f"{len(O)} family-station stacks")
    K.log("\n" + str(O.groupby("sta").agg(n=("off_s", "size"), off_med=("off_s", "median"),
                                          off_absmed=("off_s", lambda x: np.median(np.abs(x))),
                                          peak=("peak_mad", "median"), coh=("coh", "median"))))
    K.log(f"\nALL: |off| median {np.median(np.abs(O.off_s)):.2f}s  p90 "
          f"{np.quantile(np.abs(O.off_s), 0.9):.2f}s  coh median {O.coh.median():.2f}")


if __name__ == "__main__":
    main()
