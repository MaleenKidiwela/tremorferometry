#!/usr/bin/env python
"""Phase picks on the long per-station stacks.  Prefix `u_`.

Four observables per (family, station), every one of them a COLUMN, never a gate:
  t_env    envelope maximum over the 3 components  -- EXACTLY the current locator's observable
           (k_grow cuts the template 2 s before this sample, k_locate adds the 2 s back)
  t_onS    S ONSET: Maeda AIC minimum on the horizontal amplitude in a window that ends at the
           horizontal envelope maximum; this is what Chestler picks and what the envelope maximum
           is not -- the maximum of a 1-8 Hz LFE stack sits 0.3-1.5 s into the S wavetrain and
           moves with the coda, the onset does not
  t_onP    P ONSET on the vertical, searched strictly BEFORE the S onset (1.0-9.0 s earlier),
           AIC-refined about the strongest vertical envelope excursion in that window
  sp_xc    S-P from the ENVELOPE cross-correlation of the vertical against the horizontal at the
           same station, lags 1.0-9.0 s (Chestler & Creager 2017's "S minus P times by cross
           correlating S- versus P- waves at the same stations")
plus the quality columns used as WEIGHTS (never as cuts): the envelope peak height in MAD, the
onset sharpness, the P height in MAD and the S-P cross-correlation coefficient.

All times are in seconds relative to the family's NODE time (the long stack starts
U.LONG_PRE = 12 s before it), so they are directly comparable with the k_locate convention
`off + TPL_PRE`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                        # noqa: E402

SMOOTH_S = 0.4
AIC_PRE = 4.0          # AIC search starts this far before the envelope maximum
AIC_POST = 0.5
P_MIN_LEAD = 1.0       # a P must precede the S onset by at least this
P_MAX_LEAD = 9.0       # ... and at most this (S-P is 2.5-7 s over this aperture and depth range)


def envelope(x, fs=U.FS, smooth_s=SMOOTH_S):
    from scipy.signal import hilbert
    e = np.sqrt((np.abs(hilbert(x, axis=-1)) ** 2).sum(0))
    k = max(int(round(smooth_s * fs)), 1)
    return np.convolve(e, np.ones(k) / k, mode="same")


def mad_height(e, k):
    med = float(np.median(e))
    mad = float(1.4826 * np.median(np.abs(e - med)))
    return (float(e[k]) - med) / mad if mad > 0 else 0.0


def aic_onset(x, k0, k1):
    """Maeda AIC minimum of the 1-D series x on [k0, k1); returns the index or -1."""
    k0 = max(int(k0), 2)
    k1 = min(int(k1), x.size - 2)
    if k1 - k0 < 6:
        return -1
    seg = np.asarray(x[k0:k1], float)
    n = seg.size
    c1 = np.cumsum(seg)
    c2 = np.cumsum(seg ** 2)
    k = np.arange(1, n - 1)
    v1 = np.maximum(c2[k - 1] / k - (c1[k - 1] / k) ** 2, 1e-30)
    m = n - k
    v2 = np.maximum((c2[-1] - c2[k - 1]) / m - ((c1[-1] - c1[k - 1]) / m) ** 2, 1e-30)
    aic = k * np.log(v1) + m * np.log(v2)
    return int(k0 + k[int(np.argmin(aic))])


def xc_lag(a, b, lmin, lmax, fs=U.FS):
    """Lag L (s) in [lmin, lmax] maximising corr(a[t], b[t-L]); returns (lag, cc)."""
    a = np.asarray(a, float) - a.mean()
    b = np.asarray(b, float) - b.mean()
    n = a.size
    i0, i1 = int(round(lmin * fs)), int(round(lmax * fs))
    best, bl = -2.0, np.nan
    na = np.sqrt((a ** 2).sum())
    for L in range(i0, i1 + 1):
        if L <= 0 or L >= n:
            continue
        aa, bb = a[L:], b[:n - L]
        d = na * np.sqrt((bb ** 2).sum())
        if d <= 0:
            continue
        c = float((aa * bb).sum() / (np.sqrt((aa ** 2).sum()) * np.sqrt((bb ** 2).sum())))
        if c > best:
            best, bl = c, L / fs
    return bl, best


def pick_one(A, fs=U.FS):
    """A = (3, N) ZNE long stack.  Returns a dict of observables in seconds from the node time."""
    e3 = envelope(A, fs)
    eH = envelope(A[1:], fs)
    eZ = envelope(A[:1], fs)
    k3 = int(np.argmax(e3))
    kH = int(np.argmax(eH))
    h3 = mad_height(e3, k3)
    hH = mad_height(eH, kH)
    hamp = np.sqrt((A[1:] ** 2).sum(0))
    kon = aic_onset(hamp, kH - int(AIC_PRE * fs), kH + int(AIC_POST * fs))
    if kon < 0:
        kon = kH
    # onset sharpness: post/pre amplitude ratio across the onset (0.5 s either side)
    w = int(0.5 * fs)
    pre = hamp[max(kon - w, 0):kon]
    post = hamp[kon:kon + w]
    sharp = float(post.mean() / max(pre.mean(), 1e-30)) if pre.size and post.size else 1.0
    # P: strongest vertical envelope excursion strictly before the S onset
    a0 = kon - int(P_MAX_LEAD * fs)
    a1 = kon - int(P_MIN_LEAD * fs)
    tP, hP = np.nan, 0.0
    if a1 - a0 > int(0.5 * fs) and a0 >= 0:
        kp = a0 + int(np.argmax(eZ[a0:a1]))
        hP = mad_height(eZ, kp)
        kpo = aic_onset(np.abs(A[0]), kp - int(2.0 * fs), kp + int(0.5 * fs))
        if kpo < 0:
            kpo = kp
        tP = (kpo - U.I_NODE) / fs
    # S-P by envelope cross-correlation (vertical vs horizontal) over the S wavetrain
    b0 = max(kon - int(P_MAX_LEAD * fs) - int(2 * fs), 0)
    b1 = min(kH + int(6 * fs), eH.size)
    sp, spcc = np.nan, -2.0
    if b1 - b0 > int(8 * fs):
        sp, spcc = xc_lag(eH[b0:b1], eZ[b0:b1], P_MIN_LEAD, P_MAX_LEAD, fs)
    return dict(t_env=(k3 - U.I_NODE) / fs, h_env=h3,
                t_envH=(kH - U.I_NODE) / fs, h_envH=hH,
                t_onS=(kon - U.I_NODE) / fs, sharp=sharp,
                t_onP=tP, h_P=hP, sp_xc=sp, sp_cc=spcc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stacks", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    z = np.load(a.stacks, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    keys = [k for k in z.files if k != "meta"]
    U.log(f"{len(keys)} station stacks")
    rows = []
    for i, k in enumerate(keys):
        g, s = k.split("_", 1)
        d = pick_one(np.asarray(z[k], float))
        d.update(grp=int(g), sta=s, n_stack=int(meta.get(k, 0)))
        rows.append(d)
        if (i + 1) % 2000 == 0:
            U.log(f"  {i+1}/{len(keys)}")
    P = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    P.to_parquet(a.out, index=False)
    U.log(f"wrote {a.out}  {len(P)} picks, {P.grp.nunique()} families; "
          f"P found {int(np.isfinite(P.t_onP).sum())}, S-P xc cc>0.5 {int((P.sp_cc > 0.5).sum())}")


if __name__ == "__main__":
    main()
