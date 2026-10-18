#!/usr/bin/env python
"""Rational-chain equivalence check (P4-0 section 5 recipe).

For each test station-day, two 20 Hz products from the SAME raw trace:
  A  rational twin  (u_cache.cache_one_rational): bandpass 1-8 Hz at the NATIVE rate, then
     resample_poly(up, down) for the exact rational fs/20  (50 Hz -> 2/5);
  B  frozen integer chain (o_cache_month.cache_one) applied to a 100 Hz version of the trace:
     resample_poly to 100 Hz first, then bandpass 1-8 Hz at 100 Hz, then decimate by 5.
Reported per component: Pearson cc, max |A-B| normalised by the RMS of A, and the RMS ratio.
Interior samples only (5 s guard at each end) so filtfilt edge transients are not counted.
A station whose cc < 0.999 is flagged "rational, non-equivalent" for the freeze.
-> results/rational_equivalence.csv
"""
from __future__ import annotations
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "lfe_picker", "scripts")))
import d_common as D                                                   # noqa: E402
import o_cache_month as OC                                             # noqa: E402
import o_chanorder as chanorder                                        # noqa: E402

FS = 20.0
GUARD_S = 5.0
CC_BAR = 0.999


def load_raw(netsta, day):
    from obspy import read
    d = pd.Timestamp(day)
    p = os.path.join(D.TRACES, netsta, str(d.year), f"{int(d.dayofyear):03d}.mseed")
    st = read(p)
    for tr in st:
        fs_ = float(tr.stats.sampling_rate)
        fr = round(fs_)
        if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4:
            tr.stats.sampling_rate = float(fr)
    st.merge(method=1, fill_value=0)
    chans = sorted({tr.stats.channel for tr in st})
    order = [c for c in chanorder.resolve_order(chans, station=netsta, strict=True) if c]
    fs0 = float(st.select(channel=order[0])[0].stats.sampling_rate)
    npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in order)
    X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in order])
    return X, fs0, order


def chain_rational(X, fs0):
    from fractions import Fraction
    Y = OC.bp(X, fs0)
    fr = Fraction(fs0 / FS).limit_denominator(1000)
    up, down = fr.denominator, fr.numerator
    return resample_poly(Y, up, down, axis=-1) if (up, down) != (1, 1) else Y


def chain_frozen_via100(X, fs0):
    """Resample the trace to 100 Hz, then run the FROZEN integer chain (1-8 Hz at 100, /5)."""
    from fractions import Fraction
    if fs0 != 100.0:
        fr = Fraction(100.0 / fs0).limit_denominator(1000)
        X = resample_poly(X, fr.numerator, fr.denominator, axis=-1)
    Y = OC.bp(X, 100.0)
    return resample_poly(Y, 1, 5, axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", default="2010-08-18,2010-08-20")
    ap.add_argument("--xg", default="XG.BH01,XG.CL13,XG.TB01")
    ap.add_argument("--pb", default="PB.B001,PB.B003,PB.B005")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    rows = []
    for netsta in a.xg.split(",") + a.pb.split(","):
        for day in a.days.split(","):
            try:
                X, fs0, order = load_raw(netsta, day)
            except Exception as e:
                rows.append(dict(netsta=netsta, date=day, status="error",
                                 err=f"{type(e).__name__}: {e}"))
                continue
            A = chain_rational(X, fs0)
            B = chain_frozen_via100(X, fs0)
            n = min(A.shape[1], B.shape[1])
            g = int(GUARD_S * FS)
            A, B = A[:, g:n - g], B[:, g:n - g]
            for i, c in enumerate(order):
                aa, bb = A[i], B[i]
                rms = float(np.sqrt(np.mean(aa ** 2)))
                cc = float(np.corrcoef(aa, bb)[0, 1]) if rms > 0 else np.nan
                rows.append(dict(netsta=netsta, date=day, chan=c, fs_native=fs0,
                                 chain=D.chain_for(fs0), n=int(aa.size), cc=cc,
                                 max_abs_diff=float(np.max(np.abs(aa - bb))),
                                 max_diff_over_rms=float(np.max(np.abs(aa - bb)) / max(rms, 1e-30)),
                                 rms_ratio=float(np.sqrt(np.mean(bb ** 2)) / max(rms, 1e-30)),
                                 status="ok", err=""))
            D.log(f"[equiv] {netsta} {day} fs={fs0:g} cc="
                  f"{','.join(f'{r["cc"]:.6f}' for r in rows[-len(order):])}")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(D.RES, "rational_equivalence.csv"), index=False)
    ok = R[R.status == "ok"]
    V = (ok.groupby(["netsta", "chain"])
           .agg(min_cc=("cc", "min"), med_cc=("cc", "median"),
                max_diff_over_rms=("max_diff_over_rms", "max"),
                med_rms_ratio=("rms_ratio", "median")).reset_index())
    V["verdict"] = np.where(V.min_cc >= CC_BAR, "equivalent", "rational, non-equivalent")
    print(V.to_string(index=False))
    V.to_csv(os.path.join(D.RES, "rational_equivalence_verdict.csv"), index=False)


if __name__ == "__main__":
    main()
