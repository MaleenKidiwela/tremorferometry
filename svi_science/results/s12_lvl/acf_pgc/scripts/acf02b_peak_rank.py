#!/usr/bin/env python
"""S12c P1 honesty check (post-hoc, labelled): how unusual is the predicted-band peak? Slide a band of the same width
over 6-30 s lags and record (envelope-peak ratio to the +-4 s local median, years with the same argmax). Fraction of
sliding bands at least as strong (ratio >= observed and years >= observed) = chance of such a peak anywhere."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy.signal import hilbert
from acf_common import *
B = pd.read_csv(f"{WS}/tables/predicted_lag_bands.csv")
rows = []
for ns in ["CN.PGC", "CN.VGZ", "CN.SNB"]:
    for band in ["b05_2"]:
        D, H, R = load(ns, range(2018, 2022), band, pairs=["ZZ", "RR", "TT", "ZR"])
        ok = hour_mask(ns, D, R, "tremorfree")
        yrs = D.astype("datetime64[Y]").astype(int) + 1970
        for pr, ph in [("ZZ", "PP"), ("RR", "SS"), ("TT", "SS"), ("ZR-", "PS"), ("ZR+", "PS")]:
            S, n = day_stack(H[pr.rstrip("+-")], ok); S = split_pair(S, pr)
            g = np.isfinite(S[:, 0]); L = np.nanmean(S[g], 0)
            Ys = [np.nanmean(S[g & (yrs == y)], 0) for y in range(2018, 2022) if (g & (yrs == y)).sum() > 30]
            t = LAG_AUTO; env = np.abs(hilbert(L)); envy = [np.abs(hilbert(y)) for y in Ys]
            def score(lo, hi):
                mm = (t >= lo) & (t <= hi); k = np.argmax(env[mm]); tp = t[mm][k]
                nb = (np.abs(t - tp) <= 4) & (t >= 2)
                return env[mm][k] / np.median(env[nb]), sum(abs(t[mm][np.argmax(e[mm])] - tp) <= 0.3 for e in envy)
            for iface in ["top", "bottom"]:
                r = B[(B.station == ns) & (B.interface == iface)].iloc[0]
                lo, hi = r[f"{ph}_lo"] - 0.3, r[f"{ph}_hi"] + 0.3
                ro, yo = score(lo, hi)
                w = hi - lo; sl = [score(s, s + w) for s in np.arange(6.0, 30.0 - w, 0.1)]
                frac = np.mean([(a >= ro) and (b >= yo) for a, b in sl])
                rows.append(dict(station=ns, band=band, pair=pr, interface=iface, ratio=round(ro, 2), years=yo,
                                 frac_sliding_as_strong=round(frac, 3)))
T = pd.DataFrame(rows); T.to_csv(f"{WS}/tables/peak_rank_pilot.csv", index=False); print(T.to_string(index=False))
