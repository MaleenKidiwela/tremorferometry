#!/usr/bin/env python
"""Step 1c' -- shared helpers for the z_ stages.  Prefix `z_`.

Nothing here changes an engine constant: these are cut/stack utilities that read the frozen 20 Hz
cache through `w_common.DayLoader` on the ABSOLUTE time axis, plus the kept/hidden family
bookkeeping the D8 diagnostic, the oracle ceiling and the scorer all need.
"""
from __future__ import annotations

import collections
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

FS = W.FS
OUT = f"{W.RES}/step1c_prime"
os.makedirs(OUT, exist_ok=True)


def bostock_year(year=None, kept_only=True, window_days=None):
    """Bostock 2005 rows (bkey, t_origin, date, hh), hidden families dropped when kept_only."""
    import p_common as P
    year = W.STEP1_YEAR if year is None else year
    hid, _, _ = S.holdouts()
    B = P.bostock()
    B = B[B.year == year].copy()
    if kept_only:
        B = B[~B.bkey.isin(hid)]
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B["hh"] = B.time.dt.hour.astype(int)
    if window_days is not None:
        B = B[B.date.isin(set(window_days))]
    return B.reset_index(drop=True)


def busiest_family_hours(B, n, min_events=10):
    """[(bkey, date, hour, n_events)] -- the busiest (family, hour) cells."""
    g = (B.groupby(["bkey", "date", "hh"]).size().rename("n").reset_index()
         .sort_values("n", ascending=False))
    g = g[g.n >= min_events]
    return [(int(r.bkey), str(r.date), int(r.hh), int(r.n)) for r in g.head(n).itertuples()]


class Cutter:
    """Absolute-time cuts out of the frozen 20 Hz cache (one LRU day loader)."""

    def __init__(self, maxn=40):
        self.loader = W.DayLoader(maxn=maxn)

    def day(self, sta, t_abs):
        ts = pd.Timestamp(float(t_abs), unit="s")
        return self.loader.get(sta, int(ts.year), int(ts.dayofyear))

    def cut(self, sta, t_abs, len_s):
        """(3, n) float64 at the absolute start time, or None."""
        v = self.day(sta, t_abs)
        if v is None:
            return None
        t0, X = v
        n = int(round(len_s * FS))
        i0 = int(round((t_abs - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        seg = X[:, i0:i0 + n].astype(np.float64)
        if not np.isfinite(seg).all():
            return None
        return seg


def norm_rows(Sg):
    """Per-trace rms normalisation of a (N, 3, n) stack input; degenerate traces dropped."""
    e = (Sg ** 2).sum((1, 2))
    ok = e > 0
    Sg = Sg[ok]
    e = np.sqrt(e[ok] / (Sg.shape[1] * Sg.shape[2]))
    return Sg / e[:, None, None]


def stack_coherence(Sg):
    """N x <stack^2> / <trace^2>: 1.0 = incoherent, > 1 = coherent (Merlin's Q1/(3) statistic)."""
    if len(Sg) < 2:
        return np.nan
    Sg = norm_rows(np.asarray(Sg, float))
    if len(Sg) < 2:
        return np.nan
    st = Sg.mean(0)
    return float(len(Sg) * (st ** 2).mean() / max((Sg ** 2).mean(), 1e-300))


def envelope_1d(x, smooth_s=0.5, fs=FS):
    from scipy.signal import hilbert
    env = np.sqrt((np.abs(hilbert(x, axis=-1)) ** 2).sum(0))
    k = max(int(round(smooth_s * fs)), 1)
    return np.convolve(env, np.ones(k) / k, mode="same")


def mad_height(env, i, guard_s=5.0, fs=FS):
    """(peak - median) / (1.4826 * MAD) of the envelope, background outside +- guard_s of i."""
    m = np.ones(len(env), bool)
    a, b = max(i - int(guard_s * fs), 0), min(i + int(guard_s * fs), len(env))
    m[a:b] = False
    bg = env[m]
    if bg.size < 20:
        return np.nan
    med = float(np.median(bg))
    mad = float(1.4826 * np.median(np.abs(bg - med)))
    return float((env[i] - med) / mad) if mad > 0 else np.nan


def fam_tt_map(det, ST=None):
    return S.fam_tt(det, ST=(W.station_table() if ST is None else ST))


def site_median(per_sta, mem, order):
    """{site: median over its member stations} of a {sta: value} dict."""
    out = {}
    for q in order:
        v = [per_sta[s] for s in mem[q] if s in per_sta and np.isfinite(per_sta[s])]
        if v:
            out[q] = float(np.median(v))
    return out


def raw_path(sta, date):
    import glob
    ts = pd.Timestamp(date)
    hit = glob.glob(f"{W.STUDY}/traces/bed2005_polaris/*.{sta}/{ts.year}/{ts.dayofyear:03d}.mseed")
    return hit[0] if hit else None
