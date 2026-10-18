#!/usr/bin/env python
"""Reproduce-Bostock loop -- shared helpers.  Prefix `k_`.

Nothing here is frozen; this is a development chain whose only job is to find where our
implementation loses Bostock's known families.  Truth = the 100 KEPT families.  The 30 hidden
families in lfe_picker/results/holdouts.json are never read by any k_ script until the one blind
run.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)   # = lfe_south when run from lfe_south/scripts
RES = f"{STUDY}/results"
OUT = f"{RES}/repro"
PICKER = "/home/jovyan/tremorferometry/lfe_picker"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
sys.path.insert(0, HERE)
sys.path.insert(0, f"{PICKER}/scripts")
os.makedirs(OUT, exist_ok=True)

FS = 20.0
NDAY = 1728000
CACHES = [f"{STUDY}/cache_duc", f"{STUDY}/cache20hz"]
RAW2005 = f"{STUDY}/traces/bed2005_polaris"

# 2005 bed window (results/step1c_bed.json).  The bed is switchable by environment so the SAME
# literal chain can be run on the 2012 / 2010 PERMANENT-station beds without a code change:
#   K_YEAR, K_DAY0, K_DAY1  (K_OUT relocates the output directory)
BED_YEAR = int(os.environ.get("K_YEAR", "2005"))
DAY0 = os.environ.get("K_DAY0", "2005-08-29")
DAY1 = os.environ.get("K_DAY1", "2005-09-21")
OUT = os.environ.get("K_OUT", OUT)
os.makedirs(OUT, exist_ok=True)
QUIET_DAYS = [f"2005-0{m}-{d:02d}" for m, d in
              [(6, 25), (6, 26), (6, 27), (6, 28), (6, 29), (6, 30)] + [(7, d) for d in range(1, 15)]]
# 12 detection sites of the frozen bed; SHVB carries a 4.76 Hz tone (STEP1_RESULTS F1) -> excluded
# from every sum here.  JRBC/LZB/SNB/VGZ were "withheld" for the blind Gates; the reproduce loop
# uses the 11 clean detection sites for discovery statistics and may use the withheld ones only for
# TIME CORRECTION of catalogued events (never for a discovery decision).
DET_SITES = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC",
             "B003", "B005", "B006", "B007"]
CLEAN11 = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC", "B003", "B005"]
ALL2005 = ["GLBC", "JRBC", "LCBC", "MGCB", "SHVB", "TWBB", "LZB", "PGC", "SNB", "VGZ",
           "B001", "B003", "B004", "B005", "B006", "B007", "SILB", "SSIB", "TSJB", "TWKB"]
SHARED_0512 = ["PGC", "VGZ", "SNB", "B001", "B003", "B004", "B005", "B006", "B007"]


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


def bed_days():
    d = pd.date_range(DAY0, DAY1, freq="D")
    return [x.strftime("%Y-%m-%d") for x in d]


# ---------------------------------------------------------------- catalogue / geometry ---------
_C = {}


def families():
    import p_common as P
    if "F" not in _C:
        _C["F"] = P.families()
    return _C["F"]


def station_table():
    import w_common as W
    if "ST" not in _C:
        _C["ST"] = W.station_table()
    return _C["ST"]


def bostock(year=None, kept_only=None):
    """Full Bostock export on the frozen time base (p_common.bostock)."""
    import p_common as P
    if "B" not in _C:
        _C["B"] = P.bostock()
    B = _C["B"]
    if year is not None:
        B = B[B.year == year]
    if kept_only is True:
        B = B[~B.hidden]
    elif kept_only is False:
        B = B[B.hidden]
    return B.reset_index(drop=True)


def bed_catalogue(kept_only=True):
    B = bostock(BED_YEAR, kept_only=kept_only).copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B = B[B.date.isin(set(bed_days()))]
    return B.reset_index(drop=True)


def fam_tt(stas=None):
    """{(bkey, sta): predicted S delay} -- Crosson 1-D + a5 station correction (0 where absent)."""
    import s_common as S
    stas = ALL2005 if stas is None else list(stas)
    k = tuple(sorted(stas))
    if k not in _C:
        _C[k] = S.fam_tt(stas, ST=station_table())
    return _C[k]


def epi_km(b, sta):
    import s_common as S
    F, ST = families(), station_table()
    return float(S.epi_km(F.lat[b], F.lon[b], ST.lat[sta], ST.lon[sta]))


def nearest_stas(b, stas=None, n=4):
    stas = DET_SITES if stas is None else stas
    return sorted(stas, key=lambda s: epi_km(b, s))[:n]


# ---------------------------------------------------------------- 20 Hz cache ------------------
_IDX = {}


def cache_index():
    if not _IDX:
        for r in CACHES:
            for p in glob.glob(f"{r}/*/*/*.npz"):
                parts = p.split(os.sep)
                sta = parts[-3].split(".")[1]
                _IDX.setdefault((sta, int(parts[-2]), int(parts[-1][:3])), p)
    return _IDX


def load_day(sta, date):
    """(t0, (3, N) float32 ZNE) from the frozen 20 Hz 1-8 Hz cache, rotated for boreholes."""
    import s_common as S
    import w_common as W                                                     # noqa: F401
    ts = pd.Timestamp(date)
    idx = cache_index()
    key = (sta, int(ts.year), int(ts.dayofyear))
    if key not in idx:
        return None
    if sta in S.AZ12:
        v = S.load_day_zne(sta, key[1], key[2], idx)
        # some station-years are cached with a single channel (EHZ-only); a (1, N) array reached
        # the autocorrelation and crashed the 2010 bed.  Only 3-C station-days are usable here.
        return v if (v is not None and v[1].shape[0] == 3) else None
    z = np.load(idx[key])
    x = z["x"]
    # SOUTHERN FIDELITY EDIT (agent T, 2026-09-23): accept 1-C and 2-C station-days.  ~65 % of
    # Ducellier's permanent bed is vertical-only (NC EHZ/SHZ); his detector correlates PER CHANNEL
    # and averages, so those stations are full members of the network.  Requiring 3-C here would
    # run a thinned network, which is not the published method.
    if x.shape[0] < 1:
        return None
    return float(z["t0"]), np.ascontiguousarray(x, dtype=np.float32)


class DayCache:
    """{sta: (t0, X)} for ONE date, built once and reused for every event of that date."""

    def __init__(self, stas):
        self.stas = list(stas)
        self.date = None
        self.D = {}

    def day(self, date):
        if date != self.date:
            self.D = {}
            for s in self.stas:
                v = load_day(s, date)
                if v is not None:
                    self.D[s] = v
            self.date = date
        return self.D

    def cut(self, sta, t_abs, len_s, date=None):
        date = (pd.Timestamp(float(t_abs), unit="s").strftime("%Y-%m-%d") if date is None
                else date)
        D = self.day(date)
        v = D.get(sta)
        if v is None:
            return None
        t0, X = v
        n = int(round(len_s * FS))
        i0 = int(round((t_abs - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        seg = X[:, i0:i0 + n].astype(np.float64)
        return seg if np.isfinite(seg).all() else None


# ---------------------------------------------------------------- native-rate cache ------------
def raw_path(sta, date):
    ts = pd.Timestamp(date)
    hit = glob.glob(f"{RAW2005}/*.{sta}/{ts.year}/{ts.dayofyear:03d}.mseed")
    return hit[0] if hit else None


def load_day_native(sta, date, fmin=2.0, fmax=8.0, target_fs=100.0):
    """(t0, (3, N) float64 ZNE, fs) at NATIVE rate from the raw miniSEED, band-passed fmin-fmax.

    Same rotation convention as the 20 Hz cache (s_common.AZ12 for the boreholes)."""
    import obspy
    import s_common as S
    from scipy.signal import butter, sosfiltfilt
    p = raw_path(sta, date)
    if p is None:
        return None
    try:
        st = obspy.read(p)
    except Exception:
        return None
    st.merge(method=1, fill_value=0)
    if len(st) < 3:
        return None
    fs = float(st[0].stats.sampling_rate)
    if any(abs(float(tr.stats.sampling_rate) - fs) > 1e-6 for tr in st):
        return None
    ts = pd.Timestamp(date)
    t_day = ts.timestamp()
    n = int(round(86400 * fs))
    ch = {}
    for tr in st:
        c = tr.stats.channel
        x = np.zeros(n, np.float64)
        i0 = int(round((tr.stats.starttime.timestamp - t_day) * fs))
        a, b = max(i0, 0), min(i0 + tr.stats.npts, n)
        if b > a:
            x[a:b] = tr.data[a - i0:b - i0]
        ch.setdefault(c, np.zeros(n))
        ch[c] = ch[c] + x
    keys = sorted(ch)
    zc = [c for c in keys if c.endswith("Z")]
    if not zc:
        return None
    if sta in S.AZ12 and any(c.endswith("1") for c in keys):
        a1, a2 = S.AZ12[sta]
        c1 = [c for c in keys if c.endswith("1")][0]
        c2 = [c for c in keys if c.endswith("2")][0]
        r1, r2 = np.radians(a1), np.radians(a2)
        # [N; E] = inv([[cos a1, cos a2],[sin a1, sin a2]]) applied to [c1; c2]
        M = np.array([[np.cos(r1), np.cos(r2)], [np.sin(r1), np.sin(r2)]])
        N, E = M @ np.vstack([ch[c1], ch[c2]])
        X = np.vstack([ch[zc[0]], N, E])
    else:
        nc = [c for c in keys if c.endswith("N")]
        ec = [c for c in keys if c.endswith("E")]
        if not nc or not ec:
            return None
        X = np.vstack([ch[zc[0]], ch[nc[0]], ch[ec[0]]])
    X = X - X.mean(1, keepdims=True)
    sos = butter(4, [fmin / (fs / 2), fmax / (fs / 2)], btype="band", output="sos")
    X = sosfiltfilt(sos, X, axis=-1)
    return t_day, np.ascontiguousarray(X), fs


class NativeDayCache:
    def __init__(self, stas, fmin=2.0, fmax=8.0):
        self.stas, self.fmin, self.fmax = list(stas), fmin, fmax
        self.date, self.D = None, {}

    def day(self, date):
        if date != self.date:
            self.D = {}
            for s in self.stas:
                v = load_day_native(s, date, self.fmin, self.fmax)
                if v is not None:
                    self.D[s] = v
            self.date = date
        return self.D

    def cut(self, sta, t_abs, len_s, date=None):
        date = (pd.Timestamp(float(t_abs), unit="s").strftime("%Y-%m-%d") if date is None
                else date)
        v = self.day(date).get(sta)
        if v is None:
            return None
        t0, X, fs = v
        n = int(round(len_s * fs))
        i0 = int(round((t_abs - t0) * fs))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        seg = X[:, i0:i0 + n]
        return seg if np.isfinite(seg).all() else None


# ---------------------------------------------------------------- statistics -------------------
def stack_coherence(Sg):
    """N * <stack^2> / <trace^2>.  1 = incoherent."""
    Sg = np.asarray(Sg, float)
    if len(Sg) < 2:
        return np.nan
    e = (Sg ** 2).mean((1, 2))
    ok = e > 0
    Sg = Sg[ok] / np.sqrt(e[ok])[:, None, None]
    if len(Sg) < 2:
        return np.nan
    return float(len(Sg) * (Sg.mean(0) ** 2).mean() / max((Sg ** 2).mean(), 1e-300))


def norm3(x):
    """per-trace (3-component jointly) rms normalisation of a (3, n) segment."""
    r = np.sqrt((x ** 2).mean())
    return x / r if r > 0 else x
