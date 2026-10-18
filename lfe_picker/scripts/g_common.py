#!/usr/bin/env python
"""Gates — shared machinery (waveform store, travel times, stacking, band ratios).

Detector-agnostic: nothing here knows about templates, Bostock families or the matched filter.
A PROPOSAL is (a set of absolute detection times on the ORIGIN axis, a proposed source cell
lat/lon/dep).  Everything the gates need is derived from that plus a station/region config.

Measurement chain is the FROZEN one: cache20hz = 1-8 Hz Butterworth-4 zero-phase at the native
rate, integer-decimated to 20 Hz, components Z / N|1 / E|2 (scripts/o_cache_month.py).  The HF
ratio cannot be measured in that cache (Nyquist 10 Hz) and is therefore read from the NATIVE-rate
day traces, exactly as lfe_extension/scripts/a1_hf_veto.py does.
"""
from __future__ import annotations

import glob
import os
import sys
import threading
from collections import OrderedDict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PICKER = os.path.dirname(HERE)
ROOT = os.path.dirname(PICKER)
REP = f"{ROOT}/lfe_catalog_repair"
EXT = f"{ROOT}/lfe_extension"
RESULTS = f"{PICKER}/results"
GATES = f"{RESULTS}/gates"
FIGS = f"{PICKER}/figures"

FS = 20.0
HALF_W = 20.0                      # +-20 s cut around the predicted S
NW = int(round(2 * HALF_W * FS))   # 800 samples
IDX_S = int(round(HALF_W * FS))    # 400 = the predicted S sample

# read-only caches, searched in this order
CACHE_DIRS = [f"{REP}/cache20hz", f"{PICKER}/cache20hz", f"{EXT}/cache20hz"]
TRACE_DIRS = [f"{REP}/traces", f"{PICKER}/traces", f"{EXT}/traces"]

R_KM_DEG = 111.195
LAT_REF = 48.6

sys.path.insert(0, HERE)


# ----------------------------------------------------------------------------- geometry / stations
def epi_km(lat1, lon1, lat2, lon2, lat_ref=LAT_REF):
    cs = np.cos(np.radians(lat_ref))
    return np.hypot((np.asarray(lon1, float) - lon2) * R_KM_DEG * cs,
                    (np.asarray(lat1, float) - lat2) * R_KM_DEG)


def stations(path=f"{REP}/results/e5r_stations.csv"):
    return pd.read_csv(path).set_index("sta")


# The two co-located clusters collapse to ONE SITE (Merlin F10(i); a5_locate.SITE_MEMBERS).
SITE_MEMBERS = {
    "PGC": ["PGC", "B009", "B010", "B011"],
    "LRIV": ["LRIV", "B005", "B006", "B007"],
    "VGZ": ["VGZ"], "SNB": ["SNB"], "NLLB": ["NLLB"],
    "B926": ["B926"], "B001": ["B001"], "B003": ["B003"], "B004": ["B004"],
}
STA2SITE = {m: q for q, mem in SITE_MEMBERS.items() for m in mem}


def station_corrections(path=f"{EXT}/results/a5_station_corrections.csv"):
    """a5 intraslab-earthquake station terms; corr_S is ADDED to the Crosson prediction to give the
    expected OBSERVED arrival (a5_relocate subtracts it when it converts a pick into a travel time)."""
    C = pd.read_csv(path).set_index("sta")
    out = {}
    for s in C.index:
        cs, cp = float(C.loc[s, "corr_S"]), float(C.loc[s, "corr_P"])
        out[s] = (0.0 if not np.isfinite(cs) else cs, 0.0 if not np.isfinite(cp) else cp)
    return out


# ----------------------------------------------------------------------------- travel times
_TT = {}


def tt_table(cache=f"{GATES}/g_tt_table.npz", xmax=260.0, dx=1.0, zmin=4.0, zmax=70.0, dz=0.5):
    """Crosson 1-D (d3_vel1d) travel-time table, (xs, zs, TP, TS)."""
    key = "tt"
    if key in _TT:
        return _TT[key]
    if cache and os.path.exists(cache):
        d = np.load(cache)
        _TT[key] = (d["xs"], d["zs"], d["TP"], d["TS"])
        return _TT[key]
    import o_vel1d as V   # identical copy of lfe_extension d3_vel1d (Crosson 1-D)
    xs = np.arange(0.0, xmax + dx, dx)
    zs = np.arange(zmin, zmax + dz, dz)
    TP = np.empty((zs.size, xs.size))
    TS = np.empty((zs.size, xs.size))
    for i, z in enumerate(zs):
        for j, x in enumerate(xs):
            TP[i, j] = V.tt(x, z, "P")
            TS[i, j] = V.tt(x, z, "S")
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez(cache, xs=xs, zs=zs, TP=TP, TS=TS)
    _TT[key] = (xs, zs, TP, TS)
    return _TT[key]


def _interp2(xs, zs, T, x, z):
    xi = np.clip((np.asarray(x, float) - xs[0]) / (xs[1] - xs[0]), 0, xs.size - 1.001)
    zi = np.clip((np.asarray(z, float) - zs[0]) / (zs[1] - zs[0]), 0, zs.size - 1.001)
    x0 = np.floor(xi).astype(int); z0 = np.floor(zi).astype(int)
    fx = xi - x0; fz = zi - z0
    return ((1 - fz) * ((1 - fx) * T[z0, x0] + fx * T[z0, x0 + 1])
            + fz * ((1 - fx) * T[z0 + 1, x0] + fx * T[z0 + 1, x0 + 1]))


def ttSP(dist_km, dep_km):
    xs, zs, TP, TS = tt_table()
    return (float(_interp2(xs, zs, TS, dist_km, dep_km)),
            float(_interp2(xs, zs, TP, dist_km, dep_km)))


_SLAB = {}


def slab2_depth(lats, lons, path=f"{ROOT}/data/cas_slab2_input_04-18.csv"):
    from scipy.interpolate import griddata
    if "pts" not in _SLAB:
        sl = pd.read_csv(path, usecols=["lat", "lon", "depth", "etype"])
        sl = sl[~sl.etype.isin(["BA", "TO"])]
        sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna()
        _SLAB["pts"] = sl[["lon", "lat"]].values
        _SLAB["z"] = sl.depth.values
    pts = np.column_stack([np.atleast_1d(lons), np.atleast_1d(lats)])
    z = griddata(_SLAB["pts"], _SLAB["z"], pts, method="linear")
    zn = griddata(_SLAB["pts"], _SLAB["z"], pts, method="nearest")
    return np.where(np.isfinite(z), z, zn)


# ----------------------------------------------------------------------------- waveform store
class DayStore:
    """RAM-cached 20 Hz day arrays.  Thread-safe, LRU-bounded (default ~24 GB)."""

    def __init__(self, cache_dirs=None, max_days=1200):
        self.dirs = list(cache_dirs or CACHE_DIRS)
        self.max_days = max_days
        self._d = OrderedDict()
        self._lock = threading.Lock()
        self._index = {}

    def index(self, netsta, year):
        k = (netsta, year)
        if k not in self._index:
            got = {}
            for c in self.dirs:
                for f in sorted(glob.glob(f"{c}/{netsta}/{year}/*.npz")):
                    d = int(os.path.basename(f)[:3])
                    got.setdefault(d, f)
            self._index[k] = got
        return self._index[k]

    def available_days(self, netsta, year):
        return sorted(self.index(netsta, year))

    def get(self, netsta, year, doy):
        k = (netsta, year, doy)
        with self._lock:
            if k in self._d:
                self._d.move_to_end(k)
                return self._d[k]
        p = self.index(netsta, year).get(doy)
        if p is None:
            v = None
        else:
            try:
                z = np.load(p)
                v = (float(z["t0"]), z["x"])
            except Exception:
                v = None
        with self._lock:
            self._d[k] = v
            self._d.move_to_end(k)
            while len(self._d) > self.max_days:
                self._d.popitem(last=False)
        return v

    def span(self, netsta, year, doys):
        """(t_lo, t_hi) covering the cached data of those days (for random-time draws)."""
        lo, hi = np.inf, -np.inf
        for d in doys:
            v = self.get(netsta, year, d)
            if v is None:
                continue
            lo = min(lo, v[0])
            hi = max(hi, v[0] + v[1].shape[1] / FS)
        return (lo, hi) if np.isfinite(lo) else (np.nan, np.nan)

    def cut(self, netsta, starts, n=NW, normalise=True):
        """Cut n-sample 3-C windows starting at absolute times `starts`.  Returns (W, keep)."""
        starts = np.asarray(starts, float)
        if starts.size == 0:
            return None, np.zeros(0, int)
        ts = pd.to_datetime(starts, unit="s")
        yrs = ts.year.to_numpy()
        dys = ts.dayofyear.to_numpy()
        segs, keep = [], []
        order = np.argsort(yrs * 1000 + dys, kind="stable")
        for j in order:
            v = self.get(netsta, int(yrs[j]), int(dys[j]))
            if v is None:
                continue
            t0, X = v
            i0 = int(round((starts[j] - t0) * FS))
            if i0 < 0 or i0 + n > X.shape[1]:
                continue
            seg = X[:, i0:i0 + n].astype(np.float64)
            if normalise:
                r = np.sqrt((seg ** 2).mean())
                if not np.isfinite(r) or r <= 0:
                    continue
                seg = seg / r
            elif not np.isfinite(seg).all():
                continue
            segs.append(seg)
            keep.append(j)
        if not segs:
            return None, np.zeros(0, int)
        return np.stack(segs), np.asarray(keep, int)

    def stack(self, netsta, s_times, n=NW, pre=HALF_W):
        """Linear stack of rms-normalised windows aligned so that index IDX_S is the predicted S."""
        W, keep = self.cut(netsta, np.asarray(s_times, float) - pre, n=n)
        if W is None:
            return None, 0
        return W.mean(0), len(keep)


# ----------------------------------------------------------------------------- statistics
def sdb(stack, idx_s=IDX_S, sig=(-1.0, 5.0), noi=(-12.0, -6.0)):
    """10 log10 E(S+sig) / E(S+noi) on a 3-C stack.  Merlin's Gate-1 statistic."""
    a1 = idx_s + int(round(sig[0] * FS)); b1 = idx_s + int(round(sig[1] * FS))
    a0 = idx_s + int(round(noi[0] * FS)); b0 = idx_s + int(round(noi[1] * FS))
    s = stack[:, a1:b1]; nz = stack[:, a0:b0]
    return float(10 * np.log10(max((s ** 2).mean(), 1e-300) / max((nz ** 2).mean(), 1e-300)))


def envelope(st, smooth_s=0.5):
    from scipy.signal import hilbert
    e = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    k = max(2, int(round(smooth_s * FS)))
    return np.convolve(e, np.ones(k) / k, mode="same")


def peak_offset_s(stack, idx_s=IDX_S, search_s=6.0):
    env = envelope(stack)
    w = int(round(search_s * FS))
    a, b = max(idx_s - w, 0), min(idx_s + w, len(env))
    return (a + int(np.argmax(env[a:b])) - idx_s) / FS


def cc_template(stack, tpl, idx_s=IDX_S, pre_tpl=1.0, lag_s=1.0):
    """max 3-C cc of a stack against an L2-normalised (3, n) template that starts at S - pre_tpl."""
    nL = tpl.shape[1]
    i0 = idx_s - int(round(pre_tpl * FS))
    lag = int(round(lag_s * FS))
    best = -9.0
    for k in range(-lag, lag + 1):
        a = i0 + k
        if a < 0 or a + nL > stack.shape[1]:
            continue
        seg = stack[:, a:a + nL]
        nn = np.linalg.norm(seg, axis=1)
        if np.any(nn <= 0):
            continue
        best = max(best, float((seg / nn[:, None] * tpl).sum() / tpl.shape[0]))
    return best


# ----------------------------------------------------------------------------- native-rate HF ratio
class NativeStore:
    """Native-rate day traces for the HF-ratio class test (the 20 Hz cache stops at 10 Hz)."""

    def __init__(self, trace_dirs=None, max_days=60):
        self.dirs = list(trace_dirs or TRACE_DIRS)
        self.max_days = max_days
        self._d = OrderedDict()

    def get(self, netsta, year, doy):
        k = (netsta, year, doy)
        if k in self._d:
            self._d.move_to_end(k)
            return self._d[k]
        v = None
        for root in self.dirs:
            p = f"{root}/{netsta}/{year}/{doy:03d}.mseed"
            if os.path.exists(p) and os.path.getsize(p) > 0:
                try:
                    v = self._read(p, netsta, year, doy)
                except Exception:
                    v = None
                if v is not None:
                    break
        self._d[k] = v
        self._d.move_to_end(k)
        while len(self._d) > self.max_days:
            self._d.popitem(last=False)
        return v

    @staticmethod
    def _read(p, netsta, year, doy):
        from obspy import read
        import o_chanorder as chanorder
        st = read(p)
        for tr in st:
            r = float(tr.stats.sampling_rate)
            if r > 0 and abs(r - round(r)) / r < 1e-4 and r != round(r):
                tr.stats.sampling_rate = float(round(r))
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
        keep = [c for c in order if c is not None]
        if not keep:
            return None
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        fs = float(st[0].stats.sampling_rate)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        return X, fs, t0


def _bp(x, fs, f1, f2):
    from scipy.signal import butter, sosfiltfilt
    ny = fs / 2.0
    f2 = min(f2, 0.95 * ny)
    sos = butter(4, [f1 / ny, f2 / ny], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def hf_ratio_members(nstore, netsta, year_doy_times, win_s=12.0, hf=(10.0, 19.0), lf=(2.0, 8.0)):
    """Per-member E(HF)/E(LF) over the win_s window starting at the predicted S - 1 s, on the NATIVE
    trace.  Returns the array of per-member ratios (the statistic Track A calibrated: LFE core median
    ~0.25, veto >= 3.5).  The COHERENT stack cannot carry this band (alignment jitter >> 1/(2*19 Hz)),
    so the class statistic is the median over members, not a ratio measured on the stack."""
    out = []
    for (year, doy), times in year_doy_times.items():
        v = nstore.get(netsta, year, doy)
        if v is None:
            continue
        X, fs, t0 = v
        n = int(round(win_s * fs))
        # filter the whole day once per band
        try:
            Xh = _bp(X, fs, *hf)
            Xl = _bp(X, fs, *lf)
        except Exception:
            continue
        for t in times:
            i0 = int(round((t - 1.0 - t0) * fs))
            if i0 < 0 or i0 + n > X.shape[1]:
                continue
            eh = float((Xh[:, i0:i0 + n] ** 2).sum())
            el = float((Xl[:, i0:i0 + n] ** 2).sum())
            if el > 0:
                out.append(eh / el)
    return np.asarray(out, float)


def align_stack(S, iters=2, lag=30):
    """cc-aligned linear stack — VERBATIM from lfe_extension/scripts/e5r_templates.align_stack
    (rms-normalised windows (N,3,n) -> stack (3,n), aligned windows, per-event shifts in samples).
    Used ONLY for the Gate-2 pick stacks, so that the picks the a5 estimator receives are made on the
    same kind of stack a5 itself picked on.  The Gate-1 statistic never uses it: an alignment fitted
    on the data would import exactly the coherence it is supposed to test."""
    N, nc, n = S.shape
    tot = np.zeros(N, int)
    st = S.mean(0)
    ar = np.arange(n)
    for _ in range(iters):
        ref = st - st.mean(1, keepdims=True)
        ccs = np.empty((2 * lag + 1, N))
        for kk, k in enumerate(range(-lag, lag + 1)):
            ccs[kk] = np.einsum("ijk,jk->i", S, np.roll(ref, -k, axis=1))
        best = ccs.argmax(0) - lag
        idx = (ar[None, :] - best[:, None]) % n
        S = S[np.arange(N)[:, None, None], np.arange(nc)[None, :, None], idx[:, None, :]]
        tot += best
        st = S.mean(0)
    return st, S, tot
