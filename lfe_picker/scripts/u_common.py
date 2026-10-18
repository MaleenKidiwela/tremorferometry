#!/usr/bin/env python
"""P4-0 shared helpers - Puget/Olympic site-noise + registration-proxy inventory.

Everything this phase writes carries the `u_` prefix and lands in results/p4/, figures/u_*.png,
logs/, traces/, cache20hz/.  Nothing outside lfe_picker is ever written.

Geometry conventions are taken verbatim from results/p3d_puget_geometry/t_geometry2.py:
  * distances are HYPOCENTRAL (sqrt(epicentral^2 + interface depth^2));
  * stations within SITE_MERGE_KM = 1 km are ONE site (single-linkage);
  * cells are the Slab2 28-45 km interface points already tabulated in cell_geometry_by_year.csv,
    tremor_active = >= 100 PNSN epicentres within 20 km.
"""
from __future__ import annotations
import os
import json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P4 = os.path.join(HERE, "results", "p4")
GEO = os.path.join(HERE, "results", "p3d_puget_geometry")
TRACES = os.path.join(HERE, "traces")
CACHE = os.path.join(HERE, "cache20hz")
FIGS = os.path.join(HERE, "figures")
LOGS = os.path.join(HERE, "logs")
EXT_TRACES = os.path.abspath(os.path.join(HERE, "..", "lfe_extension", "traces"))
TREMOR = os.path.abspath(os.path.join(HERE, "..", "catalogs", "pnsn_tremor_cascadia_full.csv"))

R_E = 111.19
SITE_MERGE_KM = 1.0
LAT0, LAT1 = 47.0, 48.3                 # Puget / Olympic band (PLAN P4)
DMAX_KM = 60.0                          # hypocentral radius around an active cell
PROBE_Y0, PROBE_Y1 = 2023, 2025
N_TREMOR_DAYS, N_QUIET_DAYS = 4, 4
QUIET_RADIUS_KM, QUIET_PAD_D = 75.0, 3  # quiet = no epicentre within 75 km on day +- 3 d
TREMOR_RADIUS_KM, TREMOR_MIN_N = 50.0, 50
REF_STATIONS = [("CN", "PGC"), ("PB", "B011")]
PGC_KEY = "CN.PGC"
BAND_PREF = ["HH", "BH", "EH", "CH", "SH", "DP"]     # velocity seismometer bands only
FS_CACHE = 20.0
PSD_BAND = (2.0, 8.0)
E_WIN_S = 6.0                           # the incoherent S-window length (registration proxy)
MAX_CONCURRENT_FDSN = 8

for _d in (P4, FIGS, LOGS, TRACES, CACHE):
    os.makedirs(_d, exist_ok=True)


# ---------------------------------------------------------------- geometry (t_geometry2 verbatim)
def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(lambda v: np.asarray(v, float), (la1, lo1, la2, lo2))
    return np.hypot((lo1 - lo2) * R_E * np.cos(np.deg2rad((la1 + la2) / 2.0)), (la1 - la2) * R_E)


def hyp_km(clat, clon, cdep, slat, slon):
    return np.hypot(epi_km(clat, clon, slat, slon), np.asarray(cdep, float))


def azimuth(clat, clon, slat, slon):
    dx = (np.asarray(slon, float) - clon) * R_E * np.cos(np.deg2rad(clat))
    return np.degrees(np.arctan2(dx, (np.asarray(slat, float) - clat) * R_E)) % 360.0


def az_gap(az):
    if len(az) < 2:
        return 360.0
    a = np.sort(np.asarray(az))
    return float(max(np.diff(a).max(), 360.0 - (a[-1] - a[0])))


def proj(la, lo):
    la = np.asarray(la, float); lo = np.asarray(lo, float)
    return np.column_stack([lo * R_E * np.cos(np.deg2rad(la)), la * R_E])


def active_cells():
    """The tremor-active Slab2 interface cells inside the P4 band (47.0-48.3 N)."""
    G = pd.read_csv(os.path.join(GEO, "cell_geometry_by_year.csv"))
    C = G[G.network == "all"].drop_duplicates("cell_id")
    C = C[(C.lat >= LAT0) & (C.lat <= LAT1) & (C.tremor_active)]
    return C[["cell_id", "lat", "lon", "depth_km", "tremor_n20", "tremor_w"]].reset_index(drop=True)


def cluster_sites(D, lat="lat", lon="lon", key="key"):
    """Single-linkage 1 km site merge; returns a Series key -> site name."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
    D = D.reset_index(drop=True)
    if len(D) == 1:
        return pd.Series({D[key].iloc[0]: D[key].iloc[0].split(".")[-1]})
    xy = proj(D[lat].values, D[lon].values)
    cl = fcluster(linkage(pdist(xy), method="single"), t=SITE_MERGE_KM, criterion="distance")
    out = {}
    for c in np.unique(cl):
        g = D[cl == c]
        name = "+".join(sorted(g[key].str.split(".").str[-1]))
        for k in g[key]:
            out[k] = name
    return pd.Series(out)


def site_map():
    """key (NET.STA) -> site name, from the frozen P3D sites.csv (1 km single-linkage)."""
    S = pd.read_csv(os.path.join(GEO, "sites.csv"))
    out = {}
    for r in S.itertuples():
        for m in str(r.members).split(";"):
            out[m] = r.site
    return out


# ---------------------------------------------------------------- statistics
def band_median_psd(x, fs, f0, f1, nperseg=256, resp_amp=None, drop_zero=True):
    """Median-over-segments Welch PSD, then the median over the [f0, f1] band.

    x           (n_comp, N) float array (counts, band-passed by the frozen chain)
    resp_amp    (n_comp, n_freq) |R(f)| in counts/(m/s) evaluated on the rfft grid; if given the
                PSD is divided by |R|^2 so the answer is a VELOCITY PSD in (m/s)^2/Hz.
    Segments containing a gap (exact-zero sample, the fill value of the frozen merge) are dropped.
    Returns (band_median per component, psd (n_comp, n_freq), freqs, n_seg_used).
    """
    x = np.asarray(x, float)
    nc, N = x.shape
    hop = nperseg // 2
    nseg = max((N - nperseg) // hop + 1, 0)
    if nseg < 4:
        return np.full(nc, np.nan), None, None, 0
    idx = np.arange(nseg) * hop
    w = np.hanning(nperseg)
    U = (w ** 2).sum() * fs
    f = np.fft.rfftfreq(nperseg, 1.0 / fs)
    good = np.ones(nseg, bool)
    if drop_zero:
        for i, s in enumerate(idx):
            if np.any(np.abs(x[:, s:s + nperseg]).sum(0) == 0):
                good[i] = False
    if good.sum() < 4:
        return np.full(nc, np.nan), None, f, int(good.sum())
    use = idx[good]
    P = np.zeros((nc, f.size))
    for c in range(nc):
        seg = np.stack([x[c, s:s + nperseg] for s in use])
        seg = (seg - seg.mean(1, keepdims=True)) * w
        P[c] = np.median(np.abs(np.fft.rfft(seg, axis=1)) ** 2, axis=0) * 2.0 / U
    if resp_amp is not None:
        P = P / np.maximum(np.asarray(resp_amp, float) ** 2, 1e-30)
    m = (f >= f0) & (f <= f1)
    return np.median(P[:, m], axis=1), P, f, int(good.sum())


def med_win_power(x, n, drop_zero=True):
    """Median over the day of the MEAN SQUARE in non-overlapping n-sample windows, per component.

    This is the incoherent (alignment-free) S-window energy statistic used as the registration
    proxy: with n = 6 s * 20 Hz it is the same window the P3D S-excess clause uses, but with no
    stacking and no alignment, so it measures how much a station's 2-8 Hz power rises on a
    tremor day relative to a quiet day.  Windows containing a gap-fill zero are dropped.
    """
    x = np.asarray(x, float)
    nc, N = x.shape
    nw = N // n
    if nw < 10:
        return np.full(nc, np.nan), 0
    X = x[:, :nw * n].reshape(nc, nw, n)
    e = (X ** 2).mean(2)
    if drop_zero:
        good = ~np.any(np.abs(X).sum(0) == 0, axis=1)
    else:
        good = np.ones(nw, bool)
    if good.sum() < 10:
        return np.full(nc, np.nan), int(good.sum())
    return np.median(e[:, good], axis=1), int(good.sum())


def butter_band_amp(f, fs_native, f1=1.0, f2=8.0, order=4):
    """|H(f)| of the frozen chain's zero-phase Butterworth (sosfiltfilt -> |H|^2 applied twice)."""
    from scipy.signal import butter, sosfreqz
    sos = butter(order, [f1 / (fs_native / 2), f2 / (fs_native / 2)], btype="band", output="sos")
    _, h = sosfreqz(sos, worN=2 * np.pi * np.asarray(f, float) / fs_native)
    return np.abs(h) ** 2            # filtfilt = forward+reverse -> amplitude is |H|^2


def cache_path(netsta, y, doy):
    return os.path.join(CACHE, netsta, str(y), f"{doy:03d}.npz")


def trace_path(netsta, y, doy):
    return os.path.join(TRACES, netsta, str(y), f"{doy:03d}.mseed")


def jdump(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=str)
