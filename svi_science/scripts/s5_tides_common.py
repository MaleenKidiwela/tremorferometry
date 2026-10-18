"""Shared helpers for the S5 tidal-stress setup."""
import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
import pyTMD.predict as P
import pyTMD.constituents as PC

ROOT = Path('/home/jovyan/tremorferometry')
OUT = ROOT / 'svi_science/results/s5_setup'
TIDES = ROOT / 'data/tides'
MAJOR = ['q1', 'o1', 'p1', 'k1', 'n2', 'm2', 's2', 'k2', '2n2']
MODEL_CONST = MAJOR + ['m4']           # constituents taken from the ocean model
RHO_W, G0 = 1025.0, 9.81
CENTER = (48.6, -123.6)                 # SVI family cloud centre (lat, lon)
_EPOCH = np.datetime64('1992-01-01T00:00:00')


def days1992(times):
    return (np.asarray(times, dtype='datetime64[ns]') - _EPOCH) / np.timedelta64(1, 'D')


def predict(hc, times, corrections='FES', minor=True):
    """hc: dict constituent -> complex array (n,) [hc = A exp(-i G)]; returns (n, nt) real."""
    ds = xr.Dataset({k: (['station'], np.atleast_1d(np.asarray(v, dtype=complex))) for k, v in hc.items()})
    t = days1992(times)
    out = P.time_series(t, ds, corrections=corrections).values
    if minor:
        mi = P.infer_minor(t, ds, corrections=corrections)
        out = out + np.asarray(mi)
    return out


def design(times, consts, corrections='FES'):
    """LSQ design matrix: eta = sum_c f (Re hc cos th - Im hc sin th)."""
    MJD = days1992(times) + 48622.0
    pu, pf, Gd = PC.arguments(MJD, consts, corrections=corrections)
    th = np.radians(Gd) + pu
    return np.hstack([pf * np.cos(th), -pf * np.sin(th)])


def harmonic_fit(times, eta, consts, corrections='FES'):
    X = design(times, consts, corrections)
    X = np.c_[X, np.ones(len(times))]
    ok = np.isfinite(eta)
    c, *_ = np.linalg.lstsq(X[ok], eta[ok], rcond=None)
    n = len(consts)
    hc = c[:n] + 1j * c[n:2 * n]
    res = eta[ok] - X[ok] @ c
    return dict(zip(consts, hc)), np.std(res)


def gc_dist_az(lat0, lon0, lat, lon):
    """great-circle distance (km) and azimuth (deg from N) from (lat0,lon0) to points."""
    R = 6371.0
    p0, l0, p, l = map(np.radians, (lat0, lon0, lat, lon))
    dl = l - l0
    a = np.sin((p - p0) / 2) ** 2 + np.cos(p0) * np.cos(p) * np.sin(dl / 2) ** 2
    d = 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    az = np.degrees(np.arctan2(np.sin(dl) * np.cos(p), np.cos(p0) * np.sin(p) - np.sin(p0) * np.cos(p) * np.cos(dl)))
    return d, az % 360
