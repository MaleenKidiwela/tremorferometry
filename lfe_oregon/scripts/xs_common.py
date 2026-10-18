"""Shared helpers for the lfe_oregon chain (Thomas & Bostock 2015 reproduction)."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(ROOT)
CACHE = f"{ROOT}/cache40"
RES = f"{ROOT}/results"
FIGS = f"{ROOT}/figures"
FS = 40.0
NDT = 1.0 / FS
STATIONS = pd.read_csv(f"{PROJ}/notes/thomas2015_candidate_stations.csv")
STATIONS["netsta"] = STATIONS.net + "." + STATIONS.sta
COORD = {r.netsta: (r.lat, r.lon) for r in STATIONS.itertuples()}
DETECT_DAYS = ["2009-09-03", "2009-09-04", "2009-09-05"]          # T&B 2015 sec.2
ALL_DAYS = [str(d.date()) for d in pd.date_range("2009-08-03", "2009-09-17")]
XS_STATIONS = ["XA.ALP0", "XA.ALS0", "XA.BWO0"]                  # T&B 2015 sec.2


def utm_km(netsta):
    from pyproj import Proj
    la, lo = COORD[netsta]
    x, y = Proj(proj="utm", zone=10)(lo, la)
    return x / 1e3, y / 1e3


def load_day(netsta, date):
    """(t0, X (3,N) float32 ZNE, m (N,) bool) or None."""
    p = f"{CACHE}/{netsta}/{date}.npz"
    if not os.path.exists(p):
        return None
    z = np.load(p)
    return float(z["t0"]), np.ascontiguousarray(z["x"], np.float32), z["m"].astype(bool)


def patch_ab_common():
    """Point lfe_array_discovery's ab_common at THIS workspace's 1-10 Hz cache (in-process only;
    the lfe_array_discovery files are not modified)."""
    import sys
    sys.path.insert(0, f"{PROJ}/lfe_array_discovery/scripts")
    import ab_common as A
    A.cache_path = lambda sta, date: f"{CACHE}/{sta}/{date}.npz"
    A.null_remap = lambda sta, date, k=None: date
    A.BAND = (1.0, 10.0)
    assert A.FS == FS
    os.environ.pop("AB_MUTE_DIR", None)
    return A
