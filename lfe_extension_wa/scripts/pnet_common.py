#!/usr/bin/env python
"""PNET (lfe_picker P1b PhaseNet-style U-Net, trained on Bostock SVI LFEs) on the WA permanent stations.

Model input contract (lfe_picker/scripts/p_scan.cf_day, p_build_mm.py, p1b_train.py):
  60-s windows at 20 Hz, 3 components in Z, N, E order (TRUE geographic N/E: the SVI store was ROTATED with the
  StationXML azimuths), the frozen 1-8 Hz Butterworth-4 zero-phase chain, each window peak-normalised over its three
  components; overlapping windows with step = 1/3 window, only the central third kept.
The WA cache (wa_cache.py) is the same 20 Hz 1-8 Hz chain but in the RECORDED order Z / N|1 / E|2 (not rotated), so:
  * PB boreholes EH1/EH2 are rotated to N/E with the StationXML azimuths below (lfe_picker/results/p4/stationxml/PB.xml,
    identical in lfe_array_discovery/results/stationxml/PB.xml);
  * US.NLWA BHN/BHE and BH1/BH2 are at 4 / 94 deg (US.xml, 2007-02-27 onwards) -> rotated;
  * UW.HDW records EHZ + ONE horizontal (EH2 at 45 deg): a single horizontal cannot be rotated -> Z only;
  * 1-component stations (UW EHZ-only) enter as (Z, 0, 0); the SVI model never saw this input (a column: n_horiz);
  * samples outside the cache validity mask `m` are set to 0 (cf_day zeroes the CF where all channels are 0);
  * any 3-C station with a dead horizontal -> both horizontals 0 (cannot rotate).
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

PICKER = "/home/jovyan/tremorferometry/lfe_picker/scripts"
PNET_W05 = "/home/jovyan/tremorferometry/lfe_picker/results/p1/models/pnet_w05.pt"
OUT = f"{W.RES}/pnet"
LOG = f"{W.STUDY}/logs/pnet.log"

# EH1/EH2 (or BH1/BH2) azimuths, deg CW from N, dip 0 (StationXML; one epoch per station in 2008-2011)
AZ12 = {
    "PB.B001": (335.5, 65.5), "PB.B003": (239.0, 329.0), "PB.B004": (167.5, 257.5), "PB.B005": (90.5, 180.5),
    "PB.B006": (201.0, 291.0), "PB.B007": (36.5, 126.5), "PB.B009": (0.0, 90.0), "PB.B010": (228.5, 318.5),
    "PB.B011": (11.0, 101.0), "PB.B013": (31.5, 301.5), "PB.B014": (170.0, 260.0), "PB.B017": (206.0, 296.0),
    "PB.B018": (21.5, 111.5), "PB.B019": (305.0, 35.0), "PB.B941": (173.5, 263.5), "PB.B943": (301.5, 31.5),
    "PB.B926": (299.0, 29.0),
    "US.NLWA": (4.0, 94.0),
}
Z_ONLY = {"UW.HDW"}                   # EHZ + a single horizontal (EH2 at 45 deg): horizontals dropped

# Stage B: whole stations (= whole 3-km sites, each of these is its own site) never used for fine-tuning
HELD_OUT = ["UW.HDW", "UW.OSD", "PB.B018", "PB.B017", "US.NLWA", "UW.GNW"]


def rot(netsta):
    a1, a2 = np.radians(AZ12[netsta])
    A = np.array([[np.cos(a1), np.sin(a1)], [np.cos(a2), np.sin(a2)]])
    return np.linalg.inv(A)           # (c1, c2) -> (N, E); exact for any two non-parallel horizontals


def to_zne(netsta, x, m, live):
    """cache (3, N) Z / N|1 / E|2  ->  (3, N) float32 Z, N, E; invalid samples zeroed.  Returns (y, n_horiz)."""
    y = np.array(x, np.float32, copy=True)
    y[:, ~m] = 0.0
    live = np.asarray(live, bool)
    y[~live] = 0.0
    if netsta in Z_ONLY or not (live[1] and live[2]):
        y[1:] = 0.0
        return y, 0
    if netsta in AZ12:
        R = rot(netsta).astype(np.float32)
        y[1:3] = R @ y[1:3]
    return y, 2


def load_zne(netsta, date):
    p = W.cache_path(netsta, date)
    if not os.path.exists(p):
        return None
    z = np.load(p)
    y, nh = to_zne(netsta, z["x"], z["m"], z["live"])
    return dict(x=y, m=z["m"], live=z["live"], node=str(z["node"]), n_horiz=nh)


def cf_path(tag, netsta, date):
    return f"{OUT}/cf_{tag}/{netsta}/{date}.npy"


def load_cf(tag, netsta, date):
    p = cf_path(tag, netsta, date)
    if not os.path.exists(p):
        return None
    return np.load(p).astype(np.float32)


def maxfilt(v, half=10):
    """centred running max over +-half samples (0.5 s at 20 Hz)"""
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(v, size=2 * half + 1, mode="constant", cval=0.0)


def auc(pos, neg):
    """Mann-Whitney AUC (ties = 0.5), NaNs dropped"""
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos = pos[np.isfinite(pos)]; neg = neg[np.isfinite(neg)]
    if len(pos) < 5 or len(neg) < 5:
        return np.nan
    from scipy.stats import rankdata
    r = rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def label_pairs(Q):
    """Stage B label rule on the build-half template QC table (tpl_qc_build.csv):
    a (family, station node) pair is LABELLED when its stack (>= 30 build-half events) has a clear-ish arrival,
    imp >= 3 and dur < 3 s (the catalogue's clear rule is imp > 5 and dur < 3 s); label weight 1.0 if clear,
    0.5 if only clear-ish.  S time = event origin t + s_off (the stack's envelope peak)."""
    q = Q[(Q.n_events >= W.MIN_EVENTS) & Q.s_off.notna()].copy()
    q["lab"] = (q.imp >= 3) & (q.dur < 3)
    q = q[q.lab].copy()
    q["lab_w"] = np.where(q.clear.fillna(False).astype(bool), 1.0, 0.5)
    return q[["fid", "netsta", "node", "s_off", "lab_w", "imp", "dur", "clear", "n_events"]]
