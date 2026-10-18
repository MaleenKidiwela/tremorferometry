#!/usr/bin/env python
"""Shared definitions for the WA permanent-station TRIAL extension (lfe_extension_wa/).

Frozen measurement chain copied from lfe_extension (SVI, results/frozen_operating_point.json):
  20 Hz, 1-8 Hz Butterworth-4 zero-phase, component order Z / N|1 / E|2 resolved NOT rotated,
  12-s templates starting 0.5 s before the S (envelope) peak, per-channel L2-normalised,
  station cc = mean of its channels' normalised cc, weight = n_channels x SNR weight,
  site merge B (within-site w-weighted mean, site weight sqrt(sum w^2)), per family-day network MAD,
  record >= 8.0 MAD, T_ext 8.25, T_core 8.75, >= 2 sites per family-day, one-to-one +-2 s scoring,
  site-rolled scrambled control (one shift in [300, 3000] s per site per day).
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
STUDY = f"{ROOT}/lfe_extension_wa"
RES = f"{STUDY}/results"
SETUP = os.environ.get("WA_SETUP", f"{RES}/setup")   # step 4b: WA_SETUP=results/step4b/setup
CACHE = f"{STUDY}/cache20hz"
TRACES = f"{STUDY}/traces"
EXT = f"{ROOT}/data/extension_stations"
TRACES_DL = f"{STUDY}/traces_dl"   # production downloads, tiers A+B (scripts/wa_dl.py, 2026-10-02); only days missing
                                    # from the roots above, so appending it changes no existing raw_path() answer
TRACES_DL2 = f"{STUDY}/traces_dl2"   # tiers C, D, D2 + temporary networks (scripts/wa_dl2.py); appended once WA_DL2_DONE
                                     # existed (WA production queue, 2026-10-03); only days missing from the roots above
TRACE_ROOTS = [f"{ROOT}/lfe_array_discovery/traces", f"{ROOT}/lfe_extension/traces",
               f"{ROOT}/lfe_catalog_repair/traces", TRACES, TRACES_DL, TRACES_DL2]

FS = 20.0
NDAY = int(86400 * FS)          # 1,728,000 samples, midnight-aligned grid
TPL_S = 12.0                    # frozen template length
PRE = 0.5                       # template starts PRE s before the S (envelope) peak
NT = int(TPL_S * FS)            # 240
LONG_PRE, LONG_LEN = 20.0, 52.0  # stack window [t-20, t+32]
CLEAR_PRE, CLEAR_LEN = 12.0, 44.0  # the catalogue's stack window for the clear-arrival rule [t-12, t+32]
S_SEARCH = (0.0, 15.0)          # S peak searched inside the Royer detection window [t, t+15]
MIN_EVENTS = 30                 # e5r_templates.py: families with < 30 events get no template
T_REC, T_EXT, T_CORE = 8.0, 8.25, 8.75
MATCH_S = 2.0
BLOCK_S, GUARD_S = 7200.0, 60.0  # build/test halves (e4_template_test.py --half, block 7200 s, guard 60 s)
SITE_KM = 3.0                   # stations < 3 km apart are one site (catalogue n_sites rule)
ROLL_S = (300.0, 3000.0)
HF_BAND = (10.0, 19.0)
TEMP_NETS = {"XU", "XG", "ZK", "Z3", "7D", "ZU", "YW", "XD", "ZZ", "TA", "XA", "YB", "3F", "X4", "PO", "C8"}

REGIONS = ["Royer_NW_WA", "Royer_SW_WA", "Olympic_2010", "Olympic_2010_armB"]

# windows: label -> list of (start, end, kind)
WINDOWS = {
    "V2008": [("2008-05-01", "2008-05-26", "episode"), ("2008-06-16", "2008-06-18", "post")],
    "V2010": [("2010-07-14", "2010-07-18", "quiet"), ("2010-08-07", "2010-09-04", "episode"),
              ("2010-09-14", "2010-09-18", "quiet")],
    "T2012": [("2012-08-20", "2012-08-24", "quiet"), ("2012-09-05", "2012-10-13", "episode"),
              ("2012-10-18", "2012-10-22", "quiet")],
    # step 4b (2026-09-29): A2011 = array-era templating only (NW array detections, 2011-08 ETS; never scored);
    # T2013 / T2014 = the strong WA ETS of tremor_episodes.csv with the 5 nearest tremor-free days either
    # side (no window in the NW_WA or OLY2010 PNSN tremor-day tables)
    "A2011": [("2011-07-30", "2011-08-27", "episode")],
    "T2013": [("2013-08-30", "2013-09-03", "quiet"), ("2013-09-04", "2013-10-13", "episode"),
              ("2013-10-14", "2013-10-18", "quiet")],
    "T2014": [("2014-11-10", "2014-11-13", "quiet"), ("2014-11-15", "2014-11-15", "quiet"),
              ("2014-11-18", "2014-12-14", "episode"), ("2014-12-15", "2014-12-19", "quiet")],
}


def window_days(w):
    out = []
    for a, b, k in WINDOWS[w]:
        for d in pd.date_range(a, b, freq="D"):
            out.append((d.strftime("%Y-%m-%d"), k))
    return out


def all_days(windows=("V2008", "V2010", "T2012")):
    """the trial's windows by default (step 4b passes its own list)"""
    return [(w, d, k) for w in windows for d, k in window_days(w)]


def km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1) - lo2) * 111.195 * np.cos(np.radians((np.asarray(la1) + la2) / 2)),
                    (np.asarray(la1) - la2) * 111.195)


def cache_path(netsta, date):
    ts = pd.Timestamp(date)
    return f"{CACHE}/{netsta}/{ts.year}/{int(ts.dayofyear):03d}.npz"


def raw_candidates(netsta, date):
    ts = pd.Timestamp(date)
    return [f"{r}/{netsta}/{ts.year}/{int(ts.dayofyear):03d}.mseed" for r in TRACE_ROOTS]


def raw_path(netsta, date):
    for p in raw_candidates(netsta, date):
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def half_of(t):
    """0 = build, 1 = test, -1 = within GUARD_S of a block edge."""
    t = np.asarray(t, float)
    b = np.floor(t / BLOCK_S)
    r = t - b * BLOCK_S
    h = (b % 2).astype(int)
    h[(r < GUARD_S) | (r > BLOCK_S - GUARD_S)] = -1
    return h


def load_families():
    return pd.read_csv(f"{SETUP}/families.csv")


def load_stations():
    return pd.read_csv(f"{SETUP}/stations.csv")


def envelope_hilbert_sum(x):
    from scipy.signal import hilbert
    return np.abs(hilbert(x, axis=-1)).sum(0)


def env_rule(x, fs=FS):
    """catalogue clear-arrival rule (data/catalogue/scripts/11_build_families_v1_1.py env_rule)."""
    from scipy.signal import hilbert
    x = np.asarray(x, float)
    keep = [c for c in x if np.all(np.isfinite(c)) and np.any(c != 0)]
    if not keep:
        return np.nan, np.nan
    e = np.abs(hilbert(np.array(keep), axis=1)).sum(0)
    med = np.median(e)
    if not med > 0:
        return np.nan, np.nan
    return float(e.max() / med), float((e > e.max() / 2).sum() / fs)


def bp(x, fs, f1, f2, order=4):
    """Butterworth band-pass designed WITH fs = the actual sampling rate (the 20-vs-40 sps slip)."""
    from scipy.signal import butter, sosfiltfilt
    sos = butter(order, [f1, f2], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def match_one_to_one(a, b, tol=MATCH_S):
    """greedy one-to-one by |dt| (smallest first) between sorted arrays a, b -> (ia, ib) matched indices."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0, int), np.zeros(0, int)
    ob = np.argsort(b); bs = b[ob]
    pairs = []
    lo = np.searchsorted(bs, a - tol, "left"); hi = np.searchsorted(bs, a + tol, "right")
    for i in range(len(a)):
        for j in range(lo[i], hi[i]):
            pairs.append((abs(a[i] - bs[j]), i, j))
    pairs.sort()
    ua, ub, ia, ib = set(), set(), [], []
    for _, i, j in pairs:
        if i in ua or j in ub:
            continue
        ua.add(i); ub.add(j); ia.append(i); ib.append(ob[j])
    return np.array(ia, int), np.array(ib, int)


def log(path, *a):
    import time
    s = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + " ".join(str(x) for x in a)
    print(s, flush=True)
    if path:
        with open(path, "a") as f:
            f.write(s + "\n")
