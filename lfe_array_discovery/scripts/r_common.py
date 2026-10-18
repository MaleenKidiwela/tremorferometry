#!/usr/bin/env python
"""Shared helpers for agent R (confirmation-rule calibration + the two Washington
confirmations).  Prefix `r_`.

Nothing here changes any recipe constant.  It only (i) loads/saves templates in exactly the
layout `k_grow` writes, (ii) computes epicentral / hypocentral distances, and (iii) reads and
writes the CALIBRATED CONFIRMATION RULE as a small json so that the WA stages cannot silently
use a different rule from the one the north measured.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")

OUT = f"{RES}/r_confirm"
RULE = f"{OUT}/calibrated_rule.json"


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


# ---------------------------------------------------------------- templates -------------------
def load_tpl(path, grps=None):
    """{g: {sta: (tpl (3,n) float32, off_s, peak_mad)}} -- k_grow's own layout."""
    z = np.load(path, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        g = int(g)
        if grps is not None and g not in grps:
            continue
        TPL.setdefault(g, {})[s] = (np.asarray(z[k], np.float32), float(v[0]), float(v[1]))
    return TPL


def save_tpl(path, TPL):
    np.savez_compressed(
        path,
        **{f"{g}_{s}": TPL[g][s][0].astype(np.float32) for g in TPL for s in TPL[g]},
        meta=np.array(json.dumps({f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                  for g in TPL for s in TPL[g]})))


# ---------------------------------------------------------------- geometry --------------------
def epi_km(la1, lo1, la2, lo2):
    import s_common as S
    return S.epi_km(la1, lo1, la2, lo2)


def hyp_km(la1, lo1, dp1, la2, lo2, dp2):
    d = epi_km(la1, lo1, la2, lo2)
    return np.sqrt(np.asarray(d) ** 2 + (np.asarray(dp1) - np.asarray(dp2)) ** 2)


# ---------------------------------------------------------------- the rule --------------------
def write_rule(d, path=RULE):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(d, open(path, "w"), indent=1, default=float)
    return d


def read_rule(path=RULE):
    if not os.path.exists(path):
        raise SystemExit(f"missing calibrated rule {path} -- run scripts/r_north.py --stage calib")
    return json.load(open(path))


def pct_inclusive(x, frac=0.85, side="ge"):
    """The threshold that KEEPS at least `frac` of the finite values in x.

    side='ge' -- a clause of the form  value >= T  (e.g. detection days): T is the largest
                 integer/observed value such that >= frac of x satisfies it.
    side='le' -- a clause of the form  value <= T  (e.g. distance, false alarms): T is the
                 frac-quantile of x (so exactly frac of x satisfies it).
    """
    v = np.asarray([float(t) for t in x if np.isfinite(t)])
    if not len(v):
        return np.nan
    if side == "le":
        return float(np.quantile(v, frac))
    cand = np.unique(v)
    best = float(cand.min())
    for c in cand:
        if float((v >= c).mean()) >= frac:
            best = float(c)
    return best
