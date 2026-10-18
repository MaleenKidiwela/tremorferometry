#!/usr/bin/env python
"""Shared definitions for the WA PRE-PRODUCTION steps (2026-10-02, results/wa_pre/).

Read-only on everything that exists; every output goes under results/wa_pre/.  The two download queues'
partially written trace roots (lfe_extension_wa/traces_dl, traces_dl2) are NEVER read here, and neither is
lfe_array_discovery/traces_sealed (that project's sealed holdout).  `wa_common.TRACE_ROOTS` is cut accordingly
by `patch_wa()` before any frozen module (wa_cache / wa_mf / wa_hf) is used.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wa_common as W

ROOT = W.ROOT
STUDY = W.STUDY
PRE = f"{STUDY}/results/wa_pre"
LOGF = f"{STUDY}/logs/wa_pre.log"
T0, T1 = "2007-01-01", "2026-09-15"          # production span (WA_PROD_PLAN W6)
SEALED_MONTHS = {"2007-01", "2009-05"}         # one-look sealed (2007-01 headline, 2009-05 semi-independent)
EXCLUDED_MONTHS = {"2011-04"}                  # no array stations in the truth -> not fit truth

# raw roots on disk; production-visible = wa_common.TRACE_ROOTS minus the download roots
RAW_ROOTS_ALL = ["lfe_array_discovery/traces", "lfe_extension/traces", "lfe_catalog_repair/traces",
                 "lfe_extension_wa/traces", "lfe_picker/traces", "lfe_south/traces", "data/waveforms"]
NEVER = ("traces_dl", "traces_sealed")
PROD_ROOTS = [r for r in W.TRACE_ROOTS if not any(n in r for n in NEVER)]
WA_CACHE = W.CACHE
SVI_CACHES = [f"{ROOT}/lfe_extension/cache20hz", f"{ROOT}/lfe_extension/cache20hz_r1"]


def log(*a):
    os.makedirs(os.path.dirname(LOGF), exist_ok=True)
    W.log(LOGF, *a)


def core_list():
    inv = pd.read_csv(f"{STUDY}/results/wa_prod_plan/inventory_by_station.csv")
    c = inv[inv.raw_c14_26 >= 3900].netsta.tolist()
    assert len(c) == 21, c
    return c


def nodes():
    N = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]]
                  + [pd.read_csv(f"{STUDY}/results/step4b/setup/extra_nodes.csv")]).drop_duplicates("node")
    N["ns"] = N.net + "." + N.sta
    return N


def expected_days(ns, N=None):
    """station-days inside any metadata epoch (node table) of the station, within T0..T1."""
    N = nodes() if N is None else N
    g = N[N.ns == ns]
    s = set()
    for r in g.itertuples():
        a = max(pd.Timestamp(r.start), pd.Timestamp(T0)); b = min(pd.Timestamp(r.end), pd.Timestamp(T1) + pd.Timedelta(days=1))
        if b > a:
            s |= {d.strftime("%Y-%m-%d") for d in pd.date_range(a, b - pd.Timedelta(days=1), freq="D")}
    return s


def raw_files(ns, date, roots=RAW_ROOTS_ALL):
    ts = pd.Timestamp(date)
    out = []
    for r in roots:
        base = r if r.startswith("/") else f"{ROOT}/{r}"
        p = f"{base}/{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
        if os.path.exists(p) and os.path.getsize(os.path.realpath(p)) > 0 and not any(n in os.path.realpath(p) for n in NEVER):
            out.append(p)
    return out


def prod_raw_path(ns, date):
    f = raw_files(ns, date, PROD_ROOTS)
    return f[0] if f else None


def patch_wa(cache_dirs=None):
    """cut the download roots out of wa_common.TRACE_ROOTS; optionally make cache_path look in extra
    cache directories FIRST (new pre-production caches), then the WA cache."""
    W.TRACE_ROOTS[:] = PROD_ROOTS
    if cache_dirs:
        orig = W.CACHE

        def cache_path(netsta, date, _dirs=tuple(cache_dirs) + (orig,)):
            ts = pd.Timestamp(date)
            tail = f"{netsta}/{ts.year}/{int(ts.dayofyear):03d}.npz"
            for c in _dirs:
                if os.path.exists(f"{c}/{tail}"):
                    return f"{c}/{tail}"
            return f"{orig}/{tail}"
        W.cache_path = cache_path
    return W


def marker(name, text=""):
    os.makedirs(PRE, exist_ok=True)
    open(f"{PRE}/{name}", "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + text + "\n")


def tremor_days():
    out = {}
    for rg in ["NW_WA", "SW_WA", "OLY2010"]:
        t = pd.read_csv(f"{W.EXT}/{rg}_tremor_days.csv")
        out[rg] = dict(zip(pd.to_datetime(t.day).dt.strftime("%Y-%m-%d"), t.n_windows.astype(int)))
    return out


REG_TREM = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}
PNSN_START = "2010-01-06"
