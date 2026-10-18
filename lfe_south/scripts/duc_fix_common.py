#!/usr/bin/env python
"""DUC_FIX (2026-09-30) -- shared definitions for the Ducellier & Creager (2022) extension re-run
with the FIXED detector.  New files only; nothing of t6/t7/t11/sh_v2 is edited.

  his templates  catalogs/ducellier_catalog/src/templates_v1/<fam>/<STA>_<CHA>.pkl  (20 Hz, 1.5-9 Hz)
  his stations   data/Ducellier/stations_permanent.txt ; his families data/Ducellier/families_permanent.txt
  his catalogue rule (his figure / analysis scripts, e.g. data/Ducellier/cumulative_LFEs.py:42):
                 keep a detection when cc * nchannel >= threshold_perm[family] (threshold_cc.txt);
                 a family with threshold_perm == 0 is one of the SIX he eliminated.
  cache          cache_duc/<NET.STA>/<yyyy>/<doy>.npz  (x (C,N) f32, valid (C,N) bool or validp packed,
                 chans, t0, fs) -- the t6 chain; days 2008-04-01..05-10 already on disk are reused.
"""
import os, sys, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import ROOT, WS, RES, TR, FS, log, day_path

DUC = f"{ROOT}/catalogs/ducellier_catalog"
FIX = f"{RES}/duc_fix"
CACHE_DUC = f"{WS}/cache_duc"
TR2 = f"{WS}/traces_duc"          # restricted-channel downloads of this run (traces/ is shared; not written)
RESP = f"{WS}/resp"
T_START, T_END = "2008-01-01", "2026-10-01"      # scan [T_START, T_END)
PROV = {'NC': 'NCEDC', 'BK': 'NCEDC', 'PB': 'IRIS'}


def his_rows():
    rows = []
    for ln in open(f"{DUC}/data/Ducellier/stations_permanent.txt"):
        p = ln.split()
        if len(p) < 10: continue
        loc = "" if p[3] == "--" else (p[3].zfill(2) if p[3].isdigit() else p[3])   # his '0' = '00'
        rows.append(dict(sta=p[0], net=p[1], chans=p[2], loc=loc, prov=p[4],
                         lat=float(p[5]), lon=float(p[6]), on=p[7], off=p[8]))
    return pd.DataFrame(rows)


def families():
    fp = pd.read_csv(f"{DUC}/data/Ducellier/families_permanent.txt", sep=r"\s+", header=None,
                     names=["family", "stations", "duration"])
    fp["family"] = fp.family.astype(str)
    return fp


def thresholds():
    t = pd.read_csv(f"{DUC}/data/Ducellier/threshold_cc.txt", sep=r"\s+", header=None,
                    names=["family", "thr_fame", "thr_perm"])
    t["family"] = t.family.astype(str)
    return dict(zip(t.family, t.thr_perm))


def his_templates():
    """{fam: {(sta, cha): np.float64 array}} for every station his families_permanent.txt names."""
    fp = families(); T = {}
    for r in fp.itertuples():
        stas = set(r.stations.split(","))
        d = f"{DUC}/src/templates_v1/{r.family}"
        if not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".pkl"): continue
            sta, cha = fn[:-4].split("_", 1)
            if sta not in stas: continue
            tr, ori = pickle.load(open(os.path.join(d, fn), "rb"))
            T.setdefault(r.family, {})[(sta, cha)] = np.asarray(tr.data, dtype=np.float64)
    return T


def his_catalogue():
    T = pd.read_parquet(f"{RES}/ducellier_detections.parquet")
    T = T[T.cat == "net"].copy(); T["fam"] = T.fam.astype(str)
    T["t"] = T.t.astype("datetime64[ns]")
    return T


def cache_path(net, sta, d):
    d = pd.Timestamp(d)
    return f"{CACHE_DUC}/{net}.{sta}/{d.year}/{d.dayofyear:03d}.npz"


def load_cache(net, sta, d):
    """(t0, x (C,N) float64, valid (C,N) bool, chans list) or None."""
    p = cache_path(net, sta, d)
    if not os.path.exists(p): return None
    try:
        z = np.load(p)
        x = z["x"].astype(np.float64)
        if "valid" in z.files: v = z["valid"].astype(bool)
        else: v = np.unpackbits(z["validp"], axis=1)[:, :x.shape[1]].astype(bool)
        return float(z["t0"]), x, v, [str(c) for c in z["chans"]]
    except Exception:
        return None


def epochs():
    E = pd.read_csv(f"{FIX}/station_epochs.csv", dtype={"loc": str}).fillna({"loc": ""})
    E["start"] = pd.to_datetime(E.start); E["end"] = pd.to_datetime(E.end)
    return E


def epoch_lookup(E=None):
    """function (sta, day) -> epoch row (namedtuple) or None."""
    E = epochs() if E is None else E
    by = {s: list(g.sort_values("start").itertuples()) for s, g in E.groupby("sta")}
    def f(sta, d):
        d = pd.Timestamp(d)
        for r in by.get(sta, []):
            if r.start <= d < r.end: return r
        return None
    return f
