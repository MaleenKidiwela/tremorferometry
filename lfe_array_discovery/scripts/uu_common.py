#!/usr/bin/env python
"""Shared machinery for the LOCATION REBUILD (agent U).  Prefix `u_`.

Nothing in the detection or growth chain is touched.  This module only reads grown templates,
their detections and the day caches, and provides:
  * the CALIBRATION TRUTH -- Bostock's published family coordinates for the 58 blind-rediscovered
    northern families (development) and the 18 located sealed families (blind confirmation);
  * travel-time tables (Crosson 1-D via p_vel1d, with a Vp/Vs knob) on an arbitrary grid;
  * the misfit engine (L2 / robust L1 / weighted), the coarse+fine grid search, the envelope
    back-projection, the station-pair cross-correlation differential times, a double-difference
    pass and a station bootstrap.

BED ADAPTER.  `bed()` returns an object with `station_table()`, `DayCache`, `station_corr()` and
the bed day list.  For the north bed that is `k_common`; for the Washington beds it is
`o_common` (site/beam bed) registered by `o_run.py`.  No k_/n_/o_ file is edited.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
ROOT = "/home/jovyan/tremorferometry"
RES = f"{STUDY}/results"
OUT = os.environ.get("U_OUT", f"{RES}/u_locate")
PICKER = f"{ROOT}/lfe_picker"
sys.path.insert(0, HERE)
sys.path.append(f"{PICKER}/scripts")
os.makedirs(OUT, exist_ok=True)

FS = 20.0
LONG_PRE = 12.0          # the long stack starts this far before the node time (k_grow.LONG_PRE)
LONG_LEN = 44.0          # ... and is this long                                 (k_grow.LONG_LEN)
NLONG = int(round(LONG_LEN * FS))
I_NODE = int(round(LONG_PRE * FS))
MAX_STACK = 200          # k_grow.MAX_STACK
CHESTLER = (f"{ROOT}/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = (np.asarray(x, float) for x in (la1, lo1, la2, lo2))
    return np.sqrt(((la1 - la2) * 111.2) ** 2
                   + ((lo1 - lo2) * 111.2 * np.cos(np.radians(0.5 * (la1 + la2)))) ** 2)


# ================================================================= calibration truth ============
def bostock_families():
    """PUBLISHED family coordinates.  SOURCE, verified: every row of
    `lfe_catalog_repair/results/e5r_families_2012.csv` is one of the 130 distinct (lat, lon, dep)
    triples carried by the event rows of Bostock et al. (2015)'s own catalogue export
    `results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv`; `bkey` is this project's
    index into them.  Convention: geographic lat/lon, depth in km POSITIVE DOWN below sea level
    (Bostock's published convention; the stations of this bed are at 0-1.3 km elevation, which is
    a <=0.4 s S systematic and is carried in the station terms, never removed silently)."""
    F = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv").set_index("bkey")
    return F[["lat", "lon", "dep"]]


def verify_truth_source():
    """Assert the family table really is Bostock's published coordinates.  Returns a dict."""
    F = bostock_families()
    B = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/bostock/"
                    "bostock2015_sloweq_FULL_2003-2013.csv", usecols=["lat", "lon", "dep"])
    trip = set(zip(np.round(B.lat, 4), np.round(B.lon, 4), np.round(B.dep, 4)))
    ours = set(zip(np.round(F.lat, 4), np.round(F.lon, 4), np.round(F.dep, 4)))
    return dict(n_families=len(F), n_catalog_triples=len(trip),
                n_ours_in_catalog=len(ours & trip), dep_min=float(F.dep.min()),
                dep_max=float(F.dep.max()), dep_median=float(F.dep.median()))


def north_truth():
    """(dev, blind) DataFrames: grp -> bkey -> published lat/lon/dep.

    dev   = templates matching one of the 58 blind-rediscovered KEPT families (development).
    blind = templates matching one of the 18 located SEALED families (used ONCE, for confirmation;
            no variant is chosen on them)."""
    C = pd.read_csv(f"{RES}/repro/classify15.csv")
    F = bostock_families()
    out = []
    for col, tag in (("kept_bkey", "dev"), ("hid_bkey", "blind")):
        S = C[C[col] >= 0][["grp", col, "n_det", "n_days"]].rename(columns={col: "bkey"}).copy()
        S["set"] = tag
        out.append(S)
    T = pd.concat(out, ignore_index=True)
    T = T[T.bkey.isin(F.index)]
    T["t_lat"] = F.lat[T.bkey].values
    T["t_lon"] = F.lon[T.bkey].values
    T["t_dep"] = F.dep[T.bkey].values
    return T


# ================================================================= bed adapters =================
class Bed:
    def __init__(self, mod, name, grid, corr_source):
        self.m = mod
        self.name = name
        self.grid = grid
        self.corr_source = corr_source

    def station_table(self):
        return self.m.station_table()

    def DayCache(self, stas):
        if hasattr(self.m, "DayCache"):
            return self.m.DayCache(stas)
        import k_common as K       # o_common rebinds k_common's bed facts in activate()
        return K.DayCache(stas)

    def days(self):
        return self.m.bed_days()


class BedMem(Bed):
    def station_table(self):
        return self.m.member_table()

    def DayCache(self, stas):
        import o_chestler as OC
        return OC.BandCache(stas)


def bed(name="north"):
    """north = the 2005 POLARIS bed (k_common); wa = the 2010 Olympic array bed (o_common)."""
    if name == "north":
        import k_common as M
        return Bed(M, "north", dict(lat=(47.30, 49.80), lon=(-125.80, -121.90), dep=(6.0, 76.0)),
                   "a5")
    import o_common as M
    sys.modules["n_common"] = M
    if name == "wa_mem":
        # Chestler's own chain (Arm B) correlates MEMBER STATIONS, not the array beams, so its
        # coordinates come from o_common.member_table and its data from o_chestler.BandCache
        # (2-8 Hz).  n_locate crashed on exactly this (KeyError 'B005', logs/o_armb_locate.log):
        # the 739 Arm-B families have NEVER been located.
        return BedMem(M, "wa_mem",
                      dict(lat=(46.30, 49.00), lon=(-125.60, -121.60), dep=(6.0, 76.0)),
                      "intraslab")
    return Bed(M, "wa", dict(lat=(46.30, 49.00), lon=(-125.60, -121.60), dep=(6.0, 76.0)),
               "intraslab")


def station_terms_published(b, stas):
    """The station terms the CURRENT locator uses, unchanged (a5 intraslab S corrections in the
    north, x_station_corrections in Washington).  Missing -> 0.0, flagged, never imputed."""
    if b.corr_source == "a5":
        import s_common as S
        try:
            c = S.station_corr(list(stas))
        except Exception:
            c = {}
        cS = {s: float(c.get(s, 0.0) or 0.0) for s in stas}
        A = pd.read_csv(f"{ROOT}/lfe_extension/results/a5_station_corrections.csv").set_index("sta")
        cP = {s: (float(A.corr_P[s]) if s in A.index and np.isfinite(A.corr_P.get(s, np.nan))
                  else 0.0) for s in stas}
        flag = {s for s in stas if abs(cS[s]) < 1e-12}
        return cS, cP, flag
    corr, flagged = b.m.station_corrections()
    cS = {s: float(corr.get(s, 0.0)) for s in stas}
    return cS, {s: cS[s] / 1.78 for s in stas}, set(flagged)


# ================================================================= travel times =================
_TT = {}


def tt_table(vpvs=1.78, xmax=320.0, dx=0.5, zmin=4.0, zmax=80.0, dz=0.5):
    """(xs, zs, TP, TS) Crosson 1-D first-arrival tables (p_vel1d, exact flat-layer shooting)."""
    key = (vpvs, xmax, dx, zmin, zmax, dz)
    if key in _TT:
        return _TT[key]
    cache = f"{OUT}/tt_{vpvs:.3f}_{xmax:.0f}_{dx}_{zmin}_{zmax}_{dz}.npz"
    if os.path.exists(cache):
        d = np.load(cache)
        _TT[key] = (d["xs"], d["zs"], d["TP"], d["TS"])
        return _TT[key]
    import p_vel1d as V
    xs = np.arange(0.0, xmax + dx, dx)
    zs = np.arange(zmin, zmax + dz, dz)
    TP = np.empty((zs.size, xs.size), np.float32)
    TS = np.empty((zs.size, xs.size), np.float32)
    for i, z in enumerate(zs):
        for j, x in enumerate(xs):
            TP[i, j] = V.tt(float(x), float(z), "P", vpvs=vpvs) if _acc(V) else V.tt(float(x), float(z), "P")
            TS[i, j] = V.tt(float(x), float(z), "S", vpvs=vpvs) if _acc(V) else V.tt(float(x), float(z), "S")
    np.savez(cache, xs=xs, zs=zs, TP=TP, TS=TS)
    _TT[key] = (xs, zs, TP, TS)
    return _TT[key]


def _acc(V):
    import inspect
    return "vpvs" in inspect.signature(V.tt).parameters


def interp2(xs, zs, T, x, z):
    xi = np.clip((np.asarray(x, float) - xs[0]) / (xs[1] - xs[0]), 0, xs.size - 1.001)
    zi = np.clip((np.asarray(z, float) - zs[0]) / (zs[1] - zs[0]), 0, zs.size - 1.001)
    x0 = np.floor(xi).astype(int)
    z0 = np.floor(zi).astype(int)
    fx = xi - x0
    fz = zi - z0
    return ((1 - fz) * ((1 - fx) * T[z0, x0] + fx * T[z0, x0 + 1])
            + fz * ((1 - fx) * T[z0 + 1, x0] + fx * T[z0 + 1, x0 + 1])).astype(np.float32)


class Grid:
    """Regular lat/lon/depth grid with the (nnode, nsta) S and P travel-time tables."""

    def __init__(self, stas, ST, lat, lon, dep, vpvs=1.78, termS=None, termP=None):
        self.stas = list(stas)
        self.lats = np.arange(lat[0], lat[1] + 1e-9, lat[2])
        self.lons = np.arange(lon[0], lon[1] + 1e-9, lon[2])
        self.deps = np.arange(dep[0], dep[1] + 1e-9, dep[2])
        LA, LO, DP = np.meshgrid(self.lats, self.lons, self.deps, indexing="ij")
        self.la, self.lo, self.dp = LA.ravel(), LO.ravel(), DP.ravel()
        xs, zs, TP, TS = tt_table(vpvs=vpvs)
        self.TS = np.empty((self.la.size, len(stas)), np.float32)
        self.TP = np.empty_like(self.TS)
        termS = termS or {}
        termP = termP or {}
        for k, s in enumerate(self.stas):
            d = epi_km(self.la, self.lo, float(ST.lat[s]), float(ST.lon[s]))
            self.TS[:, k] = interp2(xs, zs, TS, d, self.dp) + float(termS.get(s, 0.0))
            self.TP[:, k] = interp2(xs, zs, TP, d, self.dp) + float(termP.get(s, 0.0))
        self.shape = (self.lats.size, self.lons.size, self.deps.size)

    def edge(self, la, lo, dp):
        return bool(la <= self.lats[0] + 1e-6 or la >= self.lats[-1] - 1e-6
                    or lo <= self.lons[0] + 1e-6 or lo >= self.lons[-1] - 1e-6
                    or dp <= self.deps[0] + 1e-6 or dp >= self.deps[-1] - 1e-6)


def local_grid(stas, ST, la0, lo0, dp0, half=(0.06, 0.08, 4.0), step=(0.004, 0.005, 0.25),
               vpvs=1.78, termS=None, termP=None):
    return Grid(stas, ST,
                (la0 - half[0], la0 + half[0], step[0]),
                (lo0 - half[1], lo0 + half[1], step[1]),
                (max(dp0 - half[2], 4.0), dp0 + half[2], step[2]),
                vpvs=vpvs, termS=termS, termP=termP)


# ================================================================= misfit / search ==============
def misfit(TTsel, obs, w=None, robust=False):
    """Demeaned differential misfit at every node.  TTsel (nnode, k), obs (k,).

    L2 = the current locator's rms.  ROBUST = weighted median-centred mean |residual| (L1), which
    is what stops one 8-second pick outlier from dragging the solution onto a grid face."""
    obs = np.asarray(obs, np.float32)
    if w is None:
        w = np.ones(obs.size, np.float32)
    w = np.asarray(w, np.float32)
    w = w / w.sum()
    if robust:
        R = obs[None, :] - TTsel
        c = np.median(R, axis=1, keepdims=True)
        return (np.abs(R - c) * w[None, :]).sum(1)
    R = obs[None, :] - TTsel
    R = R - (R * w[None, :]).sum(1, keepdims=True)
    return np.sqrt(((R ** 2) * w[None, :]).sum(1))


def misfit_sp(grid, idx, tS, tP, sigS=0.30, sigSP=0.50, w=None, robust=False):
    """S differential misfit plus the S-P term at every station that has a P.  Chestler's two
    observable classes.  Returns the total cost (nnode,)."""
    ok = np.isfinite(tS)
    c = misfit(grid.TS[:, idx[ok]], tS[ok], None if w is None else w[ok], robust) / sigS
    m = np.isfinite(tP) & ok
    if m.any():
        sp = (grid.TS[:, idx[m]] - grid.TP[:, idx[m]]) - (tS[m] - tP[m])[None, :]
        add = np.abs(sp).mean(1) if robust else np.sqrt((sp ** 2).mean(1))
        c = c + add / sigSP
    return c


def shuffle_control(grid, idx, obs, w=None, robust=False, n=50, seed=3):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        p = rng.permutation(idx.size)
        out.append(float(misfit(grid.TS[:, idx[p]], obs, w, robust).min()))
    return np.array(out)


def best_node(grid, cost):
    i = int(np.argmin(cost))
    return float(grid.la[i]), float(grid.lo[i]), float(grid.dp[i]), float(cost[i]), i


# ================================================================= slab2 ========================
_SL = {}


def slab2(lat, lon):
    if "f" not in _SL:
        import p_common as P
        _SL["f"] = P.slab2_interp()
    return np.asarray(_SL["f"](np.asarray(lat), np.asarray(lon)), float)
