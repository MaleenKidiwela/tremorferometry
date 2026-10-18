#!/usr/bin/env python
"""Track A-2 step 5 — shared location machinery (copied/adapted from
lfe_catalog_repair/scripts/d3_stage4_locate.py; lfe_catalog_repair is READ-ONLY).

Differences from the D3 original, all deliberate (Merlin plan-review F10 + task spec):
  * the station set is the PERMANENT-network one and the two co-located CLUSTERS COLLAPSE TO ONE
    SITE each: PGC-site = {PGC, B009, B010, B011} (one travel time, PGC coordinates) and LRIV-site =
    {LRIV, B005, B006, B007}.  D3 treated them as independent stations, which is F10(i).
  * memory-lean grid search: the (nz, nnode, ksta) tensor is never materialised inside a solve, so a
    62-km / 1-km grid x 101 depths can be re-solved thousands of times (shuffled-pick null).
  * the cost function, its sigmas (SS 0.30 s, SSP 0.50 s, SZ 5 km) and the Slab2 prior are UNCHANGED
    from d3_stage4_locate so Run B is the same estimator D3 validated.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d3_vel1d as V                                          # noqa: E402

REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
SLAB2 = "/home/jovyan/tremorferometry/data/cas_slab2_input_04-18.csv"
SS, SSP, SZ = 0.30, 0.50, 5.0        # sigma_S (s), sigma_(S-P) (s), sigma_slab2 (km) -- as D3
LAT0, LON0 = 48.50, -123.75          # grid centre (fixed, NOT the family location: no truth leak)

# --- the SITE set (Merlin F10: PGC cluster = ONE station, LRIV cluster = ONE station) -------------
SITE_MEMBERS = {
    "PGC":  ["PGC", "B009", "B010", "B011"],
    "LRIV": ["LRIV", "B005", "B006", "B007"],
    "VGZ":  ["VGZ"], "SNB": ["SNB"], "NLLB": ["NLLB"],
    "B926": ["B926"], "B001": ["B001"], "B003": ["B003"], "B004": ["B004"],
}
STA2SITE = {m: q for q, mem in SITE_MEMBERS.items() for m in mem}


def station_table():
    S = pd.read_csv(f"{REP}/results/e5r_stations.csv").set_index("sta")
    return {s: (float(S.loc[s, "lat"]), float(S.loc[s, "lon"])) for s in S.index}


def site_ll(stab=None):
    """Site -> (lat, lon): the reference member's coordinates (PGC for the PGC cluster, LRIV for the
    LRIV cluster; the members are <1 km apart, far below the 5 km bar)."""
    stab = stab or station_table()
    return {q: stab[mem[0]] for q, mem in SITE_MEMBERS.items()}


def km_e_n(lat, lon, lat0=LAT0, lon0=LON0):
    return ((lon - lon0) * 111.2 * np.cos(np.radians(lat0)), (lat - lat0) * 111.2)


def slab2_depth(lats, lons):
    from scipy.interpolate import griddata
    sl = pd.read_csv(SLAB2, usecols=["lat", "lon", "depth", "etype"])
    sl = sl[~sl.etype.isin(["BA", "TO"])]
    sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna()
    pts = np.column_stack([np.atleast_1d(lons), np.atleast_1d(lats)])
    z = griddata(sl[["lon", "lat"]].values, sl.depth.values, pts, method="linear")
    zn = griddata(sl[["lon", "lat"]].values, sl.depth.values, pts, method="nearest")
    return np.where(np.isfinite(z), z, zn)


def build_tt_table(xmax=220.0, dx=1.0, zmin=4.0, zmax=70.0, dz=0.5, vpvs=V.VPVS_DEFAULT,
                   cache=f"{os.path.dirname(HERE)}/results/a5_tt_table.npz"):
    key = f"{xmax}_{dx}_{zmin}_{zmax}_{dz}_{vpvs}"
    if cache and os.path.exists(cache):
        d = np.load(cache)
        if str(d["key"]) == key:
            return d["xs"], d["zs"], d["TP"], d["TS"]
    xs = np.arange(0.0, xmax + dx, dx)
    zs = np.arange(zmin, zmax + dz, dz)
    TP = np.empty((zs.size, xs.size)); TS = np.empty((zs.size, xs.size))
    for i, z in enumerate(zs):
        for j, x in enumerate(xs):
            TP[i, j] = V.tt(x, z, "P", vpvs); TS[i, j] = V.tt(x, z, "S", vpvs)
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez(cache, xs=xs, zs=zs, TP=TP, TS=TS, key=np.array(key))
    return xs, zs, TP, TS


def interp_tt(tab, x, z):
    xs, zs, T = tab
    xi = np.clip((x - xs[0]) / (xs[1] - xs[0]), 0, xs.size - 1.001)
    zi = np.clip((z - zs[0]) / (zs[1] - zs[0]), 0, zs.size - 1.001)
    x0, z0 = np.floor(xi).astype(int), np.floor(zi).astype(int)
    fx, fz = xi - x0, zi - z0
    return ((1 - fz) * ((1 - fx) * T[z0, x0] + fx * T[z0, x0 + 1])
            + fz * ((1 - fx) * T[z0 + 1, x0] + fx * T[z0 + 1, x0 + 1]))


def grid_axes(halfwidth=62.0, dx=1.0, zmin=10.0, zmax=60.0, dz=0.5):
    ge = np.arange(-halfwidth, halfwidth + dx, dx)
    gn = np.arange(-halfwidth, halfwidth + dx, dx)
    E, Nn = np.meshgrid(ge, gn, indexing="ij")
    return E.ravel(), Nn.ravel(), np.arange(zmin, zmax + dz, dz)


def precompute_tt(E, Nn, gz, sites, sll, tabP, tabS):
    se = np.array([km_e_n(*sll[s])[0] for s in sites])
    sn = np.array([km_e_n(*sll[s])[1] for s in sites])
    d = np.sqrt((E[:, None] - se[None, :]) ** 2 + (Nn[:, None] - sn[None, :]) ** 2)
    TTS = np.empty((len(gz), d.shape[0], d.shape[1]), np.float32)
    TTP = np.empty_like(TTS)
    for i, z in enumerate(gz):
        zz = np.full_like(d, z)
        TTS[i] = interp_tt(tabS, d, zz); TTP[i] = interp_tt(tabP, d, zz)
    return TTS, TTP


def locate(tS, tP, sel, TTS, TTP, E, Nn, gz, z_slab, use_prior=True):
    """Origin-time-free grid search; identical cost to d3_stage4_locate.locate_from_tables but
    accumulated station by station (no (nz, nnode, k) temporary).

    tS, tP  arrays over the FULL site list (NaN where not picked); sel = indices to use.
    Returns (lat, lon, z, cost, n_sta, n_SP, at_bound)."""
    sel = np.asarray(sel, int)
    if sel.size < 3:
        return (np.nan,) * 4 + (int(sel.size), 0, 0)
    tS = np.asarray(tS, float); tP = np.asarray(tP, float)
    okS = np.isfinite(tS[sel])
    sel = sel[okS]
    if sel.size < 3:
        return (np.nan,) * 4 + (int(sel.size), 0, 0)
    dS_obs = (tS[sel] - tS[sel].mean()).astype(np.float32)
    Smean = np.zeros(TTS.shape[:2], np.float32)
    for j in sel:
        Smean += TTS[:, :, j]
    Smean /= len(sel)
    cost = np.zeros(TTS.shape[:2], np.float32)
    for a, j in enumerate(sel):
        r = (TTS[:, :, j] - Smean) - dS_obs[a]
        cost += (r * r) * np.float32(1.0 / SS ** 2)
    nsp = 0
    for j in sel:
        if np.isfinite(tP[j]):
            r = (TTS[:, :, j] - TTP[:, :, j]) - np.float32(tS[j] - tP[j])
            cost += (r * r) * np.float32(1.0 / SSP ** 2)
            nsp += 1
    if use_prior and np.isfinite(z_slab):
        cost += ((gz - z_slab) ** 2 / SZ ** 2).astype(np.float32)[:, None]
    iz, ih = np.unravel_index(int(np.argmin(cost)), cost.shape)
    e, n, z = float(E[ih]), float(Nn[ih]), float(gz[iz])
    at_b = int(iz == 0 or iz == len(gz) - 1
               or abs(e) >= E.max() - 1e-6 or abs(n) >= Nn.max() - 1e-6)
    return (LAT0 + n / 111.2, LON0 + e / (111.2 * np.cos(np.radians(LAT0))), z,
            float(cost[iz, ih]), int(len(sel)), nsp, at_b)


def azimuth_gap(lat, lon, sites, sll):
    """Largest azimuthal gap (deg) seen from (lat, lon) by the given sites."""
    if len(sites) < 2:
        return 360.0
    az = sorted(np.degrees(np.arctan2(*km_e_n(*sll[s], lat, lon))) % 360 for s in sites)
    d = np.diff(az + [az[0] + 360.0])
    return float(d.max())


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.5))
    return float(np.hypot((lon1 - lon2) * 111.195 * cs, (lat1 - lat2) * 111.195))
