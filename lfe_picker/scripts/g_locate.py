#!/usr/bin/env python
"""Gate 2 location machinery — the a5 estimator (which is d3_stage4_locate's cost function,
verified identical there) with a 2-D INTERFACE mode added.

Copied and adapted from lfe_extension/scripts/a5_locate.py (read-only source).  Differences:
  * INTERFACE mode: the depth of every horizontal node is fixed to Slab2 at that node, so the
    search is over a 2-D interface surface (Merlin's P0 grid / Gate 2 "2-D interface grid").
    FREE-DEPTH mode (the 3-D search of a5) is kept as the sanity run.
  * region-agnostic: grid centre, half-width and the site list come from the caller, not module
    constants, so the same code screens Vancouver Island and, later, Puget / Oregon.
  * the misfit is returned as an RMS S-differential RESIDUAL in seconds (a truth-free statistic),
    which is what the shuffled-pick null is compared on.
Cost function and sigmas (SS 0.30 s, SSP 0.50 s, SZ 5 km) are UNCHANGED from a5/d3.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402

SS, SSP, SZ = 0.30, 0.50, 5.0


class Grid:
    """Horizontal grid + per-site travel times, built once per region."""

    def __init__(self, sites, site_ll, lat0, lon0, halfwidth=62.0, dx=1.0,
                 zmin=10.0, zmax=60.0, dz=0.5):
        self.sites = list(sites)
        self.lat0, self.lon0 = float(lat0), float(lon0)
        ge = np.arange(-halfwidth, halfwidth + dx, dx)
        E, N = np.meshgrid(ge, ge, indexing="ij")
        self.E, self.N = E.ravel(), N.ravel()
        self.gz = np.arange(zmin, zmax + dz, dz)
        self.lat = self.lat0 + self.N / 111.2
        self.lon = self.lon0 + self.E / (111.2 * np.cos(np.radians(self.lat0)))
        se = np.array([self._e(site_ll[s][0], site_ll[s][1]) for s in self.sites])
        sn = np.array([self._n(site_ll[s][0], site_ll[s][1]) for s in self.sites])
        self.dist = np.sqrt((self.E[:, None] - se[None, :]) ** 2 + (self.N[:, None] - sn[None, :]) ** 2)
        xs, zs, TP, TS = G.tt_table()
        # 3-D tensors (nz, nnode, nsite)
        self.TTS = np.empty((len(self.gz), *self.dist.shape), np.float32)
        self.TTP = np.empty_like(self.TTS)
        for i, z in enumerate(self.gz):
            zz = np.full_like(self.dist, z)
            self.TTS[i] = G._interp2(xs, zs, TS, self.dist, zz)
            self.TTP[i] = G._interp2(xs, zs, TP, self.dist, zz)
        # interface (2-D) layer: z = Slab2 at the node
        self.z_if = G.slab2_depth(self.lat, self.lon).astype(np.float32)
        zz = np.repeat(self.z_if[:, None], self.dist.shape[1], axis=1)
        self.TTS_if = G._interp2(xs, zs, TS, self.dist, zz).astype(np.float32)
        self.TTP_if = G._interp2(xs, zs, TP, self.dist, zz).astype(np.float32)

    def _e(self, lat, lon):
        return (lon - self.lon0) * 111.2 * np.cos(np.radians(self.lat0))

    def _n(self, lat, lon):
        return (lat - self.lat0) * 111.2

    def inside(self, lat, lon):
        e, n = self._e(lat, lon), self._n(lat, lon)
        return bool(abs(e) <= self.E.max() and abs(n) <= self.N.max())


def _resid_rms(tS, sel, TTS_sel):
    """RMS of the S-differential residual at every node (nnode,) for a 2-D TTS_sel (nnode, k)."""
    dS = tS[sel] - tS[sel].mean()
    pred = TTS_sel - TTS_sel.mean(1, keepdims=True)
    r = pred - dS[None, :]
    return np.sqrt((r ** 2).mean(1))


def locate_interface(grid, tS, tP, sel):
    """2-D interface search.  Returns dict(lat, lon, z, rms_s, n_sites, n_SP, at_bound)."""
    sel = np.asarray(sel, int)
    tS = np.asarray(tS, float); tP = np.asarray(tP, float)
    sel = sel[np.isfinite(tS[sel])]
    if sel.size < 3:
        return None
    TTS = grid.TTS_if[:, sel]
    dS = (tS[sel] - tS[sel].mean()).astype(np.float32)
    pred = TTS - TTS.mean(1, keepdims=True)
    r = pred - dS[None, :]
    cost = (r ** 2).sum(1) / SS ** 2
    nsp = 0
    for a, j in enumerate(sel):
        if np.isfinite(tP[j]):
            rr = (grid.TTS_if[:, j] - grid.TTP_if[:, j]) - np.float32(tS[j] - tP[j])
            cost = cost + (rr ** 2) / SSP ** 2
            nsp += 1
    ih = int(np.argmin(cost))
    rms = float(np.sqrt((r[ih] ** 2).mean()))
    at_b = int(abs(grid.E[ih]) >= grid.E.max() - 1e-6 or abs(grid.N[ih]) >= grid.N.max() - 1e-6)
    return dict(lat=float(grid.lat[ih]), lon=float(grid.lon[ih]), z=float(grid.z_if[ih]),
                rms_s=rms, cost=float(cost[ih]), n_sites=int(sel.size), n_SP=nsp, at_bound=at_b)


def locate_free(grid, tS, tP, sel, z_prior=np.nan, use_prior=False):
    """3-D free-depth search (the a5/d3 estimator, accumulated station by station)."""
    sel = np.asarray(sel, int)
    tS = np.asarray(tS, float); tP = np.asarray(tP, float)
    sel = sel[np.isfinite(tS[sel])]
    if sel.size < 3:
        return None
    dS = (tS[sel] - tS[sel].mean()).astype(np.float32)
    Smean = np.zeros(grid.TTS.shape[:2], np.float32)
    for j in sel:
        Smean += grid.TTS[:, :, j]
    Smean /= len(sel)
    cost = np.zeros(grid.TTS.shape[:2], np.float32)
    for a, j in enumerate(sel):
        r = (grid.TTS[:, :, j] - Smean) - dS[a]
        cost += (r * r) * np.float32(1.0 / SS ** 2)
    nsp = 0
    for j in sel:
        if np.isfinite(tP[j]):
            r = (grid.TTS[:, :, j] - grid.TTP[:, :, j]) - np.float32(tS[j] - tP[j])
            cost += (r * r) * np.float32(1.0 / SSP ** 2)
            nsp += 1
    if use_prior and np.isfinite(z_prior):
        cost += ((grid.gz - z_prior) ** 2 / SZ ** 2).astype(np.float32)[:, None]
    iz, ih = np.unravel_index(int(np.argmin(cost)), cost.shape)
    pred = grid.TTS[iz, ih, sel]
    r = (pred - pred.mean()) - dS
    at_b = int(iz == 0 or iz == len(grid.gz) - 1
               or abs(grid.E[ih]) >= grid.E.max() - 1e-6 or abs(grid.N[ih]) >= grid.N.max() - 1e-6)
    return dict(lat=float(grid.lat[ih]), lon=float(grid.lon[ih]), z=float(grid.gz[iz]),
                rms_s=float(np.sqrt((r ** 2).mean())), cost=float(cost[iz, ih]),
                n_sites=int(sel.size), n_SP=nsp, at_bound=at_b)


def azimuth_gap(lat, lon, sites, site_ll):
    if len(sites) < 2:
        return 360.0
    az = sorted(np.degrees(np.arctan2(
        (site_ll[s][1] - lon) * 111.2 * np.cos(np.radians(lat)),
        (site_ll[s][0] - lat) * 111.2)) % 360 for s in sites)
    d = np.diff(az + [az[0] + 360.0])
    return float(d.max())


def synthetic_floor(grid, site_ll, sites_used, lat, lon, z, tP_mask, n_real=50,
                    sigma_S=0.3, static_S=0.5, seed=0, mode="interface"):
    """Forward-model (lat, lon, z) through Crosson to the sites actually picked, add 0.3 s S pick
    noise and a per-site static U(-0.5, +0.5) s (scaled by 1/1.78 for P), relocate, and return the
    median / p90 horizontal error.  a5_synth_floor recipe, restricted to this proposal's geometry."""
    rng = np.random.default_rng(seed)
    xs, zs, TP, TS = G.tt_table()
    idx = [grid.sites.index(s) for s in sites_used]
    d = np.array([np.hypot(grid._e(lat, lon) - grid._e(*site_ll[s]),
                           grid._n(lat, lon) - grid._n(*site_ll[s])) for s in sites_used])
    ttS = G._interp2(xs, zs, TS, d, np.full_like(d, z))
    ttP = G._interp2(xs, zs, TP, d, np.full_like(d, z))
    dh = []
    for _ in range(n_real):
        st = rng.uniform(-static_S, static_S, len(sites_used))
        tS = np.full(len(grid.sites), np.nan)
        tP = np.full(len(grid.sites), np.nan)
        for a, j in enumerate(idx):
            tS[j] = ttS[a] + rng.normal(0, sigma_S) + st[a]
            if tP_mask[a]:
                tP[j] = ttP[a] + rng.normal(0, sigma_S / 1.78) + st[a] / 1.78
        r = (locate_interface(grid, tS, tP, idx) if mode == "interface"
             else locate_free(grid, tS, tP, idx))
        if r is not None:
            dh.append(G.epi_km(r["lat"], r["lon"], lat, lon))
    dh = np.asarray(dh, float)
    return (float(np.median(dh)) if dh.size else np.nan,
            float(np.quantile(dh, 0.9)) if dh.size else np.nan, dh)
