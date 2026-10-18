#!/usr/bin/env python
"""1-D layered velocity model + travel times for the D3 stage-4 family locations.

WHY A NEW MODEL FILE: the project has NO velocity model. The only wave speed anywhere in the
pipeline is the single homogeneous coda S velocity BETA = 3.5 km/s
(`fault_tomography/inversion/build_G_captured.py:13`, `fault_tomography/kernels/residence_kernel_mc.py:41`),
which is a coda-residence scalar, not a location model. Stage 4 therefore uses an explicitly stated
external 1-D model.

MODEL: the Crosson (1976) western-Washington layered P model, which is the PNSN/UW standard for
this region (the "P3/P4" location models derive from it), with Vp/Vs = 1.78.
    depth (km)   Vp (km/s)
      0 -  4        5.40
      4 -  9        6.38
      9 - 16        6.59
     16 - 20        6.73
     20 - 25        6.86
     25 - 32        7.00
      > 32          7.90     (mantle)
CAVEATS, stated once and carried into the results:
  * a 1-D model over a subduction zone biases interface depths by a few km (the real structure has a
    low-velocity overlying wedge and a high-velocity slab); Slab2 is therefore used as a weak depth
    prior, and the depth uncertainty reported is the FORMAL one from the station bootstrap, not an
    estimate of that systematic;
  * Vp/Vs is elevated (1.8-2.0) in parts of the Cascadia forearc; sensitivity to Vp/Vs = 1.75 / 1.82
    is reported alongside the preferred solution.

Travel times are computed by exact shooting for a flat-layered model (Snell ray parameter solved by
bisection on the horizontal distance), including the direct/refracted first arrival. Sources are at
depth, receivers at the surface, so times are computed for a downgoing-then-upgoing single leg
(source below receiver = one turning-free segment): the standard flat-layer integral
   x(p) = sum_i  dz_i * p*v_i / sqrt(1 - (p*v_i)^2)
   t(p) = sum_i  dz_i / (v_i * sqrt(1 - (p*v_i)^2))
over the layers between the receiver (z=0) and the source depth, plus a head-wave/refracted branch
check (t_refr = x/v_below + intercepts) so the first arrival is taken as the minimum of the two.
"""
from __future__ import annotations

import numpy as np

# layer tops (km) and Vp (km/s)
ZTOP = np.array([0.0, 4.0, 9.0, 16.0, 20.0, 25.0, 32.0])
VP = np.array([5.40, 6.38, 6.59, 6.73, 6.86, 7.00, 7.90])
VPVS_DEFAULT = 1.78


def layer_profile(z_src, vpvs=VPVS_DEFAULT, phase="S"):
    """Layer thicknesses and velocities between the surface and z_src."""
    v = VP / (vpvs if phase.upper() == "S" else 1.0)
    tops = ZTOP
    zs = float(z_src)
    dz, vv = [], []
    for i in range(len(tops)):
        top = tops[i]
        bot = tops[i + 1] if i + 1 < len(tops) else np.inf
        if zs <= top:
            break
        h = min(bot, zs) - top
        if h > 0:
            dz.append(h)
            vv.append(v[i])
    return np.asarray(dz), np.asarray(vv)


def _xt_of_p(p, dz, vv):
    s = p * vv
    s = np.clip(s, 0.0, 1.0 - 1e-12)
    c = np.sqrt(1.0 - s * s)
    x = float((dz * s / c).sum())
    t = float((dz / (vv * c)).sum())
    return x, t


def direct_time(x_km, z_src, vpvs=VPVS_DEFAULT, phase="S"):
    """First-arrival time for the direct (non-turning) ray, by bisection on the ray parameter."""
    dz, vv = layer_profile(z_src, vpvs, phase)
    if dz.size == 0:
        return np.nan
    if x_km <= 1e-6:
        return float((dz / vv).sum())
    pmax = 1.0 / vv.max() * (1 - 1e-9)
    lo, hi = 0.0, pmax
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        x, _t = _xt_of_p(mid, dz, vv)
        if x < x_km:
            lo = mid
        else:
            hi = mid
    _x, t = _xt_of_p(0.5 * (lo + hi), dz, vv)
    return t


def refracted_time(x_km, z_src, vpvs=VPVS_DEFAULT, phase="S"):
    """Head-wave branches along each interface BELOW the source (mantle refraction etc.)."""
    v = VP / (vpvs if phase.upper() == "S" else 1.0)
    best = np.inf
    dz, vv = layer_profile(z_src, vpvs, phase)
    if dz.size == 0:
        return np.inf
    nsrc = dz.size
    for j in range(nsrc, len(v)):          # refractor = top of layer j (below the source)
        vr = v[j]
        # legs: receiver side (surface->refractor) and source side (source->refractor)
        tops = list(ZTOP) + [np.inf]
        # receiver-side layers 0..j-1 full thickness
        rz, rv = [], []
        for i in range(j):
            h = tops[i + 1] - tops[i]
            rz.append(h)
            rv.append(v[i])
        # source-side layers from z_src down to the refractor
        sz, sv = [], []
        for i in range(len(v)):
            top, bot = tops[i], tops[i + 1]
            a = max(top, z_src)
            b = min(bot, tops[j])
            if b > a:
                sz.append(b - a)
                sv.append(v[i])
        rz, rv = np.asarray(rz), np.asarray(rv)
        sz, sv = np.asarray(sz), np.asarray(sv)
        if np.any(rv >= vr) or (sv.size and np.any(sv >= vr)):
            continue
        def leg(zz, vvv):
            if zz.size == 0:
                return 0.0, 0.0
            c = np.sqrt(np.maximum(1.0 - (vvv / vr) ** 2, 1e-12))
            return float((zz * (vvv / vr) / c).sum()), float((zz / (vvv * c)).sum())
        xr, tr = leg(rz, rv)
        xsq, tsq = leg(sz, sv)
        xoff = xr + xsq
        if x_km < xoff:
            continue
        t = tr + tsq + (x_km - xoff) / vr
        best = min(best, t)
    return best


def tt(x_km, z_src, phase="S", vpvs=VPVS_DEFAULT):
    """First-arrival travel time (s) for horizontal offset x_km and source depth z_src."""
    a = direct_time(x_km, z_src, vpvs, phase)
    b = refracted_time(x_km, z_src, vpvs, phase)
    return float(min(a, b))


def tt_grid(x_km, z_km, phase="S", vpvs=VPVS_DEFAULT):
    """Vectorised over arrays of (x, z) by looping (grids here are small: <1e5 points)."""
    x_km = np.atleast_1d(np.asarray(x_km, float))
    z_km = np.atleast_1d(np.asarray(z_km, float))
    out = np.empty(np.broadcast(x_km, z_km).shape)
    it = np.nditer([x_km, z_km, out], op_flags=[["readonly"], ["readonly"], ["writeonly"]])
    for a, b, o in it:
        o[...] = tt(float(a), float(b), phase, vpvs)
    return out


def slab2_depth(lats, lons):
    """Slab2 interface depth (km) at (lat, lon), the project's own source and filtering."""
    import pandas as pd
    from scipy.interpolate import griddata
    sl = pd.read_csv("data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
    sl = sl[~sl.etype.isin(["BA", "TO"])]
    sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna()
    pts = np.column_stack([np.atleast_1d(lons), np.atleast_1d(lats)])
    z = griddata(sl[["lon", "lat"]].values, sl.depth.values, pts, method="linear")
    zn = griddata(sl[["lon", "lat"]].values, sl.depth.values, pts, method="nearest")
    return np.where(np.isfinite(z), z, zn)


if __name__ == "__main__":
    print("Crosson (1976) western Washington, Vp/Vs = %.2f" % VPVS_DEFAULT)
    for z in (20, 30, 40, 50):
        for x in (0, 20, 50, 100):
            print(f"  z={z:3d} km x={x:4d} km  tP={tt(x, z, 'P'):6.2f} s  "
                  f"tS={tt(x, z, 'S'):6.2f} s  S-P={tt(x, z, 'S')-tt(x, z, 'P'):5.2f} s")
