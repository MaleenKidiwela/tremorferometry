#!/usr/bin/env python
"""STEP 2 (Merlin recipe): build the forward operator G (single-scatter ellipsoid, 2-4 s window) over the
certified (cell,station) pairs, AND the per-row on-fault CAPTURED FRACTION (F5) evaluated on a coarse 3-D
volume grid — the check that says whether the interface-only model captures the coda sensitivity (and whether
the historic 0.77 is inflated). Reads res_catalog/{pairs.csv,cells.csv}. Writes res_catalog/G.npz.
G rows = normalized single-scatter kernel over the 872 fault cells (sign-flipped: dtau/tau = -K.dbeta/beta),
matching production forward.py. Captured fraction uses ONLY the volume grid (uniform dV cancels)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "fault_tomography/kernels")
from scipy.interpolate import griddata

BETA, ELL, W1, W2 = 3.5, 40.0, 2.0, 4.0
import os
OUT = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
pairs = pd.read_csv(f"{OUT}/pairs.csv")
cells = pd.read_csv(f"{OUT}/cells.csv").sort_values("cell").reset_index(drop=True)
cidx = {c: i for i, c in enumerate(cells.cell)}
M = len(cells)

# local flat projection (km), origin = mean of cells
lat0, lon0 = cells.cell_lat.mean(), cells.cell_lon.mean()
def proj(lat, lon):
    return np.column_stack([(np.asarray(lon) - lon0) * 111 * np.cos(np.radians(lat0)),
                            (np.asarray(lat) - lat0) * 111])
CXY = proj(cells.src_lat.values, cells.src_lon.values)       # (M,2) cell scatterer horizontal = LFE family centroid
CZ = cells.depth_km.values                                    # (M,) cell depth = family-median interface depth

# ---- coarse 3-D volume grid for the captured-fraction (uniform dV) ----
sl = pd.read_csv("data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]; sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna()
HSLAB = 1.5   # interface half-thickness (km): production H_THICK=3.0 -> +/-1.5 (Merlin: f and h must match)
def build_volume(hdeg, dz):
    mlat = np.arange(cells.cell_lat.min() - 0.5, cells.cell_lat.max() + 0.5 + 1e-9, hdeg)
    mlon = np.arange(cells.cell_lon.min() - 0.5, cells.cell_lon.max() + 0.5 + 1e-9, hdeg)
    zl = np.arange(dz, 57.0 + 1e-9, dz)
    GLA, GLO = np.meshgrid(mlat, mlon, indexing="ij")
    vlat = GLA.ravel(); vlon = GLO.ravel()
    vs = griddata(sl[["lon", "lat"]].values, sl.depth.values, np.column_stack([vlon, vlat]), method="linear")
    vsn = griddata(sl[["lon", "lat"]].values, sl.depth.values, np.column_stack([vlon, vlat]), method="nearest")
    vs = np.where(np.isfinite(vs), vs, vsn)
    VXY = proj(np.repeat(vlat, len(zl)), np.repeat(vlon, len(zl)))
    VZ = np.tile(zl, len(vlat)); VSLAB = np.repeat(vs, len(zl))
    return VXY, VZ, (np.abs(VZ - VSLAB) <= HSLAB)
VXY, VZ, IFACE = build_volume(0.2, 1.0)                       # base grid: 0.2 deg x 1 km layers
print(f"cells M={M} | volume pts={len(VZ)} | interface(+/-{HSLAB}km) pts={IFACE.sum()}")


def raw_ss(pxy, pz, s_h, d, r_h):
    rs_h2 = ((pxy - s_h) ** 2).sum(1)
    Rr = np.sqrt(((pxy - r_h) ** 2).sum(1) + pz ** 2)
    direct = np.sqrt(((r_h - s_h) ** 2).sum() + d ** 2)
    lo, hi = direct + BETA * W1, direct + BETA * W2
    K = np.zeros(len(pz))
    for zsrc in (d, -d):
        Rs = np.sqrt(rs_h2 + (pz - zsrc) ** 2)
        path = Rs + Rr
        m = (path >= lo) & (path <= hi)
        K[m] += np.exp(-path[m] / ELL) / (Rs[m] ** 2 * Rr[m] ** 2)
    return K, direct


G = np.zeros((len(pairs), M), np.float32)
captured = np.zeros(len(pairs)); offset = np.zeros(len(pairs))
sxy = proj(pairs.sta_lat.values, pairs.sta_lon.values)
psrc = proj(pairs.src_lat.values, pairs.src_lon.values)       # per-pair SOURCE = that station's family centroid
for i, r in pairs.iterrows():
    ci = cidx[r.cell]
    s_h = psrc[i]; d = float(r.src_dep); r_h = sxy[i]
    kf, direct = raw_ss(CXY, CZ, s_h, d, r_h)                 # over fault cells -> G row (normalized)
    tot = kf.sum()
    if tot > 0:
        G[i] = -(kf / tot)
    kv, _ = raw_ss(VXY, VZ, s_h, d, r_h)                      # over volume -> captured fraction
    sv = kv.sum()
    captured[i] = (kv[IFACE].sum() / sv) if sv > 0 else 0.0
    offset[i] = np.sqrt(((r_h - s_h) ** 2).sum())             # source-receiver horizontal offset (km)

nz = (np.abs(G).sum(1) > 0)
print(f"G shape {G.shape} | all-zero rows: {(~nz).sum()} ({100*(~nz).mean():.1f}%)")
print(f"captured fraction (+/-{HSLAB}km, base grid): median {np.median(captured):.3f}  p25 {np.percentile(captured,25):.3f}  "
      f">5% {100*(captured>0.05).mean():.0f}%  >10% {100*(captured>0.10).mean():.0f}%")
print(f"source-receiver offset km: median {np.median(offset):.0f}  p90 {np.percentile(offset,90):.0f}  max {offset.max():.0f}")

# ---- grid-convergence check (Merlin): recompute f on a REFINED grid for a subsample; accept if within ~20% ----
VXY2, VZ2, IFACE2 = build_volume(0.1, 0.5)                    # halve both spacings
rng = np.random.RandomState(0); sub = rng.choice(len(pairs), size=min(200, len(pairs)), replace=False)
cf_base, cf_ref = [], []
for i in sub:
    s_h = psrc[i]; d = float(pairs.src_dep.iloc[i]); r_h = sxy[i]
    kv, _ = raw_ss(VXY, VZ, s_h, d, r_h);   cf_base.append(kv[IFACE].sum()/kv.sum() if kv.sum()>0 else 0)
    kv2, _ = raw_ss(VXY2, VZ2, s_h, d, r_h); cf_ref.append(kv2[IFACE2].sum()/kv2.sum() if kv2.sum()>0 else 0)
mb, mr = np.median(cf_base), np.median(cf_ref)
print(f"[convergence] median f base {mb:.3f} vs refined {mr:.3f}  -> ratio {mr/mb if mb else float('nan'):.2f} "
      f"({'PASS <20%' if mb and abs(mr/mb-1)<0.20 else 'CHECK'})")

np.savez_compressed(f"{OUT}/G.npz", G=G, captured=captured.astype(np.float32), offset=offset.astype(np.float32),
                    cell=cells.cell.values, cell_lat=cells.cell_lat.values, cell_lon=cells.cell_lon.values,
                    depth_km=CZ, frac_deep=cells.frac_deep.values, depth_mixed=cells.depth_mixed.values.astype(bool),
                    src_lat=cells.src_lat.values, src_lon=cells.src_lon.values, cxy=CXY, lat0=lat0, lon0=lon0,
                    pair_tag=pairs.tag.values, pair_cell=pairs.cell.values, pair_site=pairs.site.values,
                    pair_sigma=pairs.sigma_daily.values)
print(f"wrote {OUT}/G.npz")
