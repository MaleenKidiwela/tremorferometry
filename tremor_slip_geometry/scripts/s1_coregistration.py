#!/usr/bin/env python
"""S1 — per-event co-registration of PNSN tremor and Gualandi (2025) modeled slip, 2010-2026.

Measures, for every confirmed slip patch, WHERE the tremor sits relative to the modeled
slip (across-strike, along-strike, in depth, at the depth edges) and how much slip there
is per tremor detection; then asks whether any of those quantities drifts over 16 years.

--------------------------------------------------------------------------------------
INPUTS (all read-only, all outside this folder)
--------------------------------------------------------------------------------------
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.csv         thrust patch catalog
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz         label (T x 742, int; 0=none,
                                                                          else patch id), inc30_mm,
                                                                          days (T, 'YYYY-MM-DD'),
                                                                          patch_lat/lon/depth_km/area_km2 (742)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm_mirror.{csv,npz}   back-slip control, same structure
  fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv   per-patch tremor confirmation
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  key 'A_raw' (T x 13),
                                                                          same daily grid as `days`
  catalogs/pnsn_tremor_cascadia_full.csv                                  PNSN tremor (time, lat, lon;
                                                                          the depth column is a placeholder
                                                                          and is IGNORED)
  /home/jovyan/cascadia-sse-strain  (mesh + sse_strain package)           Ld (742 x 13) dip-slip patterns
                                                                          [m per unit amplitude];
                                                                          Slab2 depth interpolator f(lon, lat)

  Thrust slip per triangle, cumulative, metres:   raw_cum = -(Ld @ A_raw.T)   (742 x T)
  Mirror (back-slip) per triangle, metres:        raw_cum = +(Ld @ A_raw.T)   (742 x T)

  EVENT SETS
    CONFIRMED = tremor_vs_slip_patches: catalog=='thrust' & confirmed & seg in {north, central}
                & Mw >= 6.0 & not broad_pulse
    CONTROL   = tremor_vs_slip_patches: catalog=='mirror'  & Mw >= 6.0 & seg in {north, central}

--------------------------------------------------------------------------------------
EXACT DEFINITIONS (per patch p, id = pid, start/end dates mapped to indices ts, te on `days`)
--------------------------------------------------------------------------------------
  Local km frame:   x = (lon + 124) * 111 * cos(47 deg);   y = (lat - 47) * 111.

  Footprint F     = triangles whose label == pid on at least one day in [ts, te].
  acc (mm, 742)   = 1e3 * (raw_cum[:, te] - raw_cum[:, max(ts - 30, 0)]), clipped at >= 0.
                    (slip accumulated over the event, including the 30-d trailing window that
                     precedes the threshold crossing).
  Slip centroid   = area_km2 * acc weighted centroid of the triangle centroids over F (km frame).
  Slip depth      = same-weight (area*acc) weighted mean of Slab2 depth at the triangle centroids over F.
  Slip edges      = among ALL triangles whose latitude lies within F's latitude span +- 0.1 deg AND
                    whose acc >= 3 mm: up-dip edge = shallowest such triangle's depth;
                    down-dip edge = deepest such triangle's depth (km).
  Tremor during DUR = PNSN detections with date in [days[ts], days[te]] that lie within 40 km
                    (km frame, epicentral) of at least one F triangle centroid.
                    N = count.  Tremor centroid = mean of their km coordinates.
                    Tremor depth of a detection = Slab2 interface depth at its epicentre;
                    med / p10 / p90 over the finite values.
  Local frame at the slip centroid:
                    depth gradient of f by central differences of +- 0.05 deg, converted to km;
                    d = unit vector along the gradient (direction of INCREASING depth, i.e. down-dip);
                    s = d rotated 90 deg, sign chosen so that s has a positive north component.
  Q1 (km) = (tremor centroid - slip centroid) . d     (+ = tremor DOWN-DIP of slip)
  Q2 (km) = (tremor centroid - slip centroid) . s     (+ = tremor NORTH along strike)
  Q3 (km) = tremor median depth - slip-weighted depth
  Q4a(km) = tremor p10 depth - slip up-dip edge depth
  Q4b(km) = tremor p90 depth - slip down-dip edge depth
  potency (km^3) = sum over F of area_km2 * acc_mm * 1e-6
  Q5      = potency / N   (reported as log10 Q5)

  CONTROLS
    (a) mirror set: identical computation on the mirror catalog with back-slip positive.
    (b) time-shifted tremor: for each CONFIRMED patch, 100 rigid shifts of the tremor dates by a
        random integer of >= 365 and <= (span - 365) days, wrapped inside
        2010-01-01 .. 2026-05-27; the spatial 40-km mask is unchanged. Q1 and Q3 are recomputed
        per shift -> null mean/sd -> z1 = (Q1 - mean)/sd, z3 = (Q3 - mean)/sd.

  TRENDS (confirmed set only; per segment and combined)
    Theil-Sen slope of Q1, Q3, Q4a, Q4b, log10 Q5 against event start date (units per YEAR);
    p from 1000 permutations of the dates among events (two-sided on |slope|);
    Spearman rho with its own asymptotic p.

--------------------------------------------------------------------------------------
OUTPUTS (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------
  results/s1_events.csv     one row per confirmed patch: id, start, end, seg, Mw, N, potency_km3,
                            slip_depth_km, tremor_depth_med/p10/p90, slip_updip_edge_km,
                            slip_downdip_edge_km, Q1..Q5 (Q4 as Q4a/Q4b), log10Q5, z1, z3
  results/s1_controls.csv   mirror rows, same columns minus z1/z3
  results/s1_trends.csv     one row per (set, segment, quantity): n, Theil-Sen slope/yr, CI,
                            permutation p, Spearman rho/p
  figures/s1_offsets_vs_time.png   4 panels vs date (Q1; Q3; Q4a & Q4b; log10 Q5),
                            north blue / central green / mirror grey, Theil-Sen line, slope+p in titles

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s1_coregistration.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import binomtest, spearmanr, theilslopes, ttest_1samp

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
OUT_RES = ROOT / "tremor_slip_geometry/results"
OUT_FIG = ROOT / "tremor_slip_geometry/figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)
N_SHIFT = 100
N_PERM = 1000
TREMOR_RADIUS_KM = 40.0
EDGE_ACC_MM = 3.0
PRE_DAYS = 30
COS47 = np.cos(np.deg2rad(47.0))

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    """Local km frame: x east of 124W at 47N, y north of 47N."""
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def to_lonlat(x, y):
    return x / (111.0 * COS47) - 124.0, y / 111.0 + 47.0


# ----------------------------------------------------------------------------------
# 1. Slip model
# ----------------------------------------------------------------------------------
log("loading slip model (load_solution + invert, ~1 min) ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()  # Ld: (742, 13) metres per unit amplitude
slabf = slab_depth_interpolator(fault)
log(f"  Ld {Ld.shape}, fault.n_patch {fault.n_patch}")

A_raw = np.nan_to_num(
    np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"]
)
CUM_THRUST = -(Ld @ A_raw.T)  # (742, T) metres, thrust positive
CUM_MIRROR = +(Ld @ A_raw.T)  # (742, T) metres, back-slip positive
log(f"  A_raw {A_raw.shape}; cum slip {CUM_THRUST.shape} m")


# ----------------------------------------------------------------------------------
# 2. Patch label cubes and geometry
# ----------------------------------------------------------------------------------
NPZ = {
    "thrust": np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True),
    "mirror": np.load(SSE / "slip_patches_5mm_mirror.npz", allow_pickle=True),
}
LABEL = {k: v["label"] for k, v in NPZ.items()}
DAYS = np.asarray(NPZ["thrust"]["days"], dtype=str)
assert (np.asarray(NPZ["mirror"]["days"], dtype=str) == DAYS).all(), "day grids differ"
assert len(DAYS) == A_raw.shape[0], "A_raw / days length mismatch"

PLAT = NPZ["thrust"]["patch_lat"]
PLON = NPZ["thrust"]["patch_lon"]
PAREA = NPZ["thrust"]["patch_area_km2"]
PDEPTH = slabf(PLON, PLAT)  # identical to npz patch_depth_km (checked: max|diff| 1.8e-15)
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
DAY_DT = pd.to_datetime(DAYS)
log(f"  {len(DAYS)} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles")

CUM = {"thrust": CUM_THRUST, "mirror": CUM_MIRROR}


# ----------------------------------------------------------------------------------
# 3. Tremor
# ----------------------------------------------------------------------------------
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DEPTH = slabf(tr["lon"].values, tr["lat"].values)  # Slab2 depth at epicentre; NaN off-mesh
TR_DAY = (tr["t"].values.astype("datetime64[D]") - np.datetime64("2010-01-01")).astype(int)

SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1  # inclusive span, days
TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(
    f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}; "
    f"Slab2 depth finite for {np.isfinite(TR_DEPTH).mean() * 100:.1f}%; wrap span {SPAN} d"
)


# ----------------------------------------------------------------------------------
# 4. Event sets
# ----------------------------------------------------------------------------------
tvs = pd.read_csv(SSE / "tremor_vs_slip_patches.csv")
confirmed = tvs[
    (tvs.catalog == "thrust")
    & (tvs.confirmed.astype(bool))
    & (tvs.seg.isin(["north", "central"]))
    & (tvs.Mw >= 6.0)
    & (~tvs.broad_pulse.astype(bool))
].copy()
control = tvs[
    (tvs.catalog == "mirror") & (tvs.Mw >= 6.0) & (tvs.seg.isin(["north", "central"]))
].copy()
confirmed = confirmed.sort_values("start").reset_index(drop=True)
control = control.sort_values("start").reset_index(drop=True)
log(
    f"  CONFIRMED n={len(confirmed)} {confirmed.seg.value_counts().to_dict()}; "
    f"CONTROL n={len(control)} {control.seg.value_counts().to_dict()}"
)


# ----------------------------------------------------------------------------------
# 5. Per-patch measurement
# ----------------------------------------------------------------------------------
def local_frame(lon, lat):
    """Down-dip unit vector d and along-strike unit vector s (north-positive) in the km frame."""
    h = 0.05
    dzdlon = slabf(lon + h, lat) - slabf(lon - h, lat)
    dzdlat = slabf(lon, lat + h) - slabf(lon, lat - h)
    if not (np.isfinite(dzdlon) and np.isfinite(dzdlat)):
        return None, None
    dzdx = dzdlon / (2 * h * 111.0 * COS47)  # km depth per km east
    dzdy = dzdlat / (2 * h * 111.0)  # km depth per km north
    g = np.array([dzdx, dzdy], float)
    n = np.linalg.norm(g)
    if not np.isfinite(n) or n == 0:
        return None, None
    d = g / n
    s = np.array([-d[1], d[0]])  # rotate +90 deg
    if s[1] < 0:
        s = -s
    return d, s


def measure(pid, start, end, catalog, want_null=False):
    """All S1 quantities for one patch. Returns a dict (values NaN where undefined)."""
    rec = {"id": int(pid), "start": start, "end": end, "catalog": catalog}
    if start not in DAY_INDEX or end not in DAY_INDEX:
        rec["note"] = "date off day-grid"
        return rec
    ts, te = DAY_INDEX[start], DAY_INDEX[end]
    lab = LABEL[catalog][ts : te + 1]
    F = np.any(lab == pid, axis=0)
    nF = int(F.sum())
    rec["n_tri_F"] = nF
    if nF == 0:
        rec["note"] = "empty footprint"
        return rec

    acc = 1e3 * (CUM[catalog][:, te] - CUM[catalog][:, max(ts - PRE_DAYS, 0)])
    acc = np.clip(acc, 0.0, None)  # mm

    w = PAREA * acc
    wF = np.where(F, w, 0.0)
    sw = wF.sum()
    rec["potency_km3"] = float((PAREA[F] * acc[F] * 1e-6).sum())
    rec["max_acc_mm"] = float(acc[F].max())
    if sw <= 0:
        rec["note"] = "zero accumulated slip on footprint"
        return rec

    sx = float((wF * PX).sum() / sw)
    sy = float((wF * PY).sum() / sw)
    rec["slip_x_km"], rec["slip_y_km"] = sx, sy
    slon, slat = to_lonlat(sx, sy)
    rec["slip_lon"], rec["slip_lat"] = float(slon), float(slat)
    rec["slip_depth_km"] = float((wF * PDEPTH).sum() / sw)

    # slip depth edges within the footprint's latitude span +- 0.1 deg
    lat_lo, lat_hi = PLAT[F].min() - 0.1, PLAT[F].max() + 0.1
    band = (PLAT >= lat_lo) & (PLAT <= lat_hi) & (acc >= EDGE_ACC_MM)
    if band.any():
        rec["slip_updip_edge_km"] = float(PDEPTH[band].min())
        rec["slip_downdip_edge_km"] = float(PDEPTH[band].max())
        rec["n_tri_edgeband"] = int(band.sum())
    else:
        rec["slip_updip_edge_km"] = np.nan
        rec["slip_downdip_edge_km"] = np.nan
        rec["n_tri_edgeband"] = 0

    d, s = local_frame(slon, slat)
    if d is None:
        rec["note"] = "no slab gradient at slip centroid"
        return rec
    rec["ddip_x"], rec["ddip_y"] = float(d[0]), float(d[1])

    # spatial tremor mask: within 40 km of any footprint centroid
    idx = TREE.query_ball_point(np.column_stack([PX[F], PY[F]]), r=TREMOR_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(i, int) for i in idx]))
            if len(idx) else np.array([], int))
    rec["n_tremor_near_alltime"] = int(near.size)

    d0 = int((np.datetime64(start) - SHIFT_LO).astype(int))
    d1 = int((np.datetime64(end) - SHIFT_LO).astype(int))

    def q13_from(sel):
        """Q1, Q3, N and the depth percentiles for a boolean selection of `near`."""
        n = int(sel.sum())
        if n == 0:
            return np.nan, np.nan, 0, np.nan, np.nan, np.nan, np.nan, np.nan
        ii = near[sel]
        tx, ty = TR_X[ii].mean(), TR_Y[ii].mean()
        dz = TR_DEPTH[ii]
        dz = dz[np.isfinite(dz)]
        if dz.size == 0:
            med = p10 = p90 = np.nan
        else:
            med, p10, p90 = np.median(dz), np.percentile(dz, 10), np.percentile(dz, 90)
        q1 = (tx - sx) * d[0] + (ty - sy) * d[1]
        q2 = (tx - sx) * s[0] + (ty - sy) * s[1]
        q3 = med - rec["slip_depth_km"]
        return q1, q2, n, q3, med, p10, p90, tx

    if near.size == 0:
        rec["N"] = 0
        return rec
    dnear = TR_DAY[near]
    sel = (dnear >= d0) & (dnear <= d1)
    q1, q2, n, q3, med, p10, p90, tx = q13_from(sel)
    rec["N"] = n
    rec["tremor_x_km"] = tx if n else np.nan
    rec["tremor_y_km"] = float(TR_Y[near[sel]].mean()) if n else np.nan
    rec["tremor_depth_med"] = med
    rec["tremor_depth_p10"] = p10
    rec["tremor_depth_p90"] = p90
    rec["Q1_downdip_km"] = q1
    rec["Q2_strike_km"] = q2
    rec["Q3_depth_km"] = q3
    rec["Q4a_updip_km"] = p10 - rec["slip_updip_edge_km"]
    rec["Q4b_downdip_km"] = p90 - rec["slip_downdip_edge_km"]
    rec["Q5_km3_per_det"] = rec["potency_km3"] / n if n > 0 else np.nan
    rec["log10Q5"] = np.log10(rec["Q5_km3_per_det"]) if n > 0 and rec["potency_km3"] > 0 else np.nan

    if want_null:
        shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
        n1 = np.full(N_SHIFT, np.nan)
        n3 = np.full(N_SHIFT, np.nan)
        nn = np.zeros(N_SHIFT, int)
        for k, sh in enumerate(shifts):
            dsh = (dnear + int(sh)) % SPAN
            selk = (dsh >= d0) & (dsh <= d1)
            a, _, nk, c, *_ = q13_from(selk)
            n1[k], n3[k], nn[k] = a, c, nk
        rec["null_Q1_mean"] = float(np.nanmean(n1)) if np.isfinite(n1).any() else np.nan
        rec["null_Q1_sd"] = float(np.nanstd(n1, ddof=1)) if np.isfinite(n1).sum() > 1 else np.nan
        rec["null_Q3_mean"] = float(np.nanmean(n3)) if np.isfinite(n3).any() else np.nan
        rec["null_Q3_sd"] = float(np.nanstd(n3, ddof=1)) if np.isfinite(n3).sum() > 1 else np.nan
        rec["null_N_med"] = float(np.median(nn))
        rec["null_n_ok"] = int(np.isfinite(n1).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            rec["z1"] = (q1 - rec["null_Q1_mean"]) / rec["null_Q1_sd"]
            rec["z3"] = (q3 - rec["null_Q3_mean"]) / rec["null_Q3_sd"]
    return rec


log("measuring CONFIRMED patches (with 100-shift null each) ...")
rows = []
for _, r in confirmed.iterrows():
    rec = measure(r["id"], r["start"], r["end"], "thrust", want_null=True)
    rec.update({"seg": r["seg"], "Mw": float(r["Mw"]), "dur_d": int(r["dur_d"])})
    rows.append(rec)
    log(
        f"  thrust {int(r['id']):5d} {r['start']} {r['seg']:7s} "
        f"N={rec.get('N', 0):6d} Q1={rec.get('Q1_downdip_km', np.nan):7.1f} "
        f"Q3={rec.get('Q3_depth_km', np.nan):6.1f} z1={rec.get('z1', np.nan):6.2f}"
    )
events = pd.DataFrame(rows)

log("measuring MIRROR control patches ...")
crows = []
for _, r in control.iterrows():
    rec = measure(r["id"], r["start"], r["end"], "mirror", want_null=False)
    rec.update(
        {
            "seg": r["seg"],
            "Mw": float(r["Mw"]),
            "dur_d": int(r["dur_d"]),
            "confirmed": bool(r["confirmed"]),
            "broad_pulse": bool(r["broad_pulse"]),
        }
    )
    crows.append(rec)
controls = pd.DataFrame(crows)
log(f"  mirror rows {len(controls)}, with N>0: {(controls.get('N', pd.Series(dtype=float)) > 0).sum()}")


# ----------------------------------------------------------------------------------
# 6. Trends
# ----------------------------------------------------------------------------------
QUANTS = [
    ("Q1_downdip_km", "Q1 down-dip offset (km)"),
    ("Q3_depth_km", "Q3 depth offset (km)"),
    ("Q4a_updip_km", "Q4a tremor p10 - slip up-dip edge (km)"),
    ("Q4b_downdip_km", "Q4b tremor p90 - slip down-dip edge (km)"),
    ("log10Q5", "log10 Q5 (km^3 per detection)"),
]


def year_float(dates):
    dt = pd.to_datetime(pd.Series(list(dates)))
    return (dt - pd.Timestamp("2010-01-01")).dt.total_seconds().values / (365.25 * 86400.0)


def trend(df, col):
    sub = df[np.isfinite(df[col])] if col in df else df.iloc[0:0]
    n = len(sub)
    out = dict(n=n, slope=np.nan, icept=np.nan, lo=np.nan, hi=np.nan, p_perm=np.nan,
               rho=np.nan, p_spear=np.nan)
    if n < 5:
        return out
    x = year_float(sub["start"])
    y = sub[col].values.astype(float)
    sl, ic, lo, hi = theilslopes(y, x, 0.95)
    out["slope"], out["lo"], out["hi"], out["icept"] = sl, lo, hi, ic
    perm = np.empty(N_PERM)
    for i in range(N_PERM):
        perm[i] = theilslopes(y, RNG.permutation(x), 0.95)[0]
    out["p_perm"] = float((np.abs(perm) >= abs(sl)).mean())
    rho, ps = spearmanr(x, y)
    out["rho"], out["p_spear"] = float(rho), float(ps)
    return out


log("computing Theil-Sen trends + permutation p ...")
trows = []
groups = [("confirmed", "combined", events), ("confirmed", "north", events[events.seg == "north"]),
          ("confirmed", "central", events[events.seg == "central"]),
          ("mirror", "combined", controls), ("mirror", "north", controls[controls.seg == "north"]),
          ("mirror", "central", controls[controls.seg == "central"])]
for setname, seg, df in groups:
    for col, lab in QUANTS:
        t = trend(df, col)
        trows.append(dict(set=setname, seg=seg, quantity=col, label=lab, **t))
trends = pd.DataFrame(trows)


# ----------------------------------------------------------------------------------
# 7. Write outputs
# ----------------------------------------------------------------------------------
ev_cols = [
    "id", "start", "end", "seg", "Mw", "dur_d", "N", "potency_km3", "slip_depth_km",
    "tremor_depth_med", "tremor_depth_p10", "tremor_depth_p90",
    "slip_updip_edge_km", "slip_downdip_edge_km",
    "Q1_downdip_km", "Q2_strike_km", "Q3_depth_km", "Q4a_updip_km", "Q4b_downdip_km",
    "Q5_km3_per_det", "log10Q5", "z1", "z3",
    "n_tri_F", "n_tri_edgeband", "max_acc_mm", "slip_lon", "slip_lat",
    "null_Q1_mean", "null_Q1_sd", "null_Q3_mean", "null_Q3_sd", "null_N_med", "null_n_ok", "note",
]
for c in ev_cols:
    if c not in events:
        events[c] = np.nan
events[ev_cols].to_csv(OUT_RES / "s1_events.csv", index=False)

ct_cols = [c for c in ev_cols if not c.startswith(("z1", "z3", "null_"))] + ["confirmed", "broad_pulse"]
for c in ct_cols:
    if c not in controls:
        controls[c] = np.nan
controls[ct_cols].to_csv(OUT_RES / "s1_controls.csv", index=False)
trends.to_csv(OUT_RES / "s1_trends.csv", index=False)
log(f"wrote results -> {OUT_RES}")


# ----------------------------------------------------------------------------------
# 8. Figure
# ----------------------------------------------------------------------------------
COL = {"north": "#1f5fbf", "central": "#1f9e4a"}
GREY = "#9a9a9a"
ev_dt = pd.to_datetime(events["start"])
ct_dt = pd.to_datetime(controls["start"])

fig, axes = plt.subplots(4, 1, figsize=(11.5, 14), sharex=True)
panels = [
    ("Q1_downdip_km", ["Q1_downdip_km"], "Q1  tremor - slip, down-dip (km)"),
    ("Q3_depth_km", ["Q3_depth_km"], "Q3  tremor median depth - slip depth (km)"),
    ("Q4", ["Q4a_updip_km", "Q4b_downdip_km"], "Q4  tremor p10 - up-dip edge / p90 - down-dip edge (km)"),
    ("log10Q5", ["log10Q5"], "log10 Q5  (km$^3$ potency per detection)"),
]
MK = {"Q4a_updip_km": "^", "Q4b_downdip_km": "v"}

for ax, (key, cols, ylab) in zip(axes, panels):
    titles = []
    for col in cols:
        mk = MK.get(col, "o")
        # mirror control
        m = np.isfinite(controls[col])
        ax.plot(ct_dt[m], controls.loc[m, col], mk, ms=4, mfc="none", color=GREY, alpha=0.55,
                lw=0, label=("mirror " + col.split('_')[0]) if len(cols) > 1 else "mirror")
        for seg in ("north", "central"):
            m = (events.seg == seg) & np.isfinite(events[col])
            ax.plot(ev_dt[m], events.loc[m, col], mk, ms=8, color=COL[seg], lw=0,
                    label=(f"{seg} " + col.split('_')[0]) if len(cols) > 1 else seg)
        t = trends[(trends.set == "confirmed") & (trends.seg == "combined") & (trends.quantity == col)]
        if len(t) and np.isfinite(t.slope.iloc[0]):
            sl, ic = t.slope.iloc[0], t.icept.iloc[0]
            sub = events[np.isfinite(events[col])]
            xx = year_float(sub["start"])
            order = np.argsort(xx)
            ax.plot(pd.to_datetime(sub["start"]).values[order], (ic + sl * xx)[order], "-",
                    color="k" if col == cols[0] else "0.4", lw=1.6,
                    label=f"Theil-Sen {col.split('_')[0]}")
            titles.append(f"{col.split('_')[0]}: {sl:+.3f}/yr (p={t.p_perm.iloc[0]:.3f}, "
                          f"rho={t.rho.iloc[0]:+.2f}, n={int(t.n.iloc[0])})")
    ax.axhline(0, color="0.7", lw=0.8, zorder=0)
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_title("  |  ".join(titles), fontsize=9)
    ax.legend(fontsize=7, ncol=4, loc="best", framealpha=0.85)
    ax.grid(alpha=0.25)
axes[-1].set_xlabel("event start date")
fig.suptitle("S1 tremor - slip co-registration, confirmed Mw>=6 north/central patches "
             "(grey = mirror back-slip control)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(OUT_FIG / "s1_offsets_vs_time.png", dpi=140)
plt.close(fig)
log(f"wrote figure -> {OUT_FIG / 's1_offsets_vs_time.png'}")


# ----------------------------------------------------------------------------------
# 9. Print summary
# ----------------------------------------------------------------------------------
SUMQ = ["Q1_downdip_km", "Q2_strike_km", "Q3_depth_km", "Q4a_updip_km", "Q4b_downdip_km", "log10Q5"]

print("\n" + "=" * 100)
print("S1 CO-REGISTRATION SUMMARY")
print("=" * 100)


def block(df, name):
    print(f"\n--- {name} (n={len(df)}) ---")
    hdr = f"{'group':10s} {'n':>4s} " + " ".join(f"{q.split('_')[0]:>9s}" for q in SUMQ) + \
          f" {'fQ1>0':>7s} {'medN':>7s} {'medPot':>9s}"
    print(hdr)
    for seg in ["combined", "north", "central"]:
        sub = df if seg == "combined" else df[df.seg == seg]
        if len(sub) == 0:
            continue
        vals = " ".join(f"{np.nanmedian(sub[q].values.astype(float)):9.2f}"
                        if q in sub and np.isfinite(sub[q].values.astype(float)).any() else f"{'nan':>9s}"
                        for q in SUMQ)
        q1 = sub["Q1_downdip_km"].values.astype(float)
        q1 = q1[np.isfinite(q1)]
        f1 = np.mean(q1 > 0) if q1.size else np.nan
        print(f"{seg:10s} {len(sub):4d} {vals} {f1:7.2f} "
              f"{np.nanmedian(sub['N'].values.astype(float)):7.0f} "
              f"{np.nanmedian(sub['potency_km3'].values.astype(float)):9.4f}")


block(events, "CONFIRMED (thrust)")
block(controls, "MIRROR control (all)")
block(controls[controls.N > 0] if "N" in controls else controls, "MIRROR control (N>0)")

print("\n--- spread comparison (IQR of the offsets; S1 decision point) ---")
for q in ["Q1_downdip_km", "Q2_strike_km", "Q3_depth_km"]:
    a = events[q].values.astype(float)
    a = a[np.isfinite(a)]
    b = controls[q].values.astype(float)
    b = b[np.isfinite(b)]
    print(f"{q:18s} thrust med {np.median(a):7.2f} IQR {np.subtract(*np.percentile(a, [75, 25])):6.2f} (n={a.size})"
          f"   mirror med {np.median(b):7.2f} IQR {np.subtract(*np.percentile(b, [75, 25])):6.2f} (n={b.size})")

print("\n--- time-shift null z-scores (confirmed; per-event z, then the group-level test) ---")
for zc in ["z1", "z3"]:
    z = events[zc].values.astype(float)
    z = z[np.isfinite(z)]
    tt = ttest_1samp(z, 0.0)
    bt = binomtest(int((z < 0).sum()), z.size, 0.5)
    print(f"{zc}: n={z.size} median {np.median(z):+.2f} mean {z.mean():+.2f} "
          f"(SE {z.std(ddof=1) / np.sqrt(z.size):.2f})  |z|>2 in {(np.abs(z) > 2).sum()}/{z.size}"
          f"   group t p={tt.pvalue:.2g}; sign test {int((z < 0).sum())}/{z.size} negative p={bt.pvalue:.2g}")
print(f"thin nulls (fewer than 50 of {N_SHIFT} shifts gave any tremor): "
      f"{int((events.null_n_ok < 50).sum())} events -> "
      f"{sorted(events.loc[events.null_n_ok < 50, 'id'].astype(int).tolist())}")

print("\n--- TRENDS (Theil-Sen per year, permutation p over 1000 date shuffles) ---")
print(f"{'set':10s} {'seg':9s} {'quantity':18s} {'n':>3s} {'slope/yr':>10s} "
      f"{'95% CI':>20s} {'p_perm':>7s} {'rho':>6s} {'p_sp':>7s}")
for _, r in trends.iterrows():
    if not np.isfinite(r.slope):
        print(f"{r.set:10s} {r.seg:9s} {r.quantity:18s} {int(r.n):3d}   (n too small)")
        continue
    print(f"{r.set:10s} {r.seg:9s} {r.quantity:18s} {int(r.n):3d} {r.slope:10.4f} "
          f"[{r.lo:8.4f},{r.hi:8.4f}] {r.p_perm:7.3f} {r.rho:6.2f} {r.p_spear:7.3f}")

print("\nfiles:")
for p in [OUT_RES / "s1_events.csv", OUT_RES / "s1_controls.csv", OUT_RES / "s1_trends.csv",
          OUT_FIG / "s1_offsets_vs_time.png"]:
    print(f"  {p}  ({p.stat().st_size / 1024:.1f} kB)")
log("done")
