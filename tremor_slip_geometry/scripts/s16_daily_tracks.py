#!/usr/bin/env python
"""S16 PART B helper — save the per-event DAILY tracks (slip rate, tremor count, tremor energy).

Why this script exists
----------------------
S9, S11 and S16 PART A all fit the daily tremor-vs-slip-rate relation, but each of them saved only
the FITTED per-event exponents (results/s9_rate_scaling.csv, s11_energy_events.csv,
s16_fb_events.csv) — the daily series themselves were never written to disk, and the crashed
session's scratch cache was wiped with the pod.  The paper figure set (S16 PART B, PLAN Round 5)
needs the example log-log panels, so the tracks have to be regenerated.

NOTHING IS REDEFINED HERE.  The slip-model loader, the km frame, the footprint / 60-km
neighbourhood geometry, the R7 regressor, the daily measures, the eps/transform conventions and
`ts_slope_fast` are COPIED VERBATIM from scripts/s16_fb_aggregation.py (which copied them from
scripts/s11_scaling_robustness.py, which copied them from scripts/s9_within_event.py).  The ONLY
difference is that the near-tremor selection is NOT restricted to 2017+ here, so the same code
serves the 27-event count fits (S9/S11) and the 18-event amplitude fits (S16 PART A):
  * the energy/magnitude columns are populated only from 2017-01-01, so E_* stays NaN before that;
  * all 18 amplitude-era event windows start after 2017-02-14, so the count series is bit-identical
    to S16's (which pre-filtered the detections to >= 2017-01-01).
The script therefore VALIDATES itself: the per-event Theil-Sen slopes recomputed from the saved
tracks must reproduce the stored exponents of S9 / S11 / S16 exactly.

INPUTS (read-only)
  tremor_slip_geometry/results/s1_events.csv                                    27 confirmed events
  tremor_slip_geometry/results/{s9_rate_scaling,s11_energy_events,s16_fb_events}.csv  validation
  tremor_slip_geometry/results/s11_exponents.csv                                S11's eps for E
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz               label/days/patch_*
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  A_raw (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                                        PNSN tremor
  /home/jovyan/cascadia-sse-strain                                              mesh + sse_strain

DEFINITIONS (verbatim from S16 PART A / S11)
  window     [start - 10 d, end + 10 d] on the slip-patch daily grid
  footprint  triangles labelled with the patch id on >= 1 day of [start, end]
  R7 (mm km^2/d) = sum over the 60-km neighbourhood of area_km2 * max(0, 7-d centred mean of the
                   daily thrust-slip increment), with cum slip = -1e3 * (Ld @ A_raw.T)
  near tremor    PNSN detections within 60 km (km frame) of any footprint triangle centroid
  C          daily detection count (fitted as log10(C + 1), S11's eps)
  E_sum      daily sum of per-detection `energy`      (eps 0)
  E_med      daily MEDIAN per-detection `energy`      (eps 0)  <- the Frank & Brodsky aggregation
  A1         sqrt(E_med)                              (eps 0)
  fit day    R7 > 0 AND C >= 1 (and the measure finite/positive), >= 8 such days for a fit

OUTPUTS (new files only)
  tremor_slip_geometry/results/s16_daily_tracks.csv      one row per (event, day) in the window
  tremor_slip_geometry/results/s16_track_validation.csv  recomputed vs stored per-event exponents

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s16_daily_tracks.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"

PRE_A = 10
POST_A = 10
NB_RADIUS_KM = 60.0
MIN_FIT_DAYS = 8
COS47 = np.cos(np.deg2rad(47.0))
ENERGY_ERA = np.datetime64("2017-01-01")

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def ts_slope_fast(x, y, nmin=MIN_FIT_DAYS):
    """Theil-Sen slope (median of pairwise slopes) — copied from s11_scaling_robustness.py."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < nmin:
        return np.nan
    x, y = x[m], y[m]
    dx = x[:, None] - x[None, :]
    if not np.any(dx != 0):
        return np.nan
    dy = y[:, None] - y[None, :]
    k = dx != 0
    return float(np.median(dy[k] / dx[k]))


def group_quantiles(dd, v, nW, qs):
    """Per-day quantiles of v — copied from s16_fb_aggregation.py."""
    out = {q: np.full(nW, np.nan) for q in qs}
    if dd.size == 0:
        return out
    order = np.lexsort((v, dd))
    dds = dd[order]
    vs = v[order]
    counts = np.bincount(dds, minlength=nW)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    nz = counts > 0
    base = starts[nz]
    cnt = counts[nz]
    for q in qs:
        pos = (cnt - 1) * q
        lo_i = np.floor(pos).astype(np.int64)
        hi_i = np.ceil(pos).astype(np.int64)
        frac = pos - lo_i
        out[q][nz] = vs[base + lo_i] * (1.0 - frac) + vs[base + hi_i] * frac
    return out


# ---------------------------------------------------------------------------- slip model
log("loading slip model (load_solution + invert, ~1 min) ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()

A_raw = np.nan_to_num(
    np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"]
)
CUM_MM = -1e3 * (Ld @ A_raw.T)
DIFF = np.zeros_like(CUM_MM)
DIFF[:, 1:] = np.diff(CUM_MM, axis=1)
RATE7 = np.clip(uniform_filter1d(DIFF, size=7, axis=1, mode="nearest"), 0.0, None)
log(f"  Ld {Ld.shape}, A_raw {A_raw.shape}, RATE7 max {RATE7.max():.3f} mm/d")

# ---------------------------------------------------------------------------- geometry
NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
assert NT == A_raw.shape[0]
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
log(f"  {NT} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles")

# ---------------------------------------------------------------------------- tremor
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                 usecols=["time", "lat", "lon", "magnitude", "num_stas", "duration", "energy"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_ENERGY = tr["energy"].fillna(0.0).values.astype(float)
TR_MAG = tr["magnitude"].values.astype(float)
# pandas trap: to_datetime gives datetime64[us] here -> cast to ns before epoch arithmetic
TR_DAY = tr["t"].values.astype("datetime64[ns]").astype("datetime64[D]")
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(f"  {len(tr)} detections {TR_DAY.min()}..{TR_DAY.max()}")

try:
    _s11x = pd.read_csv(OUT_RES / "s11_exponents.csv")
    EPS_S11_E = float(_s11x.loc[_s11x.measure == "E", "eps"].iloc[0])
except Exception:                                            # pragma: no cover
    EPS_S11_E = 0.0
log(f"  S11 eps for E = {EPS_S11_E:.4f}")

# ---------------------------------------------------------------------------- events
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  {len(ev_in)} confirmed events; {ev_in.seg.value_counts().to_dict()}")

rows = []
val = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    seg = str(r["seg"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1

    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = (np.unique(np.concatenate([np.asarray(h, int) for h in hits])) if len(hits)
          else np.array([], int))
    nb = np.union1d(nb, np.flatnonzero(F))
    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in thits])) if len(thits)
            else np.array([], int))

    aw = PAREA[nb][:, None]
    R7 = (aw * RATE7[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)

    wdays = np.array([np.datetime64(DAYS[i]) for i in range(w0, w1 + 1)])
    dday = (TR_DAY[near] - wdays[0]).astype(int)
    inw = (dday >= 0) & (dday < nW)
    dd = dday[inw].astype(np.int64)
    pos = near[inw]

    C = np.bincount(dd, minlength=nW).astype(float)
    en = TR_ENERGY[pos]
    is17 = TR_DAY[pos] >= ENERGY_ERA
    E_sum = np.full(nW, np.nan)
    E_med = np.full(nW, np.nan)
    C17 = np.bincount(dd[is17], minlength=nW).astype(float) if is17.any() else np.zeros(nW)
    if is17.any():
        E_sum = np.bincount(dd[is17], weights=en[is17], minlength=nW)
        E_sum = np.where(C17 > 0, E_sum, np.nan)
        E_med = group_quantiles(dd[is17], en[is17], nW, (0.50,))[0.50]
    A1 = np.sqrt(E_med)

    for k in range(nW):
        rows.append(dict(id=pid, seg=seg, ev_start=r["start"], ev_end=r["end"],
                         Mw=float(r["Mw"]), date=str(wdays[k]), day_in_window=k - PRE_A,
                         R7=float(R7[k]), C=float(C[k]),
                         E_sum=float(E_sum[k]) if np.isfinite(E_sum[k]) else np.nan,
                         E_med=float(E_med[k]) if np.isfinite(E_med[k]) else np.nan,
                         A1_sqrt_Emed=float(A1[k]) if np.isfinite(A1[k]) else np.nan,
                         in_fit=bool((R7[k] > 0) and (C[k] >= 1))))

    base = (R7 > 0) & (C >= 1)
    lr = np.where(R7 > 0, np.log10(np.where(R7 > 0, R7, 1.0)), np.nan)
    def fit(v, eps):
        y = np.where(np.isfinite(v) & (v > 0), np.log10(np.where(np.isfinite(v) & (v > 0), v, 1.0)
                                                        + eps), np.nan)
        ok = base & np.isfinite(y)
        return ts_slope_fast(lr[ok], y[ok]), int(ok.sum())
    b_C, n_C = fit(C, 1.0)
    b_Esum0, _ = fit(E_sum, 0.0)
    b_EsumS11, _ = fit(E_sum, EPS_S11_E)
    b_Emed, n_E = fit(E_med, 0.0)
    b_A1, _ = fit(A1, 0.0)
    val.append(dict(id=pid, start=r["start"], seg=seg, n_fit_days_C=n_C, n_fit_days_E=n_E,
                    b_C_recomp=b_C, b_E_sum_recomp=b_Esum0, b_E_sum_epsS11_recomp=b_EsumS11,
                    b_E_med_recomp=b_Emed, b_A1_recomp=b_A1))

TRK = pd.DataFrame(rows)
TRK.to_csv(OUT_RES / "s16_daily_tracks.csv", index=False)
log(f"  wrote results/s16_daily_tracks.csv ({len(TRK)} event-days, {TRK.id.nunique()} events)")

# ---------------------------------------------------------------------------- validation
V = pd.DataFrame(val)
s9 = pd.read_csv(OUT_RES / "s9_rate_scaling.csv")
s9 = s9[s9.id > 0][["id", "exponent", "n_days_fit"]].rename(
    columns={"exponent": "b_C_s9", "n_days_fit": "n_days_fit_s9"})
s11 = pd.read_csv(OUT_RES / "s11_energy_events.csv")[
    ["id", "b_C_R7", "b_E_R7", "n_days_fit"]].rename(
    columns={"b_C_R7": "b_C_s11", "b_E_R7": "b_E_epsS11_s11", "n_days_fit": "n_days_fit_s11"})
s16 = pd.read_csv(OUT_RES / "s16_fb_events.csv")[
    ["id", "b_C", "b_E_sum", "b_E_med", "b_A1_sqrt_Emed", "n_days"]].rename(
    columns={"b_C": "b_C_s16", "b_E_sum": "b_E_sum_s16", "b_E_med": "b_E_med_s16",
             "b_A1_sqrt_Emed": "b_A1_s16", "n_days": "n_days_s16"})
V = V.merge(s9, on="id", how="left").merge(s11, on="id", how="left").merge(s16, on="id", how="left")
pairs = [("b_C_recomp", "b_C_s9"), ("b_C_recomp", "b_C_s11"), ("b_C_recomp", "b_C_s16"),
         ("b_E_sum_epsS11_recomp", "b_E_epsS11_s11"), ("b_E_sum_recomp", "b_E_sum_s16"),
         ("b_E_med_recomp", "b_E_med_s16"), ("b_A1_recomp", "b_A1_s16")]
print("\nVALIDATION — recomputed per-event exponents vs stored")
for a, b in pairs:
    d = (V[a] - V[b]).abs()
    n = int(d.notna().sum())
    print(f"  {a:24s} vs {b:18s}  n={n:2d}  max|diff| {np.nanmax(d.values):.3e}")
    V[f"d_{a}_{b}"] = V[a] - V[b]
V.to_csv(OUT_RES / "s16_track_validation.csv", index=False)
log("  wrote results/s16_track_validation.csv")
log("done")
