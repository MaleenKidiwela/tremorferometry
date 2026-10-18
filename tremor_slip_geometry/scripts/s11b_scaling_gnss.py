#!/usr/bin/env python
"""S11b — MODEL-FREE version of the tremor-vs-slip-rate scaling: raw GNSS instead of vbICA slip.

Follow-up to S11 (NOTES.md, 2026-09-10).  S11 established that the sub-linear exponent
b(C vs R7) = 0.40 / b(E vs R7) = 0.64 is not a tremor-side artefact (no count saturation, no
network-sensitivity dependence, energy agrees, the smoothing dependence is errors-in-variables
attenuation of a relation to the SMOOTH rate).  The one blocker left by the S10 literature check
is that the regressor is the vbICA slip rate: a 13-fixed-pattern inversion.  S11b removes the
slip model entirely and regresses the same tremor measures against a slip-rate proxy built
straight from GNSS displacements.

======================================================================================
INPUTS (read-only; everything written lives in tremor_slip_geometry/)
======================================================================================
  /home/jovyan/cascadia-sse-strain/data/2026-08-18/X.mat    (v7.3; read through
      sse_strain.load_solution, which is the h5py reader used by scripts/s12_representability.py)
        Xd/ts        (639 x 7169)  daily displacement, mm, 639 channels = 213 stations x [e,n,u]
        Xd/data_mask (639 x 7169)  1 = real datum, 0 = filled gap (8.5 % of entries)
        Xd/llh       station lon/lat/height (repeated 3x, decimated by sol.station_llh)
        Xd/name      channel names, e.g. ALBHe / ALBHn / ALBHu
        Xd/timeline  (7169,) decimal year, 2007-01-02 .. 2026-08-18 -- the SAME day grid as the
                     slip model, so no re-gridding is needed anywhere
  /home/jovyan/cascadia-sse-strain/data/2026-08-18/ICA.mat  U (639 x 26), S (26), V (7169 x 26)
        and ind_comps.mat -> the 13 retained components (sol.component_index)
  tremor_slip_geometry/results/s1_events.csv                27 confirmed north/central events
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz    patch labels -> the event
        FOOTPRINT only (the events are the same objects S9/S11 used; the slip model enters the
        DEFINITION of the event window and footprint, never the regressor)
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  A_raw, used ONLY
        to rebuild S11's vbICA R7 as the reference column
  catalogs/pnsn_tremor_cascadia_full.csv                    PNSN tremor (energy etc. 2017+ only)

======================================================================================
THE GNSS SLIP-RATE PROXIES  (no slip patterns, no Green's functions, no inversion)
======================================================================================
  Per channel (east and north only; the vertical is dropped):
    1. mask the filled gaps to NaN (RAW source only);
    2. remove, by one least-squares fit over the WHOLE record, a constant, a linear trend, and
       365.25-d and 182.6-d sinusoids (4 seasonal terms: sin/cos annual + sin/cos semiannual);
    3. 7-d centred mean, NaN-aware, requiring >= 4 valid days in the window (15-d and 31-d
       versions are also built: 15-d as a noise ladder, and 31-d because the vbICA amplitudes
       A_raw are TRAILING-30-DAY quantities, so the modeled R7 is effectively a ~30-d rate --
       a 7-d GNSS rate is NOT the same bandwidth and would be compared unfairly.  A CENTRED
       31-d mean is used so the bandwidth matches without importing the trailing window's lag);
    4. daily velocity = first difference of the smoothed displacement, mm/d.
  Per event, stations within 120 km of any footprint centroid (>= 5 required):
    ETS DIRECTION u = unit vector of the station-mean horizontal displacement CHANGE over
      [start, end]  (this comes out west / west-south-west in every event, i.e. the trench-ward
      ETS reversal);  a global fixed direction (the median event direction) is also used as a
      robustness variant so the per-event direction estimate cannot be doing the work.
    PROJECTED RATE  G(d) = ( mean_stations v_e(d), mean_stations v_n(d) ) . u     [mm/d]
      (a day needs >= 5 stations with a finite velocity; the positive part is used for the
      log-log fit, exactly as the modeled rate is clipped at 0 in S9/S11).
    SUM-|E| RATE    G(d) = sum_stations |v_e(d)|   [mm/d], the simple alternative asked for.
  Two data sources:
    RAW  = Xd/ts with the gap mask applied  -> completely independent of the ICA/slip model;
    ICA13 = U[:,k] S[k] V[:,k]^T summed over the 13 RETAINED components -> denoised, still uses
       NO slip patterns and no inversion, but it is NOT independent of the vbICA decomposition
       (it is the same 13 time functions), so it is reported as an intermediate, not as the
       model-free result.
  Proxy columns actually fitted:
    Graw_proj  (raw, 7-d, per-event direction)      Graw_projFix (raw, 7-d, global direction)
    Graw15_proj(raw, 15-d, per-event direction)     Graw31_proj  (raw, 31-d = BANDWIDTH-MATCHED
    Graw_sumE  (raw, 7-d, sum |v_e|)                              to the modeled rate)
    Gica_proj  (ICA13, 7-d, per-event direction)    Gica31_proj  (ICA13, 31-d, matched)
    Gica_sumE  (ICA13, 7-d, sum |v_e|)
    R7_vbICA   (S11's modeled potency rate)  -- the reference, recomputed here identically

======================================================================================
WHAT IS MEASURED  (the S11 analysis, unchanged apart from the regressor)
======================================================================================
  Window [start-10, end+10]; tremor within 60 km of the footprint (S9/S11 geometry).
  Measures: C = daily detection count (all 27 events); E = daily summed energy (18 events from
  2017, when the energy column starts).
  b = Theil-Sen slope of log10(measure + eps) on log10(proxy) over days with measure > 0 and
  proxy > 0, >= 8 days;  eps = 1 for C, the global 1st percentile of positive daily E for E.
  Pooled fit with event dummies = within-event demeaning of both axes, then OLS and Theil-Sen.
  NULL: 200 rigid >= 365-d wrapped shifts of the near-tremor dates with the proxy unchanged
  (the S4/S9/S11 control); the 2017+ measures use their own 2017-01-01.. wrap.
  Day-subsets:  all;  hisnr = days with proxy > 2 x the pre-event noise sd (sd of the proxy over
  [start-70, start-10], 60 d);  R_lower / R_upper = below / above the event's own median
  log10 proxy (the S11 saturation check).
  SNR per event = mean proxy over [start, end] divided by that pre-event 60-d sd.
  ATTENUATION: raw GNSS daily velocity is noisy (see NOTES), so the raw exponent is attenuated by
  errors in the regressor exactly as S11's b(C,R1) was.  Reported per event: the reliability
  ratio lambda = Theil-Sen slope of log10(Gica_proj) on log10(Graw_proj), and of log10(R7_vbICA)
  on log10(Graw_proj), with the implied de-attenuated exponent b / lambda.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s11b_exponents.csv   one row per (measure, proxy, subset): n events, b median/IQR,
                               90 % bootstrap CI of the median, # events with b < 1, pooled
                               OLS/Theil-Sen/r, null mean/sd, per-event z, sign test
  results/s11b_events.csv      one row per event: n stations, ETS direction, SNR, day counts,
                               and b for every proxy x measure, plus the attenuation ratios
  results/s11b_proxy_check.csv one row per event: correlations and Theil-Sen relations between
                               the proxies (Graw vs Gica vs R7_vbICA), lag-1 autocorrelations
  figures/s11b_scaling_gnss.png  6 panels

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s11b_scaling_gnss.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree
from scipy.stats import binomtest, spearmanr, theilslopes, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
SOLDIR = "/home/jovyan/cascadia-sse-strain/data/2026-08-18"
MESH = "/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt"

RNG = np.random.default_rng(20260910)

PRE_A, POST_A = 10, 10          # window [start-10, end+10]                     (S9/S11)
NB_RADIUS_KM = 60.0             # tremor / neighbourhood radius                  (S9/S11)
GNSS_RADIUS_KM = 120.0          # stations within 120 km of the footprint
MIN_STA = 5                     # >= 5 stations for an event and for a day
MIN_FIT_DAYS = 8                # >= 8 usable days for a per-event exponent      (S9/S11)
MIN_HALF_DAYS = 6
N_SHIFT = 200
N_BOOT = 2000
PRE_NOISE_LO, PRE_NOISE_HI = 70, 10     # pre-event noise window [start-70, start-10]
SNR_CUT = 2.0                   # hisnr subset: proxy > 2 x pre-event sd
COS47 = np.cos(np.deg2rad(47.0))
ENERGY_ERA = np.datetime64("2017-01-01")

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def ts_slope_fast(x, y, nmin=MIN_FIT_DAYS):
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


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan, 0
    return (float(np.median(v)), float(np.percentile(v, 25)),
            float(np.percentile(v, 75)), int(v.size))


def fmt_mi(v):
    m, q1, q3, n = med_iqr(v)
    if n == 0:
        return "      n/a          (n= 0)"
    return f"{m:+8.3f} [{q1:+7.3f},{q3:+7.3f}] (n={n:2d})"


def boot_median_ci(v, nboot=N_BOOT):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, 5)), float(np.percentile(bs, 95))


# ==================================================================================
# 1. GNSS: X.mat + ICA.mat through sse_strain (the s12_representability.py reader)
# ==================================================================================
log("loading GNSS solution (X.mat + ICA.mat) ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402

fault = load_fault_ascii(MESH)
sol = load_solution(SOLDIR, fault)

TS = np.asarray(sol.X["ts"], dtype=float)                 # (639, 7169) mm
MASK = np.asarray(sol.X["data_mask"]).astype(bool)        # 1 = real datum
TL = sol.timeline                                         # decimal year, (7169,)
GDATES = np.asarray(sol.dates, dtype="datetime64[D]")
CHAN = [str(x) for x in np.asarray(sol.X["name"], dtype=object).ravel()]
CMP = np.array([c[-1] for c in CHAN])
STA_LLH = sol.station_llh                                 # (213, 3)
NSTA_ALL = STA_LLH.shape[0]
IDX13 = sol.component_index
log(f"  ts {TS.shape} mm, mask real {MASK.mean() * 100:.1f}%, {NSTA_ALL} stations, "
    f"{TL.size} days {GDATES[0]}..{GDATES[-1]}")
log(f"  retained ICA components (0-based) {IDX13}  -> {len(IDX13)} of {sol.S.size}")

REC13 = sol.U[:, IDX13] @ np.diag(sol.S[IDX13]) @ sol.V[:, IDX13].T   # (639, 7169) mm
log(f"  rms: data {np.nanstd(np.where(MASK, TS, np.nan)):.2f} mm, "
    f"13-component reconstruction {REC13.std():.2f} mm "
    f"-> the ETS signal is a small part of the raw daily scatter")


def detrend_deseason(Y):
    """Remove constant + linear trend + 365.25-d and 182.6-d sinusoids, per channel."""
    T = TL.size
    D = np.column_stack([np.ones(T), TL - TL.mean(),
                         np.sin(2 * np.pi * TL), np.cos(2 * np.pi * TL),
                         np.sin(4 * np.pi * TL), np.cos(4 * np.pi * TL)])
    out = np.full_like(Y, np.nan)
    for i in range(Y.shape[0]):
        y = Y[i]
        m = np.isfinite(y)
        if m.sum() < 200:
            continue
        c, *_ = np.linalg.lstsq(D[m], y[m], rcond=None)
        out[i] = y - D @ c
    return out


def nan_smooth(a, w, minv):
    f = np.isfinite(a)
    z = np.where(f, a, 0.0)
    num = uniform_filter1d(z, w, axis=-1, mode="nearest") * w
    den = uniform_filter1d(f.astype(float), w, axis=-1, mode="nearest") * w
    return np.where(den >= minv, num / np.maximum(den, 1e-9), np.nan)


def velocity(disp):
    v = np.full_like(disp, np.nan)
    v[:, 1:] = np.diff(disp, axis=1)
    return v


log("detrending + deseasoning 639 channels (raw and ICA13) ...")
RAW_R = detrend_deseason(np.where(MASK, TS, np.nan))
ICA_R = detrend_deseason(REC13.copy())
D_RAW7 = nan_smooth(RAW_R, 7, 4)
D_RAW15 = nan_smooth(RAW_R, 15, 8)
D_RAW31 = nan_smooth(RAW_R, 31, 16)
D_ICA7 = nan_smooth(ICA_R, 7, 4)
D_ICA31 = nan_smooth(ICA_R, 31, 16)
V_RAW7, V_RAW15, V_RAW31 = velocity(D_RAW7), velocity(D_RAW15), velocity(D_RAW31)
V_ICA7, V_ICA31 = velocity(D_ICA7), velocity(D_ICA31)
log(f"  residual rms after detrend/deseason {np.nanstd(RAW_R):.2f} mm; "
    f"daily velocity rms raw7 {np.nanstd(V_RAW7):.3f}, raw15 {np.nanstd(V_RAW15):.3f}, "
    f"raw31 {np.nanstd(V_RAW31):.3f}, ica7 {np.nanstd(V_ICA7):.3f}, "
    f"ica31 {np.nanstd(V_ICA31):.3f} mm/d")

IE = np.flatnonzero(CMP == "e")
IN = np.flatnonzero(CMP == "n")
assert IE.size == IN.size == NSTA_ALL
SX, SY = to_km(STA_LLH[:, 0], STA_LLH[:, 1])


# ==================================================================================
# 2. vbICA reference rate R7 (rebuilt exactly as in S9/S11)
# ==================================================================================
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz",
                              allow_pickle=True)["A_raw"])
CUM_MM = -1e3 * (Ld @ A_raw.T)
DIFF = np.zeros_like(CUM_MM)
DIFF[:, 1:] = np.diff(CUM_MM, axis=1)
RATE7 = np.clip(uniform_filter1d(DIFF, size=7, axis=1, mode="nearest"), 0.0, None)

NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON, PAREA = NPZ["patch_lat"], NPZ["patch_lon"], NPZ["patch_area_km2"]
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
assert NT == TL.size, "GNSS timeline and slip-model day grid differ"
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
log(f"  vbICA reference rebuilt: {NT} model days, {len(PLAT)} triangles")


# ==================================================================================
# 3. Tremor (loader copied from s9_within_event.py / s11_scaling_robustness.py)
# ==================================================================================
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                 usecols=["time", "lat", "lon", "energy"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_ENERGY = tr["energy"].fillna(0.0).values.astype(float)
SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1
TR_DAY = tr["t"].values.astype("datetime64[D]")
TR_OFF = (TR_DAY - SHIFT_LO).astype(int)
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
SHIFT_LO17 = ENERGY_ERA
SPAN17 = int((SHIFT_HI - SHIFT_LO17).astype(int)) + 1
log(f"  {len(tr)} detections; wrap spans {SPAN} d (all) / {SPAN17} d (2017+)")


def window_subset(sorted_off, s, w0_off, nW, span):
    a = (w0_off - s) % span
    b = (a + nW - 1) % span
    if a <= b:
        idx = np.arange(np.searchsorted(sorted_off, a, "left"),
                        np.searchsorted(sorted_off, b, "right"))
    else:
        idx = np.concatenate([
            np.arange(np.searchsorted(sorted_off, a, "left"), sorted_off.size),
            np.arange(0, np.searchsorted(sorted_off, b, "right"))])
    d = (sorted_off[idx] + s) % span - w0_off
    return idx, d


# ==================================================================================
# 4. Events
# ==================================================================================
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  events {len(ev_in)}  {ev_in.seg.value_counts().to_dict()}")

PROXIES = ["Graw_proj", "Graw_projFix", "Graw15_proj", "Graw31_proj", "Graw_sumE",
           "Gica_proj", "Gica31_proj", "Gica_sumE", "R7_vbICA"]
MEASURES = ["C", "E"]
MEAS_IS17 = {"C": False, "E": True}
SUBSETS = ["all", "hisnr", "R_lower", "R_upper"]

# --- pass 0: the ETS direction of every event, to build the global fixed direction ------------
log("pass 0: per-event station sets and ETS directions ...")
EVSET = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts_i, te_i = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    F = np.any(LABEL[ts_i:te_i + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    dmin = np.hypot(SX[:, None] - fx[None, :], SY[:, None] - fy[None, :]).min(axis=1)
    sel = np.flatnonzero(dmin <= GNSS_RADIUS_KM)
    de = D_RAW7[IE[sel], te_i] - D_RAW7[IE[sel], ts_i]
    dn = D_RAW7[IN[sel], te_i] - D_RAW7[IN[sel], ts_i]
    ue, un = np.nanmean(de), np.nanmean(dn)
    nrm = np.hypot(ue, un)
    u = np.array([ue / nrm, un / nrm]) if nrm > 0 else np.array([np.nan, np.nan])
    EVSET.append(dict(pid=pid, ts=ts_i, te=te_i, sel=sel, u=u, F=F,
                      dmin_med=float(np.median(dmin[sel])) if sel.size else np.nan))
UGLOB = np.array([np.nanmedian([e["u"][0] for e in EVSET]),
                  np.nanmedian([e["u"][1] for e in EVSET])])
UGLOB = UGLOB / np.hypot(*UGLOB)
log(f"  global fixed ETS direction (median over events) = "
    f"({UGLOB[0]:+.3f}, {UGLOB[1]:+.3f})  = azimuth "
    f"{(np.degrees(np.arctan2(UGLOB[0], UGLOB[1])) % 360):.1f} deg "
    f"(west-south-west = the trench-ward ETS reversal)")


def fit_block(mm_dict, prox, noise, names):
    """Every (measure, proxy, subset) exponent for one window.  Module level on purpose: a
    nested def inside the event loop would close over the LAST event's prox/noise."""
    store = {}
    for mname in names:
        v = mm_dict[mname]
        good = np.isfinite(v) & (v > 0)
        eps = 1.0 if mname == "C" else EPS_E
        y_full = np.where(good, np.log10(np.where(good, v, 1.0) + eps), np.nan)
        for pname in PROXIES:
            P = prox[pname]
            okP = np.isfinite(P) & (P > 0)
            ok = good & okP
            key0 = (mname, pname)
            if ok.sum() < MIN_HALF_DAYS:
                for sname in SUBSETS:
                    store[key0 + (sname,)] = np.nan
                continue
            lp = np.where(okP, np.log10(np.where(okP, P, 1.0)), np.nan)
            x, y = lp[ok], y_full[ok]
            store[key0 + ("all",)] = ts_slope_fast(x, y)
            thr = SNR_CUT * noise[pname] if np.isfinite(noise[pname]) else np.nan
            hs = (ok & (P > thr)) if np.isfinite(thr) else np.zeros_like(ok)
            store[key0 + ("hisnr",)] = (ts_slope_fast(lp[hs], y_full[hs])
                                        if hs.sum() >= MIN_FIT_DAYS else np.nan)
            med = np.median(x)
            lo = x <= med
            store[key0 + ("R_lower",)] = (ts_slope_fast(x[lo], y[lo], MIN_HALF_DAYS)
                                          if lo.sum() >= MIN_HALF_DAYS else np.nan)
            store[key0 + ("R_upper",)] = (ts_slope_fast(x[~lo], y[~lo], MIN_HALF_DAYS)
                                          if (~lo).sum() >= MIN_HALF_DAYS else np.nan)
    return store


def proj_rate(V, sel, u):
    ve = np.nanmean(V[IE[sel]], axis=0)
    vn = np.nanmean(V[IN[sel]], axis=0)
    ns = np.isfinite(V[IE[sel]]).sum(axis=0)
    p = ve * u[0] + vn * u[1]
    return np.where(ns >= MIN_STA, p, np.nan), ns


def sum_abs_e(V, sel):
    ve = V[IE[sel]]
    ns = np.isfinite(ve).sum(axis=0)
    return np.where(ns >= MIN_STA, np.nansum(np.abs(ve), axis=0), np.nan)


# ==================================================================================
# 5. Pass 1 — per-event proxies, tremor measures, exponents, nulls
# ==================================================================================
log("pass 1: proxies + tremor measures + exponents + nulls ...")
OBS = {(m, p, s): {} for m in MEASURES for p in PROXIES for s in SUBSETS}
NUL = {(m, p, s): {} for m in MEASURES for p in PROXIES for s in SUBSETS}
EVROWS, POOL_ROWS, CHKROWS = [], [], []
EPS_E = None                                  # filled after pass 1 (needs the pooled E values)
EXAMPLE = {}

for es, (_, r) in zip(EVSET, ev_in.iterrows()):
    pid, ts_i, te_i, sel, u = es["pid"], es["ts"], es["te"], es["sel"], es["u"]
    w0, w1 = max(0, ts_i - PRE_A), min(NT - 1, te_i + POST_A)
    nW = w1 - w0 + 1
    W = np.arange(w0, w1 + 1)
    is17 = pd.Timestamp(r["start"]) >= pd.Timestamp("2017-01-01")

    # ---- GNSS proxies over the whole record, then sliced (the pre-event noise needs the pad) --
    p_raw, ns_raw = proj_rate(V_RAW7, sel, u)
    p_rawfix, _ = proj_rate(V_RAW7, sel, UGLOB)
    p_raw15, _ = proj_rate(V_RAW15, sel, u)
    p_raw31, _ = proj_rate(V_RAW31, sel, u)
    p_ica, _ = proj_rate(V_ICA7, sel, u)
    p_ica31, _ = proj_rate(V_ICA31, sel, u)
    s_raw = sum_abs_e(V_RAW7, sel)
    s_ica = sum_abs_e(V_ICA7, sel)

    F = es["F"]
    fx, fy = PX[F], PY[F]
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = np.unique(np.concatenate([np.asarray(h, int) for h in hits]))
    nb = np.union1d(nb, np.flatnonzero(F))
    r7 = (PAREA[nb][:, None] * RATE7[np.ix_(nb, W)]).sum(axis=0)

    full = {"Graw_proj": p_raw, "Graw_projFix": p_rawfix, "Graw15_proj": p_raw15,
            "Graw31_proj": p_raw31, "Graw_sumE": s_raw, "Gica_proj": p_ica,
            "Gica31_proj": p_ica31, "Gica_sumE": s_ica}
    prox = {k: v[W] for k, v in full.items()}
    prox["R7_vbICA"] = r7

    # pre-event noise (60 d ending 10 d before the start) and SNR, per proxy
    pre = np.arange(max(0, ts_i - PRE_NOISE_LO), ts_i - PRE_NOISE_HI)
    dur = np.arange(ts_i, te_i + 1)
    noise, snr = {}, {}
    for k in PROXIES:
        arr = full[k] if k in full else None
        if arr is None:
            v_pre = r7[:0]
            noise[k] = np.nan
            snr[k] = np.nan
            continue
        v_pre = arr[pre]
        v_pre = v_pre[np.isfinite(v_pre)]
        noise[k] = float(np.std(v_pre)) if v_pre.size >= 20 else np.nan
        v_dur = arr[dur]
        v_dur = v_dur[np.isfinite(v_dur)]
        snr[k] = (float(v_dur.mean() / noise[k])
                  if (v_dur.size and np.isfinite(noise[k]) and noise[k] > 0) else np.nan)
    # the vbICA reference gets its noise from the same rule on its own full series
    r7_full = (PAREA[nb][:, None] * RATE7[nb]).sum(axis=0)
    v_pre = r7_full[pre]
    noise["R7_vbICA"] = float(np.std(v_pre)) if v_pre.size >= 20 else np.nan
    snr["R7_vbICA"] = (float(r7_full[dur].mean() / noise["R7_vbICA"])
                       if noise["R7_vbICA"] and noise["R7_vbICA"] > 0 else np.nan)

    # ---- tremor measures over the window ------------------------------------------------
    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near = np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
    off = TR_OFF[near]
    order = np.argsort(off, kind="stable")
    idx_sorted = near[order]
    soff = off[order]
    enr = TR_ENERGY[idx_sorted]
    m17 = TR_DAY[idx_sorted] >= ENERGY_ERA
    soff17 = (TR_DAY[idx_sorted][m17] - SHIFT_LO17).astype(int)
    enr17 = enr[m17]
    w0_off = w0 - IDX_OF_SHIFT_LO
    w0_off17 = int((np.datetime64(DAYS[w0]) - SHIFT_LO17).astype(int))

    def measures(dday, e_arr, want17):
        C = (np.bincount(dday, minlength=nW).astype(float) if dday.size else np.zeros(nW))
        if want17 and dday.size:
            E = np.bincount(dday, weights=e_arr, minlength=nW)
        else:
            E = np.zeros(nW)
        return {"C": C, "E": E}

    sub, dday = window_subset(soff, 0, w0_off, nW, SPAN)
    meas = measures(dday.astype(np.int64), enr[sub], is17)

    EXAMPLE.setdefault("cand", []).append((float(r["Mw"]), pid, W.copy(),
                                           meas["C"].copy(), prox["Graw_proj"].copy(),
                                           prox["Graw31_proj"].copy(),
                                           prox["Gica_proj"].copy(),
                                           prox["R7_vbICA"].copy(),
                                           str(r["start"])))

    row = dict(id=pid, start=str(r["start"]), seg=str(r["seg"]), Mw=float(r["Mw"]),
               dur_d=int(r["dur_d"]), is17=is17, n_sta=int(sel.size),
               med_sta_dist_km=es["dmin_med"], u_e=u[0], u_n=u[1],
               azim_deg=float(np.degrees(np.arctan2(u[0], u[1])) % 360),
               n_days_W=nW, C_sum=float(meas["C"].sum()))
    for k in PROXIES:
        row[f"noise_{k}"] = noise[k]
        row[f"snr_{k}"] = snr[k]
        row[f"fracpos_{k}"] = float(np.nanmean(prox[k] > 0))
    EVROWS.append(row)

    dfr = dict(id=pid, seg=str(r["seg"]), is17=is17, C=meas["C"], E=meas["E"])
    for k in PROXIES:
        dfr[k] = prox[k]
    POOL_ROWS.append(pd.DataFrame(dfr))

    # ---- proxy inter-comparison ------------------------------------------------------------
    chk = dict(id=pid, start=str(r["start"]), lam_ica_on_raw=np.nan,
               lam_r7_on_raw=np.nan, corr_log_raw_ica=np.nan,
               corr_log_raw_r7=np.nan, corr_log_ica_r7=np.nan,
               lam_r7_on_raw31=np.nan, corr_log_raw31_r7=np.nan,
               ac1_raw31=np.nan)
    gr, gi, rr7 = prox["Graw_proj"], prox["Gica_proj"], prox["R7_vbICA"]
    g31 = prox["Graw31_proj"]
    m31 = np.isfinite(g31) & (g31 > 0) & (rr7 > 0)
    if m31.sum() >= MIN_FIT_DAYS:
        chk["lam_r7_on_raw31"] = ts_slope_fast(np.log10(g31[m31]), np.log10(rr7[m31]))
        chk["corr_log_raw31_r7"] = float(np.corrcoef(np.log10(g31[m31]),
                                                     np.log10(rr7[m31]))[0, 1])
        la = np.log10(g31[m31])
        chk["ac1_raw31"] = (float(np.corrcoef(la[:-1], la[1:])[0, 1]) if la.size > 4 else np.nan)
    mm = np.isfinite(gr) & np.isfinite(gi) & (gr > 0) & (gi > 0) & (rr7 > 0)
    if mm.sum() >= MIN_FIT_DAYS:
        chk["lam_ica_on_raw"] = ts_slope_fast(np.log10(gr[mm]), np.log10(gi[mm]))
        chk["lam_r7_on_raw"] = ts_slope_fast(np.log10(gr[mm]), np.log10(rr7[mm]))
        chk["corr_log_raw_ica"] = float(np.corrcoef(np.log10(gr[mm]), np.log10(gi[mm]))[0, 1])
        chk["corr_log_raw_r7"] = float(np.corrcoef(np.log10(gr[mm]), np.log10(rr7[mm]))[0, 1])
        chk["corr_log_ica_r7"] = float(np.corrcoef(np.log10(gi[mm]), np.log10(rr7[mm]))[0, 1])
    for k, a in (("raw", gr), ("ica", gi), ("r7", rr7)):
        f = np.isfinite(a) & (a > 0)
        la = np.log10(a[f])
        chk[f"ac1_{k}"] = (float(np.corrcoef(la[:-1], la[1:])[0, 1]) if la.size > 4 else np.nan)
    chk["n_days_cmp"] = int(mm.sum())
    CHKROWS.append(chk)

    # ---- exponents + nulls -----------------------------------------------------------------
    es["_ctx"] = dict(prox=prox, noise=noise, meas=meas, is17=is17, nW=nW,
                      soff=soff, soff17=soff17, enr=enr, enr17=enr17,
                      w0_off=w0_off, w0_off17=w0_off17)
    log(f"  {pid:5d} {r['start']} {r['seg']:7s} nsta={sel.size:3d} "
        f"az={row['azim_deg']:5.1f} SNRraw={snr['Graw_proj']:+5.2f} "
        f"SNRica={snr['Gica_proj']:+6.2f} SNRr7={snr['R7_vbICA']:+6.2f}")

POOL = pd.concat(POOL_ROWS, ignore_index=True)
pe = POOL.loc[POOL.is17 & (POOL.E > 0), "E"].values
EPS_E = float(np.percentile(pe, 1))
log(f"  pooled {len(POOL)} event-days; eps for E = {EPS_E:.4g} (1st percentile of positive daily E)")

log("pass 2: exponents and 200-draw shifted-tremor nulls ...")
for es in EVSET:
    ctx = es["_ctx"]
    pid = es["pid"]
    names = ["C", "E"] if ctx["is17"] else ["C"]
    fo = fit_block(ctx["meas"], ctx["prox"], ctx["noise"], names)
    for k, b in fo.items():
        OBS[k][pid] = b
    store = {k: [] for k in fo}
    nW = ctx["nW"]
    for sh in RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT):
        sub_n, dd_n = window_subset(ctx["soff"], int(sh), ctx["w0_off"], nW, SPAN)
        Cn = (np.bincount(dd_n.astype(np.int64), minlength=nW).astype(float)
              if dd_n.size else np.zeros(nW))
        fn = fit_block({"C": Cn, "E": np.zeros(nW)}, ctx["prox"], ctx["noise"], ["C"])
        for k, b in fn.items():
            store[k].append(b)
    if ctx["is17"] and ctx["soff17"].size:
        for sh in RNG.integers(365, max(366, SPAN17 - 365 + 1), size=N_SHIFT):
            sub_n, dd_n = window_subset(ctx["soff17"], int(sh), ctx["w0_off17"], nW, SPAN17)
            dd = dd_n.astype(np.int64)
            En = (np.bincount(dd, weights=ctx["enr17"][sub_n], minlength=nW)
                  if dd.size else np.zeros(nW))
            fn = fit_block({"C": np.zeros(nW), "E": En}, ctx["prox"], ctx["noise"], ["E"])
            for k, b in fn.items():
                store[k].append(b)
    for k, lst in store.items():
        if lst:
            NUL[k][pid] = np.asarray(lst, float)

EVD = pd.DataFrame(EVROWS)
CHK = pd.DataFrame(CHKROWS)
for m in MEASURES:
    for p in PROXIES:
        EVD[f"b_{m}_{p}"] = [OBS[(m, p, "all")].get(i, np.nan) for i in EVD.id]
EVD["lam_ica_on_raw"] = CHK.lam_ica_on_raw.values
EVD["lam_r7_on_raw"] = CHK.lam_r7_on_raw.values
EVD["b_C_Graw_deatt_ica"] = EVD.b_C_Graw_proj / EVD.lam_ica_on_raw
EVD["b_C_Graw_deatt_r7"] = EVD.b_C_Graw_proj / EVD.lam_r7_on_raw
EVD["b_E_Graw_deatt_r7"] = EVD.b_E_Graw_proj / EVD.lam_r7_on_raw
EVD["lam_r7_on_raw31"] = CHK.lam_r7_on_raw31.values
EVD["b_C_Graw31_deatt_r7"] = EVD.b_C_Graw31_proj / EVD.lam_r7_on_raw31
EVD["b_E_Graw31_deatt_r7"] = EVD.b_E_Graw31_proj / EVD.lam_r7_on_raw31
EVD.to_csv(OUT_RES / "s11b_events.csv", index=False)
CHK.to_csv(OUT_RES / "s11b_proxy_check.csv", index=False)


# ==================================================================================
# 6. Exponent table
# ==================================================================================
def pooled_fit(df, mname, pname):
    v = df[mname].values
    P = df[pname].values
    ok = np.isfinite(v) & (v > 0) & np.isfinite(P) & (P > 0)
    if ok.sum() < 20:
        return np.nan, np.nan, np.nan, int(ok.sum())
    d = df.loc[ok]
    eps = 1.0 if mname == "C" else EPS_E
    dfx = pd.DataFrame(dict(i=d["id"].values, x=np.log10(d[pname].values),
                            y=np.log10(d[mname].values + eps)))
    dfx["xc"] = dfx.x - dfx.groupby("i").x.transform("mean")
    dfx["yc"] = dfx.y - dfx.groupby("i").y.transform("mean")
    if np.ptp(dfx.xc.values) == 0:
        return np.nan, np.nan, np.nan, len(dfx)
    return (float(np.polyfit(dfx.xc.values, dfx.yc.values, 1)[0]),
            ts_slope_fast(dfx.xc.values, dfx.yc.values, 20),
            float(np.corrcoef(dfx.xc.values, dfx.yc.values)[0, 1]), len(dfx))


rows = []
for mname in MEASURES:
    dfm = POOL[POOL.is17] if MEAS_IS17[mname] else POOL
    for pname in PROXIES:
        ols, tss, rr, npool = pooled_fit(dfm, mname, pname)
        for sname in SUBSETS:
            key = (mname, pname, sname)
            bs = np.array([b for b in OBS[key].values() if np.isfinite(b)])
            m, q1, q3, n = med_iqr(bs)
            blo, bhi = boot_median_ci(bs)
            zs, ps = [], []
            for pid in OBS[key]:
                if not np.isfinite(OBS[key][pid]) or pid not in NUL[key]:
                    continue
                d = NUL[key][pid]
                d = d[np.isfinite(d)]
                if d.size < 20:
                    continue
                ps.append(float((d >= OBS[key][pid]).mean()))
                zs.append((OBS[key][pid] - d.mean()) / d.std() if d.std() > 0 else np.nan)
            zs = np.asarray(zs, float)
            zf = zs[np.isfinite(zs)]
            npos = int((zf > 0).sum())
            nl = [NUL[key][p] for p in NUL[key]]
            nl = np.concatenate(nl) if nl else np.array([])
            nl = nl[np.isfinite(nl)]
            rows.append(dict(measure=mname, proxy=pname, subset=sname, n_events=n,
                             b_median=m, b_q25=q1, b_q75=q3, b_boot_lo90=blo, b_boot_hi90=bhi,
                             n_below_1=int((bs < 1).sum()) if n else 0,
                             pooled_ols=ols, pooled_theilsen=tss, pooled_r=rr,
                             n_pooled_days=npool,
                             null_mean=float(nl.mean()) if nl.size else np.nan,
                             null_sd=float(nl.std()) if nl.size else np.nan,
                             z_median=float(np.median(zf)) if zf.size else np.nan,
                             n_events_z=int(zf.size), n_z_positive=npos,
                             sign_p=binomtest(npos, zf.size).pvalue if zf.size else np.nan,
                             n_events_p05=int(np.sum(np.asarray(ps) < 0.05)) if ps else 0))
EXP = pd.DataFrame(rows)
EXP.to_csv(OUT_RES / "s11b_exponents.csv", index=False)


# ==================================================================================
# 7. Print
# ==================================================================================
print("\n" + "=" * 112)
print("S11b — MODEL-FREE SCALING: tremor vs a GNSS slip-rate proxy (no slip patterns)")
print("=" * 112)
print(f"  {NSTA_ALL} GNSS stations, 639 channels, {TL.size} daily epochs "
      f"{GDATES[0]}..{GDATES[-1]} (same grid as the slip model).")
print(f"  detrended + deseasoned (365.25 d and 182.6 d), 7-d centred mean, differenced; "
      f"stations within {GNSS_RADIUS_KM:.0f} km of the footprint (>= {MIN_STA}).")
print(f"  data rms {np.nanstd(np.where(MASK, TS, np.nan)):.2f} mm, residual after detrend/"
      f"deseason {np.nanstd(RAW_R):.2f} mm, 13-component reconstruction {REC13.std():.2f} mm.")
print(f"  global ETS direction ({UGLOB[0]:+.3f},{UGLOB[1]:+.3f}) = azimuth "
      f"{(np.degrees(np.arctan2(UGLOB[0], UGLOB[1])) % 360):.1f} deg; per-event azimuths "
      f"{EVD.azim_deg.min():.0f}-{EVD.azim_deg.max():.0f} deg (all west-ish = trench-ward).")
print(f"  stations per event {fmt_mi(EVD.n_sta)}")

print("\n  PER-EVENT SIGNAL-TO-NOISE of the proxies "
      "(mean rate over DUR / sd of the same proxy over [start-70, start-10]):")
for k in PROXIES:
    print(f"    {k:14s} SNR {fmt_mi(EVD[f'snr_{k}'])}   noise sd "
          f"{np.nanmedian(EVD[f'noise_{k}']):.4g}   days with rate > 0 "
          f"{np.nanmedian(EVD[f'fracpos_{k}']):.2f}")

print("\n" + "=" * 112)
print("EXPONENTS  b = Theil-Sen slope of log10(measure + eps) on log10(proxy)")
print("=" * 112)
for sname in SUBSETS:
    print(f"\n  --- subset = {sname} " + "-" * 84)
    print(f"    {'measure':7s} {'proxy':14s} {'n':>3s}  {'b median [IQR]':>28s}  "
          f"{'boot90 CI':>18s}  {'pooled OLS/TS':>16s}  {'null mean+-sd':>15s} "
          f"{'sign':>7s} {'p_sign':>9s} {'<1':>4s}")
    for mname in MEASURES:
        for pname in PROXIES:
            x = EXP[(EXP.measure == mname) & (EXP.proxy == pname) & (EXP.subset == sname)]
            if x.empty:
                continue
            x = x.iloc[0]
            print(f"    {mname:7s} {pname:14s} {int(x.n_events):3d}  "
                  f"{x.b_median:+8.3f} [{x.b_q25:+7.3f},{x.b_q75:+7.3f}]  "
                  f"[{x.b_boot_lo90:+7.3f},{x.b_boot_hi90:+7.3f}]  "
                  f"{x.pooled_ols:+7.3f}/{x.pooled_theilsen:+7.3f}  "
                  f"{x.null_mean:+7.3f}+-{x.null_sd:5.3f} "
                  f"{int(x.n_z_positive):3d}/{int(x.n_events_z):3d} {x.sign_p:9.2e} "
                  f"{int(x.n_below_1):3d}")

print("\n" + "=" * 112)
print("PROXY INTER-COMPARISON  (are the GNSS proxies measuring the same rate as the model?)")
print("=" * 112)
print(f"  corr(log Graw, log Gica) {fmt_mi(CHK.corr_log_raw_ica)}")
print(f"  corr(log Graw, log R7)   {fmt_mi(CHK.corr_log_raw_r7)}")
print(f"  corr(log Gica, log R7)   {fmt_mi(CHK.corr_log_ica_r7)}")
print(f"  corr(log Graw31, log R7) {fmt_mi(CHK.corr_log_raw31_r7)}   <- the BANDWIDTH-MATCHED "
      f"proxy")
print(f"  lag-1 autocorrelation of log10 proxy: raw7 {fmt_mi(CHK.ac1_raw)}   "
      f"raw31 {fmt_mi(CHK.ac1_raw31)}   ica7 {fmt_mi(CHK.ac1_ica)}   R7 {fmt_mi(CHK.ac1_r7)}")
print("  (the vbICA A_raw is a TRAILING-30-DAY amplitude, so R7 is effectively a ~30-d rate: only "
      "the 31-d\n   GNSS proxy is the same bandwidth, and its autocorrelation is the one that "
      "matches.)")
print("\n  ATTENUATION (errors in the regressor; the S11 identity b_noisy = b_clean * lambda):")
print(f"    lambda = Theil-Sen slope of log Gica on log Graw  {fmt_mi(CHK.lam_ica_on_raw)}")
print(f"    lambda = Theil-Sen slope of log R7   on log Graw  {fmt_mi(CHK.lam_r7_on_raw)}")
print(f"    de-attenuated b(C vs Graw) using the ICA lambda   {fmt_mi(EVD.b_C_Graw_deatt_ica)}")
print(f"    de-attenuated b(C vs Graw) using the R7  lambda   {fmt_mi(EVD.b_C_Graw_deatt_r7)}")
print(f"    de-attenuated b(E vs Graw) using the R7  lambda   {fmt_mi(EVD.b_E_Graw_deatt_r7)}")
print(f"    lambda = Theil-Sen slope of log R7 on log Graw31   {fmt_mi(CHK.lam_r7_on_raw31)}")
print(f"    de-attenuated b(C vs Graw31) using that lambda     {fmt_mi(EVD.b_C_Graw31_deatt_r7)}")
print(f"    de-attenuated b(E vs Graw31) using that lambda     {fmt_mi(EVD.b_E_Graw31_deatt_r7)}")

print("\n" + "=" * 112)
print("AGREEMENT WITH THE vbICA-BASED EXPONENT")
print("=" * 112)
for mname, ref in (("C", 0.401), ("E", 0.644)):
    a = EVD[f"b_{mname}_Graw_proj"]
    a31 = EVD[f"b_{mname}_Graw31_proj"]
    c = EVD[f"b_{mname}_Gica_proj"]
    c31 = EVD[f"b_{mname}_Gica31_proj"]
    d = EVD[f"b_{mname}_R7_vbICA"]
    print(f"\n  measure {mname}:  raw GNSS  7-d {fmt_mi(a)}")
    print(f"              raw GNSS 31-d {fmt_mi(a31)}   <- BANDWIDTH-MATCHED, model-free")
    print(f"              ICA13     7-d {fmt_mi(c)}")
    print(f"              ICA13    31-d {fmt_mi(c31)}")
    print(f"              vbICA R7      {fmt_mi(d)}   (S11 reported {ref:+.3f})")
    for lab, x, y in ((f"b({mname},Graw31) vs b({mname},R7)", a31, d),
                      (f"b({mname},Gica31) vs b({mname},R7)", c31, d),
                      (f"b({mname},Gica) vs b({mname},R7)", c, d),
                      (f"b({mname},Graw) vs b({mname},R7)", a, d)):
        mm = np.isfinite(x) & np.isfinite(y)
        if mm.sum() >= 6:
            rho, pv = spearmanr(x[mm], y[mm])
            dd = (x - y)[mm]
            print(f"    {lab}: Spearman rho = {rho:+.3f} p = {pv:.4f}; "
                  f"paired difference {fmt_mi(dd)} Wilcoxon p = {wilcoxon(dd).pvalue:.4f}")

print("\n  SATURATION CHECK on the GNSS proxies (upper minus lower half of log10 proxy):")
for mname in MEASURES:
    for pname in ("Graw_proj", "Graw31_proj", "Gica31_proj", "R7_vbICA"):
        lo = np.array([OBS[(mname, pname, "R_lower")].get(i, np.nan) for i in EVD.id])
        hi = np.array([OBS[(mname, pname, "R_upper")].get(i, np.nan) for i in EVD.id])
        mm = np.isfinite(lo) & np.isfinite(hi)
        if mm.sum() >= 6:
            dd = (hi - lo)[mm]
            print(f"    {mname} vs {pname:12s} lower {np.nanmedian(lo[mm]):+.3f} upper "
                  f"{np.nanmedian(hi[mm]):+.3f}  diff {fmt_mi(dd)} "
                  f"Wilcoxon p = {wilcoxon(dd).pvalue:.4f}")


# ==================================================================================
# 8. FIGURE
# ==================================================================================
log("figure ...")
fig = plt.figure(figsize=(19, 11.5))
gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.27,
                      left=0.055, right=0.985, top=0.92, bottom=0.07)

# A: example event
ax = fig.add_subplot(gs[0, 0])
cand = sorted(EXAMPLE["cand"], key=lambda z: -z[0])[0]
_, pid_x, Wx, Cx, GRx, G31x, GIx, R7x, startx = cand
dd = np.arange(len(Wx))
ax.plot(dd, Cx / np.nanmax(Cx), "k-", lw=2.2, label="tremor count (scaled)")
ax.plot(dd, GRx / np.nanmax(np.abs(GRx)), "-", color="tab:green", lw=1.1, alpha=0.7,
        label="raw GNSS 7-d rate (scaled)")
ax.plot(dd, G31x / np.nanmax(np.abs(G31x)), "-", color="tab:olive", lw=2.2,
        label="raw GNSS 31-d rate (scaled)")
ax.plot(dd, GIx / np.nanmax(np.abs(GIx)), "-", color="tab:blue", lw=2.0,
        label="ICA13 GNSS rate (scaled)")
ax.plot(dd, R7x / np.nanmax(R7x), "-", color="tab:purple", lw=2.0, label="vbICA R7 (scaled)")
ax.axhline(0, color="grey", lw=0.7)
ax.set_xlabel("day within the window [start-10, end+10]")
ax.set_ylabel("scaled to its own maximum")
ax.set_title(f"A  event {pid_x} ({startx}), the largest Mw\nraw GNSS is the same signal buried "
             f"in daily noise", fontsize=10)
ax.legend(fontsize=7.5)

# B: exponent by proxy
ax = fig.add_subplot(gs[0, 1])
order = ["Graw_proj", "Graw_projFix", "Graw15_proj", "Graw31_proj", "Graw_sumE",
         "Gica_proj", "Gica31_proj", "Gica_sumE", "R7_vbICA"]
pos = np.arange(len(order))
dataC = [np.array([OBS[("C", p, "all")].get(i, np.nan) for i in EVD.id]) for p in order]
dataC = [d[np.isfinite(d)] for d in dataC]
bp = ax.boxplot(dataC, positions=pos - 0.17, widths=0.3, patch_artist=True,
                showfliers=False, medianprops=dict(color="k", lw=1.6))
for b in bp["boxes"]:
    b.set_facecolor("tab:blue")
    b.set_alpha(0.55)
dataE = [np.array([OBS[("E", p, "all")].get(i, np.nan) for i in EVD.id]) for p in order]
dataE = [d[np.isfinite(d)] for d in dataE]
bp2 = ax.boxplot(dataE, positions=pos + 0.17, widths=0.3, patch_artist=True,
                 showfliers=False, medianprops=dict(color="k", lw=1.6))
for b in bp2["boxes"]:
    b.set_facecolor("tab:red")
    b.set_alpha(0.55)
ax.axhline(1.0, color="k", ls="--", lw=1.2)
ax.axhline(0.0, color="grey", lw=0.7)
ax.set_xticks(pos)
ax.set_xticklabels(order, rotation=25, ha="right", fontsize=8)
ax.set_ylabel("exponent b over events")
ax.set_ylim(-1.0, 2.0)
ax.set_title("B  exponent by slip-rate proxy\nblue = count, red = energy; dashed = linear (b=1)",
             fontsize=10)

# C: per-event agreement
ax = fig.add_subplot(gs[0, 2])
for lab, col, xk in (("raw GNSS 7-d", "tab:green", "Graw_proj"),
                     ("raw GNSS 31-d (matched)", "tab:olive", "Graw31_proj"),
                     ("ICA13 31-d", "tab:blue", "Gica31_proj")):
    x = EVD[f"b_C_{xk}"].values
    y = EVD["b_C_R7_vbICA"].values
    mm = np.isfinite(x) & np.isfinite(y)
    rho = spearmanr(x[mm], y[mm])[0] if mm.sum() >= 6 else np.nan
    ax.scatter(x[mm], y[mm], s=45, c=col, edgecolor="k", lw=0.4,
               label=f"{lab}  rho = {rho:+.2f}")
lim = [-0.6, 1.3]
ax.plot(lim, lim, "k-", lw=1)
ax.axhline(0, color="grey", lw=0.6)
ax.axvline(0, color="grey", lw=0.6)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("b(C vs GNSS proxy)")
ax.set_ylabel("b(C vs vbICA R7)")
ax.set_title("C  per-event agreement with the model-based exponent", fontsize=10)
ax.legend(fontsize=8)

# D: pooled scatter
ax = fig.add_subplot(gs[1, 0])
for lab, col, xk in (("raw GNSS 7-d", "tab:green", "Graw_proj"),
                     ("raw GNSS 31-d (matched)", "tab:olive", "Graw31_proj")):
    d = POOL[(POOL.C > 0) & (POOL[xk] > 0)].copy()
    d["x"] = np.log10(d[xk])
    d["y"] = np.log10(d.C + 1.0)
    d["xc"] = d.x - d.groupby("id").x.transform("mean")
    d["yc"] = d.y - d.groupby("id").y.transform("mean")
    ax.scatter(d.xc, d.yc, s=7, c=col, alpha=0.32, lw=0)
    b = EXP[(EXP.measure == "C") & (EXP.proxy == xk) & (EXP.subset == "all")].iloc[0]
    xx = np.linspace(d.xc.quantile(0.01), d.xc.quantile(0.99), 20)
    ax.plot(xx, b.pooled_theilsen * xx, "-", color=col, lw=2.5,
            label=f"{lab}: pooled b = {b.pooled_theilsen:+.2f}")
xx = np.linspace(-1.5, 1.5, 5)
ax.plot(xx, xx, "k--", lw=1, label="b = 1 (linear)")
ax.set_xlim(-2.0, 1.6)
ax.set_xlabel("log10 GNSS rate  (within-event demeaned)")
ax.set_ylabel("log10 (count + 1)  (within-event demeaned)")
ax.set_title("D  pooled fit with event dummies, model-free regressor", fontsize=10)
ax.legend(fontsize=8)

# E: SNR vs attenuation
ax = fig.add_subplot(gs[1, 1])
ax.scatter(EVD.snr_Graw_proj, EVD.b_C_Graw_proj, s=45, c="tab:green", edgecolor="k", lw=0.4,
           label="raw GNSS")
ax.scatter(EVD.snr_Graw31_proj, EVD.b_C_Graw31_proj, s=45, c="tab:olive",
           edgecolor="k", lw=0.4, label="raw GNSS 31-d (matched)")
ax.scatter(EVD.snr_Gica_proj, EVD.b_C_Gica_proj, s=45, c="tab:blue", edgecolor="k",
           lw=0.4, label="ICA13 GNSS 7-d")
ax.scatter(EVD.snr_R7_vbICA, EVD.b_C_R7_vbICA, s=45, c="tab:purple", edgecolor="k", lw=0.4,
           label="vbICA R7")
ax.set_xscale("log")
ax.axhline(0, color="grey", lw=0.7)
ax.axhline(1.0, color="k", ls="--", lw=1)
ax.set_xlabel("per-event proxy SNR  (mean rate in DUR / pre-event daily sd)")
ax.set_ylabel("b(C vs proxy)")
ax.set_title("E  the exponent tracks the proxy's SNR\n"
             "(errors in the regressor attenuate b toward 0)", fontsize=10)
ax.legend(fontsize=8)

# F: nulls
ax = fig.add_subplot(gs[1, 2])
nl = np.concatenate([NUL[("C", "Graw31_proj", "all")][p]
                     for p in NUL[("C", "Graw31_proj", "all")]])
nl = nl[np.isfinite(nl)]
ax.hist(nl, bins=60, density=True, color="0.78", label=f"shifted-tremor null (n={nl.size})")
for lab, col, xk in (("b(C, raw 7-d)", "tab:green", "b_C_Graw_proj"),
                     ("b(C, raw 31-d)", "tab:olive", "b_C_Graw31_proj"),
                     ("b(E, raw 31-d)", "tab:red", "b_E_Graw31_proj")):
    v = EVD[xk].dropna().values
    ax.hist(v, bins=10, density=True, histtype="step", lw=2.2, color=col,
            label=f"{lab} (n={v.size})")
ax.axvline(1.0, color="k", ls="--", lw=1.2, label="b = 1")
ax.set_xlim(-1.2, 1.8)
ax.set_xlabel("exponent b")
ax.set_ylabel("density")
ax.set_title("F  observed exponents vs the shifted-tremor null", fontsize=10)
ax.legend(fontsize=7)

fig.suptitle("S11b — model-free tremor-vs-slip-rate scaling: PNSN tremor against a GNSS "
             "slip-rate proxy (27 confirmed events, 2010-2026)", fontsize=13)
fig.savefig(OUT_FIG / "s11b_scaling_gnss.png", dpi=140)
plt.close(fig)

log("done.  wrote results/s11b_exponents.csv, s11b_events.csv, s11b_proxy_check.csv and "
    "figures/s11b_scaling_gnss.png")
