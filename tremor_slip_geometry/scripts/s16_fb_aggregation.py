#!/usr/bin/env python
"""S16 PART A — match the Frank & Brodsky (2019) AGGREGATION so the exponents can be compared.

Why
---
Frank & Brodsky (2019, Sci. Adv. 5, eaaw9386, Guerrero) calibrate

        Mdot_0^geo = alpha * A^beta ,    beta = 3.17 +- 0.29 ,   alpha = 10^(13.87 +- 0.19)

where **A is the MEDIAN displacement amplitude of the individual LFEs in a time bin**.  Inverted,
that is  A ~ Mdot^(1/3.17) = Mdot^0.315 — a sub-linear exponent for a PER-EVENT (per-LFE)
amplitude.

Our S11 exponent is for a DIFFERENT aggregation: the daily SUM of per-detection radiated energy,
b(E_sum vs R7) = +0.644.  S10 (results/s10_literature_check.md, candidate 1) flagged that two
readings of that number are possible and that they disagree:
  * summed-amplitude reading  : E_sum ~ R^0.64, energy ~ amplitude^2  ->  A_sum ~ R^0.32 == 1/3.17
  * per-detection reading     : b(E/C) = 0.64 - 0.40 = 0.23  ->  A_per ~ R^0.12, a factor ~2.6 off
and told us to decide which is like-for-like BEFORE the number goes anywhere near an abstract.
This script does that: it computes the F&B statistic itself — the **median per-detection amplitude
per event-day** — on our data, with the same fitting machinery, the same shifted-tremor null and
the same 18 events that carry amplitudes.

======================================================================================
INPUTS (read-only)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv                        27 confirmed events; the 18
                                                                    with start >= 2017-01-01 are used
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz   label, days, patch_lat/lon/area
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                            PNSN tremor; energy, magnitude,
                                                                    num_stas, duration populated
                                                                    ONLY from 2017-01-01
  /home/jovyan/cascadia-sse-strain                                  mesh + sse_strain loaders

  The slip-model / geometry / tremor / km-frame / event loaders and `ts_slope_fast` are COPIED
  VERBATIM from scripts/s11_scaling_robustness.py (which copied them from s9_within_event.py).
  Cross-check printed at the end: b(C, R7) on these 18 events must reproduce
  results/s11_energy_events.csv, and b(E_sum, R7) must reproduce b_E_R7 there.

======================================================================================
THE CATALOGUE'S OWN ENERGY-MAGNITUDE RELATION (measured here, not assumed)
======================================================================================
  Over all 482,659 detections from 2017-01-01 on:
        log10(energy) = 1.781 * magnitude + 2.904 + 0.007*log10(duration),  resid sd 0.052 dex
  and within one 0.1-magnitude bin the energies span exactly 10^(0.1*1.781) = 1.51, i.e. `energy`
  is a DETERMINISTIC monotone function of an unrounded magnitude and `magnitude` is that same
  quantity rounded to 0.1.  So energy and magnitude are ONE observable here, not two, and the
  amplitude proxy one builds from them depends entirely on the magnitude convention assumed:

    A1  A = sqrt(median energy)               assumes  E ~ A^2  (energy = squared envelope over a
                                              near-constant duration: `duration` takes only the two
                                              values 150 s and 300 s).  This is the LIKE-FOR-LIKE
                                              proxy for F&B, whose A is a measured displacement
                                              amplitude.   -> log10 A = 0.891 * M + const
    A2  A = 10^(0.5 * median magnitude)       the form the task specifies.  ASSUMPTION: that
                                              log10 E = 1.0 * M, i.e. that PNSN `magnitude` is a
                                              log-ENERGY magnitude of unit slope.  The catalogue
                                              says the slope is 1.781, so A2 COMPRESSES the true
                                              amplitude range by 0.5/0.891 = 0.561 and its exponent
                                              is 0.561 x the A1 exponent BY CONSTRUCTION.
    A3  A = 10^(0.891 * median magnitude)     A1 rewritten through the fitted E-M relation; agrees
                                              with A1 to the magnitude rounding.
    A4  A = 10^(1.0 * median magnitude)       if MeL is read as a LOCAL (amplitude) magnitude,
                                              M = log10 A + const.  Also consistent with E ~ A^2
                                              at the fitted slope 1.781 ~ 2.
  All four are reported.  A1 is the primary; A2 is reported because it was asked for, with the
  factor 0.561 stated so it cannot be mistaken for an independent measurement.
  Two AGGREGATED amplitude measures are added for contrast, because S10's "summed-amplitude
  reading" of the 0.64 energy exponent needs testing on its own terms:
    Asum_sqrt_Esum  = sqrt(daily SUM of energy)   (the reading S10 floated: 0.64 / 2 = 0.32)
    Asum_lin        = daily SUM of sqrt(energy)   (the LITERAL summed amplitude; its exponent must
                                                   be b_C + b_A_per, not b_Esum / 2)

======================================================================================
REGRESSION DIRECTION (this is where the F&B comparison actually turns)
======================================================================================
  F&B fit log Mdot_geo = log alpha + beta log A, i.e. the MOMENT RATE is the dependent variable.
  A regression is not invertible unless r = 1, and our per-event r^2 of log A on log R is ~0.11,
  so the two directions differ by ~9x.  All three estimators are therefore reported:
      forward  b_fwd  = slope of log A on log R                       (our S11 convention)
      reverse  beta   = slope of log R on log A                       (F&B's convention)
      RMA      b_rma  = sign(r) * sd(log A)/sd(log R) = b_fwd_ols/|r| (direction-symmetric,
                        errors in both variables; beta_rma = 1/b_rma)
  1/b_fwd is also carried in the CSV as `naive_inv_fwd_median` so nobody re-derives it wrongly.

======================================================================================
TIME-BIN LADDER
======================================================================================
  F&B's A is a median "over a time bin" whose length is not ours.  A shorter bin means a noisier
  median (fewer detections) and a noisier rate, which ATTENUATES the forward slope, so the bin
  length is a candidate explanation for any gap.  The amplitude exponents are therefore recomputed
  on 1-, 3- and 7-day bins (median per-detection energy over all detections in the bin; R7 averaged
  over the bin).  Nulls are run for the 1-day bin only.

======================================================================================
WHAT IS COMPUTED  (per event, per event-day, 18 events from 2017)
======================================================================================
  Daily tremor measures over the near-footprint detections (60 km of the footprint, S9/S11 geometry):
    C          detection count                                     (S11's measure, eps = 1)
    E_sum      summed per-detection energy                          (S11's E)
    E_med      MEDIAN per-detection energy                          <- the F&B aggregation
    E_mean     E_sum / C                                            (S11's EperDet, for reference)
    E_p10/25/75/90   the same day's energy percentiles              (which part of the distribution
                                                                     responds to slip rate?)
    M_med      median magnitude ;  M_max  daily maximum magnitude
    NS_med     median num_stas                                      (detection-threshold diagnostic)
    E_med_ns8  median energy among num_stas >= 8 detections; C_ns8   (threshold-controlled version)
    A1..A4     the amplitude proxies above

  Exponents: b = Theil-Sen slope of log10(measure) on log10(R7) over event-days with the measure
  positive, R7 > 0 and C >= 1; >= 8 days for a per-event b.  eps = 1 for the counts (S11's
  convention, so the cross-check is exact), eps = 0 for the strictly positive energies/amplitudes,
  and the magnitudes are regressed untransformed (units: magnitude per decade of R).
  Reported per measure: median [IQR] over events, 90 % bootstrap CI of the median (2000 event
  resamples), the pooled event-dummy fit (LSDV within-estimator, OLS and Theil-Sen), and the
  200-draw shifted-tremor null (rigid >= 365-d wrapped shifts of the near-tremor dates inside
  2017-01-01..2026-05-27, slip series untouched — S11's energy-era null).

  F&B DIRECTION: beta = Theil-Sen slope of log10(R7) on log10(A) — the regression F&B actually ran
  (moment rate as the DEPENDENT variable).  Reported for every amplitude proxy, per event and
  pooled, against 3.17 +- 0.29.  Also 1/b_A (the naive inversion of our forward slope) and the OLS
  versions of both directions, because a reverse regression is NOT 1/forward unless r = 1:
  for OLS, beta_reverse = b_forward / r^2, and r^2 here is ~0.1, so the two differ by ~10x.

  CROSS-EVENT check: one point per event (event-median A vs event-mean R7 over its own days),
  n = 18, both directions — closer to F&B's "8 SSEs" granularity than the daily fit.

OUTPUT
  tremor_slip_geometry/results/s16_fb_aggregation.csv   one row per (measure, direction, fit level)

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s16_fb_aggregation.py
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

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree
from scipy.stats import theilslopes, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_RES.mkdir(parents=True, exist_ok=True)
SCRATCH = Path(os.environ.get(
    "S16_SCRATCH",
    "/tmp/claude-1000/-home-jovyan-tremorferometry/"
    "03499166-a677-42d8-9397-e2cf5e28621f/scratchpad/s16"))
SCRATCH.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)

PRE_A = 10
POST_A = 10
NB_RADIUS_KM = 60.0
MIN_FIT_DAYS = 8
N_SHIFT = 200
N_BOOT = 2000
COS47 = np.cos(np.deg2rad(47.0))
ENERGY_ERA = np.datetime64("2017-01-01")
FB_BETA = 3.17
FB_BETA_SD = 0.29

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def seg_of_lat(lat):
    lat = np.asarray(lat, float)
    out = np.full(lat.shape, "south", dtype=object)
    out[lat >= 46.0] = "north"
    out[(lat >= 44.0) & (lat < 46.0)] = "central"
    return out


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


def ols_slope(x, y, nmin=MIN_FIT_DAYS):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < nmin or np.ptp(x[m]) == 0:
        return np.nan, np.nan
    p = np.polyfit(x[m], y[m], 1)
    r = float(np.corrcoef(x[m], y[m])[0, 1])
    return float(p[0]), r


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan, 0
    return (float(np.median(v)), float(np.percentile(v, 25)),
            float(np.percentile(v, 75)), int(v.size))


def boot_median_ci(v, nboot=N_BOOT, lo=5, hi=95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


def group_quantiles(dd, v, nW, qs):
    """Per-day quantiles of v (linear interpolation, == np.percentile) without a python loop."""
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


# ==================================================================================
# 1. Slip model  (loader copied from s11_scaling_robustness.py)
# ==================================================================================
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
np.save(SCRATCH / "cum_mm.npy", CUM_MM.astype(np.float32))
np.save(SCRATCH / "rate7.npy", RATE7.astype(np.float32))
log(f"  Ld {Ld.shape}, A_raw {A_raw.shape}, RATE7 max {RATE7.max():.3f} mm/d (cached to scratch)")


# ==================================================================================
# 2. Geometry / labels
# ==================================================================================
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


# ==================================================================================
# 3. Tremor
# ==================================================================================
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                 usecols=["time", "lat", "lon", "magnitude", "num_stas", "duration", "energy"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_LAT = tr["lat"].values
TR_X, TR_Y = to_km(tr["lon"].values, TR_LAT)
TR_ENERGY = tr["energy"].fillna(0.0).values.astype(float)
TR_NSTA = tr["num_stas"].fillna(0.0).values.astype(float)
TR_MAG = tr["magnitude"].values.astype(float)
TR_DUR = tr["duration"].fillna(0.0).values.astype(float)

SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
# pandas trap: to_datetime gives datetime64[us] here -> cast to ns before epoch arithmetic
TR_DAY = tr["t"].values.astype("datetime64[ns]").astype("datetime64[D]")
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
SHIFT_LO17 = ENERGY_ERA
SPAN17 = int((SHIFT_HI - SHIFT_LO17).astype(int)) + 1
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(f"  {len(tr)} detections; 2017+ wrap span {SPAN17} d")

# --- the catalogue's own energy-magnitude relation (measured, printed, and stored) -------------
m17_all = TR_DAY >= ENERGY_ERA
e17 = TR_ENERGY[m17_all]
g17 = TR_MAG[m17_all]
d17 = TR_DUR[m17_all]
okem = np.isfinite(e17) & (e17 > 0) & np.isfinite(g17)
Xem = np.column_stack([g17[okem], np.log10(np.where(d17[okem] > 0, d17[okem], 150.0)),
                       np.ones(okem.sum())])
bem, *_ = np.linalg.lstsq(Xem, np.log10(e17[okem]), rcond=None)
resid_em = np.log10(e17[okem]) - Xem @ bem
EM_SLOPE = float(bem[0])
A_EXP_FROM_E = 0.5 * EM_SLOPE          # log10 A = A_EXP_FROM_E * M   under  E ~ A^2
log(f"  E-M relation on {okem.sum()} 2017+ detections: log10 E = {bem[0]:.4f} M "
    f"+ {bem[1]:.4f} log10(dur) + {bem[2]:.4f}; resid sd {resid_em.std():.4f} dex; "
    f"-> A = 10^({A_EXP_FROM_E:.3f} M)")

# amplitude proxies: name -> exponent on median magnitude (None = sqrt of median energy)
AMP_PROXIES = {
    "A1_sqrt_Emed": None,
    "A2_mag_0p5": 0.5,
    "A3_mag_EMfit": A_EXP_FROM_E,
    "A4_mag_1p0": 1.0,
    "Asum_sqrt_Esum": "sum",
    "Asum_lin": "sumlin",
}


# ==================================================================================
# 4. Events — the 18 with amplitudes
# ==================================================================================
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
ev_in = ev_in[pd.to_datetime(ev_in["start"]) >= pd.Timestamp("2017-01-01")].reset_index(drop=True)
log(f"  {len(ev_in)} events from 2017-01-01: {ev_in.seg.value_counts().to_dict()}")


def epoch_of(datestr: str) -> str:
    y = int(datestr[:4])
    return "2010-2014" if y <= 2014 else ("2015-2019" if y <= 2019 else "2020-2026")


def window_subset(sorted_off, s, w0_off, nW, span):
    """Detections whose SHIFTED day (o + s) mod span lands inside the window.  (S9/S11 helper.)"""
    a = (w0_off - s) % span
    b = (a + nW - 1) % span
    if a <= b:
        idx = np.arange(np.searchsorted(sorted_off, a, "left"),
                        np.searchsorted(sorted_off, b, "right"))
    else:
        idx = np.concatenate([
            np.arange(np.searchsorted(sorted_off, a, "left"), sorted_off.size),
            np.arange(0, np.searchsorted(sorted_off, b, "right")),
        ])
    d = (sorted_off[idx] + s) % span - w0_off
    return idx, d


QS = (0.10, 0.25, 0.50, 0.75, 0.90)


def daily_measures(dday, pos, nW, arr):
    """All PART-A daily measures for one window."""
    out = {}
    z = np.zeros(nW)
    nan = np.full(nW, np.nan)
    if dday.size == 0:
        for k in ("C", "E_sum", "C_ns8", "A_lin_sum"):
            out[k] = z.copy()
        for k in ("E_med", "E_p10", "E_p25", "E_p75", "E_p90", "E_mean", "M_med", "M_max",
                  "NS_med", "E_med_ns8"):
            out[k] = nan.copy()
        return out
    dd = dday.astype(np.int64)
    out["C"] = np.bincount(dd, minlength=nW).astype(float)
    en = arr["energy"][pos]
    out["E_sum"] = np.bincount(dd, weights=en, minlength=nW)
    out["A_lin_sum"] = np.bincount(dd, weights=np.sqrt(en), minlength=nW)
    qq = group_quantiles(dd, en, nW, QS)
    out["E_p10"], out["E_p25"], out["E_med"], out["E_p75"], out["E_p90"] = (
        qq[0.10], qq[0.25], qq[0.50], qq[0.75], qq[0.90])
    with np.errstate(divide="ignore", invalid="ignore"):
        out["E_mean"] = np.where(out["C"] > 0, out["E_sum"] / np.maximum(out["C"], 1e-12), np.nan)
    mg = arr["mag"][pos]
    okm = np.isfinite(mg)
    out["M_med"] = group_quantiles(dd[okm], mg[okm], nW, (0.50,))[0.50]
    mm = np.full(nW, -np.inf)
    if okm.any():
        np.maximum.at(mm, dd[okm], mg[okm])
    out["M_max"] = np.where(np.isfinite(mm), mm, np.nan)
    ns = arr["nsta"][pos]
    out["NS_med"] = group_quantiles(dd, ns, nW, (0.50,))[0.50]
    s8 = ns >= 8
    out["C_ns8"] = (np.bincount(dd[s8], minlength=nW).astype(float) if s8.any() else z.copy())
    out["E_med_ns8"] = (group_quantiles(dd[s8], en[s8], nW, (0.50,))[0.50] if s8.any()
                        else nan.copy())
    return out


def add_amplitudes(meas):
    for name, ex in AMP_PROXIES.items():
        if ex is None:
            meas[name] = np.sqrt(meas["E_med"])
        elif ex == "sum":
            meas[name] = np.sqrt(meas["E_sum"])
        elif ex == "sumlin":
            meas[name] = meas["A_lin_sum"]
        else:
            meas[name] = 10.0 ** (ex * meas["M_med"])
    meas["E_sum_epsS11"] = meas["E_sum"]          # same series, S11's eps (see MEAS_SPEC)
    return meas


# S11 regressed E with eps = the global 1st percentile of positive daily E over the same 18
# events' days.  It is read from results/s11_exponents.csv (and re-verified from this script's own
# POOL at the end) so `E_sum_epsS11` reproduces S11's b_E exactly, while `E_sum` uses eps = 0,
# which is the correct choice for a strictly positive measure.
try:
    _s11x = pd.read_csv(OUT_RES / "s11_exponents.csv")
    EPS_S11_E = float(_s11x.loc[_s11x.measure == "E", "eps"].iloc[0])
except Exception:                                            # pragma: no cover
    EPS_S11_E = 0.0
log(f"  S11 eps for E = {EPS_S11_E:.4f} (from results/s11_exponents.csv)")

# measure -> (eps, transform)   transform "log" | "none"
MEAS_SPEC = {
    "C": (1.0, "log"),
    "C_ns8": (1.0, "log"),
    "E_sum": (0.0, "log"),
    "E_sum_epsS11": (EPS_S11_E, "log"),
    "E_med": (0.0, "log"),
    "E_mean": (0.0, "log"),
    "E_p10": (0.0, "log"),
    "E_p25": (0.0, "log"),
    "E_p75": (0.0, "log"),
    "E_p90": (0.0, "log"),
    "E_med_ns8": (0.0, "log"),
    "M_med": (0.0, "none"),
    "M_max": (0.0, "none"),
    "NS_med": (0.0, "none"),
}
for name in AMP_PROXIES:
    MEAS_SPEC[name] = (0.0, "log")
MEAS_NAMES = list(MEAS_SPEC)
AMP_NAMES = list(AMP_PROXIES)


def yvals(meas, name):
    eps, trans = MEAS_SPEC[name]
    v = np.asarray(meas[name], float)
    if trans == "log":
        good = np.isfinite(v) & (v > 0)
        return np.where(good, np.log10(np.where(good, v, 1.0) + eps), np.nan)
    return np.where(np.isfinite(v), v, np.nan)


def fit_window(meas, R7):
    """Forward b (measure on rate), reverse beta (rate on amplitude), RMA slope, for one window."""
    C = meas["C"]
    base = (R7 > 0) & (C >= 1)
    lr = np.where(R7 > 0, np.log10(np.where(R7 > 0, R7, 1.0)), np.nan)
    fwd, fwd_ols, fwd_r, rev, rev_ols, rma, ndays = {}, {}, {}, {}, {}, {}, {}
    for name in MEAS_NAMES:
        y = yvals(meas, name)
        ok = base & np.isfinite(y)
        ndays[name] = int(ok.sum())
        fwd[name] = ts_slope_fast(lr[ok], y[ok])
        o, r = ols_slope(lr[ok], y[ok])
        fwd_ols[name], fwd_r[name] = o, r
        if name in AMP_NAMES:
            rev[name] = ts_slope_fast(y[ok], lr[ok])
            rev_ols[name] = ols_slope(y[ok], lr[ok])[0]
            # reduced major axis (geometric-mean) slope: direction-symmetric, errors in both vars
            if ok.sum() >= MIN_FIT_DAYS and np.std(lr[ok]) > 0:
                rma[name] = float(np.sign(r) * np.std(y[ok]) / np.std(lr[ok]))
            else:
                rma[name] = np.nan
    return fwd, fwd_ols, fwd_r, rev, rev_ols, rma, ndays


# ==================================================================================
# 5. Pass 1 — per-event daily series, observed fits, nulls
# ==================================================================================
log("pass 1: per-event daily series + observed fits + 200-draw nulls ...")
EV = []
POOL_ROWS = []
OBS_F, OBS_FO, OBS_FR, OBS_R, OBS_RO, OBS_RMA, NDAYS = {}, {}, {}, {}, {}, {}, {}
NUL_F = {n: {} for n in MEAS_NAMES}
NUL_R = {n: {} for n in AMP_NAMES}
BINS = (1, 3, 7)
BIN_FWD = {b: {n: {} for n in AMP_NAMES} for b in BINS}
BIN_REV = {b: {n: {} for n in AMP_NAMES} for b in BINS}
for n in MEAS_NAMES:
    OBS_F[n], OBS_FO[n], OBS_FR[n], NDAYS[n] = {}, {}, {}, {}
for n in AMP_NAMES:
    OBS_R[n], OBS_RO[n], OBS_RMA[n] = {}, {}, {}

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
    near = near[TR_DAY[near] >= ENERGY_ERA]          # amplitudes only exist from 2017

    aw = PAREA[nb][:, None]
    R7 = (aw * RATE7[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)

    soff = (TR_DAY[near] - SHIFT_LO17).astype(int)
    order = np.argsort(soff, kind="stable")
    idx_sorted = near[order]
    soff = soff[order]
    arr = dict(energy=TR_ENERGY[idx_sorted], mag=TR_MAG[idx_sorted], nsta=TR_NSTA[idx_sorted])
    w0_off = int((np.datetime64(DAYS[w0]) - SHIFT_LO17).astype(int))

    sub, dday = window_subset(soff, 0, w0_off, nW, SPAN17)
    meas = add_amplitudes(daily_measures(dday, sub, nW, arr))
    fwd, fwd_ols, fwd_r, rev, rev_ols, rma, ndays = fit_window(meas, R7)
    for n in MEAS_NAMES:
        OBS_F[n][pid] = fwd[n]
        OBS_FO[n][pid] = fwd_ols[n]
        OBS_FR[n][pid] = fwd_r[n]
        NDAYS[n][pid] = ndays[n]
    for n in AMP_NAMES:
        OBS_R[n][pid] = rev[n]
        OBS_RO[n][pid] = rev_ols[n]
        OBS_RMA[n][pid] = rma[n]

    # --- TIME-BIN LADDER: re-aggregate the SAME detections into 1/3/7-day bins -----------
    for nb_d in BINS:
        bid = (dday.astype(np.int64) // nb_d)
        nB = int(np.ceil(nW / nb_d))
        mb = add_amplitudes(daily_measures(bid, sub, nB, arr))
        Rb = np.array([R7[i * nb_d:(i + 1) * nb_d].mean() for i in range(nB)])
        fb, _, _, rb, _, _, _ = fit_window(mb, Rb)
        for n in AMP_NAMES:
            BIN_FWD[nb_d][n][pid] = fb[n]
            BIN_REV[nb_d][n][pid] = rb[n]

    # --- shifted-tremor null (S11 energy-era wrap) --------------------------------------
    store_f = {n: [] for n in MEAS_NAMES}
    store_r = {n: [] for n in AMP_NAMES}
    shifts = RNG.integers(365, max(366, SPAN17 - 365 + 1), size=N_SHIFT)
    for sh in shifts:
        sub_n, dd_n = window_subset(soff, int(sh), w0_off, nW, SPAN17)
        mn = add_amplitudes(daily_measures(dd_n, sub_n, nW, arr))
        f2, _, _, r2, _, _, _ = fit_window(mn, R7)
        for n in MEAS_NAMES:
            store_f[n].append(f2[n])
        for n in AMP_NAMES:
            store_r[n].append(r2[n])
    for n in MEAS_NAMES:
        NUL_F[n][pid] = np.asarray(store_f[n], float)
    for n in AMP_NAMES:
        NUL_R[n][pid] = np.asarray(store_r[n], float)

    dfr = dict(id=pid, seg=seg, epoch=epoch_of(r["start"]), day=np.arange(nW), R7=R7)
    for k in MEAS_NAMES:
        dfr[k] = meas[k]
    POOL_ROWS.append(pd.DataFrame(dfr))
    EV.append(dict(pid=pid, seg=seg, start=str(r["start"]), Mw=float(r["Mw"]),
                   epoch=epoch_of(r["start"]), nW=nW, n_near=int(near.size), meas=meas, R7=R7))
    log(f"  {pid:5d} {r['start']} {seg:7s} near={near.size:6d} nd={ndays['E_med']:3d} "
        f"b(C)={fwd['C']:+6.3f} b(Esum)={fwd['E_sum']:+6.3f} b(Emed)={fwd['E_med']:+6.3f} "
        f"b(A1)={fwd['A1_sqrt_Emed']:+6.3f} beta(A1)={rev['A1_sqrt_Emed']:+8.2f}")

POOL = pd.concat(POOL_ROWS, ignore_index=True)
log(f"  {len(EV)} events, {len(POOL)} event-days")


# ==================================================================================
# 6. Pooled event-dummy fits
# ==================================================================================
def pooled(df, name, reverse=False):
    """LSDV within-estimator: demean both axes per event.  reverse -> log R on log A."""
    y = yvals({name: df[name].values}, name)
    R = df["R7"].values
    C = df["C"].values
    ok = np.isfinite(y) & (R > 0) & (C >= 1)
    if ok.sum() < 20:
        return np.nan, np.nan, np.nan, int(ok.sum())
    x = np.log10(np.where(R > 0, R, 1.0))
    d = pd.DataFrame(dict(i=df["id"].values[ok], x=x[ok], y=y[ok]))
    if reverse:
        d = d.rename(columns={"x": "y", "y": "x"})
    d["xc"] = d.x - d.groupby("i").x.transform("mean")
    d["yc"] = d.y - d.groupby("i").y.transform("mean")
    if np.ptp(d.xc.values) == 0:
        return np.nan, np.nan, np.nan, len(d)
    o = float(np.polyfit(d.xc.values, d.yc.values, 1)[0])
    t = ts_slope_fast(d.xc.values, d.yc.values, 20)
    rr = float(np.corrcoef(d.xc.values, d.yc.values)[0, 1])
    return o, t, rr, len(d)


# ==================================================================================
# 7. Assemble the output table
# ==================================================================================
log("assembling results ...")
rows = []


def null_stats(nuld, obsd):
    """Per-event one-sided p (fraction of draws >= obs), group sign test on z, pooled null moments."""
    zs, ps = [], []
    allnull = []
    for pid, ob in obsd.items():
        nl = nuld.get(pid)
        if nl is None or not np.isfinite(ob):
            continue
        nl = nl[np.isfinite(nl)]
        if nl.size < 20:
            continue
        allnull.append(nl)
        ps.append(float((nl >= ob).mean()))
        sd = nl.std(ddof=1)
        zs.append((ob - nl.mean()) / sd if sd > 0 else np.nan)
    allnull = np.concatenate(allnull) if allnull else np.array([])
    zs = np.asarray(zs, float)
    zs_ok = zs[np.isfinite(zs)]
    npos = int((zs_ok > 0).sum())
    n = zs_ok.size
    sign_p = np.nan
    if n:
        from scipy.stats import binomtest
        sign_p = float(binomtest(npos, n, 0.5, alternative="greater").pvalue)
    return dict(null_mean=float(allnull.mean()) if allnull.size else np.nan,
                null_sd=float(allnull.std(ddof=1)) if allnull.size else np.nan,
                z_median=float(np.median(zs_ok)) if n else np.nan,
                n_events_z=n, n_z_positive=npos, sign_p=sign_p,
                n_events_p05=int(np.sum(np.asarray(ps) < 0.05)))


def add_row(measure, direction, level, bs, nuld=None, obsd=None, pooled_tuple=None, **extra):
    m, q1, q3, n = med_iqr(bs)
    lo, hi = boot_median_ci(bs)
    row = dict(measure=measure, direction=direction, level=level, n_events=n,
               b_median=m, b_q25=q1, b_q75=q3, b_boot_lo90=lo, b_boot_hi90=hi,
               n_below_1=int(np.sum(np.asarray(bs, float)[np.isfinite(bs)] < 1.0)))
    if pooled_tuple is not None:
        o, t, rr, npool = pooled_tuple
        row.update(pooled_ols=o, pooled_theilsen=t, pooled_r=rr, n_pooled_days=npool)
    if nuld is not None and obsd is not None:
        row.update(null_stats(nuld, obsd))
    row.update(extra)
    rows.append(row)


# ---- forward exponents (measure on log rate) -------------------------------------------
for name in MEAS_NAMES:
    bs = np.array([OBS_F[name][p] for p in OBS_F[name]], float)
    nd = np.array([NDAYS[name][p] for p in NDAYS[name]], float)
    add_row(name, "forward_measure_on_rate", "per_event_theilsen", bs,
            nuld=NUL_F[name], obsd=OBS_F[name],
            pooled_tuple=pooled(POOL, name),
            median_n_days=float(np.nanmedian(nd)),
            eps=MEAS_SPEC[name][0], transform=MEAS_SPEC[name][1])
    bo = np.array([OBS_FO[name][p] for p in OBS_FO[name]], float)
    br = np.array([OBS_FR[name][p] for p in OBS_FR[name]], float)
    add_row(name, "forward_measure_on_rate", "per_event_ols", bo,
            median_r=float(np.nanmedian(br)),
            median_r2=float(np.nanmedian(br ** 2)))

# ---- reverse (F&B direction): log rate on log amplitude --------------------------------
for name in AMP_NAMES:
    bs = np.array([OBS_R[name][p] for p in OBS_R[name]], float)
    fwd = np.array([OBS_F[name][p] for p in OBS_F[name]], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        naive = 1.0 / fwd
    add_row(name, "reverse_rate_on_amplitude_FB", "per_event_theilsen", bs,
            nuld=NUL_R[name], obsd=OBS_R[name],
            pooled_tuple=pooled(POOL, name, reverse=True),
            fb_beta=FB_BETA, fb_beta_sd=FB_BETA_SD,
            naive_inv_fwd_median=float(np.nanmedian(naive)))
    bo = np.array([OBS_RO[name][p] for p in OBS_RO[name]], float)
    add_row(name, "reverse_rate_on_amplitude_FB", "per_event_ols", bo, fb_beta=FB_BETA)

# ---- RMA (direction-symmetric) slopes --------------------------------------------------
for name in AMP_NAMES:
    bs = np.array([OBS_RMA[name][p] for p in OBS_RMA[name]], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / bs
    add_row(name, "rma_amp_on_rate", "per_event_rma", bs,
            fb_beta=FB_BETA, fb_beta_sd=FB_BETA_SD,
            beta_rma_median=float(np.nanmedian(inv)),
            note="reduced major axis = sign(r)*sd(logA)/sd(logR); beta_rma = 1/b_rma")

# ---- time-bin ladder -------------------------------------------------------------------
for nb_d in BINS:
    for name in AMP_NAMES:
        bs = np.array([BIN_FWD[nb_d][name][p] for p in BIN_FWD[nb_d][name]], float)
        add_row(name, "forward_measure_on_rate", f"bin{nb_d}d_theilsen", bs, bin_days=nb_d)
        rs = np.array([BIN_REV[nb_d][name][p] for p in BIN_REV[nb_d][name]], float)
        add_row(name, "reverse_rate_on_amplitude_FB", f"bin{nb_d}d_theilsen", rs,
                bin_days=nb_d, fb_beta=FB_BETA)

# ---- cross-event fits (one point per event) --------------------------------------------
ev_lvl = []
for e in EV:
    ok = (e["R7"] > 0) & (e["meas"]["C"] >= 1)
    if ok.sum() < MIN_FIT_DAYS:
        continue
    d = dict(id=e["pid"], seg=e["seg"], start=e["start"], Mw=e["Mw"], epoch=e["epoch"],
             logR_mean=float(np.mean(np.log10(e["R7"][ok]))),
             C_sum=float(np.nansum(e["meas"]["C"][ok])))
    for name in ("E_med", "E_mean", "E_sum", "M_med") + tuple(AMP_NAMES):
        v = np.asarray(e["meas"][name], float)[ok]
        v = v[np.isfinite(v)]
        d[name] = float(np.median(v)) if v.size else np.nan
    ev_lvl.append(d)
EVL = pd.DataFrame(ev_lvl)
for name in AMP_NAMES:
    x = np.log10(EVL[name].values)
    yr = EVL["logR_mean"].values
    b_f = ts_slope_fast(yr, x, 8)
    b_r = ts_slope_fast(x, yr, 8)
    o_f, r_f = ols_slope(yr, x, 8)
    o_r, _ = ols_slope(x, yr, 8)
    rows.append(dict(measure=name, direction="cross_event_amp_on_rate", level="n18_theilsen",
                     n_events=int(np.isfinite(x).sum()), b_median=b_f, pooled_ols=o_f,
                     pooled_r=r_f))
    rows.append(dict(measure=name, direction="cross_event_rate_on_amp_FB", level="n18_theilsen",
                     n_events=int(np.isfinite(x).sum()), b_median=b_r, pooled_ols=o_r,
                     pooled_r=r_f, fb_beta=FB_BETA, fb_beta_sd=FB_BETA_SD))

# ---- paired differences that matter -----------------------------------------------------
def paired(a_name, b_name, dct=OBS_F, tag=""):
    pids = [p for p in dct[a_name] if np.isfinite(dct[a_name][p]) and np.isfinite(dct[b_name][p])]
    d = np.array([dct[a_name][p] - dct[b_name][p] for p in pids], float)
    if d.size < 4:
        return
    lo, hi = boot_median_ci(d)
    try:
        wp = float(wilcoxon(d).pvalue)
    except Exception:
        wp = np.nan
    rows.append(dict(measure=f"{a_name}_minus_{b_name}{tag}", direction="paired_difference",
                     level="per_event_theilsen", n_events=d.size, b_median=float(np.median(d)),
                     b_q25=float(np.percentile(d, 25)), b_q75=float(np.percentile(d, 75)),
                     b_boot_lo90=lo, b_boot_hi90=hi, wilcoxon_p=wp,
                     n_positive=int((d > 0).sum())))


paired("E_sum", "C")
paired("E_sum_epsS11", "E_sum")
paired("Asum_lin", "Asum_sqrt_Esum")
paired("E_sum", "E_med")
paired("E_med", "E_mean")
paired("E_p90", "E_med")
paired("E_p10", "E_med")
paired("E_med_ns8", "E_med")

# ---- the E-M relation row ---------------------------------------------------------------
rows.append(dict(measure="catalog_EM_relation", direction="diagnostic", level="all_2017plus",
                 n_events=int(okem.sum()), b_median=EM_SLOPE,
                 pooled_r=float(np.corrcoef(g17[okem], np.log10(e17[okem]))[0, 1]),
                 note=(f"log10 E = {bem[0]:.4f}*M + {bem[1]:.4f}*log10(dur) + {bem[2]:.4f}; "
                       f"resid sd {resid_em.std():.4f} dex; amplitude exponent on M "
                       f"= {A_EXP_FROM_E:.4f} under E ~ A^2")))

_ep = POOL["E_sum"].values
_ep = _ep[np.isfinite(_ep) & (_ep > 0)]
EPS_S11_E_RECOMP = float(np.percentile(_ep, 1))
log(f"  S11 eps for E: from CSV {EPS_S11_E:.4f}  recomputed from this POOL "
    f"{EPS_S11_E_RECOMP:.4f}  (diff {abs(EPS_S11_E - EPS_S11_E_RECOMP):.3e})")
rows.append(dict(measure="eps_S11_E_check", direction="diagnostic", level="all_2017plus",
                 b_median=EPS_S11_E, pooled_ols=EPS_S11_E_RECOMP,
                 note="S11's eps for E, from results/s11_exponents.csv vs recomputed here"))

OUT = pd.DataFrame(rows)
OUT.to_csv(OUT_RES / "s16_fb_aggregation.csv", index=False)
log(f"  wrote results/s16_fb_aggregation.csv ({len(OUT)} rows)")

# per-event table for the figures / caption numbers
per_ev = []
for e in EV:
    okd = (e["R7"] > 0) & (e["meas"]["C"] >= 1) & np.isfinite(e["meas"]["A1_sqrt_Emed"])
    la = np.log10(e["meas"]["A1_sqrt_Emed"][okd])
    lrr = np.log10(e["R7"][okd])
    d = dict(id=e["pid"], start=e["start"], seg=e["seg"], Mw=e["Mw"], epoch=e["epoch"],
             n_near=e["n_near"], n_days=NDAYS["E_med"][e["pid"]],
             sd_logA1=float(np.std(la)) if la.size else np.nan,
             ptp_logA1=float(np.ptp(la)) if la.size else np.nan,
             sd_logR7=float(np.std(lrr)) if lrr.size else np.nan,
             ptp_logR7=float(np.ptp(lrr)) if lrr.size else np.nan)
    for n in MEAS_NAMES:
        d[f"b_{n}"] = OBS_F[n][e["pid"]]
    for n in AMP_NAMES:
        d[f"beta_{n}"] = OBS_R[n][e["pid"]]
    per_ev.append(d)
pd.DataFrame(per_ev).to_csv(OUT_RES / "s16_fb_events.csv", index=False)
log("  wrote results/s16_fb_aggregation.csv + results/s16_fb_events.csv")


# ==================================================================================
# 8. Report
# ==================================================================================
def line(name, direction="forward_measure_on_rate", level="per_event_theilsen"):
    r = OUT[(OUT.measure == name) & (OUT.direction == direction) & (OUT.level == level)]
    if r.empty:
        return f"   {name:16s}  n/a"
    r = r.iloc[0]
    return (f"   {name:16s} n={int(r.n_events):2d}  {r.b_median:+8.3f} "
            f"[{r.b_q25:+7.3f},{r.b_q75:+7.3f}]  CI90 [{r.b_boot_lo90:+7.3f},"
            f"{r.b_boot_hi90:+7.3f}]  pooled OLS/TS {r.get('pooled_ols', np.nan):+6.3f}/"
            f"{r.get('pooled_theilsen', np.nan):+6.3f}  null {r.get('null_mean', np.nan):+6.3f}"
            f"+-{r.get('null_sd', np.nan):5.3f}  sign {int(r.get('n_z_positive', 0)):2d}/"
            f"{int(r.get('n_events_z', 0)):2d} p={r.get('sign_p', np.nan):.2g}")


print("\n" + "=" * 116)
print("S16 PART A — Frank & Brodsky (2019) AGGREGATION MATCH")
print("=" * 116)
print(f"18 confirmed events from 2017-01-01 (energy/magnitude era); window [start-{PRE_A}, "
      f"end+{POST_A}]; near tremor within {NB_RADIUS_KM:.0f} km of the footprint;")
print("regressor = R7 (7-d centred mean of the area-weighted positive modeled potency rate over "
      "the 60-km neighbourhood).")
print(f"\nCatalogue E-M relation (measured on {okem.sum()} 2017+ detections): "
      f"log10 E = {bem[0]:.3f} M + {bem[2]:.3f}  (r = "
      f"{np.corrcoef(g17[okem], np.log10(e17[okem]))[0, 1]:.4f}, resid sd "
      f"{resid_em.std():.3f} dex)")
print(f"  -> `energy` is a deterministic monotone function of an unrounded magnitude; "
      f"amplitude ~ 10^({A_EXP_FROM_E:.3f} M) under E ~ A^2.")
print(f"  -> the requested proxy 10^(0.5 M) assumes log10 E = 1.0 M and therefore compresses the "
      f"amplitude range by 0.5/{A_EXP_FROM_E:.3f} = {0.5 / A_EXP_FROM_E:.3f}.")

print("\n(1) FORWARD exponents b: Theil-Sen slope of log10(measure) on log10(R7), median [IQR]")
for n in ("C", "C_ns8", "E_sum", "E_sum_epsS11", "E_med", "E_mean", "E_p10", "E_p25", "E_p75",
          "E_p90", "E_med_ns8"):
    print(line(n))
print("   -- magnitudes are regressed UNTRANSFORMED (magnitude units per decade of R):")
for n in ("M_med", "M_max", "NS_med"):
    print(line(n))
print("   -- amplitude proxies (log10 A on log10 R) -> compare with F&B's 1/3.17 = 0.315:")
for n in AMP_NAMES:
    print(line(n))

print("\n(2) F&B DIRECTION: beta = Theil-Sen slope of log10(R7) on log10(A); F&B report "
      f"{FB_BETA} +- {FB_BETA_SD}")
for n in AMP_NAMES:
    print(line(n, "reverse_rate_on_amplitude_FB"))
print("   (the naive 1/b_forward is in the CSV column naive_inv_fwd_median; a reverse regression "
      "is 1/forward only when r = 1)")
r2s = OUT[(OUT.level == "per_event_ols") & (OUT.measure == "A1_sqrt_Emed") &
          (OUT.direction == "forward_measure_on_rate")]
if not r2s.empty:
    print(f"   median per-event r^2 of log A1 on log R7 = "
          f"{float(r2s.iloc[0]['median_r2']):.4f}  -> OLS beta_reverse = b_forward / r^2 "
          f"= {float(OUT[(OUT.measure=='A1_sqrt_Emed') & (OUT.level=='per_event_ols') & (OUT.direction=='forward_measure_on_rate')].iloc[0]['b_median']) / max(float(r2s.iloc[0]['median_r2']), 1e-9):.1f}")

print("\n(2b) DIRECTION-SYMMETRIC (RMA / geometric-mean) slope, and 1/RMA vs F&B's 3.17")
for n in AMP_NAMES:
    r = OUT[(OUT.measure == n) & (OUT.direction == "rma_amp_on_rate")]
    if not r.empty:
        r = r.iloc[0]
        print(f"   {n:16s} b_rma {r.b_median:+7.3f} [{r.b_q25:+6.3f},{r.b_q75:+6.3f}] "
              f"CI90 [{r.b_boot_lo90:+6.3f},{r.b_boot_hi90:+6.3f}]  "
              f"beta_rma = 1/b_rma {float(r.beta_rma_median):+7.2f}")

print("\n(2c) TIME-BIN LADDER (forward b of log A on log R; bin = 1 / 3 / 7 d)")
for n in AMP_NAMES:
    vals = []
    for nb_d in BINS:
        r = OUT[(OUT.measure == n) & (OUT.direction == "forward_measure_on_rate") &
                (OUT.level == f"bin{nb_d}d_theilsen")]
        vals.append(float(r.iloc[0].b_median) if not r.empty else np.nan)
    revs = []
    for nb_d in BINS:
        r = OUT[(OUT.measure == n) & (OUT.direction == "reverse_rate_on_amplitude_FB") &
                (OUT.level == f"bin{nb_d}d_theilsen")]
        revs.append(float(r.iloc[0].b_median) if not r.empty else np.nan)
    print(f"   {n:16s} b {vals[0]:+7.3f} {vals[1]:+7.3f} {vals[2]:+7.3f}   |   "
          f"beta {revs[0]:+7.2f} {revs[1]:+7.2f} {revs[2]:+7.2f}")

print("\n(3) CROSS-EVENT (one point per event, n = 18)")
for n in AMP_NAMES:
    a = OUT[(OUT.measure == n) & (OUT.direction == "cross_event_amp_on_rate")]
    b = OUT[(OUT.measure == n) & (OUT.direction == "cross_event_rate_on_amp_FB")]
    if not a.empty:
        print(f"   {n:16s} b(A on R) {float(a.iloc[0].b_median):+7.3f}  "
              f"beta(R on A) {float(b.iloc[0].b_median):+8.2f}  r {float(a.iloc[0].pooled_r):+5.2f}")

print("\n(4) PAIRED DIFFERENCES")
for r in OUT[OUT.direction == "paired_difference"].itertuples():
    print(f"   {r.measure:26s} n={int(r.n_events):2d}  {r.b_median:+7.3f} "
          f"[{r.b_boot_lo90:+7.3f},{r.b_boot_hi90:+7.3f}]  Wilcoxon p {r.wilcoxon_p:.4g}  "
          f"{int(r.n_positive)} positive")

# cross-check against S11
try:
    s11e = pd.read_csv(OUT_RES / "s11_energy_events.csv")
    j = s11e.merge(pd.DataFrame(per_ev), on="id", suffixes=("_s11", "_s16"))
    print("\n(5) CROSS-CHECK against results/s11_energy_events.csv (must be identical):")
    for a, b in (("b_C_R7", "b_C"), ("b_E_R7", "b_E_sum_epsS11"), ("b_E_R7", "b_E_sum"),
                 ("b_EperDet_R7", "b_E_mean"), ("b_Mmax_R7", "b_M_max"),
                 ("b_C_ns8_R7", "b_C_ns8")):
        d = (j[a].values - j[b].values)
        print(f"   {a:16s} vs {b:12s} n={len(j)}  median S11 {np.nanmedian(j[a]):+7.4f}  "
              f"S16 {np.nanmedian(j[b]):+7.4f}  max|diff| {np.nanmax(np.abs(d)):.2e}")
except Exception as exc:                                     # pragma: no cover
    print(f"\n(5) cross-check unavailable: {exc}")

# scipy agreement on the fast Theil-Sen
e0 = EV[0]
ok0 = (e0["R7"] > 0) & (e0["meas"]["C"] >= 1) & np.isfinite(e0["meas"]["E_med"])
x0 = np.log10(e0["R7"][ok0])
y0 = np.log10(e0["meas"]["E_med"][ok0])
print(f"\n   CHECK Theil-Sen fast {ts_slope_fast(x0, y0):.10f} vs scipy "
      f"{theilslopes(y0, x0)[0]:.10f}")
PEV = pd.DataFrame(per_ev)
print("\n(6) DYNAMIC RANGE of the two axes within an event (this is what sets the attenuation):")
print(f"   sd(log10 A1_median amplitude)  median {PEV.sd_logA1.median():.4f} dex  "
      f"[{PEV.sd_logA1.quantile(.25):.4f}, {PEV.sd_logA1.quantile(.75):.4f}]  "
      f"peak-to-peak {PEV.ptp_logA1.median():.3f} dex")
print(f"   sd(log10 R7 slip rate)         median {PEV.sd_logR7.median():.4f} dex  "
      f"[{PEV.sd_logR7.quantile(.25):.4f}, {PEV.sd_logR7.quantile(.75):.4f}]  "
      f"peak-to-peak {PEV.ptp_logR7.median():.3f} dex")
print(f"   ratio sd(A)/sd(R) = b_rma      {PEV.sd_logA1.median() / PEV.sd_logR7.median():.4f}")
print("   -> a PNSN detection is a 300-s tremor WINDOW whose amplitude is a network-average "
      "envelope level;")
print(f"      its median varies only {PEV.sd_logA1.median():.2f} dex within an event, against "
      f"{PEV.sd_logR7.median():.2f} dex for the modeled slip rate.")
log("PART A done")
