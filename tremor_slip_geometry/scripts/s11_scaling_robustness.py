#!/usr/bin/env python
"""S11 — is the sub-linear tremor-vs-slip scaling (S9: C ~ R^0.40) REAL, or a measurement artefact?

PLAN.md Round 3, step S11.  One script, five tests, the same 27 confirmed north/central events
that S9 used (`results/s1_events.csv`).

Background (NOTES.md, 2026-09-10 S9 Part B and the orchestrator review that follows it):
  S9 measured b = Theil-Sen slope of log10(C+1) on log10(R7), C = daily tremor detections within
  60 km of the event footprint, R7 = 7-d centred mean of the area-weighted positive modeled slip
  rate over the 60-km neighbourhood.  Result: b = 0.401 [0.181, 0.478], 26/27 events positive
  against a shifted-tremor null (null b = -0.00 +- 0.28, sign p = 4.2e-07), the same in every
  segment and epoch.  TWO OPPOSITE BIASES were flagged and never quantified:
    (i)  R7 is 7-d SMOOTHED and mixture-limited (13 fixed vbICA patterns), so its daily dynamic
         range is compressed.  Regressing y on a shrunken x INFLATES the slope -> b biased HIGH.
    (ii) C is a DETECTION count on a catalogue whose sensitivity grows with the PNSN network, and
         PNSN detection times sit on a 2.5-min grid, so a day could in principle saturate ->
         b biased LOW at high slip rates.
  The orchestrator's own check (daily counts reach 1011, p99 695, i.e. > 576 windows/day) showed
  the hard 576-window cap is NOT reached by the COUNT, because several epicentres per window are
  catalogued.  S11 turns that observation into a measurement: it separates the count into a
  window-occupancy part (U, which CAN saturate at 576/day) and a multiplicity part (C/U), adds
  amplitude-like measures (energy, duration, max magnitude, populated 2017+), varies how the slip
  rate is smoothed, and simulates the smoothing bias directly.

======================================================================================
INPUTS (read-only; nothing outside tremor_slip_geometry/ is written)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv                        27 confirmed north/central
                                                                    patches (id, start, end, seg, Mw)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz   label (T x 742), days (T),
                                                                    patch_lat/lon/area_km2
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                            PNSN tremor.  Columns time, lat,
                                                                    lon, depth, magnitude, num_stas,
                                                                    duration, energy, id.  ENERGY /
                                                                    DURATION / NUM_STAS / MAGNITUDE
                                                                    ARE POPULATED ONLY FROM 2017-01-01
                                                                    (64.5 % of rows; exactly 0/NaN
                                                                    before, verified per year).
  /home/jovyan/cascadia-sse-strain                                  mesh + sse_strain loaders

  Loader code (slip model, mesh geometry, tremor, km frame, event table) is COPIED VERBATIM from
  scripts/s9_within_event.py, which copied it from scripts/s4_fronts.py.
    cum_mm = -1e3 * (Ld @ A_raw.T)       (742 x T, mm, thrust positive)
    x = (lon + 124) * 111 * cos(47 deg) km E,  y = (lat - 47) * 111 km N
    segments: north lat >= 46, central 44 <= lat < 46

======================================================================================
GEOMETRY (identical to S9)
======================================================================================
  Per event: window W = [start - 10, end + 10] day indices on the model day grid.
  Footprint F = triangles labelled with the patch id on any day of [ts, te].
  Neighbourhood NB = triangles within 60 km of any F centroid (F always included).
  Near tremor = detections within 60 km of any F centroid (epicentral).

======================================================================================
TREMOR MEASURES per event-day  (all from the SAME near-tremor set)
======================================================================================
  C      detection count                                       (S9's measure)
  U      # distinct 2.5-min windows with >= 1 detection.  Detection times sit exactly on 150-s
         multiples (verified: 100.0 % of rows), so a day has at most 576 windows.  U is a
         DURATION-like measure that saturates at 576 and cannot be inflated by multiple
         simultaneous epicentres.
  Uloc   # distinct (2.5-min window, 0.1-deg lat/lon cell) pairs.  Saturates only if every cell in
         the neighbourhood is active in every window; a space-time occupancy measure.
  CperU  C / U = mean detections per active window = the MULTIPLICITY that lets C exceed 576.
  E      sum of `energy` over the day        (2017+ events only)
  Dur    sum of `duration` over the day      (2017+ events only; the column takes only two values,
                                              150 s and 300 s, so Dur is close to 150*C by
                                              construction -- reported for completeness, not as an
                                              independent measure)
  Mmax   daily maximum `magnitude`           (2017+ events only; already a log quantity, so it is
                                              regressed UNTRANSFORMED -> units = magnitude/decade)
  EperDet  E / C = mean energy per detection (2017+; tests whether counts under-represent
                                              high-rate tremor)
  C_ns8  count using only detections with num_stas >= 8 (2017+; a network-sensitivity control)

======================================================================================
SLIP-RATE VARIANTS  (all = area-weighted sums over NB, units mm/d km^2 = a potency rate)
======================================================================================
  R7    clip(7-d centred moving average of d/dt cum_mm, 0) summed with area weights   <- S9's
  R3    clip(3-d centred moving average of d/dt cum_mm, 0) summed with area weights
  R1    clip(raw daily d/dt cum_mm, 0) summed with area weights                (per-triangle clip)
  Rcum  daily increment of the RAW cumulative potency P(d) = sum_NB area * cum_mm(d) (SIGNED sum,
        differenced, then the increment clipped at 0)                          (post-sum clip)
  R1 and Rcum differ only in when the positive part is taken: R1 clips each triangle, Rcum clips
  the neighbourhood total, so Rcum <= R1 and Rcum is zero on days when the neighbourhood is net
  back-slipping.

======================================================================================
TEST 1 — exponent matrix (measure x slip variant x day-subset)
======================================================================================
  b = Theil-Sen slope of log10(measure + eps) on log10(R) over event-days with measure > 0 and
  R > 0; eps = 1 for the count-like measures (C, U, Uloc, C_ns8), the measure's GLOBAL 1st
  percentile of positive values for E and Dur, 0 for the ratios CperU and EperDet (both >= their
  own floor and strictly positive), and Mmax is not transformed at all.
  >= 8 usable days for a per-event b.  Reported: median and IQR over events, a 90 % bootstrap CI
  of that median (2000 resamples OF EVENTS), and a POOLED fit with event dummies (within-event
  demeaning of both axes, i.e. the LSDV within-estimator; OLS and Theil-Sen versions).
  Day-subsets: `all`; `no_sat` (drop days with U >= 400, i.e. within 70 % of the 576 cap);
  `R_lower` / `R_upper` (days below / above the event's own median log10 R, >= 6 days each).
  NULL: 200 rigid shifts of the near-tremor DATES by >= 365 d, wrapped, with the SLIP series left
  unchanged -- exactly the S4/S9 control -- rerun through the whole measure pipeline, for every
  measure x slip-variant x subset combination.  p = fraction of draws with b_null >= b_obs.
  For the 2017+ measures (E, Dur, Mmax, EperDet, C_ns8) the null uses a SEPARATE wrap over
  2017-01-01..2026-05-27 with only 2017+ detections, so a shifted window never lands in the years
  where those columns are structurally zero.

======================================================================================
TEST 2 — saturation diagnostics
======================================================================================
  * fraction of event-days with U >= 500 and U >= 400 (the 576/day ceiling), per event and pooled;
  * b(C vs R) recomputed on the `no_sat` subset (days with U < 400) -- if saturation were bending
    the relation, removing the saturated days would RAISE b;
  * b(C vs U): the exponent of the count against the occupancy.  b > 1 means C keeps growing after
    U stops growing, i.e. MULTIPLICITY (several catalogued epicentres in the same 2.5-min window)
    carries the extra counts;
  * b in the LOWER half vs the UPPER half of log10 R.  A saturating relation must show a SMALLER
    exponent in the upper half; a genuine power law shows the same exponent in both.

======================================================================================
TEST 3 — network-sensitivity diagnostics
======================================================================================
  * b by epoch (2010-2014 / 2015-2019 / 2020-2026), Kruskal-Wallis + Spearman vs date;
  * Spearman of b against the event's mean num_stas (2017+ events);
  * b using only detections with num_stas >= 8 (2017+), paired against the same events' all-
    detection b (Wilcoxon).

======================================================================================
TEST 4 — energy-based scaling (2017+ events)
======================================================================================
  b_E (E vs R7 and vs R1) against b_C for the SAME events (paired Wilcoxon), and b_EperDet
  (energy per detection vs R).  If energy per detection RISES with slip rate (b_EperDet > 0) then
  counts under-represent high-rate tremor and the true exponent is b_C + b_EperDet = b_E.

======================================================================================
TEST 5 — model-side simulation of the smoothing bias, BOTH WAYS
======================================================================================
  The observed exponent depends strongly on which slip-rate variant it is measured against
  (b(C,R7) = 0.40, b(C,R3) = 0.20, b(C,R1) = 0.11).  Only one of two stories can produce that
  ladder, and a simulation decides which:
    (a) the counts really respond to the RAW daily rate, and smoothing the regressor destroys the
        relation  ->  fitting against R7 should UNDER-recover the exponent;
    (b) the counts respond to the SMOOTH rate, and R1 = R7 + differentiation noise, so fitting
        against R1 is a classic errors-in-variables ATTENUATION.
  So the simulation is run for EVERY generating rate, not just R1:
      lam(d) = k * (Rtruth(d) / max Rtruth)^b_true ,  k set so mean(lam) = mean(observed C),
      C_sim(d) ~ Poisson(lam(d)) ,  50 replicates per event, per b_true, per generating rate,
  refitted against log10(R1), log10(R3) and log10(R7) exactly as in Test 1.
  b_true in {1.0, 0.5, 0.4};  Rtruth in {R1, R3, R7}.
  DISCRIMINANT: the LADDER of recovered exponents normalised by the R7 one,
  (b_fitR1 / b_fitR7, b_fitR3 / b_fitR7), is compared with the OBSERVED ladder.  Only the
  generating rate whose simulated ladder matches the observed one is consistent with the data.
  ATTENUATION IDENTITY (an independent check of story (b)): if R1 is R7 plus noise in log space,
  the errors-in-variables prediction is  b(C,R1) = b(C,R7) * lambda1  with lambda1 = the Theil-Sen
  slope of log10 R7 on log10 R1 (the reliability ratio); the same for R3 with lambda3.  Both
  predictions are compared with the measured b(C,R1) and b(C,R3) per event.
  RATE CHARACTER: lag-1 autocorrelation of log10 R for each variant and the Pearson correlation
  of log10 R1 with log10 R7, to show how much of R1 is day-to-day differentiation noise.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s11_exponents.csv     one row per (measure, slip variant, subset): n events, n days,
                                b median / q25 / q75 / bootstrap 90 % CI, pooled OLS and Theil-Sen,
                                pooled r, null mean / sd / p, sign test, # events with b < 1
  results/s11_saturation.csv    per event: U statistics, saturated-day fractions, b(C vs R7) all /
                                no_sat / lower half / upper half, b(C vs U), b(CperU vs R7)
  results/s11_energy_events.csv per 2017+ event: b_C, b_E, b_Dur, b_Mmax, b_EperDet, b_C_ns8
                                against R7 and R1, mean num_stas, n days
  results/s11_smoothing_sim.csv generating rate x b_true x fitting variant: recovered b median /
                                IQR / mean, ratio to truth, ladder ratio to the R7 fit, n; plus the
                                observed ladder and the per-event attenuation identity
  figures/s11_scaling_robustness.png   6 panels (exponent heatmap; C vs U; lower/upper-half
                                exponents; energy vs R; simulation recovery; b by epoch vs null)

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s11_scaling_robustness.py
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
from scipy.stats import binomtest, kruskal, spearmanr, theilslopes, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)

PRE_A = 10               # window: start - 10 d          (S9)
POST_A = 10              # window: end + 10 d            (S9)
NB_RADIUS_KM = 60.0      # neighbourhood / tremor radius (S9)
MIN_FIT_DAYS = 8         # >= 8 usable days for a per-event exponent (S9)
MIN_HALF_DAYS = 6        # >= 6 usable days for a half-subset exponent (S9)
N_SHIFT = 200            # time-shifted null draws per event (S9)
N_BOOT = 2000            # bootstrap resamples of events
N_SIM = 50               # Test 5 Poisson replicates per event per b_true
SLOTS_PER_DAY = 576      # 86400 / 150 s
U_SAT_HI = 500           # "near the cap"
U_SAT_CUT = 400          # days dropped in the `no_sat` subset
CELL_DEG = 0.1           # Uloc cell size
RATE_SMOOTH_D = 7        # S9's smoothing
COS47 = np.cos(np.deg2rad(47.0))
ENERGY_ERA = np.datetime64("2017-01-01")   # energy/duration/num_stas/magnitude start here

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
    """Theil-Sen slope (median of pairwise slopes).  Identical to scipy.stats.theilslopes[0]
    but ~40x faster; verified against scipy at the end of the per-event pass."""
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


def boot_median_ci(v, nboot=N_BOOT, lo=5, hi=95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


# ==================================================================================
# 1. Slip model  (loader copied from scripts/s9_within_event.py)
# ==================================================================================
log("loading slip model (load_solution + invert, ~1 min) ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()  # (742, 13) metres per unit amplitude
slabf = slab_depth_interpolator(fault)

A_raw = np.nan_to_num(
    np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"]
)
CUM_MM = -1e3 * (Ld @ A_raw.T)  # (742, T) mm, thrust positive
log(f"  Ld {Ld.shape}, A_raw {A_raw.shape}, cum {CUM_MM.shape} mm")

DIFF = np.zeros_like(CUM_MM)
DIFF[:, 1:] = np.diff(CUM_MM, axis=1)                       # raw daily rate, signed
RATE7 = np.clip(uniform_filter1d(DIFF, size=7, axis=1, mode="nearest"), 0.0, None)
RATE3 = np.clip(uniform_filter1d(DIFF, size=3, axis=1, mode="nearest"), 0.0, None)
RATE1 = np.clip(DIFF, 0.0, None)
log(f"  RATE7 max {RATE7.max():.3f}  RATE3 max {RATE3.max():.3f}  RATE1 max {RATE1.max():.3f} mm/d")


# ==================================================================================
# 2. Geometry / labels  (loader copied from scripts/s9_within_event.py)
# ==================================================================================
NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]                       # (T, 742)
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
assert NT == A_raw.shape[0], "A_raw / days length mismatch"
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
log(f"  {NT} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles")


# ==================================================================================
# 3. Tremor  (loader copied from scripts/s9_within_event.py, + the 2017+ columns)
# ==================================================================================
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                 usecols=["time", "lat", "lon", "magnitude", "num_stas", "duration", "energy"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_LAT = tr["lat"].values
TR_X, TR_Y = to_km(tr["lon"].values, TR_LAT)
TR_SEG = seg_of_lat(TR_LAT)
TR_ENERGY = tr["energy"].fillna(0.0).values.astype(float)
TR_DUR = tr["duration"].fillna(0.0).values.astype(float)
TR_NSTA = tr["num_stas"].fillna(0.0).values.astype(float)
TR_MAG = tr["magnitude"].values.astype(float)          # NaN before 2017

SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1
TR_DAY = tr["t"].values.astype("datetime64[D]")
TR_OFF = (TR_DAY - SHIFT_LO).astype(int)               # 0 .. SPAN-1
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]

# 2.5-min slot of day, and 0.1-deg cell code
secs = (tr["t"].values - TR_DAY).astype("timedelta64[s]").astype(int)
TR_SLOT = (secs // 150).astype(np.int64)
assert TR_SLOT.max() < SLOTS_PER_DAY
on_grid = float((secs % 150 == 0).mean())
TR_CELLCODE = (np.floor(TR_LAT / CELL_DEG).astype(np.int64) * 100000
               + np.floor(tr["lon"].values / CELL_DEG).astype(np.int64))

TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}; "
    f"{on_grid * 100:.1f}% on the 150-s grid; wrap span {SPAN} d")

# 2017+ sub-catalogue for the energy-era null (its own wrap window)
SHIFT_LO17 = ENERGY_ERA
SPAN17 = int((SHIFT_HI - SHIFT_LO17).astype(int)) + 1
IS17 = TR_DAY >= ENERGY_ERA
log(f"  energy/duration/num_stas/magnitude populated for {IS17.mean() * 100:.1f}% of rows "
    f"(all and only >= {str(ENERGY_ERA)}); 2017+ wrap span {SPAN17} d")

# global eps for the amplitude measures (1st percentile of positive daily values is computed later,
# after the per-event daily series exist; a placeholder is filled in then)


# ==================================================================================
# 4. Events
# ==================================================================================
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  events {len(ev_in)}  {ev_in.seg.value_counts().to_dict()}")


def epoch_of(datestr: str) -> str:
    y = int(datestr[:4])
    return "2010-2014" if y <= 2014 else ("2015-2019" if y <= 2019 else "2020-2026")


def decyear(datestr: str) -> float:
    t = pd.Timestamp(datestr)
    y0 = pd.Timestamp(year=t.year, month=1, day=1)
    y1 = pd.Timestamp(year=t.year + 1, month=1, day=1)
    return t.year + (t - y0).days / (y1 - y0).days


# ==================================================================================
# 5. Machinery
# ==================================================================================
def window_subset(sorted_off, s, w0_off, nW, span):
    """Detections whose SHIFTED day (o + s) mod span lands inside the window.  (S9's helper.)"""
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


def daily_measures(dday, pos, nW, arr, ncell, want17):
    """All tremor measures for one window.

    dday : day-in-window index per detection (0..nW-1)
    pos  : positions of those detections into the event's sorted near-arrays
    arr  : dict of the event's sorted near-arrays (slot, cell, energy, dur, mag, nsta)
    """
    out = {}
    if dday.size == 0:
        z = np.zeros(nW)
        return dict(C=z, U=z.copy(), Uloc=z.copy(), E=z.copy(), Dur=z.copy(),
                    Mmax=np.full(nW, np.nan), C_ns8=z.copy())
    dd = dday.astype(np.int64)
    out["C"] = np.bincount(dd, minlength=nW).astype(float)

    slot = arr["slot"][pos]
    key = dd * SLOTS_PER_DAY + slot
    uk = np.unique(key)
    out["U"] = np.bincount((uk // SLOTS_PER_DAY).astype(np.int64), minlength=nW).astype(float)

    cell = arr["cell"][pos]
    key2 = key * ncell + cell
    uk2 = np.unique(key2)
    out["Uloc"] = np.bincount((uk2 // (SLOTS_PER_DAY * ncell)).astype(np.int64),
                              minlength=nW).astype(float)

    if want17:
        out["E"] = np.bincount(dd, weights=arr["energy"][pos], minlength=nW)
        out["Dur"] = np.bincount(dd, weights=arr["dur"][pos], minlength=nW)
        mm = np.full(nW, -np.inf)
        mg = arr["mag"][pos]
        okm = np.isfinite(mg)
        if okm.any():
            np.maximum.at(mm, dd[okm], mg[okm])
        out["Mmax"] = np.where(np.isfinite(mm), mm, np.nan)
        s8 = arr["nsta"][pos] >= 8
        out["C_ns8"] = (np.bincount(dd[s8], minlength=nW).astype(float) if s8.any()
                        else np.zeros(nW))
    else:
        z = np.zeros(nW)
        out["E"] = z
        out["Dur"] = z.copy()
        out["Mmax"] = np.full(nW, np.nan)
        out["C_ns8"] = z.copy()
    return out


# --- the measure registry -------------------------------------------------------------------
# name -> (needs 2017+ event, transform, eps-rule)
#   transform "log"  : y = log10(value + eps)
#   transform "none" : y = value               (Mmax; already a log quantity)
MEASURES = [
    ("C",       False, "log",  "one"),
    ("U",       False, "log",  "one"),
    ("Uloc",    False, "log",  "one"),
    ("CperU",   False, "log",  "zero"),
    ("E",       True,  "log",  "p1"),
    ("Dur",     True,  "log",  "p1"),
    ("Mmax",    True,  "none", "zero"),
    ("EperDet", True,  "log",  "zero"),
    ("C_ns8",   True,  "log",  "one"),
]
MEAS_NAMES = [m[0] for m in MEASURES]
MEAS_IS17 = {m[0]: m[1] for m in MEASURES}
MEAS_TRANS = {m[0]: m[2] for m in MEASURES}
MEAS_EPSRULE = {m[0]: m[3] for m in MEASURES}
SLIPVARS = ["R7", "R3", "R1", "Rcum"]
SUBSETS = ["all", "no_sat", "R_lower", "R_upper"]


def derive(meas):
    """CperU and EperDet from the primitives (NaN where undefined)."""
    C = meas["C"]
    U = meas["U"]
    with np.errstate(divide="ignore", invalid="ignore"):
        meas["CperU"] = np.where(U > 0, C / np.maximum(U, 1e-12), np.nan)
        meas["EperDet"] = np.where(C > 0, meas["E"] / np.maximum(C, 1e-12), np.nan)
    return meas


def fit_all(meas, slips, U, eps_map, nmin_all=MIN_FIT_DAYS, nmin_half=MIN_HALF_DAYS,
            names=None):
    """Every (measure, slip variant, subset) exponent for one window.  Returns dict[(m,v,s)] = b."""
    res = {}
    names = names or MEAS_NAMES
    for mname in names:
        v = meas[mname]
        trans = MEAS_TRANS[mname]
        eps = eps_map[mname]
        if trans == "log":
            good_m = np.isfinite(v) & (v > 0)
            y_full = np.where(good_m, np.log10(np.where(good_m, v, 1.0) + eps), np.nan)
        else:                                   # Mmax: already a log quantity, keep it as it is
            good_m = np.isfinite(v)
            y_full = np.where(good_m, v, np.nan)
        for vname in SLIPVARS:
            R = slips[vname]
            ok = good_m & (R > 0)
            if ok.sum() < nmin_half:
                for s in SUBSETS:
                    res[(mname, vname, s)] = np.nan
                continue
            lr_full = np.where(R > 0, np.log10(np.where(R > 0, R, 1.0)), np.nan)
            x = lr_full[ok]
            y = y_full[ok]
            res[(mname, vname, "all")] = ts_slope_fast(x, y, nmin_all)
            ns = ok & (U < U_SAT_CUT)
            res[(mname, vname, "no_sat")] = (ts_slope_fast(lr_full[ns], y_full[ns], nmin_all)
                                             if ns.sum() >= nmin_all else np.nan)
            med = np.median(x)
            lo = x <= med
            hi = ~lo
            res[(mname, vname, "R_lower")] = (ts_slope_fast(x[lo], y[lo], nmin_half)
                                              if lo.sum() >= nmin_half else np.nan)
            res[(mname, vname, "R_upper")] = (ts_slope_fast(x[hi], y[hi], nmin_half)
                                              if hi.sum() >= nmin_half else np.nan)
    return res


# ==================================================================================
# 6. Pass 1 — observed per-event daily series (needed to fix the eps values)
# ==================================================================================
log("pass 1: building per-event daily series ...")
EV = []          # per-event dicts of everything the later passes need
POOL_ROWS = []

for _, r in ev_in.iterrows():
    pid = int(r["id"])
    seg = str(r["seg"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    is17 = pd.Timestamp(r["start"]) >= pd.Timestamp("2017-01-01")

    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = np.unique(np.concatenate([np.asarray(h, int) for h in hits])) if len(hits) else np.array([], int)
    nb = np.union1d(nb, np.flatnonzero(F))

    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
            if len(thits) else np.array([], int))

    aw = PAREA[nb][:, None]
    sl = {}
    for vname, RAT in (("R7", RATE7), ("R3", RATE3), ("R1", RATE1)):
        sl[vname] = (aw * RAT[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)
    w0m = max(w0 - 1, 0)
    Pcum = (PAREA[nb][:, None] * CUM_MM[np.ix_(nb, np.arange(w0m, w1 + 1))]).sum(axis=0)
    rc = np.clip(np.diff(Pcum), 0.0, None)
    if rc.size < nW:                      # only if the window starts at day 0 of the grid
        rc = np.concatenate([np.zeros(nW - rc.size), rc])
    sl["Rcum"] = rc

    # sorted near-arrays (full catalogue) and the 2017+ restricted version
    off = TR_OFF[near]
    order = np.argsort(off, kind="stable")
    idx_sorted = near[order]
    arr = dict(slot=TR_SLOT[idx_sorted], cellraw=TR_CELLCODE[idx_sorted],
               energy=TR_ENERGY[idx_sorted], dur=TR_DUR[idx_sorted],
               mag=TR_MAG[idx_sorted], nsta=TR_NSTA[idx_sorted])
    _, arr["cell"] = np.unique(arr["cellraw"], return_inverse=True)
    ncell = int(arr["cell"].max()) + 1 if arr["cell"].size else 1
    soff = off[order]

    # 2017+ restricted arrays with their own wrap (for the energy-era null)
    m17 = TR_DAY[idx_sorted] >= ENERGY_ERA
    arr17 = {k: v[m17] for k, v in arr.items() if k != "cellraw"}
    arr17["cellraw"] = arr["cellraw"][m17]
    soff17 = (TR_DAY[idx_sorted][m17] - SHIFT_LO17).astype(int)
    ncell17 = ncell

    w0_off = w0 - IDX_OF_SHIFT_LO
    sub, dday = window_subset(soff, 0, w0_off, nW, SPAN)
    meas = derive(daily_measures(dday, sub, nW, arr, ncell, is17))

    EV.append(dict(pid=pid, seg=seg, start=str(r["start"]), Mw=float(r["Mw"]),
                   epoch=epoch_of(r["start"]), decyear=decyear(r["start"]),
                   is17=is17, nW=nW, w0_off=w0_off, nb=nb, n_near=int(near.size),
                   soff=soff, arr=arr, ncell=ncell, soff17=soff17, arr17=arr17,
                   ncell17=ncell17, slips=sl, meas=meas,
                   w0_off17=int((np.datetime64(DAYS[w0]) - SHIFT_LO17).astype(int))))

    dfr = dict(id=pid, seg=seg, epoch=epoch_of(r["start"]), is17=is17, day=np.arange(nW))
    for k in ("C", "U", "Uloc", "CperU", "E", "Dur", "Mmax", "EperDet", "C_ns8"):
        dfr[k] = meas[k]
    for vname in SLIPVARS:
        dfr[vname] = sl[vname]
    POOL_ROWS.append(pd.DataFrame(dfr))

POOL = pd.concat(POOL_ROWS, ignore_index=True)
log(f"  {len(EV)} events, {len(POOL)} event-days; "
    f"{sum(e['is17'] for e in EV)} events are 2017+")

# --- fix eps for the amplitude measures: global 1st percentile of positive daily values -------
EPS = {}
for mname, is17, trans, rule in MEASURES:
    if rule == "one":
        EPS[mname] = 1.0
    elif rule == "zero":
        EPS[mname] = 0.0
    else:
        v = POOL.loc[POOL.is17, mname].values
        v = v[np.isfinite(v) & (v > 0)]
        EPS[mname] = float(np.percentile(v, 1)) if v.size else 1.0
log("  eps per measure: " + ", ".join(f"{k}={EPS[k]:.4g}" for k in MEAS_NAMES))

print("\n" + "=" * 108)
print("S11 — SCALING ROBUSTNESS: is C ~ R^0.4 real, or a detection artefact?")
print("=" * 108)
print(f"27 confirmed north/central events; window [start-{PRE_A}, end+{POST_A}]; near tremor within "
      f"{NB_RADIUS_KM:.0f} km of the footprint; NB = triangles within {NB_RADIUS_KM:.0f} km.")
print(f"2.5-min grid: {on_grid * 100:.1f}% of detection times are exact multiples of 150 s "
      f"-> at most {SLOTS_PER_DAY} windows/day.")
print("\nDaily-measure summary over all 27 event-windows (event-days with C >= 1):")
sub = POOL[POOL.C >= 1]
for k in ("C", "U", "Uloc", "CperU"):
    v = sub[k].values
    print(f"   {k:8s} n={np.isfinite(v).sum():5d}  median {np.nanmedian(v):8.2f}  "
          f"p90 {np.nanpercentile(v, 90):8.2f}  max {np.nanmax(v):8.2f}")
sub17 = POOL[(POOL.C >= 1) & POOL.is17]
for k in ("E", "Dur", "Mmax", "EperDet", "C_ns8"):
    v = sub17[k].values
    v = v[np.isfinite(v)]
    print(f"   {k:8s} n={v.size:5d}  median {np.median(v):12.2f}  p90 {np.percentile(v, 90):12.2f}  "
          f"max {v.max():12.2f}   (2017+ events only)")
print("\nSlip-rate variants over the same event-days:")
for vname in SLIPVARS:
    v = POOL[vname].values
    print(f"   {vname:5s} frac>0 {np.mean(v > 0):.3f}  median(>0) {np.median(v[v > 0]):10.1f}  "
          f"max {v.max():10.1f}  log10 range(>0) "
          f"{np.log10(v[v > 0]).max() - np.log10(v[v > 0]).min():.2f} dex  "
          f"IQR(log10) {np.percentile(np.log10(v[v > 0]), 75) - np.percentile(np.log10(v[v > 0]), 25):.2f} dex")


# ==================================================================================
# 7. Pass 2 — observed exponents + shifted-tremor nulls
# ==================================================================================
log("pass 2: observed exponents + 200-draw shifted-tremor nulls per event ...")
OBS = {}                      # (measure, var, subset) -> {pid: b}
NUL = {}                      # (measure, var, subset) -> {pid: array of null b}
for key in [(m, v, s) for m in MEAS_NAMES for v in SLIPVARS for s in SUBSETS]:
    OBS[key] = {}
    NUL[key] = {}

for e in EV:
    names_all = [m for m in MEAS_NAMES if (not MEAS_IS17[m]) or e["is17"]]
    fo = fit_all(e["meas"], e["slips"], e["meas"]["U"], EPS, names=names_all)
    for k, b in fo.items():
        OBS[k][e["pid"]] = b

    # --- nulls -------------------------------------------------------------------------
    names_full = [m for m in names_all if not MEAS_IS17[m]]
    names_e17 = [m for m in names_all if MEAS_IS17[m]]
    store = {k: [] for k in fo}

    shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
    for sh in shifts:
        sub_n, dd_n = window_subset(e["soff"], int(sh), e["w0_off"], e["nW"], SPAN)
        mn = derive(daily_measures(dd_n, sub_n, e["nW"], e["arr"], e["ncell"], False))
        fn = fit_all(mn, e["slips"], mn["U"], EPS, names=names_full)
        for k, b in fn.items():
            store[k].append(b)

    if names_e17 and e["soff17"].size:
        shifts17 = RNG.integers(365, max(366, SPAN17 - 365 + 1), size=N_SHIFT)
        for sh in shifts17:
            sub_n, dd_n = window_subset(e["soff17"], int(sh), e["w0_off17"], e["nW"], SPAN17)
            mn = derive(daily_measures(dd_n, sub_n, e["nW"], e["arr17"], e["ncell17"], True))
            fn = fit_all(mn, e["slips"], mn["U"], EPS, names=names_e17)
            for k, b in fn.items():
                store[k].append(b)

    for k, lst in store.items():
        if lst:
            NUL[k][e["pid"]] = np.asarray(lst, float)

    log(f"  {e['pid']:5d} {e['start']} {e['seg']:7s} nb={len(e['nb']):3d} near={e['n_near']:6d} "
        f"b(C,R7)={fo.get(('C', 'R7', 'all'), np.nan):+6.3f} "
        f"b(U,R7)={fo.get(('U', 'R7', 'all'), np.nan):+6.3f} "
        f"b(C,R1)={fo.get(('C', 'R1', 'all'), np.nan):+6.3f} "
        f"b(E,R7)={fo.get(('E', 'R7', 'all'), np.nan):+6.3f}")

# --- sanity: the fast Theil-Sen must equal scipy, and b(C,R7,all) must equal S9 --------------
e0 = EV[0]
ok0 = (e0["meas"]["C"] >= 1) & (e0["slips"]["R7"] > 0)
x0 = np.log10(e0["slips"]["R7"][ok0])
y0 = np.log10(e0["meas"]["C"][ok0] + 1.0)
log(f"  CHECK Theil-Sen fast {ts_slope_fast(x0, y0):.10f} vs scipy {theilslopes(y0, x0)[0]:.10f}")
try:
    s9 = pd.read_csv(OUT_RES / "s9_rate_scaling.csv")
    s9e = s9[s9.id == e0["pid"]].exponent.values[0]
    log(f"  CHECK b(C vs R7) event {e0['pid']}: S11 {OBS[('C', 'R7', 'all')][e0['pid']]:.6f} "
        f"vs S9 {s9e:.6f}")
except Exception as exc:                                      # pragma: no cover
    log(f"  (S9 cross-check unavailable: {exc})")


# ==================================================================================
# 8. TEST 1 — the exponent matrix
# ==================================================================================
def pooled_fit(df, mname, vname):
    """LSDV within-estimator on the pooled event-days (event dummies = within-event demeaning)."""
    v = df[mname].values
    R = df[vname].values
    ok = ((np.isfinite(v) & (v > 0) & (R > 0)) if MEAS_TRANS[mname] == "log"
          else (np.isfinite(v) & (R > 0)))
    if ok.sum() < 20:
        return np.nan, np.nan, np.nan, int(ok.sum())
    d = df.loc[ok].copy()
    y = (np.log10(d[mname].values + EPS[mname]) if MEAS_TRANS[mname] == "log"
         else d[mname].values)
    x = np.log10(d[vname].values)
    ids = d["id"].values
    dfx = pd.DataFrame(dict(i=ids, x=x, y=y))
    dfx["xc"] = dfx.x - dfx.groupby("i").x.transform("mean")
    dfx["yc"] = dfx.y - dfx.groupby("i").y.transform("mean")
    if len(dfx) < 20 or np.ptp(dfx.xc.values) == 0:
        return np.nan, np.nan, np.nan, len(dfx)
    ols = float(np.polyfit(dfx.xc.values, dfx.yc.values, 1)[0])
    tss = ts_slope_fast(dfx.xc.values, dfx.yc.values, 20)
    rr = float(np.corrcoef(dfx.xc.values, dfx.yc.values)[0, 1])
    return ols, tss, rr, len(dfx)


log("Test 1: assembling the exponent matrix ...")
rows1 = []
for mname in MEAS_NAMES:
    dfm = POOL[POOL.is17] if MEAS_IS17[mname] else POOL
    for vname in SLIPVARS:
        ols, tss, rr, npool = pooled_fit(dfm, mname, vname)
        for sname in SUBSETS:
            key = (mname, vname, sname)
            bs = np.array([b for b in OBS[key].values() if np.isfinite(b)])
            m, q1, q3, n = med_iqr(bs)
            blo, bhi = boot_median_ci(bs)
            nl = [NUL[key][p] for p in NUL[key] if p in OBS[key] and np.isfinite(OBS[key][p])]
            nl = np.concatenate(nl) if nl else np.array([])
            nl = nl[np.isfinite(nl)]
            # per-event one-sided p, then the group sign test on z
            ps, zs = [], []
            for p in OBS[key]:
                if not np.isfinite(OBS[key][p]) or p not in NUL[key]:
                    continue
                d = NUL[key][p]
                d = d[np.isfinite(d)]
                if d.size < 20:
                    continue
                ps.append(float((d >= OBS[key][p]).mean()))
                zs.append((OBS[key][p] - d.mean()) / d.std() if d.std() > 0 else np.nan)
            zs = np.asarray(zs, float)
            zs_f = zs[np.isfinite(zs)]
            npos = int((zs_f > 0).sum())
            sign_p = binomtest(npos, zs_f.size).pvalue if zs_f.size else np.nan
            rows1.append(dict(
                measure=mname, slip_var=vname, subset=sname, n_events=n,
                b_median=m, b_q25=q1, b_q75=q3, b_boot_lo90=blo, b_boot_hi90=bhi,
                n_below_1=int((bs < 1).sum()) if n else 0,
                pooled_ols=ols, pooled_theilsen=tss, pooled_r=rr, n_pooled_days=npool,
                null_mean=float(nl.mean()) if nl.size else np.nan,
                null_sd=float(nl.std()) if nl.size else np.nan,
                z_median=float(np.median(zs_f)) if zs_f.size else np.nan,
                n_events_z=int(zs_f.size), n_z_positive=npos, sign_p=sign_p,
                n_events_p05=int(np.sum(np.asarray(ps) < 0.05)) if ps else 0,
                eps=EPS[mname], transform=MEAS_TRANS[mname]))
EXP = pd.DataFrame(rows1)
EXP.to_csv(OUT_RES / "s11_exponents.csv", index=False)

print("\n" + "=" * 108)
print("TEST 1 — EXPONENT MATRIX   b = Theil-Sen slope of log10(measure+eps) on log10(slip rate)")
print("=" * 108)
print("  (Mmax is NOT log-transformed: it is already a log quantity, so its b is magnitude units "
      "per decade of R.)")
for sname in SUBSETS:
    print(f"\n  --- subset = {sname} " + "-" * 80)
    print(f"    {'measure':9s} {'var':5s}  {'n':>3s}  {'b median [IQR]':>28s}  {'boot90 CI':>18s}  "
          f"{'pooled OLS/TS':>16s}  {'null mean+-sd':>16s}  {'z med':>6s} {'sign':>7s} {'p_sign':>9s}")
    for mname in MEAS_NAMES:
        for vname in SLIPVARS:
            x = EXP[(EXP.measure == mname) & (EXP.slip_var == vname) & (EXP.subset == sname)]
            if x.empty:
                continue
            x = x.iloc[0]
            print(f"    {mname:9s} {vname:5s} {int(x.n_events):3d}  "
                  f"{x.b_median:+8.3f} [{x.b_q25:+7.3f},{x.b_q75:+7.3f}]  "
                  f"[{x.b_boot_lo90:+7.3f},{x.b_boot_hi90:+7.3f}]  "
                  f"{x.pooled_ols:+7.3f}/{x.pooled_theilsen:+7.3f}  "
                  f"{x.null_mean:+7.3f}+-{x.null_sd:5.3f}  {x.z_median:+6.2f} "
                  f"{int(x.n_z_positive):3d}/{int(x.n_events_z):3d} {x.sign_p:9.2e}")

print("\n  KEY ROW (S9's combination, C vs R7, all days):")
k = EXP[(EXP.measure == "C") & (EXP.slip_var == "R7") & (EXP.subset == "all")].iloc[0]
print(f"    b = {k.b_median:+.3f} [{k.b_q25:+.3f}, {k.b_q75:+.3f}] "
      f"(median 90% CI [{k.b_boot_lo90:+.3f}, {k.b_boot_hi90:+.3f}]); "
      f"S9 reported +0.401 [+0.181, +0.478].")


# ==================================================================================
# 9. TEST 2 — saturation diagnostics
# ==================================================================================
print("\n" + "=" * 108)
print("TEST 2 — SATURATION DIAGNOSTICS")
print("=" * 108)
rows2 = []
for e in EV:
    U = e["meas"]["U"]
    C = e["meas"]["C"]
    act = C >= 1
    okR = (C >= 1) & (e["slips"]["R7"] > 0)
    lc = np.log10(C[act] + 1.0)
    lu = np.log10(U[act] + 1.0)
    rows2.append(dict(
        id=e["pid"], start=e["start"], seg=e["seg"], epoch=e["epoch"], decyear=e["decyear"],
        n_days_active=int(act.sum()), n_days_fit=int(okR.sum()),
        U_max=float(U.max()), U_median_active=float(np.median(U[act])) if act.any() else np.nan,
        C_max=float(C.max()), CperU_max=float(np.nanmax(e["meas"]["CperU"])),
        CperU_median=float(np.nanmedian(e["meas"]["CperU"][act])) if act.any() else np.nan,
        frac_U_ge500=float(np.mean(U[act] >= U_SAT_HI)) if act.any() else np.nan,
        frac_U_ge400=float(np.mean(U[act] >= U_SAT_CUT)) if act.any() else np.nan,
        b_C_R7_all=OBS[("C", "R7", "all")].get(e["pid"], np.nan),
        b_C_R7_nosat=OBS[("C", "R7", "no_sat")].get(e["pid"], np.nan),
        b_C_R7_lower=OBS[("C", "R7", "R_lower")].get(e["pid"], np.nan),
        b_C_R7_upper=OBS[("C", "R7", "R_upper")].get(e["pid"], np.nan),
        b_C_R1_lower=OBS[("C", "R1", "R_lower")].get(e["pid"], np.nan),
        b_C_R1_upper=OBS[("C", "R1", "R_upper")].get(e["pid"], np.nan),
        b_U_R7_all=OBS[("U", "R7", "all")].get(e["pid"], np.nan),
        b_U_R7_lower=OBS[("U", "R7", "R_lower")].get(e["pid"], np.nan),
        b_U_R7_upper=OBS[("U", "R7", "R_upper")].get(e["pid"], np.nan),
        b_CperU_R7_all=OBS[("CperU", "R7", "all")].get(e["pid"], np.nan),
        b_C_vs_U=ts_slope_fast(lu, lc, MIN_FIT_DAYS)))
SAT = pd.DataFrame(rows2)
SAT.to_csv(OUT_RES / "s11_saturation.csv", index=False)

allU = POOL.loc[POOL.C >= 1, "U"].values
print(f"  event-days with C >= 1: {allU.size};  U max over all event-days {allU.max():.0f} "
      f"(hard cap {SLOTS_PER_DAY})")
print(f"  fraction with U >= {U_SAT_HI}: {np.mean(allU >= U_SAT_HI):.4f}   "
      f"U >= {U_SAT_CUT}: {np.mean(allU >= U_SAT_CUT):.4f}   "
      f"U >= 300: {np.mean(allU >= 300):.4f}")
print(f"  per-event fraction of days with U >= {U_SAT_HI}: {fmt_mi(SAT.frac_U_ge500)}")
print(f"  per-event fraction of days with U >= {U_SAT_CUT}: {fmt_mi(SAT.frac_U_ge400)}")
print(f"  per-event max U: {fmt_mi(SAT.U_max)};  multiplicity C/U median {fmt_mi(SAT.CperU_median)}"
      f", max {fmt_mi(SAT.CperU_max)}")

print("\n  b(C vs R7): all days vs saturated days removed (U < 400)")
print(f"    all      {fmt_mi(SAT.b_C_R7_all)}")
print(f"    no_sat   {fmt_mi(SAT.b_C_R7_nosat)}")
d = (SAT.b_C_R7_nosat - SAT.b_C_R7_all).dropna()
if d.size >= 6:
    print(f"    paired difference (no_sat - all) {fmt_mi(d)}  Wilcoxon p = "
          f"{wilcoxon(d).pvalue:.4f}   (a positive difference = saturation was suppressing b)")

print("\n  b(C vs U) — does the count keep growing when the occupancy saturates?")
print(f"    {fmt_mi(SAT.b_C_vs_U)}   (b = 1 means C tracks U exactly; b > 1 means multiplicity adds "
      f"counts)")
print(f"    b(C/U vs R7) = {fmt_mi(SAT.b_CperU_R7_all)}   (the multiplicity's own response to slip "
      f"rate; b_C = b_U + b_CperU)")
print(f"    check: median b_U {np.nanmedian(SAT.b_U_R7_all):+.3f} + median b_CperU "
      f"{np.nanmedian(SAT.b_CperU_R7_all):+.3f} = {np.nanmedian(SAT.b_U_R7_all) + np.nanmedian(SAT.b_CperU_R7_all):+.3f} "
      f"vs median b_C {np.nanmedian(SAT.b_C_R7_all):+.3f}")

print("\n  LOWER vs UPPER half of log10 R (a saturating relation must flatten in the upper half):")
for lab, lo_c, hi_c in (("C vs R7", SAT.b_C_R7_lower, SAT.b_C_R7_upper),
                        ("C vs R1", SAT.b_C_R1_lower, SAT.b_C_R1_upper),
                        ("U vs R7", SAT.b_U_R7_lower, SAT.b_U_R7_upper)):
    mm = np.isfinite(lo_c) & np.isfinite(hi_c)
    print(f"    {lab}:  lower {fmt_mi(lo_c[mm])}   upper {fmt_mi(hi_c[mm])}")
    if mm.sum() >= 6:
        dd = (hi_c - lo_c)[mm]
        print(f"              upper - lower {fmt_mi(dd)}  Wilcoxon p = {wilcoxon(dd).pvalue:.4f}  "
              f"sign {int((dd > 0).sum())}/{int(mm.sum())} positive")


# ==================================================================================
# 10. TEST 3 — network sensitivity
# ==================================================================================
print("\n" + "=" * 108)
print("TEST 3 — NETWORK-SENSITIVITY DIAGNOSTICS")
print("=" * 108)
bC = pd.DataFrame(dict(
    id=[e["pid"] for e in EV], seg=[e["seg"] for e in EV], epoch=[e["epoch"] for e in EV],
    decyear=[e["decyear"] for e in EV], is17=[e["is17"] for e in EV],
    b_C_R7=[OBS[("C", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_C_R1=[OBS[("C", "R1", "all")].get(e["pid"], np.nan) for e in EV],
    b_U_R7=[OBS[("U", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_C_ns8_R7=[OBS[("C_ns8", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_E_R7=[OBS[("E", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_E_R1=[OBS[("E", "R1", "all")].get(e["pid"], np.nan) for e in EV],
    b_Dur_R7=[OBS[("Dur", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_Mmax_R7=[OBS[("Mmax", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_EperDet_R7=[OBS[("EperDet", "R7", "all")].get(e["pid"], np.nan) for e in EV],
    b_EperDet_R1=[OBS[("EperDet", "R1", "all")].get(e["pid"], np.nan) for e in EV],
))
mean_ns, C_sum = [], []
for e in EV:
    m = e["meas"]
    C_sum.append(float(m["C"].sum()))
    if e["is17"]:
        sub_n, dd_n = window_subset(e["soff"], 0, e["w0_off"], e["nW"], SPAN)
        ns = e["arr"]["nsta"][sub_n]
        mean_ns.append(float(ns[ns > 0].mean()) if np.any(ns > 0) else np.nan)
    else:
        mean_ns.append(np.nan)
bC["mean_num_stas"] = mean_ns
bC["C_sum"] = C_sum

for col in ("b_C_R7", "b_C_R1", "b_U_R7"):
    print(f"\n  {col} by epoch:")
    for g, s in bC.groupby("epoch"):
        print(f"    {g:10s} {fmt_mi(s[col])}")
    grps = [g[col].dropna().values for _, g in bC.groupby("epoch")]
    grps = [g for g in grps if g.size >= 3]
    if len(grps) >= 2:
        print(f"    Kruskal-Wallis across epochs p = {kruskal(*grps).pvalue:.4f}")
    rho, p = spearmanr(bC.decyear, bC[col], nan_policy="omit")
    print(f"    Spearman vs date rho = {rho:+.3f} p = {p:.4f}")

sub17 = bC[bC.is17]
print(f"\n  2017+ events (n = {len(sub17)}): mean num_stas per event "
      f"{fmt_mi(sub17.mean_num_stas)}")
for col in ("b_C_R7", "b_E_R7", "b_U_R7"):
    rho, p = spearmanr(sub17.mean_num_stas, sub17[col], nan_policy="omit")
    print(f"    Spearman({col}, mean num_stas) rho = {rho:+.3f}  p = {p:.4f}")

print("\n  num_stas >= 8 detections only (2017+):")
print(f"    b_C     (all detections) {fmt_mi(sub17.b_C_R7)}")
print(f"    b_C_ns8 (>= 8 stations)  {fmt_mi(sub17.b_C_ns8_R7)}")
d = (sub17.b_C_ns8_R7 - sub17.b_C_R7).dropna()
if d.size >= 6:
    print(f"    paired difference {fmt_mi(d)}  Wilcoxon p = {wilcoxon(d).pvalue:.4f}  "
          f"sign {int((d > 0).sum())}/{d.size} positive")


# ==================================================================================
# 11. TEST 4 — energy-based scaling
# ==================================================================================
print("\n" + "=" * 108)
print("TEST 4 — ENERGY-BASED SCALING (2017+ events)")
print("=" * 108)
ENE = bC[bC.is17].copy()
ENE["n_days_fit"] = [int(((e["meas"]["C"] >= 1) & (e["slips"]["R7"] > 0)).sum())
                     for e in EV if e["is17"]]
ENE["n_days_E"] = [int(((e["meas"]["E"] > 0) & (e["slips"]["R7"] > 0)).sum())
                   for e in EV if e["is17"]]
ENE["start"] = [e["start"] for e in EV if e["is17"]]
ENE = ENE[["id", "start", "seg", "epoch", "decyear", "n_days_fit", "n_days_E", "mean_num_stas",
           "C_sum", "b_C_R7", "b_C_R1", "b_U_R7", "b_E_R7", "b_E_R1", "b_Dur_R7", "b_Mmax_R7",
           "b_EperDet_R7", "b_EperDet_R1", "b_C_ns8_R7"]]
ENE.to_csv(OUT_RES / "s11_energy_events.csv", index=False)
print(ENE.to_string(index=False, float_format=lambda v: f"{v:9.3f}"))

for lab, a, b in (("E vs C, against R7", ENE.b_E_R7, ENE.b_C_R7),
                  ("E vs C, against R1", ENE.b_E_R1, ENE.b_C_R1)):
    mm = np.isfinite(a) & np.isfinite(b)
    print(f"\n  {lab}:  b_E {fmt_mi(a[mm])}   b_C {fmt_mi(b[mm])}")
    if mm.sum() >= 6:
        dd = (a - b)[mm]
        print(f"      b_E - b_C {fmt_mi(dd)}  Wilcoxon p = {wilcoxon(dd).pvalue:.4f}  "
              f"sign {int((dd > 0).sum())}/{int(mm.sum())} positive")
print(f"\n  energy per detection: b(E/C vs R7) {fmt_mi(ENE.b_EperDet_R7)}   "
      f"b(E/C vs R1) {fmt_mi(ENE.b_EperDet_R1)}")
v = ENE.b_EperDet_R7.dropna()
if v.size >= 6:
    print(f"    Wilcoxon vs 0: p = {wilcoxon(v).pvalue:.4f}; sign {int((v > 0).sum())}/{v.size} "
          f"positive (p = {binomtest(int((v > 0).sum()), v.size).pvalue:.4f})")
print(f"  duration: b(Dur vs R7) {fmt_mi(ENE.b_Dur_R7)}  "
      f"(the duration column takes only 150 s and 300 s -> Dur is ~150*C by construction)")
print(f"  max magnitude: b(Mmax vs R7) {fmt_mi(ENE.b_Mmax_R7)} magnitude units per decade of R")


# ==================================================================================
# 12. TEST 5 — smoothing-bias simulation, BOTH DIRECTIONS
# ==================================================================================
print("\n" + "=" * 108)
print("TEST 5 — MODEL-SIDE SIMULATION OF THE SMOOTHING BIAS (every generating rate)")
print("=" * 108)
print("  C_sim ~ Poisson(k * (Rtruth/max Rtruth)^b_true), k set so mean(lambda) = mean(observed C) "
      "over the event's\n  Rtruth > 0 days; b is then refitted exactly as in Test 1 against "
      "log10(R1), log10(R3) and log10(R7).")
print("  Running it for Rtruth = R1, R3 AND R7 is what makes it a DISCRIMINANT rather than a "
      "one-sided correction.")

FITVARS = ("R1", "R3", "R7")
rows5 = []
SIMTAB = {}                      # (truth, b_true, fitvar) -> median recovered b
for truth in FITVARS:
    for b_true in (1.0, 0.5, 0.4):
        rec = {v: [] for v in FITVARS}
        for e in EV:
            Rt = e["slips"][truth]
            pos = Rt > 0
            if pos.sum() < MIN_FIT_DAYS:
                continue
            Cobs = e["meas"]["C"]
            mean_target = max(Cobs[pos].mean(), 1.0)
            shape = np.where(pos, (Rt / Rt.max()) ** b_true, 0.0)
            k = mean_target / max(shape[pos].mean(), 1e-12)
            lam = k * shape
            for _ in range(N_SIM):
                Cs = RNG.poisson(lam).astype(float)
                for vname in FITVARS:
                    R = e["slips"][vname]
                    ok = (Cs >= 1) & (R > 0)
                    if ok.sum() >= MIN_FIT_DAYS:
                        rec[vname].append(ts_slope_fast(np.log10(R[ok]),
                                                        np.log10(Cs[ok] + 1.0), MIN_FIT_DAYS))
        med7 = np.nanmedian(np.asarray(rec["R7"], float)) if rec["R7"] else np.nan
        for vname in FITVARS:
            v = np.asarray(rec[vname], float)
            v = v[np.isfinite(v)]
            m, q1, q3, n = med_iqr(v)
            SIMTAB[(truth, b_true, vname)] = m
            rows5.append(dict(truth_var=truth, b_true=b_true, fit_against=vname,
                              b_recovered_median=m, b_recovered_q25=q1, b_recovered_q75=q3,
                              b_recovered_mean=float(v.mean()) if v.size else np.nan,
                              ratio_to_truth=m / b_true if n else np.nan,
                              ladder_vs_R7fit=m / med7 if (n and np.isfinite(med7) and med7 != 0)
                              else np.nan,
                              n_fits=n))
SIM = pd.DataFrame(rows5)

print(f"\n    {'truth':>6s} {'b_true':>7s} {'fit vs':>7s} {'b recovered [IQR]':>30s} "
      f"{'b/b_true':>9s} {'b/b(fit R7)':>12s} {'n':>7s}")
for _, x in SIM.iterrows():
    print(f"    {x.truth_var:>6s} {x.b_true:7.2f} {x.fit_against:>7s} "
          f"{x.b_recovered_median:+8.3f} [{x.b_recovered_q25:+7.3f},{x.b_recovered_q75:+7.3f}]"
          f" {x.ratio_to_truth:9.3f} {x.ladder_vs_R7fit:12.3f} {int(x.n_fits):7d}")

# --- the DISCRIMINANT: which generating rate reproduces the observed ladder? -----------------
obs_b = {v: float(np.nanmedian([OBS[("C", v, "all")].get(e["pid"], np.nan) for e in EV]))
         for v in FITVARS}
obs_ladder = {v: obs_b[v] / obs_b["R7"] for v in FITVARS}
print("\n  LADDER DISCRIMINANT — recovered b normalised by the b obtained fitting against R7:")
print(f"    {'generating rate':>18s}  {'b(fit R1)/b(fit R7)':>21s}  {'b(fit R3)/b(fit R7)':>21s}")
print(f"    {'OBSERVED DATA':>18s}  {obs_ladder['R1']:21.3f}  {obs_ladder['R3']:21.3f}"
      f"     <- b(C,R1)={obs_b['R1']:+.3f}, b(C,R3)={obs_b['R3']:+.3f}, b(C,R7)={obs_b['R7']:+.3f}")
lad_rows = []
for truth in FITVARS:
    for b_true in (1.0, 0.5, 0.4):
        l1 = SIM[(SIM.truth_var == truth) & (SIM.b_true == b_true)
                 & (SIM.fit_against == "R1")].ladder_vs_R7fit.values[0]
        l3 = SIM[(SIM.truth_var == truth) & (SIM.b_true == b_true)
                 & (SIM.fit_against == "R3")].ladder_vs_R7fit.values[0]
        print(f"    {truth + ' , b_true=' + format(b_true, '.1f'):>18s}  {l1:21.3f}  {l3:21.3f}")
        lad_rows.append(dict(truth_var=truth, b_true=b_true, fit_against="LADDER",
                             b_recovered_median=np.nan, b_recovered_q25=np.nan,
                             b_recovered_q75=np.nan, b_recovered_mean=np.nan,
                             ratio_to_truth=np.nan, ladder_vs_R7fit=np.nan, n_fits=0,
                             ladder_R1=l1, ladder_R3=l3,
                             obs_ladder_R1=obs_ladder["R1"], obs_ladder_R3=obs_ladder["R3"],
                             ladder_miss=abs(l1 - obs_ladder["R1"]) + abs(l3 - obs_ladder["R3"])))
LAD = pd.DataFrame(lad_rows)
best = LAD.loc[LAD.ladder_miss.idxmin()]
print(f"\n    -> the observed ladder is matched best by generating rate {best.truth_var} "
      f"(b_true = {best.b_true:.1f}; |miss| = {best.ladder_miss:.3f});")
print("       generating from R1 predicts a ladder > 1 (fitting against R1 would give the LARGEST "
      "exponent),\n       which is the opposite of what the data show.")

SIM = pd.concat([SIM, LAD], ignore_index=True)
SIM.to_csv(OUT_RES / "s11_smoothing_sim.csv", index=False)

# --- attenuation identity and rate character -------------------------------------------------
print("\n  ATTENUATION IDENTITY.  If R1 = R7 + noise in log space, errors-in-variables predicts")
print("  b(C,R1) = b(C,R7) * lambda1 with lambda1 = Theil-Sen slope of log10 R7 on log10 R1.")
lam1, lam3, pred1, pred3, ac1, ac3, ac7, cor17 = [], [], [], [], [], [], [], []
for e in EV:
    R1, R3, R7 = e["slips"]["R1"], e["slips"]["R3"], e["slips"]["R7"]
    m = (R1 > 0) & (R7 > 0) & (R3 > 0)
    if m.sum() < MIN_FIT_DAYS:
        lam1.append(np.nan); lam3.append(np.nan); pred1.append(np.nan); pred3.append(np.nan)
        ac1.append(np.nan); ac3.append(np.nan); ac7.append(np.nan); cor17.append(np.nan)
        continue
    x1, x3, x7 = np.log10(R1[m]), np.log10(R3[m]), np.log10(R7[m])
    l1 = ts_slope_fast(x1, x7, MIN_FIT_DAYS)
    l3 = ts_slope_fast(x3, x7, MIN_FIT_DAYS)
    lam1.append(l1)
    lam3.append(l3)
    b7 = OBS[("C", "R7", "all")].get(e["pid"], np.nan)
    pred1.append(b7 * l1)
    pred3.append(b7 * l3)
    cor17.append(float(np.corrcoef(x1, x7)[0, 1]))
    for arr, out in ((x1, ac1), (x3, ac3), (x7, ac7)):
        out.append(float(np.corrcoef(arr[:-1], arr[1:])[0, 1]) if arr.size > 3 else np.nan)
ATT = pd.DataFrame(dict(id=[e["pid"] for e in EV], lambda1=lam1, lambda3=lam3,
                        b_C_R7=[OBS[("C", "R7", "all")].get(e["pid"], np.nan) for e in EV],
                        b_C_R1_pred=pred1, b_C_R3_pred=pred3,
                        b_C_R1_obs=[OBS[("C", "R1", "all")].get(e["pid"], np.nan) for e in EV],
                        b_C_R3_obs=[OBS[("C", "R3", "all")].get(e["pid"], np.nan) for e in EV],
                        ac1=ac1, ac3=ac3, ac7=ac7, cor_logR1_logR7=cor17))
print(f"    reliability ratio lambda1 (log10 R7 on log10 R1) {fmt_mi(ATT.lambda1)}")
print(f"    reliability ratio lambda3 (log10 R7 on log10 R3) {fmt_mi(ATT.lambda3)}")
print(f"    predicted b(C,R1) = b(C,R7)*lambda1 {fmt_mi(ATT.b_C_R1_pred)}   "
      f"observed {fmt_mi(ATT.b_C_R1_obs)}")
print(f"    predicted b(C,R3) = b(C,R7)*lambda3 {fmt_mi(ATT.b_C_R3_pred)}   "
      f"observed {fmt_mi(ATT.b_C_R3_obs)}")
d1 = (ATT.b_C_R1_obs - ATT.b_C_R1_pred).dropna()
d3 = (ATT.b_C_R3_obs - ATT.b_C_R3_pred).dropna()
if d1.size >= 6:
    print(f"    observed - predicted, R1 {fmt_mi(d1)}  Wilcoxon p = {wilcoxon(d1).pvalue:.4f}")
if d3.size >= 6:
    print(f"    observed - predicted, R3 {fmt_mi(d3)}  Wilcoxon p = {wilcoxon(d3).pvalue:.4f}")
print("\n  RATE CHARACTER (how much of R1 is day-to-day differentiation noise):")
print(f"    lag-1 autocorrelation of log10 R:  R1 {fmt_mi(ATT.ac1)}   R3 {fmt_mi(ATT.ac3)}   "
      f"R7 {fmt_mi(ATT.ac7)}")
print(f"    corr(log10 R1, log10 R7) {fmt_mi(ATT.cor_logR1_logR7)}")
print(f"    peak amplitudes over the whole model: R7 {RATE7.max():.2f}, R3 {RATE3.max():.2f}, "
      f"R1 {RATE1.max():.2f} mm/d  (a factor {RATE1.max() / RATE7.max():.1f} between R1 and R7)")
print("    The vbICA amplitudes A_raw come from TRAILING-30-DAY GNSS displacements, so a genuine "
      "daily slip-rate\n    signal cannot have day-to-day swings of that size: R1 is dominated by "
      "differentiation noise, and fitting\n    against it attenuates the exponent.  R7 is the "
      "slip-rate variant the model can actually support.")


# 13. FIGURE
# ==================================================================================
log("figure ...")
fig = plt.figure(figsize=(19, 11.5))
gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.26,
                      left=0.05, right=0.985, top=0.925, bottom=0.06)

# --- A: exponent heatmap (subset = all) -----------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
M = np.full((len(MEAS_NAMES), len(SLIPVARS)), np.nan)
for i, mname in enumerate(MEAS_NAMES):
    for j, vname in enumerate(SLIPVARS):
        x = EXP[(EXP.measure == mname) & (EXP.slip_var == vname) & (EXP.subset == "all")]
        if not x.empty:
            M[i, j] = x.iloc[0].b_median
im = ax.imshow(M, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
ax.set_xticks(range(len(SLIPVARS)))
ax.set_xticklabels(SLIPVARS)
ax.set_yticks(range(len(MEAS_NAMES)))
ax.set_yticklabels([f"{m}{' *' if MEAS_IS17[m] else ''}" for m in MEAS_NAMES])
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isfinite(M[i, j]):
            ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=9,
                    color="k" if abs(M[i, j]) < 0.6 else "w")
plt.colorbar(im, ax=ax, fraction=0.045, label="median b over events")
ax.set_title("A  exponent matrix  b(measure vs slip rate)\n(* = 2017+ events only; Mmax untransformed)",
             fontsize=10)

# --- B: C vs U ------------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
p = POOL[POOL.C >= 1]
sc = ax.scatter(p.U, p.C, c=p.CperU, s=9, cmap="viridis", alpha=0.6, lw=0)
ax.axvline(SLOTS_PER_DAY, color="r", ls="--", lw=1.2, label=f"576 window/day cap")
ax.axvline(U_SAT_CUT, color="orange", ls=":", lw=1.2, label=f"U = {U_SAT_CUT} (no_sat cut)")
uu = np.logspace(0, np.log10(SLOTS_PER_DAY), 50)
ax.plot(uu, uu, "k-", lw=1, label="C = U (no multiplicity)")
bcu = np.nanmedian(SAT.b_C_vs_U)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("U = distinct 2.5-min windows with a detection")
ax.set_ylabel("C = detections")
ax.set_title(f"B  count vs occupancy — median b(C vs U) = {bcu:+.2f}\n"
             f"max U observed {p.U.max():.0f} of {SLOTS_PER_DAY}; "
             f"{np.mean(p.U >= U_SAT_HI) * 100:.2f} % of days >= {U_SAT_HI}", fontsize=10)
ax.legend(fontsize=7, loc="upper left")
plt.colorbar(sc, ax=ax, fraction=0.045, label="C/U multiplicity")

# --- C: lower vs upper half exponents -------------------------------------------------------
ax = fig.add_subplot(gs[0, 2])
mm = np.isfinite(SAT.b_C_R7_lower) & np.isfinite(SAT.b_C_R7_upper)
ax.scatter(SAT.b_C_R7_lower[mm], SAT.b_C_R7_upper[mm], s=42, c="tab:blue",
           edgecolor="k", lw=0.4, label="C vs R7", zorder=3)
mm2 = np.isfinite(SAT.b_U_R7_lower) & np.isfinite(SAT.b_U_R7_upper)
ax.scatter(SAT.b_U_R7_lower[mm2], SAT.b_U_R7_upper[mm2], s=32, c="tab:orange", marker="^",
           edgecolor="k", lw=0.4, label="U vs R7", zorder=3)
lim = [-1.6, 2.4]
ax.plot(lim, lim, "k-", lw=1)
ax.axhline(0, color="grey", lw=0.6)
ax.axvline(0, color="grey", lw=0.6)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("b on the LOWER half of log10 R")
ax.set_ylabel("b on the UPPER half of log10 R")
ax.set_title("C  saturation test: a saturating relation\nwould sit BELOW the 1:1 line", fontsize=10)
ax.legend(fontsize=8)
mC = np.nanmedian(SAT.b_C_R7_upper[mm] - SAT.b_C_R7_lower[mm])
ax.text(0.03, 0.03, f"median(upper - lower)\n  C: {mC:+.2f}\n  U: "
        f"{np.nanmedian(SAT.b_U_R7_upper[mm2] - SAT.b_U_R7_lower[mm2]):+.2f}",
        transform=ax.transAxes, fontsize=8, va="bottom",
        bbox=dict(fc="w", ec="0.7", alpha=0.9))

# --- D: energy vs R (2017+, within-event demeaned) ------------------------------------------
ax = fig.add_subplot(gs[1, 0])
d17 = POOL[(POOL.is17) & (POOL.E > 0) & (POOL.R7 > 0)].copy()
d17["x"] = np.log10(d17.R7)
d17["y"] = np.log10(d17.E + EPS["E"])
d17["xc"] = d17.x - d17.groupby("id").x.transform("mean")
d17["yc"] = d17.y - d17.groupby("id").y.transform("mean")
ax.scatter(d17.xc, d17.yc, s=10, c="tab:red", alpha=0.45, lw=0, label="energy")
dC = POOL[(POOL.is17) & (POOL.C > 0) & (POOL.R7 > 0)].copy()
dC["x"] = np.log10(dC.R7)
dC["y"] = np.log10(dC.C + 1.0)
dC["xc"] = dC.x - dC.groupby("id").x.transform("mean")
dC["yc"] = dC.y - dC.groupby("id").y.transform("mean")
ax.scatter(dC.xc, dC.yc, s=10, c="tab:blue", alpha=0.45, lw=0, label="count")
xx = np.linspace(d17.xc.min(), d17.xc.max(), 20)
bE = EXP[(EXP.measure == "E") & (EXP.slip_var == "R7") & (EXP.subset == "all")].iloc[0]
bCr = EXP[(EXP.measure == "C") & (EXP.slip_var == "R7") & (EXP.subset == "all")].iloc[0]
ax.plot(xx, bE.pooled_theilsen * xx, "r-", lw=2,
        label=f"E pooled b = {bE.pooled_theilsen:+.2f}")
ax.plot(xx, bCr.pooled_theilsen * xx, "b-", lw=2,
        label=f"C pooled b = {bCr.pooled_theilsen:+.2f}")
ax.plot(xx, 1.0 * xx, "k--", lw=1, label="b = 1 (linear)")
ax.set_xlabel("log10 R7  (within-event demeaned)")
ax.set_ylabel("log10 measure  (within-event demeaned)")
ax.set_title(f"D  energy vs count against slip rate (2017+, "
             f"{d17.id.nunique()} events)\nmedian b_E = {bE.b_median:+.2f}, "
             f"b_C = {bCr.b_median:+.2f}", fontsize=10)
ax.legend(fontsize=7.5)

# --- E: simulation LADDER discriminant --------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
cols = {"R1": "tab:green", "R3": "tab:orange", "R7": "tab:purple"}
xs = np.arange(3)
S0 = SIM[SIM.fit_against != "LADDER"]
for truth in ("R1", "R3", "R7"):
    for b_true, ls in ((0.4, "-"), (1.0, "--")):
        yy = np.array([S0[(S0.truth_var == truth) & (S0.b_true == b_true)
                          & (S0.fit_against == v)].ladder_vs_R7fit.values[0] for v in FITVARS])
        ax.plot(xs, yy, ls, marker="o", ms=6, lw=1.7, color=cols[truth], alpha=0.9,
                label=f"counts generated from {truth}, b_true={b_true:.1f}")
obs_y = np.array([np.nanmedian([OBS[("C", v, "all")].get(e["pid"], np.nan) for e in EV])
                  for v in FITVARS])
ax.plot(xs, obs_y / obs_y[2], "ks-", ms=10, lw=2.8, zorder=5,
        label=f"OBSERVED  b(C,R7)={obs_y[2]:+.2f}, R3={obs_y[1]:+.2f}, R1={obs_y[0]:+.2f}")
ax.axhline(1.0, color="grey", lw=0.8, ls=":")
ax.set_yscale("log")
ax.set_xticks(xs)
ax.set_xticklabels([f"fit vs {v}" for v in FITVARS])
ax.set_ylabel("recovered b, normalised by the b from the R7 fit")
ax.set_title("E  which slip-rate variant drives the counts?\ncounts generated from R1 peak at the "
             "R1 fit — the data fall the other way", fontsize=10)
ax.legend(fontsize=6.5, loc="lower right")

# --- F: b by epoch and vs the null ------------------------------------------------------------
ax = fig.add_subplot(gs[1, 2])
allnull = np.concatenate([NUL[("C", "R7", "all")][p] for p in NUL[("C", "R7", "all")]
                          if NUL[("C", "R7", "all")][p].size])
allnull = allnull[np.isfinite(allnull)]
ax.hist(allnull, bins=60, density=True, color="0.75", label=f"shifted-tremor null (n={allnull.size})")
obs_b = SAT.b_C_R7_all.dropna().values
ax.hist(obs_b, bins=12, density=True, alpha=0.65, color="tab:blue",
        label=f"observed b(C,R7) (n={obs_b.size})")
obsU = SAT.b_U_R7_all.dropna().values
ax.hist(obsU, bins=12, density=True, histtype="step", lw=2, color="tab:orange",
        label=f"observed b(U,R7) (n={obsU.size})")
bE_ev = ENE.b_E_R7.dropna().values
ax.hist(bE_ev, bins=8, density=True, histtype="step", lw=2, color="tab:red",
        label=f"observed b(E,R7) (n={bE_ev.size})")
ax.axvline(1.0, color="k", ls="--", lw=1.2, label="b = 1 (linear)")
ax.axvline(np.median(obs_b), color="tab:blue", lw=1.5)
ax.set_xlabel("exponent b")
ax.set_ylabel("density")
ax.set_xlim(-1.5, 2.0)
ax.set_title("F  observed exponents vs the shifted-tremor null\n"
             f"epochs: " + "  ".join(
                 f"{g} {np.nanmedian(s.b_C_R7):+.2f}" for g, s in bC.groupby('epoch')), fontsize=10)
ax.legend(fontsize=7)

fig.suptitle("S11 — robustness of the sub-linear tremor-vs-slip scaling (27 confirmed "
             "north/central events, 2010-2026)", fontsize=13)
fig.savefig(OUT_FIG / "s11_scaling_robustness.png", dpi=140)
plt.close(fig)


# ==================================================================================
# 14. VERDICT
# ==================================================================================
print("\n" + "=" * 108)
print("SUMMARY")
print("=" * 108)


def row(m, v, s="all"):
    x = EXP[(EXP.measure == m) & (EXP.slip_var == v) & (EXP.subset == s)]
    return x.iloc[0] if not x.empty else None


print(f"  {'combination':22s} {'n':>3s}  {'b median [IQR]':>28s}  {'boot90':>18s}  "
      f"{'pooled TS':>9s}  {'sign vs null':>13s}")
for m, v in (("C", "R7"), ("C", "R3"), ("C", "R1"), ("C", "Rcum"),
             ("U", "R7"), ("U", "R1"), ("Uloc", "R7"), ("CperU", "R7"),
             ("E", "R7"), ("E", "R1"), ("Dur", "R7"), ("Mmax", "R7"),
             ("EperDet", "R7"), ("C_ns8", "R7")):
    x = row(m, v)
    if x is None:
        continue
    print(f"  {m + ' vs ' + v:22s} {int(x.n_events):3d}  "
          f"{x.b_median:+8.3f} [{x.b_q25:+7.3f},{x.b_q75:+7.3f}]  "
          f"[{x.b_boot_lo90:+7.3f},{x.b_boot_hi90:+7.3f}]  {x.pooled_theilsen:+9.3f}  "
          f"{int(x.n_z_positive):3d}/{int(x.n_events_z):3d} p={x.sign_p:.1e}")

S0 = SIM[SIM.fit_against != "LADDER"]
r77 = S0[(S0.truth_var == "R7") & (S0.b_true == 0.4) & (S0.fit_against == "R7")].ratio_to_truth.values[0]
r71 = S0[(S0.truth_var == "R7") & (S0.b_true == 0.4) & (S0.fit_against == "R1")].ratio_to_truth.values[0]
r11 = S0[(S0.truth_var == "R1") & (S0.b_true == 0.4) & (S0.fit_against == "R1")].ratio_to_truth.values[0]
r17 = S0[(S0.truth_var == "R1") & (S0.b_true == 0.4) & (S0.fit_against == "R7")].ratio_to_truth.values[0]
bCR7 = row("C", "R7").b_median
bCR1 = row("C", "R1").b_median
print(f"\n  simulation at b_true = 0.4: generated from R7 -> recovered/true = {r77:.3f} (fit R7), "
      f"{r71:.3f} (fit R1)")
print(f"                              generated from R1 -> recovered/true = {r11:.3f} (fit R1), "
      f"{r17:.3f} (fit R7)")
print(f"  observed b(C,R7) = {bCR7:+.3f}, b(C,R1) = {bCR1:+.3f} -> the exponent is LARGEST against "
      f"the smoothed rate,")
print("  which only the 'counts follow the SMOOTH rate, R1 = R7 + differentiation noise' story "
      "produces.")
print(f"  de-biased for the R7 fit itself: b_true = {bCR7:+.3f} / {r77:.3f} = "
      f"{bCR7 / r77:+.3f}")
print(f"  saturation: {np.mean(allU >= U_SAT_HI) * 100:.2f} % of event-days have U >= {U_SAT_HI} "
      f"of {SLOTS_PER_DAY}; upper-half minus lower-half b(C,R7) = "
      f"{np.nanmedian(SAT.b_C_R7_upper - SAT.b_C_R7_lower):+.3f}")
print(f"  energy: b_E(R7) = {row('E', 'R7').b_median:+.3f} vs b_C(R7) on the same events = "
      f"{np.nanmedian(ENE.b_C_R7):+.3f}; energy per detection b = "
      f"{np.nanmedian(ENE.b_EperDet_R7):+.3f}")

log("done.  wrote results/s11_exponents.csv, s11_saturation.csv, s11_energy_events.csv, "
    "s11_smoothing_sim.csv and figures/s11_scaling_robustness.png")
