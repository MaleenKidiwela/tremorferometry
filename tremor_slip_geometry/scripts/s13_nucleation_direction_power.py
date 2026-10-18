#!/usr/bin/env python
"""S13 — tremor-event NUCLEATION latitude, PROPAGATION direction, and the POWER of the S1 bounds.

WHY THIS STEP EXISTS
====================
S10 (literature check, candidate C4) found that the only *verified* published claim of secular
evolution of the Cascadia ETS zone is

    Meng, Liu, Adhikari & Luo (2026), GRL e2025GL120690:
    "a long-term change in the locations where ETS events nucleate and their propagation directions"

Our stationarity results (S1 Q1/Q3 offsets, S4 lag L, S5 phase, S7/S9 edges) are *null* results about
the tremor-vs-slip RELATIVE geometry.  None of them measures the two quantities Meng et al. claim
changed.  A null is only publishable against a named claim if (a) it tests THAT claim's quantities and
(b) it comes with a power analysis.  S13 supplies both:

  PART A  measures nucleation latitude and propagation direction directly from the PNSN tremor
          catalogue (no slip model at all, so it is independent of the vbICA basis), and tests them
          for change over 2010-2026 with two completeness controls.
  PART B  injects synthetic drifts into the S1/S4 per-event quantities and reports the smallest drift
          the stationarity tests would have detected at 80 % power, expressed physically.

======================================================================================
INPUTS (read-only; everything outside tremor_slip_geometry/ is only ever read)
======================================================================================
  catalogs/pnsn_tremor_cascadia_full.csv            PNSN tremor master: time, lat, lon, num_stas
                                                    (num_stas is 0 before 2017-01-01 and real after;
                                                     the `depth` column is a placeholder and IGNORED)
  tremor_slip_geometry/results/s1_events.csv        27 CONFIRMED thrust slip patches (id, start, end,
                                                    seg, Mw, slip_lat, Q1_downdip_km, Q3_depth_km, ...)
  tremor_slip_geometry/results/s4_fronts_events.csv per-event tremor/slip front lead-lag L_days
  tremor_slip_geometry/results/s1_trends.csv        the published S1 Theil-Sen slopes + 95 % CIs
  tremor_slip_geometry/results/s5_tremor_episodes.csv   S5 per-0.5-deg-band episodes (CROSS-CHECK only:
                                                    S13 events are margin-wide space-time objects, S5
                                                    episodes are per-band peaks; we report how they map)

No slip model is loaded: PART A is deliberately model-free.

======================================================================================
PART A — EVENT-LEVEL NUCLEATION AND PROPAGATION
======================================================================================
TREMOR EVENT definition (space-time connected component; deterministic, no eps tuning).
  Grid: latitude bins of DLAT = 0.25 deg over 39.5-51.0 N; days from 2010-01-01 to the catalogue end.
  Smoothed count C[b, d] = number of detections on day d whose latitude falls within +-100 km
        (WIN_KM/2) of the centre of bin b  (a 200-km moving latitude window, per the spec).
  Active[b, d] = C[b, d] >= TH  (TH = 50 detections/day within the 200-km window for the full
        catalogue; see CATALOGUE VARIANTS for how TH is rescaled for a thinned catalogue).
  1-day holes along the time axis are closed (binary closing with a 1x3 structuring element) so a
        single quiet day inside an episode does not split it.
  EVENT = an 8-connected component of Active (8-connectivity lets a front that moves one 0.25-deg bin
        per day stay connected).  A detection belongs to the event that owns the (its own latitude
        bin, its day) cell; detections in inactive cells are unassigned background.
  KEEP events with >= MIN_DET (500) assigned detections AND >= 5 days.
  Rationale for this over ST-DBSCAN on (x_km, y_km, 8 km/d * t): both are space-time clustering, but
  DBSCAN needs an eps that is calibrated to a detection DENSITY, and the PNSN density triples across
  the record and varies 5x along strike -- exactly the non-stationarity under test.  A count
  THRESHOLD per fixed window is the same kind of rule but its sensitivity to the growing network can
  be measured and controlled directly (which is what the completeness controls below do).

PER EVENT
  nucleation_lat_2d      median latitude of the assigned detections in the first 2 days of the event
  nucleation_lat_p10     median latitude of the first 10 % of the assigned detections (time-ordered)
  slope_deg_d            Theil-Sen slope of the DAILY MEDIAN latitude vs day index, over days with
                         >= MIN_DET_DAY_TRACK (10) assigned detections; requires >= 5 such days
  direction              sign(slope_deg_d):  +1 northward, -1 southward
  bilateral              |slope_deg_d| < 0.1 deg/d  (= 11 km/d), "bilateral / stationary"
  extent_deg / _km       lat range of the assigned detections (and a robust p5-p95 version)
  peak_daily             maximum assigned detections on a single day
  dur_d, n_det           duration and size
  seg                    from nucleation_lat_2d:  north >= 46, central 44-46, south < 44

CONTAINMENT CHECK.  Slip start/end in s1_events.csv are 30-d TRAILING threshold crossings, so a slip
  start is late by 0-30 d.  Each of the 27 confirmed slip events is therefore matched against tremor
  events overlapping [start - 30, end] in time; a "spatial" match additionally requires the tremor
  event's latitude range (widened by 0.5 deg) to contain the slip centroid latitude slip_lat.

GROUP TESTS (per segment by NUCLEATION latitude, and margin-wide)
  (i)   Theil-Sen trend of nucleation latitude vs date, 95 % CI, and a two-sided permutation p over
        1000 shuffles of the dates (p = (#{|slope_perm| >= |slope_obs|} + 1) / (nperm + 1));
        plus Spearman; plus the 2010-2017 vs 2018-2026 medians with a Mann-Whitney U test.
  (ii)  fraction of northward-propagating events per epoch (2010-14 / 15-19 / 20-26), with exact
        binomial CIs, a 3x2 chi-square across epochs and a Fisher exact test first vs last epoch.
        Reported for THREE subsets, because the classification bar matters: over ALL events (sign of
        the slope, whatever its size); over events with |slope| >= 0.1 deg/d (= 11 km/d, the spec's
        "not bilateral/stationary" bar -- note this sits ABOVE the observed median |v| of 6 km/d, so
        it keeps only ~25 % of events); and over events with |slope| >= 2 km/d, a bar below the
        typical ETS migration speed that keeps ~88 % of them.
  (iii) nucleation latitude vs event size (Spearman vs n_det and vs peak_daily), margin-wide and
        per segment.  Whether event SIZE itself drifts (log10 n_det, log10 peak_daily vs date) is
        reported alongside, because a size drift would fake (iii).

  MULTIPLICITY AND AUTOCORRELATION.  (i) is run for 4 catalogue variants x 4 segments = 16 tests, so
  a Holm-Bonferroni p is reported for the whole nucleation_lat_2d family.  And because shuffling
  individual dates destroys the within-year clustering of tremor events, every (i) row also carries
  a YEAR-BLOCK permutation p (whole calendar years are the exchangeable unit: block order shuffled,
  blocks re-paired with the date-sorted x).

CATALOGUE VARIANTS (the completeness controls, control (a))
  full        every detection, TH = 50, MIN_DET = 500
  bright      only detections on days with >= 100 detections MARGIN-WIDE ("bright-only": a day that
              only exists as a tremor day because the network grew is dropped wholesale)
  full2017    the full catalogue restricted to 2017-01-01 onward -- the like-for-like time reference
              for the next variant
  sta8_2017   only detections with num_stas >= 8, 2017-01-01 onward (num_stas is not populated
              before 2017)
  For a thinned catalogue the threshold and the size floor are rescaled by
        r = (detections in the variant) / (detections in the full catalogue over the same dates),
  so that the variant's event detector sits at the same place in ITS OWN rate distribution
  (TH = round(50 r), MIN_DET = round(500 r)).  Both the scaled and the unscaled-TH runs are reported
  in the robustness sweep.

CONTROL (c) — SHORT-WINDOW TRENDS (added after the first run, which found a significant northward
  nucleation trend inside the 2017-2026 sub-window that is absent from the full record)
  SLIDING WINDOW.  The 2017+ window is 3433 d = 9.40 yr long.  The same Theil-Sen permutation trend
  test is run inside EVERY 9.40-yr window of the record (30-d steps, >= 20 events required), for the
  margin, the north and the south.  This answers "how often does a 9.4-yr window of this record
  produce a trend this large at p < 0.05, and with which sign?" -- the honest reference for an
  epoch-restricted result, because the date-permutation null does not preserve low-frequency
  structure and is therefore anti-conservative inside a short window.
  SOUTHERN SUB-ZONES.  The south (39.5-44 N) contains two distinct tremor belts.  It is split at
  42 N (Mendocino/Klamath vs southern+central Oregon) and the trend is recomputed WITHIN each zone,
  alongside the epoch-by-epoch mix of the two -- a change in the relative productivity of two fixed
  belts produces an apparent latitude drift with no belt having moved.

CONTROL (b) — the 2018-02 PNSN algorithm change (Wech 2021: +55 % epicentres with the new method)
  STEP TEST.  R(T) = mean daily detection count over [T, T+365) / mean over [T-365, T).
  Observed R(2018-02-01) is compared with R at 200 random dates drawn uniformly from the dates that
  admit a full 365 d on each side; report the observed ratio, the null mean, the null percentile and
  the two-sided empirical p.  The same is done for the number of tremor EVENTS per year.
  BEFORE/AFTER.  Nucleation latitude (Mann-Whitney) and northward fraction (Fisher) are compared
  across 2018-02-01 specifically, per segment and margin-wide.

ROBUSTNESS SWEEP.  The whole PART A pipeline is rerun with TH in {30, 50, 80}, WIN_KM in {150, 200,
  300} and DLAT in {0.25, 0.5}; the margin-wide and north nucleation trend and the northward fraction
  per epoch are recorded for every combination.

======================================================================================
PART B — POWER OF THE STATIONARITY BOUNDS
======================================================================================
For each stationarity quantity q with observed per-event values y_k at dates x_k (decimal years):
  1. Fit Theil-Sen to (x, y); residuals e_k = y_k - (intercept + slope * x_k).  This removes the
     (already ~zero) observed trend so the injected slope is the only trend present.
  2. For each injected slope s and each of NREAL realisations:
        y* = s * (x - 2018.0) + bootstrap resample of e (with replacement, at the observed dates)
     so the observed SCATTER and the observed date sampling are preserved exactly.
  3. Run the identical Theil-Sen trend test (NPERM permutations of the dates, two-sided p on |slope|)
     and record whether p < 0.05.
  4. POWER(s) = fraction of realisations detected.  MDD80 = the smallest |s| with POWER >= 0.80,
     by linear interpolation of the power curve between the bracketing grid points.
  s = 0 is included: its "power" is the empirical false-positive rate and must come out near 0.05.

  QUANTITIES
    Q1_downdip_km   s1_events.csv, n = 27, slopes {0, 0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0} km/yr
    Q3_depth_km     s1_events.csv, n = 27, same grid
    L_days          s4_fronts_events.csv, n = 27, slopes {0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3} d/yr
    nucleation_lat  PART A events (margin-wide and per segment), slopes in deg/yr
                    {0, 0.01, 0.02, 0.04, 0.07, 0.1}; reported also as km/yr (x 111) -- this is the
                    power to detect Meng et al.'s quantity in OUR event catalogue.
    frac_north      power to detect a change in the northward fraction between the first and last
                    epoch, by simulation at the observed epoch sizes (Fisher exact, alpha 0.05).

  PHYSICAL EXPRESSION of the S1/S4 bounds: for every quantity we report
    MDD80 x 16 yr        the smallest total 2010-2026 drift the test could have seen
    CI95 half-width x 16 yr   the published S1/S4 bound restated as a total drift
  and compare both with the mesh triangle size (28 km) and the PNSN epicentre location uncertainty
  (5-10 km).  A bound that is already inside one mesh cell / inside the location error cannot be
  tightened by more data of this kind -- that is the honest statement to make against Meng et al.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s13_tremor_events.csv          one row per kept tremor event (full catalogue variant)
  results/s13_nucleation_trends.csv      trend / epoch-split rows, all catalogue variants and segments
  results/s13_direction_epochs.csv       northward fraction per epoch + chi2 / Fisher rows
  results/s13_completeness_controls.csv  variant summaries, the 2018-02 step test, before/after tests,
                                         the sliding-window control, the southern sub-zone
                                         decomposition, the parameter robustness sweep, the
                                         containment check and the S5 cross-check (one `block` column
                                         names which analysis each row belongs to)
  results/s13_drift_power.csv            power curves, MDD80, and the physical restatement
  figures/s13_nucleation_direction.png   4x3: (a) nucleation vs date by segment with trends;
                                         (b) migration slope vs date; (c) northward fraction per
                                         epoch; (d) the same trend under the four catalogue variants;
                                         (e) daily rate with the 2018-02 marker; (f) step ratio vs its
                                         null; (g) power for Q1/Q3; (h) power for L; (i) power for
                                         nucleation latitude; (j) the sliding-window control;
                                         (k) event latitude-time map; (l) power for the northward
                                         fraction.

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s13_nucleation_direction_power.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.stats import (
    binomtest,
    chi2_contingency,
    fisher_exact,
    mannwhitneyu,
    spearmanr,
    theilslopes,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260910
RNG = np.random.default_rng(SEED)

# ---- PART A parameters -------------------------------------------------------------
LAT0, LAT1, DLAT = 39.5, 51.0, 0.25       # latitude grid for the event detector
WIN_KM = 200.0                            # moving latitude window (full width)
TH_FULL = 50                              # detections/day inside the window to call a cell active
MIN_DET = 500                             # minimum assigned detections per kept event
MIN_DUR = 5                               # minimum days per kept event
MIN_DET_DAY_TRACK = 10                    # detections/day needed for a day to enter the lat track
MIN_TRACK_DAYS = 5                        # days needed for a Theil-Sen migration slope
BILATERAL_DEG_D = 0.10                    # |slope| below this -> "bilateral / stationary" (spec)
DIR2_DEG_D = 2.0 / 111.0                  # secondary 'directional' bar at 2 km/d: the 11 km/d bar
                                          # above sits ABOVE the median observed |v| (6 km/d), so it
                                          # classifies most genuine ETS migrations as 'bilateral'
T0_DATE = pd.Timestamp("2010-01-01")
SPLIT_DATE = pd.Timestamp("2018-01-01")   # the (i) median-split date
ALGO_DATE = pd.Timestamp("2018-02-01")    # the PNSN algorithm change
STA8_DATE = pd.Timestamp("2017-01-01")    # num_stas becomes populated
EPOCHS = [("2010-14", "2010-01-01", "2014-12-31"),
          ("2015-19", "2015-01-01", "2019-12-31"),
          ("2020-26", "2020-01-01", "2026-12-31")]
N_PERM = 1000
N_STEP_NULL = 200

# ---- PART B parameters -------------------------------------------------------------
# (all four can be shrunk with S13_FAST=1 for a smoke test; the published run uses the defaults)
_F = float(os.environ.get("S13_FAST", "0"))
NREAL_SMALL = 50 if _F else 1000     # realisations for the n = 27 quantities
NPERM_SMALL = 50 if _F else 1000
NREAL_BIG = 20 if _F else 300        # realisations for the n ~ 190 nucleation-latitude quantity
NPERM_BIG = 20 if _F else 300
NPERM_AUX = 50 if _F else 300        # permutations inside the sliding-window and parameter sweeps
CHUNK_ELEMS = 3.0e7       # target elements per vectorised chunk
MESH_CELL_KM = 28.0
PNSN_LOC_KM = (5.0, 10.0)
YEARS_SPAN = 16.0

T_START = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T_START:7.1f}s] {msg}", flush=True)


def banner(msg: str) -> None:
    print("\n" + "=" * 100, flush=True)
    print(msg, flush=True)
    print("=" * 100, flush=True)


def seg_of(lat: float) -> str:
    if not np.isfinite(lat):
        return "nan"
    if lat >= 46.0:
        return "north"
    if lat >= 44.0:
        return "central"
    return "south"


def decyear(ts) -> np.ndarray:
    ts = pd.to_datetime(ts)
    y = ts.dt.year.values if hasattr(ts, "dt") else np.array([ts.year])
    start = pd.to_datetime([f"{v}-01-01" for v in y])
    end = pd.to_datetime([f"{v + 1}-01-01" for v in y])
    frac = (ts.values - start.values) / (end.values - start.values)
    return y + frac.astype(float)


# =====================================================================================
# generic Theil-Sen permutation trend test  (the S1/S4 test, reused verbatim)
# =====================================================================================
def theilsen_perm(x, y, nperm=N_PERM, rng=None, min_n=5):
    rng = rng or RNG
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    out = dict(n=n, slope=np.nan, icept=np.nan, lo=np.nan, hi=np.nan,
               p_perm=np.nan, rho=np.nan, p_spear=np.nan)
    if n < min_n or np.ptp(x) == 0 or np.ptp(y) == 0:
        return out
    sl, ic, lo, hi = theilslopes(y, x, 0.95)
    null = np.empty(nperm)
    for k in range(nperm):
        null[k] = theilslopes(rng.permutation(y), x)[0]
    out.update(slope=sl, icept=ic, lo=lo, hi=hi,
               p_perm=(np.sum(np.abs(null) >= abs(sl)) + 1) / (nperm + 1))
    rho, ps = spearmanr(x, y)
    out.update(rho=rho, p_spear=ps)
    return out


def yearblock_perm_p(dates, y, nperm=N_PERM, rng=None):
    """Theil-Sen permutation p with WHOLE CALENDAR YEARS as the exchangeable unit.

    Shuffling individual dates destroys the within-year clustering of tremor events, which makes
    the plain permutation p anti-conservative for a quantity with low-frequency structure.  Here
    the y-values are permuted in whole-year blocks (block order shuffled, blocks concatenated and
    re-paired with the date-sorted x), so year-to-year structure is preserved under the null.
    """
    rng = rng or RNG
    d = pd.DataFrame({"x": decyear(pd.Series(pd.to_datetime(dates))), "y": np.asarray(y, float)})
    d = d.dropna().sort_values("x").reset_index(drop=True)
    if len(d) < 8:
        return np.nan, np.nan
    x = d.x.values
    blocks = [g.y.values for _, g in d.groupby(np.floor(d.x).astype(int), sort=True)]
    sl = theilslopes(d.y.values, x)[0]
    hits = 0
    for _ in range(nperm):
        order = rng.permutation(len(blocks))
        yp = np.concatenate([blocks[i] for i in order])
        if abs(theilslopes(yp, x)[0]) >= abs(sl):
            hits += 1
    return sl, (hits + 1) / (nperm + 1)


# =====================================================================================
# 1. LOAD
# =====================================================================================
banner("S13  nucleation, propagation direction, and the power of the stationarity bounds")
log("loading the PNSN tremor catalogue ...")
TR = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                 usecols=["time", "lat", "lon", "num_stas"])
TR["t"] = pd.to_datetime(TR["time"])
TR = TR.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR["day"] = (TR["t"].dt.floor("D") - T0_DATE).dt.days.astype(int)
TR = TR[TR["day"] >= 0].reset_index(drop=True)
ND = int(TR["day"].max()) + 1
DAY_DATES = T0_DATE + pd.to_timedelta(np.arange(ND), unit="D")
DAILY_ALL = np.bincount(TR["day"].values, minlength=ND).astype(float)
log(f"  {len(TR)} detections {TR.t.min().date()}..{TR.t.max().date()};  {ND} days on the grid; "
    f"mean {DAILY_ALL.mean():.1f} det/d")

EV1 = pd.read_csv(OUT_RES / "s1_events.csv")
FR4 = pd.read_csv(OUT_RES / "s4_fronts_events.csv")
TRD1 = pd.read_csv(OUT_RES / "s1_trends.csv")
EP5 = pd.read_csv(OUT_RES / "s5_tremor_episodes.csv")
log(f"  s1_events {len(EV1)} rows; s4_fronts_events {len(FR4)}; s1_trends {len(TRD1)}; "
    f"s5_tremor_episodes {len(EP5)}")


# =====================================================================================
# 2. PART A — tremor-event detector
# =====================================================================================
def build_events(lat, day, nd, th, min_det, dlat=DLAT, win_km=WIN_KM,
                 lat0=LAT0, lat1=LAT1, min_dur=MIN_DUR):
    """Space-time connected-component tremor events.  Returns (events DataFrame, ev_id per detection)."""
    nb = int(round((lat1 - lat0) / dlat))
    bi = np.clip(np.floor((np.asarray(lat) - lat0) / dlat).astype(int), 0, nb - 1)
    day = np.asarray(day, int)
    H = np.zeros((nb, nd))
    np.add.at(H, (bi, day), 1.0)
    half = max(0, int(round((win_km / 2.0) / 111.0 / dlat)))
    ker = np.ones(2 * half + 1)
    C = ndimage.convolve1d(H, ker, axis=0, mode="constant", cval=0.0)
    A = C >= th
    A = ndimage.binary_closing(A, structure=np.ones((1, 3), bool))   # fill 1-day holes
    lab, ncomp = ndimage.label(A, structure=np.ones((3, 3), bool))   # 8-connectivity
    ev = lab[bi, day]
    rows = []
    for cid in range(1, ncomp + 1):
        m = ev == cid
        n = int(m.sum())
        if n < min_det:
            continue
        d = day[m]
        d0, d1 = int(d.min()), int(d.max())
        if d1 - d0 + 1 < min_dur:
            continue
        rows.append((cid, n, d0, d1))
    keep = {r[0] for r in rows}
    ev_out = np.where(np.isin(ev, list(keep)) if keep else False, ev, 0)
    edf = pd.DataFrame(rows, columns=["cid", "n_det", "d0", "d1"])
    if len(edf):
        edf["dur_d"] = edf.d1 - edf.d0 + 1
        edf = edf.sort_values("d0").reset_index(drop=True)
        edf.insert(0, "ev_id", np.arange(1, len(edf) + 1))
    return edf, ev_out


def event_stats(edf, lat, day, ev_out):
    """Per-event nucleation, migration slope, extent, size."""
    lat = np.asarray(lat, float)
    day = np.asarray(day, int)
    recs = []
    for _, r in edf.iterrows():
        m = ev_out == r.cid
        la, dy = lat[m], day[m]
        o = np.argsort(dy, kind="stable")
        la, dy = la[o], dy[o]
        d0 = int(r.d0)
        first2 = la[dy <= d0 + 1]
        n10 = max(1, int(round(0.10 * len(la))))
        first10 = la[:n10]
        # daily-median latitude track
        udays, inv = np.unique(dy, return_inverse=True)
        cnt = np.bincount(inv)
        med = np.array([np.median(la[inv == k]) for k in range(len(udays))])
        ok = cnt >= MIN_DET_DAY_TRACK
        slope = np.nan
        if ok.sum() >= MIN_TRACK_DAYS and np.ptp(udays[ok]) > 0:
            slope = theilslopes(med[ok], udays[ok].astype(float))[0]
        nuc2 = float(np.median(first2)) if len(first2) else np.nan
        recs.append(dict(
            ev_id=int(r.ev_id), cid=int(r.cid),
            start=str((T0_DATE + pd.Timedelta(days=d0)).date()),
            end=str((T0_DATE + pd.Timedelta(days=int(r.d1))).date()),
            dur_d=int(r.dur_d), n_det=int(r.n_det),
            peak_daily=int(cnt.max()), n_days_det=int(len(udays)), n_days_track=int(ok.sum()),
            nucleation_lat_2d=nuc2,
            nucleation_lat_p10=float(np.median(first10)),
            n_det_first2=int(len(first2)),
            slope_deg_d=slope,
            slope_km_d=slope * 111.0 if np.isfinite(slope) else np.nan,
            direction=int(np.sign(slope)) if np.isfinite(slope) and slope != 0 else 0,
            bilateral=bool(np.isfinite(slope) and abs(slope) < BILATERAL_DEG_D),
            lat_min=float(la.min()), lat_max=float(la.max()),
            extent_deg=float(la.max() - la.min()),
            extent_km=float((la.max() - la.min()) * 111.0),
            extent_p5_95_km=float((np.percentile(la, 95) - np.percentile(la, 5)) * 111.0),
            lat_med=float(np.median(la)),
            seg=seg_of(nuc2),
        ))
    out = pd.DataFrame(recs)
    if len(out):
        out["date"] = pd.to_datetime(out["start"])
        out["decyr"] = decyear(out["date"])
    return out


def make_variant(name, mask, th=None, min_det=None, scale=True):
    """Build events for a catalogue subset.  Threshold rescaled by the subset's rate ratio."""
    sub = TR[mask]
    d0, d1 = int(sub.day.min()), int(sub.day.max())
    ref = int(((TR.day >= d0) & (TR.day <= d1)).sum())
    r = len(sub) / ref if ref else 1.0
    if th is None:
        th = max(5, int(round(TH_FULL * r))) if scale else TH_FULL
    if min_det is None:
        min_det = max(50, int(round(MIN_DET * r))) if scale else MIN_DET
    edf, ev = build_events(sub.lat.values, sub.day.values, ND, th, min_det)
    st = event_stats(edf, sub.lat.values, sub.day.values, ev)
    log(f"  variant {name:12s} n_det={len(sub):7d} r={r:.3f} TH={th:3d} MIN_DET={min_det:4d} "
        f"-> {len(st)} events, {np.isin(ev, edf.cid.values).mean() * 100:.1f}% of its detections assigned")
    return dict(name=name, th=th, min_det=min_det, r=r, n_det=len(sub), events=st,
                assigned_frac=float(np.isin(ev, edf.cid.values).mean()) if len(edf) else 0.0)


banner("PART A.1 — tremor events from the catalogue (no slip model)")
VAR = {}
VAR["full"] = make_variant("full", np.ones(len(TR), bool), th=TH_FULL, min_det=MIN_DET)
bright_days = DAILY_ALL >= 100
VAR["bright"] = make_variant("bright", bright_days[TR.day.values])
m2017 = TR.t >= STA8_DATE
VAR["full2017"] = make_variant("full2017", m2017.values, th=TH_FULL, min_det=MIN_DET)
VAR["sta8_2017"] = make_variant("sta8_2017", (m2017 & (TR.num_stas >= 8)).values)

EVT = VAR["full"]["events"].copy()
print(f"\nFULL variant: {len(EVT)} tremor events, "
      f"{EVT.n_det.sum()} assigned detections ({EVT.n_det.sum() / len(TR) * 100:.1f} % of the catalogue)")
print("  duration  d : " + "  ".join(f"{k}={v:.0f}" for k, v in
                                     EVT.dur_d.describe()[["min", "25%", "50%", "75%", "max"]].items()))
print("  n_det       : " + "  ".join(f"{k}={v:.0f}" for k, v in
                                     EVT.n_det.describe()[["min", "25%", "50%", "75%", "max"]].items()))
print("  extent   km : " + "  ".join(f"{k}={v:.0f}" for k, v in
                                     EVT.extent_km.describe()[["min", "25%", "50%", "75%", "max"]].items()))
print("  by segment (nucleation latitude): " + str(EVT.seg.value_counts().to_dict()))
print("  events per year: " + str(EVT.groupby(EVT.date.dt.year).size().to_dict()))
print(f"  migration slope defined for {EVT.slope_deg_d.notna().sum()}/{len(EVT)} events; "
      f"median |slope| {np.nanmedian(np.abs(EVT.slope_deg_d)) * 111:.1f} km/d; "
      f"bilateral (|slope| < {BILATERAL_DEG_D} deg/d = {BILATERAL_DEG_D * 111:.0f} km/d): "
      f"{EVT.bilateral.mean() * 100:.1f} %")
print(f"  northward fraction (all events, sign of slope): "
      f"{(EVT.slope_deg_d > 0).mean():.3f}; "
      f"directional only: {(EVT.loc[~EVT.bilateral, 'slope_deg_d'] > 0).mean():.3f} "
      f"(n={int((~EVT.bilateral).sum())})")


# ---- containment of the 27 confirmed slip events ------------------------------------
banner("PART A.2 — are the 27 confirmed slip events each contained in one tremor event?")
EV1 = EV1.copy()
EV1["ts"] = pd.to_datetime(EV1["start"])
EV1["te"] = pd.to_datetime(EV1["end"])
EVT_start = pd.to_datetime(EVT["start"])
EVT_end = pd.to_datetime(EVT["end"])
cont_rows = []
for _, e in EV1.iterrows():
    w0, w1 = e.ts - pd.Timedelta(days=30), e.te
    ov = (EVT_start <= w1) & (EVT_end >= w0)
    ovd = np.minimum(EVT_end, w1) - np.maximum(EVT_start, w0)
    ovd = ovd.dt.days + 1
    spat = ov & (EVT.lat_min - 0.5 <= e.slip_lat) & (EVT.lat_max + 0.5 >= e.slip_lat)
    best = EVT.loc[ov].assign(od=ovd[ov]).sort_values("od", ascending=False)
    bestsp = EVT.loc[spat].assign(od=ovd[spat]).sort_values("od", ascending=False)
    cont_rows.append(dict(
        slip_id=int(e.id), start=str(e.ts.date()), end=str(e.te.date()), seg=e.seg,
        Mw=round(float(e.Mw), 2), slip_lat=round(float(e.slip_lat), 2),
        n_tremor_ev_time=int(ov.sum()), n_tremor_ev_timespace=int(spat.sum()),
        best_tremor_ev=int(best.ev_id.iloc[0]) if len(best) else -1,
        best_overlap_d=int(best.od.iloc[0]) if len(best) else 0,
        best_tremor_ev_sp=int(bestsp.ev_id.iloc[0]) if len(bestsp) else -1,
        best_nucleation_lat=float(bestsp.nucleation_lat_2d.iloc[0]) if len(bestsp) else np.nan,
        best_slope_deg_d=float(bestsp.slope_deg_d.iloc[0]) if len(bestsp) else np.nan,
    ))
CONT = pd.DataFrame(cont_rows)
print(CONT.to_string(index=False))
print(f"\n  slip events with >= 1 time-overlapping tremor event : "
      f"{(CONT.n_tremor_ev_time >= 1).sum()}/{len(CONT)}")
print(f"  slip events with EXACTLY ONE time-overlapping tremor event : "
      f"{(CONT.n_tremor_ev_time == 1).sum()}/{len(CONT)}")
print(f"  slip events with >= 1 time+space matching tremor event : "
      f"{(CONT.n_tremor_ev_timespace >= 1).sum()}/{len(CONT)}")
print(f"  slip events with EXACTLY ONE time+space match : "
      f"{(CONT.n_tremor_ev_timespace == 1).sum()}/{len(CONT)}")
print(f"  distribution of the number of time-overlapping tremor events: "
      f"{CONT.n_tremor_ev_time.value_counts().sort_index().to_dict()}")
print(f"  distribution of the number of time+space matches: "
      f"{CONT.n_tremor_ev_timespace.value_counts().sort_index().to_dict()}")
print(f"  distinct tremor events used as the best time+space match: "
      f"{CONT.loc[CONT.best_tremor_ev_sp > 0, 'best_tremor_ev_sp'].nunique()}")

# ---- S5 cross-check ------------------------------------------------------------------
EP5["peak_dt"] = pd.to_datetime(EP5["peak"])
n5_in = 0
for _, ep in EP5.iterrows():
    m = (EVT_start - pd.Timedelta(days=3) <= ep.peak_dt) & (EVT_end + pd.Timedelta(days=3) >= ep.peak_dt) \
        & (EVT.lat_min - 0.4 <= ep.band) & (EVT.lat_max + 0.4 >= ep.band)
    n5_in += int(m.any())
print(f"\n  S5 cross-check: {n5_in}/{len(EP5)} per-band S5 episodes ({n5_in / len(EP5) * 100:.1f} %) "
      f"fall inside an S13 margin-wide tremor event (peak within the event's days +-3 and its lat range +-0.4)")


# =====================================================================================
# 3. PART A group tests
# =====================================================================================
def epoch_of(dt):
    for name, a, b in EPOCHS:
        if pd.Timestamp(a) <= dt <= pd.Timestamp(b):
            return name
    return "other"


def nucleation_block(events, variant, col="nucleation_lat_2d", rng=None):
    rng = rng or RNG
    rows = []
    for grp in ["combined", "north", "central", "south"]:
        e = events if grp == "combined" else events[events.seg == grp]
        e = e.dropna(subset=[col])
        if len(e) < 5:
            rows.append(dict(variant=variant, group=grp, quantity=col, n=len(e)))
            continue
        r = theilsen_perm(e.decyr.values, e[col].values, rng=rng)
        a = e.loc[e.date < SPLIT_DATE, col].values
        b = e.loc[e.date >= SPLIT_DATE, col].values
        if len(a) >= 3 and len(b) >= 3:
            u, pu = mannwhitneyu(a, b, alternative="two-sided")
        else:
            u, pu = np.nan, np.nan
        # size relations
        rho_n, p_n = spearmanr(e[col].values, e.n_det.values)
        rho_pk, p_pk = spearmanr(e[col].values, e.peak_daily.values)
        rows.append(dict(
            variant=variant, group=grp, quantity=col, n=r["n"],
            slope_deg_yr=r["slope"], icept=r["icept"], lo=r["lo"], hi=r["hi"], p_perm=r["p_perm"],
            slope_km_yr=r["slope"] * 111.0 if np.isfinite(r["slope"]) else np.nan,
            drift16_km=r["slope"] * 111.0 * YEARS_SPAN if np.isfinite(r["slope"]) else np.nan,
            ci16_halfwidth_km=(max(abs(r["lo"]), abs(r["hi"])) * 111.0 * YEARS_SPAN)
            if np.isfinite(r["lo"]) else np.nan,
            rho=r["rho"], p_spear=r["p_spear"],
            n_pre2018=len(a), n_post2018=len(b),
            med_pre2018=float(np.median(a)) if len(a) else np.nan,
            med_post2018=float(np.median(b)) if len(b) else np.nan,
            dmed_deg=float(np.median(b) - np.median(a)) if len(a) and len(b) else np.nan,
            mw_U=u, mw_p=pu,
            rho_size_ndet=rho_n, p_size_ndet=p_n,
            rho_size_peak=rho_pk, p_size_peak=p_pk,
        ))
        rows[-1]["p_yearblock"] = yearblock_perm_p(e["date"], e[col].values,
                                                   nperm=NPERM_AUX, rng=rng)[1]
    return rows


def direction_block(events, variant):
    rows = []
    e0 = events.dropna(subset=["slope_deg_d"]).copy()
    e0["epoch"] = e0.date.map(epoch_of)
    for subset, lab in [(e0, "all"),
                        (e0[~e0.bilateral], "directional"),
                        (e0[np.abs(e0.slope_deg_d) >= DIR2_DEG_D], "dir_ge2kmd")]:
        for grp in ["combined", "north", "central", "south"]:
            e = subset if grp == "combined" else subset[subset.seg == grp]
            tab = []
            for name, _, _ in EPOCHS:
                ee = e[e.epoch == name]
                k = int((ee.slope_deg_d > 0).sum())
                n = int(len(ee))
                ci = binomtest(k, n).proportion_ci() if n else None
                rows.append(dict(variant=variant, subset=lab, group=grp, epoch=name, n=n, k_north=k,
                                 frac_north=k / n if n else np.nan,
                                 ci_lo=ci.low if ci else np.nan, ci_hi=ci.high if ci else np.nan,
                                 frac_bilateral=float(ee.bilateral.mean()) if n else np.nan,
                                 median_abs_slope_km_d=float(np.median(np.abs(ee.slope_deg_d)) * 111)
                                 if n else np.nan))
                tab.append((k, n - k))
            tab = np.array(tab)
            ok = tab.sum(axis=1) >= 3
            chi2 = p_chi2 = np.nan
            if ok.sum() >= 2 and tab[ok].sum() > 0 and (tab[ok].sum(axis=0) > 0).all():
                try:
                    chi2, p_chi2 = chi2_contingency(tab[ok])[:2]
                except ValueError:
                    pass
            odds = p_fis = np.nan
            if tab[0].sum() >= 3 and tab[-1].sum() >= 3:
                odds, p_fis = fisher_exact(np.array([tab[0], tab[-1]]))
            rows.append(dict(variant=variant, subset=lab, group=grp, epoch="TEST", n=int(tab.sum()),
                             k_north=int(tab[:, 0].sum()),
                             frac_north=tab[:, 0].sum() / tab.sum() if tab.sum() else np.nan,
                             chi2=chi2, p_chi2=p_chi2, odds_first_last=odds, p_fisher_first_last=p_fis))
    return rows


banner("PART A.3 — (i) nucleation-latitude trends")
NUC_ROWS = []
for v in ["full", "bright", "full2017", "sta8_2017"]:
    for col in (["nucleation_lat_2d", "nucleation_lat_p10"] if v == "full" else ["nucleation_lat_2d"]):
        NUC_ROWS += nucleation_block(VAR[v]["events"], v, col=col)
NUC = pd.DataFrame(NUC_ROWS)
show = ["variant", "group", "quantity", "n", "slope_deg_yr", "lo", "hi", "p_perm", "p_yearblock",
        "slope_km_yr", "rho", "p_spear", "med_pre2018", "med_post2018", "dmed_deg", "mw_p"]
with pd.option_context("display.width", 250, "display.max_columns", 50):
    print(NUC[show].round(4).to_string(index=False))

# Holm-Bonferroni over the whole nucleation_lat_2d family (4 variants x 4 segments)
_fam = NUC[(NUC.quantity == "nucleation_lat_2d") & NUC.p_perm.notna()].copy()
_p = _fam.p_perm.values
_m = len(_p)
_ord = np.argsort(_p)
_holm = np.empty(_m)
_run = 0.0
for _r, _i in enumerate(_ord):
    _run = max(_run, (_m - _r) * _p[_i])
    _holm[_i] = min(1.0, _run)
_fam["p_holm"] = _holm
NUC.loc[_fam.index, "p_holm"] = _holm
print(f"\nHolm-Bonferroni over the {_m} nucleation_lat_2d tests (4 catalogue variants x 4 segments):")
print(_fam.sort_values("p_perm")[["variant", "group", "n", "slope_km_yr", "p_perm", "p_holm",
                                  "p_yearblock"]].round(4).to_string(index=False))

banner("PART A.3 — (ii) northward-propagating fraction per epoch")
DIR_ROWS = []
for v in ["full", "bright", "full2017", "sta8_2017"]:
    DIR_ROWS += direction_block(VAR[v]["events"], v)
DIREP = pd.DataFrame(DIR_ROWS)
with pd.option_context("display.width", 250, "display.max_columns", 50):
    print(DIREP[DIREP.variant.isin(["full", "bright"])]
          [["variant", "subset", "group", "epoch", "n", "k_north", "frac_north", "ci_lo", "ci_hi",
            "frac_bilateral", "median_abs_slope_km_d", "chi2", "p_chi2", "p_fisher_first_last"]]
          .round(4).to_string(index=False))
print("\n(full table for all four catalogue variants in results/s13_direction_epochs.csv)")

banner("PART A.3 — (iii) nucleation latitude vs event size")
with pd.option_context("display.width", 250):
    print(NUC[NUC.variant == "full"][["group", "quantity", "n", "rho_size_ndet", "p_size_ndet",
                                      "rho_size_peak", "p_size_peak"]].round(4).to_string(index=False))
# also: does event SIZE itself drift?  (a completeness diagnostic that would fake (iii))
sz = theilsen_perm(EVT.decyr.values, np.log10(EVT.n_det.values))
sz_pk = theilsen_perm(EVT.decyr.values, np.log10(EVT.peak_daily.values))
print(f"\n  log10(n_det) per event vs date : slope {sz['slope']:+.4f}/yr "
      f"[{sz['lo']:+.4f}, {sz['hi']:+.4f}] p_perm {sz['p_perm']:.3f}")
print(f"  log10(peak_daily) vs date      : slope {sz_pk['slope']:+.4f}/yr "
      f"[{sz_pk['lo']:+.4f}, {sz_pk['hi']:+.4f}] p_perm {sz_pk['p_perm']:.3f}")


# =====================================================================================
# 4. PART A control (b) — the 2018-02 algorithm change
# =====================================================================================
banner("PART A.4 — control (b): the 2018-02 PNSN algorithm change")


def step_ratio(series, t_idx, win=365):
    a = series[max(0, t_idx - win):t_idx]
    b = series[t_idx:t_idx + win]
    if len(a) < win or len(b) < win:
        return np.nan
    return b.mean() / a.mean() if a.mean() > 0 else np.nan


algo_idx = int((ALGO_DATE - T0_DATE).days)
obs_ratio = step_ratio(DAILY_ALL, algo_idx)
cand = np.arange(365, ND - 365)
draw = RNG.choice(cand, size=N_STEP_NULL, replace=False)
null_ratio = np.array([step_ratio(DAILY_ALL, int(i)) for i in draw])
null_ratio = null_ratio[np.isfinite(null_ratio)]
pct = float((null_ratio < obs_ratio).mean() * 100)
p_two = float((np.sum(np.abs(np.log(null_ratio)) >= abs(np.log(obs_ratio))) + 1) / (len(null_ratio) + 1))
print(f"  daily detection rate, 365-d means: mean(after 2018-02-01) / mean(before) = {obs_ratio:.3f}")
print(f"    ({DAILY_ALL[algo_idx:algo_idx + 365].mean():.1f} vs "
      f"{DAILY_ALL[algo_idx - 365:algo_idx].mean():.1f} det/d)")
print(f"  null at {len(null_ratio)} random dates: mean {null_ratio.mean():.3f}, "
      f"sd {null_ratio.std():.3f}, p5 {np.percentile(null_ratio, 5):.3f}, "
      f"p95 {np.percentile(null_ratio, 95):.3f}")
print(f"  observed ratio sits at the {pct:.1f}th percentile of the null; "
      f"two-sided empirical p on |log ratio| = {p_two:.3f}")

# same for the number of tremor EVENTS per 365 d
ev_day = np.zeros(ND)
np.add.at(ev_day, (pd.to_datetime(EVT["start"]) - T0_DATE).dt.days.values, 1.0)
obs_ev = step_ratio(ev_day, algo_idx)
null_ev = np.array([step_ratio(ev_day, int(i)) for i in draw])
null_ev = null_ev[np.isfinite(null_ev)]
p_ev = float((np.sum(np.abs(np.log(np.clip(null_ev, 1e-6, None)))
                     >= abs(np.log(max(obs_ev, 1e-6)))) + 1) / (len(null_ev) + 1))
print(f"  tremor EVENTS per 365 d: after/before = {obs_ev:.3f} "
      f"(null mean {null_ev.mean():.3f}, p={p_ev:.3f})")

# before / after the algorithm change
step_rows = [dict(test="daily_rate_step_2018-02", stat="ratio_after_before_365d", value=obs_ratio,
                  null_mean=null_ratio.mean(), null_sd=null_ratio.std(),
                  null_p5=np.percentile(null_ratio, 5), null_p95=np.percentile(null_ratio, 95),
                  percentile=pct, p_two_sided=p_two, n_null=len(null_ratio)),
             dict(test="event_rate_step_2018-02", stat="ratio_after_before_365d", value=obs_ev,
                  null_mean=null_ev.mean(), null_sd=null_ev.std(),
                  null_p5=np.percentile(null_ev, 5), null_p95=np.percentile(null_ev, 95),
                  percentile=float((null_ev < obs_ev).mean() * 100), p_two_sided=p_ev,
                  n_null=len(null_ev))]

print("\n  nucleation / direction BEFORE vs AFTER 2018-02-01 (full catalogue):")
for grp in ["combined", "north", "central", "south"]:
    e = EVT if grp == "combined" else EVT[EVT.seg == grp]
    a = e.loc[e.date < ALGO_DATE, "nucleation_lat_2d"].dropna().values
    b = e.loc[e.date >= ALGO_DATE, "nucleation_lat_2d"].dropna().values
    if len(a) < 3 or len(b) < 3:
        continue
    u, pu = mannwhitneyu(a, b, alternative="two-sided")
    ea = e[(e.date < ALGO_DATE) & e.slope_deg_d.notna()]
    eb = e[(e.date >= ALGO_DATE) & e.slope_deg_d.notna()]
    ka, kb = int((ea.slope_deg_d > 0).sum()), int((eb.slope_deg_d > 0).sum())
    odds, pf = fisher_exact([[ka, len(ea) - ka], [kb, len(eb) - kb]])
    bila, bilb = ea.bilateral.mean(), eb.bilateral.mean()
    print(f"    {grp:9s} nucleation lat  med {np.median(a):.3f} (n={len(a)}) -> "
          f"{np.median(b):.3f} (n={len(b)})  d={np.median(b) - np.median(a):+.3f} deg "
          f"({(np.median(b) - np.median(a)) * 111:+.0f} km)  MW p={pu:.3f}")
    print(f"    {'':9s} frac north      {ka}/{len(ea)}={ka / max(len(ea), 1):.3f} -> "
          f"{kb}/{len(eb)}={kb / max(len(eb), 1):.3f}  Fisher p={pf:.3f}; "
          f"bilateral {bila:.3f} -> {bilb:.3f}")
    step_rows.append(dict(test="before_after_2018-02", stat=f"nucleation_lat_{grp}",
                          value=float(np.median(b) - np.median(a)), n_before=len(a), n_after=len(b),
                          med_before=float(np.median(a)), med_after=float(np.median(b)),
                          p_two_sided=float(pu)))
    step_rows.append(dict(test="before_after_2018-02", stat=f"frac_north_{grp}",
                          value=float(kb / max(len(eb), 1) - ka / max(len(ea), 1)),
                          n_before=len(ea), n_after=len(eb),
                          med_before=ka / max(len(ea), 1), med_after=kb / max(len(eb), 1),
                          p_two_sided=float(pf)))


# =====================================================================================
# 4b. PART A — is the 2017+ northward nucleation trend unusual for a short window?
# =====================================================================================
banner("PART A.4b — sliding-window control on the nucleation trend, and the southern sub-zones")
WLEN_D = int((TR.t.max() - STA8_DATE).days)          # 3433 d = the length of the 2017+ window
print(f"  the 2017+ window is {WLEN_D} d = {WLEN_D / 365.25:.2f} yr long; sliding it in 30-d steps "
      f"over 2010-2026 gives the null for 'a trend inside a short window'")
slide = []
for grp in ["combined", "north", "south"]:
    e_all = EVT if grp == "combined" else EVT[EVT.seg == grp]
    e_all = e_all.dropna(subset=["nucleation_lat_2d"])
    s0 = T0_DATE
    while s0 + pd.Timedelta(days=WLEN_D) <= TR.t.max():
        e = e_all[(e_all.date >= s0) & (e_all.date <= s0 + pd.Timedelta(days=WLEN_D))]
        if len(e) >= 20:
            r = theilsen_perm(e.decyr.values, e.nucleation_lat_2d.values, nperm=NPERM_AUX)
            slide.append(dict(group=grp, win_start=str(s0.date()), n=r["n"],
                              slope_deg_yr=r["slope"], p_perm=r["p_perm"]))
        s0 = s0 + pd.Timedelta(days=30)
SLIDE = pd.DataFrame(slide)
for grp in ["combined", "north", "south"]:
    d = SLIDE[SLIDE.group == grp]
    if not len(d):
        continue
    ws = pd.to_datetime(d.win_start)
    obs = d.iloc[[int(np.argmin(np.abs((ws - STA8_DATE).dt.days.values)))]]
    obs_sl = float(obs.slope_deg_yr.iloc[0])
    obs_p = float(obs.p_perm.iloc[0])
    obs_ws = obs.win_start.iloc[0]
    pct_sl = float((d.slope_deg_yr <= obs_sl).mean() * 100) if np.isfinite(obs_sl) else np.nan
    print(f"    {grp:9s} {len(d):3d} windows | slope range [{d.slope_deg_yr.min():+.3f}, "
          f"{d.slope_deg_yr.max():+.3f}] deg/yr | {int((d.p_perm < 0.05).sum())} with p<0.05 "
          f"({(d.p_perm < 0.05).mean() * 100:.0f} %), of those "
          f"{int(((d.p_perm < 0.05) & (d.slope_deg_yr > 0)).sum())} northward / "
          f"{int(((d.p_perm < 0.05) & (d.slope_deg_yr < 0)).sum())} southward "
          f"| window nearest 2017-01 ({obs_ws}) slope {obs_sl:+.3f} deg/yr p={obs_p:.3f} "
          f"= {pct_sl:.0f}th percentile of all windows")

# southern sub-zones: is the 2017+ southern drift a Mendocino/Oregon mix change?
print("\n  southern segment split into 39.5-42 N (Mendocino/Klamath) and 42-44 N (S+C Oregon):")
S = EVT[EVT.seg == "south"].dropna(subset=["nucleation_lat_2d"]).copy()
S["zone"] = np.where(S.nucleation_lat_2d < 42.0, "39.5-42", "42-44")
sub_rows = []
for name, a, b in EPOCHS + [("2017+", "2017-01-01", "2026-12-31")]:
    e = S[(S.date >= pd.Timestamp(a)) & (S.date <= pd.Timestamp(b))]
    f = float((e.zone == "42-44").mean()) if len(e) else np.nan
    print(f"    {name:8s} n={len(e):3d}  frac in 42-44 N = {f:.3f}  "
          f"median nucleation lat = {e.nucleation_lat_2d.median():.2f}")
    sub_rows.append(dict(block="south_subzones", stat=f"frac_42_44_{name}", value=f,
                         n_before=len(e), med_before=float(e.nucleation_lat_2d.median())))
for zone in ["39.5-42", "42-44"]:
    e = S[S.zone == zone]
    r_all = theilsen_perm(e.decyr.values, e.nucleation_lat_2d.values)
    e17 = e[e.date >= STA8_DATE]
    r_17 = theilsen_perm(e17.decyr.values, e17.nucleation_lat_2d.values)
    print(f"    WITHIN zone {zone}: full record n={r_all['n']} slope {r_all['slope'] * 111:+.2f} km/yr "
          f"p={r_all['p_perm']:.3f} | 2017+ n={r_17['n']} slope {r_17['slope'] * 111:+.2f} km/yr "
          f"p={r_17['p_perm']:.3f}")
    sub_rows.append(dict(block="south_subzones", stat=f"within_{zone}_slope_kmyr_full",
                         value=r_all["slope"] * 111, n_before=r_all["n"], p_two_sided=r_all["p_perm"]))
    sub_rows.append(dict(block="south_subzones", stat=f"within_{zone}_slope_kmyr_2017plus",
                         value=r_17["slope"] * 111, n_before=r_17["n"], p_two_sided=r_17["p_perm"]))
f1 = float((S.loc[S.date < STA8_DATE, "zone"] == "42-44").mean())
f2 = float((S.loc[S.date >= STA8_DATE, "zone"] == "42-44").mean())
n1s, n2s = int((S.date < STA8_DATE).sum()), int((S.date >= STA8_DATE).sum())
odds_z, p_z = fisher_exact([[int(f1 * n1s), n1s - int(f1 * n1s)], [int(f2 * n2s), n2s - int(f2 * n2s)]])
print(f"    zone MIX pre-2017 {f1:.3f} (n={n1s}) vs 2017+ {f2:.3f} (n={n2s})  Fisher p={p_z:.4f}")
sub_rows.append(dict(block="south_subzones", stat="zone_mix_pre_vs_post2017", value=f2 - f1,
                     n_before=n1s, n_after=n2s, med_before=f1, med_after=f2, p_two_sided=p_z))


# =====================================================================================
# 5. PART A — robustness sweep over the detector parameters
# =====================================================================================
banner("PART A.5 — robustness of the answers to the detector parameters")
sweep = []
for th in ([50] if _F else [30, 50, 80]):
    for wk in ([200.0] if _F else [150.0, 200.0, 300.0]):
        for dl in [0.25, 0.5]:
            edf, ev = build_events(TR.lat.values, TR.day.values, ND, th, MIN_DET, dlat=dl, win_km=wk)
            st = event_stats(edf, TR.lat.values, TR.day.values, ev)
            if len(st) < 20:
                sweep.append(dict(TH=th, WIN_KM=wk, DLAT=dl, n_events=len(st)))
                continue
            rc = theilsen_perm(st.decyr.values, st.nucleation_lat_2d.values, nperm=NPERM_AUX)
            n_ = st[st.seg == "north"]
            rn = theilsen_perm(n_.decyr.values, n_.nucleation_lat_2d.values, nperm=NPERM_AUX)
            st2 = st.dropna(subset=["slope_deg_d"]).copy()
            st2["epoch"] = st2.date.map(epoch_of)
            fr = {nm: float((st2.loc[st2.epoch == nm, "slope_deg_d"] > 0).mean()) for nm, _, _ in EPOCHS}
            tab = np.array([[int((st2.loc[st2.epoch == nm, "slope_deg_d"] > 0).sum()),
                             int((st2.loc[st2.epoch == nm, "slope_deg_d"] <= 0).sum())]
                            for nm, _, _ in EPOCHS])
            try:
                p_chi = chi2_contingency(tab)[1]
            except ValueError:
                p_chi = np.nan
            sweep.append(dict(TH=th, WIN_KM=wk, DLAT=dl, n_events=len(st),
                              med_dur=float(st.dur_d.median()),
                              nuc_slope_all=rc["slope"], nuc_p_all=rc["p_perm"],
                              nuc_slope_north=rn["slope"], nuc_p_north=rn["p_perm"], n_north=rn["n"],
                              frac_north_e1=fr["2010-14"], frac_north_e2=fr["2015-19"],
                              frac_north_e3=fr["2020-26"], p_chi2_epochs=p_chi,
                              frac_bilateral=float(st.bilateral.mean())))
SWEEP = pd.DataFrame(sweep)
with pd.option_context("display.width", 250, "display.max_columns", 50):
    print(SWEEP.round(4).to_string(index=False))
print(f"\n  nucleation trend (margin-wide) across {SWEEP.nuc_slope_all.notna().sum()} parameter sets: "
      f"slope range [{SWEEP.nuc_slope_all.min():+.4f}, {SWEEP.nuc_slope_all.max():+.4f}] deg/yr; "
      f"min p_perm {SWEEP.nuc_p_all.min():.3f}; sets with p<0.05: {int((SWEEP.nuc_p_all < 0.05).sum())}")
print(f"  nucleation trend (north)      : slope range "
      f"[{SWEEP.nuc_slope_north.min():+.4f}, {SWEEP.nuc_slope_north.max():+.4f}] deg/yr; "
      f"min p_perm {SWEEP.nuc_p_north.min():.3f}; sets with p<0.05: {int((SWEEP.nuc_p_north < 0.05).sum())}")
print(f"  epoch chi2 on the northward fraction: min p {SWEEP.p_chi2_epochs.min():.3f}; "
      f"sets with p<0.05: {int((SWEEP.p_chi2_epochs < 0.05).sum())} of {SWEEP.p_chi2_epochs.notna().sum()}")


# =====================================================================================
# 6. PART B — power of the stationarity bounds
# =====================================================================================
banner("PART B — power: what drift could these tests have detected?")
IU = {}


def pair_idx(n):
    if n not in IU:
        IU[n] = np.triu_indices(n, 1)
    return IU[n]


def power_curve(x, y, slopes, nreal, nperm, rng=None, label=""):
    """Injected-drift power for the Theil-Sen permutation trend test."""
    rng = rng or np.random.default_rng(SEED + 7)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    iu, ju = pair_idx(n)
    dx = (x[ju] - x[iu]).astype(np.float32)
    npairs = len(dx)
    sl_obs, ic_obs = theilslopes(y, x)[:2]
    resid = (y - (ic_obs + sl_obs * x)).astype(np.float32)
    chunk = max(1, int(CHUNK_ELEMS / (nperm * npairs)))
    out = []
    for s in slopes:
        det = 0
        done = 0
        while done < nreal:
            R = min(chunk, nreal - done)
            boot = resid[rng.integers(0, n, size=(R, n))]
            ys = (np.float32(s) * (x - 2018.0).astype(np.float32)[None, :] + boot)
            sl = np.median((ys[:, ju] - ys[:, iu]) / dx, axis=-1)              # (R,)
            perm = rng.permuted(np.broadcast_to(np.arange(n), (R, nperm, n)), axis=-1)
            yp = np.take_along_axis(ys[:, None, :], perm, axis=-1)
            nl = np.median((yp[..., ju] - yp[..., iu]) / dx, axis=-1)          # (R, nperm)
            p = (np.sum(np.abs(nl) >= np.abs(sl)[:, None], axis=1) + 1) / (nperm + 1)
            det += int((p < 0.05).sum())
            done += R
        out.append(dict(injected_slope=float(s), power=det / nreal, n_real=nreal))
        log(f"    {label:24s} slope {s:+7.4f} -> power {det / nreal:.3f}")
    df = pd.DataFrame(out)
    df["n"] = n
    df["obs_slope"] = sl_obs
    df["resid_sd"] = float(np.std(resid, ddof=1))
    df["resid_iqr"] = float(np.percentile(resid, 75) - np.percentile(resid, 25))
    return df


def mdd80(df):
    """Smallest |slope| reaching 80 % power, by linear interpolation of the power curve."""
    d = df.sort_values("injected_slope")
    s, p = d.injected_slope.values, d.power.values
    for k in range(1, len(s)):
        if p[k] >= 0.80:
            if p[k] == p[k - 1]:
                return float(s[k])
            return float(s[k - 1] + (0.80 - p[k - 1]) * (s[k] - s[k - 1]) / (p[k] - p[k - 1]))
    return np.nan


POWER = []

# ---- S1 Q1, Q3 ----------------------------------------------------------------------
EV1["decyr"] = decyear(EV1["ts"])
FR4["decyr"] = decyear(pd.to_datetime(FR4["start"]))
grid_km_big = [0.0, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00]
grid_km = [0.0, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 1.00]
grid_lag = [0.0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00]

specs = [
    ("Q1_downdip_km", "km", EV1.decyr.values, EV1.Q1_downdip_km.values, grid_km_big, "S1"),
    ("Q3_depth_km", "km", EV1.decyr.values, EV1.Q3_depth_km.values, grid_km, "S1"),
    ("L_days", "d", FR4.decyr.values, FR4.L_days.values.astype(float), grid_lag, "S4"),
]
for name, unit, x, y, grid, src in specs:
    log(f"  power for {name} (n={np.isfinite(y).sum()}, {NREAL_SMALL} realisations x {NPERM_SMALL} perms)")
    df = power_curve(x, y, grid, NREAL_SMALL, NPERM_SMALL, label=name)
    df["quantity"] = name
    df["unit_per_yr"] = unit
    df["source"] = src
    df["group"] = "confirmed_events"
    POWER.append(df)

# ---- nucleation latitude (Meng et al.'s quantity, in our catalogue) -------------------
grid_deg_big = [0.0, 0.02, 0.04, 0.07, 0.12, 0.18]   # n >= 100 groups (cost ~ n^2)
grid_deg = [0.0, 0.01, 0.02, 0.04, 0.07, 0.10, 0.15, 0.20]
for grp in ["combined", "north", "central", "south"]:
    e = EVT if grp == "combined" else EVT[EVT.seg == grp]
    e = e.dropna(subset=["nucleation_lat_2d"])
    if len(e) < 10:
        continue
    big = len(e) >= 100
    nr = NREAL_BIG if big else (30 if _F else 500)
    npm = NPERM_BIG if big else (30 if _F else 500)
    gd = grid_deg_big if big else grid_deg
    log(f"  power for nucleation_lat_deg [{grp}] (n={len(e)}, {nr} realisations x {npm} perms)")
    df = power_curve(e.decyr.values, e.nucleation_lat_2d.values, gd, nr, npm,
                     label=f"nucleation[{grp}]")
    df["quantity"] = "nucleation_lat_deg"
    df["unit_per_yr"] = "deg"
    df["source"] = "S13"
    df["group"] = grp
    POWER.append(df)

POW = pd.concat(POWER, ignore_index=True)

# ---- summary: MDD80 and the physical restatement -------------------------------------
sum_rows = []
for (q, g), d in POW.groupby(["quantity", "group"], sort=False):
    m = mdd80(d)
    unit = d.unit_per_yr.iloc[0]
    fp = float(d.loc[d.injected_slope == 0, "power"].iloc[0])
    row = dict(quantity=q, group=g, n=int(d.n.iloc[0]), unit_per_yr=unit,
               obs_slope=float(d.obs_slope.iloc[0]), resid_sd=float(d.resid_sd.iloc[0]),
               false_positive_rate_at_slope0=fp,
               mdd80_per_yr=m, mdd80_over_16yr=m * YEARS_SPAN if np.isfinite(m) else np.nan)
    if unit == "deg":
        row["mdd80_km_per_yr"] = m * 111.0
        row["mdd80_km_over_16yr"] = m * 111.0 * YEARS_SPAN
    elif unit == "km":
        row["mdd80_km_per_yr"] = m
        row["mdd80_km_over_16yr"] = m * YEARS_SPAN
    sum_rows.append(row)
PSUM = pd.DataFrame(sum_rows)

# published CI bounds restated physically
ci_rows = []
for _, r in TRD1[(TRD1["set"] == "confirmed") &
                 (TRD1.quantity.isin(["Q1_downdip_km", "Q3_depth_km", "Q4a_updip_km",
                                      "Q4b_downdip_km"]))].iterrows():
    hw = max(abs(r.lo), abs(r.hi))
    ci_rows.append(dict(quantity=r.quantity, group=r.seg, source="S1", unit_per_yr="km",
                        obs_slope=r.slope, ci_lo=r.lo, ci_hi=r.hi, p_perm=r.p_perm,
                        ci_halfwidth_per_yr=hw, ci_halfwidth_over_16yr=hw * YEARS_SPAN,
                        frac_of_mesh_cell=hw * YEARS_SPAN / MESH_CELL_KM,
                        frac_of_pnsn_locerr_5km=hw * YEARS_SPAN / PNSN_LOC_KM[0],
                        frac_of_pnsn_locerr_10km=hw * YEARS_SPAN / PNSN_LOC_KM[1]))
r4 = theilsen_perm(FR4.decyr.values, FR4.L_days.values.astype(float))
ci_rows.append(dict(quantity="L_days", group="combined", source="S4", unit_per_yr="d",
                    obs_slope=r4["slope"], ci_lo=r4["lo"], ci_hi=r4["hi"], p_perm=r4["p_perm"],
                    ci_halfwidth_per_yr=max(abs(r4["lo"]), abs(r4["hi"])),
                    ci_halfwidth_over_16yr=max(abs(r4["lo"]), abs(r4["hi"])) * YEARS_SPAN))
for _, r in NUC[(NUC.variant == "full") & (NUC.quantity == "nucleation_lat_2d")].iterrows():
    if not np.isfinite(r.slope_deg_yr):
        continue
    hw = max(abs(r.lo), abs(r.hi)) * 111.0
    ci_rows.append(dict(quantity="nucleation_lat_km", group=r.group, source="S13", unit_per_yr="km",
                        obs_slope=r.slope_km_yr, ci_lo=r.lo * 111.0, ci_hi=r.hi * 111.0,
                        p_perm=r.p_perm, ci_halfwidth_per_yr=hw,
                        ci_halfwidth_over_16yr=hw * YEARS_SPAN,
                        frac_of_mesh_cell=hw * YEARS_SPAN / MESH_CELL_KM,
                        frac_of_pnsn_locerr_5km=hw * YEARS_SPAN / PNSN_LOC_KM[0],
                        frac_of_pnsn_locerr_10km=hw * YEARS_SPAN / PNSN_LOC_KM[1]))
CIS = pd.DataFrame(ci_rows)

# ---- power to detect a change in the northward fraction -------------------------------
e0 = EVT.dropna(subset=["slope_deg_d"]).copy()
e0["epoch"] = e0.date.map(epoch_of)
n1 = int((e0.epoch == "2010-14").sum())
n3 = int((e0.epoch == "2020-26").sum())
p_base = float((e0.slope_deg_d > 0).mean())
rng2 = np.random.default_rng(SEED + 11)
frac_rows = []
for d in [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40]:
    p1 = np.clip(p_base - d / 2, 0.01, 0.99)
    p3 = np.clip(p_base + d / 2, 0.01, 0.99)
    hit = 0
    NS = 100 if _F else 2000
    k1 = rng2.binomial(n1, p1, NS)
    k3 = rng2.binomial(n3, p3, NS)
    for a, b in zip(k1, k3):
        if fisher_exact([[a, n1 - a], [b, n3 - b]])[1] < 0.05:
            hit += 1
    frac_rows.append(dict(delta_frac_north=d, power=hit / NS, n_epoch1=n1, n_epoch3=n3,
                          p_base=p_base))
FRACP = pd.DataFrame(frac_rows)
mdd_frac = mdd80(FRACP.rename(columns={"delta_frac_north": "injected_slope"}))

banner("PART B — results")
with pd.option_context("display.width", 250, "display.max_columns", 50):
    print("power curves:")
    print(POW[["quantity", "group", "n", "unit_per_yr", "injected_slope", "power"]]
          .round(4).to_string(index=False))
    print("\nsmallest detectable drift at 80 % power (MDD80):")
    print(PSUM.round(4).to_string(index=False))
    print("\npublished bounds restated physically (95 % CI half-width x 16 yr):")
    print(CIS.round(4).to_string(index=False))
    print(f"\npower to detect a CHANGE in the northward fraction between epoch 1 (n={n1}) and "
          f"epoch 3 (n={n3}), base p={p_base:.3f}, Fisher exact alpha 0.05:")
    print(FRACP.round(4).to_string(index=False))
    print(f"  -> smallest detectable change in the northward fraction at 80 % power: "
          f"{mdd_frac:.3f} (absolute, split symmetrically about the base rate)")

print("\nPHYSICAL READING")
for _, r in PSUM.iterrows():
    if r.quantity in ("Q1_downdip_km", "Q3_depth_km"):
        print(f"  {r.quantity:20s}: MDD80 = {r.mdd80_per_yr:.3f} km/yr = "
              f"{r.mdd80_over_16yr:.2f} km over 16 yr = "
              f"{r.mdd80_over_16yr / MESH_CELL_KM:.2f} mesh cells = "
              f"{r.mdd80_over_16yr / PNSN_LOC_KM[0]:.2f}-{r.mdd80_over_16yr / PNSN_LOC_KM[1]:.2f} "
              f"x the PNSN location error")
    elif r.quantity == "L_days":
        print(f"  {r.quantity:20s}: MDD80 = {r.mdd80_per_yr:.3f} d/yr = "
              f"{r.mdd80_over_16yr:.2f} d over 16 yr")
    elif r.quantity == "nucleation_lat_deg":
        print(f"  nucleation_lat [{r.group:8s}]: MDD80 = {r.mdd80_per_yr:.4f} deg/yr = "
              f"{r.mdd80_km_per_yr:.2f} km/yr = {r.mdd80_km_over_16yr:.1f} km over 16 yr "
              f"({r.mdd80_km_over_16yr / MESH_CELL_KM:.1f} mesh cells)")


# =====================================================================================
# 7. FIGURE
# =====================================================================================
banner("figure")
SEGC = {"north": "#1f77b4", "central": "#d62728", "south": "#2ca02c", "combined": "k"}
fig, AX = plt.subplots(4, 3, figsize=(19.5, 19.0))

# (a) nucleation latitude vs date
ax = AX[0, 0]
for grp in ["north", "central", "south"]:
    e = EVT[EVT.seg == grp]
    ax.scatter(e.date, e.nucleation_lat_2d, s=18, c=SEGC[grp], alpha=0.75, label=f"{grp} (n={len(e)})")
    rr = NUC[(NUC.variant == "full") & (NUC.quantity == "nucleation_lat_2d") & (NUC.group == grp)]
    if len(rr) and np.isfinite(rr.slope_deg_yr.iloc[0]):
        xs = np.array([e.decyr.min(), e.decyr.max()])
        ax.plot([e.date.min(), e.date.max()],
                rr.icept.iloc[0] + rr.slope_deg_yr.iloc[0] * xs, c=SEGC[grp], lw=2.2)
        ax.text(0.02, 0.30 - 0.07 * ["north", "central", "south"].index(grp),
                f"{grp}: {rr.slope_deg_yr.iloc[0] * 111:+.2f} km/yr  p={rr.p_perm.iloc[0]:.2f}",
                transform=ax.transAxes, color=SEGC[grp], fontsize=9)
rc = NUC[(NUC.variant == "full") & (NUC.quantity == "nucleation_lat_2d") & (NUC.group == "combined")]
ax.axvline(ALGO_DATE, color="gray", ls=":", lw=1.5)
ax.text(ALGO_DATE, ax.get_ylim()[1], " 2018-02 algo", rotation=90, va="top", fontsize=8, color="gray")
ax.set_ylabel("nucleation latitude (deg N)\nmedian lat of the first 2 days")
ax.set_title(f"(a) event nucleation latitude vs date  ({len(EVT)} tremor events)\n"
             f"margin-wide Theil-Sen {rc.slope_deg_yr.iloc[0] * 111:+.2f} km/yr, p={rc.p_perm.iloc[0]:.2f}",
             fontsize=10)
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3)

# (b) migration slope vs date
ax = AX[0, 1]
for grp in ["north", "central", "south"]:
    e = EVT[(EVT.seg == grp) & EVT.slope_deg_d.notna()]
    ax.scatter(e.date, e.slope_deg_d * 111, s=18, c=SEGC[grp], alpha=0.75)
ax.axhline(0, color="k", lw=1)
ax.axhspan(-BILATERAL_DEG_D * 111, BILATERAL_DEG_D * 111, color="gray", alpha=0.18,
           label=f"|v| < {BILATERAL_DEG_D * 111:.0f} km/d\n(bilateral/stationary, "
                 f"{EVT.bilateral.mean() * 100:.0f} %)")
r_sl = theilsen_perm(EVT.dropna(subset=["slope_deg_d"]).decyr.values,
                     EVT.dropna(subset=["slope_deg_d"]).slope_deg_d.values * 111)
ax.set_ylabel("along-strike migration speed (km/d, + = north)")
ax.set_title(f"(b) propagation direction/speed vs date\n"
             f"Theil-Sen {r_sl['slope']:+.3f} km/d per yr, p={r_sl['p_perm']:.2f}", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
ax.set_ylim(-40, 40)

# (c) northward fraction per epoch
ax = AX[0, 2]
w = 0.2
for i, grp in enumerate(["combined", "north", "central", "south"]):
    d = DIREP[(DIREP.variant == "full") & (DIREP.subset == "all") & (DIREP.group == grp)
              & (DIREP.epoch != "TEST")].sort_values("epoch")
    xs = np.arange(len(d)) + (i - 1.5) * w
    ax.bar(xs, d.frac_north.values, width=w, color=SEGC[grp], alpha=0.8, label=grp)
    ax.errorbar(xs, d.frac_north.values,
                yerr=[d.frac_north.values - d.ci_lo.values, d.ci_hi.values - d.frac_north.values],
                fmt="none", ecolor="k", lw=1, capsize=2)
tst = DIREP[(DIREP.variant == "full") & (DIREP.subset == "all") & (DIREP.group == "combined")
            & (DIREP.epoch == "TEST")]
ax.axhline(0.5, color="k", ls="--", lw=1)
ax.set_xticks(np.arange(3))
ax.set_xticklabels([e[0] for e in EPOCHS])
ax.set_ylabel("fraction propagating northward")
ax.set_title(f"(c) northward fraction per epoch\nmargin-wide chi2 p={tst.p_chi2.iloc[0]:.2f}, "
             f"Fisher(first vs last) p={tst.p_fisher_first_last.iloc[0]:.2f}", fontsize=10)
ax.legend(fontsize=8, ncol=2)
ax.grid(alpha=0.3, axis="y")
ax.set_ylim(0, 1)

# (d) completeness controls: nucleation slope per variant
ax = AX[1, 0]
vars_ = ["full", "bright", "full2017", "sta8_2017"]
for i, grp in enumerate(["combined", "north", "central", "south"]):
    xs, ys, lo, hi = [], [], [], []
    for j, v in enumerate(vars_):
        rr = NUC[(NUC.variant == v) & (NUC.quantity == "nucleation_lat_2d") & (NUC.group == grp)]
        if not len(rr) or not np.isfinite(rr.slope_deg_yr.iloc[0]):
            continue
        xs.append(j + (i - 1.5) * 0.18)
        ys.append(rr.slope_deg_yr.iloc[0] * 111)
        lo.append((rr.slope_deg_yr.iloc[0] - rr.lo.iloc[0]) * 111)
        hi.append((rr.hi.iloc[0] - rr.slope_deg_yr.iloc[0]) * 111)
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=5, color=SEGC[grp], capsize=3, label=grp)
ax.axhline(0, color="k", lw=1)
ax.set_xticks(np.arange(len(vars_)))
ax.set_xticklabels(vars_, fontsize=9)
ax.set_ylabel("nucleation-latitude trend (km/yr)")
ax.set_title("(d) completeness controls: same trend,\nfour catalogue variants (95 % CI)", fontsize=10)
ax.set_ylim(-45, 45)   # the central segment (n = 7-12) has CIs off this scale; see the CSV
ax.legend(fontsize=8, ncol=2)
ax.grid(alpha=0.3, axis="y")

# (e) daily detection rate with the 2018-02 marker
ax = AX[1, 1]
run365 = pd.Series(DAILY_ALL).rolling(365, center=True, min_periods=200).mean()
ax.plot(DAY_DATES, DAILY_ALL, lw=0.3, color="0.7")
ax.plot(DAY_DATES, run365, lw=2, color="C0", label="365-d centred mean")
ax.axvline(ALGO_DATE, color="r", ls="--", lw=1.6, label="2018-02 PNSN algorithm change")
ax.axhline(100, color="g", ls=":", lw=1.2, label="bright-day cut (100/d)")
ax.set_ylabel("detections / day")
ax.set_title(f"(e) PNSN daily detection rate\nafter/before ratio at 2018-02 = {obs_ratio:.2f} "
             f"(null {null_ratio.mean():.2f}+-{null_ratio.std():.2f}, p={p_two:.2f})", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
ax.set_ylim(0, 700)

# (f) step-ratio null
ax = AX[1, 2]
ax.hist(null_ratio, bins=30, color="0.7", edgecolor="k", lw=0.4)
ax.axvline(obs_ratio, color="r", lw=2.2, label=f"2018-02: {obs_ratio:.2f}\n(pct {pct:.0f}, p={p_two:.2f})")
ax.axvline(1.0, color="k", ls=":", lw=1)
ax.set_xlabel("365-d mean rate ratio after/before")
ax.set_ylabel("random dates")
ax.set_title(f"(f) step test at 2018-02 vs {len(null_ratio)} random dates", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")

# (g) power Q1, Q3
ax = AX[2, 0]
for q, c in [("Q1_downdip_km", "C0"), ("Q3_depth_km", "C1")]:
    d = POW[POW.quantity == q].sort_values("injected_slope")
    ax.plot(d.injected_slope, d.power, "o-", color=c, label=q)
    m = PSUM.loc[PSUM.quantity == q, "mdd80_per_yr"].iloc[0]
    if np.isfinite(m):
        ax.axvline(m, color=c, ls=":", lw=1.4)
        ax.text(m, 0.05, f" {m:.2f}", color=c, fontsize=9)
ax.axhline(0.80, color="k", ls="--", lw=1.2)
ax.axhline(0.05, color="0.5", ls=":", lw=1)
ax.set_xlabel("injected drift (km/yr)")
ax.set_ylabel("power (Theil-Sen perm. test, p<0.05)")
ax.set_title(f"(g) S1 offsets: power vs injected drift\n(n=27, {NREAL_SMALL} realisations x "
             f"{NPERM_SMALL} permutations)", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.02)

# (h) power lag
ax = AX[2, 1]
d = POW[POW.quantity == "L_days"].sort_values("injected_slope")
ax.plot(d.injected_slope, d.power, "o-", color="C2", label="S4 lead/lag L")
m = PSUM.loc[PSUM.quantity == "L_days", "mdd80_per_yr"].iloc[0]
if np.isfinite(m):
    ax.axvline(m, color="C2", ls=":", lw=1.4)
    ax.text(m, 0.05, f" {m:.2f} d/yr", color="C2", fontsize=9)
ax.axhline(0.80, color="k", ls="--", lw=1.2)
ax.axhline(0.05, color="0.5", ls=":", lw=1)
ax.set_xlabel("injected drift (d/yr)")
ax.set_ylabel("power")
ax.set_title("(h) S4 tremor-slip lag: power vs injected drift", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.02)

# (i) power nucleation latitude
ax = AX[2, 2]
for grp in ["combined", "north", "central", "south"]:
    d = POW[(POW.quantity == "nucleation_lat_deg") & (POW.group == grp)].sort_values("injected_slope")
    if not len(d):
        continue
    ax.plot(d.injected_slope * 111, d.power, "o-", color=SEGC[grp],
            label=f"{grp} (n={int(d.n.iloc[0])})")
    m = PSUM[(PSUM.quantity == "nucleation_lat_deg") & (PSUM.group == grp)]
    if len(m) and np.isfinite(m.mdd80_km_per_yr.iloc[0]):
        ax.axvline(m.mdd80_km_per_yr.iloc[0], color=SEGC[grp], ls=":", lw=1.2)
ax.axhline(0.80, color="k", ls="--", lw=1.2)
ax.axhline(0.05, color="0.5", ls=":", lw=1)
ax.set_xlabel("injected nucleation-latitude drift (km/yr)")
ax.set_ylabel("power")
ax.set_title("(i) nucleation latitude (Meng et al.'s quantity):\npower vs injected drift", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.02)

# (j) sliding-window nucleation trend
ax = AX[3, 0]
for grp in ["combined", "north", "south"]:
    d = SLIDE[SLIDE.group == grp]
    if not len(d):
        continue
    xs = pd.to_datetime(d.win_start)
    ax.plot(xs, d.slope_deg_yr * 111, "-", color=SEGC[grp], lw=1.6, label=grp)
    sig = d[d.p_perm < 0.05]
    ax.plot(pd.to_datetime(sig.win_start), sig.slope_deg_yr * 111, "o", ms=4, color=SEGC[grp])
ax.axhline(0, color="k", lw=1)
ax.axvline(STA8_DATE, color="r", ls="--", lw=1.4, label="2017-01 window start")
for grp in ["combined", "north", "south"]:
    rr = NUC[(NUC.variant == "full") & (NUC.quantity == "nucleation_lat_2d") & (NUC.group == grp)]
    if len(rr):
        ax.axhline(rr.slope_deg_yr.iloc[0] * 111, color=SEGC[grp], ls=":", lw=1.2)
ax.set_xlabel(f"start of the {WLEN_D / 365.25:.1f}-yr window")
ax.set_ylabel("nucleation trend inside the window (km/yr)")
ax.set_title("(j) sliding-window control: a 9.4-yr window of this\nrecord produces trends of either "
             "sign (filled = p<0.05)", fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# (k) latitude-time event map
ax = AX[3, 1]
for _, e in EVT.iterrows():
    c = "0.55" if not np.isfinite(e.slope_deg_d) else ("C3" if e.slope_deg_d > 0 else "C0")
    ax.plot([pd.Timestamp(e.start), pd.Timestamp(e.end)], [e.lat_min, e.lat_max]
            if (np.isfinite(e.slope_deg_d) and e.slope_deg_d > 0) else [e.lat_max, e.lat_min],
            "-", color=c, lw=1.0, alpha=0.8)
    ax.plot(pd.Timestamp(e.start), e.nucleation_lat_2d, ".", color="k", ms=3)
for yl in [44, 46]:
    ax.axhline(yl, color="k", ls=":", lw=1)
ax.set_ylabel("latitude (deg N)")
ax.set_title("(k) event latitude-time map: black dot = nucleation,\nred = northward, blue = southward",
             fontsize=10)
ax.grid(alpha=0.3)

# (l) power to detect a change in the northward fraction
ax = AX[3, 2]
ax.plot(FRACP.delta_frac_north, FRACP.power, "o-", color="C4")
ax.axhline(0.80, color="k", ls="--", lw=1.2)
ax.axhline(0.05, color="0.5", ls=":", lw=1)
if np.isfinite(mdd_frac):
    ax.axvline(mdd_frac, color="C4", ls=":", lw=1.4)
    ax.text(mdd_frac, 0.05, f" {mdd_frac:.2f}", color="C4", fontsize=9)
ax.set_xlabel("change in the northward fraction, epoch 1 -> epoch 3")
ax.set_ylabel("power (Fisher exact, p<0.05)")
ax.set_title(f"(l) power to see a change in propagation direction\n(n1={n1}, n3={n3}, "
             f"base p={p_base:.2f})", fontsize=10)
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.02)

fig.suptitle("S13 — tremor-event nucleation latitude, propagation direction, and the power of the "
             "Cascadia stationarity bounds (2010-2026)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.975])
figpath = OUT_FIG / "s13_nucleation_direction.png"
fig.savefig(figpath, dpi=130)
plt.close(fig)
log(f"  wrote {figpath}")


# =====================================================================================
# 8. WRITE
# =====================================================================================
banner("writing results")
EVT_out = EVT.copy()
EVT_out["slip_event_ids"] = ""
for _, c in CONT.iterrows():
    if c.best_tremor_ev_sp > 0:
        i = EVT_out.index[EVT_out.ev_id == c.best_tremor_ev_sp]
        EVT_out.loc[i, "slip_event_ids"] = (EVT_out.loc[i, "slip_event_ids"] + f"{c.slip_id};")
EVT_out.drop(columns=["date"]).to_csv(OUT_RES / "s13_tremor_events.csv", index=False)
NUC.to_csv(OUT_RES / "s13_nucleation_trends.csv", index=False)
DIREP.to_csv(OUT_RES / "s13_direction_epochs.csv", index=False)

ctrl = []
for v, d in VAR.items():
    ctrl.append(dict(block="variant_summary", variant=v, n_det=d["n_det"], rate_ratio_r=d["r"],
                     TH=d["th"], MIN_DET=d["min_det"], n_events=len(d["events"]),
                     assigned_frac=d["assigned_frac"],
                     med_dur_d=float(d["events"].dur_d.median()) if len(d["events"]) else np.nan,
                     frac_bilateral=float(d["events"].bilateral.mean()) if len(d["events"]) else np.nan,
                     frac_north=float((d["events"].slope_deg_d > 0).mean()) if len(d["events"]) else np.nan))
for r in step_rows:
    r = dict(r)
    r["block"] = "step_2018"
    ctrl.append(r)
for _, r in SWEEP.iterrows():
    r = r.to_dict()
    r["block"] = "param_sweep"
    ctrl.append(r)
for _, r in SLIDE.iterrows():
    r = r.to_dict()
    r["block"] = "sliding_window"
    ctrl.append(r)
for r in sub_rows:
    r = dict(r)
    ctrl.append(r)
ctrl.append(dict(block="containment", stat="slip_events_with_ge1_time_overlap",
                 value=float((CONT.n_tremor_ev_time >= 1).sum()), n_before=len(CONT)))
ctrl.append(dict(block="containment", stat="slip_events_with_exactly1_time_overlap",
                 value=float((CONT.n_tremor_ev_time == 1).sum()), n_before=len(CONT)))
ctrl.append(dict(block="containment", stat="slip_events_with_ge1_timespace_match",
                 value=float((CONT.n_tremor_ev_timespace >= 1).sum()), n_before=len(CONT)))
ctrl.append(dict(block="containment", stat="slip_events_with_exactly1_timespace_match",
                 value=float((CONT.n_tremor_ev_timespace == 1).sum()), n_before=len(CONT)))
ctrl.append(dict(block="s5_crosscheck", stat="s5_episodes_inside_an_s13_event",
                 value=float(n5_in), n_before=len(EP5)))
ctrl.append(dict(block="event_size_drift", stat="log10_ndet_slope_per_yr", value=sz["slope"],
                 med_before=sz["lo"], med_after=sz["hi"], p_two_sided=sz["p_perm"]))
ctrl.append(dict(block="event_size_drift", stat="log10_peakdaily_slope_per_yr", value=sz_pk["slope"],
                 med_before=sz_pk["lo"], med_after=sz_pk["hi"], p_two_sided=sz_pk["p_perm"]))
CTRL = pd.DataFrame(ctrl)
CTRL = CTRL[["block"] + [c for c in CTRL.columns if c != "block"]]
CONT2 = CONT.copy()
CONT2.insert(0, "block", "containment_per_event")
CTRL = pd.concat([CTRL, CONT2], ignore_index=True)
CTRL.to_csv(OUT_RES / "s13_completeness_controls.csv", index=False)

POW_out = POW.copy()
POW_out["row_type"] = "power_curve"
PSUM_out = PSUM.copy()
PSUM_out["row_type"] = "mdd80"
CIS_out = CIS.copy()
CIS_out["row_type"] = "published_bound_physical"
FRACP_out = FRACP.copy()
FRACP_out["row_type"] = "frac_north_power"
FRACP_out["quantity"] = "frac_north_epoch1_vs_epoch3"
FRACP_out["mdd80_per_yr"] = mdd_frac
PALL = pd.concat([POW_out, PSUM_out, CIS_out, FRACP_out], ignore_index=True)
PALL.to_csv(OUT_RES / "s13_drift_power.csv", index=False)

for f in ["s13_tremor_events.csv", "s13_nucleation_trends.csv", "s13_direction_epochs.csv",
          "s13_completeness_controls.csv", "s13_drift_power.csv"]:
    log(f"  wrote results/{f}  ({sum(1 for _ in open(OUT_RES / f)) - 1} rows)")
log(f"DONE in {time.time() - T_START:.0f} s")
