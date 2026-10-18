#!/usr/bin/env python
"""S9 — WITHIN-EVENT mechanics: tremor depth migration (Q14), tremor-rate vs slip-rate
scaling (Q13), and up-dip/down-dip edges referenced to CONCURRENT margin-wide tremor (Q2 redo).

PLAN.md Round 2, step S9.  Three parts, one script, 27 confirmed north/central events.

Why these definitions (from NOTES.md S1, S1b, S4):
  * The vbICA slip "centroid" is a mixture of 13 FIXED spatial patterns (S4: median 3 components
    carry >= 10 % of an event's slip; v_s is ~1/6 of v_t).  Its SPEED is not physical.  What IS
    meaningful is the daily RATE amplitude and its timing, which S4 showed is synchronous with the
    tremor front to within a few days (|L| <= 5 d in 63 % of events vs 18.5 % null, p = 4e-7).
    -> S9 uses the daily rate as a TIMING/AMPLITUDE quantity only.  The slip-weighted DEPTH track
    is reported alongside the tremor one but is explicitly a mixture-limited quantity.
  * S1b: raw tremor depth quantities are dominated by WHICH DAYS the 5 mm/30 d slip threshold
    happened to select (Spearman(north-wide window depth, Q3) = +0.72).  Every depth quantity here
    is therefore an ANOMALY: Part A subtracts the event's own overall median depth (so the level of
    the window drops out and only the within-event SHAPE survives); Part C subtracts the CONCURRENT
    same-segment margin-wide tremor percentile.
  * S1b also found that the 40-km band around the union of the 21 northern footprints contains
    100.0 % of all lat >= 46N detections.  So in the NORTH "near the footprint" is NOT a local
    filter, and the Part C concurrent-reference subtraction is close to vacuous there (it is a
    real filter in the central segment).  Stated as a caveat, not hidden.

======================================================================================
INPUTS (read-only; everything outside this folder is only READ)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv                       27 confirmed north/central
                                                                   patches (id, start, end, seg, Mw)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz  label (T x 742), days (T),
                                                                   patch_lat/lon/area_km2 (742)
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                           PNSN tremor (time, lat, lon;
                                                                   the depth COLUMN is a placeholder
                                                                   and is ignored)
  /home/jovyan/cascadia-sse-strain                                 mesh + sse_strain:
                                                                   load_fault_ascii / load_solution /
                                                                   component_slip_patterns -> Ld,
                                                                   megathrust.slab_depth_interpolator

  Slip:  cum_mm = -1e3 * (Ld @ A_raw.T)   (742 x T, thrust positive, mm)
         rate   = clip( 7-d CENTRED moving average of diff(cum_mm) , 0, None )   mm/d  (as in S4)
  Depth: Slab2 interface depth at triangle centroids (slip) and at tremor epicentres (tremor);
         NaN outside the mesh hull (6.4 % of detections) -> dropped from depth percentiles only.
  Frame: x = (lon + 124) * 111 * cos(47 deg) km E, y = (lat - 47) * 111 km N (as in S4).
  Segments: north lat >= 46, central 44 <= lat < 46, south lat < 44.

======================================================================================
PART A — tremor depth migration within events (Q14)
======================================================================================
  Per event: window WA = [ts - 10, te + 10] (day indices on the model day grid).
  Footprint F = triangles labelled with the patch id on any day in [ts, te].
  Neighbourhood NB = triangles within 60 km of any F centroid (F always included).
  Near tremor = detections within 60 km of any F centroid (epicentral, km frame).
  Per day d in WA:  n(d) = # near detections;  zmed(d) = median Slab2 depth of those with finite
  depth;  zp10(d) = their 10th percentile depth (the UP-DIP-most tremor of the day).
  DEFINED days = n(d) >= 5.   Phase tau(d) = (d - d_first) / (d_last - d_first) over the first and
  last DEFINED day (event needs >= 8 defined days and d_last > d_first).
  Fifths: bin k = min(int(5 * tau), 4), k = 0..4.  Bin value = MEDIAN over the defined days in the
  bin.  ANOMALY = bin value - (median of zmed over ALL defined days of that event).  Same for zp10.
  SLIP profile on the SAME phase axis: S(d) = sum over NB of area_km2 * rate(d) (mm/d km^2);
  z_s(d) = sum(area * rate * slab_depth) / S(d), defined where S(d) > 0 and d in [d_first, d_last];
  binned and anomalised identically (reference = median z_s over those days).
  STACK = median over events per bin; 90 % CI by 2000 bootstrap resamples OF EVENTS.
  HEADLINE = anomaly(bin 4) - anomaly(bin 0) = the first -> last fifth depth change (+ = deepening).
  NULL: for each event, 200 rigid shifts of the near-tremor DATES by >= 365 d (wrapped over
  2010-01-01..2026-05-27, exactly the S4 control) -> the same 60-km neighbourhood's tremor in a
  random other window of the same length; the whole phase/fifths/anomaly pipeline is rerun per draw
  (phase re-derived from that draw's own defined days), then stacked over events.  The draw INDEX
  is SHARED across events, so draw j is one coherent null realisation of the whole 27-event stack;
  a draw is kept only if >= 15 events produced a profile for it.  Null band = 5th-95th percentile
  of the stacked median per bin; p = fraction of kept draws whose stacked |first->last change| >=
  the observed |change|.  A second, per-event reading of the same 200 draws gives each event a
  z = (dz_obs - null mean) / null sd and its own p_null.
  ROBUSTNESS: the whole Part A profile is recomputed with the window = DUR ONLY (no +-10 d pads)
  and the two versions are compared, to check the result is not made by the pads.
  Dependence: segment (Mann-Whitney), Mw (Spearman), propagation direction (sign of v_t, the
  Theil-Sen slope of the daily mean along-strike position y over defined days), epoch
  (2010-2014 / 2015-2019 / 2020-2026; Kruskal-Wallis + Spearman vs date).
  Depth-vs-intensity: per event Spearman(zmed(d), n(d)) over defined days; stack = median rho with a
  sign test and a Wilcoxon signed-rank against 0.
  COUNT-BIAS CONTROL: log10 of the daily detection count is put through the identical fifths /
  anomaly pipeline.  A daily p10 estimated from fewer detections is pulled toward the daily median
  (i.e. biased DEEP), so if the counts collapse in the last fifth the p10 "deepening" could be an
  estimator artefact.  Reported: the stacked count profile, and Spearman across events between the
  first->last count change and the first->last depth change (a strongly negative rho would mean the
  depth change is riding on the count change).  The MEDIAN depth is essentially count-unbiased and
  is therefore the primary quantity; p10 is the secondary one.

======================================================================================
PART B — tremor rate vs slip rate scaling (Q13)
======================================================================================
  Same window WA, same F / NB / near-tremor masks.
  C(d) = near tremor count on day d.
  R(d) = sum over NB of area_km2 * positive rate(d)   [mm/d km^2, i.e. a POTENCY RATE].
  Per event: Theil-Sen slope of log10(C + 1) on log10(R) over days with C >= 1 and R > 0
  (>= 8 such days required) = the scaling EXPONENT.
  Rising / decaying halves: split at argmax_d R(d) within WA; rising = d <= argmax, decaying =
  d > argmax; each half needs >= 6 usable days.  Hysteresis = paired Wilcoxon(rising, decaying).
  POOLED fit with per-event offsets removed: within-event demeaning of BOTH log10(C+1) and
  log10(R) (algebraically identical to least squares with event dummies), then one OLS slope
  (= the LSDV within-estimator) and one Theil-Sen slope on the pooled demeaned points.  Reported
  combined, per segment and per epoch.
  NULL: 200 rigid tremor-date shifts (>= 365 d, same wrap) -> C from the shifted catalogue paired
  with the UNCHANGED R; per event the null exponent mean/sd and p = fraction of draws with
  null exponent >= observed exponent.
  CAVEAT (stated in the output): R is 7-d smoothed and mixture-limited (13 fixed patterns), so its
  daily dynamic range is compressed -> the exponent is biased HIGH relative to a true
  count-vs-slip-rate exponent; C is completeness-limited and grows with the PNSN network, which
  affects the LEVEL (absorbed by the per-event offset) more than the slope.

======================================================================================
PART C — edges referenced to concurrent tremor (Q2 redo)
======================================================================================
  Per event, over DUR = [ts, te]:
    near_p10 / near_p90 = 10th / 90th percentile Slab2 depth of the near-footprint tremor (60 km)
    seg_p10_conc / seg_p90_conc = the same percentiles over ALL tremor in the event's SEGMENT in
    the SAME window (concurrent, margin-wide reference)
    rel_p10 = near_p10 - seg_p10_conc      (- = the patch tremor reaches further UP-DIP than the
    rel_p90 = near_p90 - seg_p90_conc       concurrent margin-wide band)
  Slip edges (S1's definition, verbatim): acc = clip(cum_mm[:, te] - cum_mm[:, ts - 30], 0);
    band = triangles with lat inside the footprint's latitude span +- 0.1 deg AND acc >= 3 mm;
    slip_updip_edge = min Slab2 depth over band, slip_downdip_edge = max.
  Referenced to the segment's 2010-2026 ALL-TIME tremor percentiles:
    slip_updip_rel   = slip_updip_edge   - seg_p10_alltime   (- = slip edge shallower / up-dip of
                                                              the tremor up-dip edge)
    slip_downdip_rel = slip_downdip_edge - seg_p90_alltime   (+ = slip extends deeper)
  Also the direct within-window version: slip_updip_edge - near_p10 (concurrent, same patch).
  Trends: Theil-Sen slope vs decimal year with a 1000-permutation two-sided p on |slope|, per
  segment and combined, for each of the five quantities.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s9_depth_migration.csv        one row per event: dz_tremor_fifth (bin4-bin0), dz_slip_fifth,
                                        dz_p10_fifth, the DUR-only versions, rho_depth_count (+p),
                                        v_t, n defined days, the 5 binned anomalies for each of the
                                        three profiles, per-event null mean/sd/z/p
  results/s9_depth_migration_stack.csv  one row per (profile, bin): observed stacked median, bootstrap
                                        90 % CI, null median and 5-95 % band, n events
  results/s9_rate_scaling.csv           one row per event: exponent, n days, rising/decaying exponents,
                                        pooled-fit inputs, null mean/sd/p; plus pooled rows (id = -1)
  results/s9_edges.csv                  one row per event: near/segment percentiles, rel_p10, rel_p90,
                                        slip edges and their two referenced versions
  results/s9_edges_trends.csv           Theil-Sen slope, CI and permutation p for every Part C
                                        quantity x group
  figures/s9_within_event.png           9 panels (see FIGURE section at the bottom of the script)

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s9_within_event.py
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
from scipy.stats import (binomtest, kruskal, mannwhitneyu, spearmanr,
                         theilslopes, wilcoxon)

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)

PRE_A = 10              # Part A/B window: start - 10 d
POST_A = 10             # Part A/B window: end + 10 d
EDGE_PRE_DAYS = 30      # Part C slip accumulation window: start - 30 .. end   (S1 definition)
EDGE_ACC_MM = 3.0       # Part C slip edge threshold                            (S1 definition)
NB_RADIUS_KM = 60.0     # neighbourhood / tremor radius around footprint centroids
MIN_DET_DAY = 5         # a day is DEFINED for the tremor track at >= 5 detections
MIN_DEF_DAYS = 8        # an event needs >= 8 defined days for a phase profile
MIN_FIT_DAYS = 8        # Part B: >= 8 usable days for a per-event exponent
MIN_HALF_DAYS = 6       # Part B: >= 6 usable days in a half for the rising/decaying exponent
MIN_DUR_DET = 20        # Part C: >= 20 near detections with finite depth for the percentiles
N_BINS = 5              # fifths of the event
N_SHIFT = 200           # time-shifted null draws per event
N_BOOT = 2000           # bootstrap resamples of events for the stack CI
N_PERM = 1000           # permutations for the Part C trend p-values
RATE_SMOOTH_D = 7       # centred moving average of the daily slip rate
COS47 = np.cos(np.deg2rad(47.0))

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


# ==================================================================================
# 1. Slip model  (loader copied from scripts/s4_fronts.py)
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


def rate_pipeline(cum_mm: np.ndarray) -> np.ndarray:
    """diff along time -> 7-d centred moving average -> positive part.  (n_tri, T) mm/d."""
    r = np.zeros_like(cum_mm)
    r[:, 1:] = np.diff(cum_mm, axis=1)
    r = uniform_filter1d(r, size=RATE_SMOOTH_D, axis=1, mode="nearest")
    return np.clip(r, 0.0, None)


RATE = rate_pipeline(CUM_MM)  # (742, T) mm/d, positive
log(f"  RATE {RATE.shape}, max {RATE.max():.3f} mm/d, mean {RATE.mean():.4f} mm/d")


# ==================================================================================
# 2. Geometry / labels  (loader copied from scripts/s4_fronts.py)
# ==================================================================================
NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]                       # (T, 742)
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
PDEPTH = slabf(PLON, PLAT)
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
DAY_DT = pd.to_datetime(DAYS)
NT = len(DAYS)
assert NT == A_raw.shape[0], "A_raw / days length mismatch"
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
log(f"  {NT} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles")


# ==================================================================================
# 3. Tremor  (loader copied from scripts/s4_fronts.py)
# ==================================================================================
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_LAT = tr["lat"].values
TR_X, TR_Y = to_km(tr["lon"].values, TR_LAT)
TR_DEPTH = slabf(tr["lon"].values, TR_LAT)
TR_SEG = seg_of_lat(TR_LAT)
SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1
TR_OFF = (tr["t"].values.astype("datetime64[D]") - SHIFT_LO).astype(int)  # 0 .. SPAN-1
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(
    f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}; "
    f"Slab2 depth finite for {np.isfinite(TR_DEPTH).mean() * 100:.1f}%; wrap span {SPAN} d"
)

# segment ALL-TIME tremor depth percentiles (2010-01-01 .. 2026-05-27), for Part C
SEG_ALLTIME = {}
in_span = (TR_OFF >= 0) & (TR_OFF < SPAN)
for sname in ("north", "central", "south"):
    m = in_span & (TR_SEG == sname) & np.isfinite(TR_DEPTH)
    SEG_ALLTIME[sname] = (float(np.percentile(TR_DEPTH[m], 10)),
                          float(np.percentile(TR_DEPTH[m], 50)),
                          float(np.percentile(TR_DEPTH[m], 90)), int(m.sum()))
    log(f"  segment {sname:8s} all-time tremor depth p10/p50/p90 = "
        f"{SEG_ALLTIME[sname][0]:.2f} / {SEG_ALLTIME[sname][1]:.2f} / {SEG_ALLTIME[sname][2]:.2f} km "
        f"(n={SEG_ALLTIME[sname][3]})")


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
def window_subset(sorted_off, s, w0_off, nW):
    """Detections whose SHIFTED day (o + s) mod SPAN lands inside the window.

    sorted_off : detection day offsets (0..SPAN-1) SORTED ascending
    s          : rigid shift in days (0 = observed)
    w0_off     : window start expressed as an offset from SHIFT_LO
    Returns (positions into sorted_off, day-in-window index 0..nW-1).
    """
    a = (w0_off - s) % SPAN
    b = (a + nW - 1) % SPAN
    if a <= b:
        idx = np.arange(np.searchsorted(sorted_off, a, "left"),
                        np.searchsorted(sorted_off, b, "right"))
    else:
        idx = np.concatenate([
            np.arange(np.searchsorted(sorted_off, a, "left"), sorted_off.size),
            np.arange(0, np.searchsorted(sorted_off, b, "right")),
        ])
    d = (sorted_off[idx] + s) % SPAN - w0_off
    return idx, d


def daily_stats(dayidx, z, nW, want_y=None):
    """Per-day count, median depth, p10 depth (and mean y) over a window of nW days."""
    cnt = np.bincount(dayidx, minlength=nW).astype(int) if dayidx.size else np.zeros(nW, int)
    zmed = np.full(nW, np.nan)
    zp10 = np.full(nW, np.nan)
    ymean = np.full(nW, np.nan)
    if dayidx.size:
        order = np.argsort(dayidx, kind="stable")
        dd = dayidx[order]
        zz = z[order]
        yy = want_y[order] if want_y is not None else None
        bounds = np.searchsorted(dd, np.arange(nW + 1))
        for i in range(nW):
            a, b = bounds[i], bounds[i + 1]
            if b - a >= MIN_DET_DAY:
                v = zz[a:b]
                v = v[np.isfinite(v)]
                if v.size:
                    zmed[i] = np.median(v)
                    zp10[i] = np.percentile(v, 10)
                if yy is not None:
                    ymean[i] = yy[a:b].mean()
    return cnt, zmed, zp10, ymean


def phase_profile(series, cnt):
    """Fifth-binned ANOMALY profile of `series` over the days with cnt >= MIN_DET_DAY.

    Returns (anom[5], ref, n_def, first, last) or (None, ...) if the event is too thin.
    """
    dfn = np.flatnonzero(cnt >= MIN_DET_DAY)
    if dfn.size < MIN_DEF_DAYS:
        return None, np.nan, dfn.size, -1, -1
    first, last = dfn[0], dfn[-1]
    if last <= first:
        return None, np.nan, dfn.size, first, last
    vals = series[dfn]
    ok = np.isfinite(vals)
    if ok.sum() < MIN_DEF_DAYS:
        return None, np.nan, int(ok.sum()), first, last
    dfn, vals = dfn[ok], vals[ok]
    tau = (dfn - first) / (last - first)
    b = np.minimum((tau * N_BINS).astype(int), N_BINS - 1)
    ref = float(np.median(vals))
    anom = np.full(N_BINS, np.nan)
    for k in range(N_BINS):
        m = b == k
        if m.any():
            anom[k] = float(np.median(vals[m])) - ref
    return anom, ref, int(dfn.size), int(first), int(last)


def slip_phase_profile(zs, Sw, first, last):
    """Same fifth-binning for the slip-weighted depth, on the tremor phase axis."""
    if first < 0 or last <= first:
        return None, np.nan
    d = np.arange(first, last + 1)
    ok = (Sw[d] > 0) & np.isfinite(zs[d])
    if ok.sum() < MIN_DEF_DAYS:
        return None, np.nan
    d, vals = d[ok], zs[d][ok]
    tau = (d - first) / (last - first)
    b = np.minimum((tau * N_BINS).astype(int), N_BINS - 1)
    ref = float(np.median(vals))
    anom = np.full(N_BINS, np.nan)
    for k in range(N_BINS):
        m = b == k
        if m.any():
            anom[k] = float(np.median(vals[m])) - ref
    return anom, ref


def ts_slope(x, y, nmin=5):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < nmin:
        return np.nan
    if np.ptp(x[m]) == 0:
        return np.nan
    return float(theilslopes(y[m], x[m], 0.95)[0])


def ts_trend_perm(x, y, nperm=N_PERM, nmin=5):
    """Theil-Sen slope with a two-sided permutation p on |slope| and a 95 % CI."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = x.size
    if n < nmin or np.ptp(x) == 0:
        return np.nan, np.nan, np.nan, np.nan, n
    sl, ic, lo, hi = theilslopes(y, x, 0.95)
    obs = abs(sl)
    cnt = 0
    for _ in range(nperm):
        yp = RNG.permutation(y)
        cnt += abs(theilslopes(yp, x, 0.95)[0]) >= obs
    return float(sl), float(lo), float(hi), (cnt + 1) / (nperm + 1), n


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
        return "     n/a       (n=0)"
    return f"{m:+8.3f} [{q1:+7.3f},{q3:+7.3f}] (n={n:2d})"


# ==================================================================================
# 6. Per-event measurement (Parts A, B, C in one pass over events)
# ==================================================================================
log("measuring events (A: depth migration, B: rate scaling, C: edges) ...")
rowsA, rowsB, rowsC = [], [], []
PROF_OBS = {"tremor": {}, "slip": {}, "p10": {}, "logcount": {}}   # pid -> anom[5]
PROF_NULL = {"tremor": {}, "p10": {}}                 # pid -> (ndraw, 5)
POOL = []                                             # Part B pooled points

for _, r in ev_in.iterrows():
    pid = int(r["id"])
    seg = str(r["seg"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    w0_off = w0 - IDX_OF_SHIFT_LO

    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = np.unique(np.concatenate([np.asarray(h, int) for h in hits])) if len(hits) else np.array([], int)
    nb = np.union1d(nb, np.flatnonzero(F))

    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
            if len(thits) else np.array([], int))

    # ---- slip daily series over the window ----
    Rw = RATE[np.ix_(nb, np.arange(w0, w1 + 1))]              # (nb, nW) mm/d
    aw = PAREA[nb][:, None]
    Sw = (aw * Rw).sum(axis=0)                                # mm/d * km^2  (potency rate)
    zs = np.full(nW, np.nan)
    pos = Sw > 0
    if pos.any():
        zs[pos] = ((aw * Rw * PDEPTH[nb][:, None]).sum(axis=0)[pos]) / Sw[pos]

    # ---- observed tremor daily series ----
    off = TR_OFF[near]
    order = np.argsort(off, kind="stable")
    soff = off[order]
    sz = TR_DEPTH[near][order]
    sy = TR_Y[near][order]
    sub, dday = window_subset(soff, 0, w0_off, nW)
    cnt, zmed, zp10, ymean = daily_stats(dday, sz[sub], nW, want_y=sy[sub])

    # ---- PART A ----
    anom_t, ref_t, n_def, first, last = phase_profile(zmed, cnt)
    anom_p, ref_p, _, _, _ = phase_profile(zp10, cnt)
    lcnt = np.where(cnt > 0, np.log10(np.maximum(cnt, 1).astype(float)), np.nan)
    anom_c, ref_c, _, _, _ = phase_profile(lcnt, cnt)
    anom_s, ref_s = slip_phase_profile(zs, Sw, first, last)
    dfn = np.flatnonzero(cnt >= MIN_DET_DAY)
    v_t = ts_slope(dfn.astype(float), ymean[dfn]) if dfn.size else np.nan

    rho_dc, p_dc = np.nan, np.nan
    if dfn.size >= MIN_DEF_DAYS:
        mm = np.isfinite(zmed[dfn])
        if mm.sum() >= MIN_DEF_DAYS:
            rho_dc, p_dc = spearmanr(zmed[dfn][mm], cnt[dfn][mm])

    recA = dict(id=pid, start=r["start"], seg=seg, Mw=float(r["Mw"]), dur_d=int(r["dur_d"]),
                epoch=epoch_of(r["start"]), decyear=decyear(r["start"]),
                n_tri_F=int(F.sum()), n_tri_NB=int(nb.size), n_near_alltime=int(near.size),
                n_days_W=nW, n_def_days=n_def, first_def=first, last_def=last,
                event_med_depth_km=ref_t, event_med_p10_km=ref_p, event_med_slipdepth_km=ref_s,
                v_t_km_d=v_t, rho_depth_count=float(rho_dc), p_depth_count=float(p_dc))
    for k in range(N_BINS):
        recA[f"tremor_anom_b{k}"] = anom_t[k] if anom_t is not None else np.nan
        recA[f"p10_anom_b{k}"] = anom_p[k] if anom_p is not None else np.nan
        recA[f"slip_anom_b{k}"] = anom_s[k] if anom_s is not None else np.nan
        recA[f"logcount_anom_b{k}"] = anom_c[k] if anom_c is not None else np.nan
    recA["dlogcount_fifth"] = (anom_c[4] - anom_c[0]) if anom_c is not None else np.nan
    recA["dz_tremor_fifth_km"] = (anom_t[4] - anom_t[0]) if anom_t is not None else np.nan
    recA["dz_p10_fifth_km"] = (anom_p[4] - anom_p[0]) if anom_p is not None else np.nan
    recA["dz_slip_fifth_km"] = (anom_s[4] - anom_s[0]) if anom_s is not None else np.nan
    if anom_t is not None:
        PROF_OBS["tremor"][pid] = anom_t
    if anom_p is not None:
        PROF_OBS["p10"][pid] = anom_p
    if anom_s is not None:
        PROF_OBS["slip"][pid] = anom_s
    if anom_c is not None:
        PROF_OBS["logcount"][pid] = anom_c

    # ---- ROBUSTNESS: the same profile on the DUR window only (no +-10 d pads) ----
    nD = te - ts + 1
    subD, ddD = window_subset(soff, 0, ts - IDX_OF_SHIFT_LO, nD)
    cD, zD, pD, _ = daily_stats(ddD, sz[subD], nD)
    aD, _, _, _, _ = phase_profile(zD, cD)
    aDp, _, _, _, _ = phase_profile(pD, cD)
    recA["dz_tremor_fifth_DURonly_km"] = (aD[4] - aD[0]) if aD is not None else np.nan
    recA["dz_p10_fifth_DURonly_km"] = (aDp[4] - aDp[0]) if aDp is not None else np.nan

    # ---- PART B (observed) ----
    okB = (cnt >= 1) & (Sw > 0)
    lc = np.log10(cnt[okB] + 1.0)
    lr = np.log10(Sw[okB])
    expo = ts_slope(lr, lc, MIN_FIT_DAYS) if okB.sum() >= MIN_FIT_DAYS else np.nan
    imax = int(np.argmax(Sw))
    d_all = np.arange(nW)
    ris = okB & (d_all <= imax)
    dec = okB & (d_all > imax)
    exp_r = ts_slope(np.log10(Sw[ris]), np.log10(cnt[ris] + 1.0), MIN_HALF_DAYS) if ris.sum() >= MIN_HALF_DAYS else np.nan
    exp_d = ts_slope(np.log10(Sw[dec]), np.log10(cnt[dec] + 1.0), MIN_HALF_DAYS) if dec.sum() >= MIN_HALF_DAYS else np.nan
    if okB.sum() >= MIN_FIT_DAYS:
        POOL.append(pd.DataFrame(dict(id=pid, seg=seg, epoch=epoch_of(r["start"]),
                                      lc=lc, lr=lr)))
    recB = dict(id=pid, start=r["start"], seg=seg, Mw=float(r["Mw"]),
                epoch=epoch_of(r["start"]), decyear=decyear(r["start"]),
                n_days_fit=int(okB.sum()), n_rise=int(ris.sum()), n_decay=int(dec.sum()),
                peak_day_in_W=imax, R_max=float(Sw.max()), C_max=int(cnt.max()),
                C_sum=int(cnt.sum()), exponent=expo, exponent_rise=exp_r, exponent_decay=exp_d)

    # ---- PART C ----
    dsub, dd_dur = window_subset(soff, 0, ts - IDX_OF_SHIFT_LO, te - ts + 1)
    zdur = sz[dsub]
    zdur = zdur[np.isfinite(zdur)]
    near_p10 = float(np.percentile(zdur, 10)) if zdur.size >= MIN_DUR_DET else np.nan
    near_p50 = float(np.percentile(zdur, 50)) if zdur.size >= MIN_DUR_DET else np.nan
    near_p90 = float(np.percentile(zdur, 90)) if zdur.size >= MIN_DUR_DET else np.nan
    segm = (TR_SEG == seg) & (TR_OFF >= ts - IDX_OF_SHIFT_LO) & (TR_OFF <= te - IDX_OF_SHIFT_LO) \
        & np.isfinite(TR_DEPTH)
    zseg = TR_DEPTH[segm]
    seg_p10c = float(np.percentile(zseg, 10)) if zseg.size >= MIN_DUR_DET else np.nan
    seg_p50c = float(np.percentile(zseg, 50)) if zseg.size >= MIN_DUR_DET else np.nan
    seg_p90c = float(np.percentile(zseg, 90)) if zseg.size >= MIN_DUR_DET else np.nan

    acc = np.clip(CUM_MM[:, te] - CUM_MM[:, max(ts - EDGE_PRE_DAYS, 0)], 0.0, None)
    lat_lo, lat_hi = PLAT[F].min() - 0.1, PLAT[F].max() + 0.1
    band = (PLAT >= lat_lo) & (PLAT <= lat_hi) & (acc >= EDGE_ACC_MM)
    up_e = float(PDEPTH[band].min()) if band.any() else np.nan
    dn_e = float(PDEPTH[band].max()) if band.any() else np.nan
    a10, a50, a90, _ = SEG_ALLTIME[seg]

    recC = dict(id=pid, start=r["start"], seg=seg, Mw=float(r["Mw"]), dur_d=int(r["dur_d"]),
                epoch=epoch_of(r["start"]), decyear=decyear(r["start"]),
                n_near_dur=int(zdur.size), n_seg_dur=int(zseg.size),
                near_p10=near_p10, near_p50=near_p50, near_p90=near_p90,
                seg_p10_conc=seg_p10c, seg_p50_conc=seg_p50c, seg_p90_conc=seg_p90c,
                rel_p10=near_p10 - seg_p10c, rel_p50=near_p50 - seg_p50c,
                rel_p90=near_p90 - seg_p90c,
                slip_updip_edge=up_e, slip_downdip_edge=dn_e, n_tri_edgeband=int(band.sum()),
                seg_p10_alltime=a10, seg_p90_alltime=a90,
                slip_updip_rel=up_e - a10, slip_downdip_rel=dn_e - a90,
                slip_updip_minus_near_p10=up_e - near_p10,
                slip_downdip_minus_near_p90=dn_e - near_p90)

    # ---- NULL (shared by A and B): 200 rigid shifts of the near tremor ----
    shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
    null_t = np.full((N_SHIFT, N_BINS), np.nan)   # draw index j is SHARED across events
    null_p = np.full((N_SHIFT, N_BINS), np.nan)
    null_e = np.full(N_SHIFT, np.nan)
    for j, sh in enumerate(shifts):
        s = int(sh)
        subn, ddn = window_subset(soff, s, w0_off, nW)
        cn, zn, pn, _ = daily_stats(ddn, sz[subn], nW)
        at, _, _, _, _ = phase_profile(zn, cn)
        ap, _, _, _, _ = phase_profile(pn, cn)
        if at is not None:
            null_t[j] = at
        if ap is not None:
            null_p[j] = ap
        okn = (cn >= 1) & (Sw > 0)
        if okn.sum() >= MIN_FIT_DAYS:
            en = ts_slope(np.log10(Sw[okn]), np.log10(cn[okn] + 1.0), MIN_FIT_DAYS)
            if np.isfinite(en):
                null_e[j] = en
    PROF_NULL["tremor"][pid] = null_t
    PROF_NULL["p10"][pid] = null_p
    ne = null_e[np.isfinite(null_e)]

    recA["n_null_prof"] = int(np.isfinite(null_t[:, 0]).sum())
    nd = null_t[:, 4] - null_t[:, 0]
    nd = nd[np.isfinite(nd)]
    if nd.size >= 20:
        recA["null_dz_mean"] = float(nd.mean())
        recA["null_dz_sd"] = float(nd.std())
        recA["null_dz_p"] = (float((np.abs(nd) >= abs(recA["dz_tremor_fifth_km"])).mean())
                             if np.isfinite(recA["dz_tremor_fifth_km"]) else np.nan)
        recA["null_dz_z"] = ((recA["dz_tremor_fifth_km"] - nd.mean()) / nd.std()
                             if nd.std() > 0 else np.nan)
    else:
        recA["null_dz_mean"] = recA["null_dz_sd"] = np.nan
        recA["null_dz_p"] = recA["null_dz_z"] = np.nan
    recB["n_null_exp"] = int(ne.size)
    recB["null_exp_mean"] = float(ne.mean()) if ne.size else np.nan
    recB["null_exp_sd"] = float(ne.std()) if ne.size else np.nan
    recB["null_exp_p"] = (float((ne >= expo).mean())
                          if (ne.size and np.isfinite(expo)) else np.nan)

    rowsA.append(recA)
    rowsB.append(recB)
    rowsC.append(recC)
    log(f"  {pid:5d} {r['start']} {seg:7s} Mw{float(r['Mw']):.2f} nb={nb.size:3d} near={near.size:6d} "
        f"def={n_def:3d} dzT={recA['dz_tremor_fifth_km']:+6.2f} dzS={recA['dz_slip_fifth_km']:+6.2f} "
        f"rho={rho_dc:+.2f} b={expo if np.isfinite(expo) else float('nan'):+6.2f} "
        f"nullb={recB['null_exp_mean'] if np.isfinite(recB['null_exp_mean']) else float('nan'):+5.2f}")

A = pd.DataFrame(rowsA)
B = pd.DataFrame(rowsB)
C = pd.DataFrame(rowsC)
log(f"per-event pass done: A {A.shape} B {B.shape} C {C.shape}")


# ==================================================================================
# 7. PART A — stacks, nulls, dependences
# ==================================================================================
print("\n" + "=" * 100)
print("PART A — TREMOR DEPTH MIGRATION WITHIN EVENTS (Q14)")
print("=" * 100)
print(f"window [start-{PRE_A}, end+{POST_A}]; near = 60 km of footprint; day DEFINED at >= "
      f"{MIN_DET_DAY} detections; event needs >= {MIN_DEF_DAYS} defined days; fifths of the phase "
      f"between the first and last defined day; anomaly = bin median - event median.")

stack_rows = []
STACK = {}
for prof in ("tremor", "slip", "p10", "logcount"):
    ids = sorted(PROF_OBS[prof])
    M = np.vstack([PROF_OBS[prof][i] for i in ids]) if ids else np.zeros((0, N_BINS))
    obs = np.nanmedian(M, axis=0) if M.size else np.full(N_BINS, np.nan)
    # bootstrap over events
    lo = np.full(N_BINS, np.nan)
    hi = np.full(N_BINS, np.nan)
    if M.shape[0] >= 5:
        bs = np.empty((N_BOOT, N_BINS))
        for i in range(N_BOOT):
            sel = RNG.integers(0, M.shape[0], M.shape[0])
            bs[i] = np.nanmedian(M[sel], axis=0)
        lo = np.nanpercentile(bs, 5, axis=0)
        hi = np.nanpercentile(bs, 95, axis=0)
    STACK[prof] = dict(M=M, ids=ids, obs=obs, lo=lo, hi=hi)
    print(f"\n  {prof.upper()} profile: n_events = {M.shape[0]}")
    for k in range(N_BINS):
        print(f"    fifth {k + 1}  median {obs[k]:+7.3f} km  90%CI [{lo[k]:+7.3f},{hi[k]:+7.3f}]  "
              f"n={int(np.isfinite(M[:, k]).sum()) if M.size else 0}")

# ---- null stacks for tremor and p10 ----
NULLSTACK = {}
MIN_EV_PER_DRAW = 15                      # a null draw counts only if >= 15 events contributed
for prof in ("tremor", "p10"):
    ids = sorted(PROF_NULL[prof])
    cube = np.stack([PROF_NULL[prof][i] for i in ids])       # (n_ev, N_SHIFT, N_BINS)
    nev_j = np.isfinite(cube[:, :, 0]).sum(axis=0)           # events available per draw
    keep = nev_j >= MIN_EV_PER_DRAW
    nsm = np.nanmedian(cube[:, keep, :], axis=0)             # (n_kept_draws, N_BINS)
    ndraw = int(keep.sum())
    NULLSTACK[prof] = dict(nsm=nsm, n_ev=len(ids), ndraw=ndraw,
                           nev_per_draw_med=float(np.median(nev_j[keep])) if ndraw else np.nan)
    print(f"  {prof.upper()} null: {ndraw}/{N_SHIFT} draws kept (>= {MIN_EV_PER_DRAW} events); "
          f"median events per kept draw = {NULLSTACK[prof]['nev_per_draw_med']}")
    nlo = np.nanpercentile(nsm, 5, axis=0)
    nhi = np.nanpercentile(nsm, 95, axis=0)
    nmed = np.nanmedian(nsm, axis=0)
    obs = STACK[prof]["obs"]
    print(f"\n  {prof.upper()} NULL stack ({ndraw} shift draws over {len(ids)} events):")
    for k in range(N_BINS):
        print(f"    fifth {k + 1}  obs {obs[k]:+7.3f}  null med {nmed[k]:+7.3f}  "
              f"null 5-95% [{nlo[k]:+7.3f},{nhi[k]:+7.3f}]")
    ndz = nsm[:, 4] - nsm[:, 0]
    odz = obs[4] - obs[0]
    p = float((np.abs(ndz) >= abs(odz)).mean()) if np.isfinite(ndz).any() else np.nan
    print(f"    first->last fifth change: OBSERVED {odz:+.3f} km   null mean {np.nanmean(ndz):+.3f} "
          f"sd {np.nanstd(ndz):.3f}  -> p = {p:.4f}")
    NULLSTACK[prof]["p_dz"] = p
    NULLSTACK[prof]["obs_dz"] = odz
    NULLSTACK[prof]["null_dz"] = ndz
    for k in range(N_BINS):
        stack_rows.append(dict(profile=prof, bin=k + 1, phase_lo=k / N_BINS, phase_hi=(k + 1) / N_BINS,
                               obs_median=obs[k], boot_lo90=STACK[prof]["lo"][k],
                               boot_hi90=STACK[prof]["hi"][k], null_median=nmed[k],
                               null_lo90=nlo[k], null_hi90=nhi[k],
                               n_events=int(np.isfinite(STACK[prof]["M"][:, k]).sum()),
                               n_null_draws=ndraw))
# slip has no tremor-shift null (the slip series is unchanged under the shift)
for k in range(N_BINS):
    stack_rows.append(dict(profile="logcount", bin=k + 1, phase_lo=k / N_BINS,
                           phase_hi=(k + 1) / N_BINS,
                           obs_median=STACK["logcount"]["obs"][k],
                           boot_lo90=STACK["logcount"]["lo"][k],
                           boot_hi90=STACK["logcount"]["hi"][k], null_median=np.nan,
                           null_lo90=np.nan, null_hi90=np.nan,
                           n_events=int(np.isfinite(STACK["logcount"]["M"][:, k]).sum()),
                           n_null_draws=0))
for k in range(N_BINS):
    stack_rows.append(dict(profile="slip", bin=k + 1, phase_lo=k / N_BINS, phase_hi=(k + 1) / N_BINS,
                           obs_median=STACK["slip"]["obs"][k], boot_lo90=STACK["slip"]["lo"][k],
                           boot_hi90=STACK["slip"]["hi"][k], null_median=np.nan,
                           null_lo90=np.nan, null_hi90=np.nan,
                           n_events=int(np.isfinite(STACK["slip"]["M"][:, k]).sum()), n_null_draws=0))

# slip first->last, bootstrap only
sM = STACK["slip"]["M"]
sdz = sM[:, 4] - sM[:, 0]
bs = np.array([np.nanmedian(sdz[RNG.integers(0, sdz.size, sdz.size)]) for _ in range(N_BOOT)])
print(f"\n  SLIP first->last fifth change (mixture-limited): median {np.nanmedian(sdz):+.3f} km  "
      f"boot 90%CI [{np.nanpercentile(bs, 5):+.3f},{np.nanpercentile(bs, 95):+.3f}]  "
      f"sign test {int((sdz > 0).sum())}/{int(np.isfinite(sdz).sum())} positive  "
      f"p={binomtest(int((sdz > 0).sum()), int(np.isfinite(sdz).sum())).pvalue:.4f}")

print("\n  per-event first->last fifth change (km, + = deepening):")
print("    id      date        seg      Mw   ndef   dz_tremor  dz_p10  dz_slip   v_t(km/d)  rho(z,n)")
for _, x in A.iterrows():
    print(f"    {int(x.id):5d}  {x.start}  {x.seg:7s} {x.Mw:.2f} {int(x.n_def_days):4d}   "
          f"{x.dz_tremor_fifth_km:+8.2f} {x.dz_p10_fifth_km:+7.2f} {x.dz_slip_fifth_km:+8.2f}   "
          f"{x.v_t_km_d:+8.3f}   {x.rho_depth_count:+.3f}")

zz_a = A.null_dz_z.dropna()
print(f"\n  PER-EVENT null (each event vs its own 200 shifted windows):")
print(f"    z = (dz_obs - null mean)/null sd : {fmt_mi(zz_a)}")
print(f"    sign test {int((zz_a > 0).sum())}/{zz_a.size} positive -> p = "
      f"{binomtest(int((zz_a > 0).sum()), zz_a.size).pvalue:.4f};  "
      f"Wilcoxon vs 0 p = {wilcoxon(zz_a).pvalue:.4f}")
print(f"    events with per-event p_null < 0.05: {int((A.null_dz_p < 0.05).sum())}/{len(A)}")

print("\n  ROBUSTNESS — same profile with the window = DUR only (no +-10 d pads):")
print(f"    dz_tremor (DUR only) {fmt_mi(A.dz_tremor_fifth_DURonly_km)}  vs padded "
      f"{fmt_mi(A.dz_tremor_fifth_km)}")
rho_r, p_r = spearmanr(A.dz_tremor_fifth_km, A.dz_tremor_fifth_DURonly_km, nan_policy="omit")
print(f"    Spearman padded vs DUR-only rho = {rho_r:+.3f} p = {p_r:.2e}")

print("\n  DEPENDENCES of dz_tremor (first->last fifth depth change):")
for grp, col in (("segment", "seg"), ("epoch", "epoch")):
    for g, sub in A.groupby(col):
        print(f"    {grp:8s} {g:10s} {fmt_mi(sub.dz_tremor_fifth_km)}")
n_north = A[A.seg == "north"].dz_tremor_fifth_km.dropna()
n_cent = A[A.seg == "central"].dz_tremor_fifth_km.dropna()
if n_north.size >= 3 and n_cent.size >= 3:
    print(f"    Mann-Whitney north vs central: p = {mannwhitneyu(n_north, n_cent).pvalue:.4f}")
eg = [g.dz_tremor_fifth_km.dropna().values for _, g in A.groupby("epoch")]
eg = [g for g in eg if g.size >= 3]
if len(eg) >= 2:
    print(f"    Kruskal-Wallis across epochs: p = {kruskal(*eg).pvalue:.4f}")
rho, p = spearmanr(A.Mw, A.dz_tremor_fifth_km, nan_policy="omit")
print(f"    Spearman vs Mw:      rho = {rho:+.3f}  p = {p:.4f}")
rho, p = spearmanr(A.decyear, A.dz_tremor_fifth_km, nan_policy="omit")
print(f"    Spearman vs date:    rho = {rho:+.3f}  p = {p:.4f}")
A["prop_dir"] = np.where(A.v_t_km_d > 0, "north", np.where(A.v_t_km_d < 0, "south", "flat"))
for g, sub in A.groupby("prop_dir"):
    print(f"    prop_dir {g:6s} {fmt_mi(sub.dz_tremor_fifth_km)}")
pn = A[A.prop_dir == "north"].dz_tremor_fifth_km.dropna()
ps = A[A.prop_dir == "south"].dz_tremor_fifth_km.dropna()
if pn.size >= 3 and ps.size >= 3:
    print(f"    Mann-Whitney northward vs southward propagation: p = {mannwhitneyu(pn, ps).pvalue:.4f}")

print("\n  COUNT-BIAS CONTROL — stacked log10(daily count) anomaly by fifth "
      "(a collapse in the last fifth would bias the daily p10 deep):")
cM = STACK["logcount"]
for k in range(N_BINS):
    print(f"    fifth {k + 1}  median {cM['obs'][k]:+7.3f} dex  90%CI "
          f"[{cM['lo'][k]:+7.3f},{cM['hi'][k]:+7.3f}]  "
          f"(x{10 ** cM['obs'][k]:.2f} of the event median count)")
dlc = A.dlogcount_fifth
print(f"    per-event first->last count change {fmt_mi(dlc)} dex")
for lab, col in (("dz_tremor", A.dz_tremor_fifth_km), ("dz_p10", A.dz_p10_fifth_km)):
    rho_c, p_c = spearmanr(dlc, col, nan_policy="omit")
    print(f"    Spearman(count change, {lab}) rho = {rho_c:+.3f}  p = {p_c:.4f}")
print("    (a strongly NEGATIVE rho would mean the depth change rides on the count collapse; the "
      "daily MEDIAN depth is essentially count-unbiased, the daily p10 is not)")

print("\n  DEPTH vs INTENSITY within events: Spearman(daily median depth, daily count)")
rr = A.rho_depth_count.dropna()
print(f"    per-event rho {fmt_mi(rr)}")
npos = int((rr > 0).sum())
print(f"    sign test {npos}/{rr.size} positive -> p = {binomtest(npos, rr.size).pvalue:.4f}; "
      f"Wilcoxon vs 0 p = {wilcoxon(rr).pvalue:.4f}")
print(f"    nominally significant per event (p<0.05): "
      f"{int((A.p_depth_count < 0.05).sum())}/{len(A)} "
      f"({int(((A.p_depth_count < 0.05) & (A.rho_depth_count > 0)).sum())} positive, "
      f"{int(((A.p_depth_count < 0.05) & (A.rho_depth_count < 0)).sum())} negative)")
print("    (+rho = tremor is DEEPER on its most intense days; -rho = shallower when most intense)")

print("\n  UP-DIP FRONT (daily p10 depth) first->last fifth: "
      f"{fmt_mi(A.dz_p10_fifth_km)}  null p = {NULLSTACK['p10']['p_dz']:.4f}")


# ==================================================================================
# 8. PART B — rate scaling
# ==================================================================================
print("\n" + "=" * 100)
print("PART B — TREMOR RATE vs SLIP RATE SCALING (Q13)")
print("=" * 100)
print("exponent b = Theil-Sen slope of log10(C+1) on log10(R); C = daily near-footprint tremor "
      "count, R = area-weighted positive modeled slip rate over the 60-km neighbourhood "
      "[mm/d km^2]; days with C >= 1 and R > 0.")
print("\n    id      date        seg     ndays  b        b_rise   b_decay   null_b(mean+-sd)  p_null")
for _, x in B.iterrows():
    print(f"    {int(x.id):5d}  {x.start}  {x.seg:7s} {int(x.n_days_fit):4d}  "
          f"{x.exponent:+7.3f}  {x.exponent_rise:+7.3f}  {x.exponent_decay:+7.3f}   "
          f"{x.null_exp_mean:+6.3f}+-{x.null_exp_sd:5.3f}      "
          f"{x.null_exp_p if np.isfinite(x.null_exp_p) else float('nan'):.3f}")

print(f"\n  exponent, combined      {fmt_mi(B.exponent)}")
for g, sub in B.groupby("seg"):
    print(f"  exponent, seg {g:8s} {fmt_mi(sub.exponent)}")
for g, sub in B.groupby("epoch"):
    print(f"  exponent, epoch {g:10s} {fmt_mi(sub.exponent)}")
print(f"  null exponent pooled    {fmt_mi(B.null_exp_mean)}   (per-event null sd "
      f"{fmt_mi(B.null_exp_sd)})")
sig = B.null_exp_p < 0.05
print(f"  events with p_null < 0.05: {int(sig.sum())}/{int(np.isfinite(B.null_exp_p).sum())}")
zz = (B.exponent - B.null_exp_mean) / B.null_exp_sd
print(f"  per-event z vs its own null {fmt_mi(zz)};  sign test "
      f"{int((zz > 0).sum())}/{int(np.isfinite(zz).sum())} positive p="
      f"{binomtest(int((zz > 0).sum()), int(np.isfinite(zz).sum())).pvalue:.2e}")

mr = np.isfinite(B.exponent_rise) & np.isfinite(B.exponent_decay)
print(f"\n  HYSTERESIS (rising vs decaying half), n = {int(mr.sum())}")
print(f"    rising   {fmt_mi(B.exponent_rise[mr])}")
print(f"    decaying {fmt_mi(B.exponent_decay[mr])}")
if mr.sum() >= 6:
    d = (B.exponent_rise - B.exponent_decay)[mr]
    print(f"    rise - decay {fmt_mi(d)}  Wilcoxon p = {wilcoxon(d).pvalue:.4f}  "
          f"sign {int((d > 0).sum())}/{int(d.size)} positive")

POOLDF = pd.concat(POOL, ignore_index=True)
POOLDF["lc_c"] = POOLDF.lc - POOLDF.groupby("id").lc.transform("mean")
POOLDF["lr_c"] = POOLDF.lr - POOLDF.groupby("id").lr.transform("mean")


def pooled_slope(df):
    if len(df) < 20:
        return np.nan, np.nan, len(df), np.nan
    x, y = df.lr_c.values, df.lc_c.values
    ols = float(np.polyfit(x, y, 1)[0])
    tss = float(theilslopes(y, x, 0.95)[0])
    r = float(np.corrcoef(x, y)[0, 1])
    return ols, tss, len(df), r


print("\n  POOLED fit, per-event offsets removed (within transformation = event dummies):")
o, t, n, rr_ = pooled_slope(POOLDF)
print(f"    combined            OLS(LSDV) {o:+.3f}   Theil-Sen {t:+.3f}   n_days {n}   r {rr_:+.3f}")
for g, sub in POOLDF.groupby("seg"):
    o, t, n, rr_ = pooled_slope(sub)
    print(f"    seg   {g:10s}    OLS {o:+.3f}   Theil-Sen {t:+.3f}   n_days {n}   r {rr_:+.3f}")
for g, sub in POOLDF.groupby("epoch"):
    o, t, n, rr_ = pooled_slope(sub)
    print(f"    epoch {g:10s}    OLS {o:+.3f}   Theil-Sen {t:+.3f}   n_days {n}   r {rr_:+.3f}")

print("\n  CAVEAT: R is a 7-d smoothed, 13-fixed-pattern mixture -> its daily dynamic range is "
      "compressed, which biases b HIGH; C is completeness-limited (PNSN network growth), which "
      "moves the per-event LEVEL (removed by the offsets) more than the slope.")


# ==================================================================================
# 9. PART C — edges vs concurrent tremor
# ==================================================================================
print("\n" + "=" * 100)
print("PART C — EDGES REFERENCED TO CONCURRENT TREMOR (Q2 redo)")
print("=" * 100)
print("rel_p10/p90 = near-footprint tremor percentile during DUR minus the SAME-SEGMENT, "
      "SAME-WINDOW margin-wide percentile.  slip_*_rel = slip edge (3 mm over start-30..end, "
      "within the footprint latitude span) minus the segment's 2010-2026 tremor p10/p90.")
for sname in ("north", "central"):
    a10, a50, a90, an = SEG_ALLTIME[sname]
    print(f"  reference: {sname:8s} all-time tremor p10/p50/p90 = {a10:.2f}/{a50:.2f}/{a90:.2f} km (n={an})")

print("\n    id      date        seg     near_p10 seg_p10c rel_p10  near_p90 seg_p90c rel_p90  "
      "slip_up  slip_dn  up_rel   dn_rel   up-nearp10")
for _, x in C.iterrows():
    print(f"    {int(x.id):5d}  {x.start}  {x.seg:7s} {x.near_p10:7.2f} {x.seg_p10_conc:8.2f} "
          f"{x.rel_p10:+7.2f}  {x.near_p90:7.2f} {x.seg_p90_conc:8.2f} {x.rel_p90:+7.2f}  "
          f"{x.slip_updip_edge:7.2f} {x.slip_downdip_edge:8.2f} {x.slip_updip_rel:+7.2f} "
          f"{x.slip_downdip_rel:+7.2f} {x.slip_updip_minus_near_p10:+9.2f}")

QC = ["rel_p10", "rel_p50", "rel_p90", "slip_updip_rel", "slip_downdip_rel",
      "slip_updip_minus_near_p10", "slip_downdip_minus_near_p90"]
print("\n  medians [IQR]:")
for q in QC:
    print(f"    {q:28s} combined {fmt_mi(C[q])}")
    for g, sub in C.groupby("seg"):
        print(f"    {'':28s} {g:8s} {fmt_mi(sub[q])}")

print("\n  TRENDS vs date (Theil-Sen km/yr, 1000-permutation two-sided p on |slope|):")
trend_rows = []
for q in QC:
    for gname, sub in [("combined", C)] + [(g, s) for g, s in C.groupby("seg")]:
        sl, lo, hi, p, n = ts_trend_perm(sub.decyear.values, sub[q].values)
        trend_rows.append(dict(quantity=q, group=gname, slope_km_per_yr=sl, ci_lo=lo, ci_hi=hi,
                               p_perm=p, n=n))
        if n >= 5:
            print(f"    {q:28s} {gname:8s} n={n:2d}  slope {sl:+7.4f} [{lo:+7.4f},{hi:+7.4f}]  "
                  f"p_perm {p:.3f}")
TR = pd.DataFrame(trend_rows)
nt = int((TR.n >= 5).sum())
print(f"    ({nt} tests in this block -> Bonferroni FWER 0.05 needs p < {0.05 / max(nt, 1):.4f}; "
      f"min p = {TR[TR.n >= 5].p_perm.min():.3f})")

print("\n  KEY COMPARISON — is the slip up-dip edge shallower than the tremor up-dip edge?")
for gname, sub in [("combined", C)] + [(g, s) for g, s in C.groupby("seg")]:
    v = sub.slip_updip_minus_near_p10.dropna()
    if v.size >= 3:
        print(f"    {gname:8s} slip_updip - concurrent near tremor p10 = {fmt_mi(v)} km  "
              f"({int((v < 0).sum())}/{v.size} negative = slip shallower)")
    v2 = sub.slip_updip_rel.dropna()
    if v2.size >= 3:
        print(f"    {gname:8s} slip_updip - segment all-time tremor p10 = {fmt_mi(v2)} km  "
              f"({int((v2 < 0).sum())}/{v2.size} negative)")
    v3 = sub.slip_downdip_minus_near_p90.dropna()
    if v3.size >= 3:
        print(f"    {gname:8s} slip_downdip - concurrent near tremor p90 = {fmt_mi(v3)} km  "
              f"({int((v3 > 0).sum())}/{v3.size} positive = slip deeper)")

print("\n  CAVEAT: S1b showed the 40-km band around the union of northern footprints contains "
      "100 % of lat>=46N detections, so in the NORTH 'near the footprint' is essentially "
      "'all northern tremor' and rel_p10/rel_p90 are near-vacuous by construction; they are a "
      "real filter only in the central segment.")


# ==================================================================================
# 10. WRITE
# ==================================================================================
A.to_csv(OUT_RES / "s9_depth_migration.csv", index=False)
pd.DataFrame(stack_rows).to_csv(OUT_RES / "s9_depth_migration_stack.csv", index=False)

pooled_rows = []
for gname, sub in ([("combined", POOLDF)]
                   + [(f"seg_{g}", s) for g, s in POOLDF.groupby("seg")]
                   + [(f"epoch_{g}", s) for g, s in POOLDF.groupby("epoch")]):
    o, t, n, rr_ = pooled_slope(sub)
    pooled_rows.append(dict(id=-1, start="", seg=gname, Mw=np.nan, epoch="POOLED", decyear=np.nan,
                            n_days_fit=n, exponent=t, exponent_rise=np.nan, exponent_decay=np.nan,
                            null_exp_mean=np.nan, null_exp_sd=np.nan, null_exp_p=np.nan,
                            pooled_ols=o, pooled_theilsen=t, pooled_r=rr_))
Bout = pd.concat([B, pd.DataFrame(pooled_rows)], ignore_index=True)
Bout.to_csv(OUT_RES / "s9_rate_scaling.csv", index=False)

Cout = C.merge(TR.pivot_table(index="quantity", columns="group", values=["slope_km_per_yr", "p_perm"]),
               how="left", left_index=True, right_index=True) if False else C
Cout.to_csv(OUT_RES / "s9_edges.csv", index=False)
TR.to_csv(OUT_RES / "s9_edges_trends.csv", index=False)
log("wrote results/s9_depth_migration.csv, s9_depth_migration_stack.csv, s9_rate_scaling.csv, "
    "s9_edges.csv, s9_edges_trends.csv")


# ==================================================================================
# 11. FIGURE
# ==================================================================================
fig, axes = plt.subplots(3, 3, figsize=(19, 14))
xb = np.arange(1, N_BINS + 1)

# (0,0) stacked tremor depth anomaly vs phase, with null band and bootstrap CI
ax = axes[0, 0]
ns = NULLSTACK["tremor"]["nsm"]
ax.fill_between(xb, np.nanpercentile(ns, 5, axis=0), np.nanpercentile(ns, 95, axis=0),
                color="0.75", alpha=0.6, label=f"null 5-95% ({ns.shape[0]} shifts)")
ax.plot(xb, np.nanmedian(ns, axis=0), color="0.45", ls="--", lw=1.2, label="null median")
o = STACK["tremor"]
for row in o["M"]:
    ax.plot(xb, row, color="tab:blue", alpha=0.16, lw=0.8)
ax.fill_between(xb, o["lo"], o["hi"], color="tab:blue", alpha=0.25, label="boot 90% CI")
ax.plot(xb, o["obs"], "o-", color="tab:blue", lw=2.4, label=f"tremor median (n={o['M'].shape[0]})")
s = STACK["slip"]
ax.fill_between(xb, s["lo"], s["hi"], color="tab:red", alpha=0.2)
ax.plot(xb, s["obs"], "s-", color="tab:red", lw=2.0, label=f"modeled slip depth (n={s['M'].shape[0]})")
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("event phase (fifths)")
ax.set_ylabel("depth anomaly vs event median (km)  + = deeper")
ax.set_title(f"A. within-event depth migration\ntremor $\\Delta$(1st$\\to$5th) = "
             f"{NULLSTACK['tremor']['obs_dz']:+.2f} km, p={NULLSTACK['tremor']['p_dz']:.3f}")
ax.legend(fontsize=7.5, loc="upper left")
ax.set_xticks(xb)
ax2 = ax.twinx()   # the intensity collapse that accompanies the deepening (count-bias control)
ax2.plot(xb, STACK["logcount"]["obs"], "^:", color="0.35", lw=1.6, ms=7,
         label="log10 daily count (right)")
ax2.fill_between(xb, STACK["logcount"]["lo"], STACK["logcount"]["hi"], color="0.4", alpha=0.15)
ax2.set_ylabel("log10 detections/day anomaly (dex)", color="0.35")
ax2.tick_params(axis="y", colors="0.35")
ax2.legend(fontsize=7.5, loc="lower left")

# (0,1) up-dip front p10 vs phase
ax = axes[0, 1]
ns = NULLSTACK["p10"]["nsm"]
ax.fill_between(xb, np.nanpercentile(ns, 5, axis=0), np.nanpercentile(ns, 95, axis=0),
                color="0.75", alpha=0.6, label="null 5-95%")
o = STACK["p10"]
for row in o["M"]:
    ax.plot(xb, row, color="tab:green", alpha=0.16, lw=0.8)
ax.fill_between(xb, o["lo"], o["hi"], color="tab:green", alpha=0.25, label="boot 90% CI")
ax.plot(xb, o["obs"], "o-", color="tab:green", lw=2.4, label="daily tremor p10 (up-dip front)")
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("event phase (fifths)")
ax.set_ylabel("p10 depth anomaly (km)")
ax.set_title(f"A. up-dip-most tremor per day\n$\\Delta$(1st$\\to$5th) = "
             f"{NULLSTACK['p10']['obs_dz']:+.2f} km, p={NULLSTACK['p10']['p_dz']:.3f}")
ax.legend(fontsize=7.5)
ax.set_xticks(xb)

# (0,2) rho(depth, count) histogram
ax = axes[0, 2]
rr = A.rho_depth_count.dropna().values
ax.hist(rr, bins=np.linspace(-1, 1, 17), color="tab:purple", alpha=0.8, edgecolor="k")
ax.axvline(0, color="k", lw=1)
ax.axvline(np.median(rr), color="tab:red", lw=2, label=f"median {np.median(rr):+.2f}")
ax.set_xlabel("Spearman(daily median depth, daily count)")
ax.set_ylabel("events")
ax.set_title(f"A. deeper or shallower when most intense?\n{int((rr < 0).sum())}/{rr.size} negative, "
             f"Wilcoxon p={wilcoxon(rr).pvalue:.3f}")
ax.legend(fontsize=8)

# (1,0) exponent distribution vs null
ax = axes[1, 0]
ee = B.exponent.dropna().values
alln = np.concatenate([[m] for m in B.null_exp_mean.dropna().values])
bins = np.linspace(min(ee.min(), alln.min()) - 0.2, max(ee.max(), alln.max()) + 0.2, 22)
ax.hist(alln, bins=bins, color="0.7", alpha=0.85, edgecolor="k", label="per-event null mean")
ax.hist(ee, bins=bins, color="tab:orange", alpha=0.8, edgecolor="k", label="observed")
ax.axvline(np.median(ee), color="tab:red", lw=2, label=f"median b = {np.median(ee):+.2f}")
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("scaling exponent b   [log10(C+1) on log10(R)]")
ax.set_ylabel("events")
ax.set_title("B. tremor-rate vs slip-rate exponent")
ax.legend(fontsize=8)

# (1,1) exponent by segment and epoch
ax = axes[1, 1]
groups, labels = [], []
for g in ("north", "central"):
    v = B[B.seg == g].exponent.dropna().values
    if v.size:
        groups.append(v)
        labels.append(f"{g}\nn={v.size}")
for g in sorted(B.epoch.unique()):
    v = B[B.epoch == g].exponent.dropna().values
    if v.size:
        groups.append(v)
        labels.append(f"{g}\nn={v.size}")
bp = ax.boxplot(groups, tick_labels=labels, showmeans=True, patch_artist=True)
for patch, c in zip(bp["boxes"], plt.cm.tab20.colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.6)
for i, v in enumerate(groups):
    ax.plot(np.full(v.size, i + 1) + RNG.normal(0, 0.05, v.size), v, "k.", ms=5, alpha=0.7)
ax.axhline(np.median(ee), color="tab:red", ls="--", lw=1, label="all-event median")
ax.set_ylabel("exponent b")
ax.set_title("B. exponent by segment and epoch")
ax.tick_params(axis="x", labelsize=7.5)
ax.legend(fontsize=8)

# (1,2) rising vs decaying (hysteresis)
ax = axes[1, 2]
m = np.isfinite(B.exponent_rise) & np.isfinite(B.exponent_decay)
ax.plot(B.exponent_rise[m], B.exponent_decay[m], "o", color="tab:blue")
lim = [min(B.exponent_rise[m].min(), B.exponent_decay[m].min()) - 0.3,
       max(B.exponent_rise[m].max(), B.exponent_decay[m].max()) + 0.3]
ax.plot(lim, lim, "k--", lw=1)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("exponent, RISING half")
ax.set_ylabel("exponent, DECAYING half")
dd = (B.exponent_rise - B.exponent_decay)[m]
ax.set_title(f"B. hysteresis  (n={int(m.sum())})\nrise-decay median {np.median(dd):+.2f}, "
             f"Wilcoxon p={wilcoxon(dd).pvalue:.3f}")

# (2,0) edges vs date: rel_p10 / rel_p90
ax = axes[2, 0]
for q, c, mk in (("rel_p10", "tab:green", "o"), ("rel_p90", "tab:brown", "s")):
    for g, mkr in (("north", "o"), ("central", "^")):
        sub = C[C.seg == g]
        ax.plot(sub.decyear, sub[q], mkr, color=c, ms=6, alpha=0.85,
                label=f"{q} {g}" if g == "north" else None)
    mfin = np.isfinite(C.decyear.values) & np.isfinite(C[q].values)
    if mfin.sum() >= 5:
        sl, ic = theilslopes(C[q].values[mfin], C.decyear.values[mfin], 0.95)[:2]
        xx = np.linspace(C.decyear.min(), C.decyear.max(), 10)
        ax.plot(xx, ic + sl * xx, color=c, lw=1.5, ls="--")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("year")
ax.set_ylabel("km, near-patch minus concurrent segment-wide")
ax.set_title("C. tremor edges vs CONCURRENT margin-wide tremor\n(circles north, triangles central)")
ax.legend(fontsize=8)

# (2,1) slip edges relative to segment all-time tremor edges
ax = axes[2, 1]
for q, c in (("slip_updip_rel", "tab:blue"), ("slip_downdip_rel", "tab:red")):
    for g, mkr in (("north", "o"), ("central", "^")):
        sub = C[C.seg == g]
        ax.plot(sub.decyear, sub[q], mkr, color=c, ms=6, alpha=0.85,
                label=f"{q} {g}" if g == "north" else None)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("year")
ax.set_ylabel("slip edge minus segment all-time tremor p10 / p90 (km)")
ax.set_title("C. slip edges vs the tremor band\n(- = slip up-dip of tremor; + = slip deeper)")
ax.legend(fontsize=8)

# (2,2) summary text
ax = axes[2, 2]
ax.axis("off")
lines = [
    "S9 SUMMARY",
    "",
    f"A  tremor depth 1st->5th fifth : {NULLSTACK['tremor']['obs_dz']:+.2f} km "
    f"(null p {NULLSTACK['tremor']['p_dz']:.3f})",
    f"   modeled slip depth          : {np.nanmedian(sdz):+.2f} km (mixture-limited)",
    f"   up-dip front (p10)          : {NULLSTACK['p10']['obs_dz']:+.2f} km "
    f"(null p {NULLSTACK['p10']['p_dz']:.3f})",
    f"   rho(depth, count) median    : {np.median(rr):+.2f} "
    f"({int((rr < 0).sum())}/{rr.size} negative)",
    "",
    f"B  exponent b median [IQR]     : {np.median(ee):+.2f} "
    f"[{np.percentile(ee, 25):+.2f},{np.percentile(ee, 75):+.2f}]",
    f"   pooled (event dummies) OLS  : {pooled_slope(POOLDF)[0]:+.2f}  "
    f"Theil-Sen {pooled_slope(POOLDF)[1]:+.2f}",
    f"   null b (per-event mean)     : {np.nanmedian(B.null_exp_mean):+.2f}",
    f"   rise - decay median         : {np.median(dd):+.2f} (p {wilcoxon(dd).pvalue:.3f})",
    "",
    f"C  slip up-dip - tremor p10    : {np.nanmedian(C.slip_updip_minus_near_p10):+.1f} km "
    f"(concurrent, near)",
    f"   slip down-dip - tremor p90  : {np.nanmedian(C.slip_downdip_minus_near_p90):+.1f} km",
    f"   rel_p10 median              : {np.nanmedian(C.rel_p10):+.2f} km",
    f"   rel_p90 median              : {np.nanmedian(C.rel_p90):+.2f} km",
    f"   min trend p over {nt} tests   : {TR[TR.n >= 5].p_perm.min():.3f} "
    f"(Bonferroni needs {0.05 / max(nt, 1):.4f})",
]
ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=9.5,
        transform=ax.transAxes)

fig.suptitle("S9 — within-event mechanics: depth migration (A), tremor-rate vs slip-rate scaling (B), "
             "edges vs concurrent tremor (C)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.975])
fig.savefig(OUT_FIG / "s9_within_event.png", dpi=130)
plt.close(fig)
log("wrote figures/s9_within_event.png")
log("S9 done.")
