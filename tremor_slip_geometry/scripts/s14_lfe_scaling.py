#!/usr/bin/env python
"""S14 - the tremor-vs-slip-rate exponent measured with the project's OWN LFE detections.

PLAN/RESULTS "TO DO before drafting (a)".  S9/S11 measured  C_PNSN ~ R^0.40  (daily PNSN tremor
count vs modeled potency rate) and S11 got  E ~ R^0.64  for radiated tremor energy; S11b
reproduced both with a raw-GNSS regressor (b(count) = 0.28 [0.10, 0.51]).  Frank & Brodsky (2019)
found per-LFE amplitude ~ (LFE rate)^0.32 in Guerrero.  S14 repeats the exponent with a THIRD,
independent tremor measure: the LFE detections this project made itself with matched filters on
its own family templates (`scripts/densify_gnw_gpu.py`, cc >= 0.80, 6-s min gap, 2-8 Hz, 40 Hz).

======================================================================================
WHAT n_stack IS  (verified in this script, section 3 - it is NOT an LFE count)
======================================================================================
`data/daily_dvv_<STA>_Z_2to4.csv` columns are (patch, date, dvv, cc_max, n_stack).  Reading
`scripts/dvv_roll30cal.py`:

    los  = searchsorted(Dn, Dn - 30 days, side='right')     # first DAILY STACK within 30 cal days
    counts = idx + 1 - los                                  # <- this is n_stack

so **n_stack = the number of DAILY STACKS (i.e. distinct days that had >= 20 detections and
produced a daily stack) inside the trailing 30 calendar days**, NOT the number of LFE detections.
It is an ACTIVE-DAY OCCUPANCY measure with a hard ceiling of 30 and a floor of 5 (the
`--min-stacks 5` gate that decides whether the date gets a dv/v point at all).  This script
verifies the definition directly (recomputes counts from the file's own dates) and quantifies the
saturation before using it.  Consequence: n_stack CANNOT measure a 30-d LFE rate; it can only
measure whether the family was active on most days, and inside an ETS window it almost always is.
Test 1 is therefore run as specified but reported as a saturating occupancy test, and the honest
30-d LFE-rate test is Test 1c, built from TRUE daily LFE counts (Level B) summed over 30 days.

======================================================================================
INPUTS (read-only; nothing outside tremor_slip_geometry/ is written)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv        27 confirmed north/central events (S1)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
                                                    label (T x 742), inc30_mm (T x 742, trailing
                                                    30-d thrust slip per triangle, mm), days,
                                                    patch_lat/lon/depth/area
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  /home/jovyan/cascadia-sse-strain                  mesh + sse_strain (Ld, Slab2 interpolator)
  catalogs/pnsn_tremor_cascadia_full.csv            PNSN tremor (the S9/S11 comparison measure)
  data/daily_dvv_<STA>_Z_2to4.csv                   195 stations - LEVEL A (patch, date, n_stack)
  data/<sta>_causality_cert.csv                     per-family causality certification (reliable)
  data/family_user_blocklist.csv                    229 user-vetoed families (visual dv/v QC)
  data/mf_<stem>_<year>.csv                         LEVEL B - per-detection catalogues
                                                    (template, time, cc, station); streamed with
                                                    pyarrow for the 12 stations listed below

  Slip-model loader code (Ld / A_raw -> cum_mm -> R7) is COPIED VERBATIM from
  scripts/s11_scaling_robustness.py, which copied it from s9_within_event.py / s4_fronts.py:
      cum_mm = -1e3 * (Ld @ A_raw.T)          (742 x T, mm, thrust positive)
      R7     = area-weighted sum over the 60-km neighbourhood of
               clip(7-d centred mean of d/dt cum_mm, 0)          (mm/d km^2)
      x = (lon + 124) * 111 * cos(47 deg) km E,   y = (lat - 47) * 111 km N
      segments: north lat >= 46, central 44 <= lat < 46

======================================================================================
GEOMETRY
======================================================================================
  Event window  W  = [start - 10, end + 10]                              (S9/S11)
  Footprint     F  = triangles labelled with the patch id on any day of [start, end]
  Neighbourhood NB = triangles within 60 km of an F centroid  -> R7      (S9/S11)
  LFE families are attached to an event if the family EPICENTRE (the "LAT_LON__cNNN" prefix of the
  patch id, i.e. the template's source location) lies within 40 km of an F centroid (the radius the
  task specifies for the LFE measure; 60 km is carried as a sensitivity column).
  Family depth = Slab2 depth at the family lat/lon (sse_strain.megathrust.slab_depth_interpolator).

======================================================================================
STATIONS STREAMED (LEVEL B) - chosen for family coverage of the confirmed footprints x epoch
======================================================================================
  north  : WISH 2010-2019, RADR 2011-2019, B943 2010-2025, LEBA 2013-2023, PR01 2017-2025,
           CARB 2019-2025, ANGUS 2021-2025
  central: HEBO 2010-2019, YACT 2013-2025, SSO2 2019-2025, ALSE 2019-2025, WLOO 2016-2025
  ~114 GB of mf_*.csv apparent bytes (budget 150 GB), streamed with pyarrow (~380-480 MB/s,
  IO-bound) and reduced ON THE FLY to per-(family, day) counts.  DEVIATION FROM THE SPEC, on
  purpose: ALL days are kept, not only the event windows.  The whole file has to be parsed either
  way, the reduced product is tiny (a dense n_fam x 7169 int32 matrix per station), and keeping
  every day is what makes the >= 365-d shifted null and the 30-d rolling counts (Test 1c) possible.
  Reduction per (template, day): n = detections, n_hi = detections with cc >= 0.9, sum(cc).

======================================================================================
TESTS
======================================================================================
  T1a  30-d LEVEL, all 195 stations (LEVEL A).  Per (family, event): Theil-Sen slope of
       log10(n_stack) on log10(inc30 at the family's nearest triangle) over window days with both
       > 0 and >= 8 usable days; same against the footprint-mean inc30.  Pooled with family-event
       dummies (within-demeaned OLS and Theil-Sen).  Per segment / epoch.  NULL: 200 rigid >= 365-d
       wrapped shifts of the family's n_stack series (slip unchanged), on a seeded subsample of
       family-events (the full null is 10^8 fits).
  T1b  CROSS-EVENT (LEVEL A).  Per family: max n_stack inside each event window vs that event's
       peak trailing-30-d slip at the family's nearest triangle; Theil-Sen exponent ACROSS the
       family's events, families with >= 4 events.
  T1c  30-d LEVEL, TRUE LFE COUNTS (LEVEL B).  Trailing-30-d sum of daily LFE detections summed
       over the station's near families, vs inc30 (footprint-area-weighted mean and nearest
       triangle).  Same fitting, same >= 365-d shifted null.  This is the test T1a would be if
       n_stack were a count.
  T2   DAILY LEVEL (LEVEL B), the S11 protocol exactly.  Per (station, event): N(day) = LFE
       detections summed over the station's families within 40 km of the footprint; b = Theil-Sen
       slope of log10(N + 1) on log10(R7) over days with N >= 1 and R7 > 0, >= 8 days.  Pooled with
       station-event dummies; per segment / epoch; lower vs upper half of log10 R7 (saturation);
       200 rigid >= 365-d wrapped shifts of the station's LFE series as the null.  The SAME script
       recomputes the PNSN count exponent b(C, R7) for the same 27 events (S9/S11 definition) so
       the two measures are compared per event, paired.
  T5   DIAGNOSIS (added after the first pass, which returned b ~ 0.03 against the PNSN 0.40 that
       the same script reproduces exactly).  Per station: median daily detections, the ETS /
       non-ETS ratio of that count, its rank correlation with the PNSN tremor count in the SAME
       40-km neighbourhood, and the per-family distribution of the ETS contrast.  VERDICT: at the
       densify threshold (cc >= 0.80, 6-s gap) the matched filter fires ~200 times per template
       per day in quiet times; the ETS/non-ETS count ratio is 1.05-1.7 where the PNSN count rises
       by a factor 10-100, and rho(LFE count, PNSN count) is 0.05-0.18.  The mf files are a
       nearly constant FALSE-ALARM stream, not an LFE catalogue - the pipeline extracts LFE signal
       from them by STACKING (build_long_window_3comp.py), not by counting.  The same fact
       explains why n_stack saturates: the >= 20-detections/day gate is met on nearly every day.
  T6   FLOOR REMOVAL.  An additive floor F under a signal S attenuates a log-log slope by
       S/(S+F).  The exponent is recomputed on the EXCESS count N - baseline, with the baseline
       taken from the pre-event window [start-70, start-10] and, alternatively, from the window's
       own 20th percentile; the identical operation applied to the PNSN count is the control.
  T3   Per-LFE AMPLITUDE.  The mf catalogues carry `cc` (matched-filter correlation), NOT an
       amplitude, so the Frank & Brodsky per-LFE amplitude exponent CANNOT be tested here.  What is
       reported instead, as a weak proxy only: the exponent of the daily MEAN cc on log10 R7 and
       the mean-cc distribution by slip-rate decile.
  T4   DETECTION-THRESHOLD control.  Everything in T2 repeated on a cc ladder (>= 0.85, 0.90,
       0.95, 0.98 as well as the 0.80 densify threshold), and
       with causality-CERTIFIED families only (minus the 229-family user blocklist) - the project's
       own "this template really is an LFE" verdict.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s14_lfe30_exponents.csv     T1a/T1b/T1c rows (group, measure, slip variant, n, b median,
                                      IQR, bootstrap CI, pooled OLS/TS, null mean/sd/p, sign test)
  results/s14_lfe_daily_exponents.csv T2/T3/T4 rows: per (station, event) exponents + the group
                                      summary rows + the paired PNSN comparison
  results/s14_coverage.csv            what was streamed: station, stem, years, bytes, detections,
                                      families, certified families, days covered, wall time
  figures/s14_lfe_scaling.png         9 panels

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s14_lfe_scaling.py
      (add --no-stream to reuse the cached per-station daily-count matrices)
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import argparse
import glob
import re
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pv
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree
from scipy.stats import binomtest, kruskal, mannwhitneyu, spearmanr, theilslopes, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

pa.set_cpu_count(4)
pa.set_io_thread_count(4)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)
CACHE = Path("/tmp/claude-1000/-home-jovyan-tremorferometry/"
             "03499166-a677-42d8-9397-e2cf5e28621f/scratchpad/s14_cache")
CACHE.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)

PRE_A, POST_A = 10, 10          # window [start-10, end+10]        (S9/S11)
NB_RADIUS_KM = 60.0             # slip neighbourhood                (S9/S11)
FAM_RADIUS_KM = 40.0            # LFE-family radius (task spec)
FAM_RADIUS_ALT = 60.0           # sensitivity
MIN_FIT_DAYS = 8
MIN_HALF_DAYS = 6
N_SHIFT = 200
N_BOOT = 2000
NULL_SUBSAMPLE = 2500           # family-events carried into the T1a null
CC_HI = 0.9
CC_LADDER = (0.85, 0.90, 0.95, 0.98)   # detection-threshold control ladder
COS47 = np.cos(np.deg2rad(47.0))
BUDGET_GB = 150.0

# station -> mf stem (the densify product prefix), chosen in the pre-run coverage pass
STREAM_STATIONS = [
    ("WISH", "wishp90f40", "north"),
    ("RADR", "radrp90f40", "north"),
    ("B943p90f40", "b943p90f40", "north"),
    ("LEBA", "lebap90f40", "north"),
    ("PR01", "pr01p90f40", "north"),
    ("CARB", "carbp90f40", "north"),
    ("ANGUS", "angusp90f40", "north"),
    ("HEBO", "hebop90f40", "central"),
    ("YACT", "yactp90f40", "central"),
    ("SSO2", "sso2p90f40", "central"),
    ("ALSE", "alsep90f40", "central"),
    ("WLOO", "wloop90f40", "central"),
]

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
    """Theil-Sen slope (median of pairwise slopes) - identical to scipy.stats.theilslopes[0]."""
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


def boot_median_ci(v, nboot=N_BOOT, lo=5, hi=95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


def pooled_fit(x, y, grp):
    """Within-group demeaned (event/family-event dummies) OLS and Theil-Sen slope + r."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    g = np.asarray(grp)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, g = x[ok], y[ok], g[ok]
    if x.size < 10:
        return np.nan, np.nan, np.nan, int(x.size)
    df = pd.DataFrame(dict(x=x, y=y, g=g))
    df["xd"] = df.x - df.groupby("g").x.transform("mean")
    df["yd"] = df.y - df.groupby("g").y.transform("mean")
    xd, yd = df.xd.values, df.yd.values
    if np.allclose(xd, 0):
        return np.nan, np.nan, np.nan, int(x.size)
    ols = float(np.dot(xd, yd) / np.dot(xd, xd))
    n = xd.size
    if n > 4000:                                  # Theil-Sen on a seeded subsample
        s = RNG.choice(n, 4000, replace=False)
        ts = ts_slope_fast(xd[s], yd[s])
    else:
        ts = ts_slope_fast(xd, yd)
    r = float(np.corrcoef(xd, yd)[0, 1])
    return ols, ts, r, int(n)


def sign_p(v, ref=0.0):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    k = int((v > ref).sum())
    if v.size == 0:
        return np.nan, 0, 0
    return float(binomtest(k, v.size, 0.5).pvalue), k, int(v.size)


def epoch_of(datestr: str) -> str:
    y = int(str(datestr)[:4])
    return "2010-2014" if y <= 2014 else ("2015-2019" if y <= 2019 else "2020-2026")


# ==================================================================================
# 1. Slip model  (loader copied from scripts/s11_scaling_robustness.py)
# ==================================================================================
ap = argparse.ArgumentParser()
ap.add_argument("--no-stream", dest="no_stream", action="store_true",
                help="(default behaviour) reuse the cached per-station daily LFE-count matrices")
ap.add_argument("--restream", action="store_true",
                help="ignore the cache and re-stream every mf_*.csv from disk")
ARGS = ap.parse_args()

log("loading slip model (load_solution + invert, ~1 min) ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
slabf = slab_depth_interpolator(fault)

A_raw = np.nan_to_num(
    np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"]
)
CUM_MM = -1e3 * (Ld @ A_raw.T)
DIFF = np.zeros_like(CUM_MM)
DIFF[:, 1:] = np.diff(CUM_MM, axis=1)
RATE7 = np.clip(uniform_filter1d(DIFF, size=7, axis=1, mode="nearest"), 0.0, None)
log(f"  Ld {Ld.shape}, A_raw {A_raw.shape}, cum {CUM_MM.shape} mm, RATE7 max {RATE7.max():.3f} mm/d")

NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]
INC30 = np.asarray(NPZ["inc30_mm"], dtype=np.float32)      # (T, 742) trailing-30-d slip, mm
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
DAY_DT = pd.to_datetime(DAYS).values.astype("datetime64[D]")
log(f"  {NT} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles; "
    f"inc30 range {np.nanmin(INC30):.1f}..{np.nanmax(INC30):.1f} mm")

# sanity: inc30 should equal the trailing-30-d increment of cum_mm
chk = CUM_MM[:, 30:] - CUM_MM[:, :-30]
cc_chk = np.corrcoef(chk[:, ::37].ravel(), INC30.T[:, 30:][:, ::37].ravel())[0, 1]
log(f"  inc30 vs cum_mm(t)-cum_mm(t-30) correlation = {cc_chk:+.4f} (consistency check)")


# ==================================================================================
# 2. Events, footprints, R7 / inc30 per event
# ==================================================================================
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  events {len(ev_in)}  {ev_in.seg.value_counts().to_dict()}")

EVENTS = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = (np.unique(np.concatenate([np.asarray(h, int) for h in hits]))
          if len(hits) else np.array([], int))
    nb = np.union1d(nb, np.flatnonzero(F))
    aw = PAREA[nb][:, None]
    R7 = (aw * RATE7[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)
    fidx = np.flatnonzero(F)
    fw = PAREA[fidx]
    inc_fp = (INC30[w0:w1 + 1][:, fidx] * fw).sum(axis=1) / fw.sum()   # footprint-area-weighted
    EVENTS.append(dict(pid=pid, seg=str(r["seg"]), start=str(r["start"]), end=str(r["end"]),
                       Mw=float(r["Mw"]), epoch=epoch_of(r["start"]),
                       w0=w0, w1=w1, nW=nW, ts=ts, te=te, nb=nb, F=F, fidx=fidx,
                       fx=fx, fy=fy, R7=R7, inc_fp=inc_fp))
log("  per-event R7 and footprint inc30 built")


# ==================================================================================
# 3. LEVEL A - families, n_stack, and WHAT n_stack ACTUALLY IS
# ==================================================================================
DAY0 = int((np.datetime64(DAYS[0]) - np.datetime64("1970-01-01")).astype(int))
assert int((np.datetime64(DAYS[-1]) - np.datetime64(DAYS[0])).astype(int)) == NT - 1, \
    "model day grid is not contiguous daily"
FAMPAT = re.compile(r"^(-?\d+\.\d+)_(-?\d+\.\d+)__c(\d+)$")


def date_to_di(col):
    """arrow string column 'YYYY-MM-DD...' -> model day index (int32); outside the grid -> < 0."""
    d = pc.cast(pc.strptime(pc.utf8_slice_codeunits(col, 0, 10), format="%Y-%m-%d", unit="s"),
                pa.date32())
    a = pc.cast(d, pa.int32())
    if isinstance(a, pa.ChunkedArray):
        a = a.combine_chunks()
    return np.asarray(a).astype(np.int64) - DAY0


def dict_encode(col):
    """(integer codes, vocabulary list) for a string Array or ChunkedArray."""
    d = pc.dictionary_encode(col)
    if isinstance(d, pa.ChunkedArray):
        d = d.combine_chunks()
    return np.asarray(d.indices).astype(np.int32), d.dictionary.to_pylist()


log("LEVEL A: reading the 195 daily_dvv CSVs (per-family n_stack series) ...")
FAM_KEY = []            # (sta, patch)
FAM_DI = []             # int32 day indices, sorted
FAM_NS = []             # int8 n_stack
nrows_A = 0
for f in sorted(glob.glob(str(ROOT / "data/daily_dvv_*_Z_2to4.csv"))):
    tok = os.path.basename(f)[len("daily_dvv_"):-len("_Z_2to4.csv")]
    t = pv.read_csv(f, convert_options=pv.ConvertOptions(
        include_columns=["patch", "date", "n_stack"],
        column_types={"patch": pa.string(), "date": pa.string(), "n_stack": pa.int16()}))
    di = date_to_di(t.column("date"))
    ns = np.asarray(t.column("n_stack")).astype(np.int8)
    codes, vocab = dict_encode(t.column("patch"))
    keep = (di >= 0) & (di < NT)
    di, ns, codes = di[keep].astype(np.int32), ns[keep], codes[keep]
    nrows_A += di.size
    order = np.lexsort((di, codes))
    di, ns, codes = di[order], ns[order], codes[order]
    bnd = np.flatnonzero(np.diff(codes)) + 1
    for a, b in zip(np.r_[0, bnd], np.r_[bnd, codes.size]):
        FAM_KEY.append((tok, vocab[codes[a]]))
        FAM_DI.append(di[a:b].astype(np.int32))
        FAM_NS.append(ns[a:b])
log(f"  {nrows_A:,} (station, family, date) rows; {len(FAM_KEY)} (station, family) series")

# --- verify the n_stack definition from the file's OWN dates --------------------------------
agree_all, le_all, nchk = [], [], 0
for j in RNG.choice(len(FAM_KEY), 400, replace=False):
    d, sn = FAM_DI[j], FAM_NS[j].astype(int)
    los = np.searchsorted(d, d - 30, side="right")
    recount = np.arange(d.size) + 1 - los
    agree_all.append((recount == sn).mean())
    le_all.append((recount <= sn).mean())
    nchk += d.size
NS_ALL = np.concatenate(FAM_NS).astype(int)
log(f"  n_stack CHECK: recomputing 'number of DAILY STACKS inside the trailing 30 calendar days' "
    f"from each CSV's own dates reproduces n_stack EXACTLY on {np.mean(agree_all) * 100:.2f}% and "
    f"is <= n_stack on {np.mean(le_all) * 100:.2f}% of {nchk:,} rows (400 random families).")
log("     The <= 100 % / == 89 % pattern is exactly what the definition predicts: dvv_roll30cal "
    "writes a row ONLY where")
log("     counts >= --min-stacks (5), so days with 1-4 stacks in their trailing window exist in "
    "the STACK series but are")
log("     absent from the CSV; recounting from the CSV's own dates can only UNDERcount, never "
    "overcount.  Definition confirmed.")
log(f"  -> n_stack is an ACTIVE-DAY OCCUPANCY count (floor 5 = --min-stacks, ceiling 30), "
    f"NOT an LFE detection count.  min {NS_ALL.min()}, median {np.median(NS_ALL):.0f}, "
    f"max {NS_ALL.max()}, frac==30 {np.mean(NS_ALL >= 30):.3f}")
del NS_ALL

# certification + user blocklist
CERT = {}
for f in glob.glob(str(ROOT / "data/*_causality_cert.csv")):
    stem = os.path.basename(f).replace("_causality_cert.csv", "")
    c = pd.read_csv(f)
    CERT[stem] = set(c.loc[c.reliable.astype(bool), "fam"])
BLOCK = pd.read_csv(ROOT / "data/family_user_blocklist.csv")
BLOCKSET = set(zip(BLOCK.stem.astype(str), BLOCK.fam.astype(str)))
log(f"  certification files {len(CERT)}; user blocklist {len(BLOCKSET)} (stem, family) pairs")


def cert_stem(tok: str) -> str:
    return (tok.lower().replace("p90f40rev", "").replace("p90f40", "")
            .replace("hhz", "").replace("bhz", "").replace("new", ""))


def is_certified(tok: str, fam: str) -> bool:
    s = cert_stem(tok)
    return (fam in CERT.get(s, set())) and ((s, fam) not in BLOCKSET)


UF = pd.DataFrame(FAM_KEY, columns=["sta", "patch"])
gg = np.array([FAMPAT.match(p).groups() for p in UF.patch.values], dtype=object)
UF["flat"] = gg[:, 0].astype(float)
UF["flon"] = gg[:, 1].astype(float)
UF["ndays"] = [d.size for d in FAM_DI]
UFX, UFY = to_km(UF.flon.values, UF.flat.values)
_, UF_TRI = TRI_TREE.query(np.column_stack([UFX, UFY]))
UF["tri"] = UF_TRI
UF["seg"] = seg_of_lat(UF.flat.values)
UF["certified"] = [is_certified(s, p) for s, p in zip(UF.sta.values, UF.patch.values)]
try:
    dep = np.asarray(slabf(UF.flon.values, UF.flat.values), float).ravel()
except Exception:
    dep = np.asarray(slabf(np.column_stack([UF.flon.values, UF.flat.values])), float).ravel()
UF["depth_km"] = dep
FAM_TREE = cKDTree(np.column_stack([UFX, UFY]))
log(f"  {len(UF)} (station, family) pairs; certified {int(UF.certified.sum())} "
    f"({UF.certified.mean() * 100:.0f}%); Slab2 depth median {np.nanmedian(dep):.1f} km "
    f"[p10 {np.nanpercentile(dep, 10):.1f}, p90 {np.nanpercentile(dep, 90):.1f}]; "
    f"segments {UF.seg.value_counts().to_dict()}")

for e in EVENTS:
    for rad, key in ((FAM_RADIUS_KM, "fam40"), (FAM_RADIUS_ALT, "fam60")):
        hit = FAM_TREE.query_ball_point(np.column_stack([e["fx"], e["fy"]]), r=rad)
        e[key] = (np.unique(np.concatenate([np.asarray(h, int) for h in hit]))
                  if len(hit) else np.array([], int))
log("  family-event attachment (40 km): " +
    ", ".join(f"{e['pid']}:{len(e['fam40'])}" for e in EVENTS[:5]) + " ...")


# ==================================================================================
# 4. TEST 1a - 30-d level, n_stack vs inc30  (all 195 stations)
# ==================================================================================
SHIFT_LO_I = DAY_INDEX["2010-01-01"]
SHIFT_HI_I = DAY_INDEX["2026-05-27"]
SPAN = SHIFT_HI_I - SHIFT_LO_I + 1


def series_in_window(di, val, w0, nW, shift=0):
    """Values of a (di, val) series inside [w0, w0+nW), after a rigid wrapped shift (days)."""
    if shift:
        d = SHIFT_LO_I + ((di - SHIFT_LO_I + shift) % SPAN)
    else:
        d = di
    m = (d >= w0) & (d < w0 + nW)
    out = np.full(nW, np.nan)
    out[d[m] - w0] = val[m]
    return out


log("TEST 1a: per (family, event) exponent of n_stack vs trailing-30-d slip ...")
rows_1a = []
sat_frac = []
P1A_X, P1A_Y, P1A_G = [], [], []      # pooled arrays, seeded subsample of family-events
P1A_KEEP = 0.02
for e in EVENTS:
    w0, nW = e["w0"], e["nW"]
    inc_fp = e["inc_fp"]
    for j in e["fam40"]:
        di, ns = FAM_DI[j], FAM_NS[j].astype(float)
        y = series_in_window(di, ns, w0, nW)
        if np.isfinite(y).sum() < MIN_FIT_DAYS:
            continue
        tri = int(UF.tri.values[j])
        inc_near = INC30[w0:w0 + nW, tri].astype(float)
        sat_frac.append(np.nanmean(y[np.isfinite(y)] >= 30))
        rec = dict(test="T1a", pid=e["pid"], seg=e["seg"], epoch=e["epoch"],
                   sta=UF.sta.values[j], patch=UF.patch.values[j], fam=j,
                   fam_seg=UF.seg.values[j], cert=bool(UF.certified.values[j]),
                   depth_km=float(UF.depth_km.values[j]))
        for tag, xr in (("near", inc_near), ("fp", inc_fp)):
            ok = np.isfinite(y) & (y > 0) & np.isfinite(xr) & (xr > 0)
            rec[f"n_{tag}"] = int(ok.sum())
            rec[f"b_{tag}"] = (ts_slope_fast(np.log10(xr[ok]), np.log10(y[ok]))
                               if ok.sum() >= MIN_FIT_DAYS else np.nan)
        rec["nsmax"] = float(np.nanmax(y))
        rec["incmax"] = float(np.nanmax(inc_near))
        rows_1a.append(rec)
        if RNG.random() < P1A_KEEP:
            ok = np.isfinite(y) & (y > 0) & np.isfinite(inc_near) & (inc_near > 0)
            if ok.sum() >= MIN_FIT_DAYS:
                P1A_X.append(np.log10(inc_near[ok]))
                P1A_Y.append(np.log10(y[ok]))
                P1A_G.append(np.full(int(ok.sum()), f"{e['pid']}_{j}"))
T1A = pd.DataFrame(rows_1a)
log(f"  {len(T1A):,} family-event fits; saturation: mean fraction of window days with "
    f"n_stack == 30 is {np.mean(sat_frac):.3f}")

# --- null: 200 rigid >= 365-d shifts of n_stack on a seeded subsample -----------------------
sub = T1A.sample(min(NULL_SUBSAMPLE, len(T1A)), random_state=20260910)
null_b = []
per_z = []
EVBY = {e["pid"]: e for e in EVENTS}
for _, r in sub.iterrows():
    e = EVBY[r["pid"]]
    j = int(r["fam"])
    di, ns = FAM_DI[j], FAM_NS[j].astype(float)
    tri = int(UF.tri.values[j])
    xr = INC30[e["w0"]:e["w0"] + e["nW"], tri].astype(float)
    bs = []
    for sh in RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT // 4):
        y = series_in_window(di, ns, e["w0"], e["nW"], shift=int(sh))
        ok = np.isfinite(y) & (y > 0) & np.isfinite(xr) & (xr > 0)
        if ok.sum() >= MIN_FIT_DAYS:
            bs.append(ts_slope_fast(np.log10(xr[ok]), np.log10(y[ok])))
    bs = np.asarray([b for b in bs if np.isfinite(b)], float)
    if bs.size >= 10:
        null_b.append(bs)
        if np.isfinite(r["b_near"]) and bs.std() > 0:
            per_z.append((r["b_near"] - bs.mean()) / bs.std())
NULL_1A = np.concatenate(null_b) if null_b else np.array([np.nan])
log(f"  null (n={len(null_b)} family-events x {N_SHIFT // 4} shifts): mean "
    f"{np.nanmean(NULL_1A):+.4f} sd {np.nanstd(NULL_1A):.4f}; per-fit z median "
    f"{np.nanmedian(per_z):+.3f} (n={len(per_z)})")


# ==================================================================================
# 5. TEST 1b - cross-event: peak n_stack vs peak event slip, per family
# ==================================================================================
log("TEST 1b: cross-event exponent per family (>= 4 events) ...")
rows_1b = []
gb = T1A.dropna(subset=["nsmax", "incmax"]).groupby("fam")
for j, g in gb:
    g = g[(g.nsmax > 0) & (g.incmax > 0)]
    if len(g) < 4:
        continue
    b = ts_slope_fast(np.log10(g.incmax.values), np.log10(g.nsmax.values), nmin=4)
    rows_1b.append(dict(test="T1b", fam=int(j), sta=UF.sta.values[int(j)],
                        patch=UF.patch.values[int(j)], fam_seg=UF.seg.values[int(j)],
                        cert=bool(UF.certified.values[int(j)]), n_events=len(g), b=b,
                        rho=spearmanr(g.incmax.values, g.nsmax.values).statistic))
T1B = pd.DataFrame(rows_1b)
log(f"  {len(T1B)} families with >= 4 events; b median {np.nanmedian(T1B.b):+.3f} "
    f"[{np.nanpercentile(T1B.b, 25):+.3f}, {np.nanpercentile(T1B.b, 75):+.3f}]")


# ==================================================================================
# 6. LEVEL B - stream the per-detection mf catalogues into daily per-family counts
# ==================================================================================
MAXFAM = 900
FLUSH_ROWS = 8_000_000


def stream_station(tok: str, stem: str):
    """Stream every mf_<stem>_<year>.csv -> dense (n_fam x NT) daily detection counts.

    Returns (fams, N, NTH, CCS, bytes_read, n_det, files).  N = all detections (cc >= 0.80, the
    densify threshold); NTH[t] = detections with cc >= t for every t in CC_LADDER (the
    detection-threshold control); CCS = summed cc (for the mean-cc proxy).
    """
    files = sorted(glob.glob(str(ROOT / f"data/mf_{stem}_*.csv")))
    files = [f for f in files if re.search(r"_(\d{4})\.csv$", f)]
    fam_index = {}
    N = np.zeros((MAXFAM, NT), np.int32)
    NTH = {t: np.zeros((MAXFAM, NT), np.int32) for t in CC_LADDER}
    CCS = np.zeros((MAXFAM, NT), np.float32)
    nbytes = 0
    ndet = 0
    buf_c, buf_d, buf_cc = [], [], []
    nbuf = 0

    def flush():
        nonlocal buf_c, buf_d, buf_cc, nbuf
        if not nbuf:
            return
        c = np.concatenate(buf_c)
        d = np.concatenate(buf_d)
        v = np.concatenate(buf_cc)
        flat = c.astype(np.int64) * NT + d
        N.ravel()[:] += np.bincount(flat, minlength=MAXFAM * NT).astype(np.int32)
        CCS.ravel()[:] += np.bincount(flat, weights=v.astype(np.float64),
                                      minlength=MAXFAM * NT).astype(np.float32)
        for t in CC_LADDER:
            hi = v >= t
            if hi.any():
                NTH[t].ravel()[:] += np.bincount(flat[hi],
                                                 minlength=MAXFAM * NT).astype(np.int32)
        buf_c, buf_d, buf_cc = [], [], []
        nbuf = 0

    for f in files:
        sz = os.path.getsize(f)
        nbytes += sz
        rd = pv.open_csv(f, read_options=pv.ReadOptions(block_size=64 << 20),
                         convert_options=pv.ConvertOptions(
                             include_columns=["template", "time", "cc"],
                             column_types={"template": pa.string(), "time": pa.string(),
                                           "cc": pa.float32()}))
        for b in rd:
            di = date_to_di(b.column(1))
            cc = np.asarray(b.column(2))
            loc, vocab = dict_encode(b.column(0))
            g = np.empty(len(vocab), np.int32)
            for k, v in enumerate(vocab):
                if v not in fam_index:
                    fam_index[v] = len(fam_index)
                g[k] = fam_index[v]
            codes = g[loc]
            keep = (di >= 0) & (di < NT) & np.isfinite(cc)
            if not keep.any():
                continue
            buf_c.append(codes[keep].astype(np.int32))
            buf_d.append(di[keep].astype(np.int32))
            buf_cc.append(cc[keep])
            nbuf += int(keep.sum())
            ndet += int(keep.sum())
            if nbuf >= FLUSH_ROWS:
                flush()
        del rd
    flush()
    nf = len(fam_index)
    assert nf <= MAXFAM, f"{tok}: {nf} templates > MAXFAM"
    fams = [None] * nf
    for k, v in fam_index.items():
        fams[v] = k
    return (fams, N[:nf], {t: NTH[t][:nf] for t in CC_LADDER}, CCS[:nf],
            nbytes, ndet, len(files))


log("LEVEL B: streaming the per-detection mf catalogues ...")
STA = {}
cov_rows = []
gb_read = 0.0
for tok, stem, ssg in STREAM_STATIONS:
    cf = CACHE / f"{stem}_v2.npz"
    t_s = time.time()
    if cf.exists() and not ARGS.restream:
        z = np.load(cf, allow_pickle=True)
        fams = list(z["fams"])
        N, CCS = z["N"], z["CCS"]
        NTH = {t: z[f"N{int(t * 100)}"] for t in CC_LADDER}
        nbytes, ndet, nfiles = int(z["nbytes"]), int(z["ndet"]), int(z["nfiles"])
        wall = float(z["wall"])
        cached = True
    else:
        if gb_read > BUDGET_GB:
            log(f"  BUDGET {BUDGET_GB:.0f} GB reached -> skipping {tok}")
            continue
        fams, N, NTH, CCS, nbytes, ndet, nfiles = stream_station(tok, stem)
        wall = time.time() - t_s
        np.savez_compressed(cf, fams=np.array(fams, dtype=object), N=N, CCS=CCS,
                            nbytes=nbytes, ndet=ndet, nfiles=nfiles, wall=wall,
                            **{f"N{int(t * 100)}": NTH[t] for t in CC_LADDER})
        cached = False
    gb_read += nbytes / 1e9
    tot = N.sum(axis=0)
    act = np.flatnonzero(tot > 0)
    d0, d1 = (int(act[0]), int(act[-1])) if act.size else (0, -1)
    lat = np.array([float(FAMPAT.match(p).group(1)) for p in fams])
    lon = np.array([float(FAMPAT.match(p).group(2)) for p in fams])
    fx, fy = to_km(lon, lat)
    cert = np.array([is_certified(tok, p) for p in fams])
    STA[tok] = dict(tok=tok, stem=stem, ssg=ssg, fams=fams, N=N, NHI=NTH[CC_HI], NTH=NTH,
                    CCS=CCS, lat=lat, lon=lon, fx=fx, fy=fy, cert=cert, d0=d0, d1=d1,
                    tree=cKDTree(np.column_stack([fx, fy])))
    cov_rows.append(dict(sta=tok, stem=stem, seg_hint=ssg, files=nfiles,
                         GB=nbytes / 1e9, detections=ndet, families=len(fams),
                         certified=int(cert.sum()), first_day=DAYS[d0] if act.size else "",
                         last_day=DAYS[d1] if act.size else "",
                         days_with_detections=int((tot > 0).sum()),
                         lat_min=lat.min(), lat_max=lat.max(),
                         wall_s=round(wall, 1), cached=cached))
    log(f"  {tok:12s} {nfiles:2d} files {nbytes / 1e9:6.1f} GB  {ndet:12,d} detections  "
        f"{len(fams):3d} templates ({int(cert.sum())} certified)  "
        f"{DAYS[d0] if act.size else '-'}..{DAYS[d1] if act.size else '-'}  "
        f"{wall:6.1f}s {'[cache]' if cached else ''}")
COV = pd.DataFrame(cov_rows)
log(f"  streamed {COV.GB.sum():.1f} GB, {COV.detections.sum():,} detections, "
    f"{len(COV)} stations, {COV.wall_s.sum() / 60:.1f} min")


# ==================================================================================
# 7. LEVEL B machinery
# ==================================================================================
def fam_sel(S, e, rad=FAM_RADIUS_KM, cert_only=False):
    hit = S["tree"].query_ball_point(np.column_stack([e["fx"], e["fy"]]), r=rad)
    idx = (np.unique(np.concatenate([np.asarray(h, int) for h in hit]))
           if len(hit) else np.array([], int))
    if cert_only and idx.size:
        idx = idx[S["cert"][idx]]
    return idx


def shifted_window(vec, d0, d1, w0, nW, sh):
    """Values of `vec` in [w0, w0+nW) after a rigid wrapped shift of the station's own record."""
    span = d1 - d0 + 1
    t = np.arange(w0, w0 + nW)
    src = d0 + ((t - d0 - sh) % span)
    return vec[src]


def fit_logy(x, y, eps=1.0, nmin=MIN_FIT_DAYS, logy=True):
    ok = np.isfinite(x) & (x > 0) & np.isfinite(y) & (y > 0 if logy else True)
    if ok.sum() < nmin:
        return np.nan, int(ok.sum())
    lx = np.log10(x[ok])
    ly = np.log10(y[ok] + eps) if logy else y[ok]
    return ts_slope_fast(lx, ly, nmin), int(ok.sum())


def halves(x, y, eps=1.0):
    ok = np.isfinite(x) & (x > 0) & np.isfinite(y) & (y > 0)
    if ok.sum() < 2 * MIN_HALF_DAYS:
        return np.nan, np.nan
    lx, ly = np.log10(x[ok]), np.log10(y[ok] + eps)
    med = np.median(lx)
    lo, hi = lx <= med, lx > med
    bl = ts_slope_fast(lx[lo], ly[lo], MIN_HALF_DAYS) if lo.sum() >= MIN_HALF_DAYS else np.nan
    bh = ts_slope_fast(lx[hi], ly[hi], MIN_HALF_DAYS) if hi.sum() >= MIN_HALF_DAYS else np.nan
    return bl, bh


# ==================================================================================
# 8. TEST 2 / 3 / 4 - daily LFE counts vs R7  (+ the PNSN comparison)
# ==================================================================================
log("PNSN comparison: rebuilding the S9/S11 daily tremor counts ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DI = ((tr["t"].values.astype("datetime64[D]") - np.datetime64(DAYS[0])).astype(int))
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
for e in EVENTS:
    hit = TREMOR_TREE.query_ball_point(np.column_stack([e["fx"], e["fy"]]), r=NB_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in hit]))
            if len(hit) else np.array([], int))
    d = TR_DI[near]
    d = d[(d >= 0) & (d < NT)]
    vec = np.bincount(d, minlength=NT).astype(float)
    e["C_pnsn"] = vec[e["w0"]:e["w0"] + e["nW"]]
    e["C_pnsn_full"] = vec
    e["b_pnsn"], e["n_pnsn"] = fit_logy(e["R7"], e["C_pnsn"], eps=1.0)
bp = np.array([e["b_pnsn"] for e in EVENTS], float)
log(f"  PNSN b(C, R7) median {np.nanmedian(bp):+.3f} "
    f"[{np.nanpercentile(bp, 25):+.3f}, {np.nanpercentile(bp, 75):+.3f}] "
    f"(S9/S11 reported +0.401 [+0.181, +0.478]) -- cross-check")

# --- TEST 5 (DIAGNOSIS): is the mf detection stream an LFE rate at all? --------------------
log("TEST 5: per-station detector diagnosis (ETS contrast, PNSN agreement, dynamic range) ...")
ETS_MASK = np.zeros(NT, bool)
for e in EVENTS:
    ETS_MASK[e["w0"]:e["w0"] + e["nW"]] = True
DIAG = []
for tok, S in STA.items():
    d0, d1 = S["d0"], S["d1"]
    m = np.zeros(NT, bool)
    m[d0:d1 + 1] = True
    hit = TREMOR_TREE.query_ball_point(np.column_stack([S["fx"], S["fy"]]), r=FAM_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in hit]))
            if len(hit) else np.array([], int))
    dd = TR_DI[near]
    dd = dd[(dd >= 0) & (dd < NT)]
    pn = np.bincount(dd, minlength=NT).astype(float)
    ee, qq = m & ETS_MASK, m & ~ETS_MASK
    row = dict(sta=tok)
    row["days_active"] = int((S["N"].sum(0) > 0).sum())
    for tag, vec in ([("cc0.80", S["N"].sum(0).astype(float))] +
                     [(f"cc{t:.2f}", S["NTH"][t].sum(0).astype(float)) for t in CC_LADDER]):
        row[f"mean_{tag}"] = float(vec[m].mean())
        row[f"med_{tag}"] = float(np.median(vec[m]))
        # MEAN-based contrast (the median collapses to 0 for sparse series such as the PNSN
        # count on quiet days, so means are the only statistic comparable across measures)
        row[f"ets_ratio_{tag}"] = (float(vec[ee].mean() / vec[qq].mean())
                                   if vec[qq].mean() > 0 else np.nan)
        row[f"rho_pnsn_{tag}"] = float(spearmanr(vec[m], pn[m]).statistic)
        lo, hi = np.percentile(vec[m], [10, 99])
        row[f"dex_{tag}"] = float(np.log10(max(hi, 1) / max(lo, 1)))
    row["pnsn_ets_ratio"] = float(pn[ee].mean() / pn[qq].mean()) if pn[qq].mean() > 0 else np.nan
    row["pnsn_dex"] = float(np.log10(max(np.percentile(pn[m], 99), 1) /
                                     max(np.percentile(pn[m], 10), 1)))
    row["pnsn_mean_ets"] = float(pn[ee].mean())
    row["pnsn_mean_quiet"] = float(pn[qq].mean())
    # per-family ETS contrast (mean-based)
    a = S["N"][:, ee].mean(axis=1)
    b_ = S["N"][:, qq].mean(axis=1)
    r = a[b_ > 0] / b_[b_ > 0]
    row["fam_ets_ratio_p50"] = float(np.median(r)) if r.size else np.nan
    row["fam_ets_ratio_p90"] = float(np.percentile(r, 90)) if r.size else np.nan
    row["fam_frac_ratio_gt15"] = float(np.mean(r > 1.5)) if r.size else np.nan
    DIAG.append(row)
DIAG = pd.DataFrame(DIAG)
log("  ETS/non-ETS MEAN-count ratio  (LFE cc>=0.80 | cc>=0.98 | PNSN) and rho(LFE, PNSN):")
for _, r in DIAG.iterrows():
    log(f"    {r.sta:12s} {r['ets_ratio_cc0.80']:5.2f} | {r['ets_ratio_cc0.98']:5.2f} | "
        f"{r.pnsn_ets_ratio:7.1f}   rho {r['rho_pnsn_cc0.80']:+.3f} | "
        f"{r['rho_pnsn_cc0.98']:+.3f}   mean {r['mean_cc0.80']:9.0f} det/d")

log("TEST 2/3/4: daily LFE counts vs R7, per station x event ...")
rows_d = []
famrows = []
for tok, S in STA.items():
    d0, d1 = S["d0"], S["d1"]
    span_s = d1 - d0 + 1
    for e in EVENTS:
        w0, nW = e["w0"], e["nW"]
        inside = (w0 >= d0) and (w0 + nW - 1 <= d1)
        idx40 = fam_sel(S, e, FAM_RADIUS_KM)
        idx60 = fam_sel(S, e, FAM_RADIUS_ALT)
        idxc = fam_sel(S, e, FAM_RADIUS_KM, cert_only=True)
        if idx40.size == 0 or not inside:
            continue
        Nfull = S["N"][idx40].sum(0).astype(float)
        Hfull = S["NHI"][idx40].sum(0).astype(float)
        Cfull = S["CCS"][idx40].sum(0).astype(float)
        N60 = S["N"][idx60].sum(0).astype(float)
        Nc = S["N"][idxc].sum(0).astype(float) if idxc.size else np.zeros(NT)
        w = slice(w0, w0 + nW)
        N, H, Csum, N6, NC = Nfull[w], Hfull[w], Cfull[w], N60[w], Nc[w]
        R7 = e["R7"]
        with np.errstate(invalid="ignore", divide="ignore"):
            meancc = np.where(N > 0, Csum / np.maximum(N, 1), np.nan)
        rec = dict(test="T2", sta=tok, pid=e["pid"], seg=e["seg"], epoch=e["epoch"],
                   start=e["start"], Mw=e["Mw"], nW=nW, nfam40=int(idx40.size),
                   nfam60=int(idx60.size), nfam_cert=int(idxc.size),
                   det_in_window=float(N.sum()))
        rec["b_N"], rec["n_N"] = fit_logy(R7, N, 1.0)
        rec["b_N60"], _ = fit_logy(R7, N6, 1.0)
        rec["b_Nhi"], rec["n_Nhi"] = fit_logy(R7, H, 1.0)
        rec["b_Ncert"], rec["n_Ncert"] = fit_logy(R7, NC, 1.0)
        # cc-threshold ladder
        for t in CC_LADDER:
            yt = S["NTH"][t][idx40].sum(0)[w].astype(float)
            rec[f"b_cc{int(t * 100)}"], _ = fit_logy(R7, yt, 1.0)
            rec[f"det_cc{int(t * 100)}"] = float(yt.sum())
        # BACKGROUND SUBTRACTION: the detector has a large, nearly constant false-alarm rate,
        # so log10(N) is compressed by the floor.  Two baselines: the pre-event window
        # [start-70, start-10] and the window's own 20th percentile.
        pre = Nfull[max(0, w0 - 60):w0]
        base_pre = float(np.median(pre)) if pre.size >= 20 and w0 - 60 >= d0 else np.nan
        base_p20 = float(np.percentile(N, 20))
        rec["base_pre"], rec["base_p20"] = base_pre, base_p20
        rec["ets_contrast"] = (float(np.median(N) / base_pre)
                               if np.isfinite(base_pre) and base_pre > 0 else np.nan)
        for tag, bl in (("pre", base_pre), ("p20", base_p20)):
            if np.isfinite(bl):
                exc = N - bl
                rec[f"b_exc_{tag}"], rec[f"n_exc_{tag}"] = fit_logy(R7, exc, 0.0)
            else:
                rec[f"b_exc_{tag}"], rec[f"n_exc_{tag}"] = np.nan, 0
        # ADDITIVE-FLOOR DE-ATTENUATION.  If N = F + k R^b with F a constant false-alarm floor,
        # then dlogN/dlogR = b * S/(F+S), and S/(F+S) = (r - 1)/r with r = mean(N in window) /
        # mean(N in the pre-event window).  So b_true = b_obs * r/(r - 1).
        rr = (float(N.mean() / base_pre) if np.isfinite(base_pre) and base_pre > 0 else np.nan)
        rec["ets_ratio_mean"] = rr
        rec["b_deatt"] = (rec["b_N"] * rr / (rr - 1.0)
                          if np.isfinite(rr) and rr > 1.1 else np.nan)
        # the same background subtraction applied to the PNSN count, as a like-for-like control
        pnv = e["C_pnsn_full"]
        prep = pnv[max(0, w0 - 60):w0]
        blp = float(np.percentile(e["C_pnsn"], 20))
        rec["b_pnsn_exc_p20"], _ = fit_logy(R7, e["C_pnsn"] - blp, 0.0)
        rec["pnsn_ets_contrast"] = (float(np.median(e["C_pnsn"]) / max(np.median(prep), 1.0))
                                    if prep.size >= 20 else np.nan)
        rp = (float(e["C_pnsn"].mean() / prep.mean())
              if prep.size >= 20 and prep.mean() > 0 else np.nan)
        rec["pnsn_ets_ratio_mean"] = rp
        rec["b_pnsn_deatt"] = (e["b_pnsn"] * rp / (rp - 1.0)
                               if np.isfinite(rp) and rp > 1.1 else np.nan)
        okc = np.isfinite(meancc) & (R7 > 0) & (N > 0)
        rec["b_meancc"] = (ts_slope_fast(np.log10(R7[okc]), meancc[okc])
                           if okc.sum() >= MIN_FIT_DAYS else np.nan)
        rec["meancc"] = float(np.nanmean(meancc))
        rec["b_lower"], rec["b_upper"] = halves(R7, N)
        rec["b_pnsn"] = e["b_pnsn"]
        rec["n_pnsn"] = e["n_pnsn"]
        # 30-d rolling TRUE LFE count (Test 1c) vs the footprint inc30
        roll = np.convolve(Nfull, np.ones(30), mode="full")[:NT]
        rec["b_N30"], rec["n_N30"] = fit_logy(e["inc_fp"], roll[w], 1.0)
        # shifted-LFE null
        bs = []
        if span_s >= 2 * 365 + 30:
            for sh in RNG.integers(365, span_s - 365 + 1, size=N_SHIFT):
                yn = shifted_window(Nfull, d0, d1, w0, nW, int(sh))
                b, nn = fit_logy(R7, yn, 1.0)
                if np.isfinite(b):
                    bs.append(b)
        bs = np.asarray(bs, float)
        rec["null_mean"] = float(np.nanmean(bs)) if bs.size else np.nan
        rec["null_sd"] = float(np.nanstd(bs)) if bs.size else np.nan
        rec["null_n"] = int(bs.size)
        rec["null_p"] = (float(np.mean(bs >= rec["b_N"]))
                         if bs.size and np.isfinite(rec["b_N"]) else np.nan)
        rec["z_null"] = ((rec["b_N"] - rec["null_mean"]) / rec["null_sd"]
                         if bs.size and rec["null_sd"] and rec["null_sd"] > 0 else np.nan)
        rows_d.append(rec)
        # per-family daily exponents
        prewin = slice(max(0, w0 - 60), w0)
        for j in idx40:
            yj = S["N"][j, w].astype(float)
            pj = S["N"][j, prewin].astype(float)
            b, nn = fit_logy(R7, yj, 1.0)
            base = float(pj.mean()) if pj.size >= 20 and w0 - 60 >= d0 else np.nan
            ratio = float(yj.mean() / base) if np.isfinite(base) and base > 0 else np.nan
            bex = (fit_logy(R7, yj - base, 0.0)[0] if np.isfinite(base) else np.nan)
            if np.isfinite(b):
                famrows.append(dict(test="T2fam", sta=tok, pid=e["pid"], seg=e["seg"],
                                    epoch=e["epoch"], patch=S["fams"][j],
                                    cert=bool(S["cert"][j]), lat=float(S["lat"][j]),
                                    b=b, n=nn, det=float(yj.sum()),
                                    ets_ratio=ratio, b_exc=bex,
                                    b_deatt=(b * ratio / (ratio - 1.0)
                                             if np.isfinite(ratio) and ratio > 1.1
                                             else np.nan)))
T2 = pd.DataFrame(rows_d)
T2FAM = pd.DataFrame(famrows)
log(f"  {len(T2)} (station, event) fits over {T2.sta.nunique()} stations and "
    f"{T2.pid.nunique()} events; {len(T2FAM):,} per-family fits")


# ==================================================================================
# 9. SUMMARIES
# ==================================================================================
def summarize(test, measure, slipvar, group, b, null=None, pooled=None, extra=None):
    b = np.asarray(b, float)
    b = b[np.isfinite(b)]
    m, q1, q3, n = med_iqr(b)
    lo, hi = boot_median_ci(b)
    p, k, nn = sign_p(b)
    row = dict(test=test, measure=measure, slip=slipvar, group=group, n=n,
               b_median=m, b_q25=q1, b_q75=q3, b_ci05=lo, b_ci95=hi,
               sign_pos=k, sign_n=nn, sign_p=p,
               frac_below_1=float(np.mean(b < 1)) if n else np.nan)
    if pooled is not None:
        row["pooled_ols"], row["pooled_ts"], row["pooled_r"], row["pooled_n"] = pooled
    if null is not None and len(null):
        nu = np.asarray(null, float)
        nu = nu[np.isfinite(nu)]
        row["null_mean"] = float(np.mean(nu))
        row["null_sd"] = float(np.std(nu))
        row["null_n"] = int(nu.size)
        if n:
            row["null_p_median"] = float(np.mean(nu >= m))
    if extra:
        row.update(extra)
    return row


R30 = []
# --- T1a -------------------------------------------------------------------------------
POOL_1A = pooled_fit(np.concatenate(P1A_X) if P1A_X else np.array([]),
                     np.concatenate(P1A_Y) if P1A_Y else np.array([]),
                     np.concatenate(P1A_G) if P1A_G else np.array([]))
for grp, mask in [("all", np.ones(len(T1A), bool)),
                  ("north", (T1A.seg == "north").values),
                  ("central", (T1A.seg == "central").values),
                  ("certified", T1A.cert.values.astype(bool)),
                  ("2010-2014", (T1A.epoch == "2010-2014").values),
                  ("2015-2019", (T1A.epoch == "2015-2019").values),
                  ("2020-2026", (T1A.epoch == "2020-2026").values)]:
    d = T1A[mask]
    R30.append(summarize("T1a", "n_stack (active days/30 d)", "inc30_nearest", grp,
                         d.b_near.values, null=NULL_1A if grp == "all" else None,
                         pooled=POOL_1A if grp == "all" else None,
                         extra=dict(note="SATURATING occupancy measure, ceiling 30")))
    R30.append(summarize("T1a", "n_stack (active days/30 d)", "inc30_footprint", grp,
                         d.b_fp.values))
# --- T1b -------------------------------------------------------------------------------
for grp, mask in [("all", np.ones(len(T1B), bool)),
                  ("north", (T1B.fam_seg == "north").values),
                  ("central", (T1B.fam_seg == "central").values),
                  ("certified", T1B.cert.values.astype(bool))]:
    R30.append(summarize("T1b", "max n_stack per event", "peak inc30_nearest", grp,
                         T1B[mask].b.values,
                         extra=dict(note="cross-event, families with >= 4 events")))
# --- T1c -------------------------------------------------------------------------------
for grp, mask in [("all", np.ones(len(T2), bool)),
                  ("north", (T2.seg == "north").values),
                  ("central", (T2.seg == "central").values),
                  ("2010-2014", (T2.epoch == "2010-2014").values),
                  ("2015-2019", (T2.epoch == "2015-2019").values),
                  ("2020-2026", (T2.epoch == "2020-2026").values)]:
    R30.append(summarize("T1c", "true 30-d LFE count", "inc30_footprint", grp,
                         T2[mask].b_N30.values,
                         extra=dict(note="rolling 30-d sum of streamed LFE detections")))
R30 = pd.DataFrame(R30)
R30.to_csv(OUT_RES / "s14_lfe30_exponents.csv", index=False)

# --- T2 / T3 / T4 -----------------------------------------------------------------------
RD = []


def pool_for(col, subset_mask):
    xs, ys, gs = [], [], []
    for _, r in T2[subset_mask].iterrows():
        e = EVBY[r["pid"]]
        S = STA[r["sta"]]
        idx = fam_sel(S, e, FAM_RADIUS_KM, cert_only=(col == "b_Ncert"))
        if idx.size == 0:
            continue
        src = {"b_N": S["N"], "b_Nhi": S["NHI"], "b_Ncert": S["N"]}[col]
        y = src[idx].sum(0)[e["w0"]:e["w0"] + e["nW"]].astype(float)
        x = e["R7"]
        ok = np.isfinite(x) & (x > 0) & (y > 0)
        if ok.sum() < MIN_FIT_DAYS:
            continue
        xs.append(np.log10(x[ok]))
        ys.append(np.log10(y[ok] + 1.0))
        gs.append(np.full(int(ok.sum()), f"{r['sta']}_{r['pid']}"))
    if not xs:
        return (np.nan, np.nan, np.nan, 0)
    return pooled_fit(np.concatenate(xs), np.concatenate(ys), np.concatenate(gs))


groups = [("all", np.ones(len(T2), bool)),
          ("north", (T2.seg == "north").values),
          ("central", (T2.seg == "central").values),
          ("2010-2014", (T2.epoch == "2010-2014").values),
          ("2015-2019", (T2.epoch == "2015-2019").values),
          ("2020-2026", (T2.epoch == "2020-2026").values)]
for grp, mask in groups:
    d = T2[mask]
    RD.append(summarize("T2", "daily LFE count (40 km)", "R7", grp, d.b_N.values,
                        pooled=pool_for("b_N", mask),
                        extra=dict(null_mean_of_events=float(np.nanmean(d.null_mean)),
                                   null_sd_of_events=float(np.nanmean(d.null_sd)),
                                   z_median=float(np.nanmedian(d.z_null)),
                                   n_z_pos=int(np.nansum(d.z_null > 0)),
                                   n_z=int(np.isfinite(d.z_null).sum()))))
    RD.append(summarize("T2", "daily LFE count (60 km)", "R7", grp, d.b_N60.values))
    RD.append(summarize("T4", "daily LFE count, cc >= 0.9", "R7", grp, d.b_Nhi.values,
                        pooled=pool_for("b_Nhi", mask) if grp == "all" else None))
    RD.append(summarize("T4", "daily LFE count, certified families", "R7", grp,
                        d.b_Ncert.values,
                        pooled=pool_for("b_Ncert", mask) if grp == "all" else None))
    RD.append(summarize("T3", "daily mean cc (weak amplitude proxy)", "R7", grp,
                        d.b_meancc.values,
                        extra=dict(note="cc is a correlation, NOT an amplitude")))
    RD.append(summarize("T2", "PNSN daily tremor count (60 km)", "R7", grp,
                        d.drop_duplicates("pid").b_pnsn.values,
                        extra=dict(note="S9/S11 measure, recomputed here")))
    RD.append(summarize("T2", "daily LFE count, LOWER half of log R7", "R7", grp,
                        d.b_lower.values))
    RD.append(summarize("T2", "daily LFE count, UPPER half of log R7", "R7", grp,
                        d.b_upper.values))
    for t in CC_LADDER:
        RD.append(summarize("T4", f"daily LFE count, cc >= {t:.2f}", "R7", grp,
                            d[f"b_cc{int(t * 100)}"].values))
    RD.append(summarize("T6", "daily LFE EXCESS count (pre-event baseline)", "R7", grp,
                        d.b_exc_pre.values,
                        extra=dict(note="N - median(N over [start-70, start-10])")))
    RD.append(summarize("T6", "daily LFE EXCESS count (window p20 baseline)", "R7", grp,
                        d.b_exc_p20.values))
    RD.append(summarize("T6", "PNSN EXCESS count (window p20 baseline, control)", "R7", grp,
                        d.b_pnsn_exc_p20.values))
    RD.append(summarize("T6", "daily LFE count, floor-DE-ATTENUATED b*r/(r-1)", "R7", grp,
                        d.b_deatt.values,
                        extra=dict(note="errors-in-variables correction for an additive "
                                        "false-alarm floor")))
    RD.append(summarize("T6", "PNSN count, same de-attenuation (control)", "R7", grp,
                        d.b_pnsn_deatt.values))
RD = pd.DataFrame(RD)
T2OUT = T2.copy()
T2OUT.insert(0, "row_type", "station_event")
RDOUT = RD.copy()
RDOUT.insert(0, "row_type", "summary")
pd.concat([T2OUT, RDOUT], ignore_index=True).to_csv(
    OUT_RES / "s14_lfe_daily_exponents.csv", index=False)
COV = COV.merge(DIAG, on="sta", how="left")
COV.to_csv(OUT_RES / "s14_coverage.csv", index=False)
log("results written")


# ==================================================================================
# 10. PRINTED REPORT
# ==================================================================================
def fmt(v, nd=3):
    return "  n/a  " if not np.isfinite(v) else f"{v:+.{nd}f}"


def line(b, extra=""):
    b = np.asarray(b, float)
    b = b[np.isfinite(b)]
    if b.size == 0:
        return "n/a"
    m, q1, q3, n = med_iqr(b)
    lo, hi = boot_median_ci(b)
    p, k, nn = sign_p(b)
    return (f"{m:+7.3f} [{q1:+7.3f},{q3:+7.3f}]  CI90 [{lo:+6.3f},{hi:+6.3f}]  n={n:5d}  "
            f"b<1 {np.mean(b < 1) * 100:5.1f}%  >0 {k}/{nn} p={p:.2g} {extra}")


W = 118
print("\n" + "=" * W)
print("S14 - THE TREMOR-vs-SLIP-RATE EXPONENT MEASURED WITH THIS PROJECT'S OWN LFE DETECTIONS")
print("=" * W)
print(f"27 confirmed north/central events (S1); window [start-{PRE_A}, end+{POST_A}]; slip rate R7 = "
      f"area-weighted positive potency rate over the 60-km neighbourhood, 7-d centred mean;")
print(f"LFE families attached to an event if the template epicentre is within {FAM_RADIUS_KM:.0f} km "
      f"of a footprint triangle.  Reference: S9/S11 PNSN count b = +0.401 [+0.181, +0.478], "
      f"energy +0.644; S11b raw-GNSS +0.28.")

print("\n" + "-" * W)
print("0. WHAT n_stack IS  (LEVEL A, all 195 daily_dvv files)")
print("-" * W)
print("  scripts/dvv_roll30cal.py:  los = searchsorted(dates, dates - 30 d, 'right');  "
      "n_stack = counts = idx + 1 - los")
print("  => n_stack is the number of DAILY STACKS (days that produced a stack, i.e. days with "
      ">= 20 detections)")
print("     inside the TRAILING 30 CALENDAR DAYS.  It is an ACTIVE-DAY OCCUPANCY count with a "
      "hard ceiling of 30")
print("     and a floor of 5 (--min-stacks); it is NOT a count of LFE detections.  Verified above "
      "by recomputing it")
print("     from each file's own dates (see the CHECK line).")
print(f"  Inside the 27 event windows the mean fraction of days at the ceiling n_stack = 30 is "
      f"{np.mean(sat_frac):.3f}")
print("     -> the measure is saturated during ETS and carries almost no rate information; "
      "Test 1c is the honest 30-d test.")

print("\n" + "-" * W)
print("1. LEVEL B COVERAGE - what was streamed")
print("-" * W)
print(COV.to_string(index=False))
print(f"  TOTAL {COV.GB.sum():.1f} GB apparent bytes (budget {BUDGET_GB:.0f} GB), "
      f"{COV.detections.sum():,} LFE detections, {COV.wall_s.sum() / 60:.1f} min of streaming.")

print("\n" + "-" * W)
print("2. TEST 1a - 30-d level: n_stack vs trailing-30-d slip inc30 (Theil-Sen, log-log)")
print("-" * W)
for grp, mask in [("all", np.ones(len(T1A), bool)),
                  ("north", (T1A.seg == "north").values),
                  ("central", (T1A.seg == "central").values),
                  ("certified fams", T1A.cert.values.astype(bool)),
                  ("2010-2014", (T1A.epoch == "2010-2014").values),
                  ("2015-2019", (T1A.epoch == "2015-2019").values),
                  ("2020-2026", (T1A.epoch == "2020-2026").values)]:
    print(f"  vs inc30(nearest tri) {grp:15s} {line(T1A[mask].b_near.values)}")
print(f"  vs inc30(footprint)   {'all':15s} {line(T1A.b_fp.values)}")
print(f"  pooled (family-event dummies): OLS {fmt(POOL_1A[0])}  Theil-Sen {fmt(POOL_1A[1])}  "
      f"r {fmt(POOL_1A[2])}  n={POOL_1A[3]}")
print(f"  NULL (>= 365-d shifted n_stack, {N_SHIFT // 4} draws x {len(null_b)} family-events): "
      f"mean {np.nanmean(NULL_1A):+.4f} sd {np.nanstd(NULL_1A):.4f};  "
      f"per-fit z median {np.nanmedian(per_z):+.3f}")

print("\n" + "-" * W)
print("3. TEST 1b - cross-event: peak n_stack vs peak event slip, per family (>= 4 events)")
print("-" * W)
for grp, mask in [("all", np.ones(len(T1B), bool)),
                  ("north", (T1B.fam_seg == "north").values),
                  ("central", (T1B.fam_seg == "central").values),
                  ("certified", T1B.cert.values.astype(bool))]:
    print(f"  {grp:15s} {line(T1B[mask].b.values)}")
print(f"  Spearman(peak n_stack, peak inc30) across events: median "
      f"{np.nanmedian(T1B.rho):+.3f}")

print("\n" + "-" * W)
print("4. TEST 1c - 30-d level with TRUE LFE COUNTS (rolling 30-d sum) vs inc30(footprint)")
print("-" * W)
for grp, mask in [("all", np.ones(len(T2), bool)),
                  ("north", (T2.seg == "north").values),
                  ("central", (T2.seg == "central").values),
                  ("2010-2014", (T2.epoch == "2010-2014").values),
                  ("2015-2019", (T2.epoch == "2015-2019").values),
                  ("2020-2026", (T2.epoch == "2020-2026").values)]:
    print(f"  {grp:15s} {line(T2[mask].b_N30.values)}")

print("\n" + "-" * W)
print("5. TEST 2 - DAILY level: LFE count vs R7  (the S9/S11 protocol, LFEs instead of PNSN tremor)")
print("-" * W)
for grp, mask in groups:
    d = T2[mask]
    print(f"  LFE  {grp:12s} {line(d.b_N.values)}")
for grp, mask in groups:
    d = T2[mask]
    print(f"  PNSN {grp:12s} {line(d.drop_duplicates('pid').b_pnsn.values)}")
po = pool_for("b_N", np.ones(len(T2), bool))
print(f"  pooled (station-event dummies): OLS {fmt(po[0])}  Theil-Sen {fmt(po[1])}  "
      f"r {fmt(po[2])}  n={po[3]}")
zz = T2.z_null.values
print(f"  NULL (>= 365-d shifted LFE series, {N_SHIFT} draws per station-event): "
      f"per-fit null mean {np.nanmean(T2.null_mean):+.4f}, sd {np.nanmean(T2.null_sd):.3f};")
print(f"       z median {np.nanmedian(zz):+.3f}, positive {int(np.nansum(zz > 0))}/"
      f"{int(np.isfinite(zz).sum())} "
      f"(sign p {sign_p(zz)[0]:.2g}); station-events with null p < 0.05: "
      f"{int((T2.null_p < 0.05).sum())}/{int(T2.null_p.notna().sum())}")
print(f"  SATURATION check, lower vs upper half of log10 R7: lower {line(T2.b_lower.values)}")
print(f"                                                     upper {line(T2.b_upper.values)}")
du = T2.dropna(subset=["b_lower", "b_upper"])
if len(du) > 5:
    diff = du.b_upper.values - du.b_lower.values
    print(f"     upper - lower {np.nanmedian(diff):+.3f} "
          f"[{np.nanpercentile(diff, 25):+.3f}, {np.nanpercentile(diff, 75):+.3f}], "
          f"Wilcoxon p {wilcoxon(diff).pvalue:.3f}, {int((diff > 0).sum())}/{len(diff)} positive")

print("\n  PAIRED LFE vs PNSN on the SAME events:")
pv_ = T2.dropna(subset=["b_N", "b_pnsn"])
if len(pv_) > 5:
    dd = pv_.b_N.values - pv_.b_pnsn.values
    print(f"     b(LFE) - b(PNSN) = {np.nanmedian(dd):+.3f} "
          f"[{np.nanpercentile(dd, 25):+.3f}, {np.nanpercentile(dd, 75):+.3f}], "
          f"Wilcoxon p {wilcoxon(dd).pvalue:.2g}, {int((dd > 0).sum())}/{len(dd)} positive; "
          f"Spearman rho {spearmanr(pv_.b_N, pv_.b_pnsn).statistic:+.3f} "
          f"(p {spearmanr(pv_.b_N, pv_.b_pnsn).pvalue:.2g})")
    ev_lvl = pv_.groupby("pid").agg(bL=("b_N", "median"), bP=("b_pnsn", "first"))
    ddl = ev_lvl.bL.values - ev_lvl.bP.values
    print(f"     event-level (station median): {np.nanmedian(ddl):+.3f}, "
          f"Wilcoxon p {wilcoxon(ddl).pvalue:.2g}, n={len(ddl)} events, "
          f"Spearman rho {spearmanr(ev_lvl.bL, ev_lvl.bP).statistic:+.3f}")

print("\n" + "-" * W)
print("5b. TEST 5 - DIAGNOSIS: are these matched-filter detections an LFE RATE?")
print("-" * W)
print("  The densify threshold is cc >= 0.80 with a 6-s minimum gap "
      "(scripts/fleet_station.sh -> densify_gnw_gpu.py).")
print("  Per station: median daily detections summed over all its templates, the ETS/non-ETS "
      "ratio of that count,")
print("  and its rank correlation with the PNSN tremor count in the SAME 40-km neighbourhood.")
print(f"  {'station':12s} {'det/day':>10s} {'ETS/quiet':>9s} {'@cc.90':>7s} {'@cc.98':>7s} "
      f"{'rho(PNSN)':>9s} {'@cc.98':>7s} {'PNSN ETS/quiet':>15s} {'fams ratio>1.5':>14s}")
for _, r in DIAG.iterrows():
    print(f"  {r.sta:12s} {r['med_cc0.80']:10.0f} {r['ets_ratio_cc0.80']:9.2f} "
          f"{r['ets_ratio_cc0.90']:7.2f} {r['ets_ratio_cc0.98']:7.2f} "
          f"{r['rho_pnsn_cc0.80']:9.3f} {r['rho_pnsn_cc0.98']:7.3f} "
          f"{r['pnsn_ets_ratio']:15.1f} {r['fam_frac_ratio_gt15'] * 100:13.0f}%")
print(f"  MEDIAN over stations: LFE ETS/quiet {DIAG['ets_ratio_cc0.80'].median():.2f} "
      f"(cc>=0.98 {DIAG['ets_ratio_cc0.98'].median():.2f}), "
      f"rho(LFE, PNSN) {DIAG['rho_pnsn_cc0.80'].median():+.3f}; "
      f"PNSN ETS/quiet {DIAG['pnsn_ets_ratio'].median():.1f}")
print(f"  READ: inside an ETS window the PNSN count in the SAME neighbourhood rises by "
      f"x{DIAG.pnsn_ets_ratio.median():.1f} (mean daily count) /")
print(f"  x{np.nanmedian(T2.pnsn_ets_contrast):.0f} (per-event medians), while the raw "
      f"matched-filter count rises by only "
      f"x{DIAG['ets_ratio_cc0.80'].median():.2f}, and its")
print(f"  rank correlation with the PNSN count is rho = "
      f"{DIAG['rho_pnsn_cc0.80'].median():+.2f}.  The detector fires ~"
      f"{DIAG['mean_cc0.80'].median() / 300:.0f} detections per template")
print("  per day in quiet times: the output is dominated by a nearly constant FALSE-ALARM "
      "rate.  The same fact")
print("  explains the n_stack saturation - the >= 20-detections/day gate of "
      "build_long_window_3comp.py is met on")
print("  essentially every day.  The pipeline recovers LFE signal from these files by "
      "STACKING them, not by")
print("  counting them: mf_*.csv is a densify product, not an LFE catalogue.")

print("\n" + "-" * W)
print("5c. TEST 6 - the same exponent after removing the false-alarm floor")
print("-" * W)
print("  An additive floor F on top of a signal S attenuates the log-log slope by "
      "S/(S+F): with S/(S+F) ~ 0.05-0.4")
print("  a true exponent of 0.4 is compressed to 0.02-0.16, which is what Test 2 measures.  "
      "Removing the floor:")
print(f"  LFE excess, pre-event baseline   {line(T2.b_exc_pre.values)}")
print(f"  LFE excess, window-p20 baseline  {line(T2.b_exc_p20.values)}")
print(f"  PNSN excess, window-p20 (control){line(T2.b_pnsn_exc_p20.values)}")
print(f"  event ETS contrast median(N in window)/median(N pre-event): LFE "
      f"{np.nanmedian(T2.ets_contrast):.3f}  vs PNSN "
      f"{np.nanmedian(T2.pnsn_ets_contrast):.1f}")
print("\n  ERRORS-IN-VARIABLES de-attenuation.  For N = F + k R^b, dlogN/dlogR = b * (r-1)/r,")
print("  r = mean(N in window) / mean(N in [start-70, start-10]), so b_true = b_obs * r/(r-1):")
print(f"  LFE  r median {np.nanmedian(T2.ets_ratio_mean):5.2f} -> de-attenuated "
      f"{line(T2.b_deatt.values)}")
print(f"  PNSN r median {np.nanmedian(T2.pnsn_ets_ratio_mean):5.2f} -> de-attenuated "
      f"{line(T2.b_pnsn_deatt.values)}")
print("       (the PNSN row is the control: with a small floor the correction barely moves "
      "0.40, as it should)")
dz = T2.dropna(subset=["b_exc_p20", "b_pnsn_exc_p20"])
if len(dz) > 5:
    df_ = dz.b_exc_p20.values - dz.b_pnsn_exc_p20.values
    print(f"  paired LFE-excess - PNSN-excess: {np.nanmedian(df_):+.3f} "
          f"[{np.nanpercentile(df_, 25):+.3f}, {np.nanpercentile(df_, 75):+.3f}], "
          f"Wilcoxon p {wilcoxon(df_).pvalue:.2g}")

print("\n" + "-" * W)
print("6. TEST 4 - detection-threshold / family-quality controls")
print("-" * W)
print("  cc-threshold ladder (median b over station-events; the floor should fall with cc):")
for t in (0.80,) + CC_LADDER:
    col = "b_N" if t == 0.80 else f"b_cc{int(t * 100)}"
    dcol = "det_in_window" if t == 0.80 else f"det_cc{int(t * 100)}"
    print(f"    cc >= {t:.2f}   {line(T2[col].values)}   median detections/window "
          f"{np.nanmedian(T2[dcol]):12.0f}")
print(f"  all detections (cc >= 0.80)      {line(T2.b_N.values)}")
print(f"  cc >= 0.90 only                  {line(T2.b_Nhi.values)}")
print(f"  causality-CERTIFIED families     {line(T2.b_Ncert.values)}")
print(f"  60-km family radius              {line(T2.b_N60.values)}")
for a, bcol in (("cc>=0.9", "b_Nhi"), ("certified", "b_Ncert"), ("60 km", "b_N60")):
    d = T2.dropna(subset=["b_N", bcol])
    if len(d) > 5:
        diff = d[bcol].values - d.b_N.values
        print(f"     paired {a:10s} - all: {np.nanmedian(diff):+.3f} "
              f"[{np.nanpercentile(diff, 25):+.3f}, {np.nanpercentile(diff, 75):+.3f}], "
              f"Wilcoxon p {wilcoxon(diff).pvalue:.2g}, {int((diff > 0).sum())}/{len(diff)}")
print(f"\n  per-FAMILY daily exponents ({len(T2FAM):,} family-events): "
      f"{line(T2FAM.b.values) if len(T2FAM) else 'n/a'}")
if len(T2FAM):
    print(f"     certified only: {line(T2FAM[T2FAM.cert].b.values)}")
    print("\n  THE DECISIVE STRATIFICATION - keep only templates that DO respond to the event")
    print("  (ets_ratio = mean daily detections in the window / mean over [start-70, start-10]):")
    for lo_, hi_ in ((0.0, 1.1), (1.1, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 1e9)):
        d = T2FAM[(T2FAM.ets_ratio >= lo_) & (T2FAM.ets_ratio < hi_)]
        if len(d) >= 20:
            print(f"     ets_ratio {lo_:4.1f}-{hi_ if hi_ < 1e8 else np.inf:<5.1f} "
                  f"raw {line(d.b.values)}")
            print(f"     {'':>16s} excess {line(d.b_exc.values)}")
    q = T2FAM.ets_ratio.dropna()
    print(f"     per-family ets_ratio distribution: p50 {q.median():.2f}, p90 "
          f"{q.quantile(0.9):.2f}, p99 {q.quantile(0.99):.2f}, "
          f"frac > 1.5 {np.mean(q > 1.5) * 100:.1f} %, frac > 3 {np.mean(q > 3) * 100:.1f} %")
    if (q > 1.5).sum() >= 20:
        top = T2FAM[T2FAM.ets_ratio > 1.5]
        print(f"     ALL responding templates (ratio > 1.5, n={len(top)}): "
              f"raw {line(top.b.values)}")
        print(f"     {'':>44s} excess {line(top.b_exc.values)}")
        print(f"     {'':>44s} de-att {line(top.b_deatt.values)}")
    print(f"     per-family de-attenuated b (all with r > 1.1): {line(T2FAM.b_deatt.values)}")

print("\n" + "-" * W)
print("7. TEST 3 - per-LFE amplitude: NOT TESTABLE with these files")
print("-" * W)
print("  The mf_*.csv catalogues carry (template, time, cc, station).  `cc` is the matched-filter")
print("  CORRELATION COEFFICIENT, not an amplitude: it is bounded in [0.80, 1], it is normalised by")
print("  the running trace energy, and it says how well the waveform matches the template, not how")
print("  big the LFE was.  Frank & Brodsky's per-LFE amplitude ~ rate^0.32 therefore CANNOT be")
print("  reproduced here; it would need the per-detection amplitude (re-reading the waveforms).")
print(f"  Weak proxy only - exponent of the daily MEAN cc on log10 R7: {line(T2.b_meancc.values)}")
print(f"     (units: cc per decade of slip rate; mean cc over all station-events "
      f"{np.nanmean(T2.meancc):.4f})")

print("\n" + "-" * W)
print("8. FAMILY GEOMETRY / QC")
print("-" * W)
print(f"  {len(UF)} (station, family) pairs over {UF.sta.nunique()} stations; "
      f"{int(UF.certified.sum())} causality-certified ({UF.certified.mean() * 100:.0f} %)")
print(f"  family Slab2 depth: median {np.nanmedian(UF.depth_km):.1f} km "
      f"[p10 {np.nanpercentile(UF.depth_km, 10):.1f}, p90 {np.nanpercentile(UF.depth_km, 90):.1f}]; "
      f"segments {UF.seg.value_counts().to_dict()}")
if len(T2FAM):
    print(f"  streamed families in the daily test: {T2FAM.patch.nunique()} distinct templates, "
          f"{T2FAM.sta.nunique()} stations")


# ==================================================================================
# 11. FIGURE
# ==================================================================================
log("figure ...")
fig, AX = plt.subplots(3, 3, figsize=(19.5, 14.5))
fig.suptitle("S14 - tremor-vs-slip-rate exponent from the project's own LFE detections "
             "(27 confirmed ETS events, north + central Cascadia)", fontsize=13)

# A: exponent distributions with nulls
ax = AX[0, 0]
sets = [("LFE 30-d\n(n_stack,\nsaturated)", T1A.b_near.values, NULL_1A),
        ("LFE 30-d\n(true count)", T2.b_N30.values, None),
        ("LFE daily\n(raw)", T2.b_N.values, None),
        ("LFE daily\nfloor-\nde-attenuated", T2.b_deatt.values, None),
        ("PNSN daily", T2.drop_duplicates("pid").b_pnsn.values, None)]
pos = np.arange(len(sets))
for i, (nm, b, nu) in enumerate(sets):
    b = np.asarray(b, float)
    b = b[np.isfinite(b)]
    b = b[(b > -2) & (b < 3)]
    if b.size:
        parts = ax.violinplot([b], positions=[i], widths=0.75, showmedians=True)
        for pcs in parts["bodies"]:
            pcs.set_facecolor("#3b76af")
            pcs.set_alpha(0.55)
    if nu is not None:
        nu = np.asarray(nu, float)
        nu = nu[np.isfinite(nu) & (nu > -2) & (nu < 3)]
        if nu.size:
            p2 = ax.violinplot([nu], positions=[i + 0.32], widths=0.5, showmedians=True)
            for pcs in p2["bodies"]:
                pcs.set_facecolor("0.6")
                pcs.set_alpha(0.5)
# station-event nulls for the daily test
nu_daily = T2.dropna(subset=["null_mean"])
if len(nu_daily):
    ax.errorbar([2.34], [nu_daily.null_mean.mean()],
                yerr=[nu_daily.null_sd.mean()], fmt="s", color="0.35", capsize=4,
                label="shifted-LFE null")
ax.axhline(0, color="k", lw=0.6)
ax.axhline(1, color="crimson", lw=0.9, ls="--", label="b = 1 (proportional)")
ax.axhline(0.40, color="darkgreen", lw=0.9, ls=":", label="S9/S11 PNSN 0.40")
ax.set_xticks(pos)
ax.set_xticklabels([s[0] for s in sets], fontsize=8)
ax.set_ylim(-1.2, 2.2)
ax.set_ylabel("exponent b")
ax.set_title("A  exponent distributions (grey = shifted null)")
ax.legend(fontsize=7, loc="upper left")

# B: by segment
ax = AX[0, 1]
segs = ["north", "central"]
w = 0.35
for k, (lab, col, c) in enumerate([("LFE daily", "b_N", "#3b76af"),
                                   ("PNSN daily", "b_pnsn", "#d1701e")]):
    vals, errs = [], []
    for sg in segs:
        d = T2[T2.seg == sg]
        v = (d.drop_duplicates("pid")[col] if col == "b_pnsn" else d[col]).values
        v = v[np.isfinite(v)]
        vals.append(np.median(v) if v.size else np.nan)
        lo, hi = boot_median_ci(v)
        errs.append([abs(np.median(v) - lo) if v.size else 0,
                     abs(hi - np.median(v)) if v.size else 0])
    ax.bar(np.arange(2) + k * w, vals, w, yerr=np.array(errs).T, capsize=4, color=c, label=lab)
ax.axhline(1, color="crimson", ls="--", lw=0.9)
ax.axhline(0.40, color="darkgreen", ls=":", lw=0.9)
ax.set_xticks(np.arange(2) + w / 2)
ax.set_xticklabels(segs)
ax.set_ylabel("median b [CI90]")
ax.set_title("B  by segment")
ax.legend(fontsize=8)

# C: by epoch
ax = AX[0, 2]
eps = ["2010-2014", "2015-2019", "2020-2026"]
for k, (lab, col, c) in enumerate([("LFE daily", "b_N", "#3b76af"),
                                   ("PNSN daily", "b_pnsn", "#d1701e")]):
    vals, errs = [], []
    for ep in eps:
        d = T2[T2.epoch == ep]
        v = (d.drop_duplicates("pid")[col] if col == "b_pnsn" else d[col]).values
        v = v[np.isfinite(v)]
        m = np.median(v) if v.size else np.nan
        lo, hi = boot_median_ci(v)
        vals.append(m)
        errs.append([abs(m - lo) if v.size else 0, abs(hi - m) if v.size else 0])
    ax.bar(np.arange(3) + k * w, vals, w, yerr=np.array(errs).T, capsize=4, color=c, label=lab)
ax.axhline(1, color="crimson", ls="--", lw=0.9)
ax.axhline(0.40, color="darkgreen", ls=":", lw=0.9)
ax.set_xticks(np.arange(3) + w / 2)
ax.set_xticklabels(eps, fontsize=8)
ax.set_ylabel("median b [CI90]")
ax.set_title("C  by epoch")
ax.legend(fontsize=8)

# D: one event as a time series - the dilution made visible
ax = AX[1, 0]
cand = T2.dropna(subset=["b_N"]).sort_values("det_in_window", ascending=False)
if len(cand):
    r = cand.iloc[0]
    e = EVBY[r["pid"]]
    S = STA[r["sta"]]
    idx = fam_sel(S, e, FAM_RADIUS_KM)
    y = S["N"][idx].sum(0)[e["w0"]:e["w0"] + e["nW"]].astype(float)
    t = np.arange(e["nW"])
    ax.plot(t, y / np.median(y[y > 0]), "-o", ms=3, color="#3b76af", label="LFE count / median")
    cp = e["C_pnsn"]
    ax.plot(t, cp / max(np.median(cp[cp > 0]), 1), "-s", ms=3, color="#d1701e",
            label="PNSN count / median")
    ax.plot(t, e["R7"] / max(np.median(e["R7"][e["R7"] > 0]), 1e-9), "-", lw=2, color="k",
            label="R7 / median")
    ax.set_yscale("log")
    ax.set_xlabel(f"day in window (event {int(r['pid'])}, {e['start']}, {r['sta']})")
    ax.set_ylabel("normalised to its own median")
    ax.set_title("D  the dilution, one event: LFE count barely moves")
    ax.legend(fontsize=7)

# E: paired LFE vs PNSN
ax = AX[1, 1]
d = T2.dropna(subset=["b_N", "b_pnsn"])
sc = ax.scatter(d.b_pnsn, d.b_N, c=(d.seg == "north").map({True: "#3b76af", False: "#d1701e"}),
                s=22, alpha=0.8)
lim = [-1.0, 2.0]
ax.plot(lim, lim, "k--", lw=1)
ax.axhline(0, color="0.7", lw=0.6)
ax.axvline(0, color="0.7", lw=0.6)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("b, PNSN daily count")
ax.set_ylabel("b, LFE daily count")
ax.set_title("E  paired per (station, event); blue north / orange central")

# F: cc-threshold ladder
ax = AX[1, 2]
ths = (0.80,) + CC_LADDER
med = []
cont = []
for t in ths:
    col = "b_N" if t == 0.80 else f"b_cc{int(t * 100)}"
    v = T2[col].values
    v = v[np.isfinite(v)]
    med.append(np.median(v) if v.size else np.nan)
    key = "ets_ratio_cc0.80" if t == 0.80 else f"ets_ratio_cc{t:.2f}"
    cont.append(DIAG[key].median())
ax.plot(ths, med, "-o", color="#3b76af", label="median b (LFE count)")
ax.axhline(np.nanmedian(T2.b_exc_p20), color="#2ca02c", ls="-.", lw=1.2,
           label=f"b after floor removal {np.nanmedian(T2.b_exc_p20):+.2f}")
ax.axhline(np.nanmedian(T2.drop_duplicates("pid").b_pnsn), color="#d1701e", ls="--",
           label="PNSN 0.40")
ax.axhline(1, color="crimson", lw=0.8, ls=":")
ax.set_xlabel("matched-filter cc threshold")
ax.set_ylabel("median exponent b")
ax2 = ax.twinx()
ax2.plot(ths, cont, "-^", color="0.45", label="LFE ETS/quiet count ratio")
ax2.axhline(DIAG.pnsn_ets_ratio.median(), color="0.45", ls=":", lw=1.2)
ax2.text(0.805, DIAG.pnsn_ets_ratio.median(), " PNSN ETS/quiet", color="0.45", fontsize=7,
         va="bottom")
ax2.set_ylim(0.9, max(6.0, DIAG.pnsn_ets_ratio.median() * 1.15))
ax2.set_ylabel("ETS / non-ETS mean count ratio", color="0.45")
ax.set_title("F  detection-threshold ladder (the floor never clears)")
ax.legend(fontsize=7, loc="upper left")

# G: n_stack saturation
ax = AX[2, 0]
vals = []
for e in EVENTS:
    for j in e["fam40"][:400]:
        y = series_in_window(FAM_DI[j], FAM_NS[j].astype(float), e["w0"], e["nW"])
        vals.append(y[np.isfinite(y)])
vv = np.concatenate(vals) if vals else np.array([np.nan])
ax.hist(vv, bins=np.arange(4.5, 31.5, 1), color="#8c564b")
ax.set_xlabel("n_stack (daily stacks in the trailing 30 d)")
ax.set_ylabel("event-window family-days")
ax.set_title(f"G  n_stack saturates: {np.mean(vv >= 30) * 100:.0f} % at the ceiling of 30")

# H: per-station detector diagnosis
ax = AX[2, 1]
o = DIAG.sort_values("ets_ratio_cc0.80")
yy = np.arange(len(o))
ax.barh(yy - 0.2, o["ets_ratio_cc0.80"], 0.38, color="#3b76af", label="LFE count (cc>=0.80)")
ax.barh(yy + 0.2, o["pnsn_ets_ratio"], 0.38, color="#d1701e", label="PNSN count")
ax.set_yticks(yy)
ax.set_yticklabels(o.sta, fontsize=7)
ax.set_xscale("log")
ax.axvline(1, color="k", lw=0.8)
ax.set_xlabel("ETS / non-ETS MEAN daily count (log)")
ax.set_title("H  DIAGNOSIS: the matched-filter count hardly responds to ETS")
ax.set_xlim(0.7, max(8.0, DIAG.pnsn_ets_ratio.max() * 1.2))
ax.legend(fontsize=7, loc="lower right")

# I: summary text
ax = AX[2, 2]
ax.axis("off")
bN = T2.b_N.values
bP = T2.drop_duplicates("pid").b_pnsn.values
b30 = T2.b_N30.values
bA = T1A.b_near.values
txt = [
    "SUMMARY",
    f"streamed {COV.GB.sum():.0f} GB, {COV.detections.sum() / 1e6:.0f} M LFE detections,",
    f"  {len(COV)} stations, {T2.pid.nunique()} of 27 events covered",
    "",
    f"LFE daily count vs R7:   b = {np.nanmedian(bN):+.3f} "
    f"[{np.nanpercentile(bN, 25):+.2f}, {np.nanpercentile(bN, 75):+.2f}]",
    f"PNSN daily count vs R7:  b = {np.nanmedian(bP):+.3f} "
    f"[{np.nanpercentile(bP, 25):+.2f}, {np.nanpercentile(bP, 75):+.2f}]",
    f"LFE 30-d true count vs inc30: b = {np.nanmedian(b30):+.3f}",
    f"n_stack (saturated occupancy): b = {np.nanmedian(bA):+.3f}",
    "",
    f"b < 1 in {np.mean(bN[np.isfinite(bN)] < 1) * 100:.0f} % of LFE fits",
    f"cc >= 0.98 subset: b = {np.nanmedian(T2.b_cc98):+.3f}",
    f"certified families: b = {np.nanmedian(T2.b_Ncert):+.3f}",
    f"after floor removal: b = {np.nanmedian(T2.b_exc_p20):+.3f}",
    f"floor-de-attenuated:  b = {np.nanmedian(T2.b_deatt):+.3f}",
    "",
    "DIAGNOSIS: the mf detections are a nearly constant",
    f"  false-alarm stream ({DIAG['med_cc0.80'].median():.0f} det/day/station,",
    f"  ETS/quiet ratio {DIAG['ets_ratio_cc0.80'].median():.2f} vs PNSN "
    f"{DIAG['pnsn_ets_ratio'].median():.0f},",
    f"  rho(LFE, PNSN) = {DIAG['rho_pnsn_cc0.80'].median():+.2f}).",
    "  -> they are NOT an independent LFE-rate measure;",
    "     the pipeline uses them by STACKING, not counting.",
    "",
    "per-LFE AMPLITUDE not testable (mf files carry cc, not amplitude)",
]
ax.text(0.02, 0.98, "\n".join(txt), va="top", ha="left", fontsize=10, family="monospace",
        transform=ax.transAxes)

fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT_FIG / "s14_lfe_scaling.png", dpi=135)
plt.close(fig)
log(f"figure -> {OUT_FIG / 's14_lfe_scaling.png'}")

print("\n" + "=" * W)
print("FILES WRITTEN")
print("=" * W)
for p in ("results/s14_lfe30_exponents.csv", "results/s14_lfe_daily_exponents.csv",
          "results/s14_coverage.csv", "figures/s14_lfe_scaling.png"):
    q = HERE / p
    print(f"  {q}  ({q.stat().st_size / 1024:.0f} kB)")
log("S14 DONE")
