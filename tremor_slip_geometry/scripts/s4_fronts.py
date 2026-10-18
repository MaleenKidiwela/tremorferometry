#!/usr/bin/env python
"""S4 — along-strike FRONTS: daily tremor track vs daily slip-rate track, per confirmed event.

PLAN.md step S4, quantity Q6.  S1 showed that union footprints make the *event-averaged*
along-strike offset (Q2) erratic (thrust north median +23 km, IQR 52 km; mirror central +77 km,
IQR 125 km), because chained footprints put the two centroids far apart along strike even when
they co-register across strike.  S4 therefore never uses a footprint-averaged position: it builds
a DAILY track for each of tremor and modeled slip inside a fixed neighbourhood and compares the
two tracks in time (lead/lag, propagation speed, excursion, within-event depth evolution).

======================================================================================
INPUTS (all read-only, all OUTSIDE this folder except s1_events.csv)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv                      27 confirmed north/central
                                                                  patches (id, start, end, seg, Mw, ...)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz label (T x 742, int; 0=none else
                                                                  patch id), days (T, 'YYYY-MM-DD'),
                                                                  patch_lat/lon/depth_km/area_km2 (742)
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   'A_raw' (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                          PNSN tremor (time, lat, lon;
                                                                  the depth column is a placeholder
                                                                  and is IGNORED)
  /home/jovyan/cascadia-sse-strain                                mesh + sse_strain:
                                                                  load_fault_ascii / load_solution /
                                                                  component_slip_patterns -> Ld (742 x 13);
                                                                  megathrust.slab_depth_interpolator

  Slip:  raw_cum_mm = -1e3 * (Ld @ A_raw.T)          (742 x T), thrust positive, mm
         rate_raw[:, t] = raw_cum_mm[:, t] - raw_cum_mm[:, t-1]   (rate_raw[:, 0] = 0)
         RATE = clip( 7-day CENTRED moving average of rate_raw , 0, None )     mm/d, positive part

  Frame: x = (lon + 124) * 111 * cos(47 deg)  km east;  y = (lat - 47) * 111  km north.
         Along-strike coordinate = y (north positive).  Adequate for lead/lag INSIDE one event
         (the strike rotation matters for cross-margin comparisons, not for a single patch).
  Depth: Slab2 interface depth from slab_depth_interpolator, evaluated at triangle centroids
         (slip) and at tremor epicentres (tremor).  NaN outside the mesh hull (6.4 % of detections).

======================================================================================
EXACT DEFINITIONS  (per event; ts, te = day indices of start, end on the `days` grid)
======================================================================================
  Window   W = [ts - 45, te + 15] day indices, clipped to the grid.
  Footprint F = triangles whose label == id on at least one day in [ts, te].
  Neighbourhood NB = triangles whose centroid is within 60 km of AT LEAST ONE F centroid
             (F itself is always included).  NB is FIXED for the whole window.
             DIAGNOSTIC, because F is a union over [ts, te] and chains along strike:
             F_lat_min / F_lat_max / F_lat_span_deg and NB_y_span_km record how long the
             neighbourhood actually is.  For the largest events F spans most of the margin, so
             the "along-strike track" is a margin-wide centroid, not a local front.

  SLIP TRACK (one value per day t in W)
      S(t)   = sum over NB of area_km2 * RATE[:, t]                 (area-weighted positive rate)
      defined only where  S(t) >= 0.20 * max_{t in W} S(t)   and S(t) > 0
      y_s(t) = sum(area * RATE * y) / S(t)                          km, along-strike centroid
      z_s(t) = sum(area * RATE * slab_depth) / S(t)                 km, slip-weighted Slab2 depth

  TREMOR TRACK (one value per day t in W)
      detections within 60 km (epicentral, km frame) of at least one F centroid, dated on day t
      defined only on days with >= 5 such detections
      y_t(t) = mean y of those detections                           km
      z_t(t) = median Slab2 depth of those detections (finite values only)
      then 3-day CENTRED MEDIAN smoothing of y_t and z_t: the smoothed value on a defined day is
      the median of the defined values among days {t-1, t, t+1}; days that were undefined stay
      undefined (the smoother never creates new days).

  LEAD / LAG  L
      For each integer lag in -30..+30: Pearson r between y_s(t) and y_t(t + lag) over the days t
      where BOTH are defined; a lag is admissible only if >= 10 such days and both series have
      non-zero variance there.  lag_argmax = the admissible lag with the largest r.
      L = -lag_argmax  (days), so that POSITIVE L means the SLIP track LAGS the tremor track:
          if slip is a delayed copy, y_s(t) = y_t(t - D), then r(lag) peaks at lag = -D, i.e. L = +D.
      Reported: L_days, corr_at_L (= r at lag_argmax), corr_lag0, n_overlap_at_L, n_overlap_lag0.

  ONSET LAG = (first day in W where the slip track is defined)
            - (first day in W where the tremor track is defined)          days (+ = slip starts later)

  SPEEDS    v_t = Theil-Sen slope of y_t vs day index over the defined tremor days   km/d
            v_s = Theil-Sen slope of y_s vs day index over the defined slip days     km/d
            (sign = propagation direction, + = northward).  Requires >= 5 defined days.
            ratio = v_s / v_t (NaN if |v_t| < 1e-6).
            excursion_t = max(y_t) - min(y_t);  excursion_s = max(y_s) - min(y_s)     km
            (supplementary robust versions excursion_*_p10_90_km = p90 - p10 of the track, because
             the full range is set by single sparse days at the far end of a chained footprint)

  WITHIN-EVENT DEPTH EVOLUTION
            split the DEFINED days of each track into three consecutive equal-count thirds;
            dz_t = median(z_t over the last third) - median(z_t over the first third)   km (+ = deepening)
            dz_s = same for z_s.  Requires >= 6 defined days.

  MIXTURE DIAGNOSTIC (how much freedom the 13-pattern model had)
            per component k:  cum_k(mm) = -1e3 * Ld[:, k] outer A_raw[:, k]; the same
            diff -> 7-d centred mean -> positive-part pipeline; then
            S_k = sum over NB and over the days of W of  area_km2 * positive rate of component k.
            frac_k = S_k / sum_j S_j;   n_components = #{k : frac_k >= 0.10}.
            Also reported: top_component, top_component_frac.
            (Normalisation is by the SUM OF THE 13 POSITIVE PER-COMPONENT contributions, which is
             >= the positive part of their sum, so frac_k is a share of component activity, not of
             the net slip.)

  CONTROL (per event, 50 draws)
            The tremor detection dates of the SAME 60-km neighbourhood are rigidly shifted by a
            random integer of >= 365 and <= (span - 365) days, wrapped modulo the tremor span
            2010-01-01 .. 2026-05-27; the spatial mask is unchanged.  The tremor track is rebuilt
            with the identical >= 5-detections / 3-day-median rules and L is recomputed against the
            UNCHANGED slip track.  Records the null distribution of |corr at best lag| (and of L).
            p_null = fraction of successful draws with |r_null| >= |r_observed|.
            Two group-level readings of the same control:
              (i) is the CORRELATION meaningful?  compare the observed |r| distribution with the
                  pooled null |r|, and mean(p_null) (0.5 under the null);
              (ii) is the LAG meaningful?  frac_absL_le5 observed vs the pooled null, with a
                   binomial p (k of n events with |L| <= 5 d against the pooled null probability).

  AGGREGATE (combined / north / central)
            median and IQR (q25, q75) of L, onset lag, v_t, v_s, ratio, excursions, dz_t, dz_s;
            fraction of events with sign(v_t) == sign(v_s);
            Theil-Sen trend vs event date (per year) with a 1000-permutation two-sided p on |slope|
            for L, v_t, v_s, ratio; Spearman rho vs Mw for the same quantities.

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s4_fronts_events.csv    one row per event: id, start, end, seg, Mw, window/track sizes,
                                  L_days, lag_argmax, corr_at_L, corr_lag0, n_overlap*, onset_lag_d,
                                  v_t_km_d, v_s_km_d, v_ratio, excursion_t_km, excursion_s_km,
                                  excursion_*_p10_90_km, F_lat_min/max/span, NB_y_span_km,
                                  z_t first/last third + dz_t, z_s first/last third + dz_s,
                                  n_components, top_component(+frac), null_absr_mean/p95/p_null,
                                  n_null_ok, dir_agree, note
  results/s4_fronts_summary.csv   aggregate rows (group, quantity, kind, n, median/q25/q75 or
                                  slope/CI/p_perm or rho/p)
  figures/s4_fronts.png           8 largest-Mw events: y_t and y_s vs days from event start;
                                  plus histogram of L with the shifted-tremor null, v_s vs v_t
                                  with the 1:1 line, and L vs date.

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s4_fronts.py
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
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree
from scipy.stats import binomtest, spearmanr, theilslopes

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)

PRE_DAYS = 45           # window starts 45 d before the slip start
POST_DAYS = 15          # window ends 15 d after the slip end
NB_RADIUS_KM = 60.0     # neighbourhood / tremor radius around footprint centroids
RATE_FRAC = 0.20        # slip track defined where neighbourhood rate >= 20 % of its window max
MIN_DET_DAY = 5         # tremor track defined on days with >= 5 detections
SMOOTH_TREMOR = 3       # centred median smoothing of the tremor track, days
LAG_MAX = 30            # lead/lag search range, days
MIN_OVERLAP = 10        # minimum overlapping days for a lag to be admissible
MIN_TRACK_SPEED = 5     # minimum defined days for a Theil-Sen speed
MIN_TRACK_DEPTH = 6     # minimum defined days for the first/last-third depth change
COMP_FRAC = 0.10        # component share threshold for the mixture diagnostic
N_SHIFT = 50            # shifted-tremor control draws per event
N_PERM = 1000           # permutations for the trend p-values
RATE_SMOOTH_D = 7       # centred moving average applied to the daily slip rate
COS47 = np.cos(np.deg2rad(47.0))

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


# ----------------------------------------------------------------------------------
# 1. Slip model  (Ld, A_raw, Slab2)
# ----------------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------------
# 2. Geometry / labels
# ----------------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------------
# 3. Tremor
# ----------------------------------------------------------------------------------
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DEPTH = slabf(tr["lon"].values, tr["lat"].values)
SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1
# tremor day expressed on the model day grid (index into DAYS), via days since SHIFT_LO
TR_OFF = (tr["t"].values.astype("datetime64[D]") - SHIFT_LO).astype(int)  # 0 .. SPAN-1
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(
    f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}; "
    f"Slab2 depth finite for {np.isfinite(TR_DEPTH).mean() * 100:.1f}%; wrap span {SPAN} d"
)


# ----------------------------------------------------------------------------------
# 4. Events (S1 confirmed set)
# ----------------------------------------------------------------------------------
ev_in = pd.read_csv(OUT_RES / "s1_events.csv")
ev_in = ev_in.sort_values("start").reset_index(drop=True)
log(f"  events {len(ev_in)}  {ev_in.seg.value_counts().to_dict()}")


# ----------------------------------------------------------------------------------
# 5. Track machinery
# ----------------------------------------------------------------------------------
def centred_median(vals: np.ndarray, defined: np.ndarray, half: int = 1) -> np.ndarray:
    """Centred median over +-half days, using only defined neighbours; undefined days stay NaN."""
    out = np.full(vals.shape, np.nan)
    n = len(vals)
    for i in range(n):
        if not defined[i]:
            continue
        lo, hi = max(0, i - half), min(n, i + half + 1)
        w = vals[lo:hi][defined[lo:hi]]
        if w.size:
            out[i] = np.median(w)
    return out


def tremor_track(offsets: np.ndarray, y: np.ndarray, z: np.ndarray, w0: int, w1: int):
    """Daily tremor centroid track over model-day indices [w0, w1].

    offsets: day of each detection expressed as an index into DAYS (already shifted if control).
    Returns (y_t, z_t, ndet) arrays of length nW, NaN where fewer than MIN_DET_DAY detections.
    """
    nW = w1 - w0 + 1
    yt = np.full(nW, np.nan)
    zt = np.full(nW, np.nan)
    ndet = np.zeros(nW, int)
    sel = (offsets >= w0) & (offsets <= w1)
    if not sel.any():
        return yt, zt, ndet
    loc = offsets[sel] - w0
    ys, zs = y[sel], z[sel]
    order = np.argsort(loc, kind="stable")
    loc, ys, zs = loc[order], ys[order], zs[order]
    bounds = np.searchsorted(loc, np.arange(nW + 1))
    for i in range(nW):
        a, b = bounds[i], bounds[i + 1]
        ndet[i] = b - a
        if b - a >= MIN_DET_DAY:
            yt[i] = ys[a:b].mean()
            zz = zs[a:b]
            zz = zz[np.isfinite(zz)]
            if zz.size:
                zt[i] = np.median(zz)
    dfn = np.isfinite(yt)
    yt_s = centred_median(yt, dfn, SMOOTH_TREMOR // 2)
    zt_s = centred_median(zt, np.isfinite(zt), SMOOTH_TREMOR // 2)
    zt_s[~dfn] = np.nan  # depth track lives on the same days as the position track
    return yt_s, zt_s, ndet


def best_lag(ys: np.ndarray, yt: np.ndarray):
    """Lead/lag between the slip track ys(t) and the tremor track yt(t).

    Scans corr(ys[t], yt[t+lag]) for lag in -LAG_MAX..+LAG_MAX; returns
    (L, corr_at_L, n_at_L, corr_lag0, n_lag0, lag_argmax) with L = -lag_argmax so that
    POSITIVE L = slip lags tremor.  NaNs where no admissible lag exists.
    """
    n = len(ys)
    best = (np.nan, -np.inf, 0)
    c0, n0 = np.nan, 0
    for lag in range(-LAG_MAX, LAG_MAX + 1):
        ia = np.arange(n)
        ib = ia + lag
        ok = (ib >= 0) & (ib < n)
        ia, ib = ia[ok], ib[ok]
        m = np.isfinite(ys[ia]) & np.isfinite(yt[ib])
        k = int(m.sum())
        if k < MIN_OVERLAP:
            continue
        a, b = ys[ia[m]], yt[ib[m]]
        if a.std() < 1e-9 or b.std() < 1e-9:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        if lag == 0:
            c0, n0 = r, k
        if r > best[1]:
            best = (lag, r, k)
    if not np.isfinite(best[1]):
        return np.nan, np.nan, 0, c0, n0, np.nan
    return -best[0], best[1], best[2], c0, n0, best[0]


def ts_slope(x: np.ndarray, y: np.ndarray, nmin: int):
    if len(x) < nmin:
        return np.nan
    return float(theilslopes(y, x, 0.95)[0])


def third_change(z: np.ndarray):
    """median(last third) - median(first third) over the DEFINED days of z."""
    v = z[np.isfinite(z)]
    if v.size < MIN_TRACK_DEPTH:
        return np.nan, np.nan, np.nan
    k = v.size // 3
    if k < 1:
        return np.nan, np.nan, np.nan
    first, last = float(np.median(v[:k])), float(np.median(v[-k:]))
    return first, last, last - first


# ----------------------------------------------------------------------------------
# 6. Per-event measurement
# ----------------------------------------------------------------------------------
log("measuring events ...")
rows = []
TRACKS = {}
NULL_R = {}
NULL_L = {}

for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_DAYS), min(NT - 1, te + POST_DAYS)
    nW = w1 - w0 + 1
    rec = dict(id=pid, start=r["start"], end=r["end"], seg=r["seg"], Mw=float(r["Mw"]),
               dur_d=int(r["dur_d"]), N_s1=int(r["N"]), w_start=DAYS[w0], w_end=DAYS[w1],
               n_days_W=nW, note="")

    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    rec["n_tri_F"] = int(F.sum())
    rec["F_lat_min"] = float(PLAT[F].min())
    rec["F_lat_max"] = float(PLAT[F].max())
    rec["F_lat_span_deg"] = rec["F_lat_max"] - rec["F_lat_min"]
    fx, fy = PX[F], PY[F]

    # ---- neighbourhood: triangles within 60 km of any F centroid ----
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = np.unique(np.concatenate([np.asarray(h, int) for h in hits])) if len(hits) else np.array([], int)
    nb = np.union1d(nb, np.flatnonzero(F))
    rec["n_tri_NB"] = int(nb.size)
    rec["NB_y_span_km"] = float(np.ptp(PY[nb]))

    # ---- slip track ----
    Rw = RATE[np.ix_(nb, np.arange(w0, w1 + 1))]        # (nb, nW) mm/d
    aw = PAREA[nb][:, None]
    Sw = (aw * Rw).sum(axis=0)                          # area-weighted rate, km^2*mm/d
    smax = Sw.max()
    ok_s = (Sw >= RATE_FRAC * smax) & (Sw > 0) if smax > 0 else np.zeros(nW, bool)
    ys = np.full(nW, np.nan)
    zs = np.full(nW, np.nan)
    if ok_s.any():
        num_y = (aw * Rw * PY[nb][:, None]).sum(axis=0)
        num_z = (aw * Rw * PDEPTH[nb][:, None]).sum(axis=0)
        ys[ok_s] = num_y[ok_s] / Sw[ok_s]
        zs[ok_s] = num_z[ok_s] / Sw[ok_s]
    rec["n_days_slip_track"] = int(np.isfinite(ys).sum())
    rec["slip_rate_max_mm_d"] = float(Rw.max()) if Rw.size else np.nan

    # ---- tremor track ----
    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near = (np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
            if len(thits) else np.array([], int))
    rec["n_tremor_near_alltime"] = int(near.size)
    if near.size:
        # detection day as an index into DAYS
        off_idx = TR_OFF[near] + IDX_OF_SHIFT_LO
        yt, zt, ndet = tremor_track(off_idx, TR_Y[near], TR_DEPTH[near], w0, w1)
    else:
        yt = np.full(nW, np.nan); zt = np.full(nW, np.nan); ndet = np.zeros(nW, int)
    rec["n_days_tremor_track"] = int(np.isfinite(yt).sum())
    rec["n_det_in_W"] = int(ndet.sum())

    # ---- lead/lag ----
    L, cL, nL, c0, n0, lagmax = best_lag(ys, yt)
    rec.update(L_days=L, lag_argmax=lagmax, corr_at_L=cL, n_overlap_at_L=nL,
               corr_lag0=c0, n_overlap_lag0=n0)

    # ---- onset lag ----
    fs = np.flatnonzero(np.isfinite(ys))
    ft = np.flatnonzero(np.isfinite(yt))
    rec["onset_lag_d"] = float(fs[0] - ft[0]) if (fs.size and ft.size) else np.nan
    rec["slip_onset"] = DAYS[w0 + fs[0]] if fs.size else ""
    rec["tremor_onset"] = DAYS[w0 + ft[0]] if ft.size else ""

    # ---- speeds and excursions ----
    tt = np.arange(nW, dtype=float)
    rec["v_t_km_d"] = ts_slope(tt[ft], yt[ft], MIN_TRACK_SPEED) if ft.size else np.nan
    rec["v_s_km_d"] = ts_slope(tt[fs], ys[fs], MIN_TRACK_SPEED) if fs.size else np.nan
    vt, vs = rec["v_t_km_d"], rec["v_s_km_d"]
    rec["v_ratio"] = (vs / vt) if (np.isfinite(vt) and np.isfinite(vs) and abs(vt) > 1e-6) else np.nan
    rec["excursion_t_km"] = float(np.ptp(yt[ft])) if ft.size else np.nan
    rec["excursion_s_km"] = float(np.ptp(ys[fs])) if fs.size else np.nan
    rec["excursion_t_p10_90_km"] = (
        float(np.percentile(yt[ft], 90) - np.percentile(yt[ft], 10)) if ft.size else np.nan)
    rec["excursion_s_p10_90_km"] = (
        float(np.percentile(ys[fs], 90) - np.percentile(ys[fs], 10)) if fs.size else np.nan)
    rec["dir_agree"] = (
        float(np.sign(vt) == np.sign(vs)) if (np.isfinite(vt) and np.isfinite(vs)) else np.nan
    )

    # ---- within-event depth evolution ----
    rec["z_t_first3_km"], rec["z_t_last3_km"], rec["dz_t_km"] = third_change(zt)
    rec["z_s_first3_km"], rec["z_s_last3_km"], rec["dz_s_km"] = third_change(zs)

    # ---- mixture diagnostic ----
    e0, e1 = max(0, w0 - 12), min(NT - 1, w1 + 12)      # margin so the 7-d filter matches
    Sk = np.zeros(A_raw.shape[1])
    for k in range(A_raw.shape[1]):
        cum_k = -1e3 * np.outer(Ld[nb, k], A_raw[e0:e1 + 1, k])
        rk = rate_pipeline(cum_k)[:, (w0 - e0):(w1 - e0 + 1)]
        Sk[k] = float((PAREA[nb][:, None] * rk).sum())
    tot = Sk.sum()
    fk = Sk / tot if tot > 0 else np.zeros_like(Sk)
    rec["n_components"] = int((fk >= COMP_FRAC).sum())
    rec["top_component"] = int(np.argmax(fk)) + 1
    rec["top_component_frac"] = float(fk.max())
    rec["comp_fracs"] = ";".join(f"{v:.3f}" for v in fk)

    # ---- shifted-tremor control ----
    absr = []
    lnull = []
    if near.size:
        base_off = TR_OFF[near]
        shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
        for sh in shifts:
            off_idx = ((base_off + int(sh)) % SPAN) + IDX_OF_SHIFT_LO
            ytn, _, _ = tremor_track(off_idx, TR_Y[near], TR_DEPTH[near], w0, w1)
            Ln, cn, _, _, _, _ = best_lag(ys, ytn)
            if np.isfinite(cn):
                absr.append(abs(cn))
                lnull.append(Ln)
    absr = np.asarray(absr, float)
    lnull = np.asarray(lnull, float)
    NULL_R[pid] = absr
    NULL_L[pid] = lnull
    rec["n_null_ok"] = int(absr.size)
    rec["null_absr_mean"] = float(absr.mean()) if absr.size else np.nan
    rec["null_absr_p95"] = float(np.percentile(absr, 95)) if absr.size else np.nan
    rec["p_null"] = (
        float((absr >= abs(cL)).mean()) if (absr.size and np.isfinite(cL)) else np.nan
    )

    if not np.isfinite(L):
        rec["note"] = "no admissible lag (track overlap < 10 d)"

    TRACKS[pid] = dict(w0=w0, ts=ts, ys=ys, yt=yt, zs=zs, zt=zt, ndet=ndet)
    rows.append(rec)
    log(
        f"  {pid:5d} {r['start']} {r['seg']:7s} Mw{float(r['Mw']):.2f}  "
        f"nb={rec['n_tri_NB']:3d} dS={rec['n_days_slip_track']:3d} dT={rec['n_days_tremor_track']:3d} "
        f"L={L if np.isfinite(L) else float('nan'):6.1f} r={cL:6.3f} "
        f"vt={vt:7.3f} vs={vs:7.3f} nc={rec['n_components']}"
    )

events = pd.DataFrame(rows)


# ----------------------------------------------------------------------------------
# 7. Aggregate
# ----------------------------------------------------------------------------------
QUANTS = ["L_days", "onset_lag_d", "v_t_km_d", "v_s_km_d", "v_ratio",
          "excursion_t_km", "excursion_s_km",
          "excursion_t_p10_90_km", "excursion_s_p10_90_km",
          "dz_t_km", "dz_s_km", "corr_at_L", "corr_lag0",
          "n_components", "F_lat_span_deg", "NB_y_span_km"]
TRENDQ = ["L_days", "v_t_km_d", "v_s_km_d", "v_ratio"]


def year_float(dates):
    dt = pd.to_datetime(pd.Series(list(dates)))
    return (dt - pd.Timestamp("2010-01-01")).dt.total_seconds().values / (365.25 * 86400.0)


def trend(sub, col):
    s = sub[np.isfinite(sub[col])]
    if len(s) < 5:
        return None
    x = year_float(s["start"])
    y = s[col].values.astype(float)
    sl, ic, lo, hi = theilslopes(y, x, 0.95)
    perm = np.array([theilslopes(y, RNG.permutation(x), 0.95)[0] for _ in range(N_PERM)])
    rho, ps = spearmanr(x, y)
    return dict(n=len(s), slope=sl, lo=lo, hi=hi,
                p_perm=float((np.abs(perm) >= abs(sl)).mean()), rho=float(rho), p=float(ps))


srows = []
groups = [("combined", events), ("north", events[events.seg == "north"]),
          ("central", events[events.seg == "central"])]
for gname, sub in groups:
    for q in QUANTS:
        v = sub[q].values.astype(float)
        v = v[np.isfinite(v)]
        srows.append(dict(group=gname, quantity=q, kind="distribution", n=len(v),
                          median=np.median(v) if v.size else np.nan,
                          q25=np.percentile(v, 25) if v.size else np.nan,
                          q75=np.percentile(v, 75) if v.size else np.nan,
                          mean=v.mean() if v.size else np.nan))
    da = sub["dir_agree"].values.astype(float)
    da = da[np.isfinite(da)]
    srows.append(dict(group=gname, quantity="dir_agree_frac", kind="fraction", n=len(da),
                      median=da.mean() if da.size else np.nan))
    # observed vs null correlation
    obs = sub["corr_at_L"].abs().values.astype(float)
    obs = obs[np.isfinite(obs)]
    nul = np.concatenate([NULL_R[p] for p in sub["id"] if NULL_R.get(p, np.array([])).size]) \
        if len(sub) else np.array([])
    srows.append(dict(group=gname, quantity="abs_corr_at_L_observed", kind="distribution",
                      n=len(obs), median=np.median(obs) if obs.size else np.nan,
                      q25=np.percentile(obs, 25) if obs.size else np.nan,
                      q75=np.percentile(obs, 75) if obs.size else np.nan,
                      mean=obs.mean() if obs.size else np.nan))
    srows.append(dict(group=gname, quantity="abs_corr_at_L_null", kind="distribution",
                      n=len(nul), median=np.median(nul) if nul.size else np.nan,
                      q25=np.percentile(nul, 25) if nul.size else np.nan,
                      q75=np.percentile(nul, 75) if nul.size else np.nan,
                      mean=nul.mean() if nul.size else np.nan))
    pn = sub["p_null"].values.astype(float)
    pn = pn[np.isfinite(pn)]
    srows.append(dict(group=gname, quantity="frac_events_p_null_lt_0.05", kind="fraction",
                      n=len(pn), median=float((pn < 0.05).mean()) if pn.size else np.nan))
    srows.append(dict(group=gname, quantity="mean_p_null (0.5 = no signal)", kind="fraction",
                      n=len(pn), median=float(pn.mean()) if pn.size else np.nan))
    # is the LAG itself concentrated?  observed |L| <= 5 d vs the pooled shifted-tremor null
    Lo = sub["L_days"].values.astype(float); Lo = Lo[np.isfinite(Lo)]
    Ln = np.concatenate([NULL_L[q] for q in sub["id"] if NULL_L.get(q, np.array([])).size]) \
        if len(sub) else np.array([])
    Ln = Ln[np.isfinite(Ln)]
    if Lo.size and Ln.size:
        f_obs = float((np.abs(Lo) <= 5).mean())
        f_nul = float((np.abs(Ln) <= 5).mean())
        pb = binomtest(int((np.abs(Lo) <= 5).sum()), len(Lo), f_nul, alternative="greater").pvalue
        srows.append(dict(group=gname, quantity="frac_absL_le5_observed", kind="fraction",
                          n=len(Lo), median=f_obs))
        srows.append(dict(group=gname, quantity="frac_absL_le5_null", kind="fraction",
                          n=len(Ln), median=f_nul))
        srows.append(dict(group=gname, quantity="p_binom_absL_le5_obs_gt_null", kind="fraction",
                          n=len(Lo), median=float(pb)))
        srows.append(dict(group=gname, quantity="median_abs_L_observed", kind="fraction",
                          n=len(Lo), median=float(np.median(np.abs(Lo)))))
        srows.append(dict(group=gname, quantity="median_abs_L_null", kind="fraction",
                          n=len(Ln), median=float(np.median(np.abs(Ln)))))
    for q in TRENDQ:
        t = trend(sub, q)
        if t:
            srows.append(dict(group=gname, quantity=q, kind="trend_vs_date_per_yr",
                              n=t["n"], slope=t["slope"], lo=t["lo"], hi=t["hi"],
                              p_perm=t["p_perm"], rho=t["rho"], p=t["p"]))
        s2 = sub[np.isfinite(sub[q]) & np.isfinite(sub["Mw"])]
        if len(s2) >= 5:
            rho, p = spearmanr(s2["Mw"].values, s2[q].values)
            srows.append(dict(group=gname, quantity=q, kind="spearman_vs_Mw",
                              n=len(s2), rho=float(rho), p=float(p)))
summary = pd.DataFrame(srows)


# ----------------------------------------------------------------------------------
# 8. Write CSVs
# ----------------------------------------------------------------------------------
EV_COLS = ["id", "start", "end", "seg", "Mw", "dur_d", "w_start", "w_end", "n_days_W",
           "n_tri_F", "n_tri_NB", "n_tremor_near_alltime", "n_det_in_W",
           "n_days_slip_track", "n_days_tremor_track", "slip_rate_max_mm_d",
           "L_days", "lag_argmax", "corr_at_L", "n_overlap_at_L", "corr_lag0", "n_overlap_lag0",
           "onset_lag_d", "slip_onset", "tremor_onset",
           "v_t_km_d", "v_s_km_d", "v_ratio", "dir_agree",
           "excursion_t_km", "excursion_s_km",
           "excursion_t_p10_90_km", "excursion_s_p10_90_km",
           "F_lat_min", "F_lat_max", "F_lat_span_deg", "NB_y_span_km",
           "z_t_first3_km", "z_t_last3_km", "dz_t_km",
           "z_s_first3_km", "z_s_last3_km", "dz_s_km",
           "n_components", "top_component", "top_component_frac", "comp_fracs",
           "n_null_ok", "null_absr_mean", "null_absr_p95", "p_null", "N_s1", "note"]
events = events.reindex(columns=EV_COLS)
events.to_csv(OUT_RES / "s4_fronts_events.csv", index=False)
summary.to_csv(OUT_RES / "s4_fronts_summary.csv", index=False)
log(f"wrote results/s4_fronts_events.csv ({len(events)}) and s4_fronts_summary.csv ({len(summary)})")


# ----------------------------------------------------------------------------------
# 9. Figure
# ----------------------------------------------------------------------------------
top8 = events.sort_values("Mw", ascending=False).head(8).sort_values("start")
fig = plt.figure(figsize=(19, 13))
gs = fig.add_gridspec(3, 12, hspace=0.42, wspace=0.55,
                      left=0.05, right=0.985, top=0.93, bottom=0.06)
for j, (_, e) in enumerate(top8.iterrows()):
    ax = fig.add_subplot(gs[j // 4, 3 * (j % 4):3 * (j % 4) + 3])
    tk = TRACKS[int(e["id"])]
    d = np.arange(len(tk["ys"])) + (tk["w0"] - tk["ts"])   # days from event start
    ax.plot(d, tk["yt"], ".-", color="#c1272d", ms=4, lw=1.2, label="tremor $y_t$")
    ax.plot(d, tk["ys"], ".-", color="#1f5fa9", ms=4, lw=1.2, label="slip $y_s$")
    ax.axvline(0, color="k", lw=0.7, ls=":")
    ax.axvline(e["dur_d"], color="k", lw=0.7, ls=":")
    ax.set_title(f"{int(e['id'])} {e['start']} {e['seg'][:1].upper()} Mw{e['Mw']:.2f}\n"
                 f"L={e['L_days']:+.0f} d  r={e['corr_at_L']:.2f}  "
                 f"$v_t$={e['v_t_km_d']:+.2f} $v_s$={e['v_s_km_d']:+.2f} km/d", fontsize=8.5)
    ax.set_xlabel("days from slip start", fontsize=8)
    ax.set_ylabel("along-strike y (km)", fontsize=8)
    ax.tick_params(labelsize=7.5)
    if j == 0:
        ax.legend(fontsize=7.5, loc="best")

# panel A: histogram of L with null
axA = fig.add_subplot(gs[2, 0:4])
Lobs = events["L_days"].dropna().values
Lnul = np.concatenate([v for v in NULL_L.values() if v.size]) if NULL_L else np.array([])
bins = np.arange(-30.5, 31.5, 3)
if Lnul.size:
    axA.hist(Lnul, bins=bins, density=True, color="0.75", label=f"shifted-tremor null (n={Lnul.size})")
axA.hist(Lobs, bins=bins, density=True, histtype="step", lw=2.0, color="#1f5fa9",
         label=f"observed (n={Lobs.size})")
axA.axvline(0, color="k", lw=0.8, ls=":")
if Lobs.size:
    axA.axvline(np.median(Lobs), color="#c1272d", lw=1.8,
                label=f"median L = {np.median(Lobs):+.0f} d")
axA.set_xlabel("L (days); + = slip lags tremor", fontsize=9)
axA.set_ylabel("density", fontsize=9)
_fa = summary[(summary.group == "combined") & (summary.kind == "fraction")].set_index("quantity")["median"]
axA.set_title("A  lead/lag L, observed vs null\n"
              f"|L|<=5 d: {_fa.get('frac_absL_le5_observed', np.nan) * 100:.0f}% obs vs "
              f"{_fa.get('frac_absL_le5_null', np.nan) * 100:.0f}% null "
              f"(binom p={_fa.get('p_binom_absL_le5_obs_gt_null', np.nan):.1e})", fontsize=10)
axA.legend(fontsize=7.5)
axA.tick_params(labelsize=8)

# panel B: v_s vs v_t
axB = fig.add_subplot(gs[2, 4:8])
for seg, c in (("north", "#1f5fa9"), ("central", "#2e8b57")):
    s = events[events.seg == seg]
    axB.scatter(s["v_t_km_d"], s["v_s_km_d"], s=42, c=c, alpha=0.85, edgecolor="k",
                linewidth=0.4, label=f"{seg} (n={len(s)})")
lim = np.nanmax(np.abs(np.concatenate([events["v_t_km_d"].values, events["v_s_km_d"].values])))
lim = 1.1 * (lim if np.isfinite(lim) else 1.0)
axB.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="1:1")
axB.axhline(0, color="0.6", lw=0.6); axB.axvline(0, color="0.6", lw=0.6)
axB.set_xlim(-lim, lim); axB.set_ylim(-lim, lim)
axB.set_xlabel("$v_t$ tremor along-strike speed (km/d)", fontsize=9)
axB.set_ylabel("$v_s$ slip along-strike speed (km/d)", fontsize=9)
da = events["dir_agree"].dropna()
axB.set_title(f"B  propagation speeds (sign agreement {da.mean() * 100:.0f}% of {len(da)})",
              fontsize=10)
axB.legend(fontsize=7.5)
axB.tick_params(labelsize=8)

# panel C: L vs date
axC = fig.add_subplot(gs[2, 8:12])
dts = pd.to_datetime(events["start"])
for seg, c in (("north", "#1f5fa9"), ("central", "#2e8b57")):
    s = events[events.seg == seg]
    axC.scatter(pd.to_datetime(s["start"]), s["L_days"], s=42, c=c, alpha=0.85,
                edgecolor="k", linewidth=0.4, label=seg)
tr_row = summary[(summary.group == "combined") & (summary.quantity == "L_days")
                 & (summary.kind == "trend_vs_date_per_yr")]
if len(tr_row):
    sl = float(tr_row["slope"].iloc[0]); pp = float(tr_row["p_perm"].iloc[0])
    sub = events[np.isfinite(events["L_days"])]
    xf = year_float(sub["start"])
    ic = np.median(sub["L_days"].values - sl * xf)
    xx = np.linspace(xf.min(), xf.max(), 10)
    axC.plot(pd.Timestamp("2010-01-01") + pd.to_timedelta(xx * 365.25, "D"), sl * xx + ic,
             "r-", lw=1.5, label=f"Theil-Sen {sl:+.2f} d/yr (p={pp:.2f})")
axC.axhline(0, color="k", lw=0.8, ls=":")
axC.set_xlabel("event start date", fontsize=9)
axC.set_ylabel("L (days)", fontsize=9)
axC.set_title("C  lead/lag vs time", fontsize=10)
axC.legend(fontsize=7.5)
axC.tick_params(labelsize=8)

fig.suptitle("S4 — along-strike fronts: daily tremor track vs daily slip-rate track "
             "(27 confirmed north/central patches, 60 km neighbourhood)", fontsize=13)
fig.savefig(OUT_FIG / "s4_fronts.png", dpi=135)
plt.close(fig)
log("wrote figures/s4_fronts.png")


# ----------------------------------------------------------------------------------
# 10. Print tables
# ----------------------------------------------------------------------------------
pd.set_option("display.width", 220, "display.max_columns", 60)

print("\n=== PER-EVENT (S4 fronts) ===")
pe = events[["id", "start", "seg", "Mw", "n_days_slip_track", "n_days_tremor_track",
             "L_days", "corr_at_L", "corr_lag0", "p_null", "onset_lag_d",
             "v_t_km_d", "v_s_km_d", "v_ratio", "dir_agree", "excursion_t_km",
             "excursion_s_km", "dz_t_km", "dz_s_km", "n_components", "F_lat_span_deg"]].copy()
pe.columns = ["id", "start", "seg", "Mw", "dS", "dT", "L", "r@L", "r@0", "p_null", "onset",
              "v_t", "v_s", "ratio", "agree", "exc_t", "exc_s", "dz_t", "dz_s", "ncomp", "Fspan"]
print(pe.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))

print("\n=== AGGREGATE: distributions (median [q25, q75]) ===")
dist = summary[summary.kind == "distribution"]
print(f"{'group':9s} {'quantity':26s} {'n':>4s} {'median':>10s} {'q25':>10s} {'q75':>10s}")
for _, s in dist.iterrows():
    print(f"{s.group:9s} {s.quantity:26s} {int(s.n):4d} {s['median']:10.3f} "
          f"{s['q25']:10.3f} {s['q75']:10.3f}")

print("\n=== AGGREGATE: fractions ===")
for _, s in summary[summary.kind == "fraction"].iterrows():
    print(f"{s.group:9s} {s.quantity:32s} n={int(s.n):3d}  {s['median']:.3f}")

print("\n=== AGGREGATE: trends vs date (Theil-Sen per year, 1000-permutation p) ===")
print(f"{'group':9s} {'quantity':14s} {'n':>4s} {'slope/yr':>10s} {'CI95':>24s} {'p_perm':>8s} "
      f"{'rho':>7s} {'p_sp':>8s}")
for _, s in summary[summary.kind == "trend_vs_date_per_yr"].iterrows():
    print(f"{s.group:9s} {s.quantity:14s} {int(s.n):4d} {s['slope']:10.4f} "
          f"[{s['lo']:10.4f},{s['hi']:10.4f}] {s['p_perm']:8.3f} {s['rho']:7.2f} {s['p']:8.3f}")

print("\n=== AGGREGATE: Spearman vs Mw ===")
for _, s in summary[summary.kind == "spearman_vs_Mw"].iterrows():
    print(f"{s.group:9s} {s.quantity:14s} n={int(s.n):3d}  rho={s['rho']:+.3f}  p={s['p']:.3f}")

log("done.")
