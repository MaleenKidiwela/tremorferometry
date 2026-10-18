#!/usr/bin/env python
"""S9b — reconciling S9's within-event DEEPENING with Wech & Creager (2011)'s UP-DIP migration.

Coordinator follow-up to S9 Part A.  S9 found the stacked daily median Slab2 depth of the
near-footprint tremor rises (deepens) by +2.26 km from the first to the last fifth of an event
(200-shift null p = 0.020), and that the daily p10 (up-dip-most) tremor deepens +2.94 km.
Wech & Creager (2011, Nature Geoscience 4, 624) plotted median daily DOWN-DIP DISTANCE through
four ETS episodes and reported initiation far down-dip followed by UP-DIP migration during the
episode -- apparently the opposite sign.

The S9 five-bin profile already contains the likely answer: +0.44, -0.41, -1.08, +1.33, +2.70 km,
i.e. SHALLOWING over the first three fifths and DEEPENING over the last two.  S9b tests that
directly and checks every way the S9 sign could have been manufactured.

  (1) 10-bin phase profile, plus an ONSET / PLATEAU / TAIL decomposition, each with its own
      shifted-window null.  Reported: depth change onset -> peak and peak -> end.
  (2) Wech & Creager's own statistic: median daily DOWN-DIP DISTANCE (km from the mesh's
      shallowest up-dip edge, measured horizontally) vs DAYS FROM EVENT START, stacked over the
      27 events -- no phase normalisation, the way they plotted it.  Plus the two events of ours
      that can overlap their 2010-2011 Washington ETS.
  (3) Sensitivity of the SIGN to the event-start definition (first tremor day vs catalog slip
      start) and to the neighbourhood radius (40 / 60 / 100 km), each with a null.
  (4) Spatial control on "shallower when most intense": the along-strike position of the daily
      tremor centroid is regressed out of BOTH the daily depth and the daily count (rank-based
      partial Spearman), because depth and amplitude co-vary ALONG STRIKE (Wech 2021), which
      would confound a purely temporal reading.

======================================================================================
INPUTS  (identical loaders to scripts/s4_fronts.py and scripts/s9_within_event.py)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv                       27 confirmed north/central patches
  tremor_slip_geometry/results/s9_depth_migration.csv              S9 per-event table (cross-check)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz  label, days, patch_lat/lon/area
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw
  catalogs/pnsn_tremor_cascadia_full.csv                           PNSN tremor
  /home/jovyan/cascadia-sse-strain                                 mesh + slab_depth_interpolator

======================================================================================
DEFINITIONS
======================================================================================
  Window WA = [start - 10, end + 10] (as in S9).  Footprint F, neighbourhood radius R (default
  60 km), near tremor, daily count n(d) and daily median Slab2 depth zmed(d) exactly as in S9;
  a day is DEFINED at n(d) >= 5.

  DOWN-DIP DISTANCE.  The mesh's up-dip edge is the set of shallowest triangle centroids, one per
  0.1 deg latitude bin (the mesh stops at ~11 km depth, so this is a fixed up-dip reference line,
  NOT the trench -- it differs from the true trench by a latitude-dependent offset).
  ddist(detection) = horizontal distance in km from the epicentre to the nearest up-dip edge point.
  NOTE, measured here and important: depth and ddist are the same axis only WITHIN a latitude band.
  Pooled over the whole margin their Spearman is ~0, because the slab dip changes along strike (the
  up-dip edge runs from lon -124.4 at 40N to -126.8 at 49N while tremor depth goes the other way:
  ~48 km in the south, ~39 km in the north), so a detection can be far down-dip AND shallow.  The
  within-0.5-deg-band correlation is reported and is the meaningful one.

  LATITUDE DETRENDING (a confound S9 did not control).  A neighbourhood spans 220-720 km along
  strike (S4), so if the tremor simply MIGRATES ALONG STRIKE into a part of the margin where the
  interface is deeper, the daily median depth rises with no across-strike motion at all.  To
  separate the two, every detection also carries
      depth_resid = depth - median depth of ALL tremor in its own 0.25 deg latitude bin
      ddist_resid = ddist - median ddist of ALL tremor in its own 0.25 deg latitude bin
  (bin medians over the whole 2010-2026 catalogue, linearly interpolated in latitude).  The
  residual is "how far down-dip is this detection relative to the typical tremor at its own
  latitude", i.e. a genuine ACROSS-STRIKE coordinate.  The 10-bin phase profile is reported for
  the raw depth AND for depth_resid; only the residual version can be read as across-strike
  migration.

  PHASE (10 bins).  tau = (d - d_first) / (d_last - d_first) over the first and last DEFINED day;
  bin = min(int(10*tau), 9); bin value = median of the daily values in the bin; ANOMALY = bin value
  minus the event's median over all defined days.  Stack = median over events, 90 % CI by 2000
  bootstrap resamples of events.

  ONSET / PLATEAU / TAIL.  cs(d) = 3-day centred mean of n(d) over the defined days.  peak = the
  defined day maximising cs.  PLATEAU = the maximal run of consecutive defined days containing the
  peak with cs >= 0.5 * cs(peak).  ONSET = defined days strictly before the plateau; TAIL = defined
  days strictly after it.  Each segment needs >= 3 defined days to be used.
      d_onset_to_peak = median(zmed over PLATEAU) - median(zmed over ONSET)      (- = up-dip)
      d_peak_to_end   = median(zmed over TAIL)    - median(zmed over PLATEAU)    (+ = down-dip)
  NULL for both: the same 200 rigid >= 365 d shifts of the near-tremor dates used in S9, with the
  whole split recomputed on the shifted series; stacked over events per draw.

  START-DEFINITION TEST.  (a) TREMOR anchoring -- phase runs first defined day .. last defined day
  inside WA (this is S9's definition).  (b) CATALOG anchoring -- phase runs from the catalog slip
  start ts to the catalog slip end te, and only days in [ts, te] are used.  Both give a first->last
  fifth change; both get the shifted-window null.

  RADIUS TEST.  Everything in the start-definition test repeated at R = 40, 60, 100 km.  100 draws
  per cell (rather than 200) to keep the run short; the null is rebuilt per cell.

  PARTIAL SPEARMAN (spatial control).  Per event, over the defined days: rank-transform zmed, n and
  the daily mean along-strike position y; regress rank(zmed) on rank(y) and rank(n) on rank(y) by
  ordinary least squares; the partial rho is the Pearson correlation of the two residual vectors.
  Also reported: the raw Spearman(zmed, n) (S9's number), Spearman(zmed, y) (the spatial gradient
  being controlled for) and Spearman(n, y).

======================================================================================
OUTPUTS (all inside tremor_slip_geometry/)
======================================================================================
  results/s9b_phase10.csv        10 phase bins x {depth, latitude-detrended depth, down-dip
                                 distance, detrended down-dip distance}: observed stacked median,
                                 bootstrap 90 % CI, and (for the two depth series) the null band
  results/s9b_segments.csv       per event: peak day, segment lengths, median depth in
                                 onset/plateau/tail, d_onset_to_peak, d_peak_to_end and per-event
                                 null mean/sd/z; plus stacked rows (id = -1)
  results/s9b_downdip_stack.csv  day-aligned (days from event start AND days from first tremor)
                                 stacked median down-dip-distance anomaly, n events per lag
  results/s9b_sensitivity.csv    start definition x radius: stacked first->last change for depth and
                                 down-dip distance, with null mean/sd/p
  results/s9b_partial_rho.csv    per event: raw Spearman(z,n), Spearman(z,y), Spearman(n,y),
                                 partial Spearman(z,n | y), n days
  figures/s9b_depth_phase.png    6 panels (see the FIGURE section)

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s9b_depth_phase_reconcile.py
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
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"

RNG = np.random.default_rng(20260910)

PRE_A, POST_A = 10, 10
R_DEFAULT = 60.0
RADII = (40.0, 60.0, 100.0)
MIN_DET_DAY = 5
MIN_DEF_DAYS = 8
MIN_SEG_DAYS = 3            # onset / plateau / tail each need >= 3 defined days
PLATEAU_FRAC = 0.5          # plateau = smoothed count >= 0.5 x peak
N_BINS10 = 10
N_BINS5 = 5
N_SHIFT = 200
N_SHIFT_SENS = 100          # draws for the 6-cell sensitivity table
N_BOOT = 2000
LAG_LO, LAG_HI = -15, 70    # day-aligned stack range
COS47 = np.cos(np.deg2rad(47.0))

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


# ==================================================================================
# 1-3. Loaders (identical to s9_within_event.py)
# ==================================================================================
log("loading slip model ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
slabf = slab_depth_interpolator(fault)

NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PDEPTH = slabf(PLON, PLAT)
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
log(f"  {NT} days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles, "
    f"Slab2 depth {np.nanmin(PDEPTH):.1f}-{np.nanmax(PDEPTH):.1f} km")

# --- up-dip edge of the mesh: shallowest triangle per 0.1 deg latitude bin ---
lat_bins = np.arange(np.floor(PLAT.min() * 10) / 10, np.ceil(PLAT.max() * 10) / 10 + 0.1, 0.1)
edge_x, edge_y, edge_lat, edge_z = [], [], [], []
for lo in lat_bins:
    m = (PLAT >= lo) & (PLAT < lo + 0.1) & np.isfinite(PDEPTH)
    if m.sum():
        j = np.flatnonzero(m)[np.argmin(PDEPTH[m])]
        edge_x.append(PX[j]); edge_y.append(PY[j]); edge_lat.append(PLAT[j]); edge_z.append(PDEPTH[j])
EDGE = np.column_stack([edge_x, edge_y])
EDGE_TREE = cKDTree(EDGE)
log(f"  up-dip edge line: {len(EDGE)} points, lat {min(edge_lat):.2f}-{max(edge_lat):.2f}, "
    f"depth {np.min(edge_z):.1f}-{np.max(edge_z):.1f} km")

log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DEPTH = slabf(tr["lon"].values, tr["lat"].values)
TR_DDIST = EDGE_TREE.query(np.column_stack([TR_X, TR_Y]))[0]     # km from the up-dip edge
SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1
TR_OFF = (tr["t"].values.astype("datetime64[D]") - SHIFT_LO).astype(int)
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
mfin = np.isfinite(TR_DEPTH)
rr_dz = spearmanr(TR_DEPTH[mfin], TR_DDIST[mfin])
TR_LAT = tr["lat"].values
# within-band correlation: depth and ddist are the same axis only at fixed latitude
band_r = []
for lo in np.arange(40.0, 50.0, 0.5):
    m = mfin & (TR_LAT >= lo) & (TR_LAT < lo + 0.5)
    if m.sum() >= 500:
        band_r.append(spearmanr(TR_DEPTH[m], TR_DDIST[m]).statistic)
band_r = np.asarray(band_r)
log(f"  {len(tr)} detections; Slab2 depth finite {mfin.mean() * 100:.1f}%; "
    f"down-dip distance {TR_DDIST.min():.0f}-{TR_DDIST.max():.0f} km")
log(f"  Spearman(depth, down-dip distance): pooled {rr_dz.statistic:+.4f}  BUT within 0.5 deg "
    f"latitude bands median {np.median(band_r):+.4f} [{band_r.min():+.3f},{band_r.max():+.3f}] "
    f"over {band_r.size} bands -> the same axis locally, decoupled along strike")

# --- latitude detrending: subtract the margin-wide median at each detection's own latitude ---
LB = np.arange(39.5, 50.5 + 0.25, 0.25)
LBC, LBZ, LBD = [], [], []
for lo in LB[:-1]:
    m = (TR_LAT >= lo) & (TR_LAT < lo + 0.25)
    if m.sum() >= 50:
        LBC.append(lo + 0.125)
        LBZ.append(np.nanmedian(TR_DEPTH[m]) if np.isfinite(TR_DEPTH[m]).any() else np.nan)
        LBD.append(np.nanmedian(TR_DDIST[m]))
LBC, LBZ, LBD = np.asarray(LBC), np.asarray(LBZ), np.asarray(LBD)
# bins entirely outside the mesh hull have no median depth; they MUST be dropped before
# interpolating, otherwise np.interp propagates the NaN over a whole latitude range
okZ, okD = np.isfinite(LBZ), np.isfinite(LBD)
TR_DEPTH_R = TR_DEPTH - np.interp(TR_LAT, LBC[okZ], LBZ[okZ])
TR_DDIST_R = TR_DDIST - np.interp(TR_LAT, LBC[okD], LBD[okD])
log(f"  latitude detrending over {int(okZ.sum())}/{LBC.size} usable 0.25 deg bins "
    f"({int((~okZ).sum())} dropped, entirely outside the mesh hull): median depth "
    f"{np.nanmin(LBZ):.1f}-{np.nanmax(LBZ):.1f} km, median down-dip distance "
    f"{np.nanmin(LBD):.0f}-{np.nanmax(LBD):.0f} km; residual depth sd "
    f"{np.nanstd(TR_DEPTH_R):.2f} km (raw sd {np.nanstd(TR_DEPTH):.2f} km); "
    f"residual finite for {np.isfinite(TR_DEPTH_R).mean() * 100:.1f}% of detections")

ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
S9 = pd.read_csv(OUT_RES / "s9_depth_migration.csv")
log(f"  events {len(ev_in)} {ev_in.seg.value_counts().to_dict()}")


# ==================================================================================
# 4. Machinery
# ==================================================================================
def window_subset(sorted_off, s, w0_off, nW):
    a = (w0_off - s) % SPAN
    b = (a + nW - 1) % SPAN
    if a <= b:
        idx = np.arange(np.searchsorted(sorted_off, a, "left"),
                        np.searchsorted(sorted_off, b, "right"))
    else:
        idx = np.concatenate([np.arange(np.searchsorted(sorted_off, a, "left"), sorted_off.size),
                              np.arange(0, np.searchsorted(sorted_off, b, "right"))])
    return idx, (sorted_off[idx] + s) % SPAN - w0_off


def daily_stats(dayidx, z, zr, dd, ddr, y, nW):
    """Per day: count, median depth, median DETRENDED depth, median down-dip distance,
    median detrended down-dip distance, mean along-strike y."""
    cnt = np.bincount(dayidx, minlength=nW).astype(int) if dayidx.size else np.zeros(nW, int)
    out = [np.full(nW, np.nan) for _ in range(5)]
    if dayidx.size:
        o = np.argsort(dayidx, kind="stable")
        dsort = dayidx[o]
        cols = [z[o], zr[o], dd[o], ddr[o]]
        yy = y[o]
        bounds = np.searchsorted(dsort, np.arange(nW + 1))
        for i in range(nW):
            a, b = bounds[i], bounds[i + 1]
            if b - a >= MIN_DET_DAY:
                for c, col in enumerate(cols):
                    v = col[a:b]
                    v = v[np.isfinite(v)]
                    if v.size:
                        out[c][i] = np.median(v)
                out[4][i] = yy[a:b].mean()
    return (cnt, *out)


def binned_anom(series, defined, nbins):
    """Phase-binned anomaly of `series` over the boolean `defined` day mask."""
    dfn = np.flatnonzero(defined)
    if dfn.size < MIN_DEF_DAYS or dfn[-1] <= dfn[0]:
        return None
    v = series[dfn]
    ok = np.isfinite(v)
    if ok.sum() < MIN_DEF_DAYS:
        return None
    dd, v = dfn[ok], v[ok]
    tau = (dd - dfn[0]) / (dfn[-1] - dfn[0])
    b = np.minimum((tau * nbins).astype(int), nbins - 1)
    ref = float(np.median(v))
    out = np.full(nbins, np.nan)
    for k in range(nbins):
        m = b == k
        if m.any():
            out[k] = float(np.median(v[m])) - ref
    return out


def three_segments(cnt, zmed):
    """ONSET / PLATEAU / TAIL split at the smoothed-count peak.

    Returns dict with the peak position, segment sizes, segment median depths and the two changes.
    """
    dfn = np.flatnonzero(cnt >= MIN_DET_DAY)
    out = dict(n_def=int(dfn.size), peak_pos=np.nan, n_onset=0, n_plateau=0, n_tail=0,
               z_onset=np.nan, z_plateau=np.nan, z_tail=np.nan,
               d_onset_to_peak=np.nan, d_peak_to_end=np.nan)
    if dfn.size < MIN_DEF_DAYS:
        return out
    c = cnt[dfn].astype(float)
    cs = np.convolve(c, np.ones(3) / 3.0, mode="same")     # 3-day smoothing on the defined days
    ip = int(np.argmax(cs))
    thr = PLATEAU_FRAC * cs[ip]
    lo = ip
    while lo - 1 >= 0 and cs[lo - 1] >= thr:
        lo -= 1
    hi = ip
    while hi + 1 < cs.size and cs[hi + 1] >= thr:
        hi += 1
    onset = dfn[:lo]
    plateau = dfn[lo:hi + 1]
    tail = dfn[hi + 1:]
    out.update(peak_pos=float(ip / max(dfn.size - 1, 1)), n_onset=int(onset.size),
               n_plateau=int(plateau.size), n_tail=int(tail.size))

    def med(ix):
        if ix.size < MIN_SEG_DAYS:
            return np.nan
        v = zmed[ix]
        v = v[np.isfinite(v)]
        return float(np.median(v)) if v.size >= MIN_SEG_DAYS else np.nan

    out["z_onset"], out["z_plateau"], out["z_tail"] = med(onset), med(plateau), med(tail)
    out["d_onset_to_peak"] = out["z_plateau"] - out["z_onset"]
    out["d_peak_to_end"] = out["z_tail"] - out["z_plateau"]
    return out


def partial_spearman(z, n, y):
    """Spearman(z, n) with the along-strike position y regressed out of both (rank space)."""
    m = np.isfinite(z) & np.isfinite(n) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan
    rz, rn, ry = rankdata(z[m]), rankdata(n[m]), rankdata(y[m])
    if np.std(ry) < 1e-9:
        return np.nan
    ry1 = np.column_stack([np.ones(ry.size), ry])
    ez = rz - ry1 @ np.linalg.lstsq(ry1, rz, rcond=None)[0]
    en = rn - ry1 @ np.linalg.lstsq(ry1, rn, rcond=None)[0]
    if np.std(ez) < 1e-9 or np.std(en) < 1e-9:
        return np.nan
    return float(np.corrcoef(ez, en)[0, 1])


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan, 0
    return float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75)), int(v.size)


def fmt(v):
    m, a, b, n = med_iqr(v)
    return "     n/a       (n=0)" if n == 0 else f"{m:+8.3f} [{a:+7.3f},{b:+7.3f}] (n={n:2d})"


def stack_med(dd):
    """median over events of an (n_ev, nbins) matrix, plus a bootstrap 90 % CI."""
    M = np.asarray(dd, float)
    obs = np.nanmedian(M, axis=0)
    lo = np.full(M.shape[1], np.nan)
    hi = np.full(M.shape[1], np.nan)
    if M.shape[0] >= 5:
        bs = np.empty((N_BOOT, M.shape[1]))
        for i in range(N_BOOT):
            bs[i] = np.nanmedian(M[RNG.integers(0, M.shape[0], M.shape[0])], axis=0)
        lo, hi = np.nanpercentile(bs, 5, axis=0), np.nanpercentile(bs, 95, axis=0)
    return obs, lo, hi


# per-event cached geometry (footprint centroids) so the radius sweep is cheap
GEO = {}
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    GEO[pid] = dict(ts=ts, te=te, fx=PX[F], fy=PY[F], seg=str(r["seg"]), start=r["start"],
                    end=r["end"], Mw=float(r["Mw"]), dur=int(r["dur_d"]))


def near_sorted(pid, radius):
    """Detections within `radius` of any footprint centroid, sorted by day offset."""
    g = GEO[pid]
    h = TREMOR_TREE.query_ball_point(np.column_stack([g["fx"], g["fy"]]), r=radius)
    near = np.unique(np.concatenate([np.asarray(x, int) for x in h])) if len(h) else np.array([], int)
    near = near[np.argsort(TR_OFF[near], kind="stable")]
    return (near, TR_OFF[near], TR_DEPTH[near], TR_DEPTH_R[near], TR_DDIST[near],
            TR_DDIST_R[near], TR_Y[near])


# ==================================================================================
# 5. MAIN PASS at R = 60 km: 10-bin profiles, segments, partial rho, day-aligned stack
# ==================================================================================
log("main pass (R = 60 km) ...")
P10 = {"depth": [], "depth_resid": [], "ddist": [], "ddist_resid": []}   # 10-bin anomalies
NULL10 = {"depth": [], "depth_resid": []}                                # (n_ev, N_SHIFT, 10)
SEG_ROWS, PR_ROWS = [], []
NULL_SEG = []                              # (n_ev, N_SHIFT, 2)
ALIGN = {("depth_resid", "start"): [], ("depth_resid", "first"): [],
         ("ddist_resid", "start"): [], ("ddist_resid", "first"): []}
nlag = LAG_HI - LAG_LO + 1

for _, r in ev_in.iterrows():
    pid = int(r["id"])
    g = GEO[pid]
    ts, te = g["ts"], g["te"]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    w0_off = w0 - IDX_OF_SHIFT_LO
    near, soff, sz, szr, sd, sdr, sy = near_sorted(pid, R_DEFAULT)
    sub, dday = window_subset(soff, 0, w0_off, nW)
    cnt, zmed, zrmed, dmed, drmed, ymean = daily_stats(
        dday, sz[sub], szr[sub], sd[sub], sdr[sub], sy[sub], nW)
    dfn_mask = cnt >= MIN_DET_DAY

    SERIES = dict(depth=zmed, depth_resid=zrmed, ddist=dmed, ddist_resid=drmed)
    for nm, ser in SERIES.items():
        a = binned_anom(ser, dfn_mask, N_BINS10)
        if a is not None:
            P10[nm].append(a)

    seg = three_segments(cnt, zmed)
    seg.update(id=pid, start=g["start"], seg_name=g["seg"], Mw=g["Mw"], dur_d=g["dur"])

    # partial Spearman
    dfn = np.flatnonzero(dfn_mask)
    z_, n_, y_ = zmed[dfn], cnt[dfn].astype(float), ymean[dfn]
    m = np.isfinite(z_) & np.isfinite(y_)
    rho_zn = spearmanr(z_[m], n_[m]).statistic if m.sum() >= 8 else np.nan
    rho_zy = spearmanr(z_[m], y_[m]).statistic if m.sum() >= 8 else np.nan
    rho_ny = spearmanr(n_[m], y_[m]).statistic if m.sum() >= 8 else np.nan
    PR_ROWS.append(dict(id=pid, start=g["start"], seg=g["seg"], n_days=int(m.sum()),
                        rho_zn=rho_zn, rho_zy=rho_zy, rho_ny=rho_ny,
                        rho_zn_partial_y=partial_spearman(z_, n_, y_)))

    # day-aligned anomalies, two alignments, on the LATITUDE-DETRENDED series
    for nm in ("depth_resid", "ddist_resid"):
        ser = SERIES[nm]
        ref = np.nanmedian(ser[dfn_mask])
        for tag, anchor_i in (("start", ts - w0), ("first", dfn[0] if dfn.size else None)):
            a = np.full(nlag, np.nan)
            if anchor_i is not None:
                for k, lag in enumerate(range(LAG_LO, LAG_HI + 1)):
                    i = anchor_i + lag
                    if 0 <= i < nW and dfn_mask[i]:
                        a[k] = ser[i] - ref
            ALIGN[(nm, tag)].append(a)

    # ---- nulls: 200 shifts ----
    shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
    n10 = {k: np.full((N_SHIFT, N_BINS10), np.nan) for k in NULL10}
    nsg = np.full((N_SHIFT, 2), np.nan)
    for j, sh in enumerate(shifts):
        subn, ddn = window_subset(soff, int(sh), w0_off, nW)
        cn, zn, zrn, _, _, _ = daily_stats(
            ddn, sz[subn], szr[subn], sd[subn], sdr[subn], sy[subn], nW)
        for nm, ser in (("depth", zn), ("depth_resid", zrn)):
            a = binned_anom(ser, cn >= MIN_DET_DAY, N_BINS10)
            if a is not None:
                n10[nm][j] = a
        sg = three_segments(cn, zn)
        nsg[j] = (sg["d_onset_to_peak"], sg["d_peak_to_end"])
    for k in NULL10:
        NULL10[k].append(n10[k])
    NULL_SEG.append(nsg)
    v = nsg[:, 0][np.isfinite(nsg[:, 0])]
    w = nsg[:, 1][np.isfinite(nsg[:, 1])]
    seg["null_onset_mean"] = float(v.mean()) if v.size >= 20 else np.nan
    seg["null_onset_sd"] = float(v.std()) if v.size >= 20 else np.nan
    seg["null_tail_mean"] = float(w.mean()) if w.size >= 20 else np.nan
    seg["null_tail_sd"] = float(w.std()) if w.size >= 20 else np.nan
    seg["z_onset_null"] = ((seg["d_onset_to_peak"] - seg["null_onset_mean"]) / seg["null_onset_sd"]
                           if (v.size >= 20 and v.std() > 0) else np.nan)
    seg["z_tail_null"] = ((seg["d_peak_to_end"] - seg["null_tail_mean"]) / seg["null_tail_sd"]
                          if (w.size >= 20 and w.std() > 0) else np.nan)
    SEG_ROWS.append(seg)
    log(f"  {pid:5d} {g['start']} {g['seg']:7s} def={seg['n_def']:3d} peak@{seg['peak_pos']:.2f} "
        f"on/pl/tl={seg['n_onset']:3d}/{seg['n_plateau']:3d}/{seg['n_tail']:3d}  "
        f"d_on->pk={seg['d_onset_to_peak']:+6.2f} d_pk->end={seg['d_peak_to_end']:+6.2f}")

SEGDF = pd.DataFrame(SEG_ROWS)
PRDF = pd.DataFrame(PR_ROWS)
NULL10 = {k: np.stack(v) for k, v in NULL10.items()}
NULL_SEG = np.stack(NULL_SEG)
log("main pass done")


# ==================================================================================
# 6. PART 1 — 10-bin profile and the ONSET / PLATEAU / TAIL decomposition
# ==================================================================================
print("\n" + "=" * 100)
print("S9b PART 1 — FINER PHASE AXIS AND THE ONSET / PLATEAU / TAIL DECOMPOSITION")
print("=" * 100)

STK = {nm: stack_med(P10[nm]) for nm in P10}
NSTK = {}
for nm in NULL10:
    nev_j = np.isfinite(NULL10[nm][:, :, 0]).sum(axis=0)
    keep = nev_j >= 15
    nsm = np.nanmedian(NULL10[nm][:, keep, :], axis=0)
    NSTK[nm] = dict(nsm=nsm, med=np.nanmedian(nsm, axis=0),
                    lo=np.nanpercentile(nsm, 5, axis=0), hi=np.nanpercentile(nsm, 95, axis=0),
                    ndraw=int(keep.sum()))

rows10 = []
for nm, unit in (("depth", "km"), ("depth_resid", "km"), ("ddist", "km"), ("ddist_resid", "km")):
    obs, lo, hi = STK[nm]
    has_null = nm in NSTK
    print(f"\n  10-bin stacked anomaly — {nm}  (n_events {len(P10[nm])}"
          + (f"; null {NSTK[nm]['ndraw']}/{N_SHIFT} draws kept)" if has_null else "; no null)"))
    print("   bin  phase       obs " + unit + "     boot 90% CI       "
          + ("null med   null 5-95%" if has_null else ""))
    for k in range(N_BINS10):
        line = (f"   {k + 1:2d}  {k / 10:.1f}-{(k + 1) / 10:.1f}  {obs[k]:+8.3f}  "
                f"[{lo[k]:+7.3f},{hi[k]:+7.3f}]")
        if has_null:
            line += (f"   {NSTK[nm]['med'][k]:+8.3f}  [{NSTK[nm]['lo'][k]:+7.3f},"
                     f"{NSTK[nm]['hi'][k]:+7.3f}]")
        print(line)
        rows10.append(dict(series=nm, bin=k + 1, phase_lo=k / 10, phase_hi=(k + 1) / 10,
                           obs_median=obs[k], boot_lo90=lo[k], boot_hi90=hi[k],
                           null_median=NSTK[nm]["med"][k] if has_null else np.nan,
                           null_lo90=NSTK[nm]["lo"][k] if has_null else np.nan,
                           null_hi90=NSTK[nm]["hi"][k] if has_null else np.nan,
                           n_events=int(np.isfinite(np.array(P10[nm])[:, k]).sum())))
    kmin = int(np.nanargmin(obs))
    print(f"    SHAPE  shallowest bin {kmin + 1} ({obs[kmin]:+.3f})   "
          f"bin1->min {obs[kmin] - obs[0]:+.3f} (- = UP-DIP early)   "
          f"min->bin10 {obs[-1] - obs[kmin]:+.3f} (+ = DOWN-DIP late)   "
          f"net bin1->bin10 {obs[-1] - obs[0]:+.3f}")
    if has_null:
        nsm = NSTK[nm]["nsm"]
        for lab, o_ in (("bin1->min", obs[kmin] - obs[0]),
                        ("min->bin10", obs[-1] - obs[kmin]),
                        ("net", obs[-1] - obs[0])):
            if lab == "bin1->min":
                nd = np.nanmin(nsm, axis=1) - nsm[:, 0]
            elif lab == "min->bin10":
                nd = nsm[:, -1] - np.nanmin(nsm, axis=1)
            else:
                nd = nsm[:, -1] - nsm[:, 0]
            nd = nd[np.isfinite(nd)]
            pv = float((np.abs(nd) >= abs(o_)).mean()) if nd.size else np.nan
            print(f"      null for {lab:11s}: mean {nd.mean():+.3f} sd {nd.std():.3f} -> p = {pv:.4f}")

obsZ, loZ, hiZ = STK["depth"]
obsZR = STK["depth_resid"][0]
obsD = STK["ddist"][0]
nmedZ, nloZ, nhiZ = NSTK["depth"]["med"], NSTK["depth"]["lo"], NSTK["depth"]["hi"]
kmin = int(np.nanargmin(obsZ))

print("\n  ONSET / PLATEAU / TAIL (peak = smoothed daily-count maximum; plateau = >= 50 % of peak)")
print("    id      date        seg     n_on/pl/tl   z_on   z_pl   z_tl   on->pk   pk->end    z_on*  z_tl*")
for _, x in SEGDF.iterrows():
    print(f"    {int(x.id):5d}  {x.start}  {x.seg_name:7s} {int(x.n_onset):3d}/{int(x.n_plateau):3d}/"
          f"{int(x.n_tail):3d}  {x.z_onset:6.2f} {x.z_plateau:6.2f} {x.z_tail:6.2f}  "
          f"{x.d_onset_to_peak:+7.2f}  {x.d_peak_to_end:+7.2f}   {x.z_onset_null:+5.2f} {x.z_tail_null:+5.2f}")

seg_stack = {}
for lab, col, jj in (("onset->peak", "d_onset_to_peak", 0), ("peak->end", "d_peak_to_end", 1)):
    v = SEGDF[col].dropna()
    obs = float(np.median(v))
    nd = np.nanmedian(NULL_SEG[:, :, jj], axis=0)
    nd = nd[np.isfinite(nd)]
    p = float((np.abs(nd) >= abs(obs)).mean()) if nd.size else np.nan
    bs = np.array([np.median(v.values[RNG.integers(0, v.size, v.size)]) for _ in range(N_BOOT)])
    npos = int((v > 0).sum())
    print(f"\n  {lab:11s} stacked median {obs:+.3f} km  boot90% [{np.percentile(bs, 5):+.3f},"
          f"{np.percentile(bs, 95):+.3f}]  n={v.size}  {npos}/{v.size} positive "
          f"(sign p={binomtest(npos, v.size).pvalue:.4f})")
    print(f"              null median {np.median(nd):+.3f} sd {nd.std():.3f} ({nd.size} draws) -> p = {p:.4f}")
    print(f"              per-event z vs own null {fmt(SEGDF['z_onset_null' if jj == 0 else 'z_tail_null'])}")
    seg_stack[lab] = dict(obs=obs, p=p, n=int(v.size), npos=npos,
                          boot_lo=float(np.percentile(bs, 5)), boot_hi=float(np.percentile(bs, 95)),
                          null_med=float(np.median(nd)), null_sd=float(nd.std()))
print("\n  (- = UP-DIP / shallowing, + = DOWN-DIP / deepening)")


# ==================================================================================
# 7. PART 2 — Wech & Creager's statistic: down-dip distance vs DAYS FROM START
# ==================================================================================
print("\n" + "=" * 100)
print("S9b PART 2 — MEDIAN DAILY DOWN-DIP DISTANCE vs DAYS FROM EVENT START (their plot)")
print("=" * 100)
print(f"  down-dip distance = km from the mesh up-dip edge (~{np.min(edge_z):.0f} km depth line). "
      f"It equals depth as an axis only WITHIN a latitude band (median Spearman "
      f"{np.median(band_r):+.3f} over {band_r.size} 0.5 deg bands); pooled over the margin the two "
      f"decouple ({rr_dz.statistic:+.3f}) because the slab dip changes along strike, which is why "
      f"everything below is latitude-detrended.")
print("  NOTE: a fixed DAY axis mixes the early phase of an 85-day event with the tail of a 2-day "
      "one (our durations span 2-85 d), so this un-normalised view is much noisier than the "
      "phase-normalised one of Part 1; late lags rest on very few events (see n_ev).")
ALG = {k: np.array(v, float) for k, v in ALIGN.items()}
lags = np.arange(LAG_LO, LAG_HI + 1)
dn_rows = []
print("\n  All four series are LATITUDE-DETRENDED (across-strike only).")
print("   lag(d)   DOWN-DIP DIST km          DEPTH km")
print("            from slip  from 1st       from slip  from 1st      n_ev(slip/1st)")
for k, lg in enumerate(lags):
    v = {key: np.nanmedian(ALG[key][:, k]) for key in ALG}
    ns = int(np.isfinite(ALG[("ddist_resid", "start")][:, k]).sum())
    nf = int(np.isfinite(ALG[("ddist_resid", "first")][:, k]).sum())
    dn_rows.append(dict(lag_days=int(lg),
                        ddist_resid_from_slip_start=v[("ddist_resid", "start")],
                        ddist_resid_from_first_tremor=v[("ddist_resid", "first")],
                        depth_resid_from_slip_start=v[("depth_resid", "start")],
                        depth_resid_from_first_tremor=v[("depth_resid", "first")],
                        n_ev_slip_start=ns, n_ev_first_tremor=nf))
    if lg % 5 == 0 and (ns >= 5 or nf >= 5):
        print(f"   {lg:+5d}   {v[('ddist_resid', 'start')]:+9.2f} {v[('ddist_resid', 'first')]:+9.2f}"
              f"      {v[('depth_resid', 'start')]:+9.2f} {v[('depth_resid', 'first')]:+9.2f}"
              f"        {ns:3d}/{nf:3d}")


def early_late(A, lo1, hi1, lo2, hi2):
    i1 = (lags >= lo1) & (lags <= hi1)
    i2 = (lags >= lo2) & (lags <= hi2)
    a = np.nanmedian(A[:, i1], axis=1)
    b = np.nanmedian(A[:, i2], axis=1)
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.median(b[m] - a[m])), int(m.sum()), float(np.median(a[m])), float(np.median(b[m]))


for key in (("ddist_resid", "start"), ("ddist_resid", "first"),
            ("depth_resid", "start"), ("depth_resid", "first")):
    A = ALG[key]
    d1, n1, a1, b1 = early_late(A, 0, 9, 10, 24)
    d2, n2, a2, b2 = early_late(A, 10, 24, 25, 60)
    anc = "catalog slip start" if key[1] == "start" else "first tremor day"
    print(f"\n  {key[0]:12s} aligned on {anc:19s}: days 0-9 {a1:+.2f} -> 10-24 {b1:+.2f} "
          f"= {d1:+.2f} km (n={n1})")
    print(f"  {'':12s} {'':30s}  days 10-24 {a2:+.2f} -> 25-60 {b2:+.2f} = {d2:+.2f} km (n={n2})")

wc = ev_in[(ev_in.seg == "north") & (ev_in.start < "2012-01-01")]
print(f"\n  Events of ours that can overlap Wech & Creager's 2010-2011 Washington ETS "
      f"(PNSN catalogue starts {tr['t'].min().date()}): {len(wc)}")
for _, x in wc.iterrows():
    s = SEGDF[SEGDF.id == int(x.id)].iloc[0]
    print(f"    {int(x.id):5d} {x.start}..{x['end']} Mw{float(x.Mw):.2f}  "
          f"onset->peak {s.d_onset_to_peak:+.2f} km, peak->end {s.d_peak_to_end:+.2f} km, "
          f"net 1st->5th fifth {float(S9[S9.id == int(x.id)].dz_tremor_fifth_km.iloc[0]):+.2f} km")
print("  CAVEAT: their exact four-episode list could not be verified from this environment "
      "(no literature access); the two above are simply every northern confirmed event of ours "
      "before 2012, and their ETS years before 2010 are outside the PNSN file entirely.")


# ==================================================================================
# 8. PART 3 — does the SIGN depend on the start definition or the radius?
# ==================================================================================
print("\n" + "=" * 100)
print("S9b PART 3 — SENSITIVITY OF THE SIGN TO THE START DEFINITION AND THE RADIUS")
print("=" * 100)
print("  'tremor' anchoring = phase between the first and last >= 5-detection day in "
      "[start-10, end+10] (S9's rule).  'catalog' anchoring = phase between the catalog slip "
      "start and end, days outside DUR discarded.")
sens_rows = []
print("\n    anchor    R(km)  n_ev   DEPTH km  DEPTH_detr   null mean+-sd    p     DOWNDIP km")
for radius in RADII:
    NEAR = {int(r["id"]): near_sorted(int(r["id"]), radius) for _, r in ev_in.iterrows()}
    for anchor in ("tremor", "catalog"):
        # NOTE: the statistic is EXACTLY S9's -- stack the per-event fifth profiles, then take
        # bin5 - bin1 of the stacked medians (NOT the median of per-event differences), so the
        # (tremor, 60 km) cell must reproduce S9's +2.26 km / p = 0.020.
        Pz, Pzr, Pd, nullz = [], [], [], []
        for _, r in ev_in.iterrows():
            pid = int(r["id"])
            g = GEO[pid]
            ts, te = g["ts"], g["te"]
            if anchor == "tremor":
                w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
            else:
                w0, w1 = ts, te
            nW = w1 - w0 + 1
            w0_off = w0 - IDX_OF_SHIFT_LO
            near, soff, sz, szr, sd, sdr, sy = NEAR[pid]
            sub, dday = window_subset(soff, 0, w0_off, nW)
            cnt, zmed, zrmed, dmed, drmed, _ = daily_stats(
                dday, sz[sub], szr[sub], sd[sub], sdr[sub], sy[sub], nW)
            dm = cnt >= MIN_DET_DAY
            for lst, ser in ((Pz, zmed), (Pzr, zrmed), (Pd, dmed)):
                a = binned_anom(ser, dm, N_BINS5)
                if a is not None:
                    lst.append(a)
            per = np.full((N_SHIFT_SENS, N_BINS5), np.nan)
            for j, sh in enumerate(RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT_SENS)):
                subn, ddn = window_subset(soff, int(sh), w0_off, nW)
                cn, zn, _, _, _, _ = daily_stats(
                    ddn, sz[subn], szr[subn], sd[subn], sdr[subn], sy[subn], nW)
                a = binned_anom(zn, cn >= MIN_DET_DAY, N_BINS5)
                if a is not None:
                    per[j] = a
            nullz.append(per)
        nullz = np.stack(nullz)                       # (n_ev, N_SHIFT_SENS, 5)
        keep = np.isfinite(nullz[:, :, 0]).sum(axis=0) >= 15
        nsm = np.nanmedian(nullz[:, keep, :], axis=0)
        nd = nsm[:, 4] - nsm[:, 0]
        nd = nd[np.isfinite(nd)]
        oz = float(np.nanmedian(Pz, axis=0)[4] - np.nanmedian(Pz, axis=0)[0]) if Pz else np.nan
        ozr = float(np.nanmedian(Pzr, axis=0)[4] - np.nanmedian(Pzr, axis=0)[0]) if Pzr else np.nan
        od = float(np.nanmedian(Pd, axis=0)[4] - np.nanmedian(Pd, axis=0)[0]) if Pd else np.nan
        pv = float((np.abs(nd) >= abs(oz)).mean()) if nd.size else np.nan
        sens_rows.append(dict(anchor=anchor, radius_km=radius, n_events=len(Pz),
                              dz_depth_km=oz, dz_depth_resid_km=ozr, dz_ddist_km=od,
                              null_mean=float(np.mean(nd)) if nd.size else np.nan,
                              null_sd=float(np.std(nd)) if nd.size else np.nan, p_null=pv,
                              n_null_draws=int(nd.size)))
        print(f"    {anchor:8s} {radius:5.0f}  {len(Pz):4d}     {oz:+8.3f}   {ozr:+8.3f}   "
              f"{np.mean(nd) if nd.size else np.nan:+6.3f}+-{np.std(nd) if nd.size else np.nan:5.3f}"
              f"  {pv:.3f}   {od:+8.2f}")
SENS = pd.DataFrame(sens_rows)
print(f"\n  SIGN: positive (deepening) in {int((SENS.dz_depth_km > 0).sum())}/{len(SENS)} of the "
      f"6 anchor x radius cells; range {SENS.dz_depth_km.min():+.2f} to {SENS.dz_depth_km.max():+.2f} km."
      f"  Latitude-detrended: positive in {int((SENS.dz_depth_resid_km > 0).sum())}/{len(SENS)}, "
      f"range {SENS.dz_depth_resid_km.min():+.2f} to {SENS.dz_depth_resid_km.max():+.2f} km.")
_chk = SENS[(SENS.anchor == "tremor") & (SENS.radius_km == 60.0)].iloc[0]
print(f"  CROSS-CHECK vs S9 (same statistic, same rule): this run {_chk.dz_depth_km:+.3f} km "
      f"p={_chk.p_null:.3f}  |  S9 reported +2.256 km p=0.020 "
      f"(difference is only the {N_SHIFT_SENS}- vs 200-draw null).")


# ==================================================================================
# 9. PART 4 — spatial control on "shallower when most intense"
# ==================================================================================
print("\n" + "=" * 100)
print("S9b PART 4 — IS 'SHALLOWER WHEN MOST INTENSE' TEMPORAL, OR THE ALONG-STRIKE GRADIENT?")
print("=" * 100)
print("    id      date        seg     ndays  rho(z,n)  rho(z,y)  rho(n,y)  rho(z,n | y)")
for _, x in PRDF.iterrows():
    print(f"    {int(x.id):5d}  {x.start}  {x.seg:7s} {int(x.n_days):4d}   {x.rho_zn:+7.3f}  "
          f"{x.rho_zy:+7.3f}  {x.rho_ny:+7.3f}     {x.rho_zn_partial_y:+7.3f}")
raw, par = PRDF.rho_zn.dropna(), PRDF.rho_zn_partial_y.dropna()
print(f"\n  raw     Spearman(depth, count)      {fmt(raw)}  {int((raw < 0).sum())}/{raw.size} negative"
      f"  sign p={binomtest(int((raw < 0).sum()), raw.size).pvalue:.4f}  "
      f"Wilcoxon p={wilcoxon(raw).pvalue:.4f}")
print(f"  partial Spearman(depth, count | y)  {fmt(par)}  {int((par < 0).sum())}/{par.size} negative"
      f"  sign p={binomtest(int((par < 0).sum()), par.size).pvalue:.4f}  "
      f"Wilcoxon p={wilcoxon(par).pvalue:.4f}")
print(f"  the gradient being controlled for: Spearman(depth, along-strike y) {fmt(PRDF.rho_zy)}")
print(f"                                     Spearman(count, along-strike y) {fmt(PRDF.rho_ny)}")
d = (PRDF.rho_zn_partial_y - PRDF.rho_zn).dropna()
print(f"  partial - raw  {fmt(d)}  Wilcoxon p={wilcoxon(d).pvalue:.4f}")


# ==================================================================================
# 10. WRITE
# ==================================================================================
pd.DataFrame(rows10).to_csv(OUT_RES / "s9b_phase10.csv", index=False)
stack_extra = pd.DataFrame([
    dict(id=-1, start="STACK", seg_name=k, d_onset_to_peak=v["obs"] if k == "onset->peak" else np.nan,
         d_peak_to_end=v["obs"] if k == "peak->end" else np.nan, null_onset_mean=v["null_med"],
         null_onset_sd=v["null_sd"], z_onset_null=np.nan, z_tail_null=np.nan,
         n_def=v["n"], p_null=v["p"], boot_lo=v["boot_lo"], boot_hi=v["boot_hi"], n_pos=v["npos"])
    for k, v in seg_stack.items()])
pd.concat([SEGDF, stack_extra], ignore_index=True).to_csv(OUT_RES / "s9b_segments.csv", index=False)
pd.DataFrame(dn_rows).to_csv(OUT_RES / "s9b_downdip_stack.csv", index=False)
SENS.to_csv(OUT_RES / "s9b_sensitivity.csv", index=False)
PRDF.to_csv(OUT_RES / "s9b_partial_rho.csv", index=False)
log("wrote results/s9b_{phase10,segments,downdip_stack,sensitivity,partial_rho}.csv")


# ==================================================================================
# 11. FIGURE
# ==================================================================================
fig, axes = plt.subplots(2, 3, figsize=(19, 10))
xb = np.arange(1, N_BINS10 + 1)

ax = axes[0, 0]
ax.fill_between(xb, nloZ, nhiZ, color="0.75", alpha=0.6, label="null 5-95%")
ax.fill_between(xb, loZ, hiZ, color="tab:blue", alpha=0.25, label="boot 90% CI")
ax.plot(xb, obsZ, "o-", color="tab:blue", lw=2.4, label=f"tremor depth (n={len(P10['depth'])})")
ax.axhline(0, color="k", lw=0.7)
ax.axvline(kmin + 1, color="tab:green", ls=":", lw=1.5)
ax.annotate(f"shallowest\nbin {kmin + 1}", (kmin + 1, obsZ[kmin]), textcoords="offset points",
            xytext=(6, -26), color="tab:green", fontsize=9)
ax.set_xlabel("event phase (tenths)")
ax.set_ylabel("depth anomaly (km)   + = deeper / down-dip")
ax.set_title(f"1. 10-bin phase profile\nUP-DIP {obsZ[kmin] - obsZ[0]:+.2f} km then DOWN-DIP "
             f"{obsZ[-1] - obsZ[kmin]:+.2f} km (net {obsZ[-1] - obsZ[0]:+.2f})")
ax.legend(fontsize=8)
ax.set_xticks(xb)

ax = axes[0, 1]
oR, lR, hR = STK["depth_resid"]
ax.fill_between(xb, NSTK["depth_resid"]["lo"], NSTK["depth_resid"]["hi"], color="0.75",
                alpha=0.6, label="null 5-95%")
ax.fill_between(xb, lR, hR, color="tab:purple", alpha=0.25, label="boot 90% CI")
ax.plot(xb, oR, "s-", color="tab:purple", lw=2.4, label="latitude-detrended depth")
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("event phase (tenths)")
ax.set_ylabel("across-strike depth anomaly (km)")
kR = int(np.nanargmin(oR))
ax.set_title(f"1. ACROSS-STRIKE only (latitude detrended)\nUP-DIP {oR[kR] - oR[0]:+.2f} then "
             f"DOWN-DIP {oR[-1] - oR[kR]:+.2f} (net {oR[-1] - oR[0]:+.2f})")
ax.legend(fontsize=8)
ax.set_xticks(xb)

ax = axes[0, 2]
d1 = SEGDF.d_onset_to_peak.dropna()
d2 = SEGDF.d_peak_to_end.dropna()
bp = ax.boxplot([d1.values, d2.values], tick_labels=[f"onset$\\to$peak\nn={d1.size}",
                                                     f"peak$\\to$end\nn={d2.size}"],
                showmeans=True, patch_artist=True)
for pt, c in zip(bp["boxes"], ("tab:green", "tab:red")):
    pt.set_facecolor(c)
    pt.set_alpha(0.5)
for i, v in enumerate((d1.values, d2.values)):
    ax.plot(np.full(v.size, i + 1) + RNG.normal(0, 0.05, v.size), v, "k.", ms=6, alpha=0.7)
ax.axhline(0, color="k", lw=1)
ax.set_ylabel("depth change (km)")
ax.set_title(f"1. onset$\\to$peak {seg_stack['onset->peak']['obs']:+.2f} km "
             f"(p={seg_stack['onset->peak']['p']:.3f})\n"
             f"peak$\\to$end {seg_stack['peak->end']['obs']:+.2f} km "
             f"(p={seg_stack['peak->end']['p']:.3f})")

ax = axes[1, 0]
for key, c, nm in ((("ddist_resid", "start"), "tab:blue", "down-dip dist, from slip start"),
                   (("ddist_resid", "first"), "tab:orange", "down-dip dist, from 1st tremor"),
                   (("depth_resid", "first"), "tab:green", "depth (x5), from 1st tremor")):
    A = ALG[key]
    med = np.nanmedian(A, axis=0)
    nev = np.isfinite(A).sum(axis=0)
    ok = nev >= 8
    ax.plot(lags[ok], med[ok] * (5 if key[0] == "depth_resid" else 1), "-", color=c, lw=2, label=nm)
ax.axvline(0, color="k", lw=1)
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("days from event start")
ax.set_ylabel("latitude-detrended anomaly (km)")
ax.set_title("2. Wech & Creager's plot on our 27 events\n(no phase normalisation, across-strike only)")
ax.legend(fontsize=8)

ax = axes[1, 1]
for i, (anchor, mk) in enumerate((("tremor", "o"), ("catalog", "s"))):
    sub = SENS[SENS.anchor == anchor]
    ax.errorbar(sub.radius_km + i * 2 - 1, sub.dz_depth_km, yerr=sub.null_sd, fmt=mk, ms=9,
                capsize=4, label=f"{anchor} anchoring (err = null sd)")
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("neighbourhood radius (km)")
ax.set_ylabel("stacked 1st$\\to$last fifth depth change (km)")
ax.set_title(f"3. sign stability\npositive in {int((SENS.dz_depth_km > 0).sum())}/6 cells")
ax.legend(fontsize=8)
ax.set_xticks(list(RADII))

ax = axes[1, 2]
ax.plot(PRDF.rho_zn, PRDF.rho_zn_partial_y, "o", color="tab:brown")
ax.plot([-1, 1], [-1, 1], "k--", lw=1)
ax.axhline(0, color="k", lw=0.7)
ax.axvline(0, color="k", lw=0.7)
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_xlabel("raw Spearman(daily depth, daily count)")
ax.set_ylabel("partial, along-strike position removed")
ax.set_title(f"4. spatial control\nraw med {np.median(raw):+.2f} -> partial {np.median(par):+.2f}"
             f"  ({int((par < 0).sum())}/{par.size} still negative)")

fig.suptitle("S9b — reconciling the within-event depth migration with Wech & Creager (2011)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_FIG / "s9b_depth_phase.png", dpi=130)
plt.close(fig)
log("wrote figures/s9b_depth_phase.png")
log("S9b done.")
