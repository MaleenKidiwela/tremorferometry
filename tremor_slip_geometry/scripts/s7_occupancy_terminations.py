#!/usr/bin/env python
"""
S7 — occupancy classification (Q5) and along-strike terminations (Q17), 2010-2026.

PART A asks a *spatial* question that S3/S3b could not answer: is there a persistent set of
triangles that slip in the model but never tremor (silent slip), and a converse set that tremor
without modeled slip?  S3b showed that any comparison built on tremor COUNTS (and even on raw
tremor-DAYS) is dominated by the growth of the PNSN network -- the tremor-day rate grew x1.30 in
the north, x1.92 in the centre and x2.06 in the south between 2010-14 and 2020-26, while the
modeled slip rate was flat everywhere.  Absolute occupancy levels are therefore not comparable
across epochs.  The fix used here is to compare only RANKS: within each epoch and each segment,
every triangle's slip and tremor occupancy is expressed as a percentile rank among the triangles
of the same segment in that epoch.  A margin-wide or segment-wide change in detection sensitivity
is a monotone transformation of the occupancy of every triangle in that segment and therefore
leaves the ranks unchanged; only a change in the SPATIAL PATTERN can move a triangle's class.
(The residual it cannot remove is a spatially heterogeneous sensitivity change WITHIN a segment --
see the caveats printed at the end.)  The D3_nb tremor-day variant is used because S3b found it
the least sensitivity-dependent of the four (hardening the per-day detection threshold from >=1 to
>=3 shrank the spurious central/north contrast, and the neighbourhood removes the epicentre-
scatter channel).

PART B asks whether the along-strike ENDS of episodes recur at the same latitudes (barriers), and
whether the modeled slip stops where the tremor stops.

--------------------------------------------------------------------------------------
INPUTS (all read-only, all outside this folder)
--------------------------------------------------------------------------------------
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz    label (7169 x 742, int; 0=none,
                                                                     else patch id), days (7169,
                                                                     'YYYY-MM-DD'), patch_lat/lon/
                                                                     depth_km/area_km2 (742)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.csv    slip-patch catalog
                                                                     (substantive = ~minor & ~broad_pulse)
  fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv  per-patch tremor confirmation
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  'A_raw' (7169 x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                             PNSN tremor (time, lat, lon;
                                                                     the depth column is a placeholder
                                                                     and is IGNORED)
  /home/jovyan/cascadia-sse-strain (mesh + sse_strain)               Ld (742 x 13) dip-slip patterns
                                                                     [m per unit amplitude]; Slab2
                                                                     depth interpolator f(lon, lat)
  tremor_slip_geometry/results/s4_fronts_events.csv                  column v_t_km_d = Theil-Sen speed
                                                                     of the daily tremor centroid along
                                                                     strike (+ = northward); used for the
                                                                     propagation direction in Part B

  Thrust slip per triangle, cumulative, metres:  raw_cum = -(Ld @ A_raw.T)   (742 x T)
  Daily positive thrust slip rate (mm/d):        pos_rate = clip(movavg7(diff(1e3*raw_cum)), 0, None)
                                                 on the day grid days[1:]

--------------------------------------------------------------------------------------
EXACT DEFINITIONS
--------------------------------------------------------------------------------------
  Window        2010-01-06 .. 2026-05-27 (the span of tremor days on the mesh).
  Epochs        E1 2010-01-06..2014-12-31 (4.99 yr), E2 2015-01-01..2019-12-31 (5.00 yr),
                E3 2020-01-01..2026-05-27 (6.40 yr).  Both tremor and slip are clipped to these.
  Segments      north lat >= 46, central 44 <= lat < 46, south lat < 44 (triangle centroid).
  km frame      x = (lon + 124) * 111 * cos(47 deg), y = (lat - 47) * 111.

PART A
  Tremor->triangle   nearest triangle centroid in the km frame, dropped beyond 40 km.
  c[day, tri]        detections assigned to that triangle on that UTC calendar day.
  Neighbourhood      nb(tri) = {tri} u {triangles sharing >= 1 mesh VERTEX with tri} (0.25-deg grid,
                     rounded to 1e-6 deg).  c_nb[day, tri] = sum over nb(tri).
  D  (tremor-days/yr) # days in the epoch with c_nb[day, tri] >= 3, divided by the epoch length in
                     years.  This is S3b's D3_nb variant, the least sensitivity-dependent one.
  S  (mm/yr)         sum of pos_rate over the epoch days, divided by the epoch length in years.
  RANKS              rS, rD = percentile rank (scipy rankdata 'average', /n*100) of S and D among
                     the triangles OF THE SAME SEGMENT in THAT epoch.  Ties (many triangles have
                     D = 0) share the average rank; the tie fractions are printed.
  CLASSES per epoch  SLIP+TREMOR  rS >= 50 and rD >= 50
                     SLIP-ONLY    rS >= 50 and rD < 25
                     TREMOR-ONLY  rD >= 50 and rS < 25
                     QUIET        rS < 25 and rD < 25
                     MIXED        anything else (mutually exclusive by construction)
  Reported           per segment x epoch: fraction of triangles and fraction of AREA in each class;
                     mean/median Slab2 depth of each class (the along-dip structure).
  PERSISTENCE        fraction of triangles whose class is identical in E1, E2 and E3, per segment,
                     against a null in which the class labels are permuted independently among the
                     triangles of that segment in each epoch (5000 draws; this preserves each
                     epoch's class composition exactly and destroys only the spatial pairing).
                     Reported overall and per class; p = (1 + #{null >= obs}) / (ndraw + 1).
  ABSOLUTE MAP       quiet-slip candidate: S_tot >= 20 mm accumulated over the whole window AND
                     D3_nb_tot < 30 tremor-days in 16.4 yr.  Converse (tremor-without-slip):
                     D3_nb_tot >= 300 tremor-days AND S_tot < 5 mm.  THESE LITERAL SLIP THRESHOLDS
                     ARE DEGENERATE: because only the POSITIVE daily rate is accumulated, model
                     noise is rectified upward and every one of the 742 triangles accumulates
                     >= 123 mm over 16.4 yr (median 407 mm), so 742/742 pass '>= 20 mm' and 0/742
                     pass '< 5 mm'.  They are reported for completeness, and each is accompanied by
                     a PERCENTILE version that carries the intended meaning:
                     quiet-slip_p75    S_tot >= the 75th percentile of S_tot within its SEGMENT and
                                       D3_nb_tot < 30 tremor-days;
                     tremor-no-slip_p25  D3_nb_tot >= 300 and S_tot <= the 25th percentile of S_tot
                                       within its segment.
                     (D_own totals are reported alongside as a sensitivity check.)  For each
                     candidate the Slab2 depth is
                     compared with its SEGMENT's tremor depth p10/p90, where the tremor depth
                     distribution is Slab2 depth at the epicentre of every detection assigned to a
                     triangle of that segment: up-dip of band (depth < p10), in band, down-dip of
                     band (depth > p90).

PART B  (27 CONFIRMED patches: tremor_vs_slip_patches catalog=='thrust' & confirmed & seg in
        {north, central} & Mw >= 6.0 & ~broad_pulse -- the same set as S1/S4)
  Footprint F        triangles labelled with the patch id on any day in [ts, te].
  TREMOR extent      detections within 60 km (km frame, epicentral) of at least one F centroid with
                     date in [start, end];  south end = 5th percentile of their LATITUDE,
                     north end = 95th percentile.
  SLIP extent        triangles with accumulated slip acc = 1e3*(cum[:, te] - cum[:, ts-30]) clipped
                     at >= 0 of at least 3 mm, restricted to latitudes inside F's latitude span
                     +- 1 deg;  south/north end = min/max latitude of those triangles.
  Differences        (tremor end - slip end) in km = 111 * (lat difference), separately at the north
                     and the south termination.
  Direction          sign of v_t_km_d from results/s4_fronts_events.csv (+ = northward).  Events with
                     |v_t| <= 0.2 km/d (or missing) are 'ambiguous' and their direction is
                     re-derived as the sign of the Theil-Sen slope of the daily mean latitude of the
                     same 60-km tremor set (days with >= 5 detections) against day number.
  BARRIER TEST       all 2 x 27 tremor termination latitudes are histogrammed in 0.25-deg bins.
                     NULL: each event keeps its own extent LENGTH but its position is drawn
                     uniformly at random inside its segment's domain (5000 draws).  Two domains are
                     reported: 'strict' = the nominal segment bounds (north 46 .. mesh max 49.92,
                     central 44 .. 46) and 'domain' = those bounds widened to contain every observed
                     extent of that segment (the primary test; the strict domain gives a long event
                     almost no freedom, so it mostly reproduces the observation and has no power).
                     p per bin = (1 + #{null count >= observed}) / (ndraw + 1), Bonferroni over the
                     number of TESTABLE bins.  A bin the null can NEVER reach (outside every
                     segment domain, so its null count is identically 0) would otherwise return the
                     minimum possible p for a single observed termination, which is an artefact of
                     the domain, not a barrier; such bins are excluded from the test and from the
                     Bonferroni count and are flagged 'unreachable'.
  PERSISTENCE        for each barrier bin, the observed and null-expected counts of terminations
                     within each epoch (events assigned by start date), and the per-epoch p.
  SLIP COMPARISON    median |tremor end - slip end| in km at the north and south ends, and whether
                     the slip terminations concentrate in the SAME 0.25-deg bins as the tremor ones.
  DIRECTION TEST     Mann-Whitney U on north-end (and south-end) latitudes of northward- vs
                     southward-propagating events, and the fraction of each group's terminations
                     landing in the identified barrier bins.

--------------------------------------------------------------------------------------
OUTPUTS (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------
  results/s7_occupancy_triangles.csv     742 rows: tri, lat, lon, depth_km, area_km2, seg, per epoch
                                         S_mm_yr, D_days_yr, rank_S, rank_D, class; class_persistent
                                         ('varies' when not identical in all three); window totals
  results/s7_occupancy_summary.csv       per (segment, epoch, class): n_tri, frac_tri, frac_area,
                                         mean/median depth; plus persistence rows (obs vs null)
  results/s7_quiet_slip_candidates.csv   candidate quiet-slip and tremor-without-slip triangles with
                                         lat/lon/depth, totals, and their position relative to the
                                         segment tremor depth band
  results/s7_terminations_events.csv     27 rows: id, start, end, seg, Mw, tremor/slip north+south
                                         termination latitudes, differences (km), extents, direction
  results/s7_termination_bins.csv        per 0.25-deg bin: observed tremor terminations, null mean,
                                         p, Bonferroni p, per-epoch counts, slip-termination count
  figures/s7_occupancy_maps.png          3 epoch class maps + persistent-class map + class-depth
                                         panel (+ absolute quiet-slip / tremor-only map)
  figures/s7_terminations.png            termination histogram vs null (all + by epoch); tremor vs
                                         slip termination scatter; terminations by propagation
                                         direction

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s7_occupancy_terminations.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, rankdata, spearmanr, theilslopes

ROOT = Path("/home/jovyan/tremorferometry")
OUT = ROOT / "tremor_slip_geometry"
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
STRAIN = Path("/home/jovyan/cascadia-sse-strain")

SEED = 20260910
RNG = np.random.default_rng(SEED)
COS47 = np.cos(np.deg2rad(47.0))

WIN_LO, WIN_HI = "2010-01-06", "2026-05-27"
EPOCHS = [("E1", "2010-01-06", "2014-12-31"),
          ("E2", "2015-01-01", "2019-12-31"),
          ("E3", "2020-01-01", "2026-05-27")]
MAX_ASSIGN_KM = 40.0
DAY_DET_THR = 3          # >= 3 detections in the neighbourhood -> a tremor-day (D3_nb)
HI_RANK, LO_RANK = 50.0, 25.0
N_PERM_CLASS = 5000
# absolute-map thresholds
QS_SLIP_MM, QS_DAYS = 20.0, 30
TO_DAYS, TO_SLIP_MM = 300, 5.0
# Part B
EVENT_TREMOR_RADIUS_KM = 60.0
EVENT_EDGE_ACC_MM = 3.0
EVENT_PRE_DAYS = 30
EVENT_LAT_PAD_DEG = 1.0
BIN_W = 0.25
N_PERM_TERM = 5000
DIR_EPS = 0.2            # |v_t| <= this km/d -> ambiguous, re-derive

CLASSES = ["SLIP+TREMOR", "SLIP-ONLY", "TREMOR-ONLY", "QUIET", "MIXED"]
CLASS_COL = {"SLIP+TREMOR": "#b2182b", "SLIP-ONLY": "#2166ac", "TREMOR-ONLY": "#f4a582",
             "QUIET": "0.80", "MIXED": "#4d9221", "varies": "0.93"}

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def seg_of_lat(lat):
    lat = np.asarray(lat, dtype=float)
    out = np.full(lat.shape, "south", dtype=object)
    out[lat >= 46.0] = "north"
    out[(lat >= 44.0) & (lat < 46.0)] = "central"
    return out


def moving_average_centred(a, w=7, axis=-1):
    a = np.moveaxis(np.asarray(a, dtype=float), axis, -1)
    n = a.shape[-1]
    csum = np.concatenate([np.zeros(a.shape[:-1] + (1,)), np.cumsum(a, axis=-1)], axis=-1)
    half = w // 2
    idx = np.arange(n)
    lo, hi = np.clip(idx - half, 0, n), np.clip(idx + half + 1, 0, n)
    return np.moveaxis((csum[..., hi] - csum[..., lo]) / (hi - lo), -1, axis)


def pct_rank(x):
    """Percentile rank 0-100 with average ranks for ties."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    return rankdata(x, method="average") / len(x) * 100.0


def classify(rS, rD):
    out = np.full(rS.shape, "MIXED", dtype=object)
    out[(rS >= HI_RANK) & (rD >= HI_RANK)] = "SLIP+TREMOR"
    out[(rS >= HI_RANK) & (rD < LO_RANK)] = "SLIP-ONLY"
    out[(rD >= HI_RANK) & (rS < LO_RANK)] = "TREMOR-ONLY"
    out[(rS < LO_RANK) & (rD < LO_RANK)] = "QUIET"
    return out


# ======================================================================================
# 0. Load everything
# ======================================================================================
log("loading slip model (load_solution + invert, ~1 min) ...")
sys.path.insert(0, str(STRAIN / "src"))
from sse_strain import load_fault_ascii, load_solution                # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator             # noqa: E402

fault = load_fault_ascii(str(STRAIN / "data/cascadia_mesh.txt"))
vlon, vlat = np.asarray(fault.lon), np.asarray(fault.lat)
sol = load_solution(str(STRAIN / "data/2026-08-18"), fault)
sol.invert(verbose=False)
_, Ld = sol.component_slip_patterns()
slabf = slab_depth_interpolator(fault)

pz = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
PLAT, PLON = pz["patch_lat"], pz["patch_lon"]
PDEPTH, PAREA = pz["patch_depth_km"], pz["patch_area_km2"]
LABEL = pz["label"]
DAYS = np.asarray(pz["days"], dtype=str)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NTRI = len(PLAT)
PX, PY = to_km(PLON, PLAT)
TRI_SEG = seg_of_lat(PLAT)
log(f"  mesh {NTRI} triangles, lat {PLAT.min():.2f}..{PLAT.max():.2f}, "
    f"days {DAYS[0]}..{DAYS[-1]} ({len(DAYS)})")

A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz",
                              allow_pickle=True)["A_raw"])
CUM_MM = 1e3 * (-(Ld @ A_raw.T))                       # (742, T) mm, thrust positive
pos_rate = np.clip(moving_average_centred(np.diff(CUM_MM, axis=1), 7, axis=1), 0.0, None)
rate_days = pd.DatetimeIndex(pd.to_datetime(DAYS[1:]))
log(f"  cum slip {CUM_MM.shape} mm; pos daily rate {pos_rate.shape}")

# vertex-sharing adjacency
vmap = defaultdict(list)
for i in range(NTRI):
    for k in range(3):
        vmap[(round(float(vlon[i, k]), 6), round(float(vlat[i, k]), 6))].append(i)
ADJ = np.zeros((NTRI, NTRI), dtype=bool)
np.fill_diagonal(ADJ, True)
for tris in vmap.values():
    for a in tris:
        ADJ[a, tris] = True
NB_SIZE = ADJ.sum(1)
# a triangle with fewer than the full complement of vertex-neighbours sits on the mesh boundary
ON_EDGE = NB_SIZE < NB_SIZE.max()
log(f"  vertex-sharing neighbourhood size: median {int(np.median(NB_SIZE))} "
    f"(min {NB_SIZE.min()}, max {NB_SIZE.max()}); {int(ON_EDGE.sum())} boundary triangles")

# tremor
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TRX, TRY = to_km(tr["lon"].values, tr["lat"].values)
tree_tri = cKDTree(np.column_stack([PX, PY]))
dist, idx = tree_tri.query(np.column_stack([TRX, TRY]), k=1)
keep = (dist <= MAX_ASSIGN_KM) & (tr["t"].values >= np.datetime64(WIN_LO)) \
    & (tr["t"].values <= np.datetime64(WIN_HI) + np.timedelta64(1, "D"))
trk = tr.loc[keep].copy()
trk["tri"] = idx[keep]
trk["day"] = trk["t"].dt.floor("D")
TREE_DET = cKDTree(np.column_stack([TRX, TRY]))        # ALL detections, for Part B
TR_LAT_ALL = tr["lat"].values
TR_T_ALL = tr["t"].values
log(f"  {len(tr):,} detections; {len(trk):,} assigned within {MAX_ASSIGN_KM:.0f} km "
    f"and inside {WIN_LO}..{WIN_HI}")

uday = np.array(sorted(trk["day"].unique()))
uday_idx = pd.DatetimeIndex(uday)
day_pos = pd.Series(np.arange(len(uday)), index=uday_idx)
C = np.zeros((len(uday), NTRI), dtype=np.int32)
g = trk.groupby([trk["day"], trk["tri"]]).size()
C[day_pos.loc[pd.DatetimeIndex(g.index.get_level_values(0))].values,
  g.index.get_level_values(1).values] = g.values
C_nb = C @ ADJ.astype(np.int32)
log(f"  {len(uday):,} distinct tremor days on the mesh "
    f"{uday_idx[0].date()}..{uday_idx[-1].date()}")


# ======================================================================================
# PART A — occupancy classification
# ======================================================================================
print("\n" + "=" * 100)
print("PART A — occupancy classification per triangle and epoch (Q5)")
print("=" * 100)

S = np.zeros((NTRI, 3))          # mm/yr
D = np.zeros((NTRI, 3))          # tremor-days/yr (D3_nb)
Draw = np.zeros((NTRI, 3), dtype=int)
Down_raw = np.zeros((NTRI, 3), dtype=int)
yrs = []
for k, (name, a, b) in enumerate(EPOCHS):
    ta, tb = pd.Timestamp(a), pd.Timestamp(b)
    yr = (tb - ta).days / 365.25 + 1 / 365.25
    yrs.append(yr)
    dsel = (uday_idx >= ta) & (uday_idx <= tb)
    Draw[:, k] = (C_nb[dsel] >= DAY_DET_THR).sum(axis=0)
    Down_raw[:, k] = (C[dsel] >= 1).sum(axis=0)
    D[:, k] = Draw[:, k] / yr
    smask = np.asarray((rate_days >= ta) & (rate_days <= tb))
    S[:, k] = pos_rate[:, smask].sum(axis=1) / yr
    print(f"[{name}] {a}..{b}  {yr:.2f} yr | tremor days on mesh {int(dsel.sum()):4d} | "
          f"detections {int(C[dsel].sum()):,} | median D3_nb/yr (lat>=44N) "
          f"{np.median(D[PLAT >= 44, k]):.1f} | median S mm/yr (lat>=44N) "
          f"{np.median(S[PLAT >= 44, k]):.2f}")

# window totals (absolute map)
wsel = np.asarray((rate_days >= pd.Timestamp(WIN_LO)) & (rate_days <= pd.Timestamp(WIN_HI)))
S_tot = pos_rate[:, wsel].sum(axis=1)
D_tot = (C_nb >= DAY_DET_THR).sum(axis=0)
Down_tot = (C >= 1).sum(axis=0)
N_tot = C.sum(axis=0)
WIN_YR = (pd.Timestamp(WIN_HI) - pd.Timestamp(WIN_LO)).days / 365.25 + 1 / 365.25
print(f"[window] {WIN_LO}..{WIN_HI} = {WIN_YR:.2f} yr; S_tot median {np.median(S_tot):.1f} mm, "
      f"p90 {np.percentile(S_tot, 90):.1f} mm; D3_nb_tot median {np.median(D_tot):.0f} d")

# ---- ranks within (segment, epoch) -----------------------------------------------------
rS = np.full((NTRI, 3), np.nan)
rD = np.full((NTRI, 3), np.nan)
CLS = np.empty((NTRI, 3), dtype=object)
print("\n--- rank normalisation (percentile ranks within segment x epoch; ties get average rank) ---")
for seg in ("north", "central", "south"):
    m = TRI_SEG == seg
    for k, (name, _, _) in enumerate(EPOCHS):
        rS[m, k] = pct_rank(S[m, k])
        rD[m, k] = pct_rank(D[m, k])
    tieD = [(D[m, k] == 0).mean() for k in range(3)]
    tieS = [(S[m, k] == 0).mean() for k in range(3)]
    print(f"  {seg:8s} n {int(m.sum()):3d} | frac D3_nb == 0 per epoch "
          f"{tieD[0]:.2f}/{tieD[1]:.2f}/{tieD[2]:.2f} | frac S == 0 "
          f"{tieS[0]:.3f}/{tieS[1]:.3f}/{tieS[2]:.3f}")
for k in range(3):
    CLS[:, k] = classify(rS[:, k], rD[:, k])

persistent = np.where((CLS[:, 0] == CLS[:, 1]) & (CLS[:, 1] == CLS[:, 2]), CLS[:, 0], "varies")
srows = []

print("\n--- how much CONTENT is in the classes: co-location of slip and tremor ranks, and how "
      "reproducible the ranks themselves are ---")
print("    (the persistence test below compares class maps against SHUFFLED labels, an extremely "
      "weak null; these\n     Spearman numbers say how much of that persistence is simply the "
      "underlying quantities barely changing)")
for seg in ("north", "central", "south"):
    m = TRI_SEG == seg
    rs_ = [spearmanr(rS[m, k], rD[m, k]).statistic for k in range(3)]
    aS = spearmanr(rS[m, 0], rS[m, 2]).statistic
    aD = spearmanr(rD[m, 0], rD[m, 2]).statistic
    print(f"  {seg:8s} Spearman(rank_S, rank_D) per epoch "
          f"{rs_[0]:+.3f}/{rs_[1]:+.3f}/{rs_[2]:+.3f} | E1-vs-E3 reproducibility: "
          f"rank_S {aS:+.3f}, rank_D {aD:+.3f}")
    srows.append(dict(kind="rank_diagnostic", seg=seg, epoch="all", cls="",
                      n_tri=int(m.sum()), spearman_S_D_E1=float(rs_[0]),
                      spearman_S_D_E2=float(rs_[1]), spearman_S_D_E3=float(rs_[2]),
                      spearman_S_E1_E3=float(aS), spearman_D_E1_E3=float(aD)))

# ---- per segment x epoch fractions -----------------------------------------------------
print("\n--- class occupancy: fraction of triangles (and of AREA) per segment and epoch ---")
for seg in ("north", "central", "south"):
    m = TRI_SEG == seg
    atot = PAREA[m].sum()
    print(f"\n  {seg.upper()} (n = {int(m.sum())} triangles, {atot:,.0f} km2)")
    hdr = f"    {'class':13s}" + "".join(f"{name:>20s}" for name, _, _ in EPOCHS)
    print(hdr)
    for cl in CLASSES:
        cells = []
        for k, (name, _, _) in enumerate(EPOCHS):
            mm = m & (CLS[:, k] == cl)
            ft = mm.sum() / m.sum()
            fa = PAREA[mm].sum() / atot
            dep = PDEPTH[mm]
            cells.append(f"{ft * 100:5.1f}% /{fa * 100:5.1f}%a")
            srows.append(dict(kind="occupancy", seg=seg, epoch=name, cls=cl,
                              n_tri=int(mm.sum()), frac_tri=float(ft), frac_area=float(fa),
                              mean_depth_km=float(np.mean(dep)) if mm.sum() else np.nan,
                              median_depth_km=float(np.median(dep)) if mm.sum() else np.nan,
                              n_tri_seg=int(m.sum())))
        print(f"    {cl:13s}" + "".join(f"{c:>20s}" for c in cells))

print("\n--- along-dip structure: mean [median] Slab2 depth (km) of each class, per segment "
      "(pooled over the 3 epochs) ---")
for seg in ("north", "central", "south"):
    m = TRI_SEG == seg
    parts = []
    for cl in CLASSES:
        dep = np.concatenate([PDEPTH[m & (CLS[:, k] == cl)] for k in range(3)])
        parts.append(f"{cl} {np.mean(dep):.1f} [{np.median(dep):.1f}] (n{len(dep)})"
                     if len(dep) else f"{cl} --")
        srows.append(dict(kind="depth_pooled", seg=seg, epoch="all", cls=cl,
                          n_tri=len(dep), frac_tri=np.nan, frac_area=np.nan,
                          mean_depth_km=float(np.mean(dep)) if len(dep) else np.nan,
                          median_depth_km=float(np.median(dep)) if len(dep) else np.nan,
                          n_tri_seg=int(m.sum())))
    print(f"  {seg:8s} " + " | ".join(parts))
    print(f"           segment triangle depth range {PDEPTH[m].min():.1f}..{PDEPTH[m].max():.1f} km, "
          f"median {np.median(PDEPTH[m]):.1f}")

# ---- persistence vs shuffled-label null -------------------------------------------------
print("\n--- persistence: same class in all three epochs, vs independent within-segment shuffles "
      f"({N_PERM_CLASS} draws) ---")
for seg in ("north", "central", "south", "ALL(lat>=44)"):
    m = (PLAT >= 44.0) if seg.startswith("ALL") else (TRI_SEG == seg)
    idxs = np.where(m)[0]
    obs = float(np.mean(persistent[idxs] != "varies"))
    # null: shuffle each epoch's labels independently WITHIN each segment
    segs_here = np.unique(TRI_SEG[idxs])
    null = np.zeros(N_PERM_CLASS)
    perclass_null = {cl: np.zeros(N_PERM_CLASS) for cl in CLASSES}
    for it in range(N_PERM_CLASS):
        sh = np.empty((len(idxs), 3), dtype=object)
        for sg in segs_here:
            sm = TRI_SEG[idxs] == sg
            for k in range(3):
                sh[sm, k] = RNG.permutation(CLS[idxs][sm, k])
        same = (sh[:, 0] == sh[:, 1]) & (sh[:, 1] == sh[:, 2])
        null[it] = same.mean()
        for cl in CLASSES:
            perclass_null[cl][it] = (same & (sh[:, 0] == cl)).mean()
    p = (1 + np.sum(null >= obs)) / (N_PERM_CLASS + 1)
    print(f"  {seg:14s} n {len(idxs):3d}  persistent {obs * 100:5.1f}%  vs null "
          f"{null.mean() * 100:5.1f}% [{np.percentile(null, 5) * 100:.1f}, "
          f"{np.percentile(null, 95) * 100:.1f}]  p {p:.4f}")
    srows.append(dict(kind="persistence", seg=seg, epoch="all", cls="ANY",
                      n_tri=len(idxs), frac_tri=obs, frac_area=np.nan,
                      mean_depth_km=np.nan, median_depth_km=np.nan, n_tri_seg=len(idxs),
                      null_mean=float(null.mean()), null_p5=float(np.percentile(null, 5)),
                      null_p95=float(np.percentile(null, 95)), p_perm=float(p)))
    for cl in CLASSES:
        o = float(np.mean(persistent[idxs] == cl))
        nn = perclass_null[cl]
        pc = (1 + np.sum(nn >= o)) / (N_PERM_CLASS + 1)
        dep = PDEPTH[idxs][persistent[idxs] == cl]
        print(f"      {cl:13s} obs {o * 100:5.1f}% ({int(np.sum(persistent[idxs] == cl)):3d} tri)"
              f"  null {nn.mean() * 100:5.1f}%  p {pc:.4f}"
              + (f"   depth {np.mean(dep):.1f} +- {np.std(dep):.1f} km" if len(dep) else ""))
        srows.append(dict(kind="persistence", seg=seg, epoch="all", cls=cl,
                          n_tri=int(np.sum(persistent[idxs] == cl)), frac_tri=o,
                          frac_area=float(PAREA[idxs][persistent[idxs] == cl].sum()
                                          / PAREA[idxs].sum()),
                          mean_depth_km=float(np.mean(dep)) if len(dep) else np.nan,
                          median_depth_km=float(np.median(dep)) if len(dep) else np.nan,
                          n_tri_seg=len(idxs), null_mean=float(nn.mean()),
                          null_p5=float(np.percentile(nn, 5)),
                          null_p95=float(np.percentile(nn, 95)), p_perm=float(pc)))

pd.DataFrame(srows).to_csv(OUT / "results/s7_occupancy_summary.csv", index=False)
log("[write] results/s7_occupancy_summary.csv")

cols = dict(tri=np.arange(NTRI), lat=PLAT, lon=PLON, depth_km=PDEPTH, area_km2=PAREA, seg=TRI_SEG)
for k, (name, _, _) in enumerate(EPOCHS):
    cols[f"S_mm_yr_{name}"] = S[:, k]
    cols[f"D3nb_days_yr_{name}"] = D[:, k]
    cols[f"D3nb_days_{name}"] = Draw[:, k]
    cols[f"rank_S_{name}"] = rS[:, k]
    cols[f"rank_D_{name}"] = rD[:, k]
    cols[f"class_{name}"] = CLS[:, k]
cols["class_persistent"] = persistent
cols["S_tot_mm"] = S_tot
cols["D3nb_tot_days"] = D_tot
cols["D_own_tot_days"] = Down_tot
cols["N_det_tot"] = N_tot
cols["on_mesh_edge"] = ON_EDGE
tri_df = pd.DataFrame(cols)
tri_df.to_csv(OUT / "results/s7_occupancy_triangles.csv", index=False)
log(f"[write] results/s7_occupancy_triangles.csv ({NTRI} rows)")

# ---- absolute map ----------------------------------------------------------------------
print("\n--- ABSOLUTE map (no ranks): candidate quiet-slip and tremor-without-slip triangles ---")
det_seg = seg_of_lat(trk["lat"].values)
det_depth = slabf(trk["lon"].values, trk["lat"].values)
tremor_band = {}
for seg in ("north", "central", "south"):
    dd = det_depth[(det_seg == seg) & np.isfinite(det_depth)]
    if len(dd) < 50:
        tremor_band[seg] = (np.nan, np.nan, np.nan, 0)
        continue
    tremor_band[seg] = (float(np.percentile(dd, 10)), float(np.median(dd)),
                        float(np.percentile(dd, 90)), len(dd))
    print(f"  {seg:8s} tremor Slab2 depth p10 {tremor_band[seg][0]:.1f} / p50 "
          f"{tremor_band[seg][1]:.1f} / p90 {tremor_band[seg][2]:.1f} km  "
          f"(n {tremor_band[seg][3]:,} detections)")

# segment-relative slip percentiles (the literal mm thresholds are degenerate -- see below)
S_p75 = np.full(NTRI, np.nan)
S_p25 = np.full(NTRI, np.nan)
for sg in ("north", "central", "south"):
    m = TRI_SEG == sg
    S_p75[m] = np.percentile(S_tot[m], 75)
    S_p25[m] = np.percentile(S_tot[m], 25)

qs = (S_tot >= QS_SLIP_MM) & (D_tot < QS_DAYS)                      # literal spec criterion
to = (D_tot >= TO_DAYS) & (S_tot < TO_SLIP_MM)                      # literal spec criterion
qs75 = (S_tot >= S_p75) & (D_tot < QS_DAYS)                         # percentile version
to25 = (D_tot >= TO_DAYS) & (S_tot <= S_p25)                        # percentile version

print(f"\n  !! the literal slip thresholds are DEGENERATE: only the POSITIVE daily rate is "
      f"accumulated, so every\n     triangle rectifies model noise upward -- S_tot ranges "
      f"{S_tot.min():.0f}..{S_tot.max():.0f} mm (median {np.median(S_tot):.0f} mm) over "
      f"{WIN_YR:.1f} yr.\n     S_tot >= {QS_SLIP_MM:.0f} mm passes {int((S_tot >= QS_SLIP_MM).sum())}"
      f"/{NTRI} triangles and S_tot < {TO_SLIP_MM:.0f} mm passes "
      f"{int((S_tot < TO_SLIP_MM).sum())}/{NTRI}.\n     The literal criteria therefore reduce to a "
      f"pure TREMOR criterion (and to the empty set). Percentile\n     versions relative to each "
      f"segment's own S_tot distribution are reported alongside and are the ones to read.")
print(f"     segment S_tot p25 / p75 (mm): " + ", ".join(
    f"{sg} {np.percentile(S_tot[TRI_SEG == sg], 25):.0f} / "
    f"{np.percentile(S_tot[TRI_SEG == sg], 75):.0f}" for sg in ("north", "central", "south")))

for lab, mask, desc in (
        ("quiet_slip_literal", qs, f"S_tot >= {QS_SLIP_MM:.0f} mm and D3_nb_tot < {QS_DAYS} d"),
        ("quiet_slip_p75", qs75, f"S_tot >= segment p75 and D3_nb_tot < {QS_DAYS} d"),
        ("tremor_no_slip_literal", to, f"D3_nb_tot >= {TO_DAYS} d and S_tot < {TO_SLIP_MM:.0f} mm"),
        ("tremor_no_slip_p25", to25, f"D3_nb_tot >= {TO_DAYS} d and S_tot <= segment p25")):
    per = ", ".join(f"{sg} {int((mask & (TRI_SEG == sg)).sum())}"
                    for sg in ("north", "central", "south"))
    print(f"  {lab:24s} ({desc:52s}): {int(mask.sum()):3d} triangles  [{per}]")

qrows = []
for tag, mask in (("quiet_slip_literal", qs), ("quiet_slip_p75", qs75),
                  ("tremor_no_slip_literal", to), ("tremor_no_slip_p25", to25)):
    for i in np.where(mask)[0]:
        sg = TRI_SEG[i]
        p10, p50, p90, _ = tremor_band.get(sg, (np.nan,) * 3 + (0,))
        rel = ("no_band" if not np.isfinite(p10) else
               "updip_of_band" if PDEPTH[i] < p10 else
               "downdip_of_band" if PDEPTH[i] > p90 else "in_band")
        qrows.append(dict(kind=tag, tri=int(i), lat=float(PLAT[i]), lon=float(PLON[i]),
                          depth_km=float(PDEPTH[i]), area_km2=float(PAREA[i]), seg=sg,
                          S_tot_mm=float(S_tot[i]), seg_S_tot_p25=float(S_p25[i]),
                          seg_S_tot_p75=float(S_p75[i]),
                          D3nb_tot_days=int(D_tot[i]),
                          D_own_tot_days=int(Down_tot[i]), N_det_tot=int(N_tot[i]),
                          seg_tremor_p10_km=p10, seg_tremor_p50_km=p50, seg_tremor_p90_km=p90,
                          position_vs_tremor_band=rel,
                          on_mesh_edge=bool(ON_EDGE[i]),
                          class_persistent=str(persistent[i])))
qdf = pd.DataFrame(qrows)
qdf.to_csv(OUT / "results/s7_quiet_slip_candidates.csv", index=False)
log(f"[write] results/s7_quiet_slip_candidates.csv ({len(qdf)} rows)")

for tag in ("quiet_slip_literal", "quiet_slip_p75", "tremor_no_slip_literal",
            "tremor_no_slip_p25"):
    sub = qdf[qdf["kind"] == tag] if len(qdf) else qdf
    if not len(sub):
        print(f"\n  {tag}: none")
        continue
    print(f"\n  {tag}: {len(sub)} triangles; position vs the segment tremor depth band "
          f"(p10/p90 of the Slab2 depth of that segment's detections): "
          f"{sub['position_vs_tremor_band'].value_counts().to_dict()}")
    for sg in ("north", "central", "south"):
        s2 = sub[sub.seg == sg]
        if len(s2):
            print(f"     {sg:8s} n {len(s2):3d}  depth {s2.depth_km.min():.1f}.."
                  f"{s2.depth_km.max():.1f} km (median {s2.depth_km.median():.1f}); "
                  f"lat {s2.lat.min():.2f}..{s2.lat.max():.2f}; "
                  f"{s2['position_vs_tremor_band'].value_counts().to_dict()}")
    # the full list is in the CSV; print it here only for the (small) percentile sets
    if tag.endswith("_p75") or tag.endswith("_p25"):
        show = sub.sort_values(["seg", "lat"])[["tri", "lat", "lon", "depth_km", "seg",
                                                "S_tot_mm", "D3nb_tot_days", "D_own_tot_days",
                                                "position_vs_tremor_band", "class_persistent"]]
        with pd.option_context("display.width", 220, "display.max_rows", 400):
            print(show.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

# persistent SLIP-ONLY / TREMOR-ONLY triangles (the rank-based answer to the same question)
print("\n--- persistent rank-based SLIP-ONLY and TREMOR-ONLY triangles (same class in E1, E2 "
      "and E3) ---")
for cl in ("SLIP-ONLY", "TREMOR-ONLY"):
    m = persistent == cl
    if not m.any():
        print(f"  {cl}: NONE anywhere on the mesh")
        continue
    sub = pd.DataFrame(dict(tri=np.where(m)[0], lat=PLAT[m], lon=PLON[m], depth_km=PDEPTH[m],
                            seg=TRI_SEG[m], S_tot_mm=S_tot[m], D3nb_tot_days=D_tot[m],
                            D_own_tot_days=Down_tot[m], on_mesh_edge=ON_EDGE[m],
                            deg_to_mesh_lat_limit=np.minimum(PLAT[m] - PLAT.min(),
                                                             PLAT.max() - PLAT[m])))
    sub["position_vs_tremor_band"] = [
        ("no_band" if not np.isfinite(tremor_band[g][0]) else
         "updip_of_band" if d < tremor_band[g][0] else
         "downdip_of_band" if d > tremor_band[g][2] else "in_band")
        for d, g in zip(sub.depth_km, sub.seg)]
    print(f"  {cl}: {int(m.sum())} triangles, depth {PDEPTH[m].mean():.1f} +- "
          f"{PDEPTH[m].std():.1f} km; {int(ON_EDGE[m].sum())} of them on the MESH BOUNDARY; "
          f"median distance to the mesh latitude limit "
          f"{np.median(np.minimum(PLAT[m] - PLAT.min(), PLAT.max() - PLAT[m])):.2f} deg")
    with pd.option_context("display.width", 220):
        print(sub.sort_values(["seg", "lat"]).to_string(index=False,
                                                        float_format=lambda v: f"{v:.2f}"))


# ======================================================================================
# PART B — along-strike terminations
# ======================================================================================
print("\n" + "=" * 100)
print("PART B — along-strike terminations of tremor and modeled slip (Q17)")
print("=" * 100)

tvs = pd.read_csv(SSE / "tremor_vs_slip_patches.csv")
conf = tvs[(tvs.catalog == "thrust") & (tvs.confirmed.astype(bool))
           & (tvs.seg.isin(["north", "central"])) & (tvs.Mw >= 6.0)
           & (~tvs.broad_pulse.astype(bool))].sort_values("start").reset_index(drop=True)
s4 = pd.read_csv(OUT / "results/s4_fronts_events.csv").set_index("id")
print(f"  CONFIRMED events n = {len(conf)}  {conf.seg.value_counts().to_dict()}")

erows = []
for _, r in conf.iterrows():
    pid, start, end = int(r["id"]), str(r["start"]), str(r["end"])
    ts, te = DAY_INDEX[start], DAY_INDEX[end]
    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    rec = dict(id=pid, start=start, end=end, seg=r["seg"], Mw=float(r["Mw"]),
               dur_d=int(r["dur_d"]), n_tri_F=int(F.sum()),
               F_lat_min=float(PLAT[F].min()), F_lat_max=float(PLAT[F].max()))

    # ---- tremor extent
    ib = TREE_DET.query_ball_point(np.column_stack([PX[F], PY[F]]), r=EVENT_TREMOR_RADIUS_KM)
    near = np.unique(np.concatenate([np.asarray(i, int) for i in ib])) if len(ib) else \
        np.array([], int)
    intime = near[(TR_T_ALL[near] >= np.datetime64(start))
                  & (TR_T_ALL[near] <= np.datetime64(end) + np.timedelta64(1, "D"))]
    rec["n_det"] = int(intime.size)
    if intime.size >= 20:
        la = TR_LAT_ALL[intime]
        rec["trem_south_lat"] = float(np.percentile(la, 5))
        rec["trem_north_lat"] = float(np.percentile(la, 95))
        # daily centroid latitude for the fallback direction
        dnum = (TR_T_ALL[intime].astype("datetime64[D]")
                - np.datetime64(start)).astype(int)
        dfd = pd.DataFrame(dict(d=dnum, lat=la)).groupby("d")["lat"].agg(["mean", "size"])
        dfd = dfd[dfd["size"] >= 5]
        if len(dfd) >= 5:
            rec["v_t_recomp_deg_d"] = float(theilslopes(dfd["mean"].values,
                                                        dfd.index.values.astype(float))[0])
        else:
            rec["v_t_recomp_deg_d"] = np.nan
    else:
        rec["trem_south_lat"] = rec["trem_north_lat"] = rec["v_t_recomp_deg_d"] = np.nan

    # ---- slip extent
    acc = np.clip(CUM_MM[:, te] - CUM_MM[:, max(ts - EVENT_PRE_DAYS, 0)], 0.0, None)
    band = (PLAT >= PLAT[F].min() - EVENT_LAT_PAD_DEG) \
        & (PLAT <= PLAT[F].max() + EVENT_LAT_PAD_DEG) & (acc >= EVENT_EDGE_ACC_MM)
    rec["n_tri_slip_extent"] = int(band.sum())
    if band.any():
        rec["slip_south_lat"] = float(PLAT[band].min())
        rec["slip_north_lat"] = float(PLAT[band].max())
        rec["max_acc_mm"] = float(acc[F].max())
    else:
        rec["slip_south_lat"] = rec["slip_north_lat"] = rec["max_acc_mm"] = np.nan

    rec["trem_extent_deg"] = rec["trem_north_lat"] - rec["trem_south_lat"]
    rec["slip_extent_deg"] = rec["slip_north_lat"] - rec["slip_south_lat"]
    rec["d_north_km"] = 111.0 * (rec["trem_north_lat"] - rec["slip_north_lat"])
    rec["d_south_km"] = 111.0 * (rec["trem_south_lat"] - rec["slip_south_lat"])

    # ---- propagation direction
    v_t = float(s4.loc[pid, "v_t_km_d"]) if pid in s4.index else np.nan
    rec["v_t_km_d_s4"] = v_t
    if np.isfinite(v_t) and abs(v_t) > DIR_EPS:
        rec["direction"] = "north" if v_t > 0 else "south"
        rec["dir_source"] = "s4_v_t"
    elif np.isfinite(rec["v_t_recomp_deg_d"]) and rec["v_t_recomp_deg_d"] != 0:
        rec["direction"] = "north" if rec["v_t_recomp_deg_d"] > 0 else "south"
        rec["dir_source"] = "recomputed_theilsen"
    else:
        rec["direction"] = "unclear"
        rec["dir_source"] = "none"
    erows.append(rec)

ev = pd.DataFrame(erows)
ev["year"] = pd.to_datetime(ev["start"]).dt.year
ev["epoch"] = np.select([ev.year <= 2014, ev.year <= 2019], ["E1", "E2"], "E3")
ev.to_csv(OUT / "results/s7_terminations_events.csv", index=False)
log(f"[write] results/s7_terminations_events.csv ({len(ev)} rows)")

with pd.option_context("display.width", 250):
    print("\n" + ev[["id", "start", "seg", "Mw", "n_det", "trem_south_lat", "trem_north_lat",
                     "slip_south_lat", "slip_north_lat", "d_south_km", "d_north_km",
                     "direction", "dir_source", "epoch"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

good = ev.dropna(subset=["trem_south_lat", "trem_north_lat", "slip_south_lat", "slip_north_lat"])
print(f"\n  events with both extents defined: {len(good)} of {len(ev)}")
print(f"  tremor extent (deg): median {ev.trem_extent_deg.median():.2f} "
      f"[{ev.trem_extent_deg.min():.2f}, {ev.trem_extent_deg.max():.2f}];  "
      f"slip extent {ev.slip_extent_deg.median():.2f} "
      f"[{ev.slip_extent_deg.min():.2f}, {ev.slip_extent_deg.max():.2f}]")
print(f"  tremor - slip termination difference (km):  NORTH end median "
      f"{good.d_north_km.median():+.1f} (median |.| {good.d_north_km.abs().median():.1f}); "
      f"SOUTH end median {good.d_south_km.median():+.1f} "
      f"(median |.| {good.d_south_km.abs().median():.1f})")
print(f"  direction: {ev.direction.value_counts().to_dict()}  "
      f"(source {ev.dir_source.value_counts().to_dict()})")

# ---- barrier test ------------------------------------------------------------------------
term_lat = np.concatenate([good.trem_south_lat.values, good.trem_north_lat.values])
term_end = np.array(["south"] * len(good) + ["north"] * len(good))
term_epoch = np.concatenate([good.epoch.values, good.epoch.values])
term_dir = np.concatenate([good.direction.values, good.direction.values])
slip_term_lat = np.concatenate([good.slip_south_lat.values, good.slip_north_lat.values])

bin_lo = np.floor(min(term_lat.min(), slip_term_lat.min()) / BIN_W) * BIN_W
bin_hi = np.ceil(max(term_lat.max(), slip_term_lat.max()) / BIN_W) * BIN_W
edges = np.arange(bin_lo, bin_hi + BIN_W / 2, BIN_W)
centres = edges[:-1] + BIN_W / 2
NB = len(centres)
print(f"\n  0.25-deg bins: {NB} from {edges[0]:.2f} to {edges[-1]:.2f} "
      f"({2 * len(good)} tremor terminations)")

SEG_STRICT = {"north": (46.0, float(PLAT.max())), "central": (44.0, 46.0)}
SEG_DOMAIN = {}
for sg in ("north", "central"):
    sub = good[good.seg == sg]
    lo, hi = SEG_STRICT[sg]
    if len(sub):
        lo = min(lo, float(sub.trem_south_lat.min()))
        hi = max(hi, float(sub.trem_north_lat.max()))
    SEG_DOMAIN[sg] = (lo, hi)
print(f"  randomisation domains: strict {SEG_STRICT}, widened {SEG_DOMAIN}")


def null_counts(sub, domains, nperm=N_PERM_TERM, rng=None):
    """Uniformly reposition each event's tremor extent inside its segment domain."""
    rng = rng or np.random.default_rng(SEED)
    L = (sub.trem_north_lat - sub.trem_south_lat).values
    lo = np.array([domains[s][0] for s in sub.seg])
    hi = np.array([domains[s][1] for s in sub.seg])
    room = np.clip(hi - lo - L, 0.0, None)
    out = np.zeros((nperm, NB), dtype=np.int32)
    for it in range(nperm):
        a = lo + rng.random(len(sub)) * room
        h1, _ = np.histogram(a, bins=edges)
        h2, _ = np.histogram(a + L, bins=edges)
        out[it] = h1 + h2
    return out, float(np.mean(room))


obs_counts, _ = np.histogram(term_lat, bins=edges)
brows = []
sig_bins = {}
raw_bins = {}
NULLS = {}
for tag, domains in (("domain", SEG_DOMAIN), ("strict", SEG_STRICT)):
    nulls, mroom = null_counts(good, domains, rng=np.random.default_rng(SEED))
    NULLS[tag] = nulls
    # A bin the null can never reach (its null count is identically 0) is outside every segment
    # domain; a single observed termination there would score the minimum possible p, which is an
    # artefact of the domain, not a barrier.  Exclude such bins from the test and the correction.
    testable = nulls.max(axis=0) > 0
    n_test = int(testable.sum())
    pvals = np.full(NB, np.nan)
    pvals[testable] = (1 + (nulls[:, testable] >= obs_counts[None, testable]).sum(axis=0)) \
        / (N_PERM_TERM + 1)
    pbonf = np.clip(pvals * n_test, 0, 1)
    sig = np.where(np.isfinite(pbonf) & (pbonf < 0.05))[0]
    sig_bins[tag] = sig
    raw_bins[tag] = np.where(np.isfinite(pvals) & (pvals < 0.05))[0]
    print(f"\n  --- barrier test, '{tag}' domain (mean repositioning freedom "
          f"{mroom:.2f} deg; {n_test} of {NB} bins testable, "
          f"{int(np.sum(~testable & (obs_counts > 0)))} unreachable bins hold "
          f"{int(obs_counts[~testable].sum())} observed terminations) ---")
    siglist = [f"{centres[i]:.2f}N (obs {obs_counts[i]}, null {nulls[:, i].mean():.2f}, "
               f"p_bonf {pbonf[i]:.3g})" for i in sig]
    print(f"      bins with Bonferroni p < 0.05: {len(sig)}  {siglist}")
    rawlist = [f"{centres[i]:.2f}N (obs {obs_counts[i]}, null {nulls[:, i].mean():.2f}, "
               f"p {pvals[i]:.4f})" for i in raw_bins[tag]]
    print(f"      bins with RAW p < 0.05: {len(raw_bins[tag])}  {rawlist}")
    top = np.argsort(-obs_counts)[:8]
    for i in sorted(top):
        flag = "" if testable[i] else "  [UNREACHABLE under this null]"
        print(f"      bin {centres[i]:6.2f}N  obs {obs_counts[i]:2d}  null "
              f"{nulls[:, i].mean():5.2f} [p95 {np.percentile(nulls[:, i], 95):.0f}]  "
              f"p {pvals[i]:.4f}  p_bonf {pbonf[i]:.3g}{flag}")
    if tag == "domain":
        for i in range(NB):
            row = dict(bin_centre=float(centres[i]), obs_tremor_terminations=int(obs_counts[i]),
                       null_mean=float(nulls[:, i].mean()),
                       null_p95=float(np.percentile(nulls[:, i], 95)),
                       testable=bool(testable[i]),
                       p=float(pvals[i]), p_bonferroni=float(pbonf[i]),
                       n_bins_tested=n_test)
            for ep in ("E1", "E2", "E3"):
                row[f"obs_{ep}"] = int(np.sum((term_epoch == ep)
                                              & (term_lat >= edges[i]) & (term_lat < edges[i + 1])))
            for en in ("south", "north"):
                row[f"obs_{en}_end"] = int(np.sum((term_end == en) & (term_lat >= edges[i])
                                                  & (term_lat < edges[i + 1])))
            for dr in ("north", "south"):
                row[f"obs_prop_{dr}"] = int(np.sum((term_dir == dr) & (term_lat >= edges[i])
                                                   & (term_lat < edges[i + 1])))
            row["obs_slip_terminations"] = int(np.sum((slip_term_lat >= edges[i])
                                                      & (slip_term_lat < edges[i + 1])))
            brows.append(row)
pd.DataFrame(brows).to_csv(OUT / "results/s7_termination_bins.csv", index=False)
log(f"[write] results/s7_termination_bins.csv ({len(brows)} rows)")

# ---- epoch persistence of the barrier bins -------------------------------------------------
print("\n  --- do the same bins act as barriers in each epoch? ('domain' null) ---")
SIG = sig_bins["domain"]
SIG_LABEL = "Bonferroni-significant"
if len(SIG) == 0:
    # fall back to the raw-p<0.05 bins so the persistence question can still be answered;
    # these are CANDIDATE barriers only, not claims.
    SIG = raw_bins["domain"]
    SIG_LABEL = "raw-p<0.05 CANDIDATE (nothing survives Bonferroni)"
    print(f"      no bin survives Bonferroni; using the {len(SIG)} raw-p<0.05 candidate bins "
          f"instead: {[f'{centres[i]:.2f}N' for i in SIG]}")
for ep in ("E1", "E2", "E3"):
    sub = good[good.epoch == ep]
    if not len(sub):
        continue
    nulls_e, _ = null_counts(sub, SEG_DOMAIN, rng=np.random.default_rng(SEED + 1))
    lat_e = np.concatenate([sub.trem_south_lat.values, sub.trem_north_lat.values])
    obs_e, _ = np.histogram(lat_e, bins=edges)
    tot_obs = int(obs_e[SIG].sum())
    tot_null = nulls_e[:, SIG].sum(axis=1)
    p_e = (1 + np.sum(tot_null >= tot_obs)) / (N_PERM_TERM + 1)
    detail = ", ".join(f"{centres[i]:.2f}N:{obs_e[i]}" for i in SIG)
    print(f"      {ep}: {len(sub)} events / {2 * len(sub)} terminations; in the barrier bins "
          f"{tot_obs} (null {tot_null.mean():.2f}, p {p_e:.4f})   [{detail}]")

# ---- slip vs tremor terminations -------------------------------------------------------
print("\n  --- do the modeled-slip terminations stop at the same latitudes? ---")
print(f"      median |tremor - slip| : north end {good.d_north_km.abs().median():.1f} km, "
      f"south end {good.d_south_km.abs().median():.1f} km, both ends pooled "
      f"{np.median(np.abs(np.concatenate([good.d_north_km, good.d_south_km]))):.1f} km")
obs_slip, _ = np.histogram(slip_term_lat, bins=edges)
if len(SIG):
    print(f"      tremor terminations in the barrier bins: {obs_counts[SIG].sum()} of "
          f"{obs_counts.sum()} ({obs_counts[SIG].sum() / obs_counts.sum() * 100:.0f} %); "
          f"slip terminations in the SAME bins: {obs_slip[SIG].sum()} of {obs_slip.sum()} "
          f"({obs_slip[SIG].sum() / obs_slip.sum() * 100:.0f} %)")
sig_slip = np.argsort(-obs_slip)[:5]
print("      most-used slip termination bins: "
      + ", ".join(f"{centres[i]:.2f}N ({obs_slip[i]})" for i in sorted(sig_slip)))
print(f"      slip extents are systematically {'WIDER' if (good.slip_extent_deg > good.trem_extent_deg).mean() > 0.5 else 'NARROWER'} "
      f"than tremor extents in {max((good.slip_extent_deg > good.trem_extent_deg).mean(), (good.slip_extent_deg < good.trem_extent_deg).mean()) * 100:.0f} % of events "
      f"(median slip {good.slip_extent_deg.median():.2f} deg vs tremor "
      f"{good.trem_extent_deg.median():.2f} deg)")

# ---- direction --------------------------------------------------------------------------
print("\n  --- northward- vs southward-propagating events: same barriers? ---")
ev_np = good[good.direction == "north"]   # northward-propagating
ev_sp = good[good.direction == "south"]   # southward-propagating
print(f"      n north-propagating {len(ev_np)}, south-propagating {len(ev_sp)}, "
      f"unclear {int((good.direction == 'unclear').sum())}")
for end in ("north", "south"):
    a = ev_np[f"trem_{end}_lat"].values
    b = ev_sp[f"trem_{end}_lat"].values
    if len(a) >= 3 and len(b) >= 3:
        u = mannwhitneyu(a, b, alternative="two-sided")
        print(f"      {end}-end latitude: N-prop median {np.median(a):.2f} (n {len(a)}) vs "
              f"S-prop {np.median(b):.2f} (n {len(b)})  MWU U {u.statistic:.1f}  p {u.pvalue:.3f}")
if len(SIG):
    for lab, gg in (("north-prop", ev_np), ("south-prop", ev_sp)):
        la = np.concatenate([gg.trem_south_lat.values, gg.trem_north_lat.values])
        h, _ = np.histogram(la, bins=edges)
        if h.sum():
            print(f"      {lab}: {h[SIG].sum()} of {h.sum()} terminations in the barrier bins "
                  f"({h[SIG].sum() / h.sum() * 100:.0f} %)")


# ======================================================================================
# FIGURES
# ======================================================================================
log("plotting ...")
verts = [np.column_stack([vlon[i], vlat[i]]) for i in range(NTRI)]


def base_map(ax):
    ax.add_collection(PolyCollection(verts, facecolors="none", edgecolors="0.85", lw=0.25))
    ax.set_xlim(PLON.min() - 0.4, PLON.max() + 0.4)
    ax.set_ylim(PLAT.min() - 0.4, PLAT.max() + 0.4)
    ax.set_aspect(1.0 / COS47)
    for yv in (44, 46):
        ax.axhline(yv, color="k", lw=0.6, ls="--", alpha=0.55)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("lon", fontsize=8)


fig = plt.figure(figsize=(20, 11.5))
gs = fig.add_gridspec(2, 10, height_ratios=[2.7, 1.35], hspace=0.22, wspace=1.5)
for k, (name, a, b) in enumerate(EPOCHS):
    ax = fig.add_subplot(gs[0, 2 * k:2 * k + 2])
    base_map(ax)
    for cl in CLASSES:
        m = CLS[:, k] == cl
        ax.scatter(PLON[m], PLAT[m], c=CLASS_COL[cl], s=17, edgecolors="none", label=cl)
    ax.set_title(f"{name}  {a[:7]}..{b[:7]}", fontsize=9)
    if k == 0:
        ax.set_ylabel("lat", fontsize=8)
        ax.legend(fontsize=6.2, loc="lower left", framealpha=0.9, markerscale=1.3)

axp = fig.add_subplot(gs[0, 6:8])
base_map(axp)
for cl in CLASSES + ["varies"]:
    m = persistent == cl
    axp.scatter(PLON[m], PLAT[m], c=CLASS_COL[cl], s=17, edgecolors="none",
                label=f"{cl} ({int(m.sum())})")
axp.set_title("PERSISTENT class (same in E1, E2, E3)", fontsize=9)
axp.legend(fontsize=6.2, loc="lower left", framealpha=0.9)

axa = fig.add_subplot(gs[0, 8:10])
base_map(axa)
axa.scatter(PLON, PLAT, c="0.88", s=9, edgecolors="none")
axa.scatter(PLON[to25], PLAT[to25], c="#f4a582", s=34, edgecolors="k", lw=0.3,
            label=f"tremor, low slip ({int(to25.sum())})")
axa.scatter(PLON[qs75], PLAT[qs75], c="#2166ac", s=34, marker="s", edgecolors="k", lw=0.3,
            label=f"high slip, no tremor ({int(qs75.sum())})")
axa.scatter(PLON[qs & ~qs75], PLAT[qs & ~qs75], facecolors="none", edgecolors="#2166ac",
            s=22, marker="s", lw=0.6,
            label=f"D3nb < {QS_DAYS} d only ({int((qs & ~qs75).sum())})")
axa.set_title(f"ABSOLUTE (no ranks), PERCENTILE criteria:\nS_tot >= segment p75 & D3nb < "
              f"{QS_DAYS} d  vs  D3nb >= {TO_DAYS} d & S_tot <= p25\n"
              f"(the literal 20 mm / 5 mm bars are degenerate: {int((S_tot >= QS_SLIP_MM).sum())}"
              f"/{NTRI} and 0/{NTRI})", fontsize=7.5)
axa.legend(fontsize=6.5, loc="lower left", framealpha=0.9)

# depth panel
axd = fig.add_subplot(gs[1, 0:5])
positions, labels, data, colors = [], [], [], []
p = 0
for sgi, seg in enumerate(("north", "central", "south")):
    m = TRI_SEG == seg
    for cl in CLASSES:
        dep = np.concatenate([PDEPTH[m & (CLS[:, k] == cl)] for k in range(3)])
        if len(dep) < 3:
            continue
        positions.append(p)
        data.append(dep)
        colors.append(CLASS_COL[cl])
        labels.append(f"{cl}\n{seg}")
        p += 1
    p += 1
bp = axd.boxplot(data, positions=positions, widths=0.72, patch_artist=True, showfliers=False)
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.85)
for med in bp["medians"]:
    med.set_color("k")
axd.set_xticks(positions)
axd.set_xticklabels(labels, fontsize=6.0, rotation=35, ha="right")
axd.set_ylabel("Slab2 depth (km)", fontsize=8)
axd.invert_yaxis()
axd.grid(alpha=0.3, axis="y")
axd.set_title("along-dip structure: Slab2 depth of each occupancy class (all 3 epochs pooled)",
              fontsize=9)

axf = fig.add_subplot(gs[1, 5:8])
w = 0.26
xb = np.arange(len(CLASSES))
for k, (name, _, _) in enumerate(EPOCHS):
    m = PLAT >= 44.0
    fr = [np.mean(CLS[m, k] == cl) for cl in CLASSES]
    axf.bar(xb + (k - 1) * w, fr, width=w, label=name,
            color=["tab:blue", "tab:orange", "tab:red"][k])
axf.set_xticks(xb)
axf.set_xticklabels(CLASSES, fontsize=6.5, rotation=20)
axf.set_ylabel("fraction of triangles (lat >= 44N)", fontsize=8)
axf.set_title("class composition per epoch\n(fixed by construction of the ranks; the SPATIAL\n"
              "pattern is what the persistence test reads)", fontsize=8)
axf.legend(fontsize=7)
axf.grid(alpha=0.3, axis="y")

axpr = fig.add_subplot(gs[1, 8:10])
segl = ["north", "central", "south", "ALL(lat>=44)"]
pdf = pd.DataFrame(srows)
pp = pdf[(pdf.kind == "persistence") & (pdf.cls == "ANY")].set_index("seg")
obsv = [pp.loc[s, "frac_tri"] * 100 for s in segl]
nulv = [pp.loc[s, "null_mean"] * 100 for s in segl]
errl = [[pp.loc[s, "null_mean"] * 100 - pp.loc[s, "null_p5"] * 100 for s in segl],
        [pp.loc[s, "null_p95"] * 100 - pp.loc[s, "null_mean"] * 100 for s in segl]]
xx = np.arange(len(segl))
axpr.bar(xx - 0.2, obsv, width=0.4, color="tab:red", label="observed")
axpr.bar(xx + 0.2, nulv, width=0.4, color="0.7", yerr=errl, capsize=3, label="shuffled null")
for i, s in enumerate(segl):
    axpr.text(i, max(obsv[i], nulv[i]) + 1.5, f"p={pp.loc[s, 'p_perm']:.3f}", ha="center",
              fontsize=6.5)
axpr.set_xticks(xx)
axpr.set_xticklabels(segl, fontsize=7)
axpr.set_ylabel("% triangles with the same class\nin all three epochs", fontsize=7.5)
axpr.set_title("persistence vs within-segment shuffles", fontsize=9)
axpr.legend(fontsize=7)
axpr.grid(alpha=0.3, axis="y")

fig.suptitle("S7 Part A — slip / tremor occupancy classes per triangle and epoch "
             "(ranks within segment x epoch; tremor-days = D3_nb)", fontsize=13)
fig.savefig(OUT / "figures/s7_occupancy_maps.png", dpi=150, bbox_inches="tight")
plt.close(fig)
log("[write] figures/s7_occupancy_maps.png")

# ---------------------------------------------------------------- terminations figure
nulls_dom, _ = null_counts(good, SEG_DOMAIN, rng=np.random.default_rng(SEED))
fig2 = plt.figure(figsize=(19, 10.5))
gs2 = fig2.add_gridspec(2, 6, height_ratios=[1.15, 1.0], hspace=0.32, wspace=0.55)

ax1 = fig2.add_subplot(gs2[0, 0:3])
ax1.bar(centres, obs_counts, width=BIN_W * 0.9, color="tab:red", alpha=0.85,
        label=f"observed tremor terminations (n {obs_counts.sum()})")
ax1.plot(centres, nulls_dom.mean(axis=0), "k-", lw=1.3, label="null mean (repositioned extents)")
ax1.fill_between(centres, np.percentile(nulls_dom, 5, axis=0),
                 np.percentile(nulls_dom, 95, axis=0), color="0.55", alpha=0.35,
                 label="null 5-95 %")
for i in SIG:
    ax1.axvline(centres[i], color="tab:blue", lw=0.9, ls="--", alpha=0.8)
ax1.set_xlabel("latitude of termination (0.25-deg bins)")
ax1.set_ylabel("count")
ax1.set_title(f"tremor along-strike terminations, {len(good)} confirmed events 2010-2026\n"
              f"dashed = candidate barrier bins (raw p<0.05); Bonferroni survivors "
              f"{len(sig_bins['domain'])} of {NB} bins", fontsize=9)
ax1.legend(fontsize=7)
ax1.grid(alpha=0.3)

ax2 = fig2.add_subplot(gs2[0, 3:6])
for ep, c in zip(("E1", "E2", "E3"), ("tab:blue", "tab:orange", "tab:red")):
    la = term_lat[term_epoch == ep]
    if la.size:
        h, _ = np.histogram(la, bins=edges)
        ax2.step(centres, h, where="mid", color=c, lw=1.6,
                 label=f"{ep} (n {la.size} terminations)")
for i in SIG:
    ax2.axvline(centres[i], color="0.4", lw=0.9, ls="--", alpha=0.7)
ax2.set_xlabel("latitude of termination")
ax2.set_ylabel("count")
ax2.set_title("persistence of the barrier latitudes across epochs", fontsize=9)
ax2.legend(fontsize=7)
ax2.grid(alpha=0.3)

ax3 = fig2.add_subplot(gs2[1, 0:2])
ax3.scatter(good.slip_north_lat, good.trem_north_lat, s=36, c="tab:red", label="north end",
            edgecolors="k", lw=0.3)
ax3.scatter(good.slip_south_lat, good.trem_south_lat, s=36, c="tab:blue", marker="s",
            label="south end", edgecolors="k", lw=0.3)
lims = [min(good.slip_south_lat.min(), good.trem_south_lat.min()) - 0.3,
        max(good.slip_north_lat.max(), good.trem_north_lat.max()) + 0.3]
ax3.plot(lims, lims, "k:", lw=1)
ax3.set_xlim(lims)
ax3.set_ylim(lims)
ax3.set_xlabel("modeled-slip termination latitude")
ax3.set_ylabel("tremor termination latitude")
ax3.set_title(f"tremor vs slip terminations\nmedian |diff| N {good.d_north_km.abs().median():.0f} km, "
              f"S {good.d_south_km.abs().median():.0f} km", fontsize=9)
ax3.legend(fontsize=7)
ax3.grid(alpha=0.3)

ax4 = fig2.add_subplot(gs2[1, 2:4])
for lab, gg, c in (("north-propagating", ev_np, "tab:red"),
                   ("south-propagating", ev_sp, "tab:blue")):
    la = np.concatenate([gg.trem_south_lat.values, gg.trem_north_lat.values])
    if la.size:
        h, _ = np.histogram(la, bins=edges)
        ax4.step(centres, h, where="mid", color=c, lw=1.7, label=f"{lab} (n {la.size})")
for i in SIG:
    ax4.axvline(centres[i], color="0.4", lw=0.9, ls="--", alpha=0.7)
ax4.set_xlabel("latitude of termination")
ax4.set_ylabel("count")
ax4.set_title("terminations by propagation direction", fontsize=9)
ax4.legend(fontsize=7)
ax4.grid(alpha=0.3)

ax5 = fig2.add_subplot(gs2[1, 4:6])
for _, r in good.iterrows():
    c = {"north": "tab:red", "south": "tab:blue"}.get(r["direction"], "0.6")
    yv = pd.Timestamp(r["start"]).year + pd.Timestamp(r["start"]).dayofyear / 365.25
    ax5.plot([r["trem_south_lat"], r["trem_north_lat"]], [yv, yv], "-", color=c, lw=2.0, alpha=0.85)
    ax5.plot([r["slip_south_lat"], r["slip_north_lat"]], [yv + 0.12, yv + 0.12], "-",
             color="0.35", lw=1.1, alpha=0.8)
for i in SIG:
    ax5.axvline(centres[i], color="k", lw=0.8, ls="--", alpha=0.6)
ax5.set_xlabel("latitude")
ax5.set_ylabel("event start (year)")
ax5.set_title("extents through time: tremor (colour = propagation)\nvs modeled slip (grey, offset)",
              fontsize=9)
ax5.legend(loc="upper left", handles=[Line2D([], [], color="tab:red", lw=2, label="tremor, N-propagating"),
                    Line2D([], [], color="tab:blue", lw=2, label="tremor, S-propagating"),
                    Line2D([], [], color="0.6", lw=2, label="tremor, unclear"),
                    Line2D([], [], color="0.35", lw=1.1, label="modeled slip")], fontsize=6.5)
ax5.grid(alpha=0.3)

fig2.suptitle("S7 Part B — along-strike terminations of tremor and modeled slip, "
              "27 confirmed events 2010-2026", fontsize=13)
fig2.savefig(OUT / "figures/s7_terminations.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
log("[write] figures/s7_terminations.png")

print("\n" + "=" * 100)
print("CAVEATS (printed with the numbers on purpose)")
print("=" * 100)
print("""\
1. Rank normalisation removes a segment-wide monotone change in detection sensitivity, NOT a
   spatially heterogeneous one WITHIN a segment. S3b measured exactly such heterogeneity
   (detections per tremor-day fell x0.79 north and x0.69 central), so within-segment sensitivity
   drift is the residual risk for any class CHANGE; class PERSISTENCE is the robust part.
2. Ties: many triangles have D3_nb = 0 in an epoch and share the average rank, so a large tie
   block can sit above the 25th percentile and be scored MIXED rather than QUIET. The tie
   fractions are printed above; they make QUIET conservative (an undercount).
3. Slip is 13 fixed broad vbICA patterns on a 0.25-deg fallback mesh (~28 km triangles); adjacent
   triangles are NOT independent, and only the POSITIVE daily rate is accumulated, so S rectifies
   model noise upward. The south is at the model noise floor (S6) and its classes should not be
   read as fault behaviour.
4. Tremor depth is the Slab2 interface depth at the epicentre, not a measured hypocentre.
5. Part B: footprints are large (lat span 1.1-5.8 deg) and chained, so an 'event' extent is partly
   a catalogue-segmentation choice; the tremor p5/p95 is a robust but blunt edge estimator.
6. Part B null: with only 27 events and extents that fill much of their segment, the repositioning
   freedom is small (printed above), so the barrier test has limited power; the 'strict' domain has
   essentially none. Bonferroni over all bins is deliberately harsh.
7. Slip extents come from a 3-mm accumulated-slip threshold on a smooth 13-pattern model, so they
   are inherently smoother and wider than the tremor extents; the tremor-minus-slip differences
   measure that smoothness as much as any physical offset.""")
log("done")
