#!/usr/bin/env python
"""S8 — (A) did the three largest northern ETS events change what came after them at their site?
       (B) is ETS recurrence in Cascadia time-predictable, slip-predictable, or neither?

PLAN.md Round 2, step S8 (questions Q7 and Q18).  Model + PNSN tremor only; no dv/v.

======================================================================================
PART A — BEFORE vs AFTER the three largest northern/central confirmed events
======================================================================================
REFERENCE EVENTS (the three largest confirmed north/central patches of results/s1_events.csv):
    id 1317  2015-12-28 .. 2016-03-21  Mw 6.79  north
    id 1863  2017-02-24 .. 2017-05-04  Mw 6.63  north
    id 1605  2024-08-14 .. 2024-10-01  Mw 6.58  north

SITE of a reference event R = the triangles of its label footprint
    F_R = {triangles whose label == R.id on at least one day in [ts_R, te_R]}   (npz `label`, 742 cols)

SITE MEMBERSHIP.  A confirmed event E belongs to R's site when
    overlap(E, R) = |F_E & F_R| / min(|F_E|, |F_R|)  >=  OVERLAP_THRESH.
    Primary threshold 0.30 (as specified); 0.60 is also run as a sensitivity and written
    to the same CSV with column `ov_thresh`.
SPLIT.  BEFORE = E.start <  R.start ;  AFTER = E.start > R.end.  (R itself and anything
    starting inside [R.start, R.end] fall in neither, so R is excluded automatically.)

QUANTITIES compared BEFORE vs AFTER (one value per event E):
  Q1_downdip_km            s1 Q1: (tremor centroid - slip centroid) . down-dip unit vector, km
                           (+ = tremor down-dip of the slip centroid).  Read from s1_events.csv.
  Q3_depth_km_raw          s1 Q3: median Slab2 depth of the DUR tremor - slip-weighted depth, km.
                           Carried for context only; S1b showed it is largely set by WHICH DAYS the
                           5 mm/30 d threshold selected, not by local geometry.
  Q3_rel_marginwide_km     THE depth quantity S1b tells us to use: the DUR tremor depth near the
                           patch MINUS the concurrent margin-wide tremor depth,
                             = median Slab2 depth of detections within 40 km of F during [start, end]
                             - median Slab2 depth of ALL detections with lat >= 46N (north events)
                               or 44N <= lat < 46N (central events) in the SAME [start, end] window.
                           This removes the window-phase effect (S1b sec. 4b) that dominates raw Q3.
  tremordays3_per_potency  completeness-hardened efficiency: the number of DAYS in [start, end] with
                           >= 3 detections within 40 km of F, divided by the event potency (km^3).
                           Days saturate (a day is one day whether the network reports 3 or 30
                           detections), so this is far less sensitive to catalogue growth than the
                           S1 detection-count Q5 (see S3b).  Units: tremor-days per km^3.
  potency_km3              s1 potency: sum over F of area * accumulated slip, km^3.
  dur_d                    event duration (trailing-30-d threshold crossings), days.
  recur_prev_d             recurrence interval at the site: E.start minus the start of the previous
                           event of the SAME site chain (the site's events incl. R, sorted by start).
                           NaN for the first event of the chain.

TEST.  Mann-Whitney U (two-sided) BEFORE vs AFTER, per reference event and for the STACKED pool
    (each event contributes once per reference site it belongs to, so an event in 2 sites appears
    twice; this is the specified pooling).  n and medians reported for both groups.

NULL.  200 random reference dates drawn uniformly in 2012-01-01..2024-12-31 (one independent date
    per site per draw).  The SITE membership is held fixed at the real one; only the SPLIT DATE
    moves, and the same gap is applied (BEFORE: start < d ; AFTER: start > d + dur_R).  This is the
    null "the split date carries no information": it preserves the site, the event list, the values,
    the record length and the group-size imbalance, and destroys only the association with the big
    event.  Reported per quantity: sd of the null BEFORE-AFTER median difference, z = diff / sd,
    and a two-sided null p = (1 + #{|diff_null| >= |diff_obs|}) / (1 + n_valid_draws).
    Draws with fewer than 2 events on either side are discarded (n_null_ok is reported).

======================================================================================
PART B — time-predictable vs slip-predictable recurrence
======================================================================================
SITES = the 0.5-deg latitude bands of S5, restricted to 44.0 <= band < 50.0 (north + central) with
    >= 6 tremor episodes.  All 12 such bands qualify (13-20 episodes each).
EPISODES = results/s5_tremor_episodes.csv (S5 definition: local maxima of the 30-d trailing detection
    count >= the band's 90th percentile, greedy highest-first, >= 90 d apart; onset = first day in
    [peak-45, peak] above 25 % of the peak).  Recurrence uses PEAK dates, as in S5.
    For episode i of a band (sorted by peak):
        T_prev[i] = peak[i]   - peak[i-1]   (days; NaN for the first)
        T_next[i] = peak[i+1] - peak[i]     (days; NaN for the last)

SIZE PROXIES (per episode)
  count30_raw          the S5 `count30` = detections in the 30 d ending at the peak.
  count30_detr         count30 divided by a 5-YEAR RUNNING MEDIAN of the band's own count30
                       (median over the episodes of that band with peak within +- 2.5 yr, the
                       episode included; if fewer than 3 such episodes, the band median is used).
                       This removes the PNSN network-growth ramp that inflates late episodes.
  patch_Mw             Mw of the slip patch matched to the episode by S5 (episode -> nearest patch
                       start within [onset-30, onset+90]); NaN when unmatched.
  patch_potency_km3    potency of that patch (slip_patches_5mm.csv), NaN when unmatched.

MODELS
  TIME-PREDICTABLE:  size ~ T_prev   (a longer wait accumulates a bigger event)
  SLIP-PREDICTABLE:  size ~ T_next   (a bigger event needs longer to reload)
  Spearman rho per site; p from 2000 permutations of the SIZE values WITHIN that site.
  POOLED (and per segment): each site's size and T are replaced by their within-site fractional
  ranks (rank - 0.5)/n, then pooled and correlated; the permutation shuffles sizes within each site
  independently, 2000 draws.

ARTEFACT CONTROLS on the predictability result (added after the first run produced a significant
  pooled size ~ T_next correlation; a positive rho there can be manufactured by the episode
  DETECTOR rather than by the fault, because small episodes tend to be the closely-spaced secondary
  peaks of one ETS complex):
  (C1) INTERVAL CLUSTERING.  Spearman rho(T_prev, T_next) per band and pooled.  If short intervals
       follow short intervals, size ~ T_prev alone can induce an apparent size ~ T_next.
  (C2) PARTIAL correlation.  Rank-partial Spearman of size vs T_next CONTROLLING for T_prev
       (and of size vs T_prev controlling for T_next), same within-site permutation for p.
       A slip-predictable signal must survive removing the T_prev channel.
  (C3) DE-CLUSTERING SENSITIVITY.  Repeat the pooled test using only pairs with the interval
       >= 150 d, i.e. dropping the closely-spaced pairs that a stricter minimum-separation rule
       in the episode detector would have merged.  Rows appear as site = "*_sep150".
  (C4) RE-DETECT.  S5 picks peaks GREEDILY highest-first with a >= 90 d exclusion, so a large peak
       wins its slot and suppresses any neighbour within 90 d — which mechanically stretches the
       interval AFTER a large peak.  The strength of that mechanism scales with the exclusion
       radius, so the episode catalogue is rebuilt from the band trailing counts at PK_SEP =
       45 / 90 / 150 d (the 90 d rebuild is asserted to reproduce s5_tremor_episodes.csv exactly)
       and the pooled correlations recomputed.  Flat or DECREASING rho with PK_SEP clears it.
  (C5) IS rho_next BIGGER THAN rho_prev?  Bootstrap over the 12 bands (2000 draws) for the
       difference rho_next - rho_prev, with its 95 % CI.

NETWORK-GROWTH CHECK.  Spearman of count30 vs peak date, per band and pooled (per-band z-scores),
  reported before the predictability test; the detrended proxy is the control.

POWER.  For each n, the |rho| detectable with 80 % power at alpha = 0.05 two-sided,
  rho_80 = tanh((1.96 + 0.8416) / sqrt(n - 3)).  Reported per row.

======================================================================================
INPUTS (read-only, outside this folder except the two S-step tables)
======================================================================================
  tremor_slip_geometry/results/s1_events.csv           27 confirmed north/central events (S1)
  tremor_slip_geometry/results/s5_tremor_episodes.csv  359 tremor episodes (S5)
  tremor_slip_geometry/results/s5_matches.csv          501 episode-patch pairs (S5; consistency check)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.csv   patch Mw / potency
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz   label (T x 742), days, patch_lat/lon
  catalogs/pnsn_tremor_cascadia_full.csv               PNSN tremor (time, lat, lon)
  /home/jovyan/cascadia-sse-strain                     mesh + Slab2 depth interpolator
                                                       (slab_depth_interpolator only; the slip
                                                        solution is NOT needed here)

OUTPUTS (all inside tremor_slip_geometry/)
  results/s8_bigevents_before_after.csv   ref (1317/1863/1605/pooled) x ov_thresh x quantity:
                                          n_before, n_after, med_before, med_after, diff, p_mwu,
                                          null_sd, z, p_null, n_null_ok
  results/s8_bigevents_events.csv         the per-event input table actually used (site membership,
                                          side, and all seven quantities) — provenance for the above
  results/s8_predictability.csv           site (band, pooled, segment, *_sep150) x proxy x n:
                                          rho_prev, p_prev, rho_next, p_next, rho_80, and the
                                          rank-partial versions rho_prev_partial / rho_next_partial
  results/s8_predictability_controls.csv  the C1 interval-clustering diagnostic (per band + pooled),
                                          the C4 re-detection at PK_SEP 45/90/150 d, and the C5
                                          bootstrap CI on rho_next - rho_prev (column `control`)
  results/s8_growth_check.csv             per band: Spearman count30 vs date, raw and detrended
  figures/s8_bigevents.png                strip plots BEFORE/AFTER per quantity, pooled set
  figures/s8_predictability.png           rank-normalised size vs T_prev and vs T_next, N/C colours,
                                          plus the growth check, a per-site rho forest and the power curve

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s8_bigevents_predictability.py
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
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, rankdata, spearmanr

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260910)
TREMOR_RADIUS_KM = 40.0
COS47 = np.cos(np.deg2rad(47.0))
OV_THRESHES = [0.30, 0.60]          # 0.30 is the specified primary
N_NULL = 200                        # random reference dates
NULL_LO = np.datetime64("2012-01-01")
NULL_HI = np.datetime64("2024-12-31")
MIN_SIDE = 2                        # minimum events per side for a null draw to count
N_PERM = 2000                       # permutations for Part B
MIN_EPISODES = 6
DETREND_HALFWIN_D = 2.5 * 365.25    # 5-yr running median half-window
REFS = [1317, 1863, 1605]

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def rho80(n: int) -> float:
    """|rho| detectable with 80 % power, alpha = 0.05 two-sided (Fisher z)."""
    if n is None or n < 5:
        return np.nan
    return float(np.tanh((1.959964 + 0.841621) / np.sqrt(n - 3)))


def perm_p(obs: float, null: np.ndarray) -> float:
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return np.nan
    return float((1 + int((np.abs(null) >= abs(obs) - 1e-12).sum())) / (1 + null.size))


print("=" * 100)
print("S8 — big-event before/after (Q7) and time- vs slip-predictable recurrence (Q18)")
print("=" * 100)

# ======================================================================================
# 0. Shared data
# ======================================================================================
log("loading mesh + Slab2 interpolator ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
slabf = slab_depth_interpolator(fault)

NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]
DAYS = np.asarray(NPZ["days"], dtype=str)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PX, PY = to_km(PLON, PLAT)
log(f"  label {LABEL.shape}, days {DAYS[0]}..{DAYS[-1]}, {len(PLAT)} triangles")

log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_LAT = tr["lat"].values
TR_DEPTH = slabf(tr["lon"].values, tr["lat"].values)
TR_DATE = tr["t"].values.astype("datetime64[D]")
TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(
    f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}; "
    f"Slab2 depth finite for {np.isfinite(TR_DEPTH).mean() * 100:.1f}%"
)

EV = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  s1_events {len(EV)} rows ({EV.seg.value_counts().to_dict()})")

PATCH = pd.read_csv(SSE / "slip_patches_5mm.csv")

# ======================================================================================
# PART A
# ======================================================================================
print("\n" + "=" * 100)
print("PART A — before/after the three largest northern events (Q7)")
print("=" * 100)

# ---- footprints -----------------------------------------------------------------
FP = {}
for r in EV.itertuples():
    ts, te = DAY_INDEX[r.start], DAY_INDEX[r.end]
    FP[int(r.id)] = np.any(LABEL[ts : te + 1] == int(r.id), axis=0)

# ---- per-event quantities that must be recomputed from the catalogues ------------
log("computing per-event tremor quantities (near-patch depth, margin-wide depth, tremor-days) ...")
rows = []
for r in EV.itertuples():
    pid = int(r.id)
    F = FP[pid]
    d0 = np.datetime64(r.start)
    d1 = np.datetime64(r.end)

    idx = TREE.query_ball_point(np.column_stack([PX[F], PY[F]]), r=TREMOR_RADIUS_KM)
    near = (
        np.unique(np.concatenate([np.asarray(i, int) for i in idx]))
        if len(idx)
        else np.array([], int)
    )
    inwin = near[(TR_DATE[near] >= d0) & (TR_DATE[near] <= d1)]

    dz = TR_DEPTH[inwin]
    dz = dz[np.isfinite(dz)]
    near_med = float(np.median(dz)) if dz.size else np.nan
    near_set = set(inwin[np.isfinite(TR_DEPTH[inwin])].tolist())

    # concurrent margin-wide tremor depth in the same window
    if r.seg == "north":
        seg_mask = TR_LAT >= 46.0
        seg_lab = "lat>=46N"
    else:
        seg_mask = (TR_LAT >= 44.0) & (TR_LAT < 46.0)
        seg_lab = "44-46N"
    mw_sel = seg_mask & (TR_DATE >= d0) & (TR_DATE <= d1)
    mw_idx = np.flatnonzero(mw_sel & np.isfinite(TR_DEPTH))
    mdz = TR_DEPTH[mw_idx]
    margin_med = float(np.median(mdz)) if mdz.size else np.nan
    # how self-referential is the margin-wide reference?  (S1b: for the big northern footprints the
    # 40-km near-patch set IS essentially all northern tremor, so Q3_rel is near-zero by construction)
    shared = len(near_set.intersection(mw_idx.tolist())) if mw_idx.size else 0
    frac_shared = shared / mw_idx.size if mw_idx.size else np.nan

    # tremor-days with >= 3 detections in the footprint neighbourhood during DUR
    if inwin.size:
        udays, cnt = np.unique(TR_DATE[inwin], return_counts=True)
        d3 = int((cnt >= 3).sum())
        d1cnt = int(udays.size)
    else:
        d3 = d1cnt = 0

    rows.append(
        dict(
            id=pid,
            start=r.start,
            end=r.end,
            seg=r.seg,
            Mw=r.Mw,
            dur_d=r.dur_d,
            potency_km3=r.potency_km3,
            N_s1=r.N,
            n_tri_F=int(F.sum()),
            Q1_downdip_km=r.Q1_downdip_km,
            Q3_depth_km_raw=r.Q3_depth_km,
            tremor_depth_med_recomp=near_med,
            tremor_depth_med_s1=r.tremor_depth_med,
            marginwide_depth_med_km=margin_med,
            margin_ref=seg_lab,
            n_near_win=int(len(near_set)),
            n_marginwide=int(mdz.size),
            frac_marginwide_in_near=frac_shared,
            Q3_rel_marginwide_km=near_med - margin_med,
            tremordays1=d1cnt,
            tremordays3=d3,
            tremordays3_per_potency=d3 / r.potency_km3 if r.potency_km3 > 0 else np.nan,
        )
    )
Q = pd.DataFrame(rows)

# integrity check: our near-patch tremor median must reproduce S1's
dd = np.abs(Q.tremor_depth_med_recomp - Q.tremor_depth_med_s1)
log(f"  check: |recomputed - s1 tremor_depth_med| max = {np.nanmax(dd):.3e} km  (must be ~0)")
assert np.nanmax(dd) < 1e-6, "near-patch tremor depth does not reproduce s1_events.csv"

QUANTS = [
    ("Q1_downdip_km", "Q1 down-dip offset (km)"),
    ("Q3_depth_km_raw", "Q3 raw depth offset (km)"),
    ("Q3_rel_marginwide_km", "tremor depth - margin-wide (km)"),
    ("tremordays3_per_potency", "tremor-days(>=3) / potency (d/km3)"),
    ("potency_km3", "potency (km3)"),
    ("dur_d", "duration (d)"),
    ("recur_prev_d", "recurrence to previous (d)"),
]

print("\nPer-event quantities used in Part A (27 confirmed north/central events):")
show = Q[
    [
        "id",
        "start",
        "seg",
        "Mw",
        "dur_d",
        "potency_km3",
        "Q1_downdip_km",
        "Q3_depth_km_raw",
        "tremor_depth_med_recomp",
        "marginwide_depth_med_km",
        "Q3_rel_marginwide_km",
        "frac_marginwide_in_near",
        "tremordays3",
        "tremordays3_per_potency",
    ]
].copy()
with pd.option_context("display.width", 240, "display.max_columns", 40):
    print(show.to_string(index=False, float_format=lambda v: f"{v:9.4f}"))
print(
    "  NOTE `frac_marginwide_in_near` = fraction of the concurrent margin-wide detections that are\n"
    "  ALREADY inside the 40-km near-patch set.  Where it is near 1 the margin-wide reference is\n"
    "  self-referential and Q3_rel_marginwide_km is ~0 BY CONSTRUCTION (S1b sec. 4b/5)."
)
print(
    f"  median frac_marginwide_in_near: north {Q.loc[Q.seg=='north','frac_marginwide_in_near'].median():.3f}, "
    f"central {Q.loc[Q.seg=='central','frac_marginwide_in_near'].median():.3f}"
)

# ---- sites, sides, recurrence ----------------------------------------------------
REFROW = {int(r.id): r for r in EV.itertuples() if int(r.id) in REFS}
print("\nReference events:")
for rid in REFS:
    r = REFROW[rid]
    print(
        f"   {rid}  {r.start} .. {r.end}  Mw {r.Mw:.2f}  {r.seg}  "
        f"footprint {int(FP[rid].sum())} triangles  potency {r.potency_km3:.3f} km3"
    )

ev_rows = []
site_members = {}  # (rid, thresh) -> DataFrame of the site chain with recurrence
for thr in OV_THRESHES:
    for rid in REFS:
        rref = REFROW[rid]
        fr = FP[rid]
        rs, re_ = np.datetime64(rref.start), np.datetime64(rref.end)
        mem = []
        for r in EV.itertuples():
            f = FP[int(r.id)]
            ov = float((f & fr).sum()) / float(min(f.sum(), fr.sum()))
            if ov < thr:
                continue
            st = np.datetime64(r.start)
            side = "REF" if int(r.id) == rid else ("before" if st < rs else ("after" if st > re_ else "during"))
            mem.append(dict(id=int(r.id), start=r.start, ov=ov, side=side))
        m = pd.DataFrame(mem).sort_values("start").reset_index(drop=True)
        # recurrence within the site chain (previous event of the same site)
        st = pd.to_datetime(m["start"]).values.astype("datetime64[D]")
        rec = np.full(len(m), np.nan)
        rec[1:] = (st[1:] - st[:-1]).astype(int)
        m["recur_prev_d"] = rec
        site_members[(rid, thr)] = m
        for mm in m.itertuples():
            q = Q[Q.id == mm.id].iloc[0]
            row = dict(
                ref=rid,
                ov_thresh=thr,
                id=mm.id,
                start=mm.start,
                seg=q.seg,
                overlap=mm.ov,
                side=mm.side,
                recur_prev_d=mm.recur_prev_d,
            )
            for k, _ in QUANTS:
                if k != "recur_prev_d":
                    row[k] = q[k]
            ev_rows.append(row)
EVA = pd.DataFrame(ev_rows)
EVA.to_csv(OUT_RES / "s8_bigevents_events.csv", index=False)
log(f"  wrote results/s8_bigevents_events.csv ({len(EVA)} rows)")

for thr in OV_THRESHES:
    print(f"\n--- site membership at overlap >= {thr:.2f} ---")
    for rid in REFS:
        m = site_members[(rid, thr)]
        nb = int((m.side == "before").sum())
        na = int((m.side == "after").sum())
        print(f"   ref {rid}: site has {len(m)} confirmed events -> BEFORE {nb}, AFTER {na}")
        print(
            "      "
            + ", ".join(
                f"{r.id}({r.start[:7]},{r.side[0].upper()},ov{r.ov:.2f})" for r in m.itertuples()
            )
        )

# ---- observed + null tests --------------------------------------------------------
log("running Part A tests + 200-draw random-reference-date null ...")

# pre-draw the null dates: independent per site, shared across thresholds and quantities
null_dates = {rid: RNG.integers(0, int((NULL_HI - NULL_LO).astype(int)) + 1, N_NULL) for rid in REFS}


def split_stats(vals_before, vals_after):
    b = np.asarray(vals_before, float)
    a = np.asarray(vals_after, float)
    b = b[np.isfinite(b)]
    a = a[np.isfinite(a)]
    if b.size == 0 or a.size == 0:
        return b.size, a.size, np.nan, np.nan, np.nan, np.nan
    mb, ma = float(np.median(b)), float(np.median(a))
    try:
        p = float(mannwhitneyu(b, a, alternative="two-sided").pvalue)
    except ValueError:
        p = np.nan
    return b.size, a.size, mb, ma, ma - mb, p


out_rows = []
for thr in OV_THRESHES:
    # --- per reference and pooled observed sets
    obs_sets = {}
    for rid in REFS:
        sub = EVA[(EVA.ref == rid) & (EVA.ov_thresh == thr)]
        obs_sets[rid] = sub
    pooled = EVA[EVA.ov_thresh == thr]

    # --- null draws
    null_diffs = {rid: {k: [] for k, _ in QUANTS} for rid in REFS}
    null_diffs["pooled"] = {k: [] for k, _ in QUANTS}
    for j in range(N_NULL):
        per_ref_parts = {k: ([], []) for k, _ in QUANTS}
        ok_pool = True
        for rid in REFS:
            rref = REFROW[rid]
            dur = int(rref.dur_d)
            d = NULL_LO + int(null_dates[rid][j])
            sub = obs_sets[rid]
            st = pd.to_datetime(sub["start"]).values.astype("datetime64[D]")
            bmask = st < d
            amask = st > (d + np.timedelta64(dur, "D"))
            for k, _ in QUANTS:
                v = sub[k].to_numpy(float)
                vb, va = v[bmask], v[amask]
                vb, va = vb[np.isfinite(vb)], va[np.isfinite(va)]
                if vb.size >= MIN_SIDE and va.size >= MIN_SIDE:
                    null_diffs[rid][k].append(float(np.median(va) - np.median(vb)))
                else:
                    null_diffs[rid][k].append(np.nan)
                per_ref_parts[k][0].append(vb)
                per_ref_parts[k][1].append(va)
        for k, _ in QUANTS:
            vb = np.concatenate(per_ref_parts[k][0]) if per_ref_parts[k][0] else np.array([])
            va = np.concatenate(per_ref_parts[k][1]) if per_ref_parts[k][1] else np.array([])
            if vb.size >= MIN_SIDE and va.size >= MIN_SIDE:
                null_diffs["pooled"][k].append(float(np.median(va) - np.median(vb)))
            else:
                null_diffs["pooled"][k].append(np.nan)

    for label, sub in list((str(rid), obs_sets[rid]) for rid in REFS) + [("pooled", pooled)]:
        key = int(label) if label != "pooled" else "pooled"
        for k, kl in QUANTS:
            vb = sub.loc[sub.side == "before", k].to_numpy(float)
            va = sub.loc[sub.side == "after", k].to_numpy(float)
            nb, na, mb, ma, diff, p = split_stats(vb, va)
            nul = np.asarray(null_diffs[key][k], float)
            nul = nul[np.isfinite(nul)]
            sd = float(np.std(nul, ddof=1)) if nul.size > 2 else np.nan
            z = diff / sd if (np.isfinite(diff) and np.isfinite(sd) and sd > 0) else np.nan
            out_rows.append(
                dict(
                    ref=label,
                    ov_thresh=thr,
                    quantity=k,
                    label=kl,
                    n_before=nb,
                    n_after=na,
                    med_before=mb,
                    med_after=ma,
                    diff=diff,
                    p_mwu=p,
                    null_sd=sd,
                    z=z,
                    p_null=perm_p(diff, nul),
                    n_null_ok=int(nul.size),
                )
            )
BA = pd.DataFrame(out_rows)
BA.to_csv(OUT_RES / "s8_bigevents_before_after.csv", index=False)
log(f"  wrote results/s8_bigevents_before_after.csv ({len(BA)} rows)")

for thr in OV_THRESHES:
    tag = "PRIMARY" if thr == 0.30 else "sensitivity"
    print(f"\n--- PART A results, overlap >= {thr:.2f} ({tag}) ---")
    print(
        f"{'ref':>7} {'quantity':<24} {'nB':>3} {'nA':>3} {'med_before':>11} {'med_after':>11} "
        f"{'diff':>10} {'p_MWU':>8} {'null_sd':>9} {'z':>7} {'p_null':>7}"
    )
    for r in BA[BA.ov_thresh == thr].itertuples():
        print(
            f"{r.ref:>7} {r.quantity:<24} {r.n_before:3d} {r.n_after:3d} "
            f"{r.med_before:11.4f} {r.med_after:11.4f} {r.diff:10.4f} "
            f"{r.p_mwu:8.4f} {r.null_sd:9.4f} {r.z:7.2f} {r.p_null:7.4f}"
        )
print(
    "\n  READ p_null, NOT z.  The null distribution of a DIFFERENCE OF MEDIANS over 5-18 events is\n"
    "  discrete and heavy-tailed (the median jumps between data values), so its sd underestimates\n"
    "  the tail and |z| is inflated — e.g. Q3_depth_km_raw shows |z| ~ 4 with p_null ~ 0.2-0.3.\n"
    "  p_null is the rank-based two-sided p against the same 200 random-reference-date draws."
)

# ======================================================================================
# PART B
# ======================================================================================
print("\n" + "=" * 100)
print("PART B — time- vs slip-predictable recurrence (Q18)")
print("=" * 100)

E = pd.read_csv(OUT_RES / "s5_tremor_episodes.csv")
E["peak"] = pd.to_datetime(E["peak"])
E["onset"] = pd.to_datetime(E["onset"])
E = E[(E.band >= 44.0) & (E.band < 50.0)].copy()
cnt = E.groupby("band").size()
keep = cnt[cnt >= MIN_EPISODES].index
E = E[E.band.isin(keep)].sort_values(["band", "peak"]).reset_index(drop=True)
log(f"  {len(E)} episodes in {E.band.nunique()} bands (44-50N, >= {MIN_EPISODES} episodes)")

# consistency with s5_matches (episode->patch direction)
MM = pd.read_csv(OUT_RES / "s5_matches.csv")
mfe = MM[MM.from_episode].set_index("episode_key")["patch_id"]
chk = E.set_index("key")["patch_id"].dropna()
common = chk.index.intersection(mfe.index)
agree = int((chk.loc[common].astype(int) == mfe.loc[common].astype(int)).sum())
log(f"  episode->patch consistency with s5_matches: {agree}/{len(common)} agree")

PMW = PATCH.set_index("id")["Mw"]
PPOT = PATCH.set_index("id")["potency_km3"]
E["patch_Mw"] = [PMW.get(int(p), np.nan) if np.isfinite(p) else np.nan for p in E["patch_id"]]
E["patch_potency_km3"] = [
    PPOT.get(int(p), np.nan) if np.isfinite(p) else np.nan for p in E["patch_id"]
]

# --- recurrence intervals and the 5-yr running-median detrend ----------------------
E["T_prev"] = np.nan
E["T_next"] = np.nan
E["count30_detr"] = np.nan
E["run_med5"] = np.nan
for b, g in E.groupby("band"):
    g = g.sort_values("peak")
    ii = g.index.to_numpy()
    pk = g["peak"].values.astype("datetime64[D]").astype(int)
    E.loc[ii[1:], "T_prev"] = (pk[1:] - pk[:-1]).astype(float)
    E.loc[ii[:-1], "T_next"] = (pk[1:] - pk[:-1]).astype(float)
    c = g["count30"].to_numpy(float)
    bandmed = float(np.median(c))
    for a, (i, t) in enumerate(zip(ii, pk)):
        w = np.abs(pk - t) <= DETREND_HALFWIN_D
        m = float(np.median(c[w])) if w.sum() >= 3 else bandmed
        E.loc[i, "run_med5"] = m
        E.loc[i, "count30_detr"] = c[a] / m if m > 0 else np.nan

PROXIES = [
    ("count30_raw", "count30", "tremor count at peak (raw)"),
    ("count30_detr", "count30_detr", "tremor count at peak (5-yr detrended)"),
    ("patch_Mw", "patch_Mw", "matched slip-patch Mw"),
    ("patch_potency_km3", "patch_potency_km3", "matched slip-patch potency"),
]
E["count30_raw"] = E["count30"]


# --- network-growth check ----------------------------------------------------------
grow = []
for b, g in E.groupby("band"):
    x = g["peak"].values.astype("datetime64[D]").astype(float)
    for col in ["count30_raw", "count30_detr"]:
        y = g[col].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        rho, p = spearmanr(x[ok], y[ok]) if ok.sum() >= 5 else (np.nan, np.nan)
        grow.append(dict(band=b, seg=g["seg"].iloc[0], proxy=col, n=int(ok.sum()), rho=rho, p=p))
GR = pd.DataFrame(grow)
# pooled growth: per-band z-scores of the proxy vs date
poolg = []
for col in ["count30_raw", "count30_detr"]:
    xs, ys, segs = [], [], []
    for b, g in E.groupby("band"):
        y = g[col].to_numpy(float)
        if np.nanstd(y) == 0 or not np.isfinite(np.nanstd(y)):
            continue
        xs.append(g["peak"].values.astype("datetime64[D]").astype(float))
        ys.append((y - np.nanmean(y)) / np.nanstd(y))
        segs += list(g["seg"])
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    ok = np.isfinite(x) & np.isfinite(y)
    rho, p = spearmanr(x[ok], y[ok])
    poolg.append(dict(band="pooled", seg="all", proxy=col, n=int(ok.sum()), rho=rho, p=p))
GR = pd.concat([GR, pd.DataFrame(poolg)], ignore_index=True)
GR.to_csv(OUT_RES / "s8_growth_check.csv", index=False)

print("\n--- network-growth check: episode size vs date ---")
print(f"{'band':>7} {'seg':<8} {'proxy':<14} {'n':>3} {'rho':>7} {'p':>9}")
for r in GR.itertuples():
    bl = r.band if isinstance(r.band, str) else f"{r.band:.1f}"
    print(f"{bl:>7} {r.seg:<8} {r.proxy:<14} {r.n:3d} {r.rho:7.3f} {r.p:9.4f}")


# --- predictability ----------------------------------------------------------------
def partial_rho(y, x, z):
    """Rank-partial correlation of y and x controlling for z (Spearman rhos in the formula)."""
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(z)
    if ok.sum() < 5:
        return np.nan
    ryx = spearmanr(y[ok], x[ok])[0]
    ryz = spearmanr(y[ok], z[ok])[0]
    rxz = spearmanr(x[ok], z[ok])[0]
    den = np.sqrt(max(1 - ryz**2, 1e-12) * max(1 - rxz**2, 1e-12))
    return float((ryx - ryz * rxz) / den)


def site_rho(g, proxy, tcol, ocol, rng):
    """Spearman rho (+ within-site permutation p) and the partial version controlling for `ocol`."""
    y = g[proxy].to_numpy(float)
    x = g[tcol].to_numpy(float)
    z = g[ocol].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 5:
        return np.nan, np.nan, n, np.nan, np.nan
    xx, yy = x[ok], y[ok]
    rho, _ = spearmanr(xx, yy)
    nul = np.empty(N_PERM)
    for j in range(N_PERM):
        nul[j] = spearmanr(xx, rng.permutation(yy))[0]
    # partial: needs all three finite
    ok3 = ok & np.isfinite(z)
    if ok3.sum() >= 5:
        pr = partial_rho(y[ok3], x[ok3], z[ok3])
        nulp = np.empty(N_PERM)
        y3, x3, z3 = y[ok3], x[ok3], z[ok3]
        for j in range(N_PERM):
            nulp[j] = partial_rho(rng.permutation(y3), x3, z3)
        ppr = perm_p(pr, nulp)
    else:
        pr, ppr = np.nan, np.nan
    return float(rho), perm_p(rho, nul), n, pr, ppr


def frank(v):
    """within-site fractional ranks (rank-0.5)/n on the finite entries; NaN elsewhere."""
    v = np.asarray(v, float)
    out = np.full(v.shape, np.nan)
    ok = np.isfinite(v)
    if ok.sum() >= 2:
        out[ok] = (rankdata(v[ok]) - 0.5) / ok.sum()
    return out


rows = []
rng = np.random.default_rng(20260910)

# per-site
for b, g in E.groupby("band"):
    for pname, pcol, _ in PROXIES:
        rp, pp, npv, prp, ppp = site_rho(g, pcol, "T_prev", "T_next", rng)
        rn, pn, nnv, prn, ppn = site_rho(g, pcol, "T_next", "T_prev", rng)
        rows.append(
            dict(
                site=f"{b:.1f}",
                seg=g["seg"].iloc[0],
                proxy=pname,
                n_prev=npv,
                rho_prev=rp,
                p_prev=pp,
                n_next=nnv,
                rho_next=rn,
                p_next=pn,
                rho_prev_partial=prp,
                p_prev_partial=ppp,
                rho_next_partial=prn,
                p_next_partial=ppn,
                rho80_prev=rho80(npv),
                rho80_next=rho80(nnv),
            )
        )


def pooled_test(bands, pcol, min_interval=None, keep_xy=False, df=None):
    """Pooled Spearman on within-site fractional ranks, with the within-site permutation p.

    min_interval: if set, only episode pairs whose interval is >= this many days are used
                  (control C3, a proxy for a stricter minimum-separation de-clustering rule).
    Also returns the rank-partial rho of size vs T controlling for the OTHER interval (C2).
    """
    src = E if df is None else df
    Xp, Xn, Y, SITEID = [], [], [], []
    for b in bands:
        g = src[src.band == b]
        tp = g["T_prev"].to_numpy(float).copy()
        tn = g["T_next"].to_numpy(float).copy()
        y = g[pcol].to_numpy(float).copy()
        if min_interval is not None:
            tp[~(tp >= min_interval)] = np.nan
            tn[~(tn >= min_interval)] = np.nan
        Xp.append(frank(tp))
        Xn.append(frank(tn))
        Y.append(frank(y))
        SITEID.append(np.full(len(g), b))
    Xp, Xn, Y, SITEID = (np.concatenate(a) for a in (Xp, Xn, Y, SITEID))
    res = {}
    for tag, X, O in [("prev", Xp, Xn), ("next", Xn, Xp)]:
        ok = np.isfinite(X) & np.isfinite(Y)
        n = int(ok.sum())
        if n < 5:
            res[tag] = dict(rho=np.nan, p=np.nan, n=n, rho_par=np.nan, p_par=np.nan, xy=None)
            continue
        xo, yo, so, oo = X[ok], Y[ok], SITEID[ok], O[ok]
        rho, _ = spearmanr(xo, yo)
        grp = [np.where(so == b)[0] for b in np.unique(so)]
        nul = np.empty(N_PERM)
        for j in range(N_PERM):
            yp = yo.copy()
            for gi in grp:
                yp[gi] = rng.permutation(yo[gi])
            nul[j] = spearmanr(xo, yp)[0]
        ok3 = np.isfinite(oo)
        if ok3.sum() >= 5:
            y3, x3, o3, s3 = yo[ok3], xo[ok3], oo[ok3], so[ok3]
            rpar = partial_rho(y3, x3, o3)
            grp3 = [np.where(s3 == b)[0] for b in np.unique(s3)]
            nulp = np.empty(N_PERM)
            for j in range(N_PERM):
                yp = y3.copy()
                for gi in grp3:
                    yp[gi] = rng.permutation(y3[gi])
                nulp[j] = partial_rho(yp, x3, o3)
            ppar = perm_p(rpar, nulp)
        else:
            rpar, ppar = np.nan, np.nan
        res[tag] = dict(
            rho=float(rho),
            p=perm_p(rho, nul),
            n=n,
            rho_par=rpar,
            p_par=ppar,
            xy=(xo, yo, so) if keep_xy else None,
        )
    return res


# pooled / per segment, on within-site fractional ranks (+ the sep>=150 d de-clustering control)
POOL = {}
BANDSETS = [
    ("pooled", sorted(E.band.unique())),
    ("north", sorted(E.loc[E.seg == "north", "band"].unique())),
    ("central", sorted(E.loc[E.seg == "central", "band"].unique())),
]
for suffix, minint in [("", None), ("_sep150", 150.0)]:
    for grpname, bands in BANDSETS:
        for pname, pcol, _ in PROXIES:
            res = pooled_test(bands, pcol, minint, keep_xy=(suffix == "" and grpname == "pooled"))
            rows.append(
                dict(
                    site=grpname + suffix,
                    seg=grpname,
                    proxy=pname,
                    n_prev=res["prev"]["n"],
                    rho_prev=res["prev"]["rho"],
                    p_prev=res["prev"]["p"],
                    n_next=res["next"]["n"],
                    rho_next=res["next"]["rho"],
                    p_next=res["next"]["p"],
                    rho_prev_partial=res["prev"]["rho_par"],
                    p_prev_partial=res["prev"]["p_par"],
                    rho_next_partial=res["next"]["rho_par"],
                    p_next_partial=res["next"]["p_par"],
                    rho80_prev=rho80(res["prev"]["n"]),
                    rho80_next=rho80(res["next"]["n"]),
                )
            )
            if suffix == "" and grpname == "pooled":
                POOL[pname] = {k: (v["rho"], v["p"], v["n"], v["xy"]) for k, v in res.items()}

PR = pd.DataFrame(rows)
PR.to_csv(OUT_RES / "s8_predictability.csv", index=False)
log(f"  wrote results/s8_predictability.csv ({len(PR)} rows)")

print("\n--- PART B: time-predictable (size~T_prev) vs slip-predictable (size~T_next) ---")
print(
    f"{'site':>14} {'seg':<8} {'proxy':<18} {'n_pr':>4} {'rho_prev':>9} {'p_prev':>8} "
    f"{'n_nx':>4} {'rho_next':>9} {'p_next':>8} {'rho80_pr':>8} | "
    f"{'rho_pr|nx':>9} {'p':>7} {'rho_nx|pr':>9} {'p':>7}"
)
for r in PR.itertuples():
    print(
        f"{r.site:>14} {r.seg:<8} {r.proxy:<18} {r.n_prev:4d} {r.rho_prev:9.3f} {r.p_prev:8.4f} "
        f"{r.n_next:4d} {r.rho_next:9.3f} {r.p_next:8.4f} {r.rho80_prev:8.3f} | "
        f"{r.rho_prev_partial:9.3f} {r.p_prev_partial:7.3f} "
        f"{r.rho_next_partial:9.3f} {r.p_next_partial:7.3f}"
    )

# --- C1: interval clustering (does a short interval follow a short interval?) ---------
crows = []
for b, g in E.groupby("band"):
    tp = g["T_prev"].to_numpy(float)
    tn = g["T_next"].to_numpy(float)
    ok = np.isfinite(tp) & np.isfinite(tn)
    rho, p = spearmanr(tp[ok], tn[ok]) if ok.sum() >= 5 else (np.nan, np.nan)
    crows.append(dict(control="C1_interval_clustering", site=f"{b:.1f}", seg=g["seg"].iloc[0],
                      n=int(ok.sum()), rho_Tprev_Tnext=rho, p=p))
Xp, Xn, SID = [], [], []
for b in sorted(E.band.unique()):
    g = E[E.band == b]
    Xp.append(frank(g["T_prev"].to_numpy(float)))
    Xn.append(frank(g["T_next"].to_numpy(float)))
    SID.append(np.full(len(g), b))
Xp, Xn, SID = np.concatenate(Xp), np.concatenate(Xn), np.concatenate(SID)
ok = np.isfinite(Xp) & np.isfinite(Xn)
rho, _ = spearmanr(Xp[ok], Xn[ok])
grp = [np.where(SID[ok] == b)[0] for b in np.unique(SID[ok])]
nul = np.empty(N_PERM)
xo, yo = Xp[ok], Xn[ok]
for j in range(N_PERM):
    yp = yo.copy()
    for gi in grp:
        yp[gi] = rng.permutation(yo[gi])
    nul[j] = spearmanr(xo, yp)[0]
crows.append(dict(control="C1_interval_clustering", site="pooled", seg="all", n=int(ok.sum()),
                  rho_Tprev_Tnext=float(rho), p=perm_p(rho, nul)))
CT = pd.DataFrame(crows)
print("\n--- C1 interval clustering: Spearman rho(T_prev, T_next) (pooled on within-site ranks) ---")
for r in CT.itertuples():
    print(f"{r.site:>8} {r.seg:<8} n={r.n:4d}  rho={r.rho_Tprev_Tnext:+.3f}  p={r.p:.4f}")
print(
    "  A strongly POSITIVE rho here would mean short intervals cluster, which can manufacture an\n"
    "  apparent size~T_next from a size~T_prev relation (or vice versa) -> read the partial columns."
)

# --- C4: is the size~T_next relation made by the GREEDY DE-CLUSTERING? -----------------
# S5 selects peaks greedily highest-first with a >= 90 d exclusion, so a BIG peak wins its slot and
# suppresses any neighbour within 90 d -> the interval AFTER a big peak is mechanically stretched.
# The strength of that mechanism scales with the exclusion radius.  Rebuild the episode catalogue
# from the tremor counts at PK_SEP = 45 / 90 / 150 d (90 must reproduce S5 exactly) and see whether
# rho(size, T_next) scales with PK_SEP.  If it is flat, the de-clustering is not the driver.
print("\n--- C4 de-clustering sensitivity: rebuild the episode catalogue at PK_SEP = 45/90/150 d ---")
from scipy.signal import find_peaks  # noqa: E402

S5_T0, S5_T1 = pd.Timestamp("2010-01-06"), pd.Timestamp("2026-05-27")
TRAIL, PK_PCTL, ONSET_BACK, ONSET_FRAC, KD_MAX_KM = 30, 90, 45, 0.25, 40.0
BANDS_B = [round(44.0 + 0.5 * k, 1) for k in range(12)]

trb = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
trb["time"] = pd.to_datetime(trb["time"], errors="coerce")
trb = trb.dropna(subset=["time", "lat", "lon"])
trb["day"] = trb["time"].dt.floor("D")
trb = trb[(trb["day"] >= S5_T0) & (trb["day"] <= S5_T1)]
bx, by = to_km(trb["lon"].to_numpy(), trb["lat"].to_numpy())
dist, _ = cKDTree(np.column_stack([PX, PY])).query(np.column_stack([bx, by]))
trb = trb[dist <= KD_MAX_KM].copy()
trb["band"] = np.round(np.floor(trb["lat"].to_numpy() * 2.0) / 2.0, 1)
trb = trb[trb["band"].isin(BANDS_B)]

grid = pd.date_range(S5_T0, S5_T1, freq="D")
NG = len(grid)
gidx = pd.Series(np.arange(NG), index=grid)
cnts = np.zeros((NG, len(BANDS_B)))
bcol = {b: k for k, b in enumerate(BANDS_B)}
for (dd, bb), nn in trb.groupby(["day", "band"]).size().items():
    cnts[gidx[dd], bcol[round(float(bb), 1)]] = nn
csum = np.cumsum(cnts, axis=0)
trail = csum - np.vstack([np.zeros((TRAIL, cnts.shape[1])), csum[:-TRAIL]])


def detect(c, sep):
    thr = np.percentile(c, PK_PCTL)
    pk, _ = find_peaks(c)
    pk = pk[(pk >= TRAIL) & (c[pk] >= thr) & (c[pk] > 0)]
    keep = []
    for i in pk[np.argsort(-c[pk])]:
        if all(abs(i - j) >= sep for j in keep):
            keep.append(int(i))
    return sorted(keep)


c4rows = []
for sep in [45, 90, 150]:
    recs = []
    for k, b in enumerate(BANDS_B):
        c = trail[:, k]
        idx = detect(c, sep)
        pk = np.array([grid[i] for i in idx], dtype="datetime64[D]").astype(int)
        cc = np.array([c[i] for i in idx], float)
        tp = np.full(len(idx), np.nan)
        tn = np.full(len(idx), np.nan)
        if len(idx) > 1:
            tp[1:] = (pk[1:] - pk[:-1]).astype(float)
            tn[:-1] = (pk[1:] - pk[:-1]).astype(float)
        bandmed = float(np.median(cc)) if cc.size else np.nan
        det = np.full(len(idx), np.nan)
        for a, t in enumerate(pk):
            w = np.abs(pk - t) <= DETREND_HALFWIN_D
            m = float(np.median(cc[w])) if w.sum() >= 3 else bandmed
            det[a] = cc[a] / m if m > 0 else np.nan
        for a in range(len(idx)):
            recs.append(dict(band=b, seg="north" if b >= 46 else "central", peak=grid[idx[a]],
                             count30_raw=cc[a], count30_detr=det[a], T_prev=tp[a], T_next=tn[a]))
    D4 = pd.DataFrame(recs)
    if sep == 90:  # must reproduce S5
        got = set(f"T{r.band:.1f}_{pd.Timestamp(r.peak).date()}" for r in D4.itertuples())
        want = set(E["key"])
        print(f"    SEP=90 reproduces S5: {len(got & want)}/{len(want)} episode keys match "
              f"({len(got)} detected)")
        assert got == want, "SEP=90 rebuild does not reproduce s5_tremor_episodes.csv"
    for pcol in ["count30_raw", "count30_detr"]:
        res = pooled_test(BANDS_B, pcol, None, df=D4)
        c4rows.append(dict(pk_sep_d=sep, n_episodes=len(D4), proxy=pcol,
                           n_prev=res["prev"]["n"], rho_prev=res["prev"]["rho"], p_prev=res["prev"]["p"],
                           n_next=res["next"]["n"], rho_next=res["next"]["rho"], p_next=res["next"]["p"],
                           rho_next_partial=res["next"]["rho_par"], p_next_partial=res["next"]["p_par"]))
C4 = pd.DataFrame(c4rows)
CT = pd.concat(
    [
        CT,
        C4.assign(control="C4_declustering", site=lambda d: "det_sep" + d.pk_sep_d.astype(str),
                  seg="pooled"),
    ],
    ignore_index=True,
)
CT.to_csv(OUT_RES / "s8_predictability_controls.csv", index=False)
print(f"{'PK_SEP':>7} {'n_ep':>5} {'proxy':<14} {'rho_prev':>9} {'p':>7} {'rho_next':>9} {'p':>7} "
      f"{'rho_nx|pr':>9} {'p':>7}")
for r in C4.itertuples():
    print(f"{r.pk_sep_d:7d} {r.n_episodes:5d} {r.proxy:<14} {r.rho_prev:9.3f} {r.p_prev:7.4f} "
          f"{r.rho_next:9.3f} {r.p_next:7.4f} {r.rho_next_partial:9.3f} {r.p_next_partial:7.4f}")
print(
    "  If rho_next were made by the greedy exclusion it would GROW strongly with PK_SEP."
)

# --- C5: is rho_next actually LARGER than rho_prev?  bootstrap over sites --------------
print("\n--- C5 slip- vs time-predictable: bootstrap over the 12 bands, rho_next - rho_prev ---")
N_BOOT = 2000
bands_all = sorted(E.band.unique())
c5rows = []
for pname, pcol, _ in PROXIES:
    obs = PR[(PR.site == "pooled") & (PR.proxy == pname)].iloc[0]
    d_obs = obs.rho_next - obs.rho_prev
    boot = np.empty(N_BOOT)
    boot[:] = np.nan
    for j in range(N_BOOT):
        pick = rng.choice(bands_all, size=len(bands_all), replace=True)
        Xp, Xn, Y = [], [], []
        for b in pick:
            g = E[E.band == b]
            Xp.append(frank(g["T_prev"].to_numpy(float)))
            Xn.append(frank(g["T_next"].to_numpy(float)))
            Y.append(frank(g[pcol].to_numpy(float)))
        Xp, Xn, Y = np.concatenate(Xp), np.concatenate(Xn), np.concatenate(Y)
        op = np.isfinite(Xp) & np.isfinite(Y)
        on = np.isfinite(Xn) & np.isfinite(Y)
        if op.sum() < 5 or on.sum() < 5:
            continue
        boot[j] = spearmanr(Xn[on], Y[on])[0] - spearmanr(Xp[op], Y[op])[0]
    bb = boot[np.isfinite(boot)]
    lo, hi = np.percentile(bb, [2.5, 97.5])
    frac0 = float((bb <= 0).mean())
    c5rows.append(dict(control="C5_next_minus_prev", site="pooled", seg="all", proxy=pname,
                       diff_rho=float(d_obs), ci_lo=float(lo), ci_hi=float(hi),
                       frac_boot_le_0=frac0, n_boot=int(bb.size)))
    print(f"  {pname:<18} rho_next - rho_prev = {d_obs:+.3f}  boot 95 % CI [{lo:+.3f}, {hi:+.3f}]  "
          f"frac(boot <= 0) = {frac0:.3f}")
CT = pd.concat([CT, pd.DataFrame(c5rows)], ignore_index=True)
CT.to_csv(OUT_RES / "s8_predictability_controls.csv", index=False)

print("\nPower: |rho| detectable at 80 % power, alpha 0.05 two-sided")
for n in [12, 15, 17, 20, 100, 150, 195]:
    print(f"   n = {n:4d}  ->  |rho| >= {rho80(n):.3f}")

# ======================================================================================
# FIGURES
# ======================================================================================
log("figures ...")
C_BEF, C_AFT = "#2c6fbb", "#d1495b"
SEGC = {"north": "#1f77b4", "central": "#2ca02c"}

# ---- figure A ---------------------------------------------------------------------
fig, axes = plt.subplots(2, 4, figsize=(17.5, 8.4))
sub = EVA[EVA.ov_thresh == 0.30]
for ax, (k, kl) in zip(axes.ravel(), QUANTS):
    for xi, (sd, col) in enumerate([("before", C_BEF), ("after", C_AFT)]):
        v = sub.loc[sub.side == sd, k].to_numpy(float)
        v = v[np.isfinite(v)]
        jit = RNG.normal(0, 0.055, v.size)
        ax.plot(np.full(v.size, xi) + jit, v, "o", color=col, ms=6, alpha=0.65, mew=0)
        if v.size:
            ax.hlines(np.median(v), xi - 0.24, xi + 0.24, color=col, lw=2.6)
            ax.text(0.25 + 0.5 * xi, 0.985, f"n={v.size}", transform=ax.transAxes,
                    ha="center", va="top", fontsize=8, color=col)
    row = BA[(BA.ref == "pooled") & (BA.ov_thresh == 0.30) & (BA.quantity == k)].iloc[0]
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo - 0.22 * (hi - lo), hi)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["BEFORE", "AFTER"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_title(kl, fontsize=9.5, loc="left", fontweight="bold")
    ax.text(
        0.5,
        0.02,
        f"diff {row['diff']:+.3g}   MWU p={row['p_mwu']:.3f}\nnull sd {row['null_sd']:.3g}  z={row['z']:+.2f}  p_null={row['p_null']:.3f}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        bbox=dict(fc="white", ec="0.7", alpha=0.85, pad=2),
    )
    ax.grid(alpha=0.25, axis="y")
axes.ravel()[-1].axis("off")
axes.ravel()[-1].text(
    0.0,
    0.95,
    "S8 PART A — pooled (stacked) BEFORE vs AFTER\nthe three largest northern events\n"
    "(1317 2015-12, 1863 2017-02, 1605 2024-08)\n\n"
    "SITE = triangles of the reference footprint;\nan event joins when its footprint overlaps\n"
    ">= 30 % of the smaller footprint.\n"
    "Each event counts once per reference site.\n\n"
    "NULL = 200 random reference dates,\nuniform 2012-2024, same sites, same gap.\n"
    "z = observed median difference / null sd.\n\n"
    "CAVEAT: the 30 % rule makes the three sites\nnearly the whole northern margin, so the\n"
    "three 'sites' are largely the same events;\nref 1605 (2024) has only 2 AFTER events.",
    transform=axes.ravel()[-1].transAxes,
    va="top",
    ha="left",
    fontsize=9,
)
fig.suptitle(
    "S8 A — did the 2015-16 / 2017 / 2024 big northern ETS events change what followed at their site?",
    fontweight="bold",
)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_FIG / "s8_bigevents.png", dpi=145)
plt.close(fig)
log("  wrote figures/s8_bigevents.png")

# ---- figure B ---------------------------------------------------------------------
fig, axes = plt.subplots(3, 4, figsize=(18.5, 12.8))
for c, (pname, pcol, plab) in enumerate(PROXIES):
    for r, tag in enumerate(["prev", "next"]):
        ax = axes[r, c]
        res = POOL[pname][tag]
        if res[3] is None:
            ax.axis("off")
            continue
        xo, yo, so = res[3]
        for seg in ["north", "central"]:
            bands = sorted(E.loc[E.seg == seg, "band"].unique())
            m = np.isin(so, bands)
            ax.plot(xo[m], yo[m], "o", color=SEGC[seg], ms=5, alpha=0.6, mew=0, label=seg)
        # pooled trend line on ranks
        if np.isfinite(xo).sum() > 3:
            p1 = np.polyfit(xo, yo, 1)
            xs = np.linspace(0, 1, 10)
            ax.plot(xs, np.polyval(p1, xs), "-", color="0.25", lw=1.6)
        rho, p, n = res[0], res[1], res[2]
        model = "TIME-predictable" if tag == "prev" else "SLIP-predictable"
        ax.set_title(
            f"{model}: {plab}\nvs {'T_prev' if tag=='prev' else 'T_next'}  "
            f"rho={rho:+.3f} p={p:.3f} n={n} (80%:{rho80(n):.2f})",
            fontsize=8.5,
            loc="left",
        )
        ax.set_xlabel("within-site fractional rank of " + ("T_prev" if tag == "prev" else "T_next"))
        ax.set_ylabel("within-site fractional rank of size")
        ax.grid(alpha=0.25)
        if r == 0 and c == 0:
            ax.legend(fontsize=8, loc="upper left")

# row 3: growth check, per-site rho forest, power curve, text
ax = axes[2, 0]
for b, g in E.groupby("band"):
    seg = g["seg"].iloc[0]
    ax.plot(g["peak"], g["count30"], "-o", ms=3, lw=0.9, color=SEGC[seg], alpha=0.55)
gp = GR[(GR.band == "pooled") & (GR.proxy == "count30_raw")].iloc[0]
ax.set_title(
    f"network-growth check: count30 vs date\npooled (per-band z) rho={gp.rho:+.3f} p={gp.p:.2g}",
    fontsize=9,
    loc="left",
)
ax.set_ylabel("count30 (detections in 30 d at peak)")
ax.grid(alpha=0.25)

ax = axes[2, 1]
BANDSTR = {f"{b:.1f}" for b in E.band.unique()}
d = PR[(PR.proxy == "count30_detr") & (PR.site.isin(BANDSTR))].copy()
d = d.sort_values("site")
yy = np.arange(len(d))
ax.plot(d.rho_prev, yy - 0.15, "o", color="#7b3294", label="rho(size, T_prev)")
ax.plot(d.rho_next, yy + 0.15, "s", color="#e08214", label="rho(size, T_next)")
ax.axvline(0, color="0.4", lw=1)
for n_ in [d.n_prev.median()]:
    ax.axvspan(-rho80(int(n_)), rho80(int(n_)), color="0.85", alpha=0.6, zorder=0)
ax.set_yticks(yy)
ax.set_yticklabels(d.site)
ax.set_xlabel("Spearman rho")
ax.set_title(
    "per-site rho (detrended count30)\ngrey band = below 80 %-power detectability",
    fontsize=9,
    loc="left",
)
ax.legend(fontsize=8)
ax.grid(alpha=0.25, axis="x")

ax = axes[2, 2]
ns = np.arange(6, 220)
ax.plot(ns, [rho80(int(n)) for n in ns], "-", color="0.2")
for n_, lab in [(int(PR.loc[PR.site == "pooled", "n_prev"].max()), "pooled"), (16, "typical site")]:
    ax.plot([n_], [rho80(n_)], "o", color="#d1495b")
    ax.annotate(f"{lab}: n={n_}, |rho|>={rho80(n_):.2f}", (n_, rho80(n_)),
                textcoords="offset points", xytext=(8, 8), fontsize=8)
ax.set_xscale("log")
ax.set_xlabel("n pairs")
ax.set_ylabel("|rho| detectable at 80 % power")
ax.set_title("power", fontsize=9, loc="left")
ax.grid(alpha=0.25, which="both")

ax = axes[2, 3]
ax.axis("off")
SHORT = {"count30_raw": "cnt raw", "count30_detr": "cnt detr",
         "patch_Mw": "patch Mw", "patch_potency_km3": "patch pot"}
best = PR[PR.site == "pooled"]
txt = ["POOLED (within-site fractional ranks)", ""]
for r in best.itertuples():
    txt.append(
        f"{SHORT[r.proxy]:<9} T_prev {r.rho_prev:+.3f} p={r.p_prev:.3f} n={r.n_prev}\n"
        f"{'':<9} T_next {r.rho_next:+.3f} p={r.p_next:.3f} n={r.n_next}"
    )
txt += ["", "CONTROLS on the T_next result"]
for r in best.itertuples():
    if r.proxy.startswith("count30"):
        txt.append(f"C2 {SHORT[r.proxy]:<9} T_next|T_prev {r.rho_next_partial:+.3f} "
                   f"p={r.p_next_partial:.3f}")
for r in PR[PR.site == "pooled_sep150"].itertuples():
    if r.proxy.startswith("count30"):
        txt.append(f"C3 {SHORT[r.proxy]:<9} int>=150d    {r.rho_next:+.3f} p={r.p_next:.3f}")
for r in C4[C4.proxy == "count30_detr"].itertuples():
    txt.append(f"C4 re-detect PK_SEP={r.pk_sep_d:<3d}   {r.rho_next:+.3f} p={r.p_next:.4f}")
ctp = CT[(CT.control == "C1_interval_clustering") & (CT.site == "pooled")].iloc[0]
txt.append(f"C1 rho(T_prev,T_next)     {ctp.rho_Tprev_Tnext:+.3f} p={ctp.p:.3f}")
c5 = [r for r in c5rows if r["proxy"] == "count30_detr"][0]
txt.append(f"C5 rho_next-rho_prev {c5['diff_rho']:+.3f} "
           f"CI[{c5['ci_lo']:+.2f},{c5['ci_hi']:+.2f}]")
txt += [
    "",
    "TIME-predictable <=> size grows with T_prev",
    "SLIP-predictable <=> size sets T_next",
    "",
    "12 bands 44-50N, 13-20 episodes each;",
    "~14-month recurrence, ~16 yr record.",
    "Permutation within site, 2000 draws.",
    "",
    "READ: the T_next relation is significant in the",
    "TREMOR-COUNT proxy and survives C1-C4, but the",
    "MODEL-SLIP size proxy shows nothing at comparable",
    "power and C5 cannot separate rho_next from rho_prev.",
]
ax.text(0.0, 1.0, "\n".join(txt), va="top", ha="left", fontsize=7.1, family="monospace",
        transform=ax.transAxes)

fig.suptitle(
    "S8 B — is Cascadia ETS recurrence time-predictable (size ~ preceding interval) "
    "or slip-predictable (size ~ following interval)?",
    fontweight="bold",
)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig(OUT_FIG / "s8_predictability.png", dpi=140)
plt.close(fig)
log("  wrote figures/s8_predictability.png")

log("done.")
