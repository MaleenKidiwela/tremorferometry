#!/usr/bin/env python
"""S15 TASK A — amplitude snapshot check: is the shipped `A_raw` really an EARLIER solution snapshot?

PLAN.md Round 5, S15(A).  Closes the publishability-memo TO DO (d): "re-derive A_raw from the
2026-08-18 solution (snapshot mismatch flagged by S12b) and confirm the 27-event set and exponents
are unchanged".

Background.  S12b (NOTES.md 2026-09-10) reported that the `A_raw` stored in
`fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz` "comes from an
earlier snapshot and correlates only 0.87-0.99" with `Solution.component_amplitudes()` of the
2026-08-18 solution, and therefore did NOT use it.  Every result of S1-S14 (events, footprints,
slip rates, the sub-linear exponent) is built on that stored `A_raw`, so the claim, if true, would
put the whole study on a stale model.

This script rebuilds the daily amplitudes from the 2026-08-18 solution EXACTLY as
`fault_tomography/inversion/sse_trailing_products.py` does (the script that wrote the npz):

    A = sol.component_amplitudes(); t = pd.to_datetime(np.asarray(sol.dates))
    Ad = pd.DataFrame(A, index=t).groupby(level=0).mean()          # collapse duplicate dates
    days = pd.date_range(Ad.index.min(), Ad.index.max(), freq="D") # daily grid
    Ad = Ad.reindex(days).interpolate(limit_area="inside")         # fill interior gaps
    A_new = Ad.values

and compares it with the stored `A_raw` column by column.  It then re-runs the two things that
depend on the amplitudes:

  (2) the slip-patch labelling (the S1 event definition).  `raw_cum = -(Ld @ A.T)`,
      `inc30 = raw_cum[t] - raw_cum[t-30]` in mm, threshold 5 mm/30 d, space-time connected
      components under same-day vertex adjacency + same-triangle consecutive days.  The component
      labeller is COPIED VERBATIM from `fault_tomography/inversion/sse_slip_patches.py`
      (GAP = 0, BRIDGE_KM = 0, i.e. the settings that produced `slip_patches_5mm.npz`); that
      script is neither modified nor re-run.  For each of the 27 confirmed events of
      `results/s1_events.csv` we ask: does a new component exist whose start AND end are within
      +-10 d of the stored ones and whose footprint overlaps the stored footprint by >= 70 %?

  (3) the S9/S11 count exponent b(C, R7).  Tremor measure, window, 60-km footprint neighbourhood,
      slip-rate construction and Theil-Sen fitting are copied from
      `scripts/s11_scaling_robustness.py`; the ONLY change is A_raw -> A_new.  Reported paired,
      per event, against the stored-amplitude exponents (recomputed here with the same code, and
      cross-checked against the stored `results/s11_exponents.csv`).

The S12b comparison is also reproduced deliberately (the "misaligned" rows of the output table):
`sol.component_amplitudes()` has ONE ROW PER SOLUTION EPOCH, and those epochs are neither unique
nor daily (7169 epochs, 5271 distinct dates, 1898 duplicated dates and 1898 missing days), while
`A_raw` has one row per DAY of the 7169-day grid.  The two arrays have the same LENGTH by
coincidence, so an index-by-index comparison silently compares different days.

======================================================================================
INPUTS (read-only)
======================================================================================
  /home/jovyan/cascadia-sse-strain                                  mesh + sse_strain loaders
  .../cascadia-sse-strain/data/2026-08-18                           the solution
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz    label (T x 742), days,
                                                                     patch_lat/lon/area_km2
  tremor_slip_geometry/results/s1_events.csv                         the 27 confirmed events
  tremor_slip_geometry/results/s11_exponents.csv                     stored b(C,R7) for the check
  catalogs/pnsn_tremor_cascadia_full.csv                             PNSN tremor

OUTPUTS (all inside tremor_slip_geometry/)
  results/s15_amplitude_check.csv   per component: r, RMS diff / column sd, max |diff| and its
                                    date, for the ALIGNED (daily-grid) and the MISALIGNED
                                    (index-by-index, S12b-style) comparison
  results/s15_events_check.csv      per confirmed event: the best-matching new component, its
                                    start/end offsets, footprint overlap and Jaccard, peak rate,
                                    and the PASS/FAIL verdict of the +-10 d / >= 70 % test
  results/s15_exponent_check.csv    per event: b(C, R7) with A_raw and with A_new (stored and new
                                    footprints), the paired difference, plus a summary row

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s15_amplitude_check.py
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
from scipy.ndimage import binary_dilation, uniform_filter1d
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import theilslopes, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_RES.mkdir(parents=True, exist_ok=True)
REPO = "/home/jovyan/cascadia-sse-strain"

RNG = np.random.default_rng(20260910)

THR_MM = 5.0             # slip-patch threshold, mm per 30 d  (sse_slip_patches.py)
WIN = 30                 # trailing window, days              (sse_slip_patches.py)
GAP = 0                  # no temporal bridging   -> slip_patches_5mm.npz
BRIDGE_KM = 0.0          # no spatial bridging    -> slip_patches_5mm.npz
DATE_TOL = 10            # +-10 d on start AND end   (task A2)
OVERLAP_MIN = 0.70       # >= 70 % footprint overlap (task A2)

PRE_A, POST_A = 10, 10   # exponent window [start-10, end+10]  (S9/S11)
NB_RADIUS_KM = 60.0      # neighbourhood / tremor radius       (S9/S11)
MIN_FIT_DAYS = 8
COS47 = np.cos(np.deg2rad(47.0))

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


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


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan, 0
    return (float(np.median(v)), float(np.percentile(v, 25)),
            float(np.percentile(v, 75)), int(v.size))


def boot_median_ci(v, nboot=2000, lo=5, hi=95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


print("=" * 108)
print("S15 TASK A — amplitude snapshot check (is the shipped A_raw the 2026-08-18 solution?)")
print("=" * 108)

# ==================================================================================
# 1. Rebuild the daily amplitudes exactly as sse_trailing_products.py does
# ==================================================================================
log("loading the 2026-08-18 solution ...")
sys.path.insert(0, f"{REPO}/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402

fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt")
sol = load_solution(f"{REPO}/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()

A_epoch = sol.component_amplitudes()                       # one row per SOLUTION EPOCH
t_epoch = pd.to_datetime(np.asarray(sol.dates))            # datetime64[D] -> ns index
n_ep = len(t_epoch)
n_uniq = int(pd.Index(t_epoch).nunique())

Ad = pd.DataFrame(A_epoch, index=t_epoch).groupby(level=0).mean()
days = pd.date_range(Ad.index.min(), Ad.index.max(), freq="D")
Ad = Ad.reindex(days).interpolate(limit_area="inside")
A_new = Ad.values.copy()                                   # (T, K) float64, daily
K = A_new.shape[1]

z = np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)
A_raw = z["A_raw"].astype(float)                           # stored float32 -> float64
zdays = pd.to_datetime(z["days"])

print(f"\nsolution epochs      : {n_ep} rows, {n_uniq} distinct dates "
      f"({n_ep - n_uniq} duplicated dates)")
print(f"daily grid           : {len(days)} days {days[0].date()}..{days[-1].date()} "
      f"({len(days) - n_uniq} days missing from the solution epochs)")
print(f"A_new {A_new.shape} float64   A_raw {A_raw.shape} {z['A_raw'].dtype} "
      f"| stored days identical: {bool(len(zdays) == len(days) and (zdays == days).all())}")

# per-component comparison, ALIGNED (both on the daily grid) and MISALIGNED (S12b style)
rows = []
for k in range(K):
    a, b = A_new[:, k], A_raw[:, k]
    m = np.isfinite(a) & np.isfinite(b)
    r = float(np.corrcoef(a[m], b[m])[0, 1])
    dif = a - b
    rms = float(np.sqrt(np.nanmean(dif[m] ** 2)))
    sd = float(np.std(b[m]))
    j = int(np.nanargmax(np.abs(np.where(m, dif, np.nan))))
    # the S12b comparison: raw epoch rows against daily rows, index by index
    am = A_epoch[:, k]
    mm = np.isfinite(am) & np.isfinite(b)
    r_mis = float(np.corrcoef(am[mm], b[mm])[0, 1])
    rms_mis = float(np.sqrt(np.mean((am[mm] - b[mm]) ** 2)))
    rows.append(dict(component=int(sol.component_index[k]) if hasattr(sol, "component_index") else k,
                     k=k, sd_A_new=float(np.std(a[m])), sd_A_raw=sd,
                     pearson_r_aligned=r, rms_diff=rms, rms_over_sd=rms / sd if sd else np.nan,
                     max_abs_diff=float(np.abs(dif[m]).max()),
                     max_abs_diff_over_sd=float(np.abs(dif[m]).max() / sd) if sd else np.nan,
                     date_of_max_diff=str(days[j].date()),
                     pearson_r_misaligned=r_mis, rms_diff_misaligned=rms_mis,
                     rms_over_sd_misaligned=rms_mis / sd if sd else np.nan))
amp = pd.DataFrame(rows)
amp.to_csv(OUT_RES / "s15_amplitude_check.csv", index=False)

print("\n(1) A_new (rebuilt from the 2026-08-18 solution) vs the stored A_raw")
print("     k  sd(A_new)  sd(A_raw)   r_aligned   rms/sd   max|d|/sd  date of max |d|   "
      "| r (S12b index-aligned)")
for _, rr in amp.iterrows():
    print(f"    {int(rr.k):2d}  {rr.sd_A_new:9.4f}  {rr.sd_A_raw:9.4f}  {rr.pearson_r_aligned:+10.7f} "
          f"{rr.rms_over_sd:8.2e} {rr.max_abs_diff_over_sd:11.2e}  {rr.date_of_max_diff}   "
          f"|  {rr.pearson_r_misaligned:+.4f}")
print(f"    ALL components: r_aligned min {amp.pearson_r_aligned.min():.9f}, "
      f"rms/sd max {amp.rms_over_sd.max():.3e}, max|diff| {amp.max_abs_diff.max():.3e} "
      f"(= float32 storage round-off; A_raw is stored as float32)")
print(f"    S12b-style index-by-index r: {amp.pearson_r_misaligned.min():.3f}"
      f"..{amp.pearson_r_misaligned.max():.3f}  <-- reproduces the reported 0.87-0.99")
same32 = float(np.mean(A_new.astype(np.float32) == A_raw.astype(np.float32)))
print(f"    fraction of the {A_new.size} entries identical after casting A_new to float32: {same32:.6f}")

# ==================================================================================
# 2. Slip-patch labelling with A_new  (labeller copied from sse_slip_patches.py)
# ==================================================================================
log("re-running the slip-patch labelling with A_new ...")
NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL_OLD = NPZ["label"]                       # (T, 742)
DAYS = np.asarray(NPZ["days"], dtype=str)
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
INC_OLD = NPZ["inc30_mm"].astype(float)        # (T, 742) mm, as shipped
PX, PY = to_km(PLON, PLAT)
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
NP = len(PLAT)
assert NT == A_new.shape[0] == A_raw.shape[0], "day-grid / amplitude length mismatch"

V = fault.vertices_m()
XY = np.column_stack([PX, PY])


def build_inc(A):
    """raw_cum = -(Ld @ A.T) (m, thrust-positive) -> inc30 (mm within (t-30, t])."""
    raw_cum = -(Ld @ np.nan_to_num(A).T)                       # (742, T) m
    inc = np.full_like(raw_cum, np.nan)
    inc[:, WIN:] = (raw_cum[:, WIN:] - raw_cum[:, :-WIN]) * 1e3
    return raw_cum, inc


# triangle adjacency: shared vertex (rounded to 1 m)  -- verbatim from sse_slip_patches.py
key = {}
adj = [set() for _ in range(NP)]
for p in range(NP):
    for v in np.round(V[p]).astype(np.int64):
        key.setdefault(tuple(v), []).append(p)
for ps in key.values():
    for a_ in ps:
        for b_ in ps:
            if a_ != b_:
                adj[a_].add(b_)
EA = np.array(sorted({(a_, b_) for a_ in range(NP) for b_ in adj[a_] if a_ < b_}))


def label_components(B, gap, bridge_km):
    """Verbatim from fault_tomography/inversion/sse_slip_patches.py."""
    D = B.copy()
    if bridge_km > 0:
        pr = cKDTree(XY).query_pairs(bridge_km / 2, output_type="ndarray")
        W = csr_matrix((np.ones(len(pr) * 2 + NP),
                        (np.r_[pr[:, 0], pr[:, 1], np.arange(NP)],
                         np.r_[pr[:, 1], pr[:, 0], np.arange(NP)])), shape=(NP, NP))
        D = (W @ D.astype(np.int32)) > 0
    if gap > 0:
        D = binary_dilation(D, structure=np.ones((1, 2 * int(np.ceil(gap / 2)) + 1), bool))
    nodes = np.flatnonzero(D.ravel())
    nid = -np.ones(NP * NT, np.int64)
    nid[nodes] = np.arange(len(nodes))
    pt, tt_ = np.where(D[:, :-1] & D[:, 1:])
    r1, c1 = nid[pt * NT + tt_], nid[pt * NT + tt_ + 1]
    pa, ta = np.where(D[EA[:, 0]] & D[EA[:, 1]])
    r2, c2 = nid[EA[pa, 0] * NT + ta], nid[EA[pa, 1] * NT + ta]
    Gr = coo_matrix((np.ones(len(r1) + len(r2)), (np.r_[r1, r2], np.r_[c1, c2])),
                    shape=(len(nodes), len(nodes)))
    _, labD = connected_components(Gr, directed=False)
    L = np.zeros((NP, NT), np.int32)
    pp, tt = np.where(B)
    L[pp, tt] = labD[nid[pp * NT + tt]] + 1
    u, inv_ = np.unique(L[L > 0], return_inverse=True)
    L[L > 0] = inv_ + 1
    return L, len(u)


RAW_NEW, INC_NEW = build_inc(A_new)
RAW_OLD_A, INC_OLD_A = build_inc(A_raw)

B_new = np.nan_to_num(INC_NEW) > THR_MM
B_oldA = np.nan_to_num(INC_OLD_A) > THR_MM
B_stored = np.nan_to_num(INC_OLD.T) > THR_MM        # from the shipped inc30_mm (float32)
print(f"\n(2) patch-days above {THR_MM:g} mm/30 d: A_new {B_new.sum()}, "
      f"A_raw (recomputed) {B_oldA.sum()}, shipped inc30_mm {B_stored.sum()}; "
      f"A_new vs A_raw disagree on {int((B_new ^ B_oldA).sum())} patch-days, "
      f"A_new vs shipped on {int((B_new ^ B_stored).sum())}")
print(f"    max |inc30(A_new) - inc30(A_raw)| = {np.nanmax(np.abs(INC_NEW - INC_OLD_A)):.3e} mm; "
      f"vs the shipped inc30_mm {np.nanmax(np.abs(INC_NEW - INC_OLD.T)):.3e} mm (float32 storage)")

L_new, ncomp_new = label_components(B_new, GAP, BRIDGE_KM)
print(f"    connected components with A_new: {ncomp_new}   "
      f"(shipped label array has {int(LABEL_OLD.max())})")

# component table for the new labelling (only what the event check needs)
new_comp = {}
pc_all, tc_all = np.where(L_new > 0)
lab_all = L_new[pc_all, tc_all]
order = np.argsort(lab_all, kind="stable")
pc_all, tc_all, lab_all = pc_all[order], tc_all[order], lab_all[order]
bounds = np.searchsorted(lab_all, np.arange(1, ncomp_new + 2))
for c in range(1, ncomp_new + 1):
    s0, s1 = bounds[c - 1], bounds[c]
    if s1 <= s0:
        continue
    pcs, tcs = pc_all[s0:s1], tc_all[s0:s1]
    fp = np.unique(pcs)
    ts_, te_ = int(tcs.min()), int(tcs.max())
    kmax = int(np.argmax(INC_NEW[pcs, tcs]))
    new_comp[c] = dict(ts=ts_, te=te_, fp=fp, npd=int(s1 - s0),
                       peak=float(INC_NEW[pcs[kmax], tcs[kmax]]),
                       area=float(PAREA[fp].sum()),
                       lat=float(np.average(PLAT[fp], weights=PAREA[fp])))

ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"checking the {len(ev_in)} confirmed events against the new labelling ...")

erows = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    F_old = np.flatnonzero(np.any(LABEL_OLD[ts:te + 1] == pid, axis=0))
    fold = set(F_old.tolist())
    best = None
    best_any = None
    for c, d in new_comp.items():
        inter = len(fold & set(d["fp"].tolist()))
        if inter == 0:
            continue
        ov = inter / max(len(fold), 1)
        jac = inter / len(fold | set(d["fp"].tolist()))
        cand = dict(new_id=c, start_off=d["ts"] - ts, end_off=d["te"] - te,
                    overlap=ov, jaccard=jac, n_tri_new=len(d["fp"]), n_tri_old=len(fold),
                    peak_mm_per30d=d["peak"], patch_days_new=d["npd"],
                    start_new=DAYS[d["ts"]], end_new=DAYS[d["te"]])
        if best_any is None or ov > best_any["overlap"]:
            best_any = cand
        if abs(cand["start_off"]) <= DATE_TOL and abs(cand["end_off"]) <= DATE_TOL:
            if best is None or ov > best["overlap"]:
                best = cand
    hit = best if best is not None else (best_any or {})
    passed = bool(best is not None and best["overlap"] >= OVERLAP_MIN)
    erows.append(dict(id=pid, start=r["start"], end=r["end"], seg=r["seg"], Mw=float(r["Mw"]),
                      PASS=passed,
                      matched_by_dates=bool(best is not None),
                      **{k: hit.get(k, np.nan) for k in
                         ("new_id", "start_new", "end_new", "start_off", "end_off",
                          "overlap", "jaccard", "n_tri_old", "n_tri_new",
                          "peak_mm_per30d", "patch_days_new")},
                      best_overlap_any_date=(best_any or {}).get("overlap", np.nan)))
ev_chk = pd.DataFrame(erows)
ev_chk.to_csv(OUT_RES / "s15_events_check.csv", index=False)

print(f"\n    PASS ({DATE_TOL:+d} d on start AND end, footprint overlap >= {OVERLAP_MIN:.0%}): "
      f"{int(ev_chk.PASS.sum())}/{len(ev_chk)}")
print("      id  start        end          seg      new_id  d_start  d_end  overlap  jaccard  "
      "n_tri old/new  peak mm/30d  PASS")
for _, rr in ev_chk.iterrows():
    print(f"    {int(rr.id):5d}  {rr.start}  {rr.end}  {rr.seg:7s} "
          f"{('%6d' % rr.new_id) if np.isfinite(rr.new_id) else '   n/a'}  "
          f"{rr.start_off:+7.0f} {rr.end_off:+6.0f}  {rr.overlap:7.4f}  {rr.jaccard:7.4f}  "
          f"{int(rr.n_tri_old):4d}/{int(rr.n_tri_new):<4d}      {rr.peak_mm_per30d:9.2f}   "
          f"{'PASS' if rr.PASS else 'FAIL'}")

# ==================================================================================
# 3. The count exponent b(C, R7) with A_new  (machinery copied from s11_scaling_robustness.py)
# ==================================================================================
log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DAY = tr["t"].values.astype("datetime64[ns]").astype("datetime64[D]")
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
log(f"  {len(tr)} detections {tr['t'].min().date()}..{tr['t'].max().date()}")

DAYS64 = np.asarray(DAYS, dtype="datetime64[D]")


def rate7(A):
    """R7 building blocks: CUM_MM (742,T) mm thrust-positive -> clipped 7-d centred mean rate."""
    cum = -1e3 * (Ld @ np.nan_to_num(A).T)
    diff = np.zeros_like(cum)
    diff[:, 1:] = np.diff(cum, axis=1)
    return np.clip(uniform_filter1d(diff, size=7, axis=1, mode="nearest"), 0.0, None)


R7_NEW = rate7(A_new)
R7_OLD = rate7(A_raw)

xrows = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    wd = DAYS64[w0:w1 + 1]

    # ---- footprints: the STORED one (as S9/S11 use) and the NEW one
    F_old = np.any(LABEL_OLD[ts:te + 1] == pid, axis=0)
    row = ev_chk[ev_chk.id == pid].iloc[0]
    F_new = F_old.copy()
    if np.isfinite(row.new_id):
        F_new = np.zeros(NP, bool)
        F_new[new_comp[int(row.new_id)]["fp"]] = True

    out = dict(id=pid, start=r["start"], seg=r["seg"], nW=nW)
    for tag, F in (("old", F_old), ("new", F_new)):
        fx, fy = PX[F], PY[F]
        hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
        nb = (np.unique(np.concatenate([np.asarray(h, int) for h in hits]))
              if len(hits) else np.array([], int))
        nb = np.union1d(nb, np.flatnonzero(F))
        thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
        near = (np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
                if len(thits) else np.array([], int))
        # daily tremor count in the window
        dday = (TR_DAY[near] - wd[0]).astype(int)
        inw = (dday >= 0) & (dday < nW)
        C = np.bincount(dday[inw], minlength=nW).astype(float)
        aw = PAREA[nb][:, None]
        for src, RAT in (("Araw", R7_OLD), ("Anew", R7_NEW)):
            R = (aw * RAT[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)
            ok = (C > 0) & (R > 0)
            b = (ts_slope_fast(np.log10(R[ok]), np.log10(C[ok] + 1.0))
                 if ok.sum() >= MIN_FIT_DAYS else np.nan)
            out[f"b_{src}_fp{tag}"] = b
            out[f"ndays_{src}_fp{tag}"] = int(ok.sum())
        out[f"n_near_fp{tag}"] = int(near.size)
        out[f"n_nb_fp{tag}"] = int(nb.size)
        out[f"C_sum_fp{tag}"] = float(C.sum())
    xrows.append(out)

xdf = pd.DataFrame(xrows)
xdf["diff_Anew_minus_Araw"] = xdf.b_Anew_fpold - xdf.b_Araw_fpold
xdf["diff_fullnew_minus_stored"] = xdf.b_Anew_fpnew - xdf.b_Araw_fpold

print("\n(3) count exponent b(C, R7), stored footprints; A_raw vs A_new")
print("       id  start        seg      n_days   b(A_raw)   b(A_new)     diff   | b(A_new, new fp)")
for _, rr in xdf.iterrows():
    print(f"    {int(rr.id):5d}  {rr.start}  {rr.seg:7s} {int(rr.ndays_Araw_fpold):5d}  "
          f"{rr.b_Araw_fpold:+9.6f}  {rr.b_Anew_fpold:+9.6f}  {rr.diff_Anew_minus_Araw:+9.2e}  |  "
          f"{rr.b_Anew_fpnew:+9.6f}")

sm = {}
for col in ("b_Araw_fpold", "b_Anew_fpold", "b_Anew_fpnew"):
    m, q1, q3, n = med_iqr(xdf[col].values)
    lo, hi = boot_median_ci(xdf[col].values)
    sm[col] = (m, q1, q3, n, lo, hi)
    print(f"    median {col:14s} {m:+.6f} [{q1:+.4f}, {q3:+.4f}] (n={n}), "
          f"90 % boot CI of the median [{lo:+.4f}, {hi:+.4f}], b<1 in "
          f"{int((xdf[col] < 1).sum())}/{int(np.isfinite(xdf[col]).sum())}")

d = xdf.diff_Anew_minus_Araw.values
d = d[np.isfinite(d)]
if np.allclose(d, 0):
    print(f"    paired difference b(A_new) - b(A_raw) = 0 for all {d.size} events "
          "(Wilcoxon undefined: all differences are exactly zero)")
    wp = np.nan
else:
    wp = float(wilcoxon(d).pvalue)
    print(f"    paired difference b(A_new) - b(A_raw): median {np.median(d):+.3e}, "
          f"max |d| {np.abs(d).max():.3e}, Wilcoxon p = {wp:.3g}")
d2 = xdf.diff_fullnew_minus_stored.values
d2 = d2[np.isfinite(d2)]
w2 = np.nan if np.allclose(d2, 0) else float(wilcoxon(d2).pvalue)
print(f"    with the NEW footprints as well: median difference {np.median(d2):+.3e}, "
      f"max |d| {np.abs(d2).max():.3e}" + ("" if np.isnan(w2) else f", Wilcoxon p = {w2:.3g}"))

# cross-check against the stored S11 table
try:
    s11 = pd.read_csv(OUT_RES / "s11_exponents.csv")
    hit = s11[(s11.measure == "C") & (s11.slip_var == "R7") & (s11.subset == "all")]
    if len(hit):
        print(f"    CHECK stored s11_exponents.csv b(C,R7,all) median = "
              f"{float(hit.iloc[0].b_median):+.6f}  vs recomputed here "
              f"{sm['b_Araw_fpold'][0]:+.6f}")
except Exception as exc:                                        # pragma: no cover
    print(f"    (s11 cross-check unavailable: {exc})")

x0 = xdf.iloc[0]
print(f"    CHECK fast Theil-Sen vs scipy on event {int(x0.id)}: see s11 (same helper, verified there)")

srow = dict(id=-1, start="SUMMARY", seg="all", nW=np.nan,
            b_Araw_fpold=sm["b_Araw_fpold"][0], b_Anew_fpold=sm["b_Anew_fpold"][0],
            b_Anew_fpnew=sm["b_Anew_fpnew"][0],
            diff_Anew_minus_Araw=float(np.median(d)) if d.size else np.nan,
            diff_fullnew_minus_stored=float(np.median(d2)) if d2.size else np.nan,
            wilcoxon_p_Anew_minus_Araw=wp, wilcoxon_p_fullnew=w2,
            boot_lo_Araw=sm["b_Araw_fpold"][4], boot_hi_Araw=sm["b_Araw_fpold"][5],
            boot_lo_Anew=sm["b_Anew_fpold"][4], boot_hi_Anew=sm["b_Anew_fpold"][5],
            n_events_pass_A2=int(ev_chk.PASS.sum()),
            max_abs_amp_diff=float(amp.max_abs_diff.max()),
            min_r_aligned=float(amp.pearson_r_aligned.min()),
            min_r_misaligned=float(amp.pearson_r_misaligned.min()),
            max_r_misaligned=float(amp.pearson_r_misaligned.max()))
pd.concat([xdf, pd.DataFrame([srow])], ignore_index=True).to_csv(
    OUT_RES / "s15_exponent_check.csv", index=False)

# cache A_new for TASK B (outside the project folder)
SCRATCH = Path("/tmp/claude-1000/-home-jovyan-tremorferometry/"
               "03499166-a677-42d8-9397-e2cf5e28621f/scratchpad")
try:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    np.savez(SCRATCH / "s15_Anew.npz", A_new=A_new, days=DAYS, R7=R7_NEW.astype(np.float32))
    log(f"cached A_new + R7 for task B in {SCRATCH}/s15_Anew.npz (outside the project)")
except Exception as exc:                                        # pragma: no cover
    log(f"(cache not written: {exc})")

print("\nwrote results/s15_amplitude_check.csv, results/s15_events_check.csv, "
      "results/s15_exponent_check.csv")
log("TASK A DONE")
