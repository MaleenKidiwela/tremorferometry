#!/usr/bin/env python
"""S1b — why are the three northern 2022 events Q1/Q3 outliers?

S1 found that the northern confirmed patches starting 2022-03-31 (id 2015), 2022-05-26 (id 1583) and
2022-12-02 (id 1905) have Q1 = +21.1 / +19.6 / +28.9 km and Q3 = +3.5 / +4.6 / +5.4 km, while every
other northern year sits at Q1 ~ 0-10 km and Q3 ~ 0-2 km.  This script asks whether that is a SLIP-side
effect (the model put the slip somewhere unusual), a TREMOR-side effect (the tremor really moved
down-dip in 2022), a FOOTPRINT / measurement-aperture effect, or a combination.

--------------------------------------------------------------------------------------
INPUTS  (read-only; identical to s1_coregistration.py)
--------------------------------------------------------------------------------------
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz        label, days, patch_lat/lon/area
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm_mirror.npz mirror labels
  fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv  event sets (same filters as S1)
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
  catalogs/pnsn_tremor_cascadia_full.csv                                 PNSN tremor (time, lat, lon)
  /home/jovyan/cascadia-sse-strain                                       mesh, Ld (742 x 13), Slab2 depth f(lon,lat)
  tremor_slip_geometry/results/s1_events.csv, s1_controls.csv            the S1 tables (Q1..Q5 per event)

  Thrust cumulative slip per triangle, metres: raw_cum = -(Ld @ A_raw.T); mirror: +(Ld @ A_raw.T).
  Event window [ts, te] from the catalog dates on `days`; acc(mm) = 1e3*(raw_cum[:,te]-raw_cum[:,ts-30]), clip >= 0.
  Footprint F = triangles labelled with the patch id on any day in [ts, te].
  Local km frame and the down-dip / along-strike unit vectors are defined exactly as in S1.

  GROUPS (northern confirmed Mw>=6 patches, n=21):  "2022" = {2015, 1583, 1905};
  "2020-2021", "2023-2025" are the comparison groups named by the coordinator; "2010-2019" is the rest.

--------------------------------------------------------------------------------------
QUANTITIES ADDED HERE (beyond S1)
--------------------------------------------------------------------------------------
 Q1  slip-side vs tremor-side decomposition, per event:
   slip_depth_km              area*acc weighted mean Slab2 depth over F                    (S1 definition)
   fp_depth_km                AREA-weighted (slip-independent) mean Slab2 depth over F     -- the footprint itself
   slip_minus_fp_downdip_km   (slip centroid - area-weighted footprint centroid) . d       (+ = slip sits down-dip
                              inside its own footprint; - = slip concentrated up-dip)
   tremor_minus_fp_downdip_km (tremor centroid - footprint centroid) . d
   tremor_depth_med/p10/p90   Slab2 depth of the DUR tremor
   Q3_fp = tremor_depth_med - fp_depth_km    depth offset measured against the FOOTPRINT instead of the slip
   Q1_r20 / Q1_r40 / Q1_r60, Q3_r20/40/60    the same Q1/Q3 with the tremor aperture at 20 / 40 / 60 km
                              (aperture sensitivity: a small footprint + a 40 km disc reaches far down-dip)
 Q2  vbICA component mix, per event:
   per component k: pot_k = sum_F area_km2 * acc_k * 1e-6 (km^3), with acc_k = -1e3*Ld[:,k]*(A_raw[te,k]-A_raw[ts-30,k]);
   frac_k = pot_k / sum_k pot_k; the components with |frac_k| >= 0.10 are listed in `top_components`.
   comp_depth_k = area*|Ld[:,k]| weighted mean Slab2 depth of pattern k over ALL 742 triangles (a fixed property).
   mix_depth_km = sum_k (|frac_k| / sum|frac|) * comp_depth_k   -- the depth the component MIX alone implies.
   frac_k1k2 = frac_1 + frac_2  -- the two components that dominate almost every event; a LOW value flags an
   event the model built out of an unusual pattern mix.
 Q3  footprint up-dip test: fp_depth_km vs tremor_depth_med per northern event (column `Q3_fp`).
 Q4  mirror control by year: median Q1/Q3 of the mirror rows (from s1_controls.csv) per calendar year.
 Tremor-side epoch test: ALL northern tremor (lat >= 46 N, PLAN segment definition), p10/p50/p90 of the Slab2
   depth, for 2019-2021 / 2022 / 2023-2025 and per calendar year -- computed BOTH over the whole northern
   catalog and restricted to detections within 40 km of the union of the northern event footprints.

--------------------------------------------------------------------------------------
OUTPUTS (inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------
  results/s1b_2022_outliers.csv    one row per northern confirmed event, all columns above plus the S1 Q1/Q3,
                                   with a `group` label; the epoch tremor-depth table, the per-year mirror table
                                   and the per-year northern tremor table are appended as extra blocks with a
                                   `block` column.
  figures/s1b_2022_outliers.png    top: 3 maps (one per 2022 event) of the footprint triangles coloured by
                                   accumulated slip, DUR tremor epicentres, slip and tremor centroids, Slab2
                                   depth contours; bottom: Q1 vs footprint size, northern tremor-depth
                                   percentiles by year, and the mirror Q1/Q3 by year.
Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s1b_2022_outliers.py
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
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
OUT_RES = ROOT / "tremor_slip_geometry/results"
OUT_FIG = ROOT / "tremor_slip_geometry/figures"
COS47 = np.cos(np.deg2rad(47.0))
PRE_DAYS = 30
T0 = time.time()
TARGETS = [2015, 1583, 1905]


def log(m):
    print(f"[{time.time() - T0:6.1f}s] {m}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def to_lonlat(x, y):
    return x / (111.0 * COS47) - 124.0, y / 111.0 + 47.0


# ---------------------------------------------------------------- model + data
log("loading slip model ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402
from sse_strain.megathrust import slab_depth_interpolator  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
_, Ld = sol.component_slip_patterns()
slabf = slab_depth_interpolator(fault)
VLON, VLAT = np.asarray(fault.lon), np.asarray(fault.lat)  # (742, 3) triangle vertices

A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"])
CUM = {"thrust": -(Ld @ A_raw.T), "mirror": +(Ld @ A_raw.T)}

NPZ = {k: np.load(SSE / f"slip_patches_5mm{'' if k == 'thrust' else '_mirror'}.npz", allow_pickle=True)
       for k in ("thrust", "mirror")}
LABEL = {k: v["label"] for k, v in NPZ.items()}
DAYS = np.asarray(NPZ["thrust"]["days"], dtype=str)
PLAT, PLON = NPZ["thrust"]["patch_lat"], NPZ["thrust"]["patch_lon"]
PAREA = NPZ["thrust"]["patch_area_km2"]
PDEPTH = slabf(PLON, PLAT)
PX, PY = to_km(PLON, PLAT)
DIDX = {d: i for i, d in enumerate(DAYS)}

# fixed per-component depth (area * |dip-slip pattern| weighted)
CW = PAREA[:, None] * np.abs(Ld)
COMP_DEPTH = (CW * PDEPTH[:, None]).sum(0) / CW.sum(0)
log(f"  component depths (km): {np.array2string(COMP_DEPTH, precision=1)}")

log("loading tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TX, TY = to_km(tr["lon"].values, tr["lat"].values)
TDEP = slabf(tr["lon"].values, tr["lat"].values)
TDAY = (tr["t"].values.astype("datetime64[D]") - np.datetime64("2010-01-01")).astype(int)
TYEAR = tr["t"].dt.year.values
TLAT = tr["lat"].values
TREE = cKDTree(np.column_stack([TX, TY]))

s1ev = pd.read_csv(OUT_RES / "s1_events.csv")
s1ct = pd.read_csv(OUT_RES / "s1_controls.csv")
tvs = pd.read_csv(SSE / "tremor_vs_slip_patches.csv")


def local_frame(lon, lat):
    h = 0.05
    dzdx = (slabf(lon + h, lat) - slabf(lon - h, lat)) / (2 * h * 111.0 * COS47)
    dzdy = (slabf(lon, lat + h) - slabf(lon, lat - h)) / (2 * h * 111.0)
    g = np.array([dzdx, dzdy], float)
    d = g / np.linalg.norm(g)
    s = np.array([-d[1], d[0]])
    return d, (-s if s[1] < 0 else s)


# ---------------------------------------------------------------- per-event analysis
def analyse(pid, start, end, catalog="thrust"):
    ts, te = DIDX[start], DIDX[end]
    F = np.any(LABEL[catalog][ts:te + 1] == pid, axis=0)
    dA = A_raw[te] - A_raw[max(ts - PRE_DAYS, 0)]                      # (13,)
    sgn = -1.0 if catalog == "thrust" else +1.0
    acc_k = sgn * 1e3 * Ld * dA[None, :]                               # (742, 13) mm per component
    acc = np.clip(acc_k.sum(1), 0.0, None)                             # mm, as in S1

    w = np.where(F, PAREA * acc, 0.0)
    wa = np.where(F, PAREA, 0.0)                                       # area-only (footprint) weights
    sw, swa = w.sum(), wa.sum()
    r = dict(id=int(pid), start=start, end=end, n_tri_F=int(F.sum()),
             potency_km3=float((PAREA[F] * acc[F] * 1e-6).sum()),
             fp_area_km2=float(PAREA[F].sum()))
    if sw <= 0:
        return r
    sx, sy = (w * PX).sum() / sw, (w * PY).sum() / sw
    fx, fy = (wa * PX).sum() / swa, (wa * PY).sum() / swa
    r["slip_depth_km"] = float((w * PDEPTH).sum() / sw)
    r["fp_depth_km"] = float((wa * PDEPTH).sum() / swa)
    r["fp_lat_span_deg"] = float(PLAT[F].max() - PLAT[F].min())
    r["fp_strike_extent_km"] = float(r["fp_lat_span_deg"] * 111.0)
    r["fp_depth_min_km"] = float(PDEPTH[F].min())
    r["fp_depth_max_km"] = float(PDEPTH[F].max())
    slon, slat = to_lonlat(sx, sy)
    r["slip_lon"], r["slip_lat"] = float(slon), float(slat)
    r["fp_lon"], r["fp_lat"] = [float(v) for v in to_lonlat(fx, fy)]
    d, s = local_frame(slon, slat)
    r["slip_minus_fp_downdip_km"] = float((sx - fx) * d[0] + (sy - fy) * d[1])
    r["slip_minus_fp_strike_km"] = float((sx - fx) * s[0] + (sy - fy) * s[1])

    d0, d1 = (np.datetime64(start) - np.datetime64("2010-01-01")).astype(int), \
             (np.datetime64(end) - np.datetime64("2010-01-01")).astype(int)
    for R in (20.0, 40.0, 60.0):
        idx = TREE.query_ball_point(np.column_stack([PX[F], PY[F]]), r=R)
        near = np.unique(np.concatenate([np.asarray(i, int) for i in idx])) if len(idx) else np.array([], int)
        sel = near[(TDAY[near] >= d0) & (TDAY[near] <= d1)] if near.size else np.array([], int)
        tag = f"r{int(R)}"
        r[f"N_{tag}"] = int(sel.size)
        if sel.size == 0:
            continue
        tx, ty = TX[sel].mean(), TY[sel].mean()
        dz = TDEP[sel][np.isfinite(TDEP[sel])]
        med = np.median(dz) if dz.size else np.nan
        r[f"Q1_{tag}"] = float((tx - sx) * d[0] + (ty - sy) * d[1])
        r[f"Q3_{tag}"] = float(med - r["slip_depth_km"])
        if R == 40.0:
            r["tremor_lon"], r["tremor_lat"] = [float(v) for v in to_lonlat(tx, ty)]
            r["tremor_depth_med"] = float(med)
            r["tremor_depth_p10"] = float(np.percentile(dz, 10)) if dz.size else np.nan
            r["tremor_depth_p90"] = float(np.percentile(dz, 90)) if dz.size else np.nan
            r["tremor_minus_fp_downdip_km"] = float((tx - fx) * d[0] + (ty - fy) * d[1])
            r["Q3_fp"] = float(med - r["fp_depth_km"])
            r["_sel40"] = sel
    # component mix
    pot_k = (PAREA[F, None] * acc_k[F] * 1e-6).sum(0)                  # (13,) signed km^3
    tot = pot_k.sum()
    frac = pot_k / tot if tot != 0 else np.full(13, np.nan)
    aw = np.abs(frac) / np.abs(frac).sum()
    r["mix_depth_km"] = float((aw * COMP_DEPTH).sum())
    r["frac_k1k2"] = float(frac[1] + frac[2])
    r["top_components"] = "|".join(f"k{k}:{frac[k]:+.2f}@{COMP_DEPTH[k]:.0f}km"
                                   for k in np.argsort(-np.abs(frac)) if abs(frac[k]) >= 0.10)
    r["n_comp_10pct"] = int((np.abs(frac) >= 0.10).sum())
    for k in range(13):
        r[f"frac_k{k}"] = float(frac[k])
    return r


north = s1ev[s1ev.seg == "north"].copy()
north["year"] = pd.to_datetime(north.start).dt.year


def group_of(y):
    if y == 2022:
        return "2022"
    if y in (2020, 2021):
        return "2020-2021"
    if y >= 2023:
        return "2023-2025"
    return "2010-2019"


north["group"] = north.year.map(group_of)
log(f"analysing {len(north)} northern confirmed events ...")
rows, sel40 = [], {}
for _, e in north.iterrows():
    r = analyse(int(e["id"]), e["start"], e["end"])
    sel40[int(e["id"])] = r.pop("_sel40", np.array([], int))
    d0 = (np.datetime64(e["start"]) - np.datetime64("2010-01-01")).astype(int)
    d1 = (np.datetime64(e["end"]) - np.datetime64("2010-01-01")).astype(int)
    mwin = (TLAT >= 46.0) & np.isfinite(TDEP) & (TDAY >= d0) & (TDAY <= d1)
    nw = float(np.median(TDEP[mwin])) if mwin.sum() >= 50 else np.nan
    r["N_northwide"] = int(mwin.sum())
    r["tremor_depth_med_northwide"] = nw
    r["delta_vs_northwide"] = r.get("tremor_depth_med", np.nan) - nw
    r["dur_d"] = int(d1 - d0 + 1)
    r["north_rate_per_day"] = float(mwin.sum()) / float(d1 - d0 + 1)
    r.update(group=e["group"], year=int(e["year"]), Mw=float(e["Mw"]),
             Q1_s1=float(e["Q1_downdip_km"]), Q3_s1=float(e["Q3_depth_km"]), N_s1=int(e["N"]))
    rows.append(r)
ev = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
assert np.allclose(ev.Q1_s1, ev.Q1_r40, atol=1e-6), "Q1 does not reproduce S1"
assert np.allclose(ev.Q3_s1, ev.Q3_r40, atol=1e-6), "Q3 does not reproduce S1"
log("  reproduces S1 Q1/Q3 exactly")

# ---------------------------------------------------------------- tremor-side epoch test
NORTH_LAT = 46.0
allF = np.zeros(742, bool)
for _, e in north.iterrows():
    ts, te = DIDX[e["start"]], DIDX[e["end"]]
    allF |= np.any(LABEL["thrust"][ts:te + 1] == int(e["id"]), axis=0)
idx = TREE.query_ball_point(np.column_stack([PX[allF], PY[allF]]), r=40.0)
in_band = np.zeros(len(tr), bool)
in_band[np.unique(np.concatenate([np.asarray(i, int) for i in idx]))] = True
is_north = (TLAT >= NORTH_LAT) & np.isfinite(TDEP)
BAND_COVER = float(in_band[is_north].mean())

trem_rows = []
for name, mask_extra in (("north_all", np.ones(len(tr), bool)), ("north_in_event_band", in_band)):
    base = is_north & mask_extra
    for lab, yrs in (("2019-2021", (2019, 2020, 2021)), ("2022", (2022,)), ("2023-2025", (2023, 2024, 2025))):
        m = base & np.isin(TYEAR, yrs)
        z = TDEP[m]
        trem_rows.append(dict(block="tremor_epoch", subset=name, epoch=lab, n=int(m.sum()),
                              p10=float(np.percentile(z, 10)), p50=float(np.median(z)),
                              p90=float(np.percentile(z, 90)), mean=float(z.mean())))
    for y in range(2010, 2027):
        m = base & (TYEAR == y)
        if m.sum() < 50:
            continue
        z = TDEP[m]
        trem_rows.append(dict(block="tremor_year", subset=name, epoch=str(y), n=int(m.sum()),
                              p10=float(np.percentile(z, 10)), p50=float(np.median(z)),
                              p90=float(np.percentile(z, 90)), mean=float(z.mean())))
trem = pd.DataFrame(trem_rows)

# Mann-Whitney: 2022 northern tremor depth vs the flanking epochs
mw_rows = []
for name, mask_extra in (("north_all", np.ones(len(tr), bool)), ("north_in_event_band", in_band)):
    base = is_north & mask_extra
    a = TDEP[base & (TYEAR == 2022)]
    for lab, yrs in (("2019-2021", (2019, 2020, 2021)), ("2023-2025", (2023, 2024, 2025))):
        b = TDEP[base & np.isin(TYEAR, yrs)]
        u = mannwhitneyu(a, b)
        mw_rows.append(dict(block="tremor_mw", subset=name, epoch=f"2022 vs {lab}",
                            n=len(a), n2=len(b), p50=float(np.median(a) - np.median(b)),
                            mean=float(a.mean() - b.mean()), p10=float(u.pvalue)))
trem = pd.concat([trem, pd.DataFrame(mw_rows)], ignore_index=True)

# ---------------------------------------------------------------- mirror control by year
mir = s1ct[(s1ct.seg == "north") & (s1ct.N > 0)].copy()
mir["year"] = pd.to_datetime(mir.start).dt.year
mir_year = (mir.groupby("year")
            .agg(n=("id", "size"), Q1_med=("Q1_downdip_km", "median"), Q3_med=("Q3_depth_km", "median"),
                 N_med=("N", "median"))
            .reset_index())
mir_year["block"] = "mirror_year"
mir_all = s1ct[s1ct.N > 0].copy()
mir_all["year"] = pd.to_datetime(mir_all.start).dt.year
mir_2022 = mir_all[mir_all.year == 2022]
mir_other = mir_all[mir_all.year != 2022]

# ---------------------------------------------------------------- write CSV
ev_out = ev.copy()
ev_out["block"] = "event"
cols_front = ["block", "group", "id", "start", "end", "year", "Mw", "n_tri_F", "fp_area_km2", "potency_km3",
              "N_s1", "Q1_s1", "Q3_s1", "fp_lat_span_deg", "fp_strike_extent_km",
              "slip_depth_km", "fp_depth_km", "fp_depth_min_km", "fp_depth_max_km",
              "tremor_depth_med", "tremor_depth_p10", "tremor_depth_p90", "Q3_fp",
              "slip_minus_fp_downdip_km", "slip_minus_fp_strike_km", "tremor_minus_fp_downdip_km",
              "mix_depth_km", "frac_k1k2", "n_comp_10pct", "top_components",
              "N_northwide", "tremor_depth_med_northwide", "delta_vs_northwide",
              "dur_d", "north_rate_per_day",
              "N_r20", "Q1_r20", "Q3_r20", "N_r40", "Q1_r40", "Q3_r40", "N_r60", "Q1_r60", "Q3_r60",
              "slip_lon", "slip_lat", "fp_lon", "fp_lat", "tremor_lon", "tremor_lat"]
cols_front += [f"frac_k{k}" for k in range(13)]
ev_out = ev_out[[c for c in cols_front if c in ev_out]]
out = pd.concat([ev_out, trem, mir_year], ignore_index=True)
out.to_csv(OUT_RES / "s1b_2022_outliers.csv", index=False)
log(f"wrote {OUT_RES / 's1b_2022_outliers.csv'}")

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.28, wspace=0.25)

glon = np.linspace(-126.5, -120.5, 220)
glat = np.linspace(45.5, 50.2, 220)
GL, GA = np.meshgrid(glon, glat)
GD = slabf(GL, GA)

for j, pid in enumerate(TARGETS):
    ax = fig.add_subplot(gs[0, j])
    e = ev[ev.id == pid].iloc[0]
    ts, te = DIDX[e["start"]], DIDX[e["end"]]
    F = np.any(LABEL["thrust"][ts:te + 1] == pid, axis=0)
    dA = A_raw[te] - A_raw[max(ts - PRE_DAYS, 0)]
    acc = np.clip((-1e3 * Ld * dA[None, :]).sum(1), 0.0, None)
    ax.contour(GL, GA, GD, levels=np.arange(20, 61, 5), colors="0.75", linewidths=0.6, zorder=0)
    cs = ax.contour(GL, GA, GD, levels=[30, 40, 50], colors="0.45", linewidths=1.0, zorder=0)
    ax.clabel(cs, fmt="%.0f km", fontsize=6)
    polys = [np.column_stack([VLON[i], VLAT[i]]) for i in np.where(F)[0]]
    pc = PolyCollection(polys, array=acc[F], cmap="YlOrRd", edgecolor="0.5", linewidth=0.3, zorder=2, alpha=0.9)
    ax.add_collection(pc)
    plt.colorbar(pc, ax=ax, shrink=0.7, label="accumulated slip (mm)")
    s = sel40[pid]
    ax.plot(tr["lon"].values[s], tr["lat"].values[s], ".", ms=1.6, color="#1f5fbf", alpha=0.35, zorder=3,
            label=f"tremor DUR (n={len(s)})")
    ax.plot(e.slip_lon, e.slip_lat, "*", ms=20, mfc="#d62728", mec="k", zorder=5, label="slip centroid")
    ax.plot(e.tremor_lon, e.tremor_lat, "P", ms=13, mfc="#00b8d4", mec="k", zorder=5, label="tremor centroid")
    ax.plot(e.fp_lon, e.fp_lat, "s", ms=9, mfc="none", mec="k", mew=1.5, zorder=5, label="footprint centroid")
    ax.set_xlim(-126.3, -120.8)
    ax.set_ylim(45.8, 50.1)
    ax.set_aspect(1 / COS47)
    ax.set_title(f"id {pid}  {e['start']} .. {e['end']}\n"
                 f"Q1={e.Q1_s1:+.1f} km  Q3={e.Q3_s1:+.1f} km  |  slip {e.slip_depth_km:.1f} km, "
                 f"fp {e.fp_depth_km:.1f} km, tremor {e.tremor_depth_med:.1f} km", fontsize=8.5)
    ax.legend(fontsize=6.5, loc="lower left")
    ax.grid(alpha=0.2)
    ax.set_xlabel("lon", fontsize=8)
    if j == 0:
        ax.set_ylabel("lat", fontsize=8)

GC = {"2022": "#d62728", "2020-2021": "#1f5fbf", "2023-2025": "#1f9e4a", "2010-2019": "#888888"}

ax = fig.add_subplot(gs[1, 0])
for g, sub in ev.groupby("group"):
    ax.plot(sub.n_tri_F, sub.Q1_s1, "o", color=GC[g], ms=8, label=g, lw=0)
for _, r in ev.iterrows():
    ax.annotate(str(int(r.id)), (r.n_tri_F, r.Q1_s1), fontsize=5.5, xytext=(3, 3), textcoords="offset points")
rho, pv = spearmanr(ev.n_tri_F, ev.Q1_s1)
ax.axhline(0, color="0.7", lw=0.8)
ax.set_xlabel("footprint size (triangles)")
ax.set_ylabel("Q1 down-dip offset (km)")
ax.set_title(f"Q1 vs footprint size, north  (Spearman rho={rho:+.2f}, p={pv:.3f})", fontsize=9)
ax.legend(fontsize=7)
ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 1])
for subset, ls in (("north_all", "-"), ("north_in_event_band", "--")):
    t = trem[(trem.block == "tremor_year") & (trem.subset == subset)].copy()
    t["y"] = t.epoch.astype(int)
    ax.plot(t.y, t.p50, ls, marker="o", ms=4, color="#1f5fbf" if subset == "north_all" else "#00838f",
            label=f"{subset} p50")
    ax.fill_between(t.y, t.p10, t.p90, alpha=0.12,
                    color="#1f5fbf" if subset == "north_all" else "#00838f")
ax.axvspan(2021.5, 2022.5, color="#d62728", alpha=0.12, zorder=0)
ax.invert_yaxis()
ax.set_xlabel("year")
ax.set_ylabel("Slab2 depth of tremor (km)")
ax.set_title("northern (lat>=46N) tremor depth p10/p50/p90 by year\n(red band = 2022)", fontsize=9)
ax.legend(fontsize=7)
ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 2])
my = mir_year
ax.plot(my.year, my.Q1_med, "o-", color="#666666", label="mirror north Q1 median")
ax2 = ax.twinx()
ax2.plot(my.year, my.Q3_med, "s--", color="#b06000", label="mirror north Q3 median")
for a, lab in ((ax, "Q1 (km)"), (ax2, "Q3 (km)")):
    a.set_ylabel(lab)
ax.axvspan(2021.5, 2022.5, color="#d62728", alpha=0.12, zorder=0)
ax.axhline(0, color="0.8", lw=0.8)
n22 = len(mir_2022)
ax.set_title(f"MIRROR control by year (north, N>0)\n2022 rows n={n22}: "
             f"Q1 med {np.nan if n22 == 0 else np.median(mir_2022.Q1_downdip_km):.1f}, "
             f"other years med {np.median(mir_other.Q1_downdip_km):.1f} km", fontsize=9)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="best")
ax.grid(alpha=0.25)

fig.suptitle("S1b — the three northern 2022 Q1/Q3 outliers (ids 2015, 1583, 1905)", fontsize=12)
fig.savefig(OUT_FIG / "s1b_2022_outliers.png", dpi=135, bbox_inches="tight")
plt.close(fig)
log(f"wrote {OUT_FIG / 's1b_2022_outliers.png'}")

# ---------------------------------------------------------------- print
pd.set_option("display.width", 250)
print("\n" + "=" * 118)
print("S1b — 2022 NORTHERN OUTLIERS")
print("=" * 118)

print("\n(1) SLIP SIDE vs TREMOR SIDE — per northern event")
show = ["group", "id", "start", "n_tri_F", "N_s1", "Q1_s1", "Q3_s1", "slip_depth_km", "fp_depth_km",
        "tremor_depth_med", "Q3_fp", "slip_minus_fp_downdip_km", "tremor_minus_fp_downdip_km", "mix_depth_km"]
print(ev[show].to_string(index=False, float_format=lambda v: f"{v:8.2f}"))

print("\n(1b) group medians")
gm = ev.groupby("group").agg(
    n=("id", "size"), Q1=("Q1_s1", "median"), Q3=("Q3_s1", "median"), Q3_fp=("Q3_fp", "median"),
    slip_depth=("slip_depth_km", "median"), fp_depth=("fp_depth_km", "median"),
    tremor_depth=("tremor_depth_med", "median"), slip_m_fp=("slip_minus_fp_downdip_km", "median"),
    trem_m_fp=("tremor_minus_fp_downdip_km", "median"), mix_depth=("mix_depth_km", "median"),
    n_tri=("n_tri_F", "median"), N=("N_s1", "median")).reindex(["2010-2019", "2020-2021", "2022", "2023-2025"])
print(gm.to_string(float_format=lambda v: f"{v:8.2f}"))

print("\n(1c) tremor aperture sensitivity (Q1 / Q3 at 20, 40, 60 km)")
ap = ev[["group", "id", "N_r20", "Q1_r20", "Q3_r20", "N_r40", "Q1_r40", "Q3_r40", "N_r60", "Q1_r60", "Q3_r60"]]
print(ap.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))
print("\n  group medians of Q1 by aperture:")
print(ev.groupby("group")[["Q1_r20", "Q1_r40", "Q1_r60", "Q3_r20", "Q3_r40", "Q3_r60"]].median()
      .reindex(["2010-2019", "2020-2021", "2022", "2023-2025"]).to_string(float_format=lambda v: f"{v:8.2f}"))

print("\n(2) vbICA COMPONENT MIX  (components with |fraction of accumulated potency| >= 10%)")
print(f"    fixed component depths (km): " +
      "  ".join(f"k{k}={COMP_DEPTH[k]:.1f}" for k in range(13)))
for _, r in ev.iterrows():
    star = " <<<" if r.id in TARGETS else ""
    print(f"  {r.group:10s} {int(r.id):5d} {r.start}  mix_depth={r.mix_depth_km:5.1f} km  "
          f"n>=10%={int(r.n_comp_10pct)}  {r.top_components}{star}")

print("\n(2b) is the OUTLIER a component-mix effect?  frac_k1k2 = fraction of the accumulated potency carried "
      "by k1+k2\n     (the pair that dominates almost every event); Q1 against it:")
r12, p12 = spearmanr(ev.frac_k1k2, ev.Q1_s1)
r12b, p12b = spearmanr(ev.frac_k1k2, ev.Q3_s1)
print(f"     Spearman frac_k1k2 vs Q1 rho={r12:+.2f} p={p12:.4f};  vs Q3 rho={r12b:+.2f} p={p12b:.4f}")
lo = ev[ev.frac_k1k2 < 0.5]
hi = ev[ev.frac_k1k2 >= 0.5]
print(f"     unusual mix (k1+k2 < 0.5): n={len(lo)} ids {sorted(lo.id.astype(int))} -> Q1 med "
      f"{lo.Q1_s1.median():.1f} km, Q3 med {lo.Q3_s1.median():.1f} km")
print(f"     usual   mix (k1+k2 >=0.5): n={len(hi)} -> Q1 med {hi.Q1_s1.median():.1f} km, "
      f"Q3 med {hi.Q3_s1.median():.1f} km")

print("\n(3) FOOTPRINT up-dip test — footprint area-weighted depth vs tremor median depth")
print(ev.groupby("group").agg(fp_depth=("fp_depth_km", "median"), fp_shallow=("fp_depth_min_km", "median"),
                              fp_deep=("fp_depth_max_km", "median"), tremor=("tremor_depth_med", "median"),
                              Q3_fp=("Q3_fp", "median"))
      .reindex(["2010-2019", "2020-2021", "2022", "2023-2025"]).to_string(float_format=lambda v: f"{v:8.2f}"))

print("\n    per-event margin control — was the tremor deep only NEAR the patch, or everywhere in the north?")
print(ev[["group", "id", "start", "tremor_depth_med", "tremor_depth_med_northwide", "delta_vs_northwide",
          "N_s1", "N_northwide"]].to_string(index=False, float_format=lambda v: f"{v:9.2f}"))
print(f"\n    group medians of delta_vs_northwide (km): "
      + ", ".join(f"{g} {v:+.2f}" for g, v in
                  ev.groupby("group").delta_vs_northwide.median()
                  .reindex(["2010-2019", "2020-2021", "2022", "2023-2025"]).items()))

print(f"\n    NOTE: the 40-km band around the union of the 21 northern footprints contains "
      f"{BAND_COVER * 100:.1f}% of all lat>=46N detections, so the 'north_in_event_band' subset below is "
      f"{'identical to' if BAND_COVER > 0.999 else 'nearly identical to'} 'north_all'.")
print("\n    tremor-side epoch test — northern (lat>=46N) tremor Slab2 depth")
print(trem[trem.block == "tremor_epoch"][["subset", "epoch", "n", "p10", "p50", "p90", "mean"]]
      .to_string(index=False, float_format=lambda v: f"{v:9.2f}"))
print("\n    Mann-Whitney (p in the p10 column, medians/means are DIFFERENCES 2022 minus comparison):")
print(trem[trem.block == "tremor_mw"][["subset", "epoch", "n", "n2", "p50", "mean", "p10"]]
      .to_string(index=False, float_format=lambda v: f"{v:11.4g}"))

print("\n(4) MIRROR control by year (north, N>0)")
print(mir_year[["year", "n", "Q1_med", "Q3_med", "N_med"]].to_string(index=False,
                                                                     float_format=lambda v: f"{v:8.2f}"))
print(f"\n    mirror (both segments, N>0): 2022 n={len(mir_2022)} Q1 med "
      f"{np.median(mir_2022.Q1_downdip_km) if len(mir_2022) else float('nan'):.2f} "
      f"Q3 med {np.median(mir_2022.Q3_depth_km) if len(mir_2022) else float('nan'):.2f}  |  "
      f"other years n={len(mir_other)} Q1 med {np.median(mir_other.Q1_downdip_km):.2f} "
      f"Q3 med {np.median(mir_other.Q3_depth_km):.2f}")

print("\n(4b) WINDOW-PHASE control — northern tremor RATE during each event window "
      "(ETS peak = hundreds/day, inter-ETS = a few/day)")
w = ev[["group", "id", "start", "dur_d", "north_rate_per_day", "tremor_depth_med_northwide", "Q1_s1", "Q3_s1"]]
print(w.sort_values("north_rate_per_day").to_string(index=False, float_format=lambda v: f"{v:9.2f}"))
for c, lab in (("north_rate_per_day", "rate"), ("tremor_depth_med_northwide", "north-wide depth")):
    a1, b1 = spearmanr(ev[c], ev.Q1_s1)
    a3, b3 = spearmanr(ev[c], ev.Q3_s1)
    print(f"     Spearman {lab:17s} vs Q1 rho={a1:+.2f} p={b1:.4f};  vs Q3 rho={a3:+.2f} p={b3:.4f}")

print("\n(4c) CHAINED-FOOTPRINT check (PLAN 6) — footprint along-strike extent; a long footprint makes the "
      "single\n     down-dip projection at the slip centroid mix in an along-strike separation")
cf = ev[["group", "id", "fp_lat_span_deg", "fp_strike_extent_km", "n_tri_F", "Q1_s1", "Q3_s1"]].copy()
print(cf.sort_values("fp_strike_extent_km", ascending=False).head(8)
      .to_string(index=False, float_format=lambda v: f"{v:9.2f}"))
rc, pc = spearmanr(ev.fp_strike_extent_km, ev.Q1_s1)
print(f"     Spearman strike extent vs Q1 rho={rc:+.2f} p={pc:.4f}  "
      f"(2022 extents: " + ", ".join(f"{int(i)}={v:.0f} km" for i, v in
                                     zip(ev[ev.group == '2022'].id, ev[ev.group == '2022'].fp_strike_extent_km)) + ")")

print("\n(5) SIZE CONTROL — Spearman over the 21 northern events")
for c in ["n_tri_F", "N_s1", "potency_km3", "Mw"]:
    r1, p1 = spearmanr(ev[c], ev.Q1_s1)
    r3, p3 = spearmanr(ev[c], ev.Q3_s1)
    print(f"    {c:12s} vs Q1 rho={r1:+.2f} p={p1:.4f}   vs Q3 rho={r3:+.2f} p={p3:.4f}")
small = ev[ev.n_tri_F <= 90]
print(f"    among the {len(small)} SMALL northern events (<=90 triangles): "
      f"2022 Q1 med {small[small.group == '2022'].Q1_s1.median():.1f} km vs "
      f"non-2022 {small[small.group != '2022'].Q1_s1.median():.1f} km; "
      f"Q3 {small[small.group == '2022'].Q3_s1.median():.1f} vs "
      f"{small[small.group != '2022'].Q3_s1.median():.1f} km")

print("\nfiles:")
for p in [OUT_RES / "s1b_2022_outliers.csv", OUT_FIG / "s1b_2022_outliers.png"]:
    print(f"  {p} ({p.stat().st_size / 1024:.1f} kB)")
log("done")
