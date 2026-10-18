#!/usr/bin/env python
"""
S3b — detection-ROBUST tremor efficiency: slip per TREMOR-DAY instead of slip per detection.

Follow-up to S3 (scripts/s3_efficiency_maps.py).  S3 found that the margin-wide efficiency change
flips sign with the network normalisation (raw x0.85, normalised x1.17, both p < 1e-3), and that
the only shift-invariant result was a north-vs-central contrast (central x0.64 vs north x0.96,
Mann-Whitney p 1.4e-5).  That contrast is exactly what a spatially varying growth in detection
SENSITIVITY would fake: a triangle that gained a nearby station reports more detections for the
same fault behaviour, so eff = S/N falls with no change in the fault.  This script replaces the
detection COUNT by a count of tremor-ACTIVE DAYS, which saturates: a day on which a patch is
tremoring is one day whether the network reports 3 detections or 30.

INPUTS (read-only, outside this folder)
---------------------------------------
catalogs/pnsn_tremor_cascadia_full.csv                       time, lat, lon (depth is a placeholder)
fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
                                                             patch_lat/lon/depth_km/area_km2, days
fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz   A_raw (T x 13)
/home/jovyan/cascadia-sse-strain                             mesh + sse_strain (Ld, load_solution)
INPUT from this folder (written by S3, for the completeness diagnostic only):
tremor_slip_geometry/results/s3_efficiency_triangles.csv     N_E1..E3, log10_eff_ratio_E3_E1

EXACT DEFINITIONS
-----------------
Epochs           E1 2010-01-01..2014-12-31, E2 2015-01-01..2019-12-31, E3 2020-01-01..2026-05-27
                 (identical to S3; both the tremor and the slip sums are clipped to these bounds).
Tremor->triangle nearest triangle centroid in a local km frame
                 x = (lon + 124) * 111 * cos(47 deg), y = (lat - 47) * 111, dropped beyond 40 km.
c[day, tri]      number of detections assigned to that triangle on that calendar day (UTC).
Neighbourhood    nb(tri) = {tri} union {triangles sharing at least one MESH VERTEX with tri}
                 (vertices matched on the 0.25-deg grid, rounded to 1e-6 deg).  Median 12
                 triangles per neighbourhood -> a ~0.5-0.75 deg patch.
                 c_nb[day, tri] = sum of c[day, j] over j in nb(tri).
Tremor-days      D_own [tri, e]  = # days in epoch e with c[day, tri]  >= 1
                 D3_own[tri, e]  = # days in epoch e with c[day, tri]  >= 3
                 D_nb  [tri, e]  = # days in epoch e with c_nb[day, tri] >= 1
                 D3_nb [tri, e]  = # days in epoch e with c_nb[day, tri] >= 3
Daily slip rate  raw_cum_mm = -1e3 * (Ld @ A_raw.T); rate = np.diff(raw_cum_mm, axis=1) on days[1:];
                 7-day CENTRED moving average; pos_rate = np.clip(rate, 0, None).
S_mm[tri, e]     sum of pos_rate over the days of epoch e (mm), thrust-positive, as in S3.
eff_days         S_mm / D, mm of thrust slip per tremor-active day.  DEFINED only where
                 D >= 20 AND S_mm >= 5 mm; NaN otherwise.  Computed for each of the four D
                 variants above.
NO NORMALISATION is applied here, deliberately.  The S3 factor f[e] = tot[e]/mean(tot) is built
                 from DETECTION COUNTS, and detection counts are the quantity a tremor-day is
                 designed to be insensitive to.  Multiplying eff_days by f[e] would re-inject the
                 network-growth signal that the tremor-day definition removes, and (as S3 showed)
                 the constant log-offset alone decides the sign of the margin-wide median.  The
                 residual network effect on tremor-days is a threshold effect (a day with weak
                 tremor becomes detectable), which a scalar cannot correct; D3 (>=3 detections in
                 a day) is the sensitivity check for it, and the N/D diagnostic below measures it.
Log-ratios       log10(eff_days_E3 / eff_days_E1) and log10(eff_days_E2 / eff_days_E1), for
                 triangles defined in both epochs.
Tests            (a) Mann-Whitney U (two-sided) on log10(eff_days_E3/eff_days_E1),
                     north (centroid lat >= 46) vs central (44 <= lat < 46);
                 (b) Wilcoxon signed-rank (two-sided) of the margin-wide (lat >= 44) log-ratio
                     against 0.  Raw only -- see above.
Completeness     NperD[tri, e] = N[tri, e] / D_own[tri, e]  = mean detections per tremor-active
diagnostic       day, i.e. a direct per-triangle measure of detection sensitivity.
                 dlog_NperD = log10(NperD_E3 / NperD_E1), reported as a median per segment, and
                 correlated (Spearman, across triangles) with the S3 detection-based
                 log10(eff_E3/eff_E1).  If the S3 result was sensitivity, NperD must have risen
                 much more in the centre than in the north AND dlog_NperD must be strongly
                 ANTI-correlated with the S3 log-ratio (eff = S/N, so eff_ratio ~ -dlog N).
Decomposition    log10(eff_days_E3/eff_days_E1) = log10(S_E3/S_E1) - log10(D_E3/D_E1), so the
                 E3/E1 change is split into a SLIP term (physical: did the model put less slip
                 there?) and a TREMOR-DAY term (part physical, part detection threshold: did the
                 patch light up on more days?).  Reported as a per-segment median of each term
                 plus a north-vs-central Mann-Whitney on each, over the triangles defined in BOTH
                 epochs.  This is what decides whether the central drop is slip or detection.
Profiles         median eff_days (D_own) per 0.5-deg latitude bin per epoch, 90 % bootstrap over
                 the triangles in the bin (2000 resamples, seed 20260910).

OUTPUTS (inside tremor_slip_geometry/)
--------------------------------------
results/s3b_efficiency_triangles.csv  tri, lat, lon, depth_km, area_km2, seg, per epoch
                                      N, S_mm, D_own, D3_own, D_nb, D3_nb, eff_days_* for the four
                                      variants, and the E3/E1 + E2/E1 log-ratios of each variant
results/s3b_profiles.csv              lat_bin_centre, epoch, n_tri, median_eff_days, lo90, hi90
results/s3b_tests.csv                 variant, test, statistic, p, n, medians
results/s3b_completeness.csv          per-triangle N/D in E1 and E3, dlog_NperD, s3 log-ratio,
                                      and the dlog10 S / dlog10 D decomposition terms
results/s3b_decomposition.csv         per segment: median dlog10 S, dlog10 D, dlog10 eff_days
figures/s3b_efficiency_tremordays.png epoch maps of log10(eff_days), E3/E1 change map,
                                      along-strike profile, N/D diagnostic, and the
                                      dlog(N/D)-vs-S3-log-ratio scatter

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s3b_efficiency_tremordays.py
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

ROOT = Path("/home/jovyan/tremorferometry")
OUT = ROOT / "tremor_slip_geometry"
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
STRAIN = Path("/home/jovyan/cascadia-sse-strain")

MAX_ASSIGN_KM = 40.0
MIN_D = 20
MIN_S_MM = 5.0
NBOOT = 2000
SEED = 20260910
LAT_BIN = 0.5
EPOCHS = [("E1", "2010-01-01", "2014-12-31"),
          ("E2", "2015-01-01", "2019-12-31"),
          ("E3", "2020-01-01", "2026-05-27")]
VARIANTS = ["D_own", "D3_own", "D_nb", "D3_nb"]


def km_frame(lon, lat):
    x = (np.asarray(lon) + 124.0) * 111.0 * np.cos(np.deg2rad(47.0))
    y = (np.asarray(lat) - 47.0) * 111.0
    return np.column_stack([x, y])


def seg_of_lat(lat):
    lat = np.asarray(lat)
    out = np.full(lat.shape, "south", dtype=object)
    out[lat >= 46.0] = "north"
    out[(lat >= 44.0) & (lat < 46.0)] = "central"
    return out


def moving_average_centred(a, w=7, axis=-1):
    a = np.asarray(a, dtype=float)
    a = np.moveaxis(a, axis, -1)
    n = a.shape[-1]
    csum = np.concatenate([np.zeros(a.shape[:-1] + (1,)), np.cumsum(a, axis=-1)], axis=-1)
    half = w // 2
    idx = np.arange(n)
    lo = np.clip(idx - half, 0, n)
    hi = np.clip(idx + half + 1, 0, n)
    out = (csum[..., hi] - csum[..., lo]) / (hi - lo)
    return np.moveaxis(out, -1, axis)


def boot_median_ci(x, rng, nboot=NBOOT, lo=5, hi=95):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan, np.nan, np.nan, len(x)
    meds = np.median(rng.choice(x, size=(nboot, len(x)), replace=True), axis=1)
    return float(np.median(x)), float(np.percentile(meds, lo)), \
        float(np.percentile(meds, hi)), len(x)


def main():
    rng = np.random.default_rng(SEED)
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- mesh + adjacency
    sys.path.insert(0, str(STRAIN / "src"))
    from sse_strain import load_fault_ascii, load_solution

    fault = load_fault_ascii(str(STRAIN / "data/cascadia_mesh.txt"))
    vlon, vlat = np.asarray(fault.lon), np.asarray(fault.lat)

    pz = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
    plat, plon = pz["patch_lat"], pz["patch_lon"]
    pdepth, parea = pz["patch_depth_km"], pz["patch_area_km2"]
    days = pd.to_datetime(pd.Series(pz["days"]).astype(str)).values
    ntri = len(plat)
    tri_seg = seg_of_lat(plat)
    north_of_44 = plat >= 44.0

    vmap = defaultdict(list)
    for i in range(ntri):
        for k in range(3):
            vmap[(round(float(vlon[i, k]), 6), round(float(vlat[i, k]), 6))].append(i)
    ADJ = np.zeros((ntri, ntri), dtype=bool)
    np.fill_diagonal(ADJ, True)
    for tris in vmap.values():
        for a in tris:
            ADJ[a, tris] = True
    nb_sizes = ADJ.sum(axis=1)
    print(f"[mesh] {ntri} triangles; vertex-sharing neighbourhood size: "
          f"min {nb_sizes.min()}, median {int(np.median(nb_sizes))}, max {nb_sizes.max()} "
          f"(including self)")

    # ---------------------------------------------------------------- daily positive slip rate
    sol = load_solution(str(STRAIN / "data/2026-08-18"), fault)
    sol.invert(verbose=False)
    _, Ld = sol.component_slip_patterns()
    A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz",
                                  allow_pickle=True)["A_raw"])
    pos_rate = np.clip(moving_average_centred(np.diff(-1e3 * (Ld @ A_raw.T), axis=1), 7, axis=1),
                       0.0, None)
    rate_days = pd.DatetimeIndex(days[1:])

    # ---------------------------------------------------------------- tremor -> day x triangle
    tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                     usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr["time"])
    tree = cKDTree(km_frame(plon, plat))
    dist, idx = tree.query(km_frame(tr["lon"].values, tr["lat"].values), k=1)
    keep = dist <= MAX_ASSIGN_KM
    tr = tr.loc[keep].copy()
    tr["tri"] = idx[keep]
    tr["day"] = tr["time"].dt.floor("D")
    print(f"[tremor] {int(keep.sum()):,} detections assigned within {MAX_ASSIGN_KM:.0f} km")

    uday = np.array(sorted(tr["day"].unique()))
    day_pos = pd.Series(np.arange(len(uday)), index=pd.DatetimeIndex(uday))
    C = np.zeros((len(uday), ntri), dtype=np.int32)
    g = tr.groupby([tr["day"], tr["tri"]]).size()
    di = day_pos.loc[pd.DatetimeIndex(g.index.get_level_values(0))].values
    ti = g.index.get_level_values(1).values
    C[di, ti] = g.values
    uday_idx = pd.DatetimeIndex(uday)
    print(f"[tremor] {len(uday):,} distinct tremor days on the mesh, "
          f"{uday_idx[0].date()}..{uday_idx[-1].date()}; "
          f"total detections in C = {int(C.sum()):,}")

    C_nb = C.astype(np.int32) @ ADJ.astype(np.int32)   # (nday, ntri) neighbourhood day counts

    # ---------------------------------------------------------------- per-epoch tallies
    N = np.zeros((ntri, 3), dtype=int)
    S = np.zeros((ntri, 3), dtype=float)
    D = {v: np.zeros((ntri, 3), dtype=int) for v in VARIANTS}
    for k, (name, a, b) in enumerate(EPOCHS):
        dsel = (uday_idx >= pd.Timestamp(a)) & (uday_idx <= pd.Timestamp(b))
        Ce, Cnbe = C[dsel], C_nb[dsel]
        N[:, k] = Ce.sum(axis=0)
        D["D_own"][:, k] = (Ce >= 1).sum(axis=0)
        D["D3_own"][:, k] = (Ce >= 3).sum(axis=0)
        D["D_nb"][:, k] = (Cnbe >= 1).sum(axis=0)
        D["D3_nb"][:, k] = (Cnbe >= 3).sum(axis=0)
        smask = (rate_days >= pd.Timestamp(a)) & (rate_days <= pd.Timestamp(b))
        S[:, k] = pos_rate[:, np.asarray(smask)].sum(axis=1)
        nday_ep = int(np.sum(dsel))
        print(f"[epoch {name}] {a}..{b}: calendar days "
              f"{(pd.Timestamp(b) - pd.Timestamp(a)).days + 1}, tremor days on the mesh "
              f"{nday_ep}, detections {int(N[:, k].sum()):,}, slip {S[:, k].sum():.0f} mm-tri; "
              f"median per-triangle D_own {np.median(D['D_own'][north_of_44, k]):.0f} (lat>=44N)")

    eff = {}
    okv = {}
    for v in VARIANTS:
        e = np.full((ntri, 3), np.nan)
        o = (D[v] >= MIN_D) & (S >= MIN_S_MM)
        e[o] = S[o] / D[v][o]
        eff[v], okv[v] = e, o
        print(f"[defined] {v}: " + ", ".join(
            f"{EPOCHS[k][0]}={int(o[:, k].sum())}" for k in range(3)) +
            f" of {ntri} triangles (D>={MIN_D} & S>={MIN_S_MM} mm)")

    lr = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for v in VARIANTS:
            lr[(v, "E3/E1")] = np.log10(eff[v][:, 2] / eff[v][:, 0])
            lr[(v, "E2/E1")] = np.log10(eff[v][:, 1] / eff[v][:, 0])

    cols = dict(tri=np.arange(ntri), lat=plat, lon=plon, depth_km=pdepth, area_km2=parea,
                seg=tri_seg, nb_size=nb_sizes)
    for k, (name, _, _) in enumerate(EPOCHS):
        cols[f"N_{name}"] = N[:, k]
        cols[f"S_mm_{name}"] = S[:, k]
        for v in VARIANTS:
            cols[f"{v}_{name}"] = D[v][:, k]
            cols[f"eff_days_{v}_{name}"] = eff[v][:, k]
    for v in VARIANTS:
        cols[f"log10_ratio_{v}_E3_E1"] = lr[(v, "E3/E1")]
        cols[f"log10_ratio_{v}_E2_E1"] = lr[(v, "E2/E1")]
    tri_df = pd.DataFrame(cols)
    tri_df.to_csv(OUT / "results/s3b_efficiency_triangles.csv", index=False)
    print(f"\n[write] results/s3b_efficiency_triangles.csv ({ntri} rows)")

    # ---------------------------------------------------------------- tests
    print("\n=== tests (tremor-days, RAW; no network normalisation is applied -- a tremor-day "
          "already saturates the detection count) ===")
    trows = []
    for v in VARIANTS:
        for lab in ("E3/E1", "E2/E1"):
            x = lr[(v, lab)]
            nm = np.isfinite(x) & (tri_seg == "north")
            cm = np.isfinite(x) & (tri_seg == "central")
            if lab == "E3/E1" and nm.sum() >= 3 and cm.sum() >= 3:
                u = mannwhitneyu(x[nm], x[cm], alternative="two-sided")
                print(f"{v:7s} {lab} north vs central (MWU): n {int(nm.sum())} "
                      f"(median {np.median(x[nm]):+.4f}, x{10 ** np.median(x[nm]):.3f}) vs "
                      f"{int(cm.sum())} (median {np.median(x[cm]):+.4f}, "
                      f"x{10 ** np.median(x[cm]):.3f})  U {u.statistic:.1f}  p {u.pvalue:.4g}")
                trows.append(dict(variant=v, test=f"MWU north-vs-central log10({lab})",
                                  statistic=float(u.statistic), p=float(u.pvalue),
                                  n=f"{int(nm.sum())} vs {int(cm.sum())}",
                                  medians=f"north {np.median(x[nm]):+.4f}, "
                                          f"central {np.median(x[cm]):+.4f}"))
            sel = np.isfinite(x) & north_of_44
            if sel.sum() >= 6:
                w = wilcoxon(x[sel], alternative="two-sided")
                print(f"{v:7s} {lab} margin-wide (lat>=44N) vs 0 (Wilcoxon): n {int(sel.sum())}"
                      f"  median {np.median(x[sel]):+.4f} (x{10 ** np.median(x[sel]):.3f})  "
                      f"W {w.statistic:.1f}  p {w.pvalue:.4g}")
                trows.append(dict(variant=v, test=f"Wilcoxon margin-wide log10({lab}) vs 0",
                                  statistic=float(w.statistic), p=float(w.pvalue),
                                  n=int(sel.sum()),
                                  medians=f"median {np.median(x[sel]):+.4f} "
                                          f"(x{10 ** np.median(x[sel]):.3f})"))
    pd.DataFrame(trows).to_csv(OUT / "results/s3b_tests.csv", index=False)
    print("[write] results/s3b_tests.csv")

    # ---------------------------------------------------------------- completeness diagnostic
    print("\n=== completeness diagnostic: detections per tremor-day (N / D_own) ===")
    Dn = D["D_own"]
    good = (Dn[:, 0] >= MIN_D) & (Dn[:, 2] >= MIN_D)
    NperD = np.full((ntri, 3), np.nan)
    for k in range(3):
        m = Dn[:, k] >= MIN_D
        NperD[m, k] = N[m, k] / Dn[m, k]
    with np.errstate(divide="ignore", invalid="ignore"):
        dlogNperD = np.log10(NperD[:, 2] / NperD[:, 0])

    s3 = pd.read_csv(OUT / "results/s3_efficiency_triangles.csv")
    s3lr = s3["log10_eff_ratio_E3_E1"].values

    comp = pd.DataFrame(dict(tri=np.arange(ntri), lat=plat, lon=plon, seg=tri_seg,
                             D_own_E1=Dn[:, 0], D_own_E3=Dn[:, 2], N_E1=N[:, 0], N_E3=N[:, 2],
                             NperD_E1=NperD[:, 0], NperD_E3=NperD[:, 2],
                             dlog10_NperD_E3_E1=dlogNperD,
                             s3_log10_eff_ratio_E3_E1=s3lr,
                             s3b_log10_ratio_D_own_E3_E1=lr[("D_own", "E3/E1")]))
    comp.to_csv(OUT / "results/s3b_completeness.csv", index=False)
    for seg in ("north", "central", "south"):
        m = good & (tri_seg == seg)
        if m.sum() >= 3:
            print(f"{seg:8s} n {int(m.sum()):3d}  N/D  E1 {np.median(NperD[m, 0]):5.2f} -> "
                  f"E3 {np.median(NperD[m, 2]):5.2f}  median dlog10 "
                  f"{np.median(dlogNperD[m]):+.4f} (x{10 ** np.median(dlogNperD[m]):.3f})")
    mnc = good & np.isin(tri_seg, ["north", "central"])
    nn, cc = good & (tri_seg == "north"), good & (tri_seg == "central")
    if nn.sum() >= 3 and cc.sum() >= 3:
        u = mannwhitneyu(dlogNperD[nn], dlogNperD[cc], alternative="two-sided")
        print(f"north vs central dlog10(N/D) (MWU): U {u.statistic:.1f}  p {u.pvalue:.4g}")
        trows.append(dict(variant="diagnostic", test="MWU north-vs-central dlog10(N/D)",
                          statistic=float(u.statistic), p=float(u.pvalue),
                          n=f"{int(nn.sum())} vs {int(cc.sum())}",
                          medians=f"north {np.median(dlogNperD[nn]):+.4f}, "
                                  f"central {np.median(dlogNperD[cc]):+.4f}"))
    m = good & np.isfinite(s3lr) & north_of_44
    rho, prho = spearmanr(dlogNperD[m], s3lr[m])
    print(f"Spearman dlog10(N/D) vs S3 log10(eff_E3/eff_E1), lat>=44N, n {int(m.sum())}: "
          f"rho {rho:+.3f}  p {prho:.3g}")
    trows.append(dict(variant="diagnostic",
                      test="Spearman dlog10(N/D) vs S3 log10(eff_E3/eff_E1) [lat>=44N]",
                      statistic=float(rho), p=float(prho), n=int(m.sum()), medians=""))
    m2 = good & np.isfinite(lr[("D_own", "E3/E1")]) & north_of_44
    rho2, prho2 = spearmanr(dlogNperD[m2], lr[("D_own", "E3/E1")][m2])
    print(f"Spearman dlog10(N/D) vs S3b log10(eff_days_E3/E1), lat>=44N, n {int(m2.sum())}: "
          f"rho {rho2:+.3f}  p {prho2:.3g}   <- should be much weaker if tremor-days are robust")
    trows.append(dict(variant="diagnostic",
                      test="Spearman dlog10(N/D) vs S3b log10(eff_days D_own E3/E1) [lat>=44N]",
                      statistic=float(rho2), p=float(prho2), n=int(m2.sum()), medians=""))
    pd.DataFrame(trows).to_csv(OUT / "results/s3b_tests.csv", index=False)
    print("[write] results/s3b_completeness.csv (+ diagnostics appended to s3b_tests.csv)")

    # ---------------------------------------------------------------- decomposition
    # log10(eff_days E3/E1) = log10(S E3/E1) - log10(D E3/E1)
    print("\n=== decomposition of log10(eff_days E3/E1) into slip and tremor-day terms "
          "(D_own; triangles defined in BOTH epochs) ===")
    both = okv["D_own"][:, 0] & okv["D_own"][:, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        dlogS = np.log10(S[:, 2] / S[:, 0])
        dlogD = np.log10(Dn[:, 2].astype(float) / Dn[:, 0])
        dlogE = dlogS - dlogD
    drows = []
    for seg in ("north", "central", "south"):
        m = both & (tri_seg == seg)
        if m.sum() < 3:
            continue
        print(f"{seg:8s} n {int(m.sum()):3d}   dlog10 S {np.median(dlogS[m]):+.4f} "
              f"(x{10 ** np.median(dlogS[m]):.3f})   dlog10 D {np.median(dlogD[m]):+.4f} "
              f"(x{10 ** np.median(dlogD[m]):.3f})   -> dlog10 eff_days "
              f"{np.median(dlogE[m]):+.4f} (x{10 ** np.median(dlogE[m]):.3f})")
        drows.append(dict(seg=seg, n_tri=int(m.sum()),
                          median_dlog10_S=float(np.median(dlogS[m])),
                          ratio_S=float(10 ** np.median(dlogS[m])),
                          median_dlog10_D=float(np.median(dlogD[m])),
                          ratio_D=float(10 ** np.median(dlogD[m])),
                          median_dlog10_eff_days=float(np.median(dlogE[m])),
                          ratio_eff_days=float(10 ** np.median(dlogE[m]))))
    pd.DataFrame(drows).to_csv(OUT / "results/s3b_decomposition.csv", index=False)
    nb_, cb_ = both & (tri_seg == "north"), both & (tri_seg == "central")
    for lab, arr in (("dlog10 S", dlogS), ("dlog10 D", dlogD), ("dlog10 eff_days", dlogE)):
        if nb_.sum() >= 3 and cb_.sum() >= 3:
            u = mannwhitneyu(arr[nb_], arr[cb_], alternative="two-sided")
            print(f"  north vs central {lab:16s}: U {u.statistic:7.1f}  p {u.pvalue:.4g}")
            trows.append(dict(variant="decomposition",
                              test=f"MWU north-vs-central {lab} (E3/E1, D_own both-defined)",
                              statistic=float(u.statistic), p=float(u.pvalue),
                              n=f"{int(nb_.sum())} vs {int(cb_.sum())}",
                              medians=f"north {np.median(arr[nb_]):+.4f}, "
                                      f"central {np.median(arr[cb_]):+.4f}"))
    comp["dlog10_S_E3_E1"] = dlogS
    comp["dlog10_D_own_E3_E1"] = dlogD
    comp.to_csv(OUT / "results/s3b_completeness.csv", index=False)
    pd.DataFrame(trows).to_csv(OUT / "results/s3b_tests.csv", index=False)
    print("[write] results/s3b_decomposition.csv (+ terms added to s3b_completeness.csv)")

    # ---------------------------------------------------------------- profiles
    edges = np.arange(np.floor(plat.min() / LAT_BIN) * LAT_BIN,
                      np.ceil(plat.max() / LAT_BIN) * LAT_BIN + LAT_BIN, LAT_BIN)
    prows = []
    for i in range(len(edges) - 1):
        c = 0.5 * (edges[i] + edges[i + 1])
        inb = (plat >= edges[i]) & (plat < edges[i + 1])
        for k, (name, _, _) in enumerate(EPOCHS):
            m_, lo_, hi_, n_ = boot_median_ci(eff["D_own"][inb, k], rng)
            prows.append(dict(lat_bin_centre=c, epoch=name, n_tri=n_, median_eff_days=m_,
                              lo90=lo_, hi90=hi_))
    prof = pd.DataFrame(prows)
    prof.to_csv(OUT / "results/s3b_profiles.csv", index=False)
    print(f"[write] results/s3b_profiles.csv ({len(prof)} rows)")
    print("\n=== along-strike profile, median eff_days (D_own; mm slip per tremor-day) ===")
    piv = prof.pivot(index="lat_bin_centre", columns="epoch", values="median_eff_days")
    ntr = prof.pivot(index="lat_bin_centre", columns="epoch", values="n_tri")
    print(pd.concat([ntr.add_prefix("nTri_"), piv.add_prefix("effd_")], axis=1)
          .round(3).to_string())

    print("\n=== margin-wide (lat>=44N) medians of eff_days (D_own) ===")
    for k, (name, _, _) in enumerate(EPOCHS):
        sel = north_of_44 & okv["D_own"][:, k]
        print(f"{name}: n {int(sel.sum()):3d}  median {np.median(eff['D_own'][sel, k]):.4f} "
              f"mm per tremor-day")

    # ---------------------------------------------------------------- figure
    verts = [np.column_stack([vlon[i], vlat[i]]) for i in range(ntri)]
    fig = plt.figure(figsize=(19, 11))
    gs = fig.add_gridspec(2, 12, height_ratios=[2.6, 1.25], hspace=0.26, wspace=1.4)
    ev = eff["D_own"]
    lv = np.log10(ev[np.isfinite(ev)])
    vmin, vmax = np.percentile(lv, 2), np.percentile(lv, 98)
    map_axes = []
    for k, (name, a, b) in enumerate(EPOCHS):
        ax = fig.add_subplot(gs[0, 3 * k:3 * k + 3])
        ax.add_collection(PolyCollection(verts, facecolors="none", edgecolors="0.82", lw=0.3))
        m = np.isfinite(ev[:, k])
        sc = ax.scatter(plon[m], plat[m], c=np.log10(ev[m, k]), s=26, cmap="viridis",
                        vmin=vmin, vmax=vmax, edgecolors="none")
        ax.set_title(f"{name}  {a[:7]}..{b[:7]}\nlog10 eff_days (D_own), {int(m.sum())} tri",
                     fontsize=9)
        map_axes.append(ax)
    cb = fig.colorbar(sc, ax=map_axes, fraction=0.03, pad=0.01)
    cb.set_label("log10 (mm slip per tremor-day)", fontsize=8)

    axc = fig.add_subplot(gs[0, 9:12])
    axc.add_collection(PolyCollection(verts, facecolors="none", edgecolors="0.82", lw=0.3))
    x = lr[("D_own", "E3/E1")]
    m = np.isfinite(x)
    lim = float(np.nanpercentile(np.abs(x[m]), 95)) if m.sum() else 1.0
    sc2 = axc.scatter(plon[m], plat[m], c=x[m], s=26, cmap="RdBu_r", vmin=-lim, vmax=lim,
                      edgecolors="none")
    m44 = m & north_of_44
    axc.set_title(f"change  log10(eff_days E3 / E1)\n{int(m.sum())} tri, median "
                  f"{np.nanmedian(x[m]):+.3f} (lat>=44N: {np.nanmedian(x[m44]):+.3f})",
                  fontsize=9)
    fig.colorbar(sc2, ax=axc, fraction=0.03, pad=0.01).set_label("log10 ratio", fontsize=8)
    for a_ in map_axes + [axc]:
        a_.set_xlim(plon.min() - 0.4, plon.max() + 0.4)
        a_.set_ylim(plat.min() - 0.4, plat.max() + 0.4)
        a_.set_aspect(1.0 / np.cos(np.deg2rad(47.0)))
        a_.axhline(46, color="k", lw=0.6, ls="--", alpha=0.6)
        a_.axhline(44, color="k", lw=0.6, ls="--", alpha=0.6)
        a_.set_xlabel("lon", fontsize=8)
        a_.tick_params(labelsize=7)
    map_axes[0].set_ylabel("lat", fontsize=8)

    colors = {"E1": "tab:blue", "E2": "tab:orange", "E3": "tab:red"}
    axp = fig.add_subplot(gs[1, 0:4])
    for name, _, _ in EPOCHS:
        s = prof[prof.epoch == name].sort_values("lat_bin_centre")
        axp.fill_between(s.lat_bin_centre, s.lo90, s.hi90, color=colors[name], alpha=0.18)
        axp.plot(s.lat_bin_centre, s.median_eff_days, "o-", color=colors[name], ms=3.5, lw=1.6,
                 label=name)
    axp.set_yscale("log")
    for xv in (44, 46):
        axp.axvline(xv, color="k", lw=0.6, ls="--", alpha=0.6)
    axp.set_xlabel("latitude of 0.5-deg bin centre")
    axp.set_ylabel("median mm slip / tremor-day")
    axp.set_title("along-strike profile, eff_days (90 % bootstrap over triangles)", fontsize=9)
    axp.grid(alpha=0.3)
    axp.legend(fontsize=8)

    axd = fig.add_subplot(gs[1, 4:8])
    segcol = {"north": "tab:blue", "central": "tab:green", "south": "0.6"}
    for seg in ("south", "central", "north"):
        m = good & (tri_seg == seg)
        axd.scatter(NperD[m, 0], NperD[m, 2], s=16, c=segcol[seg], label=seg, alpha=0.8,
                    edgecolors="none")
    hi = np.nanmax(NperD[good]) * 1.05
    axd.plot([0, hi], [0, hi], "k:", lw=1)
    axd.set_xlim(0, hi)
    axd.set_ylim(0, hi)
    axd.set_xlabel("detections per tremor-day, E1")
    axd.set_ylabel("detections per tremor-day, E3")
    axd.set_title("detection-sensitivity diagnostic N/D\n(above the 1:1 line = more sensitive "
                  "in E3)", fontsize=9)
    axd.grid(alpha=0.3)
    axd.legend(fontsize=8)

    axs = fig.add_subplot(gs[1, 8:12])
    for seg in ("central", "north"):
        m = good & (tri_seg == seg) & np.isfinite(s3lr)
        axs.scatter(dlogNperD[m], s3lr[m], s=16, c=segcol[seg], label=f"{seg} (S3)", alpha=0.85,
                    edgecolors="none")
    for seg in ("central", "north"):
        m = good & (tri_seg == seg) & np.isfinite(lr[("D_own", "E3/E1")])
        axs.scatter(dlogNperD[m], lr[("D_own", "E3/E1")][m], s=16, marker="^",
                    facecolors="none", edgecolors=segcol[seg], label=f"{seg} (S3b days)",
                    alpha=0.85)
    axs.axhline(0, color="k", lw=0.8)
    axs.axvline(0, color="k", lw=0.8)
    axs.set_xlabel("dlog10(N / D), E3 vs E1  (detection sensitivity growth)")
    axs.set_ylabel("log10(eff E3 / E1)")
    axs.set_title(f"S3 detection-based (dots) rho {rho:+.2f}  vs\n"
                  f"S3b tremor-day (triangles) rho {rho2:+.2f}", fontsize=9)
    axs.grid(alpha=0.3)
    axs.legend(fontsize=7)

    fig.suptitle("S3b — tremor efficiency measured in TREMOR-DAYS (detection-count robust)",
                 fontsize=13)
    fig.savefig(OUT / "figures/s3b_efficiency_tremordays.png", dpi=150, bbox_inches="tight")
    print("[write] figures/s3b_efficiency_tremordays.png")


if __name__ == "__main__":
    main()
