#!/usr/bin/env python
"""
S3 — tremor-efficiency maps per epoch: how much modelled thrust slip the fault accumulates
per PNSN tremor detection, per mesh triangle, in three epochs, and whether that changed.

PLAN.md step S3 (Q5 / H3).  No dv/v anywhere in this study.

INPUTS (all read-only, outside this folder)
-------------------------------------------
catalogs/pnsn_tremor_cascadia_full.csv
    PNSN tremor master catalogue: time, lat, lon (the `depth` column is a PLACEHOLDER, unused).
fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
    patch_lat, patch_lon, patch_depth_km, patch_area_km2 (742 mesh triangles); days (T strings).
fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz
    A_raw (T x 13) vbICA component amplitudes.
/home/jovyan/cascadia-sse-strain  (mesh data/cascadia_mesh.txt; sse_strain.load_fault_ascii,
    load_solution, Solution.component_slip_patterns -> Ld (742 x 13)).

EXACT DEFINITIONS
-----------------
Epochs           E1 2010-01-01..2014-12-31 (5.0 yr), E2 2015-01-01..2019-12-31 (5.0 yr),
                 E3 2020-01-01..2026-05-27 (6.4 yr, tremor-catalogue end).  Both N and S_mm are
                 clipped to the SAME epoch bounds, so the efficiency ratio S/N is per-detection
                 and the unequal epoch length mostly cancels (see caveats in NOTES.md).
Tremor->triangle nearest triangle centroid in a local km frame
                 x = (lon + 124) * 111 * cos(47 deg), y = (lat - 47) * 111  (scipy cKDTree),
                 assignments farther than 40 km dropped.
N[tri, e]        tremor detections assigned to that triangle inside epoch e.
Daily slip rate  raw_cum_mm = -1e3 * (Ld @ A_raw.T)   (742 x T, mm, thrust-positive)
                 rate = np.diff(raw_cum_mm, axis=1) (mm/d on days[1:]), then a 7-day CENTRED
                 moving average, then pos_rate = np.clip(rate, 0, None).
S_mm[tri, e]     sum of pos_rate over the days of epoch e (mm of thrust slip).
eff[tri, e]      S_mm / N, mm of slip per tremor detection.  DEFINED only where N >= 50 and
                 S_mm >= 5; NaN otherwise.
Normalisation    tot[e] = detections assigned to triangles with centroid lat >= 44N in epoch e.
                 f[e]   = tot[e] / mean(tot[E1], tot[E2], tot[E3]).
                 N_norm = N / f[e];  eff_norm = S_mm / N_norm = eff * f[e].
                 This removes the first-order network-growth change in detection RATE (it does
                 NOT remove a spatially varying change in detectability).
Change maps      log10(eff_E3 / eff_E1) and log10(eff_E2 / eff_E1), triangles defined in both
                 epochs; the same two ratios for eff_norm.
Profiles         median eff (and eff_norm) over the defined triangles in each 0.5-deg latitude
                 bin, per epoch, with a 90 % bootstrap interval (2000 resamples of the triangles
                 in the bin, seed 20260910; bins with < 3 defined triangles are NaN).
Tests            (a) Mann-Whitney U (two-sided) on log10(eff_E3/eff_E1) between north (centroid
                     lat >= 46) and central (44 <= lat < 46) triangles;
                 (b) Wilcoxon signed-rank (two-sided) of the margin-wide (lat >= 44)
                     log10(eff_E3/eff_E1) against 0; same for E2/E1.
                 Run SEPARATELY on the raw and on the normalised log-ratio.  The Mann-Whitney
                 (a) is identical for the two (both groups get the same constant offset
                 log10(f3/f1)), but the one-sample Wilcoxon (b) is NOT shift-invariant, so its
                 statistic and p differ between raw and normalised -- both are reported.

OUTPUTS (inside tremor_slip_geometry/)
--------------------------------------
results/s3_efficiency_triangles.csv  tri, lat, lon, depth_km, area_km2, seg,
                                     N_E1..E3, S_mm_E1..E3, eff_E1..E3, eff_norm_E1..E3,
                                     log10_eff_ratio_E3_E1, log10_eff_ratio_E2_E1,
                                     log10_effnorm_ratio_E3_E1, log10_effnorm_ratio_E2_E1
results/s3_profiles.csv              lat_bin_centre, epoch, n_tri, median_eff, eff_lo90, eff_hi90,
                                     median_eff_norm, eff_norm_lo90, eff_norm_hi90
results/s3_tests.csv                 test, statistic, p, n, medians
figures/s3_efficiency_maps.png       3 epoch maps of log10(eff_norm) as coloured triangle
                                     centroids over the grey mesh, the E3/E1 change map, and the
                                     along-strike profile panels (raw and normalised).

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s3_efficiency_maps.py
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path("/home/jovyan/tremorferometry")
OUT = ROOT / "tremor_slip_geometry"
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
STRAIN = Path("/home/jovyan/cascadia-sse-strain")

MAX_ASSIGN_KM = 40.0
MIN_N = 50
MIN_S_MM = 5.0
NBOOT = 2000
SEED = 20260910
LAT_BIN = 0.5
EPOCHS = [("E1", "2010-01-01", "2014-12-31"),
          ("E2", "2015-01-01", "2019-12-31"),
          ("E3", "2020-01-01", "2026-05-27")]


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

    # ---------------------------------------------------------------- mesh
    sys.path.insert(0, str(STRAIN / "src"))
    from sse_strain import load_fault_ascii, load_solution

    fault = load_fault_ascii(str(STRAIN / "data/cascadia_mesh.txt"))
    vlon, vlat = np.asarray(fault.lon), np.asarray(fault.lat)   # (742, 3) triangle vertices

    pz = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
    plat, plon = pz["patch_lat"], pz["patch_lon"]
    pdepth, parea = pz["patch_depth_km"], pz["patch_area_km2"]
    days = pd.to_datetime(pd.Series(pz["days"]).astype(str)).values
    ntri = len(plat)
    tri_seg = seg_of_lat(plat)

    # ---------------------------------------------------------------- daily positive slip rate
    sol = load_solution(str(STRAIN / "data/2026-08-18"), fault)
    sol.invert(verbose=False)
    _, Ld = sol.component_slip_patterns()
    A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz",
                                  allow_pickle=True)["A_raw"])
    raw_cum_mm = -1e3 * (Ld @ A_raw.T)
    rate = moving_average_centred(np.diff(raw_cum_mm, axis=1), 7, axis=1)
    pos_rate = np.clip(rate, 0.0, None)
    rate_days = pd.DatetimeIndex(days[1:])
    print(f"[mesh] {ntri} triangles; [slip] pos_rate {pos_rate.shape}, "
          f"{rate_days[0].date()}..{rate_days[-1].date()}")

    # ---------------------------------------------------------------- tremor
    tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                     usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr["time"])
    tree = cKDTree(km_frame(plon, plat))
    dist, idx = tree.query(km_frame(tr["lon"].values, tr["lat"].values), k=1)
    keep = dist <= MAX_ASSIGN_KM
    tr = tr.loc[keep].copy()
    tr["tri"] = idx[keep]
    print(f"[tremor] {int(keep.sum())} / {len(keep)} detections assigned within "
          f"{MAX_ASSIGN_KM:.0f} km")

    # ---------------------------------------------------------------- per-epoch N and S_mm
    N = np.zeros((ntri, 3), dtype=int)
    S = np.zeros((ntri, 3), dtype=float)
    tot = np.zeros(3, dtype=float)
    north_of_44 = plat >= 44.0
    ttime = tr["time"].values
    ttri = tr["tri"].values
    for k, (name, a, b) in enumerate(EPOCHS):
        a64, b64 = np.datetime64(a), np.datetime64(b) + np.timedelta64(1, "D")
        m = (ttime >= a64) & (ttime < b64)
        N[:, k] = np.bincount(ttri[m], minlength=ntri)
        dm = (rate_days >= pd.Timestamp(a)) & (rate_days <= pd.Timestamp(b))
        S[:, k] = pos_rate[:, np.asarray(dm)].sum(axis=1)
        tot[k] = N[north_of_44, k].sum()
        yrs = (pd.Timestamp(b) - pd.Timestamp(a)).days / 365.25
        print(f"[epoch {name}] {a}..{b} ({yrs:.2f} yr)  detections(all mesh) "
              f"{int(N[:, k].sum()):,}  detections(lat>=44N) {int(tot[k]):,}  "
              f"slip sum {S[:, k].sum():.0f} mm-triangle  slip days {int(np.sum(dm))}")

    fnorm = tot / tot.mean()
    print("\n[normalisation] f[e] = tot[e] / mean(tot) = " +
          ", ".join(f"{EPOCHS[k][0]}={fnorm[k]:.4f}" for k in range(3)) +
          f"   (mean tot = {tot.mean():,.0f})")

    eff = np.full((ntri, 3), np.nan)
    ok = (N >= MIN_N) & (S >= MIN_S_MM)
    eff[ok] = S[ok] / N[ok]
    eff_norm = eff * fnorm[None, :]
    for k, (name, _, _) in enumerate(EPOCHS):
        print(f"[defined] {name}: {int(ok[:, k].sum())} / {ntri} triangles pass "
              f"N>={MIN_N} & S>={MIN_S_MM} mm")

    with np.errstate(divide="ignore", invalid="ignore"):
        lr31 = np.log10(eff[:, 2] / eff[:, 0])
        lr21 = np.log10(eff[:, 1] / eff[:, 0])
        lrn31 = np.log10(eff_norm[:, 2] / eff_norm[:, 0])
        lrn21 = np.log10(eff_norm[:, 1] / eff_norm[:, 0])

    tri_df = pd.DataFrame(dict(
        tri=np.arange(ntri), lat=plat, lon=plon, depth_km=pdepth, area_km2=parea, seg=tri_seg,
        N_E1=N[:, 0], N_E2=N[:, 1], N_E3=N[:, 2],
        S_mm_E1=S[:, 0], S_mm_E2=S[:, 1], S_mm_E3=S[:, 2],
        eff_E1=eff[:, 0], eff_E2=eff[:, 1], eff_E3=eff[:, 2],
        eff_norm_E1=eff_norm[:, 0], eff_norm_E2=eff_norm[:, 1], eff_norm_E3=eff_norm[:, 2],
        log10_eff_ratio_E3_E1=lr31, log10_eff_ratio_E2_E1=lr21,
        log10_effnorm_ratio_E3_E1=lrn31, log10_effnorm_ratio_E2_E1=lrn21))
    tri_df.to_csv(OUT / "results/s3_efficiency_triangles.csv", index=False)
    print(f"\n[write] results/s3_efficiency_triangles.csv ({ntri} rows; "
          f"{int(np.isfinite(lr31).sum())} triangles with an E3/E1 ratio)")

    # ---------------------------------------------------------------- along-strike profiles
    edges = np.arange(np.floor(plat.min() / LAT_BIN) * LAT_BIN,
                      np.ceil(plat.max() / LAT_BIN) * LAT_BIN + LAT_BIN, LAT_BIN)
    prows = []
    for i in range(len(edges) - 1):
        c = 0.5 * (edges[i] + edges[i + 1])
        inb = (plat >= edges[i]) & (plat < edges[i + 1])
        for k, (name, _, _) in enumerate(EPOCHS):
            m, lo, hi, n = boot_median_ci(eff[inb, k], rng)
            mn, lon_, hin, _ = boot_median_ci(eff_norm[inb, k], rng)
            prows.append(dict(lat_bin_centre=c, epoch=name, n_tri=n,
                              median_eff=m, eff_lo90=lo, eff_hi90=hi,
                              median_eff_norm=mn, eff_norm_lo90=lon_, eff_norm_hi90=hin))
    prof = pd.DataFrame(prows)
    prof.to_csv(OUT / "results/s3_profiles.csv", index=False)
    print(f"[write] results/s3_profiles.csv ({len(prof)} rows)")

    print("\n=== along-strike profile: median efficiency (mm slip per detection) ===")
    piv = prof.pivot(index="lat_bin_centre", columns="epoch", values="median_eff")
    pivn = prof.pivot(index="lat_bin_centre", columns="epoch", values="median_eff_norm")
    ntr = prof.pivot(index="lat_bin_centre", columns="epoch", values="n_tri")
    tbl = pd.concat([ntr.add_prefix("nTri_"), piv.add_prefix("eff_"),
                     pivn.add_prefix("effnorm_")], axis=1).round(4)
    print(tbl.to_string())

    print("\n=== margin-wide (lat>=44N) medians of the DEFINED triangles ===")
    for k, (name, _, _) in enumerate(EPOCHS):
        sel = north_of_44 & ok[:, k]
        print(f"{name}: n {int(sel.sum()):3d}  median eff {np.median(eff[sel, k]):.4f}  "
              f"median eff_norm {np.median(eff_norm[sel, k]):.4f} mm/detection")

    # ---------------------------------------------------------------- tests
    trows = []
    both31 = np.isfinite(lr31) & north_of_44
    n_mask = both31 & (tri_seg == "north")
    c_mask = both31 & (tri_seg == "central")
    print("\n=== tests on log10(eff_E3 / eff_E1) ===")
    if n_mask.sum() >= 3 and c_mask.sum() >= 3:
        u = mannwhitneyu(lr31[n_mask], lr31[c_mask], alternative="two-sided")
        print(f"north vs central (Mann-Whitney): n_north {int(n_mask.sum())} "
              f"(median {np.median(lr31[n_mask]):+.4f}), n_central {int(c_mask.sum())} "
              f"(median {np.median(lr31[c_mask]):+.4f})  U {u.statistic:.1f}  p {u.pvalue:.4g}")
        trows.append(dict(test="MWU north-vs-central log10(eff_E3/eff_E1) "
                               "[shift-invariant: same for raw and normalised]",
                          statistic=float(u.statistic), p=float(u.pvalue),
                          n=f"{int(n_mask.sum())} vs {int(c_mask.sum())}",
                          medians=f"north {np.median(lr31[n_mask]):+.4f}, "
                                  f"central {np.median(lr31[c_mask]):+.4f}"))
    else:
        print(f"north vs central: too few triangles ({int(n_mask.sum())}, {int(c_mask.sum())})")

    for lab, lr, lrn in (("E3/E1", lr31, lrn31), ("E2/E1", lr21, lrn21)):
        sel = np.isfinite(lr) & north_of_44
        if sel.sum() >= 6:
            # the one-sample signed-rank test is NOT shift-invariant -> test both versions
            for kind, x in (("raw", lr[sel]), ("norm", lrn[sel])):
                w = wilcoxon(x, alternative="two-sided")
                print(f"margin-wide {lab} vs 0 (Wilcoxon, {kind:4s}): n {int(sel.sum())}  "
                      f"median {np.median(x):+.4f} (x{10 ** np.median(x):.3f})  "
                      f"W {w.statistic:.1f}  p {w.pvalue:.4g}")
                trows.append(dict(test=f"Wilcoxon margin-wide log10(eff_{kind}_{lab}) vs 0",
                                  statistic=float(w.statistic), p=float(w.pvalue),
                                  n=int(sel.sum()),
                                  medians=f"median {np.median(x):+.4f} "
                                          f"(x{10 ** np.median(x):.3f})"))
    pd.DataFrame(trows).to_csv(OUT / "results/s3_tests.csv", index=False)
    print("[write] results/s3_tests.csv")

    # ---------------------------------------------------------------- figure
    verts = [np.column_stack([vlon[i], vlat[i]]) for i in range(ntri)]
    fig = plt.figure(figsize=(19, 10.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[2.6, 1.15], hspace=0.22, wspace=0.28)

    vals = eff_norm[np.isfinite(eff_norm)]
    lv = np.log10(vals)
    vmin, vmax = np.percentile(lv, 2), np.percentile(lv, 98)
    axes_map = []
    for k, (name, a, b) in enumerate(EPOCHS):
        ax = fig.add_subplot(gs[0, k])
        ax.add_collection(PolyCollection(verts, facecolors="none", edgecolors="0.82", lw=0.3))
        m = np.isfinite(eff_norm[:, k])
        sc = ax.scatter(plon[m], plat[m], c=np.log10(eff_norm[m, k]), s=26, cmap="viridis",
                        vmin=vmin, vmax=vmax, edgecolors="none")
        ax.set_title(f"{name}  {a[:7]}..{b[:7]}\nlog10 eff_norm, {int(m.sum())} tri "
                     f"(f={fnorm[k]:.3f})", fontsize=9)
        axes_map.append((ax, sc))
    cb = fig.colorbar(axes_map[-1][1], ax=[a for a, _ in axes_map], fraction=0.03, pad=0.01)
    cb.set_label("log10 (mm slip per normalised detection)", fontsize=8)

    ax = fig.add_subplot(gs[0, 3])
    ax.add_collection(PolyCollection(verts, facecolors="none", edgecolors="0.82", lw=0.3))
    m = np.isfinite(lrn31)
    lim = float(np.nanpercentile(np.abs(lrn31[m]), 95)) if m.sum() else 1.0
    sc = ax.scatter(plon[m], plat[m], c=lrn31[m], s=26, cmap="RdBu_r", vmin=-lim, vmax=lim,
                    edgecolors="none")
    m44 = m & north_of_44
    ax.set_title(f"change  log10(eff_norm E3 / E1)\n{int(m.sum())} tri, median "
                 f"{np.nanmedian(lrn31[m]):+.3f} (lat>=44N: {np.nanmedian(lrn31[m44]):+.3f})",
                 fontsize=9)
    fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01).set_label("log10 ratio", fontsize=8)

    for a_, _ in axes_map:
        pass
    for a_ in [a for a, _ in axes_map] + [ax]:
        a_.set_xlim(plon.min() - 0.4, plon.max() + 0.4)
        a_.set_ylim(plat.min() - 0.4, plat.max() + 0.4)
        a_.set_aspect(1.0 / np.cos(np.deg2rad(47.0)))
        a_.axhline(46, color="k", lw=0.6, ls="--", alpha=0.6)
        a_.axhline(44, color="k", lw=0.6, ls="--", alpha=0.6)
        a_.set_xlabel("lon", fontsize=8)
        a_.tick_params(labelsize=7)
    axes_map[0][0].set_ylabel("lat", fontsize=8)

    colors = {"E1": "tab:blue", "E2": "tab:orange", "E3": "tab:red"}
    for j, (col, lo, hi, lab) in enumerate([
            ("median_eff", "eff_lo90", "eff_hi90", "raw efficiency"),
            ("median_eff_norm", "eff_norm_lo90", "eff_norm_hi90", "network-normalised")]):
        axp = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
        for name, _, _ in EPOCHS:
            s = prof[prof.epoch == name].sort_values("lat_bin_centre")
            axp.fill_between(s.lat_bin_centre, s[lo], s[hi], color=colors[name], alpha=0.18)
            axp.plot(s.lat_bin_centre, s[col], "o-", color=colors[name], ms=3.5, lw=1.6,
                     label=name)
        axp.set_yscale("log")
        axp.axvline(46, color="k", lw=0.6, ls="--", alpha=0.6)
        axp.axvline(44, color="k", lw=0.6, ls="--", alpha=0.6)
        axp.set_xlabel("latitude of 0.5-deg bin centre")
        axp.set_ylabel("median mm slip / detection")
        axp.set_title(f"along-strike profile — {lab} (90 % bootstrap over triangles)",
                      fontsize=9)
        axp.grid(alpha=0.3)
        axp.legend(fontsize=8)

    fig.suptitle("S3 — tremor efficiency (modelled thrust slip per PNSN detection) per triangle "
                 "and epoch", fontsize=13)
    fig.savefig(OUT / "figures/s3_efficiency_maps.png", dpi=150, bbox_inches="tight")
    print("[write] figures/s3_efficiency_maps.png")


if __name__ == "__main__":
    main()
