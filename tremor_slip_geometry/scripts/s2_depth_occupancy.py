#!/usr/bin/env python
"""
S2 — depth occupancy through time: where on the interface (in Slab2 depth) do tremor and
modelled thrust slip sit, per year and per segment, 2010-2026, and does either drift?

PLAN.md step S2 (Q4 / H2).  No dv/v anywhere in this study.

INPUTS (all read-only, outside this folder)
-------------------------------------------
catalogs/pnsn_tremor_cascadia_full.csv
    PNSN tremor master catalogue: time, lat, lon (the `depth` column is a PLACEHOLDER and is
    NOT used).  748k rows, 2010-01-06 .. 2026-05-27.
fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
    patch_lat, patch_lon, patch_depth_km, patch_area_km2  (742 mesh triangles), days (T strings,
    daily 2007-01-02..2026-08-18).
fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv
    per-patch tremor confirmation (catalog, id, start, end, seg, Mw, broad_pulse, confirmed).
fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz
    A_raw (T x 13) vbICA component amplitudes.
/home/jovyan/cascadia-sse-strain  (mesh data/cascadia_mesh.txt; sse_strain.load_fault_ascii,
    load_solution, Solution.component_slip_patterns -> Ld (742 x 13);
    sse_strain.megathrust.slab_depth_interpolator).

EXACT DEFINITIONS
-----------------
Slab2 depth      f(lon, lat) from slab_depth_interpolator(fault), positive down, km.  Defined only
                 inside the 10-60 km mesh band -> ~94 % of tremor detections get a finite depth;
                 detections outside the band are counted in N but excluded from depth percentiles.
Tremor->triangle each detection is assigned to the NEAREST triangle centroid in a local km frame
                 x = (lon + 124) * 111 * cos(47 deg), y = (lat - 47) * 111  (scipy cKDTree);
                 assignments farther than 40 km are dropped.
Segment          of a triangle (and hence of a detection, via its triangle): north lat >= 46,
                 central 44 <= lat < 46, south lat < 44, by triangle CENTROID latitude.
CONFIRMED patch  tremor_vs_slip_patches.csv rows with catalog == 'thrust' & confirmed &
                 seg in {north, central} & Mw >= 6.0 & ~broad_pulse   (27 patches: 21 N, 6 C).
DUR window       [start, end] of a confirmed patch, inclusive, whole days.
Daily slip rate  raw_cum_mm = -1e3 * (Ld @ A_raw.T)          (742 x T, mm, thrust-positive)
                 rate       = np.diff(raw_cum_mm, axis=1)     (mm/d, valid on days[1:])
                 rate       = 7-day CENTRED moving average of that (edges: shorter window)
                 pos_rate   = np.clip(rate, 0, None)          (only slip in the thrust sense)
Yearly slip      acc_mm[tri, y] = sum of pos_rate over the days of year y.
Slip weight      w[tri, y] = patch_area_km2[tri] * acc_mm[tri, y]   (area x accumulated slip).
slip_depth       sum(w * depth) / sum(w) over the triangles of the segment, for that year.
slip_edge10/90   10th / 90th percentile of the depth distribution weighted by w, i.e. the depths
                 at which the cumulative w-vs-depth curve reaches 10 % / 90 % (linear interp).
completeness     N(detections, seg, year) / mean over 2010-2026 of N(detections, seg, year).
                 A crude proxy for network growth only -- 2026 is PARTIAL (tremor ends 2026-05-27)
                 and is flagged; the mean includes it, which drags the mean slightly down.
Theil-Sen trend  scipy.stats.theilslopes(value, year) on the yearly series, units km/yr.
                 Permutation p: 2000 random shuffles of the YEAR labels, two-sided
                 p = (1 + #{|slope_perm| >= |slope_obs|}) / (1 + nperm).  Spearman rho reported too.
                 Run twice: year_set '2010-2025' (complete years, PRIMARY) and '2010-2026'.
Half-record test tremor depth samples 2010-2017 vs 2018-2026 per segment; the larger half is
                 randomly subsampled to the size of the smaller (200 draws, seed 20260910) to
                 remove the completeness bias; reported p is the MEDIAN Mann-Whitney (two-sided)
                 p over the draws, alongside the full-sample p.

OUTPUTS (inside tremor_slip_geometry/)
--------------------------------------
results/s2_yearly.csv   year, seg, partial, n_tremor, n_tremor_depth, p10, p50, p90,
                        n_dur, n_dur_depth, dur_p10, dur_p50, dur_p90,
                        slip_depth, slip_edge10, slip_edge90, slip_sum_mm, completeness
results/s2_trends.csv   seg, quantity, year_set, n_years, theilsen_slope_km_per_yr, lo, hi,
                        p_perm, spearman_rho, spearman_p
figures/s2_depth_occupancy.png  per segment: yearly tremor p10/p50/p90 band + slip depth/edges
                        vs year (row 1, depth axis inverted); completeness (row 2).

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s2_depth_occupancy.py
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
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, theilslopes

ROOT = Path("/home/jovyan/tremorferometry")
OUT = ROOT / "tremor_slip_geometry"
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
STRAIN = Path("/home/jovyan/cascadia-sse-strain")

MAX_ASSIGN_KM = 40.0
NPERM = 2000
NSUB = 200
SEED = 20260910
YEARS = list(range(2010, 2027))
SEGS = ["north", "central", "south"]


def km_frame(lon, lat):
    """Local km frame used for the nearest-triangle assignment."""
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
    """Centred moving average; window shrinks at the edges (no wrap, no zero padding)."""
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


def weighted_quantile(values, weights, q):
    """q-quantile of `values` under `weights` (linear interp on the cumulative weight curve)."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if ok.sum() == 0:
        return np.nan
    v, w = values[ok], weights[ok]
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    cw = (cw - 0.5 * w) / cw[-1]  # mid-point rule
    return float(np.interp(q, cw, v))


def theilsen_with_perm(years, vals, rng, nperm=NPERM):
    years = np.asarray(years, dtype=float)
    vals = np.asarray(vals, dtype=float)
    ok = np.isfinite(years) & np.isfinite(vals)
    years, vals = years[ok], vals[ok]
    if len(years) < 4:
        return dict(n=len(years), slope=np.nan, lo=np.nan, hi=np.nan, p_perm=np.nan,
                    rho=np.nan, p_rho=np.nan)
    slope, icept, lo, hi = theilslopes(vals, years, 0.90)
    perm = np.empty(nperm)
    for i in range(nperm):
        perm[i] = theilslopes(vals, rng.permutation(years), 0.90)[0]
    p_perm = (1 + int(np.sum(np.abs(perm) >= abs(slope)))) / (1 + nperm)
    rho, p_rho = spearmanr(years, vals)
    return dict(n=len(years), slope=float(slope), lo=float(lo), hi=float(hi),
                p_perm=float(p_perm), rho=float(rho), p_rho=float(p_rho))


def main():
    rng = np.random.default_rng(SEED)
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- mesh + Slab2 depth
    sys.path.insert(0, str(STRAIN / "src"))
    from sse_strain import load_fault_ascii, load_solution
    from sse_strain.megathrust import slab_depth_interpolator

    fault = load_fault_ascii(str(STRAIN / "data/cascadia_mesh.txt"))
    fdep = slab_depth_interpolator(fault)

    pz = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
    plat, plon = pz["patch_lat"], pz["patch_lon"]
    pdepth, parea = pz["patch_depth_km"], pz["patch_area_km2"]
    days = pd.to_datetime(pd.Series(pz["days"]).astype(str)).values
    ntri = len(plat)
    tri_seg = seg_of_lat(plat)
    print(f"[mesh] {ntri} triangles, depth {pdepth.min():.1f}-{pdepth.max():.1f} km, "
          f"segments: " + ", ".join(f"{s}={int((tri_seg == s).sum())}" for s in SEGS))

    # ---------------------------------------------------------------- daily positive slip rate
    sol = load_solution(str(STRAIN / "data/2026-08-18"), fault)
    sol.invert(verbose=False)
    _, Ld = sol.component_slip_patterns()
    A_raw = np.nan_to_num(np.load(SSE / "components_at_cells_trailing30.npz",
                                  allow_pickle=True)["A_raw"])
    raw_cum_mm = -1e3 * (Ld @ A_raw.T)                    # (742, T) mm, thrust-positive
    rate = np.diff(raw_cum_mm, axis=1)                    # (742, T-1) mm/d, on days[1:]
    rate = moving_average_centred(rate, 7, axis=1)
    pos_rate = np.clip(rate, 0.0, None)
    rate_days = pd.DatetimeIndex(days[1:])
    rate_year = rate_days.year.values
    print(f"[slip] pos_rate {pos_rate.shape}, days {rate_days[0].date()}..{rate_days[-1].date()}, "
          f"margin total positive slip {pos_rate.sum():.0f} mm-triangle")

    # ---------------------------------------------------------------- tremor
    tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv",
                     usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr["time"])
    tree = cKDTree(km_frame(plon, plat))
    dist, idx = tree.query(km_frame(tr["lon"].values, tr["lat"].values), k=1)
    keep = dist <= MAX_ASSIGN_KM
    print(f"[tremor] {len(tr)} detections, {keep.sum()} ({100 * keep.mean():.1f} %) within "
          f"{MAX_ASSIGN_KM:.0f} km of a triangle centroid")
    tr = tr.loc[keep].copy()
    tr["tri"] = idx[keep]
    tr["dist_km"] = dist[keep]
    tr["seg"] = tri_seg[tr["tri"].values]
    tr["depth_km"] = fdep(tr["lon"].values, tr["lat"].values)
    tr["year"] = tr["time"].dt.year
    fin = np.isfinite(tr["depth_km"].values)
    print(f"[tremor] Slab2 depth finite for {100 * fin.mean():.1f} % of assigned detections "
          f"(mesh band 10-60 km)")

    # ---------------------------------------------------------------- DUR windows
    conf = pd.read_csv(SSE / "tremor_vs_slip_patches.csv")
    m = ((conf["catalog"] == "thrust") & conf["confirmed"].astype(bool)
         & conf["seg"].isin(["north", "central"]) & (conf["Mw"] >= 6.0)
         & (~conf["broad_pulse"].astype(bool)))
    conf = conf.loc[m].copy()
    conf["start"] = pd.to_datetime(conf["start"])
    conf["end"] = pd.to_datetime(conf["end"])
    print(f"[patches] {len(conf)} confirmed patches: " +
          ", ".join(f"{s}={int((conf['seg'] == s).sum())}" for s in SEGS))

    in_dur = np.zeros(len(tr), dtype=bool)
    tseg = tr["seg"].values
    ttime = tr["time"].values
    for _, r in conf.iterrows():
        sel = ((tseg == r["seg"]) & (ttime >= np.datetime64(r["start"]))
               & (ttime <= np.datetime64(r["end"]) + np.timedelta64(1, "D")))
        in_dur |= sel
    tr["in_dur"] = in_dur
    print(f"[patches] {in_dur.sum()} detections inside a same-segment confirmed DUR window")

    # ---------------------------------------------------------------- yearly table
    rows = []
    for seg in SEGS:
        tseg_mask = tr["seg"].values == seg
        tri_mask = tri_seg == seg
        for y in YEARS:
            ym = tseg_mask & (tr["year"].values == y)
            d = tr["depth_km"].values[ym]
            d = d[np.isfinite(d)]
            dm = ym & tr["in_dur"].values
            dd = tr["depth_km"].values[dm]
            dd = dd[np.isfinite(dd)]
            acc = pos_rate[:, rate_year == y].sum(axis=1)
            w = parea * acc
            w = np.where(tri_mask, w, 0.0)
            rows.append(dict(
                year=y, seg=seg, partial=(y == 2026),
                n_tremor=int(ym.sum()), n_tremor_depth=int(len(d)),
                p10=np.percentile(d, 10) if len(d) else np.nan,
                p50=np.percentile(d, 50) if len(d) else np.nan,
                p90=np.percentile(d, 90) if len(d) else np.nan,
                n_dur=int(dm.sum()), n_dur_depth=int(len(dd)),
                dur_p10=np.percentile(dd, 10) if len(dd) else np.nan,
                dur_p50=np.percentile(dd, 50) if len(dd) else np.nan,
                dur_p90=np.percentile(dd, 90) if len(dd) else np.nan,
                slip_depth=(float(np.sum(w * pdepth) / np.sum(w)) if np.sum(w) > 0 else np.nan),
                slip_edge10=weighted_quantile(pdepth, w, 0.10),
                slip_edge90=weighted_quantile(pdepth, w, 0.90),
                slip_sum_mm=float(acc[tri_mask].sum()),
            ))
    yearly = pd.DataFrame(rows)
    yearly["completeness"] = yearly.groupby("seg")["n_tremor"].transform(lambda s: s / s.mean())
    yearly.to_csv(OUT / "results/s2_yearly.csv", index=False)
    print(f"\n[write] results/s2_yearly.csv  ({len(yearly)} rows)")

    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print("\n=== s2_yearly (rounded) ===")
        show = yearly[["year", "seg", "n_tremor", "p10", "p50", "p90", "dur_p50",
                       "slip_depth", "slip_edge10", "slip_edge90", "completeness"]].round(2)
        for seg in SEGS:
            print(f"\n-- {seg} --")
            print(show[show.seg == seg].to_string(index=False))

    # ---------------------------------------------------------------- trends
    quantities = ["p10", "p50", "p90", "dur_p50", "slip_depth", "slip_edge10", "slip_edge90"]
    trows = []
    for seg in SEGS:
        sub = yearly[yearly.seg == seg]
        for yset, ymax in (("2010-2025", 2025), ("2010-2026", 2026)):
            s = sub[sub.year <= ymax]
            for q in quantities:
                r = theilsen_with_perm(s["year"].values, s[q].values, rng)
                trows.append(dict(seg=seg, quantity=q, year_set=yset, n_years=r["n"],
                                  theilsen_slope_km_per_yr=r["slope"], lo90=r["lo"], hi90=r["hi"],
                                  p_perm=r["p_perm"], spearman_rho=r["rho"],
                                  spearman_p=r["p_rho"]))
    trends = pd.DataFrame(trows)
    trends.to_csv(OUT / "results/s2_trends.csv", index=False)
    print(f"\n[write] results/s2_trends.csv  ({len(trends)} rows)")
    print("\n=== s2_trends (PRIMARY year_set 2010-2025; slope km/yr, + = deepening) ===")
    pr = trends[trends.year_set == "2010-2025"].copy()
    pr[["theilsen_slope_km_per_yr", "lo90", "hi90", "spearman_rho"]] = pr[
        ["theilsen_slope_km_per_yr", "lo90", "hi90", "spearman_rho"]].round(3)
    pr[["p_perm", "spearman_p"]] = pr[["p_perm", "spearman_p"]].round(4)
    print(pr.to_string(index=False))
    print("\n=== s2_trends (2010-2026, incl. partial 2026) ===")
    pr2 = trends[trends.year_set == "2010-2026"].copy()
    pr2[["theilsen_slope_km_per_yr", "spearman_rho"]] = pr2[
        ["theilsen_slope_km_per_yr", "spearman_rho"]].round(3)
    pr2[["p_perm", "spearman_p"]] = pr2[["p_perm", "spearman_p"]].round(4)
    print(pr2[["seg", "quantity", "theilsen_slope_km_per_yr", "p_perm",
               "spearman_rho", "spearman_p"]].to_string(index=False))
    nq = len(quantities)
    print(f"\n[multiple testing] {nq} quantities x 3 segments = {3 * nq} tests per year_set; "
          f"Bonferroni FWER 0.05 threshold p < {0.05 / (3 * nq):.4f}")

    # ---------------------------------------------------------------- half-record comparison
    print("\n=== half-record tremor depth: 2010-2017 vs 2018-2026 "
          "(equal-size subsampled Mann-Whitney) ===")
    hrows = []
    for seg in SEGS:
        s = tr[(tr.seg == seg) & np.isfinite(tr.depth_km.values)]
        a = s.loc[s.year <= 2017, "depth_km"].values
        b = s.loc[s.year >= 2018, "depth_km"].values
        if len(a) < 50 or len(b) < 50:
            print(f"{seg}: too few detections ({len(a)}, {len(b)})")
            continue
        p_full = mannwhitneyu(a, b, alternative="two-sided").pvalue
        n = min(len(a), len(b))
        ps, dmed = [], []
        for _ in range(NSUB):
            aa = rng.choice(a, n, replace=False) if len(a) > n else a
            bb = rng.choice(b, n, replace=False) if len(b) > n else b
            ps.append(mannwhitneyu(aa, bb, alternative="two-sided").pvalue)
            dmed.append(np.median(bb) - np.median(aa))
        hrows.append(dict(seg=seg, n_first=len(a), n_second=len(b), n_sub=n,
                          p10_first=np.percentile(a, 10), p50_first=np.percentile(a, 50),
                          p90_first=np.percentile(a, 90),
                          p10_second=np.percentile(b, 10), p50_second=np.percentile(b, 50),
                          p90_second=np.percentile(b, 90),
                          d_p50_km=np.median(b) - np.median(a),
                          d_p50_sub_med_km=float(np.median(dmed)),
                          p_full=float(p_full), p_sub_median=float(np.median(ps))))
        h = hrows[-1]
        print(f"{seg:8s} n {len(a):6d} vs {len(b):6d} (sub {n})  "
              f"p10 {h['p10_first']:5.1f}->{h['p10_second']:5.1f}  "
              f"p50 {h['p50_first']:5.1f}->{h['p50_second']:5.1f}  "
              f"p90 {h['p90_first']:5.1f}->{h['p90_second']:5.1f} km  "
              f"dp50 {h['d_p50_km']:+5.2f} km  p_full {p_full:.2e}  "
              f"p_sub(med) {np.median(ps):.2e}")
    if hrows:
        pd.DataFrame(hrows).to_csv(OUT / "results/s2_half_record.csv", index=False)
        print("[write] results/s2_half_record.csv")

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=True,
                             gridspec_kw=dict(height_ratios=[2.4, 1]))
    for j, seg in enumerate(SEGS):
        s = yearly[yearly.seg == seg].sort_values("year")
        ax = axes[0, j]
        ax.fill_between(s.year, s.p10, s.p90, color="tab:blue", alpha=0.18,
                        label="tremor p10-p90")
        ax.plot(s.year, s.p50, "o-", color="tab:blue", lw=2, ms=4, label="tremor p50")
        ax.plot(s.year, s.p10, "-", color="tab:blue", lw=0.8, alpha=0.7)
        ax.plot(s.year, s.p90, "-", color="tab:blue", lw=0.8, alpha=0.7)
        ax.plot(s.year, s.dur_p50, "s--", color="tab:cyan", lw=1.2, ms=3.5,
                label="tremor p50 (DUR only)")
        ax.fill_between(s.year, s.slip_edge10, s.slip_edge90, color="tab:red", alpha=0.15,
                        label="slip 10-90 % weight")
        ax.plot(s.year, s.slip_depth, "o-", color="tab:red", lw=2, ms=4,
                label="slip-weighted depth")
        ax.plot(s.year, s.slip_edge10, "-", color="tab:red", lw=0.8, alpha=0.7)
        ax.plot(s.year, s.slip_edge90, "-", color="tab:red", lw=0.8, alpha=0.7)
        ax.axvspan(2025.5, 2026.5, color="0.85", zorder=0)
        ax.invert_yaxis()
        ax.set_title(f"{seg}  (n tremor {int(s.n_tremor.sum()):,})")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("Slab2 depth (km, positive down)")
            ax.legend(fontsize=7, loc="lower left", framealpha=0.9)
        ax2 = axes[1, j]
        ax2.bar(s.year, s.completeness, color="0.5", width=0.7)
        ax2.axhline(1.0, color="k", lw=0.8, ls=":")
        ax2.axvspan(2025.5, 2026.5, color="0.85", zorder=0)
        ax2.set_xlabel("year")
        ax2.set_xticks(range(2010, 2027, 2))
        ax2.grid(alpha=0.3, axis="y")
        if j == 0:
            ax2.set_ylabel("completeness proxy\n(N / segment mean)")
    fig.suptitle("S2 — tremor Slab2-depth occupancy vs modelled slip depth, per year and segment "
                 "(grey band: partial 2026)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT / "figures/s2_depth_occupancy.png", dpi=150)
    print("[write] figures/s2_depth_occupancy.png")


if __name__ == "__main__":
    main()
