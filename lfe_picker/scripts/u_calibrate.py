#!/usr/bin/env python
"""P4-0a step 1b - calibrate the site proxy against MEASURED single-event registration rho.

Aggregates results/p4/north_station_day_stats.csv to one row per north station:
  noise_<C>_db   = median 2-8 Hz velocity PSD over the QUIET days, in dB relative to CN.PGC on the
                   SAME days, for the vertical (Z), the horizontals (H = median of N and E) and 3-C;
  contrast_<C>_db= 10 log10 (median over tremor days / median over quiet days) of the incoherent
                   6-s window power, same component splits;
  dist_field_km  = mean HYPOCENTRAL distance to the 130 Bostock family locations (the source field
                   rho was measured against), which is the construction the Puget side mirrors.

Then fits rho/rho_PGC against (distance, noise, contrast), exactly as Merlin's P4-0a prescribes,
and reports Spearman with n and signs plus the LEAVE-ONE-OUT prediction error in rho/rho_PGC.
Accept: Spearman >= 0.6 with the expected signs (distance negative, noise negative, contrast
positive) and LOO median error <= 0.25.  Spearman < 0.4 -> the proxy is not predictive and Puget
must be selected on measured geometry + tremor contrast alone.

Writes results/p4/north_calibration.csv, results/p4/rho_calibration.json, figures/u_rho_calibration.png
"""
from __future__ import annotations
import os
import sys
import json
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402

MODELS = {
    "prescribed_3c":   ["dist_field_km", "noise_3c_db", "contrast_3c_db"],
    "vertical_terms":  ["dist_field_km", "noise_Z_db", "contrast_Z_db"],
    "horizontal_terms": ["dist_field_km", "noise_H_db", "contrast_H_db"],
    "dist_contrast":   ["dist_field_km", "contrast_3c_db"],
    "contrast_only":   ["contrast_3c_db"],
    "dist_only":       ["dist_field_km"],
    "noise_only":      ["noise_3c_db"],
}
PRIMARY = "prescribed_3c"


def field_distances(lat, lon):
    F = pd.read_csv(os.path.join(U.GEO, "north_families.csv")).drop_duplicates("bkey")
    d = U.hyp_km(F.lat.values[None, :], F.lon.values[None, :], F.depth_km.values[None, :],
                 np.asarray(lat, float)[:, None], np.asarray(lon, float)[:, None])
    return d.mean(1), d.min(1)


def aggregate(S, PD, ref="CN.PGC"):
    dt = dict(zip(PD.date, PD.day_type))
    S = S[S.status == "ok"].copy()
    S["day_type"] = S.date.map(dt)
    ref_ = S[S.netsta == ref].set_index("date")
    rows = []
    for ns, g in S.groupby("netsta"):
        g = g.set_index("date")
        q = [d for d in g.index if dt.get(d) == "quiet"]
        t = [d for d in g.index if dt.get(d) == "tremor"]
        r = dict(netsta=ns, n_quiet=len(q), n_tremor=len(t),
                 coverage_med=float(g.coverage.median()),
                 resp_mode=g.resp_mode.mode().iloc[0] if len(g) else "")
        for lab, cols in [("Z", ["psd_Z"]), ("H", ["psd_N", "psd_E"]),
                          ("3c", ["psd_3c"])]:
            qq = [d for d in q if d in ref_.index]
            num = np.nanmedian([np.nanmedian(np.atleast_1d(g.loc[d, cols].values)) for d in qq]) if qq else np.nan
            den = np.nanmedian([np.nanmedian(np.atleast_1d(ref_.loc[d, cols].values)) for d in qq]) if qq else np.nan
            r[f"noise_{lab}_db"] = float(10 * np.log10(num / den)) if np.isfinite(num) and np.isfinite(den) else np.nan
            # all-day version, identical to the Puget probe's psd_ratio_med
            alld = [d for d in g.index if d in ref_.index]
            n2 = np.nanmedian([np.nanmedian(np.atleast_1d(g.loc[d, cols].values)) for d in alld]) if alld else np.nan
            d2 = np.nanmedian([np.nanmedian(np.atleast_1d(ref_.loc[d, cols].values)) for d in alld]) if alld else np.nan
            r[f"noise_allday_{lab}_db"] = float(10 * np.log10(n2 / d2)) if np.isfinite(n2) and np.isfinite(d2) else np.nan
        for lab, cols in [("Z", ["e6_Z"]), ("H", ["e6_N", "e6_E"]), ("3c", ["e6_3c"])]:
            tt = np.nanmedian([np.nanmedian(np.atleast_1d(g.loc[d, cols].values)) for d in t]) if t else np.nan
            qq2 = np.nanmedian([np.nanmedian(np.atleast_1d(g.loc[d, cols].values)) for d in q]) if q else np.nan
            r[f"contrast_{lab}_db"] = float(10 * np.log10(tt / qq2)) if np.isfinite(tt) and np.isfinite(qq2) and qq2 > 0 else np.nan
        # per-tremor-day contrast against the quiet median (for the ">= +3 dB on >= 3 of 4 days" bar)
        base = np.nanmedian([g.loc[d, "e6_3c"] for d in q]) if q else np.nan
        for i, d in enumerate(sorted(t)):
            r[f"contrast_day{i+1}_db"] = float(10 * np.log10(g.loc[d, "e6_3c"] / base)) if np.isfinite(base) else np.nan
        r["n_tremor_days_ge3db"] = int(sum(1 for i in range(len(t))
                                           if np.isfinite(r.get(f"contrast_day{i+1}_db", np.nan))
                                           and r[f"contrast_day{i+1}_db"] >= 3.0))
        rows.append(r)
    return pd.DataFrame(rows)


def loo_fit(D, cols, y="rho_rel"):
    X = D[cols].values.astype(float)
    yv = D[y].values.astype(float)
    ok = np.isfinite(X).all(1) & np.isfinite(yv)
    X, yv, idx = X[ok], yv[ok], np.where(ok)[0]
    n = len(yv)
    A = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(A, yv, rcond=None)
    pred_in = A @ beta
    pred_loo = np.full(n, np.nan)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        b, *_ = np.linalg.lstsq(A[m], yv[m], rcond=None)
        pred_loo[i] = A[i] @ b
    rs_in, p_in = spearmanr(pred_in, yv)
    rs_loo, p_loo = spearmanr(pred_loo, yv)
    err = np.abs(pred_loo - yv)
    return dict(cols=cols, n=int(n), beta=[float(b) for b in beta],
                spearman_in=float(rs_in), p_in=float(p_in),
                spearman_loo=float(rs_loo), p_loo=float(p_loo),
                loo_err_median=float(np.median(err)), loo_err_p90=float(np.percentile(err, 90)),
                r2_in=float(1 - np.sum((pred_in - yv) ** 2) / np.sum((yv - yv.mean()) ** 2)),
                idx=[int(i) for i in idx],
                pred_loo=[float(v) for v in pred_loo], pred_in=[float(v) for v in pred_in])


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    S = pd.read_csv(os.path.join(U.P4, "north_station_day_stats.csv"))
    PD = pd.read_csv(os.path.join(U.P4, "north_probe_days.csv"))
    P = pd.read_csv(os.path.join(U.P4, "north_stations.csv"))
    A = aggregate(S, PD).merge(P, on="netsta", how="left")
    A["dist_field_km"], A["dist_field_min_km"] = field_distances(A.lat.values, A.lon.values)
    A = A.sort_values("rho_rel", ascending=False).reset_index(drop=True)
    A.to_csv(os.path.join(U.P4, "north_calibration.csv"), index=False)
    D = A[A.primary].reset_index(drop=True)                    # n=14, VGZ (derived) held out
    print(f"[cal] {len(D)} primary stations (VGZ derived, held out of the fit)")
    print(A[["netsta", "rho_rel", "dist_field_km", "noise_Z_db", "noise_H_db", "noise_3c_db",
             "contrast_Z_db", "contrast_H_db", "contrast_3c_db", "n_tremor", "n_quiet",
             "primary"]].to_string(index=False))
    out = {"n_primary": int(len(D)), "held_out": ["CN.VGZ (rho derived, not tabulated)"],
           "univariate": {}, "models": {}}
    for c in ["dist_field_km", "dist_field_min_km", "noise_Z_db", "noise_H_db", "noise_3c_db",
              "noise_allday_3c_db", "contrast_Z_db", "contrast_H_db", "contrast_3c_db"]:
        rs, p = spearmanr(D[c], D.rho_rel)
        out["univariate"][c] = dict(spearman=float(rs), p=float(p), n=int(D[c].notna().sum()))
        print(f"  univariate {c:22s} rho_s {rs:+.3f}  p={p:.3f}")
    for name, cols in MODELS.items():
        r = loo_fit(D, cols)
        out["models"][name] = {k: v for k, v in r.items() if k not in ("pred_loo", "pred_in", "idx")}
        print(f"  model {name:18s} n={r['n']:2d} in-sample rho_s {r['spearman_in']:+.3f} "
              f"R2 {r['r2_in']:+.3f} | LOO rho_s {r['spearman_loo']:+.3f} (p={r['p_loo']:.3f}) "
              f"median |err| {r['loo_err_median']:.3f}")
        if name == PRIMARY:
            D["pred_loo"] = np.nan
            D.loc[r["idx"], "pred_loo"] = r["pred_loo"]
            D.loc[r["idx"], "pred_in"] = r["pred_in"]
    pm = out["models"][PRIMARY]
    signs_ok = (pm["beta"][1] < 0 and pm["beta"][2] < 0 and pm["beta"][3] > 0)
    out["primary"] = PRIMARY
    out["signs_expected"] = bool(signs_ok)
    out["verdict"] = ("ACCEPT" if (pm["spearman_loo"] >= 0.6 and signs_ok
                                   and pm["loo_err_median"] <= 0.25)
                      else ("NOT_PREDICTIVE" if pm["spearman_loo"] < 0.4 else "MARGINAL"))
    # VGZ as an external check of the primary model
    v = A[~A.primary]
    if len(v):
        X = np.r_[1.0, v.iloc[0][MODELS[PRIMARY]].values.astype(float)]
        out["vgz_check"] = dict(predicted=float(X @ np.array(pm["beta"])),
                                derived_rho_rel=float(v.iloc[0].rho_rel))
    U.jdump(out, os.path.join(U.P4, "rho_calibration.json"))
    D.to_csv(os.path.join(U.P4, "north_calibration_primary.csv"), index=False)
    print(f"[cal] PRIMARY model {PRIMARY}: beta={np.round(pm['beta'], 4).tolist()}  "
          f"signs_expected={signs_ok}  LOO rho_s={pm['spearman_loo']:+.3f}  "
          f"median|err|={pm['loo_err_median']:.3f}  -> {out['verdict']}")
    if "vgz_check" in out:
        print(f"[cal] VGZ external check: predicted {out['vgz_check']['predicted']:.3f} vs "
              f"derived {out['vgz_check']['derived_rho_rel']:.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for axi, (c, lab) in zip(ax[:2], [("dist_field_km", "mean hypocentral distance to the 130 families (km)"),
                                      ("contrast_3c_db", "tremor/quiet 6-s energy contrast (dB, 3-C)")]):
        axi.scatter(D[c], D.rho_rel, s=45, c="#1f77b4", edgecolor="k", linewidth=0.4)
        for r in D.itertuples():
            axi.annotate(r.netsta.split(".")[1], (getattr(r, c), r.rho_rel), fontsize=7,
                         xytext=(3, 3), textcoords="offset points")
        rs, p = spearmanr(D[c], D.rho_rel)
        axi.set_xlabel(lab)
        axi.set_ylabel(r"measured $\rho/\rho_{PGC}$")
        axi.set_title(f"Spearman {rs:+.2f} (p={p:.3f}, n={len(D)})")
        axi.grid(alpha=0.25)
    ax[2].scatter(D.pred_loo, D.rho_rel, s=45, c="#d62728", edgecolor="k", linewidth=0.4)
    for r in D.itertuples():
        ax[2].annotate(r.netsta.split(".")[1], (r.pred_loo, r.rho_rel), fontsize=7,
                       xytext=(3, 3), textcoords="offset points")
    lim = [min(D.pred_loo.min(), D.rho_rel.min()) - 0.1, max(D.pred_loo.max(), D.rho_rel.max()) + 0.1]
    ax[2].plot(lim, lim, "k--", lw=1)
    ax[2].set_xlim(lim); ax[2].set_ylim(lim)
    ax[2].set_xlabel(r"leave-one-out predicted $\rho/\rho_{PGC}$")
    ax[2].set_ylabel(r"measured $\rho/\rho_{PGC}$")
    ax[2].set_title(f"LOO Spearman {pm['spearman_loo']:+.2f}, median |err| {pm['loo_err_median']:.2f}")
    ax[2].grid(alpha=0.25)
    fig.suptitle("P4-0a: calibrating the site-noise / registration proxy on the 14 north stations "
                 "with measured single-event rho", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(U.FIGS, "u_rho_calibration.png"), dpi=150)
    print(f"[fig] -> {U.FIGS}/u_rho_calibration.png")


if __name__ == "__main__":
    main()
