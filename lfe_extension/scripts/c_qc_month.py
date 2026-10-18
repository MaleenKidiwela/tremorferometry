#!/usr/bin/env python
"""Track C step 4 — monthly QC (Merlin Step 7).

Per month, per SITE:
  quiet_rate   declustered detections per template-day on QUIET family-days (no PNSN tremor in the
               family box that day and >= 3 days from any tremor day) counted over the family-days
               that SITE contributed to, at T_ext and at T_core;
  lambda       that site's RUNNING FLOOR: the median quiet_rate over every earlier month already in
               results/trackC/qc_summary.csv (>= 2 earlier months required, else NaN);
  flag_2lambda quiet_rate > 2 x lambda  -> the site is misbehaving (Merlin's rule).
Per month, whole network:
  unobservable fraction of family-days (< 2 sites with a usable template AND data),
  EQ-flag fraction of detections (within 60 s of a catalogued earthquake),
  HF-ratio distribution (median / p90 / fraction above the frozen cut) and its shift vs the running
  median — the response-change alarm the frozen HF cut needs,
  detections per day per tier, station-day coverage and the daily noise level per station.

Writes results/trackC/qc_<YYYY-MM>.csv (per-site rows + one network row), appends to
results/trackC/qc_summary.csv, and regenerates figures/trackC_monthly_qc.png.

  python scripts/c_qc_month.py 2015-12
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

SUM = None


def month_qc(ym, out_dir=None):
    OUT = out_dir or C.TRC
    fz = C.frozen()
    T_ext = float(fz["T_ext"]); T_core = float(fz["T_core"])
    hf_cut = float(fz.get("hf_veto", {}).get("cut", np.nan))
    fdp = f"{OUT}/family_days_{ym}.parquet"
    dtp = f"{OUT}/detections_{ym}.parquet"
    if not os.path.exists(fdp):
        raise SystemExit(f"[qc] {fdp} missing — run scripts/c_run_month.py {ym} first")
    FD = pd.read_parquet(fdp)
    DET = pd.read_parquet(dtp) if os.path.exists(dtp) else pd.DataFrame()
    sump = f"{OUT}/qc_summary.csv"
    prev = pd.read_csv(sump) if os.path.exists(sump) else pd.DataFrame()
    if len(prev):
        prev = prev[prev.ym != ym]

    obs = FD[FD.observable]
    quiet = obs[obs.quiet_flag]
    rows = []
    sites = sorted({q for s in obs.sites.dropna().unique() for q in str(s).split(",")})
    for q in sites:
        sel = quiet[quiet.sites.fillna("").str.split(",").apply(lambda L: q in L)]
        nfd = len(sel)
        r_ext = float(sel.n_det_ext.sum() / nfd) if nfd else np.nan
        r_core = float(sel.n_det_core.sum() / nfd) if nfd else np.nan
        lam_e = lam_c = np.nan
        if len(prev):
            p = prev[(prev.site == q) & (prev.kind == "site") & prev.quiet_rate_ext.notna()]
            if len(p) >= 2:
                lam_e = float(p.quiet_rate_ext.median()); lam_c = float(p.quiet_rate_core.median())
        rows.append(dict(ym=ym, kind="site", site=q, n_quiet_family_days=nfd,
                         quiet_rate_ext=r_ext, quiet_rate_core=r_core,
                         lambda_ext=lam_e, lambda_core=lam_c,
                         flag_2lambda=bool(np.isfinite(lam_e) and r_ext > 2 * lam_e)))
    # ---- network row
    nfd_all = len(FD)
    unobs = float((~FD.observable).mean()) if nfd_all else np.nan
    nd = len(DET)
    net = dict(ym=ym, kind="network", site="ALL", n_family_days=nfd_all,
               frac_unobservable=unobs,
               n_quiet_family_days=int(len(quiet)),
               quiet_rate_ext=float(quiet.n_det_ext.sum() / len(quiet)) if len(quiet) else np.nan,
               quiet_rate_core=float(quiet.n_det_core.sum() / len(quiet)) if len(quiet) else np.nan,
               n_det_rec=nd)
    if nd:
        net.update(n_det_ext=int(DET.tier_ext.sum()), n_det_core=int(DET.tier_core.sum()),
                   n_det_core_hfpass=int((DET.tier_core & DET.hf_pass).sum()),
                   n_cat_ext=int(DET.cat_ext.sum()) if "cat_ext" in DET else np.nan,
                   n_cat_core=int(DET.cat_core.sum()) if "cat_core" in DET else np.nan,
                   frac_quiet_det_core=float(DET[DET.tier_core].quiet_flag.mean())
                   if "quiet_flag" in DET and DET.tier_core.any() else np.nan,
                   epoch_template=("|".join(f"{k}:{v}" for k, v in
                                            DET.epoch_template.value_counts().items())
                                   if "epoch_template" in DET else ""),
                   hf_fs=float(DET.hf_fs.median()) if "hf_fs" in DET and DET.hf_fs.notna().any() else np.nan,
                   frac_eq=float(DET.eq_flag.mean()),
                   frac_eq_core=float(DET[DET.tier_core].eq_flag.mean()) if DET.tier_core.any() else np.nan,
                   frac_tremor=float(DET.tremor_flag.mean()),
                   hfr_median=float(np.nanmedian(DET.hfr)), hfr_p90=float(np.nanpercentile(DET.hfr.dropna(), 90))
                   if DET.hfr.notna().any() else np.nan,
                   hfr_frac_above_cut=float((DET.hfr > hf_cut).mean()),
                   hfr_frac_defined=float(DET.hfr.notna().mean()),
                   median_n_sites=float(DET.n_sites.median()),
                   median_wsh_pgc=float(DET.wsh_pgc.median()))
        if len(prev):
            p = prev[(prev.kind == "network") & prev.hfr_median.notna()]
            if len(p) >= 2:
                net["hfr_median_running"] = float(p.hfr_median.median())
                net["hfr_shift"] = net["hfr_median"] / max(net["hfr_median_running"], 1e-12)
                net["flag_hfr_shift"] = bool(net["hfr_shift"] > 2 or net["hfr_shift"] < 0.5)
    # ---- station-day coverage / noise from the cache pass
    npq = f"{OUT}/daily_noise_{ym}.parquet"
    if os.path.exists(npq):
        N = pd.read_parquet(npq)
        net["n_station_days_cached"] = int((N.status != "error").sum())
        net["n_station_days_error"] = int((N.status == "error").sum())
        net["median_frac_nonzero"] = float(N.frac_nonzero.median())
        for r in N.groupby("netsta").agg(e8=("e8_Z", "median"), cov=("frac_nonzero", "median")).reset_index().itertuples():
            rows.append(dict(ym=ym, kind="station", site=r.netsta, noise_e8_Z_median=float(r.e8),
                             frac_nonzero_median=float(r.cov)))
    rows.append(net)
    Q = pd.DataFrame(rows)
    Q.to_csv(f"{OUT}/qc_{ym}.csv", index=False)
    allq = pd.concat([prev, Q], ignore_index=True) if len(prev) else Q
    allq = allq.sort_values(["ym", "kind", "site"])
    allq.to_csv(sump, index=False)
    return Q, allq


def figure(allq, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    net = allq[allq.kind == "network"].copy()
    site = allq[allq.kind == "site"].copy()
    if not len(net):
        return
    net["x"] = pd.to_datetime(net.ym + "-01")
    fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    ax[0].plot(net.x, net.get("n_det_ext", pd.Series(dtype=float)), "o-", label="T_ext")
    ax[0].plot(net.x, net.get("n_det_core", pd.Series(dtype=float)), "s-", label="T_core")
    ax[0].plot(net.x, net.get("n_det_core_hfpass", pd.Series(dtype=float)), "^-", label="T_core + HF veto")
    ax[0].set_yscale("symlog"); ax[0].set_ylabel("detections / month"); ax[0].legend(fontsize=8)
    ax[0].set_title("Track C monthly QC")
    for q, g in site.groupby("site"):
        g = g.copy(); g["x"] = pd.to_datetime(g.ym + "-01")
        ax[1].plot(g.x, g.quiet_rate_ext, "o-", ms=3, label=q)
    ax[1].set_ylabel("quiet-day det / template-day (T_ext)")
    ax[1].legend(fontsize=6, ncol=4)
    ax[2].plot(net.x, net.frac_unobservable, "o-", color="C3", label="unobservable family-days")
    ax[2].plot(net.x, net.get("frac_eq", pd.Series(dtype=float)), "s-", color="C1", label="EQ-flagged detections")
    ax[2].set_ylabel("fraction"); ax[2].legend(fontsize=8)
    ax[3].plot(net.x, net.get("hfr_median", pd.Series(dtype=float)), "o-", color="C4", label="median HF ratio")
    ax[3].plot(net.x, net.get("hfr_frac_above_cut", pd.Series(dtype=float)), "s-", color="C5",
               label="fraction above the frozen cut")
    ax[3].set_yscale("log"); ax[3].set_ylabel("HF ratio"); ax[3].legend(fontsize=8)
    ax[3].set_xlabel("month")
    for b in ax:
        b.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ym")
    ap.add_argument("--out-dir", default="")
    a = ap.parse_args()
    Q, allq = month_qc(a.ym, a.out_dir or None)
    pd.set_option("display.width", 250)
    print(Q[Q.kind == "site"][["site", "n_quiet_family_days", "quiet_rate_ext", "quiet_rate_core",
                               "lambda_ext", "flag_2lambda"]].round(4).to_string(index=False))
    net = Q[Q.kind == "network"].iloc[0]
    print(json.dumps({k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                      for k, v in net.dropna().to_dict().items() if k not in ("kind", "site")},
                     indent=1, default=str))
    os.makedirs(f"{C.EXT}/figures", exist_ok=True)
    figure(allq, f"{C.EXT}/figures/trackC_monthly_qc.png")
    print(f"[qc] -> {(a.out_dir or C.TRC)}/qc_{a.ym}.csv, qc_summary.csv, figures/trackC_monthly_qc.png")


if __name__ == "__main__":
    main()
