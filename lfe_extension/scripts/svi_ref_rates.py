#!/usr/bin/env python
"""SVI_REF R7 -- firing rate per family from the sum of P, 2014-01 .. 2026-09, with step checks.

Bostock-equivalent rate (Merlin): per family-day  rate = (sumP_same(data) - sumP_same(scrambled control)) x k(family, set)
  k(family, set) = his events / sumP_same on the truth fit days (2011 inside his window + 2012 test blocks, out-of-fold P) for
  the same station set (S15 / S14 / S13), SHRUNK toward the station-set mean in log space (empirical Bayes: between-family
  variance tau^2 by moments, sampling variance 1/his + 1/sumP), with a day-bootstrap (200 resamples of the truth days)
  for its uncertainty.  The day's set comes from the stations the production matched filter loaded that day.
  The control sum of P is subtracted AND shown.  Uncertainty = Bernoulli (sum p(1-p)) + fold-model spread + k bootstrap.
  If the inter-ETS validation failed (OUT/val/val_summary.json), quiet-day rates are given only as a control-subtracted
  bracket: [ (sumP_data - sumP_ctl) x k_lo , (N_data - N_ctl) at >= T_core of the epoch ].
Step checks at every station change (svi_ref_common.CHANGES) and every R2 flag (a station-month status change): ratios of the
  network rate in +-21 d and +-45 d windows, against the same calendar windows in the other years (z of log ratio); the window
  is marked when it contains PNSN tremor days (an SSE next to a step is not separable from the activity).
Also: the fixed-S13 count series (if the S13 production re-run exists) against the production series.
Outputs OUT/rates/: k_factors.csv, famday.parquet, monthly_family.csv, monthly_total.csv, step_checks.csv, fixedS13_compare.csv

  python scripts/svi_ref_rates.py
"""
from __future__ import annotations
import os, sys, json, glob, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

ORr = f"{R.OUT}/rates"
NB = 200


def day_sets(run):
    out = {}
    for p in glob.glob(f"{R.OUT}/feat/{run}/*.json"):
        m = json.load(open(p)); st = set(m.get("stations", []))
        out[m["date"]] = "S15" if {"LRIV", "SNB"} <= st else ("S14" if "SNB" in st else "S13")
        out[m["date"] + "_n"] = len(st)
    return out


def shrink(his, sp, sset):
    """log-space empirical-Bayes shrinkage of k = his / sumP toward the station-set mean."""
    ok = (his >= 5) & (sp >= 1)
    lk = np.log(np.where(ok, his / np.maximum(sp, 1e-9), np.nan))
    s2 = 1.0 / np.maximum(his, 1) + 1.0 / np.maximum(sp, 1e-3)
    mu = np.log(his.sum() / max(sp.sum(), 1e-9))
    v = np.nanvar(lk[ok]) if ok.sum() >= 3 else 0.0
    tau2 = max(v - np.nanmean(s2[ok]), 0.0) if ok.sum() >= 3 else 0.0
    w = np.where(ok, tau2 / (tau2 + s2), 0.0)
    return np.exp(np.where(ok, w * lk + (1 - w) * mu, mu)), mu, tau2


def kfactors(A):
    fit = A[(A.row_role == "fit") & (A.series == "data") & A.m_same.notna() & A.sset_run.isin(["S15", "S14", "S13"])]
    days = np.array(sorted(fit.date.unique()))
    B = R.bostock([y for y in (2011, 2012) if y != R.SEALED_YEAR])
    B = B[B.date.isin(days)]
    B = B[~((B.year == 2012) & (R.block_role(B.t_origin.to_numpy()) != "fit"))]
    fams = np.array(sorted(set(fit.fam.astype(int)) | set(B.fam.astype(int))))
    H = B.groupby(["fam", "date"]).size().unstack(fill_value=0).reindex(index=fams, columns=days, fill_value=0).to_numpy(float)
    SP = {ss: fit[fit.sset_run == ss].groupby(["fam", "date"]).P_same.sum().unstack(fill_value=0.0)
              .reindex(index=fams, columns=days, fill_value=0.0).to_numpy(float) for ss in ("S15", "S14", "S13")}
    def calc(w):
        his = H @ w
        return {ss: (shrink(his, SP[ss] @ w, ss), his, SP[ss] @ w) for ss in SP}
    base = calc(np.ones(len(days)))
    rng = np.random.default_rng(17); BS = {ss: [] for ss in SP}
    for _ in range(NB):
        w = np.bincount(rng.integers(0, len(days), len(days)), minlength=len(days)).astype(float)
        o = calc(w)
        for ss in o: BS[ss].append(o[ss][0][0])
    rows = []
    for ss, ((k, mu, tau2), his, sp) in base.items():
        bs = np.array(BS[ss])
        for j, f in enumerate(fams):
            rows.append(dict(fam=int(f), sset=ss, his_n=his[j], sumP=sp[j], k_raw=his[j] / sp[j] if sp[j] > 0 else np.nan,
                             k=k[j], k_lo=float(np.percentile(bs[:, j], 2.5)), k_hi=float(np.percentile(bs[:, j], 97.5)),
                             set_mean_k=float(np.exp(mu)), tau2=tau2, n_truth_days=len(days)))
    return pd.DataFrame(rows)


def main():
    warnings.filterwarnings("ignore"); R.ensure(ORr)
    cfg = json.load(open(f"{R.OUT}/conf/runs.json"))
    A = pd.read_parquet(f"{R.OUT}/conf/candidates_conf.parquet",
                        columns=["run", "kind", "sset_run", "series", "date", "year", "row_role", "fam", "t", "P_same", "P_any", "m_same", "net_mad"])
    K = kfactors(A); K.to_csv(f"{ORr}/k_factors.csv", index=False)
    Kd = {(r.fam, r.sset): (r.k, r.k_lo, r.k_hi) for r in K.itertuples()}
    P = pd.read_parquet(f"{R.OUT}/conf/P_by_candidate.parquet")
    P = P[P.run == cfg["prod"]]
    val = json.load(open(f"{R.OUT}/val/val_summary.json")) if os.path.exists(f"{R.OUT}/val/val_summary.json") else {}
    mode = val.get("inter_ets_mode", "control_subtracted_bracket")
    DS = day_sets(cfg["prod"])
    TB = C.tiers_by_epoch()
    P["Tcore"] = P.date.str[:7].map(lambda m: TB.get(C.epoch_of_month(m), (8.75, 8.25))[0])
    fsd = "P_same_WAVE_foldsd" if "P_same_WAVE_foldsd" in P else None
    g = P.groupby(["fam", "date", "series"]).agg(n=("P_same", "size"), sumP=("P_same", "sum"),
                                                 varB=("P_same", lambda p: float((p * (1 - p)).sum())),
                                                 n_core=("net_mad", lambda v: 0), quiet=("quiet", "max"))
    g = g.reset_index()
    nc = P[P.net_mad >= P.Tcore].groupby(["fam", "date", "series"]).size().rename("n_core")
    g = g.drop(columns=["n_core"]).merge(nc.reset_index(), on=["fam", "date", "series"], how="left").fillna({"n_core": 0})
    if fsd:
        sd = P.groupby(["fam", "date", "series"])[fsd].sum().rename("foldsd_sum").reset_index()
        g = g.merge(sd, on=["fam", "date", "series"], how="left")
    W = g.pivot_table(index=["fam", "date"], columns="series", values=[c for c in ["n", "sumP", "varB", "n_core", "foldsd_sum"] if c in g],
                      fill_value=0.0)
    W.columns = [f"{a}_{b}" for a, b in W.columns]
    W = W.reset_index()
    for c in ["n_data", "n_c1", "sumP_data", "sumP_c1", "varB_data", "n_core_data", "n_core_c1", "foldsd_sum_data"]:
        if c not in W: W[c] = 0.0
    W["sset"] = W.date.map(DS); W["n_sta_day"] = W.date.map(lambda d: DS.get(d + "_n", np.nan))
    W["k"] = [Kd.get((int(f), s), (np.nan,) * 3)[0] for f, s in zip(W.fam, W.sset)]
    W["k_lo"] = [Kd.get((int(f), s), (np.nan,) * 3)[1] for f, s in zip(W.fam, W.sset)]
    W["k_hi"] = [Kd.get((int(f), s), (np.nan,) * 3)[2] for f, s in zip(W.fam, W.sset)]
    W["excess"] = W.sumP_data - W.sumP_c1
    W["rate"] = W.excess * W.k
    W["rate_sd"] = np.sqrt(W.varB_data + W.foldsd_sum_data ** 2) * W.k
    q = P.groupby("date").quiet.max()
    W["quiet"] = W.date.map(q).fillna(False).astype(bool)
    W["bracket_lo"] = np.where(W.quiet, W.excess * W.k_lo, np.nan)
    W["bracket_hi"] = np.where(W.quiet, W.n_core_data - W.n_core_c1, np.nan)
    W["month"] = W.date.str[:7]
    # v1.1 counts per family-day
    v = []
    for m in R.MONTHS:
        p = f"{C.TRC}/detections_{m}.parquet"
        if os.path.exists(p):
            V = pd.read_parquet(p, columns=["fam", "date", "cat_core", "cat_ext"])
            v.append(V.groupby(["fam", "date"])[["cat_core", "cat_ext"]].sum().reset_index())
    V = pd.concat(v, ignore_index=True); V.columns = ["fam", "date", "v11_cat_core", "v11_cat_ext"]
    W = W.merge(V, on=["fam", "date"], how="left").fillna({"v11_cat_core": 0, "v11_cat_ext": 0})
    W["inter_ets_mode"] = mode
    W.to_parquet(f"{ORr}/famday.parquet", index=False)
    days_m = W.groupby("month").date.nunique()
    MF = W.groupby(["month", "fam"]).agg(sumP=("sumP_data", "sum"), sumP_ctl=("sumP_c1", "sum"), rate=("rate", "sum"),
                                          rate_var=("rate_sd", lambda s: float((s ** 2).sum())), v11_cat_core=("v11_cat_core", "sum"),
                                          v11_cat_ext=("v11_cat_ext", "sum"), quiet_days_with_candidates=("quiet", "sum")).reset_index()
    MF["days"] = MF.month.map(days_m); MF["rate_per_day"] = MF.rate / MF.days; MF["rate_per_day_sd"] = np.sqrt(MF.rate_var) / MF.days
    MF["quiet_days"] = MF.month.map(W.drop_duplicates("date").groupby("month").quiet.sum()).fillna(0).astype(int)
    MF.to_csv(f"{ORr}/monthly_family.csv", index=False)
    MT = W.groupby("month").agg(sumP=("sumP_data", "sum"), sumP_ctl=("sumP_c1", "sum"), rate=("rate", "sum"),
                                rate_var=("rate_sd", lambda s: float((s ** 2).sum())), bracket_lo=("bracket_lo", "sum"),
                                bracket_hi=("bracket_hi", "sum"), v11_cat_core=("v11_cat_core", "sum"),
                                n_quiet_famdays_with_candidates=("quiet", "sum"), sset=("sset", lambda s: ",".join(sorted(set(s.dropna()))))).reset_index()
    MT["days"] = MT.month.map(days_m); MT["rate_per_day"] = MT.rate / MT.days; MT["rate_per_day_sd"] = np.sqrt(MT.rate_var) / MT.days
    # quiet is a network-level day flag; famday holds only family-days with candidates, so count quiet days from the dates
    # (2026-10-07: the old n_quiet_famdays / quiet_days counted only family-days with a candidate; renamed *_with_candidates)
    qd = W.drop_duplicates("date").groupby("month").quiet.sum()
    MT["n_quiet_days"] = MT.month.map(qd).fillna(0).astype(int); MT["n_quiet_famdays"] = MT.n_quiet_days * W.fam.nunique()
    MT["epoch"] = MT.month.map(C.epoch_of_month); MT["inter_ets_mode"] = mode
    MT.to_csv(f"{ORr}/monthly_total.csv", index=False)
    # ---- step checks
    D = W.groupby("date").agg(rate=("rate", "sum"), excess=("excess", "sum"), v11=("v11_cat_core", "sum")).reset_index()
    D["d"] = pd.to_datetime(D.date)
    Q = pd.read_csv(C.PNSN, usecols=["time", "lat", "lon"]); Q["time"] = pd.to_datetime(Q.time, errors="coerce")
    F = pd.read_csv(R.FAMS12); lat0, lat1, lon0, lon1 = C.family_box(F)
    Q = Q[Q.lat.between(lat0, lat1) & Q.lon.between(lon0, lon1)].dropna(subset=["time"])
    tdays = set(Q.time.dt.strftime("%Y-%m-%d"))
    ch = list(R.CHANGES)
    M2 = f"{R.OUT}/r2/r2_health_monthly.csv"
    if os.path.exists(M2):
        H = pd.read_csv(M2); H = H[H.template_station].sort_values(["sta", "month"])
        H["prev"] = H.groupby("sta").status_eff.shift()
        for r in H[(H.prev.notna()) & (H.prev != H.status_eff)].itertuples():
            ch.append((f"{r.month}-01", f"R2 {r.sta} {r.prev}->{r.status_eff}"))
    rows = []
    s = D.set_index("d")
    for dd, why in ch:
        d0 = pd.Timestamp(dd)
        for w in (21, 45):
            def ratio(t0, col):
                a = s.loc[t0 - pd.Timedelta(days=w):t0 - pd.Timedelta(days=1), col].sum()
                b = s.loc[t0:t0 + pd.Timedelta(days=w - 1), col].sum()
                return b / a if a > 0 else np.nan
            for col in ("rate", "excess", "v11"):
                r0 = ratio(d0, col)
                ref = [ratio(d0.replace(year=y), col) for y in range(2014, 2027) if y != d0.year
                       and pd.Timestamp(R.PROD_T0) + pd.Timedelta(days=w) <= d0.replace(year=y) <= pd.Timestamp(R.PROD_T1) - pd.Timedelta(days=w)]
                ref = np.log(np.array([x for x in ref if np.isfinite(x) and x > 0]))
                z = (np.log(r0) - ref.mean()) / ref.std(ddof=1) if np.isfinite(r0) and r0 > 0 and len(ref) >= 4 and ref.std() > 0 else np.nan
                ntrem = sum(1 for x in pd.date_range(d0 - pd.Timedelta(days=w), d0 + pd.Timedelta(days=w)) if str(x.date()) in tdays)
                rows.append(dict(date=dd, change=why, window_d=w, series=col, ratio_after_before=r0, z_vs_other_years=z,
                                 n_ref_years=len(ref), tremor_days_in_window=ntrem, flagged=bool(np.isfinite(z) and abs(z) > 2),
                                 sse_adjacent=ntrem > 0))
    pd.DataFrame(rows).to_csv(f"{ORr}/step_checks.csv", index=False)
    # ---- fixed-S13 cross-check
    if cfg.get("prod_S13"):
        S = A[(A.run == cfg["prod_S13"]) & (A.series == "data")]
        Sc = A[(A.run == cfg["prod_S13"]) & (A.series == "c1")]
        a1 = (S.groupby(S.date.str[:7]).P_same.sum() - Sc.groupby(Sc.date.str[:7]).P_same.sum()).rename("S13_excess")
        a0 = MT.set_index("month")[["sumP", "sumP_ctl"]].eval("sumP - sumP_ctl").rename("prod_excess")
        Tq = pd.concat([a1, a0], axis=1).dropna(); Tq["ratio"] = Tq.S13_excess / Tq.prod_excess
        Tq.to_csv(f"{ORr}/fixedS13_compare.csv")
    R.log(f"[rates] mode {mode}; famdays {len(W)}; months {len(MT)}; steps flagged "
          f"{int(pd.DataFrame(rows).flagged.sum()) if rows else 0} of {len(rows)}; k-factors median by set "
          f"{K.groupby('sset').k.median().round(3).to_dict()}")


if __name__ == "__main__":
    main()
