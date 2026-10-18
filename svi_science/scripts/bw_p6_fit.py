#!/usr/bin/env python
"""bw_p6_fit -- P6 of the SVI backward extension: fit the per-set confidence models (rules: results/backward/p6/P6.md,
written before computing).

  B0 / B1 / B2  per-set refit, G1 recipe (svi_ref_fit.fit_cal: HGB + isotonic) on the G1 truth-day runs
                (backward_gates/mf/feat/g1t2011_<set>, g1t2012_<set>; 2011 role fit + 2012 TEST blocks), features WAVE_GEN
                (forward WAVE minus per-station z and has_LRIV/has_SNB).  WAVE_Z (G1 variant with z) fitted out of fold, dev only.
  era / era_noamp  the production era run p5e_2004 on Bostock's 2004 window days, BUILD blocks only; era_noamp = without rel_amp
                (applied to era days before 2003-07-31, KELB gain).
  B3 / B4       the frozen forward models (not fitted here).
Outputs p6/models/bw_models.joblib, MODELS_FROZEN_BW.json, p6/tables/dev_eval.csv, dev_reliability.csv, k_factors_bw.csv,
oof_rows.parquet.  Never reads 2003 / 2005 / 2007 / 2008 / 2009 / 2013 truth.
"""
from __future__ import annotations
import os, sys, glob, json, time, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
R, C = G.patch_forward()
import svi_ref_fit as F
import svi_ref_rates as RR
warnings.filterwarnings("ignore")

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
P6 = f"{BW}/p6"
GATES_FEAT = f"{G.OUT}/mf/feat"
STAZ_ALL = [c for c in F.WAVE if c.startswith("z_")]
WAVE_GEN = [c for c in F.WAVE if not c.startswith("z_") and c not in ("has_LRIV", "has_SNB")]
WAVE_GEN_NOAMP = [c for c in WAVE_GEN if c != "rel_amp"]
DT = pd.read_csv(f"{G.REFINE}/day_table.csv", dtype={"date": str})
ETS04 = ("2004-07-02", "2004-07-27")


def log(m): G.log(m)


def load_dir(path, dates=None):
    import pyarrow.parquet as pq
    fs = sorted(p for p in glob.glob(f"{path}/*.parquet") if not p.endswith("_fd.parquet"))
    if dates is not None:
        s = set(dates); fs = [p for p in fs if os.path.basename(p)[:10] in s]
    parts = []
    for p in fs:
        if os.path.getsize(p) == 0 or pq.read_metadata(p).num_rows == 0: continue
        x = pd.read_parquet(p); x["row"] = np.arange(len(x)); parts.append(x)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def prep(A):
    for c in F.WAVE:
        if c not in A: A[c] = np.nan
    for c in ("has_LRIV", "has_SNB", "keep_xfam8"):
        if c in A: A[c] = A[c].astype(float)
    return A


def g1_rows(s):
    L = []
    for year in (2011, 2012):
        p = f"{GATES_FEAT}/g1t{year}_{s}"
        A = load_dir(p)
        roles = DT[DT.kind == f"truth{year}"].set_index("date").role.to_dict()
        A = A[A.date.map(roles) == "fit"].copy()
        if year == 2012:
            A = A[R.block_role(A.t.to_numpy()) == "fit"]
        A["year"] = year
        B = R.bostock([year])
        d = (A.series == "data").to_numpy()
        A["m_same"] = np.nan; A["m_any"] = np.nan
        A.loc[d, "m_same"] = (R.nearest_dt(A.loc[d, "t"].to_numpy(float), A.loc[d, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                           B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
        A.loc[d, "m_any"] = (R.nearest_dt(A.loc[d, "t"].to_numpy(float), A.loc[d, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                          B.fam.to_numpy(), False) <= R.MATCH_S).astype(float)
        L.append(A)
    return prep(pd.concat(L, ignore_index=True))


def g1_truth(s):
    """his events on the fit rows' days (2012 test blocks only) -> for k factors."""
    out = []
    for year in (2011, 2012):
        days = sorted(DT[(DT.kind == f"truth{year}") & (DT.role == "fit")].date)
        B = R.bostock([year]); B = B[B.date.isin(days)]
        if year == 2012:
            B = B[R.block_role(B.t_origin.to_numpy()) == "fit"]
        out.append(B)
    return pd.concat(out, ignore_index=True)


def era_rows():
    days = [str(x.date()) for x in pd.date_range(*ETS04)]
    A = load_dir(f"{BW}/mf/feat/p5e_2004", days)
    A = A[R.block_role(A.t.to_numpy()) == "insample"].copy()
    A["year"] = 2004
    B = G.bostock_all(); B = B[B.date.isin(days)]
    Bb = B[R.block_role(B.t_origin.to_numpy()) == "insample"]
    d = (A.series == "data").to_numpy()
    A["m_same"] = np.nan; A["m_any"] = np.nan
    A.loc[d, "m_same"] = (R.nearest_dt(A.loc[d, "t"].to_numpy(float), A.loc[d, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                       B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
    A.loc[d, "m_any"] = (R.nearest_dt(A.loc[d, "t"].to_numpy(float), A.loc[d, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                      B.fam.to_numpy(), False) <= R.MATCH_S).astype(float)
    return prep(A), Bb


def oof_fit(A, feats, label):
    """5 outer folds by day (G1 seed 11): OOF P for every row (data and control), fold models, final model."""
    feats = [c for c in feats if A[c].astype(float).notna().any() and A[c].astype(float).nunique() > 1]
    X = A[feats].to_numpy(float); y = A[label].to_numpy(); dates = A.date.to_numpy()
    dl = np.array(sorted(set(dates)))
    fo = dict(zip(dl, np.random.default_rng(11).permutation(np.arange(len(dl)) % 5)))
    fold = np.array([fo[d] for d in dates])
    dm = (A.series == "data").to_numpy() & np.isfinite(y)
    P = np.full(len(A), np.nan); fms = []
    for k in range(5):
        tr = dm & (fold != k)
        m = F.fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
        P[fold == k] = m.predict_proba(X[fold == k])[:, 1]; fms.append(m)
    mfin = F.fit_cal(feats, X[dm], y[dm].astype(int), dates[dm])
    return feats, P, fms, mfin


def dev_metrics(A, P, label, tag):
    D = (A.series == "data").to_numpy() & A[label].notna().to_numpy(); Cc = (A.series == "c1").to_numpy()
    y = A.loc[D, label].to_numpy(); p = P[D]
    sl, ic, Rr = F.calib(y, p)
    r = dict(model=tag, label=label, n=int(D.sum()), n_pos=int(y.sum()), n_days=int(A.loc[D, "date"].nunique()),
             n_ctl=int(Cc.sum()), auc=F.auc(y, p), auc_net_mad=F.auc(y, A.loc[D, "net_mad"]), slope=sl,
             matched_over_sumP=float(y.sum() / p.sum()), ctl_frac_P_gt05=float(np.mean(P[Cc] > 0.5)) if Cc.any() else np.nan,
             ctl_over_data_sumP=float(np.nansum(P[Cc]) / p.sum()))
    return r, Rr.assign(model=tag, label=label)


def kf(fit, B, tag):
    """forward svi_ref_rates.kfactors logic on (data fit rows with OOF P_same, his events) for one set."""
    days = np.array(sorted(fit.date.unique()))
    B = B[B.date.isin(days)]
    fams = np.array(sorted(set(fit.fam.astype(int)) | set(B.fam.astype(int))))
    H = B.groupby(["fam", "date"]).size().unstack(fill_value=0).reindex(index=fams, columns=days, fill_value=0).to_numpy(float)
    SP = fit.groupby(["fam", "date"]).P_same.sum().unstack(fill_value=0.0).reindex(index=fams, columns=days, fill_value=0.0).to_numpy(float)
    one = np.ones(len(days))
    (k, mu, tau2) = RR.shrink(H @ one, SP @ one, tag)
    rng = np.random.default_rng(17); BS = []
    for _ in range(RR.NB):
        w = np.bincount(rng.integers(0, len(days), len(days)), minlength=len(days)).astype(float)
        BS.append(RR.shrink(H @ w, SP @ w, tag)[0])
    BS = np.array(BS); his = H @ one; sp = SP @ one
    return pd.DataFrame([dict(fam=int(f), sset=tag, his_n=his[j], sumP=sp[j], k_raw=his[j] / sp[j] if sp[j] > 0 else np.nan,
                              k=k[j], k_lo=float(np.percentile(BS[:, j], 2.5)), k_hi=float(np.percentile(BS[:, j], 97.5)),
                              set_mean_k=float(np.exp(mu)), tau2=tau2, n_truth_days=len(days)) for j, f in enumerate(fams)])


def main():
    import joblib
    t0 = time.time(); os.makedirs(f"{P6}/models", exist_ok=True); os.makedirs(f"{P6}/tables", exist_ok=True)
    mp = f"{P6}/models/bw_models.joblib"
    if os.path.exists(f"{P6}/ONELOOK_OPENED.json"):
        raise SystemExit("[p6fit] the one look is already opened: models are frozen, refusing to refit")
    MOD, FOLDS, FEATS, dev, rel, KF, OOF = {}, {}, {}, [], [], [], []
    for s in ("CN4", "B1", "B2"):
        A = g1_rows(s); log(f"[p6fit] {s}: {len(A)} rows {A.series.value_counts().to_dict()} days {A.date.nunique()}")
        for lab in ("m_same", "m_any"):
            feats, P, fms, mfin = oof_fit(A, WAVE_GEN, lab)
            key = f"{s}_{lab[2:]}"; MOD[key] = mfin; FOLDS[key] = fms; FEATS[key] = feats
            r, Rr = dev_metrics(A, P, lab, f"{s}_WAVE_GEN"); dev.append(r); rel.append(Rr)
            log(f"[p6fit] {key}: {r}")
            if lab == "m_same":
                A["P_same"] = P
                KF.append(kf(A[(A.series == "data")], g1_truth(s), s))
                OOF.append(A[["date", "series", "fam", "t", "net_mad", "m_same", "m_any", "P_same", "year"]].assign(model=s))
        feats, P, _, _ = oof_fit(A, [c for c in F.WAVE if c not in ("has_LRIV", "has_SNB")], "m_same")
        r, Rr = dev_metrics(A, P, "m_same", f"{s}_WAVE_Z(G1 variant, dev only)"); dev.append(r); rel.append(Rr)
        log(f"[p6fit] {s} WAVE_Z dev: {r}")
        if s == "CN4":                                                   # 2010 out of sample (as G1), reported only
            A0 = load_dir(f"{GATES_FEAT}/g1t2010_CN4")
            roles = DT[DT.kind == "truth2010"].set_index("date").role.to_dict()
            A0 = prep(A0[A0.date.map(roles) == "s9"].copy()); B = R.bostock([2010]); d0 = (A0.series == "data").to_numpy()
            A0["m_same"] = np.nan
            A0.loc[d0, "m_same"] = (R.nearest_dt(A0.loc[d0, "t"].to_numpy(float), A0.loc[d0, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                                 B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
            P0 = MOD["CN4_same"].predict_proba(A0[FEATS["CN4_same"]].to_numpy(float))[:, 1]
            r, Rr = dev_metrics(A0, P0, "m_same", "CN4_WAVE_GEN 2010 out of sample"); dev.append(r); rel.append(Rr)
            log(f"[p6fit] CN4 2010 oos: {r}")
    A, Bb = era_rows(); log(f"[p6fit] era: {len(A)} rows {A.series.value_counts().to_dict()} days {A.date.nunique()}; his build events {len(Bb)}")
    for tag, fl in (("ERA", WAVE_GEN), ("ERA_NOAMP", WAVE_GEN_NOAMP)):
        for lab in ("m_same", "m_any"):
            feats, P, fms, mfin = oof_fit(A, fl, lab)
            key = f"{tag}_{lab[2:]}"; MOD[key] = mfin; FOLDS[key] = fms; FEATS[key] = feats
            r, Rr = dev_metrics(A, P, lab, f"{tag}_WAVE_GEN"); dev.append(r); rel.append(Rr)
            log(f"[p6fit] {key}: {r}")
            if lab == "m_same":
                A["P_same"] = P
                KF.append(kf(A[A.series == "data"], Bb, tag))
                OOF.append(A[["date", "series", "fam", "t", "net_mad", "m_same", "m_any", "P_same", "year"]].assign(model=tag))
    joblib.dump(dict(models=MOD, folds=FOLDS, feats=FEATS), mp)
    info = dict(created_utc=pd.Timestamp.now("UTC").isoformat(), sha256=R.sha256(mp), feats=FEATS,
                forward_models=f"{G.REFINE}/conf/models.joblib", forward_sha256=R.sha256(f"{G.REFINE}/conf/models.joblib"),
                recipe="svi_ref_fit.fit_cal (HGB + isotonic CalibratedClassifierCV, GroupKFold(3) by day); 5 outer folds by day, seed 11",
                fit_rows={"CN4/B1/B2": "g1t2011_<set> role fit + g1t2012_<set> role fit, test blocks", "ERA": "p5e_2004 2004-07-02..07-27 build blocks"},
                apply={"B0": "CN4", "B1": "B1", "B2": "B2", "B3": "forward WAVE (S14)", "B4": "forward WAVE (S15)",
                       "era": "ERA (ERA_NOAMP before 2003-07-31)"})
    json.dump(info, open(f"{P6}/models/MODELS_FROZEN_BW.json", "w"), indent=1)
    pd.DataFrame(dev).to_csv(f"{P6}/tables/dev_eval.csv", index=False)
    pd.concat(rel).to_csv(f"{P6}/tables/dev_reliability.csv", index=False)
    pd.concat(KF).to_csv(f"{P6}/tables/k_factors_bw.csv", index=False)
    pd.concat(OOF).to_parquet(f"{P6}/tables/oof_rows.parquet", index=False)
    log(f"[p6fit] done {time.time() - t0:.0f}s; sha256 {info['sha256'][:16]}")
    print(pd.DataFrame(dev).round(3).to_string())


if __name__ == "__main__":
    main()
