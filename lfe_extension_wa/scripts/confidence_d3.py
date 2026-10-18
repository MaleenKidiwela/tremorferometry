#!/usr/bin/env python
"""CONFIDENCE add-on D3 (2026-09-30): the known-answer model of D2 (fit on the TEST half, 5 outer folds by day;
A2011 from a model fitted on the whole test half) but the SCRAMBLED-CONTROL candidates enter the fit as extra known
negatives, so the model cannot score a candidate high only because of its family / busy-hour context.

  rows    data candidates with truth labels (m_own / m_nbr3 / m_any) + control candidates c1, c2, c3 of the same
          family-days (label 0; together weighted like the data rows).  c4 never enters any fit (never-seen check).
  model   HistGradientBoosting (same settings as confidence_fit.hgb) with sample weights, 3 inner folds by day;
          isotonic calibration fitted on the inner out-of-fold scores of DATA rows only (so P stays calibrated for data).
  feats   NO_STEP2 (step-1 features + P1 + PNET off-the-shelf + fine-tuned); step-2 attribution exists for c1 only.
  -> results/confidence/d3_*.csv, candidates_confidence_d3.parquet
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import confidence_fit as CF

FEATS = CF.FSETS["NO_STEP2"]
TRAIN_CTL = ["c1", "c2", "c3"]


class D3Model:
    def __init__(self, models, iso):
        self.models, self.iso = models, iso

    def predict(self, X):
        s = np.mean([m.predict_proba(X)[:, 1] for m in self.models], 0)
        return self.iso.predict(s)


def fit_d3(X, y, isdata, dates, seed=0):
    from sklearn.model_selection import GroupKFold
    from sklearn.isotonic import IsotonicRegression
    nd = isdata.sum(); nc = (~isdata).sum()
    w = np.where(isdata, 1.0, nd / max(nc, 1))
    oof = np.full(len(y), np.nan)
    models = []
    for tr, te in GroupKFold(3).split(X, y, dates):
        m = CF.hgb(FEATS, seed)
        m.fit(X[tr], y[tr], sample_weight=w[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
        models.append(m)
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(oof[isdata], y[isdata])
    return D3Model(models, iso)


def main():
    t0 = time.time()
    A = CF.assemble()
    A = CF.crossfit_P1(A)
    X = A[FEATS].to_numpy(float)
    dates = A.date.to_numpy()
    tst = (A.window.isin(["V2008", "V2010"]) & (A.half == 1)).to_numpy()
    isd = (A.series == "data").to_numpy()
    lab_ok = (isd & A.truth_ok.to_numpy())
    ctl_ok = A.series.isin(TRAIN_CTL).to_numpy() & A.truth_ok.to_numpy()   # truth_ok = the family's array run scanned the day
    tdays = np.array(sorted(set(dates[tst])))
    rng = np.random.default_rng(11)                    # same outer folds as confidence_fit D2
    fo = dict(zip(tdays, rng.permutation(np.arange(len(tdays)) % 5)))
    fold = np.array([fo.get(d, -1) for d in dates])
    a11 = (A.window == "A2011").to_numpy()
    PR = {}
    for lab, lc in CF.LABELS.items():
        y = A[lc].to_numpy().astype(int)
        y = np.where(isd, y, 0)
        pred = np.full(len(A), np.nan)
        for k in range(5):
            tr = (lab_ok | ctl_ok) & tst & (fold != k)
            m = fit_d3(X[tr], y[tr], isd[tr], dates[tr])
            sc = tst & (fold == k)
            pred[sc] = m.predict(X[sc])
        tr = (lab_ok | ctl_ok) & tst
        m = fit_d3(X[tr], y[tr], isd[tr], dates[tr])
        pred[a11] = m.predict(X[a11])
        PR[f"D3_hgb_NO_STEP2_{lab}"] = pred
        CF.log(f"[d3] {lab} {time.time()-t0:.0f}s")
    O = pd.read_parquet(f"{CF.OD}/candidates_confidence.parquet")
    assert (O.fid.values == A.fid.values).all() and (O.t.values == A.t.values).all()
    for k, v in PR.items():
        O[k] = v
    O.to_parquet(f"{CF.OD}/candidates_confidence.parquet", index=False)
    # ---------------- metrics: D3 next to D2 comparators
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score
    models = ["D3_hgb_NO_STEP2", "D2_hgb_NO_STEP2", "D2_hgb_FULL"]
    D = O[(O.series == "data") & O.truth_ok]
    rows, crow = [], []
    for es in CF.EVAL:
        ge = D[CF.eval_set_mask(D, es)]
        for rg in ["ALL"] + list(CF.RS.values()):
            gr = ge if rg == "ALL" else ge[ge.reg == rg]
            if len(gr) < 50:
                continue
            for lab, lc in CF.LABELS.items():
                y = gr[lc].to_numpy().astype(int)
                r = dict(eval_set=es, region=rg, label=lab, n=len(gr), n_pos=int(y.sum()))
                for mo in models:
                    p = gr[f"{mo}_{lab}"].to_numpy()
                    r[f"auc_{mo}"], r[f"ap_{mo}"] = CF.auc_ap(y, p)
                    sl, ic, _ = CF.calib(y, p)
                    r[f"slope_{mo}"] = sl; r[f"obs_exp_{mo}"] = y.sum() / p.sum() if p.sum() > 0 else np.nan
                    # within family-day AUC (pairs weighted)
                    num = den = 0.0
                    for _, h in gr.groupby(["fid", "date"]):
                        yy = h[lc].to_numpy().astype(int)
                        if 0 < yy.sum() < len(yy):
                            wv = yy.sum() * (len(yy) - yy.sum()); num += wv * roc_auc_score(yy, h[f"{mo}_{lab}"].to_numpy()); den += wv
                    r[f"auc_withinFD_{mo}"] = num / den if den else np.nan
                for c in ["P1", "pnet_max_f"]:
                    num = den = 0.0
                    for _, h in gr.groupby(["fid", "date"]):
                        yy = h[lc].to_numpy().astype(int)
                        if 0 < yy.sum() < len(yy):
                            wv = yy.sum() * (len(yy) - yy.sum()); num += wv * roc_auc_score(yy, h[c].to_numpy()); den += wv
                    r[f"auc_withinFD_{c}"] = num / den if den else np.nan
                rows.append(r)
    pd.DataFrame(rows).to_csv(f"{CF.OD}/d3_auc_calib.csv", index=False)
    rows = []
    for es in CF.EVAL:
        ge = O[CF.eval_set_mask(O, es)]
        for rg in ["ALL"] + list(CF.RS.values()):
            gr = ge if rg == "ALL" else ge[ge.reg == rg]
            for days_ in ["all", "quiet", "episode"]:
                gk = gr if days_ == "all" else (gr[gr.quiet] if days_ == "quiet" else gr[~gr.quiet])
                d = gk[gk.series == "data"]
                if not len(d):
                    continue
                for mo in models:
                    for lab in ["own", "any"]:
                        key = f"{mo}_{lab}"
                        for cset in (["c1", "c2", "c3"], ["c4"]):
                            c = gk[gk.series.isin(cset)]
                            rows.append(dict(eval_set=es, region=rg, days=days_, model=key, controls=",".join(cset),
                                             n_data=len(d), n_ctl_per_real=len(c) / len(cset),
                                             data_fracP_gt05=(d[key] > 0.5).mean(), ctl_fracP_gt05=(c[key] > 0.5).mean(),
                                             data_sumP=d[key].sum(), ctl_sumP_per_real=c[key].sum() / len(cset),
                                             ctl_over_data_sumP=(c[key].sum() / len(cset)) / d[key].sum() if d[key].sum() > 0 else np.nan))
    pd.DataFrame(rows).to_csv(f"{CF.OD}/d3_control_quiet.csv", index=False)
    RF = pd.read_parquet(f"{W.RES}/step1_prob/rate_family_day.parquet")
    RF = RF[RF.truth_ok & (RF.kind == "episode")]
    rows = []
    Dd = O[O.series == "data"]
    for es in CF.EVAL:
        ge = Dd[CF.eval_set_mask(Dd, es)]
        w = "A2011" if es == "A2011" else es.split("_")[0]
        R = RF[RF.window == w][["fid", "date", "region", "truth_own_scored", "truth_nbr3_scored"]]
        R = R.assign(reg=R.region.map(CF.RS))
        for key, tcol in [("D3_hgb_NO_STEP2_own", "truth_own_scored"), ("D3_hgb_NO_STEP2_nbr3", "truth_nbr3_scored")]:
            s = ge.groupby(["fid", "date"])[key].sum().rename("sP").reset_index()
            RR = R.merge(s, on=["fid", "date"], how="left").fillna({"sP": 0})
            for rg in ["ALL"] + list(CF.RS.values()):
                gr = RR if rg == "ALL" else RR[RR.reg == rg]
                if len(gr) < 5 or gr[tcol].sum() == 0:
                    continue
                rows.append(dict(eval_set=es, region=rg, model=key, truth=tcol, n_famdays=len(gr), sumP=gr.sP.sum(),
                                 truth_count=int(gr[tcol].sum()), ratio=gr.sP.sum() / gr[tcol].sum(),
                                 spearman=spearmanr(gr.sP, gr[tcol]).statistic))
    pd.DataFrame(rows).to_csv(f"{CF.OD}/d3_expected_count.csv", index=False)
    CF.log(f"[d3] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
