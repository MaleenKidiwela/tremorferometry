#!/usr/bin/env python
"""Track E step 3 — train p_real on 2012, VALIDATE on the held-out year 2010.

TARGET (stated, and the only thing that is actually observable):
    y = 1  the detection is matched one-to-one (+-2 s) to a Bostock event of its OWN family
    y = 0  otherwise
Bostock's catalog exists only inside the 2010 / 2012 ETS episodes, so the model is TRAINED and
CALIBRATED on the ETS-day population, where both classes live on the same days and the model
therefore cannot cheat on a day-level property.  Quiet-day (inter-ETS) detections are added as
NEGATIVES (Track A step 3b: >= 78 % of them are noise at T_core, >= 94 % at T_ext), and the
SITE-ROLLED chance control is used as an independent negative reference, not as a training class
(the control has no rel_amp / hfr / wsh_pgc, so training on it would let the model separate the
classes on a missing-value pattern instead of on physics).

  p_real  = calibrated P(y = 1 | features)  on the ETS-day prior.  VALIDATED on 2010.
  p_real_prior = the same likelihood ratio re-based on an independently measured realness prior for
  the population the event belongs to (frozen point precision_any_ets inside tremor windows, Track A
  step 3b real fraction outside them).  Reported, never validated — Bostock cannot validate it.

Labels are never features.  Context (tremor_flag, quiet_flag, season, epoch) is never a feature.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GroupKFold
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
RES = f"{EXT}/results"
OUT = f"{RES}/trackE"
FIG = f"{EXT}/figures"

FEATS = ["net_mad", "net_cc", "rel_amp", "hfr_rel", "n_sites_obs", "wsh_pgc",
         "n_sites_above", "frac_sites_above", "z_site_med", "z_site_min", "z_ex_top",
         "conf_z", "conf_ze", "ze_site_med",
         "n_same10", "n_same60", "n_other5", "burst_z"]
MONO = {"net_mad": 1, "n_sites_above": 1, "frac_sites_above": 1, "conf_z": 1, "conf_ze": 1,
        "burst_z": 1, "hfr_rel": -1}
# feature groups for the ablation
GROUPS = {
    "physics":  ["net_mad", "net_cc", "rel_amp", "hfr_rel", "n_sites_obs", "wsh_pgc"],
    "sites":    ["n_sites_above", "frac_sites_above", "z_site_med", "z_site_min", "z_ex_top"],
    "confirm":  ["conf_z", "conf_ze", "ze_site_med"],
    "burst":    ["n_same10", "n_same60", "n_other5", "burst_z"],
}
# control-safe subset (everything the site-rolled control also has)
SAFE = [f for f in FEATS if f not in ("rel_amp", "hfr_rel")]


def mk_model(kind, feats):
    if kind == "hgb":
        return HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.06, max_leaf_nodes=15, min_samples_leaf=60,
            l2_regularization=1.0, early_stopping=False, random_state=0,
            monotonic_cst=[MONO.get(f, 0) for f in feats])
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(max_iter=2000, C=1.0))


def cv_scores(kind, feats, X, y, g, nfold=5):
    oof = np.full(len(y), np.nan)
    for tr, te in GroupKFold(n_splits=nfold).split(X, y, groups=g):
        m = mk_model(kind, feats)
        m.fit(X[feats].to_numpy()[tr], y[tr])
        oof[te] = m.predict_proba(X[feats].to_numpy()[te])[:, 1]
    return oof


def reliability(p, y, edges):
    rows = []
    for i in range(len(edges) - 1):
        last = (i == len(edges) - 2)
        m = (p >= edges[i]) & ((p <= edges[i + 1]) if last else (p < edges[i + 1]))
        if m.sum() == 0:
            rows.append(dict(lo=edges[i], hi=edges[i + 1], n=0, p_pred=np.nan, p_obs=np.nan, dev=np.nan))
            continue
        rows.append(dict(lo=edges[i], hi=edges[i + 1], n=int(m.sum()), p_pred=float(p[m].mean()),
                         p_obs=float(y[m].mean()), dev=float(y[m].mean() - p[m].mean())))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--f2012", default=f"{OUT}/feat_ref2012.parquet")
    ap.add_argument("--f2010", default=f"{OUT}/feat_ref2010.parquet")
    ap.add_argument("--min-bin-n", type=int, default=30)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    rep = {}

    A = pd.read_parquet(a.f2012)
    B = pd.read_parquet(a.f2010)
    for D in (A, B):
        D["match_own"] = D.match_own.fillna(False).astype(bool)
    # the ETS window of each year = the doy span over which Bostock's catalog has events
    def ets_span(D):
        d = D[(D.ctrl == 0) & D.match_own].doy
        return int(d.min()), int(d.max())
    e12, e10 = ets_span(A), ets_span(B)
    rep["ets_span_2012"], rep["ets_span_2010"] = e12, e10
    A12 = A[(A.ctrl == 0) & A.doy.between(*e12)].reset_index(drop=True)          # training ETS pop
    Q12 = A[(A.ctrl == 0) & (A.doy < e12[0])].reset_index(drop=True)             # quiet-day negatives
    C12 = A[A.ctrl == 1].reset_index(drop=True)                                  # rolled control
    B10 = B[(B.ctrl == 0) & B.doy.between(*e10)].reset_index(drop=True)          # VALIDATION
    C10 = B[B.ctrl == 1].reset_index(drop=True)
    P12 = A[(A.ctrl == 0) & (A.doy > e12[1])].reset_index(drop=True)             # 2012 post-SSE

    TR = pd.concat([A12.assign(src="ets"), Q12.assign(src="quiet")], ignore_index=True)
    TR["y"] = TR.match_own.astype(int)
    rep["labels"] = dict(
        train_ets_n=int(len(A12)), train_ets_pos=int(A12.match_own.sum()),
        train_ets_prior=float(A12.match_own.mean()),
        train_quiet_neg_n=int(len(Q12)), train_total=int(len(TR)),
        train_pos_frac=float(TR.y.mean()),
        control_2012_n=int(len(C12)), control_2010_n=int(len(C10)),
        valid_2010_n=int(len(B10)), valid_2010_pos=int(B10.match_own.sum()),
        valid_2010_prior=float(B10.match_own.mean()),
        post_sse_2012_n=int(len(P12)))
    print("[labels]", json.dumps(rep["labels"], indent=1), flush=True)

    y = TR.y.to_numpy(); g = TR.date.to_numpy()
    sel = {}
    for kind in ("hgb", "logreg"):
        oof = cv_scores(kind, FEATS, TR, y, g)
        m_ets = (TR.src == "ets").to_numpy()
        sel[kind] = dict(auc_all=float(roc_auc_score(y, oof)),
                         auc_ets=float(roc_auc_score(y[m_ets], oof[m_ets])))
        print(f"[cv-2012] {kind}: {sel[kind]}", flush=True)
    rep["cv_2012"] = sel
    kind = max(sel, key=lambda k: sel[k]["auc_ets"])
    rep["model_selected"] = kind
    print(f"[cv-2012] selected {kind} (by ETS-day grouped-CV AUC on 2012)", flush=True)

    # final model + isotonic calibration on the OOF scores of the ETS-day population
    oof = cv_scores(kind, FEATS, TR, y, g)
    m_ets = (TR.src == "ets").to_numpy()
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(oof[m_ets], y[m_ets])
    model = mk_model(kind, FEATS)
    model.fit(TR[FEATS].to_numpy(), y)

    def score(D):
        raw = model.predict_proba(D[FEATS].to_numpy())[:, 1]
        return raw, np.clip(iso.predict(raw), 1e-4, 1 - 1e-4)

    # ---------------- VALIDATION on 2010 ----------------
    raw10, p10 = score(B10)
    y10 = B10.match_own.to_numpy().astype(int)
    edges = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    R10 = reliability(p10, y10, edges)
    R10.to_csv(f"{OUT}/reliability_2010.csv", index=False)
    big = R10[R10.n >= a.min_bin_n]
    worst = float(np.nanmax(np.abs(big.dev))) if len(big) else np.nan
    rep["valid_2010"] = dict(
        auc=float(roc_auc_score(y10, p10)), n=int(len(B10)), prior=float(y10.mean()),
        max_abs_dev_bins_ge_n=worst, min_bin_n=a.min_bin_n,
        calibration_pass=bool(np.isfinite(worst) and worst <= 0.10),
        bins=R10.to_dict("records"))
    for thr in (0.5, 0.8, 0.9):
        m = p10 >= thr
        rep["valid_2010"][f"prec_ge_{thr}"] = float(y10[m].mean()) if m.sum() else np.nan
        rep["valid_2010"][f"n_ge_{thr}"] = int(m.sum())
        rep["valid_2010"][f"recall_ge_{thr}"] = float(y10[m].sum() / max(y10.sum(), 1))
    # in-sample 2012 (OOF) for reference
    p12 = np.clip(iso.predict(oof[m_ets]), 1e-4, 1 - 1e-4)
    R12 = reliability(p12, y[m_ets], edges)
    R12.to_csv(f"{OUT}/reliability_2012_oof.csv", index=False)
    b12 = R12[R12.n >= a.min_bin_n]
    rep["oof_2012"] = dict(auc=float(roc_auc_score(y[m_ets], p12)),
                           max_abs_dev=float(np.nanmax(np.abs(b12.dev))) if len(b12) else np.nan)
    print("[valid-2010]", json.dumps({k: v for k, v in rep["valid_2010"].items() if k != "bins"},
                                     indent=1), flush=True)

    # ---- is the 2010 miss a LEVEL shift or a SHAPE failure?  One free scalar (the population
    # log-odds anchor) is fitted on 2010 and the reliability re-measured: if the residual is inside
    # the bar, the score's ORDERING transfers and only the absolute level is network-dependent.
    from scipy.optimize import minimize_scalar
    lo10 = np.log(p10 / (1 - p10))
    def _nll(b):
        q = np.clip(1 / (1 + np.exp(-(lo10 + b))), 1e-9, 1 - 1e-9)
        return -(y10 * np.log(q) + (1 - y10) * np.log(1 - q)).sum()
    shift = float(minimize_scalar(_nll, bounds=(-4, 4), method="bounded").x)
    p10a = 1 / (1 + np.exp(-(lo10 + shift)))
    R10a = reliability(p10a, y10, edges)
    R10a.to_csv(f"{OUT}/reliability_2010_anchored.csv", index=False)
    b10a = R10a[R10a.n >= a.min_bin_n]
    prior_shift = float(np.log(y10.mean() / (1 - y10.mean())) -
                        np.log(A12.match_own.mean() / (1 - A12.match_own.mean())))
    rep["valid_2010_anchored"] = dict(
        logodds_shift=shift, prior_only_shift=prior_shift,
        max_abs_dev=float(np.nanmax(np.abs(b10a.dev))) if len(b10a) else np.nan,
        passes_bar=bool(len(b10a) and np.nanmax(np.abs(b10a.dev)) <= 0.10),
        bins=R10a.to_dict("records"))
    print("[valid-2010-anchored]", json.dumps({k: v for k, v in
          rep["valid_2010_anchored"].items() if k != "bins"}, indent=1), flush=True)

    # ---------------- ablation: what carries the score ----------------
    abl = {}
    for name, fs in list(GROUPS.items()) + [("all", FEATS), ("all_minus_burst",
                                            [f for f in FEATS if f not in GROUPS["burst"]]),
                                            ("all_minus_confirm",
                                            [f for f in FEATS if f not in GROUPS["confirm"]]),
                                            ("all_minus_sites",
                                            [f for f in FEATS if f not in GROUPS["sites"]]),
                                            ("safe_for_control", SAFE)]:
        o = cv_scores(kind, fs, TR, y, g)
        mo = mk_model(kind, fs); mo.fit(TR[fs].to_numpy(), y)
        pv = mo.predict_proba(B10[fs].to_numpy())[:, 1]
        abl[name] = dict(n_feat=len(fs), cv2012_auc_ets=float(roc_auc_score(y[m_ets], o[m_ets])),
                         valid2010_auc=float(roc_auc_score(y10, pv)))
        print(f"[abl] {name:18s} {abl[name]}", flush=True)
    rep["ablation"] = abl

    # ---------------- the site-rolled chance control ----------------
    msafe = mk_model(kind, SAFE); msafe.fit(TR[SAFE].to_numpy(), y)
    osafe = cv_scores(kind, SAFE, TR, y, g)
    isos = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(osafe[m_ets], y[m_ets])
    ctl = {}
    for nm, Dc, Dd in (("2012", C12, A12), ("2010", C10, B10)):
        if not len(Dc):
            continue
        pc = np.clip(isos.predict(msafe.predict_proba(Dc[SAFE].to_numpy())[:, 1]), 1e-4, 1 - 1e-4)
        pdat = np.clip(isos.predict(msafe.predict_proba(Dd[SAFE].to_numpy())[:, 1]), 1e-4, 1 - 1e-4)
        ctl[nm] = dict(n_ctrl=int(len(Dc)), ctrl_median_p=float(np.median(pc)),
                       ctrl_frac_ge_05=float((pc >= 0.5).mean()),
                       ctrl_frac_ge_08=float((pc >= 0.8).mean()),
                       data_median_p=float(np.median(pdat)),
                       data_frac_ge_08=float((pdat >= 0.8).mean()),
                       pos_frac_ge_08=float((pdat[Dd.match_own.to_numpy()] >= 0.8).mean()))
        print(f"[control-{nm}]", ctl[nm], flush=True)
    rep["control"] = ctl

    # ---------------- quiet-day / post-SSE behaviour ----------------
    for nm, Dq in (("quiet2012", Q12), ("postsse2012", P12)):
        if not len(Dq):
            continue
        _, pq = score(Dq)
        rep.setdefault("offets_2012", {})[nm] = dict(
            n=int(len(Dq)), median_p=float(np.median(pq)), frac_ge_05=float((pq >= 0.5).mean()),
            frac_ge_08=float((pq >= 0.8).mean()))
        print(f"[{nm}]", rep["offets_2012"][nm], flush=True)

    # ---------------- figures ----------------
    BLUE, ORANGE = "#2a78d6", "#eb6834"
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
    ax[0].plot([0, 1], [0, 1], color="#52514e", ls="--", lw=1, zorder=2)
    for R, lab, c, ls in ((R12, "2012 (out-of-fold, calibrated here)", BLUE, "-"),
                          (R10, "2010 (held-out, as delivered)", ORANGE, "-"),
                          (R10a, "2010 after ONE log-odds anchor", ORANGE, ":")):
        m = R.n >= a.min_bin_n
        ax[0].plot(R.p_pred[m], R.p_obs[m], "o" + ls, color=c, label=lab, ms=4, lw=1.8,
                   mfc="white" if ls == ":" else c, zorder=3)
        if ls == "-":
            for _, r in R[m].iterrows():
                ax[0].annotate(int(r.n), (r.p_pred, r.p_obs), fontsize=6, xytext=(2, -8),
                               textcoords="offset points", color=c)
    ax[0].fill_between([0, 1], [-0.1, 0.9], [0.1, 1.1], color="#f0efec", zorder=0,
                       label="+-0.1 calibration bar")
    ax[0].set_xlabel("predicted p_real"); ax[0].set_ylabel("observed Bostock own-match fraction")
    ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1); ax[0].legend(fontsize=8)
    ax[0].set_title(f"reliability — {kind}, 2010 AUC {rep['valid_2010']['auc']:.3f}; "
                    f"anchored max dev {rep['valid_2010_anchored']['max_abs_dev']:.3f}")
    for s_ in ("top", "right"):
        ax[0].spines[s_].set_visible(False); ax[1].spines[s_].set_visible(False)
    ax[1].hist([p10[y10 == 1], p10[y10 == 0]], bins=20, stacked=False,
               color=[BLUE, ORANGE],
               label=["2010 Bostock-matched", "2010 unmatched"], density=True)
    ax[1].set_xlabel("p_real"); ax[1].set_ylabel("density"); ax[1].legend(fontsize=8)
    ax[1].set_title("2010 score separation")
    fig.tight_layout(); fig.savefig(f"{FIG}/e_reliability_2010.png", dpi=140); plt.close(fig)

    # permutation importance on the held-out year
    from sklearn.inspection import permutation_importance
    pi = permutation_importance(model, B10[FEATS].to_numpy(), y10, n_repeats=8,
                                random_state=0, scoring="roc_auc", n_jobs=8)
    imp = pd.DataFrame(dict(feature=FEATS, imp=pi.importances_mean, sd=pi.importances_std)
                       ).sort_values("imp", ascending=False)
    imp.to_csv(f"{OUT}/feature_importance_2010.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.barh(imp.feature[::-1], imp.imp[::-1], xerr=imp.sd[::-1], color="#2a78d6")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_xlabel("permutation importance (drop in 2010 AUC)")
    ax.set_title("what carries p_real (held-out year)")
    fig.tight_layout(); fig.savefig(f"{FIG}/e_feature_importance.png", dpi=140); plt.close(fig)
    rep["importance_2010"] = imp.to_dict("records")

    import joblib
    joblib.dump(dict(model=model, iso=iso, feats=FEATS, kind=kind, shift_2010=shift,
                     model_safe=msafe, iso_safe=isos, feats_safe=SAFE,
                     prior_train=float(A12.match_own.mean())), f"{OUT}/model_v1.joblib")
    json.dump(rep, open(f"{OUT}/train_report.json", "w"), indent=1, default=float)
    print("[done] model + report written", flush=True)


if __name__ == "__main__":
    main()
