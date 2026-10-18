#!/usr/bin/env python
"""STEP 1 (Merlin, 2026-09-30) — per-repeat probability by LOCAL FALSE-DISCOVERY RATE, with the site-rolled
scrambled controls as the empirical noise model.  Every candidate keeps its row; P is a column (rank, not cut).

Method (one simple, stable density-ratio estimate):
  f(x)  = feature density of DATA candidates, f0(x) = feature density of SCRAMBLED-CONTROL candidates
          (same days, same families, same record tier >= 8.0 MAD; K = 4 independent site-rolled controls).
  A gradient-boosted classifier (sklearn HistGradientBoosting, isotonic-calibrated by 3-fold CalibratedClassifierCV
  with folds grouped by DAY) estimates p(x) = P(data | x) on the fit set.  Then
        f0(x)/f(x) = [(1-p)/p] * N_data / N_ctl.
  pi0 (fraction of data candidates that are noise) from QUIET-DAY behaviour: on quiet days (tremor-free by the
  window design: V2010 'quiet', V2008 'post') every data candidate is taken as noise, so
        r_q = N_data(quiet) / [N_ctl(quiet) / K]     (per region; global if < 500 quiet data candidates)
  calibrates the control count to the data's noise count; pi0 = r_q (N_ctl/K) / N_data.  Hence
        P(real | x) = 1 - pi0 f0/f = 1 - r_q (1-p)/(p K),  clipped to [0, 1].
  (Real LFEs on 'quiet' days are counted as noise -> conservative.)
Fit set = BUILD half (wa_common.half_of == 0) of V2008 + V2010 (all series).  Held out: V2008 / V2010 TEST half,
A2011 (independent episode, all).  Leave-one-control-out refits (K = 3) give (i) a model spread for the rate
product and (ii) P for a control realization that never entered its model.

Density features (same series only): dens_own_1h / dens_own_10m = other candidates of the same family within
+-1 h / +-10 min; dens_grp_1h = candidates of the family's <= 3-km group (incl. itself) within +-1 h (minus self).

  python scripts/step1_fit.py            -> results/step1_prob/{candidates.parquet, calibration.csv, ...}
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import wa_score as WS

FEATDIR = os.environ.get("STEP1_FEATDIR", f"{W.RES}/step1_prob/feat")   # STEP 2 group refit overrides
MF1 = os.environ.get("STEP1_MF1", "mt1_all")
OD = os.environ.get("STEP1_OUT", f"{W.RES}/step1_prob")   # smoke: STEP1_OUT=results/step1_prob/smoke
LOG = f"{W.STUDY}/logs/step1.log"
CTLS = ["c1", "c2", "c3", "c4"]
FEATS = ["net_mad", "net_cc", "n_sta", "n_sites", "wsum", "site_cc_max", "site_cc_2nd", "site_z_max", "site_z_2nd",
         "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top", "hfr_rel", "dens_own_1h", "dens_own_10m",
         "dens_grp_1h", "site_noise_mad_med", "region_code"]
# PNET refit (2026-09-30): extra feature columns present in STEP1_FEATDIR's parquet files (default: none)
FEATS += [c for c in os.environ.get("STEP1_EXTRA_FEATS", "").split(",") if c]
REGC = {r: i for i, r in enumerate(W.REGIONS)}
BINS = np.linspace(0, 1, 11)


def near(x, ref, tol=W.MATCH_S):
    x = np.asarray(x, float)
    if len(ref) == 0 or len(x) == 0:
        return np.zeros(len(x), bool)
    j = np.clip(np.searchsorted(ref, x), 1, max(len(ref) - 1, 1))
    if len(ref) == 1:
        return np.abs(ref[0] - x) <= tol
    return np.minimum(np.abs(ref[j] - x), np.abs(ref[j - 1] - x)) <= tol


def count_within(tq, tref, h):
    """number of tref (sorted) within [tq-h, tq+h]"""
    return np.searchsorted(tref, tq + h, "right") - np.searchsorted(tref, tq - h, "left")


def fit_model(X, y, groups, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.model_selection import GroupKFold
    base = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=500,
                                          l2_regularization=1.0, early_stopping=True, validation_fraction=0.1,
                                          n_iter_no_change=20, categorical_features=[FEATS.index("region_code")],
                                          random_state=seed)
    cv = list(GroupKFold(3).split(X, y, groups))
    m = CalibratedClassifierCV(base, method="isotonic", cv=cv)
    m.fit(X, y)
    return m


def p_real(m, X, rq_row, K):
    p = m.predict_proba(X)[:, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        odds_c = np.where(p > 0, (1 - p) / p, np.inf)
    return np.clip(1 - rq_row * odds_c / K, 0, 1), p


def rq_table(A, fitmask, ctls):
    K = len(ctls)
    q = A[fitmask & A.quiet & A.series.isin(["data"] + ctls)]
    out = {}
    nd_all = int((q.series == "data").sum()); nc_all = int(q.series.isin(ctls).sum())
    g = nd_all / (nc_all / K) if nc_all else 1.0
    for rg in W.REGIONS:
        qq = q[q.region == rg]
        nd = int((qq.series == "data").sum()); nc = int(qq.series.isin(ctls).sum())
        out[rg] = dict(n_data_quiet=nd, n_ctl_quiet=nc, r_q=(nd / (nc / K)) if (nd >= 500 and nc > 0) else g,
                       source="region" if (nd >= 500 and nc > 0) else "global")
    out["ALL"] = dict(n_data_quiet=nd_all, n_ctl_quiet=nc_all, r_q=g, source="global")
    return out


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(OD, exist_ok=True)
    Fm = W.load_families()
    A = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{FEATDIR}/feat_*.parquet"))], ignore_index=True)
    wk = {d: (w, k) for w in ["V2008", "V2010", "A2011"] for d, k in W.window_days(w)}
    A["window"] = A.date.map(lambda d: wk[d][0]); A["kind"] = A.date.map(lambda d: wk[d][1])
    A["quiet"] = A.kind.isin(["quiet", "post"])
    A["half"] = W.half_of(A.t.to_numpy())
    reg = dict(zip(Fm.fid, Fm.region))
    A["region"] = A.fid.map(reg); A["region_code"] = A.region.map(REGC)
    A["frac_z3"] = A.n_z3 / A.n_sites
    # tremor-day label (label only)
    tdays = {}
    for rg, tr in WS.TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        tdays[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    A["tremor_day_label"] = [d in tdays[r] for d, r in zip(A.date, A.region)]
    # ---------------- density features (same series)
    la, lo = Fm.lat_use.to_numpy(float), Fm.lon_use.to_numpy(float)
    dk = W.km(la[:, None], lo[:, None], la[None], lo[None])
    NB = {}
    for i, f in enumerate(Fm.fid):
        mm = np.nan_to_num(dk[i] <= 3.0, nan=False).astype(bool); mm[i] = True
        NB[int(f)] = [int(x) for x in Fm.fid[mm]]
    A = A.sort_values(["series", "fid", "t"]).reset_index(drop=True)
    A["dens_own_1h"] = 0; A["dens_own_10m"] = 0; A["dens_grp_1h"] = 0
    for s, gs in A.groupby("series"):
        TT = {int(f): g.t.to_numpy() for f, g in gs.groupby("fid")}
        for f, g in gs.groupby("fid"):
            tq = g.t.to_numpy()
            A.loc[g.index, "dens_own_1h"] = count_within(tq, TT[f], 3600) - 1
            A.loc[g.index, "dens_own_10m"] = count_within(tq, TT[f], 600) - 1
            A.loc[g.index, "dens_grp_1h"] = sum(count_within(tq, TT[h], 3600) for h in NB[int(f)] if h in TT) - 1
    W.log(LOG, f"[fit] {len(A)} candidates {A.series.value_counts().to_dict()}; features {time.time()-t0:.0f}s")
    # ---------------- fit set / held-out sets
    ctls = [c for c in CTLS if (A.series == c).any()]
    K = len(ctls)
    A["fit"] = A.window.isin(["V2008", "V2010"]) & (A.half == 0)
    A["eval_set"] = np.select([A.window.isin(["V2008", "V2010"]) & (A.half == 1), A.window == "A2011",
                               A.fit], [A.window + "_test", "A2011", "fit_in_sample"], "guard")
    X = A[FEATS].to_numpy(float)
    fm = A.fit.to_numpy() & A.series.isin(["data"] + ctls).to_numpy()
    y = (A.series == "data").to_numpy().astype(int)
    RQ = rq_table(A, A.fit, ctls)
    W.log(LOG, f"[fit] K={K} r_q {json.dumps(RQ)}")
    m = fit_model(X[fm], y[fm], A.date.to_numpy()[fm])
    rq_row = A.region.map({r: RQ[r]["r_q"] for r in W.REGIONS}).to_numpy(float)
    A["P"], A["p_data"] = p_real(m, X, rq_row, K)
    W.log(LOG, f"[fit] main model {time.time()-t0:.0f}s")
    # univariate separation (fit set) + permutation importance (subsample)
    from sklearn.metrics import roc_auc_score
    imp = []
    sub = np.nonzero(fm)[0]
    rs = np.random.default_rng(0)
    ss = rs.choice(sub, min(len(sub), 200000), replace=False)
    base_ll = None
    from sklearn.metrics import log_loss
    pb = m.predict_proba(X[ss])[:, 1]; base_ll = log_loss(y[ss], np.clip(pb, 1e-6, 1 - 1e-6))
    for j, fn in enumerate(FEATS):
        v = X[fm, j]; ok = np.isfinite(v)
        auc = roc_auc_score(y[fm][ok], v[ok]) if ok.sum() > 100 and len(set(y[fm][ok])) == 2 else np.nan
        Xp = X[ss].copy(); Xp[:, j] = rs.permutation(Xp[:, j])
        ll = log_loss(y[ss], np.clip(m.predict_proba(Xp)[:, 1], 1e-6, 1 - 1e-6))
        imp.append(dict(feature=fn, auc_data_vs_ctl=auc, perm_logloss_gain=ll - base_ll, frac_nan=float((~ok).mean())))
    pd.DataFrame(imp).sort_values("perm_logloss_gain", ascending=False).to_csv(f"{OD}/feature_importance.csv", index=False)
    # ---------------- leave-one-control-out refits
    for c in ctls:
        keep = [x for x in ctls if x != c]
        fmk = fm & (A.series != c).to_numpy()
        RQk = rq_table(A, A.fit, keep)
        mk = fit_model(X[fmk], y[fmk], A.date.to_numpy()[fmk], seed=1)
        rqk = A.region.map({r: RQk[r]["r_q"] for r in W.REGIONS}).to_numpy(float)
        A[f"P_loco_{c}"], _ = p_real(mk, X, rqk, len(keep))
        RQ[f"loco_{c}"] = {r: RQk[r]["r_q"] for r in RQk}
        W.log(LOG, f"[fit] loco {c} {time.time()-t0:.0f}s")
    json.dump(dict(K=K, controls=ctls, r_q=RQ, features=FEATS,
                   n_fit=dict(A[fm].series.value_counts()), method="HistGradientBoosting + isotonic (3-fold by day)"),
              open(f"{OD}/model_info.json", "w"), indent=1, default=int)
    # ---------------- truth matches
    T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    scan = pd.read_csv(f"{W.SETUP}/scan_days.csv")
    sd = {s: set(g.date) for s, g in scan.groupby("det_source")}
    src = dict(zip(Fm.fid, Fm.det_source)); vw = dict(zip(Fm.fid, Fm.val_windows))
    A["truth_ok"] = [(d in sd.get(src[f], set())) and (w in vw[f].split(",")) for f, d, w in zip(A.fid, A.date, A.window)]
    Tday = {d: g for d, g in T.groupby("date")}
    for lev in ["own", "nbr3", "any"]:
        A[f"m_{lev}"] = False
    for (d, f), g in A[A.truth_ok].groupby(["date", "fid"]):
        tg = Tday.get(d)
        if tg is None:
            continue
        tt = g.t.to_numpy()
        A.loc[g.index, "m_own"] = near(tt, np.sort(tg[tg.fid == f].t.to_numpy()))
        A.loc[g.index, "m_nbr3"] = near(tt, np.sort(tg[tg.fid.isin(NB[int(f)])].t.to_numpy()))
        A.loc[g.index, "m_any"] = near(tt, np.sort(tg.t.to_numpy()))
    A.to_parquet(f"{OD}/candidates.parquet", index=False)
    W.log(LOG, f"[fit] candidates written {time.time()-t0:.0f}s")
    # ---------------- (a) calibration
    A["pbin"] = np.clip(np.digitize(A.P, BINS[1:-1]), 0, 9)
    ev = A[A.truth_ok & A.eval_set.isin(["V2010_test", "V2008_test", "A2011", "fit_in_sample"])]
    rows = []
    for es, ge in ev.groupby("eval_set"):
        for rg in ["ALL"] + W.REGIONS:
            gr = ge if rg == "ALL" else ge[ge.region == rg]
            for s in ["data", "control"]:
                gs = gr[gr.series == "data"] if s == "data" else gr[gr.series.isin(ctls)]
                for b, gb in gs.groupby("pbin"):
                    rows.append(dict(eval_set=es, region=rg, series=s, pbin=f"{BINS[b]:.1f}-{BINS[b+1]:.1f}", n=len(gb),
                                     mean_P=gb.P.mean(), frac_own=gb.m_own.mean(), frac_nbr3=gb.m_nbr3.mean(),
                                     frac_any=gb.m_any.mean()))
    CAL = pd.DataFrame(rows)
    CAL.to_csv(f"{OD}/calibration.csv", index=False)
    srow = []
    for (es, rg, s), g in CAL.groupby(["eval_set", "region", "series"]):
        g = g[g.n >= 20]
        r = dict(eval_set=es, region=rg, series=s, n=int(g.n.sum()), sum_P=float((g.mean_P * g.n).sum()))
        for lev in ["own", "nbr3", "any"]:
            r[f"matched_{lev}"] = float((g[f"frac_{lev}"] * g.n).sum())
            if len(g) >= 3 and g.mean_P.max() - g.mean_P.min() > 0.2:
                wv = g.n.to_numpy(float); xx = g.mean_P.to_numpy(); yy = g[f"frac_{lev}"].to_numpy()
                Xd = np.stack([np.ones_like(xx), xx], 1)
                beta = np.linalg.lstsq(Xd * np.sqrt(wv)[:, None], yy * np.sqrt(wv), rcond=None)[0]
                r[f"slope_{lev}"] = beta[1]; r[f"icpt_{lev}"] = beta[0]
            else:
                r[f"slope_{lev}"] = np.nan; r[f"icpt_{lev}"] = np.nan
            r[f"obs_over_exp_{lev}"] = r[f"matched_{lev}"] / r["sum_P"] if r["sum_P"] > 0 else np.nan
        srow.append(r)
    pd.DataFrame(srow).to_csv(f"{OD}/calibration_slopes.csv", index=False)
    # ---------------- (b) control P distribution (main model on held-out halves; LOCO = never-seen realization)
    rows = []
    for es in ["V2010_test", "V2008_test", "A2011", "fit_in_sample", "ALL"]:
        ge = A if es == "ALL" else A[A.eval_set == es]
        for rg in ["ALL"] + W.REGIONS:
            gr = ge if rg == "ALL" else ge[ge.region == rg]
            for kind in ["all", "quiet", "episode"]:
                gk = gr if kind == "all" else (gr[gr.quiet] if kind == "quiet" else gr[~gr.quiet])
                d = gk[gk.series == "data"]; c = gk[gk.series.isin(ctls)]
                lo_ = np.concatenate([gk[gk.series == cc][f"P_loco_{cc}"].to_numpy() for cc in ctls]) if len(c) else np.array([])
                if not len(d) and not len(c):
                    continue
                rows.append(dict(eval_set=es, region=rg, days=kind, n_data=len(d), n_ctl_per_real=len(c) / K,
                                 data_meanP=d.P.mean(), data_fracP_gt05=(d.P > 0.5).mean(), data_fracP_gt09=(d.P > 0.9).mean(),
                                 data_sumP=d.P.sum(),
                                 ctl_meanP=c.P.mean(), ctl_fracP_gt01=(c.P > 0.1).mean(), ctl_fracP_gt05=(c.P > 0.5).mean(),
                                 ctl_sumP_per_real=c.P.sum() / K,
                                 ctl_loco_fracP_gt05=(lo_ > 0.5).mean() if len(lo_) else np.nan,
                                 ctl_loco_fracP_gt01=(lo_ > 0.1).mean() if len(lo_) else np.nan,
                                 ctl_loco_sumP_per_real=lo_.sum() / K if len(lo_) else np.nan))
    pd.DataFrame(rows).to_csv(f"{OD}/control_and_quiet.csv", index=False)
    # tremor-label split (label only) for held-out sets
    rows = []
    for es in ["V2010_test", "A2011", "V2008_test"]:
        ge = A[A.eval_set == es]
        for rg in ["ALL"] + W.REGIONS:
            gr = ge if rg == "ALL" else ge[ge.region == rg]
            for lab, gl in gr.groupby("tremor_day_label"):
                d = gl[gl.series == "data"]; c = gl[gl.series.isin(ctls)]
                nfd = d.groupby(["fid", "date"]).ngroups
                rows.append(dict(eval_set=es, region=rg, tremor_day=bool(lab), n_famdays=nfd,
                                 data_sumP_per_famday=d.P.sum() / nfd if nfd else np.nan,
                                 ctl_sumP_per_famday=c.P.sum() / K / nfd if nfd else np.nan,
                                 data_fracP_gt05=(d.P > 0.5).mean(), ctl_fracP_gt05=(c.P > 0.5).mean()))
    pd.DataFrame(rows).to_csv(f"{OD}/tremor_label_split.csv", index=False)
    # ---------------- rate product: per family and day, sum P (+ uncertainty)
    D = A[A.series == "data"].copy()
    D["scored"] = D.eval_set.isin(["V2010_test", "V2008_test", "A2011"])
    D["pp"] = D.P * (1 - D.P)
    g = D.groupby(["fid", "date"])
    R = pd.DataFrame(dict(n_cand=g.size(), sumP=g.P.sum(), var_binom=g.pp.sum(),
                          sumP_scored_half=g.apply(lambda x: x.P[x.scored].sum(), include_groups=False),
                          n_cand_scored_half=g.scored.sum()))
    for c in ctls:
        R[f"sumP_loco_{c}"] = g[f"P_loco_{c}"].sum()
    lc = R[[f"sumP_loco_{c}" for c in ctls]]
    R["sd_model"] = lc.std(axis=1, ddof=1)
    R["sd_total"] = np.sqrt(R.var_binom + R.sd_model.fillna(0) ** 2)
    Ca = A[A.series.isin(ctls)]
    Cc = Ca.groupby(["fid", "date"]).P.sum() / K
    R["sumP_ctl_per_real"] = Cc.reindex(R.index).fillna(0)
    Cs = Ca[Ca.eval_set.isin(["V2010_test", "V2008_test", "A2011"])].groupby(["fid", "date"]).P.sum() / K
    R["sumP_ctl_scored_half"] = Cs.reindex(R.index).fillna(0)
    R = R.reset_index()
    FD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{W.RES}/mf_{MF1}/famday_*.parquet"))])
    FD = FD[FD.observable][["fid", "date", "window", "kind", "n_sites"]]
    R = FD.merge(R, on=["fid", "date"], how="left")
    for c in ["n_cand", "sumP", "var_binom", "sumP_scored_half", "n_cand_scored_half", "sumP_ctl_per_real", "sumP_ctl_scored_half"]:
        R[c] = R[c].fillna(0)
    R["region"] = R.fid.map(reg); R["family_uid"] = R.fid.map(dict(zip(Fm.fid, Fm.family_uid)))
    R["tremor_day_label"] = [d in tdays[r] for d, r in zip(R.date, R.region)]
    # truth counts (scored half; only where the family's own run scanned the day)
    R["truth_ok"] = [(d in sd.get(src[f], set())) and (w in vw[f].split(",")) for f, d, w in zip(R.fid, R.date, R.window)]
    Ts = T.copy(); Ts["scored"] = (Ts.half == 1) | (Ts.window == "A2011")
    Ts = Ts[Ts.scored]
    own = Ts.groupby(["fid", "date"]).size()
    R["truth_own_scored"] = [own.get((f, d), 0) for f, d in zip(R.fid, R.date)]
    TsD = {d: g for d, g in Ts.groupby("date")}
    grp = []
    for f, d in zip(R.fid, R.date):
        tg = TsD.get(d)
        if tg is None:
            grp.append(0); continue
        tt = np.sort(tg[tg.fid.isin(NB[int(f)])].t.to_numpy())
        grp.append(int(1 + (np.diff(tt) > W.MATCH_S).sum()) if len(tt) else 0)
    R["truth_nbr3_scored"] = grp
    R.to_parquet(f"{OD}/rate_family_day.parquet", index=False)
    # ---------------- (e) expected-count check
    from scipy.stats import spearmanr
    rows = []
    E = R[R.truth_ok & (R.kind == "episode")]
    for es, ge in [("V2010_test", E[E.window == "V2010"]), ("V2008_test", E[E.window == "V2008"]), ("A2011", E[E.window == "A2011"])]:
        for rg in ["ALL"] + W.REGIONS:
            gr = ge if rg == "ALL" else ge[ge.region == rg]
            if len(gr) < 5:
                continue
            r = dict(eval_set=es, region=rg, n_famdays=len(gr), sumP=gr.sumP_scored_half.sum(),
                     truth_own=int(gr.truth_own_scored.sum()), truth_nbr3=int(gr.truth_nbr3_scored.sum()),
                     sumP_ctl=float(gr.sumP_ctl_scored_half.sum()))
            r["ratio_P_over_own"] = r["sumP"] / r["truth_own"] if r["truth_own"] else np.nan
            r["ratio_P_over_nbr3"] = r["sumP"] / r["truth_nbr3"] if r["truth_nbr3"] else np.nan
            r["spearman_own"] = spearmanr(gr.sumP_scored_half, gr.truth_own_scored).statistic
            r["spearman_nbr3"] = spearmanr(gr.sumP_scored_half, gr.truth_nbr3_scored).statistic
            r["spearman_ncand_own"] = spearmanr(gr.n_cand_scored_half, gr.truth_own_scored).statistic
            rows.append(r)
    pd.DataFrame(rows).to_csv(f"{OD}/expected_count.csv", index=False)
    W.log(LOG, f"[fit] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
