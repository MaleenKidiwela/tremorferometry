#!/usr/bin/env python
"""CONFIDENCE (2026-09-30, user-approved) — a per-repeat confidence fitted against the KNOWN ANSWER (array-era truth)
instead of data-vs-scrambled-control.

Candidates  the step-1 candidates (T1 all-templates MF, record tier >= 8.0 MAD) on V2008 / V2010 / A2011:
            results/step1_prob/candidates.parquet (series data, c1..c4 = site-rolled scrambled controls).
Features    step-1 features (21, incl. region_code) + the step-1 lfdr P ('P1'; cross-fitted by day on the build half so
            no row sees a P from a model fitted on its own day) + PNET off-the-shelf (_w) and WA fine-tuned (_f):
            pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05, pnet_rel_max + step-2 attribution where available:
            grp_size (members in the family's <= 3-km group), P_member (softmax attribution of the nearest group-MF
            candidate within +-2 s in the group frame t + delta_f; NaN if none / singleton), n_elig.
            Tremor-day labels are NOT features.  Step-2 attribution exists for data and the first control (c1) only.
Labels      data candidates on days the family's array run scanned (truth_ok), matched within +-2 s (nearest) to an
            array detection of: own family (m_own) / any family of the <= 3-km group (m_nbr3) / any family (m_any).
Model       HistGradientBoosting (shallow, strongly regularised) + isotonic calibration (CalibratedClassifierCV, 3
            folds grouped by day); logistic regression (+ isotonic) for comparison.
Designs     D1 (as specified): fit on the BUILD half (half 0) of V2008 + V2010; evaluate V2010 / V2008 TEST half and
               A2011.
            D2 (template-leak-free): the templates are stacks of BUILD-half array detections, so on the build half a
               truth-matched candidate correlates with a template that contains that very event (inflated features;
               see the leakage table).  D2 fits on the TEST half with 5 outer folds grouped by day (every test-half
               row, data and control, scored by a model that never saw its day) and evaluates A2011 with a model fitted
               on the whole test half.  A2011 is independent in both designs.
Checks      (a) AUC / AP by region vs step-1 P, PNET alone, net_mad alone; (b) reliability (WLS slope on 10 P bins);
            (c) the model applied to scrambled-control candidates; (d) quiet vs episode days; (e) sum P per
            family-day vs the true array count; (f) permutation importance + ablations; held-out-station PNET.

  WA_SETUP=results/merlin_tests/t1_setup python scripts/confidence_fit.py  -> results/confidence/*
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OD = f"{W.RES}/confidence"
LOG = f"{W.STUDY}/logs/confidence.log"
RS = {"Royer_NW_WA": "NW", "Royer_SW_WA": "SW", "Olympic_2010": "Olympic", "Olympic_2010_armB": "Arm B"}
S1F = ["net_mad", "net_cc", "n_sta", "n_sites", "wsum", "site_cc_max", "site_cc_2nd", "site_z_max", "site_z_2nd",
       "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top", "hfr_rel", "dens_own_1h", "dens_own_10m",
       "dens_grp_1h", "site_noise_mad_med", "region_code"]
PN = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_rel_max"]
PNW = [c + "_w" for c in PN]; PNF = [c + "_f" for c in PN]
PNHO = [c + "_ho_w" for c in PN] + [c + "_ho_f" for c in PN]
S2F = ["grp_size", "P_member", "n_elig"]
FSETS = {
    "FULL": S1F + ["P1"] + PNW + PNF + S2F,
    "NO_PNET": S1F + ["P1"] + S2F,
    "NO_STEP2": S1F + ["P1"] + PNW + PNF,
    "NO_P1": S1F + PNW + PNF + S2F,
    "PNET_HO": S1F + ["P1"] + PNHO + S2F,          # PNET read only at the 6 never-fine-tuned stations
    "S1_ONLY": S1F,
}
LABELS = {"own": "m_own", "nbr3": "m_nbr3", "any": "m_any"}
CTLS = ["c1", "c2", "c3", "c4"]
BINS = np.linspace(0, 1, 11)
EVAL = ["V2010_test", "V2008_test", "A2011"]


def log(s):
    W.log(LOG, s)


def near_idx(x, ref, tol=W.MATCH_S):
    x = np.asarray(x, float); ref = np.asarray(ref, float)
    out = np.full(len(x), -1)
    if len(ref) == 0 or len(x) == 0:
        return out
    j = np.clip(np.searchsorted(ref, x), 1, max(len(ref) - 1, 1))
    if len(ref) == 1:
        return np.where(np.abs(ref[0] - x) <= tol, 0, -1)
    jl = j - 1
    d1 = np.abs(ref[j] - x); d0 = np.abs(ref[jl] - x)
    k = np.where(d1 < d0, j, jl); d = np.minimum(d0, d1)
    return np.where(d <= tol, k, -1)


# ------------------------------------------------------------------ data assembly
def assemble():
    t0 = time.time()
    A = pd.read_parquet(f"{W.RES}/step1_prob/candidates.parquet")
    A = A.rename(columns={"P": "P1_orig"})
    for tag, suf in [("w05", "_w"), ("ft1", "_f")]:
        B = pd.read_parquet(f"{W.RES}/pnet/step1_{tag}/candidates.parquet", columns=["fid", "t", "series"] + PN)
        assert (B.fid.values == A.fid.values).all() and (B.t.values == A.t.values).all() and (B.series.values == A.series.values).all()
        for c in PN:
            A[c + suf] = B[c].to_numpy()
        H = pd.concat([pd.read_parquet(p) for p in sorted(__import__("glob").glob(f"{OD}/pnet_ho_{tag}/feat_*.parquet"))])
        H = H.rename(columns={c + "_ho": c + "_ho" + suf for c in PN})
        A = A.merge(H, on=["fid", "t", "series"], how="left")
    log(f"[asm] PNET merged {time.time()-t0:.0f}s, {len(A)} rows")
    # step-2 attribution
    M = pd.read_csv(f"{W.RES}/step2_groups/members.csv")
    gsz = M.groupby("gid").fid.transform("size")
    A["grp_size"] = A.fid.map(dict(zip(M.fid, gsz))).astype(float)
    dl = dict(zip(M.fid, M.delta_s)); gof = dict(zip(M.fid, M.gid))
    AP = pd.read_parquet(f"{W.RES}/step2_groups/attrib_pairs.parquet", columns=["gid", "fid", "t", "series", "n_elig", "P_member"])
    AP["series"] = AP.series.map({"data": "data", "control": "c1"})
    A["P_member"] = np.nan; A["n_elig"] = np.nan
    APg = {k: g.sort_values("t") for k, g in AP.groupby(["series", "fid"])}
    for (s, f), g in A[A.series.isin(["data", "c1"])].groupby(["series", "fid"]):
        ap = APg.get((s, f))
        if ap is None:
            continue
        tg = g.t.to_numpy() + dl[f]
        k = near_idx(tg, ap.t.to_numpy())
        ok = k >= 0
        A.loc[g.index[ok], "P_member"] = ap.P_member.to_numpy()[k[ok]]
        A.loc[g.index[ok], "n_elig"] = ap.n_elig.to_numpy()[k[ok]]
    log(f"[asm] step-2 merged {time.time()-t0:.0f}s; P_member defined data {A[A.series=='data'].P_member.notna().mean():.3f} "
        f"c1 {A[A.series=='c1'].P_member.notna().mean():.3f}")
    A["reg"] = A.region.map(RS)
    return A


def crossfit_P1(A):
    """step-1 lfdr refitted on half of the build-half DAYS, predicted on the other half (both directions).
    Build-half rows get the cross-fitted P; test half / A2011 keep the original step-1 P (fitted on the build half)."""
    import step1_fit as S1
    t0 = time.time()
    A["P1"] = A.P1_orig.to_numpy()
    ctls = CTLS
    K = len(ctls)
    days = np.array(sorted(A[A.fit].date.unique()))
    rng = np.random.default_rng(7)
    fold = dict(zip(days, rng.permutation(np.arange(len(days)) % 2)))
    A["_bfold"] = A.date.map(fold).fillna(-1).astype(int)
    X = A[S1.FEATS].to_numpy(float)
    y = (A.series == "data").to_numpy().astype(int)
    for k in (0, 1):
        fitm = (A.fit & (A._bfold == k)).to_numpy()
        prm = (A.fit & (A._bfold == 1 - k)).to_numpy()
        RQ = S1.rq_table(A, pd.Series(fitm, index=A.index), ctls)
        m = S1.fit_model(X[fitm], y[fitm], A.date.to_numpy()[fitm])
        rq = A.region.map({r: RQ[r]["r_q"] for r in W.REGIONS}).to_numpy(float)
        P, _ = S1.p_real(m, X[prm], rq[prm], K)
        A.loc[prm, "P1"] = P
        log(f"[p1cf] fold {k}: fit {fitm.sum()} predict {prm.sum()} {time.time()-t0:.0f}s")
    A.drop(columns="_bfold", inplace=True)
    return A


# ------------------------------------------------------------------ models
def hgb(feats, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=300,
                                          l2_regularization=2.0, early_stopping=True, validation_fraction=0.15,
                                          n_iter_no_change=20, categorical_features=[feats.index("region_code")],
                                          random_state=seed)


def logit_pipe(feats):
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler, OneHotEncoder, FunctionTransformer
    from sklearn.compose import ColumnTransformer
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    ir = feats.index("region_code")
    num = [i for i in range(len(feats)) if i != ir]
    cnt = {"n_sta", "n_sites", "n_z2", "n_z3", "n_z5", "dens_own_1h", "dens_own_10m", "dens_grp_1h", "grp_size", "n_elig"}
    lg = [i for i in num if feats[i] in cnt]
    oth = [i for i in num if feats[i] not in cnt]
    pre = ColumnTransformer([
        ("lg", make_pipeline(SimpleImputer(strategy="median", add_indicator=True), FunctionTransformer(np.log1p), StandardScaler()), lg),
        ("num", make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler()), oth),
        ("reg", OneHotEncoder(handle_unknown="ignore"), [ir])])
    return make_pipeline(pre, LogisticRegression(C=0.1, max_iter=2000))


def fit_cal(kind, feats, X, y, groups, seed=0):
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.model_selection import GroupKFold
    base = hgb(feats, seed) if kind == "hgb" else logit_pipe(feats)
    cv = list(GroupKFold(3).split(X, y, groups))
    m = CalibratedClassifierCV(base, method="isotonic", cv=cv)
    m.fit(X, y)
    return m


def run_design(A, design, kind, fset, lab):
    """returns predictions for every row of A that the design scores (NaN elsewhere) and the fitted model(s)"""
    feats = FSETS[fset]
    X = A[feats].to_numpy(float)
    y = A[LABELS[lab]].to_numpy().astype(int)
    lab_ok = ((A.series == "data") & A.truth_ok).to_numpy()
    vv = A.window.isin(["V2008", "V2010"]).to_numpy()
    pred = np.full(len(A), np.nan)
    dates = A.date.to_numpy()
    models = []
    if design == "D1":
        tr = lab_ok & vv & (A.half == 0).to_numpy()
        m = fit_cal(kind, feats, X[tr], y[tr], dates[tr])
        sc = (vv & (A.half == 1).to_numpy()) | (A.window == "A2011").to_numpy()
        pred[sc] = m.predict_proba(X[sc])[:, 1]
        models.append(m)
    else:
        tst = vv & (A.half == 1).to_numpy()
        tdays = np.array(sorted(set(dates[tst])))
        rng = np.random.default_rng(11)
        fo = dict(zip(tdays, rng.permutation(np.arange(len(tdays)) % 5)))
        fold = np.array([fo.get(d, -1) for d in dates])
        for k in range(5):
            tr = lab_ok & tst & (fold != k)
            m = fit_cal(kind, feats, X[tr], y[tr], dates[tr])
            sc = tst & (fold == k)
            pred[sc] = m.predict_proba(X[sc])[:, 1]
        tr = lab_ok & tst
        m = fit_cal(kind, feats, X[tr], y[tr], dates[tr])
        sc = (A.window == "A2011").to_numpy()
        pred[sc] = m.predict_proba(X[sc])[:, 1]
        models.append(m)
    return pred, models


# ------------------------------------------------------------------ metrics
def auc_ap(y, s):
    from sklearn.metrics import roc_auc_score, average_precision_score
    ok = np.isfinite(s)
    y = y[ok]; s = s[ok]
    if len(y) < 20 or y.sum() < 5 or y.sum() == len(y):
        return np.nan, np.nan
    return roc_auc_score(y, s), average_precision_score(y, s)


def calib(y, p):
    b = np.clip(np.digitize(p, BINS[1:-1]), 0, 9)
    rows = []
    for k in range(10):
        mm = b == k
        if mm.sum():
            rows.append((k, mm.sum(), p[mm].mean(), y[mm].mean()))
    R = pd.DataFrame(rows, columns=["bin", "n", "mean_P", "frac"])
    g = R[R.n >= 20]
    slope = icpt = np.nan
    if len(g) >= 3 and g.mean_P.max() - g.mean_P.min() > 0.2:
        wv = g.n.to_numpy(float); xx = g.mean_P.to_numpy(); yy = g.frac.to_numpy()
        Xd = np.stack([np.ones_like(xx), xx], 1)
        beta = np.linalg.lstsq(Xd * np.sqrt(wv)[:, None], yy * np.sqrt(wv), rcond=None)[0]
        icpt, slope = beta
    return slope, icpt, R


def eval_set_mask(A, es):
    if es == "A2011":
        return (A.window == "A2011").to_numpy()
    w = es.split("_")[0]
    return ((A.window == w) & (A.half == 1)).to_numpy()


def main():
    t0 = time.time()
    os.makedirs(OD, exist_ok=True)
    A = assemble()
    A = crossfit_P1(A)
    # nearest station of each family (for the held-out-station subset)
    Fm = W.load_families(); S = W.load_stations()
    HO = ["UW.HDW", "UW.OSD", "PB.B018", "PB.B017", "US.NLWA", "UW.GNW"]
    near_sta = {}
    for r in Fm.itertuples():
        if np.isfinite(r.lat_use) and np.isfinite(r.lon_use):
            d = W.km(r.lat_use, r.lon_use, S.lat.to_numpy(), S.lon.to_numpy())
            near_sta[int(r.fid)] = S.netsta.to_numpy()[np.argsort(d)[:2]].tolist()
    A["near_ho"] = A.fid.map(lambda f: any(s in HO for s in near_sta.get(int(f), []))).astype(bool)
    log(f"[main] assembled {time.time()-t0:.0f}s; families with a held-out station among their 2 nearest: "
        f"{sum(any(s in HO for s in v) for v in near_sta.values())}/{len(near_sta)}")
    # ---------------- leakage table: build vs test separation of raw features (own label)
    rows = []
    D = A[(A.series == "data") & A.truth_ok]
    for w in ["V2008", "V2010"]:
        for h, hn in [(0, "build"), (1, "test")]:
            for rg in ["ALL"] + list(RS.values()):
                g = D[(D.window == w) & (D.half == h)]
                g = g if rg == "ALL" else g[g.reg == rg]
                if len(g) < 50:
                    continue
                r = dict(window=w, half=hn, region=rg, n=len(g), frac_own=g.m_own.mean(), frac_any=g.m_any.mean())
                for c in ["net_mad", "n_z2", "site_cc_max", "P1_orig", "P1", "pnet_max_f", "pnet_max_w"]:
                    r[f"auc_own_{c}"] = auc_ap(g.m_own.to_numpy().astype(int), g[c].to_numpy(float))[0]
                rows.append(r)
    pd.DataFrame(rows).to_csv(f"{OD}/leakage_build_vs_test.csv", index=False)
    # ---------------- fits
    runs = [("D1", "hgb", "FULL"), ("D2", "hgb", "FULL"), ("D1", "logit", "FULL"), ("D2", "logit", "FULL"),
            ("D2", "hgb", "NO_PNET"), ("D2", "hgb", "NO_STEP2"), ("D2", "hgb", "NO_P1"), ("D2", "hgb", "PNET_HO"),
            ("D2", "hgb", "S1_ONLY"), ("D1", "hgb", "NO_PNET"), ("D1", "hgb", "NO_STEP2")]
    PRED = {}
    keep_models = {}
    for design, kind, fset in runs:
        for lab in LABELS:
            key = f"{design}_{kind}_{fset}_{lab}"
            p, ms = run_design(A, design, kind, fset, lab)
            PRED[key] = p
            if kind == "hgb" and fset in ("FULL", "NO_STEP2"):
                keep_models[key] = ms[-1]
            log(f"[fit] {key} {time.time()-t0:.0f}s")
    Pdf = pd.DataFrame(PRED)
    OUTC = ["fid", "t", "series", "date", "window", "kind", "quiet", "half", "region", "reg", "tremor_day_label",
            "truth_ok", "m_own", "m_nbr3", "m_any", "near_ho", "P1_orig", "P1", "net_mad", "pnet_max_w", "pnet_max_f",
            "pnet_max_ho_w", "pnet_max_ho_f", "grp_size", "P_member"]
    OUT = pd.concat([A[OUTC].reset_index(drop=True), Pdf], axis=1)
    OUT.to_parquet(f"{OD}/candidates_confidence.parquet", index=False)
    log(f"[main] predictions written {time.time()-t0:.0f}s")
    # ---------------- (a) AUC / AP by region, (b) calibration
    rows, crow, relrows = [], [], []
    D = OUT[(OUT.series == "data") & OUT.truth_ok]
    for es in EVAL:
        ge = D[eval_set_mask(D, es)]
        for rg in ["ALL"] + list(RS.values()):
            gr = ge if rg == "ALL" else ge[ge.reg == rg]
            for sub in ["all", "near_ho"]:
                gs = gr if sub == "all" else gr[gr.near_ho]
                if len(gs) < 50:
                    continue
                for lab, lc in LABELS.items():
                    y = gs[lc].to_numpy().astype(int)
                    r = dict(eval_set=es, region=rg, subset=sub, label=lab, n=len(gs), n_pos=int(y.sum()), base_rate=y.mean())
                    for sname, s in [("P1", gs.P1.to_numpy()), ("pnet_max_f", gs.pnet_max_f.to_numpy()),
                                     ("pnet_max_w", gs.pnet_max_w.to_numpy()), ("pnet_max_ho_f", gs.pnet_max_ho_f.to_numpy()),
                                     ("pnet_max_ho_w", gs.pnet_max_ho_w.to_numpy()), ("net_mad", gs.net_mad.to_numpy())]:
                        r[f"auc_{sname}"], r[f"ap_{sname}"] = auc_ap(y, s)
                    for design, kind, fset in runs:
                        key = f"{design}_{kind}_{fset}_{lab}"
                        p = gs[key].to_numpy()
                        r[f"auc_{design}_{kind}_{fset}"], r[f"ap_{design}_{kind}_{fset}"] = auc_ap(y, p)
                    rows.append(r)
                    if sub != "all":
                        continue
                    for design, kind, fset in runs + [("step1", "lfdr", "P1")]:
                        key = f"{design}_{kind}_{fset}_{lab}" if design != "step1" else "P1"
                        p = gs[key].to_numpy()
                        ok = np.isfinite(p)
                        if ok.sum() < 50:
                            continue
                        sl, ic, R = calib(y[ok], p[ok])
                        crow.append(dict(eval_set=es, region=rg, label=lab, model=f"{design}_{kind}_{fset}", n=int(ok.sum()),
                                         slope=sl, icpt=ic, sumP=float(p[ok].sum()), matched=int(y[ok].sum()),
                                         obs_over_exp=float(y[ok].sum() / p[ok].sum()) if p[ok].sum() > 0 else np.nan,
                                         brier=float(np.mean((p[ok] - y[ok]) ** 2))))
                        if (design, kind, fset) in [("D2", "hgb", "FULL"), ("D1", "hgb", "FULL"), ("step1", "lfdr", "P1")]:
                            R = R.assign(eval_set=es, region=rg, label=lab, model=f"{design}_{kind}_{fset}")
                            relrows.append(R)
    pd.DataFrame(rows).to_csv(f"{OD}/auc_ap.csv", index=False)
    pd.DataFrame(crow).to_csv(f"{OD}/calibration_slopes.csv", index=False)
    pd.concat(relrows).to_csv(f"{OD}/reliability.csv", index=False)
    log(f"[main] (a)(b) {time.time()-t0:.0f}s")
    # day-block bootstrap of the headline AUC difference (model - step-1 P), ALL regions
    rs = np.random.default_rng(3)
    brow = []
    from sklearn.metrics import roc_auc_score
    for es in EVAL:
        ge = D[eval_set_mask(D, es)]
        days = ge.date.unique()
        idx = {d: np.nonzero(ge.date.to_numpy() == d)[0] for d in days}
        for lab, lc in [("own", "m_own"), ("any", "m_any")]:
            for design in ["D1", "D2"]:
                key = f"{design}_hgb_FULL_{lab}"
                y = ge[lc].to_numpy().astype(int); p = ge[key].to_numpy(); q = ge.P1.to_numpy()
                diffs = []
                for b in range(200):
                    ii = np.concatenate([idx[d] for d in rs.choice(days, len(days), replace=True)])
                    if y[ii].sum() < 5:
                        continue
                    diffs.append(roc_auc_score(y[ii], p[ii]) - roc_auc_score(y[ii], q[ii]))
                brow.append(dict(eval_set=es, label=lab, design=design, n_days=len(days),
                                 dAUC_mean=np.mean(diffs), dAUC_lo=np.percentile(diffs, 2.5), dAUC_hi=np.percentile(diffs, 97.5)))
    pd.DataFrame(brow).to_csv(f"{OD}/bootstrap_dauc.csv", index=False)
    # ---------------- (c) scrambled control, (d) quiet vs episode
    rows = []
    for es in EVAL:
        ge = OUT[eval_set_mask(OUT, es)]
        for rg in ["ALL"] + list(RS.values()):
            gr = ge if rg == "ALL" else ge[ge.reg == rg]
            for days_ in ["all", "quiet", "episode"]:
                gk = gr if days_ == "all" else (gr[gr.quiet] if days_ == "quiet" else gr[~gr.quiet])
                if not len(gk):
                    continue
                for key in [f"{d}_hgb_{f}_{l}" for d in ["D1", "D2"] for f in ["FULL", "NO_STEP2"] for l in LABELS] + ["P1"]:
                    d = gk[gk.series == "data"]
                    cset = ["c1"] if "FULL" in key else CTLS
                    c = gk[gk.series.isin(cset)]
                    if not len(d):
                        continue
                    nfd = d.groupby(["fid", "date"]).ngroups
                    rows.append(dict(eval_set=es, region=rg, days=days_, model=key, controls=",".join(cset),
                                     n_data=len(d), n_ctl_per_real=len(c) / len(cset),
                                     data_meanP=d[key].mean(), data_fracP_gt05=(d[key] > 0.5).mean(),
                                     ctl_meanP=c[key].mean(), ctl_fracP_gt05=(c[key] > 0.5).mean(),
                                     data_sumP=d[key].sum(), ctl_sumP_per_real=c[key].sum() / len(cset),
                                     data_sumP_per_famday=d[key].sum() / nfd if nfd else np.nan,
                                     ctl_over_data_sumP=(c[key].sum() / len(cset)) / d[key].sum() if d[key].sum() > 0 else np.nan))
    pd.DataFrame(rows).to_csv(f"{OD}/control_quiet.csv", index=False)
    # tremor-day LABEL split (label only)
    rows = []
    for es in EVAL:
        ge = OUT[eval_set_mask(OUT, es)]
        for rg in ["ALL"] + list(RS.values()):
            gr = ge if rg == "ALL" else ge[ge.reg == rg]
            for lab_, gl in gr.groupby("tremor_day_label"):
                for key in ["D2_hgb_FULL_own", "D2_hgb_FULL_any", "D1_hgb_FULL_own", "D1_hgb_FULL_any", "P1"]:
                    d = gl[gl.series == "data"]; c = gl[gl.series == "c1"]
                    if not len(d):
                        continue
                    rows.append(dict(eval_set=es, region=rg, tremor_day=bool(lab_), model=key, n_data=len(d), n_c1=len(c),
                                     data_fracP_gt05=(d[key] > 0.5).mean(), c1_fracP_gt05=(c[key] > 0.5).mean(),
                                     data_meanP=d[key].mean(), c1_meanP=c[key].mean()))
    pd.DataFrame(rows).to_csv(f"{OD}/tremor_label_split.csv", index=False)
    log(f"[main] (c)(d) {time.time()-t0:.0f}s")
    # ---------------- (e) expected count per family-day
    from scipy.stats import spearmanr
    RF = pd.read_parquet(f"{W.RES}/step1_prob/rate_family_day.parquet")
    RF = RF[RF.truth_ok & (RF.kind == "episode")]
    rows = []
    Dd = OUT[OUT.series == "data"].copy()
    for es in EVAL:
        ge = Dd[eval_set_mask(Dd, es)]
        w = "A2011" if es == "A2011" else es.split("_")[0]
        R = RF[RF.window == w][["fid", "date", "region", "truth_own_scored", "truth_nbr3_scored", "sumP_scored_half", "n_cand_scored_half"]]
        R = R.assign(reg=R.region.map(RS))
        for key, tcol in [("D2_hgb_FULL_own", "truth_own_scored"), ("D1_hgb_FULL_own", "truth_own_scored"),
                          ("D2_hgb_FULL_nbr3", "truth_nbr3_scored"), ("D1_hgb_FULL_nbr3", "truth_nbr3_scored"),
                          ("D2_hgb_NO_PNET_own", "truth_own_scored"), ("P1", "truth_own_scored")]:
            s = ge.groupby(["fid", "date"])[key].sum().rename("sP").reset_index()
            RR = R.merge(s, on=["fid", "date"], how="left").fillna({"sP": 0})
            for rg in ["ALL"] + list(RS.values()):
                gr = RR if rg == "ALL" else RR[RR.reg == rg]
                if len(gr) < 5 or gr[tcol].sum() == 0:
                    continue
                rows.append(dict(eval_set=es, region=rg, model=key, truth=tcol, n_famdays=len(gr), sumP=gr.sP.sum(),
                                 truth_count=int(gr[tcol].sum()), ratio=gr.sP.sum() / gr[tcol].sum(),
                                 spearman=spearmanr(gr.sP, gr[tcol]).statistic,
                                 spearman_ncand=spearmanr(gr.n_cand_scored_half, gr[tcol]).statistic))
    pd.DataFrame(rows).to_csv(f"{OD}/expected_count.csv", index=False)
    log(f"[main] (e) {time.time()-t0:.0f}s")
    # ---------------- (f) permutation importance on held-out rows (D2 final model -> A2011; D1 model -> test halves + A2011)
    from sklearn.metrics import log_loss
    rows = []
    blocks = {"PNET (all)": PNW + PNF, "PNET fine-tuned": PNF, "PNET off-the-shelf": PNW, "step-1 P": ["P1"],
              "step-2 attribution": S2F}
    for key, es_list in [("D2_hgb_FULL_own", ["A2011"]), ("D2_hgb_FULL_any", ["A2011"]),
                         ("D1_hgb_FULL_own", EVAL), ("D1_hgb_FULL_any", EVAL)]:
        m = keep_models[key]
        lab = key.split("_")[-1]
        feats = FSETS["FULL"]
        msk = np.zeros(len(A), bool)
        for es in es_list:
            msk |= eval_set_mask(A, es)
        msk &= ((A.series == "data") & A.truth_ok).to_numpy()
        ii = np.nonzero(msk)[0]
        rr = np.random.default_rng(5)
        if len(ii) > 150000:
            ii = rr.choice(ii, 150000, replace=False)
        X = A[feats].to_numpy(float)[ii]; y = A[LABELS[lab]].to_numpy().astype(int)[ii]
        p0 = np.clip(m.predict_proba(X)[:, 1], 1e-6, 1 - 1e-6)
        a0 = auc_ap(y, p0)[0]; l0 = log_loss(y, p0)
        items = [(f, [f]) for f in feats] + list(blocks.items())
        for name, cols in items:
            Xp = X.copy()
            perm = rr.permutation(len(X))
            for c in cols:
                j = feats.index(c); Xp[:, j] = X[perm, j]    # permute the block jointly
            pp = np.clip(m.predict_proba(Xp)[:, 1], 1e-6, 1 - 1e-6)
            rows.append(dict(model=key, eval_rows="+".join(es_list), feature=name, is_block=name in blocks,
                             auc_drop=a0 - auc_ap(y, pp)[0], logloss_gain=log_loss(y, pp) - l0, n=len(ii)))
    pd.DataFrame(rows).to_csv(f"{OD}/feature_importance.csv", index=False)
    json.dump(dict(features=FSETS, runs=runs, labels=LABELS, held_out_stations=HO,
                   n_train=dict(D1=int(((A.series == "data") & A.truth_ok & A.window.isin(["V2008", "V2010"]) & (A.half == 0)).sum()),
                                D2=int(((A.series == "data") & A.truth_ok & A.window.isin(["V2008", "V2010"]) & (A.half == 1)).sum())),
                   model="HistGradientBoosting(max_leaf_nodes=15,min_samples_leaf=300,l2=2,lr=0.05,<=300 it, early stop) "
                         "+ isotonic CalibratedClassifierCV(GroupKFold(3) by day); logistic C=0.1 + isotonic"),
              open(f"{OD}/model_info.json", "w"), indent=1)
    log(f"[main] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
