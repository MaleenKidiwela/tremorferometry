#!/usr/bin/env python
"""WA PRODUCTION W7 + W8: per-repeat confidence P (SVI_REF WAVE standard, lfe_extension/scripts/svi_ref_fit.py pattern).

Rows      every candidate of the truth run (networks E1 / E2 / E3 / all; series data, c1, c2, rev), features from
          wa_prod_mf (+ HF from wa_prod_hf, + PNET block from wa_prod_pnet for the gate).
Labels    own (an own-family array detection within +-2 s), nbr3 (any family <= 3 km from the candidate's family), any
          (any family); defined only on days the candidate family's own array run scanned (truth exists only there).
Roles     day roles frozen in setup/day_table.csv: fit = 2007-10, 2008-05, 2009-02/04, 2010-03, 2010-07/08/09, 2011-08,
          2011-12; sealed one look = 2007-01 (headline) + 2009-05 (semi-independent); 2011-04 excluded (no array stations
          in its truth); V2008 / V2010 BUILD halves (2-h blocks) in-sample (never fitted or scored); block guards out;
          rows whose contributing template was built that day (tpl_insample) in-sample.
Models    HEADLINE WAVE = waveform features only: net_mad, net_cc, n_sta, n_sites, nch_sum, n_sta_day, nch_per_sta, wsum,
          site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, frac_z3, dom_top, hfr_rel, site_noise_mad_med, z_<E3 station>.
          Never: candidate density, family / region id, step-1 P, step-2 attribution, tremor, season, date, network id.
          FULL = WAVE + density / competition + family and region code + day-level MAD / median; NO_FAM; NO_DENSITY;
          MF_ONLY.  HistGradientBoosting + isotonic (GroupKFold by day), 5 outer folds by day -> out-of-fold P.
          Fit rows are a uniform random subsample (<= FIT_MAX rows; the base rate is kept).
Olympic   Olympic A / Arm B truth exists only in 2010: out-of-fold, and a LEAVE-REGION-OUT WAVE model (fitted without
          any Olympic row) scored on Olympic 2010 (credible if its AUC is within 0.03 of the out-of-fold AUC).
Sealed    `open_sealed` writes conf/SEALED_OPENED.json FIRST, refuses to run twice, labels 2007-01 / 2009-05 once,
          scores the frozen models, writes the bars and judges the PNET gate in the same single look.
PNET gate WITH_PNET (WAVE + PNET block) vs NO_PNET (= WAVE) fitted on the same fit rows; adopt only if dAUC(own) >= 0.01
          with a day-bootstrap 95 % CI above 0 on the sealed months and slope 0.8-1.2 (user 2026-10-03); the block made from
          the stations the fine-tune never trained on (B017, B018, GNW, HDW) is reported as a column.
Apply     every production candidate (runs prod / rich) gets every model's P (+ the WAVE fold-model spread); refit
          risk = share of rows with P_own_WAVE >= 0.75 whose 'older model' P (FULL features fitted on the V2008 + V2010
          test halves only, i.e. the truth the earlier WA confidence used) is < 0.2.

  python scripts/wa_prod_fit.py fit | open_sealed | apply --run prod | refit_risk
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.W

OC = f"{P.OUT}/conf"
NFOLD = 5
FIT_MAX = int(os.environ.get("WA_PROD_FIT_MAX", "600000"))
BINS = np.linspace(0, 1, 11)
LABELS = dict(own="m_own", nbr3="m_nbr3", any="m_any")
HO_PNET = ["PB.B017", "PB.B018", "UW.GNW", "UW.HDW"]
PNET = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_mean", "pnet_rel_max"]
PNET_HO = [c + "_ho" for c in PNET]


def feature_sets():
    E3 = P.core_sets().get("E3", [])
    STAZ = [f"z_{s.split('.')[1]}" for s in E3]
    WAVE = ["net_mad", "net_cc", "n_sta", "n_sites", "nch_sum", "n_sta_day", "nch_per_sta", "wsum", "site_cc_max",
            "site_cc_2nd", "site_z_max", "site_z_2nd", "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top",
            "hfr_rel", "site_noise_mad_med"] + STAZ
    DENS = ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]
    FAMC = ["fam_code", "reg_code"]
    DAYL = ["mad_scale", "med_scale"]
    FS = dict(WAVE=WAVE, FULL=WAVE + DENS + FAMC + DAYL, NO_FAM=WAVE + DENS + DAYL, NO_DENSITY=WAVE + FAMC + DAYL,
              MF_ONLY=["net_mad", "net_cc", "n_sites", "n_sta"])
    return FS, WAVE


def log(m):
    P.log(m)


# ------------------------------------------------------------------------------------------------ assembly
def load_run(run, dates=None, cols=None):
    parts = []
    for p in sorted(glob.glob(f"{P.OUT}/mf/{run}/det_*.parquet")):
        ds = os.path.basename(p)[4:14]
        if dates is not None and ds not in dates:
            continue
        D = pd.read_parquet(p)
        if not len(D):
            continue
        h = f"{P.OUT}/hf/{run}/hf_{ds}.parquet"
        if os.path.exists(h):
            D = D.merge(pd.read_parquet(h), on="k", how="left")
        else:
            D["hfr"] = np.nan; D["hfr_rel"] = np.nan
        q = f"{P.OUT}/pnet/{run}/pnet_{ds}.parquet"
        if os.path.exists(q):
            D = D.merge(pd.read_parquet(q), on="k", how="left")
        fd = pd.read_parquet(f"{P.OUT}/mf/{run}/famday_{ds}.parquet")
        for ser in ["data", "c1", "c2", "rev"]:
            if f"mad_{ser}" in fd:
                m = fd.set_index(["net", "fid"])[[f"mad_{ser}", f"med_{ser}"]]
                sel = D.series == ser
                if sel.any():
                    ix = pd.MultiIndex.from_arrays([D.loc[sel, "net"], D.loc[sel, "fid"]])
                    D.loc[sel, "mad_scale"] = m[f"mad_{ser}"].reindex(ix).to_numpy()
                    D.loc[sel, "med_scale"] = m[f"med_{ser}"].reindex(ix).to_numpy()
        parts.append(D if cols is None else D[[c for c in cols if c in D]])
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def context(A):
    """density / competition per (run, net, series, date) -- FULL-model features only."""
    F = P.families()
    nb = P.nbr_map()
    A = A.sort_values(["run", "net", "series", "date", "fid", "t"]).reset_index(drop=True)
    te = A.t.to_numpy(float); fam = A.fid.to_numpy().astype(int)
    out = {c: np.zeros(len(A), np.int32) for c in ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]}
    cw = lambda tq, tr, w: np.searchsorted(tr, tq + w, side="right") - np.searchsorted(tr, tq - w, side="left")
    for _, gi in A.groupby(["run", "net", "series", "date"]).indices.items():
        fg = fam[gi]
        TT = {f: np.sort(te[gi[fg == f]]) for f in np.unique(fg)}
        allt = np.sort(te[gi])
        for f in np.unique(fg):
            pos = gi[fg == f]; tq = te[pos]
            out["dens_own_1h"][pos] = cw(tq, TT[f], 3600) - 1
            out["dens_own_10m"][pos] = cw(tq, TT[f], 600) - 1
            out["dens_grp_1h"][pos] = sum(cw(tq, TT[x], 3600) for x in nb.get(f, [f]) if x in TT) - 1
            out["dens_grp_10m"][pos] = sum(cw(tq, TT[x], 600) for x in nb.get(f, [f]) if x in TT) - 1
            out["n_oth_4s"][pos] = cw(tq, allt, 4.0 - 1e-9) - cw(tq, TT[f], 4.0 - 1e-9)
    for c, v in out.items():
        A[c] = v
    return A


def add_static(A):
    F = P.families()
    A["fam_code"] = A.fid.astype(int)
    A["region"] = A.fid.map(dict(zip(F.fid, F.region)))
    A["reg"] = A.fid.map(dict(zip(F.fid, F.reg)))
    A["reg_code"] = A.reg.map({"NW": 0, "SW": 1, "OlyA": 2, "ArmB": 3}).astype(float)
    A["nch_per_sta"] = A.nch_sum / A.n_sta_day.replace(0, np.nan)
    A["month"] = A.date.str[:7]
    DT = P.day_table()
    A["role"] = A.date.map(dict(zip(DT.date, DT.role)))
    # row role: V-window build halves in-sample, guards out
    rr = A.role.astype(object).copy()
    for w, (a0, b0) in P.WINDOW_BUILD.items():
        m = (A.date >= a0) & (A.date <= b0) & (A.role == "fit")
        h = W.half_of(A.loc[m, "t"].to_numpy())
        rr.loc[m] = np.where(h == 0, "insample", np.where(h < 0, "guard", "fit"))
    if "tpl_insample" in A:
        rr[A.tpl_insample.fillna(False).astype(bool).to_numpy()] = "insample"
    A["row_role"] = rr
    return A


def label(A, sealed=False):
    F = P.families()
    T = pd.read_parquet(f"{P.OUT}/setup/sealed/truth_sealed.parquet" if sealed else f"{P.OUT}/setup/truth.parquet")
    SC = pd.read_csv(f"{P.OUT}/setup/scan_days.csv")
    scanned = set(zip(SC.det_source, SC.date))
    src = dict(zip(F.fid, F.det_source))
    nb = P.nbr_map()
    for c in LABELS.values():
        A[c] = np.nan
    lab = (A.series == "data") & np.array([(src.get(int(f)), d) in scanned for f, d in zip(A.fid, A.date)])
    if sealed:
        lab &= A.role.astype(str).str.startswith("sealed").to_numpy()
    else:
        lab &= ~A.role.astype(str).str.startswith("sealed").to_numpy()
    idx = np.flatnonzero(lab.to_numpy() if hasattr(lab, "to_numpy") else lab)
    tt = T.t.to_numpy(float); tf = T.fid.to_numpy().astype(int)
    order = np.argsort(tt); tt_s = tt[order]
    byf = {f: np.sort(tt[tf == f]) for f in np.unique(tf)}

    def near(r, q):
        if not len(r):
            return np.full(len(q), np.inf)
        j = np.clip(np.searchsorted(r, q), 1, max(len(r) - 1, 1))
        return np.minimum(np.abs(r[j - 1] - q), np.abs(r[np.minimum(j, len(r) - 1)] - q))
    t = A.t.to_numpy(float); f_ = A.fid.to_numpy().astype(int)
    m_own = np.full(len(A), np.nan); m_nb = np.full(len(A), np.nan); m_any = np.full(len(A), np.nan)
    sub = pd.DataFrame(dict(i=idx, f=f_[idx]))
    for f, g in sub.groupby("f"):
        ii = g.i.to_numpy(); q = t[ii]
        m_own[ii] = near(byf.get(f, np.zeros(0)), q) <= W.MATCH_S
        dn = np.full(len(q), np.inf)
        for x in nb.get(f, [f]):
            dn = np.minimum(dn, near(byf.get(x, np.zeros(0)), q))
        m_nb[ii] = dn <= W.MATCH_S
    m_any[idx] = near(tt_s, t[idx]) <= W.MATCH_S
    A["m_own"] = m_own; A["m_nbr3"] = m_nb; A["m_any"] = m_any
    # truth_array_frac of the family-month (carried in every truth table)
    TA = pd.read_csv(f"{P.OUT}/setup/truth_array_frac.csv")
    ta = {(int(r.fid), r.month): r.truth_array_frac for r in TA.itertuples()}
    A["truth_array_frac"] = [ta.get((int(f), m), np.nan) for f, m in zip(A.fid, A.month)]
    log(f"[label] {'SEALED ' if sealed else ''}labelled {len(idx)}; own {np.nanmean(m_own):.3f} nbr3 {np.nanmean(m_nb):.3f} any {np.nanmean(m_any):.3f}")
    return A


# ------------------------------------------------------------------------------------------------ models
def hgb(feats, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    # 2026-10-03 fix: 264 families > HGB 255-category limit -> fam_code numeric (sensitivity models only), reg_code categorical
    cat = [i for i, c in enumerate(feats) if c == "reg_code"] or None
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=200,
                                          l2_regularization=2.0, early_stopping=True, validation_fraction=0.15,
                                          n_iter_no_change=20, categorical_features=cat, random_state=seed)


def fit_cal(feats, X, y, groups, seed=0):
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.model_selection import GroupKFold
    cv = list(GroupKFold(min(3, len(np.unique(groups)))).split(X, y, groups))
    m = CalibratedClassifierCV(hgb(feats, seed), method="isotonic", cv=cv)
    m.fit(X, y)
    return m


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, float); s = np.asarray(s, float); ok = np.isfinite(y) & np.isfinite(s)
    y, s = y[ok], s[ok]
    if len(y) < 20 or y.sum() < 5 or y.sum() == len(y):
        return np.nan
    return float(roc_auc_score(y, s))


def calib(y, p):
    b = np.clip(np.digitize(p, BINS[1:-1]), 0, 9); rows = []
    for k in range(10):
        mm = b == k
        if mm.sum():
            rows.append((k, int(mm.sum()), float(p[mm].mean()), float(y[mm].mean())))
    Rr = pd.DataFrame(rows, columns=["bin", "n", "mean_P", "frac"])
    g = Rr[Rr.n >= 20]; slope = icpt = np.nan
    if len(g) >= 3 and g.mean_P.max() - g.mean_P.min() > 0.2:
        wv = g.n.to_numpy(float); xx = g.mean_P.to_numpy(); yy = g.frac.to_numpy()
        Xd = np.stack([np.ones_like(xx), xx], 1)
        icpt, slope = np.linalg.lstsq(Xd * np.sqrt(wv)[:, None], yy * np.sqrt(wv), rcond=None)[0]
    return slope, icpt, Rr


def evaluate(A, keys, tag):
    rows = []
    for (net, reg), _ in A[A.series == "data"].groupby(["net", "reg"]):
        D = A[(A.series == "data") & (A.net == net) & (A.reg == reg) & A.m_own.notna()]
        Cc = A[(A.series == "c1") & (A.net == net) & (A.reg == reg)]
        Cr = A[(A.series == "rev") & (A.net == net) & (A.reg == reg)]
        if not len(D):
            continue
        for lab, lc in LABELS.items():
            y = D[lc].to_numpy()
            base = dict(eval=tag, net=net, reg=reg, label=lab, n=len(D), n_pos=int(np.nansum(y)), n_days=D.date.nunique(),
                        auc_net_mad=auc(y, D.net_mad), med_truth_array_frac=float(np.nanmedian(D.truth_array_frac)))
            for key in [k for k in keys if k.startswith(f"P_{lab}_")]:
                p = D[key].to_numpy(float)
                ok = np.isfinite(p) & np.isfinite(y)
                if ok.sum() < 20:
                    continue
                sl, ic, _ = calib(y[ok], p[ok])
                rows.append(dict(base, model=key[len(f"P_{lab}_"):], auc=auc(y, p), slope=sl, icpt=ic, sumP=float(p[ok].sum()),
                                 matched=int(y[ok].sum()),
                                 matched_over_sumP=float(y[ok].sum() / p[ok].sum()) if p[ok].sum() > 0 else np.nan,
                                 ctl_frac_P_gt05=float(np.nanmean(Cc[key] > 0.5)) if len(Cc) and key in Cc else np.nan,
                                 rev_frac_P_gt05=float(np.nanmean(Cr[key] > 0.5)) if len(Cr) and key in Cr else np.nan,
                                 n_data_P075=int((p >= P.P_CUT).sum()),
                                 n_c1_P075=int((Cc[key] >= P.P_CUT).sum()) if len(Cc) else 0,
                                 n_rev_P075=int((Cr[key] >= P.P_CUT).sum()) if len(Cr) else 0))
    return pd.DataFrame(rows)


def by_array_frac(A, key="P_own_WAVE"):
    D = A[(A.series == "data") & A.m_own.notna() & A.truth_array_frac.notna()].copy()
    D["frac_bin"] = pd.cut(D.truth_array_frac, [-0.01, 0.45, 0.55, 0.65, 1.0], labels=["<=0.45", "0.45-0.55", "0.55-0.65", ">0.65"])
    rows = []
    for (net, reg, b), g in D.groupby(["net", "reg", "frac_bin"], observed=True):
        rows.append(dict(net=net, reg=reg, truth_array_frac=str(b), n=len(g), n_pos=int(g.m_own.sum()), auc=auc(g.m_own, g[key]),
                         auc_net_mad=auc(g.m_own, g.net_mad)))
    return pd.DataFrame(rows)


def boot(D, f1, f2, lab, nb=200, seed=3):
    from sklearn.metrics import roc_auc_score
    rs = np.random.default_rng(seed); y = D[LABELS[lab]].to_numpy().astype(int)
    p = D[f1].to_numpy(float); q = D[f2].to_numpy(float)
    dl = D.date.unique(); ix = {d: np.flatnonzero(D.date.to_numpy() == d) for d in dl}
    out = []
    for _ in range(nb):
        ii = np.concatenate([ix[d] for d in rs.choice(dl, len(dl), replace=True)])
        ok = np.isfinite(p[ii]) & np.isfinite(q[ii])
        if y[ii][ok].sum() < 5 or y[ii][ok].sum() == ok.sum():
            continue
        out.append(roc_auc_score(y[ii][ok], p[ii][ok]) - roc_auc_score(y[ii][ok], q[ii][ok]))
    return (float(np.mean(out)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (np.nan,) * 3


TRUTH_RUNS = ["truth", "truth_e12"]


def load_runs(runs, dates=None):
    parts = [load_run(r, dates=dates) for r in runs]
    parts = [x for x in parts if len(x)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def assemble_truth():
    t0 = time.time()
    A = load_runs(TRUTH_RUNS)
    A = add_static(A)
    A = context(A)
    log(f"[asm] truth rows {len(A)} {A.groupby(['net', 'series']).size().to_dict()} {time.time() - t0:.0f}s")
    return A


# ------------------------------------------------------------------------------------------------ fit
def fit(a):
    import joblib
    P.ensure(OC); t00 = time.time()
    FS, WAVE = feature_sets()
    A = assemble_truth()
    A = label(A, sealed=False)
    data = (A.series == "data").to_numpy().copy(); role = A.row_role.astype(str).to_numpy(); dates = A.date.to_numpy()
    fitm = data & (role == "fit") & A.m_own.notna().to_numpy()
    fitc = (role == "fit")
    fdays = np.array(sorted(set(dates[fitm])))
    rng = np.random.default_rng(11)
    fo = dict(zip(fdays, rng.permutation(np.arange(len(fdays)) % NFOLD)))
    fold = np.array([fo.get(d, -1) for d in dates]); fold[~fitc] = -1
    A["fold"] = fold
    # uniform subsample of fit rows (base rate kept)
    fit_idx = np.flatnonzero(fitm)
    if len(fit_idx) > FIT_MAX:
        fit_idx = np.sort(np.random.default_rng(5).choice(fit_idx, FIT_MAX, replace=False))
    fsub = np.zeros(len(A), bool); fsub[fit_idx] = True
    log(f"[fit] fit rows {int(fitm.sum())} (used {int(fsub.sum())}) on {len(fdays)} days; own pos {np.nanmean(A.m_own[fsub]):.3f}")
    MODELS, FOLDM = {}, {}
    oly = A.reg.isin(["OlyA", "ArmB"]).to_numpy().copy()
    have_pnet = "pnet_max" in A and A.loc[fitm & (A.run == "truth").to_numpy(), "pnet_max"].notna().mean() > 0.5
    pmask = A.pnet_max.notna().to_numpy() if "pnet_max" in A else np.zeros(len(A), bool)
    plan = [(fs, FS[fs], lab) for fs in FS for lab in LABELS]
    if have_pnet:
        plan += [("NO_PNET", WAVE, "own"), ("WITH_PNET", WAVE + PNET, "own"), ("WITH_PNET_HO", WAVE + PNET_HO, "own")]
    plan += [("WAVE_LRO", WAVE, lab) for lab in LABELS]           # Olympic leave-region-out
    plan += [("OLD", FS["FULL"], "own")]                           # the 'older model' for the refit-risk check
    for fs, feats, lab in plan:
        X = A[feats].to_numpy(float); y = A[LABELS[lab]].to_numpy()
        pred = np.full(len(A), np.nan)
        msk = fsub.copy()
        if fs == "WAVE_LRO":
            msk &= ~oly
        if fs == "OLD":
            msk &= A.date.between("2008-05-01", "2008-06-18").to_numpy() | A.date.between("2010-07-14", "2010-09-18").to_numpy()
        if fs in ("NO_PNET", "WITH_PNET", "WITH_PNET_HO"):
            msk &= pmask                                       # the gate compares both models on the SAME rows
        if fs not in ("WAVE_LRO", "OLD"):
            fms = []
            for k in range(NFOLD):
                tr = msk & (fold != k)
                sc = fitc & (fold == k)
                if not sc.any() or tr.sum() < 1000:
                    continue
                m = fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
                pred[sc] = m.predict_proba(X[sc])[:, 1]
                if fs == "WAVE":
                    fms.append(m)
            if fms:
                FOLDM[f"{fs}_{lab}"] = fms
        m = fit_cal(feats, X[msk], y[msk].astype(int), dates[msk])
        if fs == "WAVE_LRO":
            sc = oly.copy()
        elif fs == "OLD":
            sc = np.ones(len(A), bool)
        else:
            sc = ~fitc
        sc &= ~A.role.astype(str).str.startswith("sealed").to_numpy()
        pred[sc] = m.predict_proba(X[sc])[:, 1]
        A[f"P_{lab}_{fs}"] = pred; MODELS[f"{fs}_{lab}"] = m
        log(f"[fit] {fs} {lab}: n {int(msk.sum())} pos {int(np.nansum(y[msk]))} {time.time() - t00:.0f}s")
    joblib.dump(dict(models=MODELS, folds=FOLDM, feature_sets=FS, pnet=dict(NO_PNET=WAVE, WITH_PNET=WAVE + PNET,
                     WITH_PNET_HO=WAVE + PNET_HO)), f"{OC}/models.joblib")
    info = dict(features=FS, never=["dens_* (WAVE)", "fam_code/reg_code (WAVE)", "step-1 P", "step-2 attribution", "tremor",
                                    "season", "date", "network id"], n_fit_rows=int(fitm.sum()), n_fit_used=int(fsub.sum()),
                fit_days=len(fdays), nfold=NFOLD, labels=LABELS, have_pnet=bool(have_pnet),
                model="HGB(max_leaf_nodes=15,min_samples_leaf=200,l2=2,lr=0.05,<=300 it,early stop) + isotonic (GroupKFold(3) by day)",
                sha256_models=P.sha256(f"{OC}/models.joblib"), frozen_utc=pd.Timestamp.now("UTC").isoformat())
    json.dump(info, open(f"{OC}/MODELS_FROZEN.json", "w"), indent=1)
    keys = [k for k in A.columns if k.startswith("P_") and not k.startswith("P_own_OLD")]
    Af = A[A.row_role == "fit"]
    E = evaluate(Af, keys, "dev_oof")
    # Olympic leave-region-out vs out-of-fold on Olympic 2010
    lro = []
    for lab in LABELS:
        for reg in ("OlyA", "ArmB"):
            for net in sorted(A.net.unique()):
                D = Af[(Af.series == "data") & (Af.reg == reg) & (Af.net == net) & Af[LABELS[lab]].notna()]
                if not len(D):
                    continue
                a_oof = auc(D[LABELS[lab]], D[f"P_{lab}_WAVE"]); a_lro = auc(D[LABELS[lab]], D[f"P_{lab}_WAVE_LRO"])
                lro.append(dict(label=lab, reg=reg, net=net, n=len(D), n_pos=int(D[LABELS[lab]].sum()), auc_oof=a_oof,
                                auc_lro=a_lro, diff=a_lro - a_oof, credible=bool(np.isfinite(a_lro) and abs(a_lro - a_oof) <= 0.03)))
    pd.DataFrame(lro).to_csv(f"{OC}/olympic_lro.csv", index=False)
    E.to_csv(f"{OC}/eval_dev.csv", index=False)
    by_array_frac(Af).to_csv(f"{OC}/eval_dev_by_truth_array_frac.csv", index=False)
    ex = A[A.role == "excluded"]
    if len(ex):
        evaluate(ex, keys, "excluded_2011_04").to_csv(f"{OC}/eval_excluded_2011_04.csv", index=False)
    A.drop(columns=[c for c in A.columns if c.startswith("z_")]).to_parquet(f"{OC}/truth_rows.parquet", index=False)
    log(f"[fit] done {time.time() - t00:.0f}s; dev WAVE:\n" + E[E.model.isin(["WAVE", "FULL", "MF_ONLY"])]
        [["net", "reg", "label", "model", "n", "n_pos", "auc", "auc_net_mad", "slope", "matched_over_sumP", "ctl_frac_P_gt05",
          "rev_frac_P_gt05"]].round(3).to_string())


# ------------------------------------------------------------------------------------------------ sealed (one look)
def open_sealed(a):
    import joblib
    mk = f"{OC}/SEALED_OPENED.json"
    if os.path.exists(mk):
        log(f"[sealed] {mk} exists: the sealed months were already opened once; NOT re-scored"); return
    info = json.load(open(f"{OC}/MODELS_FROZEN.json"))
    if P.sha256(f"{OC}/models.joblib") != info["sha256_models"]:
        raise SystemExit("[sealed] models.joblib changed after freezing: refusing")
    json.dump(dict(opened_utc=pd.Timestamp.now("UTC").isoformat(), sealed=P.SEALED_MONTHS, models_sha256=info["sha256_models"],
                   rule="one look; nothing fitted or tuned on it (user 2026-10-03)"), open(mk, "w"), indent=1)
    MD = joblib.load(f"{OC}/models.joblib")
    M = MD["models"]; FS = MD["feature_sets"]; PS = MD["pnet"]
    DT = P.day_table()
    sd = set(DT[DT.role.str.startswith("sealed")].date)
    A = load_runs(TRUTH_RUNS, dates=sd)
    A = add_static(A); A = context(A)
    A = label(A, sealed=True)
    for key, m in M.items():
        fs, lab = key.rsplit("_", 1)
        feats = FS.get(fs) or PS.get(fs) or (FS["WAVE"] if fs == "WAVE_LRO" else FS["FULL"] if fs == "OLD" else None)
        if feats is None or any(c not in A for c in feats):
            continue
        A[f"P_{lab}_{fs}"] = m.predict_proba(A[feats].to_numpy(float))[:, 1]
    keys = [k for k in A.columns if k.startswith("P_")]
    rows = []
    for role in ("sealed_headline", "sealed_semi"):
        E = evaluate(A[A.role == role], keys, role)
        rows.append(E)
    E = pd.concat(rows, ignore_index=True)
    E.to_csv(f"{OC}/eval_sealed.csv", index=False)
    by_array_frac(A).to_csv(f"{OC}/eval_sealed_by_truth_array_frac.csv", index=False)
    bars = {}
    for role in ("sealed_headline", "sealed_semi"):
        for lab in LABELS:
            for net in sorted(E.net.unique()):
                r = E[(E["eval"] == role) & (E.label == lab) & (E.model == "WAVE") & (E.net == net)]
                if not len(r):
                    continue
                D = A[(A.role == role) & (A.series == "data") & (A.net == net) & A[LABELS[lab]].notna()]
                d = boot(D, f"P_{lab}_WAVE", "net_mad", lab)
                n_ = r.n.sum(); w_ = r.n / n_
                bars[f"{role}|{net}|{lab}"] = dict(auc_WAVE_by_region={k: v for k, v in zip(r.reg, r.auc)},
                                                   auc_WAVE_pooled=auc(D[LABELS[lab]], D[f"P_{lab}_WAVE"]),
                                                   auc_net_mad_pooled=auc(D[LABELS[lab]], D.net_mad), dAUC_vs_net_mad=d,
                                                   slope_by_region={k: v for k, v in zip(r.reg, r.slope)},
                                                   matched_over_sumP=float((r.matched.sum()) / r.sumP.sum()) if r.sumP.sum() > 0 else np.nan,
                                                   ctl_frac_P_gt05=float((r.ctl_frac_P_gt05 * w_).sum()),
                                                   rev_frac_P_gt05=float((r.rev_frac_P_gt05 * w_).sum()))
    json.dump(bars, open(f"{OC}/sealed_bars.json", "w"), indent=1, default=float)
    gate = dict(evaluated=False, reason="no PNET block on the sealed rows")
    if "P_own_WITH_PNET" in A and "pnet_max" in A and A.pnet_max.notna().any():
        D = A[A.role.str.startswith("sealed") & (A.series == "data") & A.m_own.notna() & A.pnet_max.notna()]
        d = boot(D, "P_own_WITH_PNET", "P_own_NO_PNET", "own")
        dho = boot(D, "P_own_WITH_PNET_HO", "P_own_NO_PNET", "own") if "P_own_WITH_PNET_HO" in D else (np.nan,) * 3
        sl, ic, _ = calib(D.m_own.to_numpy(), D.P_own_WITH_PNET.to_numpy())
        gate = dict(evaluated=True, n=len(D), auc_with=auc(D.m_own, D.P_own_WITH_PNET), auc_no=auc(D.m_own, D.P_own_NO_PNET),
                    dAUC=d, dAUC_heldout_station_block=dho, slope_with=sl,
                    rule="adopt if dAUC(own) >= 0.01 with day-bootstrap 95 % CI > 0 and slope 0.8-1.2 (user 2026-10-03); the "
                         "held-out-station block (B017, B018, GNW, HDW) is reported as a column")
        gate["adopt"] = bool(np.isfinite(d[0]) and d[0] >= 0.01 and d[1] > 0 and np.isfinite(sl) and 0.8 <= sl <= 1.2)
        gate["heldout_station_block_ge0"] = bool(np.isfinite(dho[0]) and dho[0] >= 0)      # reported, not part of the rule
    json.dump(gate, open(f"{OC}/pnet_gate.json", "w"), indent=1, default=float)
    A.drop(columns=[c for c in A.columns if c.startswith("z_")]).to_parquet(f"{OC}/sealed_rows.parquet", index=False)
    log(f"[sealed] scored once. bars: {json.dumps(bars, default=float)[:2500]}\n[sealed] PNET gate {gate}")


# ------------------------------------------------------------------------------------------------ apply
def apply(a):
    import joblib
    MD = joblib.load(f"{OC}/models.joblib")
    M = MD["models"]; FS = MD["feature_sets"]; FOLD = MD["folds"]
    gate = json.load(open(f"{OC}/pnet_gate.json")) if os.path.exists(f"{OC}/pnet_gate.json") else {}
    od = f"{OC}/P/{a.run}"
    P.ensure(od)
    dets = sorted(glob.glob(f"{P.OUT}/mf/{a.run}/det_*.parquet"))
    months = sorted({os.path.basename(p)[4:11] for p in dets})
    months = [m for j, m in enumerate(months) if j % a.nshard == a.shard]      # fixed-10 run: months split over processes
    keep = ["run", "date", "net", "series", "fid", "k", "i", "t", "net_cc", "net_mad", "n_sta", "n_sites", "n_sta_day", "hfr_rel",
            "tpl_insample", "reg", "role"]
    t0 = time.time(); ntot = 0
    for m in months:
        if os.path.exists(f"{od}/P_{m}.parquet") and not a.force:
            continue
        dd = {os.path.basename(p)[4:14] for p in dets if os.path.basename(p)[4:11] == m}
        A = load_run(a.run, dates=dd)
        if not len(A):
            continue
        A = add_static(A); A = context(A)
        out = A[[c for c in keep if c in A]].copy()
        for key, mod in M.items():
            fs, lab = key.rsplit("_", 1)
            if fs in ("NO_PNET", "WITH_PNET", "WITH_PNET_HO", "WAVE_LRO"):
                continue
            feats = FS.get(fs) or FS["FULL"]
            out[f"P_{lab}_{fs}"] = mod.predict_proba(A[feats].to_numpy(float))[:, 1].astype(np.float32)
        for lab in LABELS:
            fm = FOLD.get(f"WAVE_{lab}", [])
            if fm:
                X = A[FS["WAVE"]].to_numpy(float)
                Pf = np.stack([q.predict_proba(X)[:, 1] for q in fm])
                out[f"P_{lab}_WAVE_foldsd"] = Pf.std(0, ddof=1).astype(np.float32)
        if gate.get("adopt") and "pnet_max" in A and A.pnet_max.notna().mean() > 0.5:
            out["P_own_WAVE_PNET"] = M["WITH_PNET_own"].predict_proba(A[MD["pnet"]["WITH_PNET"]].to_numpy(float))[:, 1]
        out["P_max_WAVE"] = out[[f"P_{l}_WAVE" for l in LABELS]].max(1)
        out.to_parquet(f"{od}/P_{m}.parquet.tmp", index=False)
        os.replace(f"{od}/P_{m}.parquet.tmp", f"{od}/P_{m}.parquet")
        ntot += len(out)
        log(f"[apply {a.run}] {m}: {len(out)} rows; total {ntot}; {time.time() - t0:.0f}s")
    P.marker(f"conf/APPLY_{a.run}_DONE" + (f"_s{a.shard}" if a.nshard > 1 else ""), f"{ntot} rows")


def refit_risk(a):
    rows = []
    for p in sorted(glob.glob(f"{OC}/P/*/P_*.parquet")):
        D = pd.read_parquet(p, columns=["series", "P_own_WAVE", "P_own_OLD"])
        D = D[D.series == "data"]
        hi = D.P_own_WAVE >= P.P_CUT
        rows.append(dict(file=os.path.basename(p), n_hi=int(hi.sum()), n_hi_old_lt02=int((hi & (D.P_own_OLD < 0.2)).sum()),
                         n=len(D), n_old_hi_wave_lt02=int(((D.P_own_OLD >= P.P_CUT) & (D.P_own_WAVE < 0.2)).sum())))
    R = pd.DataFrame(rows)
    frac = R.n_hi_old_lt02.sum() / max(R.n_hi.sum(), 1)
    frac2 = R.n_old_hi_wave_lt02.sum() / max(R.n.sum(), 1)
    res = dict(n_rows_P_WAVE_ge_075=int(R.n_hi.sum()), n_of_them_old_lt_02=int(R.n_hi_old_lt02.sum()), refit_risk=float(frac),
               reverse_share_old_ge075_wave_lt02=float(frac2), floor=0.2 if frac < 0.02 else 0.1,
               rule="window floor P_max(WAVE) >= 0.2 if refit risk < 2 %, else 0.1 (user 2026-10-03)",
               older_model="FULL features fitted on the V2008 + V2010 test halves only (the truth of the earlier WA confidence)")
    json.dump(res, open(f"{OC}/refit_risk.json", "w"), indent=1)
    log(f"[refit_risk] {res}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["fit", "open_sealed", "apply", "refit_risk"])
    ap.add_argument("--run", default="prod"); ap.add_argument("--force", action="store_true")
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)   # apply only
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(fit=fit, open_sealed=open_sealed, apply=apply, refit_risk=refit_risk)[a.what](a)
