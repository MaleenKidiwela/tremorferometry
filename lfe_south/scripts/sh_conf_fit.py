#!/usr/bin/env python
"""SH_CONF stage 4 (CPU): per-repeat confidence fitted against Shelly (2025)'s catalogue (the known answer),
evaluated on independent days, applied to every forward candidate.  The WA design (lfe_extension_wa/scripts/
confidence_fit.py, 'D2'): HistGradientBoosting (shallow, regularised) + isotonic calibration (3 folds grouped by day).

Rows      every >= 8 x MAD candidate of the v3b test-day run and of the forward run (series 'data') + the site-rolled
          scrambled control of the same days (series 'c1'); features from sh_conf_mf / sh_conf_hf / sh_conf_pnet.
Context   computed here per (run, series, day): dens_own_1h/10m (other candidates of the family), dens_grp_1h/10m
          (candidates of its <= 3-km group, minus itself), n_oth_4s (other-family candidates within +-4 s);
          S1 flags: pass_9x, pass_mcc (mean cc >= 0.12), pass_floor (summed cc >= 2.75), win4s_raw (4-s cross-family
          winner among all >= 8 x MAD candidates), in_S1_9x (the Text S1 catalogue: sh_v3_common.apply_rules(k=9)).
Labels    data rows before 2024-12-01: m_same = one of his events of the SAME family within +-2 s (nearest),
          m_any = any of his events within +-2 s.  Two models (one per label).
Leakage   v3b templates are stacks of build-day events, so every build day of ANY family is in-sample: never fitted,
          never evaluated ('insample' role, P still written).  fit = v3b 2018-2023 strict test days, 5 outer folds
          grouped by day (out-of-fold P); ho2024 = v3b 2024 strict test days and fwd_ovl = forward 2023-12..2024-11
          strict days that are not fit days, both scored by the model fitted on all fit days (independent).
          Scrambled-control rows on fit days are scored by the fold model that never saw their day.
Outputs   OUT/candidates_conf.parquet (every row, every feature, P), OUT/eval_*.csv, OUT/famday_*.parquet, OUT/monthly.csv,
          OUT/forward_P_by_det_id.parquet (+ a copy results/cut/sh_ext/conf_P_by_det_id.parquet, joinable to the cut
          index on det_id).

  python scripts/sh_conf_fit.py
"""
import os, sys, json, time, glob, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_conf_common as C

MF = ["net_mad", "meancc", "ccsum", "nlive", "live_frac", "nchan", "nsta", "n_sites", "n_sites_live", "site_cc_max", "site_cc_2nd",
      "site_z_max", "site_z_2nd", "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top", "site_noise_mad_med"]
CTX = ["hfr_rel", "dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]
S1 = ["pass_9x", "pass_mcc", "pass_floor", "win4s_raw", "in_S1_9x"]
PNB = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_rel_max"]
PNW = [c + "_w" for c in PNB]; PNF = [c + "_f" for c in PNB]
FAM = ["fam_code"]
LABELS = {"same": "m_same", "any": "m_any"}
EVAL = ["fit_oof", "ho2024", "fwd_ovl"]
BINS = np.linspace(0, 1, 11)
NFOLD = 5


def log(s): C.log(s)


# ------------------------------------------------------------------------------------------------ assembly
def assemble():
    t0 = time.time()
    DT = pd.read_csv(f"{C.OUT}/day_table.csv", dtype={"date": str})
    tag = "wf" if glob.glob(f"{C.OUT}/pnet/feat_wf/*/*.parquet") else "w"
    parts = []
    for r in DT.itertuples():
        p = f"{C.OUT}/feat/{r.run}/{r.date}.parquet"
        if not os.path.exists(p): continue
        F = pd.read_parquet(p)
        if not len(F): continue
        h = f"{C.OUT}/hf/{r.run}/{r.date}.parquet"
        if os.path.exists(h):
            H = pd.read_parquet(h); F = F.merge(H.drop(columns=["fam"]), on=["series", "k"], how="left")
        q = f"{C.OUT}/pnet/feat_{tag}/{r.run}/{r.date}.parquet"
        if os.path.exists(q):
            Q = pd.read_parquet(q); F = F.merge(Q.drop(columns=["fam"]), on=["series", "k"], how="left")
        parts.append(F)
    A = pd.concat(parts, ignore_index=True)
    A = A.merge(DT[["run", "date", "role", "build_any", "year"]], on=["run", "date"], how="left")
    A["fam"] = A.fam.astype(str); A["fam_code"] = A.fam.astype(int) - 1
    A["live_frac"] = A.nlive / A.nchan.clip(lower=1)
    A["frac_z3"] = A.n_z3 / A.n_sites.clip(lower=1)
    BF = C.build_fam_days()
    A["build_own"] = [d in BF.get(f, set()) for f, d in zip(A.fam, A.date)]
    for c in PNW + PNF:
        if c not in A: A[c] = np.nan
    log(f"[asm] {len(A)} rows {A.groupby(['run', 'series']).size().to_dict()} pnet tag {tag}; {time.time() - t0:.0f}s")
    return A, tag


def count_within(tq, tref, w):
    return np.searchsorted(tref, tq + w, side="right") - np.searchsorted(tref, tq - w, side="left")


def context(A):
    from sh_v3_common import decluster, apply_rules
    t0 = time.time()
    G = C.fam_groups()
    A = A.sort_values(["run", "series", "date", "fam", "t_epoch"]).reset_index(drop=True)
    n = len(A); te = A.t_epoch.to_numpy()
    BUF = {c: np.zeros(n, np.int64) for c in ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]}
    for (run, ser, d), gi in A.groupby(["run", "series", "date"]).indices.items():
        fam = A.fam.to_numpy()[gi]
        TT = {f: te[gi[fam == f]] for f in np.unique(fam)}
        allt = np.sort(te[gi])
        for f, tq in TT.items():
            pos = gi[fam == f]
            BUF["dens_own_1h"][pos] = count_within(tq, tq, 3600) - 1
            BUF["dens_own_10m"][pos] = count_within(tq, tq, 600) - 1
            BUF["dens_grp_1h"][pos] = sum(count_within(tq, TT[x], 3600) for x in G[f] if x in TT) - 1
            BUF["dens_grp_10m"][pos] = sum(count_within(tq, TT[x], 600) for x in G[f] if x in TT) - 1
            BUF["n_oth_4s"][pos] = count_within(tq, allt, 4.0 - 1e-9) - count_within(tq, tq, 4.0 - 1e-9)
    for c, v in BUF.items():
        A[c] = v
    A["pass_9x"] = A.net_mad >= 9.0; A["pass_mcc"] = A.meancc >= 0.12; A["pass_floor"] = A.ccsum >= 2.75
    A["win4s_raw"] = False; A["in_S1_9x"] = False
    for (run, ser), g in A.groupby(["run", "series"]):
        D = g[["fam", "ccsum", "net_mad", "meancc"]].copy()
        D["t"] = pd.to_datetime(g.t_epoch, unit="s").astype("datetime64[ns]")
        A.loc[g.index, "win4s_raw"] = decluster(D.reset_index(drop=True))
        E = apply_rules(D, 9)
        A.loc[E.index, "in_S1_9x"] = True
    for c in S1: A[c] = A[c].astype(bool)
    log(f"[ctx] density / S1 flags {time.time() - t0:.0f}s; in_S1_9x {A.groupby(['run', 'series']).in_S1_9x.sum().to_dict()}")
    return A


def label(A):
    S = C.his_catalogue()
    S = S.sort_values("t_epoch")
    tall = S.t_epoch.to_numpy(); fall = S.fam.to_numpy()
    TF = {f: g.t_epoch.to_numpy() for f, g in S.groupby("fam")}
    A["m_same"] = np.nan; A["m_any"] = np.nan; A["dt_any"] = np.nan; A["his_fam_any"] = ""
    lab = (A.series == "data") & (pd.to_datetime(A.date) < C.T_HIS1)
    idx = A.index[lab]
    t = A.loc[idx, "t_epoch"].to_numpy()
    j = np.clip(np.searchsorted(tall, t), 1, len(tall) - 1)
    d0 = np.abs(tall[j - 1] - t); d1 = np.abs(tall[j] - t)
    k = np.where(d1 < d0, j, j - 1); dd = np.minimum(d0, d1)
    A.loc[idx, "m_any"] = (dd <= C.MATCH_S).astype(float); A.loc[idx, "dt_any"] = dd
    A.loc[idx, "his_fam_any"] = np.where(dd <= C.MATCH_S, fall[k], "")
    for f, gi in A.loc[idx].groupby("fam").groups.items():
        tf = TF.get(f, np.zeros(0))
        tq = A.loc[gi, "t_epoch"].to_numpy()
        if not len(tf):
            A.loc[gi, "m_same"] = 0.0; continue
        jj = np.clip(np.searchsorted(tf, tq), 1, max(len(tf) - 1, 1))
        dm = np.minimum(np.abs(tf[jj - 1] - tq), np.abs(tf[np.minimum(jj, len(tf) - 1)] - tq))
        A.loc[gi, "m_same"] = (dm <= C.MATCH_S).astype(float)
    log(f"[lab] labelled data rows {int(lab.sum())}; frac same {A.m_same.mean():.3f} any {A.m_any.mean():.3f}")
    return A


# ------------------------------------------------------------------------------------------------ models
def hgb(feats, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    cat = [feats.index("fam_code")] if "fam_code" in feats else None
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


def auc_ap(y, s):
    from sklearn.metrics import roc_auc_score, average_precision_score
    y = np.asarray(y, float); s = np.asarray(s, float); ok = np.isfinite(s) & np.isfinite(y)
    y = y[ok]; s = s[ok]
    if len(y) < 20 or y.sum() < 5 or y.sum() == len(y): return np.nan, np.nan
    return roc_auc_score(y, s), average_precision_score(y, s)


def calib(y, p):
    b = np.clip(np.digitize(p, BINS[1:-1]), 0, 9); rows = []
    for k in range(10):
        mm = b == k
        if mm.sum(): rows.append((k, int(mm.sum()), float(p[mm].mean()), float(y[mm].mean())))
    R = pd.DataFrame(rows, columns=["bin", "n", "mean_P", "frac"])
    g = R[R.n >= 20]; slope = icpt = np.nan
    if len(g) >= 3 and g.mean_P.max() - g.mean_P.min() > 0.2:
        wv = g.n.to_numpy(float); xx = g.mean_P.to_numpy(); yy = g.frac.to_numpy()
        Xd = np.stack([np.ones_like(xx), xx], 1)
        icpt, slope = np.linalg.lstsq(Xd * np.sqrt(wv)[:, None], yy * np.sqrt(wv), rcond=None)[0]
    return slope, icpt, R


def main():
    warnings.filterwarnings("ignore")
    t00 = time.time()
    A, tag = assemble()
    A = context(A)
    A = label(A)
    have_f = A.pnet_max_f.notna().any()
    FSETS = {"FULL": MF + CTX + S1 + PNW + (PNF if have_f else []) + FAM,
             "NO_PNET": MF + CTX + S1 + FAM,
             "NO_FAM": MF + CTX + S1 + PNW + (PNF if have_f else []),
             "MF_S1": MF + S1}
    if have_f:
        FSETS["PNET_W_ONLY"] = MF + CTX + S1 + PNW + FAM
    data = (A.series == "data").to_numpy()
    role = A.role.to_numpy()
    dates = A.date.to_numpy()
    fitm = data & (role == "fit")
    fdays = np.array(sorted(set(dates[fitm])))
    rng = np.random.default_rng(11)
    fo = dict(zip(fdays, rng.permutation(np.arange(len(fdays)) % NFOLD)))
    A["fold"] = [fo.get(d, -1) if r == "fit" else -1 for d, r in zip(dates, role)]
    fold = A.fold.to_numpy()
    # sensitivity design: fit days include other families' build days (own-family build days excluded)
    ownb = data & (A.run == "v3b").to_numpy() & (A.year <= 2023).to_numpy() & ~A.build_own.to_numpy()
    PRED = {}; FOLDP = {}; MODELS = {}
    runs = [(fs, lab) for fs in FSETS for lab in LABELS]
    for fs, lab in runs:
        feats = FSETS[fs]; X = A[feats].to_numpy(float); y = A[LABELS[lab]].to_numpy()
        pred = np.full(len(A), np.nan); fp = []
        for k in range(NFOLD):
            tr = fitm & (fold != k)
            m = fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
            sc = (role == "fit") & (fold == k)            # data AND control rows of the held-out fold days
            pred[sc] = m.predict_proba(X[sc])[:, 1]
            if fs in ("FULL", "NO_PNET"):
                fwd = (A.run == "fwd").to_numpy()
                q = np.full(len(A), np.nan); q[fwd] = m.predict_proba(X[fwd])[:, 1]; fp.append(q)
        m = fit_cal(feats, X[fitm], y[fitm].astype(int), dates[fitm])
        sc = role != "fit"
        pred[sc] = m.predict_proba(X[sc])[:, 1]
        PRED[f"P_{lab}_{fs}"] = pred; MODELS[f"{fs}_{lab}"] = m
        if fp: FOLDP[f"{fs}_{lab}"] = np.stack(fp)
        log(f"[fit] {fs} {lab} n_fit {int(fitm.sum())} pos {int(np.nansum(y[fitm]))} {time.time() - t00:.0f}s")
    # sensitivity: own-build-day exclusion only (FULL)
    for lab in LABELS:
        feats = FSETS["FULL"]; X = A[feats].to_numpy(float); y = A[LABELS[lab]].to_numpy()
        m = fit_cal(feats, X[ownb], y[ownb].astype(int), dates[ownb])
        sc = np.isin(role, ["ho2024", "fwd_ovl"])
        p = np.full(len(A), np.nan); p[sc] = m.predict_proba(X[sc])[:, 1]
        PRED[f"P_{lab}_FULL_ownbuild"] = p
        log(f"[fit] FULL_ownbuild {lab} n_fit {int(ownb.sum())} {time.time() - t00:.0f}s")
    for k, v in PRED.items(): A[k] = v
    A["P_same"] = A.P_same_FULL; A["P_any"] = A.P_any_FULL
    A.to_parquet(f"{C.OUT}/candidates_conf.parquet", index=False)
    json.dump(dict(features=FSETS, pnet_tag=tag, n_fit=int(fitm.sum()), fit_days=len(fdays), nfold=NFOLD,
                   model="HistGradientBoosting(max_leaf_nodes=15,min_samples_leaf=200,l2=2,lr=0.05,<=300 it,early stop)"
                         " + isotonic CalibratedClassifierCV(GroupKFold(3) by day)"),
              open(f"{C.OUT}/model_info.json", "w"), indent=1)
    import joblib
    joblib.dump(MODELS, f"{C.OUT}/models.joblib")
    log(f"[main] predictions written {time.time() - t00:.0f}s")
    evaluate(A, FSETS, MODELS)
    apply_forward(A, FOLDP)
    log(f"[main] done {time.time() - t00:.0f}s")


# ------------------------------------------------------------------------------------------------ evaluation
def eval_mask(A, es):
    r = A.role.to_numpy()
    return r == ("fit" if es == "fit_oof" else es)


def evaluate(A, FSETS, MODELS):
    from sklearn.metrics import roc_auc_score, log_loss
    t0 = time.time()
    rows, crow, rel, frows, wrows = [], [], [], [], []
    keys = [k for k in A.columns if k.startswith("P_") and k.count("_") >= 2]
    for es in EVAL:
        em = eval_mask(A, es)
        D = A[em & (A.series == "data").to_numpy()]
        Cc = A[em & (A.series == "c1").to_numpy()]
        for lab, lc in LABELS.items():
            y = D[lc].to_numpy()
            r = dict(eval_set=es, label=lab, n=len(D), n_pos=int(np.nansum(y)), base_rate=float(np.nanmean(y)) if len(D) else np.nan,
                     n_days=D.date.nunique(), n_ctl=len(Cc))
            for b in ["net_mad", "meancc", "site_z_max", "pnet_max_w", "pnet_max_f", "in_S1_9x"]:
                r[f"auc_{b}"], r[f"ap_{b}"] = auc_ap(y, D[b].to_numpy(float))
            for key in [k for k in keys if k.startswith(f"P_{lab}_")]:
                mname = key[len(f"P_{lab}_"):]
                p = D[key].to_numpy(float)
                if not np.isfinite(p).any(): continue
                r[f"auc_{mname}"], r[f"ap_{mname}"] = auc_ap(y, p)
                ok = np.isfinite(p) & np.isfinite(y)
                sl, ic, R = calib(y[ok], p[ok])
                pc = Cc[key].to_numpy(float)
                crow.append(dict(eval_set=es, label=lab, model=mname, n=int(ok.sum()), slope=sl, icpt=ic, sumP=float(p[ok].sum()),
                                 matched=int(y[ok].sum()), obs_over_exp=float(y[ok].sum() / p[ok].sum()) if p[ok].sum() > 0 else np.nan,
                                 brier=float(np.mean((p[ok] - y[ok]) ** 2)),
                                 data_fracP_gt05=float((p[ok] > 0.5).mean()), ctl_fracP_gt05=float((pc > 0.5).mean()) if len(pc) else np.nan,
                                 data_meanP=float(p[ok].mean()), ctl_meanP=float(np.nanmean(pc)) if len(pc) else np.nan,
                                 n_ctl=len(pc), ctl_sumP=float(np.nansum(pc)),
                                 ctl_over_data_sumP=float(np.nansum(pc) / p[ok].sum()) if p[ok].sum() > 0 else np.nan))
                if mname == "FULL":
                    rel.append(R.assign(eval_set=es, label=lab))
            rows.append(r)
            # within family-day ranking
            for mname in ["FULL", "net_mad"]:
                key = f"P_{lab}_FULL" if mname == "FULL" else "net_mad"
                aucs, ws = [], []
                for (f, d), g in D.groupby(["fam", "date"]):
                    yy = g[lc].to_numpy()
                    if 0 < yy.sum() < len(yy):
                        aucs.append(roc_auc_score(yy, g[key].to_numpy(float))); ws.append(len(yy))
                wrows.append(dict(eval_set=es, label=lab, score=mname, n_famdays=len(aucs),
                                  auc_within_famday=float(np.average(aucs, weights=ws)) if aucs else np.nan))
            # by family
            for f, g in D.groupby("fam"):
                yy = g[lc].to_numpy(); p = g[f"P_{lab}_FULL"].to_numpy(float)
                a1, a2 = auc_ap(yy, p) if (np.nansum(yy) >= 10 and (len(yy) - np.nansum(yy)) >= 10) else (np.nan, np.nan)
                frows.append(dict(eval_set=es, label=lab, fam=int(f), n=len(g), n_pos=int(np.nansum(yy)), auc=a1, ap=a2,
                                  auc_net_mad=auc_ap(yy, g.net_mad.to_numpy(float))[0] if np.isfinite(a1) else np.nan,
                                  sumP=float(np.nansum(p)), matched=int(np.nansum(yy)),
                                  obs_over_exp=float(np.nansum(yy) / np.nansum(p)) if np.nansum(p) > 0 else np.nan,
                                  ctl_n=int((Cc.fam == f).sum()), ctl_sumP=float(Cc[Cc.fam == f][f"P_{lab}_FULL"].sum())))
    pd.DataFrame(rows).to_csv(f"{C.OUT}/eval_auc.csv", index=False)
    pd.DataFrame(crow).to_csv(f"{C.OUT}/eval_calibration.csv", index=False)
    pd.concat(rel).to_csv(f"{C.OUT}/eval_reliability.csv", index=False)
    pd.DataFrame(frows).to_csv(f"{C.OUT}/eval_by_family.csv", index=False)
    pd.DataFrame(wrows).to_csv(f"{C.OUT}/eval_within_famday.csv", index=False)
    # day-block bootstrap of AUC(FULL) - AUC(net_mad) and AUC(FULL) - AUC(NO_PNET)
    rs = np.random.default_rng(3); brow = []
    for es in ["ho2024", "fwd_ovl"]:
        D = A[eval_mask(A, es) & (A.series == "data").to_numpy()]
        dl = D.date.unique(); ix = {d: np.nonzero(D.date.to_numpy() == d)[0] for d in dl}
        for lab, lc in LABELS.items():
            y = D[lc].to_numpy().astype(int)
            for other in ["net_mad", f"P_{lab}_NO_PNET", "in_S1_9x"]:
                p = D[f"P_{lab}_FULL"].to_numpy(); q = D[other].to_numpy(float)
                diffs = []
                for b in range(200):
                    ii = np.concatenate([ix[d] for d in rs.choice(dl, len(dl), replace=True)])
                    if y[ii].sum() < 5 or y[ii].sum() == len(ii): continue
                    diffs.append(roc_auc_score(y[ii], p[ii]) - roc_auc_score(y[ii], q[ii]))
                brow.append(dict(eval_set=es, label=lab, vs=other.replace(f"P_{lab}_", ""), n_days=len(dl),
                                 dAUC_mean=np.mean(diffs), dAUC_lo=np.percentile(diffs, 2.5), dAUC_hi=np.percentile(diffs, 97.5)))
    pd.DataFrame(brow).to_csv(f"{C.OUT}/eval_bootstrap.csv", index=False)
    # permutation importance (final FULL model; independent rows ho2024 + fwd_ovl)
    feats = FSETS["FULL"]
    blocks = {"PNET (all)": [c for c in feats if c in PNW + PNF], "PNET fine-tuned": [c for c in feats if c in PNF],
              "PNET off-the-shelf": [c for c in feats if c in PNW], "S1 flags": S1, "density / competition": [c for c in CTX if c != "hfr_rel"],
              "family id": FAM}
    msk = np.isin(A.role.to_numpy(), ["ho2024", "fwd_ovl"]) & (A.series == "data").to_numpy()
    rr = np.random.default_rng(5); prow = []
    for lab, lc in LABELS.items():
        m = MODELS[f"FULL_{lab}"]
        X = A.loc[msk, feats].to_numpy(float); y = A.loc[msk, lc].to_numpy().astype(int)
        p0 = np.clip(m.predict_proba(X)[:, 1], 1e-6, 1 - 1e-6); a0 = auc_ap(y, p0)[0]; l0 = log_loss(y, p0)
        for name, cols in [(f, [f]) for f in feats] + [(k, v) for k, v in blocks.items() if v]:
            Xp = X.copy(); perm = rr.permutation(len(X))
            for c in cols:
                j = feats.index(c); Xp[:, j] = X[perm, j]
            pp = np.clip(m.predict_proba(Xp)[:, 1], 1e-6, 1 - 1e-6)
            prow.append(dict(label=lab, feature=name, is_block=name in blocks, auc_drop=a0 - auc_ap(y, pp)[0],
                             logloss_gain=log_loss(y, pp) - l0, n=len(y)))
    pd.DataFrame(prow).to_csv(f"{C.OUT}/eval_permutation.csv", index=False)
    # single-station PNET (matched same vs unmatched), training vs held-out stations
    srow = []
    D = A[msk]
    for ns in C.STAS:
        for s in ("w", "f"):
            c = f"pn{s}_{ns}"
            if c not in D: continue
            v = D[c].to_numpy(float)
            srow.append(dict(netsta=ns, model={"w": "off-the-shelf", "f": "fine-tuned"}[s], held_out_of_finetune=ns in C.PNET_HO,
                             n=int(np.isfinite(v).sum()), auc_same=auc_ap(D.m_same.to_numpy(), v)[0], auc_any=auc_ap(D.m_any.to_numpy(), v)[0]))
    pd.DataFrame(srow).to_csv(f"{C.OUT}/eval_pnet_station.csv", index=False)
    log(f"[eval] {time.time() - t0:.0f}s")


# ------------------------------------------------------------------------------------------------ application
def apply_forward(A, FOLDP):
    t0 = time.time()
    fwd = (A.run == "fwd").to_numpy()
    rows = []
    for run in ["v3b", "fwd"]:
        R = A[A.run == run]
        for (ser, f, d), g in R.groupby(["series", "fam", "date"]):
            r = dict(run=run, series=ser, fam=f, date=d, role=g.role.iloc[0], n_cand=len(g), n_S1_9x=int(g.in_S1_9x.sum()))
            for lab in LABELS:
                p = g[f"P_{lab}"].to_numpy(float)
                r[f"sumP_{lab}"] = float(np.nansum(p)); r[f"sdB_{lab}"] = float(np.sqrt(np.nansum(p * (1 - p))))
                r[f"sumP_{lab}_NO_PNET"] = float(np.nansum(g[f"P_{lab}_NO_PNET"].to_numpy(float)))
                if run == "fwd" and f"FULL_{lab}" in FOLDP:
                    s = FOLDP[f"FULL_{lab}"][:, g.index.to_numpy()].sum(1)
                    r[f"sdM_{lab}"] = float(np.std(s, ddof=1))
                if ser == "data" and g[LABELS[lab]].notna().any():
                    r[f"matched_{lab}"] = int(np.nansum(g[LABELS[lab]].to_numpy()))
            rows.append(r)
    FD = pd.DataFrame(rows)
    # his catalogue counts per family-day
    S = C.his_catalogue(); S["date"] = S.t.dt.strftime("%Y-%m-%d")
    hc = S.groupby(["fam", "date"]).size().rename("his_n").reset_index()
    FD = FD.merge(hc, on=["fam", "date"], how="left")
    FD.loc[pd.to_datetime(FD.date) < C.T_HIS1, "his_n"] = FD.loc[pd.to_datetime(FD.date) < C.T_HIS1, "his_n"].fillna(0)
    # complete grid for data: every processed day x family (zero rows where no candidate)
    DT = pd.read_csv(f"{C.OUT}/day_table.csv", dtype={"date": str})
    DT = DT[[os.path.exists(f"{C.OUT}/feat/{r.run}/{r.date}.json") for r in DT.itertuples()]]
    grid = pd.DataFrame([(r.run, r.date, r.role, f) for r in DT.itertuples() for f in C.FAMS], columns=["run", "date", "role", "fam"])
    Dd = FD[FD.series == "data"].drop(columns=["series", "role", "his_n"])
    Cd = FD[FD.series == "c1"][["run", "fam", "date", "n_cand", "sumP_same", "sumP_any"]].rename(
        columns={"n_cand": "ctl_n", "sumP_same": "ctl_sumP_same", "sumP_any": "ctl_sumP_any"})
    G = grid.merge(Dd, on=["run", "date", "fam"], how="left").merge(Cd, on=["run", "date", "fam"], how="left").merge(hc, on=["fam", "date"], how="left")
    for c in G.columns:
        if c.startswith(("n_", "sumP", "sdB", "sdM", "ctl_", "matched")): G[c] = G[c].fillna(0)
    G.loc[pd.to_datetime(G.date) < C.T_HIS1, "his_n"] = G.loc[pd.to_datetime(G.date) < C.T_HIS1, "his_n"].fillna(0)
    G.to_parquet(f"{C.OUT}/famday.parquet", index=False)
    # monthly series (per processed day): v3b test days until 2023-11, forward from 2023-12
    G["month"] = G.date.str[:7]
    use = ((G.run == "v3b") & (G.date < "2023-12-01")) | (G.run == "fwd")
    Gm = G[use].copy()
    def agg(g):
        nd = g.date.nunique()
        out = dict(days=nd, S1_rate=g.n_S1_9x.sum() / nd, cand_rate=g.n_cand.sum() / nd,
                   his_sameday_rate=g.his_n.sum() / nd if g.his_n.notna().all() else np.nan)
        for lab in LABELS:
            out[f"sumP_{lab}_rate"] = g[f"sumP_{lab}"].sum() / nd
            v = (g[f"sdB_{lab}"] ** 2).sum() + ((g[f"sdM_{lab}"] ** 2).sum() if f"sdM_{lab}" in g else 0.0)
            out[f"sumP_{lab}_rate_sd"] = np.sqrt(v) / nd
            out[f"ctl_sumP_{lab}_rate"] = g[f"ctl_sumP_{lab}"].sum() / nd
        return pd.Series(out)
    Mo = Gm.groupby(["month", "fam"]).apply(agg).reset_index()
    Mt = Gm.groupby("month").apply(agg).reset_index(); Mt["fam"] = "all"
    Mo = pd.concat([Mo, Mt], ignore_index=True)
    S["month"] = S.t.dt.strftime("%Y-%m")
    hm = S[S.t < C.T_HIS1].groupby(["month", "fam"]).size().rename("his_n").reset_index()
    hm = pd.concat([hm, hm.groupby("month").his_n.sum().reset_index().assign(fam="all")], ignore_index=True)
    dim = {m: pd.Period(m, freq="M").days_in_month for m in sorted(set(hm.month) | set(Mo.month))}
    hm["his_rate"] = hm.his_n / hm.month.map(dim)
    Mo = Mo.merge(hm, on=["month", "fam"], how="outer").sort_values(["fam", "month"])
    Mo.to_csv(f"{C.OUT}/monthly.csv", index=False)
    # forward P keyed by the cut-window det_id
    F = A[fwd & (A.series == "data").to_numpy()][["fam", "t_epoch", "date", "seg", "role", "net_mad", "in_S1_9x", "P_same", "P_any",
                                                 "P_same_NO_PNET", "P_any_NO_PNET", "m_same", "m_any"]].copy()
    F["key"] = np.round(F.t_epoch * C.FS).astype(np.int64)
    Dc = pd.read_parquet(f"{C.RES}/cut/sh_ext/detections.parquet", columns=["det_id", "fam", "t_epoch"])
    Dc["fam"] = Dc.fam.astype(str); Dc["key"] = np.round(Dc.t_epoch * C.FS).astype(np.int64)
    J = Dc.merge(F.drop(columns=["t_epoch"]), on=["fam", "key"], how="left").drop(columns=["key"])
    J = J.rename(columns={"role": "conf_role", "m_same": "his_same_fam_2s_nearest", "m_any": "his_any_fam_2s_nearest"})
    J["in_sample"] = J.conf_role.eq("insample")
    J.to_parquet(f"{C.OUT}/forward_P_by_det_id.parquet", index=False)
    if os.path.abspath(C.OUT) == os.path.abspath(f"{C.RES}/sh_conf"):          # never from a smoke / scratch run
        J.to_parquet(f"{C.RES}/cut/sh_ext/conf_P_by_det_id.parquet", index=False)
    log(f"[apply] forward rows {len(F)}; det_id joined {int(J.P_same.notna().sum())}/{len(J)}; famday {len(G)}; {time.time() - t0:.0f}s")
    if J.P_same.isna().any():
        log(f"[apply] WARNING {int(J.P_same.isna().sum())} det_id without P")


if __name__ == "__main__":
    main()
