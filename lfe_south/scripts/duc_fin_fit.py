#!/usr/bin/env python
"""DUC_FIN step 4 (CPU): per-repeat confidence against Ducellier's catalogue (the SH_CONF design, scripts/sh_conf_fit.py).

Rows      every >= 8 x MAD candidate of the feature re-run 2008-01..2012-07 (series 'data') + the site-rolled scrambled
          control of the same days (series 'c1'); features from duc_fin_mf / duc_fin_hf / duc_fin_pnet (+ context here).
Features  MF: xmad (cc / hourly MAD), cc (channel mean), ccsum, nchannel, nlive, live_frac, n_tpl, n_new_tpl, nsta,
          n_sites, n_sites_live, site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, frac_z3, dom_top, site_noise_mad_med;
          CTX: hfr_rel, dens_own_1h/10m, dens_grp_1h/10m (<= 3 km family group), n_oth_4s, win_cross_2s (highest xmad of
          all families within +-2 s); RULE: pass_his, margin (ccsum / his family threshold), eliminated (threshold 0);
          PNET (_w off-the-shelf, _f fine-tuned on his events): pnet_max/2nd, n_sites_pnet03/05, pnet_rel_max;
          FAM: family id (categorical).  Per-station cc / z / PNET values are kept as columns (not fitted).
Labels    data rows before 2012-01-01: m_same = one of his CATALOGUE events (his rule: cc x nch >= threshold; the 60 kept
          families) of the same family within +-2 s; m_any = of any family within +-2 s.  Also m_same_raw / m_any_raw
          (his raw 8 x MAD detections, all 66 families; a FULL model each, secondary).
Roles     duc_fin_common.day_table: fit (5 outer folds grouped by day, out-of-fold P), ho_month (6 held-out months),
          ind2011a (Jan-Apr 2011), gap2011 (May-Dec 2011: his catalogue lost stations, D2011.md -> evaluated separately,
          never fitted), ext2012 (application), pnet_train (in-sample for PNET: never fitted / evaluated, P written).
Model     HistGradientBoosting(max_leaf_nodes=15, min_samples_leaf=200, l2=2, lr=0.05, <= 300 it, early stop) +
          isotonic CalibratedClassifierCV(GroupKFold(3) by day).
Outputs   OUT/candidates_conf.parquet (every row, every feature, P), OUT/eval_*.csv, OUT/famday.parquet, OUT/monthly.csv,
          OUT/P_by_det_id.parquet and CUT_OUT/conf_P_by_det_id.parquet (joinable to the cut index on det_id).
"""
import os, sys, json, time, glob, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
from sh_conf_fit import hgb, fit_cal, auc_ap, calib, count_within

MF_ = ["xmad", "cc", "ccsum", "nchannel", "nlive", "live_frac", "n_tpl", "n_new_tpl", "nsta", "n_sites", "n_sites_live",
       "site_cc_max", "site_cc_2nd", "site_z_max", "site_z_2nd", "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top",
       "site_noise_mad_med"]
CTX = ["hfr_rel", "dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s", "win_cross_2s"]
RULE = ["pass_his", "margin", "eliminated"]
PNB = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_rel_max"]
PNW = [c + "_w" for c in PNB]; PNF = [c + "_f" for c in PNB]
FAM = ["fam_code"]
LABELS = {"same": "m_same", "any": "m_any"}
RAWLAB = {"same_raw": "m_same_raw", "any_raw": "m_any_raw"}
EVAL = ["fit_oof", "ho_month", "ind2011a", "gap2011"]
NFOLD = 5
log = F.log


def assemble():
    t0 = time.time(); DT = F.day_table()
    tag = "wf" if glob.glob(f"{F.OUT}/pnet/feat_wf/*.parquet") else "w"
    parts = []
    for r in DT.itertuples():
        p = f"{F.OUT}/feat/{r.date}.parquet"
        if not os.path.exists(p) or os.path.getsize(p) == 0: continue
        A = pd.read_parquet(p)
        if not len(A): continue
        h = f"{F.OUT}/hf/{r.date}.parquet"
        if os.path.exists(h): A = A.merge(pd.read_parquet(h), on=["series", "k"], how="left")
        q = f"{F.OUT}/pnet/feat_{tag}/{r.date}.parquet"
        if os.path.exists(q): A = A.merge(pd.read_parquet(q), on=["series", "k"], how="left")
        A["date"] = r.date; parts.append(A)
    A = pd.concat(parts, ignore_index=True)
    A = A.merge(DT[["date", "role", "year", "month"]], on="date", how="left")
    A["fam"] = A.fam.astype(str)
    fams = sorted(DC.families().family.astype(str)); A["fam_code"] = A.fam.map({f: i for i, f in enumerate(fams)})
    A["xmad"] = A.cc / A.mad
    A["live_frac"] = A.nlive / A.nchannel.clip(lower=1); A["frac_z3"] = A.n_z3 / A.n_sites.clip(lower=1)
    A["margin"] = np.where(A.thr_perm > 0, A.ccsum / A.thr_perm.where(A.thr_perm > 0, np.nan), np.nan)
    A["eliminated"] = (A.thr_perm <= 0).astype(float); A["pass_his"] = A.pass_his.astype(float)
    for c in PNW + PNF + ["hfr_rel"]:
        if c not in A: A[c] = np.nan
    log(f"[asm] {len(A)} rows {A.groupby(['series']).size().to_dict()} pnet tag {tag}; {time.time() - t0:.0f}s")
    return A, tag


def context(A):
    t0 = time.time(); G = F.fam_groups()
    A = A.sort_values(["series", "date", "fam", "t_epoch"]).reset_index(drop=True)
    n = len(A); te = A.t_epoch.to_numpy(); xm = A.xmad.to_numpy()
    BUF = {c: np.zeros(n, np.int64) for c in ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]}
    win = np.zeros(n, bool)
    for (ser, d), gi in A.groupby(["series", "date"]).indices.items():
        fam = A.fam.to_numpy()[gi]
        TT = {f: te[gi[fam == f]] for f in np.unique(fam)}
        o = np.argsort(te[gi]); ts = te[gi][o]; xs = xm[gi][o]
        for f, tq in TT.items():
            pos = gi[fam == f]
            BUF["dens_own_1h"][pos] = count_within(tq, tq, 3600) - 1
            BUF["dens_own_10m"][pos] = count_within(tq, tq, 600) - 1
            BUF["dens_grp_1h"][pos] = sum(count_within(tq, TT[x], 3600) for x in G.get(f, [f]) if x in TT) - 1
            BUF["dens_grp_10m"][pos] = sum(count_within(tq, TT[x], 600) for x in G.get(f, [f]) if x in TT) - 1
            BUF["n_oth_4s"][pos] = count_within(tq, ts, 4.0 - 1e-9) - count_within(tq, tq, 4.0 - 1e-9)
        lo = np.searchsorted(ts, ts - 2.0, "left"); hi = np.searchsorted(ts, ts + 2.0, "right")
        w = np.array([xs[k] >= xs[lo[k]:hi[k]].max() for k in range(len(ts))])
        win[gi[o]] = w
    for c, v in BUF.items(): A[c] = v
    A["win_cross_2s"] = win.astype(float)
    log(f"[ctx] {time.time() - t0:.0f}s")
    return A


def label(A):
    H = F.his_catalogue()
    lab = (A.series == "data") & (pd.to_datetime(A.date) < pd.Timestamp("2012-01-01"))
    idx = A.index[lab]; t = A.loc[idx, "t_epoch"].to_numpy(); f = A.loc[idx, "fam"].to_numpy()
    for nm, Hs in (("", H[H.hp]), ("_raw", H)):
        A[f"m_same{nm}"] = np.nan; A[f"m_any{nm}"] = np.nan
        A.loc[idx, f"m_any{nm}"] = (F.nearest(t, Hs.t_epoch.to_numpy()) <= F.MATCH_S).astype(float)
        A.loc[idx, f"m_same{nm}"] = (F.nearest(t, Hs.t_epoch.to_numpy(), f, Hs.fam.to_numpy()) <= F.MATCH_S).astype(float)
    log(f"[lab] labelled data rows {int(lab.sum())}; frac same {A.m_same.mean():.3f} any {A.m_any.mean():.3f} "
        f"same_raw {A.m_same_raw.mean():.3f} any_raw {A.m_any_raw.mean():.3f}")
    return A


def main():
    warnings.filterwarnings("ignore"); t00 = time.time()
    A, tag = assemble(); A = context(A); A = label(A)
    have_f = A.pnet_max_f.notna().any()
    FSETS = {"FULL": MF_ + CTX + RULE + PNW + (PNF if have_f else []) + FAM, "NO_PNET": MF_ + CTX + RULE + FAM,
             "NO_FAM": MF_ + CTX + RULE + PNW + (PNF if have_f else []), "MF_RULE": MF_ + RULE}
    if have_f: FSETS["PNET_W_ONLY"] = MF_ + CTX + RULE + PNW + FAM
    data = (A.series == "data").to_numpy(); role = A.role.to_numpy(); dates = A.date.to_numpy()
    fitm = data & (role == "fit")
    fdays = np.array(sorted(set(dates[fitm]))); rng = np.random.default_rng(11)
    fo = dict(zip(fdays, rng.permutation(np.arange(len(fdays)) % NFOLD)))
    A["fold"] = [fo.get(d, -1) if r == "fit" else -1 for d, r in zip(dates, role)]; fold = A.fold.to_numpy()
    PRED = {}; MODELS = {}
    runs = [(fs, lab, lc) for fs in FSETS for lab, lc in LABELS.items()] + [("FULL", lab, lc) for lab, lc in RAWLAB.items()]
    for fs, lab, lc in runs:
        feats = FSETS[fs]; X = A[feats].to_numpy(float); y = A[lc].to_numpy()
        pred = np.full(len(A), np.nan)
        for k in range(NFOLD):
            tr = fitm & (fold != k); sc = (role == "fit") & (fold == k)
            if not sc.any() or len(np.unique(dates[tr])) < 2: continue
            m = fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
            pred[sc] = m.predict_proba(X[sc])[:, 1]
        m = fit_cal(feats, X[fitm], y[fitm].astype(int), dates[fitm])
        sc = role != "fit"
        if sc.any(): pred[sc] = m.predict_proba(X[sc])[:, 1]
        PRED[f"P_{lab}_{fs}"] = pred; MODELS[f"{fs}_{lab}"] = m
        log(f"[fit] {fs} {lab} n_fit {int(fitm.sum())} pos {int(np.nansum(y[fitm]))} {time.time() - t00:.0f}s")
    for k, v in PRED.items(): A[k] = v
    A["P_same"] = A.P_same_FULL; A["P_any"] = A.P_any_FULL
    A["in_sample_pnet"] = A.role.eq("pnet_train")
    A.to_parquet(f"{F.OUT}/candidates_conf.parquet", index=False)
    json.dump(dict(features=FSETS, pnet_tag=tag, n_fit=int(fitm.sum()), fit_days=len(fdays), nfold=NFOLD,
                   model="HistGradientBoosting(max_leaf_nodes=15,min_samples_leaf=200,l2=2,lr=0.05,<=300 it,early stop) + "
                         "isotonic CalibratedClassifierCV(GroupKFold(3) by day)", ho_months=F.HO_MONTHS, pnet_held_out=F.PNET_HO),
              open(f"{F.OUT}/model_info.json", "w"), indent=1)
    import joblib; joblib.dump(MODELS, f"{F.OUT}/models.joblib")
    evaluate(A, FSETS, MODELS); apply_all(A)
    log(f"[main] done {time.time() - t00:.0f}s")


def evaluate(A, FSETS, MODELS):
    from sklearn.metrics import roc_auc_score, log_loss
    t0 = time.time(); rows, crow, rel, frows, wrows = [], [], [], [], []
    keys = [k for k in A.columns if k.startswith("P_") and k.count("_") >= 2]
    for es in EVAL:
        em = A.role.to_numpy() == ("fit" if es == "fit_oof" else es)
        D = A[em & (A.series == "data").to_numpy()]; Cc = A[em & (A.series == "c1").to_numpy()]
        for lab, lc in {**LABELS, **RAWLAB}.items():
            y = D[lc].to_numpy()
            r = dict(eval_set=es, label=lab, n=len(D), n_pos=int(np.nansum(y)), base_rate=float(np.nanmean(y)) if len(D) else np.nan,
                     n_days=D.date.nunique(), n_ctl=len(Cc))
            for b in ["xmad", "cc", "ccsum", "margin", "pass_his", "site_z_max", "pnet_max_w", "pnet_max_f"]:
                r[f"auc_{b}"], r[f"ap_{b}"] = auc_ap(y, D[b].to_numpy(float))
            for key in [k for k in keys if k.startswith(f"P_{lab}_")]:
                mname = key[len(f"P_{lab}_"):]
                if lab in LABELS and mname.startswith("raw"): continue
                p = D[key].to_numpy(float)
                if not np.isfinite(p).any(): continue
                r[f"auc_{mname}"], r[f"ap_{mname}"] = auc_ap(y, p)
                ok = np.isfinite(p) & np.isfinite(y)
                if ok.sum() < 20: continue
                sl, ic, R = calib(y[ok], p[ok]); pc = Cc[key].to_numpy(float)
                crow.append(dict(eval_set=es, label=lab, model=mname, n=int(ok.sum()), slope=sl, icpt=ic, sumP=float(p[ok].sum()),
                                 matched=int(y[ok].sum()), obs_over_exp=float(y[ok].sum() / p[ok].sum()) if p[ok].sum() > 0 else np.nan,
                                 brier=float(np.mean((p[ok] - y[ok]) ** 2)), data_fracP_gt05=float((p[ok] > 0.5).mean()),
                                 ctl_fracP_gt05=float((pc > 0.5).mean()) if len(pc) else np.nan, data_meanP=float(p[ok].mean()),
                                 ctl_meanP=float(np.nanmean(pc)) if len(pc) else np.nan, n_ctl=len(pc), ctl_sumP=float(np.nansum(pc)),
                                 ctl_over_data_sumP=float(np.nansum(pc) / p[ok].sum()) if p[ok].sum() > 0 else np.nan))
                if mname == "FULL": rel.append(R.assign(eval_set=es, label=lab))
            rows.append(r)
            if lab not in LABELS: continue
            for mname in ["FULL", "xmad"]:
                key = f"P_{lab}_FULL" if mname == "FULL" else "xmad"; aucs, ws = [], []
                for (f, d), g in D.groupby(["fam", "date"]):
                    yy = g[lc].to_numpy()
                    if 0 < yy.sum() < len(yy): aucs.append(roc_auc_score(yy, g[key].to_numpy(float))); ws.append(len(yy))
                wrows.append(dict(eval_set=es, label=lab, score=mname, n_famdays=len(aucs),
                                  auc_within_famday=float(np.average(aucs, weights=ws)) if aucs else np.nan))
            for f, g in D.groupby("fam"):
                yy = g[lc].to_numpy(); p = g[f"P_{lab}_FULL"].to_numpy(float)
                a1 = auc_ap(yy, p)[0] if (np.nansum(yy) >= 10 and (len(yy) - np.nansum(yy)) >= 10) else np.nan
                frows.append(dict(eval_set=es, label=lab, fam=f, n=len(g), n_pos=int(np.nansum(yy)), auc=a1, sumP=float(np.nansum(p)),
                                  matched=int(np.nansum(yy)), obs_over_exp=float(np.nansum(yy) / np.nansum(p)) if np.nansum(p) > 0 else np.nan,
                                  ctl_n=int((Cc.fam == f).sum()), ctl_sumP=float(Cc[Cc.fam == f][f"P_{lab}_FULL"].sum())))
    pd.DataFrame(rows).to_csv(f"{F.OUT}/eval_auc.csv", index=False)
    pd.DataFrame(crow).to_csv(f"{F.OUT}/eval_calibration.csv", index=False)
    if rel: pd.concat(rel).to_csv(f"{F.OUT}/eval_reliability.csv", index=False)
    pd.DataFrame(frows).to_csv(f"{F.OUT}/eval_by_family.csv", index=False)
    pd.DataFrame(wrows).to_csv(f"{F.OUT}/eval_within_famday.csv", index=False)
    # day-block bootstrap on the independent sets
    rs = np.random.default_rng(3); brow = []
    for es in ["ho_month", "ind2011a"]:
        D = A[(A.role == es) & (A.series == "data")]
        dl = D.date.unique(); ix = {d: np.nonzero(D.date.to_numpy() == d)[0] for d in dl}
        for lab, lc in LABELS.items():
            y = D[lc].to_numpy().astype(int)
            for other in ["xmad", "pass_his", f"P_{lab}_NO_PNET", f"P_{lab}_MF_RULE"]:
                p = D[f"P_{lab}_FULL"].to_numpy(); q = D[other].to_numpy(float); diffs = []
                for b in range(200):
                    ii = np.concatenate([ix[d] for d in rs.choice(dl, len(dl), replace=True)])
                    if y[ii].sum() < 5 or y[ii].sum() == len(ii): continue
                    diffs.append(roc_auc_score(y[ii], p[ii]) - roc_auc_score(y[ii], q[ii]))
                if diffs:
                    brow.append(dict(eval_set=es, label=lab, vs=other.replace(f"P_{lab}_", ""), n_days=len(dl), dAUC_mean=np.mean(diffs),
                                     dAUC_lo=np.percentile(diffs, 2.5), dAUC_hi=np.percentile(diffs, 97.5)))
    pd.DataFrame(brow).to_csv(f"{F.OUT}/eval_bootstrap.csv", index=False)
    # permutation importance (final FULL model; independent rows)
    feats = FSETS["FULL"]
    blocks = {"PNET (all)": [c for c in feats if c in PNW + PNF], "PNET fine-tuned": [c for c in feats if c in PNF],
              "PNET off-the-shelf": [c for c in feats if c in PNW], "his-rule flags": RULE,
              "density / competition": [c for c in CTX if c != "hfr_rel"], "family id": FAM,
              "site cc / z": [c for c in feats if c.startswith(("site_", "n_z", "frac_z", "dom_"))]}
    msk = np.isin(A.role.to_numpy(), ["ho_month", "ind2011a"]) & (A.series == "data").to_numpy()
    rr = np.random.default_rng(5); prow = []
    for lab, lc in LABELS.items():
        m = MODELS[f"FULL_{lab}"]; X = A.loc[msk, feats].to_numpy(float); y = A.loc[msk, lc].to_numpy().astype(int)
        p0 = np.clip(m.predict_proba(X)[:, 1], 1e-6, 1 - 1e-6); a0 = auc_ap(y, p0)[0]; l0 = log_loss(y, p0)
        for name, cols in [(f, [f]) for f in feats] + [(k, v) for k, v in blocks.items() if v]:
            Xp = X.copy(); perm = rr.permutation(len(X))
            for c in cols: j = feats.index(c); Xp[:, j] = X[perm, j]
            pp = np.clip(m.predict_proba(Xp)[:, 1], 1e-6, 1 - 1e-6)
            prow.append(dict(label=lab, feature=name, is_block=name in blocks, auc_drop=a0 - auc_ap(y, pp)[0],
                             logloss_gain=log_loss(y, pp) - l0, n=len(y)))
    pd.DataFrame(prow).to_csv(f"{F.OUT}/eval_permutation.csv", index=False)
    # single-station PNET: off-the-shelf (before) vs fine-tuned (after), training vs held-out stations
    srow = []; D = A[msk]
    for ns in F.stations():
        for s in ("w", "f"):
            c = f"pn{s}_{ns}"
            if c not in D: continue
            v = D[c].to_numpy(float)
            srow.append(dict(netsta=ns, model={"w": "off-the-shelf", "f": "fine-tuned"}[s], held_out_of_finetune=ns in F.PNET_HO,
                             n=int(np.isfinite(v).sum()), auc_same=auc_ap(D.m_same.to_numpy(), v)[0], auc_any=auc_ap(D.m_any.to_numpy(), v)[0]))
    pd.DataFrame(srow).to_csv(f"{F.OUT}/eval_pnet_station.csv", index=False)
    log(f"[eval] {time.time() - t0:.0f}s")


def apply_all(A):
    t0 = time.time()
    rows = []
    for (ser, f, d), g in A.groupby(["series", "fam", "date"]):
        r = dict(series=ser, fam=f, date=d, role=g.role.iloc[0], n_cand=len(g), n_pass_his=int(g.pass_his.sum()))
        for lab in LABELS:
            p = g[f"P_{lab}"].to_numpy(float)
            r[f"sumP_{lab}"] = float(np.nansum(p)); r[f"sdB_{lab}"] = float(np.sqrt(np.nansum(p * (1 - p))))
            r[f"sumP_{lab}_NO_PNET"] = float(np.nansum(g[f"P_{lab}_NO_PNET"].to_numpy(float)))
            if ser == "data" and g[LABELS[lab]].notna().any(): r[f"matched_{lab}"] = int(np.nansum(g[LABELS[lab]].to_numpy()))
        rows.append(r)
    FD = pd.DataFrame(rows)
    H = F.his_catalogue(); H["date"] = H.t.dt.strftime("%Y-%m-%d")
    hc = H[H.hp].groupby(["fam", "date"]).size().rename("his_n").reset_index()
    Dd = FD[FD.series == "data"].drop(columns=["series"])
    Cd = FD[FD.series == "c1"][["fam", "date", "n_cand", "sumP_same", "sumP_any"]].rename(
        columns={"n_cand": "ctl_n", "sumP_same": "ctl_sumP_same", "sumP_any": "ctl_sumP_any"})
    DT = F.day_table(); DT = DT[[os.path.exists(f"{F.OUT}/feat/{d}.json") for d in DT.date]]
    fams = sorted(DC.families().family.astype(str))
    grid = pd.DataFrame([(r.date, r.role, f) for r in DT.itertuples() for f in fams], columns=["date", "role", "fam"])
    G = grid.merge(Dd.drop(columns=["role"]), on=["date", "fam"], how="left").merge(Cd, on=["date", "fam"], how="left").merge(hc, on=["fam", "date"], how="left")
    for c in G.columns:
        if c.startswith(("n_", "sumP", "sdB", "ctl_", "matched", "his_n")): G[c] = G[c].fillna(0)
    G.loc[G.date >= "2012-01-01", "his_n"] = np.nan
    G.to_parquet(f"{F.OUT}/famday.parquet", index=False)
    G["month"] = G.date.str[:7]
    def agg(g):
        nd = g.date.nunique()
        o = dict(days=nd, pass_his_rate=g.n_pass_his.sum() / nd, cand_rate=g.n_cand.sum() / nd, his_rate=g.his_n.sum() / nd if g.his_n.notna().all() else np.nan)
        for lab in LABELS:
            o[f"sumP_{lab}_rate"] = g[f"sumP_{lab}"].sum() / nd; o[f"sumP_{lab}_rate_sd"] = np.sqrt((g[f"sdB_{lab}"] ** 2).sum()) / nd
            o[f"ctl_sumP_{lab}_rate"] = g[f"ctl_sumP_{lab}"].sum() / nd
        return pd.Series(o)
    Mo = G.groupby(["month", "fam"]).apply(agg).reset_index()
    Mt = G.groupby("month").apply(agg).reset_index(); Mt["fam"] = "all"
    pd.concat([Mo, Mt], ignore_index=True).to_csv(f"{F.OUT}/monthly.csv", index=False)
    # P keyed by the cut-window det_id
    Dt = A[A.series == "data"][["fam", "t_epoch", "date", "role", "xmad", "pass_his", "P_same", "P_any", "P_same_NO_PNET", "P_any_NO_PNET",
                                 "P_same_raw_FULL", "P_any_raw_FULL", "m_same", "m_any", "m_same_raw", "m_any_raw", "in_sample_pnet"]].copy()
    Dt["key"] = np.round(Dt.t_epoch * DC.FS).astype(np.int64)
    Dt.drop(columns=["key"]).to_parquet(f"{F.OUT}/P_by_candidate.parquet", index=False)
    pc = f"{F.CUT_OUT}/detections.parquet"
    if os.path.exists(pc):
        Dc = pd.read_parquet(pc, columns=["det_id", "fam", "t_epoch"]); Dc["fam"] = Dc.fam.astype(str)
        Dc["key"] = np.round(Dc.t_epoch * DC.FS).astype(np.int64)
        J = Dc.merge(Dt.drop(columns=["t_epoch"]), on=["fam", "key"], how="left").drop(columns=["key"])
        J = J.rename(columns={"role": "conf_role", "m_same": "his_same_fam_2s", "m_any": "his_any_fam_2s"})
        J.to_parquet(f"{F.OUT}/P_by_det_id.parquet", index=False); J.to_parquet(f"{F.CUT_OUT}/conf_P_by_det_id.parquet", index=False)
        log(f"[apply] det_id joined {int(J.P_same.notna().sum())}/{len(J)} -> {F.CUT_OUT}/conf_P_by_det_id.parquet")
        if J.P_same.isna().mean() > F.ERR_MAX:
            log("[apply] FAIL: > 3 % of the cut detections have no P"); sys.exit(8)
    else:
        log(f"[apply] WARNING no {pc}; P not joined to det_id")
    log(f"[apply] {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
