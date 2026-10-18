#!/usr/bin/env python
"""SVI_REF R5 (+ the PNET gate) -- per-repeat confidence P fitted against Bostock's catalogue on days with truth, with the
scrambled control on every day, a one-time sealed check on 2013, and the application to every 2014-2026 candidate.

Rows      every >= 8 x MAD candidate ('data') and every scrambled-control candidate ('c1') of the truth runs
          (2011 / 2012 / 2013 with the S15, S14 and S13 station sets; 2010 with its 9 stations) and of the production run.
Day roles frozen in OUT/day_table.csv; row roles: 2012 'fit' days are split by Bostock's alternate 2-h blocks
          (test block = fit, build block = insample: never fitted or scored, block edge = guard).
Fit       data rows with row role 'fit' (2011 inside his window + 2012 test blocks), all three station sets, 5 outer folds
          GroupKFold BY DAY (out-of-fold P for every fit-day row, control rows included); final model on all fit rows.
          HistGradientBoosting (shallow, regularised) + isotonic CalibratedClassifierCV(GroupKFold(3) by day): the SH_CONF /
          DUC_FIN design.  Labels: m_same (one of his events of the SAME family within +-2 s) and m_any (any family).
Models    WAVE (HEADLINE): network / site / station-z / HF / live-station / station-set / cleaning-flag features ONLY --
          no candidate density, no family id, no low-weight-site confirmation, no day-level MAD;
          FULL = WAVE + density / competition + family id + day-level MAD; NO_FAM; NO_DENSITY; MF_ONLY (network MF only).
          The low-weight-site confirmation (conf_z, conf_ze) and every S-energy column (ez_*) are in NO model: they are the
          independent validation instrument (svi_ref_val.py).  Never tremor, season or date.
PNET gate NO_PNET (= WAVE features) vs WITH_PNET (WAVE + PNET block), both fitted on 2011 ONLY (pnet_w05 trained on
          2008/09/10/12); judged once on 2013 inside open_sealed: adopt if dAUC(same) >= 0.01 with a day-bootstrap 95 % CI
          excluding 0 and a calibration slope 0.8-1.2.
Sealed    `open_sealed` writes OUT/conf/SEALED_OPENED.json FIRST (time, model hashes), then labels the sealed year once,
          scores every frozen model, and refuses to run again.  Nothing is fitted or tuned on it.
Outputs   OUT/conf/: models.joblib, MODELS_FROZEN.json, eval_dev_*.csv, eval_sealed_*.csv, sealed_bars.json, pnet_gate.json,
          candidates_conf.parquet (every non-sealed row + features + P), P_by_candidate.parquet (production).

  python scripts/svi_ref_fit.py fit          python scripts/svi_ref_fit.py open_sealed       python scripts/svi_ref_fit.py apply
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

OC = f"{R.OUT}/conf"
NFOLD = 5
BINS = np.linspace(0, 1, 11)
STAZ = [f"z_{s}" for s in R.S15]
WAVE = (["net_mad", "net_cc", "rel_amp", "n_sites", "n_sta", "n_sta_day", "wsh_pgc", "con_pgc", "site_cc_max", "site_cc_2nd",
         "site_z_max", "site_z_2nd", "site_z_min", "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top", "hfr", "hfr_rel", "hf_fs",
         "has_LRIV", "has_SNB", "n_sta_clbad", "n_sta_lines"] + STAZ)
DENS = ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s", "keep_xfam8"]
DAYL = ["mad_scale", "med_scale", "site_noise_mad_med"]
FAM = ["fam_code"]
MF_ONLY = ["net_mad", "net_cc", "rel_amp", "n_sites", "n_sta"]
PNET = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_mean", "pnet_rel_max"]
FSETS = dict(WAVE=WAVE, FULL=WAVE + DENS + FAM + DAYL, NO_FAM=WAVE + DENS + DAYL, NO_DENSITY=WAVE + FAM + DAYL, MF_ONLY=MF_ONLY)
PSETS = dict(NO_PNET=WAVE, WITH_PNET=WAVE + PNET)
LABELS = dict(same="m_same", any="m_any")
NEVER = ["conf_z", "conf_ze", "tremor", "quiet", "date", "season"]
FIT_YEARS = [y for y in (2011, 2012) if y != R.SEALED_YEAR]          # real run: 2011 + 2012 test blocks
PNET_FIT_YEAR = 2011 if R.SEALED_YEAR != 2011 else 2012              # real run: 2011 (smoke stand-ins only differ)


def log(m): R.log(m)


def runs_cfg():
    return json.load(open(f"{OC}/runs.json"))


def run_list(cfg, sealed=False):
    """[(run, kind, year, sset)] of the truth runs (sealed-year runs only when sealed) + production runs."""
    out = []
    for y, d in cfg["truth"].items():
        y = int(y)
        if (y == R.SEALED_YEAR) != sealed:
            continue
        for ss, run in d.items():
            out.append((run, f"truth{y}", y, ss))
    if not sealed:
        out.append((cfg["prod"], "prod", 0, "eff"))
        if cfg.get("prod_S13"):
            out.append((cfg["prod_S13"], "prod", 0, "S13"))
    return out


def load_feat(run):
    fs = sorted(p for p in glob.glob(f"{R.OUT}/feat/{run}/*.parquet") if not p.endswith("_fd.parquet") and os.path.getsize(p) > 0)
    parts = []
    for p in fs:
        try:
            F = pd.read_parquet(p)
        except Exception:
            continue
        if not len(F): continue
        q = f"{R.OUT}/pnet/feat/{run}/{os.path.basename(p)}"
        if os.path.exists(q):
            Q = pd.read_parquet(q)
            F = F.merge(Q, on=["series", "k"], how="left")
        parts.append(F)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def assemble(runs):
    t0 = time.time()
    DT = R.load_day_table()
    parts = []
    for run, kind, year, ss in runs:
        F = load_feat(run)
        if not len(F):
            log(f"[asm] {run}: no rows"); continue
        F["run"] = run; F["kind"] = kind; F["sset_run"] = ss
        parts.append(F)
        log(f"[asm] {run}: {len(F)} rows {F.series.value_counts().to_dict()}")
    A = pd.concat(parts, ignore_index=True)
    A = A.merge(DT[["date", "kind", "role"]], on=["date", "kind"], how="left")
    A["year"] = A.date.str[:4].astype(int)
    A["row_role"] = A.role.astype(object)
    m12 = (A.year == 2012) & (A.role == "fit")
    A.loc[m12, "row_role"] = R.block_role(A.loc[m12, "t"].to_numpy())
    A["fam_code"] = A.fam.astype(int)
    for c in PNET + STAZ + ["n_sta_clbad", "n_sta_lines"]:
        if c not in A: A[c] = np.nan
    for c in ("has_LRIV", "has_SNB", "keep_xfam8"):
        A[c] = A[c].astype(float)
    log(f"[asm] {len(A)} rows {time.time() - t0:.0f}s; roles {A.groupby(['kind', 'row_role']).size().to_dict()}")
    return A


def fam_groups(km=3.0):
    F = pd.read_csv(R.FAMS12)
    la, lo = F.lat.to_numpy(), F.lon.to_numpy()
    d = C.epi_km(la[:, None], lo[:, None], la[None], lo[None])
    fl = F.bkey.astype(int).tolist()
    return {f: [fl[j] for j in np.nonzero(d[i] <= km)[0]] for i, f in enumerate(fl)}


def cw(tq, tref, w):
    return np.searchsorted(tref, tq + w, side="right") - np.searchsorted(tref, tq - w, side="left")


def context(A):
    """density / competition features per (run, series, date); 2012 rows count only TEST-block neighbours."""
    t0 = time.time(); G = fam_groups()
    A = A.sort_values(["run", "series", "date", "fam", "t"]).reset_index(drop=True)
    n = len(A); te = A.t.to_numpy(float); fam = A.fam.to_numpy().astype(int)
    usable = ~((A.year == 2012) & (A.row_role != "fit")).to_numpy()
    BUF = {c: np.zeros(n, np.int64) for c in ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]}
    for key, gi in A.groupby(["run", "series", "date"]).indices.items():
        fg = fam[gi]; ug = usable[gi]
        TT = {f: np.sort(te[gi[(fg == f) & ug]]) for f in np.unique(fg)}
        allt = np.sort(te[gi[ug]])
        for f in np.unique(fg):
            pos = gi[fg == f]; tq = te[pos]; selfc = usable[pos].astype(int)
            BUF["dens_own_1h"][pos] = cw(tq, TT[f], 3600) - selfc
            BUF["dens_own_10m"][pos] = cw(tq, TT[f], 600) - selfc
            BUF["dens_grp_1h"][pos] = sum(cw(tq, TT[x], 3600) for x in G.get(f, [f]) if x in TT) - selfc
            BUF["dens_grp_10m"][pos] = sum(cw(tq, TT[x], 600) for x in G.get(f, [f]) if x in TT) - selfc
            BUF["n_oth_4s"][pos] = cw(tq, allt, 4.0 - 1e-9) - cw(tq, TT[f], 4.0 - 1e-9)
    for c, v in BUF.items(): A[c] = v
    log(f"[ctx] density {time.time() - t0:.0f}s")
    return A


def quiet_table():
    """per UTC date: PNSN tremor count in the family box and the Track C quiet rule (none that day, >= 3 days from any)."""
    F = pd.read_csv(R.FAMS12); lat0, lat1, lon0, lon1 = C.family_box(F)
    tr = pd.read_csv(C.PNSN, usecols=["time", "lat", "lon"]); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"]); tr = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1)]
    cnt = tr.time.dt.strftime("%Y-%m-%d").value_counts()
    days = pd.date_range("2010-01-01", R.PROD_T1, freq="D")
    td = pd.to_datetime(pd.Index(cnt[cnt > 0].index)).sort_values()
    Q = pd.DataFrame(dict(date=[str(d.date()) for d in days]))
    Q["n_tremor_day"] = Q.date.map(cnt).fillna(0).astype(int)
    di = days.values.astype("datetime64[D]").astype(np.int64); ti = td.values.astype("datetime64[D]").astype(np.int64)
    j = np.searchsorted(ti, di)
    dl = np.abs(di - ti[np.clip(j - 1, 0, len(ti) - 1)]); dh = np.abs(ti[np.clip(j, 0, len(ti) - 1)] - di)
    Q["quiet"] = (Q.n_tremor_day == 0) & (np.minimum(dl, dh) > 3)
    return Q


def label(A, years, allow_sealed=False):
    B = R.bostock(years, allow_sealed=allow_sealed)
    A["m_same"] = np.nan; A["m_any"] = np.nan
    lab = (A.series == "data") & A.year.isin(years) & A.role.isin(["fit", "s9", "sealed"])
    idx = A.index[lab]
    t = A.loc[idx, "t"].to_numpy(float); f = A.loc[idx, "fam"].to_numpy()
    A.loc[idx, "m_same"] = (R.nearest_dt(t, f, B.t_origin.to_numpy(), B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
    A.loc[idx, "m_any"] = (R.nearest_dt(t, f, B.t_origin.to_numpy(), B.fam.to_numpy(), False) <= R.MATCH_S).astype(float)
    log(f"[lab] years {years}: labelled {int(lab.sum())}; same {A.loc[idx, 'm_same'].mean():.3f} any {A.loc[idx, 'm_any'].mean():.3f}")
    return A, B


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


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, float); s = np.asarray(s, float); ok = np.isfinite(y) & np.isfinite(s)
    y, s = y[ok], s[ok]
    if len(y) < 20 or y.sum() < 5 or y.sum() == len(y): return np.nan
    return float(roc_auc_score(y, s))


def calib(y, p):
    b = np.clip(np.digitize(p, BINS[1:-1]), 0, 9); rows = []
    for k in range(10):
        mm = b == k
        if mm.sum(): rows.append((k, int(mm.sum()), float(p[mm].mean()), float(y[mm].mean())))
    Rr = pd.DataFrame(rows, columns=["bin", "n", "mean_P", "frac"])
    g = Rr[Rr.n >= 20]; slope = icpt = np.nan
    if len(g) >= 3 and g.mean_P.max() - g.mean_P.min() > 0.2:
        wv = g.n.to_numpy(float); xx = g.mean_P.to_numpy(); yy = g.frac.to_numpy()
        Xd = np.stack([np.ones_like(xx), xx], 1)
        icpt, slope = np.linalg.lstsq(Xd * np.sqrt(wv)[:, None], yy * np.sqrt(wv), rcond=None)[0]
    return slope, icpt, Rr


def sset_of(A):
    s = A.sset_run.astype(str).copy()
    e = s == "eff"
    s[e] = np.where(A.has_LRIV[e] > 0, "S15", np.where(A.has_SNB[e] > 0, "S14", "S13"))
    return s


def evaluate(A, keys, tag, ssets=("S15", "S14", "S13", "S9")):
    """per station set x label x model on rows A (data with labels, c1 rows for the control bar)."""
    rows, rel = [], []
    A = A.assign(sset=sset_of(A))
    for ss in ssets:
        D = A[(A.series == "data") & (A.sset == ss) & A.m_same.notna()]
        Cc = A[(A.series == "c1") & (A.sset == ss)]
        if not len(D): continue
        for lab, lc in LABELS.items():
            y = D[lc].to_numpy()
            base = dict(eval=tag, sset=ss, label=lab, n=len(D), n_pos=int(np.nansum(y)), n_days=D.date.nunique(), n_ctl=len(Cc),
                        auc_net_mad=auc(y, D.net_mad))
            for key in [k for k in keys if k.startswith(f"P_{lab}_")]:
                p = D[key].to_numpy(float)
                if not np.isfinite(p).any(): continue
                ok = np.isfinite(p) & np.isfinite(y)
                sl, ic, Rr = calib(y[ok], p[ok]); pc = Cc[key].to_numpy(float)
                rows.append(dict(base, model=key[len(f"P_{lab}_"):], auc=auc(y, p), slope=sl, icpt=ic,
                                 sumP=float(p[ok].sum()), matched=int(y[ok].sum()),
                                 matched_over_sumP=float(y[ok].sum() / p[ok].sum()) if p[ok].sum() > 0 else np.nan,
                                 ctl_frac_P_gt05=float(np.nanmean(pc > 0.5)) if len(pc) else np.nan,
                                 ctl_sumP=float(np.nansum(pc)), ctl_over_data_sumP=float(np.nansum(pc) / p[ok].sum()) if p[ok].sum() > 0 else np.nan))
                rel.append(Rr.assign(eval=tag, sset=ss, label=lab, model=key))
    return pd.DataFrame(rows), (pd.concat(rel) if rel else pd.DataFrame())


def boot(D, f1, f2, lab, nb=200, seed=3):
    """day-bootstrap of AUC(f1) - AUC(f2) on data rows D."""
    from sklearn.metrics import roc_auc_score
    rs = np.random.default_rng(seed); y = D[LABELS[lab]].to_numpy().astype(int)
    p = D[f1].to_numpy(float); q = D[f2].to_numpy(float)
    dl = D.date.unique(); ix = {d: np.flatnonzero(D.date.to_numpy() == d) for d in dl}
    out = []
    for _ in range(nb):
        ii = np.concatenate([ix[d] for d in rs.choice(dl, len(dl), replace=True)])
        if y[ii].sum() < 5 or y[ii].sum() == len(ii): continue
        ok = np.isfinite(p[ii]) & np.isfinite(q[ii])
        out.append(roc_auc_score(y[ii][ok], p[ii][ok]) - roc_auc_score(y[ii][ok], q[ii][ok]))
    return (float(np.mean(out)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (np.nan,) * 3


# ------------------------------------------------------------------------------------------------ fit
def fit(a):
    import joblib
    R.ensure(OC); t00 = time.time()
    cfg = runs_cfg()
    A = assemble(run_list(cfg, sealed=False))
    A = context(A)
    Q = quiet_table(); A = A.merge(Q, on="date", how="left")
    years = [y for y in (2010, 2011, 2012) if y != R.SEALED_YEAR]
    A, B = label(A, years)
    FY = A.year.isin(FIT_YEARS).to_numpy(); PY = (A.year == PNET_FIT_YEAR).to_numpy()
    data = (A.series == "data").to_numpy(); role = A.row_role.to_numpy(); dates = A.date.to_numpy()
    fitm = data & (role == "fit") & FY & A.m_same.notna().to_numpy()
    fitc = (role == "fit") & FY                                            # data AND control rows of fit days
    fdays = np.array(sorted(set(dates[fitm])))
    rng = np.random.default_rng(11)
    fo = dict(zip(fdays, rng.permutation(np.arange(len(fdays)) % NFOLD)))
    fold = np.array([fo.get(d, -1) for d in dates])
    fold[~fitc] = -1
    A["fold"] = fold
    MODELS, PRED, FOLDM = {}, {}, {}
    plan = [(fs, FSETS[fs], lab, fitm) for fs in FSETS for lab in LABELS]
    m11 = fitm & PY
    have_pnet = A.loc[PY, "pnet_max"].notna().mean() > 0.5 if PY.any() else False
    if have_pnet:
        plan += [(ps, PSETS[ps], "same", m11) for ps in PSETS] + [(ps, PSETS[ps], "any", m11) for ps in PSETS]
    for fs, feats, lab, msk in plan:
        X = A[feats].to_numpy(float); y = A[LABELS[lab]].to_numpy()
        pred = np.full(len(A), np.nan); fms = []
        fit_scope = fitc & (PY if fs in PSETS else True)
        for k in range(NFOLD):
            tr = msk & (fold != k)
            if not (fit_scope & (fold == k)).any(): continue
            m = fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
            sc = fit_scope & (fold == k)
            if sc.any(): pred[sc] = m.predict_proba(X[sc])[:, 1]
            if fs == "WAVE": fms.append(m)
        m = fit_cal(feats, X[msk], y[msk].astype(int), dates[msk])
        sc = ~fit_scope & (A.role != "sealed").to_numpy()
        if fs in PSETS: sc &= A.pnet_max.notna().to_numpy()
        if sc.any(): pred[sc] = m.predict_proba(X[sc])[:, 1]
        PRED[f"P_{lab}_{fs}"] = pred; MODELS[f"{fs}_{lab}"] = m
        if fms: FOLDM[f"{fs}_{lab}"] = fms
        log(f"[fit] {fs} {lab}: n_fit {int(msk.sum())} pos {int(np.nansum(y[msk]))} {time.time() - t00:.0f}s")
    for k, v in PRED.items(): A[k] = v
    A["P_same"] = A.P_same_WAVE; A["P_any"] = A.P_any_WAVE
    # fold-model spread of the headline on the production rows (for the rate uncertainty)
    prod = (A.kind == "prod").to_numpy()
    for lab in LABELS:
        fm = FOLDM.get(f"WAVE_{lab}", [])
        if fm:
            X = A.loc[prod, WAVE].to_numpy(float)
            P = np.stack([m.predict_proba(X)[:, 1] for m in fm])
            A.loc[prod, f"P_{lab}_WAVE_foldsd"] = P.std(0, ddof=1)
            np.save(f"{OC}/foldP_{lab}.npy", P.astype(np.float32))
    joblib.dump(dict(models=MODELS, folds=FOLDM), f"{OC}/models.joblib")
    info = dict(features=dict(**FSETS, **PSETS), never_features=NEVER + ["ez_*", "cc_* (per station)"], n_fit=int(fitm.sum()),
                fit_days=len(fdays), nfold=NFOLD, have_pnet_2011=bool(have_pnet), labels=LABELS,
                model="HistGradientBoosting(max_leaf_nodes=15,min_samples_leaf=200,l2=2,lr=0.05,<=300 it,early stop) + "
                      "isotonic CalibratedClassifierCV(GroupKFold(3) by day)", sealed_year=R.SEALED_YEAR,
                sha256_models=R.sha256(f"{OC}/models.joblib"), frozen_utc=pd.Timestamp.now("UTC").isoformat())
    json.dump(info, open(f"{OC}/MODELS_FROZEN.json", "w"), indent=1)
    keys = [k for k in A.columns if k.startswith("P_") and k.count("_") >= 2 and not k.endswith("foldsd")]
    E, Rl = evaluate(A[(A.row_role == "fit") & A.year.isin(FIT_YEARS)], keys, "dev_oof")
    E2, Rl2 = evaluate(A[A.role == "s9"], keys, "s9_2010", ssets=("S9",))
    pd.concat([E, E2]).to_csv(f"{OC}/eval_dev.csv", index=False); pd.concat([Rl, Rl2]).to_csv(f"{OC}/eval_dev_reliability.csv", index=False)
    A.to_parquet(f"{OC}/candidates_conf.parquet", index=False)
    log(f"[fit] done {time.time() - t00:.0f}s; dev:\n" + E[E.model.isin(["WAVE", "FULL", "MF_ONLY", "NO_PNET", "WITH_PNET"])]
        [["sset", "label", "model", "auc", "auc_net_mad", "slope", "matched_over_sumP", "ctl_frac_P_gt05"]].round(3).to_string())


# ------------------------------------------------------------------------------------------------ the one-time sealed check
def open_sealed(a):
    import joblib
    mk = f"{OC}/SEALED_OPENED.json"
    if os.path.exists(mk):
        log(f"[sealed] {mk} exists: the sealed year was already opened once; NOT re-scored"); return
    info = json.load(open(f"{OC}/MODELS_FROZEN.json"))
    if R.sha256(f"{OC}/models.joblib") != info["sha256_models"]:
        raise SystemExit("[sealed] models.joblib changed after freezing: refusing")
    json.dump(dict(opened_utc=pd.Timestamp.now("UTC").isoformat(), sealed_year=R.SEALED_YEAR, models_sha256=info["sha256_models"],
                   rule="one look; nothing fitted or tuned on it; user decision 2026-10-02 (lfe_picker's sealed 2013 referee)"),
              open(mk, "w"), indent=1)
    M = joblib.load(f"{OC}/models.joblib")["models"]
    cfg = runs_cfg()
    A = assemble(run_list(cfg, sealed=True))
    A = context(A); A = A.merge(quiet_table(), on="date", how="left")
    A, B = label(A, [R.SEALED_YEAR], allow_sealed=True)
    for key, m in M.items():
        fs, lab = key.rsplit("_", 1)
        feats = FSETS.get(fs) or PSETS.get(fs)
        if fs in PSETS and A.pnet_max.notna().mean() < 0.5: continue
        A[f"P_{lab}_{fs}"] = m.predict_proba(A[feats].to_numpy(float))[:, 1]
    keys = [k for k in A.columns if k.startswith("P_") and k.count("_") >= 2]
    sc = A[A.role == "sealed"]
    E, Rl = evaluate(sc, keys, f"sealed_{R.SEALED_YEAR}", ssets=("S15", "S14", "S13", "S9"))
    E.to_csv(f"{OC}/eval_sealed.csv", index=False); Rl.to_csv(f"{OC}/eval_sealed_reliability.csv", index=False)
    D = sc[(sc.series == "data") & sc.m_same.notna()].assign(sset=lambda x: sset_of(x))
    bars = {}
    for ss in ("S15", "S14", "S13"):
        Ds = D[D.sset == ss]
        for lab in LABELS:
            r = E[(E.sset == ss) & (E.label == lab) & (E.model == "WAVE")]
            if not len(r): continue
            r = r.iloc[0]
            d = boot(Ds, f"P_{lab}_WAVE", "net_mad", lab)
            bars[f"{ss}_{lab}"] = dict(auc_WAVE=r.auc, auc_net_mad=r.auc_net_mad, dAUC_vs_net_mad=d, slope=r.slope,
                                       matched_over_sumP=r.matched_over_sumP, ctl_frac_P_gt05=r.ctl_frac_P_gt05,
                                       pass_slope=bool(0.8 <= r.slope <= 1.2) if np.isfinite(r.slope) else False,
                                       pass_matched_over_sumP=bool(0.8 <= r.matched_over_sumP <= 1.2),
                                       pass_ctl=bool(r.ctl_frac_P_gt05 <= 0.01),
                                       wave_beats_net_mad_by_0p05=bool(np.isfinite(r.auc) and r.auc >= r.auc_net_mad + 0.05))
    json.dump(bars, open(f"{OC}/sealed_bars.json", "w"), indent=1, default=float)
    gate = dict(evaluated=False)
    if "P_same_WITH_PNET" in A:
        Dp = D[D.P_same_WITH_PNET.notna()]
        d = boot(Dp, "P_same_WITH_PNET", "P_same_NO_PNET", "same")
        y = Dp.m_same.to_numpy(); p = Dp.P_same_WITH_PNET.to_numpy()
        sl, ic, _ = calib(y, p)
        gate = dict(evaluated=True, n=len(Dp), auc_with=auc(y, p), auc_no=auc(y, Dp.P_same_NO_PNET), dAUC=d, slope_with=sl,
                    rule="adopt if dAUC(same) >= 0.01, its day-bootstrap 95 % CI excludes 0, and slope 0.8-1.2")
        gate["adopt"] = bool(np.isfinite(d[0]) and d[0] >= 0.01 and d[1] > 0 and np.isfinite(sl) and 0.8 <= sl <= 1.2)
    json.dump(gate, open(f"{OC}/pnet_gate.json", "w"), indent=1, default=float)
    sc.to_parquet(f"{OC}/sealed_rows.parquet", index=False)
    log(f"[sealed] scored once: bars {json.dumps(bars, default=float)[:1500]}\n[sealed] PNET gate {gate}")


# ------------------------------------------------------------------------------------------------ application
def apply(a):
    import joblib
    A = pd.read_parquet(f"{OC}/candidates_conf.parquet")
    M = joblib.load(f"{OC}/models.joblib")["models"]
    gate = json.load(open(f"{OC}/pnet_gate.json")) if os.path.exists(f"{OC}/pnet_gate.json") else {}
    prod = A.kind == "prod"
    if gate.get("adopt"):
        cfg = runs_cfg(); parts = []
        for q in sorted(glob.glob(f"{R.OUT}/pnet/feat/{cfg['prod']}/*.parquet")):
            Q = pd.read_parquet(q); Q["date"] = os.path.basename(q)[:10]; Q["run"] = cfg["prod"]; parts.append(Q)
        if parts:
            Q = pd.concat(parts, ignore_index=True)
            A = A.drop(columns=[c for c in Q.columns if c not in ("series", "k", "date", "run") and c in A.columns])
            A = A.merge(Q, on=["run", "date", "series", "k"], how="left")
            prod = A.kind == "prod"
            for c in PNET:
                if c not in A: A[c] = np.nan
    if gate.get("adopt") and A.loc[prod, "pnet_max"].notna().mean() > 0.5:
        for lab in LABELS:
            A.loc[prod, f"P_{lab}_WAVE_PNET"] = M[f"WITH_PNET_{lab}"].predict_proba(A.loc[prod, PSETS["WITH_PNET"]].to_numpy(float))[:, 1]
    P = A[prod].copy()
    keep = ["run", "date", "series", "k", "fam", "t", "net_mad", "sset_run", "n_sta_day", "has_LRIV", "has_SNB", "quiet",
            "n_tremor_day", "conf_z", "conf_ze", "keep_xfam8", "hfr", "n_sta_clbad"] + \
           [c for c in P.columns if c.startswith("P_")]
    P = P[[c for c in keep if c in P.columns]]
    P["month"] = P.date.str[:7]
    # join key to v1.1 (results/trackC/detections_<M>.parquet row number = det_row) on (family, origin time)
    rows = []
    for m in sorted(P.month.unique()):
        p = f"{C.TRC}/detections_{m}.parquet"
        if not os.path.exists(p): continue
        V = pd.read_parquet(p, columns=["fam", "t", "cat_core", "cat_ext"]).reset_index().rename(columns={"index": "det_row"})
        V["month"] = m; rows.append(V)
    V = pd.concat(rows, ignore_index=True)
    V["key"] = np.round(V.t * R.FS).astype(np.int64); P["key"] = np.round(P.t * R.FS).astype(np.int64)
    P = P.merge(V.drop(columns=["t"]), on=["month", "fam", "key"], how="left").drop(columns=["key"])
    P["in_v11"] = P.det_row.notna() & (P.series == "data")
    P.to_parquet(f"{OC}/P_by_candidate.parquet", index=False)
    log(f"[apply] production rows {len(P)} ({P.series.value_counts().to_dict()}); joined to v1.1 {int(P.in_v11.sum())} "
        f"of {int((P.series == 'data').sum())} data rows; PNET column {'yes' if 'P_same_WAVE_PNET' in P else 'no'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["fit", "open_sealed", "apply"])
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(fit=fit, open_sealed=open_sealed, apply=apply)[a.what](a)
