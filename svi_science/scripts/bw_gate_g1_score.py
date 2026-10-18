#!/usr/bin/env python
"""bw_gate_g1_score -- G1: sensitivity of the reduced backward networks (CN4 / B1 / B2) on the forward truth days.

Scoring = lfe_extension/scripts/svi_ref_r3.py (same functions re-implemented on explicit run paths): one-to-one recall
of Bostock's events (own family, +-2 s) by every candidate >= T, with T of the reduced set chosen so that its
scrambled-control rate (control candidates per observable family-day, scored days) MATCHES the full set's rate at the
forward operating thresholds 8.25 / 8.75 MAD.  Scored rows = role 'fit' (2011 in his window; 2012 TEST blocks only) and
role 's9' (2010).  2013 is never read.
Plus (reported, SCOPE model bar): an out-of-fold WAVE confidence model per station set, exactly the forward recipe
(svi_ref_fit.fit_cal: HistGradientBoosting + isotonic, 5 folds grouped by day, 2011 + 2012 test blocks), with the full
set (S15) refitted by the same code as the reference; CN4's final model is then scored once on 2010 (out of sample).

Out: OUT/tables/g1_recall.csv, g1_families.csv, g1_observable.csv, g1_model.csv, g1_decision.json; figure g1_recall.png
"""
from __future__ import annotations
import os, sys, glob, json, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
R, C = G.patch_forward()
warnings.filterwarnings("ignore")

TAB = f"{G.OUT}/tables"
REF = {2010: ("S9", f"{G.REFINE}/feat/r3t2010_nc"), 2011: ("S15", f"{G.REFINE}/feat/r3t2011_nc"),
       2012: ("S15", f"{G.REFINE}/feat/r3t2012_nc")}
NEW = lambda y, s: f"{G.OUT}/mf/feat/g1t{y}_{s}"
SETS = {2010: ["CN4"], 2011: ["CN4", "B2", "B1"], 2012: ["CN4", "B2", "B1"]}
COLS = None


def load(path, cols=None):
    import pyarrow.parquet as pq
    fs = sorted(p for p in glob.glob(f"{path}/*.parquet") if not p.endswith("_fd.parquet"))
    fs = [p for p in fs if os.path.getsize(p) > 0 and pq.read_metadata(p).num_rows > 0]
    A = pd.concat([pd.read_parquet(p, columns=cols) for p in fs], ignore_index=True)
    fd = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{path}/*_fd.parquet")) if os.path.getsize(p) > 0],
                   ignore_index=True)
    return A, fd


DT = pd.read_csv(f"{G.REFINE}/day_table.csv", dtype={"date": str})


def scored(A, year):
    roles = DT[DT.kind == f"truth{year}"].set_index("date").role.to_dict()
    A = A[A.date.map(roles).isin(["fit", "s9"])].copy()
    if year == 2012:
        A = A[R.block_role(A.t.to_numpy()) == "fit"]
    return A


def scored_days(year):
    return sorted(DT[(DT.kind == f"truth{year}") & DT.role.isin(["fit", "s9"])].date)


def ctl_rate(A, fd, T, days):
    n = int(((A.series == "c1") & (A.net_mad >= T) & A.date.isin(days)).sum())
    m = int(((fd.series == "c1") & fd.observable.astype(bool) & fd.date.isin(days)).sum())
    return n / max(m, 1)


def matched_T(A, fd, target, days):
    grid = np.round(np.arange(8.0, 14.0, 0.01), 2)
    c = A[(A.series == "c1") & A.date.isin(days)].net_mad.to_numpy()
    m = int(((fd.series == "c1") & fd.observable.astype(bool) & fd.date.isin(days)).sum())
    for T in grid:
        if (c >= T).sum() / max(m, 1) <= target + 1e-12:
            return float(T)
    return float(grid[-1])


def truth(year, days):
    B = R.bostock([year]); B = B[B.date.isin(days)]
    if year == 2012:
        B = B[R.block_role(B.t_origin.to_numpy()) == "fit"]
    return B


def score_set(A, fd, year, T, B):
    D = A[(A.series == "data") & (A.net_mad >= T)]
    hit = R.one2one_recall(D.t.to_numpy(), D.fam.to_numpy(), B.t_origin.to_numpy(), B.fam.to_numpy())
    dt = R.nearest_dt(D.t.to_numpy(), D.fam.to_numpy(), B.t_origin.to_numpy(), B.fam.to_numpy(), True)
    prec = float((dt <= R.MATCH_S).mean()) if len(D) else np.nan
    fam = B.fam.to_numpy(); fn = pd.Series(hit).groupby(fam).sum()
    return hit, prec, len(D), fn


def recall_tables():
    rows, frows, orows = [], [], []
    for year in (2010, 2011, 2012):
        rs, rp = REF[year]
        Ar, fdr = load(rp, ["series", "fam", "date", "t", "net_mad"]); Ar = scored(Ar, year)
        days = scored_days(year); B = truth(year, days)
        for s in SETS[year]:
            if not os.path.exists(NEW(year, s)):
                G.log(f"[g1] missing run {NEW(year, s)}"); continue
            An, fdn = load(NEW(year, s), ["series", "fam", "date", "t", "net_mad"]); An = scored(An, year)
            # observability: fraction of scored days each family has >= 2 live sites (data series)
            for lab, fdx in ((s, fdn), (rs, fdr)):
                o = fdx[(fdx.series == "data") & fdx.date.isin(days)].groupby("fam").observable.mean()
                for f_ in range(130):
                    orows.append(dict(year=year, sset=lab, fam=f_, frac_days_observable=float(o.get(f_, 0.0))))
            for Tref in (8.25, 8.75):
                target = ctl_rate(Ar, fdr, Tref, days)
                Tm = matched_T(An, fdn, target, days)
                hr, pr, ndr, fr = score_set(Ar, fdr, year, Tref, B)
                hn, pn, ndn, fn_ = score_set(An, fdn, year, Tm, B)
                rows.append(dict(year=year, sset=s, ref=rs, T_ref=Tref, ctl_rate_target=target, T_matched=Tm,
                                 ctl_rate_new=ctl_rate(An, fdn, Tm, days), n_events=len(B),
                                 recall_ref=hr.mean(), recall_new=hn.mean(), recall_ratio=hn.mean() / max(hr.mean(), 1e-9),
                                 matched_ref=int(hr.sum()), matched_new=int(hn.sum()),
                                 precision_ref=pr, precision_new=pn, n_data_ref=ndr, n_data_new=ndn,
                                 n_ctl_ref=int(((Ar.series == "c1") & (Ar.net_mad >= Tref)).sum()),
                                 n_ctl_new=int(((An.series == "c1") & (An.net_mad >= Tm)).sum()),
                                 kept_of_ref=float((hr & hn).sum() / max(hr.sum(), 1)),
                                 n_fam_red_ref=int((fr >= 5).sum()), n_fam_red_new=int((fn_ >= 5).sum()),
                                 n_days=len(days)))
                for f_ in sorted(set(fr.index) | set(fn_.index)):
                    frows.append(dict(year=year, sset=s, T_ref=Tref, fam=int(f_), n_events=int((B.fam == f_).sum()),
                                      matched_ref=int(fr.get(f_, 0)), matched_new=int(fn_.get(f_, 0))))
                G.log(f"[g1] {year} {s} T{Tref}: recall {hr.mean():.4f} -> {hn.mean():.4f} (x{rows[-1]['recall_ratio']:.2f}) "
                      f"at T {Tm}; precision {pr:.3f} -> {pn:.3f}; fam red {rows[-1]['n_fam_red_ref']} -> {rows[-1]['n_fam_red_new']}")
    Rr = pd.DataFrame(rows); Fr = pd.DataFrame(frows); Or = pd.DataFrame(orows).drop_duplicates(["year", "sset", "fam"])
    Rr.to_csv(f"{TAB}/g1_recall.csv", index=False); Fr.to_csv(f"{TAB}/g1_families.csv", index=False)
    Or.to_csv(f"{TAB}/g1_observable.csv", index=False)
    return Rr, Fr, Or


# ------------------------------------------------------------------------------------------- OOF confidence model
def model_tables():
    import svi_ref_fit as F
    out = []
    parts = {}
    for s in ["S15", "CN4", "B2", "B1"]:
        L = []
        for year in (2011, 2012):
            p = REF[year][1] if s == "S15" else NEW(year, s)
            if not os.path.exists(p): continue
            A, _ = load(p); A = scored(A, year); A["year"] = year
            B = R.bostock([year])
            d = A.series == "data"
            A["m_same"] = np.nan
            A.loc[d, "m_same"] = (R.nearest_dt(A.loc[d, "t"].to_numpy(float), A.loc[d, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                               B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
            L.append(A)
        if not L: continue
        A = pd.concat(L, ignore_index=True)
        for c in F.WAVE:
            if c not in A: A[c] = np.nan
        for c in ("has_LRIV", "has_SNB"):
            A[c] = A[c].astype(float)
        parts[s] = A
        feats = [c for c in F.WAVE if A[c].astype(float).notna().any() and A[c].astype(float).nunique() > 1]   # all-NaN / constant columns of absent stations dropped
        X = A[feats].to_numpy(float); y = A.m_same.to_numpy(); dates = A.date.to_numpy()
        dl = np.array(sorted(set(dates)))
        fo = dict(zip(dl, np.random.default_rng(11).permutation(np.arange(len(dl)) % 5)))
        fold = np.array([fo[d] for d in dates])
        P = np.full(len(A), np.nan); dm = (A.series == "data").to_numpy()
        for k in range(5):
            tr = dm & (fold != k)
            m = F.fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
            P[fold == k] = m.predict_proba(X[fold == k])[:, 1]
        A["P"] = P
        D = A[dm]; Cc = A[~dm]
        sl, ic, _ = F.calib(D.m_same.to_numpy(), D.P.to_numpy())
        out.append(dict(sset=s, eval="oof_2011_2012test", n_feats=len(feats), n=len(D), n_pos=int(D.m_same.sum()), n_days=len(dl),
                        auc=F.auc(D.m_same, D.P), auc_net_mad=F.auc(D.m_same, D.net_mad), slope=sl,
                        matched_over_sumP=float(D.m_same.sum() / D.P.sum()), ctl_frac_P_gt05=float((Cc.P > 0.5).mean()),
                        ctl_over_data_sumP=float(Cc.P.sum() / D.P.sum())))
        G.log(f"[g1 model] {out[-1]}")
        if s == "CN4" and os.path.exists(NEW(2010, "CN4")):
            mfin = F.fit_cal(feats, X[dm], y[dm].astype(int), dates[dm])
            A0, _ = load(NEW(2010, "CN4")); A0 = scored(A0, 2010)
            for c in F.WAVE:
                if c not in A0: A0[c] = np.nan
            for c in ("has_LRIV", "has_SNB"):
                A0[c] = A0[c].astype(float)
            B = R.bostock([2010]); d0 = A0.series == "data"
            A0["m_same"] = np.nan
            A0.loc[d0, "m_same"] = (R.nearest_dt(A0.loc[d0, "t"].to_numpy(float), A0.loc[d0, "fam"].to_numpy(), B.t_origin.to_numpy(),
                                                 B.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
            A0["P"] = mfin.predict_proba(A0[feats].to_numpy(float))[:, 1]
            D = A0[d0]; Cc = A0[~d0]
            sl, ic, _ = F.calib(D.m_same.to_numpy(), D.P.to_numpy())
            out.append(dict(sset="CN4", eval="2010_out_of_sample", n=len(D), n_pos=int(D.m_same.sum()), n_days=D.date.nunique(),
                            auc=F.auc(D.m_same, D.P), auc_net_mad=F.auc(D.m_same, D.net_mad), slope=sl,
                            matched_over_sumP=float(D.m_same.sum() / D.P.sum()), ctl_frac_P_gt05=float((Cc.P > 0.5).mean()),
                            ctl_over_data_sumP=float(Cc.P.sum() / D.P.sum())))
            G.log(f"[g1 model] {out[-1]}")
    M = pd.DataFrame(out); M.to_csv(f"{TAB}/g1_model.csv", index=False)
    return M


def decide(Rr, Or, M):
    dec = {}
    for s, frac in (("CN4", 0.40), ("B2", 0.80), ("B1", 0.40)):
        g = Rr[(Rr.sset == s) & Rr.year.isin([2011, 2012])]
        if not len(g): continue
        pooled = {}
        for T, h in g.groupby("T_ref"):
            rr = h.matched_ref.sum() / h.n_events.sum(); rn = h.matched_new.sum() / h.n_events.sum()
            pooled[str(T)] = dict(recall_ref=rr, recall_new=rn, ratio=rn / rr, fam_red_ref=int(h.n_fam_red_ref.sum()),
                                  fam_red_new=int(h.n_fam_red_new.sum()),
                                  precision_ref=float((h.precision_ref * h.n_data_ref).sum() / h.n_data_ref.sum()),
                                  precision_new=float((h.precision_new * h.n_data_new).sum() / max(h.n_data_new.sum(), 1)))
        rec_ok = all(v["ratio"] >= frac for v in pooled.values())
        o = Or[(Or.sset == s) & Or.year.isin([2011, 2012])].groupby("fam").frac_days_observable.mean()
        n_obs = int((o >= 0.90).sum())
        fam_ok = all(v["fam_red_new"] >= 0.5 * v["fam_red_ref"] for v in pooled.values())
        d = dict(rule_recall_frac=frac, pooled=pooled, recall_ok=bool(rec_ok), n_fam_observable_ge2sites_90pct=n_obs,
                 obs_ok=bool(n_obs >= 100), fam_red_ok=bool(fam_ok))
        d["pass"] = bool(rec_ok and n_obs >= 100 and fam_ok) if s != "B2" else bool(rec_ok)
        if len(M):
            d["model"] = M[M.sset == s].to_dict("records")
        dec[s] = d
    if len(M):
        dec["S15_model_reference"] = M[M.sset == "S15"].to_dict("records")
    json.dump(dec, open(f"{TAB}/g1_decision.json", "w"), indent=1, default=float)
    G.log("[g1] decision " + json.dumps({k: v.get("pass") for k, v in dec.items() if isinstance(v, dict)}))
    return dec


def figure(Rr, Fr):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    cols = {"CN4": "#c0392b", "B1": "#e67e22", "B2": "#2980b9"}
    for s, g in Rr.groupby("sset"):
        for T, mk in ((8.25, "o"), (8.75, "s")):
            h = g[g.T_ref == T]
            ax[0].plot(h.year, h.recall_ratio, mk + "-", color=cols[s], label=f"{s}, matched to full set at {T}")
    ax[0].axhline(0.4, ls=":", c="k"); ax[0].axhline(0.8, ls="--", c="k", lw=.8)
    ax[0].set_xticks([2010, 2011, 2012]); ax[0].set_ylabel("recall / full-set recall (matched scrambled-control rate)")
    ax[0].set_ylim(0, 1.15); ax[0].legend(fontsize=7); ax[0].set_title("Recall of Bostock's events, reduced vs full network")
    h = Fr[(Fr.T_ref == 8.75) & Fr.year.isin([2011, 2012])].groupby(["sset", "fam"])[["matched_ref", "matched_new", "n_events"]].sum().reset_index()
    for s, g in h.groupby("sset"):
        ax[1].scatter(g.matched_ref / g.n_events.clip(lower=1), g.matched_new / g.n_events.clip(lower=1), s=10, color=cols[s], label=s, alpha=.7)
    ax[1].plot([0, .6], [0, .6], "k-", lw=.7); ax[1].plot([0, .6], [0, .24], "k:", lw=.7)
    ax[1].set_xlabel("per-family recall, full set (2011 + 2012, T 8.75)"); ax[1].set_ylabel("per-family recall, reduced set")
    ax[1].legend(fontsize=8); ax[1].set_title("Per family (dotted = 0.4 x full)")
    plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/g1_recall.png", dpi=130); plt.close()


if __name__ == "__main__":
    Rr, Fr, Or = recall_tables()
    M = model_tables() if "--no-model" not in sys.argv else pd.DataFrame()
    decide(Rr, Or, M)
    figure(Rr, Fr)
