#!/usr/bin/env python
"""Stage 6: catalogue flags + scoring.

Flags on every data detection (rows are never deleted, every rule is a column):
  tier_ext / tier_core   net_mad >= 8.25 / 8.75
  hf_pass                not (hfr defined AND the HF node has a cut AND hfr > cut); an undefined ratio
                         or an uncalibrated node is never vetoed (SVI rule)
  keep_xfam_*            SVI frozen cross-family rule: among families sharing >= 50 % of their template
                         stations, a detection is dropped if a station-sharing family has a higher raw
                         network cc within +-2 s (evaluated on each tier's set)
  cat_ext / cat_core     tier & hf_pass & keep_xfam
HF cut (validation mode): per HF station node, the 99th percentile of the ratio over BUILD-half
detections (>= T_ext) matched one-to-one (+-2 s) to an own-family array detection, >= 30 of them.

Validation (known answer): per family and window, on observable family-days the family's own array run
scanned, TEST half only: recall = own array detections matched one-to-one within +-2 s / own array
detections; precision_own = our catalogue detections matched to own / ours; precision_any = ours within
+-2 s of ANY WA family's array detection.  Scrambled control: site-rolled detections per family-day
(same tiers; the HF veto and cross-family rule cannot be applied to it, as in SVI) beside the data rate.

Trial (no truth): per family, observable days, live permanent sites, catalogue counts on tremor days vs
other days (tremor days = the region's PNSN tremor-day table, a LABEL only), scrambled-control rate.

  python scripts/wa_score.py --mode validation --mf val --tpl build
  python scripts/wa_score.py --mode trial --mf t2012p1 --tpl build --cuts results/hf_cuts.csv
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

LOG = f"{W.STUDY}/logs/score.log"
TREM = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}


def load_mf(mf):
    od = f"{W.RES}/mf_{mf}"
    DET = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/det_*.parquet"))], ignore_index=True)
    FD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/famday_*.parquet"))], ignore_index=True)
    H = pd.read_parquet(f"{od}/hf.parquet") if os.path.exists(f"{od}/hf.parquet") else pd.DataFrame(
        columns=["fid", "t", "hf_sta", "hf_node", "hf_fs", "hfr"])
    return DET, FD, H


def xfam_keep(D, share):
    """D: tier-set rows (fid, t, net_cc) sorted by t -> bool keep."""
    t = D.t.to_numpy(); f = D.fid.to_numpy(); c = D.net_cc.to_numpy()
    keep = np.ones(len(D), bool)
    lo = np.searchsorted(t, t - W.MATCH_S, "left"); hi = np.searchsorted(t, t + W.MATCH_S, "right")
    for i in range(len(D)):
        for j in range(lo[i], hi[i]):
            if j != i and f[j] != f[i] and c[j] > c[i] and (f[i], f[j]) in share:
                keep[i] = False
                break
    return keep


TIERS = ["ext", "core", "ext_ep", "core_ep"]
GRID = np.round(np.arange(8.0, 100.0001, 0.25), 2)


def frozen_th(windows):
    return {w: {"ext": W.T_EXT, "core": W.T_CORE, "ext_ep": W.T_EXT, "core_ep": W.T_CORE} for w in windows}


def derive_ep(D, FD, TH):
    """SVI c_tiers_by_epoch.py rule, one network epoch per window: T_core = smallest T on a 0.25-MAD grid
    >= 8.0 with quiet-day FA <= 0.1 per template-day, T_ext = smallest with quiet FA <= 1.0; FA counts
    data detections with hf_pass on the window's quiet (non-episode) observable family-days."""
    rows = []
    for w in sorted(TH):
        q = FD[(FD.window == w) & FD.observable & (FD.kind != "episode")]
        n = len(q)
        key = set(zip(q.fid, q.date))
        x = D[(D.window == w) & D.hf_pass]
        x = x[[k in key for k in zip(x.fid, x.date)]].net_mad.to_numpy()
        fa = np.array([(x >= T).sum() / n if n else np.nan for T in GRID])
        tc = GRID[np.argmax(fa <= 0.1)] if n and (fa <= 0.1).any() else np.nan
        te = GRID[np.argmax(fa <= 1.0)] if n and (fa <= 1.0).any() else np.nan
        # no quiet days -> keep the frozen pair (SVI anchor); quiet days but no T on the grid meets the
        # FA bar -> inf (nothing passes), recorded, never silently frozen
        TH[w]["core_ep"] = float(tc) if np.isfinite(tc) else (np.inf if n else W.T_CORE)
        TH[w]["ext_ep"] = float(te) if np.isfinite(te) else (np.inf if n else W.T_EXT)
        rows.append(dict(window=w, n_quiet_famdays=n, quiet_kinds=",".join(sorted(set(q.kind))),
                         T_ext_ep=TH[w]["ext_ep"], T_core_ep=TH[w]["core_ep"], derived=bool(n),
                         fa_at_8p25=float(fa[GRID == 8.25][0]) if n else np.nan,
                         fa_at_8p75=float(fa[GRID == 8.75][0]) if n else np.nan))
    return TH, pd.DataFrame(rows)


def flag(DET, H, Q, cuts, TH):
    D = DET[DET.series == "data"].copy()
    C = DET[DET.series == "control"].copy()
    D = D.merge(H[["fid", "t", "hf_node", "hfr"]], on=["fid", "t"], how="left")
    D["hf_cut"] = D.hf_node.map(dict(zip(cuts.hf_node, cuts.cut))) if len(cuts) else np.nan
    D["hf_pass"] = ~(D.hfr.notna() & D.hf_cut.notna() & (D.hfr > D.hf_cut))
    for X in (D, C):
        for tier in TIERS:
            thr = X.window.map({w: TH[w][tier] for w in TH}).astype(float)
            X[f"thr_{tier}"] = thr
            X[f"tier_{tier}"] = X.net_mad >= thr
    S = Q.groupby("fid").netsta.apply(lambda v: frozenset(v)).to_dict()
    share = set()
    fl = sorted(S)
    for i, a in enumerate(fl):
        for b in fl[i + 1:]:
            k = len(S[a] & S[b])
            if k and k >= 0.5 * min(len(S[a]), len(S[b])):
                share.add((a, b)); share.add((b, a))
    for tier in TIERS:
        m = (D[f"tier_{tier}"] & D.hf_pass).to_numpy()
        sub = D[m].sort_values("t")
        k = pd.Series(xfam_keep(sub, share), index=sub.index)
        D[f"keep_xfam_{tier}"] = False
        D.loc[k.index, f"keep_xfam_{tier}"] = k.to_numpy()
        D[f"cat_{tier}"] = D[f"tier_{tier}"] & D.hf_pass & D[f"keep_xfam_{tier}"]
    D["half"] = W.half_of(D.t.to_numpy())
    C["half"] = W.half_of(C.t.to_numpy())
    return D, C


def calibrate_cuts(D, T):
    rows = []
    B = D[(D.half == 0) & D.tier_ext & D.hfr.notna()]
    mt = []
    TB = {k: v.t.to_numpy() for k, v in T[T.half == 0].groupby(["fid", "date"])}
    for (f, ds), g in B.groupby(["fid", "date"]):
        tt = TB.get((f, ds), np.zeros(0))
        ia, _ = W.match_one_to_one(g.t.to_numpy(), tt)
        mt.append(g.iloc[ia])
    M = pd.concat(mt) if mt else B.iloc[:0]
    for nd, g in M.groupby("hf_node"):
        rows.append(dict(hf_node=nd, n_calib=len(g), cut=float(np.percentile(g.hfr, 99)) if len(g) >= 30 else np.nan,
                         median_matched=float(g.hfr.median())))
    for nd in set(B.hf_node.dropna()) - set(M.hf_node):
        rows.append(dict(hf_node=nd, n_calib=0, cut=np.nan, median_matched=np.nan))
    return pd.DataFrame(rows)


def validation(D, C, FD, T, F, scan):
    rows = []
    Tall = np.sort(T.t.to_numpy())
    E0 = D.iloc[:0]
    Dg = dict(tuple(D.groupby("fid"))); Cg = dict(tuple(C.groupby("fid")))
    Tg = dict(tuple(T.groupby("fid"))); FDg = dict(tuple(FD.groupby("fid")))
    for r in F.itertuples():
        sd = set(scan[scan.det_source == r.det_source].date)
        Df = Dg.get(r.fid, E0); Cf = Cg.get(r.fid, C.iloc[:0]); Tf = Tg.get(r.fid, T.iloc[:0])
        FDf = FDg.get(r.fid, FD.iloc[:0])
        for w in r.val_windows.split(","):
            fd = FDf[FDf.window == w]
            fdo = fd[fd.observable & fd.date.isin(sd) & (fd.kind == "episode")]
            days = set(fdo.date)
            tt = Tf[Tf.date.isin(days)]
            row = dict(fid=r.fid, family_uid=r.family_uid, region=r.region, window=w,
                       n_days_window=int((fd.kind == "episode").sum()), n_days_obs=len(days),
                       med_sites=float(fdo.n_sites.median()) if len(fdo) else 0.0,
                       n_truth_window=int(Tf.date.isin(set(fd.date)).sum()),
                       n_truth_test=int((tt.half == 1).sum()))
            for tier in TIERS:
                for lab, col in [("cat", f"cat_{tier}"), ("raw", f"tier_{tier}")]:
                    ours = Df[Df.date.isin(days) & Df[col]]
                    ot = ours[ours.half == 1]
                    nm = 0
                    for ds in days:
                        a_ = ot[ot.date == ds].t.to_numpy(); b_ = tt[(tt.date == ds) & (tt.half == 1)].t.to_numpy()
                        ia, _ = W.match_one_to_one(a_, b_)
                        nm += len(ia)
                    if len(ot):
                        x = ot.t.to_numpy(); j = np.clip(np.searchsorted(Tall, x), 1, len(Tall) - 1)
                        near = np.minimum(np.abs(Tall[j] - x), np.abs(Tall[j - 1] - x)) <= W.MATCH_S
                        nany = int(near.sum())
                    else:
                        nany = 0
                    k = f"{lab}_{tier}"
                    row[f"n_ours_test_{k}"] = len(ot)
                    row[f"n_match_{k}"] = nm
                    row[f"recall_{k}"] = nm / row["n_truth_test"] if row["n_truth_test"] else np.nan
                    row[f"prec_own_{k}"] = nm / len(ot) if len(ot) else np.nan
                    row[f"prec_any_{k}"] = nany / len(ot) if len(ot) else np.nan
                ct = Cf[Cf.date.isin(days) & Cf[f"tier_{tier}"] & (Cf.half == 1)]
                ncm = 0
                for ds in days:
                    ia, _ = W.match_one_to_one(ct[ct.date == ds].t.to_numpy(), tt[(tt.date == ds) & (tt.half == 1)].t.to_numpy())
                    ncm += len(ia)
                row[f"n_ctl_test_{tier}"] = len(ct)
                row[f"n_ctl_match_{tier}"] = ncm
                row[f"prec_own_ctl_{tier}"] = ncm / len(ct) if len(ct) else np.nan
                nctl = int((Cf.date.isin(days) & Cf[f"tier_{tier}"]).sum())
                ndat = int((Df.date.isin(days) & Df[f"tier_{tier}"]).sum())
                row[f"ctl_per_famday_{tier}"] = nctl / len(days) if days else np.nan
                row[f"data_per_famday_{tier}"] = ndat / len(days) if days else np.nan
                qd = set(fd[fd.observable & (fd.kind != "episode")].date)
                row[f"n_quiet_obs"] = len(qd)
                row[f"quiet_data_per_famday_{tier}"] = (int((Df.date.isin(qd) & Df[f"tier_{tier}"]).sum()) / len(qd)) if qd else np.nan
                row[f"quiet_ctl_per_famday_{tier}"] = (int((Cf.date.isin(qd) & Cf[f"tier_{tier}"]).sum()) / len(qd)) if qd else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_val(V):
    out = []
    for keys, g in [(("ALL", "ALL"), V)] + [((w, "ALL"), g) for w, g in V.groupby("window")] + \
            [((w, rg), g) for (w, rg), g in V.groupby(["window", "region"])]:
        g2 = g[g.n_days_obs > 0]
        row = dict(window=keys[0], region=keys[1], n_fam=len(g), n_fam_observable=len(g2),
                   n_truth_test=int(g2.n_truth_test.sum()))
        for k in [f"{l}_{t}" for t in TIERS for l in ("cat", "raw")]:
            nt, nm, no = g2.n_truth_test.sum(), g2[f"n_match_{k}"].sum(), g2[f"n_ours_test_{k}"].sum()
            row[f"recall_{k}"] = nm / nt if nt else np.nan
            row[f"prec_own_{k}"] = nm / no if no else np.nan
            row[f"prec_any_{k}"] = (g2[f"prec_any_{k}"] * g2[f"n_ours_test_{k}"]).sum() / no if no else np.nan
            row[f"med_fam_recall_{k}"] = float(g2[f"recall_{k}"].median()) if len(g2) else np.nan
        for tier in TIERS:
            nc = g2[f"n_ctl_test_{tier}"].sum()
            row[f"prec_own_ctl_{tier}"] = g2[f"n_ctl_match_{tier}"].sum() / nc if nc else np.nan
            row[f"recall_ctl_{tier}"] = g2[f"n_ctl_match_{tier}"].sum() / g2.n_truth_test.sum() if g2.n_truth_test.sum() else np.nan
            fdays = g2.n_days_obs.sum()
            row[f"data_per_famday_{tier}"] = (g2[f"data_per_famday_{tier}"] * g2.n_days_obs).sum() / fdays if fdays else np.nan
            row[f"ctl_per_famday_{tier}"] = (g2[f"ctl_per_famday_{tier}"] * g2.n_days_obs).sum() / fdays if fdays else np.nan
            q = g[g.n_quiet_obs > 0]
            qd = q.n_quiet_obs.sum()
            row[f"quiet_data_per_famday_{tier}"] = (q[f"quiet_data_per_famday_{tier}"] * q.n_quiet_obs).sum() / qd if qd else np.nan
            row[f"quiet_ctl_per_famday_{tier}"] = (q[f"quiet_ctl_per_famday_{tier}"] * q.n_quiet_obs).sum() / qd if qd else np.nan
        out.append(row)
    return pd.DataFrame(out)


def trial(D, C, FD, F):
    tdays = {}
    for rg, tr in TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        tdays[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    rows = []
    for r in F.itertuples():
        fd = FD[FD.fid == r.fid]
        fdo = fd[fd.observable]
        td = set(fdo.date) & tdays[r.region]
        od_ = set(fdo.date) - tdays[r.region]
        row = dict(fid=r.fid, family_uid=r.family_uid, region=r.region, n_days=len(fd), n_days_obs=len(fdo),
                   max_sites=int(fd.n_sites.max()) if len(fd) else 0,
                   med_sites_obs=float(fdo.n_sites.median()) if len(fdo) else 0.0,
                   ge2_sites=bool(len(fdo) > 0),
                   stations=",".join(sorted({s for v in fdo.stations for s in v.split(",") if s})),
                   n_tremor_days=len(td), n_other_days=len(od_))
        for tier in TIERS:
            mine = D[(D.fid == r.fid) & D[f"cat_{tier}"]]
            ctl = C[(C.fid == r.fid) & C[f"tier_{tier}"]]
            row[f"n_cat_{tier}"] = len(mine)
            row[f"n_cat_{tier}_tremor"] = int(mine.date.isin(td).sum())
            row[f"n_cat_{tier}_other"] = int(mine.date.isin(od_).sum())
            row[f"rate_{tier}_tremor"] = row[f"n_cat_{tier}_tremor"] / len(td) if td else np.nan
            row[f"rate_{tier}_other"] = row[f"n_cat_{tier}_other"] / len(od_) if od_ else np.nan
            row[f"ctl_rate_{tier}_tremor"] = int(ctl.date.isin(td).sum()) / len(td) if td else np.nan
            row[f"ctl_rate_{tier}_other"] = int(ctl.date.isin(od_).sum()) / len(od_) if od_ else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["validation", "trial"], required=True)
    ap.add_argument("--mf", required=True)
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--cuts", default=f"{W.RES}/hf_cuts.csv")
    a = ap.parse_args()
    F = W.load_families()
    DET, FD, H = load_mf(a.mf)
    Q = pd.concat([pd.read_csv(f"{W.RES}/templates/tpl_qc_{t}.csv") for t in a.tpl.split(",")])
    Q = Q[Q.clear.fillna(False).astype(bool)]
    od = f"{W.RES}/score_{a.mf}"
    os.makedirs(od, exist_ok=True)
    wins = sorted(set(DET.window)) if len(DET) else []
    if a.mode == "validation":
        T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
        scan = pd.read_csv(f"{W.SETUP}/scan_days.csv")
        D0, _ = flag(DET, H, Q, pd.DataFrame(columns=["hf_node", "cut"]), frozen_th(wins))
        cuts = calibrate_cuts(D0, T)
        cuts.to_csv(a.cuts, index=False)
        W.log(LOG, f"[score] HF cuts: {len(cuts)} nodes, {cuts.cut.notna().sum()} calibrated")
        D1, _ = flag(DET, H, Q, cuts, frozen_th(wins))
        TH, TT = derive_ep(D1, FD, frozen_th(wins))
        TT.to_csv(f"{od}/tiers_by_epoch.csv", index=False)
        W.log(LOG, "[score] epoch tiers\n" + TT.to_string())
        D, C = flag(DET, H, Q, cuts, TH)
        D.to_parquet(f"{od}/detections_flagged.parquet", index=False)
        V = validation(D, C, FD, T, F, scan)
        V.to_csv(f"{od}/validation_by_family.csv", index=False)
        S = summarize_val(V)
        S.to_csv(f"{od}/validation_summary.csv", index=False)
        W.log(LOG, "[score] validation summary\n" + S[["window", "region", "n_fam", "n_fam_observable", "n_truth_test",
                                                      "recall_cat_core", "prec_own_cat_core", "prec_own_ctl_core",
                                                      "recall_cat_core_ep", "prec_own_cat_core_ep", "prec_own_ctl_core_ep",
                                                      "data_per_famday_core", "ctl_per_famday_core"]].to_string())
    else:
        cuts = pd.read_csv(a.cuts) if os.path.exists(a.cuts) else pd.DataFrame(columns=["hf_node", "cut"])
        D1, _ = flag(DET, H, Q, cuts, frozen_th(wins))
        TH, TT = derive_ep(D1, FD, frozen_th(wins))
        TT.to_csv(f"{od}/tiers_by_epoch.csv", index=False)
        W.log(LOG, "[score] epoch tiers\n" + TT.to_string())
        D, C = flag(DET, H, Q, cuts, TH)
        D.to_parquet(f"{od}/detections_flagged.parquet", index=False)
        R = trial(D, C, FD, F)
        R.to_csv(f"{od}/trial_by_family.csv", index=False)
        W.log(LOG, f"[score] trial {a.mf}: families observable {int((R.n_days_obs > 0).sum())}, "
              f"cat_core {int(R.n_cat_core.sum())} (tremor {int(R.n_cat_core_tremor.sum())} / other {int(R.n_cat_core_other.sum())}); "
              f"cat_core_ep {int(R.n_cat_core_ep.sum())} (tremor {int(R.n_cat_core_ep_tremor.sum())} / other {int(R.n_cat_core_ep_other.sum())})")


if __name__ == "__main__":
    main()
