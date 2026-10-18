#!/usr/bin/env python
"""GO/NO-GO Step B scorer (2026-10-02) = merlin_t1_score.py (T1) restricted to V2010 (test half) + A2011, with:
  * --res / --mf : which matched-filter run (results/mf_mt1_all = original 42-station T1; results/wa_gonogo/mf_core21)
  * --stations   : templates (Q: HF station pool, cross-family rule) restricted to the stations of that network
  * HF cuts calibrated on the V2010 build half only (V2008 is not in this test) -- the same for both runs
  * separation vs BOTH controls: sepboth = (data+1)/(max(control, control2)+1) on label days (T1 used control 1)
  * per family: >= 30-event templates and clear sites on the 21 core stations (columns)
Everything else (tiers, flags, one-to-one +-2 s scoring, levels own / nbr3 / any) is the T1 code unchanged.

Original T1 docstring follows.
Merlin test T1 scoring — clear-only vs ALL-templates network matched filter, three matching levels.

Same flags as wa_score.py (frozen T_ext 8.25 / T_core 8.75 + per-epoch tiers by the SVI FA-only rule, HF veto
calibrated per HF node on V2008/V2010 build-half own matches, frozen cross-family rule on the template stations
actually used), same observable family-days (>= 2 sites with a template, episode days the family's own array
run scanned), same one-to-one +-2 s own-family match.  Scored half: TEST half for V2008 / V2010 (templates are
from the build half), ALL of A2011 (independent episode).  A2011 has no quiet days: its per-epoch tiers are
V2010's (nearest network epoch; the frozen pair is reported beside it).

Matching levels (columns, nothing dropped):
  own      one-to-one +-2 s with the family's own array detections (WA_TRIAL / STEP4B definition)
  nbr3     any family within 3 km (epicentral, lat_use/lon_use; non-transitive neighbourhood, so no chaining;
           unlocated families use their approximate position; the one family without a position = own only)
  any      any of the 329 families
  recall_L = own scored truth with >= 1 of OUR scored detections of a family in L(f) within +-2 s
  prec_L   = our scored detections within +-2 s of an array detection of a family in L(f)
  (own level: one-to-one, as before; nbr3 / any: nearest-within-2 s, as the trial's prec_any)
Scrambled control: the control's detections (tier only) scored by the same rules = the chance level per level.
Separation (STEP4B): per family, observable label days, raw tier counts: sep = (data+1)/(control+1), and the
same on the control = (control2+1)/(control+1).  Label = PNSN tremor days of the region (2010, 2011); 2008 has
no catalogue -> the window's episode days.  Label only, never a gate.

  WA_SETUP=results/merlin_tests/t1_setup python scripts/merlin_t1_score.py --mf mt1_all --tplmode all
  -> results/merlin_tests/t1/score_<mf>/{by_family.csv, summary.csv, tiers_by_epoch.csv, hf_cuts.csv}
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import wa_score as WS

LOG = f"{W.STUDY}/logs/gonogo_score.log"
TIERS = ["ext", "core", "ext_ep", "core_ep"]
LEVELS = ["own", "nbr3", "any"]
SEPS = [1.5, 2, 3, 5]
WINS = ["V2010", "A2011"]


def near(x, ref, tol=W.MATCH_S):
    """bool: x within tol of any element of sorted ref"""
    x = np.asarray(x, float)
    if len(ref) == 0 or len(x) == 0:
        return np.zeros(len(x), bool)
    if len(ref) == 1:
        return np.abs(ref[0] - x) <= tol
    j = np.clip(np.searchsorted(ref, x), 1, len(ref) - 1)
    d = np.minimum(np.abs(ref[j] - x), np.abs(ref[j - 1] - x))
    return d <= tol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mf", required=True)
    ap.add_argument("--tplmode", choices=["clear", "all"], required=True)
    ap.add_argument("--tpl", default="build")
    ap.add_argument("--res", required=True, help="directory holding mf_<mf>/ and templates/")
    ap.add_argument("--stations", required=True, help="stations.csv of the network (restricts Q)")
    ap.add_argument("--od", required=True)
    ap.add_argument("--core", default=f"{W.RES}/wa_gonogo/setup/stations.csv")
    ap.add_argument("--wins", default="V2010,A2011", help="supplement: V2008 (SW truth exists only 2007-08)")
    a = ap.parse_args()
    global WINS
    WINS = a.wins.split(",")
    od = a.od
    os.makedirs(od, exist_ok=True)
    W.RES = a.res
    F = W.load_families()
    T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    scan = pd.read_csv(f"{W.SETUP}/scan_days.csv")
    DET, FD, H = WS.load_mf(a.mf)
    DET = DET[DET.window.isin(WINS)]
    FD = FD[FD.window.isin(WINS)]
    Q = pd.concat([pd.read_csv(f"{W.RES}/templates/tpl_qc_{t}.csv") for t in a.tpl.split(",")])
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    Qall = Q.copy()
    Q = Q[Q.clear] if a.tplmode == "clear" else Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna() & Q.w.notna()]
    NETW = set(pd.read_csv(a.stations).netsta)
    Q = Q[Q.netsta.isin(NETW)]
    CORE = pd.read_csv(a.core)
    CSITE = dict(zip(CORE.netsta, CORE.site))
    q30 = Qall[(Qall.n_events >= W.MIN_EVENTS) & Qall.offset_s.notna() & Qall.w.notna() & Qall.netsta.isin(CSITE)]
    qcl = Qall[Qall.clear & Qall.netsta.isin(CSITE)]
    CT = dict(n_tpl30_core=q30.groupby("fid").size(),
              n_sites30_core=q30.groupby("fid").netsta.apply(lambda v: len({CSITE[x] for x in v})),
              n_clear_sites_core=qcl.groupby("fid").netsta.apply(lambda v: len({CSITE[x] for x in v})),
              tpl30_core_stations=q30.groupby("fid").netsta.apply(lambda v: ",".join(sorted(set(v)))))
    wins = sorted(set(DET.window))
    # ---- flags (wa_score): HF cut calibration on V2008/V2010 build half only
    D0, _ = WS.flag(DET, H, Q, pd.DataFrame(columns=["hf_node", "cut"]), WS.frozen_th(wins))
    CW = [w for w in WINS if w in ("V2008", "V2010")]       # build-half calibration windows present
    Tc = T[T.window.isin(CW)]
    cuts = WS.calibrate_cuts(D0[D0.window.isin(CW)], Tc)
    cuts.to_csv(f"{od}/hf_cuts.csv", index=False)
    D1, _ = WS.flag(DET, H, Q, cuts, WS.frozen_th(wins))
    TH, TT = WS.derive_ep(D1, FD, WS.frozen_th(wins))
    if "A2011" in TH and "V2010" in TH:
        TH["A2011"]["core_ep"] = TH["V2010"]["core_ep"]; TH["A2011"]["ext_ep"] = TH["V2010"]["ext_ep"]
        TT.loc[TT.window == "A2011", ["T_core_ep", "T_ext_ep"]] = [TH["V2010"]["core_ep"], TH["V2010"]["ext_ep"]]
        TT.loc[TT.window == "A2011", "quiet_kinds"] = "none: V2010 epoch tiers used"
    TT.to_csv(f"{od}/tiers_by_epoch.csv", index=False)
    W.log(LOG, f"[t1:{a.mf}] tiers\n{TT.to_string()}")
    D, C = WS.flag(DET, H, Q, cuts, TH)
    C2 = DET[DET.series == "control2"].copy()
    for tier in TIERS:
        C2[f"tier_{tier}"] = C2.net_mad >= C2.window.map({w: TH[w][tier] for w in TH}).astype(float)
    for X in (D, C, C2):
        X["half"] = W.half_of(X.t.to_numpy()) if len(X) else []
        X["scored"] = (X.half == 1) | (X.window == "A2011")
    T = T.copy()
    T["scored"] = (T.half == 1) | (T.window == "A2011")
    has_c2 = len(C2) > 0
    # ---- neighbourhoods
    la, lo = F.lat_use.to_numpy(float), F.lon_use.to_numpy(float)
    dk = W.km(la[:, None], lo[:, None], la[None], lo[None])
    NB = {}
    for i, f in enumerate(F.fid):
        m = dk[i] <= 3.0
        m[i] = True
        NB[f] = set(F.fid[np.nan_to_num(m, nan=False).astype(bool)])
    # ---- tremor-day label
    tdays = {}
    for rg, tr in WS.TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        tdays[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    # ---- per-day pools
    Tday = {d: g for d, g in T.groupby("date")}
    catcols = [f"cat_{t}" for t in TIERS] + [f"tier_{t}" for t in TIERS]
    Dday = {d: g for d, g in D.groupby("date")}
    Cday = {d: g for d, g in C.groupby("date")}
    E0 = D.iloc[:0]; C0 = C.iloc[:0]; T0 = T.iloc[:0]
    FDg = dict(tuple(FD.groupby("fid")))
    reg = dict(zip(F.fid, F.region))
    rows = []
    for r in F.itertuples():
        sd = set(scan[scan.det_source == r.det_source].date)
        fdf = FDg.get(r.fid, FD.iloc[:0])
        for w in r.val_windows.split(","):
            if w not in wins:
                continue
            fd = fdf[fdf.window == w]
            fdo = fd[fd.observable & fd.date.isin(sd) & (fd.kind == "episode")]
            days = sorted(set(fdo.date))
            alld = sorted(set(fd[fd.date.isin(sd) & (fd.kind == "episode")].date))
            row = dict(fid=r.fid, family_uid=r.family_uid, region=r.region, window=w, n_days_scan=len(alld),
                       n_days_obs=len(days), med_sites=float(fdo.n_sites.median()) if len(fdo) else 0.0,
                       n_nbr3=len(NB[r.fid]) - 1,
                       **{k: (v.get(r.fid, 0) if k != "tpl30_core_stations" else v.get(r.fid, "")) for k, v in CT.items()})
            ntt = 0
            ntall = 0
            acc = {}
            for ds in alld:
                tg = Tday.get(ds, T0)
                own_all = tg[(tg.fid == r.fid) & tg.scored].t.to_numpy()
                ntall += len(own_all)
            row["n_truth_scored_all"] = ntall
            for ds in days:
                tg = Tday.get(ds, T0)
                tsc = tg[tg.scored]
                own = np.sort(tsc[tsc.fid == r.fid].t.to_numpy())
                ntt += len(own)
                t_own_all = np.sort(tg[tg.fid == r.fid].t.to_numpy())
                t_nb = np.sort(tg[tg.fid.isin(NB[r.fid])].t.to_numpy())
                t_any = np.sort(tg.t.to_numpy())
                dg = Dday.get(ds, E0); dg = dg[dg.scored]
                cg = Cday.get(ds, C0); cg = cg[cg.scored]
                for src, pool, cols in [("d", dg, [(f"cat_{t}", f"cat_{t}") for t in TIERS] + [(f"raw_{t}", f"tier_{t}") for t in TIERS]),
                                        ("c", cg, [(f"ctl_{t}", f"tier_{t}") for t in TIERS])]:
                    for key, col in cols:
                        P = pool[pool[col]]
                        mine = np.sort(P[P.fid == r.fid].t.to_numpy())
                        pnb = np.sort(P[P.fid.isin(NB[r.fid])].t.to_numpy())
                        pany = np.sort(P.t.to_numpy())
                        ia, _ = W.match_one_to_one(mine, own)
                        k = acc.setdefault(key, dict(n_ours=0, m_own=0, p_nbr3=0, p_any=0, r_nbr3=0, r_any=0))
                        k["n_ours"] += len(mine)
                        k["m_own"] += len(ia)
                        k["p_nbr3"] += int(near(mine, t_nb).sum())
                        k["p_any"] += int(near(mine, t_any).sum())
                        k["r_nbr3"] += int(near(own, pnb).sum())
                        k["r_any"] += int(near(own, pany).sum())
            row["n_truth_test"] = ntt
            for key, k in acc.items():
                row[f"n_ours_{key}"] = k["n_ours"]
                row[f"n_match_own_{key}"] = k["m_own"]
                row[f"recall_own_{key}"] = k["m_own"] / ntt if ntt else np.nan
                row[f"recall_nbr3_{key}"] = k["r_nbr3"] / ntt if ntt else np.nan
                row[f"recall_any_{key}"] = k["r_any"] / ntt if ntt else np.nan
                row[f"n_rec_nbr3_{key}"] = k["r_nbr3"]; row[f"n_rec_any_{key}"] = k["r_any"]
                row[f"prec_own_{key}"] = k["m_own"] / k["n_ours"] if k["n_ours"] else np.nan
                row[f"prec_nbr3_{key}"] = k["p_nbr3"] / k["n_ours"] if k["n_ours"] else np.nan
                row[f"prec_any_{key}"] = k["p_any"] / k["n_ours"] if k["n_ours"] else np.nan
                row[f"n_prec_nbr3_{key}"] = k["p_nbr3"]; row[f"n_prec_any_{key}"] = k["p_any"]
            # rates (raw tiers, all halves, observable scanned episode days) + separation on label days
            ddf = D[D.fid == r.fid]; cdf = C[C.fid == r.fid]; c2f = C2[C2.fid == r.fid] if has_c2 else C2
            fo = fd[fd.observable]
            lab = set(fo[fo.kind == "episode"].date) if w == "V2008" else set(fo.date) & tdays[r.region]
            row["n_label_days"] = len(lab)
            for tier in TIERS:
                nd = int((ddf.date.isin(days) & ddf[f"tier_{tier}"]).sum())
                nc = int((cdf.date.isin(days) & cdf[f"tier_{tier}"]).sum())
                row[f"data_per_famday_{tier}"] = nd / len(days) if days else np.nan
                row[f"ctl_per_famday_{tier}"] = nc / len(days) if days else np.nan
                ld = int((ddf.date.isin(lab) & ddf[f"tier_{tier}"]).sum())
                lc = int((cdf.date.isin(lab) & cdf[f"tier_{tier}"]).sum())
                lc2 = int((c2f.date.isin(lab) & c2f[f"tier_{tier}"]).sum()) if has_c2 else np.nan
                row[f"lab_data_{tier}"] = ld; row[f"lab_ctl_{tier}"] = lc; row[f"lab_ctl2_{tier}"] = lc2
                row[f"sep_{tier}"] = (ld + 1) / (lc + 1) if lab else np.nan
                row[f"sep_ctl_{tier}"] = (lc2 + 1) / (lc + 1) if (lab and has_c2) else np.nan
                row[f"sepboth_{tier}"] = (ld + 1) / (max(lc, lc2) + 1) if (lab and has_c2) else np.nan
            rows.append(row)
    V = pd.DataFrame(rows)
    V["run"] = a.mf
    V.to_csv(f"{od}/by_family.csv", index=False)
    # ---- summary
    out = []
    groups = [(("ALL", "ALL"), V)] + [((w, "ALL"), g) for w, g in V.groupby("window")] + \
        [(("ALL", rg), g) for rg, g in V.groupby("region")] + [((w, rg), g) for (w, rg), g in V.groupby(["window", "region"])]
    for (w, rg), g in groups:
        g2 = g[g.n_days_obs > 0]
        row = dict(run=a.mf, window=w, region=rg, n_famwin=len(g), n_famwin_observable=len(g2),
                   n_fam_observable=g2.fid.nunique(), n_truth_test=int(g2.n_truth_test.sum()),
                   n_truth_scored_all=int(g.n_truth_scored_all.sum()), famdays_obs=int(g2.n_days_obs.sum()))
        nt = row["n_truth_test"]
        for key in [f"cat_{t}" for t in TIERS] + [f"raw_{t}" for t in TIERS] + [f"ctl_{t}" for t in TIERS]:
            if f"n_ours_{key}" not in g2:
                continue
            no = g2[f"n_ours_{key}"].sum()
            row[f"n_ours_{key}"] = int(no)
            for lev, cm, cp in [("own", f"n_match_own_{key}", f"n_match_own_{key}"),
                                ("nbr3", f"n_rec_nbr3_{key}", f"n_prec_nbr3_{key}"),
                                ("any", f"n_rec_any_{key}", f"n_prec_any_{key}")]:
                row[f"recall_{lev}_{key}"] = g2[cm].sum() / nt if nt else np.nan
                row[f"recall_all_{lev}_{key}"] = g2[cm].sum() / row["n_truth_scored_all"] if row["n_truth_scored_all"] else np.nan
                row[f"prec_{lev}_{key}"] = g2[cp].sum() / no if no else np.nan
        fdays = g2.n_days_obs.sum()
        for tier in TIERS:
            row[f"data_per_famday_{tier}"] = (g2[f"data_per_famday_{tier}"] * g2.n_days_obs).sum() / fdays if fdays else np.nan
            row[f"ctl_per_famday_{tier}"] = (g2[f"ctl_per_famday_{tier}"] * g2.n_days_obs).sum() / fdays if fdays else np.nan
            o = g[g.n_label_days > 0]
            for s in SEPS:
                row[f"sep{s:g}_{tier}_data"] = int((o[f"sep_{tier}"] >= s).sum())
                row[f"sep{s:g}_{tier}_ctl"] = int((o[f"sep_ctl_{tier}"] >= s).sum()) if has_c2 else np.nan
                row[f"sepboth{s:g}_{tier}_data"] = int((o[f"sepboth_{tier}"] >= s).sum()) if has_c2 else np.nan
                row[f"sepboth{s:g}_{tier}_nfam"] = int(o[o[f"sepboth_{tier}"] >= s].fid.nunique()) if has_c2 else np.nan
                row[f"sep{s:g}_{tier}_nfam"] = int(o[o[f"sep_{tier}"] >= s].fid.nunique())
                row[f"sep{s:g}_{tier}_ctl_nfam"] = int(o[o[f"sep_ctl_{tier}"] >= s].fid.nunique()) if has_c2 else np.nan
        out.append(row)
    S = pd.DataFrame(out)
    S.to_csv(f"{od}/summary.csv", index=False)
    W.log(LOG, f"[t1:{a.mf}] summary\n" + S[["window", "region", "n_fam_observable", "n_truth_test", "n_truth_scored_all",
                                                 "recall_own_cat_core", "prec_own_cat_core", "recall_nbr3_cat_core",
                                                 "prec_nbr3_cat_core", "prec_any_cat_core", "prec_own_ctl_core",
                                                 "data_per_famday_core", "ctl_per_famday_core",
                                                 "sep2_core_data", "sep2_core_ctl"]].to_string())


if __name__ == "__main__":
    main()
