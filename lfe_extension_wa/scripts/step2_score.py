#!/usr/bin/env python
"""STEP 2 (Merlin, 2026-09-30) part 2 — GROUP detection scored against the union of its members' array detections,
compared with the T1 family-level detection summed over the members.

Runs compared (same days, same truth, same rules):
  grp      the group matched filter (results/mf_step2grp: group templates, the frozen T1 chain, 42 permanent stations,
           --drop-dead 0.05, the SAME site-rolled scrambled controls as T1 — same seed, same loaded stations)
  t1union  T1 (results/mf_mt1_all) detections of the member families, each moved to the group frame (t + delta_m) and
           merged within +-2 s (start-anchored) = "family-level detection summed over members"; T1's own HF cuts and
           per-epoch tiers (results/merlin_tests/t1/score_mt1_all)
Flags (wa_score): frozen T_ext 8.25 / T_core 8.75 + per-epoch tiers (SVI FA-only rule on the run's quiet days; A2011 has
none -> V2010's), HF veto per HF node (group run: calibrated on its own build-half matches to the group truth), frozen
cross-family rule (cat) and tier-only (raw).  Truth = the group union (results/step2_groups/setup/truth.parquet).
Observable group-day: grp = >= 2 sites with a group template (wa_mf); t1union = >= 1 member observable that day.
Scored half: TEST half (group frame) for V2008 / V2010, all of A2011.  Scanned group-days: setup/scan_days.csv.
Separation (label days: V2008 = observable episode days; else observable & PNSN tremor days of the group's region; a
label only): sep = (data + 1) / (control + 1) and on the control (control2 + 1) / (control + 1), raw tier counts.
Also: recall of truth whose representative member is well aligned (ref, co-detection prior, or stack-cc ratio >= 1.2).

  WA_SETUP=results/step2_groups/setup python scripts/step2_score.py
  -> results/step2_groups/score/{by_group.csv, summary.csv, sep_curve.csv, tiers_by_epoch.csv, hf_cuts.csv, singleton_check.csv}
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import wa_score as WS

G = f"{W.RES}/step2_groups"
OD = f"{G}/score"
LOG = f"{W.STUDY}/logs/step2.log"
TIERS = ["ext", "core", "ext_ep", "core_ep"]
WINS = ["V2008", "V2010", "A2011"]
SEP_GRID = [1.0, 1.25, 1.5, 2, 3, 5, 10]


def merge_times(t, v):
    """start-anchored +-2 s merge of sorted times; returns representative times (max v in each cluster)"""
    if len(t) == 0:
        return t
    o = np.argsort(t); t = t[o]; v = v[o]
    out = []
    i = 0
    while i < len(t):
        j = i + 1
        while j < len(t) and t[j] - t[i] <= W.MATCH_S:
            j += 1
        out.append(t[i + int(np.argmax(v[i:j]))])
        i = j
    return np.array(out)


def near(x, ref, tol=W.MATCH_S):
    x = np.asarray(x, float)
    if len(ref) == 0 or len(x) == 0:
        return np.zeros(len(x), bool)
    j = np.clip(np.searchsorted(ref, x), 1, max(len(ref) - 1, 1))
    if len(ref) == 1:
        return np.abs(ref[0] - x) <= tol
    return np.minimum(np.abs(ref[j] - x), np.abs(ref[j - 1] - x)) <= tol


def pools(D, cols):
    """{col: {(gid, date): merged times}}"""
    out = {}
    for c in cols:
        x = D[D[c]]
        out[c] = {k: merge_times(g.t.to_numpy(), g.net_mad.to_numpy()) for k, g in x.groupby(["gid", "date"])}
    return out


def flags_group(mf, TG):
    DET, FD, H = WS.load_mf(mf)
    DET = DET[DET.window.isin(WINS)]; FD = FD[FD.window.isin(WINS)]
    Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_step2grp.csv")
    Q = Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna() & Q.w.notna()]
    wins = sorted(set(DET.window))
    D0, _ = WS.flag(DET, H, Q, pd.DataFrame(columns=["hf_node", "cut"]), WS.frozen_th(wins))
    Tc = TG[TG.window.isin(["V2008", "V2010"])]
    cuts = WS.calibrate_cuts(D0[D0.window.isin(["V2008", "V2010"])], Tc)
    D1, _ = WS.flag(DET, H, Q, cuts, WS.frozen_th(wins))
    TH, TT = WS.derive_ep(D1, FD, WS.frozen_th(wins))
    if "A2011" in TH and "V2010" in TH:
        TH["A2011"]["core_ep"] = TH["V2010"]["core_ep"]; TH["A2011"]["ext_ep"] = TH["V2010"]["ext_ep"]
        TT.loc[TT.window == "A2011", ["T_core_ep", "T_ext_ep"]] = [TH["V2010"]["core_ep"], TH["V2010"]["ext_ep"]]
        TT.loc[TT.window == "A2011", "quiet_kinds"] = "none: V2010 epoch tiers used"
    D, C = WS.flag(DET, H, Q, cuts, TH)
    C2 = DET[DET.series == "control2"].copy()
    for tier in TIERS:
        C2[f"tier_{tier}"] = C2.net_mad >= C2.window.map({w: TH[w][tier] for w in TH}).astype(float)
    for X in (D, C, C2):
        X["gid"] = X.fid
    obs = FD[FD.observable].assign(gid=FD.fid)[["gid", "date", "window", "kind"]].drop_duplicates()
    return D, C, C2, obs, cuts, TT


def flags_t1(M):
    DET, FD, H = WS.load_mf("mt1_all")
    DET = DET[DET.window.isin(WINS)]; FD = FD[FD.window.isin(WINS)]
    Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
    Q = Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna() & Q.w.notna()]
    sd = f"{W.RES}/merlin_tests/t1/score_mt1_all"
    cuts = pd.read_csv(f"{sd}/hf_cuts.csv")
    TT = pd.read_csv(f"{sd}/tiers_by_epoch.csv")
    TH = {r.window: {"ext": W.T_EXT, "core": W.T_CORE, "ext_ep": float(r.T_ext_ep), "core_ep": float(r.T_core_ep)}
          for r in TT.itertuples()}
    D, C = WS.flag(DET, H, Q, cuts, TH)
    C2 = DET[DET.series == "control2"].copy()
    for tier in TIERS:
        C2[f"tier_{tier}"] = C2.net_mad >= C2.window.map({w: TH[w][tier] for w in TH}).astype(float)
    g_of = dict(zip(M.fid, M.gid)); dl = dict(zip(M.fid, M.delta_s))
    for X in (D, C, C2):
        X["gid"] = X.fid.map(g_of)
        X["t"] = X.t + X.fid.map(dl)
    fo = FD[FD.observable].copy()
    fo["gid"] = fo.fid.map(g_of)
    obs = fo[["gid", "date", "window", "kind"]].drop_duplicates()
    return D, C, C2, obs


def score(run, D, C, C2, obs, GR, TG, scan, tdays, align_ok):
    for X in (D, C, C2):
        X["half"] = W.half_of(X.t.to_numpy()) if len(X) else []
        X["scored"] = (X.half == 1) | (X.window == "A2011")
    TG = TG.copy(); TG["scored"] = (TG.half == 1) | (TG.window == "A2011")
    TG["align_ok"] = TG.member.map(align_ok).fillna(False).astype(bool)
    dcols = [f"cat_{t}" for t in TIERS] + [f"tier_{t}" for t in TIERS]
    PD = pools(D, dcols); PC = pools(C, [f"tier_{t}" for t in TIERS]); PC2 = pools(C2, [f"tier_{t}" for t in TIERS]) if len(C2) else {}
    PDs = pools(D[D.scored], dcols); PCs = pools(C[C.scored], [f"tier_{t}" for t in TIERS])
    Tday = {k: g for k, g in TG.groupby(["fid", "date"])}
    Tall = np.sort(TG.t.to_numpy())
    sd = {s: set(g.date) for s, g in scan.groupby("det_source")}
    OB = {g: x for g, x in obs.groupby("gid")}
    rows = []
    E0 = np.zeros(0)
    for r in GR.itertuples():
        g = r.fid
        for w in r.val_windows.split(","):
            wd = {d for d, k in W.window_days(w) if k == "episode"}
            scd = sd.get(r.det_source, set()) & wd
            ob = OB.get(g)
            obs_w = set(ob[(ob.window == w)].date) if ob is not None else set()
            days = sorted(scd & obs_w)
            row = dict(run=run, gid=g, region=r.region, n_members=r.n_members, window=w, n_days_scan=len(scd),
                       n_days_obs=len(days))
            tru_all = sum(int(Tday[(g, d)].scored.sum()) for d in scd if (g, d) in Tday)
            row["n_truth_scored_all"] = tru_all
            nt = 0; nt_al = 0
            acc = {}
            for d in days:
                tg = Tday.get((g, d))
                own = np.sort(tg[tg.scored].t.to_numpy()) if tg is not None else E0
                own_al = np.sort(tg[tg.scored & tg.align_ok].t.to_numpy()) if tg is not None else E0
                nt += len(own); nt_al += len(own_al)
                for key, pool in [(f"cat_{t}", PDs[f"cat_{t}"]) for t in TIERS] + [(f"raw_{t}", PDs[f"tier_{t}"]) for t in TIERS] + \
                        [(f"ctl_{t}", PCs[f"tier_{t}"]) for t in TIERS]:
                    mine = pool.get((g, d), E0)
                    ia, _ = W.match_one_to_one(mine, own)
                    ib, _ = W.match_one_to_one(mine, own_al)
                    k = acc.setdefault(key, [0, 0, 0, 0])
                    k[0] += len(mine); k[1] += len(ia); k[2] += int(near(mine, Tall).sum()); k[3] += len(ib)
            row["n_truth_test"] = nt; row["n_truth_test_alignok"] = nt_al
            for key, (no, nm, na, nal) in acc.items():
                row[f"n_ours_{key}"] = no; row[f"n_match_{key}"] = nm; row[f"n_any_{key}"] = na; row[f"n_match_al_{key}"] = nal
                row[f"recall_{key}"] = nm / nt if nt else np.nan
                row[f"recall_alignok_{key}"] = nal / nt_al if nt_al else np.nan
                row[f"recall_all_{key}"] = nm / tru_all if tru_all else np.nan
                row[f"prec_{key}"] = nm / no if no else np.nan
                row[f"prec_any_{key}"] = na / no if no else np.nan
            lab = set(days) if w == "V2008" else (set(obs_w) & tdays[r.region])
            row["n_label_days"] = len(lab)
            for tier in TIERS:
                nd = sum(len(PD[f"tier_{tier}"].get((g, d), E0)) for d in days)
                nc = sum(len(PC[f"tier_{tier}"].get((g, d), E0)) for d in days)
                row[f"data_per_gday_{tier}"] = nd / len(days) if days else np.nan
                row[f"ctl_per_gday_{tier}"] = nc / len(days) if days else np.nan
                ld = sum(len(PD[f"tier_{tier}"].get((g, d), E0)) for d in lab)
                lc = sum(len(PC[f"tier_{tier}"].get((g, d), E0)) for d in lab)
                lc2 = sum(len(PC2[f"tier_{tier}"].get((g, d), E0)) for d in lab) if PC2 else np.nan
                row[f"lab_data_{tier}"] = ld; row[f"lab_ctl_{tier}"] = lc; row[f"lab_ctl2_{tier}"] = lc2
                row[f"sep_{tier}"] = (ld + 1) / (lc + 1) if lab else np.nan
                row[f"sep_ctl_{tier}"] = (lc2 + 1) / (lc + 1) if (lab and PC2) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def summarize(V):
    out = []
    V = V.copy()
    V["size_class"] = V.n_members.map(lambda k: "1" if k == 1 else ("2" if k == 2 else ("3-6" if k <= 6 else "23")))
    grp = [(("ALL", "ALL", "ALL"), V)] + [((w, "ALL", "ALL"), g) for w, g in V.groupby("window")] + \
        [(("ALL", rg, "ALL"), g) for rg, g in V.groupby("region")] + [(("ALL", "ALL", sc), g) for sc, g in V.groupby("size_class")] + \
        [((w, rg, "ALL"), g) for (w, rg), g in V.groupby(["window", "region"])]
    for run, VR in V.groupby("run"):
        for (w, rg, sc), g0 in grp:
            g = g0[g0.run == run]
            g2 = g[g.n_days_obs > 0]
            row = dict(run=run, window=w, region=rg, size_class=sc, n_groupwin=len(g), n_groupwin_obs=len(g2),
                       n_groups_obs=g2.gid.nunique(), n_families_obs=int(g2.drop_duplicates("gid").n_members.sum()),
                       n_truth_test=int(g2.n_truth_test.sum()), n_truth_scored_all=int(g.n_truth_scored_all.sum()),
                       n_truth_alignok=int(g2.n_truth_test_alignok.sum()), gdays_obs=int(g2.n_days_obs.sum()))
            nt = row["n_truth_test"]
            for key in [f"{a}_{t}" for a in ("cat", "raw", "ctl") for t in TIERS]:
                if f"n_ours_{key}" not in g2:
                    continue
                no = g2[f"n_ours_{key}"].sum(); nm = g2[f"n_match_{key}"].sum()
                row[f"n_ours_{key}"] = int(no)
                row[f"recall_{key}"] = nm / nt if nt else np.nan
                row[f"recall_all_{key}"] = nm / row["n_truth_scored_all"] if row["n_truth_scored_all"] else np.nan
                row[f"recall_alignok_{key}"] = g2[f"n_match_al_{key}"].sum() / row["n_truth_alignok"] if row["n_truth_alignok"] else np.nan
                row[f"prec_{key}"] = nm / no if no else np.nan
                row[f"prec_any_{key}"] = g2[f"n_any_{key}"].sum() / no if no else np.nan
            fd = g2.n_days_obs.sum()
            for tier in TIERS:
                row[f"data_per_gday_{tier}"] = (g2[f"data_per_gday_{tier}"] * g2.n_days_obs).sum() / fd if fd else np.nan
                row[f"ctl_per_gday_{tier}"] = (g2[f"ctl_per_gday_{tier}"] * g2.n_days_obs).sum() / fd if fd else np.nan
                o = g[g.n_label_days > 0]
                for s in SEP_GRID[2:]:
                    row[f"sep{s:g}_{tier}_data"] = int((o[f"sep_{tier}"] >= s).sum())
                    row[f"sep{s:g}_{tier}_ctl"] = int((o[f"sep_ctl_{tier}"] >= s).sum())
                row[f"n_label_groupwin_{tier}"] = len(o)
            out.append(row)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mf", default="step2grp")
    a = ap.parse_args()
    os.makedirs(OD, exist_ok=True)
    GR = pd.read_csv(f"{W.SETUP}/families.csv")
    M = pd.read_csv(f"{G}/members.csv")
    TG = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    scan = pd.read_csv(f"{W.SETUP}/scan_days.csv")
    align_ok = dict(zip(M.fid, M.is_ref | M.codet_prior.astype(bool) | (M.align_ratio >= 1.2)))
    tdays = {}
    for rg, tr in WS.TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        tdays[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    D, C, C2, obs, cuts, TT = flags_group(a.mf, TG)
    cuts.to_csv(f"{OD}/hf_cuts.csv", index=False); TT.to_csv(f"{OD}/tiers_by_epoch.csv", index=False)
    W.log(LOG, f"[score] group run flagged: data {len(D)} control {len(C)} control2 {len(C2)}\n{TT.to_string()}")
    Vg = score("grp", D, C, C2, obs, GR, TG, scan, tdays, align_ok)
    # singleton identity check vs T1 (data detections of singleton groups == the family's T1 detections)
    Dt, Ct, C2t, obst = flags_t1(M)
    single = GR[GR.n_members == 1]
    a_ = D[D.gid.isin(single.fid)][["gid", "date", "t", "net_mad"]].sort_values(["gid", "t"]).reset_index(drop=True)
    b_ = Dt[Dt.gid.isin(single.fid) & Dt.date.isin(set(D.date))][["gid", "date", "t", "net_mad"]].sort_values(["gid", "t"]).reset_index(drop=True)
    same = len(a_) == len(b_) and np.allclose(a_.t, b_.t) and np.allclose(a_.net_mad, b_.net_mad, atol=1e-4)
    pd.DataFrame([dict(n_grp=len(a_), n_t1=len(b_), identical=bool(same),
                       max_dmad=float(np.abs(a_.net_mad - b_.net_mad).max()) if len(a_) == len(b_) and len(a_) else np.nan)]).to_csv(
        f"{OD}/singleton_check.csv", index=False)
    W.log(LOG, f"[score] singleton check: group run {len(a_)} vs T1 {len(b_)} data detections, identical={same}")
    Vt = score("t1union", Dt, Ct, C2t, obst, GR, TG, scan, tdays, align_ok)
    V = pd.concat([Vg, Vt], ignore_index=True)
    V.to_csv(f"{OD}/by_group.csv", index=False)
    S = summarize(V)
    S.to_csv(f"{OD}/summary.csv", index=False)
    rows = []
    for run, g in V.groupby("run"):
        for tier in TIERS:
            o = g[g.n_label_days > 0]
            for s in SEP_GRID:
                rows.append(dict(run=run, tier=tier, sep=s, n_label=len(o), n_data=int((o[f"sep_{tier}"] >= s).sum()),
                                 n_ctl=int((o[f"sep_ctl_{tier}"] >= s).sum())))
    pd.DataFrame(rows).to_csv(f"{OD}/sep_curve.csv", index=False)
    W.log(LOG, "[score] summary\n" + S[(S.window == "ALL") & (S.size_class == "ALL")][
        ["run", "region", "n_groups_obs", "n_families_obs", "n_truth_test", "recall_raw_core", "recall_all_raw_core",
         "prec_raw_core", "prec_cat_core", "prec_ctl_core", "data_per_gday_core", "ctl_per_gday_core", "sep3_core_data", "sep3_core_ctl"]].to_string())


if __name__ == "__main__":
    main()
