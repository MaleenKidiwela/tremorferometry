#!/usr/bin/env python
"""DUC_FIN step 1 -- why do we find 1.4-2.8 x more catalogue events than Ducellier in May-Dec 2011?

  python scripts/duc_fin_d2011.py events    part A (catalogue files only): one-to-one matching by month, his channel
                                            counts, inferred stations missing from his catalogue, the extras' rule
                                            margin / clustering / co-occurrence  -> OUT/d2011/*.csv + OUT/D2011.md
  python scripts/duc_fin_d2011.py control   part B (needs OUT/feat from duc_fin_mf.py): our extras re-summed on HIS
                                            station set, per-station share of the extras, and the scrambled control's
                                            his-rule pass rate by month -> appended to OUT/D2011.md

Ours = results/duc_fix/mf1/<m>.parquet (pass 1: his templates only, no re-templated station), his = the 'net'
catalogue (results/ducellier_detections.parquet); his catalogue rule = cc x nchannel >= threshold_perm (60 kept
families).  Matching: greedy one-to-one, closest pairs first, within +-2 s.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
from t_common import log

OD = f"{F.OUT}/d2011"
M0, M1 = "2008-01", "2011-12"


def ours_mf1():
    P = []
    for m in pd.period_range(M0, M1, freq="M"):
        p = f"{F.FIX}/mf1/{m}.parquet"
        if os.path.exists(p): P.append(pd.read_parquet(p))
    D = pd.concat(P, ignore_index=True); D["fam"] = D.fam.astype(str)
    D["t"] = D.t.astype("datetime64[ns]"); D["m"] = D.t.dt.strftime("%Y-%m")
    return D


def his_month_nch(H):
    """his mode nchannel per family-month (all his detections) and the channel count of his full station list."""
    HR = DC.his_rows(); nc = {r.sta: len(r.chans.split(",")) for r in HR.itertuples()}
    fp = DC.families(); famst = {r.family: r.stations.split(",") for r in fp.itertuples()}
    G = H.groupby(["fam", "m"]).nch.agg(lambda x: int(x.value_counts().index[0])).rename("his_mode_nch").reset_index()
    G["max_nch"] = G.fam.map(lambda f: sum(nc[s] for s in famst[f]))
    return G, nc, famst


def infer_missing(G, nc, famst, months):
    """per month: greedy set of stations whose removal best reproduces his mode nchannel over all families
    (exact matches counted).  Returns DataFrame month, missing stations, families matched before/after."""
    rows = []
    allst = sorted(nc)
    for m in months:
        g = G[G.m == m]
        if not len(g): continue
        fam = g.fam.tolist(); obs = g.his_mode_nch.to_numpy(); mx = g.max_nch.to_numpy()
        memb = np.array([[s in famst[f] for s in allst] for f in fam])        # (nf, ns)
        w = np.array([nc[s] for s in allst])
        miss = np.zeros(len(allst), bool)
        def score(ms): return int((mx - (memb[:, ms] * w[ms]).sum(1) == obs).sum())
        best = score(miss); base = best
        while True:
            cand = [(score(miss | (np.arange(len(allst)) == j)), j) for j in range(len(allst)) if not miss[j]]
            sc, j = max(cand)
            if sc <= best: break
            miss[j] = True; best = sc
        rows.append(dict(month=m, n_fam=len(fam), fam_match_full_list=base, fam_match_after=best,
                         missing=",".join(np.array(allst)[miss]), n_missing_ch=int(w[miss].sum())))
    return pd.DataFrame(rows)


def events():
    os.makedirs(OD, exist_ok=True)
    thr = DC.thresholds()
    D = ours_mf1(); D = D[D.t < pd.Timestamp("2012-01-01")]
    H = F.his_catalogue(); H = H[(H.t >= pd.Timestamp("2008-01-01")) & (H.t < pd.Timestamp("2012-01-01"))].reset_index(drop=True)
    H["m"] = H.t.dt.strftime("%Y-%m")
    D["te"] = D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    D["hp"] = D.pass_his == 1
    log(f"[d2011] ours {len(D)} ({int(D.hp.sum())} pass his rule), his {len(H)} ({int(H.hp.sum())})")
    # --- raw (all 8 x MAD) one-to-one, same family
    D = D.reset_index(drop=True)
    ia, ib, dt = F.match_one_to_one(D.te.to_numpy(), D.fam.to_numpy(), H.t_epoch.to_numpy(), H.fam.to_numpy(), same_fam=True)
    D["m_raw"] = False; D.loc[ia, "m_raw"] = True
    D["his_raw_idx"] = -1; D.loc[ia, "his_raw_idx"] = ib
    H["m_raw"] = False; H.loc[ib, "m_raw"] = True
    # --- one-to-one, his rule both sides; same family and any family
    O = D[D.hp].reset_index(drop=True); Hh = H[H.hp].reset_index(drop=True)
    res = {}
    for tag, same in (("same", True), ("any", False)):
        ia, ib, dt = F.match_one_to_one(O.te.to_numpy(), O.fam.to_numpy(), Hh.t_epoch.to_numpy(), Hh.fam.to_numpy(), same_fam=same)
        O[f"m_{tag}"] = False; O.loc[ia, f"m_{tag}"] = True
        Hh[f"m_{tag}"] = False; Hh.loc[ib, f"m_{tag}"] = True
    rows = []
    for m in sorted(set(D.m) | set(H.m)):
        o = O[O.m == m]; h = Hh[Hh.m == m]; d = D[D.m == m]; hr = H[H.m == m]
        rows.append(dict(month=m, his_rule=len(h), ours_rule=len(o), ratio=len(o) / max(len(h), 1),
                         matched_same=int(o.m_same.sum()), his_missed_same=int((~h.m_same).sum()), our_extra_same=int((~o.m_same).sum()),
                         recall_same=h.m_same.mean() if len(h) else np.nan, precision_same=o.m_same.mean() if len(o) else np.nan,
                         matched_any=int(o.m_any.sum()), his_missed_any=int((~h.m_any).sum()), our_extra_any=int((~o.m_any).sum()),
                         recall_any=h.m_any.mean() if len(h) else np.nan, precision_any=o.m_any.mean() if len(o) else np.nan,
                         his_raw=len(hr), ours_raw=len(d), raw_recall_same=hr.m_raw.mean() if len(hr) else np.nan,
                         raw_precision_same=d.m_raw.mean() if len(d) else np.nan,
                         his_nch_med_rule=h.nch.median() if len(h) else np.nan, his_nch_med_raw=hr[hr.thr > 0].nch.median(),
                         ours_nch_med_rule=o.nchannel.median() if len(o) else np.nan, ours_nch_med_raw=d[d.thr_perm > 0].nchannel.median()))
    Mo = pd.DataFrame(rows); Mo.to_csv(f"{OD}/monthly_one_to_one.csv", index=False)
    # --- his channel counts and the stations his catalogue lacks
    G, nc, famst = his_month_nch(H[H.thr > 0])
    G.to_csv(f"{OD}/his_mode_nch_by_family_month.csv", index=False)
    mi = infer_missing(G, nc, famst, sorted(G.m.unique()))
    mi.to_csv(f"{OD}/his_missing_stations_by_month.csv", index=False)
    # --- the May-Dec 2011 extras (ours, his rule, no same-family his rule event within 2 s)
    late = O.m.between("2011-05", "2011-12"); early = O.m.between("2008-01", "2011-04")
    X = O[late & ~O.m_same].copy(); Mt = O[late & O.m_same].copy(); Xe = O[early & ~O.m_same].copy()
    Hall = H.sort_values("t_epoch"); tall = Hall.t_epoch.to_numpy()
    Hr = Hh.sort_values("t_epoch")
    for Z in (X, Mt, Xe):
        Z["margin"] = Z.ccsum / Z.thr_perm
        # his raw detection of the same family within 2 s (his detector fired; his rule failed?)
        dr = F.nearest(Z.te.to_numpy(), H.t_epoch.to_numpy(), Z.fam.to_numpy(), H.fam.to_numpy())
        Z["his_raw_same_2s"] = dr <= 2.0
        # his rule events of OTHER families within 2 s / 10 s; his rule events of any family within 1 h
        do = np.full(len(Z), np.inf)
        for f in Z.fam.unique():
            mm = (Z.fam == f).to_numpy(); r = np.sort(Hr.t_epoch.to_numpy()[Hr.fam.to_numpy() != f])
            if len(r): do[mm] = F.nearest(Z.te.to_numpy()[mm], r)
        Z["his_rule_other_2s"] = do <= 2.0; Z["his_rule_other_10s"] = do <= 10.0
        Z["his_rule_any_1h"] = F.nearest(Z.te.to_numpy(), Hr.t_epoch.to_numpy()) <= 3600.0
        # clustering: own extras within the same hour (bursts)
        hk = Z.t.dt.floor("h"); Z["n_same_hour"] = hk.map(hk.value_counts())
    # his values for the raw-matched extras
    X["his_idx"] = X.his_raw_idx.astype(int)
    ok = X.his_idx >= 0
    X["his_cc"] = np.nan; X["his_nch"] = np.nan
    X.loc[ok, "his_cc"] = H.cc.to_numpy()[X.his_idx[ok].to_numpy()]
    X.loc[ok, "his_nch"] = H.nch.to_numpy()[X.his_idx[ok].to_numpy()]
    X["his_margin"] = X.his_cc * X.his_nch / X.thr_perm
    X.to_parquet(f"{OD}/extras_may_dec_2011.parquet", index=False); Mt.to_parquet(f"{OD}/matched_may_dec_2011.parquet", index=False)
    def summ(Z, name):
        return dict(set=name, n=len(Z), med_nch=Z.nchannel.median(), med_cc=Z.cc.median(), med_margin=Z.margin.median(),
                    frac_his_raw_same_2s=Z.his_raw_same_2s.mean(), frac_his_rule_other_2s=Z.his_rule_other_2s.mean(),
                    frac_his_rule_other_10s=Z.his_rule_other_10s.mean(), frac_his_rule_any_1h=Z.his_rule_any_1h.mean(),
                    med_n_same_hour=Z.n_same_hour.median(), frac_in_hours_ge20=(Z.n_same_hour >= 20).mean(),
                    n_families=Z.fam.nunique(), top_family_share=Z.fam.value_counts(normalize=True).iloc[0] if len(Z) else np.nan)
    S = pd.DataFrame([summ(X, "extras 2011-05..12"), summ(Mt, "matched 2011-05..12"), summ(Xe, "extras 2008-01..2011-04")])
    S.to_csv(f"{OD}/extras_summary.csv", index=False)
    xs = dict(n=len(X), his_raw_found=int(ok.sum()), med_his_cc=float(X.his_cc.median()), med_our_cc=float(X[ok].cc.median()),
              med_his_nch=float(X.his_nch.median()), med_our_nch=float(X[ok].nchannel.median()),
              frac_his_margin_lt1=float((X[ok].his_margin < 1).mean()),
              frac_his_fail_due_to_nch=float(((X[ok].his_margin < 1) & (X[ok].his_cc * X[ok].nchannel >= X[ok].thr_perm)).mean()))
    json.dump(xs, open(f"{OD}/extras_his_raw.json", "w"), indent=1)
    write_md_a(Mo, mi, S, xs)
    log(f"[d2011] events done: {json.dumps(xs)}")


def verdict_a(Mo, mi, S, xs):
    late = Mo[Mo.month.between("2011-05", "2011-12")]; early = Mo[Mo.month.between("2008-01", "2011-04")]
    x, mt, xe = S.iloc[0], S.iloc[1], S.iloc[2]
    m5 = mi[mi.month.between("2011-05", "2011-12")]; m0 = mi[mi.month.between("2010-09", "2011-04")]
    import collections
    c5 = collections.Counter(",".join(m5.missing.fillna("")).split(",")); c0 = set(",".join(m0.missing.fillna("")).split(","))
    newmiss = sorted(k for k, v in c5.items() if k and v >= 3 and k not in c0)
    return "\n".join([
        f"* **His catalogue lost stations from May 2011; ours did not.** His median channel count for his rule events fell from "
        f"{early.his_nch_med_rule.median():.0f} (2008-01..2011-04, monthly median) to {late.his_nch_med_rule.median():.0f} (May-Dec 2011), "
        f"and for all his detections from {early.his_nch_med_raw.median():.0f} to {late.his_nch_med_raw.median():.0f}; ours stayed at "
        f"{early.ours_nch_med_rule.median():.0f} -> {late.ours_nch_med_rule.median():.0f}. With his full station list, his channel count "
        f"is explained exactly for {m0.fam_match_full_list.median():.0f} of 60 families per month before May 2011 and for "
        f"{m5.fam_match_full_list.median():.0f} of 60 after; removing a few stations explains {m5.fam_match_after.median():.0f} of 60. "
        f"Stations inferred missing in >= 3 of the May-Dec 2011 months but in none of 2010-09..2011-04: {', '.join(newmiss)} "
        f"(BK GASB / WDC / YBH changed location code to 00 on 2011-06-16 / 05-06 / 06-03); in Nov-Dec 2011 many NC analogue "
        f"stations follow (their SHZ -> EHZ switch).",
        f"* **Fewer channels make his catalogue rule harder to pass** (the rule sums the correlation over channels). Of our "
        f"{xs['n']} extra events, {xs['his_raw_found']} are also in his RAW detections of the same family within 2 s with the same "
        f"correlation (his {xs['med_his_cc']:.3f} vs ours {xs['med_our_cc']:.3f}) but {xs['med_his_nch']:.0f} instead of "
        f"{xs['med_our_nch']:.0f} channels: {xs['frac_his_fail_due_to_nch']:.0%} of them would pass his rule with his correlation and our channel count.",
        f"* **The extras look like ordinary repeats, not an artefact.** They spread over {int(x.n_families)} families (largest share "
        f"{x.top_family_share:.0%}), {x.frac_his_rule_any_1h:.0%} fall within 1 h of one of his catalogued events (matched events: "
        f"{mt.frac_his_rule_any_1h:.0%}; extras before May 2011: {xe.frac_his_rule_any_1h:.0%}), median {x.med_n_same_hour:.0f} per hour "
        f"(matched: {mt.med_n_same_hour:.0f}); only {x.frac_in_hours_ge20:.0%} sit in hours with >= 20 extras (matched {mt.frac_in_hours_ge20:.0%}). "
        f"Their rule margin (summed cc / his threshold) is {x.med_margin:.2f} (matched {mt.med_margin:.2f}): weaker events near his threshold.",
        f"* **Before May 2011 the match is one-to-one:** same-family recall {early.matched_same.sum() / early.his_rule.sum():.3f}, precision "
        f"{early.matched_same.sum() / early.ours_rule.sum():.3f} over 2008-01..2011-04 (his rule both sides, +-2 s).",
        "* **Decision (part A): his data gap, not our artefact.** Nothing is re-scanned. May-Dec 2011 is kept, evaluated separately "
        "in the confidence step (his catalogue there is incomplete) and never used to fit. Part B below re-sums our events on his "
        "remaining stations and adds the scrambled control."])


def write_md_a(Mo, mi, S, xs):
    late = Mo[Mo.month.between("2011-05", "2011-12")]; early = Mo[Mo.month.between("2008-01", "2011-04")]
    mi2 = mi[mi.month.between("2010-09", "2011-12")]
    L = ["# DUC_FIN step 1: the 2011 excess (May-Dec 2011 ours 1.4-2.8 x his under his rule)", "",
         f"_Auto-written by `scripts/duc_fin_d2011.py`, {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%M} UTC. Files in `{OD}/`._", "",
         "## Plain-language verdict (part A, catalogue files only)", "",
         verdict_a(Mo, mi, S, xs), "",
         "## 1. One-to-one matching by month (his rule both sides, +-2 s) -- `monthly_one_to_one.csv`", "",
         f"2008-01..2011-04: his {int(early.his_rule.sum())}, ours {int(early.ours_rule.sum())}, same-family recall "
         f"{early.matched_same.sum() / early.his_rule.sum():.3f}, precision {early.matched_same.sum() / early.ours_rule.sum():.3f}.",
         f"2011-05..12: his {int(late.his_rule.sum())}, ours {int(late.ours_rule.sum())}, same-family recall "
         f"{late.matched_same.sum() / max(late.his_rule.sum(), 1):.3f}, precision {late.matched_same.sum() / max(late.ours_rule.sum(), 1):.3f}; "
         f"(a) his events we miss {int(late.his_missed_same.sum())}, (b) our extra events {int(late.our_extra_same.sum())} "
         f"(any family: {int(late.our_extra_any.sum())}).", "",
         Mo[Mo.month >= "2010-06"][["month", "his_rule", "ours_rule", "ratio", "matched_same", "his_missed_same", "our_extra_same",
                                    "recall_same", "precision_same", "our_extra_any", "his_nch_med_rule", "his_nch_med_raw",
                                    "ours_nch_med_rule", "ours_nch_med_raw"]].round(3).to_markdown(index=False), "",
         "## 2. His channel counts: which stations his catalogue lost -- `his_missing_stations_by_month.csv`", "",
         "For every family and month, his most common nchannel is compared with the channel count of his own station list "
         "for that family; the stations whose removal reproduces the most families exactly are found greedily.", "",
         mi2.to_markdown(index=False), "",
         "## 3. The May-Dec 2011 extras -- `extras_summary.csv`, `extras_his_raw.json`", "",
         S.round(3).to_markdown(index=False), "",
         f"Of the {xs['n']} extras, {xs['his_raw_found']} are ALSO in his raw (8 x MAD) catalogue of the same family within 2 s, "
         f"i.e. his detector fired but his catalogue rule failed. For those: his cc {xs['med_his_cc']:.3f} vs ours "
         f"{xs['med_our_cc']:.3f} (median), his nchannel {xs['med_his_nch']:.0f} vs ours {xs['med_our_nch']:.0f}; "
         f"{xs['frac_his_margin_lt1']:.1%} fail his rule with his numbers, and "
         f"{xs['frac_his_fail_due_to_nch']:.1%} would pass with his cc and OUR channel count.", ""]
    open(f"{F.OUT}/D2011.md", "w").write("\n".join(L) + "\n")


def control():
    """part B: needs OUT/feat/<date>.parquet (duc_fin_mf.py) for 2008-01..2011-12."""
    import glob
    X = pd.read_parquet(f"{OD}/extras_may_dec_2011.parquet"); X["fam"] = X.fam.astype(str)
    mi = pd.read_csv(f"{OD}/his_missing_stations_by_month.csv")
    miss = {r.month: set(str(r.missing).split(",")) - {"nan", ""} for r in mi.itertuples()}
    no = F.netof()
    rows = []; ctl = []; share = []
    DT = F.day_table()
    for d in DT[(DT.date >= "2008-01-01") & (DT.date < "2012-01-01")].date:
        p = f"{F.OUT}/feat/{d}.parquet"
        if not os.path.exists(p): continue
        A = pd.read_parquet(p); A["fam"] = A.fam.astype(str); m = d[:7]
        c = A[A.series == "c1"]
        ctl.append(dict(month=m, date=d, n_c1=len(c), n_c1_pass=int(c.pass_his.sum()),
                        n_data=int((A.series == "data").sum()), n_data_pass=int(A[A.series == "data"].pass_his.sum())))
        if not ("2011-05" <= m <= "2011-12"): continue
        dd = A[(A.series == "data") & (A.pass_his == 1)]
        if not len(dd): continue
        key = np.round(dd.t_epoch.to_numpy() * 20).astype(np.int64)
        xk = X[X.t.dt.strftime("%Y-%m-%d") == d]
        xkey = set(zip(xk.fam, np.round(xk.te.to_numpy() * 20).astype(np.int64)))
        isx = np.array([(f, k) in xkey for f, k in zip(dd.fam, key)])
        stc = [c_[3:] for c_ in dd.columns if c_.startswith("cc_")]
        ms = {f"{no.get(s, 'NC')}.{s}" if "." not in s else s for s in miss.get(m, set())}
        for j, (_, r) in enumerate(dd.iterrows()):
            tot = 0.0; his = 0.0; per = {}
            for s in stc:
                cc_ = r.get(f"cc_{s}"); n_ = r.get(f"nc_{s}")
                if not (np.isfinite(cc_) and np.isfinite(n_) and n_ > 0): continue
                v = cc_ * n_; tot += v; per[s] = v
                if s not in ms: his += v
            nh = sum(int(r[f"nc_{s}"]) for s in stc if np.isfinite(r.get(f"nc_{s}", np.nan)) and r.get(f"nc_{s}", 0) > 0 and s not in ms)
            rows.append(dict(date=d, month=m, fam=r.fam, extra=bool(isx[j]), ccsum=r.ccsum, thr=r.thr_perm, ccsum_his_set=his,
                             nch_his_set=nh, pass_his_set=bool(r.thr_perm > 0 and his >= r.thr_perm),
                             share_missing=1 - his / tot if tot > 0 else np.nan))
            if isx[j]:
                for s, v in per.items(): share.append(dict(netsta=s, v=v, tot=tot))
    R = pd.DataFrame(rows); Cc = pd.DataFrame(ctl); Sh = pd.DataFrame(share)
    R.to_parquet(f"{OD}/rule_events_his_station_set.parquet", index=False); Cc.to_csv(f"{OD}/control_pass_by_day.csv", index=False)
    cm = Cc.groupby("month")[["n_c1", "n_c1_pass", "n_data", "n_data_pass"]].sum().reset_index()
    cm["c1_pass_over_data_pass"] = cm.n_c1_pass / cm.n_data_pass.clip(lower=1)
    cm.to_csv(f"{OD}/control_pass_by_month.csv", index=False)
    st = Sh.groupby("netsta").v.sum().div(Sh.drop_duplicates(["tot"]).tot.sum() if len(Sh) else 1).sort_values(ascending=False).rename("share_of_extra_ccsum")
    st.to_csv(f"{OD}/extras_station_share.csv")
    g = R.groupby("extra").agg(n=("ccsum", "size"), frac_pass_his_set=("pass_his_set", "mean"), med_share_missing=("share_missing", "median"),
                               med_nch_his_set=("nch_his_set", "median")).reset_index()
    g.to_csv(f"{OD}/extras_his_station_set.csv", index=False)
    pre = cm[cm.month < "2011-05"]; post = cm[cm.month >= "2011-05"]
    L = ["", "## 4. Part B (feature re-run): extras on his station set, station shares, scrambled control", "",
         "Each our rule-passing event of May-Dec 2011 re-summed over only the stations his catalogue still had that month "
         "(section 2), at the same time and with our per-station correlations -- `extras_his_station_set.csv`:", "",
         g.round(3).to_markdown(index=False), "",
         "Share of the extras' summed correlation carried by each station (top 12) -- `extras_station_share.csv`:", "",
         st.head(12).round(3).to_frame().to_markdown(), "",
         "Scrambled control (same detector, each site's correlation shifted to a different time of the same day): "
         "events passing his rule, control / data -- `control_pass_by_month.csv`:", "",
         f"2008-01..2011-04: control {int(pre.n_c1_pass.sum())} vs data {int(pre.n_data_pass.sum())} "
         f"({pre.n_c1_pass.sum() / max(pre.n_data_pass.sum(), 1):.4f}); 2011-05..12: control {int(post.n_c1_pass.sum())} vs data "
         f"{int(post.n_data_pass.sum())} ({post.n_c1_pass.sum() / max(post.n_data_pass.sum(), 1):.4f}).", "",
         cm[cm.month >= "2010-09"].round(4).to_markdown(index=False), ""]
    open(f"{F.OUT}/D2011.md", "a").write("\n".join(L) + "\n")
    log(f"[d2011] control part done: {g.to_dict('records')}; control pass ratio pre {pre.n_c1_pass.sum() / max(pre.n_data_pass.sum(), 1):.4f} "
        f"post {post.n_c1_pass.sum() / max(post.n_data_pass.sum(), 1):.4f}")


if __name__ == "__main__":
    {"events": events, "control": control}[sys.argv[1]]()
