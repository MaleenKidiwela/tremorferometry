#!/usr/bin/env python
"""E5 step 5 STEP E (+ the STEP G catalog) -- score the from-scratch 2012 network matched filter against Bostock's TEST-half events.
Pre-registered bars (results/merlin_audit_e4_2026-09-11.md, prescription 5) are REPORTED as met / not met; nothing is tuned.

Inputs (tag = the e4_template_test.py run): results/e5r_detections_<tag>.parquet (every declustered detection >= 7 MAD),
results/e5_budget_<tag>.csv (per day/L/family/T counts incl. the rolled control), results/e4_tt_events_<tag>.parquet (per Bostock event:
own-family network MAD at +-1 / +-2 s, rolled), results/e4_tt_station_<tag>.csv, results/e5r_gt_bostock2012.parquet, e5r_families_2012.csv,
e5r_template_qc_2012.csv, e5_prediction_2012_recall.csv.
Definitions: quiet day = 2012-08-15..08-29 (doy 228-242; 0 Bostock events); ETS day = >= 100 Bostock events that UTC day; FA per template
per day = mean over (family, quiet day) of the declustered detections of that family's network; operating threshold = 8 MAD unless the
quiet-day FA at 8 exceeds 1/template/day (then both 8 and the smallest T in {9, 10} meeting it are reported).  Cross-family declustering
(the catalog): a detection is kept iff its network MAD is the maximum of all detections of ANY family within +-2 s (Bostock-style max-cc
assignment; MAD-normalised so families with different station counts are comparable).  Recall = a Bostock TEST-half event has (own) a
detection of its OWN family >= T within +-2 s (like-for-like with the E5 prediction), (cat) a kept catalog detection of any family within
+-2 s.  Chance = the same with the detection times shifted by +600 s, and the rolled-network score at the events.
Outputs results/e5r_score_<tag>.json, e5r_score_<tag>_{fa,recall_mag,recall_nsta,recall_family,daily,extras}.csv,
        results/e5_catalog_2012_sVI[_L12].parquet (tag e5r_snr only; the primary product)."""
from __future__ import annotations
import argparse, json, os
import numpy as np, pandas as pd
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
MAGBINS = [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]
NSTABINS = [(1, 4, "1-4"), (5, 7, "5-7"), (8, 10, "8-10"), (11, 99, "11-15")]
QUIET = set(range(228, 243))


def min_absdiff(a, b_sorted):
    a = np.asarray(a, float)
    if len(a) == 0:
        return np.zeros(0)
    if len(b_sorted) == 0:
        return np.full(len(a), np.inf)
    j = np.searchsorted(b_sorted, a)
    lo = b_sorted[np.clip(j - 1, 0, len(b_sorted) - 1)]; hi = b_sorted[np.clip(j, 0, len(b_sorted) - 1)]
    return np.minimum(np.abs(a - lo), np.abs(a - hi))


def nearest_idx(a, b_sorted, tol):
    """index into b_sorted of the nearest element within tol, else -1."""
    a = np.asarray(a, float); j = np.searchsorted(b_sorted, a)
    lo = np.clip(j - 1, 0, len(b_sorted) - 1); hi = np.clip(j, 0, len(b_sorted) - 1)
    pick = np.where(np.abs(a - b_sorted[lo]) <= np.abs(a - b_sorted[hi]), lo, hi)
    ok = np.abs(a - b_sorted[pick]) <= tol
    return np.where(ok, pick, -1)


def local_max_keep(t, v, w=2.0):
    o = np.argsort(t, kind="stable"); ts = t[o]; vs = v[o]; n = len(ts)
    lo = np.searchsorted(ts, ts - w, "left"); hi = np.searchsorted(ts, ts + w, "right")
    keep = np.fromiter((vs[i] >= vs[lo[i]:hi[i]].max() for i in range(n)), bool, n)
    out = np.zeros(n, bool); out[o] = keep
    return out


def local_max_keep_cc(t, v, fam, share_ok, w=2.0):
    """E6 (Merlin E5r-audit F8): cross-family assignment by RAW weighted cc among families that SHARE >= 50 % of their stations
    (share_ok[fa, fb] bool matrix, bkey-indexed); a detection is dropped iff a detection of a station-sharing family within +-w s
    has a higher cc (ties: the earlier one wins).  Families on disjoint station sets never suppress one another."""
    o = np.argsort(t, kind="stable"); ts = t[o]; vs = v[o]; fs = fam[o]; n = len(ts)
    lo = np.searchsorted(ts, ts - w, "left"); hi = np.searchsorted(ts, ts + w, "right")
    keep = np.ones(n, bool)
    for i in range(n):
        fi = fs[i]; vi = vs[i]
        for j in range(lo[i], hi[i]):
            if j == i or not share_ok[fi, fs[j]]:
                continue
            if vs[j] > vi or (vs[j] == vi and j < i):
                keep[i] = False
                break
    out = np.zeros(n, bool); out[o] = keep
    return out


def doy_of(t):
    return pd.DatetimeIndex(pd.to_datetime(np.asarray(t, float), unit="s")).dayofyear.to_numpy()


def bins_of(x, bins):
    lab = np.array([None] * len(x), dtype=object)
    for lo, hi, name in bins:
        lab[(x >= lo) & (x < hi)] = name
    return lab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="e5r_snr")
    ap.add_argument("--write-catalog", action="store_true")
    ap.add_argument("--rng", type=int, default=3)
    # E6 (2026-09-11, Merlin E5r-audit F3/F8, prescription 4i-ii) -- OFF by default (E5r scoring bit-for-bit unchanged):
    #   --e6            outputs go to results/e6_score_<tag>*; operating thresholds {7.5, 8, 8.5, 9, 10}; the "rolled <= 5 %" bar uses
    #                   the SITE-rolled control (budget column n_det_roll_site; station-rolled kept as a lower bound); catalog recall is
    #                   reported CHANCE-CORRECTED, (cat - chance) / (1 - chance), and met_cat is set from the corrected value; own-family
    #                   recall is the primary number; cross-family assignment by RAW weighted cc among families sharing >= 50 % of their
    #                   stations (local_max_keep_cc) instead of network MAD; the template-less families' catalog "recall" is reported under
    #                   both assignments.  --qc-csv / --stations define the per-family station sets (lever runs use other sets).
    ap.add_argument("--e6", action="store_true")
    ap.add_argument("--t-ops", default="7.5,8,8.5,9,10")
    ap.add_argument("--qc-csv", default=f"{RES}/e5r_template_qc_2012.csv")
    ap.add_argument("--stations", default="")
    ap.add_argument("--write-T", type=float, default=8.0, help="E6: threshold of the catalog written with --write-catalog")
    a = ap.parse_args()
    rng = np.random.default_rng(a.rng)
    PFX = "e6_score" if a.e6 else "e5r_score"
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size()
    ets_days = set(per_day[per_day >= 100].index)
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv").set_index("bkey")
    qc = pd.read_csv(a.qc_csv if a.e6 else f"{RES}/e5r_template_qc_2012.csv")
    if a.e6 and a.stations:
        qc = qc[qc.sta.isin(a.stations.split(","))]
    nsta_fam = qc[qc.usable & qc.template].groupby("bkey").size()
    share_ok = None
    if a.e6:
        # E6 F8: station sets per family -> families "share" stations when |A & B| / min(|A|, |B|) >= 0.5 (bkey-indexed bool matrix)
        sets = qc[qc.usable & qc.template].groupby("bkey").sta.apply(set).to_dict()
        nb = int(max(sets) + 1) if sets else 1
        share_ok = np.zeros((nb + 1, nb + 1), bool)
        for fa_, A_ in sets.items():
            for fb_, B_ in sets.items():
                share_ok[fa_, fb_] = len(A_ & B_) / max(min(len(A_), len(B_)), 1) >= 0.5
    det = pd.read_parquet(f"{RES}/e5r_detections_{a.tag}.parquet")
    bud = pd.read_csv(f"{RES}/e5_budget_{a.tag}.csv")
    ev = pd.read_parquet(f"{RES}/e4_tt_events_{a.tag}.parquet")
    st = pd.read_csv(f"{RES}/e4_tt_station_{a.tag}.csv")
    pred = pd.read_csv(f"{RES}/e5_prediction_2012_recall.csv")
    pred = pred[(pred.set == "set15") & (pred["T"] == 8.0)]
    test = gt[(gt.half == "test") & (~gt.guard)].copy()
    test["nsta"] = test.bkey.map(nsta_fam).fillna(0).astype(int)
    test["mbin"] = bins_of(test.mag.to_numpy(), MAGBINS)
    test["nbin"] = bins_of(test.nsta.to_numpy(), NSTABINS)
    t_all_sorted = np.sort(gt.t_origin.to_numpy())
    out = dict(tag=a.tag, n_test_events=int(len(test)), n_test_templated=int((test.nsta > 0).sum()), n_families_templated=int(len(nsta_fam)),
               ets_days=sorted(int(d) for d in ets_days), quiet_days=sorted(QUIET))
    fa_rows, mag_rows, nsta_rows, fam_rows, daily_rows, ex_rows = [], [], [], [], [], []
    for L in sorted(det.L.unique()):
        bL = bud[bud.L == L]
        # ---------------- quiet-day FA vs threshold (per template per day) and the operating threshold
        q = bL[bL.doy.isin(QUIET)]; e_ = bL[bL.doy.isin(ets_days)]
        fa_T = {}
        for T in sorted(bL["T"].unique()):
            qt = q[q["T"] == T]; et = e_[e_["T"] == T]
            fa_T[T] = float(qt.n_det.mean())
            fa_rows.append(dict(L=L, T=T, quiet_fa_per_tpl_day=float(qt.n_det.mean()), quiet_fa_p90=float(qt.n_det.quantile(.9)),
                                quiet_rolled_per_tpl_day=float(qt.n_det_roll.mean()), n_quiet_family_days=int(len(qt)),
                                ets_det_per_tpl_day=float(et.n_det.mean()), ets_rolled_per_tpl_day=float(et.n_det_roll.mean()),
                                ets_rolled_frac=float(et.n_det_roll.sum() / max(et.n_det.sum(), 1)),
                                ets_precision_own=float(et.n_match_own.sum() / max(et.n_det.sum(), 1)),
                                ets_precision_any=float(et.n_match_any.sum() / max(et.n_det.sum(), 1)),
                                all_rolled_frac=float(bL[bL["T"] == T].n_det_roll.sum() / max(bL[bL["T"] == T].n_det.sum(), 1))))
        # 8 MAD is the PRE-REGISTERED operating point; 9 and 10 are reported alongside as diagnostics (never as a re-tuning).
        T_ops = [float(x) for x in a.t_ops.split(",")] if a.e6 else [8.0, 9.0, 10.0]
        T_ops = [T for T in T_ops if T in set(bL["T"].unique())]
        T_meet = [T for T in T_ops if fa_T.get(T, np.inf) <= 1.0]
        out[f"L{int(L)}_smallest_T_with_quiet_fa_le1"] = T_meet[0] if T_meet else None
        out[f"L{int(L)}_quiet_fa"] = {str(k): v for k, v in fa_T.items()}
        out[f"L{int(L)}_T_ops"] = T_ops
        dL = det[det.L == L]
        evL = ev[(ev.L == L) & (ev.half == "test") & (~ev.guard)]
        for T in T_ops:
            key = f"L{int(L)}_T{T:g}" if a.e6 else f"L{int(L)}_T{int(T)}"     # E6: half-integer thresholds must not collide (int(8.5) == 8)
            D = dL[dL.net_mad >= T].copy()
            keep_mad = local_max_keep(D.t.to_numpy(), D.net_mad.to_numpy(), 2.0)
            if a.e6:
                keep = local_max_keep_cc(D.t.to_numpy(), D.net_cc.to_numpy().astype(float), D.fam.to_numpy().astype(int), share_ok, 2.0)
                # template-less families' catalog "recall" under the OLD (max network MAD) assignment, for the F8 comparison
                Cm = D[keep_mad].sort_values("t"); tCm = np.sort(Cm.t.to_numpy())
                tl = test[test.nsta == 0]
                tl_old = float((min_absdiff(tl.t_origin.to_numpy(), tCm) <= 2.0).mean()) if len(tl) else np.nan
            else:
                keep = keep_mad
            C = D[keep].copy().sort_values("t").reset_index(drop=True)
            tC = C.t.to_numpy()
            # matched Bostock event (any family) within +-2 s
            j = nearest_idx(tC, t_all_sorted, 2.0)
            order_gt = np.argsort(gt.t_origin.to_numpy(), kind="stable")          # E6: file row of the matched event (gt is row-indexed 0..N-1)
            gts = gt.iloc[order_gt].reset_index(drop=True)
            C["matched_bostock_idx"] = np.where(j >= 0, gts.index.to_numpy()[np.clip(j, 0, len(gts) - 1)], -1)
            C["matched_bostock_id"] = np.where(j >= 0, order_gt[np.clip(j, 0, len(gts) - 1)], -1)
            C["matched_bostock_fam"] = np.where(j >= 0, gts.bkey.to_numpy()[np.clip(j, 0, len(gts) - 1)], -1)
            C["matched_bostock_mag"] = np.where(j >= 0, gts.mag.to_numpy()[np.clip(j, 0, len(gts) - 1)], np.nan)
            C["matched_same_family"] = (C.matched_bostock_fam == C.fam.astype(int))
            C["doy"] = doy_of(tC)
            # ---------------- recall on TEST-half events
            own_hit2 = np.zeros(len(test), bool); own_hit5 = np.zeros(len(test), bool)
            for fk, g in test.groupby("bkey"):
                td = np.sort(dL[(dL.fam == float(fk)) & (dL.net_mad >= T)].t.to_numpy())
                dd = min_absdiff(g.t_origin.to_numpy(), td)
                own_hit2[g.index.to_numpy() if False else test.index.get_indexer(g.index)] = dd <= 2.0
                own_hit5[test.index.get_indexer(g.index)] = dd <= 5.0
            dcat = min_absdiff(test.t_origin.to_numpy(), tC)
            cat_hit2 = dcat <= 2.0; cat_hit5 = dcat <= 5.0
            dsh = min_absdiff(test.t_origin.to_numpy(), np.sort(tC + 600.0))
            chance_cat2 = dsh <= 2.0
            evm = evL.set_index(["fam", "t_origin"]) if len(evL) else None
            # own-family rolled score at the events (chance for the own-family recall)
            roll_hit = evL[evL.net_mad_shift >= T].shape[0] / max(len(evL), 1)
            test["own2"], test["own5"], test["cat2"], test["cat5"], test["ch2"] = own_hit2, own_hit5, cat_hit2, cat_hit5, chance_cat2

            def rec(sub):
                return dict(n=int(len(sub)), recall_own_2s=float(sub.own2.mean()) if len(sub) else np.nan, recall_own_5s=float(sub.own5.mean()) if len(sub) else np.nan,
                            recall_cat_2s=float(sub.cat2.mean()) if len(sub) else np.nan, recall_cat_5s=float(sub.cat5.mean()) if len(sub) else np.nan,
                            chance_cat_2s=float(sub.ch2.mean()) if len(sub) else np.nan)
            R_all = rec(test); R_tpl = rec(test[test.nsta > 0])
            R_all.update(L=L, T=T, mbin="all"); mag_rows.append(R_all)
            r2 = rec(test[test.nsta > 0]); r2.update(L=L, T=T, mbin="all_templated_families"); mag_rows.append(r2)
            for lo, hi, name in MAGBINS:
                r = rec(test[test.mbin == name]); r.update(L=L, T=T, mbin=name); mag_rows.append(r)
            for lo, hi, name in NSTABINS:
                r = rec(test[test.nbin == name]); r.update(L=L, T=T, nbin=name); nsta_rows.append(r)
            r = rec(test[test.nsta == 0]); r.update(L=L, T=T, nbin="0 (no template)"); nsta_rows.append(r)
            for fk, g in test.groupby("bkey"):
                r = rec(g); r.update(L=L, T=T, bkey=int(fk), nsta=int(nsta_fam.get(fk, 0)), n_build=int(F.loc[fk, "n_build"]), mag_med=float(g.mag.median()))
                fam_rows.append(r)
            # ---------------- rolled / quiet / extras
            e_T = e_[e_["T"] == T]; q_T = q[q["T"] == T]
            rolled_frac_ets_sta = float(e_T.n_det_roll.sum() / max(e_T.n_det.sum(), 1))
            # E6: the "rolled <= 5 %" bar uses the SITE-rolled control when the run wrote it (F3: the per-station roll understates chance)
            rolled_frac_ets = float(e_T.n_det_roll_site.sum() / max(e_T.n_det.sum(), 1)) if (a.e6 and "n_det_roll_site" in e_T.columns) else rolled_frac_ets_sta
            quiet_fa = float(q_T.n_det.mean())
            quiet_roll_site = float(q_T.n_det_roll_site.mean()) if "n_det_roll_site" in q_T.columns else np.nan
            post_T = bL[(bL["T"] == T) & bL.doy.isin(range(270, 280))]
            post_fa = float(post_T.n_det.mean()); post_roll_site = float(post_T.n_det_roll_site.mean()) if "n_det_roll_site" in post_T.columns else np.nan
            C["is_extra"] = C.matched_bostock_idx < 0
            X = C[C.is_extra]; Mt = C[~C.is_extra]
            frac_extra_ets = float(X.doy.isin(ets_days).mean()) if len(X) else np.nan
            amp_x, amp_m = float(X.amp.median()) if len(X) else np.nan, float(Mt.amp.median()) if len(Mt) else np.nan
            amp_m_own = float(Mt[Mt.matched_same_family].amp.median()) if len(Mt) else np.nan
            # extras within 60 s of a Bostock event of a family OTHER than the assigned one, vs random times on the same days
            near60 = np.zeros(len(X), bool); rnd60 = np.zeros(len(X), bool)
            if len(X):
                for fk, g in X.groupby("fam"):
                    toth = np.sort(gt[gt.bkey != int(fk)].t_origin.to_numpy())
                    near60[X.index.get_indexer(g.index)] = min_absdiff(g.t.to_numpy(), toth) <= 60.0
                    day0 = np.floor(g.t.to_numpy() / 86400.0) * 86400.0
                    rt = day0 + rng.uniform(0, 86400.0, len(g))
                    rnd60[X.index.get_indexer(g.index)] = min_absdiff(rt, toth) <= 60.0
            # earthquake contamination
            eqn = 0; eqc = np.nan
            eqs = []
            for fn in ("eq_catalog_pgc.csv", "eq_catalog_b011.csv"):
                p = f"/home/jovyan/tremorferometry/data/{fn}"
                if os.path.exists(p):
                    e = pd.read_csv(p); e["t"] = pd.to_datetime(e.ot).dt.tz_localize(None).astype("datetime64[ns]").astype("int64") / 1e9
                    eqs.append(e[(e.t >= gt.t_origin.min() - 15 * 86400) & (e.t <= gt.t_origin.max() + 15 * 86400)])
            if eqs:
                eq_t = np.sort(pd.concat(eqs).drop_duplicates("t").t.to_numpy())
                eqn = int((min_absdiff(tC, eq_t) <= 60.0).sum())
                eqc = float((min_absdiff(np.floor(tC / 86400) * 86400 + rng.uniform(0, 86400, len(tC)), eq_t) <= 60.0).mean()) if len(tC) else np.nan
            ex_rows.append(dict(L=L, T=T, n_catalog=int(len(C)), n_matched_any=int((~C.is_extra).sum()), n_matched_own=int(C.matched_same_family.sum()),
                                n_extras=int(len(X)), frac_extras_on_ets_days=frac_extra_ets, median_amp_extras=amp_x, median_amp_matched=amp_m,
                                median_amp_matched_own=amp_m_own, extras_near60_other_family=float(near60.mean()) if len(X) else np.nan,
                                random_near60_other_family=float(rnd60.mean()) if len(X) else np.nan, rolled_frac_ets_days=rolled_frac_ets,
                                quiet_fa_per_tpl_day=quiet_fa, precision_any=float((~C.is_extra).mean()) if len(C) else np.nan,
                                precision_own=float(C.matched_same_family.mean()) if len(C) else np.nan, n_eq_within60=eqn, eq_chance_frac=eqc,
                                n_eq=int(len(eq_t)) if eqs else 0, own_rolled_recall_at_events=float(roll_hit)))
            # ---------------- daily counts ours vs his
            dc = C.groupby("doy").size(); dm = C[~C.is_extra].groupby("doy").size()
            for d in range(228, 280):
                daily_rows.append(dict(L=L, T=T, doy=d, ours=int(dc.get(d, 0)), ours_matched=int(dm.get(d, 0)), bostock=int(per_day.get(d, 0)),
                                       rolled=float(bL[(bL.doy == d) & (bL["T"] == T)].n_det_roll.sum())))
            # ---------------- PASS / FAIL (pre-registered) + prediction ratio
            mb = {r["mbin"]: r for r in mag_rows if r["L"] == L and r["T"] == T and "mbin" in r}
            bars = dict(
                recall_overall=dict(bar=">= 0.25", own=mb["all"]["recall_own_2s"], cat=mb["all"]["recall_cat_2s"], chance_cat=mb["all"]["chance_cat_2s"]),
                recall_Mge2=dict(bar=">= 0.60", own=mb["M>=2.0"]["recall_own_2s"], cat=mb["M>=2.0"]["recall_cat_2s"]),
                recall_M18_20=dict(bar=">= 0.35", own=mb["M1.8-2.0"]["recall_own_2s"], cat=mb["M1.8-2.0"]["recall_cat_2s"]),
                recall_M16_18=dict(bar=">= 0.20", own=mb["M1.6-1.8"]["recall_own_2s"], cat=mb["M1.6-1.8"]["recall_cat_2s"]),
                rolled_frac_ets=dict(bar="<= 0.05", value=rolled_frac_ets),
                quiet_fa=dict(bar="<= 1 /template/day", value=quiet_fa),
                extras_on_ets_days=dict(bar=">= 0.90 with lower median amplitude", value=frac_extra_ets, amp_extras=amp_x, amp_matched=amp_m))
            thr = dict(recall_overall=0.25, recall_Mge2=0.60, recall_M18_20=0.35, recall_M16_18=0.20)
            binname = dict(recall_overall="all", recall_Mge2="M>=2.0", recall_M18_20="M1.8-2.0", recall_M16_18="M1.6-1.8")
            for k, v in thr.items():
                bars[k]["met_own"] = bool(bars[k]["own"] >= v)
                if a.e6:
                    # E6 F3: chance-corrected catalog recall (cat - chance) / (1 - chance); met_cat from the CORRECTED value; own is primary
                    ch = mb[binname[k]]["chance_cat_2s"]; cc_ = (bars[k]["cat"] - ch) / max(1.0 - ch, 1e-9)
                    bars[k]["chance_cat"] = ch; bars[k]["cat_corrected"] = float(cc_); bars[k]["met_cat"] = bool(cc_ >= v); bars[k]["primary"] = "own"
                else:
                    bars[k]["met_cat"] = bool(bars[k]["cat"] >= v)
            bars["rolled_frac_ets"]["met"] = bool(rolled_frac_ets <= 0.05)
            bars["rolled_frac_ets"]["control"] = "site-rolled" if (a.e6 and "n_det_roll_site" in e_T.columns) else "station-rolled"
            bars["rolled_frac_ets"]["station_rolled_value"] = rolled_frac_ets_sta
            bars["quiet_fa"]["met"] = bool(quiet_fa <= 1.0)
            if a.e6:
                bars["quiet_fa"]["site_rolled"] = quiet_roll_site; bars["quiet_fa"]["site_over_data"] = quiet_roll_site / max(quiet_fa, 1e-9)
                bars["post_sse_fa"] = dict(bar="<= 1 /template/day (09-26..10-05)", value=post_fa, site_rolled=post_roll_site, site_over_data=post_roll_site / max(post_fa, 1e-9), met=bool(post_fa <= 1.0))
                bars["templateless_cat_recall"] = dict(n_events=int((test.nsta == 0).sum()), mad_assignment=tl_old,
                                                       cc_assignment=float(test[test.nsta == 0].cat2.mean()) if (test.nsta == 0).any() else np.nan)
            bars["extras_on_ets_days"]["met"] = bool(np.isfinite(frac_extra_ets) and frac_extra_ets >= 0.90 and amp_x < amp_m)
            pr = {}
            for mode in sorted(pred["mode"].unique()):
                for alpha in (0.38, 0.2):
                    rowp = pred[(pred.L == L) & (pred["mode"] == mode) & (np.isclose(pred.alpha, alpha)) & (pred.g == "g0.86")]
                    if not len(rowp):
                        continue
                    rowp = rowp.iloc[0]
                    pr[f"{mode}_alpha{alpha}"] = dict(pred_all=float(rowp.recall_all), meas_own_all=mb["all"]["recall_own_2s"], ratio_all=mb["all"]["recall_own_2s"] / rowp.recall_all,
                                                     **{f"ratio_{nm}": mb[nm]["recall_own_2s"] / float(rowp[f"recall_{nm}"]) for _, _, nm in MAGBINS},
                                                     **{f"pred_{nm}": float(rowp[f"recall_{nm}"]) for _, _, nm in MAGBINS})
            out[key] = dict(bars=bars, prediction=pr, n_catalog=int(len(C)), n_own_hits=int(own_hit2.sum()), n_cat_hits=int(cat_hit2.sum()))
            print(f"\n=== {key}: catalog {len(C)} detections ({int((~C.is_extra).sum())} matched any, {int(C.matched_same_family.sum())} own); "
                  f"recall own/cat +-2 s {mb['all']['recall_own_2s']:.3f}/{mb['all']['recall_cat_2s']:.3f} (chance {mb['all']['chance_cat_2s']:.3f}, rolled-at-events {roll_hit:.3f}); "
                  f"by mag own: " + ", ".join(f"{nm} {mb[nm]['recall_own_2s']:.3f}" for _, _, nm in MAGBINS) +
                  f"; rolled/ETS {rolled_frac_ets:.3f}; quiet FA {quiet_fa:.2f}; extras {len(X)} ({frac_extra_ets:.2f} on ETS days; amp {amp_x:.3g} vs matched {amp_m:.3g})", flush=True)
            for k, v in bars.items():
                print(f"   {k:20s} {v}")
            if a.write_catalog and ((not a.e6 and a.tag == "e5r_snr" and T == 8.0) or (a.e6 and T == a.write_T)):
                cat = pd.DataFrame(dict(
                    time=pd.to_datetime(C.t.to_numpy(), unit="s").astype("datetime64[ns]"), family_id=C.fam.astype(int).to_numpy(),
                    sb_family_id=C.fam.astype(int).map(F.sb_fam_id).to_numpy(), lat=C.fam.astype(int).map(F.lat).to_numpy(),
                    lon=C.fam.astype(int).map(F.lon).to_numpy(), depth=C.fam.astype(int).map(F.dep).to_numpy(),
                    rel_amp=C.amp.to_numpy(), magnitude=np.log10(np.clip(C.amp.to_numpy(), 1e-6, None)), n_stations=C.n_sta.to_numpy(),
                    network_mad=C.net_mad.to_numpy(), cc_best=C.net_cc.to_numpy(), template_len_s=L,
                    matched_bostock_idx=np.where(C.matched_bostock_idx >= 0, C.matched_bostock_idx, np.nan),
                    matched_bostock_mag=C.matched_bostock_mag.to_numpy(), matched_same_family=C.matched_same_family.to_numpy(), doy=C.doy.to_numpy()))
                if a.e6:
                    cat["matched_bostock_id"] = np.where(C.matched_bostock_id >= 0, C.matched_bostock_id, np.nan)     # file row of e5r_gt_bostock2012.parquet
                    cat = cat.drop(columns=["matched_bostock_idx"])
                fn = f"{RES}/e5_catalog_2012_sVI{'' if L == 8.0 else '_L12'}.parquet" if not a.e6 else \
                    f"{RES}/e6_catalog_2012_sVI_{a.tag}_T{T:g}{'' if L == 8.0 else '_L12'}.parquet"
                cat.to_parquet(fn, index=False); print(f"[write] {fn} ({len(cat)} events)")
    pd.DataFrame(fa_rows).to_csv(f"{RES}/{PFX}_{a.tag}_fa.csv", index=False)
    pd.DataFrame(mag_rows).to_csv(f"{RES}/{PFX}_{a.tag}_recall_mag.csv", index=False)
    pd.DataFrame(nsta_rows).to_csv(f"{RES}/{PFX}_{a.tag}_recall_nsta.csv", index=False)
    pd.DataFrame(fam_rows).to_csv(f"{RES}/{PFX}_{a.tag}_recall_family.csv", index=False)
    pd.DataFrame(daily_rows).to_csv(f"{RES}/{PFX}_{a.tag}_daily.csv", index=False)
    pd.DataFrame(ex_rows).to_csv(f"{RES}/{PFX}_{a.tag}_extras.csv", index=False)
    # per-station registration at test events, M1.8-2.1 (failure-signature diagnostic)
    reg = []
    e18 = ev[(ev.L == 8.0) & (ev.half == "test") & (~ev.guard) & (ev.mag >= 1.8) & (ev.mag < 2.1)]
    s8 = st[st.L == 8.0]
    for c in [c for c in ev.columns if c.startswith("sta_")]:
        s_ = c[4:]; v = e18[c].dropna()
        rb = s8[s8.sta == s_].rand3_median.median()
        reg.append(dict(sta=s_, n_ev=int(len(v)), mad_median=float(v.median()) if len(v) else np.nan, rand_median=float(rb),
                        excess=float(v.median() - rb) if len(v) else np.nan, frac_ge4=float((v >= 4).mean()) if len(v) else np.nan))
    pd.DataFrame(reg).to_csv(f"{RES}/{PFX}_{a.tag}_station_registration.csv", index=False)
    json.dump(out, open(f"{RES}/{PFX}_{a.tag}.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print("\n[fa]"); print(pd.DataFrame(fa_rows).round(3).to_string(index=False))
    print("\n[station registration at test events, M1.8-2.1, 8 s]"); print(pd.DataFrame(reg).round(2).to_string(index=False))
    print(f"[write] results/e5r_score_{a.tag}.json + csv tables")


if __name__ == "__main__":
    main()
