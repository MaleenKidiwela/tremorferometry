#!/usr/bin/env python
"""E5 step 1 analysis -- the DETECTION-BUDGET table (Merlin E4-audit F1, prescription 1).

Inputs  results/e5_budget_<tag>.csv   per (day, L, family, threshold): declustered detections, rolled-control
                                      detections, matches (+-2 s) to the same / any Bostock family, test-block split,
                                      8-MAD extras within 60 s of a NEIGHBOURING family's event vs random times
        results/e4_tt_events_<tag>.parquet   per event: net_mad_2s (max over +-2 s), half, mag, n_ch
QUIET day = <= 50 Bostock events that UTC day (all families, full catalogue); everything else = ETS day.
Outputs results/e5_budget_table_<tag>.csv, e5_budget_recall_<tag>.csv, e5_budget_bynsta_<tag>.csv
"""
from __future__ import annotations
import sys, numpy as np, pandas as pd

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
MB = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1")]


def quiet_days(year=2010):
    c = pd.read_csv(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
    c = c[c.year == year]
    t = pd.to_datetime(dict(year=c.year, month=c.month, day=c.day)).astype("datetime64[ns]")
    n = t.dt.dayofyear.value_counts()
    return n


def main(tag):
    B = pd.read_csv(f"{RES}/e5_budget_{tag}.csv")
    E = pd.read_parquet(f"{RES}/e4_tt_events_{tag}.parquet")
    nday = quiet_days()
    B["n_bostock_day"] = B.doy.map(nday).fillna(0).astype(int)
    B["quiet"] = B.n_bostock_day <= 50
    B["rate"] = B.n_det * 86400.0 / B.n_valid_s
    B["rate_unmatched_any"] = (B.n_det - B.n_match_any) * 86400.0 / B.n_valid_s
    B["rate_roll"] = B.n_det_roll * 86400.0 / B.n_valid_s
    rows = []
    for (L, T), g in B.groupby(["L", "T"]):
        q, e = g[g.quiet], g[~g.quiet]
        r = dict(L=L, T=T, n_quiet_famdays=len(q), n_ets_famdays=len(e), quiet_days=int(q.doy.nunique()), ets_days=int(e.doy.nunique()),
                 quiet_fa_mean=q.rate.mean(), quiet_fa_p50=q.rate.median(), quiet_fa_p90=q.rate.quantile(.9),
                 quiet_fa_unmatched_mean=q.rate_unmatched_any.mean(), quiet_roll_mean=q.rate_roll.mean(),
                 ets_det_mean=e.rate.mean(), ets_roll_mean=e.rate_roll.mean(),
                 ets_precision_own=e.n_match_own.sum() / max(e.n_det.sum(), 1),
                 ets_precision_any=e.n_match_any.sum() / max(e.n_det.sum(), 1),
                 ets_precision_own_testblocks=e.n_match_own_test.sum() / max(e.n_det_test.sum(), 1),
                 ets_precision_any_testblocks=e.n_match_any_test.sum() / max(e.n_det_test.sum(), 1),
                 ets_roll_matchfrac_own=e.n_roll_match_own.sum() / max(e.n_det_roll.sum(), 1),
                 ets_extras_per_tpl_day=((e.n_det - e.n_match_own) * 86400.0 / e.n_valid_s).mean(),
                 ets_total_det=int(e.n_det.sum()), ets_total_matched_own=int(e.n_match_own.sum()),
                 ets_total_roll=int(e.n_det_roll.sum()),
                 ev_recall_matched_own=e.n_match_own.sum() / max(e.n_ev_own.sum(), 1))
        # extras near a NEIGHBOURING family's event (60 s) vs random times, ETS days
        r["extras_near60_other_frac"] = e.n_extra_near60_other.sum() / max(e.n_extra.sum(), 1)
        r["random_near60_other_frac"] = e.n_rand_near60_other.sum() / max(e.n_rand.sum(), 1)
        # recall from the event table (+-2 s common lag), all events and TEST half
        h = E[E.L == L]
        r["recall_all_events"] = float((h.net_mad_2s >= T).mean())
        r["recall_test_half"] = float((h[h.half == "test"].net_mad_2s >= T).mean())
        r["recall_shift_ctrl"] = float((h.net_mad_shift >= T).mean())
        for lo, hi, lab in MB:
            k = (h.mag > lo) & (h.mag <= hi)
            r[f"recall_{lab}"] = float((h.net_mad_2s[k] >= T).mean())
            r[f"recall_test_{lab}"] = float((h.net_mad_2s[k & (h.half == 'test')] >= T).mean())
        rows.append(r)
    TB = pd.DataFrame(rows); TB.to_csv(f"{RES}/e5_budget_table_{tag}.csv", index=False)
    # recall by number of stations carrying a template (n_ch / 3), L = 8 and 12, thresholds
    bn = []
    for L in sorted(E.L.unique()):
        h = E[E.L == L].copy(); h["n_sta"] = (h.n_ch // 3).astype(int)
        for ns, g in h.groupby("n_sta"):
            if ns == 0 or len(g) < 30:
                continue
            bn.append(dict(L=L, n_sta=ns, n_ev=len(g), mag_med=g.mag.median(), raw_med=g.net_mad.median(),
                           **{f"recall{int(T)}": float((g.net_mad_2s >= T).mean()) for T in (6, 7, 8, 9, 10)}))
    BN = pd.DataFrame(bn); BN.to_csv(f"{RES}/e5_budget_bynsta_{tag}.csv", index=False)
    # precision / FA by the number of stations carrying a template (family-day level), L = 8, T = 8
    pn = []
    for (T_, ns), g in B[(B.L == 8.0) & (B["T"].isin([8.0, 10.0]))].groupby(["T", "n_sta"]):
        q, e = g[g.quiet], g[~g.quiet]
        if len(e) < 20:
            continue
        pn.append(dict(T=T_, n_sta=ns, n_famdays=len(e), n_fam=int(e.fam.nunique()),
                       quiet_fa=(q.n_det * 86400 / q.n_valid_s).mean() if len(q) else np.nan,
                       ets_det=(e.n_det * 86400 / e.n_valid_s).mean(), ets_roll=(e.n_det_roll * 86400 / e.n_valid_s).mean(),
                       precision_own=e.n_match_own.sum() / max(e.n_det.sum(), 1), precision_any=e.n_match_any.sum() / max(e.n_det.sum(), 1),
                       recall_matched_own=e.n_match_own.sum() / max(e.n_ev_own.sum(), 1)))
    PN = pd.DataFrame(pn); PN.to_csv(f"{RES}/e5_budget_bynsta_precision_{tag}.csv", index=False)
    print("precision by stations carrying a template (L=8):"); print(PN.round(3).to_string(index=False))
    # operating threshold: smallest T with quiet-day FA mean <= 1/template/day at L = 8 (and precision >= 0.85)
    t8 = TB[TB.L == 8.0].sort_values("T")
    op = t8[(t8.quiet_fa_mean <= 1.0) & (t8.ets_precision_own >= 0.85)]
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 60)
    print(f"=== {tag}: quiet days {sorted(B[B.quiet].doy.unique())} ===")
    cols = ["L", "T", "quiet_fa_mean", "quiet_fa_p90", "quiet_fa_unmatched_mean", "quiet_roll_mean", "ets_det_mean", "ets_roll_mean",
            "ets_precision_own", "ets_precision_any", "ets_precision_own_testblocks", "ets_roll_matchfrac_own",
            "recall_all_events", "recall_test_half", "recall_M<1.5", "recall_M1.5-1.8", "recall_M1.8-2.1", "recall_M>=2.1",
            "extras_near60_other_frac", "random_near60_other_frac"]
    print(TB[cols].round(3).to_string(index=False))
    print(f"operating threshold at L=8 (quiet FA <= 1/tpl/day AND precision >= 0.85): "
          f"{op['T'].min() if len(op) else 'NONE in 5-10'}")
    print(BN.round(3).to_string(index=False))


if __name__ == "__main__":
    for t in sys.argv[1:] or ["e5_budget_op"]:
        main(t)
