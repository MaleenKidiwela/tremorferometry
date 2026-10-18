#!/usr/bin/env python
"""E5 step 3 analysis -- SNR-weighted vs equal-weight network, HELD OUT (test half of alternate 2-h blocks).
Tags e5_w_{op,ring}_{eq,snr}.  Families 1 and 15 (Merlin F10: ~0 excess everywhere) reported separately.
Outputs results/e5_weighting_table.csv (per tag, L, magbin, family set), e5_weighting_corr.csv (do build-half
stack S-SNRs predict test-half single-station registration?), e5_weighting_fa.csv (quiet-day FA under weights)."""
from __future__ import annotations
import numpy as np, pandas as pd, os
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
MB = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1"), (-9, 9, "all")]
TAGS = ["e5_w_op_eq", "e5_w_op_snr", "e5_w_ring_eq", "e5_w_ring_snr"]
BLOCK = [1.0, 15.0]


def main():
    rows, fa = [], []
    for tag in TAGS:
        if not os.path.exists(f"{RES}/e4_tt_events_{tag}.parquet"):
            continue
        E = pd.read_parquet(f"{RES}/e4_tt_events_{tag}.parquet")
        NF = pd.read_csv(f"{RES}/e4_tt_netfam_{tag}.csv")
        for L in sorted(E.L.unique()):
            nb = NF[NF.L == L]
            for fams, lab in ((None, "all_fams"), (BLOCK, "no_1_15"), ("only", "only_1_15")):
                g = E[E.L == L]; nbL = nb
                if fams is None:
                    pass
                elif fams == "only":
                    g = g[g.fam.isin(BLOCK)]; nbL = nb[nb.fam.isin(BLOCK)]
                else:
                    g = g[~g.fam.isin(fams)]; nbL = nb[~nb.fam.isin(fams)]
                bgv = float(nbL.net_rand_med.median()) if len(nbL) else np.nan
                for lo, hi, mlab in MB:
                    h = g[(g.mag > lo) & (g.mag <= hi)]
                    if len(h) < 30:
                        continue
                    hv = h[h.net_mad.notna()]
                    rows.append(dict(tag=tag, L=L, fams=lab, magbin=mlab, n=len(h), n_scored=len(hv),
                                     raw_med=float(hv.net_mad.median()), rand_bg=bgv, excess=float(hv.net_mad.median() - bgv),
                                     frac_ge8=float((h.net_mad >= 8).mean()), frac_ge8_2s=float((h.net_mad_2s >= 8).mean()),
                                     frac_ge10_2s=float((h.net_mad_2s >= 10).mean()),
                                     frac_ge8_scored=float((hv.net_mad >= 8).mean()),
                                     shift_frac_ge8=float((h.net_mad_shift >= 8).mean()),
                                     net_fa8_mean=float(nbL.net_fa8.mean()), shift_fa8_mean=float(nbL.shift_fa8.mean())))
        if os.path.exists(f"{RES}/e5_budget_{tag}.csv"):
            B = pd.read_csv(f"{RES}/e5_budget_{tag}.csv")
            nd = pd.read_csv(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"); nd = nd[nd.year == 2010]
            cnt = pd.to_datetime(dict(year=nd.year, month=nd.month, day=nd.day)).dt.dayofyear.value_counts()
            B["quiet"] = B.doy.map(cnt).fillna(0) <= 50
            for (L, T), gg in B.groupby(["L", "T"]):
                q = gg[gg.quiet]; e = gg[~gg.quiet]
                fa.append(dict(tag=tag, L=L, T=T, quiet_fa_mean=(q.n_det * 86400 / q.n_valid_s).mean(),
                               quiet_fa_p90=(q.n_det * 86400 / q.n_valid_s).quantile(.9),
                               quiet_roll_mean=(q.n_det_roll * 86400 / q.n_valid_s).mean(),
                               ets_det_mean=(e.n_det * 86400 / e.n_valid_s).mean(), ets_roll_mean=(e.n_det_roll * 86400 / e.n_valid_s).mean(),
                               ets_precision_own_test=e.n_match_own_test.sum() / max(e.n_det_test.sum(), 1),
                               ets_precision_any_test=e.n_match_any_test.sum() / max(e.n_det_test.sum(), 1)))
    W = pd.DataFrame(rows); W.to_csv(f"{RES}/e5_weighting_table.csv", index=False)
    F = pd.DataFrame(fa); F.to_csv(f"{RES}/e5_weighting_fa.csv", index=False)
    # FA-MATCHED comparison: the threshold at which the weighted network has the equal-weight network's quiet-day
    # FA at 8 MAD (log-linear interpolation of FA in T between the budget thresholds), and the held-out recall there
    # (net_mad_2s >= T, M1.8-2.1 and all, families 1/15 excluded) -- the honest gain of weighting.
    fm = []
    for net in ("op", "ring"):
        te, tw = f"e5_w_{net}_eq", f"e5_w_{net}_snr"
        if not (os.path.exists(f"{RES}/e5_budget_{te}.csv") and os.path.exists(f"{RES}/e5_budget_{tw}.csv")):
            continue
        for L in (8.0, 12.0):
            fe = F[(F.tag == te) & (F.L == L)].set_index("T").quiet_fa_mean; fw = F[(F.tag == tw) & (F.L == L)].set_index("T").quiet_fa_mean
            target = float(fe.get(8.0, np.nan))
            Ts = np.array(sorted(fw.index)); vals = np.log(np.maximum(fw.loc[Ts].to_numpy(), 1e-4))
            Tm = float(np.interp(np.log(max(target, 1e-4)), vals[::-1], Ts[::-1]))     # FA decreases with T
            Ee = pd.read_parquet(f"{RES}/e4_tt_events_{te}.parquet"); Ew = pd.read_parquet(f"{RES}/e4_tt_events_{tw}.parquet")
            Ee = Ee[(Ee.L == L) & ~Ee.fam.isin(BLOCK)]; Ew = Ew[(Ew.L == L) & ~Ew.fam.isin(BLOCK)]
            k = lambda d: (d.mag > 1.8) & (d.mag <= 2.1)
            fm.append(dict(net=net, L=L, quiet_fa_equal_T8=target, T_weighted_matched=Tm,
                           quiet_fa_weighted_T8=float(fw.get(8.0, np.nan)),
                           recall_equal_T8_M18_21=float((Ee.net_mad_2s[k(Ee)] >= 8).mean()),
                           recall_weighted_T8_M18_21=float((Ew.net_mad_2s[k(Ew)] >= 8).mean()),
                           recall_weighted_Tmatched_M18_21=float((Ew.net_mad_2s[k(Ew)] >= Tm).mean()),
                           recall_equal_T8_all=float((Ee.net_mad_2s >= 8).mean()),
                           recall_weighted_T8_all=float((Ew.net_mad_2s >= 8).mean()),
                           recall_weighted_Tmatched_all=float((Ew.net_mad_2s >= Tm).mean())))
    FM = pd.DataFrame(fm); FM.to_csv(f"{RES}/e5_weighting_famatched.csv", index=False)
    print("FA-MATCHED gain of weighting:"); print(FM.round(3).to_string(index=False))
    # do build-half weights predict test-half registration?  per (fam, sta): build S-SNR vs test single-station excess
    Wt = pd.read_csv(f"{RES}/e5_weights_snr.csv")
    src = "e5_w_ring_eq" if os.path.exists(f"{RES}/e4_tt_events_e5_w_ring_eq.parquet") else TAGS[0]
    E = pd.read_parquet(f"{RES}/e4_tt_events_{src}.parquet"); ST = pd.read_csv(f"{RES}/e4_tt_station_{src}.csv")
    bg = ST[ST.L == 8.0].groupby(["fam", "sta"]).rand3_median.median()
    cc = []
    for (fam, sta), r in bg.items():
        c = f"sta_{sta}"
        if c not in E.columns:
            continue
        h = E[(E.L == 8.0) & (E.fam == fam) & E[c].notna()]
        if len(h) < 20:
            continue
        w = Wt[(Wt.fam == fam) & (Wt.sta == sta)]
        if not len(w):
            continue
        cc.append(dict(fam=fam, sta=sta, n=len(h), excess_test=float(h[c].median() - r), snr_db=float(w.snr_db.iloc[0]), w=float(w.w.iloc[0])))
    C = pd.DataFrame(cc)
    if len(C):
        C["rho"] = C.snr_db.corr(C.excess_test, method="spearman")
        C["rho_within_fam"] = np.nanmedian([g.snr_db.corr(g.excess_test, method="spearman") for _, g in C.groupby("fam") if len(g) >= 4])
    C.to_csv(f"{RES}/e5_weighting_corr.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(W[(W.L == 8.0) & W.magbin.isin(["M1.8-2.1", "all"]) & (W.fams != "only_1_15")].round(3).to_string(index=False))
    print(W[(W.L == 8.0) & (W.fams == "only_1_15") & (W.magbin == "all")].round(3).to_string(index=False))
    print(F[(F.L == 8.0) & (F["T"].isin([8.0, 10.0]))].round(3).to_string(index=False))
    if len(C):
        print(f"build S-SNR vs test excess: Spearman rho = {C.rho.iloc[0]:.3f} over {len(C)} (fam, sta) pairs; within-family median rho {C.rho_within_fam.iloc[0]:.3f}")
        print(C.groupby("sta").agg(n=("fam", "size"), snr=("snr_db", "median"), exc=("excess_test", "median"), w=("w", "median")).round(2).to_string())


if __name__ == "__main__":
    main()
