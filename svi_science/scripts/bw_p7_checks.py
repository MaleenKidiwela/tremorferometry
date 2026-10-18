#!/usr/bin/env python
"""bw_p7_checks -- P7 step checks (a)-(e), control fairness, scrambled-control ΣP, withheld S-energy, station-set invariance and
the G6 verdict per era (rules: results/backward/p7/P7.md, written before computing).
Needs: the catalogue (bw_p7_catalogue.py), p7/tables/k_chain.csv + k_links.csv (bw_p7_kchain.py), p6/tables/onelook.csv.
Out: p7/tables/step_checks.csv, step_a.json, fairness_slope.csv, step_d.csv, step_e.csv, control_sumP.csv, withheld_energy.csv,
     withheld_energy_bins.csv, invariance.csv, invariance_monthly.csv, g6.json, daily.parquet
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_p7_common as Q7
warnings.filterwarnings("ignore")

P7 = Q7.P7; T7 = f"{P7}/tables"
NB = 2000
EDGES = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0001]
BOUNDARIES = [  # (date, label, series, kind)
    ("2005-07-01", "set CN4 -> B1", "frozen", "set"), ("2005-11-02", "set B1 -> B2", "frozen", "set"),
    ("2007-11-13", "set B2 -> B3", "frozen", "set"), ("2009-08-27", "set B3 -> B4 (S15)", "frozen", "set"),
    ("2009-11-20", "LRIV first data = the real B3 -> B4 change (added after the ladder showed LRIV absent before; deviation)", "frozen", "set"),
    ("2006-07-01", "detector era -> frozen B2", "primary", "detector"), ("2014-01-01", "backward B4 -> forward production", "junction", "junction"),
    ("2006-04-20", "SNB CMG40T -> CMG3ESP/GD2", "frozen", "instrument"), ("2007-01-27", "NLLB CMG3T month", "frozen", "instrument"),
    ("2007-02-28", "VGZ CMG3ESP/GD1 -> CMG3T", "frozen", "instrument"), ("2008-12-11", "SNB GD2 -> GD1", "frozen", "instrument"),
    ("2009-07-31", "VGZ BH -> HH", "frozen", "instrument"), ("2011-03-17", "SNB BH -> HH", "frozen", "instrument"),
    ("2003-07-31", "KELB gain 6.17e7 -> 3.07e8 (+ POLARIS arrivals 07-25..08-12)", "era", "era-station"),
    ("2004-06-04", "KELB -> KLNB rename", "era", "era-station"), ("2004-08-05", "KHVB ends", "era", "era-station"),
    ("2004-09-30", "SHVB ends", "era", "era-station"), ("2005-01-11", "SHDB ends", "era", "era-station"),
    ("2005-05-10", "TWGB ends", "era", "era-station"), ("2006-05-13", "MGCB / TWKB end", "era", "era-station"),
    ("2006-06-16", "SSIB / TSJB end", "era", "era-station")]


def ratio_boot(a_num, a_den, b_num, b_den, seed=0):
    """ln((Σa_num/Σa_den) / (Σb_num/Σb_den)) with a 2000-draw day bootstrap within each side."""
    ra = a_num.sum() / max(a_den.sum(), 1e-9); rb = b_num.sum() / max(b_den.sum(), 1e-9)
    if not (ra > 0 and rb > 0) or len(a_num) < 5 or len(b_num) < 5:
        return dict(ratio=ra / rb if rb != 0 else np.nan, z=np.nan, sigma_ln=np.nan, before=rb, after=ra, n_before=len(b_num), n_after=len(a_num))
    rng = np.random.default_rng(seed); L = []
    for _ in range(NB):
        i = rng.integers(0, len(a_num), len(a_num)); j = rng.integers(0, len(b_num), len(b_num))
        x = a_num[i].sum() / max(a_den[i].sum(), 1e-9); y = b_num[j].sum() / max(b_den[j].sum(), 1e-9)
        if x > 0 and y > 0: L.append(np.log(x / y))
    s = np.std(L, ddof=1) if len(L) > 50 else np.nan
    return dict(ratio=ra / rb, z=np.log(ra / rb) / s if s > 0 else np.nan, sigma_ln=s, before=rb, after=ra, n_before=len(b_num),
                n_after=len(a_num), frac_boot_pos=len(L) / NB)


def main():
    os.makedirs(T7, exist_ok=True)
    A = pd.read_parquet(Q7.CAT, columns=["run", "date", "series", "fam", "P", "net_set", "detector", "multi_site_z2", "burst_day",
                                         "conf_ze", "quiet", "quiet_src", "primary"])
    FD = Q7.load_fd(Q7.PROD)
    D = Q7.daily(A, FD); DF = Q7.day_flags()
    D = D.join(DF, on="date")
    KC = Q7.k_chain("all"); KCm = Q7.k_chain("ms")
    Kmap = {s: KC.loc[s, "K"] for s in KC.index}
    D["K"] = D.net_set.map(Kmap); D["K_ms"] = D.net_set.map({s: KCm.loc[s, "K"] for s in KCm.index})
    D["exc"] = D.spd - D.spc; D["exc_kc"] = D.exc * D.K; D["exc_ms"] = D.spd_ms - D.spc_ms
    D.to_parquet(f"{T7}/daily.parquet", index=False)
    res = {}
    # ---------------- (a)
    a = {}
    for f, k in (("p4b_e0_decision.json", "step1_E0"), ("p4b_era_decision.json", "step2_era")):
        p = f"{Q7.BW}/p4b/tables/{f}"
        a[k] = json.load(open(p)) if os.path.exists(p) else None
    json.dump(a, open(f"{T7}/step_a.json", "w"), indent=1, default=float)
    # ---------------- (b) step checks
    fw = pd.read_parquet("/home/jovyan/tremorferometry/lfe_extension/results/refine/rates/famday.parquet", columns=["date", "excess"])
    fwd = fw.groupby("date").excess.sum()
    rows = []
    for d0, lab, ser, kind in BOUNDARIES:
        t0 = pd.Timestamp(d0)
        bef = [str(x.date()) for x in pd.date_range(t0 - pd.Timedelta(days=60), t0 - pd.Timedelta(days=1))]
        aft = [str(x.date()) for x in pd.date_range(t0, t0 + pd.Timedelta(days=59))]
        if ser == "junction":
            B = D[(D.detector == "frozen") & D.date.isin(bef)]
            import bw_p4_common as P4
            aft_between = dict(zip(aft, ~P4.bostock_near(aft, pad=3)))          # 2014: no env flag; between from the window dates
            for cl in ("between", "all"):
                b = B if cl == "all" else B[B[cl]]
                ad = [d for d in aft if d in fwd.index and (cl == "all" or aft_between[d])]
                an = fwd.reindex(ad).to_numpy(float)
                r = ratio_boot(an, np.ones(len(an)), b.exc.to_numpy(), np.ones(len(b)))
                rows.append(dict(date=d0, change=lab, kind=kind, series="network excess per day (backward B4 vs forward prod)", corrected="n/a (same model)",
                                 day_class=cl, **r))
            continue
        if ser == "frozen":
            S = D[D.detector == "frozen"]
        elif ser == "era":
            S = D[D.detector == "era"]
        else:
            S = D[D.primary]
        for cl in ("between", "env_quiet"):
            Sx = S[S[cl]]
            b = Sx[Sx.date.isin(bef)]; f = Sx[Sx.date.isin(aft)]
            for corr, col in (("k-corrected", "exc_kc"), ("uncorrected", "exc")):
                r = ratio_boot(f[col].to_numpy(), f.m_data.to_numpy(float), b[col].to_numpy(), b.m_data.to_numpy(float))
                rows.append(dict(date=d0, change=lab, kind=kind, series=ser, corrected=corr, day_class=cl, **r))
    SC = pd.DataFrame(rows)
    SC["flagged"] = SC.z.abs() > 2
    SC["near_ets"] = [bool(DF.ets.reindex([str(x.date()) for x in pd.date_range(pd.Timestamp(d) - pd.Timedelta(days=60),
                                                                               pd.Timestamp(d) + pd.Timedelta(days=59))]).fillna(False).any()) for d in SC.date]
    SC.to_csv(f"{T7}/step_checks.csv", index=False)
    # diagnostic (added after the first run, not a pass rule): the same statistic at 'placebo' boundaries = the 1st of every month
    # with no listed change within +-60 d -> how often |z| > 2 happens with no network change at all (activity alone)
    bset = [pd.Timestamp(b[0]) for b in BOUNDARIES]
    pl = []
    for ser in ("frozen", "era"):
        S = D[D.detector == ser]
        lo_, hi_ = (pd.Timestamp("2003-03-01"), pd.Timestamp("2013-11-01")) if ser == "frozen" else (pd.Timestamp("2003-04-01"), pd.Timestamp("2006-05-01"))
        for t0 in pd.date_range(lo_, hi_, freq="MS"):
            if any(abs((t0 - b).days) <= 60 for b in bset): continue
            bef = [str(x.date()) for x in pd.date_range(t0 - pd.Timedelta(days=60), t0 - pd.Timedelta(days=1))]
            aft = [str(x.date()) for x in pd.date_range(t0, t0 + pd.Timedelta(days=59))]
            for cl in ("between", "env_quiet"):
                Sx = S[S[cl]]; b = Sx[Sx.date.isin(bef)]; f = Sx[Sx.date.isin(aft)]
                r = ratio_boot(f.exc_kc.to_numpy(), f.m_data.to_numpy(float), b.exc_kc.to_numpy(), b.m_data.to_numpy(float))
                pl.append(dict(series=ser, date=str(t0.date()), day_class=cl, **r))
    PL = pd.DataFrame(pl); PL["flagged"] = PL.z.abs() > 2; PL.to_csv(f"{T7}/step_placebo.csv", index=False)
    res["b_placebo"] = PL.groupby(["series", "day_class"]).agg(n=("z", lambda z: int(z.notna().sum())), frac_flagged=("flagged", "mean"),
                                                                 median_abs_z=("z", lambda z: float(z.abs().median()))).reset_index().to_dict("records")
    key = SC[SC.kind.isin(["set", "detector", "junction"]) & (SC.day_class == "between") & SC.corrected.isin(["k-corrected", "n/a (same model)"])]
    res["b"] = dict(pass_=bool(not key.flagged.any() and key.z.notna().all()), flagged=key[key.flagged | key.z.isna()][["date", "change", "ratio", "z"]].to_dict("records"))
    # ---------------- (c) control fairness slope (env-quiet days)
    fr = []
    for (det, ns), g in D[D.env_quiet & (D.m_data > 0)].groupby(["detector", "net_set"]):
        for rw, (yd, xc) in (("all", (g.exc, g.spc)), ("ms", (g.exc_ms, g.spc_ms))):
            y = (yd / g.m_data).to_numpy(); x = (xc / g.m_data).to_numpy()
            if len(y) < 10 or np.std(x) == 0:
                continue
            sl = np.polyfit(x, y, 1)[0]; rng = np.random.default_rng(3); bs = []
            for _ in range(NB):
                i = rng.integers(0, len(y), len(y))
                if np.std(x[i]) > 0: bs.append(np.polyfit(x[i], y[i], 1)[0])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            fr.append(dict(detector=det, net_set=ns, rows=rw, n_days=len(y), slope=sl, lo95=lo, hi95=hi, pass_ci_contains_0=bool(lo <= 0 <= hi),
                           small_abs_le_025=bool(abs(sl) <= 0.25), mean_exc_per_fd=float(y.mean()), mean_ctl_per_fd=float(x.mean())))
    FR = pd.DataFrame(fr); FR.to_csv(f"{T7}/fairness_slope.csv", index=False)
    res["c"] = FR[FR.rows == "all"][["detector", "net_set", "slope", "lo95", "hi95", "pass_ci_contains_0"]].to_dict("records")
    # ---------------- (d) recall vs ladder k
    OL = pd.read_csv(f"{Q7.P6}/tables/onelook.csv"); KL = pd.read_csv(f"{T7}/k_links.csv")
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension/scripts")
    import bw_gate_common as G
    R, C = G.patch_forward()
    CC = pd.read_parquet(f"{G.REFINE}/conf/candidates_conf.parquet", columns=["series", "sset_run", "row_role", "year", "date", "t", "fam", "P_same", "m_same"])
    CC = CC[(CC.row_role == "fit") & (CC.sset_run == "S15") & (CC.series == "data")]
    DT = pd.read_csv(f"{G.REFINE}/day_table.csv", dtype={"date": str})
    Bf = []
    for y in (2011, 2012):
        days = sorted(DT[(DT.kind == f"truth{y}") & (DT.role == "fit")].date)
        b = R.bostock([y]); b = b[b.date.isin(days)]
        if y == 2012: b = b[R.block_role(b.t_origin.to_numpy()) == "fit"]
        Bf.append(b)
    Bf = pd.concat(Bf); s = CC[CC.P_same >= 0.75]
    rec_ref = float(R.one2one_recall(s.t.to_numpy(), s.fam.to_numpy(), Bf.t_origin.to_numpy(), Bf.fam.to_numpy()).mean())
    dd = []
    kl = KL[KL.rows == "all"].set_index("pair")
    for y in sorted(OL.year.unique()):
        o = OL[OL.year == y].set_index("detector")
        if "era" in o.index and "frozen" in o.index:
            fs = o.loc["frozen", "net_set"]; pair = f"era/{fs}"
            k_ets = kl.loc[pair, "k_ets"] if pair in kl.index else np.nan
            rr = o.loc["era", "recall_P075"] / o.loc["frozen", "recall_P075"] if o.loc["frozen", "recall_P075"] > 0 else np.nan
            dd.append(dict(test="d1 era vs frozen, same window", year=y, already_seen=bool(o.loc["era", "already_seen"]), pair=pair,
                           recall_ratio=rr, k=k_ets, k_kind="overlap k, ETS days", ratio_of_ratios=rr / k_ets if k_ets > 0 else np.nan,
                           pass_x15=bool(np.isfinite(rr) and k_ets > 0 and 1 / 1.5 <= rr / k_ets <= 1.5)))
        for det, r in o.iterrows():
            ns = r.net_set; Kv = KC.loc[ns, "K"] if ns in KC.index else np.nan
            rr = r.recall_P075 / rec_ref
            dd.append(dict(test="d2 chain vs S15 reference recall", year=y, already_seen=bool(r.already_seen), pair=f"{ns}/S15", recall_ratio=rr,
                           k=1 / Kv if Kv > 0 else np.nan, k_kind="1 / K_set (ladder chain, env-quiet)", recall_ref_S15_2011_12=rec_ref,
                           ratio_of_ratios=rr * Kv if Kv > 0 else np.nan, pass_x15=bool(Kv > 0 and 1 / 1.5 <= rr * Kv <= 1.5)))
    DD = pd.DataFrame(dd); DD.to_csv(f"{T7}/step_d.csv", index=False)
    res["d"] = DD[~DD.already_seen][["test", "year", "pair", "recall_ratio", "k", "ratio_of_ratios", "pass_x15"]].to_dict("records")
    # ---------------- (e) 2010-13 between rates, env-quiet vs PNSN-quiet days (B4)
    B4 = D[(D.detector == "frozen") & (D.date >= "2010-01-10") & (D.date <= "2013-12-31")]
    ee = []
    for y, g in list(B4.groupby(B4.date.str[:4])) + [("2010-13", B4)]:
        e1 = g[g.env_quiet]; e2 = g[g.pnsn_quiet_b]
        r1 = e1.exc.sum() / max(e1.m_data.sum(), 1); r2 = e2.exc.sum() / max(e2.m_data.sum(), 1)
        ee.append(dict(year=y, n_env_quiet=len(e1), n_pnsn_quiet=len(e2), n_both=int((g.env_quiet & g.pnsn_quiet_b).sum()),
                       exc_per_fd_env_quiet=r1, exc_per_fd_pnsn_quiet=r2, ratio_env_over_pnsn=r1 / r2 if r2 > 0 else np.nan))
    EE = pd.DataFrame(ee); EE["agree_08_12"] = EE.ratio_env_over_pnsn.between(0.8, 1.2); EE.to_csv(f"{T7}/step_e.csv", index=False)
    res["e"] = EE.to_dict("records")
    # ---------------- control ΣP
    cs = []
    for (det, ns), g in D.groupby(["detector", "net_set"]):
        for y, h in list(g.groupby(g.date.str[:4])) + [("all", g)]:
            for cl, hh in (("all_days", h), ("between", h[h.between])):
                cs.append(dict(detector=det, net_set=ns, year=y, days=cl, n_days=len(hh), sumP_data=hh.spd.sum(), sumP_ctl=hh.spc.sum(),
                               ctl_over_data=hh.spc.sum() / max(hh.spd.sum(), 1e-9)))
    CS = pd.DataFrame(cs)
    from scipy.stats import spearmanr
    spr = {}
    for (det, ns), g in D.groupby(["detector", "net_set"]):
        mo = g.groupby(g.date.str[:7]).agg(spc=("spc", "sum"), n=("date", "size"))
        spr[f"{det}:{ns}"] = float(spearmanr(np.arange(len(mo)), mo.spc / mo.n).correlation) if len(mo) >= 6 else np.nan
    CS["spearman_monthly_ctl_vs_time"] = [spr.get(f"{a}:{b}") for a, b in zip(CS.detector, CS.net_set)]
    CS.to_csv(f"{T7}/control_sumP.csv", index=False)
    # ---------------- withheld S-energy (quiet days of the catalogue)
    W = A[A.quiet & A.conf_ze.notna()].copy(); W["bin"] = pd.cut(W.P, EDGES, right=False)
    we, wb = [], []
    for (det, ns), g in W.groupby(["detector", "net_set"]):
        Eb = []
        for b, h in g.groupby("bin", observed=False):
            d_ = h[h.series == "data"]; c_ = h[h.series == "c1"]
            Eb.append(dict(detector=det, net_set=ns, P_bin=str(b), n_data=len(d_), n_ctl=len(c_), conf_ze_data=float(d_.conf_ze.mean()) if len(d_) else np.nan,
                           conf_ze_ctl=float(c_.conf_ze.mean()) if len(c_) else np.nan))
        Eb = pd.DataFrame(Eb); wb.append(Eb)
        hi = g[(g.series == "data") & (g.P >= 0.5)]; ct = g[g.series == "c1"]
        import svi_ref_val as V
        d = V.boot_mean_diff(hi.conf_ze.values, hi.date.values, ct.conf_ze.values, ct.date.values) if len(hi) >= 30 else (np.nan,) * 3
        ok = Eb[Eb.n_data >= 30]
        rho = float(spearmanr(np.arange(len(ok)), ok.conf_ze_data).correlation) if len(ok) >= 3 else np.nan
        we.append(dict(detector=det, net_set=ns, quiet_src=",".join(sorted(g.quiet_src.unique())), n_quiet_days=g.date.nunique(),
                       n_data=int((g.series == "data").sum()), n_ctl=len(ct), n_data_P_ge_05=len(hi), diff=d[0], lo95=d[1], hi95=d[2],
                       spearman_bins=rho, pass_=bool(np.isfinite(d[1]) and d[1] > 0 and np.isfinite(rho) and rho > 0)))
    WE = pd.DataFrame(we); WE.to_csv(f"{T7}/withheld_energy.csv", index=False); pd.concat(wb).to_csv(f"{T7}/withheld_energy_bins.csv", index=False)
    # ---------------- station-set invariance (P4-2 runs, 509 E0 PNSN-quiet days 2014-01..2017-02)
    kb = pd.read_csv(f"{Q7.P6}/tables/k_factors_bw.csv"); kf = pd.read_csv(f"{G.REFINE}/rates/k_factors.csv")
    setk = {s: float(g.set_mean_k.iloc[0]) for s, g in kb.groupby("sset")}
    setk["B3"] = float(kf[kf.sset == "S14"].set_mean_k.iloc[0]); setk["S15fix"] = float(kf[kf.sset == "S15"].set_mean_k.iloc[0])
    inv_runs = ["p4e0_CN4", "p4e0_B1", "p4e0_B2", "p4e0_B3", "p4e0_S15fix"]
    FDi = Q7.load_fd(inv_runs)
    md = FDi[FDi.observable & (FDi.series == "data")].groupby(["run", "date"]).size()
    agg = {}
    for r in inv_runs:
        P = pd.read_parquet(f"{Q7.P6}/P/{r}_P.parquet")
        g = P.groupby(["date", "series"]).P.sum().unstack(fill_value=0.0)
        g["m"] = md.reindex([(r, d) for d in g.index]).to_numpy()
        agg[r] = g
    ref = agg["p4e0_S15fix"]; iv, ivm = [], []
    for r in inv_runs[:-1]:
        s_ = r.split("_")[1]; g = agg[r]; days = sorted(set(g.index) & set(ref.index))
        a1 = g.loc[days]; a0 = ref.loc[days]
        mo = pd.Series(days).str[:7].to_numpy()
        x1 = setk[s_] * a1["data"]; x0 = setk["S15fix"] * a0["data"]
        M = pd.DataFrame(dict(set_=x1.groupby(mo).sum(), ref=x0.groupby(mo).sum())); M["ratio"] = M.set_ / M.ref
        e1 = (a1["data"] - a1.get("c1", 0)); e0 = (a0["data"] - a0.get("c1", 0))
        Me = pd.DataFrame(dict(set_=(setk[s_] * e1).groupby(mo).sum(), ref=(setk["S15fix"] * e0).groupby(mo).sum())); Me["ratio"] = Me.set_ / Me.ref
        kq = (e1.sum() / a1.m.sum()) / (e0.sum() / a0.m.sum())
        Kv = KC.loc[s_, "K"] if s_ in KC.index else np.nan
        d = dict(net_set=s_, n_days=len(days), n_months=len(M), k_set_bostock=setk[s_], k_S15_bostock=setk["S15fix"],
                 pooled_ratio=float(M.set_.sum() / M.ref.sum()), median_monthly_ratio=float(M.ratio.median()),
                 pooled_ratio_excess=float(Me.set_.sum() / Me.ref.sum()), median_monthly_ratio_excess=float(Me.ratio.median()),
                 quiet_k_excess_set_over_S15=float(kq), ladder_1_over_K=float(1 / Kv) if Kv > 0 else np.nan)
        d["pass_"] = bool(0.8 <= d["pooled_ratio"] <= 1.2 and 0.8 <= d["median_monthly_ratio"] <= 1.2)
        d["quiet_k_vs_ladder_agree_x2"] = bool(np.isfinite(d["ladder_1_over_K"]) and kq > 0 and max(kq, d["ladder_1_over_K"]) / min(kq, d["ladder_1_over_K"]) <= 2)
        iv.append(d); ivm.append(M.assign(net_set=s_).reset_index().rename(columns={"index": "month"}))
    IV = pd.DataFrame(iv); IV.to_csv(f"{T7}/invariance.csv", index=False); pd.concat(ivm).to_csv(f"{T7}/invariance_monthly.csv", index=False)
    # ---------------- G6
    g6 = {}
    csa = CS[(CS.year == "all") & (CS.days == "all_days")].set_index(["detector", "net_set"])
    for (det, ns) in csa.index:
        ctl_ok = bool(csa.loc[(det, ns), "ctl_over_data"] <= 0.05 and abs(csa.loc[(det, ns), "spearman_monthly_ctl_vs_time"] or 0) < 0.3)
        w = WE[(WE.detector == det) & (WE.net_set == ns)]
        en_ok = bool(len(w) and w.pass_.iloc[0])
        iv_ = IV[IV.net_set == ns]
        inv_ok = bool(len(iv_) and iv_.pass_.iloc[0]) if ns not in ("S15", "ERA") else False
        g6[f"{det}:{ns}"] = dict(control_sumP_pass=ctl_ok, ctl_over_data=float(csa.loc[(det, ns), "ctl_over_data"]), withheld_energy_pass=en_ok,
                                 invariance_pass=inv_ok, invariance_note="S15/B4 = reference (forward invariance failed: S13 vs S15)" if ns == "S15" else
                                 ("no era run on 2014-17 days: not testable" if ns == "ERA" else ""),
                                 inter_ets_mode="rate" if (ctl_ok and en_ok and inv_ok) else "control_subtracted_bracket")
    res["g6"] = g6
    json.dump(res, open(f"{T7}/g6.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
    print(SC[SC.day_class == "between"].round(3).to_string(index=False)); print(FR.round(3).to_string(index=False)); print(DD.round(3).to_string(index=False))
    print(EE.round(3).to_string(index=False)); print(CS[CS.year == "all"].round(3).to_string(index=False)); print(WE.round(3).to_string(index=False))
    print(IV.round(3).to_string(index=False)); print(json.dumps(g6, indent=1, default=float))


if __name__ == "__main__":
    main()
