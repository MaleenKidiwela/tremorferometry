#!/usr/bin/env python
"""bw_p7_rates -- P7 item 6: per-family daily ΣP-excess rate tables 2003-2013 in the forward format
(lfe_extension/results/refine/rates/famday.parquet, monthly_family.csv, monthly_total.csv), uncorrected and completeness-ladder
k-corrected (rules: results/backward/p7/P7.md).  Needs the catalogue, p7/tables/k_chain.csv, p7/tables/g6.json, p6/tables/k_factors_bw.csv.
Writes NEW files under data/catalogue/ (refuses to overwrite): svi_rates_v1_4_backward_famday.parquet, _monthly_family.csv,
_monthly_total.csv, _monthly_total_by_detector.csv, _k_factors.csv, _meta.json.
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_p7_common as Q7
warnings.filterwarnings("ignore")

ROOT = Q7.ROOT
OUT = f"{ROOT}/data/catalogue/svi_rates_v1_4_backward"
TLO = {"CN4": 8.79, "B1": 8.69, "B2": 8.75, "B3": 8.75, "S15": 8.75, "ERA": 8.25}
THI = {"CN4": 9.82, "B1": 9.56, "B2": 8.75, "B3": 8.75, "S15": 8.75, "ERA": 8.25}


def main():
    if os.path.exists(f"{OUT}_famday.parquet") and "--force" not in sys.argv:
        raise SystemExit("[rates] output exists; refusing to overwrite")
    A = pd.read_parquet(Q7.CAT, columns=["run", "date", "series", "fam", "P", "P_same_WAVE_foldsd", "net_mad", "net_set", "detector",
                                         "multi_site_z2", "n_sta_day"])
    A["Tlo"] = A.net_set.map(TLO); A["Thi"] = A.net_set.map(THI)
    A["varB"] = A.P * (1 - A.P); A["core"] = (A.net_mad >= A.Tlo).astype(int); A["core_hi"] = (A.net_mad >= A.Thi).astype(int)
    A["P_ms"] = A.P * A.multi_site_z2.astype(bool)
    g = A.groupby(["run", "fam", "date", "series"]).agg(n=("P", "size"), sumP=("P", "sum"), varB=("varB", "sum"), n_core=("core", "sum"),
                                                        n_core_Thi=("core_hi", "sum"), foldsd_sum=("P_same_WAVE_foldsd", "sum"),
                                                        sumP_ms=("P_ms", "sum")).reset_index()
    W = g.pivot_table(index=["run", "fam", "date"], columns="series", values=["n", "sumP", "varB", "n_core", "n_core_Thi", "foldsd_sum", "sumP_ms"],
                      fill_value=0.0)
    W.columns = [f"{a}_{b}" for a, b in W.columns]; W = W.reset_index()
    FD = Q7.load_fd(Q7.PROD)
    ob = FD.pivot_table(index=["run", "fam", "date"], columns="series", values="observable", aggfunc="max", fill_value=False)
    ob.columns = [f"m_{c}" for c in ob.columns]; ob = ob.reset_index()
    W = ob.merge(W, on=["run", "fam", "date"], how="outer")
    for c in ["n_data", "n_c1", "sumP_data", "sumP_c1", "varB_data", "n_core_data", "n_core_c1", "n_core_Thi_data", "n_core_Thi_c1",
              "foldsd_sum_data", "sumP_ms_data", "sumP_ms_c1"]:
        if c not in W: W[c] = 0.0
        W[c] = W[c].fillna(0.0)
    for c in ("m_data", "m_c1"):
        W[c] = W[c].fillna(False).astype(bool).astype(int)
    W = W.drop(columns=[c for c in ("varB_c1", "foldsd_sum_c1") if c in W])
    W["detector"] = np.where(W.run.str.startswith("p5e_"), "era", "frozen")
    W["net_set"] = np.where(W.detector == "era", "ERA", W.date.map(Q7.set_of_date))
    W["sset"] = W.net_set
    inera = (W.date >= Q7.ERA_DET[0]) & (W.date <= Q7.ERA_DET[1])
    W["primary"] = np.where(W.detector == "era", inera, ~inera)
    nsd = A.groupby(["run", "date"]).n_sta_day.max()
    W["n_sta_day"] = [nsd.get((r, d), np.nan) for r, d in zip(W.run, W.date)]
    # Bostock k per (family, set)
    kb = pd.read_csv(f"{Q7.P6}/tables/k_factors_bw.csv"); kf = pd.read_csv(f"{ROOT}/lfe_extension/results/refine/rates/k_factors.csv")
    kb = kb[kb.sset.isin(["CN4", "B1", "B2", "ERA"])]
    kf = kf[kf.sset.isin(["S14", "S15"])].assign(sset=lambda x: x.sset.map({"S14": "B3", "S15": "S15"}))
    K = pd.concat([kb, kf], ignore_index=True); K.to_csv(f"{OUT}_k_factors.csv", index=False)
    Kd = {(int(r.fam), r.sset): (r.k, r.k_lo, r.k_hi) for r in K.itertuples()}
    Km = {s: float(g_.set_mean_k.iloc[0]) for s, g_ in K.groupby("sset")}
    def kk(f, s, i):
        v = Kd.get((int(f), s))
        return v[i] if v is not None and np.isfinite(v[i]) else Km.get(s, np.nan)
    W["k"] = [kk(f, s, 0) for f, s in zip(W.fam, W.sset)]; W["k_lo"] = [kk(f, s, 1) for f, s in zip(W.fam, W.sset)]
    W["k_hi"] = [kk(f, s, 2) for f, s in zip(W.fam, W.sset)]
    W["excess"] = W.sumP_data - W.sumP_c1; W["excess_ms"] = W.sumP_ms_data - W.sumP_ms_c1
    W["rate"] = W.excess * W.k; W["rate_sd"] = np.sqrt(W.varB_data + W.foldsd_sum_data ** 2) * W.k
    DF = Q7.day_flags()
    for c in ("ets", "between", "env_quiet", "pnsn_quiet_b"):
        W[c] = W.date.map(DF[c]).fillna(False).astype(bool)
    W["quiet_env3_CN4"] = W.env_quiet; W["quiet_pnsn"] = W.pnsn_quiet_b
    W["quiet"] = np.where(W.date >= "2010-01-10", W.quiet_pnsn, W.quiet_env3_CN4)
    W["bracket_lo"] = np.where(W.quiet, W.excess * W.k_lo, np.nan)
    W["bracket_hi"] = np.where(W.quiet, W.n_core_data - W.n_core_c1, np.nan)
    W["month"] = W.date.str[:7]; W["v11_cat_core"] = 0.0; W["v11_cat_ext"] = 0.0
    bd = pd.read_parquet(Q7.CAT, columns=["run", "date", "burst_day"]).drop_duplicates(["run", "date"]).set_index(["run", "date"]).burst_day
    W["burst_day"] = [bool(bd.get((r, d), False)) for r, d in zip(W.run, W.date)]
    g6 = json.load(open(f"{Q7.P7}/tables/g6.json"))["g6"]
    W["inter_ets_mode"] = [g6.get(f"{d}:{s}", {}).get("inter_ets_mode", "control_subtracted_bracket") for d, s in zip(W.detector, W.net_set)]
    KC = Q7.k_chain("all")
    W["K_set"] = W.net_set.map(KC.K); W["K_lo95"] = W.net_set.map(KC.K_lo95); W["K_hi95"] = W.net_set.map(KC.K_hi95)
    W["excess_kc"] = W.excess * W.K_set
    W["k_S15"] = [kk(f, "S15", 0) for f in W.fam]
    W["rate_kc"] = W.excess_kc * W.k_S15
    W["era"] = W.date.map(lambda d: {"CN4": "B0", "B1": "B1", "B2": "B2", "B3": "B3", "S15": "B4"}[Q7.set_of_date(d)])
    cols = ["fam", "date", "n_data", "n_c1", "sumP_data", "sumP_c1", "varB_data", "n_core_data", "n_core_c1", "foldsd_sum_data", "sset",
            "n_sta_day", "k", "k_lo", "k_hi", "excess", "rate", "rate_sd", "quiet", "bracket_lo", "bracket_hi", "month", "v11_cat_core",
            "v11_cat_ext", "inter_ets_mode",
            "run", "detector", "primary", "era", "net_set", "m_data", "m_c1", "n_core_Thi_data", "n_core_Thi_c1", "quiet_env3_CN4", "quiet_pnsn",
            "ets", "between", "burst_day", "sumP_ms_data", "sumP_ms_c1", "excess_ms", "K_set", "K_lo95", "K_hi95", "excess_kc", "k_S15", "rate_kc"]
    W = W[cols].sort_values(["date", "detector", "fam"]).reset_index(drop=True)
    W.to_parquet(f"{OUT}_famday.parquet", index=False)
    # monthly tables (primary detector = one series per month, joinable to the forward monthly tables)
    def monthly(X, by):
        days_m = X.groupby(by[:-1] + ["month"] if len(by) > 1 else ["month"]).date.nunique() if False else None
        return days_m
    P = W[W.primary]
    days_m = P.groupby("month").date.nunique()
    MF = P.groupby(["month", "fam"]).agg(sumP=("sumP_data", "sum"), sumP_ctl=("sumP_c1", "sum"), rate=("rate", "sum"),
                                         rate_var=("rate_sd", lambda s: float((s ** 2).sum())), v11_cat_core=("v11_cat_core", "sum"),
                                         v11_cat_ext=("v11_cat_ext", "sum"), excess=("excess", "sum"), excess_kc=("excess_kc", "sum"),
                                         rate_kc=("rate_kc", "sum"), m_data=("m_data", "sum")).reset_index()
    qd = P.drop_duplicates("date").groupby("month").quiet.sum()
    MF["days"] = MF.month.map(days_m); MF["rate_per_day"] = MF.rate / MF.days; MF["rate_per_day_sd"] = np.sqrt(MF.rate_var) / MF.days
    MF["rate_kc_per_day"] = MF.rate_kc / MF.days; MF["excess_per_obs_famday"] = MF.excess / MF.m_data.replace(0, np.nan)
    MF["quiet_days"] = MF.month.map(qd).fillna(0).astype(int)
    MF.to_csv(f"{OUT}_monthly_family.csv", index=False)
    def mt(X):
        dm = X.groupby("month").date.nunique()
        T = X.groupby("month").agg(sumP=("sumP_data", "sum"), sumP_ctl=("sumP_c1", "sum"), rate=("rate", "sum"),
                                   rate_var=("rate_sd", lambda s: float((s ** 2).sum())), bracket_lo=("bracket_lo", "sum"),
                                   bracket_hi=("bracket_hi", "sum"), v11_cat_core=("v11_cat_core", "sum"), excess=("excess", "sum"),
                                   excess_ms=("excess_ms", "sum"), excess_kc=("excess_kc", "sum"), rate_kc=("rate_kc", "sum"),
                                   m_data=("m_data", "sum"), sset=("sset", lambda s: ",".join(sorted(set(s.dropna())))),
                                   detector=("detector", lambda s: ",".join(sorted(set(s))))).reset_index()
        T["days"] = T.month.map(dm); T["rate_per_day"] = T.rate / T.days; T["rate_per_day_sd"] = np.sqrt(T.rate_var) / T.days
        T["rate_kc_per_day"] = T.rate_kc / T.days; T["excess_per_obs_famday"] = T.excess / T.m_data.replace(0, np.nan)
        T["excess_kc_per_obs_famday"] = T.excess_kc / T.m_data.replace(0, np.nan)
        T["excess_noburst"] = T.month.map(X[~X.burst_day].groupby("month").excess.sum()).fillna(0.0)
        q = X.drop_duplicates("date").groupby("month").quiet.sum()
        T["n_quiet_days"] = T.month.map(q).fillna(0).astype(int); T["n_quiet_famdays"] = T.n_quiet_days * X.fam.nunique()
        T["inter_ets_mode"] = T.month.map(X.groupby("month").inter_ets_mode.agg(lambda s: ",".join(sorted(set(s)))))
        return T
    mt(P).to_csv(f"{OUT}_monthly_total.csv", index=False)
    pd.concat([mt(W[W.detector == d]).assign(series=d) for d in ("frozen", "era")]).to_csv(f"{OUT}_monthly_total_by_detector.csv", index=False)
    meta = dict(created_utc=pd.Timestamp.now("UTC").isoformat(), rows=len(W), rules="svi_science/results/backward/p7/P7.md",
                forward_format="lfe_extension/results/refine/rates/famday.parquet (+ extra columns)",
                uncorrected="excess = ΣP_data − ΣP_control; rate = excess × k(family, set) (Bostock k of the set, forward definition)",
                k_corrected="excess_kc = excess × K_set (completeness ladder chained to the S15/B4 level); rate_kc = excess_kc × k(family, S15)",
                join_to_forward="forward famday 2014-2026 has the same first 25 columns; use primary == True rows for one series per day",
                K_chain=Q7.k_chain("all")[["K", "K_lo95", "K_hi95"]].to_dict("index"))
    json.dump(meta, open(f"{OUT}_meta.json", "w"), indent=1, default=float)
    print(json.dumps(meta, indent=1, default=float)[:3000])
    print(mt(P).groupby(mt(P).month.str[:4])[["excess", "excess_kc", "rate", "rate_kc", "days"]].sum().round(1).to_string())


if __name__ == "__main__":
    main()
