#!/usr/bin/env python
"""PNET check, step 3: real-share table from tables/per_candidate.parquet (see pc_analyze.py for definitions)."""
import numpy as np, pandas as pd
OD = "/home/jovyan/tremorferometry/lfe_extension_wa/results/wa_prod/quietwave/pnet_check"; TB = f"{OD}/tables"
C = pd.read_parquet(f"{TB}/per_candidate.parquet")
RULES = {"PNET>=0.75": ["both", "pnet_only"], "WAVE>=0.75": ["both", "wave_only"], "both": ["both"],
         "PNET only (WAVE<0.75)": ["pnet_only"], "WAVE only (PNET<0.75)": ["wave_only"]}
COLS = ["out_cc", "out_coda", "out_P", "e_snr"]
SER = {"data": "data", "random": "random times", "c1": "scrambled control", "rev": "reversed control"}
rng = np.random.default_rng(11)
T = C[(C.cls == "tremor30") & (C.series == "data")]


def fam_means(reg, sels):
    h = T[(T.reg == reg) & T.sel.isin(sels)]
    cnt = h.groupby("fid").size()
    return h.groupby("fid")[COLS].mean()[cnt >= (10 if reg == "NW" else 5)]


def frac(h, col, tm, nb=500):
    h = h[h.fid.isin(tm.index) & np.isfinite(h[col])]
    if len(h) < 5:
        return np.nan, np.nan, np.nan, len(h)
    den = tm.loc[h.fid, col].to_numpy(); num = h[col].to_numpy()
    if col in HB and HB[col] is not None:            # strength-matched: x the both-set ratio in the candidate's net_cc bin
        den = den * HB[col].reindex(pd.cut(h.net_cc, BINS)).to_numpy()
    est = num.sum() / den.sum()
    days = h.date.to_numpy(); ud, inv = np.unique(days, return_inverse=True)
    nd = np.bincount(inv, num, len(ud)); dd = np.bincount(inv, den, len(ud))
    bs = []
    for _ in range(nb):
        k = rng.integers(0, len(ud), len(ud)); bs.append(nd[k].sum() / dd[k].sum())
    return est, np.quantile(bs, 0.05), np.quantile(bs, 0.95), len(h)


BINS = [0, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.25, 0.3, 1.0]
HB = {}


def netcc_curve(reg, col):
    """tremor-day BOTH rows: mean of (out / family mean) per net_cc bin = how a real event's score grows with net cc."""
    h = T[(T.reg == reg) & (T.sel == "both")]; fm = fam_means(reg, ["both"])
    h = h[h.fid.isin(fm.index) & np.isfinite(h[col])]
    r = h[col].to_numpy() / fm.loc[h.fid, col].to_numpy()
    g = pd.Series(r).groupby(pd.cut(h.net_cc, BINS).to_numpy(), observed=False).mean()
    return g


rows = []
for reg, classes in (("NW", ["quiet_cov", "unknown", "tremor30"]), ("SW", ["tremor30"]), ("ArmB", ["tremor30"])):
    for cls in classes:
        for rule, sels in RULES.items():
            same = fam_means(reg, sels); both = fam_means(reg, ["both"])
            for ser in ("data", "random", "c1", "rev"):
                h = C[(C.reg == reg) & (C.cls == cls) & (C.series == ser) & C.sel.isin(sels)]
                if not len(h):
                    continue
                d = dict(reg=reg, day_class=cls, rule=rule, series=SER[ser], n=len(h), n_days=h.date.nunique(),
                         n_families=h.fid.nunique(), net_cc_med=h.net_cc.median(), out_cc_mean=h.out_cc.mean(),
                         out_coda_mean=h.out_coda.mean(), out_P_mean=h.out_P.mean())
                for nm, tm in (("same_rule", same), ("vs_both", both), ("vs_both_netcc", both)):
                    for col in ("out_cc", "out_coda", "e_snr"):
                        HB.clear()
                        if nm == "vs_both_netcc":
                            if reg != "NW" or col == "e_snr":
                                continue
                            HB[col] = netcc_curve(reg, col)
                        e, lo, hi, nn = frac(h, col, tm)
                        tag = f"{nm}_{col}"
                        d[tag] = e; d[tag + "_p05"] = lo; d[tag + "_p95"] = hi
                    d[f"{nm}_n_matched"] = nn
                rows.append(d)
R = pd.DataFrame(rows)
R.to_csv(f"{TB}/real_share_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
show = ["reg", "day_class", "rule", "series", "n", "net_cc_med", "out_cc_mean", "out_coda_mean", "same_rule_out_cc", "same_rule_out_cc_p05",
        "same_rule_out_cc_p95", "same_rule_out_coda", "same_rule_out_coda_p05", "same_rule_out_coda_p95", "vs_both_out_cc",
        "vs_both_netcc_out_cc", "vs_both_netcc_out_cc_p05", "vs_both_netcc_out_cc_p95", "vs_both_netcc_out_coda"]
print(R[R.series.isin(["data", "random times"]) | (R.n >= 20)][show].round(3).to_string())
