#!/usr/bin/env python
"""r23b_f6_rates (2026-10-11) -- Finding 6a: ETS-episode sizes per family 2003-2026 from the Bostock-equivalent `rate`.
Rules: svi_science/results/r23b/PLAN.md §3a (pre-registered). No between-episode rates are used.

    python svi_science/scripts/r23b_f6_rates.py
"""
from __future__ import annotations
import itertools, json, time
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
TB = f"{ROOT}/svi_science/results/r23b/tables"
OUT = f"{TB}/f6"
R23 = f"{ROOT}/svi_science/results/r23/tables"
BWR = f"{ROOT}/data/catalogue/svi_rates_v1_4_backward_famday.parquet"
FWR = f"{ROOT}/lfe_extension/results/refine/rates/famday.parquet"
NB = 2000
rng = np.random.default_rng(606)
import os; os.makedirs(OUT, exist_ok=True)
LOGF = open(f"{ROOT}/svi_science/results/r23b/logs/r23b_f6_rates.log", "a")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()


def ts_slope(t, a):
    t = np.asarray(t, float); a = np.asarray(a, float)
    i, j = np.triu_indices(len(t), 1); dt = t[j] - t[i]; da = a[j] - a[i]
    m = np.abs(dt) > 1e-9
    return float(np.median(da[m] / dt[m])) if m.any() else np.nan


def table(F, W):
    rb = pd.read_parquet(BWR, columns=["fam", "date", "primary", "rate", "rate_kc", "sset"]); rb = rb[rb.primary].drop(columns="primary")
    rf = pd.read_parquet(FWR, columns=["fam", "date", "rate", "sset"]).assign(rate_kc=lambda x: x.rate)
    R = pd.concat([rb, rf[rb.columns]], ignore_index=True); R = R[R.fam.isin(F.fid)]
    V = pd.read_parquet(f"{R23}/events_2003_2026.parquet", columns=["series", "fam", "date"]); V = V[V.series == "data"]
    B = pd.read_parquet(f"{ROOT}/svi_science/results/p3_bostock/bostock_events_mapped.parquet", columns=["fid", "date"])
    B = B[B.fid.isin(F.fid)].rename(columns={"fid": "fam"})
    X = pd.read_parquet(f"{R23}/famday_tremor.parquet")
    rows = []
    for r in W.itertuples():
        dd = set(pd.date_range(r.start, r.end).strftime("%Y-%m-%d"))
        g = R[R.date.isin(dd)]
        s = g.groupby("fam").agg(S=("rate", "sum"), S_kc=("rate_kc", "sum"), frac_S15=("sset", lambda x: float((x == "S15").mean())))
        n = V[V.date.isin(dd)].groupby("fam").size().rename("n_once")
        nb = B[B.date.isin(dd)].groupby("fam").size().rename("n_bostock")
        xe = X[X.date.isin(dd) & (X.src == "env") & (X.n30 > 0)].groupby("fam").size().rename("expdays_env")
        xp = X[X.date.isin(dd) & (X.src == "pnsn") & (X.n30 > 0)].groupby("fam").size().rename("expdays_pnsn")
        d = pd.DataFrame(index=pd.Index(F.fid.astype(int), name="fam")).join([s, n, nb, xe, xp]).fillna(0).reset_index()
        d["episode"] = r.episode; rows.append(d)
    T = pd.concat(rows, ignore_index=True)
    T["S"] = T.S.clip(lower=0); T["S_kc"] = T.S_kc.clip(lower=0)
    return T


def lratio(M, early, late, fams):
    """median over families of log10(mean late / mean early)."""
    a = M.loc[fams, early].mean(1).to_numpy(); b = M.loc[fams, late].mean(1).to_numpy()
    ok = (a > 0) & (b > 0)
    return float(np.median(np.log10(b[ok] / a[ok]))) if ok.any() else np.nan


def tot_ratio(M, early, late):
    return float(np.log10(M[late].sum(0).mean() / M[early].sum(0).mean()))


def compare(T, col, early, late, nperm=10000, minev=20.0, exact=False):
    M = T.pivot_table(index="fam", columns="episode", values=col, aggfunc="sum").fillna(0)
    fams = M.index[(M[early].sum(1) >= minev) & (M[late].sum(1) >= minev)]
    obs = lratio(M, early, late, fams); tobs = tot_ratio(M, early, late)
    eps = list(early) + list(late)
    if exact:
        combos = list(itertools.combinations(range(len(eps)), len(late)))
        nul = []; tnul = []
        for c in combos:
            l2 = [eps[k] for k in c]; e2 = [eps[k] for k in range(len(eps)) if k not in c]
            nul.append(lratio(M, e2, l2, fams)); tnul.append(tot_ratio(M, e2, l2))
    else:
        nul, tnul = [], []
        for _ in range(nperm):
            p = rng.permutation(len(eps)); l2 = [eps[k] for k in p[:len(late)]]; e2 = [eps[k] for k in p[len(late):]]
            nul.append(lratio(M, e2, l2, fams)); tnul.append(tot_ratio(M, e2, l2))
    nul = np.array(nul); tnul = np.array(tnul)
    p = float((np.sum(np.abs(nul) >= abs(obs) - 1e-12)) / len(nul)) if exact else float((np.sum(np.abs(nul) >= abs(obs)) + 1) / (len(nul) + 1))
    pt = float((np.sum(np.abs(tnul) >= abs(tobs) - 1e-12)) / len(tnul)) if exact else float((np.sum(np.abs(tnul) >= abs(tobs)) + 1) / (len(tnul) + 1))
    bs, bt = [], []
    fa = np.array(fams)
    for _ in range(NB):
        e2 = list(rng.choice(early, len(early))); l2 = list(rng.choice(late, len(late)))
        f2 = rng.choice(fa, len(fa))
        Mb = M.loc[f2]
        a = np.array([Mb[e].to_numpy() for e in e2]).mean(0); b = np.array([Mb[l].to_numpy() for l in l2]).mean(0)
        ok = (a > 0) & (b > 0)
        bs.append(np.median(np.log10(b[ok] / a[ok])) if ok.any() else np.nan)
        bt.append(np.log10(M[l2].sum(0).mean() / M[e2].sum(0).mean()))
    ci = [float(np.nanquantile(bs, 0.025)), float(np.nanquantile(bs, 0.975))]
    cit = [float(np.nanquantile(bt, 0.025)), float(np.nanquantile(bt, 0.975))]
    return dict(col=col, early=list(early), late=list(late), n_families=int(len(fams)), median_log10_ratio=obs,
                ratio=10 ** obs, ci95_log10=ci, ci95_ratio=[10 ** ci[0], 10 ** ci[1]], p_perm=p,
                total_log10_ratio=tobs, total_ratio=10 ** tobs, total_ci95_ratio=[10 ** cit[0], 10 ** cit[1]], p_perm_total=pt,
                n_perm=int(len(nul)), exact=exact)


def main():
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    W = pd.read_csv(f"{TB}/windows.csv")
    T = table(F, W)
    T.to_parquet(f"{OUT}/fam_episode_sizes.parquet", index=False)
    big = W[W.big].copy()
    fs15 = T.groupby("episode").frac_S15.mean()
    big["frac_S15"] = big.episode.map(fs15)
    S15f = big[big.seg == "S15f"].episode.tolist()
    FWD_S15 = big[(big.seg == "FWD") & (big.start < "2019-01-01") & (big.frac_S15 >= 0.9)].episode.tolist()
    early, late = big[big.start < "2014-01-01"].episode.tolist(), big[big.start >= "2014-01-01"].episode.tolist()
    log(f"S15f big {S15f}; FWD-S15 big {FWD_S15} (frac S15 {big.set_index('episode').frac_S15.reindex(FWD_S15).round(2).tolist()})")
    R = {}
    R["L4L_rate"] = compare(T, "S", S15f, FWD_S15, exact=True)
    R["L4L_once"] = compare(T, "n_once", S15f, FWD_S15, exact=True, minev=5)
    R["DEC_rate"] = compare(T, "S", early, late)
    R["DEC_rate_kc"] = compare(T, "S_kc", early, late)
    R["S15f_vs_allFWD_rate"] = compare(T, "S", S15f, late)
    for k, v in R.items():
        log(f"{k}: n_fam {v['n_families']} ratio {v['ratio']:.3f} CI {np.round(v['ci95_ratio'], 3).tolist()} p {v['p_perm']:.4f} | "
            f"total ratio {v['total_ratio']:.3f} CI {np.round(v['total_ci95_ratio'], 3).tolist()} p {v['p_perm_total']:.4f}")
    # episode level
    ep = T.groupby("episode").agg(S=("S", "sum"), S_kc=("S_kc", "sum"), n_once=("n_once", "sum"), n_bostock=("n_bostock", "sum"),
                                  expdays_env=("expdays_env", "sum"), expdays_pnsn=("expdays_pnsn", "sum"))
    E = W.set_index("episode").join(ep, rsuffix="_t")
    E["S_per_pnsn_expday"] = np.where(E.expdays_pnsn > 0, E.S / E.expdays_pnsn, np.nan)
    E["S_per_env_expday"] = np.where(E.expdays_env > 0, E.S / E.expdays_env, np.nan)
    E["frac_S15"] = E.index.map(fs15)
    E.reset_index().to_csv(f"{OUT}/episode_sizes.csv", index=False)
    Bg = E[E.big]
    a, b = Bg.loc[early, "S"], Bg.loc[late, "S"]
    R["episode_level"] = dict(median_S_2003_13=float(a.median()), median_S_2014_26=float(b.median()),
                              mw_p=float(stats.mannwhitneyu(a, b).pvalue),
                              median_S_kc_2003_13=float(Bg.loc[early, "S_kc"].median()),
                              mw_p_kc=float(stats.mannwhitneyu(Bg.loc[early, "S_kc"], Bg.loc[late, "S_kc"]).pvalue))
    M = E[E.major]
    y1 = 11.0; y2 = (pd.Timestamp("2026-09-15") - pd.Timestamp("2014-01-01")).days / 365.25
    R["release_per_year_majors"] = dict(S_2003_13=float(M[M.start < "2014-01-01"].S.sum() / y1), S_2014_26=float(M[M.start >= "2014-01-01"].S.sum() / y2),
                                        S_kc_2003_13=float(M[M.start < "2014-01-01"].S_kc.sum() / y1),
                                        S15f_2010_13=float(M[M.seg == "S15f"].S.sum() / ((pd.Timestamp("2014-01-01") - pd.Timestamp("2009-08-27")).days / 365.25)),
                                        S_2014_18=float(M[(M.start >= "2014-01-01") & (M.start < "2019-01-01")].S.sum() / 5.0),
                                        big_S_2003_13=float(Bg.loc[early, "S"].sum() / y1), big_S_2014_26=float(Bg.loc[late, "S"].sum() / y2))
    # tremor-normalised (never across the env/PNSN switch)
    pn = Bg[Bg.start >= "2010-01-06"]
    s15 = pn.loc[S15f, "S_per_pnsn_expday"]; f15 = pn.loc[FWD_S15, "S_per_pnsn_expday"]; fall = pn.loc[late, "S_per_pnsn_expday"]
    R["trem_norm_pnsn"] = dict(S15f=s15.round(2).to_dict(), FWD_S15=f15.round(2).to_dict(), FWD_all_median=float(fall.median()),
                               ratio_FWDS15_S15f=float(f15.mean() / s15.mean()), ratio_FWDall_S15f=float(fall.median() / s15.median()),
                               mw_p_FWDS15=float(stats.mannwhitneyu(s15, f15).pvalue), mw_p_FWDall=float(stats.mannwhitneyu(s15, fall).pvalue))
    ev = Bg.loc[early]
    R["trem_norm_env_2003_13"] = dict(values=ev.S_per_env_expday.round(2).to_dict(),
                                      spearman_year=[float(x) for x in stats.spearmanr(pd.to_datetime(ev.start).dt.year, ev.S_per_env_expday)])
    # Bostock check: his per-family counts over his 10 windows (fixed method); our S vs his counts
    bost = W[W.n_bostock > 0].sort_values("start")
    bost = bost[bost.start < "2014-01-01"]
    Mb = T.pivot_table(index="fam", columns="episode", values="n_bostock", aggfunc="sum")[bost.episode]
    tt = (pd.to_datetime(bost.start) - pd.Timestamp("2003-01-01")).dt.days.to_numpy() / 3652.5
    good = Mb.index[(Mb >= 1).sum(1) >= 6]
    def slopes(order):
        return np.array([ts_slope(tt, np.log10(Mb.loc[f, order].to_numpy() + 1)) for f in good])
    sl = slopes(list(bost.episode)); obs = float(np.nanmean(sl))
    nul = [float(np.nanmean(slopes(list(rng.permutation(bost.episode))))) for _ in range(2000)]
    bsm = [float(np.nanmean(rng.choice(sl, len(sl)))) for _ in range(2000)]
    R["bostock_trend"] = dict(episodes=list(bost.episode), n_families=int(len(good)), mean_slope_log10_per_decade=obs,
                              ci95=[float(np.quantile(bsm, 0.025)), float(np.quantile(bsm, 0.975))],
                              p_order=float((np.sum(np.abs(nul) >= abs(obs)) + 1) / 2001),
                              his_totals=Mb.sum(0).astype(int).to_dict())
    vv = T[T.episode.isin(bost.episode) & (T.n_bostock + T.S > 0)].merge(W[["episode", "seg"]], on="episode")
    R["rate_vs_bostock_spearman"] = {sg: float(stats.spearmanr(g.S, g.n_bostock)[0]) for sg, g in vv.groupby("seg")}
    # family shares, decade vs decade
    Ms = T[T.episode.isin(early + late)].pivot_table(index="fam", columns="episode", values="S").fillna(0)
    sh = Ms / Ms.sum(0)
    rho = stats.spearmanr(sh[early].mean(1), sh[late].mean(1))
    R["family_share_spearman_decades"] = dict(rho=float(rho[0]), p=float(rho[1]))
    # verdict (PLAN §3a)
    L = R["L4L_rate"]; D_ = R["DEC_rate"]; Dk = R["DEC_rate_kc"]
    inc = lambda c: c[0] <= 1.0 <= c[1]
    if inc(L["ci95_ratio"]) and (inc(D_["ci95_ratio"]) or inc(Dk["ci95_ratio"])):
        v = "holds (no change)"
    elif (not inc(L["ci95_ratio"])) and L["p_perm"] < 0.01 and not inc(D_["ci95_ratio"]) and not inc(Dk["ci95_ratio"]) and \
            np.sign(L["median_log10_ratio"]) == np.sign(D_["median_log10_ratio"]) == np.sign(Dk["median_log10_ratio"]):
        v = "changes"
    else:
        v = "undecided"
    R["verdict_3a"] = v
    log(f"episode level {R['episode_level']}\nrelease/yr {R['release_per_year_majors']}\ntrem PNSN {R['trem_norm_pnsn']}\n"
        f"trem env {R['trem_norm_env_2003_13']}\nbostock {R['bostock_trend']}\nrate vs bostock {R['rate_vs_bostock_spearman']}\n"
        f"shares {R['family_share_spearman_decades']}\nVERDICT 3a: {v}")
    json.dump(R, open(f"{OUT}/f6_rates.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def posthoc():
    """POST-HOC (PLAN Revision 2): like-for-like without the two low-LFE ETS; all majors per tremor family-exposure-day."""
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    W = pd.read_csv(f"{TB}/windows.csv")
    T = pd.read_parquet(f"{OUT}/fam_episode_sizes.parquet")
    E = pd.read_csv(f"{OUT}/episode_sizes.csv")
    R = {}
    S15f = ["2010-08b1", "2011-08", "2012-08", "2013-09"]
    R["L4L_rate_clean"] = compare(T, "S", S15f, ["2015-12", "2018-06"], exact=True)
    R["L4L_once_clean"] = compare(T, "n_once", S15f, ["2015-12", "2018-06"], exact=True, minev=5)
    M = E[E.major].copy()
    a = M[M.seg == "S15f"]; b = M[(M.seg == "FWD") & (M.start < "2019-01-01") & (M.frac_S15 >= 0.9)]
    c = M[(M.seg == "FWD")]
    bc = b[~b.episode.isin(["2014-11", "2017-03"])]
    R["majors_per_pnsn_expday"] = dict(S15f=a.set_index("episode").S_per_pnsn_expday.round(2).to_dict(),
                                       FWD_S15=b.set_index("episode").S_per_pnsn_expday.round(2).to_dict(),
                                       median_S15f=float(a.S_per_pnsn_expday.median()), median_FWD_S15=float(b.S_per_pnsn_expday.median()),
                                       ratio_medians=float(b.S_per_pnsn_expday.median() / a.S_per_pnsn_expday.median()),
                                       mw_p=float(stats.mannwhitneyu(a.S_per_pnsn_expday, b.S_per_pnsn_expday).pvalue),
                                       ratio_medians_clean=float(bc.S_per_pnsn_expday.median() / a.S_per_pnsn_expday.median()),
                                       mw_p_clean=float(stats.mannwhitneyu(a.S_per_pnsn_expday, bc.S_per_pnsn_expday).pvalue),
                                       ratio_medians_allFWD=float(c.S_per_pnsn_expday.median() / a.S_per_pnsn_expday.median()),
                                       mw_p_allFWD=float(stats.mannwhitneyu(a.S_per_pnsn_expday, c.S_per_pnsn_expday).pvalue),
                                       once_ratio_medians=float((b.n_once / b.expdays_pnsn).median() / (a.n_once / a.expdays_pnsn).median()),
                                       once_mw_p=float(stats.mannwhitneyu(a.n_once / a.expdays_pnsn, b.n_once / b.expdays_pnsn).pvalue))
    for k, v in R.items():
        log(f"POSTHOC {k}: {json.dumps(v, default=float)}")
    json.dump(R, open(f"{OUT}/f6_rates_posthoc.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    import sys
    posthoc() if len(sys.argv) > 1 and sys.argv[1] == "posthoc" else main()
