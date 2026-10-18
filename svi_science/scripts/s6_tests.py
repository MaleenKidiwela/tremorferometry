"""S6 step 3: the five pre-registered tests + winner's curse + figures. Rules: results/s6/PLAN_S6.md.
Input results/s6/sizes.parquet (s6_measure.py). Output results/s6/tests.json, tests_table.csv, per-test csvs, figures."""
import json
import numpy as np, pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
R = f"{ROOT}/svi_science/results"; OUT = f"{R}/s6"; FIG = f"{ROOT}/svi_science/figures"
NB, NP = 2000, 2000
rng = np.random.default_rng(6)

S = pd.read_parquet(f"{OUT}/sizes.parquet")
S = S[S.kind == "event"].copy()
S["year"] = S.time_utc.dt.year + S.time_utc.dt.dayofyear / 365.25
def mkvar(S):
  return {
    "headline_complete": (S.P_same >= 0.75) & S.keep_xfam8 & S.complete.fillna(False),
    "all_rows_075": (S.P_same >= 0.75) & S.keep_xfam8,
    "P050_complete": (S.P_same >= 0.5) & S.keep_xfam8 & S.complete.fillna(False),
    "P085_complete": (S.P_same >= 0.85) & S.keep_xfam8 & S.complete.fillna(False),
    "raw_075_complete": (S.P_same >= 0.75) & S.complete.fillna(False),
  }


def boot_fam(cell_df, fam_col, val_col, nb=NB):
    """bootstrap mean of val over families (resample families, keep all their cells)"""
    g = {f: x[val_col].values for f, x in cell_df.groupby(fam_col)}
    fams = list(g)
    if not fams: return (np.nan, np.nan)
    bs = []
    for _ in range(nb):
        pick = rng.choice(len(fams), len(fams))
        v = np.concatenate([g[fams[i]] for i in pick]); bs.append(np.nanmean(v))
    return tuple(np.nanpercentile(bs, [2.5, 97.5]))


def verdict(eff, p, ci, eff_all, eff_N):
    seen = (p < 0.01) and (ci[0] > 0 or ci[1] < 0) and (np.sign(eff_all) == np.sign(eff)) and (abs(eff_N) < abs(eff))
    return "seen" if seen else "not seen"


# ---------------------------------------------------------------------------------------------- test 1: front vs behind
K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv"); K = K[K.region == "SVI"]
TR = pd.read_parquet(f"{R}/s3/tremor_svi.parquet")
FR = pd.read_csv(f"{R}/s3/families_svi_frame.csv").rename(columns={"fid": "fam"})
S = S.merge(FR[["fam", "along_strike_km"]].rename(columns={"along_strike_km": "as_s3"}), on="fam", how="left")


def arrival_map(tr):
    tr = tr.assign(b=np.floor(tr.as_km / 5) * 5)
    g = tr.groupby("b").time.apply(lambda s: s.sort_values().iloc[int(0.1 * (len(s) - 1))] if len(s) >= 20 else pd.NaT)
    return g.dropna()


maps = {e: arrival_map(TR[TR.ep == e]) for e in K.episode}
S["lag_d"] = np.nan
for e, amap in maps.items():
    ix = S.episode.astype(str) == e
    b = np.floor(S.loc[ix, "as_s3"] / 5) * 5
    S.loc[ix, "lag_d"] = (S.loc[ix, "time_utc"] - b.map(amap)).dt.total_seconds() / 86400
S["front_cls"] = np.where(S.lag_d.between(-1, 2, inclusive="left"), "front", np.where(S.lag_d >= 2, "behind", None))


def t1_effect(D, col="dM"):
    c = D[D.front_cls.notna()]
    rows = []
    for (f, e), g in c.groupby(["fam", "episode"], observed=True):
        a = g[g.front_cls == "behind"][col].dropna(); b = g[g.front_cls == "front"][col].dropna()
        if len(a) >= 3 and len(b) >= 3: rows.append((f, e, a.median() - b.median()))
    return pd.DataFrame(rows, columns=["fam", "episode", "d"])


def t1(D, col="dM", perm=True):
    C = t1_effect(D, col)
    if not len(C): return dict(n_cells=0)
    eff = C.d.mean(); ci = boot_fam(C, "fam", "d")
    out = dict(n_cells=len(C), n_fam=C.fam.nunique(), effect=eff, ci=ci)
    if perm:
        keep = D[D.front_cls.notna()].set_index(["fam", "episode"]).index.isin(C.set_index(["fam", "episode"]).index)
        X = D[D.front_cls.notna()][keep][["fam", "episode", "front_cls", col]].dropna().copy()
        grp = X.groupby(["fam", "episode"], observed=True).indices
        lab = X.front_cls.values.copy(); v = X[col].values
        null = []
        for _ in range(NP):
            L = lab.copy()
            for ix in grp.values(): L[ix] = rng.permutation(L[ix])
            ds = [np.median(v[ix][L[ix] == "behind"]) - np.median(v[ix][L[ix] == "front"]) for ix in grp.values()]
            null.append(np.mean(ds))
        null = np.array(null); out["p"] = float((np.sum(np.abs(null) >= abs(eff)) + 1) / (NP + 1))
    # column: Spearman dM vs lag within cell ranks
    c = D[D.lag_d.between(-1, 30)].dropna(subset=[col]).copy()
    c["rl"] = c.groupby(["fam", "episode"], observed=True).lag_d.rank(pct=True); c["rv"] = c.groupby(["fam", "episode"], observed=True)[col].rank(pct=True)
    out["spearman_within_cell"] = float(stats.spearmanr(c.rl, c.rv, nan_policy="omit")[0]); out["n_lag_rows"] = len(c)
    return out


# ---------------------------------------------------------------------------------------------- test 2: big vs small episodes
EP = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv")
MJ = EP[EP.major].copy(); med = MJ.family_exposure_days.median(); MJ["big"] = MJ.family_exposure_days > med
bigmap = dict(zip(MJ.episode, MJ.big)); sizemap = dict(zip(MJ.episode, MJ.family_exposure_days))
S["ep_big"] = S.episode.astype(str).map(bigmap)


def t2_eff(D, bm, col="dM"):
    c = D[D.major == True].dropna(subset=[col]).copy(); c["b"] = c.episode.astype(str).map(bm)
    rows = []
    for f, g in c.groupby("fam"):
        a = g[g.b == True][col]; b = g[g.b == False][col]
        if len(a) and len(b): rows.append((f, a.median() - b.median(), len(a), len(b)))
    return pd.DataFrame(rows, columns=["fam", "d", "na", "nb"])


def t2(D, col="dM", perm=True):
    C = t2_eff(D, bigmap, col); eff = C.d.mean(); ci = boot_fam(C, "fam", "d")
    out = dict(n_fam=len(C), effect=eff, ci=ci, effect_ge3_each=float(C[(C.na >= 3) & (C.nb >= 3)].d.mean()),
               n_fam_ge3_each=int(((C.na >= 3) & (C.nb >= 3)).sum()))
    if perm:
        eps = list(bigmap); vals = np.array([bigmap[e] for e in eps]); null = []
        for _ in range(NP):
            bm = dict(zip(eps, rng.permutation(vals))); null.append(t2_eff(D, bm, col).d.mean())
        null = np.array(null); out["p"] = float((np.sum(np.abs(null) >= abs(eff)) + 1) / (NP + 1))
    E = D[D.major == True].groupby(D.episode.astype(str))[col].agg(["median", "count"]).reset_index()
    E["size"] = E.episode.map(sizemap); E = E[E["count"] >= 10]
    out["episode_spearman_vs_logsize"] = float(stats.spearmanr(E["median"], np.log10(E["size"]))[0]) if len(E) > 3 else np.nan
    out["episode_spearman_p"] = float(stats.spearmanr(E["median"], np.log10(E["size"]))[1]) if len(E) > 3 else np.nan
    out["n_episodes"] = len(E)
    G = pd.read_csv(f"{R}/s1/episode_geodetic.csv")
    G = G[(G.SN > 3) & (G.south_only != True)][["episode", "D"]]
    EG = E.merge(G, on="episode")
    out["geodetic_spearman"] = float(stats.spearmanr(EG["median"], EG.D)[0]) if len(EG) > 3 else np.nan
    out["geodetic_n"] = len(EG)
    return out, E


# ---------------------------------------------------------------------------------------------- test 3: between vs ETS
BE = pd.read_csv(f"{R}/s7/between_events.csv", parse_dates=["time_utc"])
S["between"] = S.set_index(["fam", "time_utc"]).index.isin(BE.set_index(["fam", "time_utc"]).index)


def t3(D, col="dM", flag="between", perm=True):
    c = D.dropna(subset=[col]); rows = []; pools = {}
    for f, g in c.groupby("fam"):
        a = g[g[flag] == True][col]; b = g[g.major == True][col]
        if len(a) and len(b): rows.append((f, a.median() - b.median(), len(a))); pools[f] = (b.values, len(a))
    C = pd.DataFrame(rows, columns=["fam", "d", "n_between"])
    if not len(C): return dict(n_fam=0)
    eff = C.d.mean(); ci = boot_fam(C, "fam", "d")
    out = dict(n_fam=len(C), n_events=int(C.n_between.sum()), effect=eff, ci=ci)
    if perm:
        null = []
        for _ in range(NP):
            ds = []
            for f, (b, k) in pools.items():
                pick = rng.choice(len(b), k, replace=False) if k <= len(b) else rng.choice(len(b), k)
                ds.append(np.median(b[pick]) - np.median(b))
            null.append(np.mean(ds))
        null = np.array(null); out["p"] = float((np.sum(np.abs(null) >= abs(eff)) + 1) / (NP + 1))
    return out


# ---------------------------------------------------------------------------------------------- test 4: 12-year trend
def ts_slope(t, y):
    i, j = np.triu_indices(len(t), 1); dt = t[j] - t[i]; ok = dt != 0
    return np.median((y[j] - y[i])[ok] / dt[ok]) * 10 if ok.any() else np.nan


def t4(D, col="dM", perm=True, nperm=1000):
    c = D.dropna(subset=[col]); fams = []
    for f, g in c.groupby("fam"):
        if len(g) >= 20 and g.year.max() - g.year.min() >= 5: fams.append(f)
    rows = []; blocks = {}
    for f in fams:
        g = c[c.fam == f]
        rows.append((f, ts_slope(g.year.values, g[col].values), len(g)))
        blk = np.where(g.episode.notna(), g.episode.astype(str), g.time_utc.dt.strftime("m%Y-%m"))
        blocks[f] = (g.year.values, g[col].values, blk)
    C = pd.DataFrame(rows, columns=["fam", "slope", "n"])
    eff = C.slope.mean(); ci = boot_fam(C, "fam", "slope")
    out = dict(n_fam=len(C), effect_per_decade=eff, ci=ci, median_slope=float(C.slope.median()))
    if perm:
        null = []
        for _ in range(nperm):
            ss = []
            for f, (t, y, blk) in blocks.items():
                ub = np.unique(blk); tb = {b: np.median(t[blk == b]) for b in ub}
                newc = dict(zip(ub, rng.permutation([tb[b] for b in ub])))
                t2 = np.array([t[k] - tb[blk[k]] + newc[blk[k]] for k in range(len(t))])
                ss.append(ts_slope(t2, y))
            null.append(np.nanmean(ss))
        null = np.array(null); out["p"] = float((np.sum(np.abs(null) >= abs(eff)) + 1) / (nperm + 1)); out["nperm"] = nperm
    return out, C


# ---------------------------------------------------------------------------------------------- test 5: depth
FSZ = pd.read_csv(f"{OUT}/family_sizes.csv")


def t5():
    F = FSZ.copy(); F["sbin"] = np.floor(F.along_strike_km / 10)
    out = {}
    for lab, col in [("bostock_mag_med", "mag_med"), ("our_M_f", "M_f")]:
        x = F[["dep", col, "sbin"]].dropna()
        r = stats.spearmanr(x.dep, x[col])[0]
        null = []
        for _ in range(NP):
            d = x.groupby("sbin").dep.transform(lambda s: rng.permutation(s.values)); null.append(stats.spearmanr(d, x[col])[0])
        p = float((np.sum(np.abs(null) >= abs(r)) + 1) / (NP + 1))
        bs = []
        for _ in range(NB):
            k = rng.choice(len(x), len(x)); bs.append(stats.spearmanr(x.dep.values[k], x[col].values[k])[0])
        out[lab] = dict(n_fam=len(x), spearman=r, p=p, ci=tuple(np.nanpercentile(bs, [2.5, 97.5])))
    VAR = mkvar(S)
    H = S[VAR["headline_complete"]].dropna(subset=["dM"])
    iqr = H.groupby("fam").dM.agg(lambda v: np.subtract(*np.percentile(v, [75, 25])) if len(v) >= 20 else np.nan).dropna()
    x = F.set_index("fam").loc[iqr.index]
    out["iqr_vs_depth"] = dict(n_fam=len(iqr), spearman=float(stats.spearmanr(x.dep, iqr)[0]), p_param=float(stats.spearmanr(x.dep, iqr)[1]))
    out["M_f_vs_mag_med_spearman"] = float(stats.spearmanr(F.M_f, F.mag_med, nan_policy="omit")[0])
    out["spreading_n"] = float(F.spreading_n.iloc[0])
    return out


# ---------------------------------------------------------------------------------------------- run
def main():
    global VAR
    VAR = mkvar(S)
    res = {}
    H = S[VAR["headline_complete"]]
    A = S[VAR["all_rows_075"]]
    # winner's curse
    def wc(D):
        c = D.dropna(subset=["dM", "net_mad"]).copy()
        c["a"] = c.groupby("fam").dM.rank(pct=True); c["b"] = c.groupby("fam").net_mad.rank(pct=True)
        return float(stats.spearmanr(c.a, c.b)[0]), len(c)
    NE = pd.read_parquet(f"{OUT}/sizes.parquet"); NE = NE[NE.kind == "nonevent"]
    res["winners_curse"] = dict(all_rows_075=wc(A), headline_complete=wc(H), nonevent_reference=wc(NE),
                                frac_complete_075=float(A.complete.mean()), frac_low_sn_075=float(A.low_sn.mean()))
    # size summary
    res["sizes"] = dict(n_head_all=len(A), n_head_dM=int(A.dM.notna().sum()), n_complete=len(H),
                        dM_q=list(np.nanpercentile(A.dM, [5, 25, 50, 75, 95])),
                        dM_complete_q=list(np.nanpercentile(H.dM, [5, 25, 50, 75, 95])),
                        nonevent_dM_q=list(np.nanpercentile(NE.dM, [5, 25, 50, 75, 95])),
                        M_approx_q=list(np.nanpercentile(A.M_approx, [5, 25, 50, 75, 95])),
                        M_approx_complete_q=list(np.nanpercentile(H.M_approx, [5, 25, 50, 75, 95])),
                        family_iqr_median=float(H.groupby("fam").dM.agg(lambda v: np.subtract(*np.percentile(v, [75, 25])) if len(v) >= 20 else np.nan).median()),
                        family_ratio_p95_p5_median=float(H.groupby("fam").dM.agg(lambda v: np.subtract(*np.percentile(v, [95, 5])) if len(v) >= 20 else np.nan).median()))
    tab = []
    # test 1
    t = {v: t1(S[m], perm=(v == "headline_complete")) for v, m in VAR.items()}
    tN = t1(S[VAR["headline_complete"]], col="dN", perm=False)
    t["noise_dN"] = tN; res["t1_front_vs_behind"] = t
    h = t["headline_complete"]
    tab.append(dict(test="1 behind minus at-front", effect=h["effect"], ci_lo=h["ci"][0], ci_hi=h["ci"][1], p=h["p"], n=h["n_cells"],
                    effect_all_rows=t["all_rows_075"]["effect"], effect_noise=tN["effect"],
                    verdict=verdict(h["effect"], h["p"], h["ci"], t["all_rows_075"]["effect"], tN["effect"])))
    # front-lag curve rows for the figure
    c = S[VAR["all_rows_075"] & S.lag_d.between(-3, 25)].dropna(subset=["dM"])
    c[["fam", "episode", "lag_d", "dM", "dN", "complete"]].to_csv(f"{OUT}/t1_lag_rows.csv", index=False)
    # test 2
    t = {}; Es = None
    for v, m in VAR.items():
        o, E = t2(S[m], perm=(v == "headline_complete")); t[v] = o
        if v == "headline_complete": Es = E
    oN, _ = t2(S[VAR["headline_complete"]], col="dN", perm=False); t["noise_dN"] = oN
    res["t2_big_vs_small"] = t; Es.to_csv(f"{OUT}/t2_episode_medians.csv", index=False)
    h = t["headline_complete"]
    tab.append(dict(test="2 big minus small episodes", effect=h["effect"], ci_lo=h["ci"][0], ci_hi=h["ci"][1], p=h["p"], n=h["n_fam"],
                    effect_all_rows=t["all_rows_075"]["effect"], effect_noise=oN["effect"],
                    verdict=verdict(h["effect"], h["p"], h["ci"], t["all_rows_075"]["effect"], oN["effect"])))
    # test 3
    t = {v: t3(S[m], perm=(v == "headline_complete")) for v, m in VAR.items()}
    t["outside_all_episodes_column"] = t3(S[VAR["headline_complete"]].assign(outside=lambda d: d.in_episode == False), flag="outside", perm=True)
    t["outside_all_episodes_all_rows"] = t3(S[VAR["all_rows_075"]].assign(outside=lambda d: d.in_episode == False), flag="outside", perm=False)
    tN = t3(S[VAR["headline_complete"]], col="dN", perm=False); t["noise_dN"] = tN
    res["t3_between_vs_ETS"] = t
    h = t["headline_complete"]
    tab.append(dict(test="3 between minus ETS", effect=h.get("effect"), ci_lo=h["ci"][0], ci_hi=h["ci"][1], p=h.get("p"), n=h.get("n_events"),
                    effect_all_rows=t["all_rows_075"].get("effect"), effect_noise=tN.get("effect"),
                    verdict=verdict(h["effect"], h["p"], h["ci"], t["all_rows_075"]["effect"], tN["effect"])))
    # test 4
    t = {}
    for v, m in VAR.items():
        o, C = t4(S[m], perm=(v == "headline_complete")); t[v] = o
        if v == "headline_complete": C.to_csv(f"{OUT}/t4_family_slopes.csv", index=False)
    oN, CN = t4(S[VAR["headline_complete"]], col="dN", perm=False); t["noise_dN"] = oN
    oS, _ = t4(S[VAR["headline_complete"]], col="dM_senscorr_only", perm=False); t["sensitivity_corrected_only"] = oS
    res["t4_years"] = t
    h = t["headline_complete"]
    tab.append(dict(test="4 trend per decade", effect=h["effect_per_decade"], ci_lo=h["ci"][0], ci_hi=h["ci"][1], p=h["p"], n=h["n_fam"],
                    effect_all_rows=t["all_rows_075"]["effect_per_decade"], effect_noise=oN["effect_per_decade"],
                    verdict=verdict(h["effect_per_decade"], h["p"], h["ci"], t["all_rows_075"]["effect_per_decade"], oN["effect_per_decade"])))
    # test 5
    o = t5(); res["t5_depth"] = o
    for lab in ["bostock_mag_med", "our_M_f"]:
        h = o[lab]
        tab.append(dict(test=f"5 depth Spearman ({lab})", effect=h["spearman"], ci_lo=h["ci"][0], ci_hi=h["ci"][1], p=h["p"], n=h["n_fam"],
                        effect_all_rows=np.nan, effect_noise=np.nan,
                        verdict="seen" if (h["p"] < 0.01 and (h["ci"][0] > 0 or h["ci"][1] < 0)) else "not seen"))
    T = pd.DataFrame(tab); T.to_csv(f"{OUT}/tests_table.csv", index=False)
    json.dump(res, open(f"{OUT}/tests.json", "w"), indent=1, default=lambda x: float(x) if np.isscalar(x) else str(x))
    S.drop(columns=[c for c in S.columns if c.startswith("_")]).to_parquet(f"{OUT}/sizes_annotated.parquet")
    print(T.round(4).to_string()); print(json.dumps(res["winners_curse"], indent=1)); print(json.dumps(res["sizes"], indent=1, default=float))


if __name__ == "__main__":
    main()
