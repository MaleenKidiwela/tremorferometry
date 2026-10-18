#!/usr/bin/env python
"""s10_multi_analysis -- S10 multi steps 2-5 (multi/PLAN.md §2-5): family-episode CC^2-weighted means, cross-family
agreement (pairwise r, family vs rest, 10,000-draw within-family episode-label permutation), PGC-only / borehole-only
subsets, scrambled control through the same path, lapse check, pooled episode series with family bootstrap and
leave-one-family-out, within-episode change vs days since the front. Values in % dv/v (= -eps * 100)."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M

NPERM, NBOOT, MIN_SHARED, MIN_VALS = int(os.environ.get("S10M_NPERM", 10000)), int(os.environ.get("S10M_NBOOT", 1000)), 8, 3
log = S.Log(f"{M.LOGD}/analysis.log")
RES = {}


def load():
    D = pd.read_parquet(f"{M.TB}/daily.parquet")
    U = pd.read_csv(f"{M.TB}/episode_units.csv")
    D = D.merge(U[["episode", "cls", "start", "big"]], on="episode", how="left")
    D["v24"] = -D.eps24 * 100; D["w24"] = D.cc24.clip(lower=0) ** 2; D["s24"] = D.sig24 * 100
    D["v35"] = -D.eps35 * 100; D["w35"] = D.cc35.clip(lower=0) ** 2; D["s35"] = D.sig35 * 100
    D["vc"] = -D.c2_eps24 * 100; D["wc"] = D.c2_cc24.clip(lower=0) ** 2; D["sc"] = D.c2_sig24 * 100
    D["vc35"] = -D.c2_eps35 * 100; D["wc35"] = D.c2_cc35.clip(lower=0) ** 2; D["sc35"] = D.c2_sig35 * 100
    D["one24"] = np.where(D.v24.notna(), 1.0, np.nan)
    return D, U


GROUPS = {"core5": M.CORE5, "pgc": ["CN.PGC"], "bh": M.BOREHOLES, "all7": M.STATIONS}
VARS = {"d24": ("v24", "w24", "s24", "cc24"), "d35": ("v35", "w35", "s35", "cc35"), "c24": ("vc", "wc", "sc", "c2_cc24"),
        "c35": ("vc35", "wc35", "sc35", "c2_cc35"), "u24": ("v24", "one24", "s24", "cc24")}


def fe_table(D, var, grp, cls, min_vals=MIN_VALS):
    v, w, s, c = VARS[var]
    X = D[D.netsta.isin(GROUPS[grp]) & D.cls.isin(cls) & D[v].notna() & (D[w] > 0)]
    rows = []
    for (f, e), g in X.groupby(["fam", "episode"]):
        if len(g) < min_vals:
            continue
        ww = g[w].to_numpy(); x = g[v].to_numpy(); sg = g[s].fillna(g[s].median()).to_numpy()
        m = np.sum(ww * x) / ww.sum()
        neff = ww.sum() ** 2 / np.sum(ww ** 2)
        sb = np.sqrt(np.sum(ww ** 2 * sg ** 2)) / ww.sum()
        ss = np.sqrt(np.sum(ww * (x - m) ** 2) / ww.sum()) / np.sqrt(max(neff - 1, 1))
        rows.append(dict(fam=f, episode=e, value=m, sig_boot=sb, sig_scatter=ss, sigma=max(sb, ss), n=len(g),
                         n_days=g.date.nunique(), n_sta=g.netsta.nunique(), wsum=ww.sum(), mean_cc=float(np.mean(g[c]))))
    return pd.DataFrame(rows)


def matrix(F):
    V = F.pivot(index="fam", columns="episode", values="value")
    Wt = F.pivot(index="fam", columns="episode", values="wsum").reindex_like(V)
    return V, Wt


def pair_r(Vv, Mk):
    """Vv (nf, ne) values with 0 where missing, Mk mask -> n_shared, r matrices (Pearson over shared episodes)"""
    n = Mk @ Mk.T
    Sx = Vv @ Mk.T; Sy = Sx.T
    Sxx = (Vv ** 2) @ Mk.T; Syy = Sxx.T
    Sxy = Vv @ Vv.T
    num = n * Sxy - Sx * Sy
    den = np.sqrt(np.clip(n * Sxx - Sx ** 2, 0, None) * np.clip(n * Syy - Sy ** 2, 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    return n, r


def agreement(F, tag, perm=True, seed=M.SEED):
    if len(F) == 0:
        return {}
    V, Wt = matrix(F)
    Mk = V.notna().to_numpy().astype(float); Vv = np.nan_to_num(V.to_numpy())
    n, r = pair_r(Vv, Mk)
    iu = np.triu_indices(len(V), 1)
    elig = n[iu] >= MIN_SHARED
    rr = r[iu][elig]
    out = dict(set=tag, n_fam=len(V), n_ep=V.shape[1], n_pairs=int(elig.sum()), M=float(np.nanmedian(rr)) if elig.any() else np.nan,
               mean_r=float(np.nanmean(rr)) if elig.any() else np.nan, frac_pos=float(np.mean(rr > 0)) if elig.any() else np.nan)
    # Spearman per pair (column)
    if elig.any():
        sp = []
        A = V.to_numpy()
        for i, j in zip(iu[0][elig], iu[1][elig]):
            m = ~np.isnan(A[i]) & ~np.isnan(A[j]); sp.append(stats.spearmanr(A[i, m], A[j, m])[0])
        out["M_spearman"] = float(np.nanmedian(sp))
    # family vs rest
    A = V.to_numpy(); Wa = np.nan_to_num(Wt.to_numpy())
    fr = []
    for i in range(len(A)):
        o = np.delete(np.arange(len(A)), i)
        num = np.nansum(np.nan_to_num(A[o]) * Wa[o] * (~np.isnan(A[o])), 0); den = np.sum(Wa[o] * (~np.isnan(A[o])), 0)
        rest = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        m = ~np.isnan(A[i]) & ~np.isnan(rest)
        fr.append(dict(fam=V.index[i], n=int(m.sum()), r=stats.pearsonr(A[i, m], rest[m])[0] if m.sum() >= MIN_SHARED else np.nan))
    FR = pd.DataFrame(fr)
    out["fam_vs_rest_median_r"] = float(FR.r.median()); out["fam_vs_rest_n"] = int(FR.r.notna().sum())
    out["fam_vs_rest_frac_pos"] = float((FR.r.dropna() > 0).mean())
    if perm and elig.any():
        rng = np.random.default_rng(seed)
        idx = [np.flatnonzero(Mk[i] > 0) for i in range(len(V))]
        Ms = np.empty(NPERM)
        for b in range(NPERM):
            Vp = Vv.copy()
            for i, ii in enumerate(idx):
                Vp[i, ii] = Vv[i, rng.permutation(ii)]
            _, rp = pair_r(Vp, Mk)
            Ms[b] = np.nanmedian(rp[iu][elig])
        out.update(perm_p99=float(np.quantile(Ms, 0.99)), perm_p95=float(np.quantile(Ms, 0.95)), perm_median=float(np.median(Ms)),
                   perm_p=float((1 + np.sum(Ms >= out["M"])) / (NPERM + 1)))
        np.save(f"{M.TB}/perm_{tag}.npy", Ms)
    pairs = pd.DataFrame(dict(fam_i=V.index[iu[0]], fam_j=V.index[iu[1]], n_shared=n[iu], r=r[iu], eligible=elig))
    return out, pairs, FR


def pooled(D, var, grp, cls, fams=None):
    v, w, s, c = VARS[var]
    X = D[D.netsta.isin(GROUPS[grp]) & D.cls.isin(cls) & D[v].notna() & (D[w] > 0)]
    if fams is not None:
        X = X[X.fam.isin(fams)]
    g = X.assign(wx=X[w] * X[v]).groupby("episode")
    P = pd.DataFrame(dict(value=g.wx.sum() / g[w].sum(), n_fam=g.fam.nunique(), n_vals=g.size(), mean_cc=g[c].mean(), wsum=g[w].sum()))
    return P, X


def pooled_boot(X, var, seed):
    """family bootstrap SD of the pooled per-episode value"""
    v, w, s, c = VARS[var]
    fams = np.array(sorted(X.fam.unique()))
    agg = X.assign(wx=X[w] * X[v]).groupby(["fam", "episode"])[["wx", w]].sum()
    WX = agg.wx.unstack().reindex(fams).fillna(0).to_numpy(); WW = agg[w].unstack().reindex(fams).fillna(0).to_numpy()
    eps_ = agg.wx.unstack().columns
    rng = np.random.default_rng(seed); B = np.empty((NBOOT, WX.shape[1]))
    for b in range(NBOOT):
        k = np.bincount(rng.integers(0, len(fams), len(fams)), minlength=len(fams)).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            B[b] = (k @ WX) / (k @ WW)
    return pd.Series(np.nanstd(B, 0, ddof=1), index=eps_)


def within_episode(D, var, grp, FR):
    v, w, s, c = VARS[var]
    X = D[D.netsta.isin(GROUPS[grp]) & (D.cls == "major") & D[v].notna() & (D[w] > 0)].copy()
    X["date_t"] = pd.to_datetime(X.date)
    fr = FR.set_index(["fam", "episode"]).front
    X["front"] = [fr.get((f, e), pd.NaT) for f, e in zip(X.fam, X.episode)]
    X["x_src"] = np.where(X.front.notna(), "S5 front", "episode start")
    X["front"] = X.front.fillna(pd.to_datetime(X.start))
    X["x"] = (X.date_t - X.front).dt.days
    gk = ["fam", "netsta", "epoch", "episode"]
    X["gw"] = X.groupby(gk)[w].transform("sum"); X["gwx"] = (X[w] * X[v]).groupby([X[k] for k in gk]).transform("sum")
    X["gn"] = X.groupby(gk)[v].transform("size")
    X = X[X.gn >= 2]
    X["y"] = X[v] - X.gwx / X.gw
    X = X[(X.x >= -2) & (X.x <= 30)]
    return X


def we_stats(X, w, seed):
    def st(Y):
        b = Y.assign(wy=Y[w] * Y.y).groupby("x")
        prof = b.wy.sum() / b[w].sum()
        m = (Y.x >= 0) & (Y.x <= 20)
        Z = Y[m]
        if len(Z) < 5 or Z.x.nunique() < 3:
            return prof, np.nan, np.nan
        sl = np.polyfit(Z.x, Z.y, 1, w=np.sqrt(Z[w]))[0]
        e = Y[(Y.x >= 0) & (Y.x <= 3)]; l = Y[Y.x >= 4]
        con = np.sum(e[w] * e.y) / e[w].sum() - np.sum(l[w] * l.y) / l[w].sum()
        return prof, sl, con
    if len(X) < 5:
        return pd.Series(dtype=float), dict(slope=np.nan, slope_sig=np.nan, contrast=np.nan, contrast_sig=np.nan, n=len(X), n_fam=0)
    prof, sl, con = st(X)
    fams = np.array(sorted(X.fam.unique())); rng = np.random.default_rng(seed)
    G = {f: g for f, g in X.groupby("fam")}
    bs, bc = [], []
    for _ in range(NBOOT):
        pick = rng.choice(fams, len(fams))
        Y = pd.concat([G[f] for f in pick], ignore_index=True)
        _, a, b = st(Y); bs.append(a); bc.append(b)
    return prof, dict(slope=sl, slope_sig=float(np.nanstd(bs)), contrast=con, contrast_sig=float(np.nanstd(bc)), n=len(X),
                      n_fam=len(fams))


def es_frame(U):
    """episode stacks shaped like the daily table (one 'day' = one episode stack)"""
    E = pd.read_parquet(f"{M.TB}/episode_stacks.parquet").merge(U[["episode", "cls", "start", "big"]], on="episode")
    E["date"] = E.episode
    E["v24"] = -E.eps24 * 100; E["w24"] = E.cc24.clip(lower=0) ** 2; E["s24"] = E.sig24 * 100
    E["v35"] = -E.eps35 * 100; E["w35"] = E.cc35.clip(lower=0) ** 2; E["s35"] = E.sig35 * 100
    for k in ("24", "35"):
        if f"c2_eps{k}" not in E:
            E[f"c2_eps{k}"] = np.nan; E[f"c2_cc{k}"] = np.nan; E[f"c2_sig{k}"] = np.nan
    E["vc"] = -E.c2_eps24 * 100; E["wc"] = E.c2_cc24.clip(lower=0) ** 2; E["sc"] = E.c2_sig24 * 100
    E["vc35"] = -E.c2_eps35 * 100; E["wc35"] = E.c2_cc35.clip(lower=0) ** 2; E["sc35"] = E.c2_sig35 * 100
    E["one24"] = np.where(E.v24.notna(), 1.0, np.nan)
    return E


def wls_line(t, x, sg):
    w = 1 / sg ** 2; A = np.vstack([t - t.mean(), np.ones(len(t))]).T
    b, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], x * np.sqrt(w), rcond=None)
    cov = np.linalg.pinv(A.T @ (A * w[:, None])); chi = np.sum(w * (x - A @ b) ** 2) / max(len(t) - 2, 1)
    return float(b[0]), float(np.sqrt(cov[0, 0] * max(chi, 1)))


def series_battery(tag, D, U, F, F35):
    MAJ, ALL = ["major"], ["major", "burst"]
    PS = {}
    for k, (var, grp) in {"core5_d24": ("d24", "core5"), "pgc_d24": ("d24", "pgc"), "bh_d24": ("d24", "bh"),
                          "core5_d35": ("d35", "core5"), "core5_c24": ("c24", "core5"), "all7_d24": ("d24", "all7")}.items():
        P, X = pooled(D, var, grp, ALL)
        if len(P):
            P["sigma_famboot"] = pooled_boot(X, var, M.SEED).reindex(P.index)
        PS[k] = P
    T = U.set_index("episode")[["start", "end", "cls", "big"]]
    for k, P in PS.items():
        T = T.join(P.add_prefix(f"{k}_"), how="left")
    T = T[T.core5_d24_value.notna()]
    base = PS["core5_d24"].value
    fams = sorted(D.fam.unique()); lo = {}
    for f in fams:
        P, _ = pooled(D[D.fam != f], "d24", "core5", ALL); lo[f] = P.value
    LO = pd.DataFrame(lo).reindex(base.index)
    T["lofo_max_abs_d"] = LO.sub(base, axis=0).abs().max(axis=1)
    T["lofo_min"] = LO.min(axis=1); T["lofo_max"] = LO.max(axis=1)
    T = T.sort_values("start"); T.to_csv(f"{M.TB}/pooled_series_{tag}.csv")
    maj = T[T.cls == "major"]
    R = {}
    lofo_r = [stats.pearsonr(LO.loc[maj.index, f].fillna(base[maj.index]), base[maj.index])[0] for f in fams]
    R["lofo_min_r"] = float(np.min(lofo_r)); R["lofo_median_max_abs_d"] = float(maj.lofo_max_abs_d.median())
    R["lofo_max_max_abs_d"] = float(maj.lofo_max_abs_d.max())
    for k in ("core5_d24", "pgc_d24", "bh_d24", "core5_c24", "core5_d35"):
        if f"{k}_sigma_famboot" not in maj:
            continue
        x = maj[f"{k}_value"]; sg = maj[f"{k}_sigma_famboot"]; ok = x.notna() & sg.notna() & (maj[f"{k}_n_fam"] >= 3)
        R[f"{k}_n"] = int(ok.sum())
        if ok.sum() < 3:
            continue
        R[f"{k}_sd"] = float(x[ok].std()); R[f"{k}_rms_sigma"] = float(np.sqrt(np.mean(sg[ok] ** 2)))
        R[f"{k}_amp"] = float(np.sqrt(max(0.0, x[ok].var() - np.mean(sg[ok] ** 2))))
        wm = np.average(x[ok], weights=1 / sg[ok] ** 2)
        R[f"{k}_chi2red"] = float(np.sum(((x[ok] - wm) / sg[ok]) ** 2) / (ok.sum() - 1))
        R[f"{k}_chi2_p"] = float(stats.chi2.sf(np.sum(((x[ok] - wm) / sg[ok]) ** 2), ok.sum() - 1))
        R[f"{k}_n_beyond3sig"] = int((np.abs(x[ok] - wm) > 3 * sg[ok]).sum())
        tt = (pd.to_datetime(maj.start[ok]).dt.year + pd.to_datetime(maj.start[ok]).dt.dayofyear / 365.25).to_numpy()
        R[f"{k}_trend"], R[f"{k}_trend_sig"] = wls_line(tt, x[ok].to_numpy(), sg[ok].to_numpy())
        # seasonal (annual) sinusoid, post hoc
        ph = 2 * np.pi * (pd.to_datetime(maj.start[ok]).dt.dayofyear.to_numpy() + 10) / 365.25
        A = np.vstack([np.cos(ph), np.sin(ph), np.ones(len(ph))]).T; w = 1 / sg[ok].to_numpy() ** 2
        b, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], x[ok].to_numpy() * np.sqrt(w), rcond=None)
        chi1 = np.sum(w * (x[ok].to_numpy() - A @ b) ** 2); chi0 = np.sum(w * (x[ok].to_numpy() - wm) ** 2)
        R[f"{k}_annual_peak_doy"] = float((np.degrees(np.arctan2(b[1], b[0])) % 360) / 360 * 365.25 - 10) % 365.25
        if k == "core5_d24":
            RES[f"annual_model_{tag}"] = [float(b[0]), float(b[1])]
            xs = x[ok].to_numpy(); rng = np.random.default_rng(M.SEED); amps = np.empty(NPERM)
            for i_ in range(NPERM):
                pp = rng.permutation(len(xs))
                bb, *_ = np.linalg.lstsq(A * np.sqrt(w[pp])[:, None], xs[pp] * np.sqrt(w[pp]), rcond=None); amps[i_] = np.hypot(bb[0], bb[1])
            R[f"{k}_annual_perm_p"] = float((1 + np.sum(amps >= np.hypot(b[0], b[1]))) / (NPERM + 1))
        R[f"{k}_annual_amp"] = float(np.hypot(b[0], b[1])); R[f"{k}_annual_p"] = float(stats.f.sf(((chi0 - chi1) / 2) / (chi1 / max(ok.sum() - 3, 1)), 2, max(ok.sum() - 3, 1)))
    for a_, b_ in (("pgc_d24", "bh_d24"), ("core5_d24", "core5_c24")):
        if f"{b_}_value" in maj:
            ok = maj[f"{a_}_value"].notna() & maj[f"{b_}_value"].notna()
            if ok.sum() >= 4:
                rr, pp = stats.pearsonr(maj[f"{a_}_value"][ok], maj[f"{b_}_value"][ok])
                R[f"r_{a_}_vs_{b_}"] = float(rr); R[f"p_{a_}_vs_{b_}"] = float(pp); R[f"n_{a_}_vs_{b_}"] = int(ok.sum())
    # ICC / one-way ANOVA of family-episode values on episode
    g = F.groupby("episode").value; k_i = g.size(); a = len(k_i); N = len(F)
    msb = np.sum(k_i * (g.mean() - F.value.mean()) ** 2) / (a - 1); msw = np.sum((F.value - F.episode.map(g.mean())) ** 2) / (N - a)
    k0 = (N - np.sum(k_i ** 2) / N) / (a - 1)
    R["icc"] = float((msb - msw) / (msb + (k0 - 1) * msw)); R["anova_F"] = float(msb / msw); R["anova_p"] = float(stats.f.sf(msb / msw, a - 1, N - a))
    # lapse
    ok = maj.core5_d24_value.notna() & maj.core5_d35_value.notna() & (maj.core5_d24_n_fam >= 3)
    x, y = maj.core5_d24_value[ok], maj.core5_d35_value[ok]
    R["lapse_pooled_slope"] = float(stats.theilslopes(y, x)[0]); R["lapse_pooled_sign"] = float(np.mean(np.sign(x) == np.sign(y)))
    R["lapse_pooled_r"] = float(stats.pearsonr(x, y)[0]); R["lapse_pooled_n"] = int(ok.sum())
    R["lapse_pass"] = bool(R["lapse_pooled_slope"] >= 0.6 and R["lapse_pooled_sign"] >= 0.7)
    Fd = F.merge(F35, on=["fam", "episode"], suffixes=("", "_35"))
    R["lapse_fe_slope"] = float(stats.theilslopes(Fd.value_35, Fd.value)[0]); R["lapse_fe_sign"] = float(np.mean(np.sign(Fd.value) == np.sign(Fd.value_35)))
    R["lapse_fe_r"] = float(stats.pearsonr(Fd.value, Fd.value_35)[0])
    dd = D[D.netsta.isin(M.CORE5) & D.v24.notna() & D.v35.notna()]
    R["lapse_unit_slope"] = float(stats.theilslopes(dd.v35, dd.v24)[0]); R["lapse_unit_sign"] = float(np.mean(np.sign(dd.v24) == np.sign(dd.v35)))
    R["lapse_unit_r"] = float(stats.pearsonr(dd.v24, dd.v35)[0])
    RES[f"battery_{tag}"] = R
    log(f"battery {tag}: {json.dumps(R, default=float)}")


def main():
    D, U = load()
    log(f"daily values {len(D)}; non-railed 2-4 s {D.v24.notna().sum()}; families {D.fam.nunique()}")
    # rail / coverage summary by station
    cov = D.groupby("netsta").agg(n=("v24", "size"), ok24=("v24", lambda x: x.notna().mean()), ok35=("v35", lambda x: x.notna().mean()),
                                  okc=("vc", lambda x: x.notna().mean()), cc24=("cc24", "median"), sig24=("s24", "median"))
    cov.to_csv(f"{M.TB}/coverage_by_station.csv"); log("coverage:\n" + cov.to_string())
    MAJ, ALL, BUR = ["major"], ["major", "burst"], ["burst"]
    # family-episode tables
    sets = {"core5_d24": ("d24", "core5", MAJ), "pgc_d24": ("d24", "pgc", MAJ), "bh_d24": ("d24", "bh", MAJ),
            "core5_d35": ("d35", "core5", MAJ), "core5_c24": ("c24", "core5", MAJ), "core5_c35": ("c35", "core5", MAJ),
            "pgc_c24": ("c24", "pgc", MAJ), "bh_c24": ("c24", "bh", MAJ),
            "core5_d24_all": ("d24", "core5", ALL), "core5_d24_bursts": ("d24", "core5", BUR), "all7_d24": ("d24", "all7", MAJ),
            "core5_u24": ("u24", "core5", MAJ), "pgc_d35": ("d35", "pgc", MAJ), "bh_d35": ("d35", "bh", MAJ)}
    FE, AG, PAIRS = {}, [], {}
    for k, (var, grp, cls) in sets.items():
        F = fe_table(D, var, grp, cls); FE[k] = F
        F.to_csv(f"{M.TB}/fe_{k}.csv", index=False)
        res = agreement(F, k)
        if res:
            out, pairs, FRr = res; AG.append(out); PAIRS[k] = pairs
            pairs.to_csv(f"{M.TB}/pairs_{k}.csv", index=False); FRr.to_csv(f"{M.TB}/fam_vs_rest_{k}.csv", index=False)
            log(f"{k}: {json.dumps({a: (round(b, 4) if isinstance(b, float) else b) for a, b in out.items()})}")
    AGT = pd.DataFrame(AG); AGT.to_csv(f"{M.TB}/agreement.csv", index=False)
    ag = AGT.set_index("set")
    PASS = bool(ag.loc["core5_d24", "M"] > ag.loc["core5_d24", "perm_p99"] and ag.loc["pgc_d24", "M"] > 0 and ag.loc["bh_d24", "M"] > 0)
    CTRL_PASS = bool(ag.loc["core5_c24", "M"] > ag.loc["core5_c24", "perm_p99"]) if ("core5_c24" in ag.index and np.isfinite(ag.loc["core5_c24", "M"])) else None
    RES.update(headline_pass=PASS, control_pass=CTRL_PASS, agreement=AGT.to_dict("records"))
    # Revision 1: episode-stack variant through the same battery (post hoc)
    ES = es_frame(U)
    AG_ES, FE_ES = [], {}
    for k, (var, grp, cls) in {"es_core5_d24": ("d24", "core5", MAJ), "es_pgc_d24": ("d24", "pgc", MAJ), "es_bh_d24": ("d24", "bh", MAJ),
                               "es_core5_d35": ("d35", "core5", MAJ), "es_core5_c24": ("c24", "core5", MAJ),
                               "es_core5_c35": ("c35", "core5", MAJ), "es_core5_d24_all": ("d24", "core5", ALL),
                               "es_all7_d24": ("d24", "all7", MAJ), "es_pgc_d35": ("d35", "pgc", MAJ), "es_bh_d35": ("d35", "bh", MAJ),
                               "es_bh_c24": ("c24", "bh", MAJ)}.items():
        F = fe_table(ES, var, grp, cls, min_vals=1 if grp == "pgc" else 2); FE_ES[k] = F
        F.to_csv(f"{M.TB}/fe_{k}.csv", index=False)
        res = agreement(F, k)
        if res:
            out, pairs, FRr = res; AG_ES.append(out)
            pairs.to_csv(f"{M.TB}/pairs_{k}.csv", index=False); FRr.to_csv(f"{M.TB}/fam_vs_rest_{k}.csv", index=False)
            log(f"{k}: {json.dumps({a: (round(b, 4) if isinstance(b, float) else b) for a, b in out.items()})}")
    AGT = pd.concat([AGT, pd.DataFrame(AG_ES)], ignore_index=True); AGT.to_csv(f"{M.TB}/agreement.csv", index=False)
    ag = AGT.set_index("set")
    RES["agreement"] = AGT.to_dict("records")
    RES["es_pass_same_rule"] = bool(ag.loc["es_core5_d24", "M"] > ag.loc["es_core5_d24", "perm_p99"] and ag.loc["es_pgc_d24", "M"] > 0
                                    and ag.loc["es_bh_d24", "M"] > 0)
    RES["es_control_pass"] = bool(ag.loc["es_core5_c24", "M"] > ag.loc["es_core5_c24", "perm_p99"]) if ("es_core5_c24" in ag.index and np.isfinite(ag.loc["es_core5_c24", "M"])) else None
    cmp = FE["core5_d24"].merge(FE_ES["es_core5_d24"], on=["fam", "episode"], suffixes=("", "_es"))
    RES["daymean_vs_epstack_r"] = float(stats.pearsonr(cmp.value, cmp.value_es)[0])
    RES["daymean_vs_epstack_slope"] = float(stats.theilslopes(cmp.value_es, cmp.value)[0])
    for tag, DD, FEx in (("dm", D, FE["core5_d24"]), ("es", ES, FE_ES["es_core5_d24"])):
        series_battery(tag, DD, U, FEx, FE["core5_d35"] if tag == "dm" else FE_ES["es_core5_d35"])
    # Revision 1 item 5 (post hoc, added after the annual cycle was seen): agreement after removing the pooled annual cycle
    AG_DS = []
    for tag, FEd in (("dm", FE), ("es", FE_ES)):
        b0, b1 = RES[f"annual_model_{tag}"]
        for k in [k for k in FEd if k.endswith("core5_d24") or k.endswith("pgc_d24") or k.endswith("bh_d24") or k.endswith("core5_d35")]:
            F = FEd[k].merge(U[["episode", "start"]], on="episode")
            ph = 2 * np.pi * (pd.to_datetime(F.start).dt.dayofyear.to_numpy() + 10) / 365.25
            F["value"] = F.value - (b0 * np.cos(ph) + b1 * np.sin(ph))
            res = agreement(F, f"{k}_deseason" if k.startswith("es_") else f"dm_{k}_deseason")
            if res:
                AG_DS.append(res[0]); log(f"deseasoned {k}: M {res[0]['M']:.4f} p99 {res[0].get('perm_p99', np.nan):.4f} p {res[0].get('perm_p', np.nan):.4f}")
    AGT = pd.concat([AGT, pd.DataFrame(AG_DS)], ignore_index=True); AGT.to_csv(f"{M.TB}/agreement.csv", index=False)
    RES["agreement"] = AGT.to_dict("records")
    # within-episode
    FRt = pd.read_csv(f"{S.ROOT}/svi_science/results/s5/step3_front_arrivals.csv")
    FRt["front"] = pd.to_datetime(FRt.t_front).dt.normalize()
    WE = {}
    prof_all = {}
    for k, (var, grp) in {"core5_d24": ("d24", "core5"), "pgc_d24": ("d24", "pgc"), "bh_d24": ("d24", "bh"),
                          "core5_c24": ("c24", "core5"), "core5_d35": ("d35", "core5")}.items():
        X = within_episode(D, var, grp, FRt)
        prof, st = we_stats(X, VARS[var][1], M.SEED)
        st["n_front_S5"] = int((X.x_src == "S5 front").sum()); WE[k] = st; prof_all[k] = prof
        if k == "core5_d24":
            X.to_parquet(f"{M.TB}/within_episode_points.parquet", index=False)
        log(f"within-episode {k}: {st}")
    pd.DataFrame(prof_all).to_csv(f"{M.TB}/within_episode_profile.csv")
    pd.DataFrame(WE).T.to_csv(f"{M.TB}/within_episode.csv")
    w0 = WE["core5_d24"]
    RES["within"] = WE
    RES["within_consistent"] = bool((abs(w0["slope"]) > 3 * w0["slope_sig"] or abs(w0["contrast"]) > 3 * w0["contrast_sig"])
                                    and abs(WE["core5_c24"]["contrast"]) < 2 * WE["core5_c24"]["contrast_sig"]
                                    and np.sign(WE["pgc_d24"]["contrast"]) == np.sign(WE["bh_d24"]["contrast"]))
    json.dump(RES, open(f"{M.TB}/result.json", "w"), indent=1, default=float)
    log("RESULT " + json.dumps({k: v for k, v in RES.items() if k not in ("agreement", "within")}, default=float))
    open(f"{M.OUT}/ANALYSIS_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
