"""D. Fixed-station cross-era test for WA-NW LFE activity (rules: svi_science/results/d_fixed_station/PLAN.md).

python svi_science/scripts/d_fixed_station.py once     # once-per-source counts from the P files (parallel over months)
python svi_science/scripts/d_fixed_station.py stats    # measures 1-5, tables, figures
"""
import glob
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = "/home/jovyan/tremorferometry"
WA = f"{ROOT}/lfe_extension_wa/results"
OUT = f"{ROOT}/svi_science/results/d_fixed_station"
FIG = f"{ROOT}/svi_science/figures"
NETS = {  # name: (post dir, P dir, P net value, strict available)
    "prod": (f"{WA}/wa_prod/post_v2", f"{WA}/wa_prod/conf/P/prod", "core", False),
    "core9": (f"{WA}/wa_fixed10/post_v2_core9", f"{WA}/wa_fixed10/conf/P/prod", "core", True),
    "all10": (f"{WA}/wa_fixed10/post_v2_all10", f"{WA}/wa_fixed10/conf/P/prod", "all", True),
}
ERAS = {"E1": ("2007-01-01", "2011-12-31"), "E2": ("2012-01-01", "2014-12-31"), "E3": ("2015-01-01", "2026-09-15")}
COVER = "2009-08-06"
PAIRS = [("E1", "E3"), ("E2", "E3"), ("E1", "E2")]
NB, NP = 2000, 10000
rng = np.random.default_rng(20261007)
FAM = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv")
GRP = dict(zip(FAM.fid.astype(int), FAM.source_group.fillna("none").astype(str)))


def log(m):
    print(m, flush=True)
    with open(f"{OUT}/d.log", "a") as f:
        f.write(m + "\n")


# ------------------------------------------------------------------------------------------------ once per source
def _merge_count(t):
    t = np.sort(t)
    return 1 + int((np.diff(t) > 2.0).sum()) if len(t) else 0


def once_month(args):
    path, netval = args
    D = pd.read_parquet(path, columns=["date", "net", "series", "fid", "t", "reg", "P_own_WAVE"],
                        filters=[("reg", "==", "NW"), ("net", "==", netval), ("P_own_WAVE", ">=", 0.75)])
    if not len(D):
        return None
    D["grp"] = D.fid.astype(int).map(GRP).fillna("none")
    raw = D.groupby(["fid", "date", "series"]).size().rename("n_raw").reset_index()
    raw["grp"] = raw.fid.astype(int).map(GRP).fillna("none")
    g = D.groupby(["grp", "date", "series"]).t.agg(lambda x: _merge_count(x.to_numpy())).rename("g_merged").reset_index()
    gr = D.groupby(["grp", "date", "series"]).size().rename("g_raw").reset_index()
    raw = raw.merge(g, on=["grp", "date", "series"]).merge(gr, on=["grp", "date", "series"])
    raw["n_once"] = raw.n_raw * raw.g_merged / raw.g_raw
    return raw


def once():
    for name, (_, pdir, netval, _) in NETS.items():
        if name == "all10" and False:
            continue
        files = sorted(glob.glob(f"{pdir}/P_*.parquet"))
        with Pool(12) as p:
            parts = [x for x in p.map(once_month, [(f, netval) for f in files]) if x is not None]
        X = pd.concat(parts, ignore_index=True)
        X.to_parquet(f"{OUT}/once_{name}.parquet", index=False)
        log(f"[once] {name}: {len(X)} fid-day-series rows, raw {X.n_raw.sum():.0f}, once {X.n_once.sum():.1f}")


# ------------------------------------------------------------------------------------------------ load
def era_of(d):
    d = np.asarray(d).astype(str)
    return np.select([d <= ERAS["E1"][1], d <= ERAS["E2"][1], d <= ERAS["E3"][1]], ["E1", "E2", "E3"], "")


def load(name):
    post, _, _, has_strict = NETS[name]
    cols = ["unit_kind", "unit", "label", "date", "region", "era", "day_class", "n_trem30"] + \
        [f"n{t}_{s}" for t in ("050", "075", "085") for s in ("data", "rev", "c1")]
    R = pd.read_parquet(f"{post}/rates_daily_v2.parquet", columns=cols,
                        filters=[("label", "==", "own"), ("unit_kind", "==", "family"), ("region", "==", "NW")])
    R["fid"] = R.unit.astype(int)
    if has_strict:
        S = pd.read_parquet(f"{post}/unit_day_strict.parquet")
        S = S[S.unit_kind == "family"]
        R = R.merge(S[["unit", "date", "n_fixed", "n_fixed_live", "strict"]], on=["unit", "date"], how="left")
        R["live_frac"] = R.n_fixed_live / R.n_fixed
        R["strict"] = R.strict.fillna(False).astype(bool)
    else:
        L = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/live_core_by_day.csv")
        cs = json.load(open(f"{WA}/wa_prod/setup/core_sets.json"))
        size = {e: len(cs[e]) for e in ("E1", "E2", "E3")}
        L["era"] = era_of(L.date)
        L["live_frac"] = L.n_core_live / L.era.map(size)
        R = R.merge(L[["date", "live_frac"]], on="date", how="left")
        R["strict"] = False
    O = pd.read_parquet(f"{OUT}/once_{name}.parquet")
    O["unit"] = O.fid.astype(int).astype(str)
    Op = O.pivot_table(index=["unit", "date"], columns="series", values=["n_once", "n_raw"], fill_value=0)
    Op.columns = [f"{a}_{b}" for a, b in Op.columns]
    Op = Op.reset_index()
    R = R.merge(Op, on=["unit", "date"], how="left")
    for c in [c for c in R.columns if c.startswith(("n_once", "n_raw"))]:
        R[c] = R[c].fillna(0.0)
    for s in ("data", "rev", "c1"):
        for c in (f"n_once_{s}", f"n_raw_{s}"):
            if c not in R:
                R[c] = 0.0
    R["net"] = R.n075_data - R.n075_rev
    R["net_scr"] = R.n075_data - R.n075_c1
    R["net_once"] = R.n_once_data - R.n_once_rev
    R["net050"] = R.n050_data - R.n050_rev
    R["net085"] = R.n085_data - R.n085_rev
    R["one"] = 1.0
    dd = FAM.set_index("fid")
    md = FAM.downdip_km.median(); mz = FAM.dep.median()
    R["deep"] = R.fid.map(dd.downdip_km > md)
    R["deep_z"] = R.fid.map(dd.dep > mz)
    T = R[R.day_class == "tremor"].copy()
    log(f"[load] {name}: tremor unit-days {len(T)} (all rows {len(R)}); once raw vs n075 data "
        f"{T.n_raw_data.sum():.0f} vs {T.n075_data.sum():.0f}, rev {T.n_raw_rev.sum():.0f} vs {T.n075_rev.sum():.0f}")
    return T


def clusters(dates):
    u = np.array(sorted(set(dates)))
    dt = np.diff(pd.to_datetime(u).values.astype("datetime64[D]").astype(int))
    cid = np.concatenate([[0], np.cumsum(dt > 5)])
    return dict(zip(u, cid))


# ------------------------------------------------------------------------------------------------ stats core
def csum(T, num, den, mask):
    """per (era, cluster) sums."""
    g = T[mask].groupby(["era", "cl"]).agg(num=(num, "sum"), den=(den, "sum")).reset_index()
    return {e: g[g.era == e][["num", "den"]].to_numpy() for e in ERAS}


def rate(a):
    return a[:, 0].sum() / a[:, 1].sum() if len(a) and a[:, 1].sum() > 0 else np.nan


def boot_rate(a, n=NB):
    if not len(a):
        return np.full(n, np.nan)
    ii = rng.integers(0, len(a), (n, len(a)))
    return a[ii, 0].sum(1) / a[ii, 1].sum(1)


def perm_p(a, b, n=NP):
    if not len(a) or not len(b):
        return np.nan
    obs = abs(np.log(rate(a) / rate(b))) if rate(a) > 0 and rate(b) > 0 else np.nan
    if not np.isfinite(obs):
        return np.nan
    X = np.vstack([a, b]); na = len(a); cnt = 0
    for _ in range(n):
        p = rng.permutation(len(X))
        ra, rb = rate(X[p[:na]]), rate(X[p[na:]])
        v = abs(np.log(ra / rb)) if ra > 0 and rb > 0 else np.inf
        cnt += v >= obs - 1e-12
    return (cnt + 1) / (n + 1)


def perm_p_ds(A, B, n=NP):
    """A, B: dicts deep/shallow -> cluster arrays aligned on clusters (num_d, den_d, num_s, den_s) for eras a and b."""
    def lr(X):
        rd = X[:, 0].sum() / X[:, 1].sum(); rs = X[:, 2].sum() / X[:, 3].sum()
        return np.log(rd / rs) if rd > 0 and rs > 0 else np.nan
    obs = abs(lr(A) - lr(B))
    X = np.vstack([A, B]); na = len(A); cnt = 0
    for _ in range(n):
        p = rng.permutation(len(X))
        v = abs(lr(X[p[:na]]) - lr(X[p[na:]]))
        cnt += (not np.isfinite(v)) or v >= obs - 1e-12
    return (cnt + 1) / (n + 1), obs


# ------------------------------------------------------------------------------------------------ stats driver
VARIANTS = {"rev075": "net", "scr075": "net_scr", "once_rev075": "net_once", "rev050": "net050", "rev085": "net085"}


def stats():
    nets = {n: load(n) for n in NETS}
    cmap = clusters(np.concatenate([T.date.to_numpy() for T in nets.values()]))
    checks = {}
    rows, rat = [], []
    strata, wls, dsrows = [], [], []
    for name, T in nets.items():
        T["cl"] = T.date.map(cmap)
        post = NETS[name][0]
        daysets = {"all": np.ones(len(T), bool)}
        if NETS[name][3]:
            daysets["strict"] = T.strict.to_numpy()
        # acceptance: reproduce post_v2 summary (wave view, NW tremor, cut_0.75)
        for ds in daysets:
            S = pd.read_csv(f"{post}/{'strict/' if ds == 'strict' else ''}rates_summary_v2.csv")
            S = S[(S.view == "wave_classrules") & (S.region == "NW") & (S.day_class == "tremor") & (S.column == "cut_0.75")]
            for e in ERAS:
                m = daysets[ds] & (T.era == e).to_numpy()
                mine = T.net[m].sum() / m.sum()
                ref = float(S[S.era == e].net_per_unit_day.iloc[0])
                checks[f"{name}_{ds}_{e}"] = dict(mine=mine, post_v2=ref, unit_days=int(m.sum()),
                                                  ref_unit_days=int(S[S.era == e].unit_days.iloc[0]), ok=bool(abs(mine - ref) < 1e-9))
        post09 = (T.date >= COVER).to_numpy()
        for ds, dm in daysets.items():
            for var, col in VARIANTS.items():
                meas = {"m1_per_famday": (col, "one", dm), "m1_per_famday_post2009": (col, "one", dm & post09),
                        "m2_per_epicentre": (col, "n_trem30", dm & post09)}
                if var in ("rev075", "once_rev075"):
                    meas["m1_deep"] = (col, "one", dm & T.deep.fillna(False).to_numpy().astype(bool))
                    meas["m1_shallow"] = (col, "one", dm & (T.deep == False).to_numpy())
                for mname, (num, den, mask) in meas.items():
                    C = csum(T, num, den, mask)
                    bs = {}
                    for e in ERAS:
                        r = rate(C[e]); b = boot_rate(C[e]); bs[e] = b
                        rows.append(dict(network=name, days=ds, variant=var, measure=mname, era=e, rate=r,
                                         lo95=np.nanpercentile(b, 2.5), hi95=np.nanpercentile(b, 97.5),
                                         num=C[e][:, 0].sum(), den=C[e][:, 1].sum(), clusters=len(C[e]),
                                         data=T.loc[mask & (T.era == e).to_numpy(), "n075_data"].sum()))
                    for a, b_ in PAIRS:
                        rr = bs[a] / bs[b_]
                        lr = np.log(rr[np.isfinite(rr) & (rr > 0)])
                        rat.append(dict(network=name, days=ds, variant=var, measure=mname, pair=f"{a}/{b_}",
                                        ratio=rate(C[a]) / rate(C[b_]),
                                        lo95=np.exp(np.percentile(lr, 2.5)), hi95=np.exp(np.percentile(lr, 97.5)),
                                        lo99=np.exp(np.percentile(lr, 0.5)), hi99=np.exp(np.percentile(lr, 99.5)),
                                        p_perm=perm_p(C[a], C[b_]) if var in ("rev075", "once_rev075", "scr075") else np.nan))
            # deep / shallow permutation (rev075 and once)
            for var in ("rev075", "once_rev075"):
                col = VARIANTS[var]
                for dname, dcol in (("downdip", "deep"), ("hyp2000_depth", "deep_z")):
                    X = T[dm & T[dcol].notna().to_numpy()].copy()
                    X["nd"] = np.where(X[dcol] == True, X[col], 0.0); X["ud"] = np.where(X[dcol] == True, 1.0, 0.0)
                    X["ns"] = np.where(X[dcol] == False, X[col], 0.0); X["us"] = np.where(X[dcol] == False, 1.0, 0.0)
                    g = X.groupby(["era", "cl"])[["nd", "ud", "ns", "us"]].sum().reset_index()
                    g = g[(g.ud > 0) & (g.us > 0)]
                    A = {e: g[g.era == e][["nd", "ud", "ns", "us"]].to_numpy() for e in ERAS}
                    for e in ERAS:
                        a = A[e]; ds_ = (a[:, 0].sum() / a[:, 1].sum()) / (a[:, 2].sum() / a[:, 3].sum())
                        ii = rng.integers(0, len(a), (NB, len(a)))
                        bb = (a[ii, 0].sum(1) / a[ii, 1].sum(1)) / (a[ii, 2].sum(1) / a[ii, 3].sum(1))
                        dsrows.append(dict(network=name, days=ds, variant=var, split=dname, era=e, deep_rate=a[:, 0].sum() / a[:, 1].sum(),
                                           shallow_rate=a[:, 2].sum() / a[:, 3].sum(), deep_over_shallow=ds_,
                                           lo95=np.nanpercentile(bb, 2.5), hi95=np.nanpercentile(bb, 97.5), clusters=len(a)))
                    for a_, b_ in PAIRS:
                        p, _ = perm_p_ds(A[a_], A[b_], n=4000)
                        dsrows.append(dict(network=name, days=ds, variant=var, split=dname, era=f"{a_}vs{b_}", p_perm=p))
        # n_sta_day covariate: strata + day-level WLS (all days, measure 1, rev075)
        lf = T.live_frac.to_numpy()
        st = np.select([lf >= 0.999, lf >= 0.8], ["1.0", "0.8-1"], "<0.8")
        T["stratum"] = st
        for (e, s), g in T.groupby(["era", "stratum"]):
            strata.append(dict(network=name, era=e, live_frac_stratum=s, unit_days=len(g), net=g.net.sum(),
                               rate=g.net.sum() / len(g), net_once=g.net_once.sum() / len(g)))
        Dd = T.groupby("date").agg(net=("net", "sum"), u=("one", "sum"), lf=("live_frac", "mean"), cl=("cl", "first"),
                                   era=("era", "first")).reset_index()
        Dd = Dd[Dd.lf.notna()]
        xbar = np.average(Dd.lf, weights=Dd.u)

        def fit(D):
            X = np.column_stack([np.ones(len(D)), (D.era == "E1"), (D.era == "E2"), D.lf - xbar]).astype(float)
            w = D.u.to_numpy(); y = (D.net / D.u).to_numpy()
            bta = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
            return bta
        b0 = fit(Dd)
        cls_by_era = {e: Dd[Dd.era == e].cl.unique() for e in ERAS}
        idx = {c: np.where(Dd.cl.to_numpy() == c)[0] for c in Dd.cl.unique()}
        bb = []
        for _ in range(1000):
            take = np.concatenate([idx[c] for e in ERAS for c in rng.choice(cls_by_era[e], len(cls_by_era[e]))])
            bb.append(fit(Dd.iloc[take]))
        bb = np.array(bb)
        for k, e in ((1, "E1"), (2, "E2")):
            r = (b0[0] + b0[k]) / b0[0]; rb = (bb[:, 0] + bb[:, k]) / bb[:, 0]
            wls.append(dict(network=name, term=f"{e}/E3 at mean live fraction", value=r, lo95=np.percentile(rb, 2.5),
                            hi95=np.percentile(rb, 97.5)))
        wls.append(dict(network=name, term="slope per +0.1 live fraction (net/famday)", value=b0[3] / 10,
                        lo95=np.percentile(bb[:, 3] / 10, 2.5), hi95=np.percentile(bb[:, 3] / 10, 97.5)))
        wls.append(dict(network=name, term="pooled mean live fraction", value=xbar))
        log(f"[stats] {name} done")
    pd.DataFrame(rows).to_csv(f"{OUT}/era_rates.csv", index=False)
    pd.DataFrame(rat).to_csv(f"{OUT}/era_ratios.csv", index=False)
    pd.DataFrame(strata).to_csv(f"{OUT}/strata_live.csv", index=False)
    pd.DataFrame(wls).to_csv(f"{OUT}/wls_covariate.csv", index=False)
    pd.DataFrame(dsrows).to_csv(f"{OUT}/deep_shallow.csv", index=False)
    json.dump(checks, open(f"{OUT}/checks.json", "w"), indent=1, default=float)
    episodes(nets)


def episodes(nets):
    E = pd.read_csv(f"{WA}/episode_rates/episodes_nw_final.csv")
    rows = []
    for name, T in nets.items():
        dsets = {"all": np.ones(len(T), bool)}
        if NETS[name][3]:
            dsets["strict"] = T.strict.to_numpy()
        d = T.date.to_numpy()
        for r in E.itertuples():
            w = (d >= r.start) & (d <= r.end)
            for ds, dm in dsets.items():
                m = w & dm
                g = T[m]
                rows.append(dict(network=name, days=ds, episode=r.episode, start=r.start, end=r.end, major=r.major,
                                 kind=r.kind, era=era_of([r.start])[0], unit_days=int(m.sum()), net=g.net.sum(),
                                 net_once=g.net_once.sum(), net_scr=g.net_scr.sum(), data=g.n075_data.sum(),
                                 rev=g.n075_rev.sum(), c1=g.n075_c1.sum(), epicentres=g.n_trem30.sum(),
                                 m1=g.net.sum() / m.sum() if m.sum() else np.nan,
                                 m1_once=g.net_once.sum() / m.sum() if m.sum() else np.nan,
                                 m2=g.net.sum() / g.n_trem30.sum() if r.kind == "pnsn" and g.n_trem30.sum() > 0 else np.nan,
                                 live_frac=g.live_frac.mean()))
    S = pd.DataFrame(rows)
    S.to_csv(f"{OUT}/episode_series.csv", index=False)
    tr = []
    for (name, ds), g in S.groupby(["network", "days"]):
        for meas in ("m1", "m1_once", "m2"):
            x = g[g.major & (g.unit_days >= 100) & g[meas].notna()].copy()
            nonpos = int((x[meas] <= 0).sum())
            x = x[x[meas] > 0]
            if len(x) < 5:
                continue
            yr = pd.to_datetime(x.start).dt.year + pd.to_datetime(x.start).dt.dayofyear / 365.25
            y = np.log(x[meas].to_numpy()); w = x.unit_days.to_numpy().astype(float); t = yr.to_numpy()
            from scipy.stats import spearmanr
            rho = spearmanr(t, y).statistic
            pr = np.array([spearmanr(t, rng.permutation(y)).statistic for _ in range(NP)])
            p = (np.sum(np.abs(pr) >= abs(rho) - 1e-12) + 1) / (NP + 1)

            def slope(tt, yy, ww):
                X = np.column_stack([np.ones(len(tt)), tt - 2016])
                return np.linalg.lstsq(X * np.sqrt(ww)[:, None], yy * np.sqrt(ww), rcond=None)[0][1]
            s0 = slope(t, y, w)
            bs = []
            for _ in range(NB):
                ii = rng.integers(0, len(t), len(t)); bs.append(slope(t[ii], y[ii], w[ii]))
            tr.append(dict(network=name, days=ds, measure=meas, episodes=len(x), nonpositive_dropped=nonpos,
                           first=x.start.min(), last=x.start.max(), spearman_rho=rho, p_perm=p,
                           slope_pct_per_yr=100 * (np.exp(s0) - 1), lo95=100 * (np.exp(np.percentile(bs, 2.5)) - 1),
                           hi95=100 * (np.exp(np.percentile(bs, 97.5)) - 1)))
            # E3-only trend (one network era for production)
            x3 = x[x.era == "E3"]
            if len(x3) >= 5:
                yr3 = pd.to_datetime(x3.start).dt.year + pd.to_datetime(x3.start).dt.dayofyear / 365.25
                y3 = np.log(x3[meas].to_numpy()); t3 = yr3.to_numpy(); w3 = x3.unit_days.to_numpy().astype(float)
                rho3 = spearmanr(t3, y3).statistic
                pr3 = np.array([spearmanr(t3, rng.permutation(y3)).statistic for _ in range(NP)])
                bs3 = []
                for _ in range(NB):
                    ii = rng.integers(0, len(t3), len(t3)); bs3.append(slope(t3[ii], y3[ii], w3[ii]))
                tr.append(dict(network=name, days=ds, measure=meas + "_E3only", episodes=len(x3), nonpositive_dropped=nonpos,
                               first=x3.start.min(), last=x3.start.max(), spearman_rho=rho3,
                               p_perm=(np.sum(np.abs(pr3) >= abs(rho3) - 1e-12) + 1) / (NP + 1),
                               slope_pct_per_yr=100 * (np.exp(slope(t3, y3, w3)) - 1),
                               lo95=100 * (np.exp(np.percentile(bs3, 2.5)) - 1), hi95=100 * (np.exp(np.percentile(bs3, 97.5)) - 1)))
    pd.DataFrame(tr).to_csv(f"{OUT}/episode_trend.csv", index=False)
    log("[episodes] done")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    {"once": once, "stats": stats}[sys.argv[1]]()
