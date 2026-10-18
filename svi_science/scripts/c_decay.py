"""C: LFE activity at a patch after the slip front passes (SVI + Washington-NW). Rules: results/c_decay/PLAN_C.md.

Outputs results/c_decay/: cells_<run>.csv, stack_fits.csv, cell_fits.csv, bursts_<run>.csv, result.json, events_lag.parquet
"""
import json
import numpy as np
import pandas as pd
from scipy import optimize, special, stats

ROOT = "/home/jovyan/tremorferometry"
R = f"{ROOT}/svi_science/results"
OUT = f"{R}/c_decay"
rng = np.random.default_rng(20261007)
NSUR = 500
NPERM = 5000
W = 25.0  # fit window (days)


def to_pt(t):
    return pd.Series(pd.to_datetime(t).values).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.tz_localize(None).values


# ------------------------------------------------------------------ load
K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv")
EPS = {reg: K[K.region == reg].episode.tolist() for reg in ["SVI", "NW"]}
EA = pd.read_csv(f"{R}/s3/episodes_all.csv")
F = pd.read_csv(f"{R}/s0/families.csv")
q1, q2 = F.downdip_km.quantile([1 / 3, 2 / 3])
S = pd.read_csv(f"{R}/s3/sources_nw.csv")
n1, n2 = S.dd_km.quantile([1 / 3, 2 / 3])
CUTS = {"SVI": (q1, q2), "NW": (n1, n2)}


def dclass(dd, reg):
    a, b = CUTS[reg]
    return np.where(dd >= b, "deep", np.where(dd <= a, "shallow", "mid"))


def load_svi(series="data", pmin=0.75):
    E = pd.read_parquet(f"{R}/s3/events_svi.parquet")
    E = E[(E.series == series) & (E.P_same >= pmin) & E.counted_once & E.ep.isin(EPS["SVI"])].copy()
    E["patch"] = "G" + E.group.astype(str)
    return E[["time", "ep", "as_km", "dd_km", "patch", "P_same"]].assign(region="SVI", qualC=False)


def load_svi_s13():
    T = pd.read_parquet(f"{R}/s0/analysis_table_S13.parquet", columns=["series", "fam", "P_same", "keep_xfam8", "time_utc"])
    T = T[(T.series == "data") & (T.P_same >= 0.75) & T.keep_xfam8.astype(bool)].copy()
    ep = np.full(len(T), None, dtype=object); tv = T.time_utc.values
    for _, r in EA[(EA.region == "SVI") & EA.episode.isin(EPS["SVI"])].iterrows():
        m = (tv >= np.datetime64(pd.Timestamp(r.w0))) & (tv < np.datetime64(pd.Timestamp(r.w1)))
        ep[m] = r.episode
    T["ep"] = ep; T = T[T.ep.notna()]
    T = T.merge(F[["fid", "along_strike_km", "downdip_km", "group"]].rename(columns={"fid": "fam"}), on="fam")
    T = T.rename(columns={"time_utc": "time", "along_strike_km": "as_km", "downdip_km": "dd_km"})
    T["patch"] = "G" + T.group.astype(str)
    return T[["time", "ep", "as_km", "dd_km", "patch", "P_same"]].assign(region="SVI", qualC=False)


def load_nw():
    E = pd.read_parquet(f"{R}/s3/events_nw.parquet")
    E = E[(E.series == "data") & (E.P_same >= 0.75) & E.selected_v2 & E.counted_once & E.ep.isin(EPS["NW"])].copy()
    E["dcl"] = dclass(E.dd_km.values, "NW")
    E["patch"] = "B" + (np.floor(E.as_km / 10) * 10).astype(int).astype(str) + "_" + E.dcl
    return E[["time", "ep", "as_km", "dd_km", "patch", "P_same"]].assign(region="NW", qualC=(E.quality == "C").values)


def prep(E, reg):
    E = E.copy()
    E["time"] = pd.to_datetime(E.time).astype("datetime64[ns]")
    E["tpt"] = to_pt(E.time)
    E["dcl"] = dclass(E.dd_km.values, reg)
    return E.sort_values("time").reset_index(drop=True)


TR = {reg: pd.read_parquet(f"{R}/s3/tremor_{reg.lower()}.parquet") for reg in ["SVI", "NW"]}


def arrival_map(tr):
    tr = tr.assign(b=np.floor(tr.as_km / 5) * 5)
    g = tr.groupby("b").time.apply(lambda s: s.sort_values().iloc[int(0.1 * (len(s) - 1))] if len(s) >= 20 else pd.NaT)
    return g.dropna()


def add_lag(E, reg):
    E = E.copy(); E["lag"] = np.nan; E["arr"] = pd.NaT
    for ep in E.ep.unique():
        am = arrival_map(TR[reg][TR[reg].ep == ep]); m = E.ep == ep
        arr = (np.floor(E.loc[m, "as_km"] / 5) * 5).map(am)
        E.loc[m, "arr"] = arr.values
        E.loc[m, "lag"] = (E.loc[m, "time"] - arr).dt.total_seconds().values / 86400
    return E


# ------------------------------------------------------------------ shape fits
GRID = np.linspace(0, W, 20001)


def _nll(f, t):
    z = np.trapezoid(f(GRID), GRID)
    v = f(t)
    if not np.isfinite(z) or z <= 0 or np.any(v <= 0) or not np.all(np.isfinite(v)):
        return 1e12
    return -(np.sum(np.log(v)) - len(t) * np.log(z))


MODELS = {
    "omori": (lambda th: (lambda x: (x + np.exp(th[0])) ** (-th[1])), [[np.log(0.05), 0.5], [np.log(0.5), 1.0], [np.log(3), 2.0]],
              [(np.log(1e-3), np.log(30)), (-1, 6)]),
    "exponential": (lambda th: (lambda x: np.exp(-x / np.exp(th[0])) + np.exp(th[1])), [[np.log(1), -6], [np.log(5), -3], [np.log(0.3), -8]],
                    [(np.log(0.05), np.log(200)), (-20, 3)]),
    "plateau_stop": (lambda th: (lambda x: special.expit(-(x - th[0]) / np.exp(th[1]))), [[2, np.log(0.5)], [6, np.log(1)], [12, np.log(2)]],
                     [(-5, 50), (np.log(0.02), np.log(20))]),
}


def fit_shapes(t):
    t = np.asarray(t); out = {}
    for name, (mk, starts, bnds) in MODELS.items():
        best = None
        for s0 in starts:
            r = optimize.minimize(lambda th: _nll(mk(th), t), s0, method="L-BFGS-B", bounds=bnds)
            if best is None or r.fun < best.fun:
                best = r
        out[name] = dict(nll=float(best.fun), aic=float(2 * 2 + 2 * best.fun), th=best.x.tolist())
    aics = sorted((v["aic"], k) for k, v in out.items())
    return out, aics[0][1], aics[1][0] - aics[0][0]


def peak_shift(lags, lo=-2, hi=10):
    h, e = np.histogram(lags, bins=np.arange(lo, hi + 0.5, 0.5))
    return e[np.argmax(h)] + 0.25


# ------------------------------------------------------------------ cell metrics
def d80(t):
    return np.quantile(t, 0.9) - np.quantile(t, 0.1)


def cell_table(E, reg, run):
    rows = []
    for (ep, pa), g in E.groupby(["ep", "patch"]):
        td = (g.time - g.time.min()).dt.total_seconds().values / 86400
        r = dict(region=reg, run=run, ep=ep, patch=pa, n=len(g), dd=g.dd_km.mean(), as_km=g.as_km.mean(),
                 dcl=dclass(np.array([g.dd_km.mean()]), reg)[0], qualC_frac=float(g.qualC.mean()))
        if len(g) >= 30:
            r["D80"] = d80(td)
            r["D80_n30"] = float(np.median([d80(rng.choice(td, 30, replace=False)) for _ in range(200)]))
            arr = g.arr.dropna()
            if len(arr):
                r["T90"] = (g.time.quantile(0.9) - arr.median()).total_seconds() / 86400
            r["own_arr_minus_tremor_d"] = (g.time.quantile(0.05) - arr.median()).total_seconds() / 86400 if len(arr) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def perm_within(C, ycol, xcol="dd", strata=("ep",), nperm=NPERM, min_cells=3):
    C = C.dropna(subset=[ycol, xcol]).copy()
    groups = [g for _, g in C.groupby(list(strata)) if len(g) >= min_cells]
    if not groups:
        return np.nan, np.nan, 0
    def stat(gs):
        return np.nanmean([stats.spearmanr(g[0], g[1])[0] for g in gs])
    data = [(g[xcol].values, g[ycol].values) for g in groups]
    obs = stat(data)
    nul = np.array([stat([(rng.permutation(x), y) for x, y in data]) for _ in range(nperm)])
    p = (np.sum(np.abs(nul) >= abs(obs)) + 1) / (nperm + 1)
    return float(obs), float(p), len(groups)


# ------------------------------------------------------------------ Q2 bursts
def iets(t_s):
    d = np.diff(np.sort(t_s))
    return d[d <= 3 * 86400]


def bursts(t_s, g):
    t = np.sort(t_s); d = np.diff(t)
    sizes = []; cur = 1
    for x in d:
        if x < g:
            cur += 1
        else:
            sizes.append(cur); cur = 1
    sizes.append(cur)
    return np.array(sizes)


def cv(x):
    return np.std(x) / np.mean(x) if len(x) > 2 and np.mean(x) > 0 else np.nan


def surrogate(tpt_s, kind):
    if kind == "uni":
        return rng.uniform(tpt_s.min(), tpt_s.max(), len(tpt_s))
    u = 86400 if kind == "day" else 3600
    return np.floor(tpt_s / u) * u + rng.uniform(0, u, len(tpt_s))


def q2_cells(E, reg, run, nmin=30):
    rows = []; sur = {k: [] for k in ["uni", "day", "hour"]}; sizes = {"obs": [], "day": []}; normiet = {"obs": [], "day": [], "uni": []}
    for (ep, pa), g in E.groupby(["ep", "patch"]):
        if len(g) < nmin:
            continue
        ts = (pd.to_datetime(g.tpt) - pd.Timestamp("2014-01-01")).dt.total_seconds().values
        x = iets(ts); m = x.mean(); gg = 0.1 * m
        b = bursts(ts, gg)
        r = dict(region=reg, run=run, ep=ep, patch=pa, n=len(g), dd=g.dd_km.mean(), dcl=dclass(np.array([g.dd_km.mean()]), reg)[0],
                 cv=cv(x), mean_iet_min=m / 60, g_min=gg / 60, frac_burst=float(b[b >= 2].sum() / len(ts)),
                 frac_burst_10min=float((lambda bb: bb[bb >= 2].sum() / len(ts))(bursts(ts, 600))),
                 frac_burst_60min=float((lambda bb: bb[bb >= 2].sum() / len(ts))(bursts(ts, 3600))),
                 n_iet_lt5s=int((np.diff(np.sort(ts)) < 5).sum()), max_burst=int(b.max()))
        sizes["obs"].append(b[b >= 2]); normiet["obs"].append(x / m)
        for k in sur:
            cvs = []; fb = []
            for i in range(NSUR):
                s = surrogate(ts, k); xs = iets(s); cvs.append(cv(xs))
                bs = bursts(s, 0.1 * xs.mean()); fb.append(bs[bs >= 2].sum() / len(s))
                if k == "day" and i < 20:
                    sizes["day"].append(bs[bs >= 2])
                if k in ("day", "uni") and i < 5:
                    normiet[k].append(xs / xs.mean())
            sur[k].append(np.array(cvs))
            r[f"cv_{k}_med"] = float(np.median(cvs)); r[f"fb_{k}_med"] = float(np.median(fb))
        rows.append(r)
    C = pd.DataFrame(rows)
    res = {}
    if len(C):
        obs_med = np.nanmedian(C.cv)
        for k in sur:
            M = np.nanmedian(np.vstack(sur[k]), axis=0)
            res[k] = dict(obs_median_cv=float(obs_med), sur_median_cv=float(np.median(M)), sur_q99=float(np.quantile(M, 0.99)),
                          p=float((np.sum(M >= obs_med) + 1) / (len(M) + 1)))
    return C, res, {k: (np.concatenate(v) if v else np.array([])) for k, v in sizes.items()}, \
        {k: (np.concatenate(v) if v else np.array([])) for k, v in normiet.items()}


def size_fits(s):
    s = np.asarray(s); s = s[s >= 2]; n = len(s)
    if n < 20:
        return dict(n=n)
    nll_pl = lambda a: -(np.sum(-a * np.log(s)) - n * np.log(special.zeta(a, 2)))
    a = optimize.minimize_scalar(nll_pl, bounds=(1.01, 10), method="bounded").x
    mq = np.mean(s - 2); q = mq / (1 + mq)
    lpl = -a * np.log(s) - np.log(special.zeta(a, 2))
    lex = np.log(1 - q) + (s - 2) * np.log(q) if q > 0 else np.zeros(n)
    d = lpl - lex; Rr = d.sum(); sd = d.std() * np.sqrt(n)
    z = Rr / sd if sd > 0 else 0.0
    p = float(2 * stats.norm.sf(abs(z)))
    verdict = "power law" if (Rr > 0 and p < 0.01) else ("exponential" if (Rr < 0 and p < 0.01) else "undecided")
    return dict(n=n, alpha=float(a), q=float(q), mean=float(s.mean()), max=int(s.max()), P_ge5=float((s >= 5).mean()),
                LLR=float(Rr), z=float(z), p=p, verdict=verdict)


# ------------------------------------------------------------------ Q3 re-activation
def reactivation(days, counts):
    """days sorted ints with counts; returns (bool, first_day, gap_start, return_day)."""
    A = days[counts >= 3]
    if len(A) < 2:
        return False, None, None
    for i in range(len(A) - 1):
        if A[i + 1] - A[i] - 1 >= 3:
            between = (days > A[i]) & (days < A[i + 1])
            if counts[between].sum() == 0:
                return True, int(A[i]), int(A[i + 1])
    return False, None, None


def day_counts(tpt):
    d = (pd.to_datetime(tpt).values.astype("datetime64[D]") - np.datetime64("2014-01-01")).astype(int)
    u, c = np.unique(d, return_counts=True)
    return u, c


def q3_cells(E, reg, run, nmin=10):
    rows = []; Tr = TR[reg].copy(); Tr["tpt"] = to_pt(Tr.time)
    for ep, ge in E.groupby("ep"):
        alltimes = ge.tpt.values
        tre = Tr[Tr.ep == ep]
        for pa, g in ge.groupby("patch"):
            if len(g) < nmin:
                continue
            u, c = day_counts(g.tpt)
            ok, d0, d1 = reactivation(u, c)
            exp = []
            for _ in range(NSUR):
                s = rng.choice(alltimes, len(g), replace=False); uu, cc = day_counts(s); exp.append(reactivation(uu, cc)[0])
            r = dict(region=reg, run=run, ep=ep, patch=pa, n=len(g), dd=g.dd_km.mean(), dcl=dclass(np.array([g.dd_km.mean()]), reg)[0],
                     react=ok, react_expected=float(np.mean(exp)))
            if ok:
                pk = u[np.argmax(c)]
                r.update(main_day=int(d0), return_day=int(d1), return_after_peak_d=int(d1 - pk),
                         return_PT=str(np.datetime64("2014-01-01") + np.timedelta64(int(d1), "D")))
                lo, hi = g.as_km.min() - 5, g.as_km.max() + 5
                tt = tre[(tre.as_km >= lo) & (tre.as_km <= hi)]
                tu, tc = day_counts(tt.tpt) if len(tt) else (np.array([], int), np.array([], int))
                cmap = dict(zip(tu, tc))
                gapdays = [cmap.get(d, 0) for d in range(d0 + 1, d1)]
                r["tremor_return"] = bool(cmap.get(d1, 0) >= 3 and sum(np.array(gapdays) <= 1) >= 3)
            rows.append(r)
    C = pd.DataFrame(rows)
    if len(C) and C.react.any():
        C["sync_other_patch"] = False
        for ep, g in C[C.react].groupby("ep"):
            for i, r in g.iterrows():
                C.loc[i, "sync_other_patch"] = bool((np.abs(g.drop(index=i).return_day - r.return_day) <= 1).any())
    return C


def q3_test(C):
    C = C.copy(); C["y"] = C.react.astype(float)
    C["ct"] = pd.qcut(C.n.rank(method="first"), 3, labels=False)
    obs = stats.spearmanr(C.y, C.dd)[0]
    nul = []; idx = [g.index.values for _, g in C.reset_index(drop=True).groupby(["ep", "ct"])]
    ddv = C.dd.values; yv = C.y.values
    for _ in range(NPERM):
        dd = ddv.copy()
        for ii in idx:
            dd[ii] = rng.permutation(ddv[ii])
        nul.append(stats.spearmanr(yv, dd)[0])
    nul = np.array(nul)
    p = (np.sum(np.abs(nul) >= abs(obs)) + 1) / (NPERM + 1)
    fr = C.groupby("dcl").agg(n_cells=("y", "size"), frac=("y", "mean"), expected=("react_expected", "mean")).to_dict("index")
    return dict(rho=float(obs), p=float(p), n_cells=int(len(C)), frac_obs=float(C.y.mean()),
                frac_expected_regionwide=float(C.react_expected.mean()), by_class=fr)


# ------------------------------------------------------------------ main
def main():
    res = {"cuts_dd_km": {k: [float(v[0]), float(v[1])] for k, v in CUTS.items()}}
    runs = {("SVI", "headline"): load_svi(), ("SVI", "p050"): load_svi(pmin=0.5), ("SVI", "S13"): load_svi_s13(),
            ("NW", "headline"): load_nw(), ("SVI", "scrambled"): load_svi("c1", 0.5)}
    D = {}
    for (reg, run), E in runs.items():
        E = add_lag(prep(E, reg), reg); D[(reg, run)] = E
        res[f"{reg}_{run}_n"] = int(len(E)); res[f"{reg}_{run}_n_no_arrival"] = int(E.lag.isna().sum())
    pd.concat([D[k].assign(run=k[1]) for k in D]).to_parquet(f"{OUT}/events_lag.parquet")

    # ---- Q1 stacks + fits
    srows = []
    for reg in ["SVI", "NW"]:
        for run in (["headline", "S13"] if reg == "SVI" else ["headline"]):
            E = D[(reg, run)].dropna(subset=["lag"])
            for cl in ["all", "shallow", "mid", "deep"]:
                x = E.lag.values if cl == "all" else E[E.dcl == cl].lag.values
                tpk = peak_shift(x); t = x - tpk; t = t[(t >= 0) & (t <= W)]
                fits, win, dA = fit_shapes(t)
                srows.append(dict(region=reg, run=run, dclass=cl, n=len(t), t_peak_d=tpk, winner=win, dAIC=dA,
                                  frac_before_arrival=float((x < 0).mean()), median_lag=float(np.median(x)),
                                  lag_q90=float(np.quantile(x, 0.9)),
                                  **{f"{k}_aic": v["aic"] for k, v in fits.items()},
                                  **{f"{k}_th": json.dumps([round(z, 4) for z in v["th"]]) for k, v in fits.items()}))
    SF = pd.DataFrame(srows); SF.to_csv(f"{OUT}/stack_fits.csv", index=False)

    crow = []
    for (reg, run), E in D.items():
        if run not in ("headline", "S13"):
            continue
        for (ep, pa), g in E.groupby(["ep", "patch"]):
            td = (g.time - g.time.min()).dt.total_seconds().values / 86400
            tpk = peak_shift(td, 0, td.max() + 0.5)
            t = td - tpk; t = t[(t >= 0) & (t <= W)]
            if len(t) < 30:
                continue
            fits, win, dA = fit_shapes(t)
            crow.append(dict(region=reg, run=run, ep=ep, patch=pa, n_fit=len(t), dd=g.dd_km.mean(),
                             dcl=dclass(np.array([g.dd_km.mean()]), reg)[0], winner=win if dA >= 2 else "undecided", dAIC=dA,
                             best=win, omori_p=fits["omori"]["th"][1], exp_tau=float(np.exp(fits["exponential"]["th"][0])),
                             plateau_T=fits["plateau_stop"]["th"][0], plateau_w=float(np.exp(fits["plateau_stop"]["th"][1]))))
    CF = pd.DataFrame(crow); CF.to_csv(f"{OUT}/cell_fits.csv", index=False)
    res["shape"] = {}
    for (reg, run), g in CF.groupby(["region", "run"]):
        dec = g[g.winner != "undecided"].winner.value_counts()
        st = SF[(SF.region == reg) & (SF.run == run) & (SF.dclass == "all")].iloc[0]
        top = dec.index[0] if len(dec) else None
        cond_cells = len(dec) > 0 and dec.iloc[0] >= 0.4 * dec.sum() and (len(dec) == 1 or dec.iloc[0] >= 1.5 * dec.iloc[1])
        res["shape"][f"{reg}_{run}"] = dict(stack_winner=st.winner, stack_dAIC=float(st.dAIC), n_cells=int(len(g)),
                                            cell_winners=dec.to_dict(), cell_top=top, cells_condition=bool(cond_cells),
                                            pass_ruleS=bool(st.dAIC >= 10 and cond_cells and top == st.winner),
                                            by_dclass={c: gg[gg.winner != "undecided"].winner.value_counts().to_dict()
                                                       for c, gg in g.groupby("dcl")})

    # ---- Q1 duration
    res["duration"] = {}; CT = []
    for (reg, run), E in D.items():
        if run not in ("headline", "S13", "p050"):
            continue
        C = cell_table(E, reg, run); CT.append(C)
        rho, p, ng = perm_within(C, "D80_n30")
        rho2, p2, _ = perm_within(C, "T90")
        Cd = C.dropna(subset=["D80_n30"])
        res["duration"][f"{reg}_{run}"] = dict(rho_D80n30_dd=rho, p=p, n_episodes=ng, n_cells=int(len(Cd)),
                                               rho_T90_dd=rho2, p_T90=p2,
                                               median_D80_by_class=Cd.groupby("dcl").D80.median().to_dict(),
                                               median_D80n30_by_class=Cd.groupby("dcl").D80_n30.median().to_dict(),
                                               median_T90_by_class=Cd.groupby("dcl").T90.median().to_dict(),
                                               median_D80=float(Cd.D80.median()), median_D80n30=float(Cd.D80_n30.median()))
    CT = pd.concat(CT); CT.to_csv(f"{OUT}/cells_duration.csv", index=False)

    # ---- Q1 episode size
    L = pd.read_csv(f"{R}/s1/leakage_svi_vs_wa.csv")
    res["size"] = {}
    for reg in ["SVI", "NW"]:
        C = CT[(CT.region == reg) & (CT.run == "headline")].dropna(subset=["D80_n30"])
        ep_med = C.groupby("ep").D80_n30.median()
        k = K[K.region == reg].set_index("episode")
        out = {}
        sizes = {"n_tremor": k.n_tremor, "n_lfe": k.n_headline, "ep_duration_lfe": k.lfe_duration_d}
        if reg == "SVI":
            sizes["gnss_D"] = L.set_index("episode").D
        else:
            m = {}
            for _, r in L.iterrows():
                for e in str(r.nw_overlap).split(";"):
                    m[e] = r.D_WA
            sizes["gnss_D_WA"] = pd.Series(m)
        for nm, sz in sizes.items():
            x = pd.concat([ep_med, sz.rename("s")], axis=1, join="inner").dropna()
            if len(x) < 5:
                continue
            obs = stats.spearmanr(x.iloc[:, 0], x.s)[0]
            nul = np.array([stats.spearmanr(rng.permutation(x.iloc[:, 0].values), x.s)[0] for _ in range(10000)])
            out[nm] = dict(rho=float(obs), p=float((np.sum(np.abs(nul) >= abs(obs)) + 1) / 10001), n=int(len(x)))
        res["size"][reg] = out
        res["size"][reg]["episode_median_D80n30"] = ep_med.round(2).to_dict()

    # ---- Q2
    res["q2"] = {}; QC = []; SZ = {}; NI = {}
    for (reg, run), E in D.items():
        if run == "p050":
            continue
        nmin = 30
        if run == "scrambled":
            nmin = 30 if (E.groupby(["ep", "patch"]).size() >= 30).sum() >= 5 else 10
        C, r2, sz, ni = q2_cells(E, reg, run, nmin)
        QC.append(C); SZ[(reg, run)] = sz; NI[(reg, run)] = ni
        r2["nmin"] = nmin; r2["n_cells"] = int(len(C))
        if len(C):
            r2["median_cv_ratio_to_Sday"] = float(np.nanmedian(C.cv / C.cv_day_med))
            r2["median_frac_burst"] = float(C.frac_burst.median()); r2["median_frac_burst_Sday"] = float(C.fb_day_med.median())
            r2["median_frac_burst_uni"] = float(C.fb_uni_med.median()); r2["median_frac_burst_hour"] = float(C.fb_hour_med.median())
            r2["median_frac_burst_60min"] = float(C.frac_burst_60min.median()); r2["median_frac_burst_10min"] = float(C.frac_burst_10min.median())
            r2["median_mean_iet_min"] = float(C.mean_iet_min.median()); r2["iet_lt5s_total"] = int(C.n_iet_lt5s.sum())
            if run in ("headline", "S13"):
                r2["depth_cv"] = perm_within(C, "cv")[:2]; r2["depth_frac_burst"] = perm_within(C, "frac_burst")[:2]
                r2["cv_by_class"] = C.groupby("dcl").cv.median().to_dict()
                r2["fb_by_class"] = C.groupby("dcl").frac_burst.median().to_dict()
                r2["cv_ratio_by_class"] = (C.cv / C.cv_day_med).groupby(C.dcl).median().to_dict()
        res["q2"][f"{reg}_{run}"] = r2
    QC = pd.concat(QC); QC.to_csv(f"{OUT}/cells_bursts.csv", index=False)
    res["burst_sizes"] = {}
    for (reg, run), sz in SZ.items():
        if run not in ("headline", "S13"):
            continue
        res["burst_sizes"][f"{reg}_{run}"] = dict(obs=size_fits(sz["obs"]), Sday=size_fits(sz["day"]))
    np.savez(f"{OUT}/sizes_iets.npz", **{f"{r}_{u}_size_{k}": v for (r, u), d in SZ.items() for k, v in d.items()},
             **{f"{r}_{u}_iet_{k}": v for (r, u), d in NI.items() for k, v in d.items()})
    res["burst_sizes_by_class"] = {}
    for reg in ["SVI", "NW"]:
        pass  # per-class sizes from cells_bursts not stored per event; pooled per region only (see PLAN)

    # ---- Q3
    res["q3"] = {}; RC = []
    for (reg, run), E in D.items():
        if run not in ("headline", "S13"):
            continue
        C = q3_cells(E, reg, run); RC.append(C)
        t = q3_test(C)
        if C.react.any():
            rr = C[C.react]
            t.update(n_react=int(len(rr)), median_return_after_peak_d=float(rr.return_after_peak_d.median()),
                     frac_tremor_return=float(rr.tremor_return.mean()), frac_sync=float(rr.sync_other_patch.mean()),
                     react_by_episode=rr.ep.value_counts().to_dict())
        res["q3"][f"{reg}_{run}"] = t
    RC = pd.concat(RC); RC.to_csv(f"{OUT}/cells_reactivation.csv", index=False)
    json.dump(res, open(f"{OUT}/result.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
