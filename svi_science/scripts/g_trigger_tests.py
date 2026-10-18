#!/usr/bin/env python
"""g_trigger_tests (2026-10-11) -- Study G tests 1-3 (PLAN.md §3-§5) for SVI and WA-NW.
Reads results/g_trigger/tables/*; writes tables/{eq_list,inst_trials_*,inst_stacks,delayed_*,dose,ets_onsets}.csv
and tables/g_summary.json.  CPU only, single process, < 20 GB.

    python svi_science/scripts/g_trigger_tests.py
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger"
T = f"{G}/tables"
NS = 10 ** 9
H = 3600 * NS
DAY = 86400 * NS
NR = 2000          # random windows per trial
NSTK = 10000       # random stacks
RNG = np.random.default_rng(0)
CFG = {
    "SVI": dict(span=("2003-01-01", "2026-05-27"), segcol="seg"),
    "WA": dict(span=("2009-08-06", "2026-05-27"), segcol="era"),
    "WA_pre": dict(span=("2007-01-01", "2009-08-05"), segcol="era"),
}


def pt(ts):
    return pd.to_datetime(ts).tz_localize("UTC").tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z")


# ----------------------------------------------------------------------------------------------- earthquakes
def eq_table():
    E = pd.read_csv(f"{T}/eq_list_raw.csv", parse_dates=["t0_utc"])
    P = pd.read_csv(f"{T}/pgv_pgc.csv")
    E = E.merge(P, on="id", how="left")
    for reg in ("SVI", "WA"):
        s = E[f"{reg}_sig_est_kpa"].copy(); src = np.full(len(E), "estimate", dtype=object)
        ok = E.pgc_ok.fillna(False).astype(bool)
        if reg == "SVI":
            use = ok
            src[use] = "measured_PGC"
        else:
            use = ok & (E.WA_dkm >= 1000)
            src[use] = "PGC_proxy"
        s[use] = E.loc[use, "sig_pgc_kpa"]
        E[f"{reg}_sigma_kpa"] = s; E[f"{reg}_sigma_src"] = src
    E["clip_PGC"] = E.pgc_reason.eq("clip")
    E["t0_ns"] = E.t0_utc.astype("datetime64[ns]").astype("int64")
    return E


def trials(E, reg):
    """6-h clusters -> one trial; union window, representative = max sigma at region."""
    r = "WA" if reg.startswith("WA") else "SVI"
    E = E.sort_values("t0_ns").reset_index(drop=True)
    cl = np.zeros(len(E), int); c = 0
    for i in range(1, len(E)):
        if E.t0_ns[i] - E.t0_ns[i - 1] > 6 * H:
            c += 1
        cl[i] = c
    E = E.assign(cluster=cl)
    E["wP"] = E.t0_ns + (E[f"{r}_tP_s"] * NS).astype("int64")
    E["wL"] = E.t0_ns + (E[f"{r}_tL_s"] * NS).astype("int64")
    E["wR"] = E.t0_ns + (E[f"{r}_tR_s"] * NS).astype("int64")
    rows = []
    for k, g in E.groupby("cluster"):
        rep = g.loc[g[f"{r}_sigma_kpa"].idxmax()]
        rows.append(dict(trial=int(k), id=rep.id, n_members=len(g), members=";".join(g.id), t0_ns=int(rep.t0_ns),
                         t0_PT=pt(rep.t0_utc), mag=rep.mag, place=rep.place, ddeg=rep[f"{r}_ddeg"], dkm=rep[f"{r}_dkm"],
                         sigma_kpa=rep[f"{r}_sigma_kpa"], sigma_src=rep[f"{r}_sigma_src"],
                         sigma_est_kpa=rep[f"{r}_sig_est_kpa"], clip_PGC=bool(rep.clip_PGC),
                         inst_s=int(g.wP.min()), inst_e=int(g.wR.max() + H),
                         surf_s=int(g.wL.min()), surf_e=int(g.wR.max() + H // 2),
                         h6_s=int(g.wP.min()), h6_e=int(g.wP.min() + 6 * H)))
    Tr = pd.DataFrame(rows)
    Tr["pre_s"] = Tr.inst_s - (Tr.inst_e - Tr.inst_s); Tr["pre_e"] = Tr.inst_s
    Tr["tier"] = pd.cut(Tr.sigma_kpa, [-1, 0.1, 1, 1e9], labels=["<0.1", "0.1-1", ">=1"]).astype(str)
    Tr["regional"] = Tr.dkm < 1500
    return Tr


# ----------------------------------------------------------------------------------------------- region data
def day_ns(d):
    return pd.to_datetime(pd.Series(d)).astype("datetime64[ns]").astype("int64").to_numpy()


def load_region(reg):
    if reg == "SVI":
        Ev = pd.read_parquet(f"{T}/svi_events.parquet"); Ca = pd.read_parquet(f"{T}/svi_cand.parquet")
        O = pd.read_parquet(f"{T}/svi_obsdays.parquet").rename(columns={"seg": "grp"})
        ep = pd.read_csv(f"{T}/svi_episodes.csv"); Ev["grp"] = Ev.seg
    else:
        Ev = pd.read_parquet(f"{T}/wa_events.parquet"); Ca = pd.read_parquet(f"{T}/wa_cand.parquet")
        O = pd.read_parquet(f"{T}/wa_obsdays.parquet").rename(columns={"era": "grp"})
        ep = pd.read_csv(f"{T}/wa_episodes.csv"); Ev["grp"] = Ev.era
    t = lambda x: np.sort(x.time_utc.astype("datetime64[ns]").astype("int64").to_numpy())
    K = {"data": t(Ev[Ev.series == "data"]), "c1": t(Ev[Ev.series == "c1"]),
         "mz": t(Ev[(Ev.series == "data") & (Ev.n_z2 >= 2)]), "c1mz": t(Ev[(Ev.series == "c1") & (Ev.n_z2 >= 2)]),
         "cand_data": t(Ca[Ca.series == "data"]), "cand_c1": t(Ca[Ca.series == "c1"])}
    if reg != "SVI":
        K["rev"] = t(Ev[Ev.series == "rev"])
    # tremor epicentres within 100 km of the region reference
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    ref = (48.4837, -123.7485) if reg == "SVI" else (47.8432, -123.0282)
    dk = 111.195 * np.hypot(tr.lat - ref[0], (tr.lon - ref[1]) * np.cos(np.radians(ref[0])))
    K["trem100"] = np.sort(pd.to_datetime(tr.time[dk <= 100]).astype("datetime64[ns]").astype("int64").to_numpy())
    # day tables
    O = O[O.n_obs > 0]
    days = pd.DataFrame({"date": O.date.values, "grp": O.grp.values, "n_obs": O.n_obs.values})
    days["dns"] = day_ns(days.date)
    days = days.sort_values("dns").drop_duplicates("dns").reset_index(drop=True)
    ep["s"] = day_ns(ep.start); ep["e"] = day_ns(ep.end) + DAY  # end exclusive
    inep = np.zeros(len(days), bool); inpad = np.zeros(len(days), bool)
    for s, e in zip(ep.s, ep.e):
        inep |= (days.dns >= s) & (days.dns < e)
        inpad |= (days.dns >= s - 2 * DAY) & (days.dns < e + 2 * DAY)
    days["in_ep"] = inep; days["in_pad"] = inpad
    # daily counts for the delayed test
    for k in ("data", "c1", "mz", "c1mz"):
        dd = (K[k] // DAY) * DAY
        u, n = np.unique(dd, return_counts=True)
        days[f"n_{k}"] = pd.Series(n, index=u).reindex(days.dns).fillna(0).to_numpy()
    return K, days, ep


def count(arr, s, e):
    return np.searchsorted(arr, e, "left") - np.searchsorted(arr, s, "left")


# ----------------------------------------------------------------------------------------------- test 1
def inst(reg, Tr, K, days, allTr_ns):
    span = CFG[reg]["span"]
    s0, s1 = day_ns([span[0]])[0], day_ns([span[1]])[0] + DAY
    D = days.set_index("dns")
    near_eq = np.zeros(len(days), bool)
    for t in allTr_ns:
        near_eq |= np.abs(days.dns.to_numpy() - (t // DAY) * DAY) <= 3 * DAY
    days = days.assign(near_eq=near_eq)
    pool_ok = (~days.in_pad) & (~days.near_eq) & (days.dns >= s0) & (days.dns < s1)
    kinds = list(K)
    out, rnd, loc = [], {}, {}
    for _, r in Tr.iterrows():
        if not (s0 <= r.inst_s < s1):
            continue
        dlist = np.arange((r.inst_s // DAY) * DAY, (r.h6_e // DAY) * DAY + DAY, DAY)  # days the 6-h/inst window touch
        d_inst = np.arange((r.inst_s // DAY) * DAY, ((r.inst_e - 1) // DAY) * DAY + DAY, DAY)
        if not all(d in D.index for d in d_inst):
            status = "not_observable"
        elif any(D.loc[d, "in_ep"] for d in d_inst):
            status = "ets"
        elif any(D.loc[d, "in_pad"] for d in d_inst):
            status = "pad"
        else:
            status = "between"
        grp = D.loc[d_inst[0], "grp"] if d_inst[0] in D.index else None
        row = dict(trial=r.trial, status=status, grp=grp)
        for w in ("inst", "surf", "h6", "pre"):
            for k in kinds:
                row[f"{w}_{k}"] = int(count(K[k], r[f"{w}_s"], r[f"{w}_e"]))
        out.append(row)
        L = r.inst_e - r.inst_s
        tod = r.inst_s - (r.inst_s // DAY) * DAY
        if status == "between":
            ndays = int(np.ceil((tod + max(L, 6 * H)) / DAY))
            ok = pool_ok & (days.grp == grp)
            okd = set(days.dns[ok])
            cand = np.array([d for d in days.dns[ok] if all((d + j * DAY) in okd for j in range(ndays))])
            pick = RNG.choice(cand, NR, replace=len(cand) < NR)
            starts = pick + tod
            R = {}
            for w, Lw, off in (("inst", L, 0), ("surf", r.surf_e - r.surf_s, r.surf_s - r.inst_s),
                               ("h6", 6 * H, r.h6_s - r.inst_s)):
                for k in kinds:
                    R[f"{w}_{k}"] = count(K[k], starts + off, starts + off + Lw)
            R["n_pool"] = len(cand)
            rnd[r.trial] = R
        elif status == "ets":
            nt = int((24 * H) // L)
            ts = np.concatenate([r.inst_s - L * np.arange(1, nt + 1), r.inst_e + L * np.arange(0, nt)])
            loc[r.trial] = {k: count(K[k], ts, ts + L) for k in kinds}
    return pd.DataFrame(out), rnd, loc


def per_trial_stats(TT, rnd, loc):
    rows = []
    for _, r in TT.iterrows():
        src = rnd.get(r.trial) if r.status == "between" else (loc.get(r.trial) if r.status == "ets" else None)
        if src is None:
            rows.append({}); continue
        key = (lambda k: f"inst_{k}") if r.status == "between" else (lambda k: k)
        x = {}
        for k in ("data", "c1", "mz", "cand_data", "cand_c1", "trem100"):
            v = np.asarray(src[key(k)])
            x[f"mu_{k}"] = v.mean(); x[f"q99_{k}"] = np.quantile(v, 0.99); x[f"max_{k}"] = v.max()
            x[f"pemp_{k}"] = (1 + (v >= r[f"inst_{k}"]).sum()) / (1 + len(v))
        mu = max(x["mu_data"], 1.0 / NR)
        x["pois_p"] = stats.poisson.sf(r.inst_data - 1, mu)
        rows.append(x)
    return pd.concat([TT.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def stack_test(TT, rnd, loc, sel, w="inst", ets=False):
    S = TT[sel]
    if not len(S):
        return dict(n=0)
    src = loc if ets else rnd
    key = (lambda k: k) if ets else (lambda k: f"{w}_{k}")
    res = dict(n=len(S))
    obs = {k: S[f"{w}_{k}"].sum() for k in ("data", "c1", "mz", "cand_data", "cand_c1", "trem100")}
    obs["exc"] = obs["data"] - obs["c1"]
    picks = []
    for t in S.trial:
        a = src[t]; n = len(a[key("data")])
        picks.append(RNG.integers(0, n, NSTK))
    for k in ("data", "c1", "mz", "cand_data", "cand_c1", "trem100", "exc"):
        if k == "exc":
            M = np.sum([np.asarray(src[t][key("data")])[p] - np.asarray(src[t][key("c1")])[p] for t, p in zip(S.trial, picks)], 0)
        else:
            M = np.sum([np.asarray(src[t][key(k)])[p] for t, p in zip(S.trial, picks)], 0)
        res[f"obs_{k}"] = int(obs[k]); res[f"rnd_mean_{k}"] = float(M.mean())
        res[f"rnd_q99_{k}"] = float(np.quantile(M, 0.99)); res[f"p_{k}"] = float((1 + (M >= obs[k]).sum()) / (1 + NSTK))
        if k == "data":
            res["ub99_extra_per_trial"] = float((obs[k] - np.quantile(M, 0.01)) / len(S))
            res["rnd_q01_data"] = float(np.quantile(M, 0.01))
            res["ub99_rel_extra"] = float((obs[k] - np.quantile(M, 0.01)) / M.mean()) if M.mean() > 0 else np.nan
    res["PASS"] = bool(res["p_data"] < 0.01 and res["p_exc"] < 0.01 and res["p_mz"] < 0.01)
    return res


# ----------------------------------------------------------------------------------------------- test 2
def delayed(reg, Tr, days, allTr_ns, sel_trials, win=10, _cache={}):
    span = CFG[reg]["span"]; s0, s1 = day_ns([span[0]])[0], day_ns([span[1]])[0] + DAY
    dns = days.dns.to_numpy(); grp_arr = days.grp.to_numpy()
    btw = (~days.in_pad).to_numpy()
    idx = {int(d): i for i, d in enumerate(dns)}
    near_eq = np.zeros(len(days), bool)
    for t in allTr_ns:
        near_eq |= np.abs(dns - (t // DAY) * DAY) <= 3 * DAY
    nobs = days.n_obs.to_numpy().astype(float)
    C = {c: days[c].to_numpy().astype(float) for c in ("n_data", "n_c1", "n_mz", "n_c1mz")}
    defs = (("exc", "n_data", "n_c1"), ("data", "n_data", None), ("c1", "n_c1", None), ("mz", "n_mz", "n_c1mz"))

    def val(d0):
        A = [idx[x] for x in (d0 + j * DAY for j in range(1, win + 1)) if x in idx and btw[idx[x]]]
        B = [idx[x] for x in (d0 - j * DAY for j in range(1, 11)) if x in idx and btw[idx[x]]]
        if len(A) < min(5, win) or len(B) < 5:
            return None
        o = {}
        na, nb = nobs[A].sum(), nobs[B].sum()
        for nm, c0, c1 in defs:
            ra = (C[c0][A].sum() - (C[c1][A].sum() if c1 else 0)) / na
            rb = (C[c0][B].sum() - (C[c1][B].sum() if c1 else 0)) / nb
            o[f"after_{nm}"] = ra; o[f"diff_{nm}"] = ra - rb
        return o

    rows, pools = [], {}
    for _, r in Tr[sel_trials].iterrows():
        d0 = int((r.inst_s // DAY) * DAY)
        if not (s0 <= d0 < s1) or d0 not in idx or not btw[idx[d0]]:
            continue
        v = val(d0)
        if v is None:
            continue
        grp = grp_arr[idx[d0]]
        ck = (reg, grp, win)
        if ck not in _cache:
            cand = dns[(grp_arr == grp) & btw & ~near_eq & (dns >= s0) & (dns < s1)]
            pv = [x for x in (val(int(d)) for d in cand) if x is not None]
            _cache[ck] = pd.DataFrame(pv)
        pools[r.trial] = _cache[ck]
        rows.append(dict(trial=r.trial, grp=grp, sigma_kpa=r.sigma_kpa, **v))
    R = pd.DataFrame(rows)
    res = dict(n=len(R))
    if not len(R):
        return R, res, pools
    picks = {t: RNG.integers(0, len(pools[t]), NSTK) for t in R.trial}
    for c in [c for c in R.columns if c.startswith(("after_", "diff_"))]:
        M = np.mean([pools[t][c].to_numpy()[picks[t]] for t in R.trial], 0)
        o = R[c].mean()
        res[f"obs_{c}"] = float(o); res[f"rnd_mean_{c}"] = float(M.mean()); res[f"p_{c}"] = float((1 + (M >= o).sum()) / (1 + NSTK))
        if c == "after_exc":
            res["ub99_after_exc_minus_rnd"] = float(o - np.quantile(M, 0.01))
            res["rnd_mean_rate"] = float(M.mean())
    res["PASS"] = bool(res["p_after_exc"] < 0.01 and res["p_diff_exc"] < 0.01 and res["p_after_mz"] < 0.01
                       and res["p_diff_mz"] < 0.01 and res["p_after_c1"] >= 0.01)
    return R, res, pools


# ----------------------------------------------------------------------------------------------- test 3
def dose_inst(S):
    S = S.dropna(subset=["mu_data"])
    N = S.inst_data.to_numpy(); mu = np.maximum(S.mu_data.to_numpy(), 1.0 / NR); sg = S.sigma_kpa.to_numpy()

    def nll(b):
        lam = mu * (1 + b * sg)
        return -(stats.poisson.logpmf(N, lam)).sum()
    bs = np.concatenate([[0], np.logspace(-5, 3, 4000)])
    L = np.array([nll(b) for b in bs])
    bhat = bs[L.argmin()]
    ub = bs[(L <= L.min() + 5.412 / 2) ].max()      # one-sided 99 % profile bound
    lb = bs[(L <= L.min() + 5.412 / 2)].min()
    resp = (N + 0.5) / (mu + 0.5)
    rho = stats.spearmanr(sg, resp).statistic if len(S) > 2 else np.nan
    perm = np.array([stats.spearmanr(RNG.permutation(sg), resp).statistic for _ in range(2000)]) if len(S) > 2 else np.array([0])
    p = (1 + (perm >= rho).sum()) / (1 + len(perm))
    return dict(n=len(S), beta_hat_per_kpa=float(bhat), beta_ub99_per_kpa=float(ub), beta_lb99=float(lb),
                rho=float(rho), p_rho=float(p), scales=bool(rho > 0 and p < 0.01 and lb > 0))


def ets_onsets(reg, Tr, ep, sel):
    span = CFG[reg]["span"]; s0, s1 = day_ns([span[0]])[0], day_ns([span[1]])[0] + DAY
    on = ep.s[ep.major.astype(bool) & (ep.s >= s0) & (ep.s < s1)].to_numpy()
    tt = Tr.loc[sel, "t0_ns"].to_numpy(); tt = tt[(tt >= s0) & (tt < s1)]

    def K(tt, w):
        return int(sum(((tt <= o + DAY) & (tt >= o - w * DAY)).any() for o in on))
    res = dict(n_onsets=len(on), n_eq=len(tt))
    for w in (10, 30):
        k = K(tt, w); span_ns = s1 - s0
        sh = np.array([K(s0 + (tt - s0 + RNG.integers(0, span_ns)) % span_ns, w) for _ in range(NSTK)])
        un = np.array([K(RNG.integers(s0, s1, len(tt)), w) for _ in range(NSTK)])
        res.update({f"K{w}": k, f"K{w}_shift_mean": float(sh.mean()), f"p{w}_shift": float((1 + (sh >= k).sum()) / (1 + NSTK)),
                    f"K{w}_unif_mean": float(un.mean()), f"p{w}_unif": float((1 + (un >= k).sum()) / (1 + NSTK))})
    res["supported"] = bool(res["p10_shift"] < 0.01 and res["p10_unif"] < 0.01)
    hits = []
    for o in on:
        m = Tr.loc[sel & (Tr.t0_ns <= o + DAY) & (Tr.t0_ns >= o - 10 * DAY)]
        for _, r in m.iterrows():
            hits.append(dict(region=reg, onset=str(pd.Timestamp(o).date()), eq=r.place, eq_PT=r.t0_PT, mag=r.mag,
                             sigma_kpa=r.sigma_kpa, days_before=(o - r.t0_ns) / DAY))
    return res, hits


# ----------------------------------------------------------------------------------------------- main
def main():
    E = eq_table()
    summ, stacks, dl_rows, dose_rows, ets_rows, ets_hits = {}, [], [], [], [], []
    for reg in ("SVI", "WA", "WA_pre"):
        r = "WA" if reg.startswith("WA") else "SVI"
        Tr = trials(E, reg)
        K, days, ep = load_region(r)
        if reg == "WA_pre":
            ep = ep[ep.kind == "wecc"]
            K2, days, _ = load_region(r)
            ep2 = ep.copy()
            inep = np.zeros(len(days), bool); inpad = np.zeros(len(days), bool)
            for s, e in zip(ep2.s, ep2.e):
                inep |= (days.dns >= s) & (days.dns < e); inpad |= (days.dns >= s - 2 * DAY) & (days.dns < e + 2 * DAY)
            days["in_ep"] = inep; days["in_pad"] = inpad
        allns = Tr.t0_ns.to_numpy()
        TT, rnd, loc = inst(reg, Tr, K, days, allns)
        TT = Tr.merge(TT, on="trial")
        TT = per_trial_stats(TT, rnd, loc)
        TT["region"] = reg
        TT.to_csv(f"{T}/inst_trials_{reg}.csv", index=False)
        btw = TT.status.eq("between"); ets = TT.status.eq("ets"); strong = TT.sigma_kpa >= 1
        sets = {"strong": btw & strong, "all": btw, "0.1-1": btw & TT.tier.eq("0.1-1"), "<0.1": btw & TT.tier.eq("<0.1"),
                ">=10": btw & (TT.sigma_kpa >= 10), "regional_strong": btw & strong & TT.regional,
                "tele_strong": btw & strong & ~TT.regional}
        for g in sorted(TT.grp.dropna().unique()):
            sets[f"strong_{g}"] = btw & strong & TT.grp.eq(g)
            sets[f"all_{g}"] = btw & TT.grp.eq(g)
        for nm, sel in sets.items():
            for w in (("inst", "surf", "h6") if nm in ("strong", "all") else ("inst",)):
                stacks.append(dict(region=reg, set=nm, window=w, **stack_test(TT, rnd, loc, sel, w)))
        for nm, sel in {"ets_strong": ets & strong, "ets_all": ets}.items():
            stacks.append(dict(region=reg, set=nm, window="inst_vs_local", **stack_test(TT, rnd, loc, sel, ets=True)))
        # named individual earthquakes (rule §3)
        nT = int((btw | ets).sum())
        TT["named_abcd"] = ((TT.inst_data > TT.max_data) & (TT.pois_p < 0.01 / max(nT, 1)) & (TT.inst_mz >= 2)
                            & (TT.inst_c1 <= TT.q99_c1) & (TT.inst_cand_c1 <= TT.q99_cand_c1)).fillna(False)
        TT.to_csv(f"{T}/inst_trials_{reg}.csv", index=False)
        # test 2
        for nm, sel in (("strong", Tr.sigma_kpa >= 1), ("all", Tr.sigma_kpa > -1)):
            for win in (10, 3):
                R, res, _ = delayed(reg, Tr, days, allns, sel, win)
                R.assign(region=reg, set=nm, win=win).to_csv(f"{T}/delayed_trials_{reg}_{nm}_{win}.csv", index=False)
                dl_rows.append(dict(region=reg, set=nm, days=f"+1..+{win}", **res))
                if nm == "all" and win == 10 and len(R) > 2:
                    rr = stats.spearmanr(R.sigma_kpa, R.after_exc).statistic
                    pr = np.array([stats.spearmanr(RNG.permutation(R.sigma_kpa.values), R.after_exc).statistic for _ in range(2000)])
                    X = R.sigma_kpa.to_numpy(); Y = R.after_exc.to_numpy()
                    bs = []
                    for _ in range(2000):
                        i = RNG.integers(0, len(R), len(R))
                        bs.append(np.polyfit(X[i], Y[i], 1)[0] if np.ptp(X[i]) > 0 else np.nan)
                    base = res.get("rnd_mean_rate", np.nan)
                    dose_rows.append(dict(region=reg, test="delayed_after_exc", n=len(R), rho=rr,
                                          p_rho=(1 + (pr >= rr).sum()) / 2001,
                                          slope_per_kpa=np.polyfit(X, Y, 1)[0], slope_ub99=np.nanquantile(bs, 0.99),
                                          base_rate=base, ub99_pct_per_kpa=100 * np.nanquantile(bs, 0.99) / base if base > 0 else np.nan))
        # test 3
        for nm, sel in (("all", btw), ("ets_all", ets)):
            if nm == "ets_all":
                continue
            d = dose_inst(TT[sel])
            dose_rows.append(dict(region=reg, test=f"inst_{nm}", **d))
        for nm, sel in (("strong", Tr.sigma_kpa >= 1), ("all", Tr.sigma_kpa > -1), (">=10", Tr.sigma_kpa >= 10)):
            res, hits = ets_onsets(reg, Tr, ep, sel)
            ets_rows.append(dict(region=reg, eq_set=nm, **res)); ets_hits += [dict(h, eq_set=nm) for h in hits]
        summ[reg] = dict(n_trials=len(Tr), status=TT.status.value_counts().to_dict(),
                         n_strong_between=int((btw & strong).sum()), n_strong_ets=int((ets & strong).sum()),
                         named=TT.loc[TT.named_abcd, ["t0_PT", "place", "mag", "sigma_kpa", "status"]].to_dict("records"))
        print(reg, json.dumps(summ[reg], default=str)[:600], flush=True)
    E.to_csv(f"{T}/eq_list.csv", index=False)
    pd.DataFrame(stacks).to_csv(f"{T}/inst_stacks.csv", index=False)
    pd.DataFrame(dl_rows).to_csv(f"{T}/delayed_stacks.csv", index=False)
    pd.DataFrame(dose_rows).to_csv(f"{T}/dose.csv", index=False)
    pd.DataFrame(ets_rows).to_csv(f"{T}/ets_onsets.csv", index=False)
    pd.DataFrame(ets_hits).to_csv(f"{T}/ets_onset_hits.csv", index=False)
    json.dump(summ, open(f"{T}/g_summary.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
