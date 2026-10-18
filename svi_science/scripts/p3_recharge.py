#!/usr/bin/env python
"""p3_recharge (2026-10-07) -- P3 step 4: recharge curve of ALL deep activity vs time since the last major ETS ended,
recovery time constant per era (bootstrap over cycles), time-shift surrogates, era comparisons, and the Bostock-era
deep-share substitute. Rules frozen in svi_science/results/p3_bostock/PLAN_P3.md.
    python svi_science/scripts/p3_recharge.py
"""
import os, json, time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/p3_bostock"
S0 = f"{ROOT}/svi_science/results/s0"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
B_EP = f"{ROOT}/svi_science/results/b_ets_record/episodes_b.csv"
P0, P1 = "2014-01-01", "2026-05-27"
TMAX = 400
NB = 1000
RNG = np.random.default_rng(20261008)
BINS = [1, 8, 15, 31, 61, 91, 151, 241, 401]          # t in days since end (end day = 0 is inside the episode)
BLAB = ["0-7", "7-14", "14-30", "30-60", "60-90", "90-150", "150-240", "240-400"]
DEEP_GROUPS = [1, 4, 10, 12]
UPDIP_GROUPS = [3, 6, 8, 11]
LOG = open(f"{OUT}/p3_recharge.log", "w")
RES = {}


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


# ------------------------------------------------------------------------------------------- cycle tables
def cycle_table(days, y, n, ref, shift=0):
    """days: datetime64[D] array (all days of the period, consecutive), y counts, n exposure per day.
    ref: DataFrame start/end (reference episodes). Returns a DataFrame of usable days with cycle id and t (days since
    previous ref end, minus shift); days inside ref episodes / before the first end / beyond TMAX dropped."""
    rs = ref.start.values.astype("datetime64[D]"); re_ = ref.end.values.astype("datetime64[D]")
    inside = np.zeros(days.size, bool)
    for s, e in zip(rs, re_):
        inside |= (days >= s) & (days <= e)
    j = np.searchsorted(re_, days, side="left") - 1
    ok = (j >= 0) & ~inside
    t = np.where(ok, (days - re_[np.maximum(j, 0)]).astype(int), -1) - shift
    ok &= (t >= 1) & (t <= TMAX)
    return pd.DataFrame(dict(date=days[ok], cyc=j[ok], t=t[ok], y=y[ok], n=n[ok]))


def agg(C):
    a = C.groupby("t").agg(Y=("y", "sum"), N=("n", "sum")).reindex(range(1, TMAX + 1), fill_value=0)
    return a.Y.to_numpy(float), a.N.to_numpy(float)


def nll(p, Y, N, tt):
    r, a, tau = np.exp(p[0]), p[1], np.exp(p[2])
    mu = r * (1 - a * np.exp(-tt / tau))
    mu = np.maximum(mu, 1e-12)
    return float(np.sum(N * mu - Y * np.log(mu)))


def fit(Y, N):
    tt = np.arange(1, TMAX + 1, dtype=float)
    m = N > 0
    Y, N, tt = Y[m], N[m], tt[m]
    if Y.sum() == 0:
        return dict(r_inf=0.0, a=np.nan, tau=np.nan)
    r0 = Y.sum() / N.sum()
    best = None
    for tau0 in (3, 10, 30, 100, 300):
        for a0 in (0.5, -0.5):
            res = minimize(nll, [np.log(r0), a0, np.log(tau0)], args=(Y, N, tt), method="L-BFGS-B",
                           bounds=[(np.log(r0) - 5, np.log(r0) + 5), (-5, 1), (0, np.log(400))])
            if best is None or res.fun < best.fun:
                best = res
    return dict(r_inf=float(np.exp(best.x[0])), a=float(best.x[1]), tau=float(np.exp(best.x[2])))


def d30(Y, N):
    tt = np.arange(1, TMAX + 1)
    a, b = tt < 31, tt >= 31
    if N[a].sum() == 0 or N[b].sum() == 0 or Y[b].sum() == 0:
        return np.nan
    return float(1 - (Y[a].sum() / N[a].sum()) / (Y[b].sum() / N[b].sum()))


def curve(Y, N, r_inf):
    rows = []
    for lo, hi, lab in zip(BINS[:-1], BINS[1:], BLAB):
        s = slice(lo - 1, hi - 1)
        rows.append(dict(bin=lab, events=int(Y[s].sum()), exposure=float(N[s].sum()),
                         rate=float(Y[s].sum() / N[s].sum()) if N[s].sum() else np.nan))
    df = pd.DataFrame(rows); df["norm_rate"] = df.rate / r_inf if r_inf > 0 else np.nan
    return df


def analyse(name, C, surrogates=True):
    """C: cycle table. Fit + bootstrap over cycles + time-shift surrogates."""
    Y, N = agg(C)
    f = fit(Y, N); D = d30(Y, N)
    cyc = sorted(C.cyc.unique())
    per = {c: agg(C[C.cyc == c]) for c in cyc}
    Ys = np.array([per[c][0] for c in cyc]); Ns = np.array([per[c][1] for c in cyc])
    bt, bd, ba = [], [], []
    for _ in range(NB):
        k = RNG.integers(0, len(cyc), len(cyc))
        Yb, Nb = Ys[k].sum(0), Ns[k].sum(0)
        fb = fit(Yb, Nb); bt.append(fb["tau"]); ba.append(fb["a"]); bd.append(d30(Yb, Nb))
    bt, ba, bd = np.array(bt), np.array(ba), np.array(bd)
    out = dict(name=name, n_events=int(Y.sum()), exposure=float(N.sum()), n_cycles=len(cyc), **f, D30=D,
               tau_ci=[float(np.nanpercentile(bt, 2.5)), float(np.nanpercentile(bt, 97.5))],
               a_ci=[float(np.nanpercentile(ba, 2.5)), float(np.nanpercentile(ba, 97.5))],
               D30_ci=[float(np.nanpercentile(bd, 2.5)), float(np.nanpercentile(bd, 97.5))],
               rate_0_30_per_1000=float(1000 * Y[:30].sum() / max(N[:30].sum(), 1)),
               rate_30plus_per_1000=float(1000 * Y[30:].sum() / max(N[30:].sum(), 1)))
    if surrogates:
        sd = []
        groups = {c: g.sort_values("t") for c, g in C.groupby("cyc")}
        for _ in range(NB):
            Ysum = np.zeros(TMAX); Nsum = np.zeros(TMAX)
            for c, g in groups.items():
                L = len(g); k = RNG.integers(0, L)
                ti = g.t.to_numpy() - 1
                np.add.at(Ysum, ti, np.roll(g.y.to_numpy(), k)); np.add.at(Nsum, ti, np.roll(g.n.to_numpy(), k))
            sd.append(d30(Ysum, Nsum))
        sd = np.array(sd)
        out.update(surrogate_D30_median=float(np.nanmedian(sd)), surrogate_D30_p99=float(np.nanpercentile(sd, 99)),
                   p_surrogate=float((np.sum(sd >= D) + 1) / (NB + 1)))
        out["recovery_present"] = bool(D > out["surrogate_D30_p99"] and out["a_ci"][0] > 0)
    out["_boot_D30"] = bd; out["_boot_tau"] = bt
    log(f"{name}: " + json.dumps({k: v for k, v in out.items() if not k.startswith('_')}, default=str))
    cv = curve(Y, N, f["r_inf"]); cv.insert(0, "set", name)
    return out, cv


def compare(o1, o2, name):
    ci1, ci2 = o1["tau_ci"], o2["tau_ci"]
    overlap = not (ci1[1] < ci2[0] or ci2[1] < ci1[0])
    diff = o1["_boot_D30"] - o2["_boot_D30"]
    diff = diff[np.isfinite(diff)]
    p2 = float(2 * min(np.mean(diff <= 0), np.mean(diff >= 0)))
    lo, hi = np.percentile(diff, [2.5, 97.5])
    different = (not overlap) and p2 < 0.01
    r = dict(comparison=name, tau_1=o1["tau"], tau_1_ci=ci1, tau_2=o2["tau"], tau_2_ci=ci2, tau_ci_overlap=overlap,
             D30_1=o1["D30"], D30_2=o2["D30"], D30_diff=o1["D30"] - o2["D30"], D30_diff_ci=[float(lo), float(hi)],
             p_two_sided=p2, verdict="different" if different else "no detectable difference",
             detectable_D30_difference=float((hi - lo) / 2))
    log(f"COMPARE {json.dumps(r, default=str)}")
    return r


# ------------------------------------------------------------------------------------------- main
def main():
    F = pd.read_csv(f"{S0}/families.csv").set_index("fid", drop=False)
    deep = F.fid[F.group.isin(DEEP_GROUPS)].to_numpy()
    sets = {"deep groups 1/4/10/12": deep, "families >= 40 km": F.fid[F.dep >= 40].to_numpy(),
            "group 10": F.fid[F.group == 10].to_numpy(), "up-dip groups 3/6/8/11": F.fid[F.group.isin(UPDIP_GROUPS)].to_numpy()}
    log({k: len(v) for k, v in sets.items()})
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    maj39 = E[E.major].sort_values("start").reset_index(drop=True)
    Bp = pd.read_csv(B_EP, parse_dates=["start", "end"])
    Bs = Bp[Bp.region == "SVI"]
    big = Bs[(Bs.big == True) & (Bs.start >= "2013-01-01")].sort_values("start").reset_index(drop=True)  # incl. 2013-09 as first reset
    refs = {"39 majors": maj39, "big list": big}
    log(f"big list: {big.episode.tolist()}")
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8"])
    T = T[(T.date >= P0) & (T.date <= P1) & (T.P_same >= 0.75) & T.keep_xfam8]
    D = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable"])
    D = D[(D.run == "prod") & D.observable & (D.date >= P0) & (D.date <= P1)]
    days = np.arange(np.datetime64(P0), np.datetime64(P1) + 1)
    ds = pd.Series(days.astype(str))

    def daily(fids, series="data"):
        y = T[(T.series == series) & T.fam.isin(fids)].groupby("date").size().reindex(ds, fill_value=0).to_numpy()
        n = D[D.fid.isin(fids)].groupby("date").size().reindex(ds, fill_value=0).to_numpy()
        return y, n

    outs, curves, comps = [], [], []
    # --- primary + columns (ours, LFE)
    for rname, ref in refs.items():
        for sname, fids in sets.items():
            y, n = daily(fids)
            C = cycle_table(days, y, n, ref)
            prim = (rname == "39 majors" and sname == "deep groups 1/4/10/12")
            o, cv = analyse(f"LFE {sname} | ref {rname} | D1", C, surrogates=(sname != "group 10"))
            outs.append(o); curves.append(cv)
            if prim:
                RES["primary"] = {k: v for k, v in o.items() if not k.startswith("_")}
                Cprim = C
    # scrambled control, primary sets
    y, n = daily(deep, "c1")
    o, cv = analyse("LFE scrambled control deep groups | ref 39 majors | D1", cycle_table(days, y, n, maj39), surrogates=False)
    outs.append(o); curves.append(cv)
    # D2 (end + 7 d): t' = t - 7, t >= 8 kept
    y, n = daily(deep)
    C = cycle_table(days, y, n, maj39)
    C2 = C[C.t >= 8].copy(); C2["t"] = C2.t - 7
    o, cv = analyse("LFE deep groups | ref 39 majors | D2 end+7", C2); outs.append(o); curves.append(cv)
    # D3 family's own end: family-day table
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    Xe = X[X.n30 > 0]
    rs = maj39.start.values.astype("datetime64[D]"); re_ = maj39.end.values.astype("datetime64[D]")
    rows = []
    obs = D[D.fid.isin(deep)]
    cnt = T[(T.series == "data") & T.fam.isin(deep)].groupby(["fam", "date"]).size()
    for f in deep:
        dfx = np.sort(pd.to_datetime(Xe.date[Xe.fid == f]).values.astype("datetime64[D]"))
        le = re_.copy()
        for j, (s, e) in enumerate(zip(rs, re_)):
            k = dfx[(dfx >= s) & (dfx <= e)]
            if k.size:
                le[j] = k[-1]
        od = np.sort(pd.to_datetime(obs.date[obs.fid == f]).values.astype("datetime64[D]"))
        inside = np.zeros(od.size, bool)
        for s, e in zip(rs, re_):
            inside |= (od >= s) & (od <= e)
        j = np.searchsorted(re_, od, side="left") - 1
        ok = (j >= 0) & ~inside
        t = (od - le[np.maximum(j, 0)]).astype(int)
        ok &= (t >= 1) & (t <= TMAX)
        yy = np.array([cnt.get((f, str(d)), 0) for d in od[ok]])
        rows.append(pd.DataFrame(dict(date=od[ok], cyc=j[ok], t=t[ok], y=yy, n=1)))
    C3 = pd.concat(rows, ignore_index=True)
    o, cv = analyse("LFE deep groups | ref 39 majors | D3 family end", C3, surrogates=False); outs.append(o); curves.append(cv)
    # --- era comparison (E-a): halves of ours
    half = {}
    for lab, (a, b) in (("2014-2019", ("2014-01-01", "2019-12-31")), ("2020-2026", ("2020-01-01", P1))):
        Ch = Cprim[(Cprim.date >= np.datetime64(a)) & (Cprim.date <= np.datetime64(b))]
        o, cv = analyse(f"LFE deep groups | ref 39 majors | D1 | {lab}", Ch); outs.append(o); curves.append(cv); half[lab] = o
    comps.append(compare(half["2014-2019"], half["2020-2026"], "LFE deep, 2014-2019 vs 2020-2026"))
    # --- era comparison (E-b): PNSN deep tremor 2010-2013 vs 2014-2026
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["time"])
    tr = tr[tr.lat.between(47.8, 49.5) & tr.lon.between(-125.0, -122.5)].copy()
    dmin = np.full(len(tr), np.inf)
    for f in deep:
        dmin = np.minimum(dmin, hav_km(tr.lat.to_numpy(), tr.lon.to_numpy(), F.lat[f], F.lon[f]))
    tr = tr[dmin <= 15]
    tr["date"] = tr.time.astype("datetime64[ns]").dt.strftime("%Y-%m-%d")
    T0 = "2010-01-22"
    tdays = np.arange(np.datetime64(T0), np.datetime64(P1) + 1)
    yt = tr.groupby("date").size().reindex(pd.Series(tdays.astype(str)), fill_value=0).to_numpy()
    nt = np.ones(tdays.size)
    ref_t = pd.concat([Bs[(Bs.major == True) & (Bs.start < "2014-01-01")][["start", "end"]], maj39[["start", "end"]]]).sort_values("start").reset_index(drop=True)
    Ct = cycle_table(tdays, yt, nt, ref_t)
    era_t = {}
    for lab, (a, b) in (("2010-2013", (T0, "2013-12-31")), ("2014-2026", ("2014-01-01", P1))):
        Ce = Ct[(Ct.date >= np.datetime64(a)) & (Ct.date <= np.datetime64(b))]
        o, cv = analyse(f"PNSN deep tremor | PNSN majors | D1 | {lab}", Ce); outs.append(o); curves.append(cv); era_t[lab] = o
        Ce2 = Ce[Ce.t >= 8].copy(); Ce2["t"] = Ce2.t - 7
        o2, cv2 = analyse(f"PNSN deep tremor | PNSN majors | D2 end+7 | {lab}", Ce2); outs.append(o2); curves.append(cv2)
        era_t[lab + " D2"] = o2
    comps.append(compare(era_t["2010-2013"], era_t["2014-2026"], "PNSN deep tremor D1, 2010-2013 vs 2014-2026"))
    comps.append(compare(era_t["2010-2013 D2"], era_t["2014-2026 D2"], "PNSN deep tremor D2, 2010-2013 vs 2014-2026"))
    RES["n_deep_tremor"] = int(len(tr))

    # --- Bostock-era substitute: deep share per major episode
    M = pd.read_parquet(f"{OUT}/bostock_events_mapped.parquet")
    Eb = pd.read_csv(f"{OUT}/episodes_bostock_final.csv")
    rows = []
    prev_end = None
    for r in Eb[Eb.major].itertuples():
        m = M[(M.date >= r.start) & (M.date <= r.end)]
        rows.append(dict(era="Bostock 2003-2013", episode=r.episode, start=r.start, end=r.end, n=len(m),
                         n_deep=int(m.fid.isin(deep).sum()), interval_d=(pd.Timestamp(r.start) - pd.Timestamp(prev_end)).days if prev_end else np.nan,
                         year=pd.Timestamp(r.start).year + pd.Timestamp(r.start).dayofyear / 365.25))
        prev_end = r.end
    Td = T[T.series == "data"]
    for lab, ref in (("ours 2014-2026 big", big), ("ours 2014-2026 39 majors", maj39)):
        prev_end = "2013-10-13"
        for r in ref.itertuples():
            st, en = str(r.start.date()), str(r.end.date())
            if st < P0:
                continue
            m = Td[(Td.date >= st) & (Td.date <= en)]
            rows.append(dict(era=lab, episode=r.episode, start=st, end=en, n=len(m), n_deep=int(m.fam.isin(deep).sum()),
                             interval_d=(pd.Timestamp(st) - pd.Timestamp(prev_end)).days,
                             year=pd.Timestamp(st).year + pd.Timestamp(st).dayofyear / 365.25))
            prev_end = en
    S = pd.DataFrame(rows); S["deep_share"] = S.n_deep / S.n
    S.to_csv(f"{OUT}/deep_share_by_episode.csv", index=False)
    ds_res = {}
    for era, g in S.groupby("era"):
        g = g[g.n >= 30]
        res = dict(n_episodes=int(len(g)), mean=float(g.deep_share.mean()), cv=float(g.deep_share.std() / g.deep_share.mean()))
        for x in ("year", "interval_d"):
            h = g.dropna(subset=[x])
            rho = spearmanr(h[x], h.deep_share).statistic
            perm = np.array([spearmanr(h[x], RNG.permutation(h.deep_share.to_numpy())).statistic for _ in range(10000)])
            res[f"rho_{x}"] = float(rho); res[f"p_{x}"] = float((np.sum(np.abs(perm) >= abs(rho)) + 1) / 10001)
        ds_res[era] = res
        log(f"deep share {era}: {res}")
    RES["deep_share"] = ds_res

    O = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in o.items() if not k.startswith("_")} for o in outs])
    O.to_csv(f"{OUT}/recharge_fits.csv", index=False)
    pd.concat(curves).to_csv(f"{OUT}/recharge_curves.csv", index=False)
    pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in c.items()} for c in comps]).to_csv(f"{OUT}/recharge_era_compare.csv", index=False)
    RES["comparisons"] = comps
    json.dump(RES, open(f"{OUT}/recharge_summary.json", "w"), indent=1, default=str)
    log("done")


if __name__ == "__main__":
    main()
