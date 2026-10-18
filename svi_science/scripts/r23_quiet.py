#!/usr/bin/env python
"""r23_quiet (2026-10-11) -- R23 study 3: deep quiet after big ETS, 2003-2026 (finding 4; S7 F30 + P3 recharge).
Rules: svi_science/results/r23/PLAN.md §4. Functions cycle_table/agg/nll/fit/d30/curve/analyse copied from
svi_science/scripts/p3_recharge.py (importing it would truncate its log); only change = surrogate shifts within
(cycle, network segment) as pre-registered.     python svi_science/scripts/r23_quiet.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C

OUT = f"{C.TB}/study3"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/study3{'_union' if '--union' in sys.argv else ''}.log")
RNG = np.random.default_rng(20261011)
TMAX = 400; NB = 1000; NSH = 1000
BINS = [1, 8, 15, 31, 61, 91, 151, 241, 401]
BLAB = ["0-7", "7-14", "14-30", "30-60", "60-90", "90-150", "150-240", "240-400"]
SEGB = [("2003-01-01", "CN4f"), ("2003-02-01", "ERA"), ("2006-07-01", "B2f"), ("2007-11-13", "B3f"), ("2009-08-27", "S15f"), ("2014-01-01", "FWD")]


def seg_idx(days):
    s = np.zeros(days.size, int)
    for k, (d, _) in enumerate(SEGB):
        s[days >= np.datetime64(d)] = k
    return s


def cycle_table(days, y, n, ref, shift=0):
    rs = ref.start.values.astype("datetime64[D]"); re_ = ref.end.values.astype("datetime64[D]")
    inside = np.zeros(days.size, bool)
    for s, e in zip(rs, re_):
        inside |= (days >= s) & (days <= e)
    j = np.searchsorted(re_, days, side="left") - 1
    ok = (j >= 0) & ~inside
    t = np.where(ok, (days - re_[np.maximum(j, 0)]).astype(int), -1) - shift
    ok &= (t >= 1) & (t <= TMAX)
    return pd.DataFrame(dict(date=days[ok], cyc=j[ok], seg=seg_idx(days)[ok], t=t[ok], y=y[ok], n=n[ok]))


def agg(Cc):
    a = Cc.groupby("t").agg(Y=("y", "sum"), N=("n", "sum")).reindex(range(1, TMAX + 1), fill_value=0)
    return a.Y.to_numpy(float), a.N.to_numpy(float)


def nll(p, Y, N, tt):
    r, a, tau = np.exp(p[0]), p[1], np.exp(p[2])
    mu = np.maximum(r * (1 - a * np.exp(-tt / tau)), 1e-12)
    return float(np.sum(N * mu - Y * np.log(mu)))


def fit(Y, N):
    tt = np.arange(1, TMAX + 1, dtype=float); m = N > 0
    Y, N, tt = Y[m], N[m], tt[m]
    if Y.sum() == 0:
        return dict(r_inf=0.0, a=np.nan, tau=np.nan)
    r0 = Y.sum() / N.sum(); best = None
    for tau0 in (3, 10, 30, 100, 300):
        for a0 in (0.5, -0.5):
            res = minimize(nll, [np.log(r0), a0, np.log(tau0)], args=(Y, N, tt), method="L-BFGS-B",
                           bounds=[(np.log(r0) - 5, np.log(r0) + 5), (-5, 1), (0, np.log(400))])
            if best is None or res.fun < best.fun:
                best = res
    return dict(r_inf=float(np.exp(best.x[0])), a=float(best.x[1]), tau=float(np.exp(best.x[2])))


def d30(Y, N):
    tt = np.arange(1, TMAX + 1); a, b = tt < 31, tt >= 31
    if N[a].sum() == 0 or N[b].sum() == 0 or Y[b].sum() == 0:
        return np.nan
    return float(1 - (Y[a].sum() / N[a].sum()) / (Y[b].sum() / N[b].sum()))


def curve(Y, N, r_inf):
    rows = []
    for lo, hi, lab in zip(BINS[:-1], BINS[1:], BLAB):
        s = slice(lo - 1, hi - 1)
        rows.append(dict(bin=lab, events=int(Y[s].sum()), exposure=float(N[s].sum()), rate=float(Y[s].sum() / N[s].sum()) if N[s].sum() else np.nan))
    df = pd.DataFrame(rows); df["norm_rate"] = df.rate / r_inf if r_inf > 0 else np.nan
    return df


def analyse(name, Cc, surrogates=True):
    Y, N = agg(Cc)
    f = fit(Y, N); D = d30(Y, N)
    cyc = sorted(Cc.cyc.unique())
    per = {c: agg(Cc[Cc.cyc == c]) for c in cyc}
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
               rate_30plus_per_1000=float(1000 * Y[30:].sum() / max(N[30:].sum(), 1)),
               rate_30_60_per_1000=float(1000 * Y[30:60].sum() / max(N[30:60].sum(), 1)))
    if surrogates:
        sd = []
        groups = {c: g.sort_values("t") for c, g in Cc.groupby(["cyc", "seg"])}      # pre-registered: within (cycle, segment)
        for _ in range(NB):
            Ysum = np.zeros(TMAX); Nsum = np.zeros(TMAX)
            for c, g in groups.items():
                L = len(g); k = RNG.integers(0, L); ti = g.t.to_numpy() - 1
                np.add.at(Ysum, ti, np.roll(g.y.to_numpy(), k)); np.add.at(Nsum, ti, np.roll(g.n.to_numpy(), k))
            sd.append(d30(Ysum, Nsum))
        sd = np.array(sd)
        out.update(surrogate_D30_median=float(np.nanmedian(sd)), surrogate_D30_p99=float(np.nanpercentile(sd, 99)),
                   p_surrogate=float((np.sum(sd >= D) + 1) / (NB + 1)))
        out["recovery_present"] = bool(D > out["surrogate_D30_p99"] and out["a_ci"][0] > 0)
    log(f"{name}: " + json.dumps(out, default=str))
    cv = curve(Y, N, f["r_inf"]); cv.insert(0, "set", name)
    return out, cv


UNION = "--union" in sys.argv     # post-hoc diagnostic (PLAN Revision 2): Bostock-window majors widened to the env-(b) tremor episode


def union_majors(E):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import r23_base as B
    X = pd.read_parquet(f"{C.TB}/famday_tremor.parquet")
    q = pd.read_csv(f"{C.ROOT}/svi_science/results/backward/p4/envflag/tables/quiet_flag_env.csv", usecols=["date", "N_CN4"])
    act = set(q.date[q.N_CN4 >= 12])
    Xb = X[(X.src == "env") & (X.n30 > 0) & X.date.isin(act)].rename(columns={"fam": "fid"})
    pre = B.find_episodes(Xb, 115, "2003-01-01", "2010-01-05", B.RULE)
    E = E.copy()
    for i, r in E[E.kind == "bostock"].iterrows():
        h = pre[(pre.start <= r.end) & (pre.end >= r.start)]
        if len(h):
            E.loc[i, "start"] = min(r.start, h.start.min()); E.loc[i, "end"] = max(r.end, h.end.max())
    E["start_ts"] = pd.to_datetime(E.start); E["end_ts"] = pd.to_datetime(E.end)
    log("UNION majors: " + E[E.kind == "bostock"][["episode", "start", "end"]].to_string(index=False))
    return E


def main():
    global OUT
    F = C.families(); E = C.episodes()
    if UNION:
        OUT = f"{C.TB}/study3_union"; os.makedirs(OUT, exist_ok=True)
        E = union_majors(E)
    V = C.load_events(F, E); D = C.load_famdays(F, E)
    maj = E[E.major].reset_index(drop=True); big = E[E.big].reset_index(drop=True)
    clocks = {"all majors": maj, "big ETS": big}
    R = {"n_big_2003_13": int((big.start < "2014-01-01").sum()), "n_majors_2003_13": int((maj.start < "2014-01-01").sum()),
         "big_2003_13": big[big.start < "2014-01-01"].episode.tolist()}
    # ======================= (A) S7 quiet month: F30 of tau_prev, BE events, family+segment-preserving shuffle
    X, _ = C.tremor_famday(); Xe = X[X.n30 > 0]
    own_end = {}
    for cname, ref in clocks.items():
        rs = ref.start_ts.values.astype("datetime64[D]"); re_ = ref.end_ts.values.astype("datetime64[D]")
        for f, g in Xe.groupby("fam"):
            dfx = np.sort(pd.to_datetime(g.date).values.astype("datetime64[D]"))
            le = re_.copy()
            for j, (s, e) in enumerate(zip(rs, re_)):
                k = dfx[(dfx >= s) & (dfx <= e)]
                if k.size:
                    le[j] = k[-1]
            own_end[(cname, f)] = le

    FIDX = {f: i for i, f in enumerate(F.fid)}
    LEM = {}
    for cname, ref in clocks.items():
        re_ = ref.end_ts.values.astype("datetime64[D]")
        LEM[cname] = np.array([own_end.get((cname, f), re_) for f in F.fid])

    def taus(fams, dates, cname):
        ref = clocks[cname]
        tn, tp, _ = C.tau_clock(dates, ref.start_ts.values, ref.end_ts.values)
        d = pd.to_datetime(pd.Series(np.asarray(dates))).values.astype("datetime64[D]")
        re_ = ref.end_ts.values.astype("datetime64[D]")
        j = np.searchsorted(re_, d, side="left") - 1
        LE = LEM[cname]
        fi = np.asarray([FIDX[f] for f in fams]) if len(fams) < 50000 else pd.Series(fams).map(FIDX).to_numpy()
        ok = j >= 0
        t3 = np.full(len(d), np.nan)
        t3[ok] = (d[ok] - LE[fi[ok], j[ok]]).astype(int)
        return tp, tp - 7, t3

    rowsA = []
    poolBE = C.pools(D, D.BE.to_numpy()); poolBR = C.pools(D, D.BE_broad.to_numpy())
    for p in C.PERIODS:
        a, b = C.PERIODS[p]
        for setname, flag, pool in (("BE", "BE", poolBE), ("BE_broad", "BE_broad", poolBR)):
            for ser in ("data", "c1"):
                ev = V[(V.date >= a) & (V.date <= b) & V[flag] & (V.series == ser)]
                if len(ev) == 0:
                    continue
                keys = list(zip(ev.fam, ev.seg))
                sh = C.shuffle_dates(keys, pool, RNG, NSH) if ser == "data" else None
                for cname in clocks:
                    t1, t2, t3 = taus(ev.fam.to_numpy(), ev.date.to_numpy(), cname)
                    if sh is not None:
                        fam_rep = np.tile(ev.fam.to_numpy(), NSH)
                        s1, s2, s3 = taus(fam_rep, sh.ravel(), cname)
                        s1, s2, s3 = (x.reshape(NSH, -1) for x in (s1, s2, s3))
                    for dn, obs, shv, lo in (("D1", t1, s1 if sh is not None else None, 0), ("D2", t2, s2 if sh is not None else None, 0),
                                             ("D3", t3, s3 if sh is not None else None, 0)):
                        m = np.isfinite(obs) & (obs >= lo)
                        f30 = float(np.mean(obs[m] < 30)) if m.sum() else np.nan
                        row = dict(period=p, set=setname, series=ser, clock=cname, end_def=dn, n=int(m.sum()), F30=f30)
                        if shv is not None:
                            ms = np.isfinite(shv) & (shv >= lo)
                            f30s = np.array([np.mean(r[k] < 30) if k.sum() else np.nan for r, k in zip(shv, ms)])
                            row.update(shuffle_F30_median=float(np.nanmedian(f30s)), shuffle_F30_p01=float(np.nanpercentile(f30s, 1)),
                                       p_low=C.pct_lo(f30, f30s[np.isfinite(f30s)]))
                            # superposed table for D1
                            if dn == "D1":
                                bins = [0, 7, 14, 30, 60, 120, 240, 5000]
                                row["obs_bins"] = np.histogram(obs[m], bins)[0].tolist()
                                row["shuffle_bins_mean"] = np.round(np.mean([np.histogram(r[k], bins)[0] for r, k in zip(shv, ms)], 0), 2).tolist()
                        rowsA.append(row)
    # burst-level column (S7 deviation): each BE burst once, tau of its first event; shuffle = burst moved to a random
    # eligible day of its first family (same segment)
    for p in C.PERIODS:
        a, b = C.PERIODS[p]
        ev = V[(V.date >= a) & (V.date <= b) & V.BE & (V.series == "data")].sort_values("time_utc").copy()
        ev["burst"] = C.chain_bursts(ev.time_utc)
        fb = ev.groupby("burst").first()
        keys = list(zip(fb.fam, fb.seg))
        sh = C.shuffle_dates(keys, poolBE, RNG, NSH)
        for cname in clocks:
            t1, _, _ = taus(fb.fam.to_numpy(), fb.date.to_numpy(), cname)
            s1, _, _ = taus(np.tile(fb.fam.to_numpy(), NSH), sh.ravel(), cname); s1 = s1.reshape(NSH, -1)
            m = np.isfinite(t1); f30 = float(np.mean(t1[m] < 30))
            f30s = np.array([np.mean(r[np.isfinite(r)] < 30) for r in s1])
            rowsA.append(dict(period=p, set="BE bursts", series="data", clock=cname, end_def="D1", n=int(m.sum()), F30=f30,
                              shuffle_F30_median=float(np.median(f30s)), shuffle_F30_p01=float(np.percentile(f30s, 1)), p_low=C.pct_lo(f30, f30s)))
    A = pd.DataFrame(rowsA); A.to_csv(f"{OUT}/quiet_month_F30.csv", index=False)
    log("(A) F30 after previous end:\n" + A.drop(columns=["obs_bins", "shuffle_bins_mean"], errors="ignore").round(4).to_string(index=False))
    R["A"] = A.to_dict("records")

    # ======================= (B) P3 recharge: deep groups, daily counts / observable deep family-days
    deep = F.fid[F.deep].to_numpy()
    d0, d1 = np.datetime64("2003-01-01"), np.datetime64(C.P_END)
    days = np.arange(d0, d1 + 1); ds = pd.Series(days.astype(str))
    Dd = D[D.fam.isin(deep)]

    def daily(series="data", fids=deep, daycount=False):
        x = V[(V.series == series) & V.fam.isin(fids)]
        y = x.groupby("date").size().reindex(ds, fill_value=0).to_numpy()
        n = D[D.fam.isin(fids)].groupby("date").size().reindex(ds, fill_value=0).to_numpy()
        if daycount:
            y = (y > 0).astype(int); n = (n > 0).astype(int)
        return y, n
    outs, curves = [], []
    for p in C.PERIODS:
        a, b = C.PERIODS[p]
        sel = lambda Cc: Cc[(Cc.date >= np.datetime64(a)) & (Cc.date <= np.datetime64(b))]
        for cname in ("big ETS", "all majors"):
            ref = clocks[cname][["start_ts", "end_ts"]].rename(columns={"start_ts": "start", "end_ts": "end"})
            y, n = daily()
            Cc = sel(cycle_table(days, y, n, ref))
            o, cv = analyse(f"{p} | LFE deep groups | {cname} | D1", Cc); o.update(period=p, clock=cname, variant="D1 events")
            outs.append(o); curves.append(cv)
            if cname == "big ETS":
                C2 = Cc[Cc.t >= 8].copy(); C2["t"] = C2.t - 7
                o, cv = analyse(f"{p} | LFE deep groups | {cname} | D2 end+7", C2); o.update(period=p, clock=cname, variant="D2 events")
                outs.append(o); curves.append(cv)
                yd, nd = daily(daycount=True)
                o, cv = analyse(f"{p} | LFE deep groups | {cname} | D1 day-count", sel(cycle_table(days, yd, nd, ref)))
                o.update(period=p, clock=cname, variant="D1 day-count"); outs.append(o); curves.append(cv)
                yc, nc = daily("c1")
                o, cv = analyse(f"{p} | scrambled control deep | {cname} | D1", sel(cycle_table(days, yc, nc, ref)), surrogates=False)
                o.update(period=p, clock=cname, variant="D1 scrambled control"); outs.append(o); curves.append(cv)
                # without burst-day rows
                x = V[(V.series == "data") & V.fam.isin(deep) & (~V.burst_day_f)]
                yb = x.groupby("date").size().reindex(ds, fill_value=0).to_numpy()
                o, cv = analyse(f"{p} | LFE deep, no burst-day rows | {cname} | D1", sel(cycle_table(days, yb, n, ref)), surrogates=False)
                o.update(period=p, clock=cname, variant="D1 no burst_day rows"); outs.append(o); curves.append(cv)
    # ======================= (C) env CN4 deep tremor 2003-2013 (fixed 4-station network), PNSN 2010+ column
    TR = pd.read_parquet(f"{C.TB}/tremor_points.parquet", columns=["date", "src", "dmin_deep_km"])
    q = pd.read_csv(f"{C.ROOT}/svi_science/results/backward/p4/envflag/tables/quiet_flag_env.csv", usecols=["date", "cov_CN4"]).set_index("date")
    for src, (a, b) in (("env", ("2003-01-01", "2013-12-31")), ("pnsn", ("2010-01-22", C.P_END))):
        t = TR[(TR.src == src) & (TR.dmin_deep_km <= 15)]
        y = t.groupby("date").size().reindex(ds, fill_value=0).to_numpy()
        n = (q.cov_CN4.reindex(ds).fillna(0).to_numpy() >= 0.5).astype(int) if src == "env" else np.ones(len(ds), int)
        inr = (days >= np.datetime64(a)) & (days <= np.datetime64(b)); n = n * inr; y = y * inr
        for cname in ("big ETS", "all majors"):
            Cc = cycle_table(days, y, n, clocks[cname][["start_ts", "end_ts"]].rename(columns={"start_ts": "start", "end_ts": "end"}))
            Cc = Cc[Cc.n > 0]
            for p in ("2003-13", "2014-26"):
                pa, pb = C.PERIODS[p]
                Cp = Cc[(Cc.date >= np.datetime64(pa)) & (Cc.date <= np.datetime64(pb))]
                if Cp.y.sum() == 0:
                    continue
                o, cv = analyse(f"{p} | {src} deep tremor (<=15 km) | {cname} | D1", Cp)
                o.update(period=p, clock=cname, variant=f"{src} deep tremor"); outs.append(o); curves.append(cv)
    O = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in o.items()} for o in outs])
    O.to_csv(f"{OUT}/recharge_fits.csv", index=False)
    pd.concat(curves).to_csv(f"{OUT}/recharge_curves.csv", index=False)
    R["B"] = outs
    # ======================= verdict (2003-13), pre-registered
    a1 = A[(A.period == "2003-13") & (A.set == "BE") & (A.series == "data") & (A.clock == "all majors") & (A.end_def == "D1")].iloc[0]
    b1 = [o for o in outs if o["period"] == "2003-13" and o["clock"] == "big ETS" and o["variant"] == "D1 events"][0]
    passA = bool(a1.p_low < 0.01); passB = bool(b1["recovery_present"])
    wrongA = bool(a1.F30 >= a1.shuffle_F30_median); wrongB = bool(b1["D30"] <= 0)
    verdict = "confirmed" if passA and passB else ("partly" if passA or passB else ("rejected" if wrongA and wrongB else "undecided"))
    R["verdict_2003_13"] = dict(A_F30=float(a1.F30), A_shuffle_median=float(a1.shuffle_F30_median), A_p_low=float(a1.p_low), A_pass=passA,
                                B_D30=b1["D30"], B_p_surrogate=b1.get("p_surrogate"), B_a_ci=b1["a_ci"], B_tau=b1["tau"], B_tau_ci=b1["tau_ci"],
                                B_pass=passB, verdict=verdict)
    log(f"VERDICT 2003-13: {R['verdict_2003_13']}")
    C.dump(R, f"{OUT}/study3_summary.json")
    figure(pd.concat(curves), A)
    log("done")


def figure(CV, A):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    mids = np.array([4, 11, 22.5, 45, 75, 120, 195, 320])
    for nm, col in (("2003-13 | LFE deep groups | big ETS | D1", "#08519c"), ("2014-26 | LFE deep groups | big ETS | D1", "#fd8d3c"),
                    ("2003-13 | env deep tremor (<=15 km) | big ETS | D1", "#6baed6"), ("2014-26 | pnsn deep tremor (<=15 km) | big ETS | D1", "#fdae6b"),
                    ("2003-13 | scrambled control deep | big ETS | D1", "#bdbdbd")):
        c = CV[CV.set == nm]
        if len(c):
            ax[0].plot(mids, c.norm_rate, "o-", color=col, label=nm.replace(" | big ETS | D1", ""), lw=1.6)
    ax[0].axhline(1, color="k", lw=0.6, ls=":"); ax[0].set_xscale("log"); ax[0].set_xlabel("days since the previous big ETS ended")
    ax[0].set_ylabel("deep activity / fitted long-term level"); ax[0].legend(fontsize=7); ax[0].set_title("Deep activity after big ETS (shape)")
    a = A[(A.set == "BE") & (A.series == "data") & (A.end_def == "D1")]
    lab, obs, med, lo = [], [], [], []
    for _, r in a.iterrows():
        lab.append(f"{r.period}\n{r.clock}"); obs.append(r.F30); med.append(r.shuffle_F30_median); lo.append(r.shuffle_F30_p01)
    x = np.arange(len(lab))
    ax[1].bar(x - 0.18, obs, 0.36, color="#08519c", label="observed F30")
    ax[1].bar(x + 0.18, med, 0.36, color="#bdbdbd", label="shuffled days (median)")
    ax[1].scatter(x + 0.18, lo, marker="_", s=200, color="k", label="shuffle 1st percentile")
    ax[1].set_xticks(x); ax[1].set_xticklabels(lab, fontsize=7); ax[1].set_ylabel("share of between-episode LFEs < 30 d after the end")
    ax[1].legend(fontsize=7); ax[1].set_title("S7 quiet-month test (tremor-free between-episode LFEs)")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/r23_quiet_after_ets{'_union' if UNION else ''}.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
