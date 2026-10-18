#!/usr/bin/env python
"""q4_recharge_wa (2026-10-07) -- Washington-NW confirmatory test of the slow deep "recharge" seen in SVI (P3 step 4).
Rules frozen in svi_science/results/q4_recharge_wa/PLAN.md.  CPU only, < 20 GB.
    python svi_science/scripts/q4_recharge_wa.py build     # event table cache (P files + dedup flags, 2007-01 .. 2026-05)
    python svi_science/scripts/q4_recharge_wa.py run       # all tests
Fit / D30 / surrogate / bootstrap code copied from svi_science/scripts/p3_recharge.py (not imported: that module opens
its own log for writing at import).
"""
import os, sys, glob, json, time, zlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from multiprocessing import Pool

ROOT = "/home/jovyan/tremorferometry"
WP = f"{ROOT}/lfe_extension_wa/results/wa_prod"
OUT = f"{ROOT}/svi_science/results/q4_recharge_wa"
A = f"{ROOT}/svi_science/results/a_wa_between"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
B_EP = f"{ROOT}/svi_science/results/b_ets_record/episodes_b.csv"
S0 = f"{ROOT}/svi_science/results/s0"
P0, P1 = "2007-01-01", "2026-05-27"
ERAS = {"E1": ("2007-01-01", "2011-12-31"), "E2": ("2012-01-01", "2014-12-31"), "E3": ("2015-01-01", P1)}
TMAX = 400
NB = 1000
BINS = [1, 8, 15, 31, 61, 91, 151, 241, 401]
BLAB = ["0-7", "7-14", "14-30", "30-60", "60-90", "90-150", "150-240", "240-400"]
TT = np.arange(1, TMAX + 1, dtype=float)
LOGF = None


def log(m):
    global LOGF
    if LOGF is None:
        LOGF = open(f"{OUT}/q4.log", "a")
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


# ============================================================================================ build
def build():
    F = pd.read_csv(f"{A}/families_wa.csv")
    sg = F.set_index("fid").source_group
    cols = ["date", "series", "fid", "k", "P_own_WAVE_PNET"]
    parts = []
    for p in sorted(glob.glob(f"{WP}/conf/P/prod/P_*.parquet")):
        ym = p.split("_")[-1][:7]
        if ym > P1[:7]:
            continue
        D = pd.read_parquet(p, columns=cols + ["net", "reg"], filters=[("net", "==", "core"), ("reg", "==", "NW")])
        D = D[D.P_own_WAVE_PNET >= 0.75].drop(columns=["net", "reg"])
        if not len(D):
            continue
        fl = pd.read_parquet(f"{WP}/post_v2/dedup/flags/flags_{ym}.parquet",
                             columns=["date", "series", "fid", "k", "dup_D", "win_fid_D"])
        n0 = len(D)
        D = D.merge(fl, on=["date", "series", "fid", "k"], how="left")
        assert len(D) == n0, ym
        parts.append(D)
    E = pd.concat(parts, ignore_index=True)
    E = E[E.date <= P1]
    nmiss = int(E.dup_D.isna().sum())
    E["dup_D"] = E.dup_D.fillna(False).astype(bool)
    E["once_D"] = ~E.dup_D
    same = E.win_fid_D.map(sg).to_numpy() == E.fid.map(sg).to_numpy()
    E["once_group"] = ~(E.dup_D & same)
    E.to_parquet(f"{OUT}/events_p075.parquet")
    log(f"build: events P>=0.75 core NW {P0}..{P1}: {len(E)} (no flag row {nmiss}); by series "
        f"{E.groupby('series').size().to_dict()}; once_D {E.groupby('series').once_D.sum().to_dict()}")


# ============================================================================================ statistics (P3 code)
def nll(p, Y, N, tt):
    r, a, tau = np.exp(p[0]), p[1], np.exp(p[2])
    mu = np.maximum(r * (1 - a * np.exp(-tt / tau)), 1e-12)
    return float(np.sum(N * mu - Y * np.log(mu)))


def fit(Y, N):
    m = N > 0
    Y, N, tt = Y[m], N[m], TT[m]
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
    a, b = TT < 31, TT >= 31
    if N[a].sum() == 0 or N[b].sum() == 0 or Y[b].sum() == 0:
        return np.nan
    return float(1 - (Y[a].sum() / N[a].sum()) / (Y[b].sum() / N[b].sum()))


def beta(Y, N, t0=1):
    """Poisson rate ∝ t^beta over t >= t0 (exposure offset); profile likelihood in beta."""
    m = (N > 0) & (TT >= t0)
    Y, N, lt = Y[m], N[m], np.log(TT[m])
    if Y.sum() == 0:
        return np.nan
    sY, sYl = Y.sum(), (Y * lt).sum()
    f = lambda b: sY * np.log(np.sum(N * np.exp(b * lt))) - b * sYl
    return float(minimize_scalar(f, bounds=(-4, 4), method="bounded").x)


def agg(C):
    a = C.groupby("t").agg(Y=("y", "sum"), N=("n", "sum")).reindex(range(1, TMAX + 1), fill_value=0)
    return a.Y.to_numpy(float), a.N.to_numpy(float)


def curve(Y, N, r_ref):
    rows = []
    for lo, hi, lab in zip(BINS[:-1], BINS[1:], BLAB):
        s = slice(lo - 1, hi - 1)
        rows.append(dict(bin=lab, t_mid=(lo + hi - 1) / 2, events=float(Y[s].sum()), exposure=float(N[s].sum()),
                         rate=float(Y[s].sum() / N[s].sum()) if N[s].sum() else np.nan))
    df = pd.DataFrame(rows)
    df["norm_rate"] = df.rate / r_ref if r_ref > 0 else np.nan
    return df


def analyse(args):
    name, C, surrogates, nb = args
    rng = np.random.default_rng(zlib.crc32(name.encode()))
    Y, N = agg(C)
    f = fit(Y, N); D = d30(Y, N); B = beta(Y, N); BL = beta(Y, N, 31)
    out = dict(name=name, n_events=float(Y.sum()), exposure=float(N.sum()), n_cycles=int(C.cyc.nunique()), **f,
               D30=D, beta=B, beta_late=BL,
               rate_0_30=float(Y[:30].sum() / max(N[:30].sum(), 1e-12)),
               rate_31_400=float(Y[30:].sum() / max(N[30:].sum(), 1e-12)),
               rate_31_240=float(Y[30:240].sum() / max(N[30:240].sum(), 1e-12)),
               rate_241_400=float(Y[240:].sum() / max(N[240:].sum(), 1e-12)))
    out["mean_rate"] = float(Y.sum() / max(N.sum(), 1e-12))
    out["tau_at_limit"] = bool(np.isfinite(f["tau"]) and f["tau"] >= 395)
    if nb and Y.sum() > 0:
        cyc = sorted(C.cyc.unique())
        per = {c: agg(C[C.cyc == c]) for c in cyc}
        Ys = np.array([per[c][0] for c in cyc]); Ns = np.array([per[c][1] for c in cyc])
        bt, ba, bd, bb = [], [], [], []
        for _ in range(nb):
            k = rng.integers(0, len(cyc), len(cyc))
            Yb, Nb = Ys[k].sum(0), Ns[k].sum(0)
            fb = fit(Yb, Nb); bt.append(fb["tau"]); ba.append(fb["a"]); bd.append(d30(Yb, Nb)); bb.append(beta(Yb, Nb))
        for key, arr in (("tau", bt), ("a", ba), ("D30", bd), ("beta", bb)):
            arr = np.array(arr, float)
            out[f"{key}_ci"] = [float(np.nanpercentile(arr, 2.5)), float(np.nanpercentile(arr, 97.5))]
        out["tau_boot_share_at_limit"] = float(np.nanmean(np.array(bt) >= 395))
    if surrogates and Y.sum() > 0:
        sd, sb, sbl = [], [], []
        groups = [g.sort_values("t") for _, g in C.groupby("cyc")]
        arrs = [(g.t.to_numpy() - 1, g.y.to_numpy(float), g.n.to_numpy(float)) for g in groups]
        for _ in range(NB):
            Ysum = np.zeros(TMAX); Nsum = np.zeros(TMAX)
            for ti, yy, nn in arrs:
                k = rng.integers(0, len(ti))
                np.add.at(Ysum, ti, np.roll(yy, k)); np.add.at(Nsum, ti, np.roll(nn, k))
            sd.append(d30(Ysum, Nsum)); sb.append(beta(Ysum, Nsum)); sbl.append(beta(Ysum, Nsum, 31))
        sd, sb, sbl = np.array(sd), np.array(sb), np.array(sbl)
        out.update(sur_D30_median=float(np.nanmedian(sd)), sur_D30_p99=float(np.nanpercentile(sd, 99)),
                   p_D30=float((np.sum(sd >= D) + 1) / (NB + 1)),
                   sur_beta_median=float(np.nanmedian(sb)), sur_beta_p99=float(np.nanpercentile(sb, 99)),
                   p_beta=float((np.sum(sb >= B) + 1) / (NB + 1)),
                   p_beta_late=float((np.sum(sbl >= BL) + 1) / (NB + 1)) if np.isfinite(BL) else np.nan)
    cv = curve(Y, N, out["mean_rate"]); cv.insert(0, "set", name)
    return out, cv


# ============================================================================================ tables
def intervals(df):
    return df.start.values.astype("datetime64[D]"), df.end.values.astype("datetime64[D]")


def cycle_table(days, y, n, ref, excl, shift=0, ref_end_override=None):
    """t = days since the previous reference end; days inside any `excl` interval dropped (excl includes ref)."""
    rs, re_ = intervals(ref)
    xs, xe = intervals(excl)
    inside = np.zeros(days.size, bool)
    for s, e in zip(np.r_[rs, xs], np.r_[re_, xe]):
        inside |= (days >= s) & (days <= e)
    j = np.searchsorted(re_, days, side="left") - 1
    ok = (j >= 0) & ~inside
    t = np.where(ok, (days - re_[np.maximum(j, 0)]).astype(int), -1) - shift
    ok &= (t >= 1) & (t <= TMAX) & (n > 0)
    d = days[ok]
    era = np.where(d < np.datetime64("2012-01-01"), "E1", np.where(d < np.datetime64("2015-01-01"), "E2", "E3"))
    return pd.DataFrame(dict(date=d, cyc=j[ok], t=t[ok], y=y[ok], n=n[ok], era=era))


def era_sub(C, e):
    return C[C.era == e]


def pooled_norm(C):
    C = C.copy()
    for e, g in C.groupby("era"):
        r = g.y.sum() / g.n.sum()
        C.loc[g.index, "n"] = g.n * r if r > 0 else 0.0
    return C[C.n > 0]


# ============================================================================================ run
def run():
    t0 = time.time()
    F = pd.read_csv(f"{A}/families_wa.csv")
    ab = F[F.ab]
    q_slab = ab.slab2_dep.quantile(2 / 3); q_dd = ab.downdip_km.quantile(2 / 3); q_all = F.slab2_dep.quantile(2 / 3)
    sets = {"deep": ab.fid[ab.slab2_dep >= q_slab].to_numpy(),
            "deep_downdip": ab.fid[ab.downdip_km >= q_dd].to_numpy(),
            "deep_inclC": F.fid[F.slab2_dep >= q_all].to_numpy()}
    deep = sets["deep"]
    log(f"deep sets: thresholds slab2 {q_slab:.2f} km, downdip {q_dd:.2f} km, all-fam slab2 {q_all:.2f} km; sizes "
        f"{ {k: len(v) for k, v in sets.items()} }; deep source groups {F[F.fid.isin(deep)].source_group.nunique()}")
    Epi = pd.read_csv(f"{EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    maj = Epi[Epi.major].sort_values("start").reset_index(drop=True)
    Bp = pd.read_csv(B_EP, parse_dates=["start", "end"])
    Bn = Bp[(Bp.region == "NW")]
    big = Bn[Bn.big == True].sort_values("start").reset_index(drop=True)
    a2 = Bn[Bn.big_a2 == True].sort_values("start").reset_index(drop=True)
    log(f"majors {len(maj)}; big {len(big)}: {big.episode.tolist()}; a2 {len(a2)}")
    # big episodes must be majors (same windows)
    chk = big.merge(maj[["start", "end"]], on=["start", "end"], how="left", indicator=True)
    log(f"big windows identical to a major window: {int((chk._merge == 'both').sum())}/{len(big)}")

    E = pd.read_parquet(f"{OUT}/events_p075.parquet")
    R = pd.read_parquet(f"{WP}/post_v2/dedup/rates_daily_v2_dedup.parquet", columns=["unit_kind", "fid", "date", "region"])
    R = R[(R.region == "NW") & (R.unit_kind == "family") & (R.date >= P0) & (R.date <= P1)]
    R["fid"] = R.fid.astype(int)
    L = pd.read_csv(f"{A}/live_core_by_day.csv")
    days = np.arange(np.datetime64(P0), np.datetime64(P1) + 1)
    ds = pd.Series(days.astype(str))
    L = L.set_index("date").n_core_live.reindex(ds).to_numpy(float)
    era_d = np.where(days < np.datetime64("2012-01-01"), "E1", np.where(days < np.datetime64("2015-01-01"), "E2", "E3"))
    med = {e: np.nanmedian(L[era_d == e]) for e in ERAS}
    live_ok = np.array([L[i] >= 0.5 * med[era_d[i]] for i in range(days.size)])
    log(f"live-core median per era {med}; days failing the 50 % live rule: {int((~live_ok).sum())}")

    def daily(fids, series="data", pcut=0.75, count="once_D", live=False):
        e = E[(E.series == series) & E.fid.isin(fids) & (E.P_own_WAVE_PNET >= pcut)]
        if count != "raw":
            e = e[e[count]]
        y = e.groupby("date").size().reindex(ds, fill_value=0).to_numpy(float)
        n = R[R.fid.isin(fids)].groupby("date").size().reindex(ds, fill_value=0).to_numpy(float)
        if live:
            y = np.where(live_ok, y, 0); n = np.where(live_ok, n, 0)
        return y, n

    jobs = []

    def add(name, C, sur=True, nb=NB):
        jobs.append((name, C, sur, nb))

    y, n = daily(deep)
    log(f"deep data events (counted once, P>=0.75) per era, all days: "
        f"{ {e: int(y[era_d == e].sum()) for e in ERAS} }")
    clocks = {"big": (big, maj), "all majors": (maj, maj), "a2": (a2, maj), "big P3-as-coded": (big, big)}
    HEAD = {}
    for cname, (ref, excl) in clocks.items():
        C = cycle_table(days, y, n, ref, excl)
        for e in ERAS:
            add(f"LFE deep | clock {cname} | D1 | {e}", era_sub(C, e))
        add(f"LFE deep | clock {cname} | D1 | pooled-normalised", pooled_norm(C))
        if cname == "big":
            HEAD["C"] = C
            # D2
            C2 = C[C.t >= 8].copy(); C2["t"] = C2.t - 7
            for e in ERAS:
                add(f"LFE deep | clock big | D2 end+7 | {e}", era_sub(C2, e))
            add(f"LFE deep | clock big | D2 end+7 | pooled-normalised", pooled_norm(C2))
    # D1 variants on the big clock
    var = {"deep_downdip": dict(fids=sets["deep_downdip"]), "deep_inclC": dict(fids=sets["deep_inclC"]),
           "P>=0.85": dict(fids=deep, pcut=0.85), "source-group once": dict(fids=deep, count="once_group"),
           "raw": dict(fids=deep, count="raw"), "live-core days": dict(fids=deep, live=True)}
    for vname, kw in var.items():
        fids = kw.pop("fids")
        yv, nv = daily(fids, **kw)
        C = cycle_table(days, yv, nv, big, maj)
        for e in ERAS:
            add(f"LFE {vname} | clock big | D1 | {e}", era_sub(C, e), nb=300)
        add(f"LFE {vname} | clock big | D1 | pooled-normalised", pooled_norm(C), nb=300)
    # controls
    for ser, lab in (("rev", "time-reversed control"), ("c1", "scrambled control")):
        yc, nc = daily(deep, ser)
        C = cycle_table(days, yc, nc, big, maj)
        for e in ERAS:
            add(f"{lab} deep | clock big | D1 | {e}", era_sub(C, e), sur=True, nb=0)
        add(f"{lab} deep | clock big | D1 | all eras (not normalised)", C, sur=True, nb=0)
    # D3 family's own end (family-day table), big clock, no surrogate
    X = pd.read_parquet(f"{WP}/post/tremor_family_day.parquet")
    Xe = X[(X.n30 > 0)]
    rs, re_ = intervals(big)
    ms, me = intervals(maj)
    Ed = E[(E.series == "data") & E.fid.isin(deep) & E.once_D]
    cnt = Ed.groupby(["fid", "date"]).size()
    rows = []
    for f in deep:
        dfx = np.sort(pd.to_datetime(Xe.date[Xe.fid == f]).values.astype("datetime64[D]"))
        le = re_.copy()
        for j, (s, e) in enumerate(zip(rs, re_)):
            k = dfx[(dfx >= s) & (dfx <= e)]
            if k.size:
                le[j] = k[-1]
        od = np.sort(pd.to_datetime(R.date[R.fid == f]).values.astype("datetime64[D]"))
        inside = np.zeros(od.size, bool)
        for s, e in zip(np.r_[rs, ms], np.r_[re_, me]):
            inside |= (od >= s) & (od <= e)
        j = np.searchsorted(re_, od, side="left") - 1
        ok = (j >= 0) & ~inside
        t = (od - le[np.maximum(j, 0)]).astype(int)
        ok &= (t >= 1) & (t <= TMAX)
        dd = od[ok]
        yy = np.array([cnt.get((f, str(d)), 0) for d in dd], float)
        era = np.where(dd < np.datetime64("2012-01-01"), "E1", np.where(dd < np.datetime64("2015-01-01"), "E2", "E3"))
        rows.append(pd.DataFrame(dict(date=dd, cyc=j[ok], t=t[ok], y=yy, n=1.0, era=era)))
    C3 = pd.concat(rows, ignore_index=True)
    for e in ERAS:
        add(f"LFE deep | clock big | D3 family end | {e}", era_sub(C3, e), sur=False, nb=300)
    add(f"LFE deep | clock big | D3 family end | pooled-normalised", pooled_norm(C3), sur=False, nb=300)

    # PNSN deep tremor
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    t9 = pd.read_csv(f"{ROOT}/lfe_array_discovery/results/ab_royer/tremor/pnsn_tremor_2009.csv", skipinitialspace=True)
    t9.columns = [c.strip() for c in t9.columns]
    t9 = t9.rename(columns={"starttime": "time"})[["time", "lat", "lon"]]
    tr = pd.concat([t9, tr], ignore_index=True)
    tr["date"] = tr.time.astype(str).str.strip().str[:10]
    tr = tr[(tr.date >= "2009-08-06") & (tr.date <= P1)]
    tr = tr[tr.lat.between(46.5, 49.0) & tr.lon.between(-125.0, -121.7)].copy()
    Fd = F.set_index("fid")
    dmin = np.full(len(tr), np.inf)
    for f in deep:
        dmin = np.minimum(dmin, hav_km(tr.lat.to_numpy(float), tr.lon.to_numpy(float), Fd.lat[f], Fd.lon[f]))
    tr = tr[dmin <= 15]
    log(f"PNSN deep tremor epicentres within 15 km of a deep family 2009-08-06..{P1}: {len(tr)}")
    yt = tr.groupby("date").size().reindex(ds, fill_value=0).to_numpy(float)
    nt = (days >= np.datetime64("2009-08-06")).astype(float)
    for cname, (ref, excl) in (("big", (big, maj)), ("all majors", (maj, maj)), ("big P3-as-coded", (big, big))):
        Ct = cycle_table(days, yt, nt, ref, excl)
        for e in ERAS:
            add(f"PNSN deep tremor | clock {cname} | D1 | {e}", era_sub(Ct, e))
        add(f"PNSN deep tremor | clock {cname} | D1 | 2009-08..2026-05 pooled", Ct)
        if cname == "big":
            Ct2 = Ct[Ct.t >= 8].copy(); Ct2["t"] = Ct2.t - 7
            add(f"PNSN deep tremor | clock big | D2 end+7 | 2009-08..2026-05 pooled", Ct2, nb=300)

    # SVI re-check (P3 data, big clock, drop only big (as coded) vs drop all 39 majors)
    Fs = pd.read_csv(f"{S0}/families.csv")
    sdeep = Fs.fid[Fs.group.isin([1, 4, 10, 12])].to_numpy()
    Es = pd.read_csv(f"{EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    maj39 = Es[Es.major].sort_values("start").reset_index(drop=True)
    Bs = Bp[Bp.region == "SVI"]
    sbig = Bs[(Bs.big == True) & (Bs.start >= "2013-01-01")].sort_values("start").reset_index(drop=True)
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8"])
    T = T[(T.date >= "2014-01-01") & (T.date <= P1) & (T.P_same >= 0.75) & T.keep_xfam8 & (T.series == "data")]
    Dd = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable"])
    Dd = Dd[(Dd.run == "prod") & Dd.observable & (Dd.date >= "2014-01-01") & (Dd.date <= P1)]
    sdays = np.arange(np.datetime64("2014-01-01"), np.datetime64(P1) + 1)
    sds = pd.Series(sdays.astype(str))
    ys = T[T.fam.isin(sdeep)].groupby("date").size().reindex(sds, fill_value=0).to_numpy(float)
    ns = Dd[Dd.fid.isin(sdeep)].groupby("date").size().reindex(sds, fill_value=0).to_numpy(float)
    both = pd.concat([maj39[["start", "end"]], sbig[["start", "end"]]])
    add("SVI deep groups | clock big | D1 | P3 as coded (drop big only)", cycle_table(sdays, ys, ns, sbig, sbig))
    add("SVI deep groups | clock big | D1 | drop all majors (Q4 definition)", cycle_table(sdays, ys, ns, sbig, both))
    Cs = cycle_table(sdays, ys, ns, sbig, both)
    for lab, (a_, b_) in (("2014-2019", ("2014-01-01", "2019-12-31")), ("2020-2026", ("2020-01-01", P1))):
        add(f"SVI deep groups | clock big | D1 | drop all majors | {lab}",
            Cs[(Cs.date >= np.datetime64(a_)) & (Cs.date <= np.datetime64(b_))], nb=300)

    log(f"{len(jobs)} analyses queued (setup {time.time() - t0:.0f} s)")
    with Pool(16) as pool:
        res = pool.map(analyse, jobs, chunksize=1)
    outs = [r[0] for r in res]; curves = [r[1] for r in res]
    O = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in o.items()} for o in outs])
    O.to_csv(f"{OUT}/recharge_fits_wa.csv", index=False)
    pd.concat(curves).to_csv(f"{OUT}/recharge_curves_wa.csv", index=False)
    for o in outs:
        log(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items()}, default=str))

    # ---------------- verdict (pre-registered)
    g = {o["name"]: o for o in outs}
    h = g["LFE deep | clock big | D1 | E3"]
    e1, e2 = g["LFE deep | clock big | D1 | E1"], g["LFE deep | clock big | D1 | E2"]
    pt = g["PNSN deep tremor | clock big | D1 | 2009-08..2026-05 pooled"]
    c1 = h["D30"] > 0 and h["p_D30"] < 0.01 and h["beta"] > 0 and h["p_beta"] < 0.01
    c2 = all(x["D30"] > 0 and x["beta"] > 0 for x in (e1, e2))
    c3 = pt["D30"] > 0 and pt["beta"] > 0
    if c1 and c2 and c3:
        verdict = "SUPPORTED"
    elif (not h["D30"] > 0) or (not h["beta"] > 0) or h["p_D30"] >= 0.05:
        verdict = "NOT SUPPORTED"
    else:
        miss = [m for m, ok in (("E3 p < 0.01 (D30 and beta)", c1), ("E1/E2 same sign", c2), ("PNSN same sign", c3)) if not ok]
        verdict = "same direction, not confirmed (misses: " + "; ".join(miss) + ")"
    svi_ci = [85, 400]
    ts = {}
    for k in ("LFE deep | clock big | D1 | E3", "LFE deep | clock big | D1 | pooled-normalised"):
        ci = json.loads(json.dumps(g[k].get("tau_ci", [np.nan, np.nan])))
        ov = not (ci[1] < svi_ci[0] or ci[0] > svi_ci[1])
        ts[k] = dict(tau=g[k]["tau"], tau_ci=ci, at_limit=g[k]["tau_at_limit"], overlaps_SVI_85_400=ov)
    summ = dict(verdict=verdict, criteria=dict(E3=c1, E1E2_sign=c2, PNSN_sign=c3), time_scale=ts,
                deep_sets={k: [int(x) for x in v] for k, v in sets.items()},
                thresholds=dict(slab2=q_slab, downdip=q_dd, slab2_all=q_all), big=big.episode.tolist(),
                runtime_s=time.time() - t0)
    json.dump(summ, open(f"{OUT}/q4_summary.json", "w"), indent=1, default=str)
    log(f"VERDICT {verdict}; {json.dumps(summ['criteria'])}; time scale {json.dumps(ts, default=str)}")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    {"build": build, "run": run}[sys.argv[1]]()
