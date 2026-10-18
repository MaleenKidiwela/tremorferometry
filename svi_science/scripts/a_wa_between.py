#!/usr/bin/env python
"""a_wa_between (2026-10-07) -- Part 2 item A: Washington-NW repeat of S7 (LFEs between ETS episodes).
Rules frozen in svi_science/results/a_wa_between/PLAN_A.md before any number.  Needs a_wa_build.py outputs.
CPU only, < 15 GB.  Logic copied from s7_between.py; adaptations are listed in PLAN_A.md.

    python svi_science/scripts/a_wa_between.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/a_wa_between"
FIG = f"{ROOT}/svi_science/figures"
WP = f"{ROOT}/lfe_extension_wa/results/wa_prod"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
sys.path.insert(0, f"{ROOT}/svi_science/scripts")
import s3_build as S3B           # import only (frame projection)
FR = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["nw_frame"]
P_END = pd.Timestamp("2026-05-27")
ERAS = {"E3": ("2015-01-01", "2026-05-27"), "E2": ("2012-01-01", "2014-12-31"), "E1": ("2009-08-06", "2011-12-31")}
NSH = 1000
RNG = np.random.default_rng(20261007)
WA12 = ["CHCM", "KTBW", "P426", "ELSR", "P437", "CUSH", "P419", "PFLD", "SSHO", "SEAT", "SMAI", "P399"]
LOG = open(f"{OUT}/a.log", "w")
RES = {}


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def along(lat, lon):
    return S3B.project(FR, np.asarray(lat), np.asarray(lon))[0]


def pct_rank(obs, ctrl):
    ctrl = np.asarray(ctrl); return float((np.sum(ctrl >= obs) + 1) / (ctrl.size + 1))


def p_low(obs, ctrl):
    ctrl = np.asarray(ctrl); return float((np.sum(ctrl <= obs) + 1) / (ctrl.size + 1))


def pad_mask(dates, E, pad=2):
    d = pd.to_datetime(pd.Series(dates)).values
    m = np.zeros(len(d), bool)
    for s, e in zip(E.start.values, E.end.values):
        m |= (d >= s - np.timedelta64(pad, "D")) & (d <= e + np.timedelta64(pad, "D"))
    return m


def era_of(dates):
    d = np.asarray(dates).astype(str)
    out = np.full(d.size, "", object)
    for k, (a, b) in ERAS.items():
        out[(d >= a) & (d <= b)] = k
    return out


# ------------------------------------------------------------------------------------------- load
def load():
    F = pd.read_csv(f"{OUT}/families_wa.csv").set_index("fid", drop=False)
    E = pd.read_csv(f"{EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    T = pd.read_parquet(f"{OUT}/events_all.parquet")
    X = pd.read_parquet(f"{OUT}/exposure_famday.parquet")
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["time"])
    tr["time"] = tr.time.astype("datetime64[ns]")
    tr = tr[(tr.time >= pd.Timestamp("2009-07-01")) & tr.lat.between(45.8, 50.0) & tr.lon.between(-125.8, -120.5)]
    tr["date"] = tr.time.dt.strftime("%Y-%m-%d")
    return F, E, T, X, tr.reset_index(drop=True)


def classify(T, E):
    T = T.copy()
    T["date_ts"] = pd.to_datetime(T.date)
    T["era"] = era_of(T.date)
    T["in_pad"] = pad_mask(T.date, E, 2)
    T["in_ep_nopad"] = pad_mask(T.date, E, 0)
    T["in_major"] = pad_mask(T.date, E[E.major], 0)
    cov = T.day_class != "unknown"
    T["BE"] = ~T.in_pad & (T.day_class == "quiet")
    T["BE_broad"] = ~T.in_pad & cov
    T["out_nopad"] = ~T.in_ep_nopad & cov
    T["BE60"] = T.BE & (T.n_trem60 == 0)
    return T


# ------------------------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    log("=== A (WA-NW S7 repeat) start ===")
    F, E, T, X, tr = load()
    maj = E[E.major].sort_values("start").reset_index(drop=True)
    maj_start = maj.start.values; maj_end = maj.end.values
    log(f"episodes {len(E)} (majors {len(maj)}); families {len(F)}; primary A/B {int(F.ab.sum())}")
    n_after = int(((T.date > "2026-05-27") & (T.series == "data") & (T.P_own_WAVE_PNET >= .75) & T.once_D).sum())
    T = classify(T[(T.date >= ERAS["E1"][0]) & (T.date <= ERAS["E3"][1])], E)
    T = T[T.era != ""].copy()
    T["raw"] = True
    RES["n_after_pnsn_end_not_used"] = n_after
    log(f"P>=0.75 counted-once data events after the PNSN end (2026-05-28 ->, not used): {n_after}")
    for x in ("dep", "slab2_dep", "along_km", "downdip_km", "lat", "lon", "quality", "ab", "source_group"):
        T[x] = T.fid.map(F[x])
    T["time_PT"] = T.time_utc.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")

    # ---------------- inventory
    inv = []
    for era in ("E3", "E2", "E1"):
        for P in (0.75, 0.5, 0.85):
            for ser in ("data", "rev", "c1"):
                x = T[(T.era == era) & (T.series == ser) & (T.P_own_WAVE_PNET >= P)]
                for cnt, m in (("once_D", x.once_D), ("once_group", x.once_group), ("raw", np.ones(len(x), bool))):
                    y = x[np.asarray(m)]
                    inv.append(dict(era=era, P=P, series=ser, counting=cnt, all=len(y), BE=int(y.BE.sum()),
                                    BE_broad=int(y.BE_broad.sum()), outside_episodes_nopad=int(y.out_nopad.sum()),
                                    BE_no_tremor_60km=int(y.BE60.sum()), quiet_any=int((y.day_class == "quiet").sum())))
    INV = pd.DataFrame(inv); INV.to_csv(f"{OUT}/inventory_definitions.csv", index=False)
    log("definitions:\n" + INV[INV.counting != "once_group"].to_string(index=False))
    log("source-group once:\n" + INV[INV.counting == "once_group"].to_string(index=False))

    def sel(era, ser="data", P=0.75, cnt="once_D"):
        return T[(T.era == era) & (T.series == ser) & (T.P_own_WAVE_PNET >= P) & T[cnt]]

    # ---------------- eligible family-days per era
    X = X.copy()
    X["era"] = era_of(X.date); X = X[X.era != ""]
    X["in_pad"] = pad_mask(X.date, E, 2)
    EL_all = X[~X.in_pad & (X.day_class == "quiet")][["fid", "date", "era", "n_core_live", "n60"]].reset_index(drop=True)
    ELb_all = X[~X.in_pad & (X.day_class != "unknown")][["fid", "date", "era"]].reset_index(drop=True)
    log("eligible family-days BE / BE_broad / observable by era: " +
        str({e: (int((EL_all.era == e).sum()), int((ELb_all.era == e).sum()), int((X.era == e).sum())) for e in ERAS}))

    def pool_of(EL):
        return {f: g.date.to_numpy() for f, g in EL.groupby("fid")}

    def shuffle_dates(fams, pool, n=NSH):
        out = np.empty((n, len(fams)), dtype=object)
        for j, f in enumerate(fams):
            p = pool[f]; out[:, j] = p[RNG.integers(0, p.size, n)]
        return out

    def tau_next(dates):
        d = pd.to_datetime(pd.Series(dates)).values.astype("datetime64[D]")
        ms = maj_start.astype("datetime64[D]"); me = maj_end.astype("datetime64[D]")
        i = np.searchsorted(ms, d, side="right")
        tn = np.where(i < ms.size, (ms[np.minimum(i, ms.size - 1)] - d).astype(int), -1).astype(float)
        tn[i >= ms.size] = np.nan
        j = np.searchsorted(me, d, side="left") - 1
        tp = np.where(j >= 0, (d - me[np.maximum(j, 0)]).astype(int), -1).astype(float)
        tp[j < 0] = np.nan
        nxt = np.where(i < ms.size, i, -1)
        return tn, tp, nxt

    def chain_ids(times):
        t = times.astype("datetime64[ns]").astype("int64") / 3.6e12
        o = np.argsort(t); b = np.empty(t.size, int)
        b[o] = np.r_[0, np.cumsum(np.diff(t[o]) >= 24)]
        return b

    bins = [0, 7, 14, 30, 60, 120, 240, 5000]
    lab = ["0-7", "7-14", "14-30", "30-60", "60-120", "120-240", ">240"]

    def epoch_table(evx, ELx, pool, name, col="tau_next"):
        e = evx[np.isfinite(evx[col])]
        expo = np.histogram(ELx[col].dropna(), bins)[0]
        obs = np.histogram(e[col], bins)[0]
        shd = shuffle_dates(e.fid.to_numpy(), pool)
        dd = pd.Series(shd.ravel()).unique()
        tmap = dict(zip(dd, tau_next(dd)[0 if col == "tau_next" else 1]))
        sht = np.vectorize(tmap.get)(shd).astype(float)
        shh = np.array([np.histogram(r[np.isfinite(r)], bins)[0] for r in sht])
        tab = pd.DataFrame(dict(set=name, axis=col, bin=lab, events=obs, elig_famdays=expo,
                                rate_per_1000=1000 * obs / np.maximum(expo, 1),
                                shuffle_mean=shh.mean(0), shuffle_p01=np.percentile(shh, 1, 0),
                                shuffle_p99=np.percentile(shh, 99, 0)))
        f30 = float(np.mean(e[col] < 30)) if len(e) else np.nan
        f30s = np.nanmean(np.where(np.isfinite(sht), sht < 30, np.nan), axis=1)
        return tab, f30, f30s, sht

    pe_cache = {}
    tabs_all = []
    era_res = {}
    for era in ("E3", "E2", "E1"):
        R = {}
        A = sel(era); EL = EL_all[EL_all.era == era].copy(); ELb = ELb_all[ELb_all.era == era].copy()
        el_by_f = pool_of(EL); elb_by_f = pool_of(ELb)
        ev = A[A.BE].copy(); evb = A[A.BE_broad].copy()
        evr = sel(era, "rev")[lambda d: d.BE].copy(); evc = sel(era, "c1")[lambda d: d.BE].copy()
        ets = A[A.in_major]
        assert not (set(ev.fid) - set(el_by_f)), set(ev.fid) - set(el_by_f)
        R["n_BE"] = len(ev); R["n_BE_rev"] = len(evr); R["n_BE_c1"] = len(evc); R["n_BE_broad"] = len(evb)
        R["n_BE_broad_rev"] = int(sel(era, "rev").BE_broad.sum()); R["n_BE_broad_c1"] = int(sel(era, "c1").BE_broad.sum())
        R["control_share_rev"] = len(evr) / max(len(ev), 1); R["control_share_c1"] = len(evc) / max(len(ev), 1)
        R["n_ETS_events"] = len(ets); R["n_eligible_famdays"] = len(EL)
        R["rate_BE_per_1000_famdays"] = 1000 * len(ev) / max(len(EL), 1)
        log(f"[{era}] BE data {len(ev)}, time-reversed {len(evr)}, scrambled {len(evc)}; BE_broad {len(evb)}; "
            f"ETS(majors) events {len(ets)}; eligible family-days {len(EL)}")
        for x in (ev, evb, evr, evc):
            x["tau_next"], x["tau_prev"], x["next_major_idx"] = tau_next(x.date)
            x["next_major"] = [maj.episode.iloc[k] if k >= 0 else "" for k in x.next_major_idx]
            x["burst"] = chain_ids(x.time_utc.to_numpy()) if len(x) else []
        EL["tau_next"], EL["tau_prev"], _ = tau_next(EL.date)
        ELb["tau_next"], ELb["tau_prev"], _ = tau_next(ELb.date)

        # ---- clustering
        sizes = ev.groupby("burst").size().to_numpy()
        obs_days = ev.date.nunique()
        sh = shuffle_dates(ev.fid.to_numpy(), el_by_f)
        sh_days = np.array([len(set(r)) for r in sh])
        R["cluster"] = dict(distinct_days=int(obs_days), shuffle_days_median=float(np.median(sh_days)),
                            shuffle_days_p01=float(np.percentile(sh_days, 1)), p_fewer_days=p_low(obs_days, sh_days),
                            clustered=bool(obs_days < np.percentile(sh_days, 1)),
                            n_bursts=int(sizes.size), n_bursts_ge3=int((sizes >= 3).sum()),
                            share_in_bursts_ge3=float(sizes[sizes >= 3].sum() / max(sizes.sum(), 1)),
                            largest_bursts=sorted(sizes.tolist(), reverse=True)[:10])
        log(f"[{era}] clustering: {R['cluster']}")

        # ---- next / previous major
        out2 = {}
        for name, evx, ELx, pool in (("BE", ev, EL, el_by_f), ("BE_broad", evb, ELb, elb_by_f)):
            for col in ("tau_next", "tau_prev"):
                tab, f30, f30s, _ = epoch_table(evx, ELx, pool, name, col)
                tab.insert(0, "era", era); tabs_all.append(tab)
                out2[f"{name}_{col}"] = dict(n=int(np.isfinite(evx[col]).sum()), F30=f30,
                                             shuffle_F30_median=float(np.nanmedian(f30s)),
                                             shuffle_F30_p01=float(np.nanpercentile(f30s, 1)),
                                             shuffle_F30_p99=float(np.nanpercentile(f30s, 99)),
                                             p=pct_rank(f30, f30s), p_low=p_low(f30, f30s))
            e2 = evx[(evx.tau_prev >= 30) & np.isfinite(evx.tau_next)]
            ELx2 = ELx[ELx.tau_prev >= 30]
            pool2 = pool_of(ELx2)
            e2 = e2[e2.fid.isin(pool2.keys())]
            tab, f30, f30s, _ = epoch_table(e2, ELx2, pool2, name + "_prev>=30", "tau_next")
            tab.insert(0, "era", era); tabs_all.append(tab)
            out2[f"{name}_tau_next_prev>=30"] = dict(n=int(len(e2)), F30=f30, shuffle_F30_median=float(np.nanmedian(f30s)),
                                                     shuffle_F30_p99=float(np.nanpercentile(f30s, 99)),
                                                     p=pct_rank(f30, f30s), p_low=p_low(f30, f30s))
        R["next_ets"] = out2
        R["ramp_pass"] = bool(out2["BE_tau_next"]["F30"] > out2["BE_tau_next"]["shuffle_F30_p99"]
                              and out2["BE_tau_next_prev>=30"]["F30"] > out2["BE_tau_next_prev>=30"]["shuffle_F30_p99"])
        # pre-registered confirmatory: month after a major ends
        c = out2["BE_tau_prev"]
        if c["n"] < 20:
            verdict = "can't tell (fewer than 20 events)"
        elif c["F30"] < c["shuffle_F30_median"] and c["p_low"] < 0.01:
            verdict = "REPEATS"
        elif c["F30"] >= c["shuffle_F30_median"] or c["p_low"] >= 0.05:
            verdict = "DOES NOT REPEAT"
        else:
            verdict = "same direction, not confirmed"
        R["post_ets_month"] = dict(c, verdict=verdict)
        log(f"[{era}] next-ETS: {json.dumps(out2)}  ramp_pass={R['ramp_pass']}  post-ETS month: {verdict}")

        # ---- burst level
        bu = ev.groupby("burst").agg(fid=("fid", "first"), date=("date", "first"), tau_next=("tau_next", "first"),
                                     tau_prev=("tau_prev", "first"), n=("fid", "size"))
        bl = {"n_bursts": int(len(bu))}
        for col in ("tau_next", "tau_prev"):
            bun = bu[np.isfinite(bu[col])]
            if len(bun) == 0:
                continue
            shb = shuffle_dates(bun.fid.to_numpy(), el_by_f)
            ddb = pd.Series(shb.ravel()).unique(); tmb = dict(zip(ddb, tau_next(ddb)[0 if col == "tau_next" else 1]))
            stb = np.vectorize(tmb.get)(shb).astype(float)
            f = float(np.mean(bun[col] < 30)); fs = np.nanmean(np.where(np.isfinite(stb), stb < 30, np.nan), 1)
            bl[col] = dict(n=int(len(bun)), F30=f, shuffle_F30_median=float(np.nanmedian(fs)), p_high=pct_rank(f, fs),
                           p_low=p_low(f, fs), median_obs=float(bun[col].median()),
                           median_shuffle=float(np.median(np.nanmedian(stb, 1))))
            if col == "tau_next":
                g = float(np.mean(bun[col] > 240)); gs = np.nanmean(np.where(np.isfinite(stb), stb > 240, np.nan), 1)
                bl[col].update(F_gt240=g, shuffle_F_gt240_median=float(np.nanmedian(gs)), p_gt240_high=pct_rank(g, gs))
        R["burst_level"] = bl
        log(f"[{era}] burst level: {bl}")

        # ---- down-dip (A/B families; ETS shares from the same set)
        fam = F[["fid", "source_group", "array_group", "copy_of", "quality", "ab", "dep", "slab2_dep", "along_km",
                 "downdip_km", "lat", "lon", "erh"]].copy()
        fam["n_BE"] = fam.fid.map(ev.fid.value_counts()).fillna(0).astype(int)
        fam["n_BE_rev"] = fam.fid.map(evr.fid.value_counts()).fillna(0).astype(int)
        fam["n_BE_c1"] = fam.fid.map(evc.fid.value_counts()).fillna(0).astype(int)
        fam["n_BE_broad"] = fam.fid.map(evb.fid.value_counts()).fillna(0).astype(int)
        fam["n_ETS"] = fam.fid.map(ets.fid.value_counts()).fillna(0).astype(int)
        fam["n_all"] = fam.fid.map(A.fid.value_counts()).fillna(0).astype(int)
        fam["elig_days_BE"] = fam.fid.map(EL.fid.value_counts()).fillna(0).astype(int)
        fam["BE_days"] = fam.fid.map(ev.groupby("fid").date.nunique()).fillna(0).astype(int)
        fam["BE_over_ETS"] = (fam.n_BE + 0.5) / (fam.n_ETS + 0.5)
        fam.to_csv(f"{OUT}/families_between_{era}.csv", index=False)
        fab = fam[fam.ab & fam.dep.notna()]
        pe = (fab.n_ETS / fab.n_ETS.sum()).to_numpy()

        def depth_test(evx, name, col="dep", fset=fab, nd=10000):
            e = evx[evx.fid.isin(fset.fid)]
            n = len(e)
            if n == 0:
                return dict(set=name, metric=col, n=0)
            pe_ = (fset.n_ETS / fset.n_ETS.sum()).to_numpy()
            draws = RNG.multinomial(n, pe_, size=nd)
            md = draws @ fset[col].to_numpy() / n
            o = float(e.fid.map(F[col]).mean())
            return dict(set=name, metric=col, n=n, mean_obs=o, mean_ETSdraw=float(md.mean()), p99=float(np.percentile(md, 99)),
                        p01=float(np.percentile(md, 1)), p_deeper=pct_rank(o, md), deeper=bool(o > np.percentile(md, 99)))
        fall = fam[fam.slab2_dep.notna()]
        DT = [depth_test(ev, "BE"), depth_test(ev, "BE", "slab2_dep"), depth_test(ev, "BE", "downdip_km"),
              depth_test(ev, "BE all families", "slab2_dep", fall),
              depth_test(evb, "BE_broad"), depth_test(evr, "BE time-reversed control"), depth_test(evc, "BE scrambled control"),
              depth_test(sel(era, P=0.85)[lambda d: d.BE], "BE P>=0.85"),
              depth_test(sel(era, cnt="once_group")[lambda d: d.BE], "BE source-group once"),
              depth_test(sel(era, cnt="raw")[lambda d: d.BE], "BE raw (not de-duplicated)")]
        # one per burst: median depth of the A/B events in each burst
        evab = ev[ev.fid.isin(fab.fid)].copy(); evab["d"] = evab.fid.map(F.dep)
        bd = evab.groupby("burst").d.median()
        if len(bd):
            mdb = RNG.multinomial(len(bd), pe, size=10000) @ fab.dep.to_numpy() / len(bd)
            DT.append(dict(set="BE one per burst", metric="dep", n=int(len(bd)), mean_obs=float(bd.mean()),
                           mean_ETSdraw=float(mdb.mean()), p99=float(np.percentile(mdb, 99)),
                           p01=float(np.percentile(mdb, 1)), p_deeper=pct_rank(float(bd.mean()), mdb),
                           deeper=bool(bd.mean() > np.percentile(mdb, 99))))
        from scipy.stats import spearmanr
        ff = fab[fab.n_ETS >= 10]
        rho, prho = spearmanr(ff.BE_over_ETS, ff.dep); rho2, p2 = spearmanr(ff.BE_over_ETS, ff.slab2_dep)
        R["depth"] = dict(tests=DT, spearman_ratio_vs_depth=float(rho), p=float(prho), n_fam=int(len(ff)),
                          spearman_ratio_vs_slab2=float(rho2), p_slab2=float(p2),
                          ETS_mean_dep_AB=float(ets[ets.fid.isin(fab.fid)].fid.map(F.dep).mean()),
                          share_BE_from_C_families=float(np.mean(~ev.fid.isin(fab.fid))) if len(ev) else None,
                          share_ETS_from_C_families=float(np.mean(~ets.fid.isin(fab.fid))) if len(ets) else None,
                          share_BE_deeper40=float(np.mean(evab.d > 40)) if len(evab) else None,
                          share_ETS_deeper40=float(np.mean(ets[ets.fid.isin(fab.fid)].fid.map(F.dep) > 40)))
        log(f"[{era}] depth: {json.dumps(R['depth'], indent=1, default=str)}")
        # per-cut mean depth (data vs controls)
        cut_rows = []
        for P in (0.5, 0.75, 0.85, 0.9):
            for ser in ("data", "rev", "c1"):
                y = sel(era, ser, P)[lambda d: d.BE]; y = y[y.fid.isin(fab.fid)]
                cut_rows.append(dict(P=P, series=ser, n=len(y), mean_dep=float(y.fid.map(F.dep).mean()) if len(y) else np.nan))
        R["depth_by_cut"] = cut_rows
        # groups (source groups) + depth bins
        grp = fam.groupby("source_group").agg(n_fam=("fid", "size"), dep=("dep", "mean"), along_km=("along_km", "mean"),
                                              n_BE=("n_BE", "sum"), n_BE_rev=("n_BE_rev", "sum"), n_BE_c1=("n_BE_c1", "sum"),
                                              n_BE_broad=("n_BE_broad", "sum"), n_ETS=("n_ETS", "sum"),
                                              qual=("quality", lambda q: "".join(sorted(set(q))))).reset_index()
        grp["BE_share"] = grp.n_BE / max(grp.n_BE.sum(), 1); grp["ETS_share"] = grp.n_ETS / max(grp.n_ETS.sum(), 1)
        grp.sort_values("n_BE", ascending=False).to_csv(f"{OUT}/groups_between_{era}.csv", index=False)
        db = pd.cut(fab.dep, [0, 25, 30, 35, 40, 45, 60])
        dbt = fab.groupby(db, observed=False).agg(n_fam=("fid", "size"), n_BE=("n_BE", "sum"), n_ETS=("n_ETS", "sum"),
                                                  n_BE_broad=("n_BE_broad", "sum")).reset_index()
        dbt["BE_share"] = dbt.n_BE / max(dbt.n_BE.sum(), 1); dbt["ETS_share"] = dbt.n_ETS / max(dbt.n_ETS.sum(), 1)
        dbt["dep"] = dbt.dep.astype(str)
        dbt.to_csv(f"{OUT}/depth_bins_{era}.csv", index=False)
        R["depth_bins"] = dbt.round(3).to_dict("records")
        R["top_groups"] = grp.sort_values("n_BE", ascending=False).head(8).round(3).to_dict("records")
        R["n_families_with_BE"] = int((fam.n_BE > 0).sum()); R["n_groups_with_BE"] = int((grp.n_BE > 0).sum())
        R["top5_families_share"] = float(fam.n_BE.nlargest(5).sum() / max(len(ev), 1))
        R["top_families"] = fam.sort_values("n_BE", ascending=False).head(8)[["fid", "source_group", "quality", "dep", "along_km", "n_BE", "BE_days", "n_ETS"]].round(1).to_dict("records")
        R["events_per_year"] = ev.date_ts.dt.year.value_counts().sort_index().to_dict()
        log(f"[{era}] groups:\n" + grp.sort_values("n_BE", ascending=False).head(10).round(3).to_string(index=False))
        log(f"[{era}] depth bins:\n" + dbt.round(3).to_string(index=False))

        # ---- detectability + station changes
        R["n_core_live"] = dict(BE_events_median=float(ev.n_core_live.median()), eligible_median=float(EL.n_core_live.median()),
                                BE_events_mean=float(ev.n_core_live.mean()), eligible_mean=float(EL.n_core_live.mean()))
        era_res[era] = R
        if era == "E3":
            keepE3 = dict(ev=ev, evb=evb, evr=evr, evc=evc, EL=EL, ELb=ELb, el_by_f=el_by_f, fam=fam, ets=ets, A=A)
    EP = pd.concat(tabs_all); EP.to_csv(f"{OUT}/superposed_epoch.csv", index=False)
    RES["eras"] = era_res

    # =================================================================== E3-only tests
    ev, evb, evr, evc, EL, ELb, el_by_f, fam = (keepE3[k] for k in ("ev", "evb", "evr", "evc", "EL", "ELb", "el_by_f", "fam"))
    # bursts table
    B = ev.sort_values("time_utc").groupby("burst").agg(start_PT=("time_PT", "first"), end_PT=("time_PT", "last"),
                                                        n=("fid", "size"), n_fam=("fid", "nunique"),
                                                        n_groups=("source_group", "nunique"),
                                                        groups=("source_group", lambda g: ";".join(sorted(set(g)))),
                                                        along_km=("along_km", "median"), dep=("dep", "median"))
    B.sort_values("n", ascending=False).to_csv(f"{OUT}/bursts_between.csv", index=False)
    log("largest bursts (E3):\n" + B.sort_values("n", ascending=False).head(12).to_string(index=False))

    # ---- location vs next start
    K = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_kinematics.csv")
    K = K[K.region == "NW"][["episode", "lfe_start10", "trem_start10", "direction"]]
    fab = F[F.ab]
    lat0, lat1 = fab.lat.min() - .2, fab.lat.max() + .2; lon0, lon1 = fab.lon.min() - .2, fab.lon.max() + .2
    trb = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1)].copy()
    trb["along"] = along(trb.lat, trb.lon)
    tst = {}
    for r in maj.itertuples():
        x = trb[(trb.time >= r.start) & (trb.time < r.end + pd.Timedelta(days=1))].sort_values("time")
        if len(x) >= 20:
            k = max(1, int(np.ceil(0.1 * len(x)))); tst[r.episode] = float(np.median(x.along.to_numpy()[:k]))
    maj["trem_start_all"] = maj.episode.map(tst)
    chk = K.merge(maj[["episode", "trem_start_all"]], on="episode")
    log("tremor start check (S3 trem_start10 vs recomputed):\n" + chk.round(1).to_string(index=False))

    def loc_test(evx, starts, name, nperm=5000):
        e = evx[(evx.tau_next <= 60) & evx.next_major.isin(starts.keys()) & evx.along_km.notna()]
        if len(e) == 0:
            return dict(set=name, n=0)
        keys = list(starts.keys()); sv = np.array([starts[k] for k in keys])
        a = e.along_km.to_numpy(); own = np.array([keys.index(k) for k in e.next_major])
        obs = float(np.median(np.abs(a - sv[own])))
        perm = np.empty(nperm)
        for i in range(nperm):
            r = RNG.integers(0, len(keys) - 1, len(e)); r = np.where(r >= own, r + 1, r)
            perm[i] = np.median(np.abs(a - sv[r]))
        return dict(set=name, n=int(len(e)), n_episodes=int(e.next_major.nunique()), median_abs_km=obs,
                    perm_median_km=float(np.median(perm)), perm_p01_km=float(np.percentile(perm, 1)),
                    p=float((np.sum(perm <= obs) + 1) / (nperm + 1)),
                    per_episode=e.groupby("next_major").along_km.agg(["size", "median"]).round(1).reset_index().to_dict("records"))
    s3st = dict(zip(K.episode, K.lfe_start10))
    evab = ev[ev.ab]; evbab = evb[evb.ab]
    L = [loc_test(evab, s3st, "BE vs S3 LFE start (14 big eps)"),
         loc_test(evab, tst, "BE vs PNSN start (all majors)"),
         loc_test(evbab, s3st, "BE_broad vs S3 LFE start"),
         loc_test(evbab, tst, "BE_broad vs PNSN start (all majors)")]
    L[0]["pass"] = bool(L[0].get("n", 0) >= 20 and L[0].get("p", 1) < 0.01)
    L[0]["verdict"] = "too few to test" if L[0].get("n", 0) < 20 else ("pass" if L[0]["pass"] else "fail")
    RES["location_next_start"] = L
    log("location vs next start: " + json.dumps(L, indent=1, default=str))
    maj.to_csv(f"{OUT}/majors_tremor_start.csv", index=False)

    # ---- tremor elsewhere (30-100 km), precomputed per (fid, date)
    flat, flon = F.lat.to_dict(), F.lon.to_dict()
    fids = F.fid[F.lat.notna()].to_numpy()
    FL = np.c_[F.loc[fids, "lat"].to_numpy(), F.loc[fids, "lon"].to_numpy()]
    far_map = {}
    for d, g in tr.groupby("date"):
        dkm = hav_km(FL[:, 0:1], FL[:, 1:2], g.lat.to_numpy()[None], g.lon.to_numpy()[None])
        c = ((dkm >= 30) & (dkm <= 100)).sum(1)
        for f, v in zip(fids, c):
            if v:
                far_map[(int(f), d)] = int(v)
    far = lambda f, d: far_map.get((int(f), d), 0)
    tt = tr.time.astype("int64").to_numpy() / 3.6e12; order = np.argsort(tt); tt = tt[order]
    tla, tlo = tr.lat.to_numpy()[order], tr.lon.to_numpy()[order]

    def far6(fid, t):
        if not np.isfinite(flat.get(fid, np.nan)):
            return 0
        th = t / 3.6e12; i0, i1 = np.searchsorted(tt, [th - 6, th + 6])
        if i1 <= i0:
            return 0
        dkm = hav_km(flat[fid], flon[fid], tla[i0:i1], tlo[i0:i1]); return int(((dkm >= 30) & (dkm <= 100)).sum())
    for x in (ev, evr, evc, evb):
        x["n_trem_30_100_day"] = [far(f, d) for f, d in zip(x.fid, x.date)]
        x["n_trem_30_100_6h"] = [far6(f, t) for f, t in zip(x.fid, x.time_utc.astype("int64"))]
    evf = ev[ev.lat.notna()]
    obs_far = float(np.mean(evf.n_trem_30_100_day > 0))
    shd = shuffle_dates(evf.fid.to_numpy(), el_by_f, 1000)
    fv = evf.fid.to_numpy()
    shf = np.array([np.mean([far(f, d) > 0 for f, d in zip(fv, dr)]) for dr in shd])
    RES["tremor_30_100"] = dict(n=int(len(evf)), share_BE_with_far_tremor_day=obs_far, shuffle_median=float(np.median(shf)),
                                shuffle_p99=float(np.percentile(shf, 99)), p=pct_rank(obs_far, shf),
                                tied=bool(obs_far > np.percentile(shf, 99)),
                                share_BE_far_6h=float(np.mean(evf.n_trem_30_100_6h > 0)),
                                share_rev_far_day=float(np.mean(evr.n_trem_30_100_day > 0)) if len(evr) else None,
                                share_c1_far_day=float(np.mean(evc.n_trem_30_100_day > 0)) if len(evc) else None,
                                share_BE_broad_far_day=float(np.mean(evb.n_trem_30_100_day > 0)))
    log(f"tremor 30-100 km: {RES['tremor_30_100']}")

    # ---- earthquakes within +-60 s
    eq = pd.read_csv(f"{ROOT}/lfe_extension/results/trackC/eq_catalog.csv")
    eq["ot"] = pd.to_datetime(eq.ot, errors="coerce"); eq = eq.dropna(subset=["ot"]).sort_values("ot")
    teq = eq.ot.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9

    def eq_near(t_s):
        i = np.clip(np.searchsorted(teq, t_s), 1, teq.size - 1)
        j = np.where(np.abs(teq[i] - t_s) < np.abs(teq[i - 1] - t_s), i, i - 1)
        return np.abs(teq[j] - t_s), j
    for x in (ev, evr, evc, evb):
        dt_, j = eq_near(x.time_utc.astype("int64").to_numpy() / 1e9)
        x["eq_dt_s"] = dt_; x["eq_mag"] = eq.mag.to_numpy()[j]
        x["eq_dist_km"] = hav_km(x.fid.map(flat).to_numpy(), x.fid.map(flon).to_numpy(), eq.lat.to_numpy()[j], eq.lon.to_numpy()[j])
    rt = (pd.to_datetime(pd.Series(shd.ravel())).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
          + RNG.uniform(0, 86400, shd.size))
    RES["earthquakes_60s"] = dict(n_BE_within_60s=int((ev.eq_dt_s <= 60).sum()), share=float(np.mean(ev.eq_dt_s <= 60)),
                                  chance_share=float(np.mean(eq_near(rt)[0] <= 60)),
                                  n_rev_within_60s=int((evr.eq_dt_s <= 60).sum()), n_c1_within_60s=int((evc.eq_dt_s <= 60).sum()),
                                  share_broad=float(np.mean(evb.eq_dt_s <= 60)),
                                  events=ev[ev.eq_dt_s <= 60][["time_PT", "fid", "eq_dt_s", "eq_mag", "eq_dist_km"]].round(1).to_dict("records")[:20])
    log(f"earthquakes: {RES['earthquakes_60s']}")
    tt_s = np.sort(ev.time_utc.astype("int64").to_numpy() / 1e9)
    RES["n_BE_with_other_family_within_120s"] = int(np.sum((np.r_[np.inf, np.diff(tt_s)] <= 120) | (np.r_[np.diff(tt_s), np.inf] <= 120)))

    # ---- station changes
    C = pd.read_csv(f"{WP}/post_v2/core_set_changes_v2.csv")
    cd = np.sort(pd.to_datetime(C.date).values.astype("datetime64[D]"))

    def days_to_change(dates):
        d = pd.to_datetime(pd.Series(dates)).values.astype("datetime64[D]")
        i = np.clip(np.searchsorted(cd, d), 1, cd.size - 1)
        return np.minimum(np.abs((cd[i] - d).astype(int)), np.abs((d - cd[i - 1]).astype(int)))
    ev["days_to_core_change"] = days_to_change(ev.date)
    for x in (evr, evc, evb):
        x["days_to_core_change"] = days_to_change(x.date)
    eld = days_to_change(EL.date)
    RES["core_changes"] = dict(share_BE_within_7d=float(np.mean(ev.days_to_core_change <= 7)),
                               share_eligible_within_7d=float(np.mean(eld <= 7)),
                               share_BE_bursts_within_7d=float(np.mean(ev.groupby("burst").days_to_core_change.min() <= 7)))
    log(f"core changes: {RES['core_changes']}")

    # ---- GNSS stack
    RES["gnss"] = gnss_stack(ev, evb, E)
    log(f"gnss: {json.dumps(RES['gnss'], indent=1)}")

    # ---- added after the first run (columns only, no rule changed): burst-level versions of the location and
    #      tremor-elsewhere tests (WA BE events sit in a few large bursts), E1/E2 burst tables
    bu = ev[ev.ab].groupby("burst").agg(along_km=("along_km", "median"), next_major=("next_major", "first"),
                                        tau_next=("tau_next", "first"), fid=("fid", "first"), date=("date", "first"),
                                        n=("fid", "size"))
    bu["ab"] = True
    RES["location_next_start_burst_level"] = [loc_test(bu, s3st, "bursts vs S3 LFE start (14 big eps)"),
                                              loc_test(bu, tst, "bursts vs PNSN start (all majors)")]
    log("burst-level location: " + json.dumps(RES["location_next_start_burst_level"], default=str))
    bf = ev[ev.lat.notna()].groupby("burst").agg(fid=("fid", "first"), date=("date", "first"),
                                                 far=("n_trem_30_100_day", lambda v: bool((v > 0).any())))
    shb = shuffle_dates(bf.fid.to_numpy(), el_by_f, 1000)
    sfb = np.array([np.mean([far(f, d) > 0 for f, d in zip(bf.fid.to_numpy(), dr)]) for dr in shb])
    RES["tremor_30_100_burst_level"] = dict(n_bursts=int(len(bf)), share_far=float(bf.far.mean()),
                                            shuffle_median=float(np.median(sfb)), shuffle_p99=float(np.percentile(sfb, 99)),
                                            p=pct_rank(float(bf.far.mean()), sfb))
    log(f"burst-level tremor 30-100 km: {RES['tremor_30_100_burst_level']}")
    bt = []
    for era in ("E1", "E2"):
        y = sel(era)[lambda d: d.BE].copy()
        y["burst"] = chain_ids(y.time_utc.to_numpy())
        y["tau_next"], y["tau_prev"], _ = tau_next(y.date)
        bt.append(y.sort_values("time_utc").groupby("burst").agg(start_PT=("time_PT", "first"), end_PT=("time_PT", "last"),
                                                                 n=("fid", "size"), n_fam=("fid", "nunique"),
                                                                 groups=("source_group", lambda g: ";".join(sorted(set(g)))),
                                                                 along_km=("along_km", "median"), dep=("dep", "median"),
                                                                 tau_next=("tau_next", "first"), tau_prev=("tau_prev", "first")).assign(era=era))
    BT = pd.concat(bt).sort_values("n", ascending=False); BT.to_csv(f"{OUT}/bursts_between_E1E2.csv", index=False)
    log("largest bursts E1/E2:\n" + BT.head(10).to_string(index=False))

    # ---- save
    keep =["time_utc", "time_PT", "date", "fid", "source_group", "array_group", "quality", "dep", "slab2_dep", "along_km",
            "downdip_km", "P_own_WAVE_PNET", "dup_D", "once_group", "n_core_live", "n_trem30", "n_trem60", "tau_next",
            "tau_prev", "next_major", "burst", "n_trem_30_100_day", "n_trem_30_100_6h", "eq_dt_s", "eq_mag", "eq_dist_km",
            "days_to_core_change", "series"]
    for x in (ev, evr, evc, evb):
        x["array_group"] = x.fid.map(F.array_group)
    ev[keep].sort_values("time_utc").to_csv(f"{OUT}/between_events.csv", index=False)
    pd.concat([evr[keep], evc[keep]]).sort_values("time_utc").to_csv(f"{OUT}/between_events_control.csv", index=False)
    evb[keep].sort_values("time_utc").to_parquet(f"{OUT}/between_broad_events.parquet")
    # E1/E2 BE events too (secondary)
    sec = []
    for era in ("E1", "E2"):
        y = sel(era)[lambda d: d.BE]
        sec.append(y[["time_utc", "time_PT", "date", "fid", "source_group", "quality", "dep", "along_km", "era"]])
    pd.concat(sec).sort_values("time_utc").to_csv(f"{OUT}/between_events_E1E2.csv", index=False)
    json.dump(RES, open(f"{OUT}/a_summary.json", "w"), indent=1, default=str)
    figures(F, E, maj, T, sel, EP, RES)
    log("=== A done ===")


# ------------------------------------------------------------------------------------------- GNSS (S7 estimator)
def gnss_stack(ev, evb, E):
    sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    from sse_strain.matio import load_mat
    import slip_build as SB
    Xd = load_mat(f"{SB.SOLDIR}/X.mat")["Xd"]
    dr = SB.fixed_dates(np.asarray(Xd["timeline"], float).ravel())
    names = np.asarray(Xd["name"]).astype(str); ts = np.asarray(Xd["ts"], float); m = np.asarray(Xd["data_mask"]).astype(bool)
    llh = np.asarray(Xd["llh"], float)[::3]
    cal = np.arange(dr[0], dr[-1] + 1); ix = (dr - cal[0]).astype(int)
    Y = np.full((ts.shape[0], cal.size), np.nan); Y[:, ix] = np.where(m, ts, np.nan)
    sta = np.array([n[:-1] for n in names[::3]])
    W = np.lib.stride_tricks.sliding_window_view(np.concatenate([Y, np.full((Y.shape[0], 4), np.nan)], 1), 5, axis=1)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); M5 = np.nanmedian(W, axis=-1)
    M5[np.isfinite(W).sum(-1) < 3] = np.nan
    s1 = json.load(open(f"{ROOT}/svi_science/results/s1/tierA_summary.json"))
    az = s1["wa_azimuths"]
    rows = [int(np.flatnonzero(sta == s)[0]) for s in WA12]
    st_as = along(llh[rows, 1], llh[rows, 0])
    north = [s for s, a in zip(WA12, st_as) if a > 0]
    U = np.array([[np.sin(np.deg2rad(az[s])), np.cos(np.deg2rad(az[s]))] for s in WA12])
    re, rn = 3 * np.array(rows), 3 * np.array(rows) + 1
    n = cal.size; d = np.arange(5, n - 1)
    dE = M5[re][:, d + 1] - M5[re][:, d - 5]; dN = M5[rn][:, d + 1] - M5[rn][:, d - 5]
    proj = dE * U[:, 0:1] + dN * U[:, 1:2]
    off = np.full(n, np.nan); off[d] = np.nanmean(proj, axis=0)
    nmask = np.array([s in north for s in WA12])
    offn = np.full(n, np.nan); offn[d] = np.nanmean(proj[nmask], axis=0)
    bad = np.zeros(n, bool)
    for s, e in zip(E.start.values.astype("datetime64[D]"), E.end.values.astype("datetime64[D]")):
        a = int((s - cal[0]).astype(int)) - 2 - 6; b = int((e - cal[0]).astype(int)) + 2 + 6
        bad[max(a, 0):min(b, n - 1) + 1] = True
    lo = int((np.datetime64("2015-01-15") - cal[0]).astype(int)); hi = int((np.datetime64(P_END.date(), "D") - cal[0]).astype(int)) - 6
    allowed = np.zeros(n, bool); allowed[lo:hi] = True; allowed &= ~bad & np.isfinite(off)
    doy = (cal - cal.astype("datetime64[Y]")).astype(int)
    alw = np.flatnonzero(allowed)

    def clusters(evx):
        c = evx.groupby("date").size().sort_index()
        dd = pd.to_datetime(c.index).values.astype("datetime64[D]")
        gid = np.r_[0, np.cumsum(np.diff(dd).astype(int) >= 3)]
        out = []
        for g in np.unique(gid):
            k = np.flatnonzero(gid == g); best = k[np.argmax(c.values[k])]
            out.append(int((dd[best] - cal[0]).astype(int)))
        return np.array(out)

    def stack(idx, name, nset=1000):
        idx = idx[allowed[idx]]
        if idx.size == 0:
            return dict(set=name, n_days=0)
        obs = float(np.nanmean(off[idx])); obsn = float(np.nanmean(offn[idx]))
        Cc = np.empty(nset); Cn = np.empty(nset)
        cand = []
        for i in idx:
            dd = np.abs(doy[alw] - doy[i]); dd = np.minimum(dd, 365 - dd)
            cand.append(alw[dd <= 30])
        for j in range(nset):
            pick = np.array([c[RNG.integers(0, c.size)] for c in cand])
            Cc[j] = np.nanmean(off[pick]); Cn[j] = np.nanmean(offn[pick])
        sig = float(Cc.std()); mu = float(Cc.mean())
        return dict(set=name, n_days=int(idx.size), mean_mm=obs, ctrl_mean_mm=mu, ctrl_sigma_mm=sig, z=(obs - mu) / sig,
                    bound_3sigma_mm=3 * sig, north_mean_mm=obsn, north_ctrl_mean=float(Cn.mean()),
                    north_sigma=float(Cn.std()), north_z=float((obsn - Cn.mean()) / Cn.std()),
                    slip_like=bool((obs - mu) / sig > 3), ctrl_values=Cc.round(4).tolist())
    cb = clusters(ev); cbb = clusters(evb)
    mins = E[~E.major & (E.start >= "2015-01-01") & (E.end <= P_END)]
    mid = np.array([int(((s + (e - s) / 2).to_datetime64().astype("datetime64[D]") - cal[0]).astype(int)) for s, e in zip(mins.start, mins.end)])
    mid = mid[(mid > 6) & (mid < n - 7)]
    res = [stack(cb, "BE cluster days"), stack(cbb, "BE_broad cluster days")]
    allowed_bak = allowed.copy(); allowed[mid] = np.isfinite(off[mid])
    res.append(stack(mid, "minor-episode mid-days (reference)")); allowed[:] = allowed_bak
    res[0]["n_clusters_total"] = int(cb.size); res[1]["n_clusters_total"] = int(cbb.size)
    res[0]["north_stations"] = north
    np.savez(f"{OUT}/gnss_stack.npz", **{r["set"].split()[0].replace("-", "_"): np.array(r.get("ctrl_values", [])) for r in res},
             obs=np.array([r.get("mean_mm", np.nan) for r in res]))
    for r in res:
        r.pop("ctrl_values", None)
    return res


# ------------------------------------------------------------------------------------------- figures
def figures(F, E, maj, T, sel, EP, RES):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6})
    MUTED = "#777777"
    toPT = lambda s: s.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.tz_localize(None)
    ev = pd.concat([sel(e)[lambda d: d.BE] for e in ("E1", "E2", "E3")])
    evb = pd.concat([sel(e)[lambda d: d.BE_broad] for e in ("E1", "E2", "E3")])
    evr = pd.concat([sel(e, "rev")[lambda d: d.BE] for e in ("E1", "E2", "E3")])
    evc = pd.concat([sel(e, "c1")[lambda d: d.BE] for e in ("E1", "E2", "E3")])

    # ---- 1. timeline
    fig, ax = plt.subplots(2, 1, figsize=(13, 6.8), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]), layout="constrained")
    a = ax[0]
    E2 = E[E.end >= "2009-08-01"]
    for r in E2.itertuples():
        a.axvspan(np.datetime64(r.start), np.datetime64(r.end + pd.Timedelta(days=1)),
                  color="#b8cbe4" if r.major else "#e3ebf5", lw=0, zorder=0)
    a.scatter(toPT(evb.time_utc).to_numpy(), evb.along_km, s=4, c="#b0b0b0", lw=0, zorder=1,
              label=f"outside episodes ±2 d but tremor ≤ 30 km ({len(evb) - len(ev)})")
    ab = ev[ev.ab]; nc = ev[~ev.ab]
    sc = a.scatter(toPT(ab.time_utc).to_numpy(), ab.along_km, s=14, c=ab.dep, cmap="viridis_r", vmin=25, vmax=45, lw=0.3,
                   edgecolor="white", label=f"between-episode, no tremor ≤ 30 km ({len(ev)}; colour = depth)", zorder=3)
    a.scatter(toPT(nc.time_utc).to_numpy(), nc.along_km, s=14, facecolor="none", edgecolor="#555555", lw=0.6, zorder=3,
              label=f"  … from quality-C families, depth unknown ({len(nc)})")
    a.scatter(toPT(evr.time_utc).to_numpy(), evr.along_km, s=40, marker="x", c="#d62728", lw=1.2, zorder=4,
              label=f"time-reversed control ({len(evr)})")
    a.scatter(toPT(evc.time_utc).to_numpy(), evc.along_km, s=40, marker="+", c="#ff7f0e", lw=1.2, zorder=4,
              label=f"scrambled control ({len(evc)})")
    for d, t in (("2012-01-01", "E2"), ("2015-01-01", "E3 (primary)")):
        a.axvline(np.datetime64(d), color="#444444", lw=0.8, ls="--"); a.text(np.datetime64(d), 62, " " + t, fontsize=8)
    a.text(np.datetime64("2009-09-01"), 62, "E1", fontsize=8)
    a.set_ylabel("along strike (km, + = NW)"); a.set_ylim(-75, 70)
    a.legend(loc="lower left", fontsize=7.5, frameon=True, framealpha=0.9, ncol=3)
    a.set_title("Washington-NW LFEs between ETS episodes, 2009-08 → 2026-05, counts = all eras (P ≥ 0.75, counted once). Shaded = the 118 NW episodes (darker = 69 majors)",
                loc="left", fontsize=10)
    fig.colorbar(sc, ax=a, pad=0.01, shrink=0.9).set_label("family depth (km, quality A/B)")
    b = ax[1]
    mo = toPT(ev.time_utc).dt.to_period("M").value_counts().sort_index()
    mob = toPT(evb.time_utc).dt.to_period("M").value_counts().sort_index()
    b.bar(mob.index.to_timestamp().to_numpy(), mob.values, width=np.timedelta64(25, "D"), color="#c8c8c8", label="all outside episodes ±2 d")
    b.bar(mo.index.to_timestamp().to_numpy(), mo.values, width=np.timedelta64(25, "D"), color="#2c7fb8", label="no tremor ≤ 30 km")
    b.set_ylabel("events / month"); b.set_yscale("symlog", linthresh=5); b.legend(fontsize=8, frameon=False, loc="upper left")
    b.set_xlabel("time (Pacific)"); b.set_xlim(np.datetime64("2009-08-01"), np.datetime64("2026-06-15"))
    fig.colorbar(sc, ax=b, pad=0.01).ax.set_visible(False)
    fig.savefig(f"{FIG}/a_timeline.png", dpi=150); plt.close(fig)

    # ---- 2. superposed epoch around ETS ends (+ before next start) + SVI comparison
    svi = pd.read_csv(f"{ROOT}/svi_science/results/s7/superposed_epoch.csv")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    for a, (col, ttl) in zip(ax[:2], (("tau_prev", "days AFTER the previous major ETS ended"),
                                      ("tau_next", "days BEFORE the next major ETS starts"))):
        t = EP[(EP.era == "E3") & (EP.set == "BE") & (EP.axis == col)]
        x = np.arange(len(t))
        a.fill_between(x, t.shuffle_p01, t.shuffle_p99, color="#dddddd", step="mid", label="family-day shuffle 1–99 %")
        a.step(x, t.shuffle_mean, where="mid", color=MUTED, lw=1, label="shuffle mean")
        a.plot(x, t.events, "o-", color="#08519c", lw=2, ms=5, label="observed (WA E3)")
        a.set_xticks(x); a.set_xticklabels(t.bin, fontsize=8); a.set_xlabel(ttl)
        a.set_ylabel("between-episode events (no tremor ≤ 30 km)"); a.legend(fontsize=8, frameon=False)
    c = RES["eras"]["E3"]["next_ets"]["BE_tau_prev"]
    ax[0].set_title(f"After an ETS ends: first 30 d hold {c['F30']:.2f} vs {c['shuffle_F30_median']:.2f} expected "
                    f"(p_low {c['p_low']:.3f})", loc="left", fontsize=9)
    c2 = RES["eras"]["E3"]["next_ets"]["BE_tau_next"]
    ax[1].set_title(f"Before the next ETS: last 30 d hold {c2['F30']:.2f} vs {c2['shuffle_F30_median']:.2f} expected",
                    loc="left", fontsize=9)
    a = ax[2]
    t = EP[(EP.era == "E3") & (EP.set == "BE") & (EP.axis == "tau_prev")]
    s = svi[(svi.set == "BE") & (svi.axis == "tau_prev")]
    x = np.arange(len(t))
    a.axhline(1, color=MUTED, lw=1)
    a.plot(x, t.events / t.shuffle_mean, "o-", color="#08519c", lw=2, label="Washington-NW (E3)")
    a.plot(x, s.events.to_numpy() / s.shuffle_mean.to_numpy(), "s--", color="#e6550d", lw=1.6, label="S. Vancouver Island (S7)")
    for era, mk in (("E2", "^"), ("E1", "v")):
        te = EP[(EP.era == era) & (EP.set == "BE") & (EP.axis == "tau_prev")]
        a.plot(x, te.events / te.shuffle_mean.replace(0, np.nan), mk + ":", color="#6baed6", lw=1, ms=4, label=f"WA {era}")
    a.set_xticks(x); a.set_xticklabels(t.bin, fontsize=8); a.set_xlabel("days AFTER the previous major ETS ended")
    a.set_ylabel("observed / expected"); a.set_yscale("log"); a.legend(fontsize=8, frameon=False)
    a.set_title("Both regions side by side", loc="left", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/a_superposed_epoch.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
