#!/usr/bin/env python
"""s7_between (2026-10-07) -- SVI science study S7: LFEs between ETS episodes.
Rules frozen in svi_science/results/s7/PLAN_S7.md before any number.  CPU only, < 10 GB.

    python svi_science/scripts/s7_between.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s7"
FIG = f"{ROOT}/svi_science/figures"
S0 = f"{ROOT}/svi_science/results/s0"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
FRAME = json.load(open(f"{S0}/frame_groups.json"))
P_END = pd.Timestamp("2026-05-27")          # PNSN catalogue end
P_START = pd.Timestamp("2014-01-01")
NSH = 1000
RNG = np.random.default_rng(20261007)
SVI14 = ["ALBH", "PGC5", "SC02", "SC04", "NEAH", "PTRF", "PTAA", "NANO", "P435", "P436", "P439", "CLRS", "COUP", "SC03"]
LOG = open(f"{OUT}/s7.log", "w")
RES = {}


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def local_km(lat, lon):
    lat0, lon0 = FRAME["origin_lat"], FRAME["origin_lon"]
    x = (np.asarray(lon) - lon0) * 111.32 * np.cos(np.deg2rad(lat0))
    y = (np.asarray(lat) - lat0) * 110.574
    return x, y


def along(lat, lon):
    x, y = local_km(lat, lon); s = np.deg2rad(FRAME["strike_az"])
    return x * np.sin(s) + y * np.cos(s)


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def pct_rank(obs, ctrl):
    ctrl = np.asarray(ctrl); return float((np.sum(ctrl >= obs) + 1) / (ctrl.size + 1))


# ------------------------------------------------------------------------------------------- load
def load():
    F = pd.read_csv(f"{S0}/families.csv")
    asx = along(F.lat, F.lon)
    log(f"frame check: max |along(lat,lon) - families.along_strike_km| = {np.max(np.abs(asx - F.along_strike_km)):.2f} km")
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    cols = ["date", "series", "fam", "P_same", "keep_xfam8", "quiet", "n_tremor_day", "n_trem30", "in_episode", "episode",
            "major", "time_utc", "n_sta_day", "flag_ep_2017_03", "flag_ep_2025_10"]
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=cols)
    T13 = pd.read_parquet(f"{S0}/analysis_table_S13.parquet", columns=cols)
    D = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable", "n_sta_day"])
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["time"])
    tr["time"] = tr.time.astype("datetime64[ns]")
    tr = tr[(tr.time >= P_START - pd.Timedelta(days=30)) & tr.lat.between(46.8, 50.3) & tr.lon.between(-126.5, -121.5)]
    tr["date"] = tr.time.dt.strftime("%Y-%m-%d")
    return F, E, T, T13, D, X, tr.reset_index(drop=True)


def pad_mask(dates, E, pad=2):
    d = pd.to_datetime(pd.Series(dates)).values
    m = np.zeros(len(d), bool)
    for s, e in zip(E.start.values, E.end.values):
        m |= (d >= s - np.timedelta64(pad, "D")) & (d <= e + np.timedelta64(pad, "D"))
    return m


def classify(T, E):
    T = T.copy()
    T["date_ts"] = pd.to_datetime(T.date)
    T["in_pad"] = pad_mask(T.date, E, 2)
    T["in_ep_nopad"] = pad_mask(T.date, E, 0)
    T["BE"] = ~T.in_pad & (T.n_trem30 == 0)
    T["BE_broad"] = ~T.in_pad
    T["out_nopad"] = ~T.in_ep_nopad
    T["box_quiet_day"] = T.n_tremor_day == 0
    return T


# ------------------------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    log("=== S7 start ===")
    F, E, T, T13, D, X, tr = load()
    F = F.set_index("fid", drop=False)
    maj = E[E.major].reset_index(drop=True)
    maj_start = maj.start.values; maj_end = maj.end.values

    T = T[(pd.to_datetime(T.date) >= P_START)]
    n_after = int(((pd.to_datetime(T.date) > P_END) & (T.series == "data") & (T.P_same >= .75) & T.keep_xfam8).sum())
    T = classify(T[pd.to_datetime(T.date) <= P_END], E)
    T13 = classify(T13[pd.to_datetime(T13.date) <= P_END], E)
    log(f"P>=0.75 de-dup data events after PNSN end (2026-05-28..09-15, not used): {n_after}")

    # ---------------- inventory table of definitions
    inv = []
    for run, TT in (("prod", T), ("prod_S13", T13)):
        for P in (0.75, 0.5):
            for ser in ("data", "c1"):
                x = TT[(TT.series == ser) & (TT.P_same >= P) & TT.keep_xfam8]
                inv.append(dict(run=run, P=P, series=ser, all=len(x), BE=int(x.BE.sum()), BE_broad=int(x.BE_broad.sum()),
                                box_no_tremor_day=int(x.box_quiet_day.sum()), quiet_flag=int(x.quiet.sum()),
                                outside_episodes_nopad=int(x.out_nopad.sum()),
                                BE_raw_no_dedup=int(TT[(TT.series == ser) & (TT.P_same >= P) & TT.BE].shape[0])))
    INV = pd.DataFrame(inv); INV.to_csv(f"{OUT}/inventory_definitions.csv", index=False)
    log("definitions:\n" + INV.to_string(index=False))

    sel = lambda TT, ser="data", P=0.75: TT[(TT.series == ser) & (TT.P_same >= P) & TT.keep_xfam8]
    A = sel(T); Ac = sel(T, "c1")
    ev = A[A.BE].copy(); evc = Ac[Ac.BE].copy(); evb = A[A.BE_broad].copy()
    ets = A[A.major & A.in_episode]
    for x in (ev, evc, evb):
        x["along_km"] = F.loc[x.fam, "along_strike_km"].values
        x["downdip_km"] = F.loc[x.fam, "downdip_km"].values
        x["dep"] = F.loc[x.fam, "dep"].values
        x["group"] = F.loc[x.fam, "group"].values
        x["time_PT"] = x.time_utc.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")
    log(f"BE events data {len(ev)}, control {len(evc)}; BE_broad data {len(evb)}; ETS (39 majors) events {len(ets)}")

    # ---------------- eligible family-days
    Dp = D[(D.run == "prod") & D.observable].copy()
    Dp = Dp[(pd.to_datetime(Dp.date) >= P_START) & (pd.to_datetime(Dp.date) <= P_END)]
    Dp = Dp.merge(X.rename(columns={"fid": "fid", "n30": "n30"}), on=["fid", "date"], how="left")
    Dp["n30"] = Dp.n30.fillna(0)
    Dp["in_pad"] = pad_mask(Dp.date, E, 2)
    EL = Dp[~Dp.in_pad & (Dp.n30 == 0)][["fid", "date", "n_sta_day"]].reset_index(drop=True)
    ELb = Dp[~Dp.in_pad][["fid", "date", "n_sta_day"]].reset_index(drop=True)
    log(f"eligible family-days BE {len(EL)} (of {len(Dp)} observable); BE_broad {len(ELb)}")
    el_by_f = {f: g.date.to_numpy() for f, g in EL.groupby("fid")}
    elb_by_f = {f: g.date.to_numpy() for f, g in ELb.groupby("fid")}
    miss = set(ev.fam) - set(el_by_f)
    assert not miss, miss

    def shuffle_dates(fams, pool, n=NSH):
        out = np.empty((n, len(fams)), dtype=object)
        for j, f in enumerate(fams):
            p = pool[f]; out[:, j] = p[RNG.integers(0, p.size, n)]
        return out

    # ---------------- per-family / group table
    fam = pd.DataFrame(dict(fid=F.fid, group=F.group, dep=F.dep, along_km=F.along_strike_km, downdip_km=F.downdip_km,
                            lat=F.lat, lon=F.lon, n_carriers=F.n_carriers))
    fam["n_BE"] = fam.fid.map(ev.fam.value_counts()).fillna(0).astype(int)
    fam["n_BE_c1"] = fam.fid.map(evc.fam.value_counts()).fillna(0).astype(int)
    fam["n_BE_broad"] = fam.fid.map(evb.fam.value_counts()).fillna(0).astype(int)
    fam["n_ETS"] = fam.fid.map(ets.fam.value_counts()).fillna(0).astype(int)
    fam["n_all"] = fam.fid.map(A.fam.value_counts()).fillna(0).astype(int)
    fam["elig_days_BE"] = fam.fid.map(EL.fid.value_counts()).fillna(0).astype(int)
    fam["BE_days"] = fam.fid.map(ev.groupby("fam").date.nunique()).fillna(0).astype(int)
    fam["BE_over_ETS"] = (fam.n_BE + 0.5) / (fam.n_ETS + 0.5)
    fam.to_csv(f"{OUT}/families_between.csv", index=False)
    grp = fam.groupby("group").agg(n_fam=("fid", "size"), dep=("dep", "mean"), along_km=("along_km", "mean"),
                                   downdip_km=("downdip_km", "mean"), n_BE=("n_BE", "sum"), n_BE_c1=("n_BE_c1", "sum"),
                                   n_BE_broad=("n_BE_broad", "sum"), n_ETS=("n_ETS", "sum")).reset_index()
    grp["BE_share"] = grp.n_BE / grp.n_BE.sum(); grp["ETS_share"] = grp.n_ETS / grp.n_ETS.sum()
    grp["broad_share"] = grp.n_BE_broad / grp.n_BE_broad.sum()
    grp.to_csv(f"{OUT}/groups_between.csv", index=False)
    log("groups:\n" + grp.round(3).to_string(index=False))
    top = fam.sort_values("n_BE", ascending=False).head(12)
    log("top BE families:\n" + top[["fid", "group", "dep", "along_km", "downdip_km", "n_BE", "BE_days", "n_ETS", "n_BE_broad"]].to_string(index=False))
    RES["n_BE"] = int(len(ev)); RES["n_BE_c1"] = int(len(evc)); RES["n_BE_broad"] = int(len(evb))
    RES["n_BE_broad_c1"] = int(Ac.BE_broad.sum())
    RES["n_families_with_BE"] = int((fam.n_BE > 0).sum())
    RES["top5_families_share"] = float(fam.n_BE.nlargest(5).sum() / max(len(ev), 1))
    RES["events_per_year"] = ev.date_ts.dt.year.value_counts().sort_index().to_dict()

    # ---------------- time clustering
    def chain(times, gap_h=24):
        t = np.sort(times.astype("datetime64[ns]").astype("int64")) / 3.6e12
        b = np.r_[0, np.cumsum(np.diff(t) >= gap_h)]
        return pd.Series(b).value_counts().to_numpy()
    sizes = chain(ev.time_utc.to_numpy())
    obs_days = ev.date.nunique()
    sh = shuffle_dates(ev.fam.to_numpy(), el_by_f)
    sh_days = np.array([len(set(r)) for r in sh])
    RES["cluster"] = dict(distinct_days=int(obs_days), shuffle_days_median=float(np.median(sh_days)),
                          shuffle_days_p01=float(np.percentile(sh_days, 1)),
                          p_fewer_days=float((np.sum(sh_days <= obs_days) + 1) / (NSH + 1)),
                          n_bursts=int(sizes.size), n_bursts_ge3=int((sizes >= 3).sum()),
                          share_in_bursts_ge3=float(sizes[sizes >= 3].sum() / sizes.sum()),
                          largest_bursts=sorted(sizes.tolist(), reverse=True)[:8])
    RES["cluster"]["clustered"] = bool(obs_days < np.percentile(sh_days, 1))
    # bursts table
    evs = ev.sort_values("time_utc").reset_index(drop=True)
    t_h = evs.time_utc.astype("int64").to_numpy() / 3.6e12
    evs["burst"] = np.r_[0, np.cumsum(np.diff(t_h) >= 24)]
    B = evs.groupby("burst").agg(start_PT=("time_PT", "first"), end_PT=("time_PT", "last"), n=("fam", "size"),
                                 n_fam=("fam", "nunique"), groups=("group", lambda g: ";".join(map(str, sorted(set(g))))),
                                 along_km=("along_km", "median"), dep=("dep", "median"))
    B.sort_values("n", ascending=False).to_csv(f"{OUT}/bursts_between.csv", index=False)
    log(f"clustering: {RES['cluster']}")
    log("largest bursts:\n" + B.sort_values("n", ascending=False).head(10).to_string(index=False))

    # ---------------- next / previous major
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
    for x in (ev, evc, evb):
        x["tau_next"], x["tau_prev"], x["next_major_idx"] = tau_next(x.date)
        x["next_major"] = [maj.episode.iloc[k] if k >= 0 else "" for k in x.next_major_idx]
    bins = [0, 7, 14, 30, 60, 120, 240, 2000]
    lab = ["0-7", "7-14", "14-30", "30-60", "60-120", "120-240", ">240"]
    ELt = EL.copy(); ELt["tau_next"], ELt["tau_prev"], _ = tau_next(ELt.date)
    ELbt = ELb.copy(); ELbt["tau_next"], ELbt["tau_prev"], _ = tau_next(ELbt.date)

    def epoch_table(evx, ELx, pool, name, col="tau_next"):
        e = evx[np.isfinite(evx[col])]
        expo = np.histogram(ELx[col].dropna(), bins)[0]
        obs = np.histogram(e[col], bins)[0]
        shd = shuffle_dates(e.fam.to_numpy(), pool)
        dd = pd.Series(shd.ravel()).unique()
        tmap = dict(zip(dd, tau_next(dd)[0 if col == "tau_next" else 1]))
        sht = np.vectorize(tmap.get)(shd).astype(float)
        shh = np.array([np.histogram(r[np.isfinite(r)], bins)[0] for r in sht])
        tab = pd.DataFrame(dict(set=name, axis=col, bin=lab, events=obs, elig_famdays=expo,
                                rate_per_1000=1000 * obs / np.maximum(expo, 1),
                                shuffle_mean=shh.mean(0), shuffle_p01=np.percentile(shh, 1, 0),
                                shuffle_p99=np.percentile(shh, 99, 0)))
        f30 = float(np.mean(e[col] < 30)); f30s = np.mean(sht < 30, axis=1)
        return tab, f30, f30s, sht

    tabs = []
    out2 = {}
    for name, evx, ELx, pool in (("BE", ev, ELt, el_by_f), ("BE_broad", evb, ELbt, elb_by_f)):
        for col in ("tau_next", "tau_prev"):
            tab, f30, f30s, sht = epoch_table(evx, ELx, pool, name, col)
            tabs.append(tab)
            out2[f"{name}_{col}"] = dict(n=int(np.isfinite(evx[col]).sum()), F30=f30, shuffle_F30_median=float(np.median(f30s)),
                                         shuffle_F30_p99=float(np.percentile(f30s, 99)), p=pct_rank(f30, f30s),
                                         p_low=float((np.sum(f30s <= f30) + 1) / (f30s.size + 1)))
        # pre-ETS excess not explained by the tail: events >= 30 d after previous end
        e2 = evx[(evx.tau_prev >= 30) & np.isfinite(evx.tau_next)]
        ELx2 = ELx[ELx.tau_prev >= 30]
        pool2 = {f: g.date.to_numpy() for f, g in ELx2.groupby("fid")}
        e2 = e2[e2.fam.isin(pool2.keys())]
        tab, f30, f30s, _ = epoch_table(e2, ELx2, pool2, name + "_prev>=30", "tau_next")
        tabs.append(tab)
        out2[f"{name}_tau_next_prev>=30"] = dict(n=int(len(e2)), F30=f30, shuffle_F30_median=float(np.median(f30s)),
                                                 shuffle_F30_p99=float(np.percentile(f30s, 99)), p=pct_rank(f30, f30s),
                                                 p_low=float((np.sum(f30s <= f30) + 1) / (f30s.size + 1)))
    EP = pd.concat(tabs); EP.to_csv(f"{OUT}/superposed_epoch.csv", index=False)
    RES["next_ets"] = out2
    RES["ramp_pass"] = bool(out2["BE_tau_next"]["F30"] > out2["BE_tau_next"]["shuffle_F30_p99"]
                            and out2["BE_tau_next_prev>=30"]["F30"] > out2["BE_tau_next_prev>=30"]["shuffle_F30_p99"])
    log("superposed epoch:\n" + EP.round(2).to_string(index=False))
    log(f"next-ETS stats: {json.dumps(out2, indent=1)}  ramp_pass={RES['ramp_pass']}")
    RES["control_tau_next"] = evc.tau_next.tolist()

    # ---------------- location vs next start (S3 starts; tremor starts for all majors)
    K = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_kinematics.csv")
    K = K[K.region == "SVI"][["episode", "lfe_start10", "trem_start10", "direction"]]
    lat0, lat1 = F.lat.min() - .2, F.lat.max() + .2; lon0, lon1 = F.lon.min() - .2, F.lon.max() + .2
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
        e = evx[(evx.tau_next <= 60) & evx.next_major.isin(starts.keys())]
        if len(e) == 0:
            return dict(set=name, n=0)
        keys = list(starts.keys()); sv = np.array([starts[k] for k in keys])
        a = e.along_km.to_numpy(); own = np.array([keys.index(k) for k in e.next_major])
        obs = float(np.median(np.abs(a - sv[own])))
        perm = np.empty(nperm)
        for i in range(nperm):
            r = RNG.integers(0, len(keys) - 1, len(e)); r = np.where(r >= own, r + 1, r)
            perm[i] = np.median(np.abs(a - sv[r]))
        p = float((np.sum(perm <= obs) + 1) / (nperm + 1))
        return dict(set=name, n=int(len(e)), n_episodes=int(e.next_major.nunique()), median_abs_km=obs,
                    perm_median_km=float(np.median(perm)), perm_p01_km=float(np.percentile(perm, 1)), p=p,
                    per_episode=e.groupby("next_major").along_km.agg(["size", "median"]).round(1).reset_index().to_dict("records"))
    s3st = dict(zip(K.episode, K.lfe_start10))
    L = [loc_test(ev, s3st, "BE vs S3 LFE start (10 eps)"),
         loc_test(ev, {k: v for k, v in tst.items()}, "BE vs PNSN start (all majors)"),
         loc_test(evb, s3st, "BE_broad vs S3 LFE start"),
         loc_test(evb, {k: v for k, v in tst.items()}, "BE_broad vs PNSN start (all majors)")]
    for l in L:
        if l.get("n", 0) and l["set"].startswith("BE vs S3"):
            l["pass"] = bool(l["n"] >= 20 and l["p"] < 0.01)
    RES["location_next_start"] = L
    log("location vs next start: " + json.dumps(L, indent=1, default=str))
    maj.to_csv(f"{OUT}/majors_tremor_start.csv", index=False)

    # ---------------- tremor elsewhere (30-100 km)
    tr_by_day = {d: g for d, g in tr.groupby("date")}
    flat, flon = F.lat.to_dict(), F.lon.to_dict()

    def has_far(fid, date):
        g = tr_by_day.get(date)
        if g is None:
            return 0
        dkm = hav_km(flat[fid], flon[fid], g.lat.to_numpy(), g.lon.to_numpy())
        return int(((dkm >= 30) & (dkm <= 100)).sum())
    cache = {}
    def far(fid, date):
        k = (fid, date)
        if k not in cache:
            cache[k] = has_far(fid, date)
        return cache[k]
    for x in (ev, evc, evb):
        x["n_trem_30_100_day"] = [far(f, d) for f, d in zip(x.fam, x.date)]
    # +-6 h
    tt = tr.time.astype("int64").to_numpy() / 3.6e12; order = np.argsort(tt); tt = tt[order]
    tla, tlo = tr.lat.to_numpy()[order], tr.lon.to_numpy()[order]

    def far6(fid, t):
        th = t / 3.6e12; i0, i1 = np.searchsorted(tt, [th - 6, th + 6])
        if i1 <= i0:
            return 0
        dkm = hav_km(flat[fid], flon[fid], tla[i0:i1], tlo[i0:i1]); return int(((dkm >= 30) & (dkm <= 100)).sum())
    for x in (ev, evc, evb):
        x["n_trem_30_100_6h"] = [far6(f, t) for f, t in zip(x.fam, x.time_utc.astype("int64"))]
    obs_far = float(np.mean(ev.n_trem_30_100_day > 0))
    shd = shuffle_dates(ev.fam.to_numpy(), el_by_f, 1000)
    fams = np.tile(ev.fam.to_numpy(), (1000, 1))
    shf = np.array([[far(f, d) > 0 for f, d in zip(fr, dr)] for fr, dr in zip(fams, shd)]).mean(1)
    RES["tremor_30_100"] = dict(share_BE_with_far_tremor_day=obs_far, shuffle_median=float(np.median(shf)),
                                shuffle_p99=float(np.percentile(shf, 99)), p=pct_rank(obs_far, shf),
                                share_BE_far_6h=float(np.mean(ev.n_trem_30_100_6h > 0)),
                                share_control_far_day=float(np.mean(evc.n_trem_30_100_day > 0)) if len(evc) else None,
                                share_BE_broad_far_day=float(np.mean(evb.n_trem_30_100_day > 0)))
    RES["tremor_30_100"]["tied"] = bool(obs_far > np.percentile(shf, 99))
    log(f"tremor 30-100 km: {RES['tremor_30_100']}")

    # ---------------- added check (after first run): catalogued earthquakes within +-60 s
    eq = pd.read_csv(f"{ROOT}/lfe_extension/results/trackC/eq_catalog.csv")
    eq["ot"] = pd.to_datetime(eq.ot, errors="coerce"); eq = eq.dropna(subset=["ot"]).sort_values("ot")
    teq = eq.ot.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9

    def eq_near(t_s):
        i = np.clip(np.searchsorted(teq, t_s), 1, teq.size - 1)
        j = np.where(np.abs(teq[i] - t_s) < np.abs(teq[i - 1] - t_s), i, i - 1)
        return np.abs(teq[j] - t_s), j
    for x in (ev, evc, evb):
        dt_, j = eq_near(x.time_utc.astype("int64").to_numpy() / 1e9)
        x["eq_dt_s"] = dt_; x["eq_mag"] = eq.mag.to_numpy()[j]
        x["eq_dist_km"] = hav_km(x.fam.map(flat).to_numpy(), x.fam.map(flon).to_numpy(), eq.lat.to_numpy()[j], eq.lon.to_numpy()[j])
    rt = (pd.to_datetime(pd.Series(shd.ravel())).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
          + RNG.uniform(0, 86400, shd.size))
    chance = float(np.mean(eq_near(rt)[0] <= 60))
    RES["earthquakes_60s"] = dict(n_BE_within_60s=int((ev.eq_dt_s <= 60).sum()), share=float(np.mean(ev.eq_dt_s <= 60)),
                                  chance_share=chance, n_control_within_60s=int((evc.eq_dt_s <= 60).sum()),
                                  share_broad=float(np.mean(evb.eq_dt_s <= 60)),
                                  events=ev[ev.eq_dt_s <= 60][["time_PT", "fam", "eq_dt_s", "eq_mag", "eq_dist_km"]].round(1).to_dict("records"))
    log(f"earthquakes: {RES['earthquakes_60s']}")
    # same-minute multi-family coincidences (several families within 120 s)
    tt_s = np.sort(ev.time_utc.astype("int64").to_numpy() / 1e9)
    RES["n_BE_with_other_family_within_120s"] = int(np.sum((np.r_[np.inf, np.diff(tt_s)] <= 120) | (np.r_[np.diff(tt_s), np.inf] <= 120)))

    # ---------------- added check: burst-level (each burst = one unit) for next-ETS and depth
    evs2 = ev.sort_values("time_utc").reset_index(drop=True)
    th2 = evs2.time_utc.astype("int64").to_numpy() / 3.6e12
    evs2["burst"] = np.r_[0, np.cumsum(np.diff(th2) >= 24)]
    bu = evs2.groupby("burst").agg(fam=("fam", "first"), date=("date", "first"), tau_next=("tau_next", "first"),
                                   tau_prev=("tau_prev", "first"), dep=("dep", "median"), n=("fam", "size"))
    bun = bu[np.isfinite(bu.tau_next)]
    shb = shuffle_dates(bun.fam.to_numpy(), el_by_f)
    ddb = pd.Series(shb.ravel()).unique(); tmb = dict(zip(ddb, tau_next(ddb)[0]))
    stb = np.vectorize(tmb.get)(shb).astype(float)
    f30b = float(np.mean(bun.tau_next < 30)); f30bs = np.mean(stb < 30, 1)
    f240b = float(np.mean(bun.tau_next > 240)); f240bs = np.mean(stb > 240, 1)
    pe = (fam.n_ETS / fam.n_ETS.sum()).to_numpy(); depv = fam.dep.to_numpy()
    mdb = RNG.multinomial(len(bu), pe, size=10000) @ depv / len(bu)
    RES["burst_level"] = dict(n_bursts=int(len(bu)), n_with_next=int(len(bun)), F30=f30b, shuffle_F30_median=float(np.median(f30bs)),
                              p_high=pct_rank(f30b, f30bs), p_low=float((np.sum(f30bs <= f30b) + 1) / (NSH + 1)),
                              F_gt240=f240b, shuffle_F_gt240_median=float(np.median(f240bs)), p_gt240_high=pct_rank(f240b, f240bs),
                              mean_dep_obs=float(bu.dep.mean()), mean_dep_ETSdraw=float(mdb.mean()),
                              dep_p99=float(np.percentile(mdb, 99)), p_deeper=pct_rank(float(bu.dep.mean()), mdb),
                              median_tau_next_obs=float(bun.tau_next.median()), median_tau_next_shuffle=float(np.median(np.nanmedian(stb, 1))))
    log(f"burst level: {RES['burst_level']}")

    # ---------------- detectability: n_sta_day of BE events vs eligible days
    RES["n_sta_day"] = dict(BE_events_median=float(ev.n_sta_day.median()), eligible_median=float(EL.n_sta_day.median()),
                            BE_events_mean=float(ev.n_sta_day.mean()), eligible_mean=float(EL.n_sta_day.mean()))
    log(f"n_sta_day: {RES['n_sta_day']}")

    # ---------------- down-dip
    pe = (fam.n_ETS / fam.n_ETS.sum()).to_numpy(); depv = fam.dep.to_numpy(); ddv = fam.downdip_km.to_numpy()

    def depth_test(evx, name, nd=10000):
        n = len(evx)
        if n == 0:
            return dict(set=name, n=0)
        draws = RNG.multinomial(n, pe, size=nd)
        md = draws @ depv / n; mdd = draws @ ddv / n
        o = float(evx.dep.mean()); odd = float(evx.downdip_km.mean())
        return dict(set=name, n=n, mean_dep_obs=o, mean_dep_ETSdraw=float(md.mean()), dep_p99=float(np.percentile(md, 99)),
                    dep_p01=float(np.percentile(md, 1)), p_deeper=pct_rank(o, md),
                    mean_downdip_obs=odd, mean_downdip_ETSdraw=float(mdd.mean()), downdip_p99=float(np.percentile(mdd, 99)),
                    p_downdip=pct_rank(odd, mdd), deeper=bool(o > np.percentile(md, 99)))
    DT = [depth_test(ev, "BE"), depth_test(evb, "BE_broad"), depth_test(evc, "BE scrambled control")]
    from scipy.stats import spearmanr
    ff = fam[fam.n_ETS >= 10]
    rho, prho = spearmanr(ff.BE_over_ETS, ff.dep); rho2, p2 = spearmanr(ff.BE_over_ETS, ff.downdip_km)
    RES["depth"] = dict(tests=DT, spearman_ratio_vs_depth=float(rho), p=float(prho), n_fam=int(len(ff)),
                        spearman_ratio_vs_downdip=float(rho2), p_downdip=float(p2),
                        ETS_mean_dep=float((ets.fam.map(F.dep)).mean()))
    # S13 check
    A13 = sel(T13); ev13 = A13[A13.BE].copy(); ev13["dep"] = F.loc[ev13.fam, "dep"].values
    ev13["downdip_km"] = F.loc[ev13.fam, "downdip_km"].values
    ets13 = A13[A13.major & A13.in_episode]
    pe13 = ets13.fam.value_counts().reindex(fam.fid).fillna(0).to_numpy(); pe13 = pe13 / pe13.sum()
    d13 = RNG.multinomial(len(ev13), pe13, size=10000) @ depv / max(len(ev13), 1)
    ev_p21 = ev[ev.date_ts < "2021-04-01"]
    RES["S13_check"] = dict(n_BE_S13=int(len(ev13)), n_BE_c1_S13=int(sel(T13, "c1").BE.sum()),
                            n_BE_prod_same_period=int(len(ev_p21)),
                            mean_dep_obs=float(ev13.dep.mean()), mean_dep_ETSdraw=float(d13.mean()),
                            dep_p99=float(np.percentile(d13, 99)), p_deeper=pct_rank(float(ev13.dep.mean()), d13),
                            top_fams_S13=ev13.fam.value_counts().head(8).to_dict(),
                            top_fams_prod_same_period=ev_p21.fam.value_counts().head(8).to_dict())
    log(f"depth: {json.dumps(RES['depth'], indent=1)}\nS13: {RES['S13_check']}")

    # ---------------- GNSS stack
    RES["gnss"] = gnss_stack(ev, evb, E, maj)
    log(f"gnss: {json.dumps(RES['gnss'], indent=1)}")

    # ---------------- save event tables
    keep = ["time_utc", "time_PT", "date", "fam", "group", "dep", "along_km", "downdip_km", "P_same", "n_sta_day",
            "n_tremor_day", "quiet", "n_trem30", "tau_next", "tau_prev", "next_major", "n_trem_30_100_day",
            "n_trem_30_100_6h", "eq_dt_s", "eq_mag", "eq_dist_km", "flag_ep_2017_03", "flag_ep_2025_10"]
    ev[keep].sort_values("time_utc").to_csv(f"{OUT}/between_events.csv", index=False)
    evc[keep].sort_values("time_utc").to_csv(f"{OUT}/between_events_control.csv", index=False)
    evb[keep + ["series"]].sort_values("time_utc").to_parquet(f"{OUT}/between_broad_events.parquet")
    json.dump(RES, open(f"{OUT}/s7_summary.json", "w"), indent=1, default=str)
    figures(F, E, maj, ev, evc, evb, fam, EP, RES, K)
    log("=== S7 done ===")


# ------------------------------------------------------------------------------------------- GNSS
def gnss_stack(ev, evb, E, maj):
    sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    from sse_strain.matio import load_mat
    import slip_build as SB
    Xd = load_mat(f"{SB.SOLDIR}/X.mat")["Xd"]
    dr = SB.fixed_dates(np.asarray(Xd["timeline"], float).ravel())
    names = np.asarray(Xd["name"]).astype(str); ts = np.asarray(Xd["ts"], float); m = np.asarray(Xd["data_mask"]).astype(bool)
    cal = np.arange(dr[0], dr[-1] + 1); ix = (dr - cal[0]).astype(int)
    Y = np.full((ts.shape[0], cal.size), np.nan); Y[:, ix] = np.where(m, ts, np.nan)
    sta = np.array([n[:-1] for n in names[::3]])
    W = np.lib.stride_tricks.sliding_window_view(np.concatenate([Y, np.full((Y.shape[0], 4), np.nan)], 1), 5, axis=1)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); M5 = np.nanmedian(W, axis=-1)
    M5[np.isfinite(W).sum(-1) < 3] = np.nan
    s1 = json.load(open(f"{ROOT}/svi_science/results/s1/tierA_summary.json"))
    az = s1["own_azimuths"]; north = s1.get("north")
    rows = [int(np.flatnonzero(sta == s)[0]) for s in SVI14]
    U = np.array([[np.sin(np.deg2rad(az[s])), np.cos(np.deg2rad(az[s]))] for s in SVI14])
    re, rn = 3 * np.array(rows), 3 * np.array(rows) + 1
    # per-day projected offset M5[d+1] - M5[d-5] for every calendar day
    n = cal.size; d = np.arange(5, n - 1)
    dE = M5[re][:, d + 1] - M5[re][:, d - 5]; dN = M5[rn][:, d + 1] - M5[rn][:, d - 5]
    proj = dE * U[:, 0:1] + dN * U[:, 1:2]                       # (n_sta, n_days)
    off = np.full(n, np.nan); off[d] = np.nanmean(proj, axis=0)
    nmask = np.array([s in (north or []) for s in SVI14])
    offn = np.full(n, np.nan)
    if nmask.any():
        offn[d] = np.nanmean(proj[nmask], axis=0)
    # allowed days
    bad = np.zeros(n, bool)
    for s, e in zip(E.start.values.astype("datetime64[D]"), E.end.values.astype("datetime64[D]")):
        a = int((s - cal[0]).astype(int)) - 2 - 6; b = int((e - cal[0]).astype(int)) + 2 + 6
        bad[max(a, 0):min(b, n - 1) + 1] = True
    lo = int((np.datetime64("2014-01-15") - cal[0]).astype(int)); hi = int((np.datetime64(P_END.date(), "D") - cal[0]).astype(int)) - 6
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
        C = np.empty(nset); Cn = np.empty(nset)
        cand = []
        for i in idx:
            dd = np.abs(doy[alw] - doy[i]); dd = np.minimum(dd, 365 - dd)
            cand.append(alw[dd <= 30])
        for j in range(nset):
            pick = np.array([c[RNG.integers(0, c.size)] for c in cand])
            C[j] = np.nanmean(off[pick]); Cn[j] = np.nanmean(offn[pick])
        sig = float(C.std()); mu = float(C.mean())
        return dict(set=name, n_days=int(idx.size), mean_mm=obs, ctrl_mean_mm=mu, ctrl_sigma_mm=sig, z=(obs - mu) / sig,
                    bound_3sigma_mm=mu + 3 * sig - mu, north_mean_mm=obsn, north_ctrl_mean=float(Cn.mean()),
                    north_sigma=float(Cn.std()), north_z=float((obsn - Cn.mean()) / Cn.std()),
                    slip_like=bool((obs - mu) / sig > 3), ctrl_values=C.round(4).tolist())
    cb = clusters(ev); cbb = clusters(evb)
    # minor episodes: centre day, own exclusion relaxed (they ARE episodes) -> evaluate directly against the control
    mins = E[~E.major]
    mid = np.array([int(((s + (e - s) / 2).to_datetime64().astype("datetime64[D]") - cal[0]).astype(int)) for s, e in zip(mins.start, mins.end)])
    mid = mid[(mid > 6) & (mid < n - 7)]
    res = [stack(cb, "BE cluster days"), stack(cbb, "BE_broad cluster days")]
    # minor reference: observed at minor mid-days, control as for BE days
    allowed_bak = allowed.copy(); allowed[mid] = np.isfinite(off[mid])
    r = stack(mid, "minor-episode mid-days (reference)"); allowed[:] = allowed_bak
    res.append(r)
    res[0]["n_clusters_total"] = int(cb.size); res[1]["n_clusters_total"] = int(cbb.size)
    np.savez(f"{OUT}/gnss_stack.npz", **{r["set"].split()[0].replace("-", "_"): np.array(r.get("ctrl_values", [])) for r in res},
             obs=np.array([r.get("mean_mm", np.nan) for r in res]))
    for r in res:
        r.pop("ctrl_values", None)
    return res


# ------------------------------------------------------------------------------------------- figures
def figures(F, E, maj, ev, evc, evb, fam, EP, RES, K):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6})
    INK, MUTED = "#222222", "#777777"
    toPT = lambda s: s.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.tz_localize(None)

    # ---- 1. timeline
    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]),
                           layout="constrained")
    a = ax[0]
    for r in E.itertuples():
        a.axvspan(np.datetime64(r.start), np.datetime64(r.end + pd.Timedelta(days=1)),
                  color="#b8cbe4" if r.major else "#e3ebf5", lw=0, zorder=0)
    tb = toPT(evb.time_utc).to_numpy(); te = toPT(ev.time_utc).to_numpy()
    a.scatter(tb, evb.along_km, s=5, c="#b0b0b0", lw=0, zorder=1,
              label=f"outside episodes ±2 d but with tremor ≤ 30 km ({len(evb) - len(ev)})")
    sc = a.scatter(te, ev.along_km, s=22, c=ev.dep, cmap="viridis_r", vmin=27, vmax=46, lw=0.4,
                   edgecolor="white", label=f"between-episode, no tremor ≤ 30 km ({len(ev)})", zorder=3)
    if len(evc):
        a.scatter(toPT(evc.time_utc).to_numpy(), evc.along_km, s=40, marker="x", c="#d62728", lw=1.2,
                  label=f"scrambled control, same definition ({len(evc)})", zorder=4)
    a.set_ylabel("along strike (km, + = NW)"); a.set_ylim(-64, 50)
    a.legend(loc="lower left", fontsize=8, frameon=True, framealpha=0.9, ncol=3)
    a.set_title("SVI LFEs between ETS episodes (P ≥ 0.75, counted once). Shaded = the 79 episodes (darker = 39 majors)",
                loc="left", fontsize=10)
    fig.colorbar(sc, ax=a, pad=0.01, shrink=0.9).set_label("family depth (km)")
    b = ax[1]
    mo = toPT(ev.time_utc).dt.to_period("M").value_counts().sort_index()
    mob = toPT(evb.time_utc).dt.to_period("M").value_counts().sort_index()
    b.bar(mob.index.to_timestamp().to_numpy(), mob.values, width=np.timedelta64(25, "D"), color="#c8c8c8", label="all outside episodes ±2 d")
    b.bar(mo.index.to_timestamp().to_numpy(), mo.values, width=np.timedelta64(25, "D"), color="#2c7fb8", label="no tremor ≤ 30 km")
    b.set_ylabel("events / month"); b.set_yscale("symlog", linthresh=5); b.legend(fontsize=8, frameon=False, loc="upper left")
    b.set_xlabel("time (Pacific)")
    fig.colorbar(sc, ax=b, pad=0.01).ax.set_visible(False)
    fig.savefig(f"{FIG}/s7_timeline.png", dpi=150); plt.close(fig)

    # ---- 2. map + depth
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    a.scatter(fam.lon, fam.lat, s=10, c="#cccccc", lw=0, zorder=1)
    m = fam.n_BE > 0
    sc = a.scatter(fam.lon[m], fam.lat[m], s=20 + 12 * fam.n_BE[m], c=fam.dep[m], cmap="viridis_r", vmin=27, vmax=46,
                   edgecolor=INK, lw=0.5, zorder=3)
    for r in fam[fam.n_BE >= 5].itertuples():
        a.annotate(f"{r.fid} ({r.n_BE})", (r.lon, r.lat), xytext=(5, 4), textcoords="offset points", fontsize=7, color=INK)
    g = fam.groupby("group")[["lon", "lat"]].mean()
    for k, r in g.iterrows():
        a.text(r.lon, r.lat - 0.035, f"G{k}", fontsize=7, color=MUTED, ha="center")
    a.set_xlabel("longitude"); a.set_ylabel("latitude"); a.set_aspect(1 / np.cos(np.deg2rad(48.5)))
    a.set_title("Families firing between episodes (size = events, colour = depth; grey = none)", loc="left", fontsize=9)
    fig.colorbar(sc, ax=a, pad=0.01).set_label("family depth (km)")
    b = ax[1]
    bins = np.arange(26, 48, 2)
    w_ets = np.histogram(fam.dep, bins, weights=fam.n_ETS)[0]; w_be = np.histogram(fam.dep, bins, weights=fam.n_BE)[0]
    w_bb = np.histogram(fam.dep, bins, weights=fam.n_BE_broad)[0]
    c = 0.5 * (bins[1:] + bins[:-1])
    b.plot(c, w_ets / w_ets.sum(), "-", color=MUTED, lw=2, label="inside the 39 major ETS")
    b.plot(c, w_bb / w_bb.sum(), "-", color="#9ecae1", lw=2, label="outside episodes ±2 d (all)")
    b.plot(c, w_be / max(w_be.sum(), 1), "-", color="#08519c", lw=2, label="between-episode, no local tremor")
    b.set_xlabel("family depth (km, Bostock)"); b.set_ylabel("share of events"); b.legend(fontsize=8, frameon=False)
    d = RES["depth"]["tests"][0]
    b.set_title(f"Between-episode mean depth {d['mean_dep_obs']:.1f} km; ETS shares predict {d['mean_dep_ETSdraw']:.1f} km",
                loc="left", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/s7_map_depth.png", dpi=150); plt.close(fig)

    # ---- 3. superposed epoch + GNSS
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, (setn, col, ttl) in zip(ax[:2], (("BE", "tau_next", "days BEFORE the next major ETS starts"),
                                             ("BE", "tau_prev", "days AFTER the previous major ETS ended"))):
        t = EP[(EP.set == setn) & (EP.axis == col)]
        x = np.arange(len(t))
        a.fill_between(x, t.shuffle_p01, t.shuffle_p99, color="#dddddd", step="mid", label="family-day shuffle 1–99 %")
        a.step(x, t.shuffle_mean, where="mid", color=MUTED, lw=1)
        a.plot(x, t.events, "o-", color="#08519c", lw=2, ms=5, label="observed events")
        tb = EP[(EP.set == "BE_broad") & (EP.axis == col)]
        a.set_xticks(x); a.set_xticklabels(t.bin, rotation=0, fontsize=8); a.set_xlabel(ttl)
        a.set_ylabel("between-episode events (no local tremor)"); a.legend(fontsize=8, frameon=False)
    ax[0].set_title("Do they pile up before the next ETS?", loc="left", fontsize=9)
    ax[1].set_title("…or right after the last one?", loc="left", fontsize=9)
    a = ax[2]
    z = np.load(f"{OUT}/gnss_stack.npz")
    g0 = RES["gnss"][0]
    if "BE" in z.files and z["BE"].size:
        a.hist(z["BE"], 40, color="#cccccc", label="season-matched random days")
        a.axvline(g0["mean_mm"], color="#08519c", lw=2, label=f"LFE days (n={g0['n_days']})")
        a.axvline(g0["ctrl_mean_mm"] + 3 * g0["ctrl_sigma_mm"], color=MUTED, ls="--", lw=1, label="3σ")
    a.set_xlabel("mean GNSS offset over ~6 d, ETS direction (mm)"); a.set_ylabel("random sets")
    a.set_title("GNSS: any slip on between-episode LFE days?", loc="left", fontsize=9); a.legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/s7_superposed_epoch.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
