#!/usr/bin/env python
"""q5_bursts (2026-10-07) -- are the deep between-ETS LFE bursts small slow slips?
Rules frozen in svi_science/results/q5_bursts/PLAN.md before any number.  CPU only, < 10 GB.

    python svi_science/scripts/q5_bursts.py
"""
from __future__ import annotations
import os, sys, json, time, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q5_bursts"
FIG = f"{ROOT}/svi_science/figures"
S0 = f"{ROOT}/svi_science/results/s0"
S7 = f"{ROOT}/svi_science/results/s7"
WAD = f"{ROOT}/svi_science/results/a_wa_between"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
P_END = pd.Timestamp("2026-05-27")
NSH = 1000
NPERM = 5000
RNG = np.random.default_rng(20261007)
SVI14 = ["ALBH", "PGC5", "SC02", "SC04", "NEAH", "PTRF", "PTAA", "NANO", "P435", "P436", "P439", "CLRS", "COUP", "SC03"]
WA12 = ["CHCM", "KTBW", "P426", "ELSR", "P437", "CUSH", "P419", "PFLD", "SSHO", "SEAT", "SMAI", "P399"]
BINS = [(0, 30), (30, 60), (60, 100)]
MU = 40e9
LOG = open(f"{OUT}/q5.log", "w")
RES = {}
sys.path.insert(0, f"{ROOT}/svi_science/scripts")
import s3_build as S3B  # noqa: E402  (frame projection only)


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def p_hi(obs, ctrl):
    ctrl = np.asarray(ctrl); return float((np.sum(ctrl >= obs) + 1) / (ctrl.size + 1))


def p_lo(obs, ctrl):
    ctrl = np.asarray(ctrl); return float((np.sum(ctrl <= obs) + 1) / (ctrl.size + 1))


def pad_mask(dates, E, pad=2):
    d = pd.to_datetime(pd.Series(dates)).values
    m = np.zeros(len(d), bool)
    for s, e in zip(E.start.values, E.end.values):
        m |= (d >= s - np.timedelta64(pad, "D")) & (d <= e + np.timedelta64(pad, "D"))
    return m


def to_pt(t):
    return pd.to_datetime(pd.Series(t)).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M").to_numpy()


# --------------------------------------------------------------------------------------------- frames
FR_SVI = json.load(open(f"{S0}/frame_groups.json"))
FR_SVI = dict(origin_lat=FR_SVI["origin_lat"], origin_lon=FR_SVI["origin_lon"], strike_az=FR_SVI["strike_az"],
              downdip_az=FR_SVI["mean_downdip_az"])
FR_WA = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["nw_frame"]


def proj(fr, lat, lon):
    a, d, _, _ = S3B.project(fr, lat, lon); return a, d


# --------------------------------------------------------------------------------------------- load
def load_tremor():
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["time"])
    tr["time"] = tr.time.astype("datetime64[ns]")
    tr = tr[(tr.time >= pd.Timestamp("2009-07-01")) & tr.lat.between(45.5, 51.0) & tr.lon.between(-127.5, -120.0)]
    tr = tr.sort_values("time").reset_index(drop=True)
    return tr


def load_svi():
    F = pd.read_csv(f"{S0}/families.csv").set_index("fid", drop=False)
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    ev = pd.read_csv(f"{S7}/between_events.csv", parse_dates=["time_utc"]).rename(columns={"fam": "fid"})
    evb = pd.read_parquet(f"{S7}/between_broad_events.parquet").rename(columns={"fam": "fid"})
    evb = evb[evb.series == "data"].copy()
    D = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable"])
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    Dp = D[(D.run == "prod") & D.observable].copy()
    Dp = Dp[(Dp.date >= "2014-01-01") & (Dp.date <= "2026-05-27")]
    Dp = Dp.merge(X, on=["fid", "date"], how="left"); Dp["n30"] = Dp.n30.fillna(0)
    Dp["in_pad"] = pad_mask(Dp.date, E, 2)
    pool = {f: g.date.to_numpy().astype("datetime64[D]") for f, g in Dp[~Dp.in_pad & (Dp.n30 == 0)].groupby("fid")}
    poolb = {f: g.date.to_numpy().astype("datetime64[D]") for f, g in Dp[~Dp.in_pad].groupby("fid")}
    for x in (ev, evb):
        x["group"] = F.loc[x.fid, "group"].astype(str).values
        x["lat"] = F.loc[x.fid, "lat"].values; x["lon"] = F.loc[x.fid, "lon"].values
        x["along_km"] = F.loc[x.fid, "along_strike_km"].values; x["downdip_km"] = F.loc[x.fid, "downdip_km"].values
        x["dep"] = F.loc[x.fid, "dep"].values
    # minor episodes + their LFEs (for test 3)
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8", "time_utc"])
    T = T[(T.series == "data") & (T.P_same >= 0.75) & T.keep_xfam8 & (T.date >= "2014-01-01") & (T.date <= "2026-05-27")]
    T = T.rename(columns={"fam": "fid"})[["date", "fid", "time_utc"]]
    return dict(name="SVI", F=F, E=E, ev=ev, evb=evb, pool=pool, poolb=poolb, T=T, fr=FR_SVI, top_group="10",
                stations=SVI14, start="2014-01-01", gnss_lo="2014-01-15",
                famF=F.assign(along=F.along_strike_km, downdip=F.downdip_km))


def load_wa():
    F = pd.read_csv(f"{WAD}/families_wa.csv").set_index("fid", drop=False)
    E = pd.read_csv(f"{EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    ev = pd.read_csv(f"{WAD}/between_events.csv", parse_dates=["time_utc"])
    evb = pd.read_parquet(f"{WAD}/between_broad_events.parquet"); evb = evb[evb.series == "data"].copy()
    e12 = pd.read_csv(f"{WAD}/between_events_E1E2.csv", parse_dates=["time_utc"])
    X = pd.read_parquet(f"{WAD}/exposure_famday.parquet")
    X = X[X.era.isin(["E3"]) & (X.date >= "2015-01-01") & (X.date <= "2026-05-27")].copy()
    X["in_pad"] = pad_mask(X.date, E, 2)
    pool = {f: g.date.to_numpy().astype("datetime64[D]") for f, g in X[~X.in_pad & (X.day_class == "quiet")].groupby("fid")}
    poolb = {f: g.date.to_numpy().astype("datetime64[D]") for f, g in X[~X.in_pad & (X.day_class != "unknown")].groupby("fid")}
    X12 = pd.read_parquet(f"{WAD}/exposure_famday.parquet")
    X12 = X12[X12.era.isin(["E1", "E2"])].copy(); X12["in_pad"] = pad_mask(X12.date, E, 2)
    pool12 = {f: g.date.to_numpy().astype("datetime64[D]") for f, g in X12[~X12.in_pad & (X12.day_class == "quiet")].groupby("fid")}
    for x in (ev, evb, e12):
        x["group"] = F.loc[x.fid, "source_group"].astype(str).values
        x["lat"] = F.loc[x.fid, "lat"].values; x["lon"] = F.loc[x.fid, "lon"].values
        x["along_km"] = F.loc[x.fid, "along_km"].values; x["downdip_km"] = F.loc[x.fid, "downdip_km"].values
        x["dep"] = F.loc[x.fid, "dep"].values
    T = pd.read_parquet(f"{WAD}/events_all.parquet", columns=["date", "series", "fid", "P_own_WAVE_PNET", "once_D", "time_utc"])
    T = T[(T.series == "data") & (T.P_own_WAVE_PNET >= 0.75) & T.once_D & (T.date >= "2015-01-01") & (T.date <= "2026-05-27")]
    T = T[["date", "fid", "time_utc"]]
    return dict(name="WA", F=F, E=E, ev=ev, evb=evb, e12=e12, pool=pool, poolb=poolb, pool12=pool12, T=T, fr=FR_WA,
                top_group="G019", stations=WA12, start="2015-01-01", gnss_lo="2015-01-15",
                famF=F.assign(along=F.along_km, downdip=F.downdip_km))


# --------------------------------------------------------------------------------------------- tremor counts per family-day
def family_day_counts(F, tr, day0):
    """-> dict fid -> (n_days, 3) int array of epicentre counts per distance bin, day index from day0."""
    di = (tr.time.values.astype("datetime64[D]") - day0).astype("timedelta64[D]").astype(int)
    nd = int(di.max()) + 40
    assert nd < 8000, nd
    tla, tlo = tr.lat.to_numpy(), tr.lon.to_numpy()
    out = {}
    for f, r in F.iterrows():
        dk = hav_km(r.lat, r.lon, tla, tlo)
        C = np.zeros((nd, 3), np.int32)
        for b, (lo, hi) in enumerate(BINS):
            m = (dk >= lo) & (dk < hi) if b < 2 else (dk >= lo) & (dk <= hi)
            np.add.at(C[:, b], di[m], 1)
        out[f] = C
    return out


def bursts_of(ev):
    ev = ev.sort_values("time_utc").reset_index(drop=True)
    t = ev.time_utc.astype("datetime64[ns]").astype("int64").to_numpy() / 3.6e12
    ev["burst"] = np.r_[0, np.cumsum(np.diff(t) >= 24)]
    return ev


# --------------------------------------------------------------------------------------------- test 1
def test1(R, ev, pool, C, day0, tr, setname, bins_tested):
    ev = bursts_of(ev)
    day = ev.time_utc.values.astype("datetime64[D]"); di = (day - day0).astype(int)
    fids = ev.fid.to_numpy()
    cnt = np.array([C[f][d] for f, d in zip(fids, di)])                    # (n_ev, 3)
    # +-6 h column
    tt = tr.time.values.astype("int64") / 3.6e12; tla, tlo = tr.lat.to_numpy(), tr.lon.to_numpy()
    c6 = np.zeros((len(ev), 3), int)
    for i, (t, la, lo) in enumerate(zip(ev.time_utc.values.astype("int64") / 3.6e12, ev.lat, ev.lon)):
        i0, i1 = np.searchsorted(tt, [t - 6, t + 6])
        if i1 > i0:
            dk = hav_km(la, lo, tla[i0:i1], tlo[i0:i1])
            for b, (l, h) in enumerate(BINS):
                c6[i, b] = int(((dk >= l) & (dk < h if b < 2 else dk <= h)).sum())
    for b in range(3):
        ev[f"n_day_{BINS[b][0]}_{BINS[b][1]}"] = cnt[:, b]; ev[f"n_6h_{BINS[b][0]}_{BINS[b][1]}"] = c6[:, b]
    # bursts
    g = ev.groupby("burst")
    bid = ev.burst.to_numpy(); nb = bid.max() + 1
    first_f = g.fid.first().to_numpy(); first_d = g.time_utc.first().values.astype("datetime64[D]")
    off = (day - first_d[bid]).astype(int)                                 # event day offset inside its burst
    obs_b = np.zeros((nb, 3), bool)
    for b in range(3):
        obs_b[:, b] = pd.Series(cnt[:, b] > 0).groupby(bid).any().to_numpy()
    # stacked count array for vectorised look-ups
    fl = sorted(C); row = {f: i for i, f in enumerate(fl)}; ARR = np.stack([C[f] for f in fl]) > 0
    fr_ = np.array([row[f] for f in fids]); nd = ARR.shape[1]
    starts = np.flatnonzero(np.r_[True, np.diff(bid) != 0])
    for f in set(first_f) | set(fids):
        assert f in pool, f"family {f} has no pool days"
    # random-day control, burst level (whole burst shifted to a random pool day of its first family)
    ctl_b = np.zeros((NSH, 3))
    for j in range(NSH):
        start = np.array([pool[f][RNG.integers(0, pool[f].size)] for f in first_f])
        dd = np.clip((start[bid] - day0).astype(int) + off, 0, nd - 1)
        cc = ARR[fr_, dd].astype(int)
        ctl_b[j] = (np.add.reduceat(cc, starts, axis=0) > 0).mean(0)
    # event level control
    pick = np.array([pool[f][RNG.integers(0, pool[f].size, NSH)] for f in fids]).T   # (NSH, n_ev)
    pdi = np.clip((pick - day0).astype(int), 0, nd - 1)
    ctl_e = ARR[fr_[None, :], pdi].mean(1)                                    # (NSH, 3)
    rows = []
    for b in range(3):
        lab = f"{BINS[b][0]}-{BINS[b][1]} km"
        ob = float(obs_b[:, b].mean()); oe = float((cnt[:, b] > 0).mean())
        rows.append(dict(region=R, set=setname, bin=lab, tested=lab in bins_tested, n_events=len(ev), n_bursts=int(nb),
                         burst_share=ob, burst_ctrl_median=float(np.median(ctl_b[:, b])),
                         burst_ctrl_p01=float(np.percentile(ctl_b[:, b], 1)), burst_ctrl_p99=float(np.percentile(ctl_b[:, b], 99)),
                         burst_p=p_hi(ob, ctl_b[:, b]), event_share=oe, event_ctrl_median=float(np.median(ctl_e[:, b])),
                         event_ctrl_p99=float(np.percentile(ctl_e[:, b], 99)), event_p=p_hi(oe, ctl_e[:, b]),
                         event_share_6h=float((c6[:, b] > 0).mean()),
                         tied=bool(lab in bins_tested and ob > np.percentile(ctl_b[:, b], 99))))
    return pd.DataFrame(rows), ev


# --------------------------------------------------------------------------------------------- test 2
def ts_slope(t, y, cross=False):
    i, j = np.triu_indices(t.size, 1)
    dt = t[j] - t[i]; ok = (dt != 0) & ((y[j] != y[i]) if cross else True)
    return float(np.median((y[j] - y[i])[ok] / dt[ok])) if ok.any() else np.nan


def perm_slope(t, y, n=NPERM, cross=False):
    i, j = np.triu_indices(t.size, 1)
    out = np.empty(n)
    for k in range(n):
        tp = RNG.permutation(t); dt = tp[j] - tp[i]; ok = (dt != 0) & ((y[j] != y[i]) if cross else True)
        out[k] = np.median((y[j] - y[i])[ok] / dt[ok])
    return out


def burst_table(D, ev, tr, label):
    """Per-burst descriptors (tests 2, 3, 4)."""
    ev = bursts_of(ev)
    fr = D["fr"]
    trt = tr.time.values.astype("int64") / 3.6e12
    rows = []
    for b, g in ev.groupby("burst"):
        t = g.time_utc.values.astype("int64") / 3.6e12
        t0, t1 = t.min(), t.max()
        fpos = g.drop_duplicates("fid")[["fid", "lat", "lon", "along_km", "downdip_km", "group"]]
        la, lo = fpos.lat.to_numpy(), fpos.lon.to_numpy()
        sep = float(max(hav_km(la[i], lo[i], la, lo).max() for i in range(len(la)))) if len(la) > 1 else 0.0
        clat, clon = float(g.lat.mean()), float(g.lon.mean())
        i0, i1 = np.searchsorted(trt, [t0 - 6, t1 + 6])
        tw = tr.iloc[i0:i1]; dk = hav_km(clat, clon, tw.lat.to_numpy(), tw.lon.to_numpy()); tw = tw[dk <= 100]
        r = dict(set=label, burst=int(b), start_PT=to_pt([g.time_utc.iloc[0]])[0], end_PT=to_pt([g.time_utc.iloc[-1]])[0],
                 start_utc=str(g.time_utc.iloc[0]), n=len(g), n_fam=int(g.fid.nunique()),
                 groups=";".join(sorted(set(g.group))), duration_h=float(t1 - t0), max_fam_sep_km=sep,
                 lat=clat, lon=clon, along_km=float(g.along_km.mean()), dep=float(g.dep.mean()),
                 n_trem_100km_pm6h=int(len(tw)), fids=";".join(map(str, sorted(set(g.fid)))),
                 main_day=str(pd.Series(g.time_utc.values.astype("datetime64[D]")).value_counts().idxmax())[:10])
        if len(tw):
            ta, _ = proj(fr, tw.lat.to_numpy(), tw.lon.to_numpy())
            r["trem_along_span_km"] = float(np.percentile(ta, 95) - np.percentile(ta, 5))
        # migration (LFE positions)
        a = g.along_km.to_numpy(float); d = g.downdip_km.to_numpy(float)
        npos_ok = (len(g) >= 4) and (sep >= 5) and np.isfinite(a).all()
        r["migr_eligible"] = bool(npos_ok)
        if npos_ok:
            for ax, y in (("along", a), ("downdip", d)):
                s = ts_slope(t, y); ps = perm_slope(t, y)
                r[f"{ax}_kmh"] = s; r[f"{ax}_p"] = p_hi(abs(s), np.abs(ps))
                r[f"{ax}_range_km"] = float(np.ptp(y))
            r["migrates"] = bool(min(r["along_p"], r["downdip_p"]) < 0.025)
            # added after the first run: same-family pairs carry no position change and force the plain
            # Theil-Sen slope to 0 when one family dominates -> cross-family pairs only (column + same pass rule)
            for ax, y in (("along", a), ("downdip", d)):
                if np.ptp(y) > 0:
                    s = ts_slope(t, y, True); ps = perm_slope(t, y, NPERM, True)
                    r[f"{ax}_x_kmh"] = s; r[f"{ax}_x_p"] = p_hi(abs(s), np.abs(ps))
                else:
                    r[f"{ax}_x_kmh"] = 0.0; r[f"{ax}_x_p"] = 1.0
            r["migrates_x"] = bool(min(r["along_x_p"], r["downdip_x_p"]) < 0.025)
        # migration (tremor) column
        if len(tw) >= 5:
            ta, td = proj(fr, tw.lat.to_numpy(), tw.lon.to_numpy())
            tt = tw.time.values.astype("int64") / 3.6e12
            for ax, y in (("trem_along", ta), ("trem_downdip", td)):
                s = ts_slope(tt, y); ps = perm_slope(tt, y, 2000)
                r[f"{ax}_kmh"] = s; r[f"{ax}_p"] = p_hi(abs(s), np.abs(ps))
        rows.append(r)
    return pd.DataFrame(rows), ev


def minor_table(D, tr, C, day0):
    E = D["E"]; T = D["T"]; F = D["famF"]; fr = D["fr"]
    m = E[~E.major & (E.start >= D["start"]) & (E.end <= P_END)].copy()
    lat0, lat1 = F.lat.min() - .3, F.lat.max() + .3; lon0, lon1 = F.lon.min() - .3, F.lon.max() + .3
    trb = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1)]
    rows = []
    for r in m.itertuples():
        x = T[(T.date >= str(r.start.date())) & (T.date <= str(r.end.date()))]
        fa = F.loc[x.fid.unique(), "along"] if len(x) else pd.Series(dtype=float)
        tw = trb[(trb.time >= r.start) & (trb.time < r.end + pd.Timedelta(days=1))]
        ta = proj(fr, tw.lat.to_numpy(), tw.lon.to_numpy())[0] if len(tw) else np.array([])
        rows.append(dict(region=D["name"], episode=r.episode, start=str(r.start.date()), end=str(r.end.date()),
                         span_days=int(r.span_days), n_lfe=len(x), n_fam=int(x.fid.nunique()),
                         fam_along_span_km=float(np.ptp(fa)) if len(fa) > 1 else 0.0, n_trem=len(tw),
                         trem_along_span_km=float(np.percentile(ta, 95) - np.percentile(ta, 5)) if len(ta) >= 5 else np.nan,
                         lfe_duration_h=float((x.time_utc.max() - x.time_utc.min()).total_seconds() / 3600) if len(x) > 1 else 0.0))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------- test 4
def cv(x):
    x = np.diff(np.sort(x)); return float(x.std() / x.mean()) if x.size > 1 and x.mean() > 0 else np.nan


def test4(D, B, ev, pool):
    first = ev.groupby("burst").fid.first().to_numpy()
    t = pd.to_datetime(B.start_utc, format="ISO8601").values.astype("datetime64[D]").astype(int).astype(float)
    top = B.groups.str.split(";").apply(lambda s: D["top_group"] in s).to_numpy()
    out = {}
    for name, m in (("all_bursts", np.ones(len(B), bool)), (f"top_group_{D['top_group']}", top),
                    ("big_bursts_ge3", B.n.to_numpy() >= 3)):
        if m.sum() < 3:
            out[name] = dict(n=int(m.sum())); continue
        o = cv(t[m]); ctl = np.empty(NSH)
        for j in range(NSH):
            ctl[j] = cv(np.array([pool[f][RNG.integers(0, pool[f].size)] for f in first[m]]).astype(int).astype(float))
        out[name] = dict(n=int(m.sum()), cv=o, ctrl_cv_median=float(np.median(ctl)), ctrl_p01=float(np.percentile(ctl, 1)),
                         ctrl_p99=float(np.percentile(ctl, 99)), p_regular=p_lo(o, ctl), p_clustered=p_hi(o, ctl),
                         median_interval_d=float(np.median(np.diff(np.sort(t[m])))),
                         verdict="regular" if o < np.percentile(ctl, 1) else ("clustered" if o > np.percentile(ctl, 99)
                                                                                else "consistent with random timing"),
                         years=sorted(set(pd.to_datetime(B.start_utc[m], format="ISO8601").dt.year.tolist())))
    # group participation
    gs = B.groups.str.split(";").explode()
    gp = gs.value_counts().rename("n_bursts").to_frame()
    gp["share_of_bursts"] = gp.n_bursts / len(B)
    gp["n_years"] = [B[B.groups.str.split(";").apply(lambda s: g in s)].start_utc.str[:4].nunique() for g in gp.index]
    gp["n_events"] = [int(ev[ev.group == g].shape[0]) for g in gp.index]
    gp["dep_mean"] = [float(ev[ev.group == g].dep.mean()) for g in gp.index]
    out["groups"] = gp.reset_index().rename(columns={"index": "group"}).head(10).to_dict("records")
    return out


# --------------------------------------------------------------------------------------------- test 5
def gnss_setup(D, E, lo_date):
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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); M5 = np.nanmedian(W, axis=-1)
    M5[np.isfinite(W).sum(-1) < 3] = np.nan
    s1 = json.load(open(f"{ROOT}/svi_science/results/s1/tierA_summary.json"))
    az = s1["own_azimuths"] if D["name"] == "SVI" else s1["wa_azimuths"]
    stl = D["stations"]
    rows = [int(np.flatnonzero(sta == s)[0]) for s in stl]
    U = np.array([[np.sin(np.deg2rad(az[s])), np.cos(np.deg2rad(az[s]))] for s in stl])
    re, rn = 3 * np.array(rows), 3 * np.array(rows) + 1
    n = cal.size; d = np.arange(5, n - 1)
    dE = M5[re][:, d + 1] - M5[re][:, d - 5]; dN = M5[rn][:, d + 1] - M5[rn][:, d - 5]
    proj_ = dE * U[:, 0:1] + dN * U[:, 1:2]
    off = np.full(n, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); off[d] = np.nanmean(proj_, axis=0)
    bad = np.zeros(n, bool)
    for s, e in zip(E.start.values.astype("datetime64[D]"), E.end.values.astype("datetime64[D]")):
        a = int((s - cal[0]).astype(int)) - 2 - 6; b = int((e - cal[0]).astype(int)) + 2 + 6
        bad[max(a, 0):min(b, n - 1) + 1] = True
    lo = int((np.datetime64(lo_date) - cal[0]).astype(int)); hi = int((np.datetime64(P_END.date(), "D") - cal[0]).astype(int)) - 6
    allowed = np.zeros(n, bool); allowed[lo:hi] = True; allowed &= ~bad & np.isfinite(off)
    doy = (cal - cal.astype("datetime64[Y]")).astype(int)
    return dict(cal=cal, off=off, allowed=allowed, doy=doy, alw=np.flatnonzero(allowed), az=az)


def gstack(G, idx, nset=NSH):
    cal, off, allowed, doy, alw = G["cal"], G["off"], G["allowed"], G["doy"], G["alw"]
    keep = np.isfinite(off[idx]); idx = idx[keep]
    obs = float(np.nanmean(off[idx]))
    cand = []
    for i in idx:
        dd = np.abs(doy[alw] - doy[i]); dd = np.minimum(dd, 365 - dd); cand.append(alw[dd <= 30])
    C = np.array([np.nanmean(off[np.array([c[RNG.integers(0, c.size)] for c in cand])]) for _ in range(nset)])
    mu, sig = float(C.mean()), float(C.std())
    return dict(n_days=int(idx.size), mean_mm=obs, ctrl_mean_mm=mu, ctrl_sigma_mm=sig, z=(obs - mu) / sig,
                slip_like=bool((obs - mu) / sig > 3), bound_mm=max(obs - mu, 0) + 3 * sig), C, keep


_FAULT = {}


def mm_per_Nm(D, lat, lon, R_km):
    """Mean projected offset (mm) at the region's stack stations per N.m of pure-thrust slip on mesh patches within R."""
    if not _FAULT:
        sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src"); sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
        from sse_strain import load_fault_ascii
        import slip_build as SB
        z = np.load(f"{ROOT}/lfe_extension/results/slip_compare/slip_cache.npz")
        _FAULT.update(fault=load_fault_ascii(SB.MESH), G=z["G"], names=[str(n)[:-1] for n in z["station_names"]], SB=SB)
    fault, Gm, names, SB = _FAULT["fault"], _FAULT["G"], _FAULT["names"], _FAULT["SB"]
    p, ok, zmesh, P = SB.locate_on_mesh(fault, lat, lon)
    c = fault.centroid_xyz_km
    dist = np.sqrt(((c - P[None, :]) ** 2).sum(1))
    sel = np.flatnonzero(dist <= R_km)
    if sel.size == 0:
        sel = np.array([p])
    npch = fault.n_patch
    area = np.asarray(fault.area)[sel] * 1e6                               # km2 -> m2
    az = D["az"]
    u = []
    for s in D["stations"]:
        i = names.index(s)
        gE = -Gm[3 * i, npch + sel]; gN = -Gm[3 * i + 1, npch + sel]          # thrust = -dip-slip on this mesh
        u.append((gE * np.sin(np.deg2rad(az[s])) + gN * np.cos(np.deg2rad(az[s]))).sum())
    m0 = MU * area.sum() * 1.0                                             # 1 m uniform slip
    return float(np.mean(u) * 1000.0 / m0), dict(patch=p, inside=bool(ok), mesh_depth_km=zmesh, n_patch=int(sel.size),
                                                 area_km2=float(area.sum() / 1e6))


def mw(m0):
    return (2 / 3) * (np.log10(m0) - 9.1)


# --------------------------------------------------------------------------------------------- main
def main():
    log("=== Q5 start ===")
    tr = load_tremor()
    day0 = np.datetime64("2009-07-01")
    regs = [load_svi(), load_wa()]
    T1, B_all, MIN_all, T4 = [], [], [], {}
    for D in regs:
        R = D["name"]
        C = family_day_counts(D["F"], tr, day0)
        D["C"] = C
        # ---- test 1
        t1a, evx = test1(R, D["ev"], D["pool"], C, day0, tr, "BE (no tremor <=30 km that day)", ["30-60 km", "60-100 km"])
        t1b, _ = test1(R, D["evb"], D["poolb"], C, day0, tr, "BE_broad (any local tremor)", ["0-30 km"])
        T1 += [t1a, t1b]
        if R == "WA":
            t1c, _ = test1(R, D["e12"], D["pool12"], C, day0, tr, "BE E1/E2 2009-14 (secondary)", [])
            T1.append(t1c)
        evx.to_csv(f"{OUT}/events_{R}_tremor_bins.csv", index=False)
        # coverage column: share of major-episode family-days with an epicentre <= 30 km
        maj = D["E"][D["E"].major & (D["E"].start >= D["start"]) & (D["E"].end <= P_END)]
        days = np.concatenate([np.arange(s, e + np.timedelta64(1, "D")) for s, e in
                               zip(maj.start.values.astype("datetime64[D]"), maj.end.values.astype("datetime64[D]"))])
        di = (days - day0).astype(int)
        cov = float(np.mean([np.mean(C[f][di, 0] > 0) for f in D["F"].fid]))
        RES[f"{R}_coverage_major_famdays_with_trem30"] = cov
        log(f"[{R}] test1:\n" + pd.concat(T1[-3 if R == 'WA' else -2:]).round(3).to_string(index=False))
        log(f"[{R}] PNSN coverage: share of major-ETS family-days with an epicentre <= 30 km = {cov:.2f}")
        # ---- tests 2-4
        B, evb_ = burst_table(D, D["ev"], tr, f"{R} primary")
        B_all.append(B.assign(region=R))
        T4[R] = test4(D, B, evb_, D["pool"])
        if R == "WA":
            B12, ev12 = burst_table(D, D["e12"], tr, "WA E1/E2 secondary")
            B_all.append(B12.assign(region=R))
            T4["WA_E1E2"] = test4(D, B12, ev12, D["pool12"])
        MIN_all.append(minor_table(D, tr, C, day0))
        D["B"] = B; D["evB"] = evb_
    T1 = pd.concat(T1); T1.to_csv(f"{OUT}/test1_tremor_by_distance.csv", index=False)
    B = pd.concat(B_all, ignore_index=True); B.to_csv(f"{OUT}/bursts.csv", index=False)
    MIN = pd.concat(MIN_all, ignore_index=True); MIN.to_csv(f"{OUT}/minor_episodes.csv", index=False)
    RES["test1"] = T1.to_dict("records")

    # ---- test 2 verdicts
    from scipy.stats import binom
    t2 = {}
    for lab, g in B.groupby("set"):
        e = g[g.migr_eligible]
        k = int(e.migrates.sum()) if len(e) else 0
        thr = int(binom.ppf(0.99, max(len(e), 1), 0.05))
        kx = int(e.migrates_x.sum()) if len(e) else 0
        t2[lab] = dict(n_bursts=len(g), n_eligible=len(e), n_migrating=k, binom99=thr, bursts_migrate=bool(k > thr),
                       n_migrating_crossfam=kx, bursts_migrate_crossfam=bool(kx > thr),
                       migrating=e[e.migrates][["start_PT", "n", "n_fam", "along_kmh", "along_p", "downdip_kmh", "downdip_p",
                                                "along_range_km", "downdip_range_km", "duration_h"]].round(3).to_dict("records"),
                       eligible=e[["start_PT", "n", "n_fam", "max_fam_sep_km", "along_kmh", "along_p", "downdip_kmh",
                                   "downdip_p", "along_x_kmh", "along_x_p", "downdip_x_kmh", "downdip_x_p", "duration_h"]].round(3).to_dict("records"))
        tm = g[g.trem_along_p.notna()] if "trem_along_p" in g else g.iloc[:0]
        t2[lab]["tremor_column"] = dict(n=len(tm), n_p_lt_0025=int(((tm.trem_along_p < .025) | (tm.trem_downdip_p < .025)).sum()) if len(tm) else 0,
                                        rows=tm[["start_PT", "n_trem_100km_pm6h", "trem_along_kmh", "trem_along_p",
                                                 "trem_downdip_kmh", "trem_downdip_p"]].round(3).to_dict("records") if len(tm) else [])
    RES["test2"] = t2
    log("test2: " + json.dumps(t2, indent=1, default=str))
    # ---- test 3 summary
    t3 = {}
    for lab, g in B.groupby("set"):
        big = g[g.n >= 3]
        t3[lab] = dict(n_big=len(big), duration_h_median=float(big.duration_h.median()),
                       duration_h_range=[float(big.duration_h.min()), float(big.duration_h.max())],
                       n_fam_median=float(big.n_fam.median()), sep_km_median=float(big.max_fam_sep_km.median()),
                       sep_km_max=float(big.max_fam_sep_km.max()),
                       share_single_position=float((big.max_fam_sep_km < 5).mean()),
                       share_with_trem100=float((big.n_trem_100km_pm6h > 0).mean()))
    for R, g in MIN.groupby("region"):
        t3[f"{R} minor episodes"] = dict(n=len(g), span_days_median=float(g.span_days.median()),
                                         span_days_range=[int(g.span_days.min()), int(g.span_days.max())],
                                         n_fam_median=float(g.n_fam.median()), fam_along_span_km_median=float(g.fam_along_span_km.median()),
                                         trem_along_span_km_median=float(g.trem_along_span_km.median()),
                                         n_trem_median=float(g.n_trem.median()), lfe_duration_h_median=float(g.lfe_duration_h.median()))
    RES["test3"] = t3
    log("test3: " + json.dumps(t3, indent=1))
    RES["test4"] = T4
    log("test4: " + json.dumps(T4, indent=1, default=str))

    # ---- test 5 GNSS
    t5 = {}
    est = []
    for D in regs:
        R = D["name"]
        G = gnss_setup(D, D["E"], D["gnss_lo"]); D.update(az=G["az"])
        Bb = D["B"][D["B"].n >= 3].sort_values("n", ascending=False)
        idx = np.array([int((np.datetime64(d) - G["cal"][0]).astype(int)) for d in Bb.main_day])
        okm = G["allowed"][idx]
        r, Cc, keep = gstack(G, idx[okm])
        Bk = Bb[okm].iloc[np.flatnonzero(keep)]
        k = {}
        for Rk in (5, 10, 20):
            vals = [mm_per_Nm(D, la, lo, Rk) for la, lo in zip(Bk.lat, Bk.lon)]
            k[Rk] = float(np.mean([v[0] for v in vals]))
            if Rk == 10:
                r["mesh"] = [v[1] for v in vals]
        r["mm_per_1e17Nm"] = {Rk: v * 1e17 for Rk, v in k.items()}
        r["bound_M0_Nm"] = {Rk: r["bound_mm"] / v for Rk, v in k.items()}
        r["bound_Mw"] = {Rk: float(mw(r["bound_M0_Nm"][Rk])) for Rk in k}
        r["bursts_used"] = Bk.start_PT.tolist(); r["excluded_not_allowed"] = Bb[~okm].start_PT.tolist()
        r["durations_h"] = Bk.duration_h.round(1).tolist()
        Tdur = Bk.duration_h.clip(lower=1).to_numpy() * 3600
        r["ide_Mw_median"] = float(np.median(mw(Tdur * 10 ** 12.5)))
        r["ide_Mw_range_pm05"] = [float(np.median(mw(Tdur * 10 ** 12))), float(np.median(mw(Tdur * 10 ** 13)))]
        r["ide_pred_offset_mm"] = float(np.median(Tdur * 10 ** 12.5) * k[10])
        est.append(((r["mean_mm"] - r["ctrl_mean_mm"]) / k[10], r["ctrl_sigma_mm"] / k[10]))
        # single biggest two
        singles = []
        for d, row in zip(idx[okm][:2], Bb[okm].head(2).itertuples()):
            s1, _, _ = gstack(G, np.array([d]))
            v, _ = mm_per_Nm(D, row.lat, row.lon, 10)
            s1.update(start_PT=row.start_PT, n=row.n, bound_Mw_R10=float(mw(s1["bound_mm"] / v)))
            singles.append(s1)
        r["single"] = singles
        t5[R] = r
        np.save(f"{OUT}/gnss_ctrl_{R}.npy", Cc)
        if R == "WA":   # secondary: E1/E2 big bursts added (allowed from 2009-08-15)
            G2 = gnss_setup(D, D["E"], "2009-08-15")
            B2 = pd.concat([D["B"], [b for b in B_all if (b.set == "WA E1/E2 secondary").any()][0]])
            B2 = B2[B2.n >= 3]
            idx2 = np.array([int((np.datetime64(d) - G2["cal"][0]).astype(int)) for d in B2.main_day])
            r2, _, _ = gstack(G2, idx2[G2["allowed"][idx2]])
            r2["bound_Mw_R10"] = float(mw(r2["bound_mm"] / k[10]))
            t5["WA_2009_2026_secondary"] = r2
    a = np.array([e[0] for e in est]); s = np.array([e[1] for e in est])
    w = 1 / s ** 2; am = float((a * w).sum() / w.sum()); sm = float(1 / np.sqrt(w.sum()))
    t5["combined_R10"] = dict(M0_est_Nm=am, M0_sigma_Nm=sm, z=am / sm, bound_M0_Nm=max(am, 0) + 3 * sm,
                              bound_Mw=float(mw(max(am, 0) + 3 * sm)))
    RES["test5"] = t5
    log("test5: " + json.dumps(t5, indent=1, default=str))
    json.dump(RES, open(f"{OUT}/q5_summary.json", "w"), indent=1, default=str)
    figures(regs, tr, T1, B, t5)
    log("=== Q5 done ===")


# --------------------------------------------------------------------------------------------- figures
BLUE, ORANGE, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#9a9893", "#0b0b0b", "#52514e"
SER = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]


def coast(ax, box):
    try:
        import cartopy.io.shapereader as shp
        rd = shp.Reader(os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_10m_coastline.shp"))
        for geom in rd.geometries():
            b = geom.bounds
            if b[2] < box[0] or b[0] > box[1] or b[3] < box[2] or b[1] > box[3]:
                continue
            for line in getattr(geom, "geoms", [geom]):
                x, y = line.xy; ax.plot(x, y, color=GREY, lw=0.5, zorder=0)
    except Exception as e:  # noqa
        log(f"coastline skipped: {e}")


def figures(regs, tr, T1, B, t5):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
    # ---- 1. map of the biggest bursts + tremor
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
    for ax, D in zip(axs, regs):
        F = D["F"]; Bb = D["B"].sort_values("n", ascending=False).head(4)
        box = [F.lon.min() - 1.0, F.lon.max() + 1.0, F.lat.min() - 0.7, F.lat.max() + 0.7]
        coast(ax, box)
        ax.scatter(F.lon, F.lat, s=8, c="#d6d4cf", edgecolor="none", label="all families", zorder=1)
        trt = tr.time.values.astype("int64") / 3.6e12
        for c, r in zip(SER, Bb.itertuples()):
            t0 = pd.Timestamp(r.start_utc).value / 3.6e12
            g = D["evB"][D["evB"].burst == r.burst]; t1 = g.time_utc.values.astype("int64").max() / 3.6e12
            i0, i1 = np.searchsorted(trt, [t0 - 12, t1 + 12]); tw = tr.iloc[i0:i1]
            tw = tw[hav_km(r.lat, r.lon, tw.lat.to_numpy(), tw.lon.to_numpy()) <= 100]
            ax.scatter(tw.lon, tw.lat, s=10, marker="x", c=c, lw=0.8, alpha=0.7, zorder=2)
            fp = F.loc[sorted(set(g.fid))]
            ax.scatter(fp.lon, fp.lat, s=60, c=c, edgecolor="white", lw=1.2, zorder=3,
                       label=f"{r.start_PT[:10]} PT · {r.n} LFEs, {r.n_fam} fam · {r.n_trem_100km_pm6h} tremor ≤100 km")
        ax.set_xlim(box[0], box[1]); ax.set_ylim(box[2], box[3])
        ax.set_aspect(1 / np.cos(np.deg2rad(np.mean(box[2:]))))
        ax.set_title(f"{'Southern Vancouver Island' if D['name'] == 'SVI' else 'Washington-NW'}: 4 biggest between-ETS bursts",
                     loc="left", fontsize=9.5)
        ax.legend(fontsize=7, loc="lower left", frameon=True, framealpha=0.9)
        ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    fig.text(0.01, 0.01, "dots = LFE families in the burst; × = PNSN tremor epicentres ≤ 100 km from the burst, 12 h before to 12 h after it",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{FIG}/q5_map_bursts.png", dpi=160); plt.close(fig)

    # ---- 2. along-strike vs time, 2 biggest bursts per region, with tremor
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True)
    for row, D in zip(axs, regs):
        Bb = D["B"].sort_values("n", ascending=False).head(2)
        trt = tr.time.values.astype("int64") / 3.6e12
        for ax, r in zip(row, Bb.itertuples()):
            g = D["evB"][D["evB"].burst == r.burst]
            t0 = pd.Timestamp(r.start_utc).value / 3.6e12; t1 = g.time_utc.values.astype("int64").max() / 3.6e12
            i0, i1 = np.searchsorted(trt, [t0 - 24, t1 + 24]); tw = tr.iloc[i0:i1]
            dk = hav_km(r.lat, r.lon, tw.lat.to_numpy(), tw.lon.to_numpy()); tw = tw[dk <= 100]
            ta, _ = proj(D["fr"], tw.lat.to_numpy(), tw.lon.to_numpy())
            ax.scatter((tw.time.values.astype("int64") / 3.6e12) - t0, ta, s=12, marker="x", c=GREY, lw=0.8,
                       label=f"PNSN tremor ≤ 100 km ({len(tw)})")
            sc = ax.scatter(g.time_utc.values.astype("int64") / 3.6e12 - t0, g.along_km, s=40, c=g.dep, cmap="viridis_r",
                            vmin=30, vmax=48, edgecolor="white", lw=0.8, zorder=3, label=f"LFEs ({len(g)}, {r.n_fam} families)")
            ax.axvspan(0, t1 - t0, color=BLUE, alpha=0.06, lw=0)
            ax.set_xlabel(f"hours from burst start ({r.start_PT} PT)"); ax.set_ylabel("along-strike (km, + = NW)")
            ax.set_title(f"{D['name']} burst {r.start_PT[:10]} PT · groups {r.groups}", loc="left", fontsize=9)
            ax.legend(fontsize=7.5, frameon=False, loc="best")
            ax.set_xlim(-24, t1 - t0 + 24)
            lo_, hi_ = ax.get_ylim(); mid_ = (lo_ + hi_) / 2
            if hi_ - lo_ < 20:
                ax.set_ylim(mid_ - 10, mid_ + 10)
    cb = fig.colorbar(sc, ax=axs, shrink=0.6, pad=0.01); cb.set_label("family depth (km)")
    fig.suptitle("Biggest between-ETS bursts: LFEs (dots) and PNSN tremor ≤ 100 km (×), ± 24 h; shaded = burst", fontsize=9.5)
    fig.savefig(f"{FIG}/q5_along_strike_time.png", dpi=160, bbox_inches="tight"); plt.close(fig)

    # ---- 3. tremor coincidence by distance, both regions (burst level)
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4), sharey=True)
    for ax, R in zip(axs, ["SVI", "WA"]):
        x = T1[(T1.region == R) & T1.tested].copy()
        x["o"] = x.bin.map({"0-30 km": 0, "30-60 km": 1, "60-100 km": 2}); x = x.sort_values("o")
        lab = [f"{b}\n({'broad' if 'broad' in s else 'BE'})" for b, s in zip(x.bin, x.set)]
        xi = np.arange(len(x))
        ax.bar(xi - 0.18, x.burst_share, 0.34, color=BLUE, label="between-ETS bursts")
        ax.bar(xi + 0.18, x.burst_ctrl_median, 0.34, color="#c9c7c1", label="random-day control (median)")
        ax.errorbar(xi + 0.18, x.burst_ctrl_median, yerr=[x.burst_ctrl_median - x.burst_ctrl_p01, x.burst_ctrl_p99 - x.burst_ctrl_median],
                    fmt="none", ecolor=INK2, capsize=3, lw=1)
        for i, r in enumerate(x.itertuples()):
            ax.text(i - 0.18, r.burst_share + 0.015, f"{100 * r.burst_share:.0f} %\np={r.burst_p:.3f}", ha="center", fontsize=7.5, color=INK)
            ax.plot(i - 0.18, r.event_share, marker="D", ms=5, color=INK, zorder=4, ls="none",
                    label="share of events" if i == 0 else None)
        ax.set_xticks(xi); ax.set_xticklabels(lab)
        n = int(x.n_bursts.iloc[0]) if len(x) else 0
        ax.set_title(f"{R}: bursts with PNSN tremor the same day, by distance", loc="left", fontsize=9)
        ax.set_ylabel("share of bursts with ≥ 1 epicentre")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.text(0.01, 0.01, "BE = no tremor ≤ 30 km that day (by definition), so 0–30 km is tested on the broad set (any local tremor). "
             "Whiskers: 1st–99th percentile of 1000 random-day draws.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(f"{FIG}/q5_tremor_by_distance.png", dpi=160); plt.close(fig)

    # ---- 4. GNSS stack
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, R in zip(axs, ["SVI", "WA"]):
        C = np.load(f"{OUT}/gnss_ctrl_{R}.npy"); r = t5[R]
        ax.hist(C, 40, color="#c9c7c1")
        ax.axvline(r["mean_mm"], color=BLUE, lw=2, label=f"biggest bursts ({r['n_days']}): {r['mean_mm']:+.2f} mm")
        ax.axvline(r["ctrl_mean_mm"] + 3 * r["ctrl_sigma_mm"], color=INK2, ls="--", lw=1,
                   label=f"control mean + 3σ  (bound {r['bound_mm']:.2f} mm ≈ Mw {r['bound_Mw'][10]:.1f})")
        ax.set_xlabel("mean GNSS offset over ~6 d, ETS direction (mm)"); ax.set_ylabel("random-day sets")
        ax.set_title(f"{R}: GNSS on the biggest bursts", loc="left", fontsize=9); ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/q5_gnss_bursts.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()
