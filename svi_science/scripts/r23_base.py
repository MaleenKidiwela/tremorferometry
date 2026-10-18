#!/usr/bin/env python
"""r23_base (2026-10-11) -- shared foundation for the R23 reruns on the 23-year SVI LFE record.
Rules: svi_science/results/r23/PLAN.md §0-1 (pre-registered). Read-only on every input; writes svi_science/results/r23/tables.

    python svi_science/scripts/r23_base.py
"""
from __future__ import annotations
import os, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/r23/tables"
S0 = f"{ROOT}/svi_science/results/s0"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
BW = f"{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet"
BWR = f"{ROOT}/data/catalogue/svi_rates_v1_4_backward_famday.parquet"
FWR = f"{ROOT}/lfe_extension/results/refine/rates/famday.parquet"
ENV = f"{ROOT}/svi_science/results/backward/p4/envflag"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
RULE = dict(frac=0.05, gap=5, minlen=5, major=0.15)       # episode_rates RULE (PLAN Revision 1), unchanged
ENV_C, ENV_K = 0.42, 12                                    # frozen env flag (ENVFLAG.md)
PNSN_T0 = "2010-01-06"; PNSN_T1 = "2026-05-27"
SEGS = [("CN4f", "2003-01-01", "2003-01-31"), ("ERA", "2003-02-01", "2006-06-30"), ("B2f", "2006-07-01", "2007-11-12"),
        ("B3f", "2007-11-13", "2009-08-26"), ("S15f", "2009-08-27", "2013-12-31"), ("FWD", "2014-01-01", "2026-12-31")]
DEEP_GROUPS = (1, 4, 10, 12)
LOGF = None
RES = {}


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()


def seg_of(dates):
    d = pd.Series(np.asarray(dates, dtype=str))
    out = np.full(len(d), "", dtype=object)
    for nm, a, b in SEGS:
        out[((d >= a) & (d <= b)).to_numpy()] = nm
    return out


def hav_km(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def find_episodes(Xe, nfam, t0, t1, rule):
    """Copied unchanged from lfe_extension_wa/scripts/episode_rates.py (kind label set by caller)."""
    days = pd.date_range(t0, t1, freq="D")
    ds = days.strftime("%Y-%m-%d")
    nexp = Xe[Xe.n30 > 0].groupby("date").fid.nunique().reindex(ds, fill_value=0).to_numpy()
    q = nexp >= rule["frac"] * nfam
    idx = np.flatnonzero(q)
    eps = []
    if len(idx):
        s = idx[0]; prev = idx[0]
        for i in idx[1:]:
            if i - prev - 1 <= rule["gap"]:
                prev = i; continue
            eps.append((s, prev)); s = i; prev = i
        eps.append((s, prev))
    rows = []
    for s, e in eps:
        L = e - s + 1
        if L < rule["minlen"]:
            continue
        fed = int(nexp[s:e + 1].sum())
        rows.append(dict(start=ds[s], end=ds[e], span_days=L, qual_days=int(q[s:e + 1].sum()), family_exposure_days=fed,
                         max_families_exposed=int(nexp[s:e + 1].max()), major_threshold=rule["major"] * nfam * 14,
                         major=fed >= rule["major"] * nfam * 14))
    return pd.DataFrame(rows)


def two_class_split(v):
    """q6_build.two_class_split, unchanged."""
    v = np.sort(np.asarray(v, float)); best = (np.inf, None)
    for i in range(1, len(v)):
        a, b = v[:i], v[i:]
        w = ((a - a.mean()) ** 2).sum() + ((b - b.mean()) ** 2).sum()
        if w < best[0]:
            best = (w, (v[i - 1] + v[i]) / 2)
    return best[1]


def overlaps(a0, a1, B):
    return ((B.start <= a1) & (B.end >= a0)).any()


# ------------------------------------------------------------------------------------------- events
def events(F):
    cols = ["date", "series", "fam", "t", "P", "detector", "net_set", "era", "primary", "keep_xfam8", "burst_day", "n_z2",
            "n_sta_day", "quiet"]
    B = pd.read_parquet(BW, columns=cols)
    # acceptance against RESULTS_BACKWARD (all families)
    sel = B.primary & B.keep_xfam8 & (B.P >= np.where(B.detector == "era", 0.90, 0.75))
    yr = B.date.str[:4]
    acc_d = B[sel & (B.series == "data")].groupby(yr).size()
    acc_c = B[sel & (B.series == "c1")].groupby(yr).size()
    RES["acc_backward_by_year_data"] = acc_d.to_dict(); RES["acc_backward_by_year_c1"] = acc_c.to_dict()
    RES["acc_backward_total"] = [int(acc_d.sum()), int(acc_c.sum())]
    # Revision 1: RESULTS_BACKWARD's 51,539 is the all-P>=0.75 total; under the mandated era P*>=0.90 rule the reference is
    # its "data at era P* 0.90" row for 2003-2005 (3,730 / 10,891 / 6,620; 2003 + 3 CN4 frozen Jan rows) and the data row 2007-13.
    ref = {"2004": 10891, "2005": 6620, "2007": 452, "2008": 2097, "2009": 1403, "2010": 1689, "2011": 1709, "2012": 3100, "2013": 1540}
    RES["acc_backward_pass"] = bool(all(int(acc_d.get(y, -1)) == v for y, v in ref.items()) and int(acc_d.get("2003", 0)) in (3730, 3733))
    log(f"backward counted-once P-rule totals data {acc_d.sum()} / control {acc_c.sum()} (per-year reference = P* row): "
        f"{'PASS' if RES['acc_backward_pass'] else 'FAIL'}")
    log(f"per year data {acc_d.to_dict()}")
    B = B[sel].copy()
    RES["backward_extra_families_events"] = int((~B.fam.isin(F.fid) & (B.series == "data")).sum())
    B = B[B.fam.isin(F.fid)]
    B["time_utc"] = pd.to_datetime(B.t, unit="s").astype("datetime64[ns]")
    B["src"] = "bw"
    B = B.drop(columns=["t", "primary", "keep_xfam8"])
    fc = ["date", "series", "fam", "P_same", "keep_xfam8", "time_utc", "n_sta_day", "quiet", "run"]
    Fw = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=fc, filters=[("P_same", ">=", 0.75)])
    Fw = Fw[(Fw.run.astype(str) == "prod") & Fw.keep_xfam8].drop(columns=["run", "keep_xfam8"]).rename(columns={"P_same": "P"})
    Fw["series"] = Fw.series.astype(str)
    Fw["detector"] = "forward"; Fw["net_set"] = "FWD"; Fw["era"] = "FWD"; Fw["src"] = "fw"
    Fw["burst_day"] = pd.NA; Fw["n_z2"] = pd.NA
    Fw["time_utc"] = Fw.time_utc.astype("datetime64[ns]")
    E = pd.concat([B, Fw[B.columns]], ignore_index=True)
    E["seg"] = seg_of(E.date)
    E["burst_day"] = E.burst_day.astype("boolean"); E["n_z2"] = E.n_z2.astype("Int16")
    E = E.sort_values("time_utc").reset_index(drop=True)
    # Bostock-equivalent weight per event: k of the (family, date) row of the rate tables
    kb = pd.read_parquet(BWR, columns=["fam", "date", "primary", "k"]); kb = kb[kb.primary].drop(columns="primary")
    kf = pd.read_parquet(FWR, columns=["fam", "date", "k"])
    K = pd.concat([kb, kf], ignore_index=True).drop_duplicates(["fam", "date"])
    E = E.merge(K, on=["fam", "date"], how="left")
    RES["events_missing_k"] = int(E.k.isna().sum())
    g = E.groupby(["seg", "series"]).size().unstack(fill_value=0)
    log("events (115 families, single-repeat rule) by segment:\n" + g.to_string())
    E.to_parquet(f"{OUT}/events_2003_2026.parquet", index=False)
    return E


# ------------------------------------------------------------------------------------------- observable family-days
def famdays(F):
    b = pd.read_parquet(BWR, columns=["fam", "date", "primary", "m_data", "k", "rate", "excess", "n_sta_day", "burst_day"])
    b = b[b.primary & (b.m_data > 0) & b.fam.isin(F.fid)].drop(columns=["primary", "m_data"])
    f = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable", "n_sta_day"])
    f = f[(f.run == "prod") & f.observable].rename(columns={"fid": "fam"}).drop(columns=["run", "observable"])
    D = pd.concat([b.assign(src="bw"), f.assign(src="fw")], ignore_index=True)
    D["seg"] = seg_of(D.date)
    D.to_parquet(f"{OUT}/famdays_observable.parquet", index=False)
    log(f"observable family-days: backward {len(b)}, forward {len(f)}")
    return D


# ------------------------------------------------------------------------------------------- tremor
def tremor(F):
    parts = []
    for y in range(2003, 2014):
        w = pd.read_parquet(f"{ENV}/win/env_{y}.parquet", columns=["date", "series", "sset", "win", "coh", "s3", "lat", "lon", "dfam_km"])
        w = w[(w.sset == "CN4") & (w.series == "data") & (w.coh >= ENV_C) & (w.s3 >= ENV_C) & (w.dfam_km <= 30)]
        w = w.assign(time=pd.to_datetime(w.date) + pd.to_timedelta(w.win.astype(int) * 150 + 150, unit="s"), src="env")
        parts.append(w[["time", "date", "lat", "lon", "src"]])
        # scrambled-control windows kept for a count check only
    env = pd.concat(parts, ignore_index=True)
    T = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce"); T = T.dropna(subset=["time"])
    T["time"] = T.time.astype("datetime64[ns]")
    T["date"] = T.time.dt.strftime("%Y-%m-%d")
    T = T[(T.date >= PNSN_T0) & (T.date <= PNSN_T1)]
    lat0, lat1, lon0, lon1 = F.lat.min() - 1.0, F.lat.max() + 1.0, F.lon.min() - 1.4, F.lon.max() + 1.4
    T = T[T.lat.between(lat0, lat1) & T.lon.between(lon0, lon1)].assign(src="pnsn")
    TR = pd.concat([env, T[["time", "date", "lat", "lon", "src"]]], ignore_index=True)
    TR["lat"] = TR.lat.astype(float); TR["lon"] = TR.lon.astype(float)
    # nearest-family distance and per-family counts
    la, lo = TR.lat.to_numpy(), TR.lon.to_numpy()
    dmin = np.full(len(TR), np.inf); dmin_deep = np.full(len(TR), np.inf)
    rows = []
    deep = set(F.fid[F.group.isin(DEEP_GROUPS)])
    for r in F.itertuples():
        d = hav_km(la, lo, r.lat, r.lon)
        dmin = np.minimum(dmin, d)
        if r.fid in deep:
            dmin_deep = np.minimum(dmin_deep, d)
        for src in ("env", "pnsn"):
            m = (TR.src.to_numpy() == src)
            sub = TR.date.to_numpy()
            c30 = pd.Series(sub[m & (d <= 30)]).value_counts()
            c100 = pd.Series(sub[m & (d > 30) & (d <= 100)]).value_counts()
            c15 = pd.Series(sub[m & (d <= 15)]).value_counts()
            idx = c30.index.union(c100.index)
            for dd in idx:
                rows.append((int(r.fid), dd, src, int(c30.get(dd, 0)), int(c100.get(dd, 0)), int(c15.get(dd, 0))))
    TR["dmin_fam_km"] = dmin; TR["dmin_deep_km"] = dmin_deep
    TR.to_parquet(f"{OUT}/tremor_points.parquet", index=False)
    X = pd.DataFrame(rows, columns=["fam", "date", "src", "n30", "n30_100", "n15"])
    X.to_parquet(f"{OUT}/famday_tremor.parquet", index=False)
    # acceptance: PNSN n30 must equal exposure_svi (the S7 input)
    XS = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    M = X[(X.src == "pnsn") & (X.n30 > 0)].merge(XS.rename(columns={"fid": "fam", "n30": "n30_ref"}), on=["fam", "date"], how="outer")
    M = M[(M.date >= PNSN_T0) & (M.date <= PNSN_T1)]
    mism = int((M.n30.fillna(0) != M.n30_ref.fillna(0)).sum())
    RES["acc_pnsn_n30_mismatch_rows"] = mism; RES["acc_pnsn_n30_rows"] = int(len(M))
    log(f"tremor points: env {int((TR.src=='env').sum())}, pnsn {int((TR.src=='pnsn').sum())}; PNSN n30 vs exposure_svi "
        f"mismatching family-days {mism} of {len(M)}")
    return TR, X


# ------------------------------------------------------------------------------------------- episodes
def episodes(F, X, E):
    nf = len(F)
    q = pd.read_csv(f"{ENV}/tables/quiet_flag_env.csv", usecols=["date", "N_CN4"])
    act = set(q.date[q.N_CN4 >= ENV_K])
    Xe = X[(X.src == "env") & (X.n30 > 0)].rename(columns={"fam": "fid"})[["fid", "date", "n30"]]
    Xb = Xe[Xe.date.isin(act)]
    Eb = pd.read_csv(f"{ROOT}/svi_science/results/b_ets_record/episodes_b.csv")
    P10 = Eb[(Eb.region == "SVI") & (Eb.start < "2014-01-01")].copy()       # PNSN extension 2010-2013
    P10["major"] = P10.major.astype(bool)
    bost = pd.read_csv(f"{ROOT}/svi_science/results/p3_bostock/episodes_bostock_final.csv")
    val = {}
    lists = {}
    for var, XX in (("a", Xe), ("b", Xb)):
        v = find_episodes(XX, nf, PNSN_T0, "2013-12-31", RULE)
        Mv = v[v.major]; Mp = P10[P10.major]
        rec = sum(overlaps(r.start, r.end, Mv) for r in Mp.itertuples())
        prec = np.mean([overlaps(r.start, r.end, P10) for r in Mv.itertuples()]) if len(Mv) else np.nan
        pre = find_episodes(XX, nf, "2003-01-01", "2010-01-05", RULE)
        bw = bost[bost.start < "2010-01-01"]
        nb = sum(overlaps(r.start, r.end, pre[pre.major]) for r in bw.itertuples())
        dd = []
        for r in Mp.itertuples():
            h = Mv[(Mv.start <= r.end) & (Mv.end >= r.start)]
            if len(h):
                dd.append(dict(pnsn=r.episode, d_start=(pd.Timestamp(h.start.iloc[0]) - pd.Timestamp(r.start)).days,
                               d_end=(pd.Timestamp(h.end.iloc[-1]) - pd.Timestamp(r.end)).days, n_env_majors=len(h)))
        n_se = sum(abs(x["d_start"]) <= 7 and abs(x["d_end"]) <= 7 for x in dd)   # Revision 1 (start/end criterion)
        ok = bool(rec >= 6 and prec >= 0.8 and nb == 6 and n_se >= 6)
        val[var] = dict(env_majors_2010_13=int(len(Mv)), env_episodes_2010_13=int(len(v)), pnsn_majors=int(len(Mp)),
                        recall=int(rec), precision=float(prec), bostock_windows_hit=int(nb), n_start_end_within_7d=int(n_se), pass_=ok, start_end_diff=dd,
                        pre2010_episodes=int(len(pre)), pre2010_majors=int(pre.major.sum()))
        lists[var] = pre
        v.assign(variant=var).to_csv(f"{OUT}/env_episodes_2010_2013_variant_{var}.csv", index=False)
        log(f"env variant {var}: {json.dumps({k: x for k, x in val[var].items() if k != 'start_end_diff'})}")
        log(f"   start/end diffs vs PNSN majors: {dd}")
    RES["env_episode_validation"] = val
    if val["a"]["pass_"]:
        choice = "a"
    elif val["b"]["pass_"]:
        choice = "b"
    else:
        choice = "fallback"
    RES["env_episode_choice"] = choice
    log(f"pre-2010 episode list choice: {choice}")
    if choice in ("a", "b"):
        pre = lists[choice].assign(kind=f"env_{choice}", validated=True)
    else:
        bw = bost[bost.start < "2010-01-01"][["start", "end"]].assign(major=True, kind="bostock", validated=False)
        mb = lists["b"]
        mb = mb[[not overlaps(r.start, r.end, bw) for r in mb.itertuples()]].assign(major=False, kind="env_b", validated=False)
        pre = pd.concat([bw, mb], ignore_index=True).sort_values("start")
    pre["episode"] = pre.start.str[:7]
    P10 = P10.assign(kind="pnsn", validated=True)
    Ef = pd.read_csv(f"{EPR}/episodes_svi_final.csv").assign(kind="pnsn", validated=True)
    cols = ["episode", "start", "end", "span_days", "family_exposure_days", "major", "kind", "validated"]
    EP = pd.concat([pre.reindex(columns=cols), P10.reindex(columns=cols), Ef.reindex(columns=cols)], ignore_index=True)
    EP = EP.sort_values("start").reset_index(drop=True)
    # unique episode labels
    cnt = {}
    lab = []
    for e in EP.episode:
        cnt[e] = cnt.get(e, 0) + 1
        lab.append(e if cnt[e] == 1 else f"{e}b{cnt[e]-1}")
    EP["episode"] = lab
    EP["seg"] = seg_of(EP.start)
    # sizes: S_rate (115 families, primary rate), rate_kc, Bostock counts, forward S3/B columns
    rb = pd.read_parquet(BWR, columns=["fam", "date", "primary", "rate", "rate_kc"]); rb = rb[rb.primary & rb.fam.isin(F.fid)]
    rf = pd.read_parquet(FWR, columns=["fam", "date", "rate"]); rf = rf[rf.fam.isin(F.fid)].assign(rate_kc=lambda x: x.rate)
    R = pd.concat([rb.drop(columns="primary"), rf], ignore_index=True)
    Rd = R.groupby("date")[["rate", "rate_kc"]].sum()
    ev = E[E.series == "data"]
    evd = ev.groupby("date").size()
    bostev = pd.read_parquet(f"{ROOT}/svi_science/results/p3_bostock/bostock_events_mapped.parquet", columns=["fid", "date"])
    bostev = bostev[bostev.fid.isin(F.fid)].groupby("date").size()
    s_rate, s_kc, n_once, n_bost, t5 = [], [], [], [], []
    for r in EP.itertuples():
        dd = pd.date_range(r.start, r.end).strftime("%Y-%m-%d")
        s_rate.append(float(Rd.rate.reindex(dd).fillna(0).sum())); s_kc.append(float(Rd.rate_kc.reindex(dd).fillna(0).sum()))
        n_once.append(int(evd.reindex(dd).fillna(0).sum())); n_bost.append(int(bostev.reindex(dd).fillna(0).sum()))
        x = ev[(ev.date >= r.start) & (ev.date <= r.end)].time_utc
        t5.append(x.quantile(0.05) if len(x) else pd.NaT)
    EP["S_rate"] = s_rate; EP["S_rate_kc"] = s_kc; EP["n_once"] = n_once; EP["n_bostock"] = n_bost
    EP["onset_t5_utc"] = t5
    EP["onset_t5_PT"] = pd.to_datetime(EP.onset_t5_utc).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")
    Bb = Eb[Eb.region == "SVI"][["episode", "start", "big", "s3_n_headline", "gnss_D_mm", "gnss_SN"]].rename(columns={"big": "big_B", "episode": "episode_B"})
    EP = EP.merge(Bb, on="start", how="left")
    s3big = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["big"]["SVI"]
    EP["big_S3"] = EP.episode.isin(s3big) & (EP.start >= "2014-01-01")
    M = EP[EP.major]
    thr = two_class_split(np.log10(np.maximum(M.S_rate.to_numpy(), 1.0)))
    EP["big"] = EP.major & (np.log10(np.maximum(EP.S_rate, 1.0)) > thr)
    RES["big_split_log10_S_rate"] = float(thr); RES["big_split_S_rate"] = float(10 ** thr)
    RES["n_episodes"] = int(len(EP)); RES["n_majors"] = int(EP.major.sum()); RES["n_big"] = int(EP.big.sum())
    RES["n_big_2003_13"] = int(EP[EP.start < "2014-01-01"].big.sum()); RES["n_big_2014_26"] = int(EP[EP.start >= "2014-01-01"].big.sum())
    EP.to_csv(f"{OUT}/episodes_2003_2026.csv", index=False)
    log(f"episodes 2003-2026: {len(EP)} ({EP.major.sum()} majors, {EP.big.sum()} big; split S_rate {10**thr:.0f})")
    log("majors:\n" + EP[EP.major][["episode", "start", "end", "kind", "seg", "family_exposure_days", "S_rate", "n_once", "n_bostock",
                                     "big", "big_B", "big_S3", "gnss_D_mm"]].round(1).to_string(index=False))
    return EP


def main():
    global LOGF
    os.makedirs(OUT, exist_ok=True)
    LOGF = open(f"{OUT}/../r23_base.log", "w")
    F = pd.read_csv(f"{S0}/families.csv")
    E = events(F)
    D = famdays(F)
    TR, X = tremor(F)
    EP = episodes(F, X, E)
    json.dump(RES, open(f"{OUT}/base_summary.json", "w"), indent=1, default=str)
    log("done")


if __name__ == "__main__":
    main()
