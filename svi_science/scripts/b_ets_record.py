#!/usr/bin/env python
"""B: ETS over the long record, southern Vancouver Island (SVI) and Washington-NW.
Frozen rules: svi_science/results/b_ets_record/PLAN_B.md.  Writes only under results/b_ets_record and figures/b_*.png.
Never deletes.  Times in Pacific time (PT).

  python b_ets_record.py build     # episodes + footprint tremor + LFE/GNSS columns -> episodes_b.csv
  python b_ets_record.py tests     # Q1-Q5 -> intervals.csv, pairs.csv, seasonal_monthly.csv, tests.json
  python b_ets_record.py figures
"""
from __future__ import annotations
import os, sys, json, glob, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/b_ets_record"
FIG = f"{ROOT}/svi_science/figures"
sys.path.insert(0, f"{ROOT}/lfe_extension_wa/scripts")
import episode_rates as ER  # noqa: E402

TZ = "America/Los_Angeles"
KM = 111.195
PAD = 3.0
NPERM = 10000
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
COV = {"SVI": ("2010-01-06", "2026-05-27"), "NW": ("2009-08-06", "2026-05-27")}
NFAM = {"SVI": 115, "NW": 132}
CUT_ERA = pd.Timestamp("2018-01-01")
COL = {"SVI": "#2a78d6", "NW": "#eb6834"}


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(f"{OUT}/b.log", "a") as f:
        f.write(pd.Timestamp.now(tz=TZ).strftime("%Y-%m-%d %H:%M PT ") + s + "\n")


def to_pt(t):
    t = pd.to_datetime(pd.Series(t)).astype("datetime64[ns]")
    return t.dt.tz_localize("UTC").dt.tz_convert(TZ).dt.tz_localize(None)


def padded_windows(E):           # identical to s3_build.padded_windows
    E = E.sort_values("start").reset_index(drop=True).copy()
    s = pd.to_datetime(E.start); e = pd.to_datetime(E.end) + pd.Timedelta(days=1)
    lo, hi = [], []
    for i in range(len(E)):
        gp = (s[i] - e[i - 1]).total_seconds() / 86400 if i > 0 else 99
        gn = (s[i + 1] - e[i]).total_seconds() / 86400 if i < len(E) - 1 else 99
        lo.append(s[i] - pd.Timedelta(days=min(PAD, max(gp, 0) / 2)))
        hi.append(e[i] + pd.Timedelta(days=min(PAD, max(gn, 0) / 2)))
    E["w0"] = lo; E["w1"] = hi
    return E


def name_eps(E):
    E["episode"] = E.start.str[:7]
    E["episode"] = E.groupby("episode").cumcount().map(lambda k: "" if k == 0 else f"b{k}").radd(E.episode)
    return E


# ------------------------------------------------------------------------------------------------ build
def load_tremor():
    P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon", "depth"])
    P9 = pd.read_csv(f"{ROOT}/lfe_array_discovery/results/ab_royer/tremor/pnsn_tremor_2009.csv", skipinitialspace=True)
    P9.columns = [c.strip() for c in P9.columns]
    P9 = P9.rename(columns={"starttime": "time"})[["time", "lat", "lon", "depth"]]
    P9["time"] = P9.time.astype(str).str.strip()
    P9 = P9[P9.time < "2010-01-06"]
    P = pd.concat([P9, P], ignore_index=True)
    P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
    for c in ("lat", "lon", "depth"):
        P[c] = pd.to_numeric(P[c], errors="coerce")
    return P.sort_values("time").reset_index(drop=True)


def min_dist(lat, lon, F):
    d = np.full(len(lat), np.inf)
    for la, lo in zip(F.lat.values, F.lon.values):
        d = np.minimum(d, np.hypot((lon - lo) * KM * np.cos(np.radians(la)), (lat - la) * KM))
    return d


def split2(x):
    """optimal 1-D two-class split (minimum within-class SS, exhaustive) -> lower edge of the upper class."""
    xs = np.sort(x); best = (np.inf, None)
    for k in range(1, len(xs)):
        sse = ((xs[:k] - xs[:k].mean()) ** 2).sum() + ((xs[k:] - xs[k:].mean()) ** 2).sum()
        if sse < best[0]:
            best = (sse, k)
    return xs[best[1]]


def build():
    info = {}
    Fs = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")[["fid", "lat", "lon"]]
    Fw = ER.wa_families()
    # --- SVI extension 2010-2013 with the identical rule
    Xs = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    nsv = len(ER.svi_families())
    Ext = ER.find_episodes(Xs, nsv, COV["SVI"][0], COV["SVI"][1], ER.RULE)
    Efs = pd.read_csv(f"{EPR}/episodes_svi_final.csv")
    late = Ext[Ext.start >= "2014-01-01"].reset_index(drop=True)
    same = (len(late) == len(Efs) and (late.start.values == Efs.start.values).all() and (late.end.values == Efs.end.values).all()
            and (late.family_exposure_days.values == Efs.family_exposure_days.values).all())
    cross = Ext[(Ext.start < "2014-01-01") & (Ext.end >= "2014-01-01")]
    info["svi_ext_check_pass"] = bool(same)   # Revision 1: a crossing episode is allowed (no double listing)
    info["svi_ext_crossing"] = cross[["start", "end", "family_exposure_days"]].to_dict("records")
    info["svi_n_fam_rule"] = nsv
    log(f"SVI extension reproduces 2014+ final list: {same}; episodes crossing 2014-01-01: {len(cross)}")
    early = Ext[Ext.start < "2014-01-01"].copy()
    early.insert(0, "region", "SVI")
    early["in_frozen_major"] = np.nan
    early["source_list"] = "svi_extension_2010_2013"
    Efs["source_list"] = "episodes_svi_final"
    Es = pd.concat([early, Efs], ignore_index=True) if info["svi_ext_check_pass"] else Efs
    Es = name_eps(Es.drop(columns="episode", errors="ignore").sort_values("start").reset_index(drop=True))
    Ew = pd.read_csv(f"{EPR}/episodes_nw_final.csv"); Ew["source_list"] = "episodes_nw_final"
    info["n_eps"] = {"SVI": int(len(Es)), "NW": int(len(Ew)), "SVI_major": int(Es.major.sum()), "NW_major": int(Ew.major.sum())}

    # --- tremor footprints
    P = load_tremor()
    box = (P.lat > 46.5) & (P.lat < 49.8) & (P.lon > -126) & (P.lon < -121)
    P = P[box].reset_index(drop=True)
    ds = min_dist(P.lat.values, P.lon.values, Fs); dw = min_dist(P.lat.values, P.lon.values, Fw)
    P["d_svi"] = ds; P["d_nw"] = dw
    for r in (10, 30):
        reg = np.where((ds <= r) | (dw <= r), np.where(ds <= dw, "SVI", "NW"), "")
        P[f"reg{r}"] = reg
    P["both10"] = (ds <= 10) & (dw <= 10)
    info["tremor_rows_box"] = int(len(P))
    info["tremor_both10"] = int(P.both10.sum())
    info["tremor_reg10"] = P.reg10.value_counts().to_dict()
    P["tpt"] = to_pt(P.time).values
    P.to_parquet(f"{OUT}/tremor_footprints.parquet")

    rows = []
    for reg, E in (("SVI", Es), ("NW", Ew)):
        E = padded_windows(E)
        T10 = P[P.reg10 == reg].time.values; T30 = P[P.reg30 == reg].time.values
        for r in E.itertuples():
            a, b = np.datetime64(r.w0), np.datetime64(r.w1)
            t = T10[(T10 >= a) & (T10 < b)]
            n30 = int(((T30 >= a) & (T30 < b)).sum())
            d = dict(region=reg, episode=r.episode, kind=r.kind, start=r.start, end=r.end, span_days=r.span_days,
                     family_exposure_days=r.family_exposure_days, max_families_exposed=r.max_families_exposed,
                     major=bool(r.major), in_frozen_major=r.in_frozen_major, source_list=r.source_list,
                     w0=r.w0, w1=r.w1, n_trem10=len(t) if r.kind != "wecc" else np.nan, n_trem30=n30 if r.kind != "wecc" else np.nan)
            if len(t) >= 5 and r.kind != "wecc":
                tq = pd.Series(np.sort(t))
                q = lambda p: tq.iloc[int(np.floor(p * (len(tq) - 1)))]
                d.update(t5=q(0.05), t10=q(0.10), t_first=tq.iloc[0], t95=q(0.95))
            else:
                st = pd.Timestamp(r.start) + pd.Timedelta(hours=8)    # WECC/no-tremor: list start day, 00:00 PT approx
                d.update(t5=st, t10=st, t_first=st, t95=pd.Timestamp(r.end) + pd.Timedelta(hours=8 + 24))
            rows.append(d)
    B = pd.DataFrame(rows)
    for c in ("t5", "t10", "t_first", "t95"):
        B[c + "_PT"] = to_pt(B[c]).values
    B["onset_PT"] = B.t5_PT; B["end_PT"] = B.t95_PT
    B["duration_d"] = (B.end_PT - B.onset_PT).dt.total_seconds() / 86400
    B["onset_year"] = B.onset_PT.dt.year + (B.onset_PT.dt.dayofyear - 1) / 365.25
    B["onset_doy"] = B.onset_PT.dt.dayofyear

    # --- LFE columns: SVI from S3 (majors), NW recomputed for every era
    S3 = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_all.csv")
    s3 = S3.set_index(["region", "episode"])
    for c in ("n_headline", "n_c1_075", "n_raw075", "n_once050"):
        B[f"s3_{c}"] = [s3[c].get((g, e), np.nan) for g, e in zip(B.region, B.episode)]
    B["s3_big"] = [bool(s3["big"].get((g, e), False)) for g, e in zip(B.region, B.episode)]
    B["lfe_svi_once075"] = np.where(B.region == "SVI", B.s3_n_headline, np.nan)
    B["lfe_svi_c1_075"] = np.where(B.region == "SVI", B.s3_n_c1_075, np.nan)

    C = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet",
                        columns=["det_id", "time", "region", "family", "family_uid", "selected_v2"],
                        filters=[("region", "==", "NW"), ("selected_v2", "==", True)])
    ids = set(C.det_id.values)
    fl = []
    for f in sorted(glob.glob(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_20*.parquet")):
        q = pd.read_parquet(f, columns=["det_id", "series", "dup_D"])
        q = q[(q.series == "data")]
        fl.append(q[q.det_id.isin(C.det_id)][["det_id", "dup_D"]])
    fl = pd.concat(fl).drop_duplicates("det_id")
    C = C.merge(fl, on="det_id", how="left")
    info["nw_selected_rows"] = int(len(C)); info["nw_selected_without_flag"] = int(C.dup_D.isna().sum())
    C["once"] = ~C.dup_D.fillna(False).astype(bool)
    V = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_3.csv", usecols=["family_uid", "nw_source_group"], low_memory=False)
    C["src"] = C.family_uid.map(V.drop_duplicates("family_uid").set_index("family_uid").nw_source_group)
    C = C.sort_values("time")
    tC = C.time.values.astype("datetime64[ns]")
    n_once, n_raw, n_src = [], [], []
    for r in B.itertuples():
        if r.region != "NW":
            n_once.append(np.nan); n_raw.append(np.nan); n_src.append(np.nan); continue
        i0, i1 = np.searchsorted(tC, np.datetime64(r.w0)), np.searchsorted(tC, np.datetime64(r.w1))
        sub = C.iloc[i0:i1]
        n_raw.append(len(sub)); so = sub[sub.once]; n_once.append(len(so))
        vc = so.src.value_counts(); n_src.append(int((vc >= 3).sum()))
    B["lfe_nw_once075"] = n_once; B["lfe_nw_raw075"] = n_raw; B["lfe_nw_n_sources"] = n_src
    B["nw_era"] = [None if g != "NW" else ER.wa_era(str(pd.Timestamp(s).date())) for g, s in zip(B.region, B.start)]
    chk = B[(B.region == "NW") & B.s3_n_headline.notna()]
    rel = (chk.lfe_nw_once075 - chk.s3_n_headline).abs() / chk.s3_n_headline.clip(lower=1)
    info["nw_lfe_vs_s3_max_rel_diff"] = float(rel.max()) if len(rel) else None
    info["nw_lfe_vs_s3_pass_1pct"] = bool((rel <= 0.01).all())
    log(f"NW LFE counted-once vs S3 n_headline: max rel diff {rel.max():.4f} over {len(rel)} episodes")

    # --- GNSS (SVI, S1)
    G = pd.read_csv(f"{ROOT}/svi_science/results/s1/episode_geodetic.csv").set_index("episode")
    for c, nc in (("D", "gnss_D_mm"), ("sigma", "gnss_sigma_mm"), ("SN", "gnss_SN"), ("tierB_slip_mm", "gnss_tierB_slip_mm")):
        B[nc] = [G[c].get(e, np.nan) if g == "SVI" else np.nan for g, e in zip(B.region, B.episode)]

    # --- big-ETS definitions
    B["big"] = False; B["big_split_all"] = False; B["big_a1"] = False; B["big_a2"] = False; B["big_a3"] = B.s3_big
    info["big_split"] = {}
    for reg in ("SVI", "NW"):
        m = (B.region == reg) & (B.kind == "pnsn")
        x = np.log10(B.loc[m, "n_trem10"].values + 1)
        thr_all = split2(x)
        B.loc[m, "big_split_all"] = x >= thr_all
        mm = m & B.major.astype(bool)                       # Revision 1: split among MAJOR episodes only
        xm = np.log10(B.loc[mm, "n_trem10"].values + 1)
        thr = split2(xm)
        B.loc[mm, "big"] = xm >= thr
        info["big_split"][reg] = dict(log10_threshold=float(thr), n_trem10_threshold=float(10 ** thr - 1),
                                      n_big=int((xm >= thr).sum()), n_major_pnsn=int(mm.sum()), n_pnsn=int(m.sum()),
                                      all_episode_split_threshold=float(10 ** thr_all - 1), n_big_split_all=int((x >= thr_all).sum()))
        B.loc[m, "big_a1"] = B.loc[m, "family_exposure_days"] >= 0.5 * NFAM[reg] * 14
        c0, c1 = pd.Timestamp(COV[reg][0]), pd.Timestamp(COV[reg][1])
        K = int((c1 - c0).days // 365.25)
        top = B.loc[m, "n_trem10"].nlargest(K).index
        B.loc[top, "big_a2"] = True
        info["big_split"][reg]["a2_K"] = K
    w = B.kind == "wecc"
    B.loc[w, ["big", "big_split_all", "big_a1", "big_a2"]] = True
    B = B.sort_values(["region", "onset_PT"]).reset_index(drop=True)
    B.to_csv(f"{OUT}/episodes_b.csv", index=False)
    json.dump(info, open(f"{OUT}/build_info.json", "w"), indent=1, default=str)
    log("build done", json.dumps(info, default=str))


# ------------------------------------------------------------------------------------------------ tests
def big_list(B, reg, defn, with_wecc=True):
    E = B[B.region == reg].copy()
    col = {"primary": "big", "a1": "big_a1", "a2": "big_a2", "a3": "big_a3", "a4": "big"}[defn]
    E = E[E[col].astype(bool)].sort_values("onset_PT").reset_index(drop=True)
    if not with_wecc:
        E = E[E.kind != "wecc"].reset_index(drop=True)
    if defn == "a4":          # cycles: merge onsets < 120 d apart
        grp, g = [], 0
        for i in range(len(E)):
            if i > 0 and (E.onset_PT[i] - E.onset_PT[i - 1]).days >= 120:
                g += 1
            grp.append(g)
        E["cyc"] = grp
        agg = dict(episode=lambda s: "+".join(s), kind="first", onset_PT="first", end_PT="max", onset_year="first",
                   onset_doy="first", n_trem10=lambda s: s.sum(min_count=1), n_trem30=lambda s: s.sum(min_count=1),
                   family_exposure_days=lambda s: s.sum(min_count=1), lfe_svi_once075=lambda s: s.sum(min_count=1),
                   lfe_nw_once075=lambda s: s.sum(min_count=1), gnss_D_mm=lambda s: s.sum(min_count=1),
                   nw_era="first", t10_PT="first", start="first")
        E = E.groupby("cyc").agg(agg).reset_index(drop=True)
        E["duration_d"] = (E.end_PT - E.onset_PT).dt.total_seconds() / 86400
    return E


def perm_spearman(x, y, rng, n=NPERM):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return dict(n=int(len(x)), rho=np.nan, p=np.nan)
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    rx = (rx - rx.mean()) / rx.std(); ry = (ry - ry.mean()) / ry.std()
    r0 = float(np.mean(rx * ry))
    Y = np.array([rng.permutation(ry) for _ in range(n)])
    rs = Y @ rx / len(rx)
    return dict(n=int(len(x)), rho=float(r0), p=float((np.abs(rs) >= abs(r0) - 1e-12).mean()))


def rank_detrend(v, yr):
    v, yr = np.asarray(v, float), np.asarray(yr, float)
    out = np.full(len(v), np.nan); ok = np.isfinite(v) & np.isfinite(yr)
    if ok.sum() < 4:
        return out
    r = stats.rankdata(v[ok]); b = np.polyfit(yr[ok], r, 1)
    out[ok] = r - np.polyval(b, yr[ok])
    return out


def rayleigh(theta):
    n = len(theta); C = np.cos(theta).sum(); S = np.sin(theta).sum()
    R = np.hypot(C, S) / n; z = n * R ** 2
    p = np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n) - (24 * z - 132 * z ** 2 + 76 * z ** 3 - 9 * z ** 4) / (288 * n ** 2))
    mean_doy = (np.degrees(np.arctan2(S, C)) % 360) / 360 * 365.25
    return float(R), float(min(max(p, 0), 1)), float(mean_doy)


def harmonic(rate):
    ph = 2 * np.pi * np.arange(12) / 12
    c = np.sum(rate * np.exp(-1j * ph)) * 2 / 12
    return float(np.abs(c) / np.mean(rate)), float((-np.angle(c)) % (2 * np.pi) / (2 * np.pi) * 12 + 1)   # amp/mean, peak month (1..12.99)


def tests():
    rng = np.random.default_rng(0)
    B = pd.read_csv(f"{OUT}/episodes_b.csv", parse_dates=["onset_PT", "end_PT", "t10_PT", "w0", "w1"])
    J = {"Q1": {}, "Q2": {}, "Q3": {}, "Q4": {}, "Q5": {}}
    IV = []
    defs = ["primary", "a1", "a2", "a3", "a4"]
    lists = {}
    for reg in ("SVI", "NW"):
        for d in defs:
            for ww in ((True, False) if reg == "NW" and d != "a3" else (True,)):
                key = f"{reg}|{d}|{'wecc' if (reg == 'NW' and ww and d != 'a3') else 'pnsn'}"
                E = big_list(B, reg, d, ww)
                lists[key] = E
                on = E.onset_PT.values
                iv = np.diff(on).astype("timedelta64[s]").astype(float) / 86400
                mid = (E.onset_year.values[:-1] + E.onset_year.values[1:]) / 2
                for i in range(len(iv)):
                    IV.append(dict(key=key, region=reg, defn=d, from_ep=E.episode[i], to_ep=E.episode[i + 1],
                                   from_onset_PT=E.onset_PT[i], to_onset_PT=E.onset_PT[i + 1], interval_d=iv[i],
                                   interval_mo=iv[i] / 30.44, mid_year=mid[i]))
                res = dict(n_big=int(len(E)), n_int=int(len(iv)))
                if len(iv) >= 3:
                    cv = iv.std(ddof=1) / iv.mean()
                    span = (on[-1] - on[0]).astype("timedelta64[s]").astype(float) / 86400
                    cvc = []
                    for _ in range(NPERM):
                        u = np.sort(np.r_[0, rng.uniform(0, span, len(on) - 2), span]); di = np.diff(u)
                        cvc.append(di.std(ddof=1) / di.mean())
                    cvc = np.array(cvc)
                    before = iv[E.onset_PT.values[:-1] < np.datetime64(CUT_ERA)]; after = iv[E.onset_PT.values[:-1] >= np.datetime64(CUT_ERA)]
                    mw = stats.mannwhitneyu(before, after) if len(before) >= 2 and len(after) >= 2 else None
                    res.update(mean_d=float(iv.mean()), mean_mo=float(iv.mean() / 30.44), median_d=float(np.median(iv)),
                               sd_d=float(iv.std(ddof=1)), cv=float(cv), min_d=float(iv.min()), max_d=float(iv.max()),
                               cv_control_p01=float(np.percentile(cvc, 1)), cv_control_median=float(np.median(cvc)),
                               p_cv_below_control=float((cvc <= cv).mean()),
                               drift=perm_spearman(iv, mid, rng),
                               mean_before_2018_d=float(before.mean()) if len(before) else None, n_before=int(len(before)),
                               mean_after_2018_d=float(after.mean()) if len(after) else None, n_after=int(len(after)),
                               mw_p=float(mw.pvalue) if mw is not None else None)
                    # Q2
                    q2 = {}
                    sizes = ["n_trem10", "n_trem30", "duration_d", "family_exposure_days",
                             "lfe_svi_once075" if reg == "SVI" else "lfe_nw_once075"] + (["gnss_D_mm"] if reg == "SVI" else [])
                    prev = np.r_[np.nan, iv]; nxt = np.r_[iv, np.nan]
                    for s in sizes:
                        v = E[s].values.astype(float).copy()
                        if s == "lfe_nw_once075":
                            v[E.nw_era.values != "E3"] = np.nan        # LFE sizes compared only within E3
                        vd = rank_detrend(v, E.onset_year.values)
                        q2[s] = dict(trend=perm_spearman(v, E.onset_year.values, rng),
                                     slip_pred=perm_spearman(v, prev, rng), time_pred=perm_spearman(v, nxt, rng),
                                     slip_pred_detr=perm_spearman(vd, prev, rng), time_pred_detr=perm_spearman(vd, nxt, rng))
                    J["Q2"][key] = q2
                    # Q3a onsets
                    th = 2 * np.pi * (E.onset_PT.dt.dayofyear.values - 1) / 365.25
                    R, p, md = rayleigh(th)
                    Rc = []
                    for _ in range(NPERM):
                        o = rng.uniform(0, 365.25) + np.r_[0, np.cumsum(rng.permutation(iv))]
                        Rc.append(rayleigh(2 * np.pi * (o % 365.25) / 365.25)[0])
                    J["Q3"].setdefault("onsets", {})[key] = dict(n=int(len(th)), R=R, rayleigh_p=p, mean_doy=md,
                                                                  mean_month=float(md / 30.44 + 1),
                                                                  control_p=float((np.array(Rc) >= R).mean()),
                                                                  onset_months=E.onset_PT.dt.month.tolist())
                J["Q1"][key] = res
    IV = pd.DataFrame(IV); IV.to_csv(f"{OUT}/intervals.csv", index=False)

    # ---------------- Q3b between-episode activity
    P = pd.read_parquet(f"{OUT}/tremor_footprints.parquet", columns=["tpt", "reg10", "reg30", "d_svi", "d_nw"])
    P["date"] = P.tpt.dt.normalize()
    S7 = pd.read_csv(f"{ROOT}/svi_science/results/s7/between_events.csv")
    S7["date"] = pd.to_datetime(S7.time_PT).dt.normalize()
    MON = []
    for reg in ("SVI", "NW"):
        days = pd.date_range(COV[reg][0], COV[reg][1], freq="D")
        inep = np.zeros(len(days), bool)
        E = B[B.region == reg]
        for r in E.itertuples():
            a = pd.Timestamp(r.w0).tz_localize("UTC").tz_convert(TZ).tz_localize(None)
            b = pd.Timestamp(r.w1).tz_localize("UTC").tz_convert(TZ).tz_localize(None)
            inep |= (days >= a.normalize()) & (days <= b.normalize())
        bdays = days[~inep]
        cnt = P[P.reg10 == reg].groupby("date").size().reindex(bdays, fill_value=0)
        D = pd.DataFrame({"date": bdays, "n": cnt.values})
        D["any"] = (D.n > 0).astype(int)
        if reg == "SVI":
            lfe = S7.groupby("date").size().reindex(bdays, fill_value=0).values
            D["lfe"] = lfe
            D.loc[D.date < "2014-01-01", "lfe"] = np.nan
        D["y"] = D.date.dt.year; D["m"] = D.date.dt.month
        # big-ETS days (detection diagnostic): days inside [onset, end] of primary big
        Eb = E[E.big & (E.kind == "pnsn")]
        alld = P[P.reg10 == reg].groupby("date").size().reindex(days, fill_value=0)
        bigmask = np.zeros(len(days), bool)
        for r in Eb.itertuples():
            bigmask |= (days >= r.onset_PT.normalize()) & (days <= r.end_PT.normalize())
        DB = pd.DataFrame({"date": days[bigmask], "n": alld.values[bigmask]})
        DB["y"] = DB.date.dt.year; DB["m"] = DB.date.dt.month
        res = {}
        for var, DD in (("tremor_rate", D), ("tremor_anyday", D), ("lfe_svi_s7", D), ("bigETS_day_rate", DB)):
            if var == "lfe_svi_s7" and reg != "SVI":
                continue
            col = {"tremor_rate": "n", "tremor_anyday": "any", "lfe_svi_s7": "lfe", "bigETS_day_rate": "n"}[var]
            DD2 = DD.dropna(subset=[col])
            g = DD2.groupby(["y", "m"]).agg(s=(col, "sum"), nd=(col, "size")).reset_index()
            yrs = sorted(g.y.unique())
            Smat = np.zeros((len(yrs), 12)); Nmat = np.zeros((len(yrs), 12))
            for r in g.itertuples():
                Smat[yrs.index(r.y), r.m - 1] = r.s; Nmat[yrs.index(r.y), r.m - 1] = r.nd
            rate = Smat.sum(0) / np.maximum(Nmat.sum(0), 1)
            amp, pk = harmonic(rate)
            ctrl = []; prof = []
            for _ in range(NPERM):
                sh = rng.integers(0, 12, len(yrs))
                S2 = np.array([np.roll(Smat[i], sh[i]) for i in range(len(yrs))]); N2 = np.array([np.roll(Nmat[i], sh[i]) for i in range(len(yrs))])
                r2 = S2.sum(0) / np.maximum(N2.sum(0), 1)
                ctrl.append(harmonic(r2)[0]); prof.append(r2)
            prof = np.array(prof)
            res[var] = dict(amp_over_mean=amp, peak_month=pk, control_p=float((np.array(ctrl) >= amp).mean()),
                            control_amp_p95=float(np.percentile(ctrl, 95)), n_days=int(Nmat.sum()), total=float(Smat.sum()))
            for m in range(12):
                MON.append(dict(region=reg, var=var, month=m + 1, rate=rate[m], days=Nmat.sum(0)[m],
                                ctrl_p05=np.percentile(prof[:, m], 5), ctrl_p95=np.percentile(prof[:, m], 95)))
        J["Q3"].setdefault("between", {})[reg] = res
        J["Q3"]["between"][reg]["n_between_days"] = int(len(bdays)); J["Q3"]["between"][reg]["n_days"] = int(len(days))
    # catalogue-wide diagnostic: all epicentres outside both 30 km footprints, all days 2010-01-06..2026-05-27
    Pall = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    Pall["tpt"] = to_pt(Pall.time).values
    Pall = Pall[(Pall.lat < 46.5) | (Pall.lat > 49.8) | (Pall.lon < -126) | (Pall.lon > -121)]   # outside the B box entirely
    Pall["y"] = Pall.tpt.dt.year; Pall["m"] = Pall.tpt.dt.month
    days = pd.date_range(COV["SVI"][0], COV["SVI"][1], freq="D")
    nd = pd.Series(1, index=days).groupby([days.year, days.month]).size()
    s = Pall.groupby(["y", "m"]).size().reindex(nd.index, fill_value=0)
    rate = s.groupby(level=1).sum().values / nd.groupby(level=1).sum().values
    a, pk = harmonic(rate)
    J["Q3"]["catalogue_outside_box_diag"] = dict(amp_over_mean=a, peak_month=pk, note="all PNSN epicentres outside 46.5-49.8N/126-121W (mostly Oregon/N. California), all days; dominated by southern ETS, diagnostic only")
    for m in range(12):
        MON.append(dict(region="outside", var="catalogue_rate", month=m + 1, rate=rate[m], days=nd.groupby(level=1).sum().values[m]))
    pd.DataFrame(MON).to_csv(f"{OUT}/seasonal_monthly.csv", index=False)

    # ---------------- Q4/Q5 cross-border
    K = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_kinematics.csv").set_index(["region", "episode"])
    XO = pd.read_csv(f"{ROOT}/svi_science/results/s3/cross_region_overlaps.csv")
    PR = []
    c0, c1 = pd.Timestamp(COV["SVI"][0]), pd.Timestamp(COV["SVI"][1])
    for d in ("primary", "a1", "a2", "a4"):
        Es = lists[f"SVI|{d}|pnsn"]; En = lists[f"NW|{d}|{'wecc'}"]
        En = En[(En.kind != "wecc") & (En.onset_PT >= c0)].reset_index(drop=True)
        cand = []
        for i, a in Es.iterrows():
            for j, b in En.iterrows():
                dd = (a.onset_PT - b.onset_PT).total_seconds() / 86400
                ov = (a.onset_PT <= b.end_PT) and (b.onset_PT <= a.end_PT)
                if ov or abs(dd) <= 30:
                    cand.append((abs(dd), i, j, dd))
        cand.sort(); ui, uj = set(), set(); pairs = []
        for _, i, j, dd in cand:
            if i in ui or j in uj:
                continue
            ui.add(i); uj.add(j)
            a, b = Es.loc[i], En.loc[j]
            d10 = (a.t10_PT - b.t10_PT).total_seconds() / 86400
            row = dict(defn=d, svi_episode=a.episode, nw_episode=b.episode, svi_onset_PT=a.onset_PT, nw_onset_PT=b.onset_PT,
                       delay_d=dd, delay_t10_d=d10, year=a.onset_year, svi_n_trem10=a.n_trem10, nw_n_trem10=b.n_trem10)
            x = XO[(XO.svi_episode == a.episode) & (XO.nw_episode == b.episode)]
            row["s3_lfe_nw_lead_d"] = float(x.nw_lead_days.iloc[0]) if len(x) else np.nan
            row["s3_svi_front"] = K["front_class"].get(("SVI", a.episode), None) if d == "primary" else None
            row["s3_nw_front"] = K["front_class"].get(("NW", b.episode), None) if d == "primary" else None
            row["s3_svi_dir"] = K["direction"].get(("SVI", a.episode), None) if d == "primary" else None
            row["s3_nw_dir"] = K["direction"].get(("NW", b.episode), None) if d == "primary" else None
            pairs.append(row)
        pdf = pd.DataFrame(pairs).sort_values("svi_onset_PT")
        PR.append(pdf)
        dl = pdf.delay_d.values
        npos, nneg = int((dl > 0).sum()), int((dl < 0).sum())
        bt = [np.median(rng.choice(dl, len(dl))) for _ in range(NPERM)]
        frac = max(npos, nneg) / max(npos + nneg, 1)
        # coincidence control
        span = (c1 - c0).total_seconds() / 86400
        so = (Es.onset_PT - c0).dt.total_seconds().values / 86400; no = (En.onset_PT - c0).dt.total_seconds().values / 86400

        def npair(sv):
            return int(sum(np.min(np.abs(no - x)) <= 30 for x in sv))
        obs = npair(so)
        ctrl = np.array([npair((so + rng.uniform(0, span)) % span) for _ in range(NPERM)])
        q4 = dict(n_pairs=int(len(pdf)), n_svi_big=int(len(Es)), n_nw_big=int(len(En)),
                  unpaired_svi=sorted(set(Es.episode) - set(pdf.svi_episode)), unpaired_nw=sorted(set(En.episode) - set(pdf.nw_episode)),
                  median_delay_d=float(np.median(dl)), delay_ci90=[float(np.percentile(bt, 5)), float(np.percentile(bt, 95))],
                  mean_abs_delay_d=float(np.abs(dl).mean()), n_wa_first=npos, n_svi_first=nneg, frac_same_sign=float(frac),
                  binom_p=float(stats.binomtest(npos, npos + nneg).pvalue) if npos + nneg else None,
                  n_wa_first_t10=int((pdf.delay_t10_d > 0).sum()), n_svi_first_t10=int((pdf.delay_t10_d < 0).sum()),
                  coincident_within30=obs, coincident_control_mean=float(ctrl.mean()), coincident_p=float((ctrl >= obs).mean()))
        J["Q4"][d] = q4
        # Q5
        ivs = IV[(IV.key == f"SVI|{d}|pnsn")].interval_d.values
        ivn = IV[(IV.key == f"NW|{d}|wecc") & (IV.from_onset_PT >= c0)].interval_d.values
        J["Q5"][d] = dict(delay_vs_year=perm_spearman(pdf.delay_d.values, pdf.year.values, rng),
                          absdelay_vs_year=perm_spearman(np.abs(pdf.delay_d.values), pdf.year.values, rng),
                          n_big_common_span={"SVI": int(len(Es)), "NW": int(len(En))},
                          mean_interval_d={"SVI": float(ivs.mean()), "NW": float(ivn.mean())},
                          mw_intervals_p=float(stats.mannwhitneyu(ivs, ivn).pvalue))
    pd.concat(PR).to_csv(f"{OUT}/pairs.csv", index=False)
    json.dump(J, open(f"{OUT}/tests.json", "w"), indent=1, default=str)
    log("tests done")


# ------------------------------------------------------------------------------------------------ figures
def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#52514e",
                         "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    B = pd.read_csv(f"{OUT}/episodes_b.csv", parse_dates=["onset_PT", "end_PT", "t10_PT"])
    IV = pd.read_csv(f"{OUT}/intervals.csv", parse_dates=["from_onset_PT", "to_onset_PT"])
    PRS = pd.read_csv(f"{OUT}/pairs.csv", parse_dates=["svi_onset_PT", "nw_onset_PT"])
    MON = pd.read_csv(f"{OUT}/seasonal_monthly.csv")
    J = json.load(open(f"{OUT}/tests.json"))
    grey = "#b5b3ad"

    # ---- timeline
    fig, ax = plt.subplots(4, 1, figsize=(11, 9), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 2, 2]))
    for k, reg in enumerate(("SVI", "NW")):
        a = ax[k]; E = B[B.region == reg]
        for r in E.itertuples():
            if r.kind == "wecc":
                a.axvspan(r.onset_PT, r.end_PT, color=COL[reg], alpha=0.25, hatch="//", lw=0)
                a.text(r.onset_PT, 1.02, "WECC", transform=a.get_xaxis_transform(), fontsize=7, color="#52514e")
                continue
            w = max(r.duration_d, 2)
            a.bar(r.onset_PT, r.n_trem10, width=w, align="edge", color=COL[reg] if r.big else grey, edgecolor="none")
        a.set_ylabel(f"{reg}\nPNSN tremor in footprint")
        big = E[E.big]
        for r in big.itertuples():
            if r.kind != "wecc":
                a.text(r.onset_PT, r.n_trem10 * 1.02, r.episode[2:], fontsize=6.5, rotation=90, va="bottom", ha="left", color="#52514e")
        a.set_ylim(0, E.n_trem10.max() * 1.35)
    a = ax[2]
    for reg, c in (("SVI", "lfe_svi_once075"), ("NW", "lfe_nw_once075")):
        E = B[(B.region == reg) & B[c].notna()]
        if reg == "NW":
            E = E[E.nw_era == "E3"]
        a.scatter(E.onset_PT, E[c], s=np.where(E.big, 28, 10), color=COL[reg], label=f"{reg} LFEs (counted once, P≥0.75)",
                  edgecolor="white", lw=0.5, zorder=3)
    a.set_yscale("log"); a.set_ylabel("LFE count\nper episode"); a.legend(fontsize=7, frameon=False, loc="upper left")
    a = ax[3]
    E = B[(B.region == "SVI") & B.gnss_D_mm.notna()]
    a.errorbar(E.onset_PT, E.gnss_D_mm, yerr=E.gnss_sigma_mm, fmt="o", ms=4, color=COL["SVI"], ecolor=grey, lw=1)
    a.axhline(0, color=grey, lw=0.8); a.set_ylabel("SVI GNSS offset\n(mm, Tier A)")
    a.set_xlim(pd.Timestamp("2006-10-01"), pd.Timestamp("2026-10-01"))
    ax[0].set_title("Episodes 2007–2026: coloured = big ETS (primary rule), grey = smaller listed episodes; bar width = duration", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/b_timeline.png", dpi=150); plt.close(fig)

    # ---- interval vs size
    fig, ax = plt.subplots(2, 3, figsize=(12, 7))
    for k, reg in enumerate(("SVI", "NW")):
        key = f"{reg}|primary|{'wecc' if reg == 'NW' else 'pnsn'}"
        iv = IV[IV.key == key]
        a = ax[k, 0]
        a.plot(iv.to_onset_PT, iv.interval_mo, "o-", color=COL[reg], ms=5, lw=1.2)
        for d, mk in (("a1", "s"), ("a2", "^")):
            v = IV[IV.key == f"{reg}|{d}|{'wecc' if reg == 'NW' else 'pnsn'}"]
            a.plot(v.to_onset_PT, v.interval_mo, mk, mfc="none", color=grey, ms=5, label=f"alt {d}")
        a.axhline(14, color="#52514e", ls=":", lw=1); a.text(a.get_xlim()[0], 14.3, " 14 mo (literature)", fontsize=7, color="#52514e")
        a.set_ylabel(f"{reg}: interval to this ETS (months)"); a.legend(fontsize=7, frameon=False)
        a.xaxis.set_major_locator(matplotlib.dates.YearLocator(4))
        E = big_list(B, reg, "primary", True)
        on = E.onset_PT.values
        ivd = np.diff(on).astype("timedelta64[s]").astype(float) / 86400 / 30.44
        q2 = J["Q2"][key]["n_trem10"]
        for j, (x, lab, tst) in enumerate(((np.r_[np.nan, ivd], "interval BEFORE (months)", "slip_pred"),
                                          (np.r_[ivd, np.nan], "interval AFTER (months)", "time_pred"))):
            a = ax[k, j + 1]
            a.scatter(x, E.n_trem10, color=COL[reg], s=30, edgecolor="white", lw=0.5, zorder=3)
            if reg == "SVI":
                for xx, yy, g in zip(x, E.n_trem10, E.gnss_D_mm):
                    if np.isfinite(g) and np.isfinite(xx):
                        a.annotate(f"{g:.1f} mm", (xx, yy), fontsize=6, color="#52514e", xytext=(3, 3), textcoords="offset points")
            r = q2[tst]
            a.set_title(f"Spearman ρ = {r['rho']:.2f}, p = {r['p']:.2f} (n = {r['n']})", fontsize=8, loc="left")
            a.set_xlabel(lab); a.set_ylabel(f"{reg} tremor count of the ETS")
    fig.suptitle("Wait vs size (primary big-ETS list; SVI labels = Tier-A GNSS offset)", fontsize=10, x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/b_interval_size.png", dpi=150); plt.close(fig)

    # ---- seasonal
    fig = plt.figure(figsize=(12, 7))
    for k, reg in enumerate(("SVI", "NW")):
        a = fig.add_subplot(2, 3, 3 * k + 1, projection="polar")
        a.set_theta_zero_location("N"); a.set_theta_direction(-1)
        E = big_list(B, reg, "primary", True)
        th = 2 * np.pi * (E.onset_PT.dt.dayofyear.values - 1) / 365.25
        rr = 1 + 0.6 * np.arange(len(th)) / max(len(th) - 1, 1)     # radius = order in time (inner = oldest)
        a.scatter(th, rr, color=COL[reg], s=30, edgecolor="white", lw=0.5, zorder=3)
        o = J["Q3"]["onsets"][f"{reg}|primary|{'wecc' if reg == 'NW' else 'pnsn'}"]
        a.annotate("", xy=(2 * np.pi * o["mean_doy"] / 365.25, o["R"] * 1.6), xytext=(0, 0),
                   arrowprops=dict(arrowstyle="->", color="#0b0b0b", lw=1.5))
        a.set_xticks(2 * np.pi * np.arange(12) / 12); a.set_xticklabels(list("JFMAMJJASOND")); a.set_yticks([]); a.set_ylim(0, 1.7)
        a.set_title(f"{reg} big-ETS onsets (n={o['n']})\nRayleigh p={o['rayleigh_p']:.2f}, control p={o['control_p']:.2f}", fontsize=8)
        for j, (var, lab) in enumerate((("tremor_rate", "between-episode tremor / day"), ("bigETS_day_rate", "tremor / day on big-ETS days (detection check)"))):
            a = fig.add_subplot(2, 3, 3 * k + 2 + j)
            M = MON[(MON.region == reg) & (MON["var"] == var)]
            a.fill_between(M.month, M.ctrl_p05, M.ctrl_p95, color=grey, alpha=0.4, lw=0, label="scrambled control 5–95 %")
            a.bar(M.month, M.rate, color=COL[reg], width=0.7)
            r = J["Q3"]["between"][reg][var]
            a.set_title(f"{reg}: {lab}\nannual amp/mean {r['amp_over_mean']:.2f}, peak month {r['peak_month']:.1f}, control p={r['control_p']:.3f}", fontsize=8, loc="left")
            a.set_xticks(range(1, 13)); a.set_xticklabels(list("JFMAMJJASOND"))
            if j == 0 and k == 0:
                a.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/b_seasonal.png", dpi=150); plt.close(fig)

    # ---- cross-border
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    P0 = PRS[PRS.defn == "primary"]
    a = ax[0]
    a.axhline(0, color=grey, lw=0.8)
    a.scatter(P0.svi_onset_PT, P0.delay_d, color="#0b0b0b", s=30, zorder=3, label="5th-pct tremor onset")
    a.scatter(P0.svi_onset_PT, P0.delay_t10_d, color=grey, s=18, marker="s", zorder=2, label="10th-pct onset")
    for r in P0.itertuples():
        a.annotate(r.svi_episode[2:], (r.svi_onset_PT, r.delay_d), fontsize=6.5, xytext=(3, 3), textcoords="offset points", color="#52514e")
    a.set_ylabel("SVI onset − WA-NW onset (days)\n> 0: Washington first")
    a.legend(fontsize=7, frameon=False); a.set_title("Paired big ETS: who starts first", fontsize=9, loc="left")
    a = ax[1]
    for reg in ("SVI", "NW"):
        E = big_list(B, reg, "primary", True)
        E = E[E.onset_PT >= pd.Timestamp(COV["SVI"][0])]
        a.step(E.onset_PT, np.arange(1, len(E) + 1), where="post", color=COL[reg], lw=2, label=f"{reg} primary")
        E2 = big_list(B, reg, "a2", True); E2 = E2[E2.onset_PT >= pd.Timestamp(COV["SVI"][0])]
        a.step(E2.onset_PT, np.arange(1, len(E2) + 1), where="post", color=COL[reg], lw=1.2, ls="--", label=f"{reg} a2 (16 largest)")
    a.set_ylabel("cumulative number of big ETS since 2010-01-06"); a.legend(frameon=False, fontsize=8)
    a.set_title("Recurrence clocks (common span 2010-01 → 2026-05)", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/b_crossborder.png", dpi=150); plt.close(fig)
    log("figures done")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    {"build": build, "tests": tests, "figures": figures}[sys.argv[1]]()
