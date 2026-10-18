#!/usr/bin/env python
"""WA pre-production W2: station health of the 21 core stations 2007-01-01 .. 2026-09-15 from data on disk, in the
SVI_REF R2 / KMR-check style (lfe_extension/scripts/svi_ref_r2.py; lfe_south/results/SH_CHECKS_RESULTS.md).

  epochs   channel epochs (sensor, sample rate, sensitivity, azimuth) from the LOCAL StationXML files + the SVI gain
           table (B926) -> w2/epochs.csv, w2/gains.csv; epoch-change list -> w2/epoch_changes.csv
  raw      per channel-day raw metrics (svi_ref_r2.raw_day, copied unchanged): gap fraction, DC level and drift, std,
           flat fraction, longest constant run, rms 2-8 Hz and 0.1-0.5 Hz (counts and nm/s), 2-9 Hz line prominence.
           The SVI R2 daily table (lfe_extension/results/refine/r2/r2_health_daily.parquet) already holds the 13
           shared stations 2014-01-01 .. 2026-09-15 from the same files and the same code -> reused, not recomputed.
           Raw roots: production-visible roots first, then any other root (never traces_dl*, never traces_sealed).
  days     the stack-test days: per year the strongest tremor days of the WA regions (PNSN tables, 2010+) or the array
           truth-scan days (2007-2009), never in the sealed months 2007-01 / 2009-05, never 2011-04
  cache    WA-chain 20 Hz cache (wa_cache.cache_one, frozen) for those days where the WA cache lacks them, written to
           results/wa_pre/w2/cache20hz (the WA cache itself is never written)
  (mf)     the frozen wa_mf.py on those days, network = 21 core, templates build + step4b_a (node-specific), data only
  lo       per data candidate >= 8.75 MAD: per-site cc at the candidate and at 1000 random origins (step1_features
           method) -> LEAVE-ONE-SITE-OUT network z for every site (the network without that site, median / MAD from
           the random origins, same unscaled-MAD convention as net_mad)
  stack    per station-month (selected days only): "confident repeats" for station s = candidates with leave-s's-site-out
           z >= 8.75 on >= 2 other sites AND >= 2 OTHER sites with their own site z >= 3 (so one loud site cannot make
           the event; if s's site is not in the family's network that day: net_mad >= 8.75 and >= 2 sites at z >= 3),
           declustered over families at +-2 s (best lo-z kept).  Station 3-C segment at the family's placement
           (-2..+14 s around the template start), rms-normalised, stacked per family; template = the family's template
           at the station's CURRENT node, else the latest earlier node's (carry-over, flagged).  Scrambled control =
           the same station read on a different selected day (>= 3 days away, same month if possible else <= 60 days),
           same time of day, no confident repeat of that family within +-1 h there; 10 draws (SVI R2 rule).
  summary  monthly table, effective station set per month, exclusions, swap checks at every epoch / node change,
           bad-for-dv/v windows, positive control (LRIV after 2019-05-30).
Guard: a data-reading stage fails (w2/W2_<stage>_FAILED, exit 5) if > 3 % of its station-days error or are missing.

  python scripts/pre_w2.py epochs | raw --workers 20 | days | cache --workers 16 | lo --workers 12 | stack --workers 16 | summary
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse, warnings, itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

O2 = f"{P.PRE}/w2"
CACHE2 = f"{O2}/cache20hz"
SETUP2 = f"{O2}/setup"
MF2 = f"{O2}/res/mf_w2"
FS = W.FS
PRE_S, LEN_S = 2.0, 16.0
NCTL = 10
NMIN_EV = 5
LAGMAX = int(2.0 * FS)
TCONF = W.T_CORE
SVI_DAILY = f"{P.ROOT}/lfe_extension/results/refine/r2/r2_health_daily.parquet"
SVI_MONTHLY = f"{P.ROOT}/lfe_extension/results/refine/r2/r2_health_monthly.csv"
SVI_GAINS = f"{P.ROOT}/lfe_extension/results/refine/r2/gains.csv"
SVI_T0, SVI_T1 = "2014-01-01", "2026-09-15"
XMLS = [f"{P.ROOT}/lfe_array_discovery/results/stationxml/{n}.xml" for n in ("CN", "PB", "UW")] + \
       [f"{P.ROOT}/lfe_picker/results/p4/stationxml/{n}.xml" for n in ("PB", "UW", "CN")]
# node / sensor changes named in WA_PROD_PLAN W2 (+ the node tables)
PLAN_CHANGES = [("CN.PGC", "2015-01-12"), ("CN.PGC", "2017-08-10"), ("CN.VGZ", "2017-03-03"), ("UW.DOSE", "2019-05-07"),
                ("UW.GNW", "2019-05-07"), ("UW.LRIV", "2019-05-30"), ("UW.HDW", "2020-12-04"), ("UW.HDW", "2022-09-01")]
GAP_CHECKS = [("PB.B004", "2020"), ("PB.B014", "2023"), ("PB.B014", "2024"), ("PB.B926", "2025")]


def ensure(*d):
    for x in d:
        os.makedirs(x, exist_ok=True)


def guard_fail(stage, n_bad, n_tot, msg=""):
    fr = n_bad / max(n_tot, 1)
    P.log(f"[w2:{stage}] guard: {n_bad} of {n_tot} station-days error/missing ({fr:.2%}) {msg}")
    if fr > 0.03:
        P.marker(f"w2/W2_{stage.upper()}_FAILED", f"{n_bad} of {n_tot} ({fr:.2%}) {msg}")
        sys.exit(5)


# ------------------------------------------------------------------------------------------------ epochs / gains
def epochs(a):
    ensure(O2)
    from obspy import read_inventory
    core = set(P.core_list())
    rows = []
    for p in XMLS:
        if not os.path.exists(p):
            continue
        inv = read_inventory(p)
        for net in inv:
            for st in net:
                ns = f"{net.code}.{st.code}"
                if ns not in core:
                    continue
                for ch in st:
                    if ch.code[:2] not in ("HH", "BH", "EH"):
                        continue
                    sens = ch.response.instrument_sensitivity if ch.response is not None else None
                    rows.append(dict(netsta=ns, sta=st.code, loc=ch.location_code, cha=ch.code, start=str(ch.start_date)[:10],
                                     end=str(ch.end_date)[:10] if ch.end_date else "2599-12-31", fs=float(ch.sample_rate),
                                     sensor=(ch.sensor.description if ch.sensor else "") or "",
                                     gain=float(sens.value) if sens and sens.value else np.nan,
                                     az=ch.azimuth, dip=ch.dip, src=os.path.relpath(p, P.ROOT)))
    E = pd.DataFrame(rows).drop_duplicates(["netsta", "loc", "cha", "start"], keep="first")
    G0 = pd.read_csv(SVI_GAINS, dtype={"loc": str})
    b9 = G0[G0.sta == "B926"].copy()
    for r in b9.itertuples():
        E = pd.concat([E, pd.DataFrame([dict(netsta="PB.B926", sta="B926", loc="", cha=r.cha, start=r.start, end=r.end, fs=100.0,
                                             sensor="", gain=r.gain, az=np.nan, dip=np.nan, src=r.src)])])
    E = E.sort_values(["netsta", "cha", "start"])
    E.to_csv(f"{O2}/epochs.csv", index=False)
    G = E[["sta", "cha", "start", "end", "gain", "src"]].copy()
    G.to_csv(f"{O2}/gains.csv", index=False)
    # change list: every epoch start inside 2007-01-01 .. 2026-09-15 with what changed
    ch = []
    for ns, g in E.groupby("netsta"):
        g = g[g.cha.str[:2].isin(["HH", "BH", "EH"])]
        starts = sorted({s for s in g.start if P.T0 < s <= P.T1})
        for s in starts:
            bef = g[(g.start < s) & (g.end >= s)]
            aft = g[g.start == s]
            what = []
            if set(bef.cha) != set(aft.cha):
                what.append(f"channels {','.join(sorted(set(bef.cha)))} -> {','.join(sorted(set(aft.cha)))}")
            if set(bef.sensor) != set(aft.sensor):
                what.append(f"sensor {';'.join(sorted(set(bef.sensor)))[:60]} -> {';'.join(sorted(set(aft.sensor)))[:60]}")
            gb = bef.groupby("cha").gain.first(); ga = aft.groupby("cha").gain.first()
            com = sorted(set(gb.index) & set(ga.index))
            if com and any(abs(ga[c] / gb[c] - 1) > 0.02 for c in com if gb[c] > 0):
                what.append("gain " + ";".join(f"{c} x{ga[c]/gb[c]:.3f}" for c in com if gb[c] > 0))
            ch.append(dict(netsta=ns, date=s, change="; ".join(what) or "metadata epoch only (no channel/sensor/gain change)",
                           source="StationXML"))
    for ns, d in PLAN_CHANGES:
        ch.append(dict(netsta=ns, date=d, change="node change (WA_PROD_PLAN W2 / node tables)", source="plan"))
    C = pd.DataFrame(ch).sort_values(["netsta", "date"])
    C.to_csv(f"{O2}/epoch_changes.csv", index=False)
    print(C.to_string())


def gain_of(G, sta, cha, date):
    g = G[(G.sta == sta) & (G.cha == cha) & (G.start <= date) & (G.end > date)]
    return float(g.gain.iloc[0]) if len(g) else np.nan


# ------------------------------------------------------------------------------------------------ raw metrics
def raw_day(path, sta, date, G):
    """svi_ref_r2.raw_day, copied unchanged."""
    from obspy import read
    from scipy.signal import butter, sosfiltfilt, welch, resample_poly
    from scipy.ndimage import median_filter
    st = read(path)
    for tr in st:
        fr = round(float(tr.stats.sampling_rate))
        if fr > 0 and float(tr.stats.sampling_rate) != fr and abs(float(tr.stats.sampling_rate) - fr) / fr < 1e-4:
            tr.stats.sampling_rate = float(fr)
    st.merge(method=1, fill_value=None)
    rows = []
    for tr in st:
        fs = float(tr.stats.sampling_rate)
        if fs < 1: continue
        x = tr.data
        msk = np.ma.getmaskarray(x) if np.ma.isMaskedArray(x) else np.zeros(len(x), bool)
        xv = np.asarray(np.ma.filled(x, 0) if np.ma.isMaskedArray(x) else x, np.float64)
        nexp = int(round(86400 * fs))
        gap = 1.0 - (len(xv) - msk.sum()) / nexp
        good = ~msk
        r = dict(date=date, sta=sta, seed=tr.id, cha=tr.stats.channel, loc=tr.stats.location, fs=fs, npts=len(xv), gap_frac=float(gap))
        if good.sum() < 0.1 * nexp:
            rows.append(r); continue
        xg = xv[good]
        r["dc"] = float(np.mean(xg)); r["std"] = float(np.std(xg))
        nh = int(3600 * fs)
        r["drift"] = float(np.mean(xg[-nh:]) - np.mean(xg[:nh])) if len(xg) > 2 * nh else np.nan
        eq = np.zeros(len(xv), bool); eq[1:] = (xv[1:] == xv[:-1]) & good[1:] & good[:-1]
        r["flat_frac"] = float(eq[good].mean())
        if eq.any():
            dd = np.diff(np.concatenate([[0], eq.astype(np.int8), [0]]))
            a_ = np.flatnonzero(dd == 1); b_ = np.flatnonzero(dd == -1)
            r["longest_const_s"] = float((b_ - a_).max() + 1) / fs
        else:
            r["longest_const_s"] = 0.0
        y = np.where(good, xv - r["dc"], 0.0)
        hb = sosfiltfilt(butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos"), y)
        r["rms_2_8"] = float(np.sqrt(np.mean(hb[good] ** 2)))
        q = int(round(fs / 10.0))
        y10 = resample_poly(y, 1, q) if q > 1 else y
        lb = sosfiltfilt(butter(4, [0.1 / 5.0, 0.5 / 5.0], btype="band", output="sos"), y10)
        r["rms_01_05"] = float(np.sqrt(np.mean(lb[int(600 * 10):-int(600 * 10)] ** 2))) if len(lb) > 24000 else np.nan
        gn = gain_of(G, sta, tr.stats.channel, date)
        r["gain"] = gn
        r["nms_2_8"] = r["rms_2_8"] / gn * 1e9 if np.isfinite(gn) else np.nan
        r["nms_01_05"] = r["rms_01_05"] / gn * 1e9 if np.isfinite(gn) else np.nan
        try:
            f, p = welch(y, fs=fs, nperseg=int(100 * fs))
            bg = median_filter(p, size=51, mode="nearest")
            pr = p / np.maximum(bg, 1e-300); band = (f >= 2.0) & (f <= 9.0)
            r["line_prom_max"] = float(pr[band].max()); r["n_lines4"] = int((pr[band] >= 4).sum())
            r["line_f_max"] = float(f[band][np.argmax(pr[band])])
        except Exception:
            pass
        rows.append(r)
    return rows


_G = {}


def _init_raw(G):
    warnings.filterwarnings("ignore"); _G["G"] = G


def raw_job(args):
    ns, ym, dates = args
    sta = ns.split(".")[1]
    out, stat = [], []
    t0 = time.time()
    for ds in dates:
        p = P.prod_raw_path(ns, ds) or (P.raw_files(ns, ds) or [None])[0]
        if p is None:
            stat.append((ds, ns, "missing", "")); continue
        try:
            out += raw_day(p, sta, ds, _G["G"]); stat.append((ds, ns, "ok", os.path.relpath(p, P.ROOT)))
        except Exception as e:
            stat.append((ds, ns, "error", f"{type(e).__name__}: {str(e)[:60]}"))
    return ns, ym, out, stat, time.time() - t0


def raw(a):
    ensure(f"{O2}/daily", f"{O2}/jobs")
    G = pd.read_csv(f"{O2}/gains.csv")
    core = P.core_list(); N = P.nodes()
    SVI = set(pd.read_parquet(SVI_DAILY, columns=["sta"]).sta.unique())
    done = {os.path.basename(p)[:-4] for p in glob.glob(f"{O2}/jobs/*.ok")}
    jobs = []
    for ns in core:
        exp = sorted(P.expected_days(ns, N))
        sta = ns.split(".")[1]
        if sta in SVI:
            exp = [d for d in exp if not (SVI_T0 <= d <= SVI_T1)]
        by = {}
        for d in exp:
            by.setdefault(d[:7], []).append(d)
        for ym, ds in sorted(by.items()):
            if f"{ns}_{ym}" not in done:
                jobs.append((ns, ym, ds))
    if a.stations:
        jobs = [j for j in jobs if j[0] in a.stations.split(",")]
    if a.months:
        jobs = [j for j in jobs if j[1] in a.months.split(",")]
    P.log(f"[w2:raw] {len(jobs)} station-month jobs ({sum(len(j[2]) for j in jobs)} station-days), {a.workers} workers; "
          f"SVI-reused stations 2014-01..2026-09-15: {sorted(SVI & {c.split('.')[1] for c in core})}")
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init_raw, initargs=(G,)) as ex:
        fut = [ex.submit(raw_job, j) for j in jobs]
        for n, fu in enumerate(as_completed(fut), 1):
            ns, ym, out, stat, secs = fu.result()
            pd.DataFrame(out).to_parquet(f"{O2}/daily/{ns}_{ym}.parquet", index=False)
            pd.DataFrame(stat, columns=["date", "netsta", "status", "msg"]).to_csv(f"{O2}/jobs/{ns}_{ym}_raw.csv", index=False)
            open(f"{O2}/jobs/{ns}_{ym}.ok", "w").write(f"{secs:.1f}\n")
            if n % 50 == 0 or n == len(fut):
                P.log(f"[w2:raw] {n}/{len(fut)} {time.time()-t0:.0f}s")
    S = pd.concat([pd.read_csv(f, dtype=str) for f in glob.glob(f"{O2}/jobs/*_raw.csv")], ignore_index=True)
    S.to_csv(f"{O2}/raw_status.csv", index=False)
    # the guard counts unreadable among present files; missing (expected, no file) are W1's gap list (tier A)
    pres = S[S.status != "missing"]
    guard_fail("raw", int((pres.status == "error").sum()), len(pres), f"(missing, no file: {int((S.status == 'missing').sum())})")


# ------------------------------------------------------------------------------------------------ stack-test days
def days(a):
    ensure(O2)
    td = P.tremor_days()
    ep = pd.read_csv(f"{W.EXT}/tremor_episodes.csv")
    rows = []
    for y in range(2010, 2027):
        for rg, k in (("NW_WA", 10), ("SW_WA", 6), ("OLY2010", 4)):
            s = pd.Series(td[rg]); s = s[(s.index.str[:4] == str(y)) & (s.index <= P.T1)]
            s = s[~s.index.str[:7].isin(P.SEALED_MONTHS | P.EXCLUDED_MONTHS)]
            for d, n in s.sort_values(ascending=False).head(k).items():
                if n > 0:
                    rows.append(dict(date=d, why=f"{rg} top-{k} tremor day of {y}", n_windows=int(n)))
    # 2007-2009: no PNSN catalogue -> days the array-truth runs scanned in the fit months (never 2007-01, 2009-05)
    scan = pd.read_csv(f"{P.STUDY}/results/merlin_tests/t1_setup/scan_days.csv")
    sd = sorted(set(scan.date))
    for m in ["2007-10", "2008-05", "2009-02", "2009-04", "2009-08"]:
        dd = [d for d in sd if d[:7] == m]
        if not dd:
            e = ep[(ep.region == "NW_WA") & (ep.start.str[:7] == m)]
            for r in e.itertuples():
                dd += [x.strftime("%Y-%m-%d") for x in pd.date_range(r.start, r.end)]
        step = max(1, len(dd) // 10)
        for d in dd[::step][:10]:
            rows.append(dict(date=d, why=f"array-truth scan day {m}", n_windows=np.nan))
    D = pd.DataFrame(rows).groupby("date").agg(why=("why", lambda v: "; ".join(v)), n_windows=("n_windows", "max")).reset_index()
    D = D[~D.date.str[:7].isin(P.SEALED_MONTHS | P.EXCLUDED_MONTHS)].sort_values("date")
    out = f"{O2}/stack_days.csv"
    if not os.path.exists(out):
        D.to_csv(out, index=False)
    D = pd.read_csv(out)
    print(len(D), D.date.str[:4].value_counts().sort_index().to_dict())
    # setup for the MF: core stations (step 4b setup rows, which include B926), step 4b families
    ensure(SETUP2)
    S = pd.read_csv(f"{P.STUDY}/results/step4b/setup/stations.csv")
    core = P.core_list()
    S[S.netsta.isin(core)].to_csv(f"{SETUP2}/stations.csv", index=False)
    import shutil
    shutil.copy(f"{P.STUDY}/results/step4b/setup/families.csv", f"{SETUP2}/families.csv")
    ensure(f"{O2}/res")
    if not os.path.exists(f"{O2}/res/templates"):
        os.symlink(f"{P.STUDY}/results/templates", f"{O2}/res/templates")


# ------------------------------------------------------------------------------------------------ cache
def _cache_job(j):
    import wa_cache
    try:
        return wa_cache.cache_one(j)
    except Exception as e:
        return dict(netsta=j[0], date=j[1], status=f"crash:{type(e).__name__}:{str(e)[:60]}")


def cache(a):
    ensure(CACHE2)
    os.environ["WA_EXTRA_NODES"] = f"{P.STUDY}/results/step4b/setup/extra_nodes.csv"
    W.TRACE_ROOTS[:] = P.PROD_ROOTS
    W.CACHE = CACHE2                         # wa_common.cache_path reads the module global -> writes go to w2/cache20hz
    D = pd.read_csv(f"{O2}/stack_days.csv")
    if a.months:
        D = D[D.date.str[:7].isin(a.months.split(","))]
    core = P.core_list()
    jobs, rows0, skipped = [], [], 0
    for d in D.date:
        for ns in core:
            ts = pd.Timestamp(d); tail = f"{ns}/{ts.year}/{int(ts.dayofyear):03d}.npz"
            if os.path.exists(f"{P.WA_CACHE}/{tail}"):
                skipped += 1; rows0.append(dict(netsta=ns, date=d, status="wa_cache")); continue
            p = P.prod_raw_path(ns, d)
            if p is None:
                anyp = P.raw_files(ns, d)
                rows0.append(dict(netsta=ns, date=d, status="raw_not_prod_visible" if anyp else "no_raw")); continue
            jobs.append((ns, d, p))
    P.log(f"[w2:cache] {len(jobs)} station-days to cache ({skipped} already in the WA cache), {a.workers} workers")
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(_cache_job, jobs, chunksize=4)):
            rows.append(r)
            if (i + 1) % 500 == 0:
                P.log(f"[w2:cache] {i+1}/{len(jobs)}")
    M = pd.DataFrame(rows0 + rows)
    M.to_csv(f"{O2}/cache_manifest{'_smoke' if a.months else ''}.csv", index=False)
    print(M.status.value_counts().to_string())
    att = M[~M.status.isin(["wa_cache", "no_raw", "raw_not_prod_visible"])]
    bad = att[~att.status.isin(["ok", "cached"])]
    guard_fail("cache", len(bad), len(att), str(bad.status.value_counts().to_dict()))
    # day list for the MF window injection
    D.assign(kind="episode")[["date", "kind"]].to_csv(f"{O2}/mf_days.csv", index=False)


# ------------------------------------------------------------------------------------------------ leave-one-site-out z
NRAND = 1000
_L = {}


def _init_lo():
    warnings.filterwarnings("ignore")
    os.environ["WA_SETUP"] = SETUP2; W.SETUP = SETUP2
    P.patch_wa([CACHE2])
    import wa_mf
    T, Q = wa_mf.load_tpl(["build", "step4b_a"], all_templates=True)
    S = pd.read_csv(f"{SETUP2}/stations.csv")
    _L.update(T=T, SITE=dict(zip(S.netsta, S.site)))


def lo_day(ds):
    import step1_features as S1
    T, SITE = _L["T"], _L["SITE"]
    n = W.NT; NV = W.NDAY - n + 1
    D0 = pd.Timestamp(ds).timestamp()
    J = json.load(open(f"{MF2}/run_{ds}.json")); nodes = J["nodes"]
    C = pd.read_parquet(f"{MF2}/det_{ds}.parquet")
    C = C[(C.series == "data") & (C.net_mad >= TCONF)].copy()
    out = f"{O2}/lo/lo_{ds}.parquet"
    if not len(C):
        pd.DataFrame(columns=["fid", "t", "net_mad", "site", "lo_z", "n_other_sites", "in_net"]).to_parquet(out); return ds, 0
    C["i"] = np.round((C.t.to_numpy() - D0) * FS).astype(np.int64)
    FD = pd.read_parquet(f"{MF2}/famday_{ds}.parquet"); FD = FD[FD.observable]
    fam_sta = {int(f): s.split(",") for f, s in zip(FD.fid, FD.stations)}
    rng = np.random.default_rng(1000 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    R = rng.integers(0, NV, NRAND)
    Cg = {f: g for f, g in C.groupby("fid")}
    CC, RC, NCH = {}, {}, {}
    sta_fams = {}
    for f, ss in fam_sta.items():
        for ns in ss:
            sta_fams.setdefault(ns, []).append(f)
    for ns, fl in sta_fams.items():
        p = W.cache_path(ns, ds)
        if not os.path.exists(p):
            continue
        z = np.load(p); node = str(z["node"])
        if node != nodes.get(ns):
            continue
        x = z["x"]; m = z["m"]; live = z["live"]
        bad, ew = S1.guard_np(x, live, m)
        efl = 1e-10 * np.median(ew[live])
        for f in fl:
            key = (f, ns, node)
            if key not in T or f not in Cg:
                continue
            tp = T[key]; chi = np.nonzero(tp["live"] & live)[0]
            if not len(chi):
                continue
            d = int(round(tp["off"] * FS))
            j = Cg[f].i.to_numpy() + d
            allj = np.concatenate([j, R + d]); okj = (allj >= 0) & (allj < NV); jj = np.clip(allj, 0, NV - 1)
            idx = jj[:, None] + np.arange(n)[None]
            num = np.einsum("cnk,ck->cn", x[chi[:, None, None], idx[None]], tp["tpl"][chi])
            E = np.sqrt(np.maximum(ew[chi[:, None], jj[None]], efl)).astype(np.float32)
            cc = np.where(okj & ~bad[jj], (num / E).mean(0), np.nan)
            CC.setdefault(f, {})[ns] = cc[:len(j)]; RC.setdefault(f, {})[ns] = cc[len(j):]; NCH.setdefault(f, {})[ns] = len(chi)
        del x, ew, bad, z
    rows = []
    for f, g in Cg.items():
        if f not in CC:
            continue
        sts = [ns for ns in fam_sta.get(f, []) if ns in CC[f]]
        wm = {ns: float(NCH[f][ns]) * T[(f, ns, nodes[ns])]["w"] for ns in sts}
        groups = {}
        for ns in sts:
            groups.setdefault(SITE[ns], []).append(ns)
        sites = sorted(groups)
        if len(sites) < 2:
            continue
        Ws = np.array([np.sqrt(sum(wm[s] ** 2 for s in groups[q])) for q in sites])
        SC = np.zeros((len(g), len(sites))); SR = np.zeros((NRAND, len(sites)))
        for k, q in enumerate(sites):
            mem = groups[q]; ws = np.array([wm[s] for s in mem])
            SC[:, k] = np.nan_to_num((np.stack([CC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum())
            SR[:, k] = (np.stack([RC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum()
        okr = np.isfinite(SR).all(1)
        SR = SR[okr]
        smed = np.median(SR, 0); smad = 1.4826 * np.median(np.abs(SR - smed), 0)
        Zs = (SC - smed) / np.where(smad > 0, smad, np.nan)
        Z3 = np.nan_to_num(Zs, nan=-99) >= 3
        nz3_all = Z3.sum(1)
        for k, q in enumerate(sites):
            keep = np.arange(len(sites)) != k
            if keep.sum() < 2:
                lz = np.full(len(g), np.nan)
            else:
                wk = Ws[keep] / Ws[keep].sum()
                net = SC[:, keep] @ wk; nr = SR[:, keep] @ wk
                med = np.median(nr); mad = max(np.median(np.abs(nr - med)), 1e-6)
                lz = (net - med) / mad
            rows.append(pd.DataFrame(dict(fid=f, t=g.t.to_numpy(), net_mad=g.net_mad.to_numpy(), site=q, lo_z=lz,
                                          n_other_sites=int(keep.sum()), in_net=True, site_z=Zs[:, k],
                                          n_z3_other=nz3_all - Z3[:, k], n_z3_all=nz3_all)))
    F = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    F["date"] = ds
    F.to_parquet(out + ".tmp", index=False); os.replace(out + ".tmp", out)
    return ds, len(F)


def lo(a):
    ensure(f"{O2}/lo")
    ds = sorted(os.path.basename(p)[4:14] for p in glob.glob(f"{MF2}/det_*.parquet"))
    ds = [d for d in ds if not os.path.exists(f"{O2}/lo/lo_{d}.parquet")]
    if a.months:
        ds = [d for d in ds if d[:7] in a.months.split(",")]
    P.log(f"[w2:lo] {len(ds)} days, {a.workers} workers")
    with ProcessPoolExecutor(a.workers, initializer=_init_lo) as ex:
        for fu in as_completed([ex.submit(lo_day, d) for d in ds]):
            d, nr = fu.result()
            P.log(f"[w2:lo] {d}: {nr} rows")


# ------------------------------------------------------------------------------------------------ stacks
def seg_of(x, t0, t, o):
    j = int(round((t + o - PRE_S - t0) * FS)); m = int(round(LEN_S * FS))
    if j < 0 or j + m > x.shape[1]:
        return None
    s = x[:, j:j + m].astype(np.float64)
    r = np.sqrt(np.mean(s * s))
    if not np.isfinite(r) or r <= 0:
        return None
    return s / r


def _init_stack():
    warnings.filterwarnings("ignore")
    os.environ["WA_SETUP"] = SETUP2; W.SETUP = SETUP2
    P.patch_wa([CACHE2])
    import wa_mf
    T, Q = wa_mf.load_tpl(["build", "step4b_a"], all_templates=True)
    S = pd.read_csv(f"{SETUP2}/stations.csv")
    nd = {}
    for (f, ns, node) in T:
        nd.setdefault(ns, set()).add(node)
    _L.update(T=T, SITE=dict(zip(S.netsta, S.site)), NODES_T={k: sorted(v, key=lambda s: s.split("@")[1]) for k, v in nd.items()})


def ref_template(T, NODES_T, f, ns, node):
    """template of (f, ns) at the station's current node, else the latest earlier node's (carry-over)."""
    if (f, ns, node) in T:
        return T[(f, ns, node)], False, node
    d = node.split("@")[1] if "@" in node else "9999"
    cand = [n for n in NODES_T.get(ns, []) if n.split("@")[1] <= d and (f, ns, n) in T]
    if not cand:
        cand = [n for n in NODES_T.get(ns, []) if (f, ns, n) in T]
    if not cand:
        return None, False, ""
    n = cand[-1]
    return T[(f, ns, n)], True, n


def confident_for_station(LOd, MFdet, site, fams_ok):
    """confident repeats for a station on one day: leave-its-site-out z >= TCONF on >= 2 other sites if its site is in
    the family's network, else net_mad >= TCONF (the station is not in that family's network that day); only families
    with a template at the station (any node); declustered over families at +-2 s (best z kept)."""
    import bisect
    a = LOd[(LOd.site == site) & (LOd.lo_z >= TCONF) & (LOd.n_other_sites >= 2) & (LOd.n_z3_other >= 2)][["fid", "t", "lo_z"]]
    fam_in = set(LOd[LOd.site == site].fid)
    nz = LOd.drop_duplicates(["fid", "t"])[["fid", "t", "n_z3_all"]]
    b = MFdet[~MFdet.fid.isin(fam_in)][["fid", "t", "net_mad"]].merge(nz, on=["fid", "t"], how="inner")
    b = b[b.n_z3_all >= 2][["fid", "t", "net_mad"]].rename(columns={"net_mad": "lo_z"})
    c = pd.concat([a, b], ignore_index=True)
    c = c[c.fid.isin(fams_ok)].sort_values("lo_z", ascending=False)
    keep, tk = [], []
    for k, t in enumerate(c.t.to_numpy()):
        j = bisect.bisect_left(tk, t)
        if (j < len(tk) and tk[j] - t <= W.MATCH_S) or (j > 0 and t - tk[j - 1] <= W.MATCH_S):
            continue
        tk.insert(j, t); keep.append(k)
    return c.iloc[keep].reset_index(drop=True)


def stack(a):
    ensure(f"{O2}/stacks", f"{O2}/jobs_stack")
    D = pd.read_csv(f"{O2}/stack_days.csv")
    if a.months:
        D = D[D.date.str[:7].isin(a.months.split(","))]
    S = pd.read_csv(f"{SETUP2}/stations.csv"); SITE = dict(zip(S.netsta, S.site))
    core = P.core_list()
    LO = {os.path.basename(p)[3:13]: pd.read_parquet(p) for p in glob.glob(f"{O2}/lo/lo_*.parquet")}
    MD = {}
    for ds in LO:
        x = pd.read_parquet(f"{MF2}/det_{ds}.parquet")
        MD[ds] = x[(x.series == "data") & (x.net_mad >= TCONF)]
    import wa_mf
    T, _ = wa_mf.load_tpl(["build", "step4b_a"], all_templates=True)
    fam_tpl = {}
    for (f, ns_, node) in T:
        fam_tpl.setdefault(ns_, set()).add(f)
    jobs = []
    by = {}
    for d in D.date:
        by.setdefault(d[:7], []).append(d)
    memo = {}
    for ns in core:
        site = SITE.get(ns, ns.split(".")[1])
        for ym, dd in sorted(by.items()):
            # control days: same month if >= 2 other days >= 3 d away, else the selected days within +-60 d
            pool = sorted(set(d2 for d2 in D.date if abs((pd.Timestamp(d2) - pd.Timestamp(dd[0])).days) <= 60))
            conf = {}
            for d2 in pool:
                if d2 in LO:
                    if (ns, d2) not in memo:
                        memo[(ns, d2)] = confident_for_station(LO[d2], MD[d2], site, fam_tpl.get(ns, set()))
                    conf[d2] = memo[(ns, d2)]
            conf_m = {d2: (v if d2 in dd else v.iloc[:0]) for d2, v in conf.items()}
            jobs.append((ns, ym, pool, conf_m, conf))
    P.log(f"[w2:stack] {len(jobs)} station-month jobs, {a.workers} workers")
    rows = []
    with ProcessPoolExecutor(a.workers, initializer=_init_stack) as ex:
        fut = {ex.submit(stack_job2, j): j for j in jobs}
        for fu in as_completed(fut):
            ns, ym, st, secs = fu.result()
            rows += st
    St = pd.DataFrame(rows, columns=["date", "netsta", "status"]).drop_duplicates()
    St.to_csv(f"{O2}/stack_cache_status.csv", index=False)
    # guard: expected = core station-days of the selected days inside metadata epochs with a raw file on disk
    if not os.path.exists(f"{P.PRE}/w1/w1_station_days.parquet"):
        P.log("[w2:stack] W1 table not there yet: guard skipped (smoke)"); return
    W1 = pd.read_parquet(f"{P.PRE}/w1/w1_station_days.parquet", columns=["netsta", "date", "expected", "readable", "prod_visible"])
    W1 = W1[W1.date.isin(set(D.date)) & W1.expected & W1.readable & W1.prod_visible]
    St1 = St[St.date.isin(set(D.date))]
    k = W1.merge(St1, on=["netsta", "date"], how="left")
    bad = k[k.status.fillna("nocache") != "ok"]
    guard_fail("stack", len(bad), len(k), "(selected-day station-days with readable raw but no usable cache)")


def stack_job2(args):
    ns, ym, pool, conf_m, conf_all = args
    # repeats from conf_m (the month's days); control exclusion from conf_all
    T, SITE, NT_ = _L["T"], _L["SITE"], _L["NODES_T"]
    t00 = time.time()
    X, nodes, st = {}, {}, []
    for ds in pool:
        p = W.cache_path(ns, ds)
        if not os.path.exists(p):
            st.append((ds, ns, "nocache")); continue
        try:
            z = np.load(p); X[ds] = (z["x"].astype(np.float32), float(z["t0"])); nodes[ds] = str(z["node"]); st.append((ds, ns, "ok"))
        except Exception as e:
            st.append((ds, ns, f"error:{type(e).__name__}"))
    m = int(round(LEN_S * FS))
    SD, ND, SC, NC, TP, CARRY = {}, {}, {}, {}, {}, {}
    rng = np.random.default_rng(int(ym.replace("-", "")) * 7 + sum(map(ord, ns)))
    dl = sorted(X)
    for ds in dl:
        c = conf_m.get(ds)
        if c is None or not len(c):
            continue
        same_m = [d2 for d2 in dl if d2[:7] == ym and abs((pd.Timestamp(d2) - pd.Timestamp(ds)).days) >= 3]
        cand0 = same_m if len(same_m) >= 2 else [d2 for d2 in dl if abs((pd.Timestamp(d2) - pd.Timestamp(ds)).days) >= 3]
        for r in c.itertuples():
            f = int(r.fid)
            tp, carry, tnode = ref_template(T, NT_, f, ns, nodes[ds])
            if tp is None:
                continue
            s = seg_of(X[ds][0], X[ds][1], r.t, tp["off"])
            if s is None:
                continue
            if f not in SD:
                SD[f] = np.zeros((3, m)); ND[f] = 0; SC[f] = np.zeros((NCTL, 3, m)); NC[f] = np.zeros(NCTL, int)
                TP[f] = tp["tpl"].astype(np.float64) * tp["live"][:, None]; CARRY[f] = (carry, tnode)
            SD[f] += s; ND[f] += 1
            tod = r.t - pd.Timestamp(ds).timestamp()
            for kk in range(NCTL):
                if not cand0:
                    break
                for _try in range(4):
                    d2 = cand0[rng.integers(len(cand0))]
                    t2 = pd.Timestamp(d2).timestamp() + tod
                    cb = conf_all.get(d2)
                    if cb is not None and len(cb):
                        tt = cb.t[cb.fid == f].to_numpy()
                        if len(tt) and np.min(np.abs(tt - t2)) < 3600:
                            continue
                    s2 = seg_of(X[d2][0], X[d2][1], t2, tp["off"])
                    if s2 is None:
                        continue
                    SC[f][kk] += s2; NC[f][kk] += 1; break
    fams = sorted(SD)
    np.savez_compressed(f"{O2}/stacks/{ns}_{ym}.npz", fams=np.array(fams, int),
                        SD=np.stack([SD[f] for f in fams]) if fams else np.zeros((0, 3, m)),
                        ND=np.array([ND[f] for f in fams], int), SC=np.stack([SC[f] for f in fams]) if fams else np.zeros((0, NCTL, 3, m)),
                        NC=np.stack([NC[f] for f in fams]) if fams else np.zeros((0, NCTL), int),
                        TP=np.stack([TP[f] for f in fams]) if fams else np.zeros((0, 3, W.NT)),
                        carry=np.array([CARRY[f][0] for f in fams], bool), tnode=np.array([CARRY[f][1] for f in fams]),
                        nodes=np.array(sorted({nodes[d] for d in dl if d[:7] == ym})))
    return ns, ym, st, time.time() - t00


# ------------------------------------------------------------------------------------------------ summary (svi_ref_r2 functions, copied)
def lagcc(S, Tp):
    from numpy.lib.stride_tricks import sliding_window_view
    L = Tp.shape[1]; i0 = int(PRE_S * FS)
    Wv = sliding_window_view(S, L, axis=1)[:, i0 - LAGMAX:i0 + LAGMAX + 1]
    Wc = Wv - Wv.mean(-1, keepdims=True); Tc = Tp - Tp.mean(1, keepdims=True)
    num = (Wc * Tc[:, None, :]).sum(-1)
    wn = np.sqrt((Wc ** 2).sum(-1)); tn = np.sqrt((Tc ** 2).sum(-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        per = num / (wn * tn[:, None])
        per = np.where(tn[:, None] > 0, per, np.nan)            # channels the template does not have
        mean3 = np.nanmean(np.where(np.isfinite(per), per, np.nan), 0)
        summed = num.sum(0) / np.sqrt((wn ** 2).sum(0) * (tn ** 2).sum())
    return mean3, summed


LAGS = np.arange(-LAGMAX, LAGMAX + 1)
PLACE = np.abs(LAGS) <= int(round(0.25 * FS))


def fam_stats(S, Tp):
    m3, sm = lagcc(S, Tp)
    if not np.isfinite(m3).any():
        return np.nan, np.nan, np.nan, np.nan
    j = int(np.nanargmax(m3)); ja = int(np.nanargmax(np.abs(sm)))
    return float(np.nanmax(m3[PLACE])), float(m3[j]), LAGS[j] / FS, float(sm[ja])


TRANSFORMS = [(p, s) for p in itertools.permutations(range(3)) for s in itertools.product((1, -1), repeat=3)]


def transform(S, p, s):
    return np.stack([s[i] * S[p[i]] for i in range(3)])


def rotate(S, th):
    c, s = np.cos(th), np.sin(th)
    return np.stack([S[0], c * S[1] + s * S[2], -s * S[1] + c * S[2]])


def set_stats(stacks, tpls):
    v = np.array([fam_stats(S, Tp) for S, Tp in zip(stacks, tpls)], float).reshape(-1, 4)
    if not len(v):
        return dict(cc_place=np.nan, cc_best=np.nan, lag_med=np.nan, lagcons=np.nan, cc_a2=np.nan)
    return dict(cc_place=float(np.nanmedian(v[:, 0])), cc_best=float(np.nanmedian(v[:, 1])), lag_med=float(np.nanmedian(v[:, 2])),
                lagcons=float(np.nanmean(np.abs(v[:, 2]) <= 0.25 + 1e-9)), cc_a2=float(np.nanmedian(v[:, 3])))


def station_month(path):
    z = np.load(path)
    fams, SD, ND, SC, NC, TP, carry = z["fams"], z["SD"], z["ND"], z["SC"], z["NC"], z["TP"], z["carry"]
    keep = [i for i in range(len(fams)) if ND[i] >= NMIN_EV]
    out = dict(n_fam=len(keep), n_ev=int(ND[keep].sum()) if keep else 0, n_fam_any=len(fams),
               carry_frac=float(carry[keep].mean()) if keep else np.nan, nodes_month=";".join(z["nodes"]),
               tnodes=";".join(sorted(set(z["tnode"][keep]))) if keep else "")
    if len(keep) < 3 or out["n_ev"] < 30:
        out["tested"] = False; return out
    out["tested"] = True
    tp = [TP[i] for i in keep]
    st = [SD[i] / ND[i] for i in keep]
    out.update(set_stats(st, tp))
    C_ = []
    for k in range(SC.shape[1]):
        ii = [j for j, i in enumerate(keep) if NC[i, k] >= NMIN_EV]
        if len(ii) >= 3:
            C_.append(set_stats([SC[keep[j], k] / NC[keep[j], k] for j in ii], [tp[j] for j in ii]))
    for c in ("cc_place", "cc_best", "lagcons", "cc_a2"):
        v = [x[c] for x in C_ if np.isfinite(x[c])]
        out[f"ctl_{c}_med"] = float(np.median(v)) if v else np.nan
        out[f"ctl_{c}_p95"] = float(np.quantile(v, 0.95)) if v else np.nan
    out["n_ctl_draws"] = len(C_)
    best = (-np.inf, None)
    for pp, sg in TRANSFORMS:
        v = float(np.nanmedian([fam_stats(transform(S, pp, sg), Tp)[0] for S, Tp in zip(st, tp)]))
        if v > best[0]:
            best = (v, (pp, sg))
    out["perm_best_cc"] = best[0]; out["perm_best"] = f"{best[1][0]}{best[1][1]}"
    out["perm_is_identity"] = best[1] == ((0, 1, 2), (1, 1, 1))
    rb = (-np.inf, 0.0)
    for th in np.radians(np.arange(0, 360, 10)):
        v = float(np.nanmedian([fam_stats(rotate(S, th), Tp)[0] for S, Tp in zip(st, tp)]))
        if v > rb[0]:
            rb = (v, th)
    out["rot_best_cc"] = rb[0]; out["rot_best_deg"] = float(np.degrees(rb[1]))
    out["at_control"] = bool((out["cc_place"] <= out.get("ctl_cc_place_p95", np.nan)) or (out["lagcons"] <= out.get("ctl_lagcons_p95", np.nan)))
    return out


def _sm(path):
    warnings.filterwarnings("ignore")
    b = os.path.basename(path)[:-4]
    ns, ym = b.rsplit("_", 1)
    r = station_month(path); r.update(netsta=ns, month=ym); return r


def daily_flags(Dy):
    """svi_ref_r2.daily_flags, copied unchanged."""
    Dy = Dy.sort_values(["seed", "date"]).copy()
    for c in ("rms_2_8", "rms_01_05", "dc"):
        if c not in Dy: Dy[c] = np.nan
        Dy[f"{c}_rm"] = Dy.groupby("seed")[c].transform(lambda v: v.rolling(31, center=True, min_periods=5).median())
    Dy["dc_rmad"] = Dy.groupby("seed").dc.transform(lambda v: (v - v.rolling(31, center=True, min_periods=5).median()).abs()
                                                      .rolling(31, center=True, min_periods=5).median())
    for c in ("longest_const_s", "flat_frac", "gap_frac", "std", "line_prom_max"):
        if c not in Dy: Dy[c] = np.nan
    Dy["f_dead"] = (Dy.rms_2_8 <= 1e-3 * Dy.rms_2_8_rm) | (Dy["std"] <= 0) | (Dy.longest_const_s >= 60) | Dy.rms_2_8.isna()
    Dy["f_flat"] = Dy.flat_frac >= 0.01
    Dy["f_noisy"] = Dy.rms_2_8 >= 10 * Dy.rms_2_8_rm
    Dy["f_lf_noisy"] = Dy.rms_01_05 >= 10 * Dy.rms_01_05_rm
    Dy["f_dc_jump"] = (Dy.dc - Dy.dc_rm).abs() > 20 * Dy.dc_rmad.clip(lower=1e-9)
    Dy["f_gap"] = Dy.gap_frac > 0.5
    Dy["f_line"] = Dy.line_prom_max >= 20
    return Dy


def summary(a):
    core = P.core_list(); NS = {c.split(".")[1]: c for c in core}
    fs = sorted(glob.glob(f"{O2}/daily/*.parquet"))
    own = pd.concat([pd.read_parquet(f) for f in fs if os.path.getsize(f) > 0], ignore_index=True)
    own["src"] = "wa_pre"
    svi = pd.read_parquet(SVI_DAILY)
    svi = svi[svi.sta.isin(NS) & (svi.date >= SVI_T0) & (svi.date <= SVI_T1)]
    svi = svi[[c for c in svi.columns if not (c.endswith("_rm") or c.endswith("_rmad") or c.startswith("f_"))]].assign(src="svi_r2")
    Dy = pd.concat([own, svi], ignore_index=True)
    Dy = daily_flags(Dy)
    Dy["netsta"] = Dy.sta.map(NS)
    Dy.to_parquet(f"{O2}/r2_health_daily.parquet", index=False)
    fl = [c for c in Dy.columns if c.startswith("f_")]
    # band-matched channels only for the monthly levels (the cache band of that day)
    SDy = Dy.groupby(["netsta", "date"]).agg(**{c: (c, "max") for c in fl}, seeds=("seed", lambda v: ",".join(sorted(set(v)))),
                                              fsset=("fs", lambda v: ",".join(str(x) for x in sorted(set(v))))).reset_index()
    SDy = SDy.sort_values(["netsta", "date"])
    SDy["f_meta_change"] = (SDy.groupby("netsta").seeds.shift() != SDy.seeds) | (SDy.groupby("netsta").fsset.shift() != SDy.fsset)
    SDy.loc[SDy.groupby("netsta").head(1).index, "f_meta_change"] = False
    SDy["month"] = SDy.date.str[:7]
    SDy.to_parquet(f"{O2}/r2_health_stationday.parquet", index=False)
    # noise level per station-month: median over days of the median over the day's 2-8 Hz channels (nm/s), Z separately
    Dy["month"] = Dy.date.str[:7]
    Dy["is_z"] = Dy.cha.str[-1].eq("Z")
    lev = Dy.groupby(["netsta", "month"]).agg(nms_2_8_med=("nms_2_8", "median"), rms_2_8_med=("rms_2_8", "median"),
                                             nms_01_05_med=("nms_01_05", "median"), line_prom_med=("line_prom_max", "median"),
                                             gain_med=("gain", "median")).reset_index()
    levz = Dy[Dy.is_z].groupby(["netsta", "month"]).agg(nms_2_8_z=("nms_2_8", "median"), rms_2_8_z=("rms_2_8", "median")).reset_index()
    lev = lev.merge(levz, on=["netsta", "month"], how="left").sort_values(["netsta", "month"])
    # gain / offset jump: month-to-month step of the Z 2-8 Hz level in counts > x2 with no change in nm/s gain metadata
    lev["step_counts"] = lev.groupby("netsta").rms_2_8_z.transform(lambda v: v / v.shift())
    lev["step_nms"] = lev.groupby("netsta").nms_2_8_z.transform(lambda v: v / v.shift())
    lev["f_level_step"] = (lev.step_counts > 2) | (lev.step_counts < 0.5)
    lev["f_level_step_nms"] = (lev.step_nms > 2) | (lev.step_nms < 0.5)
    # stack test
    sp = sorted(glob.glob(f"{O2}/stacks/*.npz"))
    import multiprocessing as mp
    with mp.Pool(a.workers) as pool:
        SMv = pool.map(_sm, sp, chunksize=2)
    SM = pd.DataFrame(SMv)
    SM.to_csv(f"{O2}/r2_stack_tests.csv", index=False)
    # month grid
    N = P.nodes()
    rows = []
    for ns in core:
        e = P.expected_days(ns, N)
        months = sorted({d[:7] for d in e})
        for ym in months:
            rows.append(dict(netsta=ns, month=ym, n_alive_days=sum(1 for d in e if d[:7] == ym)))
    M = pd.DataFrame(rows)
    g = SDy.groupby(["netsta", "month"])
    agg = g.agg(n_raw_days=("date", "size"), **{f"n_{c[2:]}": (c, "sum") for c in fl + ["f_meta_change"]}).reset_index()
    M = M.merge(agg, on=["netsta", "month"], how="left").merge(lev, on=["netsta", "month"], how="left")
    M = M.merge(SM, on=["netsta", "month"], how="left")
    M["n_raw_days"] = M.n_raw_days.fillna(0).astype(int)
    for c in [c for c in M.columns if c.startswith("n_") and c not in ("n_alive_days",)]:
        M[c] = M[c].fillna(0)
    M["tested"] = M.tested.fillna(False).astype(bool)
    M["status"] = "ok"
    M.loc[M.n_raw_days == 0, "status"] = "nodata"
    M.loc[(M.n_raw_days > 0) & (M.n_dead > 0.5 * M.n_raw_days), "status"] = "dead"
    ctl = M.tested & M.at_control.fillna(False).astype(bool)
    M.loc[ctl & (M.status == "ok"), "status"] = "control_level"
    M.loc[~M.tested & (M.status == "ok"), "status"] = "untested"
    # instrument epoch = node (band/fs/date) from the node tables; untested months inherit the nearest tested month
    def node_at(ns, ym):
        d = ym + "-15"
        g = N[(N.ns == ns) & (N.start <= d) & (N.end > d)]
        return g.sort_values("start").node.iloc[-1] if len(g) else ""
    M["inst_epoch"] = [node_at(s, m) for s, m in zip(M.netsta, M.month)]
    M["status_eff"] = M.status
    for (ns, ie), gg in M.groupby(["netsta", "inst_epoch"]):
        t = gg[gg.tested]
        if not len(t):
            continue
        tm = pd.PeriodIndex(t.month, freq="M")
        for i in gg.index[gg.status == "untested"]:
            d = np.abs((tm - pd.Period(M.at[i, "month"], freq="M")).map(lambda x: x.n))
            j = t.index[int(np.argmin(d))]
            M.at[i, "status_eff"] = "inherit_" + M.at[j, "status"]
    M["absent"] = M.status_eff.isin(["nodata", "dead", "control_level", "inherit_control_level", "inherit_dead"])
    # SVI R2 status for the shared stations (SVI families) as a reference column
    sv = pd.read_csv(SVI_MONTHLY)
    sv["netsta"] = sv.sta.map(NS)
    M = M.merge(sv[["netsta", "month", "status_eff"]].rename(columns={"status_eff": "svi_r2_status"}), on=["netsta", "month"], how="left")
    M.to_csv(f"{O2}/r2_health_monthly.csv", index=False)
    E = M.groupby("month").apply(lambda gg: pd.Series(dict(
        n_effective=int((~gg.absent).sum()), effective=",".join(gg[~gg.absent].netsta),
        absent=",".join(f"{s}:{st}" for s, st in zip(gg[gg.absent].netsta, gg[gg.absent].status_eff))))).reset_index()
    E.to_csv(f"{O2}/r2_effective_sets.csv", index=False)
    X = M[M.absent & ~M.status_eff.eq("nodata")].copy()
    X["reason"] = X.status_eff + " cc_place=" + X.cc_place.round(3).astype(str) + " ctl_p95=" + X.ctl_cc_place_p95.round(3).astype(str)
    X[["netsta", "month", "reason"]].to_csv(f"{O2}/r2_exclusions.csv", index=False)
    # swap checks at every change
    C = pd.read_csv(f"{O2}/epoch_changes.csv")
    sw = []
    for r in C.drop_duplicates(["netsta", "date"]).itertuples():
        gg = M[(M.netsta == r.netsta) & M.tested]
        b = gg[gg.month < r.date[:7]].tail(2); af = gg[gg.month > r.date[:7]].head(2)
        sw.append(dict(netsta=r.netsta, date=r.date, change=r.change, months_before=",".join(b.month), months_after=",".join(af.month),
                       cc_before=b.cc_place.median() if len(b) else np.nan, ctl95_before=b.ctl_cc_place_p95.median() if len(b) else np.nan,
                       cc_after=af.cc_place.median() if len(af) else np.nan, ctl95_after=af.ctl_cc_place_p95.median() if len(af) else np.nan,
                       lag_after=af.lag_med.median() if len(af) else np.nan, carry_after=af.carry_frac.median() if len(af) else np.nan,
                       perm_after=";".join(af.perm_best.astype(str)) if len(af) else "",
                       perm_identity_after=bool(af.perm_is_identity.fillna(False).all()) if len(af) else None,
                       rot_after=";".join(af.rot_best_deg.round(0).astype(str)) if len(af) else "",
                       status_after=";".join(af.status) if len(af) else "untested"))
    for ns, yr in GAP_CHECKS:
        gg = M[(M.netsta == ns) & (M.month.str[:4] == yr)]
        sw.append(dict(netsta=ns, date=yr, change="gap year check", status_after=";".join(sorted(set(gg.status_eff))),
                       cc_after=gg.cc_place.median(), ctl95_after=gg.ctl_cc_place_p95.median(),
                       months_after=",".join(gg[gg.tested].month)))
    SW = pd.DataFrame(sw)
    SW["verdict"] = np.where(SW.cc_after.notna() & SW.ctl95_after.notna(),
                             np.where(SW.cc_after > SW.ctl95_after, "above control (carry-over OK)" , "AT CONTROL LEVEL"), "untested")
    SW.to_csv(f"{O2}/r2_swap_checks.csv", index=False)
    # positive control
    L = M[(M.netsta == "UW.LRIV") & (M.month >= "2019-06") & M.tested]
    L0 = M[(M.netsta == "UW.LRIV") & (M.month < "2019-05") & M.tested]
    pc = dict(n_tested_post=len(L), n_control_level_post=int((L.status == "control_level").sum()),
              cc_place_post_median=float(L.cc_place.median()) if len(L) else np.nan,
              ctl_p95_post_median=float(L.ctl_cc_place_p95.median()) if len(L) else np.nan,
              n_tested_pre=len(L0), n_control_level_pre=int((L0.status == "control_level").sum()),
              cc_place_pre_median=float(L0.cc_place.median()) if len(L0) else np.nan,
              ctl_p95_pre_median=float(L0.ctl_cc_place_p95.median()) if len(L0) else np.nan,
              svi_reference="SVI R2: LRIV 2021-10 0.17 vs control p95 0.09 (weak, above control); A-2 -0.216 on its own event set")
    pc["flags_lriv_post2019"] = bool(len(L) > 0 and pc["n_control_level_post"] >= 0.5 * len(L))
    json.dump(pc, open(f"{O2}/r2_positive_control.json", "w"), indent=1, default=float)
    # bad for dv/v
    bad = []
    for ns, gg in SDy.groupby("netsta"):
        for c in ["f_dead", "f_flat", "f_noisy", "f_gap", "f_dc_jump", "f_meta_change", "f_line"]:
            v = gg[c].to_numpy(bool); dd = gg.date.to_numpy()
            if not v.any():
                continue
            e = np.diff(np.concatenate([[0], v.astype(np.int8), [0]]))
            for i, j in zip(np.flatnonzero(e == 1), np.flatnonzero(e == -1)):
                bad.append(dict(netsta=ns, start=dd[i], end=dd[j - 1], n_days=int(j - i), reason=c[2:]))
    for r in M[M.absent & ~M.status_eff.eq("nodata")].itertuples():
        bad.append(dict(netsta=r.netsta, start=r.month + "-01", end=str(pd.Period(r.month).end_time.date()), n_days=int(r.n_alive_days),
                        reason=f"month {r.status_eff} (stack cc {r.cc_place})"))
    for r in lev[lev.f_level_step].itertuples():
        bad.append(dict(netsta=r.netsta, start=r.month + "-01", end=r.month + "-01", n_days=0,
                        reason=f"2-8 Hz level step x{r.step_counts:.2f} counts / x{r.step_nms:.2f} nm/s vs previous month"))
    B = pd.DataFrame(bad)
    B.to_csv(f"{O2}/r2_bad_for_dvv.csv", index=False)
    P.log(f"[w2:summary] {len(M)} station-months; status {M.status_eff.value_counts().to_dict()}; exclusions {len(X)}; pc {pc}")
    print(M.status_eff.value_counts().to_string()); print(SW.to_string()); print(json.dumps(pc, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["epochs", "raw", "days", "cache", "lo", "stack", "summary"])
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--stations", default=""); ap.add_argument("--months", default="")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(epochs=epochs, raw=raw, days=days, cache=cache, lo=lo, stack=stack, summary=summary)[a.what](a)
