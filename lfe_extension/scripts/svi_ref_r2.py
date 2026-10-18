#!/usr/bin/env python
"""SVI_REF R2 -- station health 2014-2026 from the RAW counts (KMR method, lfe_south SH_KMRFIX / SH_CHECKS) plus a
per-station monthly STACK test against the frozen 2012 template with a scrambled control.

One job per (station, month), CPU workers:
  raw (native rate, every channel of the file)   per channel-day: fs, cha, loc, gap fraction, DC level, drift
        (last-hour minus first-hour mean), std, flat fraction (sample == previous), longest constant run (s),
        rms 2-8 Hz and 0.1-0.5 Hz in counts and nm/s (sensitivity from LOCAL StationXML / the A-2 audit table; no
        query), 2-9 Hz line prominence (Welch / running median: max and number >= 4).
  stack (20 Hz cache, the detector's own chain)  for every v1.1 catalogue repeat of the month inside a PNSN tremor
        window (cat_core & tremor_flag: mostly real LFEs, so the test has power), the station's 3-C segment at its
        template placement (-2 s .. +14 s around the template start), rms-normalised and stacked per family;
        SCRAMBLED CONTROL = the same station read on a DIFFERENT day of the month (>= 3 days away, same time of day,
        no T_core detection of that family within +-1 h there), 10 independent draws.  Every station alive is
        tested, including stations the detector excludes (LRIV after 2019-05-30: the positive control).
Aggregation (`summary`): per (station, month) signed 3-C cc of each family stack against its 2012 template, best of
+-2 s; station-month value = median over families with >= 5 repeats; control p95 = 95th percentile over the 10
control draws of the same median; best lag; best of the 48 channel permutation / sign transforms; best horizontal
rotation (10-deg steps).  Status: nodata / dead (> 50 % dead days) / control_level (tested and data median <= control
p95: COUNTS AS ABSENT) / ok / untested (inherits the nearest tested month of the same instrument epoch, else present).
Outputs OUT/r2/: r2_health_daily.parquet, r2_health_monthly.csv, r2_effective_sets.csv, r2_exclusions.csv
(the station-months the production matched filter leaves out), r2_swap_checks.csv, r2_bad_for_dvv.csv,
r2_positive_control.json.  Guard: > 3 % of expected raw or cache station-days erroring / missing -> exit 5.

  python scripts/svi_ref_r2.py run --workers 20 [--months 2021-10] [--stations LRIV,B926]
  python scripts/svi_ref_r2.py summary
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings, traceback, itertools
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C
import c_chanorder as chanorder

FS = R.FS
O2 = f"{R.OUT}/r2"
PRE_S, LEN_S = 2.0, 16.0
NCTL = 10
NMIN_EV = 5
LAGMAX = int(2.0 * FS)
TRC_ALL = ["PGC", "VGZ", "NLLB", "SNB", "LZB", "YOUB", "PFB", "SOKB", "LRIV", "B001", "B003", "B004", "B005", "B006",
           "B007", "B009", "B010", "B011", "B926"]


# ------------------------------------------------------------------------------------------------ sensitivity
def gain_table():
    p = f"{O2}/gains.csv"
    if os.path.exists(p):
        return pd.read_csv(p, dtype={"loc": str})
    rows = []
    A = pd.read_csv(f"{R.RES}/a4_stationxml_audit.csv", dtype={"loc": str})
    for r in A.itertuples():
        rows.append(dict(sta=r.sta, cha=r.chan, start=str(r.start)[:10], end=("2100-01-01" if str(r.end) == "open" else str(r.end)[:10]),
                         gain=float(r.sens), src="results/a4_stationxml_audit.csv"))
    from obspy import read_inventory
    for path in [f"{R.PICK}/results/p4/stationxml/PB.xml", f"{R.PICK}/results/p4/stationxml/CN.xml",
                 f"{R.ROOT}/lfe_array_discovery/results/stationxml/CN.xml", f"{R.ROOT}/lfe_array_discovery/results/stationxml/UW.xml"]:
        if not os.path.exists(path): continue
        try: inv = read_inventory(path)
        except Exception: continue
        for net in inv:
            for st in net:
                if st.code not in TRC_ALL: continue
                for ch in st:
                    if ch.code[:2] not in ("HH", "BH", "EH"): continue
                    sens = ch.response.instrument_sensitivity if ch.response is not None else None
                    if sens is None or not sens.value: continue
                    rows.append(dict(sta=st.code, cha=ch.code, start=str(ch.start_date)[:10],
                                     end=str(ch.end_date)[:10] if ch.end_date else "2100-01-01", gain=float(sens.value), src=path))
    G = pd.DataFrame(rows).drop_duplicates(["sta", "cha", "start"], keep="first")
    R.ensure(O2); G.to_csv(p, index=False)
    return G


def gain_of(G, sta, cha, date):
    g = G[(G.sta == sta) & (G.cha == cha) & (G.start <= date) & (G.end > date)]
    return float(g.gain.iloc[0]) if len(g) else np.nan


# ------------------------------------------------------------------------------------------------ raw metrics
def raw_day(path, sta, date, G):
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
    d0 = pd.Timestamp(date).timestamp()
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


# ------------------------------------------------------------------------------------------------ stacks
_G = {}


def init_tables():
    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    off = C.placement(tpl, sp, F, mv, C.station_table(), L)
    T = {k: (tpl[k][L].astype(np.float64), off[k]) for k in tpl if k in off}
    return dict(T=T, L=L)


def month_dets(ym):
    p = f"{C.TRC}/detections_{ym}.parquet"
    if not os.path.exists(p):
        return pd.DataFrame(columns=["fam", "t", "date", "cat_core", "tremor_flag", "tier_core"])
    return pd.read_parquet(p, columns=["fam", "t", "date", "cat_core", "tremor_flag", "tier_core"])


def seg_of(x, t0, t, o):
    j = int(round((t + o - PRE_S - t0) * FS)); m = int(round(LEN_S * FS))
    if j < 0 or j + m > x.shape[1]:
        return None
    s = x[:, j:j + m].astype(np.float64)
    r = np.sqrt(np.mean(s * s))
    if not np.isfinite(r) or r <= 0:
        return None
    return s / r


def job(args):
    sta, ns, ym, dates = args
    t00 = time.time()
    G = _G["gains"]; T = _G["T"]
    daily, status = [], []
    man = _G["man"]
    for ds in dates:
        mrow = man.get((ns, ds))
        p = C.trace_path(ns, int(ds[:4]), pd.Timestamp(ds).dayofyear, mrow[2] if mrow else "")
        expect = mrow is not None and mrow[0] in ("data", "present_elsewhere")
        if p is None:
            status.append((ds, ns, "missing" if expect else "nodata", "raw: no file")); continue
        try:
            daily += raw_day(p, sta, ds, G); status.append((ds, ns, "ok", ""))
        except Exception as e:
            status.append((ds, ns, "error", f"raw: {type(e).__name__}: {str(e)[:60]}"))
    # stacks
    D = _G["dets"].get(ym)
    fams = sorted({f for (f, s) in T if s == sta})
    cst = []
    X = {}
    for ds in dates:
        d = pd.Timestamp(ds)
        q = next((f"{c}/{ns}/{d.year}/{d.dayofyear:03d}.npz" for c in (R.CACHE_R1, C.CACHE)
                  if os.path.exists(f"{c}/{ns}/{d.year}/{d.dayofyear:03d}.npz")), None)
        if q is None:
            continue
        try:
            z = np.load(q); X[ds] = (z["x"].astype(np.float32), float(z["t0"])); cst.append((ds, ns, "ok", ""))
        except Exception as e:
            cst.append((ds, ns, "error", f"cache: {str(e)[:60]}"))
    m = int(round(LEN_S * FS))
    SD = {f: np.zeros((3, m)) for f in fams}; ND = {f: 0 for f in fams}
    SC = {f: np.zeros((NCTL, 3, m)) for f in fams}; NC = {f: np.zeros(NCTL, int) for f in fams}
    if D is not None and len(D) and X:
        rng = np.random.default_rng(int(ym.replace("-", "")) * 7 + sum(map(ord, sta)))
        Dc = D[D.cat_core & D.tremor_flag]
        core_by = {f: np.sort(g.t.to_numpy()) for f, g in D[D.tier_core].groupby("fam")}
        dlist = sorted(X)
        for r in Dc.itertuples():
            f = float(r.fam)
            if (f, sta) not in T or r.date not in X: continue
            tpl_, o = T[(f, sta)]
            x, t0 = X[r.date]
            s = seg_of(x, t0, r.t, o)
            if s is None: continue
            SD[f] += s; ND[f] += 1
            tod = r.t - pd.Timestamp(r.date).timestamp()
            cand = [d2 for d2 in dlist if abs((pd.Timestamp(d2) - pd.Timestamp(r.date)).days) >= 3]
            for kk in range(NCTL):
                if not cand: break
                for _try in range(4):
                    d2 = cand[rng.integers(len(cand))]
                    t2 = pd.Timestamp(d2).timestamp() + tod
                    cb = core_by.get(r.fam, np.zeros(0))
                    if len(cb) and np.min(np.abs(cb - t2)) < 3600: continue
                    s2 = seg_of(X[d2][0], X[d2][1], t2, o)
                    if s2 is None: continue
                    SC[f][kk] += s2; NC[f][kk] += 1; break
    R.ensure(f"{O2}/stacks")
    np.savez_compressed(f"{O2}/stacks/{ns}_{ym}.npz", fams=np.array(fams), SD=np.stack([SD[f] for f in fams]) if fams else np.zeros((0, 3, m)),
                        ND=np.array([ND[f] for f in fams]), SC=np.stack([SC[f] for f in fams]) if fams else np.zeros((0, NCTL, 3, m)),
                        NC=np.stack([NC[f] for f in fams]) if fams else np.zeros((0, NCTL)))
    return sta, ns, ym, daily, status, cst, time.time() - t00


def _init(G):
    warnings.filterwarnings("ignore"); _G.update(G)


def run(a):
    R.ensure(O2, f"{R.OUT}/status")
    ep = C.epochs(); man = C.manifest()
    months = [m for m in a.months.split(",") if m] or R.MONTHS
    stas = [s for s in a.stations.split(",") if s] or TRC_ALL
    nsof = R.netsta_of()
    done = {os.path.basename(p)[:-4] for p in glob.glob(f"{O2}/jobs/*.ok")}
    jobs = []
    for ym in months:
        for sta in stas:
            if sta not in nsof: continue
            ds = [str(d.date()) for d in C.month_days(ym) if len(C.alive(ep[ep.sta == sta], d)) and str(d.date()) < R.PROD_T1]
            if not ds: continue
            ns = nsof[sta]
            if f"{ns}_{ym}" in done: continue
            jobs.append((sta, ns, ym, ds))
    tabs = init_tables()
    dets = {}
    for ym in sorted({j[2] for j in jobs}):
        dets[ym] = month_dets(ym)
    G = dict(gains=gain_table(), T=tabs["T"], dets=dets,
             man={(r.netsta, r.date): (r.status, r.coverage_frac, r.note if isinstance(r.note, str) else "") for r in man.itertuples()})
    R.log(f"[r2] {len(jobs)} station-month jobs ({len(done)} done), {a.workers} workers")
    import multiprocessing as mp
    R.ensure(f"{O2}/jobs", f"{O2}/daily")
    t0 = time.time()
    with mp.Pool(a.workers, initializer=_init, initargs=(G,), maxtasksperchild=20) as pool:
        for n, (sta, ns, ym, daily, st, cst, secs) in enumerate(pool.imap_unordered(job, jobs, chunksize=1), 1):
            pd.DataFrame(daily).to_parquet(f"{O2}/daily/{ns}_{ym}.parquet", index=False)
            pd.DataFrame(st, columns=["date", "netsta", "status", "msg"]).to_csv(f"{O2}/jobs/{ns}_{ym}_raw.csv", index=False)
            pd.DataFrame(cst, columns=["date", "netsta", "status", "msg"]).to_csv(f"{O2}/jobs/{ns}_{ym}_cache.csv", index=False)
            open(f"{O2}/jobs/{ns}_{ym}.ok", "w").write(f"{secs:.1f}\n")
            if n % 50 == 0 or n == len(jobs):
                R.log(f"  {n}/{len(jobs)} {time.time() - t0:.0f}s (last {ns} {ym} {secs:.0f}s)")
    # status for the guard: raw read + cache read
    for kind in ("raw", "cache"):
        fs = sorted(glob.glob(f"{O2}/jobs/*_{kind}.csv"))
        S = pd.concat([pd.read_csv(f, dtype=str) for f in fs], ignore_index=True) if fs else pd.DataFrame(columns=["date", "netsta", "status", "msg"])
        R.write_status(f"r2{kind}", "all", S.values.tolist())


# ------------------------------------------------------------------------------------------------ summary
def lagcc(S, Tp):
    """stack S (3, m) vs template Tp (3, L) at lags -LAGMAX..LAGMAX around the template placement ->
    (mean of the per-channel cc (3-C, our convention), the A-2 component-summed cc) per lag."""
    from numpy.lib.stride_tricks import sliding_window_view
    L = Tp.shape[1]; i0 = int(PRE_S * FS)
    W = sliding_window_view(S, L, axis=1)[:, i0 - LAGMAX:i0 + LAGMAX + 1]
    Wc = W - W.mean(-1, keepdims=True); Tc = Tp - Tp.mean(1, keepdims=True)
    num = (Wc * Tc[:, None, :]).sum(-1)
    wn = np.sqrt((Wc ** 2).sum(-1)); tn = np.sqrt((Tc ** 2).sum(-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        per = num / (wn * tn[:, None])
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


def station_month(ns, sta, ym, T):
    """the stack test of one station-month (statistics described in the module docstring)."""
    p = f"{O2}/stacks/{ns}_{ym}.npz"
    if not os.path.exists(p):
        return None
    z = np.load(p)
    fams, SD, ND, SC, NC = z["fams"], z["SD"], z["ND"], z["SC"], z["NC"]
    keep = [i for i in range(len(fams)) if ND[i] >= NMIN_EV and (float(fams[i]), sta) in T]
    out = dict(n_fam=len(keep), n_ev=int(ND[keep].sum()) if keep else 0)
    if len(keep) < 3 or out["n_ev"] < 30:
        out["tested"] = False; return out
    out["tested"] = True
    tp = [T[(float(fams[i]), sta)][0] for i in keep]
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
        if v > best[0]: best = (v, (pp, sg))
    out["perm_best_cc"] = best[0]; out["perm_best"] = f"{best[1][0]}{best[1][1]}"
    out["perm_is_identity"] = best[1] == ((0, 1, 2), (1, 1, 1))
    rb = (-np.inf, 0.0)
    for th in np.radians(np.arange(0, 360, 10)):
        v = float(np.nanmedian([fam_stats(rotate(S, th), Tp)[0] for S, Tp in zip(st, tp)]))
        if v > rb[0]: rb = (v, th)
    out["rot_best_cc"] = rb[0]; out["rot_best_deg"] = float(np.degrees(rb[1]))
    # the station is at CONTROL LEVEL when either what the detector sees (cc at its fixed placement, +-0.25 s)
    # or the timing consistency (share of families whose best lag is within +-0.25 s) is at or below the
    # 95th percentile of the scrambled control (the same station read on other days)
    out["at_control"] = bool((out["cc_place"] <= out["ctl_cc_place_p95"]) or (out["lagcons"] <= out["ctl_lagcons_p95"]))
    return out


def _sm_job(k):
    return station_month(k[0], k[1], k[2], _G["T"])


def daily_flags(Dy):
    """flag channel-days against the channel's own centred 31-day running median."""
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
    T = init_tables()["T"]; nsof = R.netsta_of(); ep = C.epochs()
    fs = sorted(glob.glob(f"{O2}/daily/*.parquet"))
    Dy = pd.concat([pd.read_parquet(f) for f in fs if os.path.getsize(f) > 0], ignore_index=True) if fs else pd.DataFrame()
    Dy = daily_flags(Dy) if len(Dy) else Dy
    Dy.to_parquet(f"{O2}/r2_health_daily.parquet", index=False)
    fl = [c for c in Dy.columns if c.startswith("f_")]
    # station-day flags (any channel), meta changes (channel / loc / fs set vs previous day)
    if len(Dy):
        SDy = Dy.groupby(["sta", "date"]).agg(**{c: (c, "max") for c in fl}, seeds=("seed", lambda v: ",".join(sorted(set(v)))),
                                               fsset=("fs", lambda v: ",".join(str(x) for x in sorted(set(v))))).reset_index()
        SDy = SDy.sort_values(["sta", "date"])
        SDy["f_meta_change"] = (SDy.groupby("sta").seeds.shift() != SDy.seeds) | (SDy.groupby("sta").fsset.shift() != SDy.fsset)
        SDy.loc[SDy.groupby("sta").head(1).index, "f_meta_change"] = False
        SDy["month"] = SDy.date.str[:7]
    else:
        SDy = pd.DataFrame(columns=["sta", "date", "month"] + fl + ["f_meta_change"])
    SDy.to_parquet(f"{O2}/r2_health_stationday.parquet", index=False)
    alive_days = {}
    for sta in TRC_ALL:
        e = ep[ep.sta == sta]
        for r in e.itertuples():
            for d in pd.date_range(max(r.start, pd.Timestamp(R.PROD_T0)), min(r.end, pd.Timestamp(R.PROD_T1)), freq="D", inclusive="left"):
                alive_days.setdefault((sta, str(d.date())[:7]), set()).add(str(d.date()))
    keys = [(nsof[sta], sta, ym) for ym in R.MONTHS for sta in TRC_ALL if sta in nsof and sta in R.S15
            and (sta, ym) in alive_days and os.path.exists(f"{O2}/stacks/{nsof[sta]}_{ym}.npz")]
    import multiprocessing as mp
    with mp.Pool(a.workers, initializer=_init, initargs=(dict(T=T),)) as pool:
        SMv = pool.map(_sm_job, keys, chunksize=4)
    SM = {k[2] + k[1]: v for k, v in zip(keys, SMv)}
    Dyg = {k: g for k, g in Dy.groupby(["sta", Dy.date.str[:7]])} if len(Dy) else {}
    SDg = {k: g for k, g in SDy.groupby(["sta", "month"])} if len(SDy) else {}
    rows = []
    for ym in R.MONTHS:
        for sta in TRC_ALL:
            if sta not in nsof or (sta, ym) not in alive_days: continue
            ns = nsof[sta]
            r = dict(month=ym, sta=sta, netsta=ns, n_alive_days=len(alive_days[(sta, ym)]), template_station=sta in R.S15)
            g = SDg.get((sta, ym), SDy.iloc[0:0])
            r["n_raw_days"] = len(g)
            for c in fl + ["f_meta_change"]:
                r[f"n_{c[2:]}"] = int(g[c].sum()) if len(g) else 0
            gd = Dyg.get((sta, ym))
            for c in ("nms_2_8", "nms_01_05", "rms_2_8", "line_prom_max"):
                r[f"{c}_med"] = float(gd[c].median()) if gd is not None and c in gd else np.nan
            sm = SM.get(ym + sta)
            if sm: r.update(sm)
            rows.append(r)
    M = pd.DataFrame(rows)
    for c in ("tested", "cc_place", "ctl_cc_place_p95", "lagcons", "ctl_lagcons_p95", "perm_is_identity", "at_control",
              "perm_best", "rot_best_deg", "lag_med", "cc_a2", "cc_best"):
        if c not in M: M[c] = np.nan
    M["tested"] = M.tested.fillna(False).astype(bool)
    M["cc_med"] = M.cc_place; M["ctl_p95"] = M.ctl_cc_place_p95
    # status
    M["status"] = "ok"
    M.loc[M.n_raw_days == 0, "status"] = "nodata"
    M.loc[(M.n_raw_days > 0) & (M.n_dead > 0.5 * M.n_raw_days), "status"] = "dead"
    ctl = M.tested & M.at_control.fillna(False).astype(bool)
    M.loc[ctl & (M.status == "ok"), "status"] = "control_level"
    M.loc[~M.tested & (M.status == "ok") & M.template_station, "status"] = "untested"
    # untested: inherit the nearest tested month of the same instrument epoch (band epoch + post-2012 boundary)
    M["inst_epoch"] = [C.eff_epoch_start(s, pd.Timestamp(m + "-15"), ep) for s, m in zip(M.sta, M.month)]
    M["status_eff"] = M.status
    for (sta, ie), g in M.groupby(["sta", "inst_epoch"]):
        t = g[g.tested]
        if not len(t): continue
        tm = pd.PeriodIndex(t.month, freq="M")
        for i in g.index[g.status == "untested"]:
            d = np.abs((tm - pd.Period(M.at[i, "month"], freq="M")).map(lambda x: x.n))
            j = t.index[int(np.argmin(d))]
            M.at[i, "status_eff"] = "inherit_" + M.at[j, "status"]
    absent = M.status_eff.isin(["nodata", "dead", "control_level", "inherit_control_level", "inherit_dead"])
    M["absent"] = absent
    M.to_csv(f"{O2}/r2_health_monthly.csv", index=False)
    E = M[M.template_station].groupby("month").apply(lambda g: pd.Series(dict(
        effective=",".join(g[~g.absent].sta), absent=",".join(f"{s}:{st}" for s, st in zip(g[g.absent].sta, g[g.absent].status_eff)),
        n_effective=int((~g.absent).sum())))).reset_index()
    E.to_csv(f"{O2}/r2_effective_sets.csv", index=False)
    # exclusions for the production matched filter: template stations ABSENT for a health reason (not plain nodata,
    # not the epoch-policy LRIV drop, which the detector applies itself)
    X = M[M.template_station & M.absent & ~M.status_eff.eq("nodata")].copy()
    X = X[~((X.sta == "LRIV") & (X.month >= "2019-06"))]
    X["reason"] = X.status_eff + " cc_med=" + X.cc_med.round(3).astype(str) + " ctl_p95=" + X.ctl_p95.round(3).astype(str)
    X[["sta", "month", "reason"]].to_csv(f"{O2}/r2_exclusions.csv", index=False)
    # swap checks
    sw = []
    for sta, d in R.SWAPS:
        g = M[(M.sta == sta) & M.tested]
        b = g[g.month < d[:7]].tail(3); af = g[g.month > d[:7]].head(3)
        sw.append(dict(sta=sta, date=d, months_before=",".join(b.month), months_after=",".join(af.month),
                       cc_before=float(b.cc_med.median()) if len(b) else np.nan, ctl95_before=float(b.ctl_p95.median()) if len(b) else np.nan,
                       cc_after=float(af.cc_med.median()) if len(af) else np.nan, ctl95_after=float(af.ctl_p95.median()) if len(af) else np.nan,
                       lag_before=float(b.lag_med.median()) if len(b) else np.nan, lag_after=float(af.lag_med.median()) if len(af) else np.nan,
                       perm_after=";".join(af.perm_best.astype(str)) if len(af) else "", perm_identity_after=bool(af.perm_is_identity.all()) if len(af) else None,
                       rot_after=";".join(af.rot_best_deg.astype(str)) if len(af) else "",
                       status_after=";".join(af.status) if len(af) else ""))
    pd.DataFrame(sw).to_csv(f"{O2}/r2_swap_checks.csv", index=False)
    # positive control: LRIV after 2019-05-30 (A-2 step 4: -0.216 vs +0.542 for its matched control) + 2021-10
    L = M[(M.sta == "LRIV") & (M.month >= "2019-06") & M.tested]
    L21 = M[(M.sta == "LRIV") & (M.month == "2021-10")]
    pc = dict(n_tested_months=len(L), n_control_level=int((L.status == "control_level").sum()),
              cc_med_median=float(L.cc_med.median()) if len(L) else np.nan,
              ctl_p95_median=float(L.ctl_p95.median()) if len(L) else np.nan,
              cc_2021_10=float(L21.cc_med.iloc[0]) if len(L21) and "cc_med" in L21 and np.isfinite(L21.cc_med.iloc[0]) else np.nan,
              ctl95_2021_10=float(L21.ctl_p95.iloc[0]) if len(L21) and np.isfinite(L21.ctl_p95.iloc[0]) else np.nan,
              status_2021_10=str(L21.status.iloc[0]) if len(L21) else "", a2_reference=-0.216)
    pc["lagcons_median"] = float(L.lagcons.median()) if len(L) else np.nan
    pc["ctl_lagcons_p95_median"] = float(L.ctl_lagcons_p95.median()) if len(L) else np.nan
    pc["r2_flags_lriv_post2019"] = bool(len(L) > 0 and pc["n_control_level"] >= 0.5 * len(L))
    pa = json.load(open(f"{O2}/pc_a2.json")) if os.path.exists(f"{O2}/pc_a2.json") else {}
    pc["a2_reproduction"] = pa
    pc["pass"] = bool(pa.get("a2_inputs_LRIV", {}).get("reproduced", False))
    B9 = M[(M.sta == "B926") & M.tested]
    pc["b926_no_change_control"] = dict(n_tested=len(B9), n_ok=int((B9.status == "ok").sum()), cc_med_median=float(B9.cc_med.median()) if len(B9) else np.nan)
    json.dump(pc, open(f"{O2}/r2_positive_control.json", "w"), indent=1)
    # bad_for_dvv (KMR style): runs of station-day flags + absent months
    bad = []
    for sta, g in SDy.groupby("sta"):
        for c in ["f_dead", "f_flat", "f_noisy", "f_gap", "f_dc_jump", "f_meta_change"]:
            v = g[c].to_numpy(bool); dd = g.date.to_numpy()
            if not v.any(): continue
            e = np.diff(np.concatenate([[0], v.astype(np.int8), [0]]))
            for i, j in zip(np.flatnonzero(e == 1), np.flatnonzero(e == -1)):
                bad.append(dict(netsta=nsof.get(sta, sta), start=dd[i], end=dd[j - 1], n_days=int(j - i), reason=c[2:]))
    for r in M[M.absent & ~M.status_eff.eq("nodata")].itertuples():
        bad.append(dict(netsta=r.netsta, start=r.month + "-01", end=str((pd.Period(r.month) .end_time).date()), n_days=int(r.n_alive_days),
                        reason=f"month {r.status_eff} (stack cc {r.cc_med if hasattr(r, 'cc_med') else np.nan})"))
    pd.DataFrame(bad).to_csv(f"{O2}/r2_bad_for_dvv.csv", index=False)
    R.log(f"[r2 summary] {len(M)} station-months; status {M.status_eff.value_counts().to_dict()}; exclusions {len(X)}; "
          f"positive control {pc}")
    if a.require_pc and not pc["pass"]:
        R.log("[r2 summary] FAIL: the LRIV positive control is not reproduced"); sys.exit(6)


# ------------------------------------------------------------------------------------------------ A-2 reproduction
def pc_a2(a):
    """Positive control, Merlin R2: Track A-2 step 4 found LRIV's post-2019-05-30 stacks at cc median -0.216 against its
    2012 stacks (results/a4_epoch_identity_LRIV_new.csv).  Re-run A-2's exact procedure (a4_epoch_stack.py functions:
    predicted-S cut, align_stack, best_cc = signed 3-C cc at max |cc| over +-2 s vs the 2012 stack over S-1..S+7)
    (i) on A-2's own withheld-station detections -> must reproduce -0.216, and (ii) on the R2 event set
    (v1.1 cat_core inside tremor windows, 2021 doy 276-302) for LRIV and for B926 (no-change control)."""
    import e5r_templates as T
    import a4_epoch_stack as A4
    import d3_vel1d as V
    REPR = f"{R.REP}/results"
    stab = pd.read_csv(f"{REPR}/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{REPR}/e5r_families_2012.csv").set_index("bkey")
    Z = np.load(f"{REPR}/e5r_templates_2012.npz", allow_pickle=True)
    QC0 = pd.read_csv(f"{REPR}/e5r_template_qc_2012.csv")
    def days_of(ns, year, d0, d1):
        D = {}
        for c in (C.CACHE, f"{R.REP}/cache20hz"):
            for d in range(d0, d1 + 1):
                f = f"{c}/{ns}/{year}/{d:03d}.npz"
                if d in D or not os.path.exists(f): continue
                z = np.load(f); D[d] = (float(z["t0"]), z["x"])
        return D
    def leg(sta, ns, Dd, tmad):
        QC = QC0[(QC0.sta == sta) & QC0.usable & QC0.template].set_index("bkey")
        days = days_of(ns, 2021, 276, 302)
        nwin = int(round(2 * A4.HALF_W * FS)); rows = []
        for bkey in sorted(QC.index):
            if bkey not in F.index: continue
            g = Dd[(Dd.fam == float(bkey)) & (Dd.net_mad >= tmad)]
            if len(g) < 20: continue
            g = g.sort_values("net_mad", ascending=False).head(400)
            fr = F.loc[bkey]
            dx = T.epi_km(float(fr.lat), float(fr.lon), float(stab.loc[sta, "lat"]), float(stab.loc[sta, "lon"]))
            ttS = float(V.tt(dx, float(fr.dep), "S"))
            S, keep = T.cut_windows(days, g.t.to_numpy() + ttS - A4.HALF_W, nwin)
            if S is None or len(S) < 20: continue
            st, Sal, sh = T.align_stack(S, iters=2)
            ref = np.asarray(Z[f"{sta}|{int(bkey)}"], float)
            i_ref = A4.IDX_NOM + int(round(float(QC.loc[bkey, "s_pick_s"]) * FS)) - int(round(A4.PRE_TPL * FS))
            cc, dS = A4.best_cc(ref, st, i_ref)
            rows.append(dict(sta=sta, family=int(bkey), n_events=len(S), cc=cc, dS=dS))
        return pd.DataFrame(rows)
    out = {}
    A2 = pd.read_parquet(f"{R.RES}/e5r_detections_a4_LRIV_new.parquet")
    if "L" in A2: A2 = A2[A2.L == 12.0]
    ref_csv = pd.read_csv(f"{R.RES}/a4_epoch_identity_LRIV_new.csv"); ref_csv = ref_csv[ref_csv.status == "ok"]
    r1 = leg("LRIV", "UW.LRIV", A2, 8.75)
    out["a2_inputs_LRIV"] = dict(n_fam=len(r1), cc_median=float(r1.cc.median()) if len(r1) else np.nan,
                                 stored_a2_median=float(ref_csv.cc.median()), stored_n_fam=len(ref_csv))
    out["a2_inputs_LRIV"]["reproduced"] = bool(len(r1) and abs(out["a2_inputs_LRIV"]["cc_median"] - out["a2_inputs_LRIV"]["stored_a2_median"]) <= 0.02)
    E = pd.concat([month_dets(m) for m in ("2021-10", "2021-11")], ignore_index=True)
    E = E[E.cat_core & E.tremor_flag]
    E = E[(pd.to_datetime(E.date).dt.dayofyear >= 276) & (pd.to_datetime(E.date).dt.dayofyear <= 302)].assign(net_mad=99.0)
    r2 = leg("LRIV", "UW.LRIV", E, 0); r3 = leg("B926", "PB.B926", E, 0)
    out["r2_events_LRIV"] = dict(n_fam=len(r2), cc_median=float(r2.cc.median()) if len(r2) else np.nan)
    out["r2_events_B926"] = dict(n_fam=len(r3), cc_median=float(r3.cc.median()) if len(r3) else np.nan)
    pd.concat([r1.assign(leg="a2_inputs"), r2.assign(leg="r2_events"), r3.assign(leg="r2_events")]).to_csv(f"{O2}/pc_a2_families.csv", index=False)
    json.dump(out, open(f"{O2}/pc_a2.json", "w"), indent=1)
    R.log(f"[r2 pc_a2] {out}")
    if a.require_pc and not out["a2_inputs_LRIV"]["reproduced"]:
        R.log("[r2 pc_a2] FAIL: A-2's LRIV identity number is not reproduced from its own inputs"); sys.exit(6)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "summary", "pc_a2"])
    ap.add_argument("--workers", type=int, default=20); ap.add_argument("--months", default=""); ap.add_argument("--stations", default="")
    ap.add_argument("--require-pc", type=int, default=1)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    if a.what == "run": run(a)
    elif a.what == "pc_a2": pc_a2(a)
    else: summary(a)
