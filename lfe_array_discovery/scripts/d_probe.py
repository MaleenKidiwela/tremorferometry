#!/usr/bin/env python
"""Step 4 - array-era site probe (P4-0 statistic, NATIVE rate).

Statistic REUSED verbatim: `lfe_picker/scripts/u_fartremor.hf_one` - native sample rate,
4 disjoint 1.5 h blocks (00/06/12/18 UTC), zero-phase Butterworth-4 in 2-8 Hz and in the
10-19 Hz control, median over the day of the mean square in 6 s windows, per component, 3-C
median.  `u_common.TRACES` is pointed at this study's traces/ so nothing is re-downloaded and
no lfe_picker file is touched.  The contrast is a per-station RATIO of the same statistic on two
sets of days, so it is scale-invariant: the instrument response and the digitiser gain cancel
exactly (they are the same channel epoch on both day sets - asserted per station).

2010 (N-Olympic): 6 tremor days (PNSN >= 100 epicentres within 40 km of the array centroid)
vs 6 matched quiet days drawn from the frozen T2B quiet set, plus the S-WA-only "far" days
2010-07-25 and 2010-09-03..09-07.  Admit at >= 4/6 days >= +3 dB AND local-far >= +3 dB.

2008 (Puget/S-Olympic): NO tremor catalog exists (PNSN starts 2010-01-06), so the 6 tremor days
come from the array's OWN hourly-energy flag (Merlin Q3 fallback), computed on the 20 Hz cache:
per station the hourly median 6 s window power in the cache band, divided by that station's median
quiet-day hourly value; the network median of that ratio per hour; an hour is "hot" when it
exceeds the p95 of the same network-median series over the literature-quiet days; a day is a
tremor day at >= 6 hot hours.  It is an ENERGY FLAG, not a catalog, and is recorded as such.
There are no far-tremor days in 2008 (no catalog) -> the local-far clause is UNTESTABLE there and
admission is reported on the contrast clause alone.

  python scripts/d_probe.py flag2008 [--workers 24]
  python scripts/d_probe.py probe2010 [--workers 20]
  python scripts/d_probe.py probe2008 [--workers 20]
"""
from __future__ import annotations
import os
import sys
import json
import time
import argparse
import warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "lfe_picker", "scripts")))
import d_common as D                                                   # noqa: E402
import u_common as U                                                   # noqa: E402
U.TRACES = D.TRACES                                                    # study-local traces
U.CACHE = D.CACHE
import u_fartremor as FT                                               # noqa: E402

FAR_DAYS_2010 = ["2010-07-25", "2010-09-03", "2010-09-04", "2010-09-05", "2010-09-06",
                 "2010-09-07"]
TREMOR_R_KM, TREMOR_MIN_N = 40.0, 100
N_TREMOR, N_QUIET = 6, 6
CONTRAST_DB, MIN_DAYS, FAR_MARGIN_DB = 3.0, 4, 3.0
W2010 = ("2010-07-20", "2010-09-10")
W2008 = ("2008-04-15", "2008-06-10")
NOLY = (47.8, 48.3)


def _init():
    U.TRACES = D.TRACES
    U.CACHE = D.CACHE


def db(x):
    x = np.asarray(x, float)
    return 10.0 * np.log10(np.where(x > 0, x, np.nan))


def delivered(win=None, days=None):
    M = pd.read_csv(D.MANIFEST, keep_default_na=False, low_memory=False)
    M = M[(M.status == "data") & (M.outroot == D.TRACES)]
    M = M.drop_duplicates(["net", "sta", "day"], keep="last")
    if win:
        M = M[(M.day >= win[0]) & (M.day <= win[1])]
    if days is not None:
        M = M[M.day.isin(days)]
    return M


# ------------------------------------------------------------------ 2010 day selection
def days_2010():
    F = D.station_flags()
    S = pd.read_csv(os.path.join(D.RES, "sites_array_era_final.csv"))
    xg = S[S.net3c == "XG"]
    clat, clon = float(xg.lat.mean()), float(xg.lon.mean())
    T = D.tremor_catalog()
    T = T[(T.time >= "2010-06-01") & (T.time < "2010-10-01")]
    T["d"] = D.epi_km(clat, clon, T.lat.values, T.lon.values)
    T["day"] = T.time.dt.strftime("%Y-%m-%d")
    n40 = T[T.d <= TREMOR_R_KM].groupby("day").size()
    cand = n40[n40 >= TREMOR_MIN_N]
    cand = cand[(cand.index >= W2010[0]) & (cand.index <= W2010[1])]
    tdays = sorted(cand.sort_values(ascending=False).head(N_TREMOR).index)
    fz = json.load(open(os.path.join(D.RES, "freeze_windows.json")))
    qpool = fz["quiet_days"]["y2010"]["days"]
    Mq = delivered(days=qpool)
    nsta = Mq.groupby("day").size()
    qdays = sorted(nsta.sort_values(ascending=False).head(N_QUIET).index)
    far = [d for d in FAR_DAYS_2010 if len(delivered(days=[d]))]
    return dict(tremor=tdays, quiet=qdays, far=far, centroid=[clat, clon],
                n40=n40.reindex(sorted(set(tdays + far))).to_dict(),
                rule=f"PNSN >= {TREMOR_MIN_N} epicentres within {TREMOR_R_KM:g} km of the XG "
                     f"array centroid {clat:.3f}/{clon:.3f}")


# ------------------------------------------------------------------ 2008 hourly energy flag
def hour_energy(job):
    netsta, y, doy, date_s = job
    try:
        z = np.load(os.path.join(D.CACHE, netsta, str(y), f"{doy:03d}.npz"))
        X = z["x"].astype(np.float64)
        fs = float(z["fs"])
        nh = int(round(3600 * fs))
        n6 = int(round(6.0 * fs))
        out = []
        for h in range(24):
            B = X[:, h * nh:(h + 1) * nh]
            if B.shape[1] < nh // 2:
                out.append(np.nan)
                continue
            e, nw = U.med_win_power(B, n6)
            out.append(float(np.nanmedian(e)) if np.isfinite(e).all() else np.nan)
        return dict(netsta=netsta, date=date_s, status="ok",
                    **{f"h{h:02d}": out[h] for h in range(24)})
    except Exception as e:
        return dict(netsta=netsta, date=date_s, status="error", err=f"{type(e).__name__}: {e}")


def flag2008(workers):
    fz = json.load(open(os.path.join(D.RES, "freeze_windows.json")))
    q = fz["quiet_days"]["y2008"]["days"]
    M = delivered(win=W2008)
    Mq = delivered(days=q)
    A = pd.concat([M, Mq], ignore_index=True).drop_duplicates(["net", "sta", "day"])
    F = D.station_flags()
    ok3 = set((F[F.three_c].net + "." + F[F.three_c].sta))
    A["netsta"] = A.net + "." + A.sta
    A = A[A.netsta.isin(ok3)]
    jobs = []
    for r in A.itertuples():
        d = pd.Timestamp(r.day)
        p = os.path.join(D.CACHE, r.netsta, str(d.year), f"{int(d.dayofyear):03d}.npz")
        if os.path.exists(p):
            jobs.append((r.netsta, d.year, int(d.dayofyear), r.day))
    D.log(f"[flag2008] {len(jobs)} cached station-days, {workers} workers")
    import multiprocessing as mp
    rows, t0 = [], time.time()
    with mp.Pool(workers, initializer=_init, maxtasksperchild=50) as pool:
        for i, r in enumerate(pool.imap_unordered(hour_energy, jobs, chunksize=4), 1):
            rows.append(r)
            if i % 1000 == 0:
                D.log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    H = pd.DataFrame(rows)
    H = H[H.status == "ok"]
    hc = [f"h{h:02d}" for h in range(24)]
    base = H[H.date.isin(q)].groupby("netsta")[hc].median().median(axis=1)   # quiet-day median
    R = H.copy()
    R[hc] = H[hc].div(H.netsta.map(base), axis=0)
    net_med = R.groupby("date")[hc].median()
    qser = net_med.loc[net_med.index.isin(q)].values.ravel()
    p95 = float(np.nanpercentile(qser, 95))
    hot = (net_med > p95).sum(axis=1)
    F2 = pd.DataFrame(dict(date=net_med.index, n_hot_hours=hot.values,
                           net_med_max=net_med.max(axis=1).values,
                           is_quiet=[d in q for d in net_med.index]))
    F2["tremor_flag"] = (F2.n_hot_hours >= 6) & (~F2.is_quiet)
    F2 = F2.sort_values("date")
    F2.to_csv(os.path.join(D.RES, "hour_energy_flag_2008.csv"), index=False)
    H.to_csv(os.path.join(D.RES, "hour_energy_2008_raw.csv"), index=False)
    td = F2[F2.tremor_flag].sort_values("n_hot_hours", ascending=False).head(N_TREMOR)
    D.log(f"[flag2008] p95 of the quiet network-median hourly ratio = {p95:.3f}; "
          f"{int(F2.tremor_flag.sum())} flagged tremor days; top {N_TREMOR}: "
          f"{sorted(td.date)}")
    qd = sorted(F2[F2.is_quiet].sort_values(["n_hot_hours", "net_med_max"]).head(N_QUIET).date)
    # ZK ran only 2008-03..05 and Z3/YW are partial, so the literature-quiet pools (2008-02/03,
    # 2008-07) give those stations NO baseline at all.  Add the 6 quietest IN-WINDOW days by the
    # array's own flag (0 hot hours, lowest network-median ratio) as a fallback baseline, and
    # record per station which baseline it actually used.
    inwin = F2[(~F2.is_quiet) & (F2.n_hot_hours == 0)].sort_values("net_med_max")
    qd2 = sorted(inwin.head(N_QUIET).date)
    D.jdump(dict(tremor=sorted(td.date), quiet=qd, quiet_inwindow=qd2, p95=p95,
                 rule="ARRAY HOURLY-ENERGY FLAG (no tremor catalog before 2010-01-06): hourly "
                      "median 6 s window power in the 20 Hz cache band, per station divided by "
                      "that station's median quiet-day hourly value, network median per hour, "
                      "hour hot when > p95 of the same series on the literature-quiet days, "
                      "tremor day at >= 6 hot hours",
                 n_flagged=int(F2.tremor_flag.sum())),
            os.path.join(D.RES, "probe_days_2008.json"))


# ------------------------------------------------------------------ the probe itself
def run_probe(tag, tdays, qdays, fdays, stations, zsign, workers, qdays2=()):
    qdays2 = list(qdays2)
    days = sorted(set(tdays) | set(qdays) | set(fdays) | set(qdays2))
    M = delivered(days=days)
    M["netsta"] = M.net + "." + M.sta
    M = M[M.netsta.isin(stations)]
    jobs = []
    for r in M.itertuples():
        d = pd.Timestamp(r.day)
        p = os.path.join(D.TRACES, r.netsta, str(d.year), f"{int(d.dayofyear):03d}.mseed")
        if os.path.exists(p):
            jobs.append((r.netsta, r.day, d.year, int(d.dayofyear),
                         float(zsign.get(r.netsta, 1.0))))
    D.log(f"[probe {tag}] {len(jobs)} station-days ({len(set(j[0] for j in jobs))} stations x "
          f"{len(days)} days), {workers} workers")
    import multiprocessing as mp
    rows, t0 = [], time.time()
    with mp.Pool(workers, initializer=_init, maxtasksperchild=20) as pool:
        for i, r in enumerate(pool.imap_unordered(FT.hf_one, jobs, chunksize=1), 1):
            rows.append(r)
            if i % 250 == 0:
                D.log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    R = pd.DataFrame(rows)
    R["day_type"] = np.where(R.date.isin(tdays), "tremor",
                             np.where(R.date.isin(qdays), "quiet",
                                      np.where(R.date.isin(qdays2), "quiet_inwindow",
                                               np.where(R.date.isin(fdays), "far", "?"))))
    R.to_csv(os.path.join(D.RES, f"probe_station_days_{tag}.csv"), index=False)
    ok = R[(R.status == "ok") & np.isfinite(R.get("e6_lf_3c", np.nan))]
    out = []
    for ns, g in ok.groupby("netsta"):
        qlit = g[g.day_type == "quiet"].e6_lf_3c
        qwin = g[g.day_type == "quiet_inwindow"].e6_lf_3c
        q = qlit if len(qlit) >= 3 else (qwin if len(qwin) >= 3 else pd.concat([qlit, qwin]))
        baseline = ("literature_quiet" if len(qlit) >= 3 else
                    ("inwindow_quiet" if len(qwin) >= 3 else "mixed_or_short"))
        t = g[g.day_type == "tremor"].e6_lf_3c
        f = g[g.day_type == "far"].e6_lf_3c
        qb = float(np.nanmedian(q)) if len(q) else np.nan
        cd = db(np.nanmedian(t) / qb) if (len(t) and np.isfinite(qb)) else np.nan
        per_day = db(t.values / qb) if (len(t) and np.isfinite(qb)) else np.array([])
        fd = db(np.nanmedian(f) / qb) if (len(f) and np.isfinite(qb)) else np.nan
        gq = g[g.day_type.isin(["quiet", "quiet_inwindow"])]
        qh = gq.get("e6_hf_3c", pd.Series(dtype=float))
        th = g[g.day_type == "tremor"].get("e6_hf_3c", pd.Series(dtype=float))
        hf = (db(np.nanmedian(th) / np.nanmedian(qh))
              if len(th) and len(qh) and np.nanmedian(qh) > 0 else np.nan)
        out.append(dict(netsta=ns, net=ns.split(".")[0], sta=ns.split(".")[1],
                        n_quiet=int(len(q)), n_quiet_lit=int(len(qlit)),
                        n_quiet_inwindow=int(len(qwin)), baseline=baseline,
                        n_tremor=int(len(t)), n_far=int(len(f)),
                        fs_native=float(np.nanmedian(g.fs_native)),
                        contrast_db=cd, n_days_ge3db=int(np.nansum(per_day >= CONTRAST_DB)),
                        far_contrast_db=fd, local_minus_far_db=cd - fd,
                        hf_contrast_db=hf, hf_below_lf=bool(np.isfinite(hf) and hf < cd),
                        per_day_db=";".join(f"{v:.2f}" for v in per_day)))
    A = pd.DataFrame(out)
    A["pass_contrast"] = (A.n_days_ge3db >= MIN_DAYS) & (A.n_quiet >= 3)
    A["baseline_ok"] = A.n_quiet >= 3
    if fdays:
        A["pass_far_margin"] = A.local_minus_far_db >= FAR_MARGIN_DB
        A["admit"] = A.pass_contrast & A.pass_far_margin
    else:
        A["pass_far_margin"] = np.nan
        A["admit"] = A.pass_contrast          # local-far UNTESTABLE without a catalog
    A = A.sort_values(["net", "sta"])
    A.to_csv(os.path.join(D.RES, f"site_probe_{tag}.csv"), index=False)
    return A


def station_set(win):
    F = D.station_flags()
    F["netsta"] = F.net + "." + F.sta
    ok3 = set(F[F.three_c].netsta)
    M = delivered(win=win)
    M["netsta"] = M.net + "." + M.sta
    have = set(M.netsta) & ok3
    z = {r.netsta: float(r.zsign) for r in F.itertuples()}
    return have, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["flag2008", "probe2010", "probe2008"])
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    if a.mode == "flag2008":
        flag2008(a.workers)
        return
    if a.mode == "probe2010":
        d = days_2010()
        D.jdump(d, os.path.join(D.RES, "probe_days_2010.json"))
        D.log(f"[probe2010] tremor {d['tremor']}  quiet {d['quiet']}  far {d['far']}")
        stns, z = station_set(W2010)
        A = run_probe("2010", d["tremor"], d["quiet"], d["far"], stns, z, a.workers)
    else:
        d = json.load(open(os.path.join(D.RES, "probe_days_2008.json")))
        D.log(f"[probe2008] tremor(flag) {d['tremor']}  quiet {d['quiet']}")
        stns, z = station_set(W2008)
        A = run_probe("2008", d["tremor"], d["quiet"], [], stns, z, a.workers,
                      qdays2=d.get("quiet_inwindow", []))
    D.log(f"[probe] {int(A.admit.sum())}/{len(A)} stations admitted")
    print(A.groupby("net").agg(n=("admit", "size"), admitted=("admit", "sum"),
                               med_contrast_db=("contrast_db", "median")).to_string())


if __name__ == "__main__":
    main()
