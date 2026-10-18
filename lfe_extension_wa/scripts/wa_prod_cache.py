#!/usr/bin/env python
"""WA PRODUCTION W3: 20 Hz caches on the WA chain (wa_cache.cache_one UNCHANGED: merge, band HH > BH > EH, Z / N|1 / E|2
not rotated, 1-8 Hz Butterworth-4 at the native rate, resample to 20 Hz, midnight grid, validity mask `m` dilated 5 s).
New files go to CACHE_NEW only (cache20hz_prod/); existing WA caches (cache20hz/, wa_pre/w2/cache20hz/) are read, never
written.  WA_EXTRA_NODES = the step-4b node table (as when the trial caches were built).

  plan --phase 1   every roster station-day with raw in the production roots (not traces_dl2), 2007-01-01..2026-09-15,
                   not yet cached; truth months first, then oldest first.  XU / XG before 2012 are the truth array: never.
  plan --phase 2   the traces_dl2 station-days (only once WA_DL2_DONE exists); months >= 2012-01
  run --phase N    month by month; per month status file + marker OUT/cache/done_p<N>_<YYYY-MM>; FAILS (exit 5) when > 3 %
                   of a month's station-days error (unreadable raw, crash) -- 'empty' / 'no band' = no usable data
  determinism      rebuild 20 already-cached WA station-days from raw into a scratch dir; arrays must be identical
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, shutil, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("WA_EXTRA_NODES", "/home/jovyan/tremorferometry/lfe_extension_wa/results/step4b/setup/extra_nodes.csv")
import wa_prod_common as P
W = P.patch_wa()
import wa_cache as WC

OC = f"{P.OUT}/cache"


def plan(a):
    P.ensure(OC)
    S = P.stations()
    ok_sta = set(S.netsta)
    roots = P.prod_roots()
    if a.phase == 1:
        roots = [r for r in roots if r != P.DL2]
    else:
        if not P.dl2_ready():
            raise SystemExit("phase 2 needs WA_DL2_DONE")
        roots = [P.DL2]
    from wa_prod_setup import _list_station
    from concurrent.futures import ThreadPoolExecutor
    years = [str(y) for y in range(2007, 2027)]
    jobs = [(r, ns, years) for r in roots for ns in ok_sta if os.path.isdir(f"{r}/{ns}")]
    rows = []
    with ThreadPoolExecutor(16) as ex:
        for o in ex.map(_list_station, jobs):
            rows += o
    R = pd.DataFrame(rows, columns=["netsta", "date", "root"])
    R = R[(R.date >= P.T0) & (R.date <= P.T1)]
    R = R[~(R.netsta.str.split(".").str[0].isin(["XU", "XG"]) & (R.date < "2012-01-01"))]
    dl = [f"{P.OUT}/setup/days_truth.json", f"{P.OUT}/setup/days_prod.json"]
    if P.SMOKE and all(os.path.exists(x) for x in dl):          # smoke: only the smoke days
        keep = set(json.load(open(dl[0]))) | set(json.load(open(dl[1])))
        R = R[R.date.isin(keep)]
    rk = {r: i for i, r in enumerate(P.prod_roots())}
    R["rank"] = R.root.map(rk)
    R = R.sort_values("rank").drop_duplicates(["netsta", "date"]).drop(columns=["rank"])
    if a.phase == 2:          # only days no earlier root holds (raw_path order = the production root order)
        p1 = pd.read_parquet(f"{OC}/jobs_p1_all.parquet") if os.path.exists(f"{OC}/jobs_p1_all.parquet") else pd.DataFrame(columns=["netsta", "date"])
        R = R[~R.set_index(["netsta", "date"]).index.isin(p1.set_index(["netsta", "date"]).index)]
    R.to_parquet(f"{OC}/jobs_p{a.phase}_all.parquet", index=False)
    cached = np.array([P.find_cache(ns, d) is not None for ns, d in zip(R.netsta, R.date)])
    J = R[~cached].copy()
    J["month"] = J.date.str[:7]
    J["prio"] = np.where(J.month.isin(P.TRUTH_MONTHS), 0, 1)
    J = J.sort_values(["prio", "date", "netsta"])
    J.to_parquet(f"{OC}/jobs_p{a.phase}.parquet", index=False)
    P.log(f"[cache plan p{a.phase}] {len(R)} station-days with raw ({R.netsta.nunique()} stations); already cached {int(cached.sum())}; "
          f"to build {len(J)} ({J.month.nunique()} months)")


def _job(j):
    ns, date, raw = j
    try:
        r = WC.cache_one((ns, date, raw))
    except Exception as e:
        r = dict(netsta=ns, date=date, status=f"crash:{type(e).__name__}:{str(e)[:80]}")
    return r


def run(a):
    J = pd.read_parquet(f"{OC}/jobs_p{a.phase}.parquet")
    months = list(dict.fromkeys(J.month.tolist()))
    allm = sorted({d[:7] for d in P.days()})
    if a.phase == 2:
        allm = [m for m in allm if m >= "2012-01"]
    if a.months:
        months = [m for m in months if m in set(a.months.split(","))]
    P.log(f"[cache p{a.phase}] {len(J)} station-days in {len(months)} months, {a.workers} workers")
    t0 = time.time(); ntot = 0
    with ProcessPoolExecutor(a.workers) as ex:
        for m in months:
            if os.path.exists(f"{OC}/done_p{a.phase}_{m}"):
                continue
            g = J[J.month == m]
            if a.max_days:
                g = g[g.date.isin(sorted(g.date.unique())[:a.max_days])]
            rows = []
            fut = [ex.submit(_job, (r.netsta, r.date, f"{r.root}/{r.netsta}/{r.date[:4]}/{pd.Timestamp(r.date).dayofyear:03d}.mseed"))
                   for r in g.itertuples()]
            for f in as_completed(fut):
                rows.append(f.result())
            M = pd.DataFrame(rows)
            M.to_csv(f"{OC}/manifest_p{a.phase}_{m}.csv", index=False)
            st = []
            for r in M.itertuples():
                s = str(r.status)
                if s in ("ok", "cached"):
                    st.append((r.date, r.netsta, "ok", ""))
                elif s.startswith("unreadable") or s.startswith("crash"):
                    st.append((r.date, r.netsta, "error", s))
                else:
                    st.append((r.date, r.netsta, "nodata", s))
            P.write_status(f"cache_p{a.phase}", m, st)
            nb = sum(1 for x in st if x[2] == "error")
            fr = nb / max(len(st), 1)
            ntot += len(st)
            if fr > P.ERR_MAX:
                P.log(f"[cache p{a.phase}] {m}: {nb}/{len(st)} errors ({fr:.1%}) > 3 % -> FAIL")
                open(f"{OC}/FAILED_p{a.phase}_{m}", "w").write(f"{nb}/{len(st)}\n")
                sys.exit(5)
            P.marker(f"cache/done_p{a.phase}_{m}", f"{len(st)} station-days, errors {nb}")
            P.log(f"[cache p{a.phase}] {m}: {len(st)} built ({M.status.value_counts().to_dict()}); {ntot} total {time.time() - t0:.0f}s")
    for m in allm:                                  # months with nothing to build are complete as well
        if not os.path.exists(f"{OC}/done_p{a.phase}_{m}") and not a.months and not a.max_days:
            P.marker(f"cache/done_p{a.phase}_{m}", "no station-day to build")
    P.log(f"[cache p{a.phase}] done {ntot} station-days {time.time() - t0:.0f}s")


def determinism(a):
    """rebuild 20 WA-cache station-days (already cached, frozen chain) from raw into a scratch dir: identical arrays?"""
    M = pd.read_csv(f"{P.RES}/cache_manifest.csv")
    M = M[M.status.isin(["ok", "cached"])]
    rng = np.random.default_rng(20261003)
    pick = M.iloc[rng.choice(len(M), min(60, len(M)), replace=False)]
    tmp = tempfile.mkdtemp(prefix="wa_det_", dir=os.environ.get("TMPDIR", "/tmp"))
    rows = []
    orig = W.cache_path
    for r in pick.itertuples():
        if len(rows) >= 20:
            break
        ref = f"{P.STUDY}/cache20hz/{r.netsta}/{r.date[:4]}/{pd.Timestamp(r.date).dayofyear:03d}.npz"
        raw = P.raw_path(r.netsta, r.date)
        if not os.path.exists(ref) or raw is None:
            continue
        W.cache_path = lambda ns, d, _t=tmp: f"{_t}/{ns}/{pd.Timestamp(d).year}/{pd.Timestamp(d).dayofyear:03d}.npz"
        out = WC.cache_one((r.netsta, r.date, raw))
        W.cache_path = orig
        p = f"{tmp}/{r.netsta}/{r.date[:4]}/{pd.Timestamp(r.date).dayofyear:03d}.npz"
        if out.get("status") != "ok" or not os.path.exists(p):
            rows.append(dict(netsta=r.netsta, date=r.date, identical=False, note=str(out.get("status")))); continue
        z0, z1 = np.load(ref), np.load(p)
        same = all(np.array_equal(z0[k], z1[k]) for k in ("x", "m", "live")) and str(z0["node"]) == str(z1["node"])
        rows.append(dict(netsta=r.netsta, date=r.date, identical=bool(same), raw=raw,
                         max_abs_dx=float(np.max(np.abs(z0["x"].astype(np.float64) - z1["x"]))), node0=str(z0["node"]), node1=str(z1["node"])))
    shutil.rmtree(tmp, ignore_errors=True)
    D = pd.DataFrame(rows)
    P.ensure(OC)
    D.to_csv(f"{OC}/determinism.csv", index=False)
    ok = len(D) >= 20 and bool(D.identical.all())
    P.log(f"[determinism] {int(D.identical.sum())}/{len(D)} station-days rebuilt identical (max |dx| {D.max_abs_dx.max() if 'max_abs_dx' in D else 'nan'})")
    open(f"{OC}/DETERMINISM.txt", "w").write(f"{'PASS' if ok else 'FAIL'} {int(D.identical.sum())}/{len(D)} identical\n")
    sys.exit(0 if ok else 5)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["plan", "run", "determinism"])
    ap.add_argument("--phase", type=int, default=1)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--months", default="")
    ap.add_argument("--max-days", type=int, default=0)
    a = ap.parse_args()
    dict(plan=plan, run=run, determinism=determinism)[a.what](a)
