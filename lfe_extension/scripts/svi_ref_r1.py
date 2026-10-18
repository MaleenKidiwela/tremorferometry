#!/usr/bin/env python
"""SVI_REF R1 -- data-completeness audit of the v1.1 chain, 2014-01 .. 2026-09 (pattern: lfe_south duc_fix/DIAGNOSIS.md,
where 69 % of station-days were silently lost to a GPU cuFFT plan-cache OOM inside a stage that 'never kills the month').

Per station-day (every station alive in results/station_epochs.csv, end exclusive):
  expected   manifest row 'data' (coverage >= 0.5) or 'present_elsewhere'          (results/download_manifest.csv)
  raw        a readable-path raw day file exists (c_common.trace_path)
  cache      the 20 Hz cache file exists in cache20hz/ (or a previous R1 rebuild in cache20hz_r1/), its zip directory
             lists x / t0 (no full read), and daily_noise_<M>.parquet does not record an error for it
  mf         template stations only (the 15 with 2012 templates), not excluded by the epoch policy (LRIV >= 2019-05-30):
             the station appears in that date's family_days_<M>.parquet network string (= the matched filter loaded it)
  masked     per station-MONTH: windows masked by the robust energy floor (run_<M>.json) / (loaded days x windows per day)
A stage-month FAILS when > 3 % of its expected station-days are lost at that stage.  Failing CACHE days with a raw file are
REBUILT first (c_cache_month.cache_one, the frozen chain, written to cache20hz_r1/: nothing in cache20hz/ is touched) and
re-audited.  MF-stage losses are re-run anyway by the SVI_REF feature pass (its own > 3 % guard).  Energy-floor masking is a
property of the data (dead / gappy channels), not of the processing path: it is listed and handed to R2, never 'rebuilt'.
Outputs OUT/r1/: r1_stationdays.parquet, r1_completeness.csv (station x month x stage), r1_failing_days.csv,
r1_rebuild.csv, r1_summary.json.  Exit 5 when a cache/raw stage-month still fails after the rebuild.

  python scripts/svi_ref_r1.py --workers 16
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings, zipfile
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

O1 = f"{R.OUT}/r1"
WIN_PER_DAY = int(86400 * R.FS) - int(12 * R.FS)


def cache_ok(p):
    try:
        with zipfile.ZipFile(p) as z:
            n = set(z.namelist())
        return {"x.npy", "t0.npy"} <= n
    except Exception:
        return False


_A = {}


def _audit_station(args):
    """one station epoch row -> station-day rows (file checks are NFS-bound, so stations run in parallel)."""
    sta, ns, start, end = args
    G = _A
    rows = []
    for d in pd.date_range(max(start, pd.Timestamp(R.PROD_T0)), min(end, pd.Timestamp(R.PROD_T1)), freq="D", inclusive="left"):
        ds = str(d.date())
        st, cov, note = G["mm"].get((ns, ds), ("absent", np.nan, ""))
        try: c = float(cov)
        except (TypeError, ValueError): c = np.nan
        expected = st in ("data", "present_elsewhere") and (not np.isfinite(c) or c >= 0.5)
        raw = C.trace_path(ns, d.year, d.dayofyear, note) is not None
        cp = next((q for q in (f"{R.CACHE_R1}/{ns}/{d.year}/{d.dayofyear:03d}.npz", C.cache_path(ns, d.year, d.dayofyear))
                   if os.path.exists(q)), None)
        cok = cp is not None and cache_ok(cp) and (cp.startswith(R.CACHE_R1) or G["dn"].get((ns, ds), ("", ""))[0] != "error")
        tmpl = sta in R.S15
        es = C.eff_epoch_start(sta, d, G["ep"]) if tmpl else ""
        dropped = bool(G["pol"].get((sta, es), {}).get("drop", False)) if tmpl else False
        rows.append(dict(netsta=ns, sta=sta, date=ds, month=ds[:7], manifest=st, coverage=c, expected=expected, raw=raw,
                         cache=cok, cache_path=cp or "", cache_r1=bool(cp and cp.startswith(R.CACHE_R1)),
                         template=tmpl, policy_drop=dropped, mf_loaded=sta in G["loaded"].get(ds, set()),
                         dn_err=G["dn"].get((ns, ds), ("", ""))[1] or ""))
    return rows


def _ainit(G):
    warnings.filterwarnings("ignore"); _A.update(G)


def audit(workers=16):
    ep = C.epochs(); man = C.manifest()
    mm = {(r.netsta, r.date): (r.status, r.coverage_frac, r.note if isinstance(r.note, str) else "") for r in man.itertuples()}
    loaded = {}
    for m in R.MONTHS:
        p = f"{C.TRC}/family_days_{m}.parquet"
        if not os.path.exists(p): continue
        F = pd.read_parquet(p, columns=["date", "epoch"]).drop_duplicates()
        for r in F.itertuples():
            loaded.setdefault(r.date, set()).update(x.split(":")[0] for x in str(r.epoch).split("|") if x)
    dn = {}
    for m in R.MONTHS:
        p = f"{C.TRC}/daily_noise_{m}.parquet"
        if not os.path.exists(p): continue
        X = pd.read_parquet(p, columns=["netsta", "date", "status", "err"])
        for r in X.itertuples(): dn[(r.netsta, r.date)] = (r.status, r.err)
    G = dict(mm=mm, loaded=loaded, dn=dn, ep=ep, pol=C.epoch_policy())
    jobs = []
    for r in ep.itertuples():           # split each epoch into years so the pool stays busy
        for y in range(max(r.start.year, 2014), min(r.end.year, 2026) + 1):
            a_ = max(r.start, pd.Timestamp(f"{y}-01-01")); b_ = min(r.end, pd.Timestamp(f"{y + 1}-01-01"))
            if a_ < b_: jobs.append((r.sta, r.netsta, a_, b_))
    import multiprocessing as mp
    with mp.Pool(workers, initializer=_ainit, initargs=(G,)) as pool:
        out = pool.map(_audit_station, jobs, chunksize=1)
    D = pd.DataFrame([x for o in out for x in o]).drop_duplicates(["netsta", "date"], keep="last")
    return D


def stage_table(D):
    out = []
    masked = {}
    for m in R.MONTHS:
        p = f"{C.TRC}/run_{m}.json"
        if os.path.exists(p):
            masked[m] = json.load(open(p)).get("n_windows_masked_by_energy_floor", {}) or {}
    for (ns, m), g in D.groupby(["netsta", "month"]):
        e = g[g.expected]
        r = dict(netsta=ns, sta=g.sta.iloc[0], month=m, n_alive=len(g), n_expected=len(e))
        r["raw_lost"] = int((~e.raw).sum())
        r["cache_lost"] = int((e.raw & ~e.cache).sum())
        f = e[e.template & ~e.policy_drop & e.cache]
        r["mf_expected"] = len(f); r["mf_lost"] = int((~f.mf_loaded).sum())
        nl = int(g.mf_loaded.sum())
        mw = float(masked.get(m, {}).get(g.sta.iloc[0], 0))
        r["masked_windows"] = mw; r["masked_frac"] = mw / (nl * WIN_PER_DAY) if nl else np.nan
        for s_ in ("raw", "cache"):
            r[f"{s_}_frac"] = r[f"{s_}_lost"] / max(r["n_expected"], 1)
        r["mf_frac"] = r["mf_lost"] / max(r["mf_expected"], 1)
        r["fail_raw"] = r["raw_frac"] > R.ERR_MAX; r["fail_cache"] = r["cache_frac"] > R.ERR_MAX
        r["fail_mf"] = r["mf_frac"] > R.ERR_MAX; r["masked_gt3"] = bool(np.isfinite(r["masked_frac"]) and r["masked_frac"] > R.ERR_MAX)
        out.append(r)
    return pd.DataFrame(out)


def rebuild(D, workers):
    import c_cache_month as CM
    man = C.manifest()
    note = {(r.netsta, r.date): (r.note if isinstance(r.note, str) else "") for r in man.itertuples()}
    B = D[D.expected & D.raw & ~D.cache]
    jobs = []
    for r in B.itertuples():
        d = pd.Timestamp(r.date)
        src = C.trace_path(r.netsta, d.year, d.dayofyear, note.get((r.netsta, r.date), ""))
        cov = r.coverage if np.isfinite(r.coverage) else 1.0
        jobs.append((r.netsta, d.year, d.dayofyear, src, R.CACHE_R1, False, cov))
    R.log(f"[r1] rebuilding {len(jobs)} cache station-days into {R.CACHE_R1}")
    res = []
    if jobs:
        import multiprocessing as mp
        with mp.Pool(workers, maxtasksperchild=40) as pool:
            res = list(pool.imap_unordered(CM.cache_one, jobs, chunksize=1))
    # data-path fix for 'Can not merge traces with same ids but differing sampling rates' (a stray low-rate fragment,
    # e.g. a 1-Hz BHE piece in a 40-Hz day): keep only the traces at the channel's dominant rate in a copy under OUT/r1/fixed_raw
    # (the original file is untouched) and rebuild from that copy.  Every such fix is recorded in r1_rebuild.csv.
    fixed = []
    for i, r in enumerate(res):
        if r["status"] != "error" or "differing sampling rates" not in str(r["err"]): continue
        j = next(jj for jj in jobs if jj[0] == r["netsta"] and jj[1] == r["year"] and jj[2] == r["doy"])
        from obspy import read, Stream
        st = read(j[3])
        keep = Stream(); dropped = []
        for cid in sorted({tr.id for tr in st}):
            trs = [tr for tr in st if tr.id == cid]
            rates = pd.Series([round(float(tr.stats.sampling_rate)) for tr in trs]); w = pd.Series([tr.stats.npts for tr in trs])
            dom = w.groupby(rates).sum().idxmax()
            for tr in trs:
                (keep.append(tr) if round(float(tr.stats.sampling_rate)) == dom else dropped.append(f"{tr.id}@{tr.stats.sampling_rate}Hz:{tr.stats.npts}"))
        q = f"{O1}/fixed_raw/{r['netsta']}/{r['year']}/{r['doy']:03d}.mseed"; os.makedirs(os.path.dirname(q), exist_ok=True)
        keep.write(q, format="MSEED")
        r2 = CM.cache_one((j[0], j[1], j[2], q, j[4], True, j[6]))
        r2["fix"] = "dropped " + ";".join(dropped); r2["src_fixed"] = q
        res[i] = r2; fixed.append(r2)
    Rb = pd.DataFrame(res)
    Rb.to_csv(f"{O1}/r1_rebuild.csv", index=False)
    return Rb


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=16); ap.add_argument("--no-rebuild", action="store_true")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    R.ensure(O1)
    t0 = time.time()
    D = audit(a.workers)
    S = stage_table(D)
    S.to_csv(f"{O1}/r1_completeness_before.csv", index=False)
    nb = int((D.expected & D.raw & ~D.cache).sum())
    R.log(f"[r1] {len(D)} station-days audited ({time.time() - t0:.0f}s); expected {int(D.expected.sum())}; raw missing "
          f"{int((D.expected & ~D.raw).sum())}; cache missing/unreadable {nb}; mf not loaded "
          f"{int((D.expected & D.template & ~D.policy_drop & D.cache & ~D.mf_loaded).sum())}")
    Rb = pd.DataFrame()
    if nb and not a.no_rebuild:
        Rb = rebuild(D, a.workers)
        D = audit(a.workers); S = stage_table(D)
    D.to_parquet(f"{O1}/r1_stationdays.parquet", index=False)
    S.to_csv(f"{O1}/r1_completeness.csv", index=False)
    lost = D[(D.expected & ~D.raw) | (D.expected & D.raw & ~D.cache) | (D.expected & D.template & ~D.policy_drop & D.cache & ~D.mf_loaded)].copy()
    lost["stage"] = np.where(~lost.raw, "raw", np.where(~lost.cache, "cache", "mf"))
    fm = S[S.fail_raw | S.fail_cache | S.fail_mf | S.masked_gt3]
    lost = lost.merge(fm[["netsta", "month"]].assign(month_fails=True), on=["netsta", "month"], how="left")
    lost["month_fails"] = lost.month_fails.fillna(False).astype(bool)
    lost.to_csv(f"{O1}/r1_failing_days.csv", index=False)
    summ = dict(n_station_days=len(D), n_expected=int(D.expected.sum()),
                raw_missing=int((D.expected & ~D.raw).sum()), cache_missing_after=int((D.expected & D.raw & ~D.cache).sum()),
                cache_rebuilt_ok=int((Rb.status == "built").sum()) if len(Rb) else 0,
                cache_rebuild_errors=int((Rb.status == "error").sum()) if len(Rb) else 0,
                mf_not_loaded=int((D.expected & D.template & ~D.policy_drop & D.cache & ~D.mf_loaded).sum()),
                months_fail_raw=int(S.fail_raw.sum()), months_fail_cache=int(S.fail_cache.sum()), months_fail_mf=int(S.fail_mf.sum()),
                station_months_masked_gt3=int(S.masked_gt3.sum()),
                masked_gt3_list=S[S.masked_gt3][["sta", "month", "masked_frac"]].round(3).values.tolist()[:200],
                fail_raw_list=S[S.fail_raw][["sta", "month", "raw_frac"]].round(3).values.tolist(),
                fail_cache_list=S[S.fail_cache][["sta", "month", "cache_frac"]].round(3).values.tolist(),
                fail_mf_list=S[S.fail_mf][["sta", "month", "mf_frac"]].round(3).values.tolist(),
                note="cache/raw stage failures stop the queue; mf-stage losses are re-run by the SVI_REF feature pass (own guard); "
                     "energy-floor masking is a data property, listed for R2")
    json.dump(summ, open(f"{O1}/r1_summary.json", "w"), indent=1, default=str)
    R.log(f"[r1] {json.dumps({k: v for k, v in summ.items() if not k.endswith('_list')})}")
    if summ["months_fail_cache"] or summ["months_fail_raw"]:
        R.log("[r1] FAIL: station-months still lose > 3 % of their expected station-days at the raw / cache stage")
        sys.exit(5)


if __name__ == "__main__":
    main()
