#!/usr/bin/env python
"""P4-0 step 3 - cache the probe days to the FROZEN 20 Hz / 1-8 Hz chain.

Integer-decimation stations (fs = 40/100/200 Hz) go through `o_cache_month.cache_one` VERBATIM -
merge(method=1, fill_value=0) -> component order via o_chanorder (Z, N|1, E|2) -> Butterworth-4
zero-phase 1-8 Hz at the NATIVE rate -> resample_poly(1, fs/20).  Bit-identical to every earlier
run.

The 22 candidates at 50 / 125 Hz cannot be decimated by an integer, and the frozen `cache_one`
raises on them.  `cache_one_rational` is the same function with the single change that the
decimation is resample_poly(up, down) for the exact rational fs/20 (50 Hz -> 2/5, 125 Hz -> 8/50);
the filter, the merge, the component order and the output format are untouched.  Every row carries
`chain` so nothing built this way is ever mistaken for a frozen-chain product.

DIP is honoured downstream, not here: the cache stores RAW traces in Z/N/E order (as it always has),
and results/p4/candidates.csv carries `zsign` (-1 for the 4 stations whose vertical is declared at
dip +90, i.e. positive DOWN) which u_stats applies on read - the same contract as s_common.zsign.

Output: cache20hz/<NET.STA>/<YYYY>/<DOY>.npz, log results/p4/cache_probe.parquet
"""
from __future__ import annotations
import os
import sys
import time
import argparse
import warnings
from fractions import Fraction
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
import o_cache_month as OC                                             # noqa: E402
import o_chanorder as chanorder                                        # noqa: E402

FS = OC.FS


def cache_one_rational(job):
    """o_cache_month.cache_one with a rational (up/down) decimation to 20 Hz."""
    netsta, year, doy, src, out, force, cov_man = job
    t0w = time.time()
    dst = f"{out}/{netsta}/{year}/{doy:03d}.npz"
    try:
        if os.path.exists(dst) and not force:
            z = np.load(dst)
            X = z["x"].astype(np.float64)
            me = OC.med_win_energy(X, int(round(OC.WIN_S * FS)))
            return dict(netsta=netsta, year=year, doy=doy, status="cached", n=X.shape[1],
                        t0=float(z["t0"]), order=",".join(list(z["order"])), fs_native=np.nan,
                        e8_Z=me[0], e8_N=me[1], e8_E=me[2],
                        frac_nonzero=float((np.abs(X).sum(0) > 0).mean()),
                        cov_manifest=cov_man, secs=time.time() - t0w, err="")
        from obspy import read
        st = read(src)
        n_snap = 0
        for tr in st:
            fs_ = float(tr.stats.sampling_rate)
            fr = round(fs_)
            if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
                n_snap += 1
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
        keep = [c for c in order if c is not None]
        fs0 = float(st.select(channel=keep[0])[0].stats.sampling_rate)
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        X = OC.bp(X, fs0)
        fr = Fraction(fs0 / FS).limit_denominator(1000)
        up, down = fr.denominator, fr.numerator
        if abs(float(fr) - fs0 / FS) > 1e-9:
            raise RuntimeError(f"{netsta}: fs={fs0} not rational vs {FS}")
        if (up, down) != (1, 1):
            X = OC.resample_poly(X, up, down, axis=-1)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        os.makedirs(f"{out}/{netsta}/{year}", exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=X.astype(np.float32), t0=np.float64(t0), order=np.array(keep),
                 fs=np.float64(FS))
        os.replace(tmp, dst)
        me = OC.med_win_energy(X, int(round(OC.WIN_S * FS)))
        return dict(netsta=netsta, year=year, doy=doy, status="built", n=X.shape[1], t0=t0,
                    order=",".join(keep), fs_native=fs0, n_fs_snapped=n_snap,
                    e8_Z=me[0], e8_N=me[1], e8_E=me[2],
                    frac_nonzero=float((np.abs(X).sum(0) > 0).mean()), cov_manifest=cov_man,
                    secs=time.time() - t0w, err="")
    except Exception as e:
        return dict(netsta=netsta, year=year, doy=doy, status="error", n=0, t0=np.nan, order="",
                    fs_native=np.nan, e8_Z=np.nan, e8_N=np.nan, e8_E=np.nan, frac_nonzero=np.nan,
                    cov_manifest=cov_man, secs=time.time() - t0w, err=f"{type(e).__name__}: {e}")


def dispatch(job):
    netsta, year, doy, src, out, force, cov_man, chain = job
    j = job[:7]
    r = OC.cache_one(j) if chain == "frozen" else cache_one_rational(j)
    r["chain"] = chain
    return r


def build_jobs(force, min_cov, stations=None):
    M = pd.read_csv(os.path.join(U.P4, "download_manifest.csv"), keep_default_na=False)
    M = M[M.status == "data"].drop_duplicates(["net", "sta", "date"], keep="last")
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"))
    ref = os.path.join(U.P4, "reference_stations.csv")
    if os.path.exists(ref):
        C = pd.concat([C, pd.read_csv(ref)], ignore_index=True)
    chain = {(r.net, r.sta): r.chain for r in C.itertuples()}
    jobs, skipped = [], []
    for r in M.itertuples():
        netsta = f"{r.net}.{r.sta}"
        if stations and netsta not in stations:
            continue
        try:
            cov = float(r.coverage_frac or 0)
        except ValueError:
            cov = 0.0
        if cov < min_cov:
            skipped.append((netsta, r.date, f"coverage {cov:.3f} < {min_cov}"))
            continue
        d = pd.Timestamp(r.date)
        src = U.trace_path(netsta, d.year, int(d.dayofyear))
        if not os.path.exists(src):
            skipped.append((netsta, r.date, "no trace file"))
            continue
        ch = chain.get((r.net, r.sta), "frozen")
        fs = float(r.fs or 0)
        if fs > 0:                                      # trust the file over the metadata
            q = fs / FS
            ch = "frozen" if abs(q - round(q)) < 1e-9 else "rational"
        jobs.append((netsta, d.year, int(d.dayofyear), src, U.CACHE, force, cov, ch))
    return sorted(set(jobs)), skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--min-cov", type=float, default=0.2)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--stations", default="")
    ap.add_argument("--out", default=os.path.join(U.P4, "cache_probe.parquet"))
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    stations = set(a.stations.split(",")) if a.stations else None
    jobs, skipped = build_jobs(a.force, a.min_cov, stations)
    print(f"[cache] {len(jobs)} station-days ({len(skipped)} skipped), {a.workers} workers",
          flush=True)
    t0 = time.time()
    rows = []
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.workers, maxtasksperchild=20) as pool:
            for i, r in enumerate(pool.imap_unordered(dispatch, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 50 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    D = pd.DataFrame(rows)
    if len(D):
        D["date"] = (pd.to_datetime(D.year.astype(str), format="%Y")
                     + pd.to_timedelta(D.doy - 1, "D")).dt.strftime("%Y-%m-%d")
        D = D.sort_values(["netsta", "date"])
        D.to_parquet(a.out, index=False)
        print(D.groupby(["chain", "status"]).size().to_string())
        for r in D[D.status == "error"].head(25).itertuples():
            print(f"  ERROR {r.netsta} {r.date}: {r.err}", flush=True)
        print(f"[cache] {len(D)} days in {time.time()-t0:.0f}s -> {a.out}")
    for s in skipped[:20]:
        print(f"  SKIP {s[0]} {s[1]}: {s[2]}")
    if len(skipped) > 20:
        print(f"  ... {len(skipped)-20} more skipped")


if __name__ == "__main__":
    main()
