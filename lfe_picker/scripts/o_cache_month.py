#!/usr/bin/env python
"""Track C step 2 — build the 20 Hz cache for one month, EXACTLY as e4_cache20.py does.

Chain (FROZEN, Merlin F6):  merge(method=1, fill_value=0) -> component order Z / N|1 / E|2 resolved by
c_chanorder (borehole EH1/EH2 -> N/E) -> Butterworth-4 zero-phase 1-8 Hz at the NATIVE rate ->
integer decimation to 20 Hz with resample_poly.  One file per station-day:
    cache20hz/<NET.STA>/<YYYY>/<DOY>.npz   x (3,N) float32, t0 float64, order (3,), fs 20.0
The only differences from e4_cache20.cache_day are (i) the input path may be our own traces/ or the
read-only lfe_catalog_repair copy named in the manifest, and (ii) per station-day noise statistics are
returned for results/trackC/daily_noise_<YYYY-MM>.parquet (Merlin Step 6 component 2: the daily
sensitivity term needs each station's daily median 8-s window energy per component).

  python scripts/c_cache_month.py 2015-12 --workers 10
"""
from __future__ import annotations
import os, sys, time, argparse, warnings
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, resample_poly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_common as C
import o_chanorder as chanorder
import o_month_ready as MR

FS = 20.0
WIN_S = 8.0


def bp(x, fs, f1=1.0, f2=8.0):
    sos = butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def med_win_energy(X, n):
    """median over the day of the n-sample window energy, per component (float64 prefix sums)."""
    c = np.concatenate([np.zeros((X.shape[0], 1)), np.cumsum(X.astype(np.float64) ** 2, axis=1)], 1)
    e = c[:, n:] - c[:, :-n]
    return np.median(e, axis=1)


def cache_one(job):
    netsta, year, doy, src, out, force, cov_man = job
    t0w = time.time()
    dst = f"{out}/{netsta}/{year}/{doy:03d}.npz"
    try:
        if os.path.exists(dst) and not force:
            z = np.load(dst)
            X = z["x"].astype(np.float64)
            me = med_win_energy(X, int(round(WIN_S * FS)))
            return dict(netsta=netsta, year=year, doy=doy, status="cached", n=X.shape[1],
                        t0=float(z["t0"]), order=",".join(list(z["order"])),
                        e8_Z=me[0], e8_N=me[1], e8_E=me[2],
                        frac_nonzero=float((np.abs(X).sum(0) > 0).mean()),
                        cov_manifest=cov_man, secs=time.time() - t0w, err="")
        from obspy import read
        st = read(src)
        # OPERATIONAL FIX (Track A-2, 2026-09-16): CN fragments carry per-fragment sampling rates like
        # 99.99999 / 100.00001 Hz, and obspy refuses to merge traces "with same ids but differing
        # sampling rates" -> 41 of 504 CN station-days in 2021 were lost outright.  Snap every rate to
        # the nearest integer when it is within 1e-4 relative (i.e. < 0.01 Hz at 100 Hz, far below one
        # sample over a day); anything further out is left alone so a genuinely different rate still
        # raises instead of being silently glued together.
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
        fs0 = float(st[0].stats.sampling_rate)
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        X = bp(X, fs0)
        q = fs0 / FS
        if abs(q - round(q)) > 1e-9:
            raise RuntimeError(f"{netsta} fs={fs0} not an integer multiple of {FS}")
        q = int(round(q))
        if q > 1:
            X = resample_poly(X, 1, q, axis=-1)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        os.makedirs(f"{out}/{netsta}/{year}", exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=X.astype(np.float32), t0=np.float64(t0), order=np.array(keep), fs=np.float64(FS))
        os.replace(tmp, dst)
        me = med_win_energy(X, int(round(WIN_S * FS)))
        return dict(netsta=netsta, year=year, doy=doy, status="built", n=X.shape[1], t0=t0,
                    order=",".join(keep), fs_native=fs0, n_fs_snapped=n_snap, e8_Z=me[0], e8_N=me[1], e8_E=me[2],
                    frac_nonzero=float((np.abs(X).sum(0) > 0).mean()), cov_manifest=cov_man,
                    secs=time.time() - t0w, err="")
    except Exception as e:                                          # never kill the month for one day
        return dict(netsta=netsta, year=year, doy=doy, status="error", n=0, t0=np.nan, order="",
                    e8_Z=np.nan, e8_N=np.nan, e8_E=np.nan, frac_nonzero=np.nan, cov_manifest=cov_man,
                    secs=time.time() - t0w, err=f"{type(e).__name__}: {e}")


def build_jobs(ym, man, ep, min_cov, force, stations=None):
    st = MR.month_state(ym, man, ep, min_cov)
    sub = man[man.date.str[:7] == ym]
    note = {(r.netsta, r.date): (r.note if isinstance(r.note, str) else "") for r in sub.itertuples()}
    jobs, skipped = [], []
    for sta, ds, band, epoch_id, cov in st["usable"]:
        if stations and sta not in stations:
            continue
        d = pd.Timestamp(ds)
        netsta = ep[ep.sta == sta].netsta.iloc[0]
        p = C.trace_path(netsta, d.year, int(d.dayofyear), note.get((netsta, ds), ""))
        if p is None:
            skipped.append((netsta, ds, "no trace file"))
            continue
        jobs.append((netsta, d.year, int(d.dayofyear), p, C.CACHE, force, cov))
    return jobs, skipped, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ym")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--stations", default="", help="restrict to these station codes (debug)")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    OUTD = f"{C.PRES}/trackO"
    os.makedirs(OUTD, exist_ok=True)
    man = C.manifest(); ep = C.epochs()
    jobs, skipped, st = build_jobs(a.ym, man, ep, a.min_cov, a.force,
                                   set(a.stations.split(",")) if a.stations else None)
    print(f"[cache {a.ym}] {len(jobs)} station-days to cache ({len(skipped)} with no trace file), "
          f"{a.workers} workers, runnable={st['runnable']}", flush=True)
    t0 = time.time()
    import multiprocessing as mp
    rows = []
    if jobs:
        with mp.Pool(a.workers, maxtasksperchild=40) as pool:
            for i, r in enumerate(pool.imap_unordered(cache_one, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 50 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    D = pd.DataFrame(rows)
    if len(D):
        D["date"] = pd.to_datetime(D.year.astype(str), format="%Y") + pd.to_timedelta(D.doy - 1, "D")
        D["date"] = D.date.dt.strftime("%Y-%m-%d")
        D = D.sort_values(["netsta", "doy"])
        D.to_parquet(f"{OUTD}/daily_noise_{a.ym}.parquet", index=False)
        bad = D[D.status == "error"]
        print(D.groupby("status").size().to_string())
        for r in bad.itertuples():
            print(f"  ERROR {r.netsta} {r.year}.{r.doy:03d}: {r.err}", flush=True)
        print(f"[cache {a.ym}] done in {time.time()-t0:.0f}s -> {OUTD}/daily_noise_{a.ym}.parquet")
    else:
        print(f"[cache {a.ym}] nothing to do")
    for s in skipped[:10]:
        print(f"  SKIP {s[0]} {s[1]}: {s[2]}")


if __name__ == "__main__":
    main()
