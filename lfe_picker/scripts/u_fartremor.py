#!/usr/bin/env python
"""P4-0a step 2 - the FAR-TREMOR control and the 10-19 Hz control for the Puget candidates.

Merlin's bar: a station whose 2-8 Hz "tremor contrast" is as large for distant tremor as for local
tremor is responding to something other than the local source (tremor-season weather, regional
microseism, a scheduled cultural cycle), so far-tremor contrast must be <= +1 dB; and the contrast
in 10-19 Hz must be LOWER than in 2-8 Hz, because tremor is a low-frequency source.

FAR-TREMOR DAYS - the specified construction (reuse north tremor days on which Puget had none) is
geometrically IMPOSSIBLE: the Bostock family field sits at 48.50 N, ~25 km north of the top of the
P4 band, so every north tremor day has tremor within 50 km of a Puget active cell (checked over
2023-01-01..2026-05-27: ZERO such days).  The control is therefore built from genuinely distant
tremor at Merlin's own 100-200 km scale: days with >= 100 PNSN epicentres 120-350 km from the
nearest Puget active cell and NOT ONE within 120 km.  Both azimuths are represented (north
Vancouver Island and Oregon), so a station cannot pass by being blind to one direction.

The 2-8 Hz far contrast is computed with the SAME chain as everything else (download -> frozen
20 Hz cache -> u_stats.one) against the SAME quiet-day baseline.  The 10-19 Hz statistic cannot
come from that cache (it is band-limited to 1-8 Hz), so it is measured on the RAW native-rate
traces of the 4 local tremor and 4 quiet days already on disk, with the 2-8 Hz statistic recomputed
in the same pass as an internal check.

Writes results/p4/far_tremor_days.csv, far_station_day_stats.csv, hf_station_day_stats.csv,
control_summary.csv
"""
from __future__ import annotations
import os
import sys
import csv
import time
import argparse
import warnings
from datetime import date, datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
import u_cache as UC                                                   # noqa: E402
import u_stats as US                                                   # noqa: E402
import o_chanorder as chanorder                                        # noqa: E402
from r_download_day import fetch_day                                   # noqa: E402
from u_download import link_from_extension                             # noqa: E402

FAR_MIN, FAR_MAX = 120.0, 350.0
FAR_MIN_N = 100
N_FAR = 4
HF_BAND = (10.0, 19.0)
MANIFEST = os.path.join(U.P4, "far_download_manifest.csv")
COLS = ["net", "sta", "band", "loc", "date", "doy", "src", "status", "bytes", "n_channels",
        "fs", "coverage_frac", "n_frag", "updated_utc", "note"]


def log(m):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def pick_far_days():
    C = U.active_cells()
    T = pd.read_csv(U.TREMOR, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time)
    T = T[T.time >= f"{U.PROBE_Y0}-01-01"]
    T["day"] = T.time.dt.floor("D")
    T["d"] = cKDTree(U.proj(C.lat.values, C.lon.values)).query(
        U.proj(T.lat.values, T.lon.values))[0]
    days = pd.date_range(f"{U.PROBE_Y0}-01-01", "2026-05-27", freq="D")
    band = T[(T.d >= FAR_MIN) & (T.d <= FAR_MAX)].groupby("day").size().reindex(days, fill_value=0)
    near = T[T.d < FAR_MIN].groupby("day").size().reindex(days, fill_value=0)
    ok = band[(band >= FAR_MIN_N) & (near.values == 0)].sort_values(ascending=False)
    rows, north, south = [], 0, 0
    for d in ok.index:
        g = T[(T.day == d) & (T.d >= FAR_MIN) & (T.d <= FAR_MAX)]
        lat = float(g.lat.median())
        side = "north" if lat > 48.3 else "south"
        if side == "north" and north >= 2:
            continue
        if side == "south" and south >= 2:
            continue
        north += side == "north"
        south += side == "south"
        rows.append(dict(date=str(d.date()), day_type="far_tremor", n_far=int(ok[d]),
                         lat_med=lat, lon_med=float(g.lon.median()),
                         dist_med_km=float(g.d.median()), side=side))
        if len(rows) >= N_FAR:
            break
    return pd.DataFrame(rows)


def download(Q, workers=8):
    if not os.path.exists(MANIFEST):
        with open(MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
        done = set()
    else:
        done = {(r["net"], r["sta"], r["date"])
                for r in csv.DictReader(open(MANIFEST)) if r.get("status") in ("data", "nodata")}
    jobs = [r for r in Q if (r["net"], r["sta"], r["date"]) not in done]
    log(f"[far] {len(jobs)} station-days to fetch ({len(done)} terminal)")
    lock = __import__("threading").Lock()
    cnt = {"n": 0}

    def one(r):
        d = date.fromisoformat(r["date"])
        try:
            src = link_from_extension(r["net"], r["sta"], r["band"], d)
        except Exception as e:
            src = f"link_err:{type(e).__name__}"
        rec = fetch_day(r["net"], r["sta"], r["band"], d, r["loc"], U.TRACES,
                        prov="IRIS", tries=3, timeout=180)
        rec.update(loc=r["loc"], src=src,
                   updated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"))
        with lock:
            cnt["n"] += 1
            with open(MANIFEST, "a", newline="") as f:
                csv.DictWriter(f, COLS, extrasaction="ignore").writerow(rec)
            if cnt["n"] % 50 == 0:
                log(f"  {cnt['n']}/{len(jobs)}")
        return rec
    if jobs:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(one, jobs))
    return pd.read_csv(MANIFEST, keep_default_na=False)


# ------------------------------------------------------------------ the native-rate HF pass
def hf_one(job):
    netsta, date_s, y, doy, zsign = job
    rec = dict(netsta=netsta, date=date_s, status="error", err="")
    try:
        from obspy import read
        from scipy.signal import butter, sosfiltfilt
        src = U.trace_path(netsta, y, doy)
        st = read(src)
        for tr in st:
            fs_ = float(tr.stats.sampling_rate)
            fr = round(fs_)
            if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = [c for c in chanorder.resolve_order(chans, station=netsta, strict=True) if c]
        fs0 = float(st.select(channel=order[0])[0].stats.sampling_rate)
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in order)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in order])
        X[0] *= zsign
        rec["fs_native"] = fs0
        rec["nyq_ok"] = bool(fs0 / 2.0 > HF_BAND[1] * 1.02)
        n = int(round(U.E_WIN_S * fs0))
        # four disjoint 1.5 h blocks spread over the day (06/12/18/00 UTC): 6 h of data keeps the
        # day-wide sampling of an episodic source while cutting the native-rate filtfilt cost 4x.
        blk = int(round(1.5 * 3600 * fs0))
        starts = [int(round(h * 3600 * fs0)) for h in (0, 6, 12, 18)]
        blocks = [X[:, s0:s0 + blk] for s0 in starts if s0 + blk <= X.shape[1]]
        if not blocks:
            blocks = [X]
        rec["n_blocks"] = len(blocks)
        for lab, (f1, f2) in [("lf", (2.0, 8.0)), ("hf", HF_BAND)]:
            if f2 >= fs0 / 2.0:
                rec[f"e6_{lab}_3c"] = np.nan
                continue
            sos = butter(4, [f1 / (fs0 / 2), f2 / (fs0 / 2)], btype="band", output="sos")
            ee, nn = [], 0
            for B in blocks:
                Y = sosfiltfilt(sos, B, axis=-1)
                e, nw = U.med_win_power(Y, n)
                if np.isfinite(e).all():
                    ee.append(e)
                    nn += nw
            rec[f"e6_{lab}_3c"] = float(np.nanmedian(np.nanmedian(np.array(ee), axis=0))) if ee else np.nan
            rec[f"n_win_{lab}"] = nn
        rec["status"] = "ok"
    except Exception as e:
        rec["err"] = f"{type(e).__name__}: {e}"
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--fdsn", type=int, default=8)
    ap.add_argument("--skip-hf", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t00 = time.time()
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"), keep_default_na=False,
                    dtype={"loc": str})
    C["loc"] = C["loc"].fillna("")
    FD = pick_far_days()
    FD.to_csv(os.path.join(U.P4, "far_tremor_days.csv"), index=False)
    for r in FD.itertuples():
        log(f"  [far-day] {r.date} n_far={r.n_far:4d} median {r.dist_med_km:.0f} km "
            f"({r.side}, {r.lat_med:.2f} N)")

    Q = [dict(net=r.net, sta=r.sta, band=r.band, loc=r["loc"], date=d)
         for _, r in C.iterrows() for d in FD.date
         if str(r.first_3c) <= d <= str(r.last_3c)]
    M = download(Q, a.fdsn)
    M = M[M.status == "data"].drop_duplicates(["net", "sta", "date"], keep="last")
    log(f"[far] delivered {len(M)} of {len(Q)} station-days")

    jobs = []
    for r in M.itertuples():
        ns = f"{r.net}.{r.sta}"
        d = pd.Timestamp(r.date)
        src = U.trace_path(ns, d.year, int(d.dayofyear))
        if not os.path.exists(src):
            continue
        try:
            fs = float(r.fs)
        except (TypeError, ValueError):
            fs = 100.0
        ch = "frozen" if abs(fs / U.FS_CACHE - round(fs / U.FS_CACHE)) < 1e-9 else "rational"
        jobs.append((ns, d.year, int(d.dayofyear), src, U.CACHE, False,
                     float(r.coverage_frac or 0), ch))
    jobs = sorted(set(jobs))
    log(f"[far] caching {len(jobs)} station-days")
    import multiprocessing as mp
    with mp.Pool(a.workers, maxtasksperchild=20) as pool:
        cr = pd.DataFrame(list(pool.imap_unordered(UC.dispatch, jobs, chunksize=1)))
    log("[far] cache " + str(cr.groupby("status").size().to_dict()))

    need = M[["net", "sta", "date"]].merge(C[["net", "sta", "loc", "band"]], on=["net", "sta"])
    RESP, MODE, freqs = US.build_response_table(need, log)
    meta = {(r.net, r.sta): r for r in C.itertuples()}
    sjobs = []
    for r in M.itertuples():
        ns = f"{r.net}.{r.sta}"
        d = pd.Timestamp(r.date)
        if not os.path.exists(U.cache_path(ns, d.year, int(d.dayofyear))):
            continue
        m = meta[(r.net, r.sta)]
        try:
            fs = float(r.fs)
        except (TypeError, ValueError):
            fs = float(m.fs)
        ch = "frozen" if abs(fs / U.FS_CACHE - round(fs / U.FS_CACHE)) < 1e-9 else "rational"
        sjobs.append((ns, r.date, d.year, int(d.dayofyear), m.band, m._asdict()["loc"], fs, ch,
                      float(m.zsign), "far_tremor", "candidate", m.site, m.sensor_type))
    log(f"[far] statistics on {len(sjobs)} station-days")
    with mp.Pool(a.workers, initializer=US._init, initargs=(RESP, {"freqs": freqs}),
                 maxtasksperchild=40) as pool:
        S = pd.DataFrame(list(pool.imap_unordered(US.one, sjobs, chunksize=1)))
    S["resp_mode"] = [MODE.get((r.netsta, f"{r.band}Z", r.date), "none") for r in S.itertuples()]
    S.sort_values(["netsta", "date"]).to_csv(
        os.path.join(U.P4, "far_station_day_stats.csv"), index=False)
    log("[far] stats " + str(S.groupby("status").size().to_dict()))

    # ------------------------------------------------------------------ HF pass
    if not a.skip_hf:
        PD = pd.read_csv(os.path.join(U.P4, "probe_days.csv"))
        hj = []
        for r in C.itertuples():
            ns = f"{r.net}.{r.sta}"
            for d in PD.date:
                dd = pd.Timestamp(d)
                if os.path.exists(U.trace_path(ns, dd.year, int(dd.dayofyear))):
                    hj.append((ns, d, dd.year, int(dd.dayofyear), float(r.zsign)))
        log(f"[hf] native-rate 10-19 Hz pass on {len(hj)} station-days")
        with mp.Pool(a.workers, maxtasksperchild=8) as pool:
            H = pd.DataFrame(list(pool.imap_unordered(hf_one, hj, chunksize=1)))
        H["day_type"] = H.date.map(dict(zip(PD.date, PD.day_type)))
        H.sort_values(["netsta", "date"]).to_csv(
            os.path.join(U.P4, "hf_station_day_stats.csv"), index=False)
        log("[hf] " + str(H.groupby("status").size().to_dict()))
    log(f"[far] DONE in {(time.time()-t00)/60:.1f} min")


if __name__ == "__main__":
    main()
