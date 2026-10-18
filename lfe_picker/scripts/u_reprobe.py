#!/usr/bin/env python
"""P4-0a re-probe - Merlin's ORIGINAL clause: 6 tremor + 6 matched quiet days, ">= +3 dB on >= 4
of 6", on the geometrically decisive sites.

The 4-day / ">= 3 of 4" reading used in the first pass was our shortcut; this restores the
pre-registered one.  Nothing about the chain, the statistics, the far-tremor days or the
calibration changes - only the number of tremor days the contrast bar is evaluated on.

Order (enforced): the pre-registration (day list + station list + clause + SHA-256 of every input)
is written to results/p4/reprobe_prereg.json BEFORE the first byte is fetched.  The two extra
tremor days are chosen by PNSN epicentre count alone (>= 100 within 40 km of a P4-band active
cell), one from 2023-03/04 and one from 2024-08 or 2025-09/10, at least 3 days from any day already
used; the two extra quiet days by the same quiet rule as the first pass (no epicentre within 75 km
on the day or +- 3 d), season-matched.  No station response is looked at in either choice.

Station set (as specified): every candidate SITE that contributes to >= 10 cells of any sub-region
at <= 45 km hypocentral, plus every site currently scoring 1-2 of 4 on the native 2-8 Hz contrast.

Writes probe_days_v2.csv, probe_queue_v2.csv, station_day_stats_v2.csv, hf_station_day_stats_v2.csv,
reprobe_download_manifest.csv, reprobe_prereg.json.
"""
from __future__ import annotations
import os
import sys
import csv
import json
import time
import hashlib
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
from u_fartremor import hf_one                                         # noqa: E402
from r_download_day import fetch_day                                   # noqa: E402
from u_download import link_from_extension                             # noqa: E402

SUBREGIONS = {"R1": (47.0, 47.5), "R2": (48.0, 48.3), "R3": (47.5, 48.0)}
TREMOR_R_KM, TREMOR_MIN_N = 40.0, 100
QUIET_R_KM, QUIET_PAD_D = 75.0, 3
MIN_SEP_D = 3
CELLS_AT45_MIN = 10
WINDOWS = [("2023-03-01", "2023-04-30"), None]      # the 2nd is chosen from either late window
WINDOW_B = [("2024-08-01", "2024-08-31"), ("2025-09-01", "2025-10-31")]
MANIFEST = os.path.join(U.P4, "reprobe_download_manifest.csv")
COLS = ["net", "sta", "band", "loc", "date", "doy", "src", "status", "bytes", "n_channels",
        "fs", "coverage_frac", "n_frag", "updated_utc", "note"]


def log(m):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pick_days(used_t, used_q):
    C = U.active_cells()
    T = pd.read_csv(U.TREMOR, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time)
    T = T[T.time >= "2022-11-01"]
    T["day"] = T.time.dt.floor("D")
    T["d"] = cKDTree(U.proj(C.lat.values, C.lon.values)).query(
        U.proj(T.lat.values, T.lon.values))[0]
    days = pd.date_range("2023-01-01", "2026-05-27", freq="D")
    n40 = T[T.d <= TREMOR_R_KM].groupby("day").size().reindex(days, fill_value=0)
    n75 = T[T.d <= QUIET_R_KM].groupby("day").size().reindex(days, fill_value=0)
    qok = n75.rolling(2 * QUIET_PAD_D + 1, center=True,
                      min_periods=2 * QUIET_PAD_D + 1).max() == 0
    ut = [pd.Timestamp(d) for d in used_t]

    def far_enough(d):
        return all(abs((d - u).days) >= MIN_SEP_D for u in ut)
    new_t = []
    for win in [WINDOWS[0]] + [None]:
        if win is not None:
            cand = n40[(n40.index >= win[0]) & (n40.index <= win[1])]
        else:
            cand = pd.concat([n40[(n40.index >= a) & (n40.index <= b)] for a, b in WINDOW_B])
        cand = cand[(cand >= TREMOR_MIN_N)]
        cand = cand[[far_enough(d) and str(d.date()) not in used_t for d in cand.index]]
        if not len(cand):
            continue
        d = cand.idxmax()
        new_t.append((d, int(cand[d])))
        ut.append(d)
    qpool = [d for d in qok[qok].index if str(d.date()) not in used_q]
    new_q = []
    for d, _ in new_t:
        same = [q for q in qpool if q.year == d.year and q not in [x for x, _ in new_q]]
        pool = same if same else [q for q in qpool if q not in [x for x, _ in new_q]]
        if not pool:
            continue
        q = min(pool, key=lambda q: (min(abs((q - d).days % 365), 365 - abs((q - d).days % 365)),
                                     abs((q - d).days)))
        new_q.append((q, d))
    rows = [dict(date=str(d.date()), day_type="tremor", n_tremor_40km=n,
                 n_tremor_50km=np.nan, n_tremor_75km=int(n75[d]), source="reprobe")
            for d, n in new_t]
    rows += [dict(date=str(q.date()), day_type="quiet", n_tremor_40km=int(n40[q]),
                  n_tremor_50km=np.nan, n_tremor_75km=int(n75[q]), season_match=str(d.date()),
                  source="reprobe") for q, d in new_q]
    return pd.DataFrame(rows)


def pick_stations(C):
    """Geometrically decisive sites + the sites currently at 1-2 of 4 on the native contrast."""
    CELLS = U.active_cells()
    d = U.hyp_km(CELLS.lat.values[None, :], CELLS.lon.values[None, :],
                 CELLS.depth_km.values[None, :], C.lat.values[:, None], C.lon.values[:, None])
    keep_geo = set()
    for R, (a0, a1) in SUBREGIONS.items():
        m = ((CELLS.lat >= a0) & (CELLS.lat <= a1)).values
        n45 = (d[:, m] <= 45.0).sum(1)
        keep_geo |= set(C.site.values[n45 >= CELLS_AT45_MIN])
    V = pd.read_csv(os.path.join(U.P4, "shortlist_v1.csv")).drop_duplicates("station")
    marg = set(V[(V.n_tremor_days_ge3db_native >= 1)
                 & (V.n_tremor_days_ge3db_native <= 2)].site)
    sites = keep_geo | marg
    sel = C[C.site.isin(sites)].copy()
    return sel, sorted(keep_geo), sorted(marg)


def download(Q, workers=8):
    if not os.path.exists(MANIFEST):
        with open(MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
        done = set()
    else:
        done = {(r["net"], r["sta"], r["date"])
                for r in csv.DictReader(open(MANIFEST)) if r.get("status") in ("data", "nodata")}
    jobs = [r for r in Q if (r["net"], r["sta"], r["date"]) not in done]
    log(f"[reprobe] {len(jobs)} station-days to fetch ({len(done)} terminal)")
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
            if cnt["n"] % 40 == 0:
                log(f"  {cnt['n']}/{len(jobs)}")
        return rec
    if jobs:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(one, jobs))
    return pd.read_csv(MANIFEST, keep_default_na=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--fdsn", type=int, default=8)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t00 = time.time()
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"), keep_default_na=False,
                    dtype={"loc": str})
    PD = pd.read_csv(os.path.join(U.P4, "probe_days.csv"))
    used_t = list(PD[PD.day_type == "tremor"].date)
    used_q = list(PD[PD.day_type == "quiet"].date)
    ND = pick_days(used_t, used_q)
    sel, geo, marg = pick_stations(C)
    for r in ND.itertuples():
        log(f"  [new day] {r.date} {r.day_type:6s} n40={r.n_tremor_40km} n75={r.n_tremor_75km}")
    log(f"[reprobe] {len(sel)} stations / {sel.site.nunique()} sites "
        f"({len(geo)} geometrically decisive, {len(marg)} at 1-2 of 4)")

    # ---------------- PRE-REGISTRATION, written before anything is fetched ----------------
    prereg = dict(
        written_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        purpose="P4-0a re-probe to Merlin's original clause: 6 tremor + 6 matched quiet days, "
                "contrast >= +3 dB on >= 4 of 6 (native 2-8 Hz).",
        clause=dict(contrast_bar_db=3.0, min_days=4, n_tremor_days=6,
                    statistic="native-rate 2-8 Hz median 6-s window power, tremor day vs the "
                              "median over the station's quiet days",
                    unchanged=["predicted rho >= 0.5 from the frozen north fit",
                               "far-tremor <= +1 dB (specified) and local-far >= 3 dB (intent)",
                               "10-19 Hz contrast < 2-8 Hz contrast",
                               "delivery >= 90 % of the days in the station's own epoch",
                               "sub-region rule: >= 4 admitted sites within 45 km hypocentral of "
                               ">= 50 % of cells, >= 2 further admitted sites withheld "
                               "(withheld = site ranks 3-4 by predicted rho)"],
                    far_tremor_days="UNCHANGED from the first pass"),
        day_selection_rule=dict(tremor=">= 100 PNSN epicentres within 40 km of a 47.0-48.3 N "
                                       "tremor-active Slab2 cell; one day from 2023-03/04 and one "
                                       "from 2024-08 or 2025-09/10; >= 3 days from any day already "
                                       "used; highest epicentre count wins; no station response "
                                       "consulted",
                                quiet="no PNSN epicentre within 75 km of any active cell on the "
                                      "day or +- 3 d; same-year, nearest day-of-year to its "
                                      "tremor day"),
        station_selection_rule=f"every candidate SITE contributing to >= {CELLS_AT45_MIN} cells of "
                               f"any sub-region at <= 45 km hypocentral, plus every site currently "
                               f"at 1-2 of 4 on the native 2-8 Hz contrast",
        existing_days=dict(tremor=used_t, quiet=used_q),
        new_days=ND.to_dict("records"),
        n_stations=int(len(sel)), n_sites=int(sel.site.nunique()),
        stations=sorted((sel.net + "." + sel.sta).tolist()),
        sites_geometrically_decisive=geo, sites_marginal_1_or_2_of_4=marg,
        input_sha256={}, script_sha256={})
    for f in ["candidates.csv", "probe_days.csv", "probe_queue.csv", "station_day_stats.csv",
              "hf_station_day_stats.csv", "far_tremor_days.csv", "far_station_day_stats.csv",
              "rho_calibration.json", "north_calibration.csv", "shortlist_v1.csv",
              "subregion_summary.csv"]:
        p = os.path.join(U.P4, f)
        if os.path.exists(p):
            prereg["input_sha256"][f] = sha(p)
    for f in ["u_common.py", "u_stats.py", "u_cache.py", "u_calibrate.py", "u_fartremor.py",
              "u_shortlist_v1.py", "u_reprobe.py"]:
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
        if os.path.exists(p):
            prereg["script_sha256"][f] = sha(p)
    U.jdump(prereg, os.path.join(U.P4, "reprobe_prereg.json"))
    log(f"[reprobe] PRE-REGISTERED -> {U.P4}/reprobe_prereg.json "
        f"({len(prereg['input_sha256'])} inputs, {len(prereg['script_sha256'])} scripts hashed)")

    # ---------------- fetch ----------------
    Q = [dict(net=r.net, sta=r.sta, band=r.band, loc=r._asdict()["loc"], date=d,
              role="reprobe", site=r.site)
         for r in sel.itertuples() for d in ND.date
         if str(r.first_3c) <= d <= str(r.last_3c)]
    M = download(Q, a.fdsn)
    M = M[M.status == "data"].drop_duplicates(["net", "sta", "date"], keep="last")
    log(f"[reprobe] delivered {len(M)} of {len(Q)} station-days")

    # ---------------- cache ----------------
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
    log(f"[reprobe] caching {len(jobs)} station-days")
    import multiprocessing as mp
    with mp.Pool(a.workers, maxtasksperchild=20) as pool:
        cr = pd.DataFrame(list(pool.imap_unordered(UC.dispatch, jobs, chunksize=1)))
    log("[reprobe] cache " + str(cr.groupby("status").size().to_dict()))

    # ---------------- statistics (cache band) ----------------
    need = M[["net", "sta", "date"]].merge(sel[["net", "sta", "loc", "band"]], on=["net", "sta"])
    RESP, MODE, freqs = US.build_response_table(need, log)
    dtype = dict(zip(ND.date, ND.day_type))
    meta = {(r.net, r.sta): r for r in sel.itertuples()}
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
                      float(m.zsign), dtype[r.date], "reprobe", m.site, m.sensor_type))
    log(f"[reprobe] statistics on {len(sjobs)} station-days")
    with mp.Pool(a.workers, initializer=US._init, initargs=(RESP, {"freqs": freqs}),
                 maxtasksperchild=40) as pool:
        S = pd.DataFrame(list(pool.imap_unordered(US.one, sjobs, chunksize=1)))
    S["resp_mode"] = [MODE.get((r.netsta, f"{r.band}Z", r.date), "none") for r in S.itertuples()]
    log("[reprobe] stats " + str(S.groupby("status").size().to_dict()))

    # ---------------- native-rate 2-8 / 10-19 Hz pass ----------------
    hj = []
    for r in sel.itertuples():
        ns = f"{r.net}.{r.sta}"
        for d in ND.date:
            dd = pd.Timestamp(d)
            if os.path.exists(U.trace_path(ns, dd.year, int(dd.dayofyear))):
                hj.append((ns, d, dd.year, int(dd.dayofyear), float(r.zsign)))
    log(f"[reprobe] native-rate pass on {len(hj)} station-days")
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        H = pd.DataFrame(list(pool.imap_unordered(hf_one, hj, chunksize=1)))
    H["day_type"] = H.date.map(dtype)
    log("[reprobe] native " + str(H.groupby("status").size().to_dict()))

    # ---------------- merged v2 products ----------------
    PD2 = pd.concat([PD.assign(source="p4-0"), ND], ignore_index=True)
    PD2.to_csv(os.path.join(U.P4, "probe_days_v2.csv"), index=False)
    Q0 = pd.read_csv(os.path.join(U.P4, "probe_queue.csv"), keep_default_na=False)
    Q2 = pd.concat([Q0, pd.DataFrame(Q)[["net", "sta", "band", "loc", "date", "role", "site"]]],
                   ignore_index=True).drop_duplicates(["net", "sta", "date"])
    Q2.to_csv(os.path.join(U.P4, "probe_queue_v2.csv"), index=False)
    S0 = pd.read_csv(os.path.join(U.P4, "station_day_stats.csv"))
    pd.concat([S0, S], ignore_index=True).drop_duplicates(
        ["netsta", "date"], keep="last").sort_values(["netsta", "date"]).to_csv(
        os.path.join(U.P4, "station_day_stats_v2.csv"), index=False)
    H0 = pd.read_csv(os.path.join(U.P4, "hf_station_day_stats.csv"))
    pd.concat([H0, H], ignore_index=True).drop_duplicates(
        ["netsta", "date"], keep="last").sort_values(["netsta", "date"]).to_csv(
        os.path.join(U.P4, "hf_station_day_stats_v2.csv"), index=False)
    log(f"[reprobe] DONE in {(time.time()-t00)/60:.1f} min")


if __name__ == "__main__":
    main()
