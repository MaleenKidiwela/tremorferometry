#!/usr/bin/env python
"""P4-0a step 1 - run the P4-0 probe UNCHANGED on the 14 north stations that have a MEASURED
single-event registration rho (P3D_RESULTS.md 3.2 and A3.2).

Merlin F5: the site-noise / registration proxy is unvalidated, and the P3D class+distance model
missed SYMB by 3.1 dB.  So the proxy is calibrated here, against rho, before it is allowed to rank
a single Puget station.  Everything is computed with the SAME code as the Puget probe:
u_download.fetch_day -> u_cache.dispatch (frozen 20 Hz / 1-8 Hz chain) -> u_stats.one
(2-8 Hz velocity PSD relative to CN.PGC on the same day; incoherent 6-s energy contrast).

Days: north tremor days inside S18 / S21 / S22 with >= 100 PNSN epicentres within 40 km of the
Bostock family centroid, and matched quiet days within +-30 d with NO epicentre within 150 km
(day and +-1 d).  The SAME days for all 15 stations, so the PGC ratio is exact.

Writes results/p4/north_{stations,probe_days,download_manifest,station_day_stats}.csv and
results/p4/stationxml/NORTH_<NET>.xml.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
import u_cache as UC                                                   # noqa: E402
import u_stats as US                                                   # noqa: E402
from r_download_day import fetch_day                                   # noqa: E402
from u_download import link_from_extension                             # noqa: E402

# P3D_RESULTS.md 3.2 (V1 roster) and A3.2 (PQ).  rho_rel = rho_new / rho_PGC, same events, same N.
RHO = {
    "CN.PGC":   dict(rho_rel=1.000, n_fam=None, src="3.2 reference", ssnr=None),
    "PQ.LHLYB": dict(rho_rel=1.086, n_fam=8,  src="A3.2", ssnr=10.29),
    "PQ.DAOB":  dict(rho_rel=0.826, n_fam=33, src="A3.2", ssnr=11.12),
    "CN.SYMB":  dict(rho_rel=0.783, n_fam=64, src="3.2", ssnr=13.28),
    "CN.CLRS":  dict(rho_rel=0.773, n_fam=64, src="3.2", ssnr=12.18),
    "CN.GOBB":  dict(rho_rel=0.571, n_fam=54, src="3.2", ssnr=8.03),
    "UW.TURTL": dict(rho_rel=0.559, n_fam=33, src="3.2", ssnr=6.59),
    "UW.RNWY":  dict(rho_rel=0.488, n_fam=18, src="3.2", ssnr=7.27),
    "UW.LOPEZ": dict(rho_rel=0.398, n_fam=5,  src="3.2", ssnr=3.65),
    "UW.LRIV":  dict(rho_rel=0.366, n_fam=64, src="3.2", ssnr=5.61),
    "UW.BUCKS": dict(rho_rel=0.296, n_fam=33, src="3.2", ssnr=3.41),
    "CN.PTRF":  dict(rho_rel=0.275, n_fam=56, src="3.2", ssnr=4.07),
    "UW.OSQM":  dict(rho_rel=0.168, n_fam=18, src="3.2", ssnr=2.14),
    "UW.ODUC":  dict(rho_rel=0.158, n_fam=33, src="3.2", ssnr=2.96),
}
# VGZ is not tabulated against PGC directly; 3.2 gives rho/rho_VGZ per row, so rho_VGZ/rho_PGC is
# recoverable as the median of (rho_rel / (rho/rho_VGZ)) over the 11 V1 rows -> 0.59.  DERIVED, and
# reported separately; it is never in the primary n=14 fit.
VGZ_DERIVED = 0.590
EPISODES = ["S18", "S21", "S22"]
N_TREMOR, N_QUIET = 6, 6
TREMOR_R_KM, TREMOR_MIN = 40.0, 100
QUIET_R_KM, QUIET_PAD_D = 150.0, 1
MATCH_WINDOW_D = 30
XMLDIR = os.path.join(U.P4, "stationxml")
MANIFEST = os.path.join(U.P4, "north_download_manifest.csv")
COLS = ["net", "sta", "band", "loc", "date", "doy", "src", "status", "bytes", "n_channels",
        "fs", "coverage_frac", "n_frag", "updated_utc", "note"]


def log(m):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def family_centroid():
    F = pd.read_csv(os.path.join(U.GEO, "north_families.csv")).drop_duplicates("bkey")
    return float(F.lat.median()), float(F.lon.median()), float(F.depth_km.median()), len(F)


def fetch_meta(stations):
    """IRIS StationXML level=response for the north roster -> NORTH_<NET>.xml + the band table."""
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime, read_inventory
    os.makedirs(XMLDIR, exist_ok=True)
    by_net = {}
    for ns in stations:
        n, s = ns.split(".")
        by_net.setdefault(n, []).append(s)
    rows = []
    for net, stas in by_net.items():
        p = os.path.join(XMLDIR, f"NORTH_{net}.xml")
        if os.path.exists(p):
            inv = read_inventory(p)
        else:
            inv = Client("IRIS", timeout=300).get_stations(
                network=net, station=",".join(sorted(stas)), level="response",
                starttime=UTCDateTime("2024-01-01"), endtime=UTCDateTime("2026-06-30"))
            inv.write(p, format="STATIONXML")
        for n in inv:
            for s in n:
                for ch in s:
                    if ch.code[-1] not in "ZNE12" or ch.code[:2] not in U.BAND_PREF:
                        continue
                    if (ch.sample_rate or 0) < U.FS_CACHE:
                        continue
                    try:
                        isens = ch.response.instrument_sensitivity
                        sens, iu = float(isens.value), str(isens.input_units or "")
                    except Exception:
                        sens, iu = np.nan, ""
                    rows.append(dict(net=n.code, sta=s.code, loc=ch.location_code or "",
                                     chan=ch.code, band=ch.code[:2], fs=float(ch.sample_rate),
                                     az=float(ch.azimuth if ch.azimuth is not None else 0.0),
                                     dip=float(ch.dip if ch.dip is not None else 0.0),
                                     lat=float(ch.latitude), lon=float(ch.longitude),
                                     elev=float(ch.elevation), sensitivity=sens, input_units=iu,
                                     start=str(ch.start_date.date) if ch.start_date else "2024-01-01",
                                     end=str(ch.end_date.date) if ch.end_date else "2026-06-30"))
    R = pd.DataFrame(rows)
    pick = []
    for (net, sta), g in R.groupby(["net", "sta"]):
        best = None
        for (loc, band), gg in g.groupby(["loc", "band"]):
            slots = {c[-1] for c in gg.chan}
            if not ({"Z"} <= slots and ({"N", "1"} & slots) and ({"E", "2"} & slots)):
                continue
            if len({"N", "1"} & slots) > 1 or len({"E", "2"} & slots) > 1:
                continue
            key = (U.BAND_PREF.index(band), 0 if loc in ("", "--") else 1, loc)
            if best is None or key < best[0]:
                best = (key, loc, band, gg)
        if best is None:
            log(f"  [band] {net}.{sta}: no usable 3-C velocity band -> dropped")
            continue
        _, loc, band, gg = best
        latest = gg.sort_values("start").groupby("chan").tail(1)
        zrow = latest[latest.chan.str[-1] == "Z"]
        zdip = float(zrow.dip.iloc[0]) if len(zrow) else -90.0
        fs = float(gg.fs.max())
        pick.append(dict(net=net, sta=sta, netsta=f"{net}.{sta}", loc=loc, band=band, fs=fs,
                         chain="frozen" if abs(fs / U.FS_CACHE - round(fs / U.FS_CACHE)) < 1e-9 else "rational",
                         zsign=-1.0 if zdip > 0 else 1.0,
                         chan_az_dip=";".join(f"{r.chan}:{r.az:.1f}/{r.dip:.1f}"
                                              for r in latest.itertuples()),
                         sensitivity_med=float(np.nanmedian(latest.sensitivity.values)),
                         input_units=latest.input_units.mode().iloc[0],
                         lat=float(gg.lat.iloc[0]), lon=float(gg.lon.iloc[0]),
                         elev=float(gg.elev.iloc[0]),
                         start=str(gg.start.min()), end=str(gg.end.max())))
    return pd.DataFrame(pick), R


def pick_days(clat, clon):
    W = pd.read_csv(os.path.join(U.P4, "..", "p3d", "sse_windows.csv")) \
        if os.path.exists(os.path.join(U.P4, "..", "p3d", "sse_windows.csv")) else None
    W = W[W.id.isin(EPISODES)]
    T = pd.read_csv(U.TREMOR, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time)
    T = T[T.time >= "2024-01-01"]
    T["d_km"] = U.epi_km(T.lat.values, T.lon.values, clat, clon)
    T["day"] = T.time.dt.floor("D")
    days = pd.date_range("2024-06-01", "2026-06-30", freq="D")
    n40 = T[T.d_km <= TREMOR_R_KM].groupby("day").size().reindex(days, fill_value=0)
    n150 = T[T.d_km <= QUIET_R_KM].groupby("day").size().reindex(days, fill_value=0)
    inep = pd.Series(False, index=days)
    for r in W.itertuples():
        inep |= (days >= r.core_start) & (days <= r.core_end)
    cand = n40[(n40 >= TREMOR_MIN) & inep.values]
    tdays = list(cand.sort_values(ascending=False).index[:N_TREMOR * 3])
    # spread over the three episodes: at most 2 per episode, strongest first
    per, tsel = {}, []
    for d in sorted(tdays, key=lambda d: -n40[d]):
        ep = [r.id for r in W.itertuples() if r.core_start <= str(d.date()) <= r.core_end]
        ep = ep[0] if ep else "?"
        if per.get(ep, 0) >= 2:
            continue
        per[ep] = per.get(ep, 0) + 1
        tsel.append((d, ep))
        if len(tsel) >= N_TREMOR:
            break
    qok = n150.rolling(2 * QUIET_PAD_D + 1, center=True, min_periods=2 * QUIET_PAD_D + 1).max() == 0
    qpool = list(qok[qok].index)
    qsel, used = [], set()
    for d, ep in tsel:
        cands = [q for q in qpool if q not in used and abs((q - d).days) <= MATCH_WINDOW_D]
        if not cands:
            cands = [q for q in qpool if q not in used]
        if not cands:
            continue
        q = min(cands, key=lambda q: abs((q - d).days))
        qsel.append((q, ep, d))
        used.add(q)
    rows = [dict(date=str(d.date()), day_type="tremor", episode=ep, n_tremor_40km=int(n40[d]),
                 n_tremor_150km=int(n150[d])) for d, ep in tsel]
    rows += [dict(date=str(q.date()), day_type="quiet", episode=ep, n_tremor_40km=int(n40[q]),
                  n_tremor_150km=int(n150[q]), season_match=str(d.date())) for q, ep, d in qsel]
    return pd.DataFrame(rows).sort_values(["day_type", "date"]).reset_index(drop=True)


def download(P, PD, workers=8):
    if not os.path.exists(MANIFEST):
        with open(MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
        done = set()
    else:
        done = {(r["net"], r["sta"], r["date"])
                for r in csv.DictReader(open(MANIFEST)) if r.get("status") in ("data", "nodata")}
    jobs = [(r, d) for r in P.itertuples() for d in PD.date
            if (r.net, r.sta, d) not in done and r.start <= d <= r.end]
    log(f"[north] {len(jobs)} station-days to fetch ({len(done)} terminal), {workers} threads")
    lock = __import__("threading").Lock()
    cnt = {"n": 0}

    def one(job):
        r, ds = job
        d = date.fromisoformat(ds)
        try:
            src = link_from_extension(r.net, r.sta, r.band, d)
        except Exception as e:
            src = f"link_err:{type(e).__name__}"
        rec = fetch_day(r.net, r.sta, r.band, d, r.loc, U.TRACES, prov="IRIS", tries=3, timeout=180)
        rec.update(loc=r.loc, src=src,
                   updated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"))
        with lock:
            cnt["n"] += 1
            with open(MANIFEST, "a", newline="") as f:
                csv.DictWriter(f, COLS, extrasaction="ignore").writerow(rec)
            if cnt["n"] % 25 == 0:
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
    clat, clon, cdep, nfam = family_centroid()
    log(f"[north] Bostock family centroid {clat:.3f} {clon:.3f} depth {cdep:.1f} km "
        f"({nfam} distinct families)")
    P, R = fetch_meta(list(RHO) + ["CN.VGZ"])
    P["rho_rel"] = [RHO.get(ns, {}).get("rho_rel", VGZ_DERIVED) for ns in P.netsta]
    P["rho_src"] = [RHO.get(ns, {}).get("src", "derived from 3.2 rho/rho_VGZ column") for ns in P.netsta]
    P["rho_nfam"] = [RHO.get(ns, {}).get("n_fam", np.nan) for ns in P.netsta]
    P["ssnr_p3d"] = [RHO.get(ns, {}).get("ssnr", np.nan) for ns in P.netsta]
    P["dist_centroid_km"] = U.hyp_km(clat, clon, cdep, P.lat.values, P.lon.values)
    P["primary"] = P.netsta != "CN.VGZ"
    P.to_csv(os.path.join(U.P4, "north_stations.csv"), index=False)
    log(f"[north] {len(P)} stations, bands {P.band.value_counts().to_dict()}, "
        f"chain {P.chain.value_counts().to_dict()}, dip-inverted {int((P.zsign < 0).sum())}")
    PD = pick_days(clat, clon)
    PD.to_csv(os.path.join(U.P4, "north_probe_days.csv"), index=False)
    for r in PD.itertuples():
        log(f"  [day] {r.date} {r.day_type:6s} {r.episode} n40={r.n_tremor_40km:5d} "
            f"n150={r.n_tremor_150km:5d}")
    M = download(P, PD, a.fdsn)
    M = M[M.status == "data"].drop_duplicates(["net", "sta", "date"], keep="last")
    log(f"[north] delivered {len(M)} of {len(P) * len(PD)} station-days")

    # ---- cache (the same frozen chain) ----
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
        jobs.append((ns, d.year, int(d.dayofyear), src, U.CACHE, False, float(r.coverage_frac or 0), ch))
    jobs = sorted(set(jobs))
    log(f"[north] caching {len(jobs)} station-days")
    import multiprocessing as mp
    with mp.Pool(a.workers, maxtasksperchild=20) as pool:
        cr = pd.DataFrame(list(pool.imap_unordered(UC.dispatch, jobs, chunksize=1)))
    log("[north] cache " + str(cr.groupby("status").size().to_dict()))
    for r in cr[cr.status == "error"].head(10).itertuples():
        log(f"  CACHE ERROR {r.netsta} {r.year}.{r.doy}: {r.err}")

    # ---- statistics (the same u_stats.one) ----
    need = M[["net", "sta", "date"]].merge(P[["net", "sta", "loc", "band"]], on=["net", "sta"])
    RESP, MODE, freqs = US.build_response_table(need, log)
    dtype = dict(zip(PD.date, PD.day_type))
    meta = {r.netsta: r for r in P.itertuples()}
    sjobs = []
    for r in M.itertuples():
        ns = f"{r.net}.{r.sta}"
        d = pd.Timestamp(r.date)
        if not os.path.exists(U.cache_path(ns, d.year, int(d.dayofyear))):
            continue
        m = meta[ns]
        try:
            fs = float(r.fs)
        except (TypeError, ValueError):
            fs = float(m.fs)
        ch = "frozen" if abs(fs / U.FS_CACHE - round(fs / U.FS_CACHE)) < 1e-9 else "rational"
        sjobs.append((ns, r.date, d.year, int(d.dayofyear), m.band, m.loc, fs, ch, float(m.zsign),
                      dtype.get(r.date, "?"), "north", ns.split(".")[1], "north"))
    log(f"[north] statistics on {len(sjobs)} station-days")
    with mp.Pool(a.workers, initializer=US._init, initargs=(RESP, {"freqs": freqs}),
                 maxtasksperchild=40) as pool:
        S = pd.DataFrame(list(pool.imap_unordered(US.one, sjobs, chunksize=1)))
    S["resp_mode"] = [MODE.get((r.netsta, f"{r.band}Z", r.date), "none") for r in S.itertuples()]
    S = S.sort_values(["netsta", "date"])
    S.to_csv(os.path.join(U.P4, "north_station_day_stats.csv"), index=False)
    log("[north] stats " + str(S.groupby("status").size().to_dict()))
    log(f"[north] DONE in {(time.time()-t00)/60:.1f} min")


if __name__ == "__main__":
    main()
