#!/usr/bin/env python
"""P4-0 step 1 - who to probe and on which 8 days.

(a) CANDIDATES: every detection-capable 3-C station (all networks, Z-only excluded, strong-motion
    and citizen-science stations flagged out by the P3D `detect_capable` rule) whose HYPOCENTRAL
    distance to at least one tremor-active Slab2 cell in 47.0-48.3 N is <= 60 km, alive in 2023+.
    Co-located stations (<1 km) are tagged with the frozen P3D site name; they are probed
    individually but count once in the geometry.
(b) METADATA: IRIS StationXML level=response for every candidate -> band/loc choice (3-C velocity
    band, highest preference, >=20 Hz), channel azimuth/DIP, overall sensitivity, epochs.
(c) DAYS: 4 strong-tremor days (>=50 PNSN epicentres within 50 km of an active cell) from 4
    different episodes of 2023-2025, and 4 quiet days (NO epicentre within 75 km of any active
    cell on the day or +-3 d), each season-matched to one tremor day.  The SAME 8 days are used
    for every station; a station that has no epoch on a day is simply dropped for that day.

Writes results/p4/{candidates.csv, chan_response.csv, station_epochs.csv, probe_days.csv,
probe_queue.csv} and results/p4/stationxml/<NET>.xml.
"""
from __future__ import annotations
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402

XMLDIR = os.path.join(U.P4, "stationxml")
LOCPREF = ["", "--", "00", "01", "10", "20"]


# ------------------------------------------------------------------ (a) candidates
def candidates():
    C = U.active_cells()
    S = pd.read_csv(os.path.join(U.GEO, "stations_all_2019plus.csv"))
    D = S[S.detect_capable & S.three_c & ~S.z_only].copy()
    d = U.hyp_km(C.lat.values[None, :], C.lon.values[None, :], C.depth_km.values[None, :],
                 D.lat.values[:, None], D.lon.values[:, None])
    D["dmin_cell_km"] = d.min(1)
    D["nearest_cell"] = C.cell_id.values[d.argmin(1)]
    D["n_cells60"] = (d <= U.DMAX_KM).sum(1)
    D["n_cells45"] = (d <= 45.0).sum(1)
    D["tw_cells60"] = (C.tremor_w.values[None, :] * (d <= U.DMAX_KM)).sum(1)
    D = D[(D.dmin_cell_km <= U.DMAX_KM) & (D.last_3c >= f"{U.PROBE_Y0}-01-01")].copy()
    sm = U.site_map()
    D["site"] = [sm.get(k, k.split(".")[-1]) for k in D.key]
    return C, D.sort_values("dmin_cell_km").reset_index(drop=True)


# ------------------------------------------------------------------ (b) metadata
def fetch_meta(D, log):
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime, read_inventory
    os.makedirs(XMLDIR, exist_ok=True)
    cli = Client("IRIS", timeout=300)
    t0, t1 = UTCDateTime(f"{U.PROBE_Y0}-01-01"), UTCDateTime(f"{U.PROBE_Y1}-12-31T23:59:59")
    rows = []
    for net, g in D.groupby("net"):
        p = os.path.join(XMLDIR, f"{net}.xml")
        if os.path.exists(p):
            inv = read_inventory(p)
            log(f"  [meta] {net}: cached {p}")
        else:
            inv = cli.get_stations(network=net, station=",".join(sorted(set(g.sta))),
                                   level="response", starttime=t0, endtime=t1)
            inv.write(p, format="STATIONXML")
            log(f"  [meta] {net}: {sum(len(s) for s in inv)} stations -> {p}")
        for n in inv:
            for s in n:
                for ch in s:
                    sens, iu = np.nan, ""
                    try:
                        isens = ch.response.instrument_sensitivity
                        sens, iu = float(isens.value), str(isens.input_units or "")
                    except Exception:
                        pass
                    rows.append(dict(
                        net=n.code, sta=s.code, loc=ch.location_code or "", chan=ch.code,
                        lat=float(ch.latitude), lon=float(ch.longitude),
                        elev=float(ch.elevation), depth_m=float(ch.depth or 0.0),
                        fs=float(ch.sample_rate or 0.0),
                        az=float(ch.azimuth) if ch.azimuth is not None else np.nan,
                        dip=float(ch.dip) if ch.dip is not None else np.nan,
                        start=str(ch.start_date.date) if ch.start_date else f"{U.PROBE_Y0}-01-01",
                        end=str(ch.end_date.date) if ch.end_date else f"{U.PROBE_Y1}-12-31",
                        sensitivity=sens, input_units=iu,
                        sensor=str(getattr(ch.sensor, "description", "") or "")[:80]))
    return pd.DataFrame(rows)


def choose_band(R, D, log):
    """Per station pick the (loc, band) with >=3 components, best band, >=20 Hz."""
    R = R[R.chan.str[-1].isin(list("ZNE12"))].copy()
    R["band"] = R.chan.str[:2]
    R = R[R.band.isin(U.BAND_PREF) & (R.fs >= U.FS_CACHE)]
    pick, epochs = [], []
    for (net, sta), g in R.groupby(["net", "sta"]):
        best = None
        for (loc, band), gg in g.groupby(["loc", "band"]):
            slots = {c[-1] for c in gg.chan}
            if not ({"Z"} <= slots and (len({"N", "1"} & slots) >= 1) and (len({"E", "2"} & slots) >= 1)):
                continue
            if len({"N", "1"} & slots) > 1 or len({"E", "2"} & slots) > 1:
                continue                                    # ambiguous slot -> chanorder would stop
            key = (U.BAND_PREF.index(band),
                   LOCPREF.index(loc) if loc in LOCPREF else 99, loc)
            if best is None or key < best[0]:
                best = (key, loc, band, gg)
        if best is None:
            log(f"  [band] {net}.{sta}: NO usable 3-C velocity band >= {U.FS_CACHE:g} Hz -> dropped")
            continue
        _, loc, band, gg = best
        fs = float(gg.fs.max())
        q = fs / U.FS_CACHE
        chain = "frozen" if abs(q - round(q)) < 1e-9 else "rational"
        latest = gg.sort_values("start").groupby("chan").tail(1)
        azd = ";".join(f"{r.chan}:{0.0 if np.isnan(r.az) else r.az:.1f}/"
                       f"{0.0 if np.isnan(r.dip) else r.dip:.1f}"
                       for r in latest.sort_values(
                           "chan", key=lambda s: s.str[-1].map({"Z": 0, "N": 1, "1": 1, "E": 2, "2": 2})).itertuples())
        zrow = latest[latest.chan.str[-1] == "Z"]
        zdip = float(zrow.dip.iloc[0]) if len(zrow) and not np.isnan(zrow.dip.iloc[0]) else -90.0
        sens = float(np.nanmedian(latest.sensitivity.values)) if len(latest) else np.nan
        pick.append(dict(net=net, sta=sta, loc=loc, band=band, fs=fs, chain=chain,
                         chan_az_dip=azd, zsign=-1.0 if zdip > 0 else 1.0,
                         sensitivity_med=sens,
                         input_units=(latest.input_units.mode().iloc[0] if len(latest) else ""),
                         sensor_desc=(latest.sensor.mode().iloc[0] if len(latest) else ""),
                         n_comp=int(latest.chan.nunique())))
        for (loc2, band2), g2 in gg.groupby(["loc", "band"]):
            epochs.append(dict(net=net, sta=sta, band=band2, fs=fs, loc=loc2,
                               start=str(g2.start.min()), end=str(g2.end.max()),
                               lat=float(g2.lat.iloc[0]), lon=float(g2.lon.iloc[0]),
                               elev=float(g2.elev.iloc[0]), n_comp=int(g2.chan.nunique()),
                               chan_az_dip=azd, mode="3C"))
    return pd.DataFrame(pick), pd.DataFrame(epochs)


def alive_on(R, net, sta, loc, band, day):
    g = R[(R.net == net) & (R.sta == sta) & (R["loc"] == loc) & (R.chan.str[:2] == band)]
    g = g[(g.start <= day) & (g.end >= day)]
    return g.chan.nunique() >= 3


# ------------------------------------------------------------------ (c) days
def pick_days(C, R, P, log):
    T = pd.read_csv(U.TREMOR, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time)
    T = T[(T.time >= f"{U.PROBE_Y0 - 1}-11-01") & (T.time < f"{U.PROBE_Y1 + 1}-01-01")]
    from scipy.spatial import cKDTree
    tree = cKDTree(U.proj(C.lat.values, C.lon.values))
    T["dmin"] = tree.query(U.proj(T.lat.values, T.lon.values))[0]
    T["day"] = T.time.dt.floor("D")
    days = pd.date_range(f"{U.PROBE_Y0}-01-01", f"{U.PROBE_Y1}-12-31", freq="D")
    n50 = T[T.dmin <= U.TREMOR_RADIUS_KM].groupby("day").size().reindex(days, fill_value=0)
    n75 = T[T.dmin <= U.QUIET_RADIUS_KM].groupby("day").size().reindex(days, fill_value=0)
    quiet_ok = n75.rolling(2 * U.QUIET_PAD_D + 1, center=True, min_periods=2 * U.QUIET_PAD_D + 1).max() == 0

    def n_alive(day):
        ds = str(day.date())
        return sum(alive_on(R, r.net, r.sta, r.loc, r.band, ds) for r in P.itertuples())

    strong = n50[n50 >= U.TREMOR_MIN_N]
    epi, prev = [], None
    for d in strong.index:
        if prev is None or (d - prev).days > 7:
            epi.append([])
        epi[-1].append(d)
        prev = d
    epi = sorted(epi, key=lambda e: -max(n50[d] for d in e))[:U.N_TREMOR_DAYS]
    tdays = []
    for e in sorted(epi, key=lambda e: e[0]):
        top = sorted(e, key=lambda d: -n50[d])[:3]
        tdays.append(max(top, key=lambda d: (n_alive(d), n50[d])))
    qpool = list(quiet_ok[quiet_ok].index)
    qdays, used = [], set()
    for td in tdays:            # match the SAME YEAR first (station lifetimes), then day-of-year
        free = [d for d in qpool if d not in used]
        same = [d for d in free if d.year == td.year]
        pool = same if same else free
        cand = sorted(pool, key=lambda d: (min(abs((d - td).days % 365),
                                               365 - abs((d - td).days % 365)),
                                           abs((d - td).days)))[:5]
        if not cand:
            continue
        pickd = max(cand, key=lambda d: (n_alive(d), -abs((d - td).days)))
        qdays.append(pickd)
        used.add(pickd)
    rows = []
    for d in tdays:
        rows.append(dict(date=str(d.date()), day_type="tremor", n_tremor_50km=int(n50[d]),
                         n_tremor_75km=int(n75[d]), n_stations_alive=n_alive(d)))
    for d, td in zip(qdays, tdays):
        rows.append(dict(date=str(d.date()), day_type="quiet", n_tremor_50km=int(n50[d]),
                         n_tremor_75km=int(n75[d]), n_stations_alive=n_alive(d),
                         season_match=str(td.date())))
    PD = pd.DataFrame(rows).sort_values(["day_type", "date"]).reset_index(drop=True)
    for r in PD.itertuples():
        log(f"  [day] {r.date} {r.day_type:6s} n50={r.n_tremor_50km:4d} n75={r.n_tremor_75km:4d} "
            f"alive={r.n_stations_alive}")
    return PD


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=U.P4)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    lines = []

    def log(m):
        print(m, flush=True)
        lines.append(m)

    C, D = candidates()
    log(f"[p4-0] active cells 47.0-48.3 N: {len(C)}  |  candidate stations <= {U.DMAX_KM:g} km "
        f"hypocentral, alive {U.PROBE_Y0}+: {len(D)} in {D.net.nunique()} networks, "
        f"{D.site.nunique()} sites")
    log("      " + str(D.groupby(['net']).size().to_dict()))
    log("      " + str(D.groupby(['sensor_type']).size().to_dict()))
    R = fetch_meta(D, log)
    R.to_csv(os.path.join(a.out, "chan_response.csv"), index=False)
    P, E = choose_band(R, D, log)
    P = P.merge(D[["net", "sta", "key", "lat", "lon", "elev", "sensor_type", "site",
                   "dmin_cell_km", "nearest_cell", "n_cells45", "n_cells60", "tw_cells60",
                   "in_roster183", "first_3c", "last_3c"]], on=["net", "sta"], how="inner")
    P["zdip_inverted"] = P.zsign < 0
    P.to_csv(os.path.join(a.out, "candidates.csv"), index=False)
    E.to_csv(os.path.join(a.out, "station_epochs.csv"), index=False)
    log(f"[p4-0] usable candidates: {len(P)}  bands={P.band.value_counts().to_dict()}  "
        f"chain={P.chain.value_counts().to_dict()}  "
        f"dip-inverted Z: {int(P.zdip_inverted.sum())}  "
        f"no sensitivity: {int(P.sensitivity_med.isna().sum())}")
    PD = pick_days(C, R, P, log)
    PD.to_csv(os.path.join(a.out, "probe_days.csv"), index=False)

    # --------- the queue: every candidate x every day it has an epoch for, plus the references
    q = []
    for r in P.itertuples():
        for d in PD.date:
            if alive_on(R, r.net, r.sta, r.loc, r.band, d):
                q.append(dict(net=r.net, sta=r.sta, band=r.band, loc=r.loc, date=d,
                              role="candidate", site=r.site))
    ref = pd.DataFrame(REF_META(log))
    for r in ref.itertuples():
        for d in PD.date:
            q.append(dict(net=r.net, sta=r.sta, band=r.band, loc=r.loc, date=d,
                          role="reference", site=r.sta))
    Q = pd.DataFrame(q)
    refset = {(r["net"], r["sta"]) for r in ref.to_dict("records")}
    Q["role"] = [("reference" if (r.net, r.sta) in refset else r.role) for r in Q.itertuples()]
    Q["loc"] = Q["loc"].fillna("")
    Q = Q.drop_duplicates(["net", "sta", "date"]).reset_index(drop=True)
    Q.to_csv(os.path.join(a.out, "probe_queue.csv"), index=False)
    ref.to_csv(os.path.join(a.out, "reference_stations.csv"), index=False)
    log(f"[p4-0] queue: {len(Q)} station-days ({Q.role.value_counts().to_dict()}), "
        f"{Q.sta.nunique()} stations x {PD.date.nunique()} days")
    with open(os.path.join(a.out, "u_select.log"), "w") as f:
        f.write("\n".join(lines) + "\n")


def REF_META(log):
    """CN.PGC + PB.B011 - the reference pair; band/loc/sensitivity from IRIS StationXML."""
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime, read_inventory
    os.makedirs(XMLDIR, exist_ok=True)
    out = []
    for net, sta in U.REF_STATIONS:
        p = os.path.join(XMLDIR, f"REF_{net}.{sta}.xml")
        if os.path.exists(p):
            inv = read_inventory(p)
        else:
            inv = Client("IRIS", timeout=300).get_stations(
                network=net, station=sta, level="response",
                starttime=UTCDateTime(f"{U.PROBE_Y0}-01-01"),
                endtime=UTCDateTime(f"{U.PROBE_Y1}-12-31T23:59:59"))
            inv.write(p, format="STATIONXML")
        rows = []
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
                                     az=float(ch.azimuth or 0.0), dip=float(ch.dip or 0.0),
                                     sensitivity=sens, input_units=iu,
                                     lat=float(ch.latitude), lon=float(ch.longitude)))
        RR = pd.DataFrame(rows)
        best = None
        for (loc, band), g in RR.groupby(["loc", "band"]):
            slots = {c[-1] for c in g.chan}
            if not ({"Z"} <= slots and ({"N", "1"} & slots) and ({"E", "2"} & slots)):
                continue
            key = (U.BAND_PREF.index(band), LOCPREF.index(loc) if loc in LOCPREF else 99)
            if best is None or key < best[0]:
                best = (key, loc, band, g)
        _, loc, band, g = best
        latest = g.sort_values("chan").groupby("chan").tail(1)
        zrow = latest[latest.chan.str[-1] == "Z"]
        out.append(dict(net=net, sta=sta, loc=loc, band=band, fs=float(g.fs.max()),
                        chain="frozen" if abs(g.fs.max() / U.FS_CACHE - round(g.fs.max() / U.FS_CACHE)) < 1e-9 else "rational",
                        sensitivity_med=float(np.nanmedian(latest.sensitivity.values)),
                        input_units=latest.input_units.mode().iloc[0],
                        zsign=-1.0 if (len(zrow) and zrow.dip.iloc[0] > 0) else 1.0,
                        chan_az_dip=";".join(f"{r.chan}:{r.az:.1f}/{r.dip:.1f}" for r in latest.itertuples()),
                        lat=float(g.lat.iloc[0]), lon=float(g.lon.iloc[0])))
        log(f"  [ref] {net}.{sta}: band {band} loc'{loc}' fs={g.fs.max():g} "
            f"sens={out[-1]['sensitivity_med']:.4g} {out[-1]['input_units']}")
    return out


if __name__ == "__main__":
    main()
