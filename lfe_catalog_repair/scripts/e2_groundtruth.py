#!/usr/bin/env python
"""E2 — ground-truth preparation: station coordinates, the Bostock/Lin event tables with
PREDICTED arrival times at each station, and the catalogue-time-convention calibration.

Pieces
  1. `results/e2_stations.csv`  — net/sta/lat/lon/elev from FDSN (cached).
  2. `results/e2_stack_calib.csv` — per (Bostock family, station) S-P differential time read from
     the SeisBench `LFEStacksCascadiaBostock2015` stacks (the stacks are aligned with the P arrival
     at sample 800, 20 Hz), next to the Crosson-1D model prediction.  This is the pick-time offset
     calibration asked for in the brief and it also validates the 1-D model.
  3. `results/e2_gt_bostock.parquet` / `e2_gt_lin.parquet` — event tables restricted to the pilot
     windows with (a) the catalogue time, (b) hypocentre, (c) predicted P and S arrival at every
     pilot station from the Crosson 1-D model.
  4. The Bostock catalogue-time convention.  REMOVED 2026-09-11 (E3 step 5): the former estimator
     C = median(Lin OT - Bostock time) over co-detected events is chance-dominated (16,026 of
     18,114 Bostock events have a Lin event within 20 s BY CHANCE at 500 Lin events/day) and
     reported C = +0.086 s.  The time base is now measured per FAMILY from the residual histogram
     of our own S picks by `e3_bostock_timebase.py` -> `e3_bostock_family_offsets.csv`.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa: E402  (Crosson 1976 1-D model + exact flat-layer travel times)

RES = "lfe_catalog_repair/results"
EXT = f"{RES}/external"
STACK_MD = f"{EXT}/other/seisbench_lfestacks_cascadia/metadata.csv"
BOSTOCK = f"{EXT}/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
LIN = f"{EXT}/lin_catalog/EQloc_001_0.1_3_S.csv"

PILOT_STATIONS = ["CN.PGC", "CN.SNB", "CN.VGZ", "CN.LZB", "CN.PFB", "CN.YOUB", "CN.NLLB",
                  "PB.B011", "PB.B926", "UW.LRIV", "CN.SYMB", "CN.CLRS", "CN.PTRF",
                  "UW.WISH", "PB.B943", "UW.YACT", "UW.OCP", "UW.STW", "CN.GOBB", "PB.B926"]
WINDOWS = {"pgc2010": ("2010-05-01", "2010-11-01"),
           "vi2010": ("2010-08-10", "2010-09-14"),
           "vi2017": ("2017-02-20", "2017-04-21"),
           "offgt": ("2017-01-01", "2019-01-01")}
R_EARTH_KM_DEG = 111.195


def stations() -> pd.DataFrame:
    out = f"{RES}/e2_stations.csv"
    if os.path.exists(out):
        return pd.read_csv(out)
    from obspy.clients.fdsn import Client
    rows = []
    for prov in ("IRIS",):
        c = Client(prov, timeout=120)
        for s in sorted(set(PILOT_STATIONS)):
            net, sta = s.split(".")
            try:
                inv = c.get_stations(network=net, station=sta, level="station")
            except Exception as e:
                print(f"[warn] {s}: {e}")
                continue
            st = inv[0][0]
            rows.append(dict(net=net, sta=sta, lat=st.latitude, lon=st.longitude,
                             elev_m=st.elevation))
    df = pd.DataFrame(rows).drop_duplicates(["net", "sta"])
    df.to_csv(out, index=False)
    return df


def epi_km(lat1, lon1, lat2, lon2):
    lat0 = np.radians((np.asarray(lat1) + np.asarray(lat2)) / 2.0)
    dx = (np.asarray(lon2) - np.asarray(lon1)) * R_EARTH_KM_DEG * np.cos(lat0)
    dy = (np.asarray(lat2) - np.asarray(lat1)) * R_EARTH_KM_DEG
    return np.hypot(dx, dy)


_TT_CACHE = f"{RES}/e2_ttable_crosson.npz"
_TT = {}


def _tt_table():
    """Cached (x, z) travel-time table for the Crosson 1-D model; 0.59 ms per exact call is far too
    slow for 10^6 events x 20 stations, so everything is bilinearly interpolated off a 0.5 x 0.25 km
    grid (interpolation error < 1 ms, checked against the exact solver)."""
    global _TT
    if _TT:
        return _TT
    if os.path.exists(_TT_CACHE):
        z = np.load(_TT_CACHE)
        _TT = {k: z[k] for k in z.files}
        return _TT
    x = np.arange(0.0, 400.01, 0.5)
    d = np.arange(0.0, 80.01, 0.25)
    P = np.empty((len(x), len(d)))
    S = np.empty((len(x), len(d)))
    for i, xi in enumerate(x):
        for j, zj in enumerate(d):
            P[i, j] = v1d.tt(float(xi), float(zj), "P")
            S[i, j] = v1d.tt(float(xi), float(zj), "S")
    _TT = dict(x=x, z=d, P=P, S=S)
    np.savez_compressed(_TT_CACHE, **_TT)
    return _TT


def tt_vec(x_km, z_km, phase):
    t = _tt_table()
    x = np.clip(np.atleast_1d(np.asarray(x_km, float)), t["x"][0], t["x"][-1])
    z = np.clip(np.atleast_1d(np.asarray(z_km, float)), t["z"][0], t["z"][-1])
    x, z = np.broadcast_arrays(x, z)
    G = t[phase]
    dx, dz = t["x"][1] - t["x"][0], t["z"][1] - t["z"][0]
    ix = np.clip(((x - t["x"][0]) / dx).astype(int), 0, len(t["x"]) - 2)
    iz = np.clip(((z - t["z"][0]) / dz).astype(int), 0, len(t["z"]) - 2)
    fx = (x - t["x"][ix]) / dx
    fz = (z - t["z"][iz]) / dz
    return ((1 - fx) * (1 - fz) * G[ix, iz] + fx * (1 - fz) * G[ix + 1, iz]
            + (1 - fx) * fz * G[ix, iz + 1] + fx * fz * G[ix + 1, iz + 1])


def bostock_families():
    """Bostock family hypocentres + the SeisBench stack S-P picks, joined by nearest neighbour."""
    from scipy.spatial import cKDTree
    md = pd.read_csv(STACK_MD)
    fam = (md.groupby("source_family_id")[["source_latitude_deg", "source_longitude_deg",
                                           "source_depth_km"]].first().reset_index())
    bo = pd.read_csv(BOSTOCK)
    key = bo.groupby(["lat", "lon", "dep"]).size().reset_index(name="n_events")
    cs = np.cos(np.radians(48.6))
    A = np.c_[fam.source_latitude_deg * R_EARTH_KM_DEG,
              fam.source_longitude_deg * R_EARTH_KM_DEG * cs, fam.source_depth_km]
    B = np.c_[key.lat * R_EARTH_KM_DEG, key.lon * R_EARTH_KM_DEG * cs, key.dep]
    d, i = cKDTree(B).query(A)
    fam["bkey"] = i
    fam["join_km"] = d
    key = key.reset_index().rename(columns={"index": "bkey"})
    fam = fam.merge(key, on="bkey", how="left")
    return md, fam


def stack_calibration(md, fam, stas):
    """S-P differential time per (family, station) from the stacks vs the 1-D model."""
    m = md[md.trace_S_arrival_sample.notna()].copy()
    # the stacks put the P arrival (observed where available, else predicted) at sample 800
    p_ref = m.trace_P_arrival_sample.fillna(m.trace_P_predicted_arrival_sample)
    m["sp_obs_s"] = (m.trace_S_arrival_sample - p_ref) / m.trace_sampling_rate_hz
    m = m[p_ref.notna()]
    m = m.merge(fam[["source_family_id", "bkey", "join_km", "n_events"]], on="source_family_id")
    m = m.merge(stas.rename(columns={"sta": "station_code", "lat": "sta_lat", "lon": "sta_lon"}),
                on="station_code", how="left")
    # station coordinates carried in the stack metadata itself (all 21 stations)
    m["sta_lat"] = m.sta_lat.fillna(m.station_latitude_deg)
    m["sta_lon"] = m.sta_lon.fillna(m.station_longitude_deg)
    m["epi_km"] = epi_km(m.source_latitude_deg, m.source_longitude_deg, m.sta_lat, m.sta_lon)
    m["sp_mod_s"] = (tt_vec(m.epi_km, m.source_depth_km, "S")
                     - tt_vec(m.epi_km, m.source_depth_km, "P"))
    m["sp_resid_s"] = m.sp_obs_s - m.sp_mod_s
    cols = ["source_family_id", "station_code", "station_network_code", "bkey", "join_km",
            "n_events", "source_latitude_deg", "source_longitude_deg", "source_depth_km",
            "epi_km", "sp_obs_s", "sp_mod_s", "sp_resid_s", "trace_S_snr_db", "trace_stack_count"]
    return m[cols].sort_values(["source_family_id", "station_code"])


def predicted_arrivals(ev, stas, latc="lat", lonc="lon", depc="dep"):
    """Wide table of predicted P and S arrival OFFSETS (s after the catalogue time) per station."""
    out = {}
    for r in stas.itertuples():
        x = epi_km(ev[latc].values, ev[lonc].values, r.lat, r.lon)
        z = np.abs(ev[depc].values)
        out[f"ttP_{r.sta}"] = tt_vec(x, z, "P")
        out[f"ttS_{r.sta}"] = tt_vec(x, z, "S")
        out[f"epi_{r.sta}"] = x
    return pd.DataFrame(out, index=ev.index)


def main():
    stas = stations()
    print(stas.to_string(index=False))
    md, fam = bostock_families()
    print(f"[fam] {len(fam)} stack families joined to Bostock keys, "
          f"join dist median {fam.join_km.median():.3f} km max {fam.join_km.max():.3f} km")

    cal = stack_calibration(md, fam, stas)
    cal.to_csv(f"{RES}/e2_stack_calib.csv", index=False)
    r = cal.sp_resid_s
    print(f"[calib] S-P from stacks vs Crosson-1D over {len(cal)} (family,station) pairs: "
          f"median {r.median():+.2f} s  MAD {1.4826*np.median(np.abs(r-r.median())):.2f} s  "
          f"p10 {r.quantile(.1):+.2f}  p90 {r.quantile(.9):+.2f}")
    pgc = cal[cal.station_code == "PGC"]
    print(f"[calib] PGC only, {len(pgc)} families: S-P obs median {pgc.sp_obs_s.median():.2f} s, "
          f"model median {pgc.sp_mod_s.median():.2f} s, resid median {pgc.sp_resid_s.median():+.2f} s")

    # ---- Bostock events, whole pilot span
    bo = pd.read_csv(BOSTOCK)
    t = pd.to_datetime(dict(year=bo.year, month=bo.month, day=bo.day, hour=bo.hour,
                            minute=bo["min"], second=bo.sec.astype(int))) \
        + pd.to_timedelta((bo.sec % 1.0), unit="s")
    bo["time"] = t.astype("datetime64[ns]")
    cs = np.cos(np.radians(48.6))
    from scipy.spatial import cKDTree
    key = bo.groupby(["lat", "lon", "dep"]).size().reset_index(name="n").reset_index()
    B = np.c_[key.lat * R_EARTH_KM_DEG, key.lon * R_EARTH_KM_DEG * cs, key.dep]
    bo = bo.merge(key[["index", "lat", "lon", "dep"]].rename(columns={"index": "bkey"}),
                  on=["lat", "lon", "dep"], how="left")
    fam_by_bkey = fam.set_index("bkey")
    bo["fam_id"] = bo.bkey.map(fam_by_bkey.source_family_id)
    s0, s1 = WINDOWS["pgc2010"]
    sel = bo[(bo.time >= s0) & (bo.time < s1)].copy()
    pa = predicted_arrivals(sel, stas)
    sel = pd.concat([sel, pa], axis=1)
    sel.to_parquet(f"{RES}/e2_gt_bostock.parquet", index=False)
    print(f"[bostock] {len(sel)} events in {s0}..{s1}, {sel.bkey.nunique()} families")

    # ---- Lin events
    lin = pd.read_csv(LIN)
    lin["OT"] = pd.to_datetime(lin.OT, format="mixed").dt.tz_localize(None).astype("datetime64[ns]")
    lin["starttime"] = pd.to_datetime(lin.starttime, format="mixed").dt.tz_localize(None).astype("datetime64[ns]")
    keep = []
    for name in ("pgc2010", "vi2010", "vi2017"):
        a, b = WINDOWS[name]
        k = lin[(lin.OT >= a) & (lin.OT < b)].copy()
        k["window"] = name
        keep.append(k)
    lw = pd.concat(keep, ignore_index=True).drop_duplicates(["starttime", "OT", "lon", "lat"])
    lw["dep"] = np.abs(lw.depth)
    pa = predicted_arrivals(lw, stas)
    lw = pd.concat([lw, pa], axis=1)
    lw.to_parquet(f"{RES}/e2_gt_lin.parquet", index=False)
    print(f"[lin] {len(lw)} events kept; per window: "
          f"{lw.groupby('window').size().to_dict()}")

    # ---- catalogue-time convention: RETIRED (E3 step 5, 2026-09-11)
    # The original estimator here was C = median(Lin OT - Bostock t) over Bostock events with a Lin
    # event within +-20 s.  It is CHANCE-DOMINATED and must not be used: at ~500 Lin events/day,
    # 16,026 of the 18,114 Bostock events in this window have such a neighbour BY CHANCE, so the
    # median is pulled to ~0 whatever the truth is (it reported C = +0.086 s).  The downstream
    # scripts stopped using its value on 2026-09-10, but the file it wrote,
    # `e2_bostock_lin_dt.csv`, was still readable as a calibration -- Merlin F7 asked for it to be
    # fixed or deleted, so it is no longer written.
    #
    # The Bostock/Slow-Earthquake-DB time base is instead measured EMPIRICALLY and PER FAMILY by
    # `e3_bostock_timebase.py` -> `results/e3_bostock_family_offsets.csv` (equivalent origin time =
    # DB time + offset[family]).  Measured at PGC over 2010-05..11: median +11.75 s, sd 4.14 s,
    # range +2.5..+20.5 s over 90 of 100 families, cross-validated against 3+ other stations in an
    # independent window at r = 0.96 (PGC-minus-others sd 1.19 s).  A SINGLE constant (E2's
    # +12.46 s) is wrong by up to 9 s for individual families.
    # The E2 output `e2_bostock_lin_dt.csv` is NOT deleted or renamed (E3 rule: never touch an E2
    # output); a sibling warning file is written next to it instead so it cannot be picked up as a
    # calibration by mistake.
    _stale = f"{RES}/e2_bostock_lin_dt.csv"
    _warn = f"{RES}/e2_bostock_lin_dt.WARNING_chance_dominated.txt"
    with open(_warn, "w") as _f:
        _f.write(
            "e2_bostock_lin_dt.csv IS NOT A CALIBRATION -- DO NOT USE IT.\n\n"
            "It holds (Lin OT - Bostock time) for Bostock events with a Lin event within 20 s.\n"
            "At ~500 Lin events/day, 16,026 of the 18,114 Bostock events in the pgc2010 window\n"
            "have such a neighbour BY CHANCE, so the median is pulled to ~0 whatever the truth is\n"
            "(it reported C = +0.086 s).  Merlin audit 2026-09-11, flaw F7; E2 bug 3.\n\n"
            "The Bostock / Slow-Earthquake-DB time base is measured empirically and PER FAMILY in\n"
            "results/e3_bostock_family_offsets.csv (E3 step 5, e3_bostock_timebase.py):\n"
            "  equivalent origin time = DB time + offset[family]\n"
            "  90 of 100 families accepted; median +11.75 s, s.d. 4.14 s, range +2.50..+20.50 s;\n"
            "  cross-validated at r = 0.96 against >=3 other stations in an independent window\n"
            "  (PGC-minus-others s.d. 1.19 s).  A single constant is wrong by up to 9 s.\n"
            "The convention itself is still unidentified and needs the paper / DB maintainers.\n")
    print(f"[convention] chance-dominated C estimator REMOVED from this script; "
          f"wrote {_warn} next to the stale {os.path.basename(_stale)}")



if __name__ == "__main__":
    main()
