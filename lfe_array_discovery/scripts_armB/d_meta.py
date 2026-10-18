#!/usr/bin/env python
"""Tranche 1 - channel-level StationXML for every network in the box, and the site table.

For each of XU XG ZK Z3 YW TA PB UW CN CC US: FDSN `level=response` over the box, 2006-01-01..
2012-12-31 (permanents to 2026-12-31 so tranche 7 can choose a band per epoch) ->
  results/stationxml/<NET>.xml
  results/channels_all.csv          every channel epoch returned
  results/channels_2006_2012.csv    net, sta, loc, chan, fs, azimuth, dip, sensor, sensitivity,
                                    start, end (+ lat/lon/elev/band) for epochs overlapping the
                                    array era
  results/station_flags.csv         ehz_only | rational_rate | short_period | zsign per station
  results/sites_array_era.csv       site_id, members, lat, lon  (1 km single-linkage)
Writes results/markers/DL_META_DONE when finished.

The dip sign contract is u_select.py's: zsign = -1 when the vertical is declared at dip +90
(positive DOWN), +1 otherwise; the cache stores raw Z/N/E and the reader applies zsign.
"""
from __future__ import annotations
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402


def fetch_net(cli, net, t0, t1, force=False):
    from obspy import UTCDateTime, read_inventory
    p = os.path.join(D.XMLDIR, f"{net}.xml")
    if os.path.exists(p) and not force:
        D.log(f"  [meta] {net}: cached {p}")
        return read_inventory(p)
    inv = cli.get_stations(network=net, station="*", location="*", channel="*",
                           level="response",
                           minlatitude=D.LAT0, maxlatitude=D.LAT1,
                           minlongitude=D.LON0, maxlongitude=D.LON1,
                           starttime=UTCDateTime(t0), endtime=UTCDateTime(t1))
    inv.write(p, format="STATIONXML")
    D.log(f"  [meta] {net}: {sum(len(s) for s in inv)} station epochs -> {p}")
    return inv


def inv_rows(inv, net_default=""):
    rows = []
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
                    band=ch.code[:2], fs=float(ch.sample_rate or 0.0),
                    azimuth=float(ch.azimuth) if ch.azimuth is not None else np.nan,
                    dip=float(ch.dip) if ch.dip is not None else np.nan,
                    lat=float(ch.latitude), lon=float(ch.longitude),
                    elev=float(ch.elevation), depth_m=float(ch.depth or 0.0),
                    sensor=str(getattr(ch.sensor, "description", "") or "")[:90],
                    sensitivity=sens, input_units=iu,
                    start=str(ch.start_date.date) if ch.start_date else "1990-01-01",
                    end=str(ch.end_date.date) if ch.end_date else "2100-01-01"))
    return rows


def is_short_period(band, sensor):
    s = (sensor or "").upper()
    return bool(band in D.SHORT_PERIOD_BANDS or any(k in s for k in D.SHORT_PERIOD_PAT))


def build_flags(C):
    """One row per (net, sta), computed over the ARRAY-ERA epochs only."""
    out = []
    for (net, sta), g in C.groupby(["net", "sta"]):
        gv = g[g.chan.str[-1].isin(list("ZNE12")) & g.band.isin(D.BAND_PREF)]
        best, epochs3c = None, []
        for (loc, band), gg in gv.groupby(["loc", "band"]):
            slots = {c[-1] for c in gg.chan}
            three = ("Z" in slots and ({"N", "1"} & slots) and ({"E", "2"} & slots)
                     and len({"N", "1"} & slots) == 1 and len({"E", "2"} & slots) == 1)
            if not three or float(gg.fs.max()) < D.FS_MIN:
                continue
            key = (D.BAND_PREF.index(band),
                   D.LOCPREF.index(loc) if loc in D.LOCPREF else 99, loc)
            epochs3c.append((str(gg.start.min()), str(gg.end.max())))
            if best is None or key < best[0]:
                best = (key, loc, band, float(gg.fs.max()), gg)
        # sensor / dip must come from the SEISMOMETER channels we would actually use, not from
        # the station's strainmeter or state-of-health channels (PB EH? vs GLADWIN strainmeter).
        gs = gv if len(gv) else g
        vert = gs[gs.chan.str[-1] == "Z"]
        zdip = float(vert.dip.dropna().iloc[0]) if len(vert.dip.dropna()) else -90.0
        sensor = (gs.sensor.mode().iloc[0] if len(gs.sensor.mode()) else "")
        if best is None:
            zb = gs[gs.chan.str[-1] == "Z"]
            zloc = str(zb.loc[:, "loc"].iloc[0]) if len(zb) else ""
            zband = str(zb.band.iloc[0]) if len(zb) else str(g.band.mode().iloc[0])
            zfs = float(zb.fs.max()) if len(zb) else float(g.fs.max())
            out.append(dict(net=net, sta=sta, loc=zloc, band=zband, fs=zfs,
                            n_comp=int(gs[gs.chan.str[-1].isin(list("ZNE12"))].chan.nunique()),
                            three_c=False, ehz_only=True, rational_rate=False,
                            short_period=is_short_period(zband, sensor),
                            zsign=-1.0 if zdip > 0 else 1.0, zdip=zdip,
                            chain=D.chain_for(zfs) if zfs else "",
                            lat=float(g.lat.iloc[0]), lon=float(g.lon.iloc[0]),
                            elev=float(g.elev.iloc[0]), sensor=sensor,
                            sensitivity=float(np.nanmedian(gs.sensitivity.values)) if len(gs) else np.nan,
                            start=str(gs.start.min()) if len(gs) else str(g.start.min()),
                            end=str(gs.end.max()) if len(gs) else str(g.end.max()),
                            chans=",".join(sorted(gs.chan.unique())) if len(gs) else ",".join(sorted(g.chan.unique()))))
            continue
        _, loc, band, fs, gg = best
        out.append(dict(net=net, sta=sta, loc=loc, band=band, fs=fs,
                        n_comp=int(gg.chan.nunique()), three_c=True, ehz_only=False,
                        rational_rate=(D.chain_for(fs) == "rational"),
                        short_period=is_short_period(band, sensor),
                        zsign=-1.0 if zdip > 0 else 1.0, zdip=zdip, chain=D.chain_for(fs),
                        lat=float(gg.lat.iloc[0]), lon=float(gg.lon.iloc[0]),
                        elev=float(gg.elev.iloc[0]), sensor=sensor,
                        sensitivity=float(np.nanmedian(gg.sensitivity.values)),
                        sensor_band=band,
                        start=min(e[0] for e in epochs3c), end=max(e[1] for e in epochs3c),
                        chans=",".join(sorted(gg.chan.unique()))))
    return pd.DataFrame(out).sort_values(["net", "sta"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--prov", default="IRIS")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    from obspy.clients.fdsn import Client
    cli = Client(a.prov, timeout=300)
    rows = []
    for net in D.NETS:
        t1 = D.PERM_T1 if net in D.PERM_NETS else D.META_T1
        try:
            inv = fetch_net(cli, net, D.META_T0, t1, a.force)
        except Exception as e:
            D.log(f"  [meta] {net}: FAILED {type(e).__name__}: {str(e)[:150]}")
            continue
        rows += inv_rows(inv)
    A = pd.DataFrame(rows)
    A = A[D.in_box(A.lat.values, A.lon.values)].drop_duplicates()
    A = A.sort_values(["net", "sta", "loc", "chan", "start"]).reset_index(drop=True)
    A.to_csv(os.path.join(D.RES, "channels_all.csv"), index=False)
    C = A[(A.start <= D.META_T1) & (A.end >= D.META_T0)].reset_index(drop=True)
    cols = ["net", "sta", "loc", "chan", "fs", "azimuth", "dip", "sensor", "sensitivity",
            "start", "end", "band", "lat", "lon", "elev", "depth_m", "input_units"]
    C[cols].to_csv(os.path.join(D.RES, "channels_2006_2012.csv"), index=False)
    F = build_flags(C)
    F.to_csv(os.path.join(D.RES, "station_flags.csv"), index=False)
    # ---- data-driven site table (1 km single-linkage) over every array-era station in the box
    S = F[["net", "sta", "lat", "lon"]].copy()
    S["key"] = S.net + "." + S.sta
    sm = D.cluster_sites(S, key="key")
    S["site"] = [sm.get(k, k.split(".")[-1]) for k in S.key]
    T = (S.groupby("site")
           .agg(members=("key", lambda v: ";".join(sorted(v))), n=("key", "size"),
                lat=("lat", "mean"), lon=("lon", "mean"))
           .reset_index().sort_values(["lat", "lon"]).reset_index(drop=True))
    T.insert(0, "site_id", [f"S{i:03d}" for i in range(len(T))])
    T = T.rename(columns={"site": "site_name"})
    T.to_csv(os.path.join(D.RES, "sites_array_era.csv"), index=False)
    D.log(f"[meta] {len(A)} channel epochs, {len(C)} in 2006-2012, "
          f"{len(F)} stations, {len(T)} sites (1 km)")
    D.log("[meta] per network (stations / 3-C / EHZ-only / rational / short-period):")
    for net, g in F.groupby("net"):
        D.log(f"    {net:3s} {len(g):4d}  3C={int(g.three_c.sum()):4d}  "
              f"ehz_only={int(g.ehz_only.sum()):3d}  rational={int(g.rational_rate.sum()):3d}  "
              f"sp={int(g.short_period.sum()):3d}  fs={sorted(set(g.fs))}")
    D.marker("DL_META_DONE", f"{len(F)} stations, {len(T)} sites")


if __name__ == "__main__":
    main()
