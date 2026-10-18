#!/usr/bin/env python
"""Step 1c -- delivery probe for the 2005-09 POLARIS bed.  METADATA IS NOT DELIVERY.

Tries every (network, band) for each candidate station on four days (three inside the 2005-09 ETS
window, one quiet day) at BOTH providers, and writes results/step1c_stations.csv with the winning
(net, band, provider) plus the geometry (epicentral / hypocentral distance to the 130 Bostock
family hypocentres, families within 25 km).  Family COORDINATES only -- no event time is read.
"""
from __future__ import annotations
import json, os, sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np, pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

DAYS = ["2005-09-05", "2005-09-10", "2005-09-15", "2005-07-05"]
NETS = ["C8", "PO", "CN", "PB", "UW"]
BANDS = ["HH", "BH", "EH"]
CAND = ["SOKB", "GLBC", "SHVB", "TWKB", "MGCB", "TSJB", "JRBC", "KELB", "KHVB", "TWBB",
        "TWGB", "SILB", "SSIB", "LCBC", "GOWB",                      # POLARIS dense
        "PGC", "VGZ", "SNB", "LZB", "PFB", "YOUB", "BPCB", "GOBB",   # CN permanents
        "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011",   # PB boreholes
        "SQM", "STW", "OBC"]                                          # UW


def probe(sta):
    provs = {"IRIS": "IRIS", "NRCAN": "https://earthquakescanada.nrcan.gc.ca"}
    best = None
    lat = lon = np.nan
    for pn, url in provs.items():
        try:
            c = Client(url, timeout=120)
        except Exception:
            continue
        for net in NETS:
            for band in BANDS:
                ok, fs = 0, np.nan
                for d in DAYS:
                    try:
                        st = c.get_waveforms(net, sta, "*", band + "?", UTCDateTime(d),
                                             UTCDateTime(d) + 1800)
                        if len(st) >= 3:
                            ok += 1
                            fs = float(st[0].stats.sampling_rate)
                    except Exception:
                        pass
                if ok:
                    try:
                        inv = c.get_stations(network=net, station=sta, level="station",
                                             starttime=UTCDateTime(DAYS[0]))
                        lat = float(inv[0][0].latitude); lon = float(inv[0][0].longitude)
                    except Exception:
                        pass
                    cand = dict(sta=sta, net=net, band=band, prov=pn, n_days=ok,
                                n_probe=len(DAYS), fs=fs, lat=lat, lon=lon)
                    if best is None or ok > best["n_days"]:
                        best = cand
                    if ok == len(DAYS):
                        return best
    return best or dict(sta=sta, net="", band="", prov="", n_days=0, n_probe=len(DAYS),
                        fs=np.nan, lat=np.nan, lon=np.nan)


def main():
    F = S.families()
    with ThreadPoolExecutor(max_workers=6) as ex:
        rows = list(ex.map(probe, CAND))
    D = pd.DataFrame(rows)
    ST = W.station_table()
    for i, r in D.iterrows():
        la, lo = r.lat, r.lon
        if not np.isfinite(la) and r.sta in ST.index:
            la, lo = float(ST.lat[r.sta]), float(ST.lon[r.sta])
        D.loc[i, "lat"], D.loc[i, "lon"] = la, lo
    ok = np.isfinite(D.lat)
    D["epi_med"] = np.nan; D["hyp_med"] = np.nan; D["n_fam_epi25"] = 0
    for i, r in D[ok].iterrows():
        e = S.epi_km(F.lat.values, F.lon.values, r.lat, r.lon)
        h = S.hyp_km(F.lat.values, F.lon.values, F.dep.values, r.lat, r.lon)
        D.loc[i, "epi_med"] = round(float(np.median(e)), 1)
        D.loc[i, "hyp_med"] = round(float(np.median(h)), 1)
        D.loc[i, "n_fam_epi25"] = int((e <= 25).sum())
    D = D.sort_values(["epi_med"]).reset_index(drop=True)
    D.to_csv(f"{W.RES}/step1c_probe.csv", index=False)
    G = D[(D.n_days > 0) & np.isfinite(D.epi_med) & (D.epi_med <= 80)]
    G[["net", "sta", "lat", "lon", "band", "prov", "n_days", "epi_med", "hyp_med",
       "n_fam_epi25"]].to_csv(f"{W.RES}/step1c_stations.csv", index=False)
    print(D.to_string(index=False))
    print(f"\nDELIVERING within 80 km: {len(G)} -> {W.RES}/step1c_stations.csv")


if __name__ == "__main__":
    main()
