#!/usr/bin/env python
"""rf00_events: teleseismic event list (USGS ComCat, M>=5.8, 1993-01-01..now) and per-station P arrival table
(30-90 deg, iasp91 P time + ray parameter), restricted to the station's BH/HH channel epochs.

Outputs: data/rf_teleseis/events_M58_1993_2026.csv, data/rf_teleseis/arrivals_<NET>.<STA>.csv,
         data/rf_teleseis/stations.csv
"""
import io, os, sys, time, urllib.request
import numpy as np, pandas as pd
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/data/rf_teleseis"
STATIONS = [  # core (Gosselin 2020) first, then extension
    ("CN", "PGC"), ("CN", "SNB"), ("CN", "VGZ"), ("CN", "LZB"), ("CN", "PFB"),
    ("CN", "NLLB"), ("CN", "YOUB"), ("CN", "MGB"), ("CN", "GOBB"),
    ("UW", "GNW"), ("UW", "SP2"), ("UW", "LON"), ("UW", "OFR"), ("UW", "WISH"), ("US", "NLWA"),
]


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # noqa
            print("retry", k, e, flush=True); time.sleep(5 * (k + 1))
    raise RuntimeError(url)


def events():
    f = f"{OUT}/events_M58_1993_2026.csv"
    if os.path.exists(f):
        return pd.read_csv(f)
    rows = []
    for y in range(1993, 2027):
        u = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minmagnitude=5.8&orderby=time-asc"
             f"&starttime={y}-01-01&endtime={y + 1}-01-01")
        t = get(u)
        d = pd.read_csv(io.StringIO(t))
        print(y, len(d), flush=True)
        rows.append(d)
    E = pd.concat(rows, ignore_index=True)
    E = E[["time", "latitude", "longitude", "depth", "mag", "magType", "id"]].rename(columns={"id": "evid"})
    E.to_csv(f, index=False)
    return E


def channel_epochs(net, sta):
    u = (f"https://service.iris.edu/fdsnws/station/1/query?net={net}&sta={sta}&cha=BH?,HH?&level=channel"
         "&format=text&nodata=404")
    d = pd.read_csv(io.StringIO(get(u)), sep="|")
    d.columns = [c.strip().lstrip("#").strip() for c in d.columns]
    d["Location"] = d["Location"].fillna("").astype(str).str.strip()
    d["EndTime"] = d["EndTime"].fillna("2100-01-01T00:00:00")
    return d


def main():
    E = events()
    E["t"] = pd.to_datetime(E["time"], utc=True).dt.tz_localize(None).astype("datetime64[ns]")
    model = TauPyModel("iasp91")
    # P time / ray parameter tables on a (distance, depth) grid
    dist = np.arange(29.5, 90.6, 0.25); dep = np.arange(0, 701, 25.0)
    tP = np.full((len(dist), len(dep)), np.nan); pP = np.full_like(tP, np.nan)
    tab = f"{OUT}/ptable_iasp91.npz"
    if os.path.exists(tab):
        z = np.load(tab); tP, pP = z["tP"], z["pP"]
    else:
        for i, D in enumerate(dist):
            for j, h in enumerate(dep):
                a = model.get_travel_times(source_depth_in_km=h, distance_in_degree=D, phase_list=["P"])
                if a:
                    tP[i, j] = a[0].time; pP[i, j] = a[0].ray_param_sec_degree / 111.19  # s/km
        np.savez(tab, dist=dist, dep=dep, tP=tP, pP=pP)
    from scipy.interpolate import RegularGridInterpolator as RGI
    fT = RGI((dist, dep), tP); fp = RGI((dist, dep), pP)
    srows = []
    for net, sta in STATIONS:
        C = channel_epochs(net, sta)
        lat, lon = float(C.Latitude.iloc[0]), float(C.Longitude.iloc[0])
        srows.append(dict(net=net, sta=sta, lat=lat, lon=lon, elev=float(C.Elevation.iloc[0])))
        D = locations2degrees(lat, lon, E.latitude.values, E.longitude.values)
        m = (D >= 30) & (D <= 90) & E.depth.between(0, 700).values
        A = E[m].copy(); A["gcarc"] = D[m]
        A["baz"] = [gps2dist_azimuth(lat, lon, a, b)[1] for a, b in zip(A.latitude, A.longitude)]  # [1] = station->event azimuth = back-azimuth (fixed 2026-10-13; [2] was wrong)
        pts = np.c_[A.gcarc.values, A.depth.clip(0, 700).values]
        A["ttP"] = fT(pts); A["slow"] = fp(pts)
        A["tP"] = A.t + pd.to_timedelta(A.ttP, unit="s")
        # channel epoch choice: prefer BH, else HH; need Z + two horizontals at P time
        ch = []
        for tp in A.tP:
            s = tp.strftime("%Y-%m-%dT%H:%M:%S")
            act = C[(C.StartTime <= s) & (C.EndTime > s)]
            pick = ""
            for band in ("BH", "HH"):
                b = act[act.Channel.str.startswith(band)]
                if b.Channel.str[-1].isin(["Z"]).any() and b.Channel.size >= 3:
                    loc = b.Location.iloc[0]
                    pick = f"{loc}|{band}?"; break
            ch.append(pick)
        A["chan"] = ch
        A = A[A.chan != ""]
        A.to_csv(f"{OUT}/arrivals_{net}.{sta}.csv", index=False)
        print(net, sta, len(A), "arrivals", A.tP.min(), A.tP.max(), flush=True)
    pd.DataFrame(srows).to_csv(f"{OUT}/stations.csv", index=False)


if __name__ == "__main__":
    main()
