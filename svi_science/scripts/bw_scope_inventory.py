#!/usr/bin/env python
"""bw_scope_inventory -- SVI backward-extension scope (2003-01..2013-12), step 1: channel inventory.

Queries FDSN station metadata (level=channel) at NRCan/CNDC and EarthScope for every station inside a box
around the 115 SVI families, keeps stations whose epicentral distance to the NEAREST family is <= 60 km
(plus every forward template station regardless), and writes one row per (centre, net, sta, loc, cha, epoch)
overlapping 2003-01-01..2014-01-01.  Metadata only; no waveforms.

Out: svi_science/results/backward_scope/inventory_channels.csv, inventory_stations.csv
"""
import io, sys, time
import numpy as np, pandas as pd, requests

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
CENTRES = {"NRCAN": "https://www.earthquakescanada.nrcan.gc.ca/fdsnws",
           "EARTHSCOPE": "https://service.earthscope.org/fdsnws"}
NETS_KEEP = {"CN", "C8", "PO", "PB", "UW", "TA", "US", "XG", "4P", "XU", "X6", "YW"}   # SY synthetic, GD GPS, GS/NP strong motion
BANDS = ("HH", "BH", "EH", "SH", "HN"[:0] or "EL")   # seismic velocity bands we could template on
TEMPLATE_STAS = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
                 "B009", "B010", "B011", "B926"]
MAXD = 60.0


def hav(lat1, lon1, lat2, lon2):
    r = np.pi / 180
    a = np.sin((lat2 - lat1) * r / 2) ** 2 + np.cos(lat1 * r) * np.cos(lat2 * r) * np.sin((lon2 - lon1) * r / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def get(url, tries=4):
    for k in range(tries):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code in (200, 204):
                return r
        except Exception as e:
            print("err", e, file=sys.stderr)
        time.sleep(3 * (k + 1))
    return None


def main():
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    rows = []
    for c, base in CENTRES.items():
        url = (f"{base}/station/1/query?minlat=47.6&maxlat=49.6&minlon=-125.4&maxlon=-122.2"
               f"&starttime=2003-01-01&endtime=2014-01-01&level=channel&format=text")
        r = get(url)
        if r is None or r.status_code != 200:
            print("no metadata from", c); continue
        d = pd.read_csv(io.StringIO(r.text), sep="|", dtype=str)
        d.columns = [x.strip().lstrip("#").strip() for x in d.columns]
        d = d.apply(lambda s: s.str.strip())
        d["centre"] = c
        rows.append(d)
    D = pd.concat(rows, ignore_index=True)
    D = D.rename(columns={"Network": "net", "Station": "sta", "Location": "loc", "Channel": "cha",
                          "Latitude": "lat", "Longitude": "lon", "SampleRate": "fs", "StartTime": "start",
                          "EndTime": "end", "SensorDescription": "sensor", "Scale": "scale"})
    D["loc"] = D["loc"].fillna("")
    D = D[D.net.isin(NETS_KEEP)]
    D = D[D.cha.str[:2].isin(["HH", "BH", "EH", "SH", "EL", "HN", "EN"])]
    D["lat"] = D.lat.astype(float); D["lon"] = D.lon.astype(float); D["fs"] = D.fs.astype(float)
    dmin = [hav(la, lo, F.lat.values, F.lon.values).min() for la, lo in zip(D.lat, D.lon)]
    D["dist_min_km"] = np.round(dmin, 1)
    dmed = [np.median(hav(la, lo, F.lat.values, F.lon.values)) for la, lo in zip(D.lat, D.lon)]
    D["dist_med_km"] = np.round(dmed, 1)
    D = D[(D.dist_min_km <= MAXD) | D.sta.isin(TEMPLATE_STAS)]
    D["end"] = D["end"].fillna("").replace("", "2599-01-01T00:00:00")
    D = D[["centre", "net", "sta", "loc", "cha", "lat", "lon", "fs", "start", "end", "sensor", "scale",
           "dist_min_km", "dist_med_km"]].sort_values(["net", "sta", "cha", "start", "centre"])
    D.to_csv(f"{OUT}/inventory_channels.csv", index=False)
    S = (D.groupby(["net", "sta"]).agg(lat=("lat", "first"), lon=("lon", "first"),
                                        dist_min_km=("dist_min_km", "min"), dist_med_km=("dist_med_km", "min"),
                                        centres=("centre", lambda s: ";".join(sorted(set(s)))),
                                        bands=("cha", lambda s: ";".join(sorted(set(x[:2] for x in s)))),
                                        first=("start", "min"), last=("end", "max")).reset_index())
    S["template_station"] = S.sta.isin(TEMPLATE_STAS)
    S.to_csv(f"{OUT}/inventory_stations.csv", index=False)
    print(S.to_string())


if __name__ == "__main__":
    main()
