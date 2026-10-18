"""T-study (revised brief) step 1: ALL seismic stations with data 2019-01-01 -> today inside
46.5-49.5 N / -125.5..-121.5 E, ANY network, from the IRIS/EarthScope station service
(level=channel; CN is served by IRIS).  NCEDC is queried too for completeness (it is out of box).
No waveforms are downloaded.

Output: results/p3d_puget_geometry/{stations_all_2019plus.csv, z_only_2019plus.csv, chan_raw.csv}
"""
import json, os
import numpy as np, pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
SCRATCH = "/tmp/claude-1000/-home-jovyan-tremorferometry/d6a923e6-f9f1-4ee1-8d37-955b6b2e0f35/scratchpad"
BOX = dict(minlatitude=46.5, maxlatitude=49.5, minlongitude=-125.5, maxlongitude=-121.5)
T0, T1 = UTCDateTime("2019-01-01"), UTCDateTime("2026-09-19")
FS_MIN = 20.0

rows = []
for name in ["IRIS", "NCEDC"]:
    cache = f"{SCRATCH}/chan_{name}.csv"
    if os.path.exists(cache):
        rows.append(pd.read_csv(cache)); print("cache", name, len(rows[-1])); continue
    try:
        cli = Client(name, timeout=300)
        inv = cli.get_stations(network="*", station="*", location="*", channel="*",
                               level="channel", starttime=T0, endtime=T1,
                               **BOX)
    except Exception as ex:
        print("FAIL", name, ex); continue
    r = []
    for n in inv:
        for s in n:
            for c in s:
                r.append(dict(src=name, net=n.code, sta=s.code, loc=c.location_code or "",
                              chan=c.code, lat=c.latitude, lon=c.longitude, elev=c.elevation,
                              depth=c.depth, fs=float(c.sample_rate or 0),
                              start=str(c.start_date)[:10],
                              end=(str(c.end_date)[:10] if c.end_date else "2026-09-19"),
                              sensor=(c.sensor.description if c.sensor else ""),
                              restricted=str(c.restricted_status), az=c.azimuth, dip=c.dip))
    d = pd.DataFrame(r); d.to_csv(cache, index=False); rows.append(d)
    print(name, len(d), "channel epochs,", d.groupby(['net','sta']).ngroups, "stations")

CH = pd.concat(rows, ignore_index=True)
CH["loc"] = CH["loc"].fillna("").astype(str)      # a cached empty location code reads back as NaN,
CH["sensor"] = CH.sensor.fillna("")
CH["restricted"] = CH.get("restricted", pd.Series("", index=CH.index)).fillna("")               # and groupby would then silently DROP those rows
CH = CH[CH.chan.str.len() == 3].copy()
CH["band"], CH["inst"], CH["comp"] = CH.chan.str[0], CH.chan.str[1], CH.chan.str[2]
CH = CH[CH.inst.isin(["H", "N", "L"]) & CH.band.isin(list("BHES"))]     # seismometer or accelerometer
CH = CH[CH.fs >= FS_MIN]
CH["end"] = CH.end.clip(upper="2026-09-19")
CH = CH[CH.end >= "2019-01-01"]
CH.to_csv(f"{OUT}/chan_raw.csv", index=False)
print("channel epochs kept:", len(CH), "stations:", CH.groupby(['net','sta']).ngroups)

# ---- roster membership ----
RO = pd.read_csv(f"{ROOT}/codex-LFE/results/target_network_from_map/stations.csv")
ROSTER = set(RO.network + "." + RO.station)
F1 = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/f1_stations.csv")
F1["key"] = F1.net + "." + F1.sta
BH_F1 = set(F1[F1.is_borehole_any].key)
# temporary/experiment network codes (first char a digit, or the X/Y/Z 2-char series)
def is_temp(net):
    return bool(net[0].isdigit() or net[0] in "XYZ")

def intervals_3c(g):
    """3-C epochs per (loc, band group): start = max over the three comps, end = min."""
    out = []
    for (loc, band), gg in g.groupby(["loc", "band"]):
        comps = {c: (min(x.start), max(x.end)) for c, x in gg.groupby("comp")}
        z = comps.get("Z")
        hor = [c for c in ("N", "E", "1", "2", "3") if c in comps]
        if z is None or len(hor) < 2:
            continue
        s = max([z[0]] + [comps[c][0] for c in hor[:2]])
        e = min([z[1]] + [comps[c][1] for c in hor[:2]])
        if s < e and e >= "2019-01-01":
            out.append([max(s, "2019-01-01"), e])
    return sorted(out)

def merge_iv(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

rows = []
for (net, sta), g in CH.groupby(["net", "sta"]):
    key = f"{net}.{sta}"
    seis = g[g.inst.isin(["H", "L"])]
    acc = g[g.inst == "N"]
    iv_seis = intervals_3c(seis) if len(seis) else []
    iv_acc = intervals_3c(acc) if len(acc) else []
    only_acc = (not iv_seis) and bool(iv_acc)
    iv = merge_iv(iv_seis if iv_seis else iv_acc)
    src = "seismometer" if iv_seis else ("accelerometer_only" if iv_acc else "")
    prim = seis if len(seis) else acc
    depth = float(np.nanmax(prim.depth.values)) if len(prim) else np.nan
    bands = set(prim.band)
    borehole = (depth >= 20) or (key in BH_F1) or net == "PB"
    if borehole:
        stype = "borehole"
    elif only_acc:
        stype = "strong_motion"
    elif bands & {"B", "H"}:
        stype = "surface_broadband"
    else:
        stype = "surface_short_period"
    zonly = (not iv) and bool(len(g[g.comp == "Z"]))
    rows.append(dict(key=key, net=net, sta=sta,
                     lat=round(float(prim.lat.median()), 5), lon=round(float(prim.lon.median()), 5),
                     elev=float(prim.elev.median()), depth_m=depth,
                     sensor_type=stype, chan_src=src, accel_only=only_acc,
                     bands="/".join(sorted(bands)), fs_max=float(prim.fs.max()),
                     sensor=str(prim.sensor.dropna().iloc[0])[:60] if prim.sensor.notna().any() else "",
                     restricted=("open" if (prim.restricted == "open").any() else
                                 ";".join(sorted(set(prim.restricted.astype(str))))),
                     nets_temporary=is_temp(net), in_roster183=key in ROSTER,
                     three_c=bool(iv), z_only=zonly, n_iv=len(iv),
                     first_3c=(iv[0][0] if iv else ""), last_3c=(iv[-1][1] if iv else ""),
                     days_3c=round(sum((np.datetime64(e) - np.datetime64(s)) / np.timedelta64(1, "D")
                                       for s, e in iv), 1),
                     intervals=json.dumps(iv)))
S = pd.DataFrame(rows).sort_values(["net", "sta"])
# detection-capable = a real 3-C seismometer set (not accelerometer-only) on a non-citizen network
S["detect_capable"] = S.three_c & (S.chan_src == "seismometer") & (S.net != "AM")
S[S.three_c].to_csv(f"{OUT}/stations_all_2019plus.csv", index=False)
S[S.z_only].to_csv(f"{OUT}/z_only_2019plus.csv", index=False)
print("3-C stations:", int(S.three_c.sum()), "| Z-only:", int(S.z_only.sum()),
      "| detection-capable:", int(S.detect_capable.sum()))
print(S[S.three_c].sensor_type.value_counts().to_dict())
print("temporary-network 3-C:", int(S[S.three_c].nets_temporary.sum()),
      "| accel-only 3-C:", int(S[S.three_c].accel_only.sum()),
      "| in roster183:", int(S[S.three_c].in_roster183.sum()))
print(S[S.three_c].groupby("net").size().to_dict())
