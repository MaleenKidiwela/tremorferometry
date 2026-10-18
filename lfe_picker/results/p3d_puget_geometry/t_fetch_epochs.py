"""T-study step 1: build the post-2020 station/epoch table for the Puget-Oregon geometry study.

Read-only inputs:
  ../lfe_catalog_repair/results/f1_stations.csv      (superset of the 183-station roster and the
                                                      post-2020 candidate list; role/roster/klass)
  ../codex-LFE/results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv   (3-C epochs)
  IRIS StationXML level=channel for every pool station with no triplet epoch (no waveforms).
Output: results/p3d_puget_geometry/stations_epochs.csv  (one row per station, 3-C interval list).
"""
import json, os, sys
import numpy as np, pandas as pd
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
CACHE = "/tmp/claude-1000/-home-jovyan-tremorferometry/d6a923e6-f9f1-4ee1-8d37-955b6b2e0f35/scratchpad/iris_chan_cache.csv"
BOX = dict(minlat=42.6, maxlat=50.0, minlon=-126.8, maxlon=-120.3)
T0, T1 = "2019-06-01", "2026-09-19"
VEL_INST = set("HL")          # channel instrument code: H = seismometer, L = low gain; drop N (accel)

# ---------------- pool ----------------
F = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/f1_stations.csv")
F["key"] = F.net + "." + F.sta
pool = F[(F.lat.between(BOX["minlat"], BOX["maxlat"])) &
         (F.lon.between(BOX["minlon"], BOX["maxlon"]))].copy()
print("pool in box:", len(pool), pool.roster.value_counts().to_dict())

# ---------------- triplet epochs (project roster) ----------------
E = pd.read_csv(f"{ROOT}/codex-LFE/results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv")
E["key"] = E.network + "." + E.station
trip = {}
for k, g in E.groupby("key"):
    trip[k] = [(str(a)[:10], str(b)[:10]) for a, b in zip(g.start_utc, g.end_utc_exclusive)]
print("stations with triplet epochs:", len(trip), "; of pool:", pool.key.isin(trip).sum())

# ---------------- IRIS channel metadata for everything else ----------------
need = sorted(set(pool.key))          # query all: gives channel depth + Z-only flag + true epochs
if os.path.exists(CACHE):
    CH = pd.read_csv(CACHE)
    print("cache hit:", len(CH), "channel epochs")
else:
    cli = Client("IRIS", timeout=120)
    rows = []
    nets = sorted({k.split(".")[0] for k in need})
    for net in nets:
        stas = sorted({k.split(".")[1] for k in need if k.startswith(net + ".")})
        for i in range(0, len(stas), 60):
            chunk = ",".join(stas[i:i + 60])
            try:
                inv = cli.get_stations(network=net, station=chunk, level="channel",
                                       starttime=UTCDateTime(T0), endtime=UTCDateTime(T1))
            except Exception as ex:
                print("  FAIL", net, i, ex); continue
            for n in inv:
                for s in n:
                    for c in s:
                        rows.append(dict(net=n.code, sta=s.code, loc=c.location_code or "",
                                         chan=c.code, lat=c.latitude, lon=c.longitude,
                                         elev=c.elevation, depth=c.depth,
                                         start=str(c.start_date)[:10],
                                         end=(str(c.end_date)[:10] if c.end_date else "2026-09-19"),
                                         fs=float(c.sample_rate or 0)))
            print("  ", net, i, len(rows))
    CH = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    CH.to_csv(CACHE, index=False)
    print("IRIS channel epochs:", len(CH))

CH["key"] = CH.net + "." + CH.sta
CH = CH[CH.chan.str.len() == 3]
CH["band"] = CH.chan.str[0]; CH["inst"] = CH.chan.str[1]; CH["comp"] = CH.chan.str[2]
CH = CH[CH.inst.isin(VEL_INST) & CH.band.isin(list("BHES"))]

def three_c_intervals(g):
    """3-C intervals per (loc, band): [max over comps of comp-start, min over comps of comp-end]."""
    out = []
    for (loc, band), gg in g.groupby(["loc", "band"]):
        comps = {}
        for c, ggg in gg.groupby("comp"):
            comps[c] = (min(ggg.start), max(ggg.end))
        z = comps.get("Z")
        hor = [c for c in ("N", "E", "1", "2") if c in comps]
        if z is None or len(hor) < 2:
            continue
        s = max([z[0]] + [comps[c][0] for c in hor[:2]])
        e = min([z[1]] + [comps[c][1] for c in hor[:2]])
        if s < e:
            out.append((s, e))
    return out

rows = []
for r in pool.itertuples():
    g = CH[CH.key == r.key]
    iv = three_c_intervals(g) if len(g) else []
    src = "iris"
    if not iv and r.key in trip:
        iv = trip[r.key]; src = "triplet"
    elif r.key in trip:                       # union both sources (IRIS is authoritative, triplet fills)
        iv = iv + trip[r.key]
    has_z = bool(len(g[g.comp == "Z"]))
    depths = g.depth.dropna().values
    is_bh = bool(r.is_borehole_any) or (len(depths) and float(np.nanmax(depths)) >= 20.0) or r.net == "PB"
    bands = sorted(set(g.band)) if len(g) else []
    bb = any(b in ("B", "H") for b in bands) or str(r.klass).startswith("broadband")
    rows.append(dict(key=r.key, net=r.net, sta=r.sta, lat=r.lat, lon=r.lon, elev=r.elev,
                     roster=r.roster, role=r.role, klass=r.klass, is_borehole=is_bh,
                     is_broadband=bool(bb), has_z=has_z, n_3c_iv=len(iv), epoch_src=src,
                     intervals=json.dumps(sorted(iv))))
S = pd.DataFrame(rows)
S["three_c"] = S.n_3c_iv > 0
S["z_only"] = S.has_z & ~S.three_c
S.to_csv(f"{OUT}/stations_epochs.csv", index=False)
print(S.groupby(["roster"])[["three_c", "z_only"]].sum())
print("total", len(S), "3-C", S.three_c.sum(), "Z-only", S.z_only.sum(),
      "no metadata", int((~S.has_z & ~S.three_c).sum()))
