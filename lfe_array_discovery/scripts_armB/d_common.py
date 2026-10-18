#!/usr/bin/env python
"""lfe_array_discovery — shared constants and helpers for the download/metadata queue (agent D).

Geometry conventions are the project's standing ones (lfe_picker/scripts/u_common.py):
  * 1 km single-linkage site merge (co-located sensors are ONE site);
  * hypocentral distances;
  * the 20 Hz / 1-8 Hz cache chain is the frozen one (o_cache_month.cache_one) with the
    documented rational twin (u_cache.cache_one_rational) for non-integer decimation.
Nothing outside lfe_array_discovery/ is ever written.
"""
from __future__ import annotations
import os
import sys
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json                                                            # noqa: E402
import hashlib                                                         # noqa: E402
from datetime import datetime, timezone                                # noqa: E402
import numpy as np                                                     # noqa: E402
import pandas as pd                                                    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(HERE)
SCRIPTS = os.path.join(HERE, "scripts")
RES = os.path.join(HERE, "results")
QUEUES = os.path.join(RES, "queues")
MARKERS = os.path.join(RES, "markers")
XMLDIR = os.path.join(RES, "stationxml")
TRACES = os.path.join(HERE, "traces")
TRACES_SEALED = os.path.join(HERE, "traces_sealed")
TRACES_EQ = os.path.join(HERE, "traces_eq")
CACHE = os.path.join(HERE, "cache20hz")
LOGS = os.path.join(HERE, "logs")
FIGS = os.path.join(HERE, "figures")
PICKER = os.path.join(REPO, "lfe_picker")
EXT = os.path.join(REPO, "lfe_extension")
# read-only trees we SYMLINK from instead of re-fetching (never modified)
LINK_ROOTS = [os.path.join(EXT, "traces"), os.path.join(PICKER, "traces"),
              os.path.join(REPO, "lfe_catalog_repair", "traces")]
TREMOR = os.path.join(REPO, "catalogs", "pnsn_tremor_cascadia_full.csv")
SLAB2 = os.path.join(REPO, "data", "cas_slab2_input_04-18.csv")
INVENTORY = os.path.join(RES, "fdsn_stations_40-50N_2003-2014.txt")
MANIFEST = os.path.join(RES, "download_manifest.csv")

for _d in (RES, QUEUES, MARKERS, XMLDIR, TRACES, TRACES_SEALED, TRACES_EQ, CACHE, LOGS, FIGS):
    os.makedirs(_d, exist_ok=True)

# ------------------------------------------------------------------ frozen study constants
LAT0, LAT1, LON0, LON1 = 46.0, 48.6, -125.0, -121.5          # the box (task spec)
NETS = ["XU", "XG", "ZK", "Z3", "YW", "TA", "PB", "UW", "CN", "CC", "US"]
PERM_NETS = ["PB", "UW", "CN", "CC", "US", "TA"]             # metadata pulled to 2026 as well
META_T0, META_T1 = "2006-01-01", "2012-12-31"
PERM_T1 = "2026-12-31"
BAND_PREF = ["HH", "BH", "EH", "CH", "SH", "DP"]             # velocity seismometers only
LOCPREF = ["", "--", "00", "01", "10", "20"]
FS_CACHE = 20.0
FS_MIN = 20.0
MAX_CONCURRENT_FDSN = 8
SITE_MERGE_KM = 1.0
R_E = 111.19
CELL_KM = 5.0                                                # active-cell grid spacing
CELL_DEPTH = (28.0, 45.0)
CELL_TREMOR_R_KM, CELL_TREMOR_N = 8.0, 30                    # active = >=30 epicentres <=8 km
QUIET_RADIUS_KM, QUIET_PAD_D = 100.0, 3                      # task spec: 100 km, day +- 3 d
SHORT_PERIOD_BANDS = ("EH", "SH", "EP", "ES", "DP", "SP")
SHORT_PERIOD_PAT = ("L-4", "L4", "L22", "L-22", "S-13", "S13", "SS-1", "HS-1", "MARK",
                    "GEOPHONE", "4.5 HZ", "2 HZ", "1 HZ", "SHORT PERIOD", "SHORT-PERIOD",
                    "RANGER", "SERCEL", "GS-13", "GS13", "L28")

# extension roster (tranche 7): the N-Olympic permanents
EXT_ROSTER = [("PB", s) for s in ("B001", "B003", "B004", "B005", "B006", "B007", "B013",
                                  "B943", "B014", "B017", "B018")] + \
             [("UW", s) for s in ("HDW", "GNW", "LRIV", "DOSE")] + \
             [("CN", s) for s in ("VGZ", "PFB")]
UW_ARRAY_ERA = ["HDW", "GNW", "LRIV", "DOSE"]
CN_ARRAY_ERA = ["VGZ", "PFB"]

WINDOWS = {
    "t2_nolympic2010": ("2010-07-20", "2010-09-10"),
    "t3_sealed2011":   ("2011-07-20", "2011-09-10"),
    "t4_puget2008":    ("2008-04-15", "2008-06-10"),
    "t5_cafe2007":     ("2007-01-01", "2007-02-28"),
    "t6_intraslab":    ("2006-06-01", "2012-06-30"),
    "t7_extension":    ("2007-01-01", "2026-05-31"),
}


def log(m):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def marker(name, text=""):
    p = os.path.join(MARKERS, name)
    with open(p, "w") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}\n{text}\n")
    return p


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def script_hashes():
    out = {}
    for fn in sorted(os.listdir(SCRIPTS)):
        if fn.startswith("d_"):
            out[fn] = sha256_file(os.path.join(SCRIPTS, fn))
    return out


def jdump(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=str)


# ------------------------------------------------------------------ geometry
def proj(la, lo):
    la = np.asarray(la, float); lo = np.asarray(lo, float)
    return np.column_stack([lo * R_E * np.cos(np.deg2rad(la)), la * R_E])


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(lambda v: np.asarray(v, float), (la1, lo1, la2, lo2))
    return np.hypot((lo1 - lo2) * R_E * np.cos(np.deg2rad((la1 + la2) / 2.0)), (la1 - la2) * R_E)


def hyp_km(clat, clon, cdep, slat, slon):
    return np.hypot(epi_km(clat, clon, slat, slon), np.asarray(cdep, float))


def in_box(lat, lon):
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    return (lat >= LAT0) & (lat <= LAT1) & (lon >= LON0) & (lon <= LON1)


def cluster_sites(D, lat="lat", lon="lon", key="key", km=SITE_MERGE_KM):
    """1 km single-linkage site merge -> Series key -> site name (u_common.cluster_sites)."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
    D = D.reset_index(drop=True)
    if len(D) == 1:
        return pd.Series({D[key].iloc[0]: D[key].iloc[0].split(".")[-1]})
    xy = proj(D[lat].values, D[lon].values)
    cl = fcluster(linkage(pdist(xy), method="single"), t=km, criterion="distance")
    out = {}
    for c in np.unique(cl):
        g = D[cl == c]
        name = "+".join(sorted(g[key].str.split(".").str[-1]))
        for k in g[key]:
            out[k] = name
    return pd.Series(out)


def active_cells():
    """Slab2 interface cells (5 km grid, 28-45 km) in the box that are PNSN tremor-active."""
    sys.path.insert(0, os.path.join(PICKER, "scripts"))
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    sl = pd.read_csv(SLAB2, usecols=["lat", "lon", "depth", "etype"], low_memory=False)
    sl = sl[~sl.etype.isin(["BA", "TO"])]
    sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
    L = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
    N = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
    cs = np.cos(np.radians(0.5 * (LAT0 + LAT1)))
    lats = np.arange(LAT0, LAT1 + 1e-9, CELL_KM / R_E)
    lons = np.arange(LON0, LON1 + 1e-9, CELL_KM / (R_E * cs))
    LA, LO = np.meshgrid(lats, lons, indexing="ij")
    xy = np.column_stack([LO.ravel(), LA.ravel()])
    dep = L(xy)
    dep = np.where(np.isfinite(dep), dep, N(xy))
    G = pd.DataFrame(dict(cell_id=np.arange(xy.shape[0]), lat=LA.ravel(), lon=LO.ravel(),
                          depth_km=dep))
    G = G[(G.depth_km >= CELL_DEPTH[0]) & (G.depth_km <= CELL_DEPTH[1])].reset_index(drop=True)
    T = tremor_catalog()
    T = T[(T.time >= "2010-01-01") & (T.time < "2012-01-01")]
    from scipy.spatial import cKDTree
    tree = cKDTree(proj(T.lat.values, T.lon.values))
    n = np.array([len(tree.query_ball_point(p, CELL_TREMOR_R_KM))
                  for p in proj(G.lat.values, G.lon.values)])
    G["tremor_n8"] = n
    G["tremor_active"] = n >= CELL_TREMOR_N
    return G


def tremor_catalog():
    T = pd.read_csv(TREMOR, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time).astype("datetime64[ns]")        # memory: datetime64[us] trap
    return T


# ------------------------------------------------------------------ metadata tables
def channels(all_epochs=False):
    p = os.path.join(RES, "channels_all.csv" if all_epochs else "channels_2006_2012.csv")
    C = pd.read_csv(p, keep_default_na=False, dtype={"loc": str})
    C["loc"] = C["loc"].fillna("").astype(str).replace("nan", "")
    for c in ("fs", "azimuth", "dip", "sensitivity", "lat", "lon", "elev"):
        if c in C:
            C[c] = pd.to_numeric(C[c], errors="coerce")
    return C


def station_flags():
    return pd.read_csv(os.path.join(RES, "station_flags.csv"), keep_default_na=False,
                       dtype={"loc": str})


def sites():
    return pd.read_csv(os.path.join(RES, "sites_array_era.csv"))


def chain_for(fs):
    q = float(fs) / FS_CACHE
    return "frozen" if abs(q - round(q)) < 1e-9 else "rational"


def best_band_on(C, net, sta, day, allow_z_only=False):
    """Best 3-C (loc, band, fs) for one station on one ISO day, or None.

    C = channel table (any epoch coverage).  A (loc, band) qualifies when Z plus one of N|1 plus
    one of E|2 are declared on that day, unambiguously, at fs >= 20 Hz.  With allow_z_only the
    fallback is a vertical-only band (used ONLY for the explicitly named permanent stations -
    HDW is EHZ/EH2 through 2022 - and always tagged three_c=False in the manifest).
    """
    g = C[(C.net == net) & (C.sta == sta) & (C.start <= day) & (C.end >= day)]
    if not len(g):
        return None
    g = g[g.chan.str[-1].isin(list("ZNE12")) & g.band.isin(BAND_PREF) & (g.fs >= FS_MIN)]
    best = None
    for (loc, band), gg in g.groupby(["loc", "band"]):
        slots = {c[-1] for c in gg.chan}
        if "Z" not in slots or not ({"N", "1"} & slots) or not ({"E", "2"} & slots):
            continue
        if len({"N", "1"} & slots) > 1 or len({"E", "2"} & slots) > 1:
            continue
        key = (BAND_PREF.index(band), LOCPREF.index(loc) if loc in LOCPREF else 99, loc)
        if best is None or key < best[0]:
            best = (key, loc, band, float(gg.fs.max()))
    if best is None:
        if not allow_z_only:
            return None
        gz = g[g.chan.str[-1] == "Z"]
        if not len(gz):
            return None
        for (loc, band), gg in gz.groupby(["loc", "band"]):
            key = (BAND_PREF.index(band), LOCPREF.index(loc) if loc in LOCPREF else 99, loc)
            if best is None or key < best[0]:
                best = (key, loc, band, float(gg.fs.max()))
        return dict(loc=best[1], band=best[2], fs=best[3], chain=chain_for(best[3]),
                    three_c=False)
    return dict(loc=best[1], band=best[2], fs=best[3], chain=chain_for(best[3]), three_c=True)
