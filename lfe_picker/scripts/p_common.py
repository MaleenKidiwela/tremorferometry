#!/usr/bin/env python
"""P1/P2 shared infrastructure for the learned LFE picker.

Everything the learned models and the physics back-projection need that is NOT already frozen in the
O-stage scripts:

  * the S15 station roster, its ZNE rotation azimuths (StationXML), the site merge;
  * a memory-mapped float16 store of every cached station-day at the frozen 20 Hz / 1-8 Hz chain,
    ROTATED TO ZNE (the template path in o_* is untouched and still uses the raw Z/1/2 order);
  * the Bostock event table on the frozen time base (DB time + frozen 2012 family offset);
  * Crosson 1-D S travel times + the a5 station corrections, and the Slab2 back-projection grid.

Read-only inputs: ../lfe_extension, ../lfe_catalog_repair, ../codex-LFE, lfe_picker/results/o*.
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
CODEX = "/home/jovyan/tremorferometry/codex-LFE"
IN = f"{REP}/results"
DATA = "/home/jovyan/tremorferometry/data"
BOST = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
CACHES = [f"{PICK}/cache20hz", f"{REP}/cache20hz", f"{EXT}/cache20hz"]
MM = f"{RES}/mm"

FS = 20.0
NDAY = 1728000                       # samples in one 20 Hz day
R = 111.195
S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
       "B009", "B010", "B011", "B926"]
NETSTA = {"PGC": "CN.PGC", "VGZ": "CN.VGZ", "NLLB": "CN.NLLB", "SNB": "CN.SNB", "LRIV": "UW.LRIV",
          "B001": "PB.B001", "B003": "PB.B003", "B004": "PB.B004", "B005": "PB.B005",
          "B006": "PB.B006", "B007": "PB.B007", "B009": "PB.B009", "B010": "PB.B010",
          "B011": "PB.B011", "B926": "PB.B926"}
# co-located clusters merged into one SITE (same convention as o3_oracle / o_common.SITES)
SITE_CLUSTERS = {"PGCq": ["PGC", "B009", "B010", "B011"], "LRIVq": ["LRIV", "B005", "B006", "B007"]}
HELD_OUT_STATIONS = ["B926", "NLLB", "LRIV"]
TRAIN_YEARS = [2008, 2009, 2010, 2012]
DEV_YEARS = [2010, 2012]
GATE_YEAR = 2011

# EH1/EH2 azimuths (deg CW from N, dip 0) from StationXML.  CN/UW surface stations are BHN/BHE or
# HHN/HHE at 0/90 exactly (a4_stationxml_audit.csv), so only the boreholes need rotating.
AZ12 = {
    "B001": (335.5, 65.5), "B003": (239.0, 329.0), "B004": (167.5, 257.5),
    "B005": (90.5, 180.5), "B006": (201.0, 291.0), "B007": (36.5, 126.5),
    "B009": (0.0, 90.0), "B010": (228.5, 318.5), "B011": (11.0, 101.0),
    "B926": (299.0, 29.0),
    "PGC": (0.0, 90.0), "VGZ": (0.0, 90.0), "NLLB": (0.0, 90.0), "SNB": (0.0, 90.0),
    "LRIV": (0.0, 90.0),
}


def site_map():
    sm = {s: s for s in S15}
    for q, mem in SITE_CLUSTERS.items():
        for s in mem:
            sm[s] = q
    return sm


SITE_OF = site_map()
SITES = sorted(set(SITE_OF.values()))


def station_table():
    return pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")


def epi_km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1) - np.asarray(lo2)) * R * np.cos(np.radians(48.6)),
                    (np.asarray(la1) - np.asarray(la2)) * R)


def rot_matrix(sta):
    """2x2 that maps the recorded (c1, c2) onto (N, E).

    A component with azimuth a records  c = N cos a + E sin a.  With two orthogonal horizontals the
    2x2 [[cos a1, sin a1],[cos a2, sin a2]] is orthonormal, so its transpose is the inverse."""
    a1, a2 = np.radians(AZ12[sta])
    A = np.array([[np.cos(a1), np.sin(a1)], [np.cos(a2), np.sin(a2)]])
    return np.linalg.inv(A)                      # exact, also for the 0/90 identity case


# ------------------------------------------------------------------ day store ------------------
def cache_index(years=None, stas=None):
    """{(sta, year, doy): path} over every 20 Hz cache root, first root wins."""
    idx = {}
    for root in CACHES:
        for p in glob.glob(f"{root}/*/*/*.npz"):
            parts = p.split(os.sep)
            netsta, yr, fn = parts[-3], int(parts[-2]), parts[-1]
            sta = netsta.split(".")[1]
            if sta not in S15:
                continue
            if stas is not None and sta not in stas:
                continue
            if years is not None and yr not in years:
                continue
            idx.setdefault((sta, yr, int(fn[:3])), p)
    return idx


def day_index(tag=""):
    return pd.read_parquet(f"{MM}/day_index{tag}.parquet")


def open_days(mode="r", tag=""):
    """(index DataFrame, memmap of shape (Nrows, 3, NDAY) float16, ZNE, per-day scaled).

    `di.row` is the positional row in the memmap (failed station-days leave a hole)."""
    di = day_index(tag)
    f = f"{MM}/days_zne_f16{tag}.dat"
    nrow = int(os.path.getsize(f) // (3 * NDAY * 2))
    X = np.memmap(f, dtype=np.float16, mode=mode, shape=(nrow, 3, NDAY))
    return di, X


class DayStore:
    """Random access to (sta, year, doy) -> (3, NDAY) float16 ZNE, plus absolute-time indexing."""

    def __init__(self, mode="r", tag=""):
        self.di, self.X = open_days(mode, tag)
        self.key = {(r.sta, r.year, r.doy): i for i, r in enumerate(self.di.itertuples())}
        self.row = self.di.row.to_numpy()
        self.t0 = self.di.t0.to_numpy()
        self.nvalid = self.di.n.to_numpy()

    def has(self, sta, year, doy):
        return (sta, year, doy) in self.key

    def idx(self, sta, year, doy):
        """index into di / t0 / nvalid (NOT the memmap row)."""
        return self.key[(sta, year, doy)]

    def get(self, sta, year, doy):
        return self.X[self.row[self.key[(sta, year, doy)]]]

    def get_i(self, i):
        return self.X[self.row[i]]

    def sample_of(self, i, t_abs):
        """sample index inside the day of station-day `i` (index into di) for absolute time."""
        return int(round((float(t_abs) - self.t0[i]) * FS))


# ------------------------------------------------------------------ catalog -------------------
def families():
    return pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")


def holdouts():
    ho = json.load(open(f"{RES}/holdouts.json"))
    return {int(h["bkey"]) for h in ho["hidden_families"]}, set(ho["held_out_stations"])


def bostock():
    """Bostock events on the FROZEN time base: t_origin = DB time + frozen 2012 family offset."""
    FAM = families()
    OFF = pd.read_csv(f"{IN}/e5r_family_offsets_2012.csv")
    omap = {int(r.bkey): float(r.offset_s) for r in OFF.itertuples() if bool(r.accept)}
    B = pd.read_csv(BOST)
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour,
                                    minute=B["min"], second=B.sec.astype(int)), errors="coerce") \
        + pd.to_timedelta(B.sec % 1.0, unit="s")
    key = {(round(FAM.lat[b], 4), round(FAM.lon[b], 4), round(FAM.dep[b], 4)): int(b)
           for b in FAM.index}
    B["bkey"] = [key.get((round(la, 4), round(lo, 4), round(dp, 4)), -1)
                 for la, lo, dp in zip(B.lat, B.lon, B.dep)]
    B = B[B.bkey >= 0].copy()
    B["t_origin"] = B.time.astype("datetime64[ns]").astype("int64") / 1e9 + B.bkey.map(omap)
    B = B.dropna(subset=["t_origin"]).reset_index(drop=True)
    B["doy"] = B.time.dt.dayofyear
    hid, _ = holdouts()
    B["hidden"] = B.bkey.isin(hid)
    return B


# ------------------------------------------------------------------ travel times ---------------
def station_corr():
    c = pd.read_csv(f"{EXT}/results/a5_station_corrections.csv").set_index("sta")
    return {s: (0.0 if not np.isfinite(c.corr_S.get(s, np.nan)) else float(c.corr_S[s]))
            for s in S15}


def slab2_interp():
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    sl = pd.read_csv(f"{DATA}/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
    sl = sl[~sl.etype.isin(["BA", "TO"])]
    sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
    L = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
    N = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)

    def f(lat, lon):
        xy = np.column_stack([np.asarray(lon, float).ravel(), np.asarray(lat, float).ravel()])
        d = L(xy)
        d = np.where(np.isfinite(d), d, N(xy))
        return np.clip(d, 12.0, 55.0)
    return f


# Merlin P2 grid: 2-D cells of 3 km over 48.0-49.0 N / -124.5..-123.0 E, depth fixed on Slab2.
GRID_LAT = (48.0, 49.0)
GRID_LON = (-124.5, -123.0)
GRID_KM = 3.0


def build_grid(dlat=GRID_LAT, dlon=GRID_LON, km=GRID_KM):
    cs = np.cos(np.radians(48.5))
    nlat = int(round((dlat[1] - dlat[0]) * R / km)) + 1
    nlon = int(round((dlon[1] - dlon[0]) * R * cs / km)) + 1
    lats = np.linspace(dlat[0], dlat[1], nlat)
    lons = np.linspace(dlon[0], dlon[1], nlon)
    LA, LO = np.meshgrid(lats, lons, indexing="ij")
    dep = slab2_interp()(LA.ravel(), LO.ravel()).reshape(LA.shape)
    g = pd.DataFrame(dict(cell=np.arange(LA.size), ilat=np.repeat(np.arange(nlat), nlon),
                          ilon=np.tile(np.arange(nlon), nlat),
                          lat=LA.ravel(), lon=LO.ravel(), dep=dep.ravel()))
    g.attrs["shape"] = (nlat, nlon)
    return g, (nlat, nlon), lats, lons


def tt_table(grid, stas=S15):
    """(ncell, nsta) predicted S arrival delay after origin, Crosson + a5 station correction."""
    sys.path.insert(0, HERE)
    import p_vel1d as V
    ST = station_table()
    CR = station_corr()
    out = np.zeros((len(grid), len(stas)), np.float64)
    # Crosson tt is a slow 1-D shoot: tabulate on (distance, depth) and interpolate bilinearly.
    dgrid = np.arange(0.0, 220.0, 1.0)
    zgrid = np.arange(10.0, 60.001, 1.0)
    T = np.empty((dgrid.size, zgrid.size))
    for i, d in enumerate(dgrid):
        for j, z in enumerate(zgrid):
            T[i, j] = V.tt(float(d), float(z), "S")
    from scipy.interpolate import RegularGridInterpolator
    F = RegularGridInterpolator((dgrid, zgrid), T, bounds_error=False, fill_value=None)
    for k, s in enumerate(stas):
        d = epi_km(grid.lat.values, grid.lon.values, ST.lat[s], ST.lon[s])
        out[:, k] = F(np.column_stack([np.clip(d, 0, 219), np.clip(grid.dep.values, 10, 60)])) + CR[s]
    return out


def fam_tt(stas=S15):
    """{(bkey, sta): predicted S delay} at the 130 Bostock family hypocentres, same recipe."""
    sys.path.insert(0, HERE)
    import p_vel1d as V
    FAM, ST, CR = families(), station_table(), station_corr()
    out = {}
    for b in FAM.index:
        for s in stas:
            d = epi_km(FAM.lat[b], FAM.lon[b], ST.lat[s], ST.lon[s])
            out[(int(b), s)] = float(V.tt(float(d), float(FAM.dep[b]), "S")) + CR[s]
    return out


def fam_tt_labels(stas=S15):
    """{(bkey, sta): (S delay after origin, provenance)} EXACTLY as o_labels.py built labels_v0:
    the a5 stack S pick with its station correction added back where one exists, else Crosson ttS
    plus the same correction.  Needed to recover event-level (not row-level) label times."""
    PK = pd.read_csv(f"{EXT}/results/a5_picks.csv")
    CORR = station_corr()
    QC = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    HASTPL = {(int(r.bkey), str(r.sta)) for r in QC.itertuples() if r.usable and r.template}
    FAM, ST = families(), station_table()
    sys.path.insert(0, HERE)
    import p_vel1d as V
    tt = {}
    for r in PK.itertuples():
        s, b = str(r.sta), int(r.bkey)
        if s not in stas:
            continue
        c = CORR.get(s, 0.0)
        tt[(b, s)] = ((float(r.t_S) + c, "picked") if np.isfinite(r.t_S)
                      else (float(r.ttS) + c, "predicted"))
    for (b, s) in HASTPL:
        if (b, s) in tt or s not in stas or b not in FAM.index:
            continue
        d = epi_km(FAM.lat[b], FAM.lon[b], ST.lat[s], ST.lon[s])
        tt[(b, s)] = (float(V.tt(d, float(FAM.dep[b]), "S")) + CORR.get(s, 0.0), "predicted")
    return tt
