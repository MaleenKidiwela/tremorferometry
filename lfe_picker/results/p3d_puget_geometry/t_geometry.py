"""T-study step 2: hypocentral network density over the Cascadia LFE source zones, 2020-2026.

Cells = Slab2 interface points on a ~5 km grid with interface depth 28-45 km, in three regions:
  North  48.0-49.0 N (the zone where discovery was actually attempted)
  Puget  47.0-48.0 N
  Oregon 43.5-46.0 N
For each cell and year: 3-C SITES (co-located stations within 1 km merged) alive that year, counted
within 45/55/60 km HYPOCENTRAL distance, nearest-site hypocentral distance, azimuthal gap among
sites <= 60 km, and the surface-broadband / borehole split.
North reference: the same statistics per Bostock family (e5r_families_2012.csv) for the S15 network
and for S15+SYMB+BUCKS -- the network where discovery failed.
Outputs: cells.csv, north_families.csv, summary.csv (+ SUMMARY.md, figures/t_puget_geometry.png).
"""
import json
import numpy as np, pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
R_E, CS = 111.19, np.cos(np.deg2rad(46.0))
YEARS = list(range(2020, 2027))
YEND = {y: f"{y}-12-31" for y in YEARS}; YEND[2026] = "2026-09-19"
MIN_DAYS = 90.0                      # a station counts for a year only if alive >= 90 d in it
SITE_MERGE_KM = 1.0
REGIONS = {"North_48_49": (48.0, 49.0), "Puget_47_48": (47.0, 48.0), "Oregon_43.5_46": (43.5, 46.0)}
DMIN, DMAX, STEP_KM = 28.0, 45.0, 5.0


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(lambda v: np.asarray(v, float), (la1, lo1, la2, lo2))
    return np.hypot((lo1 - lo2) * R_E * np.cos(np.deg2rad((la1 + la2) / 2.0)), (la1 - la2) * R_E)


def hyp_km(clat, clon, cdep, slat, slon):
    return np.hypot(epi_km(clat, clon, slat, slon), float(cdep))


def azimuth(clat, clon, slat, slon):
    dx = (np.asarray(slon, float) - clon) * R_E * np.cos(np.deg2rad(clat))
    dy = (np.asarray(slat, float) - clat) * R_E
    return np.degrees(np.arctan2(dx, dy)) % 360.0


def az_gap(az):
    if len(az) == 0:
        return 360.0
    if len(az) == 1:
        return 360.0
    a = np.sort(np.asarray(az))
    d = np.diff(a)
    return float(max(d.max(), 360.0 - (a[-1] - a[0])))


# ---------------- Slab2 interface (same recipe as fault_tomography/inversion/assemble_res_catalog.py)
sl = pd.read_csv(f"{ROOT}/data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]
sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
_SLL = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
_SLN = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)


def slab_depth(clat, clon):
    xy = np.column_stack([np.asarray(clon, float), np.asarray(clat, float)])
    d = _SLL(xy); dn = _SLN(xy)
    return np.clip(np.where(np.isfinite(d), d, dn), 12.0, 55.0)


def make_cells():
    rows = []
    dlat = STEP_KM / R_E
    for reg, (la0, la1) in REGIONS.items():
        lats = np.arange(la0, la1 + 1e-9, dlat)
        for la in lats:
            dlon = STEP_KM / (R_E * np.cos(np.deg2rad(la)))
            lons = np.arange(-125.6, -120.9, dlon)
            dep = slab_depth(np.full_like(lons, la), lons)
            k = (dep >= DMIN) & (dep <= DMAX)
            for lo, dd in zip(lons[k], dep[k]):
                rows.append(dict(region=reg, lat=round(float(la), 4), lon=round(float(lo), 4),
                                 depth_km=round(float(dd), 2)))
    C = pd.DataFrame(rows)
    C.insert(0, "cell_id", [f"{r.region.split('_')[0][:3].upper()}{i:04d}" for i, r in enumerate(C.itertuples())])
    return C


# ---------------- sites ----------------
def build_sites(ST):
    """Merge stations within SITE_MERGE_KM into one site. ST needs lat/lon/is_borehole/is_broadband."""
    ST = ST.reset_index(drop=True).copy()
    xy = np.column_stack([ST.lon.values * R_E * np.cos(np.deg2rad(ST.lat.values)), ST.lat.values * R_E])
    if len(ST) > 1:
        lab = fcluster(linkage(pdist(xy), method="single"), t=SITE_MERGE_KM, criterion="distance")
    else:
        lab = np.array([1])
    ST["site_id"] = lab
    rows = []
    for sid, g in ST.groupby("site_id"):
        surf_bb = bool((~g.is_borehole & g.is_broadband).any())
        bh = bool(g.is_borehole.any())
        kind = "surface_bb" if surf_bb else ("borehole" if bh else "surface_sp")
        rows.append(dict(site=("+".join(sorted(g.sta))), lat=float(g.lat.mean()), lon=float(g.lon.mean()),
                         kind=kind, n_members=len(g), members=";".join(sorted(g.key)),
                         nets=";".join(sorted(set(g.net))),
                         intervals=json.dumps(sum([json.loads(i) for i in g.intervals], []))))
    return pd.DataFrame(rows)


def alive_days(intervals, y):
    a, b = np.datetime64(f"{y}-01-01"), np.datetime64(YEND[y])
    tot = 0.0
    for s, e in intervals:
        s, e = np.datetime64(s), np.datetime64(e)
        ov = (min(e, b) - max(s, a)) / np.timedelta64(1, "D")
        tot += max(0.0, float(ov))
    return tot


def site_year_mask(SI, y):
    return np.array([alive_days(json.loads(i), y) >= MIN_DAYS for i in SI.intervals])


def metrics(clat, clon, cdep, SI, alive):
    S = SI[alive]
    if len(S) == 0:
        return dict(n45=0, n55=0, n60=0, nearest_km=np.nan, gap_deg=360.0,
                    n45_surface_bb=0, n45_borehole=0, n45_surface_sp=0,
                    n60_surface_bb=0, n60_borehole=0, nearest_site="")
    d = hyp_km(clat, clon, cdep, S.lat.values, S.lon.values)
    k60 = d <= 60; k45 = d <= 45
    az = azimuth(clat, clon, S.lat.values[k60], S.lon.values[k60])
    kind = S.kind.values
    j = int(np.argmin(d))
    return dict(n45=int(k45.sum()), n55=int((d <= 55).sum()), n60=int(k60.sum()),
                nearest_km=round(float(d[j]), 2), gap_deg=round(az_gap(az), 1),
                n45_surface_bb=int((k45 & (kind == "surface_bb")).sum()),
                n45_borehole=int((k45 & (kind == "borehole")).sum()),
                n45_surface_sp=int((k45 & (kind == "surface_sp")).sum()),
                n60_surface_bb=int((k60 & (kind == "surface_bb")).sum()),
                n60_borehole=int((k60 & (kind == "borehole")).sum()),
                nearest_site=str(S.site.values[j]))


# ---------------- station pools ----------------
ST = pd.read_csv(f"{OUT}/stations_epochs.csv")
ST3 = ST[ST.three_c].copy()
POOLS = {"all": ST3, "project": ST3[ST3.roster == "project"]}
SITES = {k: build_sites(v) for k, v in POOLS.items()}
ZONLY = ST[ST.z_only].copy()
for k, v in SITES.items():
    print(f"pool {k}: {len(POOLS[k])} 3-C stations -> {len(v)} sites", v.kind.value_counts().to_dict())
SITES["all"].to_csv(f"{OUT}/sites_all.csv", index=False)
SITES["project"].to_csv(f"{OUT}/sites_project.csv", index=False)
ZONLY[["key", "net", "sta", "lat", "lon", "klass", "roster"]].to_csv(f"{OUT}/z_only_stations.csv", index=False)

CELLS = make_cells()
# tremor density: PNSN epicentres (2010-2026) within 20 km of the cell -- flags the cells that are in
# a real LFE/tremor source zone, as opposed to a merely geometric slice of the interface.
TR = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["lat", "lon"])
from scipy.spatial import cKDTree
proj = lambda la, lo: np.column_stack([np.asarray(lo) * R_E * np.cos(np.deg2rad(np.asarray(la))),
                                       np.asarray(la) * R_E])
tree = cKDTree(proj(TR.lat.values, TR.lon.values))
CELLS["tremor_n20"] = [len(i) for i in tree.query_ball_point(proj(CELLS.lat.values, CELLS.lon.values), 20.0)]
CELLS["tremor_active"] = CELLS.tremor_n20 >= 100
print("cells:", CELLS.region.value_counts().to_dict(),
      "tremor-active:", CELLS.groupby("region").tremor_active.sum().to_dict())

rows = []
for pool, SI in SITES.items():
    masks = {y: site_year_mask(SI, y) for y in YEARS}
    zmask = {y: np.array([alive_days(json.loads(i), y) >= MIN_DAYS for i in ZONLY.intervals])
             for y in YEARS} if len(ZONLY) else {}
    for c in CELLS.itertuples():
        for y in YEARS:
            m = metrics(c.lat, c.lon, c.depth_km, SI, masks[y])
            nz = 0
            if len(ZONLY):
                Z = ZONLY[zmask[y]]
                if len(Z):
                    nz = int((hyp_km(c.lat, c.lon, c.depth_km, Z.lat.values, Z.lon.values) <= 45).sum())
            rows.append(dict(pool=pool, cell_id=c.cell_id, region=c.region, lat=c.lat, lon=c.lon,
                             depth_km=c.depth_km, tremor_n20=c.tremor_n20,
                             tremor_active=c.tremor_active, year=y, n45_zonly=nz, **m))
CC = pd.DataFrame(rows)
CC.to_csv(f"{OUT}/cells.csv", index=False)
print("cells.csv rows:", len(CC))

# ---------------- North reference: Bostock families with S15 (+SYMB/BUCKS) ----------------
FAM = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv")
S15 = pd.read_csv(f"{ROOT}/lfe_extension/results/station_epochs.csv")
S15 = S15.groupby(["net", "sta"]).agg(lat=("lat", "first"), lon=("lon", "first"),
                                      start=("start", "min"), end=("end", "max")).reset_index()
extra = pd.DataFrame([dict(net="CN", sta="SYMB", lat=48.5593, lon=-123.7989, start="2016-09-19", end="2026-09-19"),
                      dict(net="UW", sta="BUCKS", lat=48.15541, lon=-123.83423, start="2023-09-17", end="2026-09-19")])


def prep(df):
    d = df.copy()
    d["key"] = d.net + "." + d.sta
    d["is_borehole"] = d.net == "PB"
    d["is_broadband"] = d.net != "PB"
    d["intervals"] = [json.dumps([[s, e]]) for s, e in zip(d.start, d.end)]
    return d


# the P3D frozen 15 (o_common.STAS_ORDER) -- the network P3 discovery actually ran on; kept as a
# separate variant because lfe_extension/station_epochs.csv also carries four CN/C8 stations
# (SOKB, LZB, PFB, YOUB) that were dead by 2019.
P15 = pd.read_csv(f"{ROOT}/lfe_picker/results/p3d/stations_p3d.csv")
P15 = P15[P15.role == "frozen"][["net", "sta", "lat", "lon"]].copy()
P15["start"], P15["end"] = "2003-01-01", "2026-09-19"
NETS = {"S15": prep(S15), "S15+SYMB+BUCKS": prep(pd.concat([S15, extra], ignore_index=True)),
        "P3Dfrozen15": prep(P15),
        "P3Dfrozen15+SYMB+BUCKS": prep(pd.concat([P15, extra], ignore_index=True))}
frows = []
for name, tab in NETS.items():
    SI = build_sites(tab)
    print(f"north {name}: {len(tab)} stations -> {len(SI)} sites")
    masks = {y: site_year_mask(SI, y) for y in YEARS}
    masks[0] = np.ones(len(SI), bool)          # year 0 = epoch-blind (the discovery-era geometry:
    for f in FAM.itertuples():                 # most of S15 was dead by 2020, so 2020+ masking
        for y in [0] + YEARS:                  # would flatter nothing -- both views are reported)
            frows.append(dict(network=name, bkey=int(f.bkey), lat=f.lat, lon=f.lon, depth_km=f.dep,
                              year=y, **metrics(f.lat, f.lon, f.dep, SI, masks[y])))
FF = pd.DataFrame(frows)
FF.to_csv(f"{OUT}/north_families.csv", index=False)

# ---------------- summary ----------------
STATS = ["n45", "n55", "n60", "nearest_km", "gap_deg", "n45_surface_bb", "n45_borehole"]


def agg(g):
    out = dict(n_units=len(g))
    for s in STATS:
        v = g[s].astype(float)
        out[f"{s}_p10"] = round(float(np.nanpercentile(v, 10)), 2)
        out[f"{s}_med"] = round(float(np.nanmedian(v)), 2)
        out[f"{s}_p90"] = round(float(np.nanpercentile(v, 90)), 2)
    out["n45_mean"] = round(float(np.nanmean(g.n45)), 2)
    out["frac_n45_ge3"] = round(float((g.n45 >= 3).mean()), 3)
    out["frac_n45_ge5"] = round(float((g.n45 >= 5).mean()), 3)
    return pd.Series(out)


srows = []
for (pool, reg, y), g in CC.groupby(["pool", "region", "year"]):
    srows.append(pd.concat([pd.Series(dict(kind="cells", set=pool, region=reg, year=y)), agg(g)]))
for (nm, y), g in FF.groupby(["network", "year"]):
    srows.append(pd.concat([pd.Series(dict(kind="bostock_families", set=nm, region="North_ref_48_49", year=y)), agg(g)]))
for (pool, reg, y), g in CC[CC.tremor_active].groupby(["pool", "region", "year"]):
    srows.append(pd.concat([pd.Series(dict(kind="cells_tremor_active", set=pool, region=reg, year=y)), agg(g)]))
# per 0.5-degree latitude band (the "where" question)
CC["band"] = (np.floor(CC.lat * 2) / 2).round(1)
for (pool, reg, band, y), g in CC[CC.year.isin([2024])].groupby(["pool", "region", "band", "year"]):
    srows.append(pd.concat([pd.Series(dict(kind="cells_latband", set=pool,
                                           region=f"{reg}|{band:.1f}N", year=y)), agg(g)]))
SU = pd.DataFrame(srows)
SU.to_csv(f"{OUT}/summary.csv", index=False)
print(SU[(SU.year == 2024)][["kind", "set", "region", "n_units", "n45_med", "n45_p90", "n55_med",
                             "n60_med", "nearest_km_med", "gap_deg_med"]].to_string())
