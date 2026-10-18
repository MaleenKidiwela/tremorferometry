"""T-study (revised brief) step 2: hypocentral geometry of the FULL 2019+ network in the box
46.5-49.5 N / -125.5..-121.5 E, against (1) the 130 Bostock families and (2) Slab2 interface cells,
for 2019 / 2021 / 2023 / 2025, plus (3) a ranked list of stations outside S15+SYMB/BUCKS.

Outputs: sites.csv, family_geometry_by_year.csv, cell_geometry_by_year.csv, station_ranking.csv
"""
import json
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
R_E = 111.19
YEARS = [2019, 2021, 2023, 2025]
YEND = {y: f"{y}-12-31" for y in YEARS}
MIN_DAYS, SITE_MERGE_KM = 90.0, 1.0
BOX = (46.5, 49.5, -125.5, -121.5)
DMIN, DMAX, STEP_KM = 28.0, 45.0, 5.0
S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
       "B009", "B010", "B011", "B926"]
BASE_PLUS = S15 + ["SYMB", "BUCKS"]
PUGET_BAND = (46.5, 48.0)          # "Olympics / Puget" for the ranking tally


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(lambda v: np.asarray(v, float), (la1, lo1, la2, lo2))
    return np.hypot((lo1 - lo2) * R_E * np.cos(np.deg2rad((la1 + la2) / 2.0)), (la1 - la2) * R_E)


def hyp_km(clat, clon, cdep, slat, slon):
    return np.hypot(epi_km(clat, clon, slat, slon), np.asarray(cdep, float))


def az_gap(az):
    if len(az) < 2:
        return 360.0
    a = np.sort(np.asarray(az))
    return float(max(np.diff(a).max(), 360.0 - (a[-1] - a[0])))


def azimuth(clat, clon, slat, slon):
    dx = (np.asarray(slon, float) - clon) * R_E * np.cos(np.deg2rad(clat))
    return np.degrees(np.arctan2(dx, (np.asarray(slat, float) - clat) * R_E)) % 360.0


# ---------------- sites ----------------
ST = pd.read_csv(f"{OUT}/stations_all_2019plus.csv")
DET = ST[ST.detect_capable].reset_index(drop=True).copy()
ACC = ST[~ST.detect_capable].reset_index(drop=True).copy()     # strong-motion / citizen, flagged only
xy = np.column_stack([DET.lon * R_E * np.cos(np.deg2rad(DET.lat)), DET.lat * R_E])
DET["cl"] = fcluster(linkage(pdist(xy), method="single"), t=SITE_MERGE_KM, criterion="distance")
rows = []
for cl, g in DET.groupby("cl"):
    prio = {"surface_broadband": 0, "borehole": 1, "surface_short_period": 2}
    g = g.sort_values("sensor_type", key=lambda s: s.map(prio))
    iv = sorted(sum([json.loads(i) for i in g.intervals], []))
    rows.append(dict(site="+".join(g.sta), lat=float(g.lat.mean()), lon=float(g.lon.mean()),
                     sensor_type=g.sensor_type.iloc[0], n_members=len(g),
                     members=";".join(g.key), nets=";".join(sorted(set(g.net))),
                     any_roster183=bool(g.in_roster183.any()),
                     any_temporary=bool(g.nets_temporary.any()),
                     is_S15=bool(g.sta.isin(S15).any()), is_base_plus=bool(g.sta.isin(BASE_PLUS).any()),
                     first_3c=min(g.first_3c), last_3c=max(g.last_3c), intervals=json.dumps(iv)))
SI = pd.DataFrame(rows).sort_values(["lat"]).reset_index(drop=True)
SI.to_csv(f"{OUT}/sites.csv", index=False)
print("sites:", len(SI), SI.sensor_type.value_counts().to_dict(),
      "| S15 sites:", int(SI.is_S15.sum()), "| +SYMB/BUCKS:", int(SI.is_base_plus.sum()))


def alive(intervals, y):
    a, b = np.datetime64(f"{y}-01-01"), np.datetime64(YEND[y])
    tot = 0.0
    for s, e in intervals:
        tot += max(0.0, float((min(np.datetime64(e), b) - max(np.datetime64(s), a)) / np.timedelta64(1, "D")))
    return tot >= MIN_DAYS


IV = [json.loads(i) for i in SI.intervals]
ALIVE = {y: np.array([alive(v, y) for v in IV]) for y in YEARS}
ACC_IV = [json.loads(i) for i in ACC.intervals]
ACC_ALIVE = {y: np.array([alive(v, y) for v in ACC_IV]) for y in YEARS}
for y in YEARS:
    print(f"  {y}: {ALIVE[y].sum()} sites alive "
          f"(S15 {int((ALIVE[y] & SI.is_S15.values).sum())}, "
          f"+SYMB/BUCKS {int((ALIVE[y] & SI.is_base_plus.values).sum())}), "
          f"{ACC_ALIVE[y].sum()} strong-motion/other stations")

NETWORKS = {"all": np.ones(len(SI), bool), "S15": SI.is_S15.values, "S15+SYMB+BUCKS": SI.is_base_plus.values}


def stats(clat, clon, cdep, mask):
    S = SI[mask]
    if len(S) == 0:
        return dict(n45=0, n55=0, n60=0, nearest_km=np.nan, nearest_site="", gap_deg=360.0,
                    n45_surface_bb=0, n45_borehole=0, n45_short_period=0, sites45="")
    d = hyp_km(clat, clon, cdep, S.lat.values, S.lon.values)
    k45, k60 = d <= 45, d <= 60
    t = S.sensor_type.values
    j = int(np.argmin(d))
    return dict(n45=int(k45.sum()), n55=int((d <= 55).sum()), n60=int(k60.sum()),
                nearest_km=round(float(d[j]), 2), nearest_site=str(S.site.values[j]),
                gap_deg=round(az_gap(azimuth(clat, clon, S.lat.values[k60], S.lon.values[k60])), 1),
                n45_surface_bb=int((k45 & (t == "surface_broadband")).sum()),
                n45_borehole=int((k45 & (t == "borehole")).sum()),
                n45_short_period=int((k45 & (t == "surface_short_period")).sum()),
                sites45="+".join(S.site.values[k45][np.argsort(d[k45])][:8]))


# ---------------- (1) Bostock families ----------------
FAM = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv")
frows = []
for net, base in NETWORKS.items():
    for y in YEARS:
        m = base & ALIVE[y]
        for f in FAM.itertuples():
            frows.append(dict(network=net, year=y, bkey=int(f.bkey), lat=f.lat, lon=f.lon,
                              depth_km=f.dep, n_events=int(f.n_all),
                              **stats(f.lat, f.lon, f.dep, m)))
FG = pd.DataFrame(frows)
FG.to_csv(f"{OUT}/family_geometry_by_year.csv", index=False)

# ---------------- (2) Slab2 cells over the whole box ----------------
sl = pd.read_csv(f"{ROOT}/data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]
sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
_L = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
_N = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
def slab_depth(la, lo):
    xy = np.column_stack([np.asarray(lo, float), np.asarray(la, float)])
    d = _L(xy)
    return np.clip(np.where(np.isfinite(d), d, _N(xy)), 12.0, 55.0)

crows = []
for la in np.arange(BOX[0], BOX[1] + 1e-9, STEP_KM / R_E):
    lons = np.arange(BOX[2], BOX[3] + 1e-9, STEP_KM / (R_E * np.cos(np.deg2rad(la))))
    dep = slab_depth(np.full_like(lons, la), lons)
    k = (dep >= DMIN) & (dep <= DMAX)
    for lo, dd in zip(lons[k], dep[k]):
        crows.append(dict(lat=round(float(la), 4), lon=round(float(lo), 4), depth_km=round(float(dd), 2)))
CELLS = pd.DataFrame(crows)
CELLS.insert(0, "cell_id", [f"C{i:04d}" for i in range(len(CELLS))])
TR = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
TR = TR[TR.time >= "2019-01-01"]
proj = lambda la, lo: np.column_stack([np.asarray(lo) * R_E * np.cos(np.deg2rad(np.asarray(la))),
                                       np.asarray(la) * R_E])
tree = cKDTree(proj(TR.lat.values, TR.lon.values))
CELLS["tremor_n20"] = [len(i) for i in tree.query_ball_point(proj(CELLS.lat.values, CELLS.lon.values), 20.0)]
CELLS["tremor_w"] = CELLS.tremor_n20 / max(CELLS.tremor_n20.max(), 1)
CELLS["tremor_active"] = CELLS.tremor_n20 >= 100
print("cells:", len(CELLS), "tremor-active:", int(CELLS.tremor_active.sum()),
      "| PNSN 2019+ epicentres in region:", len(TR))
crows = []
for net, base in NETWORKS.items():
    for y in YEARS:
        m = base & ALIVE[y]
        nacc = ACC_ALIVE[y]
        for c in CELLS.itertuples():
            s = stats(c.lat, c.lon, c.depth_km, m)
            if net == "all":
                d = hyp_km(c.lat, c.lon, c.depth_km, ACC[nacc].lat.values, ACC[nacc].lon.values)
                s["n45_strongmotion"] = int((d <= 45).sum())
            crows.append(dict(network=net, year=y, cell_id=c.cell_id, lat=c.lat, lon=c.lon,
                              depth_km=c.depth_km, tremor_n20=c.tremor_n20,
                              tremor_w=round(c.tremor_w, 4), tremor_active=c.tremor_active, **s))
CG = pd.DataFrame(crows)
CG.to_csv(f"{OUT}/cell_geometry_by_year.csv", index=False)

# ---------------- (3) ranking of sites outside S15+SYMB/BUCKS ----------------
Y = 2025
base_mask = SI.is_base_plus.values & ALIVE[Y]
cand = np.where((~SI.is_base_plus.values) & ALIVE[Y])[0]
FD = np.array([hyp_km(FAM.lat.values, FAM.lon.values, FAM.dep.values, SI.lat[i], SI.lon[i])
               for i in range(len(SI))])                      # site x family
CD = np.array([hyp_km(CELLS.lat.values, CELLS.lon.values, CELLS.depth_km.values, SI.lat[i], SI.lon[i])
               for i in range(len(SI))])                      # site x cell
W = CELLS.tremor_w.values * CELLS.tremor_active.values
PUG = (CELLS.lat.values >= PUGET_BAND[0]) & (CELLS.lat.values < PUGET_BAND[1])
base_f45 = (FD[base_mask] <= 45).sum(0); base_f55 = (FD[base_mask] <= 55).sum(0)
base_c45 = (CD[base_mask] <= 45).sum(0)


def greedy(D, thr, target, weight, pool, base_counts, nsel=15):
    """Add sites one at a time, each time the one that lifts the most (weighted) units to >= target
    sites within thr km."""
    cnt = base_counts.copy(); chosen = []; avail = list(pool)
    for _ in range(nsel):
        best, bg = None, 0.0
        for i in avail:
            g = float(weight[(cnt < target) & ((cnt + (D[i] <= thr)) >= target)].sum())
            if g > bg:
                best, bg = i, g
        if best is None:
            break
        cnt = cnt + (D[best] <= thr); chosen.append((best, round(bg, 2))); avail.remove(best)
    return chosen


gn = greedy(FD, 45.0, 3, np.ones(len(FAM)), cand, base_f45.astype(int))
gp = greedy(CD, 45.0, 3, W * PUG, cand, base_c45.astype(int))
grank_n = {i: (r + 1, g) for r, (i, g) in enumerate(gn)}
grank_p = {i: (r + 1, g) for r, (i, g) in enumerate(gp)}

rows = []
for i in cand:
    f45, f55 = FD[i] <= 45, FD[i] <= 55
    c45 = CD[i] <= 45
    rows.append(dict(site=SI.site[i], members=SI.members[i], nets=SI.nets[i],
                     lat=round(SI.lat[i], 4), lon=round(SI.lon[i], 4),
                     sensor_type=SI.sensor_type[i], in_roster183=SI.any_roster183[i],
                     temporary_net=SI.any_temporary[i],
                     first_3c=SI.first_3c[i], last_3c=SI.last_3c[i],
                     years_alive="/".join(str(y) for y in YEARS if ALIVE[y][i]),
                     fam45=int(f45.sum()), fam55=int(f55.sum()),
                     fam45_uncovered=int((f45 & (base_f45 == 0)).sum()),
                     fam45_lift_to3=int((f45 & (base_f45 == 2)).sum()),
                     fam55_uncovered=int((f55 & (base_f55 == 0)).sum()),
                     cells45=int(c45.sum()),
                     cells45_tw=round(float(W[c45].sum()), 2),
                     cells45_tw_puget=round(float((W * PUG)[c45].sum()), 2),
                     cells45_uncovered_tw=round(float(W[c45 & (base_c45 == 0)].sum()), 2),
                     greedy_rank_north=grank_n.get(i, (0, 0))[0],
                     greedy_gain_north=grank_n.get(i, (0, 0))[1],
                     greedy_rank_puget=grank_p.get(i, (0, 0))[0],
                     greedy_gain_puget=grank_p.get(i, (0, 0))[1]))
RK = pd.DataFrame(rows)
RK["score_north"] = RK.fam45 + 2 * RK.fam45_uncovered
RK = RK.sort_values(["greedy_rank_north", "fam45"], ascending=[True, False])
RK.to_csv(f"{OUT}/station_ranking.csv", index=False)

# ---------------- (4) Bostock's own 2003-2005 template array, for scale ----------------
BM = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia/metadata.csv",
                 usecols=["station_network_code", "station_code", "station_latitude_deg",
                          "station_longitude_deg"])
BA = (BM.groupby(["station_network_code", "station_code"])
        .agg(lat=("station_latitude_deg", "first"), lon=("station_longitude_deg", "first"))
        .reset_index().rename(columns={"station_code": "sta"}))
BA["permanent"] = BA.sta.isin(["PGC", "VGZ", "SNB", "LZB", "PFB", "YOUB"])
BA.to_csv(f"{OUT}/bostock_array_stations.csv", index=False)
bxy = np.column_stack([BA.lon * R_E * np.cos(np.deg2rad(BA.lat)), BA.lat * R_E])
BA["cl"] = fcluster(linkage(pdist(bxy), method="single"), t=SITE_MERGE_KM, criterion="distance")
BS = BA.groupby("cl").agg(lat=("lat", "mean"), lon=("lon", "mean")).reset_index()
BD = np.array([hyp_km(FAM.lat.values, FAM.lon.values, FAM.dep.values, BS.lat[i], BS.lon[i])
               for i in range(len(BS))])
pq = SI.site.isin(["JRDRB", "LHLYB", "DAOB", "ALBHB"]).values
SCEN = {"Bostock_2003_2005_array": (BD, np.ones(len(BS), bool)),
        "S15": (D_ := FD, SI.is_S15.values & ALIVE[2025]),
        "S15+SYMB+BUCKS": (FD, base_mask),
        "S15+SYMB+BUCKS+PQ4": (FD, base_mask | (pq & ALIVE[2025]))}
rows, agg = [], []
for name, (Dm, m) in SCEN.items():
    c45 = (Dm[m] <= 45).sum(0); c55 = (Dm[m] <= 55).sum(0); near = Dm[m].min(0)
    for j, b in enumerate(FAM.bkey):
        rows.append(dict(network=name, bkey=int(b), n45=int(c45[j]), n55=int(c55[j]),
                         nearest_km=round(float(near[j]), 2)))
    agg.append(dict(network=name, n_sites=int(m.sum()), n45_med=float(np.median(c45)),
                    n45_mean=round(float(c45.mean()), 2), n45_p90=float(np.percentile(c45, 90)),
                    frac_ge1=round(float((c45 >= 1).mean()), 3), frac_ge2=round(float((c45 >= 2).mean()), 3),
                    frac_ge3=round(float((c45 >= 3).mean()), 3), n55_med=float(np.median(c55)),
                    nearest_med=round(float(np.median(near)), 2)))
pd.DataFrame(rows).to_csv(f"{OUT}/bostock_array_comparison.csv", index=False)
AG = pd.DataFrame(agg)
AG.to_csv(f"{OUT}/bostock_array_comparison_summary.csv", index=False)
print("\n== Bostock array vs our sets (per family, sites <=45 km):")
print(AG.to_string(index=False))
print("\nTOP north (greedy, families to >=3 sites <=45 km):")
print(RK[RK.greedy_rank_north > 0].sort_values("greedy_rank_north")
        [["greedy_rank_north", "greedy_gain_north", "site", "nets", "sensor_type", "fam45",
          "fam45_uncovered", "years_alive"]].head(12).to_string(index=False))
print("\nTOP Puget/Olympics (greedy, tremor-weighted cells to >=3 sites <=45 km):")
print(RK[RK.greedy_rank_puget > 0].sort_values("greedy_rank_puget")
        [["greedy_rank_puget", "greedy_gain_puget", "site", "nets", "sensor_type", "cells45",
          "cells45_tw_puget", "years_alive"]].head(12).to_string(index=False))
