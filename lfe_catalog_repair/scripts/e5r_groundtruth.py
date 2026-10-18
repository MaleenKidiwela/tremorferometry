#!/usr/bin/env python
"""E5 step 5 STEP B -- Bostock's Sept-2012 events as ground truth for the from-scratch replication.

* events of results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv in 2012-08-15..2012-10-05 (the download window);
* Bostock family = unique (lat, lon, dep); mapped to the SeisBench family id by NEAREST NEIGHBOUR in (x, y, z) km, exactly as
  e2_groundtruth.bostock_families() did for 2010 (join distance reported; > 1 km -> unmatched, family kept with fid = 1000 + k);
* per-family time-base offset: provisional = the E3 2010 offset of the same family (accepted ones), else the median accepted
  offset; the 2012 offsets are re-derived on the BUILD half at PGC by e5r_offsets.py and this table is then REBUILT with them
  (--offsets csv) so the split is defined on t_origin = DB time + offset, as E4/E5 did;
* hold-out = ALTERNATE 2-h blocks of t_origin with a 60-s guard at every edge (never by event parity).
Outputs results/e5r_gt_bostock2012.parquet (time [ns], bkey, fam_id, lat, lon, dep, mag, offset_s, t_origin, half, guard, blk)
        results/e5r_families_2012.csv (per family: bkey, fam_id, lat, lon, dep, n_all, n_build, n_test, n_guard, mag_med, join_km)
"""
import argparse, numpy as np, pandas as pd
from scipy.spatial import cKDTree
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
R = 111.195

ap = argparse.ArgumentParser()
ap.add_argument("--offsets", default="", help="csv with bkey, offset_s, accept (2012 re-derived); default = E3 2010 offsets")
ap.add_argument("--t0", default="2012-08-15"); ap.add_argument("--t1", default="2012-10-06")
ap.add_argument("--block-s", type=float, default=7200.0); ap.add_argument("--guard-s", type=float, default=60.0)
a = ap.parse_args()

bo = pd.read_csv(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
t = pd.to_datetime(dict(year=bo.year, month=bo.month, day=bo.day, hour=bo.hour, minute=bo["min"],
                        second=bo.sec.astype(int))) + pd.to_timedelta(bo.sec % 1.0, unit="s")
bo["time"] = t.astype("datetime64[ns]")
bo = bo[(bo.time >= a.t0) & (bo.time < a.t1)].copy()
key = bo.groupby(["lat", "lon", "dep"]).size().reset_index(name="n_all").reset_index().rename(columns={"index": "bkey"})
# nearest SeisBench family
md = pd.read_csv(f"{SD}/metadata.csv")
fam = md.groupby("source_family_id")[["source_latitude_deg", "source_longitude_deg", "source_depth_km"]].first().reset_index()
cs = np.cos(np.radians(48.6))
A = np.c_[fam.source_latitude_deg * R, fam.source_longitude_deg * R * cs, fam.source_depth_km]
B = np.c_[key.lat * R, key.lon * R * cs, key.dep]
d, i = cKDTree(A).query(B)
key["sb_fam_id"] = fam.source_family_id.to_numpy()[i].astype(float)
key["join_km"] = d
bad = key.join_km > 1.0
key.loc[bad, "sb_fam_id"] = np.nan
key["fam_id"] = key.bkey.astype(float)   # the correlator keys families by fam_id = bkey (from scratch)
dup = key.sb_fam_id.duplicated(keep=False) & ~bad
if dup.any():
    print("[warn] several Bostock keys map to one SeisBench family:\n", key[dup])
print(f"[fam] {len(key)} Bostock 2012 families; {int((~bad).sum())} matched to SeisBench ids within 1 km "
      f"(join median {key.join_km.median():.3f} km, max matched {key.join_km[~bad].max():.3f}); unmatched: "
      f"{key.loc[bad, ['lat','lon','dep','n_all']].to_dict('records')}")
bo = bo.merge(key[["lat", "lon", "dep", "bkey", "fam_id", "sb_fam_id", "join_km"]], on=["lat", "lon", "dep"])
# offsets
fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
acc = fo[fo.accept.astype(bool)]
off2010 = dict(zip(acc.fam_id.astype(float), acc.offset_s.astype(float)))
med2010 = float(np.median(list(off2010.values())))
key["offset_2010"] = key.sb_fam_id.map(off2010)
if a.offsets:
    o = pd.read_csv(a.offsets)
    o = o[o.accept.astype(bool)]
    omap = dict(zip(o.bkey.astype(int), o.offset_s.astype(float)))
    key["offset_s"] = key.bkey.map(omap)
    key["offset_src"] = np.where(key.offset_s.notna(), "2012_build", np.where(key.offset_2010.notna(), "2010", "median"))
    key["offset_s"] = key.offset_s.fillna(key.offset_2010).fillna(med2010)
else:
    key["offset_src"] = np.where(key.offset_2010.notna(), "2010", "median2010")
    key["offset_s"] = key.offset_2010.fillna(med2010)
bo = bo.merge(key[["bkey", "offset_s", "offset_src"]], on="bkey")
bo["t_origin"] = bo.time.astype("datetime64[ns]").astype("int64") / 1e9 + bo.offset_s
blk = np.floor(bo.t_origin / a.block_s).astype("int64")
pos = bo.t_origin - blk * a.block_s
bo["blk"] = blk
bo["half"] = np.where(blk % 2 == 0, "build", "test")
bo["guard"] = (pos < a.guard_s) | (pos > a.block_s - a.guard_s)
bo = bo.sort_values("time").reset_index(drop=True)
cols = ["time", "bkey", "fam_id", "sb_fam_id", "lat", "lon", "dep", "mag", "offset_s", "offset_src", "t_origin", "blk", "half", "guard", "join_km"]
bo[cols].to_parquet(f"{RES}/e5r_gt_bostock2012.parquet", index=False)
g = bo.groupby("bkey")
F = key.merge(pd.DataFrame(dict(
    n_build=g.apply(lambda x: int(((x.half == "build") & ~x.guard).sum())),
    n_test=g.apply(lambda x: int(((x.half == "test") & ~x.guard).sum())),
    n_guard=g.guard.sum().astype(int), mag_med=g.mag.median())).reset_index(), on="bkey")
F.to_csv(f"{RES}/e5r_families_2012.csv", index=False)
print(f"[evt] {len(bo)} events {bo.time.min()} .. {bo.time.max()}; build {int(((bo.half=='build')&~bo.guard).sum())} / "
      f"test {int(((bo.half=='test')&~bo.guard).sum())} / guard {int(bo.guard.sum())}; offsets: {key.offset_src.value_counts().to_dict()}")
for k in (30, 100, 200):
    print(f"[fam] families with >= {k} build-half events: {int((F.n_build >= k).sum())} of {len(F)}")
print(F.n_build.describe().round(1).to_string())
