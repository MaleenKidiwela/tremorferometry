#!/usr/bin/env python
"""STEP 1 of the resolution-through-time study (Merlin-corrected recipe).
Assemble the certified (cell, station) pair catalog from the FINAL roster:
  - one daily_dvv_<TAG>_Z_2to4.csv per station (exclude *rev*, *MIRROR*)
  - reliable families only: <stem>_causality_cert.csv (reliable==True), fallback <stem>_fwd_vs_rev_coda.csv (ratio>1.5)
  - cc_max>0.7; snap families to 0.10 deg cells; a (cell,station) datum = station's reliable families in that cell
  - Slab2 interface depth per cell (clip 12-55 km); station coords; distinct SITES (cluster stations <2 km)
  - per-(cell,station): sigma_daily = robust MAD of the 90-day-detrended daily series (RAW -> upper bound on noise;
    only 7 boreholes have a MIRROR product, not fleet-wide), and per-month observation-day counts.
Outputs (fault_tomography/inversion/res_catalog/):
  cells.csv  pairs.csv  pair_months.parquet   + prints Merlin's accept checks.
"""
import os, glob, json
import numpy as np, pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

CC_MIN = 0.70
CELL_DEG = float(os.environ.get("CELL_DEG", "0.10"))
# monthly dv/v aggregation: "mean" = average of the ~30 daily rolling-stack values in the month (current default,
# double-smooths + adjacent months correlated); "sample" = take the daily value at the month END = one trailing
# 30-day stack per month -> non-overlapping, near-independent months.
MONTHLY_MODE = os.environ.get("MONTHLY_MODE", "mean")
SITE_KM = 2.0
OUT = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
os.makedirs(OUT, exist_ok=True)

# ---------- station coordinates ----------
fleet = pd.read_csv("data/broadband_fleet_order.csv")
FLEET_XY = {r.sta: (r.lat, r.lon) for _, r in fleet.iterrows()}
cand = pd.read_csv("data/candidate_stations_post2020.csv")
CAND_XY = {r.sta: (r.lat, r.lon) for _, r in cand.iterrows()}
ANCHOR = {"pgc": (48.6498, -123.4521), "shb": (49.5990, -123.8770), "clrs": (48.8203, -124.1309)}


def coords_for(tag, stem):
    for pre, ll in ANCHOR.items():
        if stem.startswith(pre):
            return ll
    if tag in FLEET_XY:
        return FLEET_XY[tag]
    up = stem.upper()
    if up in FLEET_XY:
        return FLEET_XY[up]
    if up in CAND_XY:
        return CAND_XY[up]
    return None


# ---------- reliable families per station ----------
def cert_fams(stem):
    f1 = f"data/{stem}_causality_cert.csv"
    if os.path.exists(f1):
        c = pd.read_csv(f1)
        return set(c[c.reliable].fam)
    f2 = f"data/{stem}_fwd_vs_rev_coda.csv"
    if os.path.exists(f2):
        c = pd.read_csv(f2)
        return set(c[c.ratio > 1.5].fam)
    return set()


# ---------- Slab2 interface depth interpolator ----------
sl = pd.read_csv("data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]
sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
SL_PTS = sl[["lon", "lat"]].values
SL_D = sl["depth"].values
_SLL = LinearNDInterpolator(SL_PTS, SL_D)      # built once (slab_depth is now called per cell)
_SLN = NearestNDInterpolator(SL_PTS, SL_D)


def slab_depth(clat, clon):
    xy = np.column_stack([np.asarray(clon, float), np.asarray(clat, float)])
    d = _SLL(xy); dn = _SLN(xy)
    d = np.where(np.isfinite(d), d, dn)
    return np.clip(d, 12.0, 55.0)


def famll(pid):
    a = pid.split("__")[0].split("_")
    return float(a[0]), float(a[1])


def snap(v):
    return round(round(v / CELL_DEG) * CELL_DEG, 2)


# ---------- iterate daily files ----------
daily_files = sorted(glob.glob("data/daily_dvv_*_Z_2to4.csv"))
pair_rows, month_rows = [], []
skipped_nocoord, skipped_nocert, used_fwdrev = [], [], []
n_stations = 0
for fp in daily_files:
    tag = os.path.basename(fp).split("daily_dvv_")[1].split("_Z_2to4")[0]
    if "rev" in tag.lower() or "MIRROR" in tag:
        continue
    stem = tag.lower().replace("p90f40", "")
    fams = cert_fams(stem)
    if not fams:
        skipped_nocert.append(tag); continue
    if not os.path.exists(f"data/{stem}_causality_cert.csv"):
        used_fwdrev.append(tag)
    ll = coords_for(tag, stem)
    if ll is None:
        skipped_nocoord.append(tag); continue
    sla, slo = ll
    d = pd.read_csv(fp, usecols=["patch", "date", "dvv", "cc_max"])
    d = d[(d.cc_max > CC_MIN) & (d.patch.isin(fams))]
    if not len(d):
        continue
    d["date"] = pd.to_datetime(d["date"], errors="coerce")
    d = d.dropna(subset=["date"])
    ll2 = d.patch.map(famll)
    d["flat"] = [a for a, _ in ll2]          # raw LFE family position (the actual source location)
    d["flon"] = [b for _, b in ll2]
    d["clat"] = [snap(a) for a, _ in ll2]
    d["clon"] = [snap(b) for _, b in ll2]
    d["cell"] = d.clat.map("{:.2f}".format) + "_" + d.clon.map("{:.2f}".format)
    n_stations += 1
    # (cell, day) mean dv/v across the station's families in that cell -> the (cell,station) daily series
    for cell, g in d.groupby("cell"):
        clat, clon = g.clat.iloc[0], g.clon.iloc[0]
        ser = g.groupby(g.date.dt.floor("D")).dvv.mean().sort_index() * 100.0  # percent
        if len(ser) < 5:
            continue
        # 90-day rolling-median detrend, then robust MAD
        det = ser - ser.rolling(90, center=True, min_periods=10).median()
        det = det.dropna()
        sig = float(1.4826 * np.median(np.abs(det - np.median(det)))) if len(det) >= 10 else float(ser.std())
        if not np.isfinite(sig) or sig <= 0:
            sig = float(np.nanstd(ser.values)) or 0.2
        # SOURCE location = centroid of THIS station's distinct families in the cell (not the snapped grid center).
        fll = g.groupby("patch")[["flat", "flon"]].first()
        src_lat = float(fll.flat.median()); src_lon = float(fll.flon.median())
        fdep = slab_depth(fll.flat.values, fll.flon.values)
        src_dep = float(np.median(fdep))                       # family-median interface depth (kernel source depth)
        pair_rows.append(dict(tag=tag, stem=stem, cell=cell, cell_lat=clat, cell_lon=clon,
                              src_lat=round(src_lat, 4), src_lon=round(src_lon, 4), src_dep=round(src_dep, 2),
                              n_fam=int(len(fll)), sta_lat=sla, sta_lon=slo, n_days=int(len(ser)), sigma_daily=round(sig, 4)))
        for per, sub in ser.groupby(ser.index.to_period("M")):
            val = float(sub.iloc[-1]) if MONTHLY_MODE == "sample" else float(sub.mean())  # last day = trailing-30d stack
            month_rows.append((tag, cell, str(per), int(sub.size), round(val, 4)))

pairs = pd.DataFrame(pair_rows)
pm = pd.DataFrame(month_rows, columns=["tag", "cell", "ym", "n_days", "dvv_month"])
print(f"stations used: {n_stations} | (cell,station) pairs: {len(pairs)} | pair-months: {len(pm)}")

# ---------- distinct SITES (cluster stations within SITE_KM) ----------
sta = pairs[["tag", "sta_lat", "sta_lon"]].drop_duplicates("tag").reset_index(drop=True)
lat0, lon0 = sta.sta_lat.mean(), sta.sta_lon.mean()
sx = (sta.sta_lon - lon0) * 111 * np.cos(np.radians(lat0))
sy = (sta.sta_lat - lat0) * 111
site_id = np.arange(len(sta))
for i in range(len(sta)):
    for j in range(i):
        if (sx[i] - sx[j]) ** 2 + (sy[i] - sy[j]) ** 2 < SITE_KM ** 2:
            site_id[site_id == site_id[i]] = site_id[j]
sta["site"] = site_id
tag2site = dict(zip(sta.tag, sta.site))
pairs["site"] = pairs.tag.map(tag2site)
print(f"distinct stations: {sta.tag.nunique()} -> distinct SITES (<{SITE_KM} km): {sta.site.nunique()}")

# ---------- cells table: FAMILY-CENTROID position + FAMILY-MEDIAN interface depth ----------
# Each cell is represented by where its LFE families actually sit (median over the per-station family
# centroids), NOT the snapped geometric grid center. depth_km = family-median Slab2 depth -> used for both
# the kernel source depth AND the DEEP (>30 km) classification. frac_deep/depth_mixed flag near-boundary cells.
cells = pairs.groupby("cell").agg(
    cell_lat=("cell_lat", "first"), cell_lon=("cell_lon", "first"),     # snapped grid center (display mesh only)
    src_lat=("src_lat", "median"), src_lon=("src_lon", "median"),       # family centroid (kernel horizontal source)
    depth_km=("src_dep", "median"),                                     # family-median depth (kernel + classification)
    frac_deep=("src_dep", lambda s: float((s > 30).mean())),
    n_pairs=("tag", "size"), n_sites=("site", "nunique"), n_fam=("n_fam", "sum")).reset_index()
cells["depth_center"] = slab_depth(cells.cell_lat.values, cells.cell_lon.values)  # DIAGNOSTIC: old center-based depth
cells["depth_mixed"] = (cells.frac_deep > 0.3) & (cells.frac_deep < 0.7)          # ambiguous deep/shallow (fuzzy ~30 km)
print(f"cells: {len(cells)} | with >=3 sites: {(cells.n_sites>=3).sum()} | "
      f"deep(>30km) family-median: {(cells.depth_km>30).sum()} (center-based was {(cells.depth_center>30).sum()}) | "
      f"depth-mixed near boundary: {int(cells.depth_mixed.sum())}")

pairs.to_csv(f"{OUT}/pairs.csv", index=False)
cells.to_csv(f"{OUT}/cells.csv", index=False)
pm.to_parquet(f"{OUT}/pair_months.parquet", index=False)

# ---------- Merlin accept checks ----------
print("\n=== ACCEPT CHECKS ===")
print(f"[rev excluded]        {'PASS' if not any('rev' in t.lower() for t in pairs.tag.unique()) else 'FAIL'}")
print(f"[all pairs have coord] {'PASS' if pairs[['sta_lat','sta_lon']].notna().all().all() else 'FAIL'}")
print(f"[sigma_daily median]   {pairs.sigma_daily.median():.3f}%  (expect ~0.2%)")
print(f"[fwd_vs_rev fallbacks]  {used_fwdrev}")
print(f"[skipped no-coord]      {skipped_nocoord}")
print(f"[skipped no-cert count] {len(skipped_nocert)}")
print(f"[date span]             {pm.ym.min()} -> {pm.ym.max()}")
json.dump(dict(stations=int(n_stations), pairs=int(len(pairs)), cells=int(len(cells)),
               sites=int(sta.site.nunique()), cc_min=CC_MIN, cell_deg=CELL_DEG),
          open(f"{OUT}/meta.json", "w"), indent=2)
print(f"\nwrote {OUT}/pairs.csv, cells.csv, pair_months.parquet")
