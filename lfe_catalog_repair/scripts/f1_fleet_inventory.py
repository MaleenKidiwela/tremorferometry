#!/usr/bin/env python
"""F1 — FLEET INVENTORY + association-region design for a fleet-wide LFE picker run.

No waveforms are downloaded.  One station-level FDSN metadata request per network
(cached as StationXML under $F1_CACHE) supplies channels / sample rates / epochs.

Stages
  A  roster      : every station this project has ever used (broadband fleet order,
                   post-2020 candidates incl. PB boreholes, 3 anchors, + everything with a
                   daily_dvv product or an inversion pair) plus the *extra* 3-C FDSN stations
                   in-region that are NOT on the project list (they cost nothing to add and
                   they are what buys association).
  B  metadata    : per (net,sta) channel epochs -> per station-year ZNE-at->=20-Hz availability
                   (2009-2026; the running year is truncated at TODAY so volumes are not overstated).
  C  source grid : 0.25 deg cells on the Slab2 interface, 30-45 km depth, 40-50 N (the ETS band);
                   per cell-year station counts at 60/100/150 km, nearest distance, azimuthal
                   gap within 100 km, borehole fraction.
  D  regions     : 7 along-strike association regions; per region-year median 3-C station count
                   and the fraction of cells with >=3 / >=4 stations ("locatable" fraction).
  E  cost        : download volume (native + decimated), GPU picker minutes, associate/locate CPU.
  F  ground truth: which regions have Lin / Bostock / Ducellier / Shelly catalogs for
                   validation and fine-tuning.

Outputs (lfe_catalog_repair/):
  results/f1_stations.csv  f1_station_years.csv  f1_cell_density.csv  f1_regions.csv
  results/f1_region_stations.csv  f1_cost_estimate.csv  f1_groundtruth.csv
  figures/f1_fleet_inventory.png
"""
from __future__ import annotations

import glob
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_catalog_repair/results"
FIG = f"{ROOT}/lfe_catalog_repair/figures"
CACHE = os.environ.get(
    "F1_CACHE",
    "/tmp/claude-1000/-home-jovyan-tremorferometry/dbe5feb4-daae-48b0-9ca9-47aa075fc3b6/scratchpad/f1_stationxml",
)
os.makedirs(CACHE, exist_ok=True)

NOW = pd.Timestamp("2026-09-10")         # run date: the open-ended epochs stop here
YEARS = list(range(2009, 2027))          # 2009..2026 inclusive (project fleet starts 2009-01-01)
MIN_FS = 20.0                            # LFEDetect input rate; ZNE must exist at >= this
BBOX = dict(minlatitude=39.4, maxlatitude=50.6, minlongitude=-128.5, maxlongitude=-118.8)
CELL_DEG = 0.25
DEPTH_BAND = (30.0, 45.0)                # ETS / LFE depth band on the interface
LAT_BAND = (40.0, 50.0)
RADII = (60.0, 100.0, 150.0)
BOREHOLE_DEPTH_M = 20.0                  # channel emplacement depth -> borehole

# networks: the project's own + the regional networks that association can lean on
NETS_PROJECT = ["UW", "UO", "CC", "CN", "BK", "NC", "PB", "US"]
NETS_EXTRA = ["PO", "C8", "NV", "TA", "IU", "OO"]
NCEDC_NETS = {"BK", "NC"}

# along-strike association regions (lat_min, lat_max, name)
REGIONS = [
    ("R1_N_VancIsl", 48.70, 50.00),
    ("R2_S_VancIsl_JdF", 47.90, 48.70),
    ("R3_Puget_Olympic", 46.60, 47.90),
    ("R4_N_Oregon", 45.20, 46.60),
    ("R5_C_Oregon", 43.60, 45.20),
    ("R6_S_Oregon", 42.00, 43.60),
    ("R7_NorCal_MTJ", 40.00, 42.00),
]

# ---------------------------------------------------------------- geometry helpers
R_EARTH = 6371.0


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = (np.radians(np.asarray(x, float)) for x in (lat1, lon1, lat2, lon2))
    d = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(d, 0, 1)))


def azimuth_deg(lat1, lon1, lat2, lon2):
    """Azimuth from point 1 (cell) to point 2 (station), degrees clockwise from north."""
    la1, lo1, la2, lo2 = (np.radians(np.asarray(x, float)) for x in (lat1, lon1, lat2, lon2))
    y = np.sin(lo2 - lo1) * np.cos(la2)
    x = np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    return np.degrees(np.arctan2(y, x)) % 360.0


def max_az_gap(az):
    """Largest gap (deg) between consecutive station azimuths seen from the cell."""
    az = np.sort(np.mod(np.asarray(az, float), 360.0))
    if az.size == 0:
        return 360.0
    if az.size == 1:
        return 360.0
    d = np.diff(az)
    return float(max(d.max(), 360.0 - (az[-1] - az[0])))


# ================================================================ STAGE A — roster
def build_roster():
    rows = []

    fleet = pd.read_csv(f"{ROOT}/data/broadband_fleet_order.csv")
    for r in fleet.itertuples():
        rows.append(dict(net=r.net, sta=r.sta, lat_tab=r.lat, lon_tab=r.lon,
                         role="fleet_broadband", region_tab=r.region,
                         fleet_order=int(r.order), tab_chans=r.chans))

    cand = pd.read_csv(f"{ROOT}/data/candidate_stations_post2020.csv")
    have = {(r["net"], r["sta"]) for r in rows}
    for r in cand.itertuples():
        if (r.net, r.sta) in have:
            continue
        rows.append(dict(net=r.net, sta=r.sta, lat_tab=r.lat, lon_tab=r.lon,
                         role="borehole" if r.net == "PB" else "candidate",
                         region_tab=np.nan, fleet_order=np.nan, tab_chans=r.chans))
        have.add((r.net, r.sta))

    for sta, (la, lo) in {"PGC": (48.6498, -123.4521), "SHB": (49.5990, -123.8770),
                          "CLRS": (48.8203, -124.1309)}.items():
        if ("CN", sta) not in have:
            rows.append(dict(net="CN", sta=sta, lat_tab=la, lon_tab=lo, role="anchor",
                             region_tab=np.nan, fleet_order=np.nan, tab_chans="BH?/HH?"))
            have.add(("CN", sta))
    df = pd.DataFrame(rows)
    df.loc[df.sta.isin(["PGC", "SHB", "CLRS"]), "role"] = "anchor"
    return df


def norm_tag(tag):
    """daily_dvv_<TAG> / pairs.csv tag -> bare station code."""
    t = tag
    for suf in ("rev", "p90f40"):
        while t.lower().endswith(suf):
            t = t[: -len(suf)]
    for suf in ("bhz", "hhz", "ehz", "new"):
        if t.lower().endswith(suf) and len(t) > len(suf):
            t = t[: -len(suf)]
    return t.upper()


def mark_usage(df):
    dvv = sorted({norm_tag(os.path.basename(f).split("daily_dvv_")[1].split("_Z_2to4")[0])
                  for f in glob.glob(f"{ROOT}/data/daily_dvv_*_Z_2to4.csv")})
    pairs = pd.read_csv(f"{ROOT}/fault_tomography/inversion/res_catalog/pairs.csv")
    inv = sorted({norm_tag(t) for t in pairs.tag.unique()})
    cert = sorted({os.path.basename(f).split("_causality_cert.csv")[0].replace("p90f40", "").upper()
                   for f in glob.glob(f"{ROOT}/data/*_causality_cert.csv")})
    df["has_dvv"] = df.sta.str.upper().isin(dvv)
    df["in_inversion"] = df.sta.str.upper().isin(inv)
    df["has_cert"] = df.sta.str.upper().isin(cert)
    known = set(df.sta.str.upper())
    missing = [s for s in set(dvv) | set(inv) if s not in known]
    return df, sorted(missing), len(dvv), len(inv)


# ============================================================= STAGE B — metadata
def fetch_network(net):
    """One station-level (level=channel) FDSN request per network, cached to StationXML."""
    from obspy import read_inventory
    from obspy.clients.fdsn import Client

    fp = f"{CACHE}/{net}.xml"
    if os.path.exists(fp) and os.path.getsize(fp) > 500:
        return read_inventory(fp)
    provs = (["NCEDC"] if net in NCEDC_NETS else ["EARTHSCOPE", "IRIS"])
    for prov in provs:
        try:
            inv = Client(prov, timeout=180).get_stations(
                network=net, station="*", channel="?H?", level="channel", **BBOX)
            inv.write(fp, format="STATIONXML")
            print(f"  [fdsn] {net:3s} via {prov}: "
                  f"{sum(len(n) for n in inv)} stations", flush=True)
            return inv
        except Exception as e:                              # noqa: BLE001
            print(f"  [fdsn] {net:3s} {prov} failed: {type(e).__name__}: {e}", flush=True)
    return None


def channel_records(nets):
    """-> DataFrame of one row per channel epoch (high-gain seismometers only)."""
    recs = []
    for net in nets:
        inv = fetch_network(net)
        if inv is None:
            continue
        for nw in inv:
            for st in nw:
                for ch in st.channels:
                    code = ch.code
                    if len(code) != 3 or code[0] not in "BHES" or code[1] != "H":
                        continue                            # high-gain seismometer bands only
                    if ch.sample_rate is None or ch.sample_rate < 1:
                        continue
                    e = ch.end_date
                    recs.append(dict(
                        net=nw.code, sta=st.code, chan=code, loc=(ch.location_code or "").strip(),
                        band=code[:2], comp=code[2], fs=float(ch.sample_rate),
                        lat=float(ch.latitude), lon=float(ch.longitude),
                        elev=float(ch.elevation), depth_m=float(ch.depth or 0.0),
                        y0=ch.start_date.year, y1=(e.year if e is not None else 2026),
                        t0=ch.start_date.datetime, t1=(e.datetime if e is not None else None)))
    return pd.DataFrame(recs)


HORIZ = ({"N", "E"}, {"1", "2"}, {"X", "Y"})


def station_year_table(ch):
    """For each (net,sta,year): does a full ZNE (or Z12) triplet at >= MIN_FS exist?"""
    import datetime as dt

    out = []
    for (net, sta), g in ch.groupby(["net", "sta"], sort=True):
        for y in YEARS:
            ys, ye = dt.datetime(y, 1, 1), dt.datetime(y + 1, 1, 1)
            act = g[(g.y0 <= y) & (g.y1 >= y)]
            # exact epoch overlap with the calendar year
            keep = []
            for r in act.itertuples():
                # NB pd.NaT is TRUTHY, so an open-ended epoch must be tested with isna()
                end = min(ye, NOW) if pd.isna(r.t1) else min(r.t1, ye, NOW)
                a = max(r.t0, ys)
                b = end
                if b > a:
                    keep.append((r, (b - a).days))
            if not keep:
                out.append(dict(net=net, sta=sta, year=y, has_3c=False, z_only=False,
                                fs=np.nan, fs_lo=np.nan, band="", loc="", chans="",
                                orient="", days=0, is_borehole=False))
                continue
            best = None
            cands = []
            z_any = False
            for (band, loc), idx in defaultdict(list, {
                    (r.band, r.loc): None for r, _ in keep}).items():
                sub = [(r, d) for r, d in keep if r.band == band and r.loc == loc]
                comps = {r.comp for r, _ in sub}
                fs = max(r.fs for r, _ in sub)
                if fs < MIN_FS:
                    continue
                z_only_here = False
                if "Z" in comps:
                    z_only_here = True
                    for hset in HORIZ:
                        if hset <= comps:
                            z_only_here = False
                            orient = "ZNE" if hset == {"N", "E"} else (
                                "Z12" if hset == {"1", "2"} else "ZXY")
                            days = int(np.median([d for r, d in sub if r.comp in hset | {"Z"}]))
                            dep = max(r.depth_m for r, _ in sub)
                            cand = dict(band=band, loc=loc, fs=fs, orient=orient, days=days,
                                        chans=",".join(sorted(r.chan for r, _ in sub)),
                                        depth_m=dep,
                                        score=(2 if band in ("BH", "HH") else 1) * 1000 + days)
                            cands.append(cand)
                            if best is None or cand["score"] > best["score"]:
                                best = cand
                            break
                    z_any = z_any or z_only_here
            if best is None:
                out.append(dict(net=net, sta=sta, year=y, has_3c=False, z_only=bool(z_any),
                                fs=np.nan, fs_lo=np.nan, band="", loc="", chans="", orient="",
                                days=0, is_borehole=False))
            else:
                out.append(dict(net=net, sta=sta, year=y, has_3c=True, z_only=False,
                                fs=best["fs"], fs_lo=float(min(c["fs"] for c in cands)),
                                band=best["band"], loc=best["loc"],
                                chans=best["chans"], orient=best["orient"], days=best["days"],
                                is_borehole=best["depth_m"] >= BOREHOLE_DEPTH_M))
    return pd.DataFrame(out)


def classify(ch_sta, sy_sta):
    """Instrument class for a station, from its channel bands and emplacement depth."""
    bands = set(ch_sta.band[ch_sta.fs >= MIN_FS])
    deep = float(ch_sta.depth_m.max() or 0.0) >= BOREHOLE_DEPTH_M
    bb = bool(bands & {"BH", "HH"})
    if bb and deep:
        return "broadband_borehole"
    if bb:
        return "broadband"
    if bands & {"EH", "SH"}:
        return "borehole_short_period" if deep else "short_period"
    return "other"


# ============================================================ STAGE C — source grid
def slab_interp():
    sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
    from sse_strain.fault import load_fault_ascii
    from sse_strain.megathrust import slab_depth_interpolator
    f = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
    return slab_depth_interpolator(f)


def source_grid():
    itp = slab_interp()
    lats = np.arange(LAT_BAND[0] + CELL_DEG / 2, LAT_BAND[1], CELL_DEG)
    lons = np.arange(-128.0 + CELL_DEG / 2, -119.0, CELL_DEG)
    LO, LA = np.meshgrid(lons, lats)
    D = itp(LO.ravel(), LA.ravel())
    ok = np.isfinite(D) & (D >= DEPTH_BAND[0]) & (D <= DEPTH_BAND[1])
    g = pd.DataFrame(dict(cell_lat=np.round(LA.ravel()[ok], 3),
                          cell_lon=np.round(LO.ravel()[ok], 3),
                          slab_km=np.round(D[ok], 1)))
    g["cell"] = g.cell_lat.map("{:.3f}".format) + "_" + g.cell_lon.map("{:.3f}".format)
    g["region"] = ""
    for name, a, b in REGIONS:
        g.loc[(g.cell_lat >= a) & (g.cell_lat < b), "region"] = name
    return g.reset_index(drop=True)


def cell_density(grid, sta, sy):
    """cell x year station-geometry metrics, for the project roster and for project+extra."""
    slat, slon = sta.lat.values, sta.lon.values
    D = haversine_km(grid.cell_lat.values[:, None], grid.cell_lon.values[:, None],
                     slat[None, :], slon[None, :])
    AZ = azimuth_deg(grid.cell_lat.values[:, None], grid.cell_lon.values[:, None],
                     slat[None, :], slon[None, :])
    proj = sta.roster.values == "project"
    bore = sta.is_borehole_any.values
    bband = np.isin(sta.klass.values, ["broadband", "broadband_borehole"])
    key = {(n, s): i for i, (n, s) in enumerate(zip(sta.net, sta.sta))}
    rows = []
    for y in YEARS:
        yy = sy[(sy.year == y) & sy.has_3c]
        idx = np.array([key[(n, s)] for n, s in zip(yy.net, yy.sta) if (n, s) in key], int)
        m3 = np.zeros(len(sta), bool)
        m3[idx] = True
        for tag, mask in (("all", m3), ("fleet", m3 & proj)):
            if not mask.any():
                sub = np.zeros((len(grid), 0))
                subaz = sub
                subb = np.zeros(0, bool)
            else:
                sub = D[:, mask]
                subaz = AZ[:, mask]
                subb = bore[mask]
            subbb = bband[mask] if mask.any() else np.zeros(0, bool)
            n60 = (sub <= 60).sum(1) if sub.shape[1] else np.zeros(len(grid), int)
            n100 = (sub <= 100).sum(1) if sub.shape[1] else np.zeros(len(grid), int)
            n150 = (sub <= 150).sum(1) if sub.shape[1] else np.zeros(len(grid), int)
            if sub.shape[1]:
                n60_bb = (sub[:, subbb] <= 60).sum(1) if subbb.any() else np.zeros(len(grid), int)
                n100_bb = (sub[:, subbb] <= 100).sum(1) if subbb.any() else np.zeros(len(grid), int)
            else:
                n60_bb = n100_bb = np.zeros(len(grid), int)
            near = sub.min(1) if sub.shape[1] else np.full(len(grid), np.nan)
            gaps = np.full(len(grid), 360.0)
            gaps60 = np.full(len(grid), 360.0)
            bfrac = np.full(len(grid), np.nan)
            bfrac60 = np.full(len(grid), np.nan)
            if sub.shape[1]:
                for i in range(len(grid)):
                    sel = sub[i] <= 100
                    if sel.any():
                        gaps[i] = max_az_gap(subaz[i][sel])
                        bfrac[i] = float(subb[sel].mean())
                    s60 = sub[i] <= 60
                    if s60.any():
                        gaps60[i] = max_az_gap(subaz[i][s60])
                        bfrac60[i] = float(subb[s60].mean())
            rows.append(pd.DataFrame(dict(
                cell=grid.cell, cell_lat=grid.cell_lat, cell_lon=grid.cell_lon,
                slab_km=grid.slab_km, region=grid.region, year=y, roster=tag,
                n60=n60, n100=n100, n150=n150, n60_bb=n60_bb, n100_bb=n100_bb,
                nearest_km=np.round(near, 1), az_gap_deg=np.round(gaps, 0),
                az_gap_60_deg=np.round(gaps60, 0),
                borehole_frac_100=np.round(bfrac, 3),
                borehole_frac_60=np.round(bfrac60, 3))))
    return pd.concat(rows, ignore_index=True), D


# ============================================================== STAGE D — regions
def region_tables(grid, sta, sy, D, dens):
    """Region membership = every station within 150 km of ANY cell of the region."""
    reg_sta, reg_year = [], []
    sy3 = sy[sy.has_3c]
    yrs_by_sta = sy3.groupby(["net", "sta"]).year.apply(lambda s: sorted(s)).to_dict()
    for name, a, b in REGIONS:
        gi = np.where(grid.region.values == name)[0]
        if not len(gi):
            continue
        dmin = D[gi].min(0)
        mem = np.where(dmin <= 150.0)[0]
        for i in mem:
            r = sta.iloc[i]
            ys = yrs_by_sta.get((r.net, r.sta), [])
            reg_sta.append(dict(region=name, net=r.net, sta=r.sta, lat=r.lat, lon=r.lon,
                                roster=r.roster, klass=r.klass, role=r.role,
                                min_cell_dist_km=round(float(dmin[i]), 1),
                                n_3c_years=len(ys),
                                y0_3c=(ys[0] if ys else None), y1_3c=(ys[-1] if ys else None),
                                years_3c=";".join(str(v) for v in ys)))
        for tag in ("all", "fleet"):
            for y in YEARS:
                d = dens[(dens.region == name) & (dens.year == y) & (dens.roster == tag)]
                if not len(d):
                    continue
                nm = sy3[sy3.year == y][["net", "sta"]].merge(
                    sta.iloc[mem][["net", "sta", "roster"]], on=["net", "sta"])
                if tag == "fleet":
                    nm = nm[nm.roster == "project"]
                reg_year.append(dict(
                    region=name, lat_min=a, lat_max=b, year=y, roster=tag,
                    n_cells=len(d), n_member_stations=len(nm),
                    med_n100=float(np.median(d.n100)), med_n60=float(np.median(d.n60)),
                    med_n150=float(np.median(d.n150)),
                    med_n60_bb=float(np.median(d.n60_bb)),
                    med_n100_bb=float(np.median(d.n100_bb)),
                    med_nearest_km=float(np.nanmedian(d.nearest_km)),
                    med_az_gap=float(np.nanmedian(d.az_gap_deg)),
                    med_az_gap_60=float(np.nanmedian(d.az_gap_60_deg)),
                    med_borehole_frac=float(np.nanmedian(d.borehole_frac_100)),
                    med_borehole_frac_60=float(np.nanmedian(d.borehole_frac_60)),
                    frac_ge3=float((d.n100 >= 3).mean()), frac_ge4=float((d.n100 >= 4).mean()),
                    frac_ge6=float((d.n100 >= 6).mean()),
                    frac_ge4_gap180=float(((d.n100 >= 4) & (d.az_gap_deg <= 180)).mean()),
                    frac_ge3_60=float((d.n60 >= 3).mean()),
                    frac_ge4_60=float((d.n60 >= 4).mean()),
                    frac_ge6_60=float((d.n60 >= 6).mean()),
                    frac_ge4_60_gap180=float(((d.n60 >= 4)
                                              & (d.az_gap_60_deg <= 180)).mean()),
                    frac_ge3_60_bb=float((d.n60_bb >= 3).mean()),
                    frac_ge4_60_bb=float((d.n60_bb >= 4).mean()),
                    frac_ge6_60_bb=float((d.n60_bb >= 6).mean())))
    return pd.DataFrame(reg_sta), pd.DataFrame(reg_year)


# ================================================================= STAGE E — cost
BYTES_PER_SAMPLE = 1.1        # measured on this project's own miniSEED day files (see NOTES)


def measure_bytes_per_sample(nmax=600, seed=0):
    """Empirical compressed-miniSEED bytes/sample from the day files already on disk.

    Two corrections over a naive size/(3*fs*86400):
      * the sample count comes from the file's own record headers (obspy headonly), not from
        the metadata sample rate — the day files hold whichever band was downloaded (mostly
        40-Hz BH), so scaling by the metadata max fs (100 Hz HH) underestimates by ~2x;
      * bytes/sample DEPENDS ON THE SAMPLE RATE (STEIM2 encodes differences: 1.43 B/sample at
        40 Hz vs 0.91 at 100 Hz here).  We therefore fit b(fs) = c0 + c1/fs across the rates
        present on disk and evaluate it per station-year.  77 % of the fleet's 3-C
        station-years are 100-Hz HH, so using a pooled median would overstate the volume ~30 %.
    Returns (model, per_rate_medians, pooled_median, n_files).
    """
    from obspy import read

    files = sorted(glob.glob(f"{ROOT}/lfe_catalog_repair/traces/*/*/*.mseed"))
    files = [f for f in files if os.path.getsize(f) >= 1e5]
    if len(files) > nmax:
        idx = np.random.default_rng(seed).choice(len(files), nmax, replace=False)
        files = [files[i] for i in sorted(idx)]
    per = defaultdict(list)
    for f in files:
        try:
            st = read(f, headonly=True)
        except Exception:                                       # noqa: BLE001
            continue
        npts = int(sum(tr.stats.npts for tr in st))
        if npts <= 0:
            continue
        fs = float(np.median([tr.stats.sampling_rate for tr in st]))
        per[round(fs)].append(os.path.getsize(f) / npts)
    if not per:
        return (lambda fs: BYTES_PER_SAMPLE), {}, BYTES_PER_SAMPLE, 0
    pooled = float(np.median(np.concatenate([np.asarray(v) for v in per.values()])))
    rates = {k: float(np.median(v)) for k, v in per.items() if len(v) >= 10}
    n = int(sum(len(v) for v in per.values()))
    if len(rates) >= 2:
        x = np.array([1.0 / k for k in rates])
        y = np.array([rates[k] for k in rates])
        c1, c0 = np.polyfit(x, y, 1)

        def model(fs):
            return float(max(0.2, c0 + c1 / float(fs)))
    else:
        def model(fs):
            return pooled
    return model, {k: round(v, 3) for k, v in sorted(rates.items())}, pooled, n


def cost_table(sy, bmodel, gpu_min_per_sta_year, gpu_src, bps_desc=""):
    rows = []
    bps20 = bmodel(MIN_FS)
    for y in YEARS:
        d = sy[(sy.year == y) & sy.has_3c]
        if not len(d):
            continue
        ns = len(d)
        frac = d.days.clip(0, 366) / 365.0
        bps = d.fs.map(bmodel)
        bps_lo = d.fs_lo.map(bmodel)
        raw_gb = float((3 * d.fs * 86400 * 365 * frac * bps).sum() / 1e9)
        lo_gb = float((3 * d.fs_lo * 86400 * 365 * frac * bps_lo).sum() / 1e9)
        dec_ms_gb = float((3 * MIN_FS * 86400 * 365 * frac * bps20).sum() / 1e9)
        dec_f32_gb = float((3 * MIN_FS * 86400 * 365 * frac * 4).sum() / 1e9)
        sy_eff = float(frac.sum())
        rows.append(dict(year=y, n_stations_3c=ns, station_years_effective=round(sy_eff, 1),
                         raw_download_GB_native=round(raw_gb, 1),
                         raw_download_GB_lowest_band=round(lo_gb, 1),
                         decim20_mseed_GB=round(dec_ms_gb, 1),
                         decim20_float32_GB=round(dec_f32_gb, 1),
                         gpu_pick_hours=round(sy_eff * gpu_min_per_sta_year / 60, 2)))
    t = pd.DataFrame(rows)
    tot = dict(year="TOTAL", n_stations_3c=int(sy[sy.has_3c].groupby(["net", "sta"]).ngroups),
               station_years_effective=round(t.station_years_effective.sum(), 1),
               raw_download_GB_native=round(t.raw_download_GB_native.sum(), 1),
               raw_download_GB_lowest_band=round(t.raw_download_GB_lowest_band.sum(), 1),
               decim20_mseed_GB=round(t.decim20_mseed_GB.sum(), 1),
               decim20_float32_GB=round(t.decim20_float32_GB.sum(), 1),
               gpu_pick_hours=round(t.gpu_pick_hours.sum(), 2))
    t = pd.concat([t, pd.DataFrame([tot])], ignore_index=True)
    t["note"] = ""
    t.loc[t.year == "TOTAL", "note"] = (
        f"bytes/sample model b(fs)=c0+c1/fs fitted to REAL day files (record-header sample "
        f"counts, compressed miniSEED STEIM2): {bps_desc}; b(20 Hz)={bmodel(MIN_FS):.2f} is an "
        f"EXTRAPOLATION. native = the band the station-year scores best on (often 100-Hz HH), "
        f"lowest_band = the cheapest band that still gives 3-C >=20 Hz (what --target-fs 40 "
        f"downloads). gpu_pick_hours are QUEUE-hours (divide by the number of concurrent "
        f"queues for wall-clock); gpu={gpu_min_per_sta_year:.2f} min/station-year — {gpu_src}")
    return t


def associate_locate_rows(sy, scope):
    """Order-of-magnitude PyOcto + NonLinLoc CPU cost, with the arithmetic in the note."""
    fp = f"{OUT}/e2_throughput.csv"
    picks_per_day = 1800.0
    src = "ASSUMED 1800 picks/station-day"
    if os.path.exists(fp):
        t = pd.read_csv(fp)
        picks_per_day = float(np.median(t.n_picks / t.n_days))
        src = (f"measured median {picks_per_day:.0f} S picks/station-day over {len(t)} E2 pilot "
               f"runs (LFEDetect threshold as run in E2; association discards most of them)")
    styr = float(sy[sy.has_3c].days.clip(0, 366).sum() / 365.0)
    n_picks = styr * 365 * picks_per_day
    PYOCTO_PICKS_PER_CPU_MIN = 2.0e5      # assumption: PyOcto 4-D octree, single core
    NLL_S_PER_EVENT = 0.03                # assumption: EDT search on a precomputed TT grid
    picks_per_event = 1500.0              # assumption: ~1 associated event per 1500 raw picks
    n_ev = n_picks / picks_per_event
    return pd.DataFrame([
        dict(year="ASSOC", scope=scope, n_stations_3c=np.nan,
             station_years_effective=round(styr, 1), raw_download_GB_native=np.nan,
             raw_download_GB_lowest_band=np.nan,
             decim20_mseed_GB=np.nan, decim20_float32_GB=np.nan,
             gpu_pick_hours=np.nan,
             note=f"picker output {n_picks/1e9:.2f}e9 S picks ({src}); PyOcto at "
                  f"{PYOCTO_PICKS_PER_CPU_MIN:.0e} picks/CPU-min -> "
                  f"{n_picks/PYOCTO_PICKS_PER_CPU_MIN/60:.0f} CPU-hours "
                  f"(embarrassingly parallel by region x day)"),
        dict(year="LOCATE", scope=scope, n_stations_3c=np.nan,
             station_years_effective=round(styr, 1), raw_download_GB_native=np.nan,
             raw_download_GB_lowest_band=np.nan,
             decim20_mseed_GB=np.nan, decim20_float32_GB=np.nan, gpu_pick_hours=np.nan,
             note=f"assume 1 associated event per {picks_per_event:.0f} raw picks -> "
                  f"{n_ev/1e6:.1f}e6 events; NonLinLoc at {NLL_S_PER_EVENT*1000:.0f} ms/event -> "
                  f"{n_ev*NLL_S_PER_EVENT/3600:.0f} CPU-hours"),
    ])


# ========================================================== STAGE F — ground truth
def groundtruth_table():
    ext = f"{ROOT}/lfe_catalog_repair/results/external"
    rows = []

    def bin_by_region(lat, label, source, note, n_total):
        for name, a, b in REGIONS:
            n = int(((lat >= a) & (lat < b)).sum())
            rows.append(dict(catalog=label, source=source, region=name,
                             n_events_or_families=n, note=note, n_total=n_total))

    fp = f"{ext}/lin_catalog/EQloc_001_0.1_3_S.csv"
    if os.path.exists(fp):
        lat = pd.read_csv(fp, usecols=["lat"]).lat.values
        bin_by_region(lat, "Lin2024_events", "Zenodo 10.5281/zenodo.10016020",
                      "1.06M located LFEs 2005-2022, no per-station picks", len(lat))
    fp = f"{ext}/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
    if os.path.exists(fp):
        b = pd.read_csv(fp, usecols=["lat", "lon"])
        bin_by_region(b.lat.values, "Bostock2015_events", "Slow Earthquake Database",
                      "269k events / 130 families 2003-2013 (origin times)", len(b))
        fam = b.drop_duplicates(["lat", "lon"])
        bin_by_region(fam.lat.values, "Bostock2015_families", "Slow Earthquake Database",
                      "family hypocentres", len(fam))
    fp = f"{ext}/other/seisbench_lfestacks_cascadia/metadata.csv"
    if os.path.exists(fp):
        m = pd.read_csv(fp, usecols=["source_latitude_deg", "source_family_id"], low_memory=False)
        m = m.dropna().drop_duplicates("source_family_id")
        bin_by_region(m.source_latitude_deg.values, "SeisBench_LFEstacks_families",
                      "SeisBench LFEStacksCascadiaBostock2015",
                      "PER-STATION P/S PICKS -> fine-tuning substrate", len(m))
    fp = f"{ROOT}/catalogs/ducellier_catalog/data/Plourde_2015/templates_list.txt"
    if os.path.exists(fp):
        t = pd.read_csv(fp, sep=r"\s+", skiprows=1,
                        names=["temp", "fam", "lat", "lon", "dep", "eh", "ez", "n"])
        f2 = t.drop_duplicates("fam")
        bin_by_region(f2.lat.values, "Ducellier2022_SouthernCascadia_families",
                      "catalogs/ducellier_catalog (Zenodo 6374239) + Plourde2015 locations",
                      "~478k LFEs, 66 templates / 37 located families, 2004-2011", len(f2))
    fp = f"{ext}/other/LFEcat_MTJ_2018-2024_shelly2025.csv"
    if os.path.exists(fp):
        s = pd.read_csv(fp)
        nf = s.ID.nunique() if "ID" in s else np.nan
        for name, a, b in REGIONS:
            n = int(nf) if name == "R7_NorCal_MTJ" else 0
            rows.append(dict(catalog="Shelly2025_MTJ_families", source="USGS 10.5066/P1TCKK7G (CC0)",
                             region=name, n_events_or_families=n,
                             note="61,441 detections / 27 families 2018-2024; NO lat/lon in the "
                                  "release -> region assigned from the paper (Mendocino/Gorda)",
                             n_total=len(s)))
    fp = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
    if os.path.exists(fp):
        p = pd.read_csv(fp, usecols=lambda c: c.lower() in ("lat", "latitude"))
        col = p.columns[0]
        bin_by_region(p[col].values, "PNSN_tremor", "catalogs/pnsn_tremor_cascadia_full.csv",
                      "tremor epicentres - weak validator everywhere, not LFE ground truth",
                      len(p))
    return pd.DataFrame(rows)


# ==================================================================== STAGE G — figure
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
INK, INK2, MUT = "#0b0b0b", "#52514e", "#a8a7a1"


def make_figure(sta, sy, dens, regy, grid):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.gridspec import GridSpec

    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUT, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})
    fig = plt.figure(figsize=(19.5, 11.2))
    gs = GridSpec(2, 7, figure=fig, height_ratios=[1.55, 1.0], hspace=0.30, wspace=0.42,
                  left=0.045, right=0.985, top=0.945, bottom=0.085)

    n3 = sy[sy.has_3c].groupby(["net", "sta"]).year.nunique()
    sta = sta.copy()
    sta["n3y"] = [int(n3.get((n, s), 0)) for n, s in zip(sta.net, sta.sta)]

    # ---- (a) station map ------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0:2])
    ax.plot(grid.cell_lon, grid.cell_lat, "s", ms=2.0, color="#e8e7e3", zorder=1)
    ex = sta[sta.roster == "extra"]
    ax.plot(ex.lon, ex.lat, ".", ms=2.6, color=MUT, alpha=0.8, zorder=3,
            label=f"non-project FDSN 3-C, <=150 km ({len(ex)})")
    pj = sta[sta.roster == "project"]
    cmap = ListedColormap(SEQ)
    norm = BoundaryNorm([0, 1, 4, 8, 11, 14, 16, 18], cmap.N)
    MARK = {"broadband": ("o", 30), "broadband_borehole": ("P", 40),
            "borehole_short_period": ("^", 36), "short_period": ("s", 22), "other": ("x", 22)}
    for k, (m, s) in MARK.items():
        d = pj[(pj.klass == k) & (pj.n3y > 0)]
        if len(d):
            ax.scatter(d.lon, d.lat, marker=m, s=s, c=d.n3y, cmap=cmap, norm=norm,
                       edgecolors="white", linewidths=0.5, zorder=4,
                       label=f"{k} ({len(d)})")
    z = pj[pj.n3y == 0]
    if len(z):
        ax.scatter(z.lon, z.lat, marker="x", s=26, c="#e34948", linewidths=1.0, zorder=5,
                   label=f"project, never 3-C >=20 Hz ({len(z)})")
    for name, a, b in REGIONS:
        ax.axhline(a, color=MUT, lw=0.7, ls="--", zorder=2)
        ax.text(-118.95, (a + b) / 2, name.split("_", 1)[1].replace("_", " "), fontsize=7,
                color=INK2, ha="right", va="center", zorder=6)
    ax.axhline(REGIONS[0][2], color=MUT, lw=0.7, ls="--", zorder=2)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.015,
                      fraction=0.045, ticks=[0, 1, 4, 8, 11, 14, 16, 18])
    cb.set_label(f"years with 3-C >=20 Hz ({YEARS[0]}-{YEARS[-1]})", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    ax.set_xlim(-128.3, -118.9); ax.set_ylim(39.6, 50.4)
    ax.set_title("a  Station inventory: class and 3-C availability", loc="left", fontsize=9.5)
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ax.legend(fontsize=6.6, loc="lower left", frameon=True, framealpha=0.92, borderpad=0.4)
    ax.grid(alpha=0.18, lw=0.5)

    # ---- (b,c,d) cell density snapshots --------------------------------------
    dcm = ListedColormap(SEQ)
    dnorm = BoundaryNorm([0, 1, 2, 3, 4, 6, 9, 40], dcm.N)
    for j, yr in enumerate((2012, 2018, 2024)):
        ax = fig.add_subplot(gs[0, 2 + j])
        d = dens[(dens.year == yr) & (dens.roster == "fleet")]
        ax.scatter(d.cell_lon, d.cell_lat, marker="s", s=13, c=d.n60, cmap=dcm, norm=dnorm,
                   linewidths=0)
        ss = sy[(sy.year == yr) & sy.has_3c & (sy.roster == "project")].merge(
            sta[["net", "sta", "lat", "lon"]], on=["net", "sta"])
        ax.plot(ss.lon, ss.lat, "^", ms=2.4, color=INK2, alpha=0.75, zorder=3)
        for name, a, b in REGIONS:
            ax.axhline(a, color=MUT, lw=0.6, ls="--")
        ax.axhline(REGIONS[0][2], color=MUT, lw=0.6, ls="--")
        ax.set_xlim(-126.4, -119.9); ax.set_ylim(39.6, 50.4)
        ax.set_title(f"{'bcd'[j]}  project-fleet 3-C stations <=60 km, {yr}",
                     loc="left", fontsize=8.6)
        ax.set_xlabel("longitude")
        if j == 0:
            ax.set_ylabel("latitude")
        else:
            ax.set_yticklabels([])
        ax.grid(alpha=0.18, lw=0.5)
        if j == 2:
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=dnorm, cmap=dcm), ax=ax, pad=0.02,
                              fraction=0.05, ticks=[0, 1, 2, 3, 4, 6, 9, 40])
            cb.set_label("N stations within 60 km", fontsize=7.5)
            cb.ax.tick_params(labelsize=7)

    # ---- (e) locatable fraction, small multiples -----------------------------
    ra = regy[regy.roster == "fleet"]
    rall = regy[regy.roster == "all"]
    for i, (name, a, b) in enumerate(REGIONS):
        ax = fig.add_subplot(gs[1, i])
        d = ra[ra.region == name].sort_values("year")
        e = rall[rall.region == name].sort_values("year")
        ax.bar(d.year, d.frac_ge4_60, width=0.74, color=CAT[0], edgecolor="none", zorder=2,
               label="fleet: >=4 3-C stations <=60 km")
        ax.plot(d.year, d.frac_ge6_60, "o", ms=4.2, color=CAT[1], mec="white", mew=0.6, zorder=4,
                label="fleet: >=6 <=60 km")
        ax.plot(d.year, d.frac_ge4_60_bb, "v", ms=4.0, color=CAT[3], mec="white", mew=0.5,
                zorder=5, label="fleet: >=4 BROADBAND <=60 km")
        ax.plot(e.year, e.frac_ge4_60, "s", ms=3.4, mfc="none", color=MUT, mew=0.9, zorder=3,
                label="all FDSN in region: >=4 <=60 km")
        ax.set_ylim(0, 1.05); ax.set_xlim(2008.3, 2026.7)
        ax.set_title(f"{name.split('_',1)[0]} {name.split('_',1)[1].replace('_',' ')}\n"
                     f"{a:.1f}-{b:.1f} N, {int(d.n_cells.iloc[0])} cells",
                     loc="left", fontsize=8)
        ax.set_xticks([2010, 2014, 2018, 2022, 2026])
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.set_ylabel("fraction of cells locatable")
        else:
            ax.set_yticklabels([])
        ax.grid(axis="y", alpha=0.2, lw=0.5)
        ax.set_axisbelow(True)
        if i == 6:
            ax.legend(fontsize=6.4, loc="upper left", frameon=True, framealpha=0.94)
    fig.text(0.045, 0.045,
             "e  Expected locatable fraction per association region and year — cells on the Slab2 "
             "interface, 30-45 km depth, 0.25 deg, station radius 60 km. Bars/markers = project "
             "fleet; open squares = every 3-C FDSN station in the region.", fontsize=8.6,
             color=INK)
    fig.text(0.045, 0.020,
             "   At a 100 km radius every cell has >=4 3-C stations in every year for both "
             "rosters (fraction 1.00 throughout): geometry at that radius is never the limit — "
             "station PROXIMITY and site class are. Metadata only; no data-quality check.",
             fontsize=8.6, color=INK2)
    fig.suptitle("F1 — fleet inventory and association-region design for a fleet-wide LFE picker "
                 "run (no waveforms downloaded)", fontsize=11.5, x=0.045, ha="left", color=INK)
    fp = f"{FIG}/f1_fleet_inventory.png"
    fig.savefig(fp, dpi=135)
    print(f"[fig] {fp}")


# ======================================================================== main
def main():
    print("=" * 100)
    print("F1 FLEET INVENTORY")
    print("=" * 100)

    # ---- A roster
    roster = build_roster()
    roster, missing, n_dvv, n_inv = mark_usage(roster)
    print(f"[A] project lists: {len(roster)} (net,sta) rows | "
          f"dv/v products {n_dvv} | inversion tags {n_inv} | unmapped tags {missing}")

    # ---- B metadata
    nets = sorted(set(NETS_PROJECT) | set(roster.net) | set(NETS_EXTRA))
    print(f"[B] FDSN station-level metadata, one request per network: {nets}")
    ch = channel_records(nets)
    print(f"    {len(ch)} high-gain seismometer channel epochs, "
          f"{ch.groupby(['net','sta']).ngroups} distinct stations in bbox")

    sy = station_year_table(ch)

    # station master: coords from FDSN, fall back to project table
    cs = ch.groupby(["net", "sta"]).agg(lat=("lat", "median"), lon=("lon", "median"),
                                        elev=("elev", "median"), depth_m=("depth_m", "max"),
                                        fs_max=("fs", "max")).reset_index()
    cs["klass"] = [classify(ch[(ch.net == n) & (ch.sta == s)], None)
                   for n, s in zip(cs.net, cs.sta)]
    cs["bands"] = [",".join(sorted(set(ch[(ch.net == n) & (ch.sta == s)]
                                       .query("fs >= @MIN_FS").band)))
                   for n, s in zip(cs.net, cs.sta)]
    sta = roster.merge(cs, on=["net", "sta"], how="outer", indicator=True)
    sta["roster"] = np.where(sta._merge == "right_only", "extra", "project")
    sta["fdsn_metadata"] = sta._merge != "left_only"
    sta["lat"] = sta.lat.fillna(sta.lat_tab)
    sta["lon"] = sta.lon.fillna(sta.lon_tab)
    sta["role"] = sta.role.fillna("fdsn_extra")
    sta["klass"] = sta.klass.fillna("unknown")
    sta = sta.drop(columns=["_merge"])
    sta = sta[sta.lat.notna() & sta.lon.notna()].reset_index(drop=True)

    # only keep 'extra' stations that are actually useful: 3-C at some point AND in-region
    g3 = sy[sy.has_3c].groupby(["net", "sta"]).size()
    sta["n_3c_years"] = [int(sy[(sy.net == n) & (sy.sta == s) & sy.has_3c].year.nunique())
                         for n, s in zip(sta.net, sta.sta)]
    sta["z_only_ever"] = [bool(sy[(sy.net == n) & (sy.sta == s)].z_only.any())
                          for n, s in zip(sta.net, sta.sta)]
    sta["is_borehole_any"] = sta.depth_m.fillna(0) >= BOREHOLE_DEPTH_M
    grid = source_grid()
    dcell = haversine_km(grid.cell_lat.values[:, None], grid.cell_lon.values[:, None],
                         sta.lat.values[None, :], sta.lon.values[None, :]).min(0)
    sta["min_cell_dist_km"] = np.round(dcell, 1)
    drop = (sta.roster == "extra") & ((sta.n_3c_years == 0) | (sta.min_cell_dist_km > 150))
    print(f"    dropping {int(drop.sum())} non-project stations with no 3-C data or >150 km "
          f"from every source cell")
    sta = sta[~drop].reset_index(drop=True)
    sy = sy.merge(sta[["net", "sta", "roster", "klass", "role"]], on=["net", "sta"], how="inner")

    print(f"[B] station master: {len(sta)} stations "
          f"({int((sta.roster=='project').sum())} project + "
          f"{int((sta.roster=='extra').sum())} extra FDSN)")
    print("    class breakdown (project only):")
    print(sta[sta.roster == "project"].klass.value_counts().to_string())
    nomd = sta[(sta.roster == "project") & (~sta.fdsn_metadata)]
    print(f"    project stations with NO FDSN metadata returned: {len(nomd)} "
          f"{sorted(nomd.net + '.' + nomd.sta)[:20]}")
    z_only = sta[(sta.roster == "project") & (sta.n_3c_years == 0)]
    print(f"    project stations that NEVER have 3-C >=20 Hz: {len(z_only)}")
    print("      " + ", ".join(sorted((z_only.net + "." + z_only.sta).tolist())))
    orient = sy[sy.has_3c].groupby("orient").size()
    print(f"    3-C orientation codes across station-years: {orient.to_dict()}")

    # ---- C grid + density
    print(f"[C] source grid: {len(grid)} cells, {CELL_DEG} deg, Slab2 {DEPTH_BAND[0]}-"
          f"{DEPTH_BAND[1]} km, {LAT_BAND[0]}-{LAT_BAND[1]} N")
    print(grid.groupby("region").size().to_string())
    dens, D = cell_density(grid, sta, sy)

    # ---- D regions
    reg_sta, regy = region_tables(grid, sta, sy, D, dens)

    # ---- E cost
    bmodel, bps_rates, bps_pooled, nmeas = measure_bytes_per_sample()
    bps_desc = ("; ".join(f"{k} Hz -> {v:.3f} B/sample" for k, v in bps_rates.items())
                + f"; pooled median {bps_pooled:.3f}")
    fp = f"{OUT}/e2_throughput.csv"
    if os.path.exists(fp):
        t = pd.read_csv(fp)
        # aggregate over ALL pilot runs: total station-days / total wall-seconds
        spm = float(t.n_days.sum() / (t.seconds.sum() / 60.0))
        lo = float(t.station_days_per_min.min())
        hi = float(t.station_days_per_min.max())
        gpu_min = 365.0 / spm
        gsrc = (f"E2 pilot measured: {spm:.1f} station-days/min aggregate over {len(t)} runs "
                f"(per-run {lo:.0f}-{hi:.0f}) on {t.gpu.iloc[0]}, "
                f"weights {sorted(set(t.model))}, {int(t.n_days.sum())} station-days total")
    else:
        gpu_min, gsrc = 3.0, "ASSUMED 1 station-year per 3 GPU-minutes (no e2_throughput.csv)"
    print(f"[E] bytes/sample measured on {nmeas} day-files (actual record-header sample "
          f"counts): {bps_desc}; model at 20 Hz {bmodel(MIN_FS):.3f} (extrapolated), "
          f"at 40 Hz {bmodel(40):.3f}, at 100 Hz {bmodel(100):.3f}")
    print(f"[E] GPU picker: {gpu_min:.2f} min/station-year — {gsrc}")
    cost = cost_table(sy[sy.roster == "project"], bmodel, gpu_min, gsrc, bps_desc)
    cost_all = cost_table(sy, bmodel, gpu_min, gsrc, bps_desc)
    cost["scope"] = "project_fleet"
    cost_all["scope"] = "project_plus_extra"
    cost = pd.concat([cost, cost_all,
                      associate_locate_rows(sy[sy.roster == "project"], "project_fleet"),
                      associate_locate_rows(sy, "project_plus_extra")], ignore_index=True)

    # ---- F ground truth
    gt = groundtruth_table()

    # ---- write
    for name, df in [("f1_stations", sta), ("f1_station_years", sy),
                     ("f1_cell_density", dens), ("f1_regions", regy),
                     ("f1_region_stations", reg_sta), ("f1_cost_estimate", cost),
                     ("f1_groundtruth", gt)]:
        p = f"{OUT}/{name}.csv"
        df.to_csv(p, index=False)
        print(f"[out] {p}  ({len(df)} rows)")

    make_figure(sta, sy, dens, regy, grid)

    # ---- print key numbers
    print("\n" + "=" * 100)
    print("KEY NUMBERS")
    print("=" * 100)
    print("\n-- station-years with 3-C >=20 Hz, by year --")
    pv = sy[sy.has_3c].pivot_table(index="year", columns="roster", values="sta",
                                   aggfunc="count").fillna(0).astype(int)
    print(pv.to_string())
    for tag, lbl in (("fleet", "PROJECT FLEET only"), ("all", "every 3-C FDSN station in region")):
        r = regy[regy.roster == tag]
        print(f"\n-- median 3-C stations within 100 km per cell [{lbl}] --")
        print(r.pivot_table(index="year", columns="region", values="med_n100")
              .round(1).to_string())
        for col, desc in (("frac_ge3", ">=3 stations <=100 km"),
                          ("frac_ge4", ">=4 stations <=100 km"),
                          ("frac_ge3_60", ">=3 stations <=60 km"),
                          ("frac_ge4_60", ">=4 stations <=60 km"),
                          ("frac_ge6_60", ">=6 stations <=60 km"),
                          ("frac_ge4_60_bb", ">=4 BROADBAND stations <=60 km")):
            print(f"\n-- locatable fraction, cells with {desc} [{lbl}] --")
            print(r.pivot_table(index="year", columns="region", values=col)
                  .round(2).to_string())
        print(f"\n-- region summary, median over 2012-2024 [{lbl}] --")
        s = regy[(regy.roster == tag) & regy.year.between(2012, 2024)].groupby("region").agg(
            n_cells=("n_cells", "median"), n_stations=("n_member_stations", "median"),
            med_n60=("med_n60", "median"), med_n100=("med_n100", "median"),
            med_n60_bb=("med_n60_bb", "median"),
            med_nearest_km=("med_nearest_km", "median"),
            med_az_gap=("med_az_gap", "median"), med_az_gap_60=("med_az_gap_60", "median"),
            bore_frac=("med_borehole_frac", "median"),
            bore_frac_60=("med_borehole_frac_60", "median"),
            frac_ge4=("frac_ge4", "median"),
            frac_ge3_60=("frac_ge3_60", "median"), frac_ge4_60=("frac_ge4_60", "median"),
            frac_ge6_60=("frac_ge6_60", "median"),
            frac_ge4_60_gap180=("frac_ge4_60_gap180", "median"),
            frac_ge4_60_bb=("frac_ge4_60_bb", "median"))
        print(s.round(2).to_string())
    print("\n-- cost --")
    with pd.option_context("display.max_colwidth", 250, "display.width", 260):
        print(cost[cost.year.isin(["TOTAL", "ASSOC", "LOCATE"])].to_string(index=False))
    print("\n-- ground truth by region (events or families) --")
    print(gt.pivot_table(index="region", columns="catalog", values="n_events_or_families",
                         aggfunc="sum").fillna(0).astype(int).to_string())
    print("\nDONE")


if __name__ == "__main__":
    main()
