#!/usr/bin/env python
"""P3D (dense-north 2019-2026) shared infrastructure.  Prefix `s_`.

Everything the P3D steps need that is not already frozen elsewhere:

  * the P3D roster = the frozen S15 (in STAS_ORDER) + the ten new 3-C stations + the two Z-only
    helpers, with the CO-LOCATED SITE MERGE corrected (CLRS -> the B926 site, OSQM -> the B001
    site) so that no site is double counted in any network sum (Merlin P3D flaw 3);
  * geometry against the 130 Bostock families (hypocentral, family depths from e5r_families_2012);
  * a cache index that is NOT filtered to S15 (the override of `p_common.cache_index` Merlin asks
    for) and a ZNE day reader for the new stations;
  * the hp event table (Track E realness_v1_1, p_real >= 0.8) restricted to SSE windows, with the
    pre-registered episode split;
  * copies (never edits) of the frozen stacking primitives from `e5r_templates.py`.

Read-only inputs: ../lfe_extension, ../lfe_catalog_repair, lfe_picker/results/{o*,p*,q*,g*}.
Nothing here writes outside lfe_picker/results/p3d and lfe_picker/figures.
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
P3D = f"{RES}/p3d"
FIG = f"{PICK}/figures"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
CACHES = [f"{PICK}/cache20hz", f"{REP}/cache20hz", f"{EXT}/cache20hz"]
MM = f"{RES}/mm"
sys.path.insert(0, HERE)

FS = 20.0
NDAY = 1728000
R = 111.195
CS = np.cos(np.radians(48.6))

# ---------------------------------------------------------------- roster -----------------------
# frozen order (o_common.STAS_ORDER) -- NEVER reordered; new stations are APPENDED after it.
STAS_FROZEN = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
               "B009", "B010", "B011", "B926"]
NEW_3C = ["SYMB", "CLRS", "PTRF", "GOBB", "RNWY", "BUCKS", "TURTL", "LOPEZ", "OSQM", "ODUC"]
HELPERS_Z = ["STW", "UWFH"]                    # Z-only: flagged, never summed into a site statistic
# ---- COORDINATOR SCOPE AMENDMENT 2026-09-19 (pre-analysis, before any PQ stack was opened) ------
# The all-network inventory (results/p3d_puget_geometry/station_ranking.csv) found four PQ-network
# broadbands on southern Vancouver Island inside 45 km of 36-66 Bostock families each -- more
# geometry than SYMB alone.  Two pre-registered NETWORK VARIANTS are carried from here on:
#   V1 = the frozen P3D roster (everything already hashed in freeze.json)
#   V2 = V1 + the four PQ stations (+ UW.OHOH / UW.ODEER if they deliver)
# The Step-3 DECISION is taken on V2; V1 is reported alongside.
PQ_3C = ["LHLYB", "JRDRB", "ALBHB", "DAOB"]
EXTRA_3C = ["OHOH", "ODEER"]
V2_3C = PQ_3C + EXTRA_3C
STAS_P3D = STAS_FROZEN + NEW_3C + V2_3C + HELPERS_Z

NETSTA_NEW = {"SYMB": "CN.SYMB", "CLRS": "CN.CLRS", "PTRF": "CN.PTRF", "GOBB": "CN.GOBB",
              "RNWY": "UW.RNWY", "BUCKS": "UW.BUCKS", "TURTL": "UW.TURTL", "LOPEZ": "UW.LOPEZ",
              "OSQM": "UW.OSQM", "ODUC": "UW.ODUC", "STW": "UW.STW", "UWFH": "UW.UWFH"}
NETSTA_FROZEN = {"PGC": "CN.PGC", "VGZ": "CN.VGZ", "NLLB": "CN.NLLB", "SNB": "CN.SNB",
                 "LRIV": "UW.LRIV", "B001": "PB.B001", "B003": "PB.B003", "B004": "PB.B004",
                 "B005": "PB.B005", "B006": "PB.B006", "B007": "PB.B007", "B009": "PB.B009",
                 "B010": "PB.B010", "B011": "PB.B011", "B926": "PB.B926"}
NETSTA_V2 = {"LHLYB": "PQ.LHLYB", "JRDRB": "PQ.JRDRB", "ALBHB": "PQ.ALBHB", "DAOB": "PQ.DAOB",
             "OHOH": "UW.OHOH", "ODEER": "UW.ODEER"}
NETSTA = {**NETSTA_FROZEN, **NETSTA_NEW, **NETSTA_V2}

# ---- SITES.  The frozen merge (o_common.SITES / p_common.SITE_CLUSTERS) plus the two co-located
# pairs Merlin's flaw 3 requires.  CLRS is 25 m from B926 and OSQM 17 m from B001: they are the SAME
# SITE and must never both enter a site sum, or the E5 site-coherent false-alarm mechanism recurs.
SITE_CLUSTERS = {
    "PGCq":  ["PGC", "B009", "B010", "B011"],
    "LRIVq": ["LRIV", "B005", "B006", "B007"],
    "B926q": ["B926", "CLRS"],                 # co-located, 25 m  -- WITHHELD site
    "B001q": ["B001", "OSQM"],                 # co-located, 17 m
}
# withheld SITE set for P3D (Merlin Step 0): never in any network sum, used as Gate-1 confirmers.
WITHHELD_SITES = ["B926q", "GOBB", "LOPEZ", "TURTL"]

# ZNE rotation azimuths.  Boreholes from p_common.AZ12 (StationXML); every CN/UW surface station in
# the P3D roster is 0/90 exactly (results/p3d/station_epochs.csv chan_az_dip column, verified).
AZ12 = {
    "B001": (335.5, 65.5), "B003": (239.0, 329.0), "B004": (167.5, 257.5),
    "B005": (90.5, 180.5), "B006": (201.0, 291.0), "B007": (36.5, 126.5),
    "B009": (0.0, 90.0), "B010": (228.5, 318.5), "B011": (11.0, 101.0), "B926": (299.0, 29.0),
}
for _s in ["PGC", "VGZ", "NLLB", "SNB", "LRIV"] + NEW_3C + V2_3C + HELPERS_Z:
    AZ12[_s] = (0.0, 90.0)

# ---- COMPONENT DIP (downloader correction, 2026-09-19) -----------------------------------------
# SEED dip is degrees from horizontal, POSITIVE DOWN.  Every CN/UW/PB station in this project
# declares its vertical at dip -90 (positive UP); ALL FOUR PQ STATIONS declare HHZ at dip +90
# (positive DOWN), so their raw Z is polarity-INVERTED relative to the rest of the network.  The
# frozen 20 Hz cache stores raw traces, so the sign has to be applied on read.  ZSIGN[sta] = -1
# turns a dip +90 vertical into the network's positive-up convention; every other station gets +1
# and is therefore bit-identical to every earlier run.
ZSIGN = {}


def _load_chan_meta():
    """Parse chan_az_dip in results/p3d/station_epochs.csv -> per-station (az1, az2) and Z sign."""
    p = f"{P3D}/station_epochs.csv"
    if not os.path.exists(p):
        return
    E = pd.read_csv(p)
    for r in E.itertuples():
        spec = str(getattr(r, "chan_az_dip", "") or "")
        az = {}
        for part in spec.split(";"):
            if ":" not in part:
                continue
            ch, ad = part.split(":", 1)
            try:
                a, d = (float(x) for x in ad.split("/"))
            except ValueError:
                continue
            az[ch[-1]] = (a, d)
        if "Z" in az:
            ZSIGN[r.sta] = -1.0 if az["Z"][1] > 0 else 1.0
        if "N" in az and "E" in az:
            AZ12[r.sta] = (az["N"][0], az["E"][0])


_load_chan_meta()
for _s in list(AZ12):
    ZSIGN.setdefault(_s, 1.0)
# S_P3D_NO_DIP=1 reproduces the PRE-CORRECTION behaviour (raw Z for every station) so the
# before/after of the dip fix can be measured with one code path.
NO_DIP = os.environ.get("S_P3D_NO_DIP", "") == "1"


def zsign(sta):
    """+1 for a dip -90 (positive-up) vertical, -1 for a dip +90 (positive-down) one."""
    return 1.0 if NO_DIP else float(ZSIGN.get(sta, 1.0))

# pre-registered episode split (Merlin Q4).  Stations alive before 2023 calibrate on S01-S13 and are
# tested on S14-S22; BUCKS/TURTL/ODUC build on S18 and are tested on S21/S22 only.
CAL_EPISODES = [f"S{i:02d}" for i in range(1, 14)]
TEST_EPISODES = ["S14", "S18", "S21", "S22"]
LATE_BUILD_EPISODE = "S18"
LATE_STATIONS = ["BUCKS", "TURTL", "ODUC", "LHLYB", "DAOB"]
LATE_TEST_EPISODES = ["S21", "S22"]
# Stations whose MEASURED data extent (station_epochs.csv first_data_month, added by the downloader
# 2026-09-19) starts after S18: they can only build on S21 and be tested on S22.
LATE2_STATIONS = ["JRDRB", "ODEER"]
LATE2_BUILD_EPISODE = "S21"
LATE2_TEST_EPISODES = ["S22"]
# PQ time split (coordinator amendment): JRDRB/ALBHB calibrate on S08-S13 and test on S14-S22;
# LHLYB/DAOB build on S18 and test on S21/S22 only.  OHOH (2022-04) calibrates on S12-S13.
# Calibration episodes narrowed to the MEASURED extents: ALBHB delivers from 2022-09 (so only S12's
# tail and S13 are available, not S08-S13 as the amendment assumed) and OHOH from 2022-05.
CAL_OVERRIDE = {"ALBHB": ["S12", "S13"], "OHOH": ["S12", "S13"]}


def cal_episodes_for(sta):
    if sta in LATE2_STATIONS:
        return [LATE2_BUILD_EPISODE]
    if sta in LATE_STATIONS:
        return [LATE_BUILD_EPISODE]
    return CAL_OVERRIDE.get(sta, CAL_EPISODES)


def test_episodes_for(sta):
    if sta in LATE2_STATIONS:
        return LATE2_TEST_EPISODES
    return LATE_TEST_EPISODES if sta in LATE_STATIONS else TEST_EPISODES
# network epochs for the FA-only tier rule (Merlin Step 5)
NET_EPOCHS = [("D1", "2019-01-01", "2021-03-31"), ("D2", "2021-04-01", "2023-08-31"),
              ("D3", "2023-09-01", "2026-12-31")]
QUIET_MONTHS = ["2023-03", "2024-12", "2025-11", "2026-06"]
GEOM_RULE_KM = 60.0          # a priori site-sum rule: median hypocentral distance to the family field


# ---------------------------------------------------------------- geometry ---------------------
def epi_km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1, float) - np.asarray(lo2, float)) * R * CS,
                    (np.asarray(la1, float) - np.asarray(la2, float)) * R)


def hyp_km(fam_lat, fam_lon, fam_dep, sta_lat, sta_lon):
    return np.hypot(epi_km(fam_lat, fam_lon, sta_lat, sta_lon), np.asarray(fam_dep, float))


def families():
    return pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")


def holdouts():
    ho = json.load(open(f"{RES}/holdouts.json"))
    hid = {int(h["bkey"]) for h in ho["hidden_families"]}
    iso = {int(h["bkey"]) for h in ho["hidden_families"] if h["stratum"] == "isolated_ge5km"}
    return hid, iso, ho


def station_table():
    """The P3D roster: frozen S15 + ten new 3-C + two Z-only helpers, with sites and roles.
    Written to results/p3d/stations_p3d.csv by s_step0_roster.py; e5r_stations.csv is NEVER edited."""
    p = f"{P3D}/stations_p3d.csv"
    if os.path.exists(p):
        return pd.read_csv(p).set_index("sta")
    return build_station_table()


def build_station_table():
    old = pd.read_csv(f"{IN}/e5r_stations.csv")
    new = pd.read_csv(f"{P3D}/station_epochs.csv")
    site_of = {}
    for q, mem in SITE_CLUSTERS.items():
        for s in mem:
            site_of[s] = q
    rows = []
    for i, s in enumerate(STAS_FROZEN):
        r = old[old.sta == s].iloc[0]
        rows.append(dict(sta=s, net=r.net, netsta=r.netsta, band=str(r.chan)[:2], fs_native=r.fs_native,
                         lat=r.lat, lon=r.lon, elev_m=r.elev_m, klass=r.klass, mode="3C",
                         role="frozen", site=site_of.get(s, s), order=i,
                         start="2003-01-01", end="2026-12-31",
                         az1=AZ12[s][0], az2=AZ12[s][1]))
    for j, s in enumerate(NEW_3C + HELPERS_Z):
        r = new[new.sta == s].iloc[0]
        rows.append(dict(sta=s, net=r.net, netsta=f"{r.net}.{r.sta}", band=r.band, fs_native=r.fs,
                         lat=r.lat, lon=r.lon, elev_m=r.elev, klass="surface_bb",
                         mode=("3C" if s in NEW_3C else "Z"),
                         role=("new" if s in NEW_3C else "helper_z"),
                         site=site_of.get(s, s), order=len(STAS_FROZEN) + j,
                         start=r.start, end=r.end, az1=AZ12[s][0], az2=AZ12[s][1]))
    T = pd.DataFrame(rows)
    # LRIV upgraded BH -> HH on 2019-05-30 (A-2 step 4); its 2012 templates FAILED transfer, so in
    # P3D it is a transfer TARGET like the new stations, not a frozen template station.
    T.loc[T.sta == "LRIV", "band_p3d"] = "HH"
    T["band_p3d"] = T.band_p3d.fillna(T.band)
    T["transfer_target"] = T.sta.isin(NEW_3C + ["LRIV"])
    T["withheld"] = T.site.isin(WITHHELD_SITES)
    T["summable"] = T["mode"].eq("3C")          # Z-only helpers are never summed
    return T.set_index("sta")


def site_table(ST=None):
    """One row per SITE: representative position (mean of members), members, median hypocentral
    distance to the 130-family field, and the a priori <=60 km geometry flag."""
    ST = station_table() if ST is None else ST
    F = families()
    rows = []
    for site, g in ST[ST.summable].groupby("site"):
        la, lo = float(g.lat.mean()), float(g.lon.mean())
        d = hyp_km(F.lat.values, F.lon.values, F.dep.values, la, lo)
        rows.append(dict(site=site, members="+".join(sorted(g.index)), n_members=len(g),
                         lat=la, lon=lo, order=int(g.order.min()),
                         role=("new" if g.role.astype(str).str.startswith("new").all() else
                               ("mixed" if g.role.astype(str).str.startswith("new").any()
                                else "frozen")),
                         med_hyp_km=float(np.median(d)), n_fam_45=int((d <= 45).sum()),
                         n_fam_55=int((d <= 55).sum()), n_fam_60=int((d <= 60).sum()),
                         withheld=bool(g.withheld.any())))
    S = pd.DataFrame(rows).sort_values("order").reset_index(drop=True)
    S["geom_le60"] = S.med_hyp_km <= GEOM_RULE_KM
    return S


def dense_sites(S=None, drop_withheld=True, at_date=None):
    """The a priori dense site set: median hypocentral distance to the family field <= 60 km."""
    S = site_table() if S is None else S
    k = S[S.geom_le60]
    if drop_withheld:
        k = k[~k.withheld]
    if at_date is not None:
        ST = station_table()
        alive = {s for s in ST.index if str(ST.start[s]) <= str(at_date) <= str(ST.end[s])}
        k = k[[any(m in alive for m in r.split("+")) for r in k.members]]
    return list(k.site)


def old7_sites():
    """The seven P3 discovery sites (S15 sites minus the withheld SNB and B926)."""
    return ["PGCq", "VGZ", "NLLB", "LRIVq", "B001q", "B003", "B004"]


# ---------------------------------------------------------------- SSE windows ------------------
def sse_windows():
    W = pd.read_csv(f"{P3D}/sse_windows.csv")
    for c in ("core_start", "core_end", "start", "end"):
        W[c] = pd.to_datetime(W[c]).astype("datetime64[ns]")
    return W


def sse_days(ids=None):
    """{episode id: DatetimeIndex of padded-window days}."""
    W = sse_windows()
    out = {}
    for r in W.itertuples():
        if ids is not None and r.id not in ids:
            continue
        out[r.id] = pd.date_range(r.start, r.end, freq="D")
    return out


def episode_of(ts):
    """episode id for an array of pandas Timestamps (padded windows), else ''."""
    W = sse_windows()
    t = pd.to_datetime(pd.Series(ts)).values.astype("datetime64[ns]")
    out = np.array([""] * len(t), dtype=object)
    for r in W.itertuples():
        m = (t >= np.datetime64(r.start)) & (t <= np.datetime64(r.end) + np.timedelta64(1, "D"))
        out[m] = r.id
    return out


# ---------------------------------------------------------------- hp events --------------------
def hp_events(min_p=0.8, episodes=None, column="hp"):
    """Track E realness v1.1 high-purity events (p_real >= 0.8) on the frozen origin-time axis.

    `t` in realness_v1_1 is the DETECTION time on the same axis Track C emits (family origin axis;
    the frozen family offsets are already in it, cf. p_common.bostock)."""
    cols = ["det_uid", "family_id", "t", "date", "p_real", "hp", "cat_core", "net_mad", "tier_epoch"]
    D = pd.read_parquet(f"{EXT}/results/trackE/realness_v1_1.parquet", columns=cols)
    D = D.rename(columns={"family_id": "bkey"})
    D["dt"] = pd.to_datetime(D.date).astype("datetime64[ns]")
    if column == "hp":
        D = D[D.p_real >= min_p]
    else:
        D = D[D.cat_core & (D.p_real >= min_p)]
    D["episode"] = episode_of(D.dt)
    if episodes is not None:
        D = D[D.episode.isin(episodes)]
    return D.reset_index(drop=True)


# ---------------------------------------------------------------- caches / days ----------------
_MASK = None


def masked_days():
    """{(sta, year, doy)} excluded by the Step-0 acceptance rule (coverage_frac < 0.9)."""
    global _MASK
    if _MASK is None:
        p = f"{P3D}/masked_station_days.csv"
        _MASK = set()
        if os.path.exists(p):
            m = pd.read_csv(p)
            if len(m):
                d = pd.to_datetime(m.date)
                _MASK = set(zip(m.sta, d.dt.year.astype(int), d.dt.dayofyear.astype(int)))
    return _MASK


def cache_index(years=None, stas=None, roots=None, use_mask=True):
    """{(sta, year, doy): path}.  OVERRIDE of p_common.cache_index: no S15 filter, so the ten new
    stations and the Z-only helpers are visible (p_common silently skips them)."""
    idx = {}
    known = set(STAS_P3D)
    # REFEREE GUARD: 2011 and 2013 are sealed for P3D.  Nothing in this workspace may open them.
    if years is not None and (set(years) & {2011, 2013}):
        raise RuntimeError("P3D may not read 2011/2013 (sealed referee years)")
    SEALED = {2011, 2013}
    for root in (CACHES if roots is None else roots):
        for p in glob.glob(f"{root}/*/*/*.npz"):
            parts = p.split(os.sep)
            netsta, yr, fn = parts[-3], int(parts[-2]), parts[-1]
            sta = netsta.split(".")[1]
            if sta not in known:
                continue
            if stas is not None and sta not in stas:
                continue
            if yr in SEALED:
                continue
            if use_mask and (sta, yr, int(fn[:3])) in masked_days():
                continue
            if years is not None and yr not in years:
                continue
            idx.setdefault((sta, yr, int(fn[:3])), p)
    return idx


def rot_matrix(sta):
    a1, a2 = np.radians(AZ12[sta])
    A = np.array([[np.cos(a1), np.sin(a1)], [np.cos(a2), np.sin(a2)]])
    return np.linalg.inv(A)


def load_day_zne(sta, year, doy, idx=None):
    """(t0, (3, n) float64 ZNE) from the frozen 20 Hz cache, rotated with the StationXML azimuths.
    Z-only helpers come back as (1, n) with a flag."""
    idx = cache_index(years=[year], stas=[sta]) if idx is None else idx
    p = idx.get((sta, year, doy))
    if p is None:
        return None
    z = np.load(p)
    x = z["x"]
    t0 = float(z["t0"])
    zs = zsign(sta)
    if x.shape[0] == 3:
        x = np.vstack([zs * x[0:1], (rot_matrix(sta) @ x[1:3].astype(np.float64))]).astype(np.float32)
    elif zs != 1.0:
        x = (zs * x).astype(np.float32)
    return t0, np.ascontiguousarray(x, dtype=np.float32)


class DayCache:
    """LRU over load_day_zne, keyed (sta, year, doy).  Keeps the current family/station pass hot."""

    def __init__(self, maxn=80, years=None, stas=None):
        from collections import OrderedDict
        self.idx = cache_index(years=years, stas=stas)
        self.d = OrderedDict()
        self.maxn = maxn

    def get(self, sta, year, doy):
        k = (sta, year, doy)
        if k in self.d:
            self.d.move_to_end(k)
            return self.d[k]
        v = load_day_zne(sta, year, doy, self.idx)
        self.d[k] = v
        if len(self.d) > self.maxn:
            self.d.popitem(last=False)
        return v

    def has_any(self, sta):
        return any(k[0] == sta for k in self.idx)


# ---------------------------------------------------------------- frozen primitives -------------
# COPIES of e5r_templates.cut_windows / align_stack / envelope / pick_s / s_snr_db.  Copied, not
# imported, because e5r_templates.py hard-codes the lfe_catalog_repair cache root at import time.
def yd_of(t):
    """(year, doy) arrays for epoch seconds, without touching pandas (hot path)."""
    dt = np.asarray(t, float).astype("datetime64[s]")
    yr = dt.astype("datetime64[Y]")
    return yr.astype(int) + 1970, (dt.astype("datetime64[D]") - yr).astype(int) + 1


def cut_fast(dc, sta, starts, n, norm=True):
    """(N,3,n) float32 windows at absolute start times + kept indices.  Day-sorted, pandas-free.

    Same cut and same rms normalisation as e5r_templates.cut_windows; only the bookkeeping differs."""
    t = np.asarray(starts, float)
    yr, doy = yd_of(t)
    order = np.lexsort((doy, yr))
    segs, keep = [], []
    cur, X, t0 = None, None, None
    for j in order:
        k = (int(yr[j]), int(doy[j]))
        if k != cur:
            v = dc.get(sta, k[0], k[1])
            cur = k
            X, t0 = (None, None) if v is None else (v[1], v[0])
        if X is None:
            continue
        i0 = int(round((t[j] - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            continue
        seg = X[:, i0:i0 + n]
        if norm:
            r = float(np.sqrt((seg.astype(np.float64) ** 2).mean()))
            if not np.isfinite(r) or r <= 0:
                continue
            seg = seg / r
        segs.append(np.asarray(seg, np.float32))
        keep.append(j)
    if not segs:
        return None, np.zeros(0, int)
    o = np.argsort(keep)
    return np.stack(segs)[o], np.asarray(keep, int)[o]


def cut_windows_abs(dc, sta, starts, n):
    """(N,3,n) rms-normalised windows at absolute start times, plus the kept indices."""
    segs, keep = [], []
    ts_ser = pd.to_datetime(np.asarray(starts, float), unit="s")
    for j, ts in enumerate(ts_ser):
        v = dc.get(sta, int(ts.year), int(ts.dayofyear))
        if v is None:
            continue
        t0, X = v
        i0 = int(round((starts[j] - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            continue
        seg = X[:, i0:i0 + n]
        r = np.sqrt((seg ** 2).mean())
        if not np.isfinite(r) or r <= 0:
            continue
        segs.append(seg / r)
        keep.append(j)
    return (np.stack(segs), np.asarray(keep)) if segs else (None, np.zeros(0, int))


def align_stack(S, iters=2, lag=30):
    """FROZEN (e5r_templates.align_stack, verbatim).  lag=30 samples = +-1.5 s at 20 Hz."""
    N, nc, n = S.shape
    tot = np.zeros(N, int)
    st = S.mean(0)
    ar = np.arange(n)
    for _ in range(iters):
        ref = st - st.mean(1, keepdims=True)
        ccs = np.empty((2 * lag + 1, N))
        for kk, k in enumerate(range(-lag, lag + 1)):
            ccs[kk] = np.einsum("ijk,jk->i", S, np.roll(ref, -k, axis=1))
        best = ccs.argmax(0) - lag
        idx = (ar[None, :] - best[:, None]) % n
        S = S[np.arange(N)[:, None, None], np.arange(nc)[None, :, None], idx[:, None, :]]
        tot += best
        st = S.mean(0)
    return st, S, tot


def envelope(st, smooth_s=0.5):
    from scipy.signal import hilbert
    env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    k = int(round(smooth_s * FS))
    return np.convolve(env, np.ones(k) / k, mode="same")


def pick_s(st, nominal_idx, search_s=4.0):
    env = envelope(st)
    w = int(round(search_s * FS))
    a, b = max(nominal_idx - w, 0), min(nominal_idx + w, len(env))
    return int(a + np.argmax(env[a:b]))


def s_snr_db(st, s_idx, win_s=4.0, noise_off_s=15.0):
    n = int(round(win_s * FS))
    a = max(s_idx - int(0.5 * FS), 0)
    sig = st[:, a:a + n]
    b = max(a - int(noise_off_s * FS), 0)
    noi = st[:, b:b + n]
    if sig.size == 0 or noi.size == 0 or (noi ** 2).mean() <= 0:
        return np.nan
    return float(10 * np.log10((sig ** 2).mean() / (noi ** 2).mean()))


def s_excess_db(st, s_idx, sig=(-1.0, 5.0), noi=(-12.0, -6.0)):
    """Gate-1 statistic: 10log10 E(S-1..S+5) / E(S-12..S-6) on a 3-C stack, NO alignment."""
    def w(a, b):
        i, j = s_idx + int(round(a * FS)), s_idx + int(round(b * FS))
        i, j = max(i, 0), min(j, st.shape[1])
        return st[:, i:j]
    a, b = w(*sig), w(*noi)
    if a.size == 0 or b.size == 0 or (b ** 2).mean() <= 0:
        return np.nan
    return float(10 * np.log10((a ** 2).mean() / (b ** 2).mean()))


# ---------------------------------------------------------------- travel times ------------------
def station_corr(stas):
    """a5 S-wave station corrections; NEW stations have none -> 0.0, and the +-1.5 s constant-offset
    search of the acceptance test is charged identically to every null (Merlin flaw 8)."""
    c = pd.read_csv(f"{EXT}/results/a5_station_corrections.csv").set_index("sta")
    out = {}
    for s in stas:
        v = c.corr_S.get(s, np.nan)
        out[s] = 0.0 if not np.isfinite(v) else float(v)
    return out


def fam_tt(stas, ST=None):
    """{(bkey, sta): predicted S delay after origin} = Crosson 1-D + a5 correction (0 for new)."""
    import p_vel1d as V
    ST = station_table() if ST is None else ST
    F = families()
    CR = station_corr(stas)
    out = {}
    for b in F.index:
        for s in stas:
            d = float(epi_km(F.lat[b], F.lon[b], ST.lat[s], ST.lon[s]))
            out[(int(b), s)] = float(V.tt(d, float(F.dep[b]), "S")) + CR[s]
    return out


def sha256_of(obj):
    import hashlib
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()
