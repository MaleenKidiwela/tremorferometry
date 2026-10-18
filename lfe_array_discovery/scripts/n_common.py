#!/usr/bin/env python
"""Puget / S-Olympic 2008 (and the 2007-01 confirmation) bed adapter for agent K's LITERAL
Bostock chain.  Prefix `n_`.

Same design as `l_common.py` (the N-Olympic 2010 adapter): NOTHING in `scripts/k_*.py` is edited.
Every `k_` module reads its bed facts through `k_common` at CALL time, so this module rebinds them
IN PROCESS and the k_ stages then run unchanged on the 2008 / 2007 data.  The recipe constants
(6 s windows, 1-8 Hz, plain per-component mean, 8 x MAD, 4 s decluster, Bungum-Husebye
realignment, 12 s templates, no merge) live in k_auto / k_group / k_grow and are NOT touched.

WHAT IS DIFFERENT FROM THE 2010 BED, AND WHY

 1. THERE IS NO TREMOR CATALOGUE.  PNSN starts 2010-01-06, so the 2010 chain's "PNSN candidate
    hour" rule cannot be evaluated at all.  The hours come from `n_hours.py`'s per-station
    tremor-hour audit, written to results/n_puget2008/tremor_hours_<bed>.csv, and this module only
    READS that file.  (The 2008 network-energy hour flag that `d_probe.py` used is 86 % daytime --
    `results/WA_V2_INPUTS.md` sec.2.2 -- and is NOT used for discovery here.)
 2. THE SUB-REGION is Puget/S-WA + S-Olympic (lat 46.2-47.8) instead of N-Olympic.  Active cells
    are still Merlin's definition (Slab2 28-45 km, >= 30 PNSN 2010-2011 epicentres within 8 km) and
    are therefore 2010-2011 cells used for a 2008 bed -- stated, not hidden (d_coverage.py sec.0).
 3. QUIET DAYS are real calendar days outside the window, not hours inside it:
      2008 -> the 20 literature-quiet days of results/hour_energy_flag_2008.csv (2008-02/03 and
              2008-07), which are cached (tranche T4 spans 2008-02-01..2008-07-31);
      2007 -> the 15 out-of-episode days downloaded with tranche T5 (2006-11-xx, 2007-04-xx).
    They are the reference for the TEMPLATE-LEVEL quiet-day false-alarm clause, which is the
    screen the 2010/2011 experience says is the one that matters (PLAN.md sec.15 RESULT).
 4. THE SITE-SHIFT NULL shifts days WITHIN a partition (window days among window days, quiet days
    among quiet days), because the quiet days are calendar-remote from the window and a single
    modulo list would swap the two populations.

Environment
  N_BED=2008|2007   which bed (default 2008)
  N_NULL_SHIFT=1    activate the day-shift null (offsets from N_NULL_SEED, default 11)
  N_OUT             output directory
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
CACHE = f"{STUDY}/cache20hz"
ROOT = "/home/jovyan/tremorferometry"
PICKER = f"{ROOT}/lfe_picker"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
CHESTLER = (f"{ROOT}/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")
sys.path.insert(0, HERE)
sys.path.insert(0, f"{PICKER}/scripts")

BED = os.environ.get("N_BED", "2008")
OUTROOT = os.environ.get("N_OUTROOT", f"{RES}/n_puget2008")

# ------------------------------------------------------------------ the two beds ---------------
SUBREGION = "Puget/S-Olympic"
LAT0, LAT1, LON0, LON1 = 46.2, 47.8, -125.0, -121.5     # sub-region box for the active cells
CELL_DEP = (28.0, 45.0)
CELL_MIN_EPI = 30
CELL_KM = 8.0
CELL_YEARS = (2010, 2011)
PROBE = f"{RES}/site_probe_2008.csv"                     # the 2008 array-era site probe
# 2008-03-11 is flagged is_quiet by the array-energy flag but carries 7 hot hours -- "probably
# real tremor on a day the literature calls quiet", and results/P_STEP2_SITES.md sec.3.2 already
# excluded it from the 2008 probe baseline.  The same exclusion is carried here, for the same
# reason and with no new rule: a day that may be a tremor day is neither a baseline day nor a
# false-alarm reference day.
QUIET_EXCLUDE = {"2008-03-11"}
VETO_STATION_MAX = 0.30                                  # WA_V2_INPUTS sec.3

if BED == "2008":
    YEAR, DAY0, DAY1 = 2008, "2008-04-15", "2008-06-10"
    TRANCHE = "T4"
    VETO_YEARS = (2008,)
elif BED == "2007":
    YEAR, DAY0, DAY1 = 2007, "2007-01-01", "2007-02-28"
    TRANCHE = "T5"
    VETO_YEARS = (2006, 2007)
else:
    raise SystemExit(f"N_BED must be 2008 or 2007, not {BED!r}")

OUT = os.environ.get("N_OUT", f"{OUTROOT}/{BED}")
HOURS = f"{OUTROOT}/tremor_hours_{BED}.csv"
QUIET_DAYS_TXT = f"{OUTROOT}/quiet_dates_{BED}.txt"

_S = {}


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------ days -----------------------
def window_days():
    return [d.strftime("%Y-%m-%d") for d in pd.date_range(DAY0, DAY1)]


def quiet_days():
    """Literature-quiet calendar days for this bed, intersected with what is CACHED."""
    if "QD" in _S:
        return _S["QD"]
    if BED == "2008":
        F = pd.read_csv(f"{RES}/hour_energy_flag_2008.csv")
        q = sorted(F[F.is_quiet].date)
    else:
        C = pd.read_parquet(f"{RES}/cache_{TRANCHE}.parquet")
        C = C[C.status == "built"]
        w = set(window_days())
        q = sorted(d for d in set(C.date) if d not in w)
    have = cached_days()
    _S["QD"] = [d for d in q if have.get(d, 0) >= 4 and d not in QUIET_EXCLUDE]
    return _S["QD"]


def cached_days():
    """{date: n cached stations} over the ADMITTED station set (all years on disk)."""
    if "CD" in _S:
        return _S["CD"]
    n = {}
    for (sta, yr, doy) in cache_index():
        d = (pd.Timestamp(year=yr, month=1, day=1)
             + pd.Timedelta(days=doy - 1)).strftime("%Y-%m-%d")
        n[d] = n.get(d, 0) + 1
    _S["CD"] = n
    return n


def bed_days():
    """Every day this bed may read: the window plus the literature-quiet days."""
    return sorted(set(window_days()) | set(quiet_days()))


# ------------------------------------------------------------------ stations -------------------
def _veto_frac():
    f = {}
    for y in VETO_YEARS:
        p = f"{RES}/wa_station_hour_veto_{y}.csv"
        if os.path.exists(p):
            V = pd.read_csv(p, usecols=["netsta", "veto"])
            g = V.groupby("netsta").veto.agg(["sum", "size"])
            for ns, r in g.iterrows():
                a, b = f.get(ns, (0, 0))
                f[ns] = (a + int(r["sum"]), b + int(r["size"]))
    return {k: v[0] / max(v[1], 1) for k, v in f.items()}


def cached_netsta():
    """netsta codes with at least one cached day inside this bed's day set."""
    if "CN" in _S:
        return _S["CN"]
    days = set(window_days())
    out = set()
    for p in glob.glob(f"{CACHE}/*/*/*.npz"):
        parts = p.split(os.sep)
        yr, doy = int(parts[-2]), int(parts[-1][:3])
        d = (pd.Timestamp(year=yr, month=1, day=1)
             + pd.Timedelta(days=doy - 1)).strftime("%Y-%m-%d")
        if d in days:
            out.add(parts[-3])
    _S["CN"] = out
    return out


def admitted():
    """DataFrame(net, sta, netsta, veto_frac, in_bed, health_drop) of the 2008 probe-admitted
    stations.  For the 2007 bed the SAME admitted set is used -- no station is re-probed and none
    is substituted; a station simply drops out if it has no 2007 data or fails the 2007 health
    bar.  (The 2008 probe is the only array-era probe there is for this array; re-probing on the
    confirmation bed would be tuning on the confirmation data.)"""
    if "A" in _S:
        return _S["A"]
    A = pd.read_csv(PROBE)
    A = A[A.admit].copy()
    # OPTIONAL station-set restriction (agent P's quiet subset, 2026-09-21).  Unset -> no effect,
    # so every agent-N number stays reproducible; set to a comma list of NET.STA (or a file path
    # holding one) it restricts the bed to that subset and NOTHING ELSE changes -- same probe,
    # same veto, same audited hours, same quiet days, same recipe constants.
    sub = os.environ.get("N_STATION_SUBSET", "").strip()
    if sub:
        if os.path.exists(sub):
            sub = open(sub).read()
        want = {x.strip() for x in sub.replace("\n", ",").split(",") if x.strip()}
        miss = want - set(A.netsta)
        if miss:
            raise SystemExit(f"N_STATION_SUBSET names non-admitted stations: {sorted(miss)}")
        A = A[A.netsta.isin(want)].copy()
    f = _veto_frac()
    A["veto_frac"] = [float(f.get(n, 0.0)) for n in A.netsta]
    have = cached_netsta()
    A["in_bed"] = [n in have for n in A.netsta]
    A["health_drop"] = A.veto_frac > VETO_STATION_MAX
    _S["A"] = A.sort_values("netsta").reset_index(drop=True)
    return _S["A"]


def keepers():
    A = admitted()
    return A[A.in_bed & ~A.health_drop]


def stations():
    return list(keepers().sta)


def netsta():
    return {r.sta: r.netsta for r in keepers().itertuples()}


# ------------------------------------------------------------------ station table ---------------
def station_table():
    if "ST" in _S:
        return _S["ST"]
    C = pd.read_csv(f"{RES}/channels_2006_2012.csv").dropna(subset=["lat", "lon"])
    ns = netsta()
    C["netsta"] = C.net + "." + C.sta
    C = C[C.netsta.isin(set(ns.values()))]
    g = (C.groupby(["net", "sta"], as_index=False)
           .agg(lat=("lat", "median"), lon=("lon", "median"), elev_m=("elev", "median"),
                fs_native=("fs", "max"), band=("band", "first")))
    g["netsta"] = g.net + "." + g.sta
    miss = sorted(set(ns) - set(g.sta))
    if miss:
        raise SystemExit(f"no channel metadata for {miss}")
    _S["ST"] = g.set_index("sta")
    return _S["ST"]


def az12_wa():
    """{sta: (az1, az2)} for stations cached as (?HZ, ?H1, ?H2) -- the PB boreholes."""
    C = pd.read_csv(f"{RES}/channels_2006_2012.csv")
    ns = netsta()
    C["netsta"] = C.net + "." + C.sta
    C = C[C.netsta.isin(set(ns.values()))]
    out = {}
    for s in ns:
        c = C[(C.sta == s) & (C.chan.str.len() == 3)]
        a1 = c[c.chan.str.endswith("1")]
        a2 = c[c.chan.str.endswith("2")]
        ah = a1[a1.chan.str[:2].isin(["EH", "BH", "HH", "SH"])]
        bh = a2[a2.chan.str[:2].isin(["EH", "BH", "HH", "SH"])]
        if len(ah) and len(bh):
            out[s] = (float(ah.azimuth.iloc[0]), float(bh.azimuth.iloc[0]))
    return out


# ------------------------------------------------------------------ cache index -----------------
def cache_index():
    """{(sta, year, doy): path} built from the EXPLICIT (net, sta) list and this study's cache
    root only (the WA cache holds C8.PA02 beside XG.PA02)."""
    if "IDX" in _S:
        return _S["IDX"]
    idx = {}
    for s, ns in netsta().items():
        for p in glob.glob(f"{CACHE}/{ns}/*/*.npz"):
            parts = p.split(os.sep)
            idx[(s, int(parts[-2]), int(parts[-1][:3]))] = p
    _S["IDX"] = idx
    return idx


# ------------------------------------------------------------------ health veto -----------------
def health_veto():
    """{(date, hour): set of vetoed bare station codes} over the admitted set."""
    if "V" in _S:
        return _S["V"]
    ns = set(netsta().values())
    d = {}
    for y in VETO_YEARS:
        p = f"{RES}/wa_station_hour_veto_{y}.csv"
        if not os.path.exists(p):
            continue
        V = pd.read_csv(p, usecols=["netsta", "date", "hour", "veto"])
        V = V[V.veto & V.netsta.isin(ns)]
        for r in V.itertuples():
            d.setdefault((r.date, int(r.hour)), set()).add(r.netsta.split(".")[1])
    _S["V"] = d
    return d


# ------------------------------------------------------------------ hours ----------------------
def candidate_hours():
    """The AUDITED tremor hours (n_hours.py).  This file is the only hour source for discovery."""
    if not os.path.exists(HOURS):
        raise SystemExit(f"missing {HOURS} -- run scripts/n_hours.py first")
    H = pd.read_csv(HOURS)
    return H[["date", "hour"]].copy()


def tremor_hour_set():
    return {(str(r.date), int(r.hour)) for r in candidate_hours().itertuples()}


def tremor_dates():
    return sorted(set(candidate_hours().date))


def quiet_hours():
    """Every hour of every literature-quiet day (the quiet CONTRAST for the pair stage uses the
    hour-matched set of results/wa_quiet_hours_2008.csv when it exists for this bed)."""
    p = f"{RES}/wa_quiet_hours_{YEAR}.csv"
    if os.path.exists(p):
        Q = pd.read_csv(p)[["date", "hour"]]
        have = set(quiet_days())
        Q = Q[Q.date.isin(have)]
        if len(Q):
            return Q.reset_index(drop=True)
    rows = [dict(date=d, hour=h) for d in quiet_days() for h in range(0, 24, 4)]
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ geometry -------------------
def active_cells():
    """Merlin's active cells inside the Puget/S-Olympic box.  NOTE: the epicentre count is from
    the 2010-2011 PNSN catalogue because no 2008 catalogue exists -- if the 2008 ETS lit up
    different cells these are the wrong cells, which is stated wherever they are used."""
    if "AC" in _S:
        return _S["AC"]
    import g_common as G
    import s_common as S
    tr = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce").astype("datetime64[ns]")
    tr = tr.dropna(subset=["time"])
    tr = tr[(tr.time >= pd.Timestamp(f"{CELL_YEARS[0]}-01-01"))
            & (tr.time < pd.Timestamp(f"{CELL_YEARS[1] + 1}-01-01"))]
    step = 0.05
    lat = np.arange(LAT0, LAT1 + 1e-9, step)
    lon = np.arange(LON0, LON1 + 1e-9, step)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    LA, LO = LA.ravel(), LO.ravel()
    Z = G.slab2_depth(LA, LO)
    ok = np.isfinite(Z) & (Z >= CELL_DEP[0]) & (Z <= CELL_DEP[1])
    LA, LO, Z = LA[ok], LO[ok], Z[ok]
    tla, tlo = tr.lat.to_numpy(), tr.lon.to_numpy()
    cnt = np.zeros(len(LA), int)
    for i in range(0, len(LA), 200):
        sl = slice(i, min(i + 200, len(LA)))
        d = S.epi_km(LA[sl][:, None], LO[sl][:, None], tla[None, :], tlo[None, :])
        cnt[sl] = (d <= CELL_KM).sum(1)
    keep = cnt >= CELL_MIN_EPI
    _S["AC"] = pd.DataFrame(dict(lat=LA[keep], lon=LO[keep], dep=Z[keep], n_epi=cnt[keep]))
    return _S["AC"]


def near_stations(km=60.0):
    """Admitted stations within `km` EPICENTRAL of any active cell -- the audit's station set."""
    import s_common as S
    ST = station_table()
    AC = active_cells()
    cla, clo = AC.lat.to_numpy(), AC.lon.to_numpy()
    out = []
    for s in stations():
        d = float(np.min(S.epi_km(float(ST.lat[s]), float(ST.lon[s]), cla, clo)))
        if d <= km:
            out.append(s)
    return sorted(out)


def station_corrections():
    """Intraslab-earthquake S station terms (results/x_station_corrections.csv).  A station
    without a usable term gets 0.0 and is FLAGGED -- never silently imputed."""
    p = f"{RES}/x_station_corrections.csv"
    corr, flagged = {}, []
    if os.path.exists(p):
        C = pd.read_csv(p).set_index("sta")
        for s in stations():
            if s in C.index and bool(C.loc[s, "usable"]):
                v = float(C.loc[s, "corr_S"])
                corr[s] = v if np.isfinite(v) else 0.0
            else:
                corr[s] = 0.0
                flagged.append(s)
    else:
        corr = {s: 0.0 for s in stations()}
        flagged = list(stations())
    return corr, flagged


# ------------------------------------------------------------------ the null -------------------
def null_offsets(stas, seed=None):
    seed = int(os.environ.get("N_NULL_SEED", "11")) if seed is None else seed
    rng = np.random.default_rng(seed)
    nw = len(window_days())
    return {s: int(rng.integers(1, max(nw, 2))) for s in sorted(stas)}


# ------------------------------------------------------------------ activation -----------------
def activate(null=None):
    """Rebind k_common onto this bed.  Idempotent.  Returns the k_common module."""
    import k_common as K
    if _S.get("ACT"):
        return K
    import s_common as S
    ST = station_table()
    K._C["ST"] = ST
    S.station_table = lambda: ST
    S.AZ12.update(az12_wa())
    for s in list(S.AZ12):
        S.ZSIGN.setdefault(s, 1.0)
    K._IDX.clear()
    K._IDX.update(cache_index())
    K.BED_YEAR, K.DAY0, K.DAY1 = YEAR, DAY0, DAY1
    K.bed_days = bed_days
    K.CLEAN11 = stations()
    K.DET_SITES = stations()

    null = (os.environ.get("N_NULL_SHIFT", "") == "1") if null is None else null
    if null:
        wd, qd = window_days(), quiet_days()
        part = {}
        for i, d in enumerate(wd):
            part[d] = ("w", i)
        for i, d in enumerate(qd):
            part[d] = ("q", i)
        lists = {"w": wd, "q": qd}
        off = null_offsets(stations())
        _S["OFF"] = off
        base = K.load_day

        def _shifted(sta, date):
            """The station's data from another day of the SAME partition, restamped with day d's
            t0 -- every downstream stage is byte-identical code on a network whose inter-station
            coherence has been destroyed."""
            if date not in part or sta not in off:
                return base(sta, date)
            kind, i = part[date]
            L = lists[kind]
            if len(L) < 2:
                return base(sta, date)
            src = L[(i + off[sta]) % len(L)]
            v = base(sta, src)
            if v is None:
                return None
            return float(pd.Timestamp(date).timestamp()), v[1]

        K.load_day = _shifted
        log(f"[null] day-shift active (within-partition), "
            f"{json.dumps({k: off[k] for k in list(off)[:6]})} ...")
    _S["ACT"] = True
    return K


def bed_summary():
    A = admitted()
    keep = keepers()
    return dict(bed=BED, subregion=SUBREGION, day0=DAY0, day1=DAY1, tranche=TRANCHE,
                n_probe_admitted=int(len(A)),
                n_absent=int((~A.in_bed).sum()), absent=sorted(A[~A.in_bed].netsta),
                n_health_dropped=int((A.in_bed & A.health_drop).sum()),
                health_dropped=sorted(A[A.in_bed & A.health_drop].netsta),
                n_stations=int(len(keep)),
                by_net={k: int(v) for k, v in keep.groupby("net").size().items()},
                stations={k: sorted(v) for k, v in keep.groupby("net").sta.apply(list).items()},
                n_window_days=len(window_days()),
                n_quiet_days=len(quiet_days()), quiet_days=quiet_days(),
                n_active_cells=int(len(active_cells())),
                n_near_stations=len(near_stations()))
