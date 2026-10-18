#!/usr/bin/env python
"""lfe_array_discovery -- WASHINGTON bed adapter for the frozen `w_` array-discovery engine.

Prefix `x_`.  NOTHING in `scripts/w_*.py` and nothing in lfe_picker is edited: every `w_` script
imports `w_common as W` and reads `W.<attr>` at CALL time, so this module rebinds W's bed-level
attributes IN PROCESS (exactly the way `w_common.patch_v` rebinds `v_common`) and then the `w_`
stage modules run unchanged on the Washington data.

WHY AN ADAPTER AND NOT A FIFTH BED MODE IN w_common
  `results/freeze_step1{,b,c}.json` hash `w_common.py`, and the Step-1b / Step-1c queues re-import
  it at every stage while they are still running.  Editing that file would (a) break their frozen
  hashes and (b) change the behaviour of a run in flight.  The engine constants are therefore
  BYTE-IDENTICAL to `freeze_step1c.json` -- this module changes only the BED (year, window, days,
  hours, sites, withheld, grid centre, active-cell depth band) and the two things Merlin's Q4/Q5
  make region-specific (the location minimum-site bar and the HF-veto station).

  The only engine constants that differ from Step 1c, each with its written reason:
    GRID_LAT0/GRID_LON0    48.50/-123.75 -> 48.00/-123.20   Merlin Q4: the WA grid centre.
    ACTIVE_CELL_DEP        (25,50) -> (28,45)               Merlin F2/Step 0: the WA active-cell
                                                            definition (Slab2 28-45 km).
    LOC_MIN_SITES          4 -> 5                           Merlin Q4 "located" = >= 5 picked
                                                            sites <= 60 km (w_common's own comment
                                                            says "the WA bar is >= 5").
  Everything else -- AC_WIN_S, AC_HOP_S, AC_MIN_LAG_S, AC_DECLUSTER_S, the site-z clauses, MIN_SITES
  4, T_ENG_FLOOR 5.0, the FA rules, the seed-cluster rule, survival, the scramble ratio, dedupe,
  the location controls, the null design -- is copied unchanged and verified byte-for-byte against
  freeze_step1c.json by `x_freeze.check_constants()`.

BED FACTS THAT FORCED A WRITTEN CHOICE (all recorded in results/freeze_wa_v1.json)
  * PB boreholes are cached as (EHZ, EH1, EH2): the horizontals are NOT rotated to N/E.  Every
    statistic this chain computes is rotation-INVARIANT in the horizontal plane -- the per-station
    autocorrelation is a 3-C mean cc of one station against itself (the same rotation on both
    sides), the growth/Gate-1 `sdb` and the HF ratio sum the squares of all three components, and
    the Gate-2 S pick is made on sqrt(H1^2 + H2^2).  `w_common.load_day` nevertheless REJECTS any
    file whose channel order is not exactly Z,N,E, which would silently drop all 19 PB boreholes.
    `load_day` is therefore replaced with one that accepts (vertical, horizontal, horizontal) in
    that order and records the convention.  No rotation is applied and none is needed.
  * The site id produced by the 1-km merge is a data-driven tie-break and is NOT stable across
    cache roots (the 2026-09-20 Step-1b bug).  The WA bed therefore carries an EXPLICIT site table
    (site key -> member stations) in results/step5_bed.json; nothing here matches on a merge id.
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
# X_CACHE lets the Step-5b INJECTION run point the whole chain (autocorrelation, growth, gates)
# at the injected cache root without touching the real one.  Unset everywhere else.
CACHE = os.environ.get("X_CACHE") or f"{STUDY}/cache20hz"
TRACES = f"{STUDY}/traces"
TRACES_EQ = f"{STUDY}/traces_eq"
OUT = f"{RES}/step5"
BED_JSON = f"{RES}/step5_bed.json"
CORR_CSV = f"{RES}/x_station_corrections.csv"
FREEZE = f"{RES}/freeze_wa_v1.json"
PICKER = "/home/jovyan/tremorferometry/lfe_picker"
ROOT = "/home/jovyan/tremorferometry"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
CHESTLER = (f"{ROOT}/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")

sys.path.insert(0, HERE)
sys.path.insert(0, f"{PICKER}/scripts")

# --------------------------------------------------------------------- WA sub-region and window --
SUBREGION = "N-Olympic"
LAT0, LAT1 = 47.8, 48.3
LON0, LON1 = -123.7, -122.7
YEAR = 2010
DAY0, DAY1 = "2010-07-20", "2010-09-10"
GRID_LAT0, GRID_LON0, GRID_HALFWIDTH = 48.00, -123.20, 62.0
CELL_DEP = (28.0, 45.0)                 # Merlin's WA active-cell depth band
CELL_MIN_EPI = 30                       # >= 30 PNSN epicentres within 8 km ...
CELL_KM = 8.0
CELL_YEARS = (2010, 2011)               # ... in 2010-2011
LOC_MIN_SITES_WA = 5                    # Merlin Q4: "located" needs >= 5 picked sites <= 60 km
HF_STATION_CANDIDATES = ("B003", "B004")   # WA 100-Hz boreholes; the probe-admitted one is used

# Withheld from EVERY autocorrelation sum, growth sum and stack, for the whole run
# (P_STEP2_SITES.md sec.5 Option B, adopted in the PLAN sec.10 addendum).
WITHHELD_STATIONS = ["B001", "B005", "B006", "B007", "LRIV", "VGZ"]
WITHHELD_SITES = {"B001": ["B001"], "B005q": ["B005", "B006", "B007", "LRIV"], "VGZ": ["VGZ"]}

_CFG = [None]
_ACT = [False]


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%Y-%m-%d %H:%M:%S')}]", *a, flush=True)


# =================================================================== the WA station table ========
_SSTA = [None]


def wa_station_table():
    """One row per station of the WA box with a 2006-2012 epoch: lat, lon, netsta, net, band, fs.

    Built from `results/channels_2006_2012.csv` (the tranche-1 channel-level StationXML), which is
    the only station metadata this study uses.  Indexed by `sta`, because every w_/v_/g_ helper
    addresses stations by bare code."""
    if _SSTA[0] is None:
        C = pd.read_csv(f"{RES}/channels_2006_2012.csv")
        C = C.dropna(subset=["lat", "lon"])
        g = (C.groupby(["net", "sta"], as_index=False)
               .agg(lat=("lat", "median"), lon=("lon", "median"), elev_m=("elev", "median"),
                    fs_native=("fs", "max"), band=("band", "first")))
        g["netsta"] = g.net + "." + g.sta
        g["klass"] = "wa"
        _SSTA[0] = g.drop_duplicates("sta", keep="first").set_index("sta")
    return _SSTA[0]


def cfg():
    if _CFG[0] is None:
        if not os.path.exists(BED_JSON):
            raise SystemExit(f"missing {BED_JSON} -- run scripts/x_setup.py first")
        _CFG[0] = json.load(open(BED_JSON))
    return _CFG[0]


def station_corrections():
    """{sta: (corr_S, corr_P)} from `x_corrections.py` (intraslab EQ station terms, a5 recipe).

    A station with fewer than the frozen minimum accepted picks, or a MAD above the frozen bar,
    gets 0.0 and is FLAGGED (`usable=False`) -- never silently imputed."""
    out = {}
    if os.path.exists(CORR_CSV):
        C = pd.read_csv(CORR_CSV).set_index("sta")
        for s in C.index:
            if not bool(C.loc[s, "usable"]):
                out[str(s)] = (0.0, 0.0)
                continue
            cs, cp = float(C.loc[s, "corr_S"]), float(C.loc[s, "corr_P"])
            out[str(s)] = (0.0 if not np.isfinite(cs) else cs, 0.0 if not np.isfinite(cp) else cp)
    return out


# =================================================================== active cells / hours ========
def active_cells(sub=True):
    """Merlin's WA active cells: Slab2 interface points inside the sub-region box, depth 28-45 km,
    with >= 30 PNSN 2010-2011 epicentres within 8 km.  Catalogue-driven only."""
    import g_common as G
    import s_common as S
    tr = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce").astype("datetime64[ns]")
    tr = tr.dropna(subset=["time"])
    tr = tr[(tr.time >= pd.Timestamp(f"{CELL_YEARS[0]}-01-01"))
            & (tr.time < pd.Timestamp(f"{CELL_YEARS[1] + 1}-01-01"))]
    step = 0.05
    if sub:
        lat = np.arange(LAT0, LAT1 + 1e-9, step)
        lon = np.arange(LON0, LON1 + 1e-9, step)
    else:
        half = GRID_HALFWIDTH / 111.195
        lat = np.arange(GRID_LAT0 - half, GRID_LAT0 + half + 1e-9, step)
        lon = np.arange(GRID_LON0 - half / np.cos(np.radians(GRID_LAT0)),
                        GRID_LON0 + half / np.cos(np.radians(GRID_LAT0)) + 1e-9, step)
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
    return pd.DataFrame(dict(lat=LA[keep], lon=LO[keep], dep=Z[keep], n_epi=cnt[keep]))


def tremor_days():
    """PNSN-corrected N-Olympic tremor days inside the window: days with >= 5 epicentres within
    AC_HOUR_EPI_KM of an active cell (the same clause the hour rule uses, aggregated to days)."""
    import w_common as W
    H = W.candidate_hours()
    return sorted(set(H.date)) if len(H) else []


# =================================================================== the bed patch ===============
def _wa_load_day(sta, year, doy, idx=None):
    """(t0, (3, n) float32) for a WA station-day, accepting UNROTATED borehole horizontals.

    `w_common.load_day` demands order == [*Z, *N, *E]; the PB boreholes are cached as
    (EHZ, EH1, EH2) and would be dropped from every hour.  Every statistic downstream is invariant
    to a rotation about the vertical (see the module docstring), so the requirement here is only
    "vertical first, two horizontals after"."""
    import w_common as W
    idx = W.cache_idx() if idx is None else idx
    p = idx.get((sta, int(year), int(doy)))
    if p is None:
        return None
    try:
        z = np.load(p)
    except Exception:                                                       # noqa: BLE001
        return None
    x = z["x"]
    if x.shape[0] != 3:
        return None
    order = [str(o) for o in z["order"]]
    if order[0][-1] != "Z" or any(o[-1] == "Z" for o in order[1:]):
        return None
    return float(z["t0"]), np.ascontiguousarray(x, dtype=np.float32)


def activate():
    """Rebind `w_common` onto the Washington bed.  Idempotent."""
    if _ACT[0]:
        import w_common as W
        return W
    C = cfg()
    import w_common as W
    import s_common as S

    # ---- bed identity
    W.BED_MODE = "w"
    W.STEP1_YEAR = int(C["year"])
    W.STEP1_DAY0, W.STEP1_DAY1 = C["day0"], C["day1"]
    W.WITHHELD = tuple(C["withheld_stations"])
    W.BED = CACHE
    W.STUDY_CACHE = CACHE
    # ---- the three region-specific constants (see the module docstring)
    W.GRID_LAT0, W.GRID_LON0, W.GRID_HALFWIDTH = GRID_LAT0, GRID_LON0, GRID_HALFWIDTH
    W.ACTIVE_CELL_DEP = tuple(CELL_DEP)
    W.ACTIVE_CELL_MIN_EPI = CELL_MIN_EPI
    W.ACTIVE_CELL_KM = CELL_KM
    W.LOC_MIN_SITES = LOC_MIN_SITES_WA

    # ---- station table and corrections
    ST = wa_station_table()
    W._ST[0] = ST
    S.station_table = lambda: ST
    corr = station_corrections()
    S.station_corr = lambda stas: {s: float(corr.get(s, (0.0, 0.0))[0]) for s in stas}

    # ---- bed roots / reader
    W.bed_roots = lambda: [CACHE]
    W.station_table = lambda: ST
    W.extra_stations = lambda: ST.reset_index()[["net", "sta", "lat", "lon"]].assign(band="",
                                                                                    prov="")
    W.load_day = _wa_load_day

    # ---- the explicit site table (never a merge id; the 2026-09-20 instability bug)
    mem = {str(k): list(v) for k, v in C["site_members"].items()}
    order = list(C["site_order"])
    of = {s: q for q, mm in mem.items() for s in mm}
    det = [s for q in order for s in mem[q]]
    W._SITES[0] = (mem, of, order, det)
    W.sites = lambda year=None: W._SITES[0]

    # ---- cells / hours / quiet days
    W.active_cells = active_cells
    quiet = set(C["quiet_days"])

    def _quiet_days(year=None):
        A = W.availability(int(C["year"]))
        return pd.DataFrame([dict(date=r.date, n_sites=r.n_sites, days_from_tremor=np.nan,
                                  quiet=bool(r.date in quiet and r.usable))
                             for r in A.itertuples()])

    W.quiet_days = _quiet_days
    _ACT[0] = True
    return W


# =================================================================== gates patch =================
def gate_context(accepted, tpl, ccp95, withheld=None, controls=None, hf_station=None,
                 hf_veto=None):
    """A `g_gates.GateContext` bound to the WASHINGTON bed.

    `g_gates`/`g_common` are region-configurable through their module attributes only, so the
    Washington run rebinds them in process, exactly as `w_common.patch_v` does for `v_common`:

      G.CACHE_DIRS / G.TRACE_DIRS  -> this study's cache20hz / traces (Gate 1 stacks, HF ratio)
      G.SITE_MEMBERS / G.STA2SITE  -> the WA site table
      G.stations()                 -> the WA station table
      G.station_corrections()      -> the WA intraslab-EQ station terms
      ctx.accepted / ctx.TPL       -> the WA candidate list and its own iterated templates
      ctx.ccp95                    -> the WA bank's own cross-family cc p95 (w_dedupe writes it)
      cfg.g0_hf_station / g0_hf_veto -> the WA 100-Hz borehole and the cut derived on WA quiet days

    NOTE (Merlin Q4, accepted and recorded): `g_common.epi_km` and the Gate-1 moveout null bind
    cos(48.6) as a DEFAULT ARGUMENT at definition time, so the ~1 % east-west distance error at
    48.0 N cannot be patched away without editing lfe_picker.  It is left in place -- the same
    error is applied to the data and to every null -- and reported.
    """
    import g_common as G
    import g_gates as GG
    C = cfg()
    ST = wa_station_table()
    G.CACHE_DIRS = [CACHE]
    G.TRACE_DIRS = [TRACES]
    mem = {str(k): list(v) for k, v in C["site_members"].items()}
    allsites = dict(mem)
    for q, mm in WITHHELD_SITES.items():
        allsites[q] = list(mm)
    G.SITE_MEMBERS = allsites
    G.STA2SITE = {m: q for q, mm in allsites.items() for m in mm}
    G.stations = lambda path=None: ST
    corr = station_corrections()
    G.station_corrections = lambda path=None: corr
    wh = tuple(withheld if withheld is not None else C["withheld_stations"])
    cfg_g = GG.GateConfig(
        sites=tuple(list(mem) + [q for q in WITHHELD_SITES]),
        withheld=wh,
        grid_lat0=GRID_LAT0, grid_lon0=GRID_LON0, grid_halfwidth=GRID_HALFWIDTH,
        g0_hf_station=tuple([hf_station or C["hf_station"]]),
        g0_hf_veto=float(hf_veto if hf_veto is not None else C["hf_veto"]),
        g2_min_sites_60km=LOC_MIN_SITES_WA)
    ctx = GG.GateContext.__new__(GG.GateContext)
    ctx.cfg = cfg_g
    ctx.verbose = False
    ctx.S = ST
    ctx.corr = corr
    ctx.store = G.DayStore(cache_dirs=[CACHE])
    ctx.nstore = G.NativeStore(trace_dirs=[TRACES])
    ctx.site_ll = {q: (float(ST.lat[mm[0]]), float(ST.lon[mm[0]]))
                   for q, mm in allsites.items() if mm[0] in ST.index}
    ctx.site_members = {q: [m for m in mm if m in ST.index] for q, mm in allsites.items()}
    ctx.accepted = accepted
    ctx.TPL = tpl
    ctx.ccp95 = ccp95
    ctx._grid = None
    ctx._tremor = None
    ctx._controls = controls
    return ctx


# =================================================================== misc ========================
def sha256_file(p):
    import hashlib
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def x_scripts():
    return sorted(glob.glob(f"{HERE}/x_*.py")) + sorted(glob.glob(f"{HERE}/x_*.sh"))


def w_scripts():
    return sorted(glob.glob(f"{HERE}/w_*.py")) + sorted(glob.glob(f"{HERE}/w_*.sh"))


def input_tables():
    return [f"{RES}/{n}" for n in
            ("sites_array_era_final.csv", "site_probe_2010.csv", "coverage_by_subregion.csv",
             "withheld_candidates.csv", "withheld_proposal.json", "delivery_low_stations.csv",
             "channels_2006_2012.csv", "intraslab_eq_2006_2012.csv", "freeze_windows.json",
             "freeze_step1c.json")]
