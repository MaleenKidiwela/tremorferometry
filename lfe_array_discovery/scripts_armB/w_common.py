#!/usr/bin/env python
"""lfe_array_discovery -- ARRAY-DISCOVERY ENGINE, shared machinery.  Prefix `w_`.

Arm A of PLAN.md sec.10 (Merlin review 2026-09-20, `results/merlin_review_plan_2026-09-20.md`,
answer Q1/(a) and prescription Steps 0-1):

    network autocorrelation -> clusters -> align_stack handoff -> matched-filter growth ->
    survival (+ template-scramble null) -> iterations -> dedupe -> locate -> Gates 0-2,
    with a WHOLE-CHAIN site-day-shift null.

STEP 1 is the NORTH POSITIVE CONTROL of that engine: it runs on the frozen 20 Hz
`lfe_catalog_repair/cache20hz` bed (2012 south-Vancouver-Island ETS, 15 stations on disk),
with VGZ/SNB/B926 withheld exactly as in the Gates suite, and is scored by the relative-clause
rediscovery scorer (`q_rediscover`) against the 10 ISOLATED hidden Bostock families.

NOTHING in lfe_picker is edited.  `v_common` / `v_engine` / `v_iterate` / `v_dedupe` are imported
and their module-level constants are PATCHED IN PROCESS (patch_v below) to carry Merlin's four
engine prerequisites:

  P1  NO T = 8.0 floor.  FA-only thresholds with an ENGINEERING floor of 5.0, and a report of
      where the rule binds (v_common.T_FLOOR -> 5.0; every call passes floor=T_ENG_FLOOR).
  P2  DATA-DRIVEN SITE_MEMBERS: 1-km single-linkage merge of the station table (not a hard-coded
      7-site dict).
  P3  station_z median/MAD over VALID SAMPLES ONLY (the C-8 gap bug: the frozen station_z takes
      the median and the lower-half MAD over all NV samples, gaps included).
  P4  MIN_SITES 3 -> 4, and per-day availability NEVER intersected.

All constants below are written verbatim into results/freeze_step1.json by w_freeze.py.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
FIG = f"{STUDY}/figures"
LOGS = f"{STUDY}/logs"
WORK = f"{RES}/step1"
PICKER = "/home/jovyan/tremorferometry/lfe_picker"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
BED = f"{REP}/cache20hz"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
sys.path.insert(0, f"{PICKER}/scripts")
for _d in (RES, FIG, LOGS, WORK):
    os.makedirs(_d, exist_ok=True)

import s_common as S                                                        # noqa: E402

_ORIG_LOAD_DAY_ZNE = S.load_day_zne

FS = 20.0
NDAY = 1728000

# ------------------------------------------------------------------ Step-1 bed ------------------
STEP1_YEAR = 2012
STEP1_DAY0, STEP1_DAY1 = "2012-08-15", "2012-10-05"      # the 52 day-files on the bed
# ---- BED MODE.  "a" = the north PERMANENT bed (Step 1a, 6 sites at 46-137 km hypocentral).
# "b" = Step 1b, the same 2012 window with the CLOSE stations that exist then (C8 SOKB, CN LZB,
# the six C8 PA0x co-located with PGC, C8 BPCB, CN PFB, CN YOUB).  Selected by a written rule:
# the maximal subset of available 3-C sites whose MEDIAN S-moveout spread over the 130 families is
# <= AC_WIN_S (the frozen 10-s common window) -- 8.8 s for the six close sites, against 13.6 s for
# the Step-1a permanent six and 15.6 s for all eleven together.  NO engine constant is retuned.
BED_MODE = os.environ.get("W_BED", "a").lower()
WITHHELD_A = ("VGZ", "SNB", "B926")                      # the Gates-suite withheld set
# + one DENSE station withheld, chosen by coverage: dropping YOUB costs the fewest families with
# >= 5 close sites within 60 km epicentral (63/130, tied with PFB) and YOUB is the farther /
# least-covering of the two (epi_med 59.8 km, 6 families <= 25 km vs PFB 55.8 km / 9).
WITHHELD_B = ("VGZ", "SNB", "B926", "YOUB")
WITHHELD = WITHHELD_A if BED_MODE == "a" else WITHHELD_B
STUDY_CACHE = f"{STUDY}/cache20hz"
STEP1B_STATIONS_CSV = f"{RES}/step1b_stations.csv"
# ---- bed "c": the 2005-09 POLARIS bed.  Everything about it (year, window, quiet days, stations,
# sites, withheld) is decided ONCE by w_setup1c.py and frozen into results/step1c_bed.json before
# any byte is scored; this module only reads it.  NO engine constant changes.
STEP1C_BED = f"{RES}/step1c_bed.json"
STEP1C_CFG = {}
STEP1C_SITE_KEYS = []
STEP1C_QUIET_DAYS = []
STEP1C_ENERGY_HOURS = f"{RES}/step1c_energy_hours.csv"
# SELECTED BY MEMBER STATION, NEVER BY SITE ID.  The site id is a data-driven tie-break (the
# member with the most cached files), so it CHANGES when a cache root is added: the
# {B005,B006,B007,LRIV} cluster is "LRIVq" on the Step-1a roots and "B006q" once the study cache
# joins them.  Matching on the id silently DROPPED that site from Step 1b -- 5 sites became 4,
# MIN_SITES 4 could then never be cleared when BPCB was absent, and the first Step-1b run returned
# 0 pairs over 460 tremor hours.  Membership is stable; the id is not.
# After the 1-km merge and the withheld removal.  RULE (written before any Step-1b scoring):
# the maximal subset of DELIVERABLE 3-C sites, added in order of increasing median epicentral
# distance, whose MEDIAN S-moveout spread over the 130 families stays <= AC_WIN_S (10 s, the
# frozen common window).  LZB 28.0 km + PGCq 37.8 + BPCB 53.0 + LRIVq 54.8 + B003 60.3 = 9.7 s;
# adding B004 would make it 10.6 s.  Step 1a's six permanent sites were 13.6 s.
STEP1B_SITE_KEYS = ["B003", "BPCB", "LRIV", "LZB", "PGC"]
SITE_MERGE_KM = 1.0                                      # P2: data-driven single-linkage merge
MIN_SITES = 4                                            # P4
GRID_LAT0, GRID_LON0, GRID_HALFWIDTH = 48.50, -123.75, 62.0
ISOLATED_HIDDEN = [8, 10, 12, 15, 29, 62, 87, 113, 115, 119]   # holdouts.json stratum isolated_ge5km

if BED_MODE == "c" and os.path.exists(STEP1C_BED):
    STEP1C_CFG = json.load(open(STEP1C_BED))
    STEP1_YEAR = int(STEP1C_CFG["year"])
    STEP1_DAY0, STEP1_DAY1 = STEP1C_CFG["day0"], STEP1C_CFG["day1"]
    WITHHELD = tuple(STEP1C_CFG["withheld"])
    STEP1C_SITE_KEYS = list(STEP1C_CFG.get("site_keys") or [])
    STEP1C_QUIET_DAYS = list(STEP1C_CFG["quiet_days"])

# ------------------------------------------------------------------ autocorrelation stage -------
AC_WIN_S = 10.0             # window length, cut at the SAME ABSOLUTE TIME at every station
AC_HOP_S = 0.5
AC_BLOCK_S = 3600.0         # 1-hour blocks
AC_MIN_LAG_S = 10.0         # |lag| < 10 s excluded
AC_DECLUSTER_S = 10.0       # peaks declustered at 10 s (NOT 30 s: LFEs in tremor are seconds apart)
AC_SITE_Z_COINC = 3.0       # >= 2 sites at this per-site z
AC_SITE_Z_CLIP = 8.0
AC_MIN_SITE_Z_SITES = 2
AC_MAX_PEAKS_PER_WINDOW = 40
AC_T_SCAN = 4.0             # peaks are RECORDED above this; T is applied afterwards
AC_Z_SUB = 4                # robust median/MAD of the per-row cc and of net_mad estimated on
                            # every AC_Z_SUB-th VALID lag (18 k samples/row); set 1 to use all
AC_HOUR_MIN_EPI = 5         # candidate hour: >= 5 PNSN epicentres ...
AC_HOUR_EPI_KM = 50.0       # ... within 50 km of any ACTIVE CELL of the region
ACTIVE_CELL_DEG = 0.05      # active-cell grid step
ACTIVE_CELL_MIN_EPI = 30    # >= 30 PNSN epicentres within ACTIVE_CELL_KM
ACTIVE_CELL_KM = 8.0
ACTIVE_CELL_DEP = (25.0, 50.0)
AC_HOUR_BUDGET = 0          # 0 = every candidate hour; >0 = top-N by epicentre count (per freeze)

# ------------------------------------------------------------------ threshold rule (P1) ---------
T_ENG_FLOOR = 5.0           # ENGINEERING floor; the FA rule must be able to bind above it
T_GRID = 0.25
T_MAX = 60.0
QUIET_FA_PAIRS_PER_HOUR = 1.0     # <= 1 false pair per quiet hour ...
QUIET_MAX_COMPONENT = 5           # ... AND zero quiet-hour components of size >= 5
QUIET_TREMOR_KM = 100.0           # no PNSN tremor within 100 km ...
QUIET_MIN_DAYS_FROM_TREMOR = 3    # ... and >= 3 days from any tremor day
QUIET_MIN_HOURS = 40              # >= 40 quiet hours (Merlin: >= 20 quiet days x 2 hours)
QUIET_HOURS_PER_DAY = 2           # raised only if the bed has < 20 quiet days (reported)
FA_PER_TPL_DAY = 0.1              # growth-stage FA rule (c_tiers_by_epoch)

# ------------------------------------------------------------------ clustering ------------------
NODE_MERGE_S = 2.0
SEED_MIN_MEMBERS = 5              # >= 5 members within the hour ...
SEED_MIN_MEMBERS_PM1H = 8         # ... or >= 8 across the hour +- 1 h
SEED_BIN_S = 600.0                # burst guard: members must span >= 2 distinct 10-minute bins
SEED_MIN_BINS = 2
SEED_BUDGET = int(os.environ.get("W_SEED_BUDGET", "600") or 600)   # clusters carried into growth

# ------------------------------------------------------------------ Step 1c' REPAIR FLAGS ------
# Merlin's Step-1c' prescription (results/merlin_review_step1_2026-09-20.md, defects F1/F2/F3/F6)
# is installed behind FLAGS whose defaults are the FROZEN 1c behaviour, so results/freeze_step1c.json
# still reproduces the 1c run byte-for-byte.  The 1c' values are set by scripts/z_step1cp_queue.sh
# and recorded in results/freeze_step1c_prime.json.
#
#   W_HEALTH_VETO       F1  path to the station-hour health veto table ("" = no veto)
#   W_WIN_TRANSIENT     --  network-median window-energy factor above which a 10-s window is
#                           excluded as row i AND as lag j (0 = off)
#   W_DROP_ZERO_MEDIAN  F6  a station whose hour energy median is 0 is DROPPED from the hour, and
#                           no NaN can reach the site sum
#   W_SEED_RULE         F3  "frozen" = >= 5 members in an hour + >= 2 ten-minute bins;
#                           "pairs2" = every component with >= 2 nodes, ranked by max pair net_mad
#   W_SEED_BUDGET       F3  clusters carried into growth (1c 600 -> 1c' 2000)
#   W_HANDOFF_MIN_SEG   F3  minimum stacked segments per station at handoff (1c 3 -> 1c' 2)
#   W_NEAR_SEED_MODE    F3  "median" = frozen (+-10 min of the cluster median);
#                           "either"  = +-10 min of EITHER end of the seed window
#   W_QUIET_HOURS_CSV   F2  explicit quiet-hour list (hour-matched + earthquake-screened)
HEALTH_VETO_CSV = os.environ.get("W_HEALTH_VETO", "")
WIN_TRANSIENT_FACTOR = float(os.environ.get("W_WIN_TRANSIENT", "0") or 0.0)
DROP_ZERO_MEDIAN = os.environ.get("W_DROP_ZERO_MEDIAN", "0") == "1"
SEED_RULE = os.environ.get("W_SEED_RULE", "frozen")
HANDOFF_MIN_SEGMENTS = int(os.environ.get("W_HANDOFF_MIN_SEG", "3") or 3)
NEAR_SEED_MODE = os.environ.get("W_NEAR_SEED_MODE", "median")
QUIET_HOURS_CSV = os.environ.get("W_QUIET_HOURS_CSV", "")


def repair_flags():
    """Every 1c' flag and its ACTIVE value -- written verbatim into the freeze."""
    return dict(W_HEALTH_VETO=HEALTH_VETO_CSV, W_WIN_TRANSIENT=WIN_TRANSIENT_FACTOR,
                W_DROP_ZERO_MEDIAN=DROP_ZERO_MEDIAN, W_SEED_RULE=SEED_RULE,
                W_SEED_BUDGET=SEED_BUDGET, W_HANDOFF_MIN_SEG=HANDOFF_MIN_SEGMENTS,
                W_NEAR_SEED_MODE=NEAR_SEED_MODE, W_QUIET_HOURS_CSV=QUIET_HOURS_CSV)


_VETO = [None]


def health_veto():
    """{(sta, date, hour)} from the station-hour health table (F1), or an empty set."""
    if _VETO[0] is None:
        V = set()
        if HEALTH_VETO_CSV and os.path.exists(HEALTH_VETO_CSV):
            D = pd.read_csv(HEALTH_VETO_CSV)
            D = D[D.veto.astype(bool)] if "veto" in D.columns else D
            V = {(str(r.sta), str(r.date), int(r.hour)) for r in D.itertuples()}
            log(f"[veto] {len(V)} vetoed station-hours from {HEALTH_VETO_CSV}")
        _VETO[0] = V
    return _VETO[0]


def is_vetoed(sta, date, hour):
    V = health_veto()
    return bool(V) and (str(sta), str(date), int(hour)) in V


def vetoed_at(sta, t_abs):
    """Veto lookup from an ABSOLUTE time (handoff / stacking)."""
    V = health_veto()
    if not V:
        return False
    ts = pd.Timestamp(float(t_abs), unit="s")
    return (str(sta), ts.strftime("%Y-%m-%d"), int(ts.hour)) in V


# ------------------------------------------------------------------ handoff / growth ------------
ALIGN_ITERS, ALIGN_LAG = 2, 30    # align_stack, +- 1.5 s at 20 Hz
HANDOFF_CUT_PRE_S = 8.0           # stack cut from cluster time - 8 s ...
HANDOFF_CUT_LEN_S = 36.0          # ... 36 s long, onset searched inside it
HANDOFF_ONSET_SEARCH = (-4.0, 16.0)   # energy-onset search window (s) relative to cluster time
TPL0_LEN_S = 20.0                 # pass-0 template: 20 s around the stack's energy onset
TPL0_PRE_S = 2.0                  # cut from onset - 2 s
TPL1_LEN_S = 12.0                 # iterations: 12 s from S - 2 s
TPL1_PRE_S = 2.0
GROW_DECLUSTER_S = 30.0
INDEP_S = 60.0
SURV_MIN_INDEP = 20
SURV_MIN_DAYS = 3
SURV_MAX_NEAR_SEED = 0.20
NEAR_SEED_S = 600.0
SCRAMBLE_RATIO = 3.0              # real independent hits >= 3 x template-scramble hits
N_ITERATIONS = 2                  # 2-3 iterations; 2 run by the queue, the 3rd only if it moves
MAX_STACK = 400

# ------------------------------------------------------------------ dedupe / locate / gates -----
DEDUPE_JACCARD = 0.50             # hit-set Jaccard >= 0.5 ...
DEDUPE_CC_MIN_SITES = 3           # ... or inter-template cc above the cross-family p95 at >= 3 sites
DEDUPE_TOL_S = 2.0
DEDUPE_LAG_S = 1.0
LOC_PICK_SEARCH_S = 3.0           # S on the horizontal envelope +- 3 s of the prediction
LOC_MIN_PICK_SNR = 3.0
LOC_N_PICK_NULL = 200
LOC_N_SYNTH = 50
LOC_SYNTH_BAR_KM = 3.0
LOC_LOSO_BAR_KM = 5.0
LOC_DEPTH_BAR_KM = 10.0
LOC_NULL_RATIO = 0.5
LOC_MIN_SITES = 4                 # >= 4 picked sites here (the bed has 6); WA bar is >= 5
LOC_MAX_GAP_DEG = 180.0
GATES_OFFSET_SEARCH_S = 30.0      # the tested value; offset_fit_s is read back
GATES = (0, 1, 2)                 # Gate 3 is UNTESTABLE inside one ETS window (Merlin F6)

# ------------------------------------------------------------------ nulls -----------------------
N_SITESHIFT_REALISATIONS = 3
SITESHIFT_MIN_DAYS = 5            # each site's whole-day circular shift, members shift together
SITESHIFT_SEED = 20260920

# ------------------------------------------------------------------ scoring ---------------------
SCORE_WINDOWS = "2012-08-15,2012-10-05"
SCORE_TOL_S = 2.0
SCORE_R_KM = 5.0
SCORE_RECALL_FRAC = 0.5           # relative clauses (q_rediscover)
SCORE_PURITY_FRAC = 0.8
SCORE_MAX_EXPLAINED = 0.50
SCORE_SPLIT_FRAC = 0.20
ACCEPT_REDISCOVERED = 3           # >= 3 of the 10 isolated hidden families
ACCEPT_SPLIT_RATE = 0.10
ACCEPT_NULL_GATED = 1             # <= 1 null "gated family" per realisation
CLUSTER_TRUTH_KM = 3.0            # scoring-truth-colocated-families memory: cluster truth <= 3 km


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%Y-%m-%d %H:%M:%S')}]", *a, flush=True)


# ================================================================== geometry / sites =============
_ST = [None]


def station_table():
    """The P3D roster, EXTENDED with the Step-1b station coordinates (lat/lon only)."""
    if _ST[0] is None:
        T = S.station_table().copy()
        E = extra_stations()
        for r in E.itertuples():
            if r.sta not in T.index:
                T.loc[r.sta] = {c: np.nan for c in T.columns}
            T.loc[r.sta, "lat"] = float(r.lat)
            T.loc[r.sta, "lon"] = float(r.lon)
        _ST[0] = T
    return _ST[0]


def bed_roots():
    return [BED] if BED_MODE == "a" else [BED, STUDY_CACHE]


def bed_cfg():
    return STEP1C_CFG


def extra_stations():
    """Step-1b and Step-1c stations with their FDSN coordinates (net, sta, lat, lon, band, prov)."""
    fr = []
    for p in (STEP1B_STATIONS_CSV, f"{RES}/step1c_stations.csv"):
        if os.path.exists(p) and os.path.getsize(p) > 0:
            fr.append(pd.read_csv(p))
    if not fr:
        return pd.DataFrame(columns=["net", "sta", "lat", "lon", "band", "prov"])
    D = pd.concat(fr, ignore_index=True)
    for c in ("band", "prov"):
        if c not in D.columns:
            D[c] = ""
    return D.drop_duplicates("sta", keep="last")


def bed_stations(year=None):
    """Stations with day files on the active bed."""
    year = (STEP1_YEAR if year is None else year)
    out = set()
    for r in bed_roots():
        for p in glob.glob(f"{r}/*/{year}/*.npz"):
            out.add(p.split(os.sep)[-3].split(".")[1])
    return sorted(out)


def netsta_of(sta):
    E = extra_stations()
    hit = E[E.sta == sta]
    if len(hit):
        return f"{hit.net.iloc[0]}.{sta}"
    return S.NETSTA.get(sta, f"XX.{sta}")


def merge_sites(stas, ST=None, km=SITE_MERGE_KM):
    """P2: DATA-DRIVEN single-linkage merge of the station table at `km`.

    Returns (SITE_MEMBERS dict, SITE_OF dict, ordered site list).  The site id is the member with
    the most day files, falling back to the alphabetically first."""
    ST = station_table() if ST is None else ST
    stas = [s for s in stas if s in ST.index]
    n = len(stas)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            d = float(S.epi_km(ST.lat[stas[i]], ST.lon[stas[i]], ST.lat[stas[j]], ST.lon[stas[j]]))
            if d <= km:
                parent[find(i)] = find(j)
    grp = {}
    for i in range(n):
        grp.setdefault(find(i), []).append(stas[i])
    members, of = {}, {}
    for g in grp.values():
        g = sorted(g)
        q = sorted(g, key=lambda s: (-sum(len(glob.glob(f"{r}/*.{s}/*/*.npz"))
                                          for r in bed_roots()), s))[0]
        q = q + "q" if len(g) > 1 else q
        members[q] = g
        for s in g:
            of[s] = q
    order = sorted(members)
    return members, of, order


_SITES = [None]


def sites(year=None):
    """(SITE_MEMBERS, SITE_OF, SITE_ORDER, DET_STAS) for the Step-1 bed, withheld removed."""
    year = (STEP1_YEAR if year is None else year)
    if _SITES[0] is None:
        stas = [s for s in bed_stations(year) if s not in WITHHELD]
        mem, of, order = merge_sites(stas)
        keys = (STEP1B_SITE_KEYS if BED_MODE == "b"
                else (STEP1C_SITE_KEYS if BED_MODE == "c" else None))
        if keys:
            order = [q for q in order if any(k in mem[q] for k in keys)]
            missing = [k for k in keys if not any(k in mem[q] for q in order)]
            if missing:
                log(f"[sites] WARNING: selected member(s) {missing} are not on this bed")
            mem = {q: mem[q] for q in order}
            of = {s: q for q, mm in mem.items() for s in mm}
        det = [s for q in order for s in mem[q]]
        _SITES[0] = (mem, of, order, det)
    return _SITES[0]


# ================================================================== bed index / availability =====
_IDX = [None]


def cache_idx(years=None, stas=None):
    """{(sta, year, doy): path} over the active bed roots.

    Own glob index (NOT `s_common.cache_index`): that function filters to the P3D roster, so the
    Step-1b stations (SOKB, LZB, PA0x, BPCB, PFB, YOUB) would be silently invisible."""
    years = (STEP1_YEAR,) if years is None else years
    if _IDX[0] is None:
        _, _, _, det = sites()
        want = set(det) | set(WITHHELD)
        idx = {}
        for r in bed_roots():
            for p in glob.glob(f"{r}/*/*/*.npz"):
                parts = p.split(os.sep)
                sta = parts[-3].split(".")[1]
                yr = int(parts[-2])
                if sta not in want or yr not in set(years):
                    continue
                idx.setdefault((sta, yr, int(parts[-1][:3])), p)
        _IDX[0] = idx
    if stas is None:
        return _IDX[0]
    return {k: v for k, v in _IDX[0].items() if k[0] in set(stas)}


def load_day(sta, year, doy, idx=None):
    """(t0, (3, n) float32 ZNE) for ANY station of the active bed.

    Boreholes and the P3D roster go through `s_common.load_day_zne` (StationXML azimuth rotation);
    the Step-1b surface broadbands are already ZNE (HHZ/HHN/HHE, BHZ/BHN/BHE, dip -90 vertical) and
    are read straight from the frozen cache."""
    idx = cache_idx() if idx is None else idx
    if sta in S.AZ12:
        return _ORIG_LOAD_DAY_ZNE(sta, int(year), int(doy), idx)
    p = idx.get((sta, int(year), int(doy)))
    if p is None:
        return None
    z = np.load(p)
    x = z["x"]
    if x.shape[0] != 3:
        return None
    order = [str(o) for o in z["order"]]
    if [o[-1] for o in order] != ["Z", "N", "E"]:
        return None
    return float(z["t0"]), np.ascontiguousarray(x, dtype=np.float32)


def bed_days(year=None):
    year = (STEP1_YEAR if year is None else year)
    idx = cache_idx()
    ds = sorted({pd.Timestamp(year=k[1], month=1, day=1) + pd.Timedelta(days=k[2] - 1)
                 for k in idx if k[1] == year})
    return ds


def availability(year=None):
    """PER-DAY availability, NEVER intersected (P4)."""
    year = (STEP1_YEAR if year is None else year)
    mem, of, order, det = sites(year)
    idx = cache_idx()
    rows = []
    for d in bed_days(year):
        ss = sorted(s for s in det if (s, int(d.year), int(d.dayofyear)) in idx)
        qq = sorted({of[s] for s in ss})
        rows.append(dict(date=d.strftime("%Y-%m-%d"), doy=int(d.dayofyear), n_sta=len(ss),
                         stas=",".join(ss), sites=",".join(qq), n_sites=len(qq),
                         usable=bool(len(qq) >= MIN_SITES)))
    return pd.DataFrame(rows)


# ================================================================== tremor catalogue =============
_TR = [None]


def tremor(year=None):
    year = (STEP1_YEAR if year is None else year)
    if _TR[0] is None:
        tr = pd.read_csv(PNSN, usecols=["time", "lat", "lon", "depth"])
        tr["time"] = pd.to_datetime(tr.time, errors="coerce").astype("datetime64[ns]")
        tr = tr.dropna(subset=["time"])
        _TR[0] = tr
    tr = _TR[0]
    return tr[tr.time.dt.year == year].reset_index(drop=True) if year else tr


def active_cells():
    """Slab2 interface cells of the north box with >= ACTIVE_CELL_MIN_EPI PNSN epicentres within
    ACTIVE_CELL_KM (2010-2013), depth inside ACTIVE_CELL_DEP.  Catalogue-driven: no family
    coordinate, hidden or kept, is used anywhere in hour or quiet-day selection."""
    import g_common as G
    tr = tremor(None)
    tr = tr[(tr.time >= pd.Timestamp("2010-01-01")) & (tr.time < pd.Timestamp("2014-01-01"))]
    half = GRID_HALFWIDTH / 111.195
    lat = np.arange(GRID_LAT0 - half, GRID_LAT0 + half + 1e-9, ACTIVE_CELL_DEG)
    lon = np.arange(GRID_LON0 - half / np.cos(np.radians(GRID_LAT0)),
                    GRID_LON0 + half / np.cos(np.radians(GRID_LAT0)) + 1e-9, ACTIVE_CELL_DEG)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    LA, LO = LA.ravel(), LO.ravel()
    Z = G.slab2_depth(LA, LO)
    ok = np.isfinite(Z) & (Z >= ACTIVE_CELL_DEP[0]) & (Z <= ACTIVE_CELL_DEP[1])
    LA, LO, Z = LA[ok], LO[ok], Z[ok]
    tla, tlo = tr.lat.to_numpy(), tr.lon.to_numpy()
    cnt = np.zeros(len(LA), int)
    for i in range(0, len(LA), 200):
        sl = slice(i, min(i + 200, len(LA)))
        d = S.epi_km(LA[sl][:, None], LO[sl][:, None], tla[None, :], tlo[None, :])
        cnt[sl] = (d <= ACTIVE_CELL_KM).sum(1)
    keep = cnt >= ACTIVE_CELL_MIN_EPI
    return pd.DataFrame(dict(lat=LA[keep], lon=LO[keep], dep=Z[keep], n_epi=cnt[keep]))


def candidate_hours(year=None, cells=None):
    """Hours of the bed with >= AC_HOUR_MIN_EPI PNSN epicentres within AC_HOUR_EPI_KM of any
    active cell.  Returns DataFrame date, hour, n_epi.

    BED "c": no tremor catalogue exists for 2005, so hours come from the ARRAY-ENERGY FLAG
    (hourly network-median band energy against the quiet-day distribution, w_energy_hours.py).
    Merlin's Q3 ruling: array energy is legitimate for HOUR SELECTION ONLY -- P4-0b's negative was
    about energy as an EVENT SEED, not as an hour flag.  No Bostock event time is used."""
    year = (STEP1_YEAR if year is None else year)
    if BED_MODE == "c":
        if not os.path.exists(STEP1C_ENERGY_HOURS):
            raise SystemExit(f"bed c needs {STEP1C_ENERGY_HOURS} (run w_energy_hours.py)")
        H = pd.read_csv(STEP1C_ENERGY_HOURS)
        H = H[H.candidate.astype(bool)][["date", "hour"]].copy()
        H["n_epi"] = 0
        A = availability(year)
        ok = set(A[A.usable].date)
        return H[H.date.isin(ok)].sort_values(["date", "hour"]).reset_index(drop=True)
    cells = active_cells() if cells is None else cells
    tr = tremor(year)
    days = {d.strftime("%Y-%m-%d") for d in bed_days(year)}
    tr = tr[tr.time.dt.strftime("%Y-%m-%d").isin(days)]
    if not len(tr):
        return pd.DataFrame(columns=["date", "hour", "n_epi"])
    tla, tlo = tr.lat.to_numpy(), tr.lon.to_numpy()
    near = np.zeros(len(tr), bool)
    cla, clo = cells.lat.to_numpy(), cells.lon.to_numpy()
    for i in range(0, len(cla), 200):
        sl = slice(i, min(i + 200, len(cla)))
        d = S.epi_km(cla[sl][:, None], clo[sl][:, None], tla[None, :], tlo[None, :])
        near |= (d <= AC_HOUR_EPI_KM).any(0)
    tr = tr[near]
    key = pd.DataFrame(dict(date=tr.time.dt.strftime("%Y-%m-%d").to_numpy(),
                            hour=tr.time.dt.hour.to_numpy()))
    H = key.groupby(["date", "hour"], as_index=False).size().rename(columns={"size": "n_epi"})
    H = H[H.n_epi >= AC_HOUR_MIN_EPI].sort_values(["date", "hour"]).reset_index(drop=True)
    A = availability(year)
    ok = set(A[A.usable].date)
    return H[H.date.isin(ok)].reset_index(drop=True)


def quiet_days(year=None):
    """Bed days with NO PNSN tremor within QUIET_TREMOR_KM of the region centre and >=
    QUIET_MIN_DAYS_FROM_TREMOR days from any such tremor day.

    BED "c" (2005) has no tremor catalogue -- the PNSN product starts 2009-08 and the Slow
    Earthquake Database's WTD catalogue could not be retrieved (its only public interface is a POST
    endpoint whose parameter shape returns [] even for a catalogue inside its own date range).  The
    quiet days are therefore LITERATURE-QUIET: >= 30 days from the published 2005-09 ETS, on the
    same stations, frozen in step1c_bed.json.  No Bostock event time is used."""
    year = (STEP1_YEAR if year is None else year)
    if BED_MODE == "c":
        A = availability(year)
        qs = set(STEP1C_QUIET_DAYS)
        return pd.DataFrame([dict(date=r.date, n_sites=r.n_sites, days_from_tremor=np.nan,
                                  quiet=bool(r.date in qs and r.usable))
                             for r in A.itertuples()])
    tr = tremor(None)
    d = S.epi_km(tr.lat.to_numpy(), tr.lon.to_numpy(), GRID_LAT0, GRID_LON0)
    td = sorted({t.date() for t in tr.time[d <= QUIET_TREMOR_KM]})
    tds = np.array([pd.Timestamp(x).value for x in td], dtype="int64")
    A = availability(year)
    rows = []
    for r in A[A.usable].itertuples():
        t = pd.Timestamp(r.date)
        gap = np.abs(tds - t.value).min() / 86400e9 if len(tds) else 1e9
        rows.append(dict(date=r.date, n_sites=r.n_sites, days_from_tremor=float(gap),
                         quiet=bool(gap >= QUIET_MIN_DAYS_FROM_TREMOR)))
    return pd.DataFrame(rows)


def quiet_hours(year=None, per_day=None, min_hours=QUIET_MIN_HOURS):
    """>= min_hours quiet hours on the quiet days of the bed, spread evenly over the 24 h clock.

    DOCUMENTED DEVIATION (bed-forced).  Merlin's rule is ">= 40 quiet hours (>= 20 quiet days x
    2 hours)".  The 2012 bed carries only 52 day-files, of which few are >= 3 days from a tremor
    day, so the DAY count cannot reach 20; the binding clause (>= 40 quiet hours on the same sites)
    is met by taking more hours per quiet day.  Both counts are reported and frozen."""
    year = (STEP1_YEAR if year is None else year)
    if QUIET_HOURS_CSV and os.path.exists(QUIET_HOURS_CSV):
        # F2: the hour-matched, earthquake-screened quiet set (scripts/z_quiet.py).  The list is
        # frozen BEFORE any 1c' autocorrelation runs and is read verbatim here.
        H = pd.read_csv(QUIET_HOURS_CSV)
        H = H[["date", "hour"]].copy()
        H["hour"] = H.hour.astype(int)
        return H.reset_index(drop=True), int(H.date.nunique())
    Q = quiet_days(year)
    Q = Q[Q.quiet].sort_values("date").reset_index(drop=True)
    nd = len(Q)
    if nd == 0:
        return pd.DataFrame(columns=["date", "hour"]), nd
    per = per_day if per_day else max(QUIET_HOURS_PER_DAY, int(np.ceil(min_hours / nd)))
    per = min(per, 24)
    hrs = np.linspace(0, 23, per).round().astype(int) if per > 1 else np.array([12])
    hrs = sorted(set(int(h) for h in hrs))
    rows = [dict(date=r.date, hour=int(h)) for r in Q.itertuples() for h in hrs]
    return pd.DataFrame(rows), nd


# ================================================================== day / hour loading ===========
class DayLoader:
    """LRU over s_common.load_day_zne on the Step-1 bed, with the OPTIONAL whole-chain
    SITE-DAY-SHIFT null: each SITE's day series is circularly shifted by a site-specific whole-day
    offset (members shift together), while the returned t0 stays on the NOMINAL day, so
    site-coherent noise is preserved and true inter-site coherence destroyed."""

    def __init__(self, shift=None, maxn=24, year=None):
        from collections import OrderedDict
        year = (STEP1_YEAR if year is None else year)
        self.idx = cache_idx()
        self.d = OrderedDict()
        self.maxn = maxn
        self.year = year
        self.days = bed_days(year)
        self.doys = [int(x.dayofyear) for x in self.days]
        self.shift = dict(shift or {})
        self.last_residual = {}
        self.of = sites(year)[1]

    def _src(self, sta, year, doy):
        q = self.of.get(sta)
        k = int(self.shift.get(q, 0))
        if k == 0 or doy not in self.doys:
            return year, doy
        i = self.doys.index(doy)
        return year, self.doys[(i + k) % len(self.doys)]

    def get(self, sta, year, doy):
        """(t0_nominal, (3, N) float32 ZNE) or None."""
        key = (sta, year, doy)
        if key in self.d:
            self.d.move_to_end(key)
            return self.d[key]
        sy, sd = self._src(sta, year, doy)
        v = load_day(sta, sy, sd, self.idx)
        if v is not None and (sy, sd) != (year, doy):
            t0_src = pd.Timestamp(year=sy, month=1, day=1) + pd.Timedelta(days=sd - 1)
            t0_nom = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)
            v = (float(v[0]) + (t0_nom.timestamp() - t0_src.timestamp()), v[1])
        self.d[key] = v
        if len(self.d) > self.maxn:
            self.d.popitem(last=False)
        return v

    def hour(self, sta, date, hour, pad_s=0.0):
        """(t_start_abs, (3, n)) for one hour of one station, cut at the SAME ABSOLUTE TIME at
        every station, or None.

        BUGFIX 2026-09-20.  This used to index from each FILE's own start (`t0 + a/FS`), and
        `w_autocorr.autocorr_hour` then required every station's returned start to match the first
        station's to 1e-6 s.  The frozen cache does NOT store a common t0: `o_cache_month` takes
        t0 from the first trace's `starttime`, which on this bed is midnight + 0.0000 s for the
        CN/C8 broadbands, +0.0084 s for the PB boreholes and +0.0150 s for UW.LRIV -- all far
        below one 20-Hz sample (0.05 s), but not equal.  The equality test therefore silently
        dropped every station whose offset differed from the first one's: Step 1a ran on 5 sites /
        9 stations instead of 6 / 12 (NLLB, PGC and LRIV dropped from every hour), and Step 1b on
        3 sites, below MIN_SITES, which is why it returned 0 pairs over 460 tremor hours.

        The cut is now indexed from the ABSOLUTE target time, per station, rounded to the nearest
        sample; a station whose sample grid is off by more than a quarter sample is rejected rather
        than silently misaligned."""
        t = pd.Timestamp(date)
        v = self.get(sta, int(t.year), int(t.dayofyear))
        if v is None:
            return None
        t0, x = v
        t_target = float(t.timestamp()) + hour * 3600.0 - pad_s
        i0f = (t_target - t0) * FS
        i0 = int(round(i0f))
        # Any sample-aligned cut of a continuous target time leaves a residual of at most half a
        # sample (25 ms at 20 Hz); that is the floor of the method, not a defect.  The residual is
        # CONSTANT per station over the hour, so it cancels exactly inside each station's own
        # autocorrelation and can only smear the cross-station sum by < 1 sample.  Recorded, not
        # rejected -- rejecting at a quarter sample threw away UW.LRIV (0.30 samples).
        if abs(i0f - i0) > 0.5 + 1e-9:
            return None
        self.last_residual[sta] = float(i0f - i0) / FS
        n = int(round((AC_BLOCK_S + 2 * pad_s) * FS))
        if i0 < 0:
            return None
        if i0 + n > x.shape[1]:
            n = x.shape[1] - i0
        if n < int(0.5 * AC_BLOCK_S * FS):
            return None
        return t0 + i0 / FS, np.ascontiguousarray(x[:, i0:i0 + n])


def site_shift(realisation, order=None, ndays=None, seed=SITESHIFT_SEED):
    """{site: whole-day circular offset} for one null realisation (members shift together)."""
    order = sites()[2] if order is None else order
    ndays = len(bed_days()) if ndays is None else ndays
    rng = np.random.default_rng(seed + 1000 * int(realisation))
    lo, hi = SITESHIFT_MIN_DAYS, ndays - SITESHIFT_MIN_DAYS
    return {q: int(rng.integers(lo, hi + 1)) for q in order}


# ================================================================== v_ patching ==================
def patch_v(shift=None, year=None):
    """Install the four engine prerequisites into the imported v_ modules, in process.

    NOTHING on disk is edited.  Returns (v_common, v_engine) after patching."""
    import v_common as V
    import v_engine as E
    import torch

    year = (STEP1_YEAR if year is None else year)
    mem, of, order, det = sites(year)
    V.SITE_MEMBERS = dict(mem)
    V.SITE_ORDER = list(order)
    V.SITE_OF = dict(of)
    V.DET_STAS = list(det)
    V.WITHHELD_STAS = list(WITHHELD)
    V.MIN_SITES = MIN_SITES                          # P4
    V.T_FLOOR = T_ENG_FLOOR                          # P1 (engineering floor only)
    V.T_GRID, V.T_MAX = T_GRID, T_MAX
    V.FA_PER_TPL_DAY = FA_PER_TPL_DAY
    V.SITE_Z_CLIP = AC_SITE_Z_CLIP
    V.SITE_Z_COINC = AC_SITE_Z_COINC
    V.DECLUSTER_S = GROW_DECLUSTER_S
    V.INDEP_S = INDEP_S
    V.TPL0_LEN_S, V.TPL0_PRE_S = TPL0_LEN_S, TPL0_PRE_S
    V.TPL1_LEN_S, V.TPL1_PRE_S = TPL1_LEN_S, TPL1_PRE_S
    V.ALIGN_ITERS, V.ALIGN_LAG = ALIGN_ITERS, ALIGN_LAG
    V.SURV_MIN_INDEP, V.SURV_MIN_DAYS = SURV_MIN_INDEP, SURV_MIN_DAYS
    V.SURV_MIN_EPISODES = 0
    V.SURV_MAX_NEAR_SEED, V.NEAR_SEED_S = SURV_MAX_NEAR_SEED, NEAR_SEED_S
    V.DEDUPE_OVERLAP, V.DEDUPE_TOL_S = DEDUPE_JACCARD, DEDUPE_TOL_S
    V.cache_idx = lambda stas=None: cache_idx(stas=stas)

    loader = DayLoader(shift=shift, maxn=200, year=year)
    # Every reader inside v_bank / v_iterate goes through the bed loader: it dispatches the
    # Step-1b surface broadbands (not in the P3D roster) AND carries the site-day shift when the
    # whole-chain null is running.  The gates' own stacking is NOT shifted -- it is scored on the
    # nominal axis against its own DayStore.
    S.load_day_zne = lambda sta, y, doy, idx=None: loader.get(sta, int(y), int(doy))

    class WDayData(E.DayData):
        """P3: median / lower-half MAD of the per-station cc over VALID SAMPLES ONLY.

        The frozen station_z takes the median and the 25th percentile over all NV samples, so a
        gappy day (whose cc is identically ~0 inside the gap) biases both and inflates z."""

        def __init__(self, date, stas, n, dev="cuda", idx=None):
            self.dev = dev
            self.n = int(n)
            self.stas, self.t0, self.DH, self.E, self.valid, self.NV = [], {}, {}, {}, {}, {}
            self.vidx, self.NFs = {}, {}
            y, doy = int(date.year), int(date.dayofyear)
            for s in stas:
                v = loader.get(s, y, doy)
                if v is None:
                    continue
                t0, x = v
                N = x.shape[1]
                if N < 20 * self.n:
                    continue
                NV = N - self.n
                NF = int(2 ** np.ceil(np.log2(N + self.n)))
                X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
                self.DH[s] = torch.fft.rfft(X, n=NF)
                xd = X.double()
                cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                             xd * xd], 1), 1)
                ew = (cs[:, self.n:N + 1] - cs[:, :NV + 1])[:, :NV]
                e3 = ew.sum(0)
                self.valid[s] = ((e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0))
                vi = torch.nonzero(self.valid[s], as_tuple=False).squeeze(1)
                self.vidx[s] = vi if vi.numel() >= 1000 else torch.arange(NV, device=dev)
                floor = 1e-10 * ew.median()
                self.E[s] = ew.clamp_min(floor).sqrt().float()
                self.NV[s] = NV
                self.t0[s] = t0
                self.NFs[s] = NF
                self.NF = NF
                self.stas.append(s)
                del X, xd, cs, ew, e3
            self.sites = sorted({V.SITE_OF[s] for s in self.stas if s in V.SITE_OF})

        def station_z(self, s, W):
            # BUGFIX vs v_engine: NF is PER STATION.  v_engine keeps a single self.NF (the last
            # station's), which is harmless only while every station-day has exactly the same
            # sample count; a short/partial day makes the rFFT lengths disagree.
            NF, NV = self.NFs[s], self.NV[s]
            WH = torch.fft.rfft(W, n=NF)
            num = torch.fft.irfft(WH.conj() * self.DH[s][None], n=NF)[:, :, :NV]
            cc = (num / self.E[s][None]).mean(1)
            del WH, num
            vi = self.vidx[s]
            ccv = cc.index_select(1, vi)                       # P3: VALID samples only
            m = ccv.median(dim=-1, keepdim=True).values
            q = torch.kthvalue(ccv, max(1, ccv.shape[-1] // 4), dim=-1, keepdim=True).values
            sc = (1.4826 * (m - q)).clamp_min(1e-9)
            del ccv
            return (cc - m) / sc

        def free(self):
            self.DH.clear(); self.E.clear(); self.valid.clear(); self.vidx.clear()
            self.NFs.clear()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    E.DayData = WDayData
    E._W_LOADER = loader
    return V, E


# ================================================================== misc =========================
def independent(times, sep_s=INDEP_S):
    t = np.asarray(times, float)
    o = np.argsort(t)
    keep = np.zeros(len(t), bool)
    last = -1e18
    for j in o:
        if t[j] - last >= sep_s:
            keep[j] = True
            last = t[j]
    return keep


def read_csv_safe(path, cols=None):
    """pd.read_csv that returns an EMPTY DataFrame (with `cols`) for a missing / header-less /
    zero-byte file.  Every stage of this engine must survive "nothing survived" without an
    exception -- a null result is a RESULT, not a crash."""
    cols = list(cols or [])
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=cols)
    try:
        d = pd.read_csv(path)
    except Exception:                                                       # noqa: BLE001
        return pd.DataFrame(columns=cols)
    for c in cols:
        if c not in d.columns:
            d[c] = pd.Series(dtype=object)
    return d


SURV_COLS = ["tpl", "tid", "n_hits", "n_indep", "n_days", "frac_near_seed", "n_scramble",
             "scramble_ratio", "survives"]
THR_COLS = ["tpl", "tid", "T_seed", "n_quiet_hits", "quiet_days", "fa_at_T", "floor_binds",
            "saturated"]
HIT_COLS = ["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]
DEDUPE_COLS = ["tid", "group", "kept", "n_indep"]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_of(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def w_scripts():
    return sorted(glob.glob(f"{HERE}/w_*.py")) + sorted(glob.glob(f"{HERE}/w_*.sh"))
