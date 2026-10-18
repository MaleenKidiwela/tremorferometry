#!/usr/bin/env python
"""S15 TASK B — a third tremor measure: the project's OWN picker-ranked envelope-peak candidates.

PLAN.md Round 5, S15(B).  S14 showed that the project's cc >= 0.80 densify detections
(`data/mf_*.csv`) are false-alarm-dominated stacking material (ETS/non-ETS ratio x1.24, rho 0.09
with PNSN) and cannot deliver the exponent's VALUE.  This step tries the OTHER product of the
pipeline: `data/<sta>_cand_filtered.parquet` — the candidates that the discovery stage found and
that the LFE picker then ranked.

WHAT A CANDIDATE IS (read `scripts/discover_nllb_pnsn_driven.py` and
`lfe_features/score_candidates.py`; this fixes the measure's independence):
  * PNSN sets WHEN to look.  For every PNSN tremor window inside the station's bbox
    (+-0.9 deg lat, +-1.35 deg lon), a `window_s` segment of that station's 2-8 Hz record is cut.
  * The STATION sets the COUNT.  Inside the window the Hilbert envelope is computed, the
    background is its median, and every separated peak above `snr_threshold` x background is
    emitted as one candidate (`peak_sep_s` minimum separation).  So a PNSN window contributes
    0, 1, 2, ... candidates depending on what that station's envelope actually did.
  * `time`  = the envelope-peak time AT THE STATION (station-side information)
  * `lat/lon` = the PNSN window's own epicentre (the SOURCE CELL; not the station, not relocated)
  * `snr`   = envelope peak / median envelope in the window
  * `p_lfe` = the tremor-picker's P(LFE) for that peak's waveform features
  * the `_cand_filtered` file keeps the top-ranked candidates (`scripts/adaptive_cand_threshold.py`
    takes the P threshold that yields ~30,000 per station, floored at 0.05, so selection is
    RANK-based; absolute P does not transfer across instrument classes)
=> **THE MEASURE IS SEMI-INDEPENDENT, and the write-up must say so.**  Candidates exist ONLY inside
PNSN tremor windows, so C_cand can never see tremor PNSN missed; what is independent of PNSN is the
MULTIPLICITY (how many envelope peaks each window produced at each station) and the picker rank.
The script therefore measures the decomposition explicitly:
    b(C_cand, R7)  =  b(C_pnsn, R7)  +  b(C_cand / C_pnsn, R7)
and reports the multiplicity exponent as the genuinely new number.

======================================================================================
METHOD (events, window, footprint, fitting, nulls IDENTICAL to S9/S11; only the tremor
measure changes.  Loaders copied from scripts/s11_scaling_robustness.py.)
======================================================================================
Base product: daily candidate counts per STATION and per 0.25-deg SOURCE-CELL latitude band
(band from the candidate lat, i.e. the PNSN window latitude), over 2010-01-01..2026-05-27.

Per event (the 27 confirmed north/central patches of `results/s1_events.csv`):
  window   W  = [start - 10, end + 10]
  footprint F = mesh triangles labelled with the patch id on any day of [start, end]
                (`slip_patches_5mm.npz`)
  candidates kept: source cell within 60 km of any F centroid AND station within 150 km of the
                area-weighted F centroid
  C_cand   = daily candidate count summed over the kept stations
  W_cand   = daily sum of p_lfe over the same candidates (picker-weighted count)
  C_pnsn   = daily PNSN detection count within 60 km of F  (S9/S11's measure, the reference)
  M_cand   = C_cand / C_pnsn                                (station-side multiplicity)
  C_q75    = daily count of the candidates in the event's TOP-QUARTILE p_lfe (picker rank control)
  C_s10    = daily count using only the 10 nearest contributing stations (fixed-network control)
  R7       = area-weighted sum over the 60-km triangle neighbourhood of the clipped 7-d centred
             mean potency rate, built from A_new = the amplitudes rebuilt from the 2026-08-18
             solution in scripts/s15_amplitude_check.py (TASK A: identical to the shipped A_raw
             to float32 round-off, so this is also S9/S11's R7)

Exponent b = Theil-Sen slope of log10(measure + eps) on log10(R7) over event-days with
measure > 0 and R7 > 0, >= 8 days; eps = 1 for the counts, the global 1st percentile of positive
daily values for W_cand, 0 for the ratio M_cand.  Median + IQR over events, 90 % bootstrap CI of
the median (2000 resamples of events), pooled LSDV fit with event dummies (OLS and Theil-Sen).

Controls:
  * SHIFTED-CANDIDATE NULL, 200 draws, two flavours:
      pooled      — one rigid >= 365-d wrapped shift of all near-candidate dates over
                    2010-01-01..2026-05-27, slip series unchanged (exactly S9/S11's null);
      per-station — each station's dates are shifted rigidly and wrapped WITHIN ITS OWN record
                    span (>= 3 yr, window inside the span), so the null cannot be satisfied by
                    swapping one station's active years for another's.  This is the stricter one,
                    because the candidate network grows x8 over the record.
  * lower vs upper half of log10 R7 (a saturating relation must flatten at high rate);
  * the 300-s WINDOW CEILING: a PNSN window is 300 s long and peaks are separated by >= 6 s, so a
    (station, window) pair can hold at most 50 candidates.  Each candidate is matched back to its
    own PNSN window (merge_asof on the exact source cell, 300-s tolerance) and the multiplicity
    distribution and ceiling occupancy are reported;
  * network growth: b by epoch (Kruskal-Wallis), Spearman of b against the event's contributing
    station count and date, and a FIXED-station-count variant using only the 10 nearest
    contributing stations of every event (the network grows from 18 to 70 usable stations);
  * ETS/quiet ratio of the daily candidate count per station (mean count inside the 27 event
    windows / mean outside, over the station's own record), paired against the SAME ratio for the
    PNSN detections inside that station's discovery bbox — the S14 diagnostic that condemned the
    mf_*.csv detections (x1.24 there);
  * paired per-event difference b(C_cand) - b(C_pnsn), Wilcoxon.

======================================================================================
INPUTS (read-only)
======================================================================================
  data/<sta>_cand_filtered.parquet      209 stations, 6.71 M candidates (time, lat, lon, snr,
                                        station, p_lfe)
  data/broadband_fleet_order.csv | data/station_slab2_depth.csv |
  data/candidate_stations_post2020.csv  station coordinates (merged; all 209 covered)
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz      labels + mesh geometry
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  A_raw fallback
  /home/jovyan/cascadia-sse-strain                                     mesh + sse_strain loaders
  catalogs/pnsn_tremor_cascadia_full.csv                               PNSN tremor
  tremor_slip_geometry/results/s1_events.csv                           the 27 confirmed events

OUTPUTS (all inside tremor_slip_geometry/)
  results/s15_cand_exponents.csv   per event x measure: b, n days, null mean/sd/p (both nulls),
                                   lower/upper-half b; plus summary rows per measure with the
                                   median, IQR, bootstrap CI, pooled fits, sign tests and the
                                   paired comparison against PNSN
  results/s15_cand_coverage.csv    per (event, station) rows: distance to the footprint centroid,
                                   candidates in the window, p_lfe stats, record span; plus
                                   per-station rows (whole-record counts, latitude-band spread,
                                   ETS/quiet ratio) and per-event summary rows
  figures/s15_candidate_scaling.png  A exponent distributions PNSN vs candidates vs nulls,
                                   B-D example event tracks, E paired per-event comparison,
                                   F coverage (stations and candidates per event over time),
                                   G lower vs upper half, H per-station ETS/quiet ratio,
                                   I the multiplicity decomposition

Run:  OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s15_candidate_scaling.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import glob
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree
from scipy.stats import binomtest, kruskal, spearmanr, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
HERE = ROOT / "tremor_slip_geometry"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)
SCRATCH = Path("/tmp/claude-1000/-home-jovyan-tremorferometry/"
               "03499166-a677-42d8-9397-e2cf5e28621f/scratchpad")

RNG = np.random.default_rng(20260910)

PRE_A, POST_A = 10, 10        # window [start-10, end+10]                  (S9/S11)
NB_RADIUS_KM = 60.0           # source-cell / tremor / triangle radius     (S9/S11)
STA_RADIUS_KM = 150.0         # station within 150 km of the F centroid    (task B)
BAND_DEG = 0.25               # source-cell latitude band                  (task B)
MIN_FIT_DAYS = 8
MIN_HALF_DAYS = 6
N_SHIFT = 200
N_BOOT = 2000
MIN_STA_SPAN_D = 3 * 365      # per-station null needs >= 3 yr of record
N_FIX_STA = 10                # fixed-network control: the 10 nearest contributing stations
PNSN_WINDOW_S = 300.0         # discover_nllb_pnsn_driven.py --window-s
PEAK_SEP_S = 6.0              # discover_nllb_pnsn_driven.py --peak-sep-s
CEIL_PER_WINDOW = int(PNSN_WINDOW_S / PEAK_SEP_S)     # 50 candidates per (station, window)
BBOX_DLAT, BBOX_DLON = 0.9, 1.35                      # fleet_station.sh discovery bbox
COS47 = np.cos(np.deg2rad(47.0))

SHIFT_LO = np.datetime64("2010-01-01")
SHIFT_HI = np.datetime64("2026-05-27")
SPAN = int((SHIFT_HI - SHIFT_LO).astype(int)) + 1

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon) + 124.0) * 111.0 * COS47, (np.asarray(lat) - 47.0) * 111.0


def ts_slope_fast(x, y, nmin=MIN_FIT_DAYS):
    """Theil-Sen slope — copied from s11_scaling_robustness.py."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < nmin:
        return np.nan
    x, y = x[m], y[m]
    dx = x[:, None] - x[None, :]
    if not np.any(dx != 0):
        return np.nan
    dy = y[:, None] - y[None, :]
    k = dx != 0
    return float(np.median(dy[k] / dx[k]))


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan, 0
    return (float(np.median(v)), float(np.percentile(v, 25)),
            float(np.percentile(v, 75)), int(v.size))


def boot_median_ci(v, nboot=N_BOOT, lo=5, hi=95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return np.nan, np.nan
    bs = np.median(v[RNG.integers(0, v.size, (nboot, v.size))], axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


def window_subset(sorted_off, s, w0_off, nW, span, lo=0):
    """Indices (into sorted_off) whose SHIFTED day lands in the window, and the day-in-window.

    Shift model: o -> lo + ((o - lo + s) mod span).  With lo = 0 and span = the global span this
    is S9/S11's helper verbatim; with lo/span = a station's own record span it wraps inside that
    station's record.
    """
    a = (w0_off - lo - s) % span
    b = (a + nW - 1) % span
    a += lo
    b += lo
    if a <= b:
        i0 = np.searchsorted(sorted_off, a, "left")
        i1 = np.searchsorted(sorted_off, b, "right")
        idx = np.arange(i0, i1)
    else:
        idx = np.concatenate([
            np.arange(np.searchsorted(sorted_off, a, "left"), sorted_off.size),
            np.arange(0, np.searchsorted(sorted_off, b, "right")),
        ])
    d = (sorted_off[idx] - lo + s) % span + lo - w0_off
    return idx, d


print("=" * 112)
print("S15 TASK B — picker-ranked candidate counts as a third tremor measure vs the slip rate")
print("=" * 112)

# ==================================================================================
# 1. Slip model  (A_new from TASK A; loaders copied from s11_scaling_robustness.py)
# ==================================================================================
log("loading the slip model ...")
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402

cache = SCRATCH / "s15_Anew.npz"
if cache.exists():
    zc = np.load(cache, allow_pickle=True)
    A_new = zc["A_new"]
    RATE7 = zc["R7"].astype(float)
    log(f"  A_new {A_new.shape} + R7 from the TASK A cache {cache}")
else:                                                            # rebuild (same code as TASK A)
    fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
    sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
    sol.invert(verbose=False)
    _, Ld = sol.component_slip_patterns()
    A = sol.component_amplitudes()
    t = pd.to_datetime(np.asarray(sol.dates))
    Ad = pd.DataFrame(A, index=t).groupby(level=0).mean()
    dd = pd.date_range(Ad.index.min(), Ad.index.max(), freq="D")
    A_new = Ad.reindex(dd).interpolate(limit_area="inside").values
    cum = -1e3 * (Ld @ np.nan_to_num(A_new).T)
    diff = np.zeros_like(cum)
    diff[:, 1:] = np.diff(cum, axis=1)
    RATE7 = np.clip(uniform_filter1d(diff, size=7, axis=1, mode="nearest"), 0.0, None)
    log(f"  A_new {A_new.shape} rebuilt from the 2026-08-18 solution")

# ==================================================================================
# 2. Geometry / labels
# ==================================================================================
NPZ = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
LABEL = NPZ["label"]                       # (T, 742)
DAYS = np.asarray(NPZ["days"], dtype=str)
DAYS64 = np.asarray(DAYS, dtype="datetime64[D]")
PLAT, PLON = NPZ["patch_lat"], NPZ["patch_lon"]
PAREA = NPZ["patch_area_km2"]
PX, PY = to_km(PLON, PLAT)
TRI_TREE = cKDTree(np.column_stack([PX, PY]))
DAY_INDEX = {d: i for i, d in enumerate(DAYS)}
NT = len(DAYS)
IDX_OF_SHIFT_LO = DAY_INDEX[str(SHIFT_LO)]
log(f"  {NT} model days {DAYS[0]}..{DAYS[-1]}; {len(PLAT)} triangles; "
    f"shift span {SPAN} d {str(SHIFT_LO)}..{str(SHIFT_HI)}")

# ==================================================================================
# 3. PNSN tremor (the reference measure)
# ==================================================================================
log("loading PNSN tremor ...")
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"])
tr = tr.dropna(subset=["lat", "lon", "t"]).reset_index(drop=True)
TR_X, TR_Y = to_km(tr["lon"].values, tr["lat"].values)
TR_DAY = tr["t"].values.astype("datetime64[ns]").astype("datetime64[D]")
TR_OFF = (TR_DAY - SHIFT_LO).astype(int)
TR_NS = tr["t"].values.astype("datetime64[ns]").astype("int64")     # cast to ns before int64
TR_LAT = tr["lat"].values.astype(float)
TR_LON = tr["lon"].values.astype(float)
TREMOR_TREE = cKDTree(np.column_stack([TR_X, TR_Y]))
log(f"  {len(tr)} PNSN detections {tr['t'].min().date()}..{tr['t'].max().date()}")

# ==================================================================================
# 4. Candidates: load every station, build the daily (station x latitude-band) counts
# ==================================================================================
log("loading candidate parquets ...")
files = sorted(glob.glob(str(ROOT / "data/*_cand_filtered.parquet")))
frames = []
for p in files:
    d = pd.read_parquet(p, columns=["time", "lat", "lon", "snr", "station", "p_lfe"])
    frames.append(d)
cd = pd.concat(frames, ignore_index=True)
del frames
cd["sta"] = cd["station"].astype(str).str.upper()
# datetime64[us] -> [ns] before any integer arithmetic (project trap)
CD_DAY = cd["time"].values.astype("datetime64[ns]").astype("datetime64[D]")
keep = (CD_DAY >= SHIFT_LO) & (CD_DAY <= SHIFT_HI)
cd = cd.loc[keep].reset_index(drop=True)
CD_DAY = CD_DAY[keep]
CD_OFF = (CD_DAY - SHIFT_LO).astype(np.int32)
CD_NS = cd["time"].values.astype("datetime64[ns]").astype("int64")  # cast to ns before int64
CD_LAT = cd["lat"].values.astype(float)
CD_LON = cd["lon"].values.astype(float)
CD_P = cd["p_lfe"].values.astype(np.float32)
CD_SNR = cd["snr"].values.astype(np.float32)
STA_NAMES, CD_STA = np.unique(cd["sta"].values.astype(str), return_inverse=True)
CD_STA = CD_STA.astype(np.int16)
NSTA = len(STA_NAMES)
log(f"  {len(files)} files -> {len(cd)} candidates, {NSTA} stations, "
    f"{CD_DAY.min()}..{CD_DAY.max()}")

# station coordinates (three tables merged)
srows = []
for f, c in (("data/broadband_fleet_order.csv", "sta"),
             ("data/station_slab2_depth.csv", "station"),
             ("data/candidate_stations_post2020.csv", "sta")):
    dd = pd.read_csv(ROOT / f).rename(columns={c: "sta"})
    srows.append(dd[["sta", "lat", "lon"]])
sc = pd.concat(srows, ignore_index=True).dropna()
sc["sta"] = sc["sta"].astype(str).str.upper()
sc = sc.groupby("sta")[["lat", "lon"]].median()
missing = [s for s in STA_NAMES if s not in sc.index]
assert not missing, f"station coordinates missing for {missing}"
STA_LAT = sc.loc[STA_NAMES, "lat"].values
STA_LON = sc.loc[STA_NAMES, "lon"].values
STA_X, STA_Y = to_km(STA_LON, STA_LAT)
log(f"  station coordinates for all {NSTA} stations "
    f"(lat {STA_LAT.min():.2f}..{STA_LAT.max():.2f})")

# --- the base product: daily counts per station and per 0.25-deg source-cell latitude band -----
BAND_LO = np.floor(CD_LAT.min() / BAND_DEG) * BAND_DEG
CD_BAND = np.floor((CD_LAT - BAND_LO) / BAND_DEG).astype(np.int32)
NBAND = int(CD_BAND.max()) + 1
NDAY = SPAN
key = (CD_STA.astype(np.int64) * NBAND + CD_BAND) * NDAY + CD_OFF
ukey, upos = np.unique(key, return_inverse=True)
CUBE_N = np.bincount(upos, minlength=ukey.size).astype(np.int32)
CUBE_W = np.bincount(upos, weights=CD_P.astype(float), minlength=ukey.size)
CUBE_DAY = (ukey % NDAY).astype(np.int32)
CUBE_BAND = ((ukey // NDAY) % NBAND).astype(np.int32)
CUBE_STA = (ukey // (NDAY * NBAND)).astype(np.int16)
log(f"  daily (station x {BAND_DEG:g}-deg band) cube: {NSTA} stations x {NBAND} bands x {NDAY} d "
    f"-> {ukey.size} occupied cells ({100 * ukey.size / (NSTA * NBAND * NDAY):.3f} % dense), "
    f"bands {BAND_LO:.2f}..{BAND_LO + NBAND * BAND_DEG:.2f} N")

# per-station daily totals (any band) and per-band daily totals — the two marginals
STA_DAILY = np.zeros((NSTA, NDAY), np.int32)
np.add.at(STA_DAILY, (CUBE_STA.astype(np.int64), CUBE_DAY.astype(np.int64)), CUBE_N)
BAND_DAILY = np.zeros((NBAND, NDAY), np.int32)
np.add.at(BAND_DAILY, (CUBE_BAND.astype(np.int64), CUBE_DAY.astype(np.int64)), CUBE_N)
STA_FIRST = np.full(NSTA, -1, np.int64)
STA_LAST = np.full(NSTA, -1, np.int64)
for s in range(NSTA):
    nz = np.flatnonzero(STA_DAILY[s] > 0)
    if nz.size:
        STA_FIRST[s], STA_LAST[s] = nz[0], nz[-1]
log(f"  station records: median span {np.median((STA_LAST - STA_FIRST + 1)[STA_FIRST >= 0]) / 365.25:.1f} yr; "
    f"stations active per year "
    f"{ {int(str(SHIFT_LO)[:4]) + i: int(((STA_FIRST <= (np.datetime64('%d-07-01' % (2010 + i)) - SHIFT_LO).astype(int)) & (STA_LAST >= (np.datetime64('%d-07-01' % (2010 + i)) - SHIFT_LO).astype(int))).sum()) for i in (0, 5, 10, 15)} }")

# ==================================================================================
# 5. Events, per-event selections
# ==================================================================================
ev_in = pd.read_csv(OUT_RES / "s1_events.csv").sort_values("start").reset_index(drop=True)
log(f"  events {len(ev_in)}  {ev_in.seg.value_counts().to_dict()}")

# unique source cells -> one KD-tree query per event instead of 6.7 M
CD_XY = np.column_stack(to_km(CD_LON, CD_LAT))
uloc, CD_LOC = np.unique(np.column_stack([CD_LAT, CD_LON]), axis=0, return_inverse=True)
CD_LOC = CD_LOC.astype(np.int32)
ULX, ULY = to_km(uloc[:, 1], uloc[:, 0])
LOC_TREE = cKDTree(np.column_stack([ULX, ULY]))
log(f"  {len(uloc)} distinct source cells among the candidates")

# ETS-day mask (union of the 27 windows) for the per-station ETS/quiet ratio
ETS_DAY = np.zeros(NDAY, bool)
for _, r in ev_in.iterrows():
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    a = max(0, ts - PRE_A - IDX_OF_SHIFT_LO)
    b = min(NDAY - 1, te + POST_A - IDX_OF_SHIFT_LO)
    if b >= a:
        ETS_DAY[a:b + 1] = True
log(f"  ETS days (union of the 27 windows, +-10 d): {ETS_DAY.sum()} of {NDAY} "
    f"({100 * ETS_DAY.mean():.1f} %)")

# PNSN daily counts inside each station's discovery bbox — the fair reference for the ETS contrast
PN_DAILY_ALL = np.zeros(NDAY)
_pok = (TR_OFF >= 0) & (TR_OFF < NDAY)
np.add.at(PN_DAILY_ALL, TR_OFF[_pok], 1.0)


def _pnsn_bbox_ratio(s, sl):
    """ETS/quiet ratio of the PNSN count inside station s's discovery bbox, over the same days."""
    m = (np.abs(TR_LAT - STA_LAT[s]) <= BBOX_DLAT) & (np.abs(TR_LON - STA_LON[s]) <= BBOX_DLON)
    m &= _pok
    d = np.zeros(NDAY)
    if m.any():
        np.add.at(d, TR_OFF[m], 1.0)
    dd = d[sl]
    e = ETS_DAY[sl]
    a = float(dd[e].mean()) if e.any() else np.nan
    b = float(dd[~e].mean()) if (~e).any() else np.nan
    return dict(pnsn_in_bbox=int(m.sum()), pnsn_mean_daily_ets=a, pnsn_mean_daily_quiet=b,
                pnsn_ets_quiet_ratio=(a / b) if (b and np.isfinite(b) and b > 0) else np.nan)


# per-station ETS/quiet ratio over the station's own record
sta_rows = []
for s in range(NSTA):
    if STA_FIRST[s] < 0:
        continue
    sl = slice(STA_FIRST[s], STA_LAST[s] + 1)
    d = STA_DAILY[s, sl].astype(float)
    e = ETS_DAY[sl]
    mets = float(d[e].mean()) if e.any() else np.nan
    mq = float(d[~e].mean()) if (~e).any() else np.nan
    bands = CUBE_BAND[CUBE_STA == s]
    sta_rows.append(dict(kind="station", sta=STA_NAMES[s], lat=STA_LAT[s], lon=STA_LON[s],
                         n_cand=int(STA_DAILY[s].sum()),
                         first_day=str(SHIFT_LO + np.timedelta64(int(STA_FIRST[s]), "D")),
                         last_day=str(SHIFT_LO + np.timedelta64(int(STA_LAST[s]), "D")),
                         span_yr=float((STA_LAST[s] - STA_FIRST[s] + 1) / 365.25),
                         days_with_cand=int((STA_DAILY[s] > 0).sum()),
                         n_bands=int(np.unique(bands).size),
                         band_lat_min=float(BAND_LO + bands.min() * BAND_DEG) if bands.size else np.nan,
                         band_lat_max=float(BAND_LO + (bands.max() + 1) * BAND_DEG) if bands.size else np.nan,
                         mean_daily_ets=mets, mean_daily_quiet=mq,
                         ets_quiet_ratio=mets / mq if (mq and np.isfinite(mq) and mq > 0) else np.nan,
                         p_lfe_median=float(np.median(CD_P[CD_STA == s])),
                         snr_median=float(np.median(CD_SNR[CD_STA == s])),
                         **_pnsn_bbox_ratio(s, sl)))
sta_df = pd.DataFrame(sta_rows)
sta_df["ratio_cand_over_pnsn"] = sta_df.ets_quiet_ratio / sta_df.pnsn_ets_quiet_ratio
r_med = float(np.nanmedian(sta_df.ets_quiet_ratio))
r_pnsn_med = float(np.nanmedian(sta_df.pnsn_ets_quiet_ratio))
_d = (sta_df.ets_quiet_ratio - sta_df.pnsn_ets_quiet_ratio).dropna().values
log(f"  per-station ETS/quiet ratio of the candidate count: median {r_med:.2f} "
    f"[{np.nanpercentile(sta_df.ets_quiet_ratio, 10):.2f}, "
    f"{np.nanpercentile(sta_df.ets_quiet_ratio, 90):.2f}], "
    f">1 for {int((sta_df.ets_quiet_ratio > 1).sum())}/{int(sta_df.ets_quiet_ratio.notna().sum())} stations")
log(f"  paired against the PNSN count in the SAME station bbox: PNSN median {r_pnsn_med:.2f}, "
    f"candidate - PNSN = {np.median(_d):+.2f}, Wilcoxon p = "
    f"{wilcoxon(_d).pvalue:.3g}, candidate higher in {int((_d > 0).sum())}/{_d.size} stations")

# ==================================================================================
# 6. Per-event series
# ==================================================================================
MEASURES = ["C_cand", "W_cand", "C_pnsn", "M_cand", "C_q75", "C_s10"]

# exact source-cell code shared by candidates and PNSN (candidate lat/lon ARE the PNSN values)
CODE_CD = (np.round(CD_LAT * 1000).astype(np.int64) * 10_000_000
           + np.round((CD_LON + 180.0) * 1000).astype(np.int64))
CODE_TR = (np.round(TR_LAT * 1000).astype(np.int64) * 10_000_000
           + np.round((TR_LON + 180.0) * 1000).astype(np.int64))
TOL_NS = int(PNSN_WINDOW_S * 1e9)

log("building per-event series ...")
EV = []
cov_rows = []
for _, r in ev_in.iterrows():
    pid = int(r["id"])
    ts, te = DAY_INDEX[r["start"]], DAY_INDEX[r["end"]]
    w0, w1 = max(0, ts - PRE_A), min(NT - 1, te + POST_A)
    nW = w1 - w0 + 1
    w0_off = w0 - IDX_OF_SHIFT_LO

    F = np.any(LABEL[ts:te + 1] == pid, axis=0)
    fx, fy = PX[F], PY[F]
    cx = float(np.average(fx, weights=PAREA[F]))
    cy = float(np.average(fy, weights=PAREA[F]))
    clat = float(np.average(PLAT[F], weights=PAREA[F]))
    clon = float(np.average(PLON[F], weights=PAREA[F]))

    # triangle neighbourhood -> R7
    hits = TRI_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    nb = (np.unique(np.concatenate([np.asarray(h, int) for h in hits]))
          if len(hits) else np.array([], int))
    nb = np.union1d(nb, np.flatnonzero(F))
    R7 = (PAREA[nb][:, None] * RATE7[np.ix_(nb, np.arange(w0, w1 + 1))]).sum(axis=0)

    # PNSN near set (S9/S11)
    thits = TREMOR_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    near_p = (np.unique(np.concatenate([np.asarray(h, int) for h in thits]))
              if len(thits) else np.array([], int))
    near_p = near_p[np.argsort(TR_OFF[near_p], kind="stable")]
    off_p = TR_OFF[near_p]
    sub, dday = window_subset(off_p, 0, w0_off, nW, SPAN)
    C_pnsn = np.bincount(dday.astype(np.int64), minlength=nW).astype(float)
    pn_win = near_p[sub]                                    # PNSN windows inside the event window

    # candidate near set: source cell within 60 km of F AND station within 150 km of the centroid
    lhits = LOC_TREE.query_ball_point(np.column_stack([fx, fy]), r=NB_RADIUS_KM)
    locs = (np.unique(np.concatenate([np.asarray(h, int) for h in lhits]))
            if len(lhits) else np.array([], int))
    loc_ok = np.zeros(len(uloc), bool)
    loc_ok[locs] = True
    dsta = np.hypot(STA_X - cx, STA_Y - cy)
    sta_ok = dsta <= STA_RADIUS_KM
    sel = loc_ok[CD_LOC] & sta_ok[CD_STA]

    order = np.argsort(CD_OFF[sel], kind="stable")
    isel = np.flatnonzero(sel)[order]                      # original candidate rows, time-sorted
    off_c = CD_OFF[isel]
    p_c = CD_P[isel].astype(float)
    s_c = CD_STA[isel]
    sub_c, dd_c = window_subset(off_c, 0, w0_off, nW, SPAN)
    ddi = dd_c.astype(np.int64)
    C_cand = np.bincount(ddi, minlength=nW).astype(float)
    W_cand = np.bincount(ddi, weights=p_c[sub_c], minlength=nW)
    with np.errstate(divide="ignore", invalid="ignore"):
        M_cand = np.where(C_pnsn > 0, C_cand / np.maximum(C_pnsn, 1e-12), np.nan)

    # picker-rank control: the event's top-quartile p_lfe candidates
    wp = p_c[sub_c]
    thr_q75 = float(np.percentile(wp, 75)) if wp.size else np.nan
    hq = wp >= thr_q75
    C_q75 = np.bincount(ddi[hq], minlength=nW).astype(float) if hq.any() else np.zeros(nW)

    # fixed-network control: only the N_FIX_STA nearest contributing stations
    ws = s_c[sub_c]
    stas_in_window = np.unique(ws)
    near_order = stas_in_window[np.argsort(dsta[stas_in_window])]
    fix_set = set(near_order[:N_FIX_STA].tolist())
    fm = np.isin(ws, list(fix_set))
    C_s10 = np.bincount(ddi[fm], minlength=nW).astype(float) if fm.any() else np.zeros(nW)

    # --- 300-s window ceiling: match every candidate back to its own PNSN window --------------
    mult_max = mult_p99 = mult_mean = frac_ceil = np.nan
    n_matched = 0
    if sub_c.size and pn_win.size:
        cdf = pd.DataFrame(dict(t=CD_NS[isel[sub_c]], code=CODE_CD[isel[sub_c]],
                                sta=ws.astype(np.int32))).sort_values("t", kind="stable")
        pdf = pd.DataFrame(dict(t=TR_NS[pn_win], code=CODE_TR[pn_win],
                                wid=pn_win.astype(np.int64))).sort_values("t", kind="stable")
        pdf = pdf.drop_duplicates(subset=["code", "t"], keep="first")
        mm = pd.merge_asof(cdf, pdf, on="t", by="code", direction="backward",
                           tolerance=TOL_NS)
        mm = mm.dropna(subset=["wid"])
        n_matched = len(mm)
        if n_matched:
            g = mm.groupby(["sta", "wid"], sort=False).size().values
            mult_max = float(g.max())
            mult_p99 = float(np.percentile(g, 99))
            mult_mean = float(g.mean())
            frac_ceil = float((g >= CEIL_PER_WINDOW).mean())

    EV.append(dict(pid=pid, seg=str(r["seg"]), start=str(r["start"]), end=str(r["end"]),
                   Mw=float(r["Mw"]), nW=nW, w0_off=w0_off, R7=R7,
                   meas=dict(C_cand=C_cand, W_cand=W_cand, C_pnsn=C_pnsn, M_cand=M_cand,
                             C_q75=C_q75, C_s10=C_s10),
                   off_c=off_c, p_c=p_c, s_c=s_c, n_near_cand=int(off_c.size),
                   n_near_pnsn=int(off_p.size), stas=stas_in_window,
                   n_sta_ok=int(sta_ok.sum()), clat=clat, clon=clon,
                   n_tri_F=int(F.sum()), n_nb=int(nb.size),
                   thr_q75=thr_q75, n_fix_sta=len(fix_set),
                   mult_max=mult_max, mult_p99=mult_p99, mult_mean=mult_mean,
                   frac_ceil=frac_ceil, n_matched=n_matched,
                   days=DAYS64[w0:w1 + 1]))

    # coverage rows, per (event, station)
    wsnr = CD_SNR[isel[sub_c]].astype(float)
    for s in stas_in_window:
        m = ws == s
        cov_rows.append(dict(kind="event_station", id=pid, start=str(r["start"]),
                             seg=str(r["seg"]), sta=STA_NAMES[s],
                             dist_to_centroid_km=float(dsta[s]),
                             n_cand_window=int(m.sum()),
                             n_days_with_cand=int(np.unique(ddi[m]).size),
                             p_lfe_median=float(np.median(wp[m])),
                             snr_median=float(np.median(wsnr[m])),
                             in_fixed10=bool(s in fix_set),
                             sta_first=str(SHIFT_LO + np.timedelta64(int(STA_FIRST[s]), "D")),
                             sta_last=str(SHIFT_LO + np.timedelta64(int(STA_LAST[s]), "D"))))

log(f"  {len(EV)} events; candidates in window: median "
    f"{int(np.median([e['meas']['C_cand'].sum() for e in EV]))}, "
    f"stations contributing: median {int(np.median([len(e['stas']) for e in EV]))}")

# eps for W_cand: global 1st percentile of positive daily values (S11's p1 rule)
allW = np.concatenate([e["meas"]["W_cand"] for e in EV])
allW = allW[np.isfinite(allW) & (allW > 0)]
EPS = {"C_cand": 1.0, "C_pnsn": 1.0, "M_cand": 0.0, "C_q75": 1.0, "C_s10": 1.0,
       "W_cand": float(np.percentile(allW, 1))}
log("  eps: " + ", ".join(f"{k}={v:.4g}" for k, v in EPS.items()))
mm_ = np.array([e["mult_mean"] for e in EV], float)
fc_ = np.array([e["frac_ceil"] for e in EV], float)
log(f"  300-s window ceiling ({CEIL_PER_WINDOW} candidates per (station, window)): "
    f"mean multiplicity {np.nanmedian(mm_):.2f} (max over events "
    f"{np.nanmax([e['mult_max'] for e in EV]):.0f}), p99 "
    f"{np.nanmedian([e['mult_p99'] for e in EV]):.0f}, fraction of (station, window) pairs at the "
    f"ceiling {np.nanmax(fc_):.2e} (max over events); "
    f"{100 * np.median([e['n_matched'] / max(e['meas']['C_cand'].sum(), 1) for e in EV]):.1f} % of "
    f"window candidates matched back to a PNSN window")

print(f"\nBASE PRODUCT — daily candidate counts per station and per {BAND_DEG:g}-deg source-cell "
      f"latitude band")
print(f"  {NSTA} stations x {NBAND} bands ({BAND_LO:.2f}-{BAND_LO + NBAND * BAND_DEG:.2f} N) "
      f"x {NDAY} days; {ukey.size} occupied (station, band, day) cells; {len(cd)} candidates "
      f"{str(CD_DAY.min())}..{str(CD_DAY.max())}")
print(f"  candidates per year: " + ", ".join(
    f"{y}:{int(((CD_DAY >= np.datetime64('%d-01-01' % y)) & (CD_DAY < np.datetime64('%d-01-01' % (y + 1)))).sum() / 1000)}k"
    for y in range(2010, 2027)))
print(f"  per-station ETS/quiet ratio of the daily count: median {r_med:.2f}, "
      f"10-90 % [{np.nanpercentile(sta_df.ets_quiet_ratio, 10):.2f}, "
      f"{np.nanpercentile(sta_df.ets_quiet_ratio, 90):.2f}]  "
      f"(S14's mf_*.csv detections: 1.24; PNSN in the same neighbourhoods: computed below)")


# ==================================================================================
# 7. Fitting + nulls
# ==================================================================================
def fit_one(y_raw, x_log, eps, trans_ok, nmin=MIN_FIT_DAYS):
    ok = np.isfinite(y_raw) & (y_raw > 0) & np.isfinite(x_log) & trans_ok
    if ok.sum() < nmin:
        return np.nan, int(ok.sum())
    return ts_slope_fast(x_log[ok], np.log10(y_raw[ok] + eps), nmin), int(ok.sum())


log("fitting observed exponents, halves, and 200-draw nulls (pooled + per-station) ...")
rows = []
OBS = {m: {} for m in MEASURES}
NULP = {m: {} for m in MEASURES}
NULS = {m: {} for m in MEASURES}

for e in EV:
    R = e["R7"]
    xl = np.where(R > 0, np.log10(np.where(R > 0, R, 1.0)), np.nan)
    rok = R > 0
    out = dict(kind="event", id=e["pid"], start=e["start"], seg=e["seg"], Mw=e["Mw"],
               nW=e["nW"], n_sta=len(e["stas"]), n_sta_within_150km=e["n_sta_ok"],
               n_cand_window=int(e["meas"]["C_cand"].sum()),
               n_pnsn_window=int(e["meas"]["C_pnsn"].sum()),
               n_near_cand_allspan=e["n_near_cand"], n_tri_F=e["n_tri_F"],
               lat_centroid=e["clat"], lon_centroid=e["clon"],
               epoch=("2010-2014" if int(e["start"][:4]) <= 2014 else
                      ("2015-2019" if int(e["start"][:4]) <= 2019 else "2020-2026")),
               year=int(e["start"][:4]), thr_q75=e["thr_q75"], n_fix_sta=e["n_fix_sta"],
               mult_mean=e["mult_mean"], mult_p99=e["mult_p99"], mult_max=e["mult_max"],
               frac_at_ceiling=e["frac_ceil"], n_matched_to_window=e["n_matched"])
    for m in MEASURES:
        b, n = fit_one(e["meas"][m], xl, EPS[m], rok)
        OBS[m][e["pid"]] = b
        out[f"b_{m}"] = b
        out[f"ndays_{m}"] = n
        # lower / upper half of log10 R7
        ok = np.isfinite(e["meas"][m]) & (e["meas"][m] > 0) & rok
        if ok.sum() >= 2 * MIN_HALF_DAYS:
            xx = xl[ok]
            yy = np.log10(e["meas"][m][ok] + EPS[m])
            med = np.median(xx)
            lo, hi = xx <= med, xx > med
            out[f"b_{m}_lower"] = (ts_slope_fast(xx[lo], yy[lo], MIN_HALF_DAYS)
                                   if lo.sum() >= MIN_HALF_DAYS else np.nan)
            out[f"b_{m}_upper"] = (ts_slope_fast(xx[hi], yy[hi], MIN_HALF_DAYS)
                                   if hi.sum() >= MIN_HALF_DAYS else np.nan)
        else:
            out[f"b_{m}_lower"] = out[f"b_{m}_upper"] = np.nan

    # --- nulls, candidate measures only (PNSN nulls are S11's, recomputed here for the pair) ---
    nW, w0_off = e["nW"], e["w0_off"]
    # (a) pooled rigid shift over the whole span (S9/S11's null)
    bn = {m: [] for m in ("C_cand", "W_cand", "C_pnsn")}
    shifts = RNG.integers(365, SPAN - 365 + 1, size=N_SHIFT)
    off_p_e = np.sort(TR_OFF[np.array([], int)]) if False else None
    for sh in shifts:
        s_, d_ = window_subset(e["off_c"], int(sh), w0_off, nW, SPAN)
        Cn = np.bincount(d_.astype(np.int64), minlength=nW).astype(float)
        Wn = np.bincount(d_.astype(np.int64), weights=e["p_c"][s_], minlength=nW)
        b, _ = fit_one(Cn, xl, EPS["C_cand"], rok)
        bn["C_cand"].append(b)
        b, _ = fit_one(Wn, xl, EPS["W_cand"], rok)
        bn["W_cand"].append(b)
    for m in ("C_cand", "W_cand"):
        NULP[m][e["pid"]] = np.asarray(bn[m], float)

    # (b) per-station shift, wrapped inside each station's own record span
    use = []
    for s in e["stas"]:
        f0, f1 = STA_FIRST[s], STA_LAST[s]
        span_s = int(f1 - f0 + 1)
        if span_s < MIN_STA_SPAN_D:
            continue
        if not (f0 <= w0_off and w0_off + nW - 1 <= f1):
            continue
        if span_s < 2 * 365 + nW:
            continue
        m = e["s_c"] == s
        use.append((int(f0), span_s, e["off_c"][m], e["p_c"][m]))
    n_used = sum(u[2].size for u in use)
    out["null_perstation_frac_cand"] = (n_used / max(e["n_near_cand"], 1))
    out["null_perstation_n_sta"] = len(use)
    if use:
        bs_c, bs_w = [], []
        for _ in range(N_SHIFT):
            Cn = np.zeros(nW)
            Wn = np.zeros(nW)
            for f0, span_s, offs, ps in use:
                sh = int(RNG.integers(365, span_s - 365 + 1))
                s_, d_ = window_subset(offs, sh, w0_off, nW, span_s, lo=f0)
                if s_.size:
                    dd = d_.astype(np.int64)
                    inw = (dd >= 0) & (dd < nW)
                    if inw.any():
                        Cn += np.bincount(dd[inw], minlength=nW)
                        Wn += np.bincount(dd[inw], weights=ps[s_][inw], minlength=nW)
            b, _ = fit_one(Cn, xl, EPS["C_cand"], rok)
            bs_c.append(b)
            b, _ = fit_one(Wn, xl, EPS["W_cand"], rok)
            bs_w.append(b)
        NULS["C_cand"][e["pid"]] = np.asarray(bs_c, float)
        NULS["W_cand"][e["pid"]] = np.asarray(bs_w, float)

    for m in ("C_cand", "W_cand"):
        for tag, store in (("pooled", NULP), ("persta", NULS)):
            v = store[m].get(e["pid"], np.array([]))
            v = v[np.isfinite(v)]
            out[f"null_{tag}_mean_{m}"] = float(v.mean()) if v.size else np.nan
            out[f"null_{tag}_sd_{m}"] = float(v.std()) if v.size else np.nan
            out[f"null_{tag}_p_{m}"] = (float((v >= OBS[m][e["pid"]]).mean())
                                        if v.size and np.isfinite(OBS[m][e["pid"]]) else np.nan)
    rows.append(out)
    log(f"  {e['pid']:5d} {e['start']} {e['seg']:7s} nsta={len(e['stas']):3d} "
        f"cand={int(e['meas']['C_cand'].sum()):7d} pnsn={int(e['meas']['C_pnsn'].sum()):6d} "
        f"b(cand)={out['b_C_cand']:+6.3f} b(w)={out['b_W_cand']:+6.3f} "
        f"b(pnsn)={out['b_C_pnsn']:+6.3f} b(mult)={out['b_M_cand']:+6.3f}")

ev_df = pd.DataFrame(rows)


# ==================================================================================
# 8. Pooled fits, summaries, paired comparison
# ==================================================================================
def pooled_fit(measure):
    xs, ys, ids = [], [], []
    for e in EV:
        R = e["R7"]
        v = e["meas"][measure]
        ok = np.isfinite(v) & (v > 0) & (R > 0)
        if ok.sum() < 4:
            continue
        xs.append(np.log10(R[ok]))
        ys.append(np.log10(v[ok] + EPS[measure]))
        ids.append(np.full(ok.sum(), e["pid"]))
    if not xs:
        return np.nan, np.nan, np.nan, 0
    d = pd.DataFrame(dict(i=np.concatenate(ids), x=np.concatenate(xs), y=np.concatenate(ys)))
    d["xc"] = d.x - d.groupby("i").x.transform("mean")
    d["yc"] = d.y - d.groupby("i").y.transform("mean")
    ols = float(np.polyfit(d.xc.values, d.yc.values, 1)[0])
    tss = ts_slope_fast(d.xc.values, d.yc.values, 20)
    rr = float(np.corrcoef(d.xc.values, d.yc.values)[0, 1])
    return ols, tss, rr, len(d)


print("\n" + "=" * 112)
print("EXPONENTS  b = Theil-Sen slope of log10(measure + eps) on log10(R7), per event")
print("=" * 112)
print(f"{'measure':10s} {'n':>3s}  {'b median [IQR]':28s} {'90% boot CI':22s} "
      f"{'pooled OLS/TS':16s} {'b<1':>7s}  null(pooled) p / sign   null(per-station) p / sign")
sum_rows = []
for m in MEASURES:
    v = np.array([OBS[m][e["pid"]] for e in EV], float)
    md, q1, q3, n = med_iqr(v)
    lo, hi = boot_median_ci(v)
    ols, tss, rr, npool = pooled_fit(m)
    nb1 = int(np.sum(v[np.isfinite(v)] < 1))
    line = (f"{m:10s} {n:3d}  {md:+8.3f} [{q1:+7.3f},{q3:+7.3f}]     "
            f"[{lo:+7.3f},{hi:+7.3f}]   {ols:+6.3f}/{tss:+6.3f}  {nb1:3d}/{np.isfinite(v).sum():<3d} ")
    srow = dict(kind="summary", measure=m, n_events=n, b_median=md, b_q25=q1, b_q75=q3,
                b_boot_lo90=lo, b_boot_hi90=hi, n_below_1=nb1,
                pooled_ols=ols, pooled_theilsen=tss, pooled_r=rr, n_pooled_days=npool,
                eps=EPS[m])
    for tag, store in (("pooled", NULP), ("persta", NULS)):
        if m in ("C_cand", "W_cand") and store[m]:
            zs, pos, sig = [], 0, 0
            nl_all = []
            for e in EV:
                nl = store[m].get(e["pid"], np.array([]))
                nl = nl[np.isfinite(nl)]
                b = OBS[m][e["pid"]]
                if nl.size < 10 or not np.isfinite(b):
                    continue
                nl_all.append(nl)
                z = (b - nl.mean()) / (nl.std() if nl.std() > 0 else np.nan)
                zs.append(z)
                pos += int(b > nl.mean())
                sig += int((nl >= b).mean() < 0.05)
            nz = len(zs)
            sp = float(binomtest(pos, nz, 0.5, alternative="greater").pvalue) if nz else np.nan
            nl_all = np.concatenate(nl_all) if nl_all else np.array([])
            srow.update({f"null_{tag}_mean": float(nl_all.mean()) if nl_all.size else np.nan,
                         f"null_{tag}_sd": float(nl_all.std()) if nl_all.size else np.nan,
                         f"null_{tag}_z_median": float(np.nanmedian(zs)) if nz else np.nan,
                         f"null_{tag}_n_events": nz, f"null_{tag}_n_positive": pos,
                         f"null_{tag}_sign_p": sp, f"null_{tag}_n_p05": sig})
            line += f" | {pos:2d}/{nz:<2d} p={sp:.1e}"
        else:
            line += " |    -        "
    # halves
    dl = ev_df[f"b_{m}_upper"] - ev_df[f"b_{m}_lower"]
    dl = dl[np.isfinite(dl)]
    if dl.size >= 6:
        srow["upper_minus_lower_median"] = float(np.median(dl))
        srow["upper_minus_lower_p"] = float(wilcoxon(dl).pvalue)
        srow["upper_minus_lower_n"] = int(dl.size)
        srow["upper_minus_lower_npos"] = int((dl > 0).sum())
    sum_rows.append(srow)
    print(line)

# paired comparison, candidates vs PNSN
print("\nPAIRED per-event comparison against the PNSN count on the same events")
pair_rows = []
for m in ("C_cand", "W_cand", "M_cand"):
    d = np.array([OBS[m][e["pid"]] - OBS["C_pnsn"][e["pid"]] for e in EV], float)
    d = d[np.isfinite(d)]
    w = float(wilcoxon(d).pvalue) if d.size >= 6 and not np.allclose(d, 0) else np.nan
    rho, rp = spearmanr(np.array([OBS[m][e["pid"]] for e in EV]),
                        np.array([OBS["C_pnsn"][e["pid"]] for e in EV]),
                        nan_policy="omit")
    lo, hi = boot_median_ci(d)
    print(f"  b({m}) - b(C_pnsn) = {np.median(d):+.3f} [{lo:+.3f}, {hi:+.3f}], "
          f"Wilcoxon p = {w:.3g}, {int((d > 0).sum())}/{d.size} positive; "
          f"Spearman(b) rho = {rho:+.3f} (p {rp:.3g})")
    pair_rows.append(dict(kind="paired", measure=m, vs="C_pnsn", n=int(d.size),
                          diff_median=float(np.median(d)), diff_boot_lo90=lo, diff_boot_hi90=hi,
                          wilcoxon_p=w, n_positive=int((d > 0).sum()),
                          spearman_rho=float(rho), spearman_p=float(rp)))

# candidate count against the PNSN count itself (how much is inherited)
bcp = []
for e in EV:
    C, P = e["meas"]["C_cand"], e["meas"]["C_pnsn"]
    ok = (C > 0) & (P > 0)
    bcp.append(ts_slope_fast(np.log10(P[ok] + 1.0), np.log10(C[ok] + 1.0))
               if ok.sum() >= MIN_FIT_DAYS else np.nan)
bcp = np.array(bcp, float)
md, q1, q3, n = med_iqr(bcp)
lo, hi = boot_median_ci(bcp)
print(f"  b(C_cand vs C_pnsn) = {md:+.3f} [{q1:+.3f}, {q3:+.3f}], boot CI [{lo:+.3f}, {hi:+.3f}] "
      f"(n={n})  -> how nearly the candidate count is a monotone function of the PNSN count")
pair_rows.append(dict(kind="paired", measure="C_cand", vs="C_pnsn_regressor", n=n,
                      diff_median=md, diff_boot_lo90=lo, diff_boot_hi90=hi))
rho_daily = []
for e in EV:
    C, P = e["meas"]["C_cand"], e["meas"]["C_pnsn"]
    if np.isfinite(C).sum() > 6:
        rho_daily.append(spearmanr(C, P).statistic)
print(f"  Spearman(daily C_cand, daily C_pnsn) within events: median "
      f"{np.nanmedian(rho_daily):+.3f} (S14's mf_*.csv value was +0.09)")

# the PNSN ETS/quiet ratio, margin-wide and paired per station inside the station's own bbox
pn_r = PN_DAILY_ALL[ETS_DAY].mean() / PN_DAILY_ALL[~ETS_DAY].mean()
print(f"  ETS/quiet contrast of the daily count: candidates {r_med:.2f} (median over "
      f"{int(sta_df.ets_quiet_ratio.notna().sum())} stations), PNSN in the SAME station bbox "
      f"{r_pnsn_med:.2f} (paired difference {np.median(_d):+.2f}, Wilcoxon p "
      f"{wilcoxon(_d).pvalue:.2g}), margin-wide PNSN {pn_r:.2f}, S14 mf_*.csv 1.24")

# ---- network growth / saturation controls -------------------------------------------------
print("\nCONTROLS")
gr_rows = []
for m in ("C_cand", "W_cand", "C_s10", "C_q75", "C_pnsn"):
    v = ev_df[f"b_{m}"].values
    ns_ = ev_df.n_sta.values.astype(float)
    yr_ = ev_df.year.values.astype(float)
    groups = [v[(ev_df.epoch == g).values] for g in ("2010-2014", "2015-2019", "2020-2026")]
    groups = [g[np.isfinite(g)] for g in groups]
    kp = float(kruskal(*[g for g in groups if g.size]).pvalue) if sum(g.size > 0 for g in groups) >= 2 else np.nan
    rs, ps = spearmanr(v, ns_, nan_policy="omit")
    ry, py = spearmanr(v, yr_, nan_policy="omit")
    print(f"  {m:8s} b by epoch {np.median(groups[0]):+.3f} / {np.median(groups[1]):+.3f} / "
          f"{np.median(groups[2]):+.3f} (Kruskal p {kp:.2f});  Spearman(b, n_sta) "
          f"{rs:+.3f} (p {ps:.2g});  Spearman(b, year) {ry:+.3f} (p {py:.2g})")
    gr_rows.append(dict(kind="growth", measure=m, b_2010_2014=float(np.median(groups[0])),
                        b_2015_2019=float(np.median(groups[1])),
                        b_2020_2026=float(np.median(groups[2])), kruskal_p=kp,
                        spearman_b_nsta=float(rs), spearman_b_nsta_p=float(ps),
                        spearman_b_year=float(ry), spearman_b_year_p=float(py)))
for m in ("C_cand", "W_cand", "C_q75", "C_s10", "C_pnsn"):
    d = (ev_df[f"b_{m}_upper"] - ev_df[f"b_{m}_lower"]).dropna().values
    pw = float(wilcoxon(d).pvalue)
    verdict = ("no significant change of slope across the rate range"
               if pw > 0.05 else
               ("steeper at high rate" if np.median(d) > 0 else "FLATTENS at high rate"))
    print(f"  {m:8s} upper minus lower half of log10 R7: {np.median(d):+.3f} "
          f"(Wilcoxon p {pw:.2g}, {int((d > 0).sum())}/{d.size} positive) -> {verdict}")
d10 = (ev_df.b_C_s10 - ev_df.b_C_cand).dropna().values
print(f"  fixed-network control (only the {N_FIX_STA} nearest contributing stations): "
      f"b = {np.nanmedian(ev_df.b_C_s10):+.3f} vs {np.nanmedian(ev_df.b_C_cand):+.3f} "
      f"(paired {np.median(d10):+.3f}, Wilcoxon p {wilcoxon(d10).pvalue:.2g})")
dq = (ev_df.b_C_q75 - ev_df.b_C_cand).dropna().values
print(f"  picker-rank control (top-quartile p_lfe only): b = {np.nanmedian(ev_df.b_C_q75):+.3f} "
      f"(paired {np.median(dq):+.3f}, Wilcoxon p {wilcoxon(dq).pvalue:.2g})")

# ---- is the higher exponent just the COUNT MAGNITUDE (the +1 offset / thinning)? -------------
# The candidate count is ~10x the PNSN count in the same window, so log10(C+1) compresses the
# low-count PNSN days more.  Two controls: (a) drop the +1 (eps = 0; the fit already keeps only
# days with count > 0), (b) binomially THIN the candidate count down to the PNSN mean per event
# and refit with the same eps = 1.
b0 = {m: [] for m in ("C_cand", "C_pnsn")}
b_thin = []
N_THIN = 20
for e in EV:
    R = e["R7"]
    xl = np.where(R > 0, np.log10(np.where(R > 0, R, 1.0)), np.nan)
    rok = R > 0
    for m in ("C_cand", "C_pnsn"):
        b, _ = fit_one(e["meas"][m], xl, 0.0, rok)
        b0[m].append(b)
    C, P = e["meas"]["C_cand"], e["meas"]["C_pnsn"]
    p_thin = min(1.0, P.sum() / max(C.sum(), 1.0))
    bs = []
    for _ in range(N_THIN):
        Ct = RNG.binomial(C.astype(np.int64), p_thin).astype(float)
        b, _ = fit_one(Ct, xl, 1.0, rok)
        bs.append(b)
    b_thin.append(float(np.nanmedian(bs)))
b_thin = np.array(b_thin, float)
d0 = np.array(b0["C_cand"], float) - np.array(b0["C_pnsn"], float)
d0 = d0[np.isfinite(d0)]
dth = b_thin - ev_df.b_C_pnsn.values
dth = dth[np.isfinite(dth)]
lo_t, hi_t = boot_median_ci(b_thin)
print(f"  count-magnitude control (a) eps = 0 instead of +1: b(C_cand) "
      f"{np.nanmedian(b0['C_cand']):+.3f}, b(C_pnsn) {np.nanmedian(b0['C_pnsn']):+.3f}, "
      f"paired difference {np.median(d0):+.3f} (Wilcoxon p {wilcoxon(d0).pvalue:.2g})")
print(f"  count-magnitude control (b) candidate count binomially THINNED to the PNSN mean "
      f"({N_THIN} draws): b = {np.nanmedian(b_thin):+.3f} [{lo_t:+.3f}, {hi_t:+.3f}], "
      f"vs PNSN paired {np.median(dth):+.3f} (Wilcoxon p {wilcoxon(dth).pvalue:.2g}, "
      f"{int((dth > 0).sum())}/{dth.size} positive)")
gr_rows.append(dict(kind="magnitude_control",
                    b_C_cand_eps0=float(np.nanmedian(b0["C_cand"])),
                    b_C_pnsn_eps0=float(np.nanmedian(b0["C_pnsn"])),
                    diff_eps0_median=float(np.median(d0)),
                    diff_eps0_wilcoxon_p=float(wilcoxon(d0).pvalue),
                    b_C_cand_thinned=float(np.nanmedian(b_thin)),
                    b_thinned_boot_lo90=lo_t, b_thinned_boot_hi90=hi_t,
                    diff_thinned_median=float(np.median(dth)),
                    diff_thinned_wilcoxon_p=float(wilcoxon(dth).pvalue),
                    n_thin_draws=N_THIN))

print(f"  300-s window ceiling: mean multiplicity per (station, window) "
      f"{np.nanmedian(mm_):.2f}, p99 {np.nanmedian([e['mult_p99'] for e in EV]):.0f}, "
      f"max {np.nanmax([e['mult_max'] for e in EV]):.0f} of the {CEIL_PER_WINDOW} ceiling; "
      f"fraction of (station, window) pairs AT the ceiling <= {np.nanmax(fc_):.1e} "
      "-> the count cannot be saturating (and saturation would LOWER b, not raise it)")

# ==================================================================================
# 9. Write the tables
# ==================================================================================
ev_out = pd.concat([ev_df, pd.DataFrame(sum_rows), pd.DataFrame(pair_rows),
                    pd.DataFrame(gr_rows)], ignore_index=True)
ev_out.to_csv(OUT_RES / "s15_cand_exponents.csv", index=False)

cov_ev = []
for e in EV:
    yrs = sorted({int(str(d)[:4]) for d in e["days"]})
    cov_ev.append(dict(kind="event", id=e["pid"], start=e["start"], end=e["end"], seg=e["seg"],
                       years=",".join(str(y) for y in yrs),
                       n_sta=len(e["stas"]), n_sta_within_150km=e["n_sta_ok"],
                       n_cand_window=int(e["meas"]["C_cand"].sum()),
                       n_cand_nearset_allspan=e["n_near_cand"],
                       n_pnsn_window=int(e["meas"]["C_pnsn"].sum()),
                       n_pnsn_nearset_allspan=e["n_near_pnsn"],
                       cand_per_pnsn=(e["meas"]["C_cand"].sum()
                                      / max(e["meas"]["C_pnsn"].sum(), 1)),
                       lat_centroid=e["clat"], lon_centroid=e["clon"],
                       n_tri_F=e["n_tri_F"], nW=e["nW"]))
cov = pd.concat([pd.DataFrame(cov_ev), pd.DataFrame(cov_rows), sta_df], ignore_index=True)
cov.to_csv(OUT_RES / "s15_cand_coverage.csv", index=False)

print("\nCOVERAGE per event")
print("     id  start        seg      years        n_sta  n_cand   n_pnsn   cand/pnsn")
for c in cov_ev:
    print(f"  {c['id']:5d}  {c['start']}  {c['seg']:7s} {c['years']:11s}  {c['n_sta']:4d}  "
          f"{c['n_cand_window']:7d}  {c['n_pnsn_window']:6d}   {c['cand_per_pnsn']:8.2f}")
ns = np.array([c["n_sta"] for c in cov_ev])
nc = np.array([c["n_cand_window"] for c in cov_ev])
yr = np.array([int(c["start"][:4]) for c in cov_ev])
print(f"  stations per event: median {np.median(ns):.0f} [{ns.min()}, {ns.max()}]; "
      f"candidates per event: median {np.median(nc):.0f} [{nc.min()}, {nc.max()}]; "
      f"Spearman(n_sta, year) = {spearmanr(ns, yr).statistic:+.2f} "
      f"(p {spearmanr(ns, yr).pvalue:.2g})")

# ==================================================================================
# 10. Figure
# ==================================================================================
log("figure ...")
fig = plt.figure(figsize=(17.5, 12.5))
gs = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.26)

# A: exponent distributions
ax = fig.add_subplot(gs[0, 0])
bins = np.linspace(-0.6, 1.4, 33)
nl = np.concatenate([v for v in NULP["C_cand"].values()])
nl = nl[np.isfinite(nl)]
ax.hist(nl, bins=bins, weights=np.full(nl.size, 1 / max(nl.size, 1)), color="0.75",
        label=f"shifted-candidate null (pooled, n={nl.size})")
nl2 = np.concatenate([v for v in NULS["C_cand"].values()]) if NULS["C_cand"] else np.array([])
nl2 = nl2[np.isfinite(nl2)]
if nl2.size:
    ax.hist(nl2, bins=bins, weights=np.full(nl2.size, 1 / nl2.size), histtype="step",
            color="k", lw=1.2, label=f"null, per-station wrap (n={nl2.size})")
for m, c in (("C_cand", "tab:red"), ("C_pnsn", "tab:blue")):
    v = np.array([OBS[m][e["pid"]] for e in EV], float)
    ax.hist(v, bins=bins, weights=np.full(v.size, 1 / v.size), histtype="step", lw=2.2, color=c,
            label=f"{m}: median {np.nanmedian(v):+.3f}")
ax.axvline(1.0, color="0.3", ls=":", lw=1)
ax.axvline(0.0, color="0.3", lw=0.6)
ax.set_xlabel("exponent b (log count vs log R7)", fontsize=8)
ax.set_ylabel("fraction", fontsize=8)
ax.set_title("A  exponent distributions: PNSN vs picker-ranked candidates", fontsize=9)
ax.legend(fontsize=6.5, loc="upper right")
ax.tick_params(labelsize=7)

# B-D: example event tracks
examples = sorted(EV, key=lambda e: -e["meas"]["C_cand"].sum())[:2]
central = [e for e in EV if e["seg"] == "central"]
examples.append(max(central, key=lambda e: e["meas"]["C_cand"].sum()))
for j, e in enumerate(examples):
    ax = fig.add_subplot(gs[0, 1] if j == 0 else (gs[0, 2] if j == 1 else gs[1, 0]))
    d = e["days"].astype("datetime64[D]").astype("O")
    ax.plot(d, e["meas"]["C_pnsn"] + 1, color="tab:blue", lw=1.4, label="PNSN count")
    ax.plot(d, e["meas"]["C_cand"] + 1, color="tab:red", lw=1.4, label="candidate count")
    ax.plot(d, e["meas"]["W_cand"] + 1, color="tab:orange", lw=1.0, ls="--",
            label="p_lfe-weighted")
    ax.set_yscale("log")
    ax.tick_params(labelsize=6.5)
    ax2 = ax.twinx()
    ax2.plot(d, e["R7"], color="k", lw=1.6, alpha=0.6)
    ax2.set_ylabel("R7 (mm/d km$^2$)", fontsize=7, color="k")
    ax2.tick_params(labelsize=6.5)
    ax.set_title(f"{'BCD'[j]}  event {e['pid']} {e['start']} ({e['seg']}), "
                 f"{len(e['stas'])} stations\nb(cand) {OBS['C_cand'][e['pid']]:+.3f}  "
                 f"b(PNSN) {OBS['C_pnsn'][e['pid']]:+.3f}", fontsize=8.5)
    if j == 0:
        ax.legend(fontsize=6, loc="upper left")
    for lb in ax.get_xticklabels():
        lb.set_rotation(25)

# E: paired per-event comparison
ax = fig.add_subplot(gs[1, 1])
bp = np.array([OBS["C_pnsn"][e["pid"]] for e in EV])
bc = np.array([OBS["C_cand"][e["pid"]] for e in EV])
cols = {"north": "tab:blue", "central": "tab:green"}
for e, x, y in zip(EV, bp, bc):
    ax.plot([x], [y], "o", ms=5, color=cols[e["seg"]])
lim = [min(bp.min(), bc.min()) - 0.1, max(bp.max(), bc.max()) + 0.1]
ax.plot(lim, lim, "k:", lw=1)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("b, PNSN count", fontsize=8)
ax.set_ylabel("b, candidate count", fontsize=8)
ax.set_title(f"E  paired per event (n={len(EV)})\nmedian difference "
             f"{np.median(bc - bp):+.3f}, Wilcoxon p "
             f"{wilcoxon(bc - bp).pvalue:.2g}", fontsize=9)
ax.tick_params(labelsize=7)

# F: coverage
ax = fig.add_subplot(gs[1, 2])
yy = [pd.Timestamp(e["start"]) for e in EV]
ax.plot(yy, [len(e["stas"]) for e in EV], "o-", color="tab:red", label="stations contributing")
ax.set_ylabel("stations within 150 km with candidates", fontsize=8, color="tab:red")
ax.tick_params(labelsize=7)
ax2 = ax.twinx()
ax2.semilogy(yy, [e["meas"]["C_cand"].sum() for e in EV], "s--", color="0.35",
             label="candidates in window")
ax2.semilogy(yy, [e["meas"]["C_pnsn"].sum() for e in EV], "^--", color="tab:blue",
             label="PNSN detections")
ax2.set_ylabel("count in window", fontsize=8)
ax2.tick_params(labelsize=7)
ax.set_title("F  coverage per event (the network grows x8)", fontsize=9)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=6.5, loc="upper left")
for lb in ax.get_xticklabels():
    lb.set_rotation(25)

# G: lower vs upper half
ax = fig.add_subplot(gs[2, 0])
for m, c in (("C_cand", "tab:red"), ("C_pnsn", "tab:blue")):
    lo_ = ev_df[f"b_{m}_lower"].values
    hi_ = ev_df[f"b_{m}_upper"].values
    ax.plot(lo_, hi_, "o", ms=5, color=c, label=f"{m}: median upper-lower "
            f"{np.nanmedian(hi_ - lo_):+.3f}")
lim = [-1.2, 1.6]
ax.plot(lim, lim, "k:", lw=1)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("b, lower half of log10 R7", fontsize=8)
ax.set_ylabel("b, upper half", fontsize=8)
ax.set_title("G  saturation check (points above the line = no flattening)", fontsize=9)
ax.legend(fontsize=6.5, loc="upper left")
ax.tick_params(labelsize=7)

# H: per-station ETS/quiet ratio
ax = fig.add_subplot(gs[2, 1])
v = sta_df.ets_quiet_ratio.values
v = v[np.isfinite(v)]
ax.hist(v, bins=np.linspace(0, max(6, np.nanpercentile(v, 99)), 40), color="tab:red", alpha=0.8)
ax.axvline(1.0, color="k", lw=1)
ax.axvline(np.median(v), color="tab:red", ls="--", lw=1.5,
           label=f"median {np.median(v):.2f}")
ax.axvline(1.24, color="0.4", ls=":", lw=1.5, label="S14 mf_*.csv 1.24")
ax.axvline(r_pnsn_med, color="tab:blue", ls="-.", lw=1.5,
           label=f"PNSN, same station bbox {r_pnsn_med:.2f}")
ax.set_xlabel("mean daily candidates, ETS windows / outside", fontsize=8)
ax.set_ylabel("stations", fontsize=8)
ax.set_title(f"H  per-station ETS/quiet contrast (n={v.size} stations)", fontsize=9)
ax.legend(fontsize=6.5)
ax.tick_params(labelsize=7)

# I: the decomposition
ax = fig.add_subplot(gs[2, 2])
labels = ["C_pnsn\n(S9/S11)", "C_cand\n(candidates)", "W_cand\n(p_lfe-wt)",
          "M_cand\n(multiplicity)", "C_s10\n(10 stns)"]
vals = [np.array([OBS[m][e["pid"]] for e in EV], float) for m in
        ("C_pnsn", "C_cand", "W_cand", "M_cand", "C_s10")]
bpx = ax.boxplot(vals, tick_labels=labels, showmeans=False, widths=0.6, patch_artist=True)
for patch, c in zip(bpx["boxes"], ("tab:blue", "tab:red", "tab:orange", "0.6", "tab:purple")):
    patch.set_facecolor(c)
    patch.set_alpha(0.6)
ax.axhline(0, color="k", lw=0.6)
ax.axhline(1, color="0.3", ls=":", lw=1)
ax.set_ylabel("exponent b vs log10 R7", fontsize=8)
ax.set_title("I  decomposition: b(C_cand) = b(C_pnsn) + b(multiplicity)", fontsize=9)
ax.tick_params(labelsize=6.5)

fig.suptitle("S15B — picker-ranked candidate counts (data/<sta>_cand_filtered.parquet) vs the "
             "modeled slip rate R7, 27 confirmed events\n"
             "CAVEAT: candidates exist only inside PNSN tremor windows (PNSN sets WHEN to look; "
             "the station's envelope peaks set the count) -> SEMI-independent",
             fontsize=10.5, y=0.985)
fig.savefig(OUT_FIG / "s15_candidate_scaling.png", dpi=115, bbox_inches="tight")
print("\nwrote results/s15_cand_exponents.csv, results/s15_cand_coverage.csv, "
      "figures/s15_candidate_scaling.png")
log("TASK B DONE")
