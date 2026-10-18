#!/usr/bin/env python3
"""D1 - separability diagnostic: can an LFE event catalogue be recovered from the STORED densify
matched-filter hits by thresholding on cc, or must the detection be redone from waveforms?

Background (tremor_slip_geometry S14, 2026-09-10).  `data/mf_<sta>p90f40_<year>.csv`
(columns template,time,cc,station) are the densify matched-filter hits at the FIXED pipeline
threshold cc >= 0.80 with a 6-s minimum gap.  At that setting the filter fires ~100 times per
template per day at rest, the ETS / non-ETS ratio of the raw daily count is only ~1.24 and its
Spearman correlation with the PNSN daily tremor count in the same 40-km neighbourhood is +0.09.
As stored the files are therefore STACKING material, not an event list: real LFEs are a ~20 %
excess riding on a near-constant false-alarm floor.

WHAT COUNTS AS "NOISE" - the central subtlety (user amendment, 2026-09-10)
--------------------------------------------------------------------------
LFEs also occur BETWEEN ETS episodes.  Low-tremor ("quiet") weeks are therefore noise PLUS
inter-ETS LFEs, not pure noise, so:
  * an ETS-minus-QUIET excess UNDER-states the real LFE population, and
  * a "precision" computed against the quiet histogram is biased low by however many real
    inter-ETS LFEs the quiet reference contains.
Every separability number below is therefore reported against FOUR reference floors, and each
number says which floor it uses:
  Q25     quiet-day histogram, local 30-d tremor < 25th percentile  (the original spec; NOT pure)
  Q10     quiet-day histogram, local 30-d tremor < 10th percentile  (purer, still not pure)
  NF-exp  NOISE-ONLY: exponential fit to log(rate) over cc 0.80-0.85 ONLY - the lowest bins,
          where false alarms outnumber any plausible LFE population by orders of magnitude -
          extrapolated upward.  Assumes an exponential tail but never that quiet weeks are pure.
  NF-emp  NOISE-ONLY: the pooled cc-histogram SHAPE of the UNRESPONSIVE templates (ETS/Q25 < 1.1,
          i.e. templates with no tremor response, presumed pure noise), taken per station,
          normalised on 0.80-0.85 and rescaled to each template's own 0.80-0.85 rate.  Makes no
          functional-form assumption.
The per-template threshold recommended for D2 is derived from the NOISE-ONLY floor at a target
false-alarm rate (`thr_d2` = NF-emp 1/day), NOT from the ETS-minus-quiet excess.

METHOD
------
1. Day grid 2009-01-01 .. 2026-12-31.  PNSN tremor (`catalogs/pnsn_tremor_cascadia_full.csv`) is
   binned to daily counts; for every template SOURCE CELL (the "LAT_LON" prefix of the template
   name) the LOCAL tremor series is the daily count of PNSN detections within 40 km of the cell,
   trailing-30-day summed.  Percentiles are taken over 2010-2026.
     ETS days = local 30-d count > p90;  Q25 = < p25;  Q10 = < p10  (<= if the percentile is 0)
   Days between p25 and p90 are excluded.  Q10 is a subset of Q25.
2. Every `mf_*.csv` of the 12 stations already streamed in S14 is streamed once with pyarrow and
   reduced on the fly to
     H [template, year, condition, cc-bin]  cc bins of 0.005 over 0.80-1.00 (40 bins)
     D [template, day]                      total detections per template per day
     Dc[template, day, coarse cc bin]       0.01-wide bins over 0.85-1.00 (+1 catch-all below),
                                            i.e. daily counts above any 0.01-quantised threshold -
                                            this is what makes the burstiness test possible
   Detections with cc >= 0.90 on ETS or quiet days are buffered WITH THEIR TIMES for the
   network-coincidence pilot.  Every rate has as its denominator the number of days on which the
   STATION produced detections at all, restricted to the template's own first..last active day.
3. Two estimators of the ETS excess: POOLED (whole record) and YEAR-MATCHED (the ETS-minus-quiet
   difference formed WITHIN each calendar year with >= 15 ETS and >= 15 quiet days, then averaged),
   because the quiet floor is not perfectly stationary.
4. Per template: rates; the responsiveness ratio; precision/recall vs threshold against each of the
   four floors; both noise-only floor models and their target-false-alarm thresholds; the
   INTER-ETS LFE CANDIDATE population (quiet-day detections above the noise-only floor) and
   whether it is bursty and tremor-correlated (as real LFEs would be) or Poisson-flat (as false
   alarms would be); and the stationarity of the floor in time.
5. Aggregates per responsiveness class, station and segment (source-cell latitude: north >= 46,
   central 44-46, south < 44); a network-coincidence pilot with a LAG SCAN (inter-station S-wave
   travel-time differences can exceed the +-6 s coincidence window, so the window is centred on
   the empirically measured peak lag as well as on zero); and the figure.

DEVIATIONS FROM THE TASK SPEC, deliberate and flagged:
  * the stream uses pyarrow's incremental CSV reader (64 MB blocks) rather than
    `pandas.read_csv(chunksize=2e6)`; it is the code S14 validated on exactly these files and is
    ~3x faster on the same three columns.
  * the year-matched estimator, the softer separability diagnostics (max attainable precision, the
    threshold at precision >= 0.5) and the coincidence lag scan are additions, so that a negative
    verdict is quantitative rather than a column of NaNs.
  * stage 1 is cached to the session scratchpad (`--use-cache`) so a change to the analysis costs
    no re-stream.

Outputs (all inside lfe_catalog_repair/):
  results/d1_template_table.csv    one row per (station, template)
  results/d1_stacked_hist.csv      cc-bin x condition x responsiveness class
  results/d1_floor_shapes.csv      the per-station empirical noise-floor shapes
  results/d1_coincidence.csv       network-coincidence pilot incl. the lag scan
  results/d1_station_summary.csv   per station, per segment, per responsiveness class
  figures/d1_cc_separability.png   15 panels
"""
from __future__ import annotations

import glob
import os
import pickle
import re
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pv
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)
pa.set_cpu_count(4)
pa.set_io_thread_count(4)

ROOT = Path("/home/jovyan/tremorferometry")
HERE = ROOT / "lfe_catalog_repair"
OUT_RES = HERE / "results"
OUT_FIG = HERE / "figures"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)
SCRATCH = Path("/tmp/claude-1000/-home-jovyan-tremorferometry/"
               "03499166-a677-42d8-9397-e2cf5e28621f/scratchpad")
CACHE = SCRATCH / "s14_cache"            # S14's template vocabularies (read-only here)
D1CACHE = SCRATCH / "d1_cache"           # this script's stage-1 cache
D1CACHE.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- configuration
DAY0_STR, DAY1_STR = "2009-01-01", "2026-12-31"
YEAR0, YEAR1 = 2009, 2026
PCT_FROM = "2010-01-01"
TREMOR_RADIUS_KM = 40.0
ROLL = 30
P_ETS, P_Q25, P_Q10 = 90.0, 25.0, 10.0
MIN_COND_DAYS = 20
MIN_YM_DAYS = 15
NBIN = 40
BIN_LO, BIN_W = 0.80, 0.005
BIN_EDGE = BIN_LO + BIN_W * np.arange(NBIN)
BIN_CTR = BIN_EDGE + BIN_W / 2
THR = np.round(np.arange(0.800, 0.9801, 0.005), 4)      # precision/recall grid (spec)
THR_BIN = np.rint((THR - BIN_LO) / BIN_W).astype(int)
THRX = np.round(np.arange(0.800, 0.9951, 0.005), 4)     # extended grid for FA thresholds
THRX_BIN = np.rint((THRX - BIN_LO) / BIN_W).astype(int)
PREC_TARGET, PREC_SOFT = 0.80, 0.50
RESP_HI, RESP_LO = 1.5, 1.1
FIT_LO_HI = 0.85          # NF-exp / NF-emp anchor band: cc 0.80-0.85 (noise-dominated)
FIT_Q_HI = 0.90           # the original spec's 0.80-0.90 quiet-tail fit, kept for comparison
MIN_UNRESP_FOR_SHAPE = 20
NCB = 16                  # coarse daily grid: bin 0 = cc<0.85, bins 1..15 = 0.85..0.99 by 0.01
CTHR = np.round(np.arange(0.85, 0.9951, 0.01), 4)          # 15 usable coarse thresholds
CC_BUF = 0.90
COINC_S = 6.0
LAG_MAX_S = 60.0
SHIFT_S = 600.0
N_CELLS_PILOT = 5
MIN_BOTH_DAYS = 20
FLUSH_ROWS = 8_000_000
COS47 = np.cos(np.deg2rad(47.0))
BUDGET_GB = 150.0

STREAM_STATIONS = [
    ("WISH", "wishp90f40"), ("RADR", "radrp90f40"), ("B943", "b943p90f40"),
    ("LEBA", "lebap90f40"), ("PR01", "pr01p90f40"), ("CARB", "carbp90f40"),
    ("ANGUS", "angusp90f40"), ("HEBO", "hebop90f40"), ("YACT", "yactp90f40"),
    ("SSO2", "sso2p90f40"), ("ALSE", "alsep90f40"), ("WLOO", "wloop90f40"),
]
if os.environ.get("D1_STATIONS"):
    _keep = set(os.environ["D1_STATIONS"].split(","))
    STREAM_STATIONS = [s for s in STREAM_STATIONS if s[0] in _keep]
USE_CACHE = "--use-cache" in sys.argv     # reuse the stage-1 cache, stream nothing

FAMPAT = re.compile(r"^(-?\d+\.\d+)_(-?\d+\.\d+)__c(\d+)$")
T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def to_km(lon, lat):
    return (np.asarray(lon, float) + 124.0) * 111.0 * COS47, (np.asarray(lat, float) - 47.0) * 111.0


def seg_of_lat(lat: float) -> str:
    return "north" if lat >= 46.0 else ("central" if lat >= 44.0 else "south")


# ================================================================ 1. day grid
DAYS = pd.date_range(DAY0_STR, DAY1_STR, freq="D")
NT = len(DAYS)
NYEAR = YEAR1 - YEAR0 + 1
EPOCH0 = int((np.datetime64(DAY0_STR) - np.datetime64("1970-01-01")).astype(int))
PCT_I0 = int((np.datetime64(PCT_FROM) - np.datetime64(DAY0_STR)).astype(int))
YEAR = DAYS.year.values.astype(np.int32)
YI = (YEAR - YEAR0).astype(np.int64)
log(f"day grid {DAY0_STR}..{DAY1_STR} ({NT} days, {NYEAR} years)")


# ================================================================ 2. PNSN tremor -> conditions
tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["di"] = (pd.to_datetime(tr.time).values.astype("datetime64[D]").astype(int)) - EPOCH0
tr = tr[(tr.di >= 0) & (tr.di < NT) & np.isfinite(tr.lat) & np.isfinite(tr.lon)]
TRX, TRY = to_km(tr.lon.values, tr.lat.values)
TRDI = tr.di.values.astype(np.int64)
TR_TREE = cKDTree(np.column_stack([TRX, TRY]))
log(f"PNSN tremor: {len(tr):,} detections {DAYS[TRDI.min()].date()}..{DAYS[TRDI.max()].date()}")

_CELL_CACHE: dict[tuple[float, float], tuple] = {}


def cell_condition(clat: float, clon: float):
    """(cond[NT], daily local tremor count, trailing-30-d count, info) for one source cell.

    cond: 0 excluded / 1 Q25-but-not-Q10 / 2 ETS / 3 Q10.  Q25 = (cond==1)|(cond==3).
    """
    key = (round(clat, 4), round(clon, 4))
    if key in _CELL_CACHE:
        return _CELL_CACHE[key]
    cx, cy = to_km(clon, clat)
    idx = TR_TREE.query_ball_point([float(cx), float(cy)], TREMOR_RADIUS_KM)
    daily = (np.bincount(TRDI[np.asarray(idx, dtype=np.int64)], minlength=NT).astype(np.float32)
             if len(idx) else np.zeros(NT, np.float32))
    cs = np.cumsum(daily.astype(np.float64))
    t30 = (cs - np.concatenate([np.zeros(ROLL), cs[:-ROLL]])).astype(np.float32)
    ref = t30[PCT_I0:]
    p_ets = float(np.percentile(ref, P_ETS))
    p25 = float(np.percentile(ref, P_Q25))
    p10 = float(np.percentile(ref, P_Q10))
    cond = np.zeros(NT, np.int8)
    if p_ets > 0:
        cond[t30 > p_ets] = 2
        cond[(t30 <= p25) if p25 <= 0 else (t30 < p25)] = 1
        cond[(t30 <= p10) if p10 <= 0 else (t30 < p10)] = 3
    info = dict(n_tremor=int(len(idx)), p_ets=p_ets, p_quiet=p25, p_quiet10=p10,
                tremor_ok=bool(p_ets > 0))
    _CELL_CACHE[key] = (cond, daily, t30, info)
    return _CELL_CACHE[key]


# ================================================================ 3. stream
def date_to_di(col):
    d = pc.cast(pc.strptime(pc.utf8_slice_codeunits(col, 0, 10), format="%Y-%m-%d", unit="s"),
                pa.date32())
    a = pc.cast(d, pa.int32())
    if isinstance(a, pa.ChunkedArray):
        a = a.combine_chunks()
    return np.asarray(a).astype(np.int64) - EPOCH0


def dict_encode(col):
    d = pc.dictionary_encode(col)
    if isinstance(d, pa.ChunkedArray):
        d = d.combine_chunks()
    return np.asarray(d.indices).astype(np.int32), d.dictionary.to_pylist()


def sec_of_day(col):
    hh = np.asarray(pc.cast(pc.utf8_slice_codeunits(col, 11, 13), pa.float64()))
    mm = np.asarray(pc.cast(pc.utf8_slice_codeunits(col, 14, 16), pa.float64()))
    ss = np.asarray(pc.cast(pc.utf8_slice_codeunits(col, 17, 29), pa.float64()))
    return hh * 3600.0 + mm * 60.0 + ss


def prefill_vocab(stem: str):
    f = CACHE / f"{stem}_v2.npz"
    if not f.exists():
        return []
    try:
        return [str(x) for x in np.load(f, allow_pickle=True)["fams"]]
    except Exception:                                          # noqa: BLE001
        return []


def stream_station(sta: str, stem: str):
    """Single pass over mf_<stem>_<year>.csv -> H, D, Dc, coincidence buffer, template metadata."""
    files = sorted(f for f in glob.glob(str(ROOT / f"data/mf_{stem}_*.csv"))
                   if re.search(r"_(\d{4})\.csv$", f))
    vocab0 = prefill_vocab(stem)
    NF = max(100, len(vocab0) + 25)          # sized from the S14 vocabulary, with slack
    fam_index: dict[str, int] = {}
    fams: list[str] = []
    COND = np.zeros((NF, NT), np.int8)
    cinfo: list[dict] = []
    HS = NYEAR * 4 * NBIN
    H = np.zeros(NF * HS, np.int64)
    D = np.zeros((NF, NT), np.int32)
    Dc = np.zeros(NF * NT * NCB, np.int32)

    def add_template(name: str) -> int:
        j = len(fam_index)
        if j >= NF:
            raise RuntimeError(f"{sta}: more than {NF} templates")
        m = FAMPAT.match(name)
        if m is None:
            raise RuntimeError(f"{sta}: template name {name!r} does not parse")
        clat, clon = float(m.group(1)), float(m.group(2))
        cond, _daily, _t30, info = cell_condition(clat, clon)
        COND[j] = cond
        cinfo.append(dict(template=name, cell_lat=clat, cell_lon=clon, **info))
        fam_index[name] = j
        fams.append(name)
        return j

    for name in vocab0:
        if name not in fam_index:
            add_template(name)
    n_pre = len(fam_index)

    buf_c, buf_d, buf_cc = [], [], []
    nbuf = 0
    bc, bd, bt, bcc = [], [], [], []
    nbytes = ndet = 0

    def flush():
        nonlocal buf_c, buf_d, buf_cc, nbuf
        if not nbuf:
            return
        c = np.concatenate(buf_c).astype(np.int64)
        d = np.concatenate(buf_d).astype(np.int64)
        v = np.concatenate(buf_cc)
        flat = c * NT + d
        D.ravel()[:] += np.bincount(flat, minlength=NF * NT).astype(np.int32)
        cd = COND.ravel()[flat].astype(np.int64)
        b = np.clip(((v - BIN_LO) / BIN_W).astype(np.int64), 0, NBIN - 1)
        H[:] += np.bincount(c * HS + YI[d] * (4 * NBIN) + cd * NBIN + b, minlength=NF * HS)
        cb = np.clip(np.floor((v - 0.85) / 0.01).astype(np.int64) + 1, 0, NCB - 1)
        Dc[:] += np.bincount(flat * NCB + cb, minlength=NF * NT * NCB).astype(np.int32)
        buf_c, buf_d, buf_cc = [], [], []
        nbuf = 0

    for f in files:
        nbytes += os.path.getsize(f)
        rd = pv.open_csv(f, read_options=pv.ReadOptions(block_size=64 << 20),
                         convert_options=pv.ConvertOptions(
                             include_columns=["template", "time", "cc"],
                             column_types={"template": pa.string(), "time": pa.string(),
                                           "cc": pa.float32()}))
        for b in rd:
            di = date_to_di(b.column(1))
            cc = np.asarray(b.column(2))
            loc, vocab = dict_encode(b.column(0))
            g = np.empty(len(vocab), np.int32)
            for k, name in enumerate(vocab):
                g[k] = fam_index[name] if name in fam_index else add_template(name)
            codes = g[loc]
            keep = (di >= 0) & (di < NT) & np.isfinite(cc)
            if not keep.any():
                continue
            ck, dk, vk = codes[keep], di[keep].astype(np.int32), cc[keep]
            buf_c.append(ck)
            buf_d.append(dk)
            buf_cc.append(vk)
            nbuf += ck.size
            ndet += int(ck.size)
            selb = (vk >= CC_BUF) & (COND[ck, dk] > 0)
            if selb.any():
                wi = np.flatnonzero(keep)[selb]
                bc.append(ck[selb].astype(np.int16))
                bd.append(dk[selb])
                bt.append(dk[selb].astype(np.float64) * 86400.0
                          + sec_of_day(b.column(1).take(pa.array(wi))))
                bcc.append(vk[selb])
            if nbuf >= FLUSH_ROWS:
                flush()
        del rd
    flush()
    nf = len(fam_index)
    BUF = dict(code=np.concatenate(bc) if bc else np.zeros(0, np.int16),
               day=np.concatenate(bd) if bd else np.zeros(0, np.int32),
               t=np.concatenate(bt) if bt else np.zeros(0),
               cc=np.concatenate(bcc) if bcc else np.zeros(0, np.float32))
    Hn = H.reshape(NF, NYEAR, 4, NBIN)[:nf]
    Dcn = Dc.reshape(NF, NT, NCB)[:nf]
    assert int(Hn.sum()) == int(D.sum()) == int(Dcn.sum()), f"{sta}: count mismatch"
    return (fams, pd.DataFrame(cinfo), COND[:nf], Hn, D[:nf], Dcn, BUF,
            nbytes, ndet, len(files), n_pre)


# ================================================================ 4. helpers
def hist_median(h: np.ndarray) -> float:
    tot = np.nansum(h)
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    return float(BIN_CTR[min(int(np.searchsorted(np.cumsum(h), 0.5 * tot)), NBIN - 1)])


def cum_median(exc: np.ndarray, tot: float) -> float:
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    k = np.flatnonzero(np.cumsum(exc) >= 0.5 * tot)
    return float(BIN_CTR[k[0]]) if k.size else np.nan


def tailsum(h: np.ndarray, ext: bool = False) -> np.ndarray:
    """Sum of h over bins with left edge >= t, for every t in THR (or THRX if ext)."""
    cs = np.concatenate([np.cumsum(h[::-1])[::-1], [0.0]])
    return cs[THRX_BIN] if ext else cs[THR_BIN]


def separability(h_sig: np.ndarray, h_floor: np.ndarray) -> dict:
    """Precision / recall of (signal - floor) against threshold, for ANY choice of floor.

    precision(t) = (sig - floor) above t / sig above t   -> the fraction of the detections kept
                                                            at that threshold that are real
    recall(t)    = (sig - floor) above t / (sig - floor) total
    """
    exc = h_sig - h_floor
    tot = float(np.nansum(exc))
    exc_above, sig_above = tailsum(np.nan_to_num(exc)), tailsum(np.nan_to_num(h_sig))
    with np.errstate(divide="ignore", invalid="ignore"):
        recall = exc_above / tot if tot > 0 else np.full(THR.size, np.nan)
        prec = np.where(sig_above > 0, exc_above / sig_above, np.nan)
    out = dict(exc=exc, exc_tot=tot, exc_above=exc_above, sig_above=sig_above,
               recall=recall, prec=prec,
               prec_max=float(np.nanmax(prec)) if np.isfinite(prec).any() else np.nan)
    for tag, target in (("p80", PREC_TARGET), ("p50", PREC_SOFT)):
        ok = np.isfinite(prec) & (prec >= target) & np.isfinite(recall)
        if ok.any():
            k = int(np.flatnonzero(ok)[np.argmax(recall[ok])])
            out[f"k_{tag}"], out[f"thr_{tag}"] = k, float(THR[k])
            out[f"rec_{tag}"], out[f"prc_{tag}"] = float(recall[k]), float(prec[k])
        else:
            out[f"k_{tag}"] = -1
            out[f"thr_{tag}"] = out[f"rec_{tag}"] = out[f"prc_{tag}"] = np.nan
    return out


def exp_tail_fit(h: np.ndarray, hi: float):
    """Exponential fit of log(rate per bin) vs cc over bins with centre < hi.

    Returns (slope, intercept, R2, predicted rate per bin over ALL bins).
    """
    m = (BIN_CTR < hi) & (h > 0) & np.isfinite(h)
    if m.sum() < 6:
        return np.nan, np.nan, np.nan, np.full(NBIN, np.nan)
    x, y = BIN_CTR[m], np.log(h[m])
    a, b = np.polyfit(x, y, 1)
    yy = a * x + b
    r2 = 1.0 - np.sum((y - yy) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-30)
    return float(a), float(b), float(r2), np.exp(np.clip(a * BIN_CTR + b, -700, 700))


ANCHOR = BIN_EDGE < FIT_LO_HI          # the 0.80-0.85 anchor band (10 bins)


def first_below(rate_above: np.ndarray, target: float, grid: np.ndarray = THRX):
    """Smallest threshold on `grid` whose predicted above-threshold rate is < target.

    Returns (nan, -1) when the modelled floor never falls below `target` anywhere on the grid -
    which is itself a result: no cc threshold reaches that false-alarm rate for that template.
    """
    i = np.flatnonzero(np.isfinite(rate_above) & (rate_above < target))
    return (float(grid[i[0]]), int(i[0])) if i.size else (np.nan, -1)


# ================================================================ 5. stage 1 (stream or cache)
STAGE1 = D1CACHE / ("stage1_" + "_".join(s for s, _ in STREAM_STATIONS) + ".pkl")
if USE_CACHE and STAGE1.exists():
    log(f"loading stage-1 cache {STAGE1.name} ({STAGE1.stat().st_size / 1e9:.2f} GB, no streaming)")
    with open(STAGE1, "rb") as fh:
        S1 = pickle.load(fh)
else:
    log("streaming ...")
    S1 = dict(cov=[], hist={}, cond={}, act={}, buf={}, code={})
    gb = 0.0
    for sta, stem in STREAM_STATIONS:
        if gb > BUDGET_GB:
            log(f"  BUDGET {BUDGET_GB:.0f} GB reached -> skipping {sta}")
            continue
        t_s = time.time()
        fams, CI, COND, H, D, Dc, BUF, nbytes, ndet, nfiles, n_pre = stream_station(sta, stem)
        gb += nbytes / 1e9
        wall = time.time() - t_s
        act = D.sum(axis=0) > 0
        aidx = np.flatnonzero(act)
        log(f"  {sta:6s} {nfiles:2d} files {nbytes / 1e9:6.2f} GB {ndet:13,d} det  "
            f"{len(fams):3d} templates ({n_pre} pre-filled)  {act.sum():5d} active days  "
            f"buf {BUF['code'].size:,}  {wall:6.1f}s")
        S1["cov"].append(dict(sta=sta, stem=stem, files=nfiles, GB=nbytes / 1e9, detections=ndet,
                              templates=len(fams), active_days=int(act.sum()),
                              first_day=str(DAYS[aidx[0]].date()) if aidx.size else "",
                              last_day=str(DAYS[aidx[-1]].date()) if aidx.size else "",
                              buffered=int(BUF["code"].size), wall_s=round(wall, 1)))
        S1["cond"][sta], S1["act"][sta], S1["buf"][sta] = COND, act, BUF
        S1["code"][sta] = {n: j for j, n in enumerate(fams)}
        Dcum = np.cumsum(Dc[:, :, ::-1], axis=2)[:, :, ::-1]     # counts above each coarse edge
        for j, name in enumerate(fams):
            ci = CI.iloc[j]
            d = D[j]
            nz = np.flatnonzero(d > 0)
            if nz.size == 0 or not bool(ci.tremor_ok):
                continue
            span = np.zeros(NT, bool)
            span[nz[0]:nz[-1] + 1] = True
            use = act & span
            cond = COND[j]
            m_e = use & (cond == 2)
            m_q25 = use & ((cond == 1) | (cond == 3))
            m_q10 = use & (cond == 3)
            n_e, n_q25, n_q10 = int(m_e.sum()), int(m_q25.sum()), int(m_q10.sum())
            if n_e < MIN_COND_DAYS or n_q25 < MIN_COND_DAYS:
                continue
            Hj = H[j].astype(np.float64)
            _c, _daily, _t30, _i = cell_condition(float(ci.cell_lat), float(ci.cell_lon))
            S1["hist"][(sta, name)] = dict(
                ets=Hj[:, 2].sum(axis=0), q25=(Hj[:, 1] + Hj[:, 3]).sum(axis=0),
                q10=Hj[:, 3].sum(axis=0),
                ets_y=Hj[:, 2], q25_y=Hj[:, 1] + Hj[:, 3],
                ne_y=np.bincount(YI[m_e], minlength=NYEAR),
                nq_y=np.bincount(YI[m_q25], minlength=NYEAR),
                n_e=n_e, n_q25=n_q25, n_q10=n_q10,
                dq=Dcum[j][m_q25, 1:].astype(np.int32),     # daily counts above each coarse thr
                de=Dcum[j][m_e, 1:].astype(np.int32),
                trem_q_daily=_daily[m_q25].copy(), trem_q_30d=_t30[m_q25].copy(),
                cell_lat=float(ci.cell_lat), cell_lon=float(ci.cell_lon),
                n_tremor=int(ci.n_tremor), p_ets=float(ci.p_ets),
                p_q25=float(ci.p_quiet), p_q10=float(ci.p_quiet10))
        del H, D, Dc, Dcum
    with open(STAGE1, "wb") as fh:
        pickle.dump(S1, fh, protocol=4)
    log(f"stage-1 cache written: {STAGE1} ({STAGE1.stat().st_size / 1e9:.2f} GB)")

COV = pd.DataFrame(S1["cov"])
HISTD = S1["hist"]
STA_COND, STA_ACT, STA_BUF, CODE_OF = S1["cond"], S1["act"], S1["buf"], S1["code"]
log(f"streamed/loaded {COV.GB.sum():.1f} GB / {COV.detections.sum():,} detections / "
    f"{COV.templates.sum()} templates; {len(HISTD)} templates pass the day gate")
if not HISTD:
    sys.exit("no templates survived the day gate")


# ================================================================ 6. pass A: responsiveness
base = []
for (sta, name), h in HISTD.items():
    rate_e = h["ets"].sum() / h["n_e"]
    rate_q = h["q25"].sum() / h["n_q25"]
    rate_q10 = h["q10"].sum() / h["n_q10"] if h["n_q10"] >= MIN_COND_DAYS else np.nan
    base.append(dict(station=sta, template=name, rate_ets=rate_e, rate_quiet=rate_q,
                     rate_quiet10=rate_q10,
                     ratio=(rate_e / rate_q if rate_q > 0 else np.nan),
                     ratio_q10=(rate_e / rate_q10
                                if np.isfinite(rate_q10) and rate_q10 > 0 else np.nan)))
BASE = pd.DataFrame(base).set_index(["station", "template"])
log(f"pass A: responsiveness ETS/Q25 median {BASE.ratio.median():.3f}; "
    f"ETS/Q10 median {BASE.ratio_q10.median():.3f}")


def pooled_shape(keys):
    """POOLED quiet-day cc histogram of a set of templates, normalised on the 0.80-0.85 anchor.

    Pooling the COUNTS (rather than averaging per-template normalised shapes) is what the
    amendment specifies and is robust: a handful of templates with freak tails cannot dominate.
    """
    if not keys:
        return None, 0
    tot = np.zeros(NBIN)
    for k in keys:
        tot += HISTD[k]["q25"] / HISTD[k]["n_q25"]      # per-day rates, so long records do not
    a = tot[ANCHOR].sum()                               # simply outweigh short ones
    return (tot / a, len(keys)) if a > 0 else (None, len(keys))


SHAPE: dict[str, np.ndarray] = {}
shape_rows = []
unresp = BASE[BASE.ratio < RESP_LO]
UKEYS = list(unresp.index)
SHAPE["__global__"], n_glob = pooled_shape(UKEYS)
n_fallback = 0
for sta in sorted({k[0] for k in HISTD}):
    kk = [k for k in UKEYS if k[0] == sta]
    sh, n = pooled_shape(kk)
    fb = (n < MIN_UNRESP_FOR_SHAPE) or (sh is None)
    n_fallback += int(fb)
    SHAPE[sta] = SHAPE["__global__"] if fb else sh
    shape_rows.append(dict(station=sta, n_unresponsive_used=n, used_global=fb,
                           **{f"shape_{BIN_EDGE[k]:.3f}": SHAPE[sta][k] for k in range(NBIN)}))
shape_rows.append(dict(station="__global__", n_unresponsive_used=n_glob, used_global=False,
                       **{f"shape_{BIN_EDGE[k]:.3f}": SHAPE["__global__"][k]
                          for k in range(NBIN)}))
pd.DataFrame(shape_rows).to_csv(OUT_RES / "d1_floor_shapes.csv", index=False)
log(f"pass A: NF-emp shapes POOLED from {n_glob} unresponsive templates "
    f"({n_fallback} stations fell back to the global shape)")


# ================================================================ 7. pass B: per-template metrics
ROWS = []
CURVES: dict[tuple[str, str], dict] = {}
for (sta, name), h in HISTD.items():
    b = BASE.loc[(sta, name)]
    n_e, n_q25, n_q10 = h["n_e"], h["n_q25"], h["n_q10"]
    h_e = h["ets"] / n_e
    h_q = h["q25"] / n_q25
    h_q10 = (h["q10"] / n_q10) if n_q10 >= MIN_COND_DAYS else np.full(NBIN, np.nan)
    # ---------- floors.  The noise-only floors are scaled SEPARATELY for ETS days and for quiet
    # days, using that condition's own 0.80-0.85 rate as the noise anchor: the background noise
    # level is not the same in the two conditions, and scaling both from the quiet anchor would
    # charge any change of noise level to "signal".
    aq, bq, r2_q90, fit_q90 = exp_tail_fit(h_q, FIT_Q_HI)        # spec's 0.80-0.90 quiet fit
    a_nf, b_nf, r2_nf, fl_exp_q = exp_tail_fit(h_q, FIT_LO_HI)   # NF-exp on quiet days
    a_nfe, b_nfe, r2_nfe, fl_exp_e = exp_tail_fit(h_e, FIT_LO_HI)   # NF-exp on ETS days
    sh = SHAPE.get(sta)
    if sh is not None:
        fl_emp_q = sh * h_q[ANCHOR].sum()
        fl_emp_e = sh * h_e[ANCHOR].sum()
    else:
        fl_emp_q = fl_emp_e = np.full(NBIN, np.nan)
    # ---------- separability of the ETS-day histogram against each floor
    S_q25 = separability(h_e, h_q)
    S_q10 = separability(h_e, h_q10) if np.isfinite(h_q10).all() else None
    S_nfe = separability(h_e, fl_exp_e)
    S_nfm = separability(h_e, fl_emp_e)
    # ---------- year-matched (Q25 floor)
    gy = np.flatnonzero((h["ne_y"] >= MIN_YM_DAYS) & (h["nq_y"] >= MIN_YM_DAYS))
    if gy.size:
        he_ym = (h["ets_y"][gy] / h["ne_y"][gy, None]).mean(axis=0)
        hq_ym = (h["q25_y"][gy] / h["nq_y"][gy, None]).mean(axis=0)
        SY = separability(he_ym, hq_ym)
    else:
        he_ym = hq_ym = np.full(NBIN, np.nan)
        SY = None
    # ---------- target-false-alarm thresholds from each floor
    fa_q90, fa_nfe, fa_nfm = (tailsum(fit_q90, True), tailsum(fl_exp_q, True),
                              tailsum(fl_emp_q, True))
    thr_q1, k_q1 = first_below(fa_q90, 1.0)
    thr_q01, k_q01 = first_below(fa_q90, 0.1)
    thr_e1, k_e1 = first_below(fa_nfe, 1.0)
    thr_e01, k_e01 = first_below(fa_nfe, 0.1)
    thr_m1, k_m1 = first_below(fa_nfm, 1.0)
    thr_m01, k_m01 = first_below(fa_nfm, 0.1)
    obs_e, obs_q = tailsum(h_e, True), tailsum(h_q, True)
    i090 = int(np.flatnonzero(THRX == 0.90)[0])
    # ---------- inter-ETS LFE candidates: quiet-day detections above the NOISE-ONLY floor.
    # thr_d2 comes from NF-emp ONLY.  There is deliberately NO fallback to NF-exp: NF-exp is
    # shown below to under-predict the measured noise floor above cc 0.90 by orders of magnitude
    # (the false-alarm tail is heavy, not exponential), so a threshold taken from it would be
    # wildly over-optimistic.  When NF-emp never reaches 1 FA/day even at cc 0.995, thr_d2 is
    # left undefined - that is the answer for that template.
    thr_d2 = thr_m1

    def cand_stats(kk: int, fa: np.ndarray):
        """(floor, quiet excess, ETS excess, burstiness and tremor-link diagnostics) at THRX[kk]."""
        if kk < 0 or not np.isfinite(fa[kk]):
            return dict.fromkeys(
                ["floor", "cq", "ce", "mu", "mu_e", "disp", "f0", "f3", "p0", "rho_d", "rho_30",
                 "r_hi", "r_lo", "tb"], np.nan)
        ci = int(np.argmin(np.abs(CTHR - THRX[kk])))
        nd = h["dq"][:, ci].astype(float)
        mu = float(nd.mean())
        tq_d, tq_30 = h["trem_q_daily"].astype(float), h["trem_q_30d"].astype(float)
        med30 = np.median(tq_30)
        hi_m, lo_m = tq_30 > med30, tq_30 <= med30
        r_hi = float(nd[hi_m].mean()) if hi_m.sum() >= 10 else np.nan
        r_lo = float(nd[lo_m].mean()) if lo_m.sum() >= 10 else np.nan
        return dict(
            floor=float(fa[kk]), cq=float(obs_q[kk] - fa[kk]), ce=float(obs_e[kk] - fa[kk]),
            mu=mu, mu_e=float(h["de"][:, ci].astype(float).mean()),
            disp=(float(nd.var() / mu) if mu > 0 else np.nan),
            f0=float((nd == 0).mean()), f3=float((nd >= 3).mean()),
            p0=(float(np.exp(-mu)) if np.isfinite(mu) else np.nan),
            rho_d=(spearmanr(nd, tq_d).statistic
                   if nd.size > 20 and nd.std() > 0 and tq_d.std() > 0 else np.nan),
            rho_30=(spearmanr(nd, tq_30).statistic
                    if nd.size > 20 and nd.std() > 0 and tq_30.std() > 0 else np.nan),
            r_hi=r_hi, r_lo=r_lo,
            tb=((r_hi / r_lo) if (np.isfinite(r_lo) and r_lo > 0) else np.nan))

    C1 = cand_stats(k_m1, fa_nfm)                       # at the NF-emp 1 FA/day threshold
    k95 = int(np.flatnonzero(THRX == 0.95)[0])          # and at a FIXED cc >= 0.95, always defined
    C95 = cand_stats(k95, fa_nfm)
    # ---------- floor stationarity
    qy = np.flatnonzero(h["nq_y"] >= MIN_COND_DAYS)
    if qy.size >= 4:
        rq = h["q25_y"][qy].sum(axis=1) / h["nq_y"][qy]
        sr = spearmanr(YEAR0 + qy, rq)
        trend, trend_p = float(sr.statistic), float(sr.pvalue)
        spread = float(np.max(rq) / max(np.min(rq), 1e-9))
        cv = float(np.std(rq) / max(np.mean(rq), 1e-12))
    else:
        trend = trend_p = spread = cv = np.nan
    CURVES[(sta, name)] = dict(h_e=h_e, h_q=h_q, h_q10=h_q10,
                               fl_exp=fl_exp_q, fl_emp=fl_emp_q,
                               fl_exp_e=fl_exp_e, fl_emp_e=fl_emp_e,
                               he_ym=he_ym, hq_ym=hq_ym,
                               prec_q25=S_q25["prec"], rec_q25=S_q25["recall"],
                               prec_nfm=S_nfm["prec"], rec_nfm=S_nfm["recall"],
                               prec_nfe=S_nfe["prec"], rec_nfe=S_nfe["recall"])
    ROWS.append(dict(
        station=sta, template=name, cell_lat=h["cell_lat"], cell_lon=h["cell_lon"],
        segment=seg_of_lat(h["cell_lat"]), n_tremor_40km=h["n_tremor"],
        p_ets_30d=h["p_ets"], p_quiet_30d=h["p_q25"], p_quiet10_30d=h["p_q10"],
        n_ets_days=n_e, n_quiet_days=n_q25, n_quiet10_days=n_q10,
        n_det_ets=int(h["ets"].sum()), n_det_quiet=int(h["q25"].sum()),
        n_det_quiet10=int(h["q10"].sum()),
        rate_ets=b.rate_ets, rate_quiet=b.rate_quiet, rate_quiet10=b.rate_quiet10,
        ratio=b.ratio, ratio_q10=b.ratio_q10,
        excess_total=S_q25["exc_tot"],
        excess_total_q10=(S_q10["exc_tot"] if S_q10 else np.nan),
        excess_total_nfe=S_nfe["exc_tot"], excess_total_nfm=S_nfm["exc_tot"],
        median_cc_quiet=hist_median(h_q), median_cc_quiet10=hist_median(h_q10),
        median_cc_ets=hist_median(h_e),
        median_cc_excess=cum_median(S_q25["exc"], S_q25["exc_tot"]),
        median_cc_excess_nfm=cum_median(np.nan_to_num(S_nfm["exc"]), S_nfm["exc_tot"]),
        q_frac_gt_090=float(h_q[BIN_EDGE >= 0.90].sum() / max(h_q.sum(), 1e-12)),
        q_frac_gt_095=float(h_q[BIN_EDGE >= 0.95].sum() / max(h_q.sum(), 1e-12)),
        e_frac_gt_090=float(h_e[BIN_EDGE >= 0.90].sum() / max(h_e.sum(), 1e-12)),
        e_frac_gt_095=float(h_e[BIN_EDGE >= 0.95].sum() / max(h_e.sum(), 1e-12)),
        thr_p80=S_q25["thr_p80"], recall_p80=S_q25["rec_p80"], prec_p80=S_q25["prc_p80"],
        thr_p50=S_q25["thr_p50"], recall_p50=S_q25["rec_p50"], prec_max=S_q25["prec_max"],
        thr_p80_q10=(S_q10["thr_p80"] if S_q10 else np.nan),
        recall_p80_q10=(S_q10["rec_p80"] if S_q10 else np.nan),
        prec_max_q10=(S_q10["prec_max"] if S_q10 else np.nan),
        thr_p80_nfe=S_nfe["thr_p80"], recall_p80_nfe=S_nfe["rec_p80"],
        prec_max_nfe=S_nfe["prec_max"],
        thr_p80_nfm=S_nfm["thr_p80"], recall_p80_nfm=S_nfm["rec_p80"],
        prec_max_nfm=S_nfm["prec_max"],
        n_years_ym=int(gy.size),
        ratio_ym=(float(he_ym.sum() / hq_ym.sum()) if gy.size and hq_ym.sum() > 0 else np.nan),
        excess_total_ym=(SY["exc_tot"] if SY else np.nan),
        thr_p80_ym=(SY["thr_p80"] if SY else np.nan),
        recall_p80_ym=(SY["rec_p80"] if SY else np.nan),
        median_cc_excess_ym=(cum_median(SY["exc"], SY["exc_tot"]) if SY else np.nan),
        median_cc_quiet_ym=(hist_median(hq_ym) if gy.size else np.nan),
        tail_slope=aq, tail_r2=r2_q90, nf_exp_slope=a_nf, nf_exp_r2=r2_nf,
        nf_exp_slope_ets=a_nfe, nf_exp_r2_ets=r2_nfe,
        thr_fa1perday=thr_q1, thr_fa0p1perday=thr_q01,
        thr_nf1_exp=thr_e1, thr_nf01_exp=thr_e01,
        thr_nf1_emp=thr_m1, thr_nf01_emp=thr_m01, thr_d2=thr_d2,
        floor_exp_above_090=float(fa_nfe[i090]), floor_emp_above_090=float(fa_nfm[i090]),
        obs_quiet_above_090=float(obs_q[i090]), obs_ets_above_090=float(obs_e[i090]),
        snr_at_thr=(float(obs_e[k_q1] / fa_q90[k_q1])
                    if k_q1 >= 0 and fa_q90[k_q1] > 0 else np.nan),
        snr_at_thr0p1=(float(obs_e[k_q01] / fa_q90[k_q01])
                       if k_q01 >= 0 and fa_q90[k_q01] > 0 else np.nan),
        cand_rate_quiet=C1["cq"], cand_rate_ets=C1["ce"], floor_rate_at_thr=C1["floor"],
        cand_obs_mean_quiet=C1["mu"], cand_obs_mean_ets=C1["mu_e"],
        cand_dispersion=C1["disp"], cand_frac_days_zero=C1["f0"],
        cand_frac_days_ge3=C1["f3"], cand_poisson_frac_zero=C1["p0"],
        cand_rho_daily_tremor=C1["rho_d"], cand_rho_30d_tremor=C1["rho_30"],
        cand_rate_hi_tremor=C1["r_hi"], cand_rate_lo_tremor=C1["r_lo"],
        cand_topbot_ratio=C1["tb"],
        # the same diagnostics at a FIXED cc >= 0.95, so they exist for every template
        cand95_rate_quiet=C95["cq"], cand95_rate_ets=C95["ce"],
        cand95_floor=C95["floor"], cand95_obs_mean_quiet=C95["mu"],
        cand95_obs_mean_ets=C95["mu_e"], cand95_dispersion=C95["disp"],
        cand95_frac_days_zero=C95["f0"], cand95_frac_days_ge3=C95["f3"],
        cand95_poisson_frac_zero=C95["p0"], cand95_rho_daily_tremor=C95["rho_d"],
        cand95_rho_30d_tremor=C95["rho_30"], cand95_topbot_ratio=C95["tb"],
        floor_trend=trend, floor_trend_p=trend_p, floor_year_spread=spread, floor_cv=cv,
        n_quiet_years=int(qy.size)))

TT = pd.DataFrame(ROWS)
TT["rclass"] = np.where(TT.ratio >= RESP_HI, "responsive",
                        np.where(TT.ratio < RESP_LO, "unresponsive", "mid"))
TT["responsive"] = TT.rclass == "responsive"
TT["unresponsive"] = TT.rclass == "unresponsive"
TT["rclass_ym"] = np.where(TT.ratio_ym >= RESP_HI, "responsive",
                           np.where(TT.ratio_ym < RESP_LO, "unresponsive", "mid"))
for tag in ("", "_q10", "_nfe", "_nfm", "_ym"):
    TT[f"repairable{tag}"] = np.isfinite(TT[f"thr_p80{tag}"]) & (TT[f"recall_p80{tag}"] >= 0.5)
TT["dcc"] = TT.median_cc_excess - TT.median_cc_quiet
TT["dcc_nfm"] = TT.median_cc_excess_nfm - TT.median_cc_quiet
TT["dcc_ym"] = TT.median_cc_excess_ym - TT.median_cc_quiet_ym
TT["floor_model_ratio_090"] = TT.floor_emp_above_090 / TT.floor_exp_above_090
TT["quiet_over_floor_090"] = TT.obs_quiet_above_090 / TT.floor_emp_above_090
# how heavy the tail is relative to an EXPONENTIAL noise model.  Measured on the UNRESPONSIVE
# templates this is a pure-noise number (they have no ETS response, so no LFE population): if it
# is as large there as for responsive templates, the heavy high-cc tail is FALSE ALARMS, and any
# threshold derived from an exponential extrapolation is invalid.
TT["quiet_over_expfloor_090"] = TT.obs_quiet_above_090 / TT.floor_exp_above_090
TT["thr_d2_defined"] = np.isfinite(TT.thr_d2)
TT.sort_values(["station", "template"]).to_csv(OUT_RES / "d1_template_table.csv", index=False)
log(f"pass B: {len(TT)} template rows")


# ================================================================ 8. stacked histograms
srows = []
for cls in ("responsive", "mid", "unresponsive", "all"):
    sub = TT if cls == "all" else TT[TT.rclass == cls]
    if sub.empty:
        continue
    keys = [(r.station, r.template) for r in sub.itertuples()]
    series = {nm: np.stack([CURVES[k][fld] for k in keys])
              for nm, fld in (("ets", "h_e"), ("quiet", "h_q"), ("quiet10", "h_q10"),
                              ("floor_nf_exp", "fl_exp"), ("floor_nf_emp", "fl_emp"),
                              ("floor_nf_exp_ets", "fl_exp_e"),
                              ("floor_nf_emp_ets", "fl_emp_e"),
                              ("ets_ym", "he_ym"), ("quiet_ym", "hq_ym"))}
    series["excess"] = series["ets"] - series["quiet"]
    series["excess_q10"] = series["ets"] - series["quiet10"]
    series["excess_nf_emp"] = series["ets"] - series["floor_nf_emp_ets"]
    series["excess_ym"] = series["ets_ym"] - series["quiet_ym"]
    series["quiet_minus_floor"] = series["quiet"] - series["floor_nf_emp"]
    for cname, M in series.items():
        for k in range(NBIN):
            srows.append(dict(bin_lo=BIN_EDGE[k], bin_hi=BIN_EDGE[k] + BIN_W, bin_ctr=BIN_CTR[k],
                              condition=cname, rclass=cls, n_templates=len(sub),
                              mean_rate_per_day=float(np.nanmean(M[:, k])),
                              sum_rate_per_day=float(np.nansum(M[:, k])),
                              median_rate_per_day=float(np.nanmedian(M[:, k]))))
SH = pd.DataFrame(srows)
SH.to_csv(OUT_RES / "d1_stacked_hist.csv", index=False)


# ================================================================ 9. coincidence pilot + lag scan
TT["cell"] = TT.cell_lat.round(3).astype(str) + "_" + TT.cell_lon.round(3).astype(str)
TT["thr_use"] = np.maximum(np.where(np.isfinite(TT.thr_d2), TT.thr_d2, 0.95), CC_BUF)
TT["thr_source"] = np.where(np.isfinite(TT.thr_d2), "nf_emp_1perday", "fallback0.95")


def pick_cells(pool: pd.DataFrame):
    """Cells with eligible templates at >= 2 stations, ranked by simultaneously-recording pairs."""
    if pool.empty:
        return pd.DataFrame()
    g = (pool.groupby("cell")
         .agg(n_sta=("station", "nunique"), n_tpl=("template", "size"),
              lat=("cell_lat", "first"), lon=("cell_lon", "first"),
              mean_ratio=("ratio", "mean")).query("n_sta >= 2"))
    npairs = []
    for cell in g.index:
        ss = sorted(pool[pool.cell == cell].station.unique())
        npairs.append(sum(int((STA_ACT[ss[i]] & STA_ACT[ss[j]]).sum()) >= 200
                          for i in range(len(ss)) for j in range(i + 1, len(ss))))
    g["n_overlap_pairs"] = npairs
    return g.sort_values(["n_overlap_pairs", "n_sta", "n_tpl"], ascending=False)


POOL = TT[TT.responsive]
cc_cnt = pick_cells(POOL)
pilot_pool = "responsive"
if len(cc_cnt) < N_CELLS_PILOT:
    thr_r = TT.ratio.quantile(0.75)
    POOL = TT[TT.ratio >= thr_r]
    cc_cnt = pick_cells(POOL)
    pilot_pool = f"top-quartile ratio (>= {thr_r:.3f})"
log(f"coincidence pilot: pool = {pilot_pool}; {len(cc_cnt)} cells with >= 2 stations")


def cell_times(cell: str, sta: str, cond_sel, daymask: np.ndarray):
    sub = POOL[(POOL.cell == cell) & (POOL.station == sta)]
    B = STA_BUF[sta]
    out = []
    for r in sub.itertuples():
        j = CODE_OF[sta][r.template]
        m = ((B["code"] == j) & (B["cc"] >= r.thr_use)
             & np.isin(STA_COND[sta][j][B["day"]], cond_sel) & daymask[B["day"]])
        if m.any():
            out.append(B["t"][m])
    return np.sort(np.concatenate(out)) if out else np.zeros(0)


def nn_lag(ta: np.ndarray, tb: np.ndarray):
    """Signed nearest-neighbour lag (tb - ta) for every ta, kept where |lag| <= LAG_MAX_S."""
    if ta.size == 0 or tb.size == 0:
        return np.zeros(0)
    i = np.clip(np.searchsorted(tb, ta), 1, tb.size - 1)
    cand = np.stack([tb[i - 1] - ta, tb[i] - ta])
    lag = cand[np.argmin(np.abs(cand), axis=0), np.arange(ta.size)]
    return lag[np.abs(lag) <= LAG_MAX_S]


def matched_frac(ta: np.ndarray, tb: np.ndarray, centre: float = 0.0, shift: float = 0.0):
    if ta.size == 0 or tb.size == 0:
        return np.nan, 0, int(ta.size)
    tbs = np.sort(tb + shift - centre)
    lo = np.searchsorted(tbs, ta - COINC_S, "left")
    hi = np.searchsorted(tbs, ta + COINC_S, "right")
    return float((hi > lo).sum()) / ta.size, int((hi > lo).sum()), int(ta.size)


COND_SEL = {"ets": (2,), "quiet": (1, 3)}
crows = []
for cell, cr in cc_cnt.head(N_CELLS_PILOT).iterrows():
    stas = sorted(POOL[POOL.cell == cell].station.unique())
    for i in range(len(stas)):
        for jj in range(i + 1, len(stas)):
            A, B = stas[i], stas[jj]
            sa = POOL[(POOL.cell == cell) & (POOL.station == A)]
            sb = POOL[(POOL.cell == cell) & (POOL.station == B)]
            both = STA_ACT[A] & STA_ACT[B]
            jA = CODE_OF[A][sa.template.iloc[0]]
            rec = dict(cell=cell, lat=cr.lat, lon=cr.lon, seg=seg_of_lat(cr.lat),
                       sta_a=A, sta_b=B, n_sta_cell=int(cr.n_sta), n_tpl_cell=int(cr.n_tpl),
                       n_tpl_a=len(sa), n_tpl_b=len(sb), n_days_both=int(both.sum()),
                       thr_a=float(sa.thr_use.min()), thr_b=float(sb.thr_use.min()))
            T = {}
            for cname, csel in COND_SEL.items():
                nd = int((both & np.isin(STA_COND[A][jA], csel)).sum())
                ta = tb = np.zeros(0)
                if nd >= MIN_BOTH_DAYS:
                    ta, tb = cell_times(cell, A, csel, both), cell_times(cell, B, csel, both)
                T[cname] = (ta, tb, nd)
            lg = nn_lag(*T["ets"][:2])
            if lg.size < 20:
                lg = nn_lag(*T["quiet"][:2])
            if lg.size:
                hh, ee = np.histogram(lg, bins=np.arange(-LAG_MAX_S, LAG_MAX_S + 1, 1.0))
                kk = int(np.argmax(hh))
                peak = float(0.5 * (ee[kk] + ee[kk + 1]))
                rec.update(lag_n=int(lg.size), lag_peak_s=peak, lag_peak_count=int(hh.max()),
                           lag_median_s=float(np.median(lg)),
                           lag_frac_within6_of_peak=float(np.mean(np.abs(lg - peak) <= COINC_S)),
                           lag_frac_within6_of_zero=float(np.mean(np.abs(lg) <= COINC_S)))
            else:
                peak = 0.0
                rec.update(lag_n=0, lag_peak_s=np.nan, lag_peak_count=0, lag_median_s=np.nan,
                           lag_frac_within6_of_peak=np.nan, lag_frac_within6_of_zero=np.nan)
            for cname in ("ets", "quiet"):
                ta, tb, nd = T[cname]
                fr0, np0, na = matched_frac(ta, tb, 0.0)
                frp, npp, _ = matched_frac(ta, tb, peak)
                rec[f"ndays_both_{cname}"] = nd
                rec[f"n_a_{cname}"] = na
                rec[f"n_b_{cname}"] = int(tb.size)
                rec[f"frac_{cname}"] = fr0
                rec[f"npair_{cname}"] = np0
                rec[f"frac_{cname}_atpeak"] = frp
                rec[f"npair_{cname}_atpeak"] = npp
                # the FRACTION's denominator changes between conditions, so the coincident pair
                # RATE per day is the statistic that can be compared across ETS and quiet
                rec[f"pairrate_{cname}"] = (np0 / nd) if nd > 0 else np.nan
                rec[f"pairrate_{cname}_atpeak"] = (npp / nd) if nd > 0 else np.nan
                rec[f"detrate_a_{cname}"] = (na / nd) if nd > 0 else np.nan
                rec[f"detrate_b_{cname}"] = (tb.size / nd) if nd > 0 else np.nan
                if cname == "ets":
                    frs, nps, _ = matched_frac(ta, tb, 0.0, SHIFT_S)
                    rec["frac_shift"] = frs
                    rec["npair_shift"] = nps
                    rec["pairrate_shift"] = (nps / nd) if nd > 0 else np.nan
                    rec["frac_poisson_exp"] = ((1 - np.exp(-tb.size / (nd * 86400.0) * 2 * COINC_S))
                                               if nd > 0 else np.nan)
            crows.append(rec)
CO = pd.DataFrame(crows)
CO.to_csv(OUT_RES / "d1_coincidence.csv", index=False)


# ================================================================ 10. station / segment summary
def summarise(sub: pd.DataFrame, key: str, val: str) -> dict:
    r = dict(group=key, value=val, n_templates=len(sub))
    if sub.empty:
        return r
    r.update(
        n_responsive=int(sub.responsive.sum()), n_unresponsive=int(sub.unresponsive.sum()),
        frac_responsive=float(sub.responsive.mean()),
        med_ratio=float(sub.ratio.median()), med_ratio_q10=float(sub.ratio_q10.median()),
        med_ratio_ym=float(sub.ratio_ym.median()),
        agg_ratio=float(sub.rate_ets.sum() / max(sub.rate_quiet.sum(), 1e-12)),
        med_rate_ets=float(sub.rate_ets.median()), med_rate_quiet=float(sub.rate_quiet.median()),
        med_excess_total=float(sub.excess_total.median()),
        frac_excess_positive=float((sub.excess_total > 0).mean()),
        frac_excess_positive_ym=float((sub.excess_total_ym > 0).mean()),
        med_median_cc_quiet=float(sub.median_cc_quiet.median()),
        med_median_cc_excess=float(sub.median_cc_excess.median()),
        med_dcc=float(sub.dcc.median()), med_dcc_nfm=float(sub.dcc_nfm.median()),
        frac_dcc_positive=float((sub.dcc > 0).mean()),
        med_prec_max=float(sub.prec_max.median()),
        med_prec_max_q10=float(sub.prec_max_q10.median()),
        med_prec_max_nfe=float(sub.prec_max_nfe.median()),
        med_prec_max_nfm=float(sub.prec_max_nfm.median()),
        frac_p80_reachable=float(np.isfinite(sub.thr_p80).mean()),
        frac_p80_reachable_nfm=float(np.isfinite(sub.thr_p80_nfm).mean()),
        med_thr_p80=float(sub.thr_p80.median()), med_recall_p80=float(sub.recall_p80.median()),
        med_thr_p80_nfm=float(sub.thr_p80_nfm.median()),
        med_recall_p80_nfm=float(sub.recall_p80_nfm.median()),
        frac_repairable=float(sub.repairable.mean()),
        frac_repairable_q10=float(sub.repairable_q10.mean()),
        frac_repairable_nfe=float(sub.repairable_nfe.mean()),
        frac_repairable_nfm=float(sub.repairable_nfm.mean()),
        frac_repairable_ym=float(sub.repairable_ym.mean()),
        med_thr_fa1=float(sub.thr_fa1perday.median()),
        med_thr_nf1_exp=float(sub.thr_nf1_exp.median()),
        med_thr_nf1_emp=float(sub.thr_nf1_emp.median()),
        med_thr_d2=float(sub.thr_d2.median()),
        med_tail_r2=float(sub.tail_r2.median()), med_nf_exp_r2=float(sub.nf_exp_r2.median()),
        med_floor_model_ratio_090=float(sub.floor_model_ratio_090.median()),
        med_quiet_over_floor_090=float(sub.quiet_over_floor_090.median()),
        med_quiet_over_expfloor_090=float(sub.quiet_over_expfloor_090.median()),
        frac_thr_d2_defined=float(sub.thr_d2_defined.mean()),
        med_cand95_rate_quiet=float(sub.cand95_rate_quiet.median()),
        med_cand95_rate_ets=float(sub.cand95_rate_ets.median()),
        frac_cand95_positive=float((sub.cand95_rate_quiet > 0).mean()),
        med_cand95_dispersion=float(sub.cand95_dispersion.median()),
        frac_cand95_bursty=float((sub.cand95_dispersion > 1.5).mean()),
        med_cand95_frac_days_zero=float(sub.cand95_frac_days_zero.median()),
        med_cand95_poisson_frac_zero=float(sub.cand95_poisson_frac_zero.median()),
        med_cand95_rho_30d=float(sub.cand95_rho_30d_tremor.median()),
        frac_cand95_rho30_pos=float((sub.cand95_rho_30d_tremor > 0).mean()),
        med_cand95_topbot_ratio=float(sub.cand95_topbot_ratio.median()),
        frac_cand95_topbot_gt1=float((sub.cand95_topbot_ratio > 1).mean()),
        med_snr_at_thr=float(sub.snr_at_thr.median()),
        med_cand_rate_quiet=float(sub.cand_rate_quiet.median()),
        med_cand_rate_ets=float(sub.cand_rate_ets.median()),
        frac_cand_positive=float((sub.cand_rate_quiet > 0).mean()),
        med_cand_dispersion=float(sub.cand_dispersion.median()),
        frac_cand_bursty=float((sub.cand_dispersion > 1.5).mean()),
        med_cand_frac_days_zero=float(sub.cand_frac_days_zero.median()),
        med_cand_poisson_frac_zero=float(sub.cand_poisson_frac_zero.median()),
        med_cand_rho_daily=float(sub.cand_rho_daily_tremor.median()),
        med_cand_rho_30d=float(sub.cand_rho_30d_tremor.median()),
        frac_cand_rho30_pos=float((sub.cand_rho_30d_tremor > 0).mean()),
        med_cand_topbot_ratio=float(sub.cand_topbot_ratio.median()),
        med_floor_trend=float(sub.floor_trend.median()),
        frac_floor_drift_p05=float((sub.floor_trend_p < 0.05).mean()),
        med_floor_year_spread=float(sub.floor_year_spread.median()),
        med_floor_cv=float(sub.floor_cv.median()),
        med_n_ets_days=float(sub.n_ets_days.median()),
        med_n_quiet_days=float(sub.n_quiet_days.median()),
        med_n_quiet10_days=float(sub.n_quiet10_days.median()))
    return r


SS = pd.DataFrame([summarise(TT, "all", "all")]
                  + [summarise(TT[TT.station == s], "station", s)
                     for s in sorted(TT.station.unique())]
                  + [summarise(TT[TT.segment == s], "segment", s)
                     for s in sorted(TT.segment.unique())]
                  + [summarise(TT[TT.rclass == c], "rclass", c)
                     for c in ("responsive", "mid", "unresponsive")])
SS = SS.merge(COV.rename(columns={"sta": "value"})[["value", "GB", "detections", "files",
                                                    "active_days"]], on="value", how="left")
SS.to_csv(OUT_RES / "d1_station_summary.csv", index=False)


# ================================================================ 11. figure
fig, ax = plt.subplots(5, 3, figsize=(17.5, 22.5))
fig.suptitle(
    "D1 - can an LFE catalogue be recovered from the stored densify cc values?\n"
    f"{len(COV)} stations, {COV.detections.sum() / 1e9:.2f} G detections, {len(TT)} templates."
    "   Floors: Q25/Q10 = quiet-day histograms (NOT pure - they contain inter-ETS LFEs);"
    "   NF-exp / NF-emp = noise-only models",
    fontsize=11, y=0.998, va="top")
CQ, CE, CX, CF = "#2166ac", "#b2182b", "#1a9850", "#8c510a"


def hist_panel(a, cls, title, fields):
    sub = SH[SH.rclass == cls]
    if sub.empty:
        a.set_title(f"{title} (none)", fontsize=9.5)
        return
    n = int(sub.n_templates.iloc[0])
    for cname, col, lw, ls, lab in fields:
        s = sub[sub.condition == cname].sort_values("bin_ctr")
        if s.empty:
            continue
        y = s.mean_rate_per_day.values
        a.step(s.bin_ctr.values, np.clip(y, 1e-7, None), where="mid", color=col, lw=lw, ls=ls,
               label=lab)
        if (y < 0).any():
            a.step(s.bin_ctr.values, np.clip(-y, 1e-7, None), where="mid", color=col,
                   lw=max(lw - 0.6, 0.7), ls=":")
    a.set_yscale("log")
    a.set_ylim(1e-5, None)
    a.set_xlabel("matched-filter cc")
    a.set_ylabel("det / template / day / 0.005 bin")
    a.set_title(f"{title}  (n={n})", fontsize=9.5)
    a.legend(fontsize=6.5)
    a.grid(alpha=0.3)


hist_panel(ax[0, 0], "responsive", f"A  responsive (ETS/Q25 >= {RESP_HI})",
           [("ets", CE, 1.8, "-", "ETS days"), ("quiet", CQ, 1.8, "-", "QUIET Q25"),
            ("quiet10", "#67a9cf", 1.2, "--", "QUIET Q10"),
            ("floor_nf_emp", CF, 1.6, "-", "NF-emp floor (quiet-scaled)"),
            ("floor_nf_emp_ets", "#bf812d", 1.6, "--", "NF-emp floor (ETS-scaled)"),
            ("floor_nf_exp", "#d8b365", 1.2, "--", "NF-exp floor (quiet-scaled)"),
            ("excess", CX, 1.4, "-", "ETS - Q25")])
hist_panel(ax[0, 1], "unresponsive", f"B  unresponsive (ETS/Q25 < {RESP_LO})",
           [("ets", CE, 1.8, "-", "ETS days"), ("quiet", CQ, 1.8, "-", "QUIET Q25"),
            ("floor_nf_emp", CF, 1.6, "-", "NF-emp noise floor"),
            ("excess", CX, 1.4, "-", "ETS - Q25")])
hist_panel(ax[0, 2], "responsive", "C  responsive: what sits ABOVE the noise-only floor",
           [("excess", CX, 1.6, "-", "ETS - Q25 (ETS excess)"),
            ("excess_nf_emp", CE, 1.8, "-", "ETS - NF-emp (ETS-scaled)"),
            ("quiet_minus_floor", CQ, 1.8, "-", "QUIET - NF-emp = inter-ETS candidates"),
            ("excess_ym", "#999999", 1.2, "--", "ETS - Q25, year-matched")])


def pr_panel(a, pf, rf, title):
    for cls, col in (("responsive", CE), ("mid", "#777777"), ("unresponsive", CQ)):
        sub = TT[TT.rclass == cls]
        if sub.empty:
            continue
        keys = [(r.station, r.template) for r in sub.itertuples()]
        P = np.stack([CURVES[k][pf] for k in keys])
        R = np.stack([CURVES[k][rf] for k in keys])
        a.plot(THR, np.nanmedian(P, axis=0), color=col, lw=2,
               label=f"{cls} precision (n={len(sub)})")
        a.plot(THR, np.nanmedian(R, axis=0), color=col, lw=1.2, ls="--", label=f"{cls} recall")
    a.axhline(PREC_TARGET, color="k", ls=":", lw=1)
    a.axhline(PREC_SOFT, color="k", ls="-.", lw=0.8)
    a.set_xlabel("cc threshold")
    a.set_ylabel("median precision (solid) / recall (dashed)")
    a.set_ylim(-0.05, 1.05)
    a.set_title(title, fontsize=9.5)
    a.legend(fontsize=6.5)
    a.grid(alpha=0.3)


pr_panel(ax[1, 0], "prec_q25", "rec_q25", "D  precision/recall vs Q25 quiet floor (spec)")
pr_panel(ax[1, 1], "prec_nfm", "rec_nfm", "E  precision/recall vs NF-emp NOISE-ONLY floor")
pr_panel(ax[1, 2], "prec_nfe", "rec_nfe", "F  precision/recall vs NF-exp NOISE-ONLY floor")

a = ax[2, 0]
for col, lab, v in ((CQ, "Q25 floor", TT.prec_max), ("#67a9cf", "Q10 floor", TT.prec_max_q10),
                    (CF, "NF-emp floor", TT.prec_max_nfm),
                    ("#d8b365", "NF-exp floor", TT.prec_max_nfe)):
    vv = v.dropna()
    if vv.size:
        a.hist(vv.clip(-0.5, 1), bins=np.linspace(-0.5, 1, 61), histtype="step", lw=1.8,
               color=col, label=f"{lab} (median {vv.median():.2f})")
a.axvline(PREC_TARGET, color="k", ls=":", lw=1)
a.axvline(PREC_SOFT, color="k", ls="-.", lw=0.8)
a.set_xlabel("max attainable precision over cc 0.80-0.98")
a.set_ylabel("templates")
a.set_title("G  best-case precision, by floor model", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[2, 1]
for col, lab, v in ((CQ, "Q25", TT.thr_p80), (CF, "NF-emp", TT.thr_p80_nfm),
                    ("#d8b365", "NF-exp", TT.thr_p80_nfe)):
    vv = v.dropna()
    if vv.size:
        a.hist(vv, bins=np.arange(0.7975, 0.9851, 0.005), histtype="step", lw=1.8, color=col,
               label=f"{lab} ({100 * np.isfinite(v).mean():.0f} % reach p80)")
a.set_xlabel(f"threshold reaching precision >= {PREC_TARGET}")
a.set_ylabel("templates")
a.set_title("H  thr_p80 by floor model", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[2, 2]
for col, lab, v in ((CQ, "Q25", TT.recall_p80), (CF, "NF-emp", TT.recall_p80_nfm),
                    ("#d8b365", "NF-exp", TT.recall_p80_nfe)):
    vv = v.dropna()
    if vv.size:
        a.hist(vv, bins=np.linspace(0, 1, 26), histtype="step", lw=1.8, color=col,
               label=f"{lab} (median {vv.median():.2f})")
a.axvline(0.5, color="k", ls=":", lw=1)
a.set_xlabel("recall at thr_p80")
a.set_ylabel("templates")
a.set_title(f"I  recall at 80 % precision; repairable: Q25 {100 * TT.repairable.mean():.1f} %, "
            f"NF-emp {100 * TT.repairable_nfm.mean():.1f} %", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[3, 0]
s = TT.quiet_over_floor_090.replace([np.inf, -np.inf], np.nan).dropna()
if s.size:
    a.hist(np.log10(s.clip(1e-3, 1e4)), bins=60, color=CQ, alpha=0.8,
           label="observed QUIET / NF-emp floor")
s2 = TT.floor_model_ratio_090.replace([np.inf, -np.inf], np.nan).dropna()
if s2.size:
    a.hist(np.log10(s2.clip(1e-3, 1e4)), bins=60, histtype="step", lw=1.8, color=CF,
           label="NF-emp / NF-exp floor")
a.axvline(0, color="k", ls=":", lw=1)
a.set_xlabel("log10 ratio of the rate above cc 0.90")
a.set_ylabel("templates")
a.set_title(f"J  floor models vs the observed quiet tail\n(quiet/floor median "
            f"{TT.quiet_over_floor_090.median():.2f}x, emp/exp median "
            f"{TT.floor_model_ratio_090.median():.2f}x)", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[3, 1]
v = TT.cand95_dispersion.replace([np.inf, -np.inf], np.nan).dropna()
if v.size:
    a.hist(np.log10(v.clip(1e-2, 1e3)), bins=60, color=CX, alpha=0.8)
a.axvline(0, color="k", lw=1.2, label="Poisson (var/mean = 1)")
a.set_xlabel("log10 dispersion index, daily quiet-day count above cc 0.95")
a.set_ylabel("templates")
a.set_title(f"K  are inter-ETS candidates BURSTY?  median var/mean {v.median():.2f} "
            f"({100 * (TT.cand95_dispersion > 1.5).mean():.0f} % > 1.5)", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[3, 2]
for col, lab, v in ((CX, "vs local 30-d tremor", TT.cand95_rho_30d_tremor),
                    ("#984ea3", "vs local daily tremor", TT.cand95_rho_daily_tremor)):
    vv = v.dropna()
    if vv.size:
        a.hist(vv, bins=np.linspace(-0.6, 0.6, 49), histtype="step", lw=1.8, color=col,
               label=f"{lab} (median {vv.median():+.3f})")
a.axvline(0, color="k", ls=":", lw=1)
a.set_xlabel("Spearman(daily count above cc 0.95, local tremor), QUIET days")
a.set_ylabel("templates")
a.set_title(f"L  are inter-ETS candidates TREMOR-LINKED?  "
            f"{100 * (TT.cand95_rho_30d_tremor > 0).mean():.0f} % positive", fontsize=9.5)
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[4, 0]
s = TT.dropna(subset=["median_cc_excess"])
a.scatter(s.median_cc_quiet, s.median_cc_excess, s=8, alpha=0.45,
          c=np.where(s.responsive, CE, np.where(s.unresponsive, CQ, "#999999")))
a.plot([0.80, 1.0], [0.80, 1.0], "k:", lw=1)
a.set_xlabel("median cc of QUIET (Q25) detections")
a.set_ylabel("median cc of the ETS excess")
a.set_title(f"M  where the excess sits  (median shift {np.nanmedian(TT.dcc):+.4f}, "
            f"{100 * (TT.dcc > 0).mean():.0f} % > 0)", fontsize=9.5)
a.grid(alpha=0.3)

a = ax[4, 1]
if not CO.empty and np.isfinite(CO.pairrate_ets).any():
    lab = ["ETS\nlag 0", "QUIET\nlag 0", "ETS\nat peak lag", "QUIET\nat peak lag",
           "ETS shifted\n+600 s"]
    vals = [np.nanmean(CO.pairrate_ets), np.nanmean(CO.pairrate_quiet),
            np.nanmean(CO.pairrate_ets_atpeak), np.nanmean(CO.pairrate_quiet_atpeak),
            np.nanmean(CO.pairrate_shift)]
    a.bar(lab, vals, color=[CE, CQ, "#e08214", "#67a9cf", "#999999"])
    for i, v in enumerate(vals):
        if np.isfinite(v):
            a.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    a.set_ylabel(f"coincident pairs per day (+-{COINC_S:.0f} s)")
    a.set_title(f"N  network coincidence RATE, {CO.cell.nunique()} cells / {len(CO)} pairs; "
                f"peak lag median {np.nanmedian(CO.lag_peak_s):+.1f} s", fontsize=9.5)
    a.tick_params(axis="x", labelsize=7)
else:
    a.set_title("N  network coincidence (no eligible pair)", fontsize=9.5)
a.grid(alpha=0.3, axis="y")

a = ax[4, 2]
g = SS[SS.group == "station"].sort_values("value")
xx = np.arange(len(g))
a.bar(xx - 0.3, g.frac_responsive, width=0.2, color=CE, label="frac responsive")
a.bar(xx - 0.1, g.frac_repairable, width=0.2, color=CQ, label="repairable (Q25 floor)")
a.bar(xx + 0.1, g.frac_repairable_nfm, width=0.2, color=CF, label="repairable (NF-emp floor)")
a.bar(xx + 0.3, g.med_prec_max_nfm, width=0.2, color="#e08214",
      label="median max precision (NF-emp)")
a.set_xticks(xx)
a.set_xticklabels(g.value, rotation=60, fontsize=8)
a.set_ylabel("fraction / precision")
a.set_title("O  per station", fontsize=9.5)
a.legend(fontsize=6.5)
a.grid(alpha=0.3, axis="y")

fig.tight_layout(rect=(0, 0, 1, 0.975))
fig.savefig(OUT_FIG / "d1_cc_separability.png", dpi=118)
plt.close(fig)


# ================================================================ 12. printed summary
def pr(*a_):
    print(*a_, flush=True)


def q(v, f="{:.3f}"):
    v = pd.Series(v).replace([np.inf, -np.inf], np.nan).dropna()
    if v.empty:
        return "n/a"
    return (f + " [" + f + ", " + f + "]").format(v.median(), np.percentile(v, 25),
                                                  np.percentile(v, 75))


pr("\n" + "=" * 108)
pr("D1 SEPARABILITY - SUMMARY")
pr("=" * 108)
pr(f"streamed {COV.GB.sum():.1f} GB, {int(COV.files.sum())} station-year files, "
   f"{COV.detections.sum():,} detections, {COV.templates.sum()} templates, "
   f"{COV.wall_s.sum() / 60:.1f} min of streaming")
pr(COV.to_string(index=False))
pr(f"\ntemplates analysed (>= {MIN_COND_DAYS} ETS and >= {MIN_COND_DAYS} Q25 recording days): "
   f"{len(TT)} of {COV.templates.sum()}")
pr(f"ETS days/template {q(TT.n_ets_days, '{:.0f}')}; Q25 days {q(TT.n_quiet_days, '{:.0f}')}; "
   f"Q10 days {q(TT.n_quiet10_days, '{:.0f}')}")
pr(f"rate_quiet(Q25) {q(TT.rate_quiet, '{:.1f}')} det/template/day; "
   f"rate_quiet(Q10) {q(TT.rate_quiet10, '{:.1f}')}; rate_ets {q(TT.rate_ets, '{:.1f}')}")
pr(f"responsiveness ETS/Q25 {q(TT.ratio)} (responsive {TT.responsive.sum()} = "
   f"{100 * TT.responsive.mean():.1f} %, unresponsive {TT.unresponsive.sum()} = "
   f"{100 * TT.unresponsive.mean():.1f} %)")
pr(f"responsiveness ETS/Q10 {q(TT.ratio_q10)};  year-matched ETS/Q25 {q(TT.ratio_ym)}")
pr(f"aggregate (detection-weighted) ETS/Q25 ratio "
   f"{TT.rate_ets.sum() / TT.rate_quiet.sum():.3f}  [S14 reported x1.24 on the same files]")
pr(f"excess_total ETS-Q25 {q(TT.excess_total, '{:+.2f}')} det/day, positive in "
   f"{100 * (TT.excess_total > 0).mean():.1f} %;  year-matched "
   f"{q(TT.excess_total_ym, '{:+.2f}')}, positive in "
   f"{100 * (TT.excess_total_ym > 0).mean():.1f} %")
w = TT.excess_total_ym.dropna()
if len(w) > 20:
    pr(f"  Wilcoxon(excess_total_ym vs 0) p = {wilcoxon(w).pvalue:.3g}")

pr("\n--- (a) WHERE the excess sits in cc ---")
for cls in ("responsive", "mid", "unresponsive", "all"):
    s = TT if cls == "all" else TT[TT.rclass == cls]
    if s.empty:
        continue
    pr(f"  {cls:12s} n={len(s):4d}  cc median: QUIET {s.median_cc_quiet.median():.4f}  "
       f"ETS {s.median_cc_ets.median():.4f}  EXCESS(Q25) {s.median_cc_excess.median():.4f} "
       f"shift {s.dcc.median():+.4f} (>0 in {100 * (s.dcc > 0).mean():.0f} %)  |  "
       f"EXCESS(NF-emp) {s.median_cc_excess_nfm.median():.4f} shift {s.dcc_nfm.median():+.4f}")
pr(f"  detections above cc 0.90: QUIET {100 * TT.q_frac_gt_090.median():.2f} %, "
   f"ETS {100 * TT.e_frac_gt_090.median():.2f} %;  above 0.95: QUIET "
   f"{100 * TT.q_frac_gt_095.median():.3f} %, ETS {100 * TT.e_frac_gt_095.median():.3f} %")
pr(f"  quiet-tail exp fit 0.80-0.90: slope {q(TT.tail_slope, '{:+.1f}')} /unit cc, "
   f"R2 {q(TT.tail_r2)} ({100 * (TT.tail_r2 > 0.9).mean():.0f} % R2>0.9)")
pr(f"  NF-exp fit 0.80-0.85 only:    slope {q(TT.nf_exp_slope, '{:+.1f}')} /unit cc, "
   f"R2 {q(TT.nf_exp_r2)} ({100 * (TT.nf_exp_r2 > 0.9).mean():.0f} % R2>0.9)")

pr("\n--- (a2) THE NOISE-ONLY FLOOR MODELS, and how much of the quiet tail is NOT noise ---")
pr(f"  rate above cc 0.90, per template per day:  NF-exp floor {q(TT.floor_exp_above_090)}   "
   f"NF-emp floor {q(TT.floor_emp_above_090)}")
pr(f"  observed QUIET above 0.90 {q(TT.obs_quiet_above_090)};  observed ETS above 0.90 "
   f"{q(TT.obs_ets_above_090)}")
pr(f"  NF-emp / NF-exp floor ratio {q(TT.floor_model_ratio_090, '{:.2f}')}x  "
   f"-> DISagreement between the two noise-only models")
pr("  IS THE HEAVY HIGH-cc TAIL NOISE OR LFEs?  observed QUIET rate above 0.90 divided by the "
   "EXPONENTIAL extrapolation, per responsiveness class:")
for _cls in ("responsive", "mid", "unresponsive"):
    _s = TT[TT.rclass == _cls]
    if _s.empty:
        continue
    pr(f"    {_cls:12s} {q(_s.quiet_over_expfloor_090, '{:.1f}')}x")
pr("    the UNRESPONSIVE templates have no ETS response and therefore no LFE population, so "
   "their factor is a PURE-NOISE measurement:")
pr("    if it is comparable to the responsive one, the heavy tail is FALSE ALARMS and NF-exp "
   "is invalid as a floor -> use NF-emp.")
pr(f"  observed QUIET / NF-emp floor {q(TT.quiet_over_floor_090, '{:.2f}')}x  "
   f"-> the quiet tail exceeds the noise floor by this factor (candidate inter-ETS LFEs)")
pr(f"  target-FA thresholds (median): Q25-fit 1/day {TT.thr_fa1perday.median():.3f}; "
   f"NF-exp 1/day {TT.thr_nf1_exp.median():.3f} / 0.1/day {TT.thr_nf01_exp.median():.3f}; "
   f"NF-emp 1/day {TT.thr_nf1_emp.median():.3f} / 0.1/day {TT.thr_nf01_emp.median():.3f}")
pr(f"  ** RECOMMENDED D2 per-template threshold (NOISE-ONLY, NF-emp at 1 FA/day): "
   f"{q(TT.thr_d2)}; it exists for "
   f"{100 * TT.thr_d2_defined.mean():.1f} % of templates (for the rest the modelled noise floor "
   f"never falls below 1 false alarm/day anywhere up to cc 0.995) **")
pr(f"  (NF-emp 0.1 FA/day threshold exists for "
   f"{100 * np.isfinite(TT.thr_nf01_emp).mean():.1f} % of templates)")

pr("\n--- (b) REPAIRABILITY by thresholding alone, BY FLOOR MODEL ---")
pr(f"  {'class':12s} {'n':>5s} | {'floor':7s} {'maxPrec':>8s} {'p>=0.5':>7s} {'p>=0.8':>7s} "
   f"{'thr_p80':>8s} {'recall':>7s} {'REPAIR':>7s}")
for cls in ("responsive", "mid", "unresponsive", "all"):
    s = TT if cls == "all" else TT[TT.rclass == cls]
    if s.empty:
        continue
    for lab, pm, t50, tp, rp, rep in (
            ("Q25", s.prec_max, s.thr_p50, s.thr_p80, s.recall_p80, s.repairable),
            ("Q10", s.prec_max_q10, s.thr_p80_q10, s.thr_p80_q10, s.recall_p80_q10,
             s.repairable_q10),
            ("NF-emp", s.prec_max_nfm, s.thr_p80_nfm, s.thr_p80_nfm, s.recall_p80_nfm,
             s.repairable_nfm),
            ("NF-exp", s.prec_max_nfe, s.thr_p80_nfe, s.thr_p80_nfe, s.recall_p80_nfe,
             s.repairable_nfe)):
        pr(f"  {cls:12s} {len(s):5d} | {lab:7s} {pm.median():8.3f} "
           f"{100 * (pm >= PREC_SOFT).mean():6.1f}% {100 * np.isfinite(tp).mean():6.1f}% "
           f"{tp.median():8.3f} {rp.median():7.3f} {100 * rep.mean():6.1f}%")
pr(f"  (year-matched, Q25 floor) repairable {100 * TT.repairable_ym.mean():.1f} %")
pr(f"  observed ETS rate / extrapolated FA rate at the Q25 1/day threshold "
   f"{q(TT.snr_at_thr, '{:.2f}')}; at 0.1/day {q(TT.snr_at_thr0p1, '{:.2f}')}")

pr("\n--- (b2) INTER-ETS LFE CANDIDATES (quiet-day detections above the NOISE-ONLY floor) ---")
pr(f"  candidate rate on QUIET days {q(TT.cand_rate_quiet, '{:.3f}')} /template/day "
   f"(positive for {100 * (TT.cand_rate_quiet > 0).mean():.1f} % of templates)")
pr(f"  candidate rate on ETS days   {q(TT.cand_rate_ets, '{:.3f}')} /template/day")
pr(f"  NF floor rate at that threshold {q(TT.floor_rate_at_thr, '{:.3f}')} /day; "
   f"observed mean daily count above it: quiet {q(TT.cand_obs_mean_quiet, '{:.3f}')}, "
   f"ETS {q(TT.cand_obs_mean_ets, '{:.3f}')}")
pr(f"  summed candidate rate over templates: QUIET {TT.cand_rate_quiet.sum():+.1f}, "
   f"ETS {TT.cand_rate_ets.sum():+.1f} det/day; ETS > QUIET for "
   f"{100 * (TT.cand_rate_ets > TT.cand_rate_quiet).mean():.1f} % of templates")
pr(f"  BURSTINESS dispersion var/mean of the daily quiet-day count above threshold "
   f"{q(TT.cand_dispersion, '{:.2f}')} "
   f"({100 * (TT.cand_dispersion > 1.5).mean():.1f} % > 1.5, "
   f"{100 * (TT.cand_dispersion > 3).mean():.1f} % > 3)")
pr(f"    frac quiet days with ZERO above threshold {q(TT.cand_frac_days_zero, '{:.3f}')} "
   f"vs Poisson {q(TT.cand_poisson_frac_zero, '{:.3f}')}; frac with >= 3 "
   f"{q(TT.cand_frac_days_ge3, '{:.3f}')}")
pr(f"  TREMOR LINK Spearman(daily count, local 30-d tremor) on quiet days "
   f"{q(TT.cand_rho_30d_tremor, '{:+.3f}')} "
   f"({100 * (TT.cand_rho_30d_tremor > 0).mean():.1f} % positive); vs local DAILY tremor "
   f"{q(TT.cand_rho_daily_tremor, '{:+.3f}')} "
   f"({100 * (TT.cand_rho_daily_tremor > 0).mean():.1f} % positive)")
pr(f"    candidate rate top vs bottom half of quiet-day tremor: ratio "
   f"{q(TT.cand_topbot_ratio, '{:.3f}')} "
   f"({100 * (TT.cand_topbot_ratio > 1).mean():.1f} % > 1)")
pr("\n  SAME DIAGNOSTICS AT A FIXED cc >= 0.95 (defined for every template, so not conditioned "
   "on a threshold existing):")
pr(f"    candidate rate quiet {q(TT.cand95_rate_quiet, '{:.4f}')} /template/day "
   f"(positive for {100 * (TT.cand95_rate_quiet > 0).mean():.1f} %); ETS "
   f"{q(TT.cand95_rate_ets, '{:.4f}')};  NF-emp floor {q(TT.cand95_floor, '{:.4f}')}")
pr(f"    summed candidate rate at 0.95: QUIET {TT.cand95_rate_quiet.sum():+.2f}, "
   f"ETS {TT.cand95_rate_ets.sum():+.2f} det/day; ETS > QUIET for "
   f"{100 * (TT.cand95_rate_ets > TT.cand95_rate_quiet).mean():.1f} % of templates")
pr(f"    BURSTINESS dispersion {q(TT.cand95_dispersion, '{:.2f}')} "
   f"({100 * (TT.cand95_dispersion > 1.5).mean():.1f} % > 1.5); frac quiet days with ZERO "
   f"{q(TT.cand95_frac_days_zero, '{:.3f}')} vs Poisson "
   f"{q(TT.cand95_poisson_frac_zero, '{:.3f}')}; frac >= 3 "
   f"{q(TT.cand95_frac_days_ge3, '{:.3f}')}")
pr(f"    TREMOR LINK rho(daily count, local 30-d tremor) "
   f"{q(TT.cand95_rho_30d_tremor, '{:+.3f}')} "
   f"({100 * (TT.cand95_rho_30d_tremor > 0).mean():.1f} % positive); top/bottom-half tremor "
   f"rate ratio {q(TT.cand95_topbot_ratio, '{:.3f}')} "
   f"({100 * (TT.cand95_topbot_ratio > 1).mean():.1f} % > 1)")

pr("\n--- (c) NETWORK COINCIDENCE pilot (with lag scan) ---")
pr(f"  pool = {pilot_pool}; {len(cc_cnt)} cells with >= 2 stations, "
   f"{min(len(cc_cnt), N_CELLS_PILOT)} used")
if not CO.empty:
    show = ["cell", "sta_a", "sta_b", "n_days_both", "thr_a", "thr_b", "lag_n", "lag_peak_s",
            "lag_frac_within6_of_zero", "lag_frac_within6_of_peak", "ndays_both_ets",
            "detrate_a_ets", "detrate_b_ets", "pairrate_ets", "pairrate_ets_atpeak",
            "ndays_both_quiet", "pairrate_quiet", "pairrate_quiet_atpeak", "pairrate_shift"]
    pr(CO[show].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    pr(f"  peak nearest-neighbour lag: median {np.nanmedian(CO.lag_peak_s):+.1f} s; |peak| > "
       f"{COINC_S:.0f} s for {100 * (CO.lag_peak_s.abs() > COINC_S).mean():.0f} % of pairs "
       f"-> a +-6 s window centred on ZERO misses those")
    pr(f"  mean coincident PAIR RATE /day:  ETS {np.nanmean(CO.pairrate_ets):.3f}  "
       f"QUIET {np.nanmean(CO.pairrate_quiet):.3f}  ETS shifted +600 s "
       f"{np.nanmean(CO.pairrate_shift):.3f}")
    pr(f"  at the peak lag:                 ETS {np.nanmean(CO.pairrate_ets_atpeak):.3f}  "
       f"QUIET {np.nanmean(CO.pairrate_quiet_atpeak):.3f}")
    pr(f"  mean coincidence FRACTION:       ETS {np.nanmean(CO.frac_ets):.4f}  QUIET "
       f"{np.nanmean(CO.frac_quiet):.4f}  shifted {np.nanmean(CO.frac_shift):.4f}  Poisson "
       f"{np.nanmean(CO.frac_poisson_exp):.4f}")
    pr(f"  ETS/shifted pair rate = "
       f"{np.nanmean(CO.pairrate_ets) / max(np.nanmean(CO.pairrate_shift), 1e-12):.2f}; "
       f"ETS/QUIET pair rate = "
       f"{np.nanmean(CO.pairrate_ets) / max(np.nanmean(CO.pairrate_quiet), 1e-12):.2f}")
else:
    pr("  no source cell had eligible templates at >= 2 stations")

pr("\n--- (d) FLOOR STATIONARITY ---")
pr(f"  Spearman(year, quiet-day rate) {q(TT.floor_trend, '{:+.3f}')}; |rho|>0.5 in "
   f"{100 * (TT.floor_trend.abs() > 0.5).mean():.1f} %; p<0.05 in "
   f"{100 * (TT.floor_trend_p < 0.05).mean():.1f} % of "
   f"{int(np.isfinite(TT.floor_trend).sum())} templates with >= 4 quiet years")
pr(f"  max/min quiet-year rate {q(TT.floor_year_spread, '{:.2f}')}x; CV {q(TT.floor_cv, '{:.2f}')}")

cols = ["value", "n_templates", "frac_responsive", "med_ratio", "agg_ratio", "med_ratio_q10",
        "med_rate_quiet", "med_dcc", "med_prec_max", "med_prec_max_nfm", "med_thr_d2",
        "frac_thr_d2_defined", "frac_repairable", "frac_repairable_nfm",
        "med_quiet_over_floor_090", "med_quiet_over_expfloor_090", "med_cand95_rate_quiet",
        "med_cand95_dispersion", "frac_cand95_bursty", "med_cand95_rho_30d",
        "frac_cand95_rho30_pos", "med_cand95_topbot_ratio", "med_tail_r2", "med_nf_exp_r2",
        "med_floor_trend", "med_floor_year_spread"]
for grp, ttl in (("station", "per station"), ("segment", "per segment (source-cell latitude)"),
                 ("rclass", "per responsiveness class")):
    pr(f"\n--- {ttl} ---")
    pr(SS[SS.group == grp][cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

pr("\nfiles written:")
for f in ("results/d1_template_table.csv", "results/d1_stacked_hist.csv",
          "results/d1_floor_shapes.csv", "results/d1_coincidence.csv",
          "results/d1_station_summary.csv", "figures/d1_cc_separability.png"):
    p = HERE / f
    pr(f"  {p}  ({p.stat().st_size / 1e3:.1f} kB)")
pr(f"\ntotal wall {(time.time() - T0) / 60:.1f} min")
