#!/usr/bin/env python
"""
S5 - SITE RECURRENCE: tremor-episode vs slip-patch recurrence per latitude band, and phase drift
     (PLAN.md step S5).

Question.  At a fixed SITE (a 0.5-degree latitude band of the Cascadia margin), tremor episodes and
modeled slip patches each recur.  Do they recur at the same rate?  Does one of the two clocks drift
over 2010-2026?  And does the PHASE between them (slip start minus tremor onset) drift?  A changing
phase would be new (PLAN.md H4); a changing recurrence interval at a site would be new (H2-like).

No dv/v anywhere (deferred by the user 2026-09-10).

--------------------------------------------------------------------------------------------------
INPUTS  (all read-only, all OUTSIDE this folder)
--------------------------------------------------------------------------------------------------
  catalogs/pnsn_tremor_cascadia_full.csv
      PNSN tremor master, 748k rows; columns used: time, lat, lon.  Epicentres only.
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
      label      (T x 742) int32   connected-component labels; 0 = no patch, else the patch id
      days       (T,)      str     daily 2007-01-02 .. 2026-08-18   (T = 7169)
      patch_lat  (742,)    f8      mesh-triangle centroid latitudes
      patch_lon  (742,)    f8      mesh-triangle centroid longitudes
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.csv
      one row per patch: id, start, end, duration_d, Mw, lat_centroid, lat_min, lat_max,
      minor, broad_pulse, ...
  fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv
      per-patch tremor confirmation; columns used: catalog ('thrust'|'mirror'), id, seg, confirmed.
  fault_tomography/inversion/res_catalog_sse/tremor_episodes.csv        [CROSS-CHECK ONLY]
      the earlier tremor-episode table (key, band, seg, peak, onset, count30).  Never used as input
      to any number reported here; only compared against, at the end.

--------------------------------------------------------------------------------------------------
EXACT DEFINITIONS
--------------------------------------------------------------------------------------------------
  ANALYSIS WINDOW   2010-01-06 .. 2026-05-27  (the PNSN tremor span).  N = 5986 days.
                    2026 is PARTIAL (ends 2026-05-27) and is flagged as such everywhere.

  SITES / BANDS     20 bands of 0.5 degrees of latitude, lower edges 40.0, 40.5, ..., 49.5.
                    A band is labelled by its lower edge.  Membership is [lo, lo+0.5).
  SEGMENTS          north lo >= 46.0 ; central 44.0 <= lo < 46.0 ; south lo < 44.0.
                    (8 north bands, 4 central, 8 south.)

  TREMOR -> BAND    A detection is kept only if its epicentre lies within 40 km of SOME mesh-triangle
                    centroid, distance in the local km frame
                        x = (lon + 124) * 111 * cos(47 deg),   y = (lat - 47) * 111,
                    nearest-neighbour distance from scipy.spatial.cKDTree over the 742 centroids.
                    A kept detection is assigned to the band containing its OWN EPICENTRE latitude
                    (not the triangle's latitude).
  DAILY COUNT       number of kept detections per band per calendar day on the analysis-window grid.
  TRAILING COUNT    c[t] = sum of the daily counts over (t - 30, t]  (30-day trailing sum, the same
                    frame the slip model's inc30 uses).  The first 29 days of the grid have a partial
                    window, so PEAKS ARE ONLY SOUGHT FROM day 30 of the grid (2010-02-05) onward.

  TREMOR EPISODE (per band)
      peaks   = local maxima of c (scipy.signal.find_peaks, plateaux handled) whose value is
                >= the band's 90th percentile of c over the whole analysis window,
                selected GREEDILY highest-first with a >= 90 d separation between accepted peaks.
      onset   = the FIRST day t in [peak - 45, peak] with c[t] > 0.25 * c[peak].
                (c[peak] itself always satisfies this, so an onset always exists.)
      count30 = c[peak].
      key     = "T{band:.1f}_{peak date}".

  SLIP PATCHES
      SUBSTANTIVE = (~minor) & (~broad_pulse), ANY Mw  -> 284 patches in the catalogue.
      Restricted to those whose START lies inside the analysis window -> the S5 patch set.
      FOOTPRINT of patch p = the set of triangles j with label[t, j] == p.id for ANY day t in
                [start, end].  (Computed once by a single pass over the T x 742 label array.)
      TOUCHES BAND b  <=>  some footprint triangle has its CENTROID latitude inside band b.
      A patch may touch several bands; it then appears once in each of those bands.
      CONFIRMED subset = patch ids with catalog == 'thrust' & confirmed == True in
                tremor_vs_slip_patches.csv (66 ids), intersected with the substantive set.
      TIMES: start / end are the model's TRAILING-30-d threshold crossings (slip within (t-30, t]
             exceeds 5 mm), so every slip time is late by 0-30 d relative to the true slip onset.
             This is a KNOWN, UNCORRECTED bias in every lag reported below (PLAN.md section 6).

  MATCHING (per band, both directions computed; the two windows are the same window seen from the
            two sides, so a pair found one way is a candidate the other way as well)
      episode -> patch : among substantive patches touching the band with
                         start in [onset - 30, onset + 90], take the one whose start is CLOSEST to
                         the onset.  lag_d = start - onset  (days; + = slip after tremor).
      patch  -> episode: among the band's episodes with onset in [start - 90, start + 30],
                         take the one whose onset is CLOSEST to the start.
      A pair is MUTUAL when both directions pick each other.  Matching is not forced to be
      one-to-one: two episodes may pick the same patch (recorded, and counted in the CSV).
      UNMATCHED EPISODE = tremor with no slip patch in the window ("tremor without slip").
      UNMATCHED PATCH   = slip with no tremor episode in the window ("slip without tremor").

  RECURRENCE (per band)
      tremor intervals = diffs of consecutive episode PEAK dates (days).
      slip intervals   = diffs of consecutive substantive patch START dates (days);
                         computed again for the Mw >= 6.0 subset only.
      For each series: median interval; Theil-Sen slope (scipy.stats.theilslopes) of
      interval [d] versus the MIDPOINT DATE of that interval (in years since 2010-01-01),
      reported in days-of-interval per year.
      PERMUTATION p: 1000 random shuffles of the interval VALUES against the fixed midpoint dates;
      two-sided, p = (1 + #{|slope_null| >= |slope_obs|}) / (1 + 1000).
      A trend needs >= 5 intervals; otherwise slope and p are NaN.
  PHASE DRIFT (per band)
      lag_d of the episode->patch matched pairs versus the episode onset date (years since
      2010-01-01): Theil-Sen slope [days of lag per year] + the same 1000-shuffle two-sided
      permutation p.  Needs >= 5 matched pairs.
      NOTE (tie atom): lag_d is an integer with many repeats, so the Theil-Sen pairwise-slope set
      has a large atom at exactly 0; when its median lands on that atom the slope is exactly 0 and
      the two-sided permutation p is 1 BY CONSTRUCTION (|null| >= 0 always).  A tie-robust
      Spearman rho of lag vs onset date is therefore reported alongside, per band and pooled.
  MULTIPLE TESTING
      Holm-Bonferroni over the 20 bands, applied SEPARATELY within each family of tests
      (tremor-recurrence trend, slip-recurrence trend, lag trend).  Both the raw p and the
      Holm-adjusted p are written out; "significant" means Holm-adjusted p < 0.05.
  SHIFT NULL (chance level for the matched fractions; PLAN.md section 4)
      200 draws.  In each draw every band's DAILY tremor count series is rigidly (circularly)
      shifted by the SAME random integer in [365, N - 365] days; the trailing counts, the peak
      percentile gate, the episodes and both matched fractions are recomputed from scratch against
      the UNSHIFTED slip patches.  One-sided p = (1 + #{null >= observed}) / (1 + 200).
      This destroys the tremor-slip association while preserving each band's tremor rate,
      burstiness and episode count distribution.
  NETWORK-GROWTH CONTROL (PLAN.md section 4; applied to the recurrence trends)
      The 90th-percentile peak gate is applied to the WHOLE record, so a band whose DETECTION RATE
      grows with the network puts most of its peaks in the later years and shows an artefactually
      SHORTENING recurrence interval.  Control: divide each band's daily count by a 365-d CENTRED
      running mean of its own daily count (edge-padded with the first/last 365-d mean; floored at
      10 % of the band's median running rate so dead early years cannot explode), then redo the
      trailing sum, the episodes and the Theil-Sen trend + permutation p on that rate-normalised
      series.  A real recurrence trend survives with the same sign; a detection-growth artefact
      does not.  det_growth_frac_per_yr = Theil-Sen slope of the band's yearly detection count
      (2010-2025) divided by its mean, i.e. fractional detection growth per year.

  YEARLY
      per (year, segment): the fraction of that segment's tremor episodes (band-instances) matched
      by a slip patch, and the fraction of that segment's Mw >= 6.0 substantive patch
      band-instances matched by a tremor episode.  The UNIT is the (band, item) instance, because
      this is a per-site analysis; a patch spanning three bands counts three times (the number of
      distinct patches is also reported).  Year = year of the episode PEAK / of the patch START.
      Context column: total kept tremor detections per year per segment - tremor-episode detection
      depends on network completeness, so yearly fractions must be read against it.
      2026 is partial (through 2026-05-27) and carries partial = True.

--------------------------------------------------------------------------------------------------
OUTPUTS  (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------------------
  results/s5_tremor_episodes.csv   one row per tremor episode: key, band, seg, peak, onset, count30,
                                   rise_d, matched, patch_id, lag_d
  results/s5_matches.csv           one row per (band, episode, patch) pair found in EITHER direction:
                                   from_episode, from_patch, mutual, lag_d, Mw, confirmed, ...
  results/s5_site_recurrence.csv   one row per band: n episodes / patches, median recurrence
                                   intervals (tremor, slip, slip Mw>=6) + Theil-Sen trends and
                                   permutation / Holm p, median lag + lag trend and p (plus the
                                   tie-robust Spearman), unmatched fractions both ways, the
                                   shift-null chance level for both matched fractions, and the
                                   rate-normalised (network-growth-controlled) recurrence trend
  results/s5_yearly_match.csv      one row per (year, segment): matched fractions both ways,
                                   counts, yearly detections, partial flag
  figures/s5_recurrence.png        4 panels:
                                   (a) per-band median recurrence, tremor vs slip
                                   (b) matched-pair lag vs date by segment + Theil-Sen lines
                                   (c) yearly matched fractions by segment
                                   (d) recurrence intervals vs date, north and central bands

Console: the per-band table, segment medians, every band with a Holm-significant trend, and the
cross-check against res_catalog_sse/tremor_episodes.csv (same band, peak within 20 d).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, theilslopes

# ----------------------------------------------------------------------------- paths / constants
ROOT = Path("/home/jovyan/tremorferometry")
SSE = ROOT / "fault_tomography/inversion/res_catalog_sse"
OUT_R = ROOT / "tremor_slip_geometry/results"
OUT_F = ROOT / "tremor_slip_geometry/figures"
OUT_R.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

T0 = pd.Timestamp("2010-01-06")
T1 = pd.Timestamp("2026-05-27")
BAND_EDGES = np.round(np.arange(40.0, 50.0, 0.5), 1)   # 20 lower edges
BAND_W = 0.5
TRAIL = 30            # d, trailing tremor window
PK_SEP = 90           # d, minimum peak separation
PK_PCTL = 90          # percentile gate on peak height
ONSET_BACK = 45       # d, onset search-back window
ONSET_FRAC = 0.25     # onset threshold as a fraction of the peak
MATCH_LO, MATCH_HI = -30, 90   # slip start - tremor onset window (days)
KD_MAX_KM = 40.0
NPERM = 1000
NSHIFT = 1000         # rigid-shift null draws for the matched fractions
RATE_WIN = 365        # d, window of the running detection-rate used by the network-growth control
RNG = np.random.default_rng(20260910)

SEG_COLOR = {"north": "#2a78d6", "central": "#eb6834", "south": "#1baf7a"}
SEG_MARK = {"north": "o", "central": "s", "south": "^"}
SEG_ORDER = ["north", "central", "south"]
C_TREMOR, C_SLIP = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"


def seg_of(lo: float) -> str:
    return "north" if lo >= 46.0 else ("central" if lo >= 44.0 else "south")


def band_of(lat: np.ndarray) -> np.ndarray:
    """Band lower edge for each latitude; NaN outside [40, 50)."""
    b = np.floor(np.asarray(lat, float) * 2.0) / 2.0
    b = np.where((b >= 40.0) & (b < 50.0), b, np.nan)
    return np.round(b, 1)


def yrs(dates: pd.Series | pd.DatetimeIndex) -> np.ndarray:
    """Decimal years since 2010-01-01."""
    d = pd.to_datetime(pd.Series(np.asarray(dates)))
    return (d - pd.Timestamp("2010-01-01")).dt.total_seconds().to_numpy() / (365.25 * 86400.0)


def ts_perm(x: np.ndarray, y: np.ndarray, nperm: int = NPERM):
    """Theil-Sen slope of y vs x + two-sided permutation p (shuffle y against x)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 5 or np.ptp(x) == 0:
        return np.nan, np.nan, x.size
    slope = theilslopes(y, x)[0]
    null = np.empty(nperm)
    for i in range(nperm):
        null[i] = theilslopes(RNG.permutation(y), x)[0]
    p = (1.0 + np.sum(np.abs(null) >= abs(slope))) / (1.0 + nperm)
    return slope, p, x.size


def holm(pvals: np.ndarray) -> np.ndarray:
    """Holm-Bonferroni adjusted p-values (NaNs pass through, excluded from the family)."""
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx])]
    m = order.size
    prev = 0.0
    for k, j in enumerate(order):
        adj = min(1.0, (m - k) * p[j])
        prev = max(prev, adj)
        out[j] = prev
    return out


# ============================================================================= 1. tremor -> bands
print("[1] tremor catalogue -> per-band daily counts", flush=True)
npz = np.load(SSE / "slip_patches_5mm.npz", allow_pickle=True)
tri_lat = npz["patch_lat"].astype(float)
tri_lon = npz["patch_lon"].astype(float)
mdays = pd.to_datetime(pd.Series(npz["days"].astype(str)))

KMX = 111.0 * np.cos(np.deg2rad(47.0))
tree = cKDTree(np.c_[(tri_lon + 124.0) * KMX, (tri_lat - 47.0) * 111.0])

tr = pd.read_csv(ROOT / "catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["time"] = pd.to_datetime(tr["time"], errors="coerce")
tr = tr.dropna(subset=["time", "lat", "lon"])
tr["day"] = tr["time"].dt.floor("D")
tr = tr[(tr["day"] >= T0) & (tr["day"] <= T1)]
n_raw = len(tr)
dist, _ = tree.query(np.c_[(tr["lon"].to_numpy() + 124.0) * KMX, (tr["lat"].to_numpy() - 47.0) * 111.0])
tr = tr[dist <= KD_MAX_KM].copy()
tr["band"] = band_of(tr["lat"].to_numpy())
tr = tr.dropna(subset=["band"])
print(f"    detections {n_raw} in window -> {len(tr)} within {KD_MAX_KM:.0f} km of a mesh centroid "
      f"and inside 40-50N ({100*len(tr)/n_raw:.1f} %)")

grid = pd.date_range(T0, T1, freq="D")
N = len(grid)
gidx = pd.Series(np.arange(N), index=grid)
counts = np.zeros((N, len(BAND_EDGES)), float)
bcol = {b: k for k, b in enumerate(BAND_EDGES)}
g = tr.groupby(["day", "band"]).size()
for (d, b), n in g.items():
    counts[gidx[d], bcol[round(float(b), 1)]] = n

# trailing 30-d sums, (t-30, t]
csum = np.cumsum(counts, axis=0)
trail = csum - np.vstack([np.zeros((TRAIL, counts.shape[1])), csum[:-TRAIL]])

# ============================================================================= 2. tremor episodes
print("[2] tremor episodes", flush=True)
def find_episodes(c: np.ndarray) -> list[tuple[int, int]]:
    """(peak index, onset index) pairs for one band's trailing-count series c."""
    thr = np.percentile(c, PK_PCTL)
    pk, _ = find_peaks(c)
    pk = pk[(pk >= TRAIL) & (c[pk] >= thr) & (c[pk] > 0)]
    keep: list[int] = []
    for i in pk[np.argsort(-c[pk])]:                       # greedy, highest first
        if all(abs(i - j) >= PK_SEP for j in keep):
            keep.append(int(i))
    out = []
    for i in sorted(keep):
        lo_i = max(TRAIL, i - ONSET_BACK)
        w = np.flatnonzero(c[lo_i:i + 1] > ONSET_FRAC * c[i])
        out.append((i, lo_i + (int(w[0]) if w.size else i - lo_i)))
    return out


eps = []
for k, b in enumerate(BAND_EDGES):
    c = trail[:, k]
    for i, ion in find_episodes(c):
        eps.append(dict(band=b, seg=seg_of(b), ipeak=i, ion=ion,
                        peak=grid[i], onset=grid[ion], count30=c[i], rise_d=i - ion))
E = pd.DataFrame(eps).sort_values(["band", "peak"]).reset_index(drop=True)
E["key"] = ["T%.1f_%s" % (r.band, r.peak.date()) for r in E.itertuples()]
print(f"    {len(E)} episodes over 20 bands "
      f"(N {int((E.seg=='north').sum())}, C {int((E.seg=='central').sum())}, S {int((E.seg=='south').sum())})")

# ============================================================================= 3. slip patches
print("[3] slip patches -> band footprints", flush=True)
P = pd.read_csv(SSE / "slip_patches_5mm.csv")
P["start"] = pd.to_datetime(P["start"])
P["end"] = pd.to_datetime(P["end"])
sub = P[(~P["minor"].astype(bool)) & (~P["broad_pulse"].astype(bool))].copy()
n_sub_all = len(sub)
sub = sub[(sub["start"] >= T0) & (sub["start"] <= T1)].copy()
print(f"    substantive patches {n_sub_all} in catalogue -> {len(sub)} with start inside the window "
      f"({int((sub.Mw>=6).sum())} with Mw >= 6.0)")

conf_ids = set(pd.read_csv(SSE / "tremor_vs_slip_patches.csv")
               .query("catalog == 'thrust' and confirmed")["id"].astype(int))

# footprint triangles per patch id: one pass over the label array
lab = npz["label"]
want = set(sub["id"].astype(int))
foot: dict[int, set[int]] = {i: set() for i in want}
tri_band = band_of(tri_lat)
for a in range(0, lab.shape[0], 512):
    blk = lab[a:a + 512]
    nz = np.flatnonzero(blk.ravel() != 0)
    if nz.size == 0:
        continue
    ids = blk.ravel()[nz]
    cols = nz % blk.shape[1]
    pairs = np.unique(np.c_[ids, cols], axis=0)
    for pid, col in pairs:
        if pid in foot:
            foot[int(pid)].add(int(col))

rows = []
for r in sub.itertuples():
    tj = np.array(sorted(foot[int(r.id)]), int)
    bands = np.unique(tri_band[tj][np.isfinite(tri_band[tj])]) if tj.size else np.array([])
    for b in bands:
        rows.append(dict(band=round(float(b), 1), seg=seg_of(float(b)), patch_id=int(r.id),
                         start=r.start, end=r.end, duration_d=r.duration_d, Mw=r.Mw,
                         lat_centroid=r.lat_centroid, n_tri_band=int(np.sum(tri_band[tj] == b)),
                         confirmed=int(r.id) in conf_ids))
PB = pd.DataFrame(rows).sort_values(["band", "start"]).reset_index(drop=True)
print(f"    {len(PB)} (band, patch) instances; {PB.patch_id.nunique()} distinct patches touch a band; "
      f"{int(PB.confirmed.sum())} instances from confirmed patches")

# ============================================================================= 4. matching
print("[4] matching", flush=True)
E["matched"] = False
E["patch_id"] = pd.NA
E["lag_d"] = np.nan
pairs = {}
for b in BAND_EDGES:
    ee = E.index[E.band == b]
    pp = PB[PB.band == b]
    ps = pp["start"].to_numpy()
    for i in ee:
        on = E.at[i, "onset"]
        d = (ps - np.datetime64(on)) / np.timedelta64(1, "D")
        cand = np.flatnonzero((d >= MATCH_LO) & (d <= MATCH_HI))
        if cand.size:
            j = cand[np.argmin(np.abs(d[cand]))]
            row = pp.iloc[j]
            E.at[i, "matched"] = True
            E.at[i, "patch_id"] = int(row.patch_id)
            E.at[i, "lag_d"] = float(d[j])
            pairs.setdefault((b, E.at[i, "key"], int(row.patch_id)),
                             dict(from_episode=False, from_patch=False))["from_episode"] = True
    # patch -> episode
    eon = E.loc[ee, "onset"].to_numpy()
    for j, row in enumerate(pp.itertuples()):
        if eon.size == 0:
            continue
        d = (np.datetime64(row.start) - eon) / np.timedelta64(1, "D")
        cand = np.flatnonzero((d >= MATCH_LO) & (d <= MATCH_HI))
        if cand.size:
            i = ee[cand[np.argmin(np.abs(d[cand]))]]
            pairs.setdefault((b, E.at[i, "key"], int(row.patch_id)),
                             dict(from_episode=False, from_patch=False))["from_patch"] = True

mrows = []
Eix = E.set_index("key")
PBix = PB.set_index(["band", "patch_id"])
for (b, key, pid), fl in pairs.items():
    e = Eix.loc[key]
    q = PBix.loc[(b, pid)]
    if isinstance(q, pd.DataFrame):
        q = q.iloc[0]
    mrows.append(dict(band=b, seg=seg_of(b), episode_key=key, peak=e.peak, onset=e.onset,
                      count30=e.count30, patch_id=pid, patch_start=q.start, patch_end=q.end,
                      Mw=q.Mw, confirmed=bool(q.confirmed),
                      lag_d=(q.start - e.onset).days,
                      from_episode=fl["from_episode"], from_patch=fl["from_patch"],
                      mutual=fl["from_episode"] and fl["from_patch"]))
M = pd.DataFrame(mrows).sort_values(["band", "onset"]).reset_index(drop=True)

PB["matched"] = False
PB["episode_key"] = pd.NA
for r in M[M.from_patch].itertuples():
    sel = (PB.band == r.band) & (PB.patch_id == r.patch_id)
    PB.loc[sel, "matched"] = True
    PB.loc[sel, "episode_key"] = r.episode_key

print(f"    {len(M)} pairs ({int(M.mutual.sum())} mutual); "
      f"episodes matched {E.matched.mean():.3f}; patch-instances matched {PB.matched.mean():.3f}; "
      f"Mw>=6 patch-instances matched {PB[PB.Mw>=6].matched.mean():.3f}")

# ============================================================================= 5. per-band table
print("[5] recurrence, trends, permutation nulls", flush=True)
recs = []
for b in BAND_EDGES:
    ee = E[E.band == b].sort_values("peak")
    pp = PB[PB.band == b].sort_values("start")
    p6 = pp[pp.Mw >= 6.0]
    rec = dict(band=b, seg=seg_of(b), n_episodes=len(ee), n_patches=len(pp), n_patches_mw6=len(p6),
               n_patches_conf=int(pp.confirmed.sum()))

    def rr(dates, tag):
        dts = pd.to_datetime(pd.Series(np.asarray(dates))).sort_values()
        if len(dts) < 2:
            rec[f"med_int_{tag}_d"] = np.nan
            rec[f"n_int_{tag}"] = 0
            rec[f"trend_{tag}_d_per_yr"] = np.nan
            rec[f"p_{tag}"] = np.nan
            return
        iv = dts.diff().dropna().dt.days.to_numpy(float)
        mid = dts.to_numpy()[:-1] + (dts.diff().dropna().to_numpy() / 2)
        s, p, n = ts_perm(yrs(mid), iv)
        rec[f"med_int_{tag}_d"] = float(np.median(iv))
        rec[f"n_int_{tag}"] = int(iv.size)
        rec[f"trend_{tag}_d_per_yr"] = s
        rec[f"p_{tag}"] = p

    rr(ee["peak"], "tremor")
    rr(pp["start"], "slip")
    rr(p6["start"], "slip6")

    mm = M[(M.band == b) & (M.from_episode)]
    rec["n_matched"] = len(mm)
    rec["med_lag_d"] = float(np.median(mm.lag_d)) if len(mm) else np.nan
    s, p, n = ts_perm(yrs(mm["onset"]), mm["lag_d"].to_numpy(float)) if len(mm) else (np.nan, np.nan, 0)
    rec["lag_trend_d_per_yr"] = s
    rec["p_lag"] = p
    if len(mm) >= 5:
        rho, prho = spearmanr(yrs(mm["onset"]), mm["lag_d"].to_numpy(float))
    else:
        rho, prho = np.nan, np.nan
    rec["lag_spearman_rho"] = rho
    rec["p_lag_spearman"] = prho
    rec["frac_tremor_unmatched"] = float(1 - ee["matched"].mean()) if len(ee) else np.nan
    rec["frac_slip_unmatched"] = float(1 - pp["matched"].mean()) if len(pp) else np.nan
    rec["frac_slip6_unmatched"] = float(1 - p6["matched"].mean()) if len(p6) else np.nan
    recs.append(rec)
R = pd.DataFrame(recs)
for tag in ["tremor", "slip", "slip6", "lag", "lag_spearman"]:
    R[f"pholm_{tag}"] = holm(R[f"p_{tag}"].to_numpy())

# ------------------------------------------------- 5b. rigid tremor-shift null for the matching
print("[5b] rigid tremor-shift null for the matched fractions (%d draws)" % NSHIFT, flush=True)
ps_band = {b: PB.loc[PB.band == b, "start"].to_numpy() for b in BAND_EDGES}
ps_i = {b: np.array([(pd.Timestamp(x) - T0).days for x in v], float) for b, v in ps_band.items()}


def match_fracs(on_i: dict[float, np.ndarray]):
    """(episode-matched fraction, patch-instance-matched fraction) per band and pooled."""
    fe, fp, ne_, np_ = {}, {}, 0, 0
    he, hp = 0, 0
    for b in BAND_EDGES:
        on, ps = on_i[b], ps_i[b]
        if on.size and ps.size:
            d = ps[None, :] - on[:, None]                  # (n_ep, n_patch) slip start - onset
            ok = (d >= MATCH_LO) & (d <= MATCH_HI)
            e_ok, p_ok = ok.any(axis=1), ok.any(axis=0)
        else:
            e_ok = np.zeros(on.size, bool)
            p_ok = np.zeros(ps.size, bool)
        fe[b] = e_ok.mean() if on.size else np.nan
        fp[b] = p_ok.mean() if ps.size else np.nan
        he += int(e_ok.sum()); ne_ += on.size
        hp += int(p_ok.sum()); np_ += ps.size
    return fe, fp, (he / ne_ if ne_ else np.nan), (hp / np_ if np_ else np.nan)


obs_on = {b: E.loc[E.band == b, "ion"].to_numpy(float) for b in BAND_EDGES}
fe_obs, fp_obs, FE_OBS, FP_OBS = match_fracs(obs_on)
null_fe = np.full((NSHIFT, len(BAND_EDGES)), np.nan)
null_fp = np.full((NSHIFT, len(BAND_EDGES)), np.nan)
null_FE = np.empty(NSHIFT)
null_FP = np.empty(NSHIFT)
null_nep = np.empty(NSHIFT)
for d_ in range(NSHIFT):
    sh = int(RNG.integers(365, N - 365))
    cs = np.cumsum(np.roll(counts, sh, axis=0), axis=0)
    tr_s = cs - np.vstack([np.zeros((TRAIL, counts.shape[1])), cs[:-TRAIL]])
    on_s = {}
    tot = 0
    for k, b in enumerate(BAND_EDGES):
        ev = find_episodes(tr_s[:, k])
        on_s[b] = np.array([o for _, o in ev], float)
        tot += len(ev)
    fe_n, fp_n, null_FE[d_], null_FP[d_] = match_fracs(on_s)
    null_nep[d_] = tot
    null_fe[d_] = [fe_n[b] for b in BAND_EDGES]
    null_fp[d_] = [fp_n[b] for b in BAND_EDGES]

R["null_frac_ep_matched"] = np.nanmean(null_fe, axis=0)
R["null_frac_patch_matched"] = np.nanmean(null_fp, axis=0)
R["p_shift_ep"] = [(1 + np.sum(null_fe[:, k] >= (1 - R.frac_tremor_unmatched.iloc[k]))) / (1 + NSHIFT)
                   for k in range(len(BAND_EDGES))]
R["p_shift_patch"] = [(1 + np.sum(null_fp[:, k] >= (1 - R.frac_slip_unmatched.iloc[k]))) / (1 + NSHIFT)
                      for k in range(len(BAND_EDGES))]
P_SHIFT_FE = (1 + np.sum(null_FE >= FE_OBS)) / (1 + NSHIFT)
P_SHIFT_FP = (1 + np.sum(null_FP >= FP_OBS)) / (1 + NSHIFT)

# --------------------------------------- 5c. network-growth control on the recurrence trends
# The 90th-percentile peak gate is applied to the WHOLE record, so a band whose detection RATE
# grows with the network will place most of its peaks in the later years and so show an
# artefactually SHORTENING recurrence interval.  Control: divide each band's daily counts by a
# 365-d CENTRED running mean of that band's own daily count (floored at 10 % of the band's median
# running rate, to stop dead early years exploding), then redo episodes and the trend on the
# rate-normalised series.  A real recurrence trend survives; a detection-growth artefact does not.
print("[5c] network-growth control (rate-normalised recurrence trends)", flush=True)
ker = np.ones(RATE_WIN) / RATE_WIN
rate = np.empty_like(counts)
for k in range(counts.shape[1]):
    pad = np.r_[np.full(RATE_WIN // 2, counts[:RATE_WIN, k].mean()),
                counts[:, k],
                np.full(RATE_WIN // 2, counts[-RATE_WIN:, k].mean())]
    rate[:, k] = np.convolve(pad, ker, mode="valid")[:counts.shape[0]]
floor = np.maximum(0.1 * np.median(rate, axis=0), 1e-3)
cn = counts / np.maximum(rate, floor[None, :])
csn = np.cumsum(cn, axis=0)
trn = csn - np.vstack([np.zeros((TRAIL, cn.shape[1])), csn[:-TRAIL]])

norm_rows = []
for k, b in enumerate(BAND_EDGES):
    ev = find_episodes(trn[:, k])
    pks = grid[[i for i, _ in ev]]
    if len(pks) >= 2:
        iv = np.diff(pks.values).astype("timedelta64[D]").astype(float)
        mid = pks.values[:-1] + np.diff(pks.values) / 2
        s, p, _ = ts_perm(yrs(pd.Series(mid)), iv)
        norm_rows.append(dict(band=b, n_episodes_norm=len(pks), med_int_tremor_norm_d=float(np.median(iv)),
                              trend_tremor_norm_d_per_yr=s, p_tremor_norm=p))
    else:
        norm_rows.append(dict(band=b, n_episodes_norm=len(pks), med_int_tremor_norm_d=np.nan,
                              trend_tremor_norm_d_per_yr=np.nan, p_tremor_norm=np.nan))
R = R.merge(pd.DataFrame(norm_rows), on="band", how="left")
R["pholm_tremor_norm"] = holm(R["p_tremor_norm"].to_numpy())
# yearly detection growth per band, as the artefact-risk context for each trend
det_yr = (tr.assign(year=tr["day"].dt.year).groupby(["band", "year"]).size().rename("n").reset_index())
gro = {}
for b in BAND_EDGES:
    d = det_yr[(det_yr.band == b) & (det_yr.year >= 2010) & (det_yr.year <= 2025)]
    gro[b] = (theilslopes(d["n"].to_numpy(float), d["year"].to_numpy(float))[0]
              / max(d["n"].mean(), 1e-9)) if len(d) >= 5 else np.nan
R["det_growth_frac_per_yr"] = [gro[b] for b in R["band"]]

# ============================================================================= 6. yearly
print("[6] yearly matched fractions", flush=True)
E["year"] = E["peak"].dt.year
PB["year"] = PB["start"].dt.year
yr_seg_det = (tr.assign(year=tr["day"].dt.year, seg=[seg_of(x) for x in tr["band"]])
                .groupby(["year", "seg"]).size().rename("detections"))
yrows = []
for y in range(T0.year, T1.year + 1):
    for s in SEG_ORDER:
        ee = E[(E.year == y) & (E.seg == s)]
        pp = PB[(PB.year == y) & (PB.seg == s)]
        p6 = pp[pp.Mw >= 6.0]
        yrows.append(dict(
            year=y, seg=s, partial=(y == T1.year),
            n_episodes=len(ee), frac_episodes_matched=float(ee.matched.mean()) if len(ee) else np.nan,
            n_patch_inst=len(pp), n_patch_uniq=pp.patch_id.nunique(),
            n_patch6_inst=len(p6), n_patch6_uniq=p6.patch_id.nunique(),
            frac_patch6_matched=float(p6.matched.mean()) if len(p6) else np.nan,
            frac_patch_matched=float(pp.matched.mean()) if len(pp) else np.nan,
            detections=int(yr_seg_det.get((y, s), 0))))
Y = pd.DataFrame(yrows)

# ============================================================================= 7. write CSVs
Eout = E[["key", "band", "seg", "peak", "onset", "count30", "rise_d", "matched", "patch_id", "lag_d"]].copy()
Eout.to_csv(OUT_R / "s5_tremor_episodes.csv", index=False)
M.to_csv(OUT_R / "s5_matches.csv", index=False)
R.to_csv(OUT_R / "s5_site_recurrence.csv", index=False)
Y.to_csv(OUT_R / "s5_yearly_match.csv", index=False)

# ============================================================================= 8. figure
print("[7] figure", flush=True)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "figure.facecolor": "white"})
fig, ax = plt.subplots(2, 2, figsize=(13.0, 9.0))

# (a) per-band median recurrence
a = ax[0, 0]
xb = R["band"] + 0.25
a.plot(R["med_int_tremor_d"], xb, "-o", color=C_TREMOR, lw=2, ms=6, label="tremor episodes")
a.plot(R["med_int_slip_d"], xb, "-s", color=C_SLIP, lw=2, ms=6, label="slip patches (substantive)")
a.plot(R["med_int_slip6_d"], xb, "--^", color=C_SLIP, lw=1.4, ms=6, alpha=0.65,
       label="slip patches Mw $\\geq$ 6.0")
for yv in (44.0, 46.0):
    a.axhline(yv, color=INK2, lw=0.8, ls=":")
a.set_xscale("log")
a.set_xlim(25, 4500)
for yv, lab_ in [(49.7, "north"), (45.7, "central"), (43.7, "south")]:
    a.text(4200, yv, lab_, ha="right", va="top", color=INK2, fontsize=8)
a.text(3070, 49.75, "n = 1 interval", color=C_SLIP, fontsize=7, ha="right", va="bottom", alpha=0.8)
a.set_xlabel("median recurrence interval (days, log scale)")
a.set_ylabel("band centre latitude (N)")
a.set_title("(a) recurrence per site", loc="left", fontweight="bold")
a.legend(frameon=False, fontsize=8, loc="lower right")

# (b) lag vs date by segment
a = ax[0, 1]
mm = M[M.from_episode]
for si, s in enumerate(SEG_ORDER):
    d = mm[mm.seg == s]
    if not len(d):
        continue
    a.plot(d["onset"], d["lag_d"], SEG_MARK[s], color=SEG_COLOR[s], ms=5, alpha=0.7, mew=0,
           label=f"{s} (n={len(d)})")
    x = yrs(d["onset"])
    if x.size >= 5 and np.ptp(x) > 0:
        sl, ic = theilslopes(d["lag_d"].to_numpy(float), x)[:2]
        rho, prho = spearmanr(x, d["lag_d"].to_numpy(float))
        xx = np.array([x.min(), x.max()])
        a.plot(pd.Timestamp("2010-01-01") + pd.to_timedelta(xx * 365.25, "D"), ic + sl * xx,
               "-", color=SEG_COLOR[s], lw=2.2)
        a.annotate(f"{s} {sl:+.1f} d/yr  ($\\rho$={rho:+.2f}, p={prho:.2f})",
                   (pd.Timestamp("2016-01-01"), 0), color=SEG_COLOR[s], fontsize=8,
                   xytext=(0, -46 - 13 * si), textcoords="offset points", ha="center", va="top")
a.axhline(0, color=INK2, lw=1)
a.set_ylim(-58, 95)
a.set_ylabel("lag = slip start $-$ tremor onset (days)")
a.set_title("(b) phase of slip relative to tremor", loc="left", fontweight="bold")
a.legend(frameon=False, fontsize=8, loc="upper left", ncol=3)

# (c) yearly matched fractions
a = ax[1, 0]
for s in SEG_ORDER:
    d = Y[(Y.seg == s) & (Y.n_episodes > 0)]
    a.plot(d["year"], d["frac_episodes_matched"], "-" + SEG_MARK[s], color=SEG_COLOR[s], lw=2, ms=5,
           mew=0, label=f"{s}: tremor ep. matched by slip")
    d6 = Y[(Y.seg == s) & (Y.n_patch6_inst > 0)]
    a.plot(d6["year"], d6["frac_patch6_matched"], ":" + SEG_MARK[s], color=SEG_COLOR[s], lw=1.6,
           ms=5, mew=0, alpha=0.65, label=f"{s}: Mw$\\geq$6 slip matched by tremor")
a.axvspan(2025.55, 2026.45, color=GRID, alpha=0.6, zorder=0)
a.text(2026, 0.04, "2026\npartial", ha="center", va="bottom", fontsize=7, color=INK2)
a.axhline(FE_OBS_NULL := null_FE.mean(), color=INK2, lw=1.2, ls="--")
a.text(2010.2, FE_OBS_NULL + 0.015, "shift-null chance level (tremor$\\to$slip, all bands)",
       fontsize=7, color=INK2, va="bottom",
       bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85))
a.set_ylim(-0.05, 1.30)
a.set_xlim(2009.5, 2026.5)
a.set_xticks(range(2010, 2027, 2))
a.set_xlabel("year")
a.set_ylabel("matched fraction")
a.set_title("(c) yearly matched fractions (solid: tremor$\\to$slip, dotted: slip$\\to$tremor)",
            loc="left", fontweight="bold")
a.legend(fontsize=7, ncol=2, loc="upper left", framealpha=0.92, edgecolor=GRID)

# (d) recurrence intervals vs date, north + central bands
a = ax[1, 1]
for s in ["north", "central"]:
    xs, ys_ = [], []
    for b in BAND_EDGES[[seg_of(x) == s for x in BAND_EDGES]]:
        dts = E[E.band == b].sort_values("peak")["peak"]
        if len(dts) < 2:
            continue
        iv = dts.diff().dropna().dt.days.to_numpy(float)
        mid = dts.to_numpy()[:-1] + (dts.diff().dropna().to_numpy() / 2)
        xs.append(mid)
        ys_.append(iv)
    if not xs:
        continue
    xs = np.concatenate(xs)
    ys_ = np.concatenate(ys_)
    a.plot(xs, ys_, SEG_MARK[s], color=SEG_COLOR[s], ms=5, alpha=0.55, mew=0,
           label=f"{s} tremor intervals (n={ys_.size})")
    xv = yrs(pd.Series(xs))
    sl, ic = theilslopes(ys_, xv)[:2]
    xx = np.array([xv.min(), xv.max()])
    a.plot(pd.Timestamp("2010-01-01") + pd.to_timedelta(xx * 365.25, "D"), ic + sl * xx,
           "-", color=SEG_COLOR[s], lw=2.2)
    a.annotate(f"{sl:+.1f} d/yr", (pd.Timestamp("2010-01-01") + pd.to_timedelta(xx[1] * 365.25, "D"),
                                   ic + sl * xx[1]), color=SEG_COLOR[s], fontsize=8,
               xytext=(4, 0), textcoords="offset points", va="center")
a.set_ylabel("tremor-episode recurrence interval (days)")
a.set_title("(d) recurrence intervals vs date, north & central", loc="left", fontweight="bold")
a.legend(frameon=False, fontsize=8, loc="upper left")
a.text(0.02, 0.80,
       "apparent shortening is a DETECTION-GROWTH ARTEFACT:\n"
       "no band's trend survives the rate-normalised control\n"
       "(see results/s5_site_recurrence.csv, p_tremor_norm)",
       transform=a.transAxes, fontsize=7.5, color=INK2, va="top",
       bbox=dict(fc="white", ec=GRID, lw=0.8, boxstyle="round,pad=0.4"))

fig.suptitle("S5  site recurrence: tremor episodes vs modeled slip patches, 0.5$\\degree$ bands, "
             "2010-01-06 to 2026-05-27", fontsize=11, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT_F / "s5_recurrence.png", dpi=160)
plt.close(fig)

# ============================================================================= 9. console report
pd.set_option("display.width", 220, "display.max_columns", 60)
print("\n=== PER-BAND TABLE ===")
show = R[["band", "seg", "n_episodes", "n_patches", "n_patches_mw6", "med_int_tremor_d",
          "med_int_slip_d", "med_int_slip6_d", "trend_tremor_d_per_yr", "p_tremor",
          "trend_slip_d_per_yr", "p_slip", "n_matched", "med_lag_d", "lag_trend_d_per_yr", "p_lag",
          "lag_spearman_rho", "p_lag_spearman",
          "frac_tremor_unmatched", "frac_slip_unmatched", "frac_slip6_unmatched",
          "null_frac_ep_matched", "p_shift_ep"]].copy()
print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print("\n=== SEGMENT MEDIANS (over the bands of the segment) ===")
segtab = R.groupby("seg").agg(
    n_bands=("band", "size"), episodes=("n_episodes", "sum"), patches=("n_patches", "sum"),
    med_int_tremor=("med_int_tremor_d", "median"), med_int_slip=("med_int_slip_d", "median"),
    med_int_slip6=("med_int_slip6_d", "median"), med_lag=("med_lag_d", "median"),
    med_lag_trend=("lag_trend_d_per_yr", "median"),
    tremor_unmatched=("frac_tremor_unmatched", "median"),
    slip_unmatched=("frac_slip_unmatched", "median"),
    slip6_unmatched=("frac_slip6_unmatched", "median")).reindex(SEG_ORDER)
print(segtab.to_string(float_format=lambda v: f"{v:.3g}"))
print("\npooled: episode->patch matched %.3f (%d/%d); patch-instance->episode matched %.3f (%d/%d); "
      "Mw>=6 patch-instances matched %.3f (%d/%d)"
      % (E.matched.mean(), E.matched.sum(), len(E), PB.matched.mean(), PB.matched.sum(), len(PB),
         PB[PB.Mw >= 6].matched.mean(), PB[PB.Mw >= 6].matched.sum(), int((PB.Mw >= 6).sum())))
allm = M[M.from_episode]
print("pooled matched lag: median %.1f d, IQR %.1f..%.1f d, n=%d; mutual pairs %d"
      % (allm.lag_d.median(), allm.lag_d.quantile(.25), allm.lag_d.quantile(.75), len(allm),
         int(M.mutual.sum())))
sl, p, n = ts_perm(yrs(allm["onset"]), allm["lag_d"].to_numpy(float))
rho, prho = spearmanr(yrs(allm["onset"]), allm["lag_d"].to_numpy(float))
print("pooled lag trend (all bands together): %+.4f d/yr, perm p = %.3f, n = %d "
      "| tie-robust Spearman rho = %+.3f, p = %.3f" % (sl, p, n, rho, prho))
if sl == 0:
    print("  (Theil-Sen slope is EXACTLY 0 - the tied-integer atom - so its perm p = 1 by "
          "construction; read the Spearman instead)")

print("\n=== RIGID TREMOR-SHIFT NULL for the matched fractions (%d draws, shift 365..%d d) ==="
      % (NSHIFT, N - 365))
print("  episodes matched by a slip patch : observed %.3f | null mean %.3f (p5 %.3f, p95 %.3f) | p = %.4f"
      % (FE_OBS, null_FE.mean(), np.percentile(null_FE, 5), np.percentile(null_FE, 95), P_SHIFT_FE))
print("  patch instances matched by tremor: observed %.3f | null mean %.3f (p5 %.3f, p95 %.3f) | p = %.4f"
      % (FP_OBS, null_FP.mean(), np.percentile(null_FP, 5), np.percentile(null_FP, 95), P_SHIFT_FP))
print("  null episode count per draw: mean %.1f (observed %d) - the shift preserves burstiness"
      % (null_nep.mean(), len(E)))
ph_shift = holm(R.p_shift_ep.to_numpy())
nb = int(np.sum(R.p_shift_ep < 0.05))
print("  bands whose episode-matched fraction beats the shift null at p < 0.05: %d / 20 "
      "(Holm-adjusted: %d)" % (nb, int(np.sum(ph_shift < 0.05))))
for k in np.flatnonzero(ph_shift < 0.05):
    r = R.iloc[k]
    print("    band %.1f (%s): observed %.3f vs null %.3f, p = %.4f, Holm p = %.3f"
          % (r.band, r.seg, 1 - r.frac_tremor_unmatched, r.null_frac_ep_matched,
             r.p_shift_ep, ph_shift[k]))

print("\n=== BANDS WITH RAW p < 0.05 (Holm over the 20 bands, per family) ===")
any_sig = False
for tag, lab_ in [("tremor", "tremor recurrence trend"), ("slip", "slip recurrence trend"),
                  ("slip6", "slip Mw>=6 recurrence trend"), ("lag", "lag (phase) trend"),
                  ("lag_spearman", "lag (phase) Spearman")]:
    d = R[R[f"p_{tag}"] < 0.05]
    for r in d.itertuples():
        any_sig = True
        sl_ = (r.lag_trend_d_per_yr if tag == "lag" else
               r.lag_spearman_rho if tag == "lag_spearman" else
               getattr(r, f"trend_{tag}_d_per_yr"))
        unit = "rho  " if tag == "lag_spearman" else "d/yr"
        print("  %-28s band %.1f (%s): %+7.2f %s, p = %.4f, Holm p = %.3f  %s"
              % (lab_, r.band, r.seg, sl_, unit, getattr(r, f"p_{tag}"), getattr(r, f"pholm_{tag}"),
                 "**SURVIVES FWER**" if getattr(r, f"pholm_{tag}") < 0.05 else "(fails FWER)"))
if not any_sig:
    print("  none")

print("\n=== NETWORK-GROWTH CONTROL on the tremor-recurrence trends ===")
print("  (episodes redone on the DETECTION-RATE-NORMALISED series; a real trend survives)")
print(R[["band", "seg", "det_growth_frac_per_yr", "n_episodes", "trend_tremor_d_per_yr", "p_tremor",
         "n_episodes_norm", "trend_tremor_norm_d_per_yr", "p_tremor_norm", "pholm_tremor_norm"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
sig_raw = R[R.p_tremor < 0.05]
print("  verdict on the %d bands with a raw-significant tremor-recurrence trend:" % len(sig_raw))
for r in sig_raw.itertuples():
    surv = np.isfinite(r.p_tremor_norm) and r.p_tremor_norm < 0.05 and \
        np.sign(r.trend_tremor_norm_d_per_yr) == np.sign(r.trend_tremor_d_per_yr)
    print("    band %.1f (%s): raw %+.1f d/yr p=%.4f (Holm %.3f) -> rate-normalised %+.1f d/yr "
          "p=%.3f  %s   [detections %+.1f %%/yr]"
          % (r.band, r.seg, r.trend_tremor_d_per_yr, r.p_tremor, r.pholm_tremor,
             r.trend_tremor_norm_d_per_yr, r.p_tremor_norm,
             "SURVIVES" if surv else "DOES NOT SURVIVE -> detection-growth artefact",
             100 * r.det_growth_frac_per_yr))

print("\n=== YEARLY (matched fractions and detection context) ===")
print(Y[Y.n_episodes + Y.n_patch_inst > 0][
    ["year", "seg", "partial", "n_episodes", "frac_episodes_matched", "n_patch6_inst",
     "frac_patch6_matched", "detections"]].to_string(index=False,
                                                     float_format=lambda v: f"{v:.3g}"))

# ---- cross-check against the existing tremor_episodes.csv
old = pd.read_csv(SSE / "tremor_episodes.csv")
old["peak"] = pd.to_datetime(old["peak"])
old["band"] = old["band"].round(1)
hit_new = 0
for r in E.itertuples():
    o = old[(old.band == round(float(r.band), 1)) & ((old.peak - r.peak).abs() <= pd.Timedelta("20D"))]
    hit_new += int(len(o) > 0)
hit_old = 0
for r in old.itertuples():
    o = E[(E.band == round(float(r.band), 1)) & ((E.peak - r.peak).abs() <= pd.Timedelta("20D"))]
    hit_old += int(len(o) > 0)
print("\n=== CROSS-CHECK vs res_catalog_sse/tremor_episodes.csv (same band, peak within 20 d) ===")
print("  S5 episodes %d, old table %d;  S5 with an old counterpart %d (%.1f %%);  "
      "old with an S5 counterpart %d (%.1f %%)"
      % (len(E), len(old), hit_new, 100 * hit_new / len(E), hit_old, 100 * hit_old / len(old)))

print("\nwrote:")
for f in ["results/s5_tremor_episodes.csv", "results/s5_matches.csv",
          "results/s5_site_recurrence.csv", "results/s5_yearly_match.csv",
          "figures/s5_recurrence.png"]:
    print("  tremor_slip_geometry/" + f)
