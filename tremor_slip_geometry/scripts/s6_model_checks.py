#!/usr/bin/env python
"""
S6 - model-only checks on the Gualandi (2025) vbICA slip model (PLAN.md step S6).

Two independent questions, both answered with the slip model + the PNSN tremor catalogue only
(no dv/v; deferred by the user 2026-09-10):

  PART A  INTER-ETS SLIP.  Does the model place thrust slip BETWEEN ETS events, and does tremor
          independently confirm those inter-event slip bursts?  Run under THREE event-day
          definitions (A0 / A1 / A2, below) because the answer turned out to depend entirely on
          which days are called "ETS days".
  PART B  THE SOUTHERN +-5 mm OSCILLATION.  What is it - a seasonal cycle, a projection of a few
          vbICA components, or unstructured (station-like) model noise?

--------------------------------------------------------------------------------------------------
INPUTS (all read-only, all outside this folder)
--------------------------------------------------------------------------------------------------
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz
      label        (T x 742) int  connected-component labels; 0 = no patch, else the patch id
      inc30_mm     (T x 742) f32  modeled THRUST slip within the trailing 30 d, (t-30, t], in mm
                                  (NaN for the first 30 days)
      days         (T,)      str  daily 2007-01-02 .. 2026-08-18  (T = 7169)
      patch_lat/patch_lon/patch_depth_km/patch_area_km2  (742,)  mesh-triangle centroids + areas
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.csv     id,start,end,Mw,lat_centroid,minor,broad_pulse
  fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv
      per-patch tremor confirmation; catalog in {'thrust','mirror'}, id, seg, Mw, broad_pulse, confirmed
  fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz
      A_raw        (T x 13) f32  raw vbICA temporal amplitudes
  catalogs/pnsn_tremor_cascadia_full.csv   PNSN tremor master (time, lat, lon), ~748k rows 2010-01..2026-05
  /home/jovyan/cascadia-sse-strain  (src/sse_strain; data/cascadia_mesh.txt; data/2026-08-18)
      -> Ld (742 x 13) down-dip slip patterns; cumulative thrust slip [mm] = -1e3 * (Ld @ A_raw.T)

--------------------------------------------------------------------------------------------------
EXACT DEFINITIONS
--------------------------------------------------------------------------------------------------
  Analysis window        2010-01-06 .. 2026-05-27 (the PNSN tremor span).  N = 5986 days.
  Latitude bands         10 one-degree bands by TRIANGLE centroid latitude: [40,41), ..., [49,50).
                         (2 triangles with lat < 40 are outside every band and are dropped.)
  Segments               north lat >= 46, central 44 <= lat < 46, south lat < 44.
                         "north/central" = bands 44-50; "south" = bands 40-44.
  Daily slip rate        raw_cum_mm = -1e3 * (Ld @ A_raw.T)          (742 x T, cumulative thrust mm)
                         rate = d/dt along time (mm/d), then a 7-d CENTRED moving average.
                         Band rate = area-weighted mean of the band's triangles.
  Band 30-d slip         area-weighted mean of inc30_mm over the band's triangles (mm per 30 d).
  Tremor -> triangle     nearest mesh-triangle centroid by cKDTree in a local km frame
                         x = (lon + 124) * 111 * cos(47 deg), y = (lat - 47) * 111; kept if <= 40 km.
                         Band daily count = detections assigned to the band's triangles that day.
  Tremor 30-d count      trailing 30-d sum of the band daily counts (same (t-30, t] frame as inc30_mm).
  Patch footprint        the set of triangles j with label[t, j] == id for ANY day t; a patch "touches"
                         a band if its footprint contains a triangle of that band.

  EVENT DAYS - THREE DEFINITIONS (a band's days are EVENT or, complementarily, INTER-EVENT):
    A0  CONFIRMED patches (the first-run definition, kept for traceability).
        tremor_vs_slip_patches: catalog=='thrust' & confirmed & seg in ('north','central') &
        Mw >= 6.0 & ~broad_pulse   -> 27 patches (21 N, 6 C).
        Event day = within [start - 30, end + 30] of any such patch touching the band.
        KNOWN DEFECT: no confirmed patch is southern, so bands 40-42 have ZERO event days, and many
        real ETS episodes (Dec 2014, Jan-Feb 2018, Dec 2019, Feb-Mar 2021, ...) are not excluded.
    A1  SUBSTANTIVE THRUST patches: slip_patches_5mm.csv with ~minor & ~broad_pulse, ANY Mw, ANY
        segment  -> 284 patches.  Event day = within [start - 30, end + 30] of any such patch
        touching the band.  This is the model's own full ETS calendar.
    A2  TREMOR-BASED ETS calendar (model-independent): per band, take the runs of >= 5 consecutive
        days on which the band's 30-d tremor count exceeds that band's 90th percentile (over the
        analysis window); event day = within +- 30 d of any such run.

  Inter-event slip burst contiguous runs of INTER-EVENT days with band 30-d slip > 3 mm, length >= 5 d.
                         peak day = day of the maximum band 30-d slip inside the run.
                         tremor-confirmed  <=>  tremor 30-d count at the peak day >= 3 x (band
                         inter-event MEDIAN tremor 30-d count)  AND  >= 10 detections.
  Inter-event tremor burst contiguous runs of INTER-EVENT days with tremor 30-d count >= 3 x median
                         AND >= 30 detections, length >= 5 d; peak day = maximum tremor 30-d count.
                         slip-confirmed <=> band 30-d slip at the peak day > 3 mm.
  Shift null             200 draws; the band's DAILY tremor-count series is rigidly (circularly)
                         shifted by a random integer in [365, T-365] days and the 30-d counts, the
                         band inter-event median, and the confirmation fractions are recomputed.
                         One shift per draw is shared by all bands / groups / definitions (paired).
                         The EVENT-DAY mask is held FIXED at its observed value (for A2 that means
                         the tremor-defined calendar is not re-derived from the shifted series -
                         the null asks only whether tremor happens to sit on the burst peaks).
                         p = (1 + #{null >= observed}) / (1 + 200).
  Welch spectrum         scipy.signal.welch on the band-mean inc30_mm series, fs = 1 /day,
                         nperseg = 1461 (4 years), Hann, detrended linearly.
                         peak period = period of maximum power with 60 d <= period <= 1096 d (3 yr).
                         ann_ratio = mean power over periods 330-400 d divided by the MEDIAN power
                         over the neighbouring periods 180-300 d and 430-900 d.
  Seasonal amplitude     month-of-year mean of the band inc30 series; amplitude = max - min over the
                         12 months.  Null: 500 CIRCULAR shifts of the daily series by a random
                         integer in [30, 700] days; p = (1 + #{null >= observed}) / (1 + 500).
                         (This null preserves the series' own spectrum, so a shift that happens to be
                         near a multiple of 365 d keeps the annual phase -> the test is conservative.)
  Component variance     for vbICA component k and the southern triangles S:
                         var_cum[k]   = sum_{i in S} Var_t( -1e3 * Ld[i,k] * A_raw[:,k] )     (mm^2)
                         var_inc30[k] = same on the trailing-30-d increment of that series.
                         Reported as a fraction of the sum over the 13 components.

--------------------------------------------------------------------------------------------------
OUTPUTS (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------------------
  results/s6_interets_bands.csv          per DEFINITION (column 'definition' = A0/A1/A2) and band:
                                         positive-slip and tremor fractions on inter-event days, mean
                                         positive rate event vs inter-event, burst counts + confirmation
                                         rates with their shift nulls.
  results/s6_interets_bursts.csv         one row per inter-event burst, per definition (slip and tremor bursts).
  results/s6_interets_conf.csv           the 3 definitions x 2 groups x 2 burst kinds confirmation table
                                         with the shift-null mean / p95 / p.
  results/s6_south_oscillation_bands.csv per band: peak period, annual power ratio, seasonal amplitude + p,
                                         frac(inc30 > 5 mm), frac(inc30 < -5 mm), corr with sum of the
                                         13 raw amplitudes' trailing-30-d means.
  results/s6_sign_symmetry_year.csv      per segment per year: fraction of triangle-days |inc30| above +-5 mm.
  results/s6_south_components.csv        the 13 vbICA components ranked by southern variance.
  results/s6_band_corr.csv               10 x 10 Pearson correlation matrix of the band-mean inc30 series.
  figures/s6_model_checks.png            9 panels: A1a/A1b/A1c per-band fractions and rate ratio for the three
                                         definitions, A2a/A2b burst confirmation vs null, B2a spectra,
                                         B2b month-of-year cycle, B3 band correlations, B4 stationarity.
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.signal import welch

ROOT = "/home/jovyan/tremorferometry"
SSE = os.path.join(ROOT, "fault_tomography/inversion/res_catalog_sse")
HERE = os.path.join(ROOT, "tremor_slip_geometry")
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

WIN0, WIN1 = "2010-01-06", "2026-05-27"
BAND_EDGES = np.arange(40, 51)  # 10 bands [40,41) ... [49,50)
BAND_LABELS = [f"{lo}-{lo+1}" for lo in BAND_EDGES[:-1]]
N_NULL_A = 200
N_NULL_SEAS = 500
DEF_NAMES = ["A0", "A1", "A2"]
DEF_DESC = {
    "A0": "CONFIRMED patches (thrust & confirmed & N/C & Mw>=6 & ~broad_pulse), n=27",
    "A1": "SUBSTANTIVE thrust patches (~minor & ~broad_pulse, any Mw/segment)",
    "A2": "TREMOR ETS calendar (runs >=5 d of band 30-d tremor > band p90, +-30 d)",
}

RNG_A = np.random.default_rng(20260910)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 80)


def hdr(s):
    print("\n" + "=" * 100)
    print(s)
    print("=" * 100)


def runs_of_true(mask):
    """Contiguous runs of True -> list of (i0, i1) inclusive index pairs."""
    m = np.asarray(mask, bool).astype(np.int8)
    d = np.diff(np.concatenate(([0], m, [0])))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1))


def roll30(x):
    """Trailing 30-day rolling sum (min_periods=30 -> NaN before day 30)."""
    return pd.Series(x).rolling(30, min_periods=30).sum().to_numpy()


def nanmed(x):
    return np.nanmedian(x) if np.isfinite(x).any() else np.nan


# ==================================================================================================
# LOAD
# ==================================================================================================
hdr("S6 model-only checks - loading")

z = np.load(os.path.join(SSE, "slip_patches_5mm.npz"), allow_pickle=True)
label = z["label"]                        # (T, 742)
inc30 = z["inc30_mm"].astype(np.float64)  # (T, 742)
days = pd.to_datetime([str(d) for d in z["days"]])
plat, plon = z["patch_lat"], z["patch_lon"]
parea = z["patch_area_km2"]
T, NT = inc30.shape
print(f"slip_patches_5mm.npz: T={T} days {days[0].date()}..{days[-1].date()}, {NT} triangles")

patch_csv = pd.read_csv(os.path.join(SSE, "slip_patches_5mm.csv"))
tvs = pd.read_csv(os.path.join(SSE, "tremor_vs_slip_patches.csv"))
conf = tvs[
    (tvs.catalog == "thrust")
    & (tvs.confirmed.astype(bool))
    & (tvs.seg.isin(["north", "central"]))
    & (tvs.Mw >= 6.0)
    & (~tvs.broad_pulse.astype(bool))
].copy()
subst = patch_csv[(~patch_csv.minor.astype(bool)) & (~patch_csv.broad_pulse.astype(bool))].copy()
print(f"A0 CONFIRMED patches: n={len(conf)}  ({conf.seg.value_counts().to_dict()})")
print(f"A1 SUBSTANTIVE thrust patches (~minor & ~broad_pulse): n={len(subst)}  "
      f"Mw {subst.Mw.min():.2f}-{subst.Mw.max():.2f}, median duration {subst.duration_d.median():.0f} d")

# --- daily thrust slip rate per triangle -----------------------------------------------------
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution  # noqa: E402

fault = load_fault_ascii("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt")
sol = load_solution("/home/jovyan/cascadia-sse-strain/data/2026-08-18", fault)
sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()  # (742, 13)
A_raw = np.nan_to_num(
    np.load(os.path.join(SSE, "components_at_cells_trailing30.npz"), allow_pickle=True)["A_raw"]
).astype(np.float64)  # (T, 13)
raw_cum_mm = -1e3 * (Ld @ A_raw.T)  # (742, T) cumulative thrust slip, mm
rate = np.full_like(raw_cum_mm, np.nan)
rate[:, 1:] = np.diff(raw_cum_mm, axis=1)
rate = pd.DataFrame(rate.T).rolling(7, center=True, min_periods=1).mean().to_numpy().T
print(f"Ld {Ld.shape}, A_raw {A_raw.shape}, rate {rate.shape} mm/d")

# --- tremor -> triangles ---------------------------------------------------------------------
tre = pd.read_csv(
    os.path.join(ROOT, "catalogs/pnsn_tremor_cascadia_full.csv"),
    usecols=["time", "lat", "lon"],
    parse_dates=["time"],
).dropna(subset=["lat", "lon", "time"])
KX = 111.0 * np.cos(np.deg2rad(47.0))
tri_xy = np.column_stack([(plon + 124.0) * KX, (plat - 47.0) * 111.0])
tre_xy = np.column_stack([(tre.lon.to_numpy() + 124.0) * KX, (tre.lat.to_numpy() - 47.0) * 111.0])
dist, idx = cKDTree(tri_xy).query(tre_xy, k=1)
keep = dist <= 40.0
print(f"tremor: {len(tre)} detections {tre.time.min().date()}..{tre.time.max().date()}; "
      f"{keep.sum()} ({100*keep.mean():.1f}%) within 40 km of a triangle centroid")
tre = tre.loc[keep].copy()
tre["tri"] = idx[keep]
tre["day"] = tre.time.dt.floor("D")
day_pos = pd.Series(np.arange(T), index=days)
tre = tre[tre.day.isin(day_pos.index)]
tre["ti"] = day_pos.loc[tre.day].to_numpy()

# --- bands -------------------------------------------------------------------------------------
band_of_tri = np.full(NT, -1, int)
for b, lo in enumerate(BAND_EDGES[:-1]):
    band_of_tri[(plat >= lo) & (plat < lo + 1)] = b
print(f"triangles per band: {[int((band_of_tri==b).sum()) for b in range(10)]} "
      f"(+{int((band_of_tri<0).sum())} triangles with lat<40 dropped)")

tre_band = band_of_tri[tre.tri.to_numpy()]
ok = tre_band >= 0
cnt_day = np.zeros((T, 10))
np.add.at(cnt_day, (tre.ti.to_numpy()[ok], tre_band[ok]), 1.0)

band_rate = np.full((T, 10), np.nan)
band_inc30 = np.full((T, 10), np.nan)
for b in range(10):
    sel = band_of_tri == b
    w = parea[sel] / parea[sel].sum()
    band_rate[:, b] = rate[sel, :].T @ w
    band_inc30[:, b] = inc30[:, sel] @ w

wmask = (days >= pd.Timestamp(WIN0)) & (days <= pd.Timestamp(WIN1))
wi = np.flatnonzero(wmask)
NW = len(wi)
wdays = days[wi]
print(f"analysis window {WIN0}..{WIN1}: {NW} days")

tre30_full = np.column_stack([roll30(cnt_day[:, b]) for b in range(10)])
br_w = band_rate[wi, :]
bi_w = band_inc30[wi, :]
cnt_w = cnt_day[wi, :]
tre30_w = tre30_full[wi, :]
pos = np.clip(br_w, 0, None)

# ==================================================================================================
# PART A - inter-ETS slip, three event-day definitions
# ==================================================================================================
hdr("PART A - inter-ETS slip: three EVENT-DAY definitions (A0 confirmed / A1 substantive / A2 tremor)")

# --- patch footprints -> bands touched, in one pass over the label array ----------------------
lab_flat = label.ravel()
nzi = np.flatnonzero(lab_flat)
key = np.unique(lab_flat[nzi].astype(np.int64) * NT + (nzi % NT))
pid_u, tri_u = (key // NT).astype(int), (key % NT).astype(int)
order = np.argsort(pid_u, kind="stable")
pid_s, tri_s, bnd_s = pid_u[order], tri_u[order], band_of_tri[tri_u[order]]
cuts = np.flatnonzero(np.diff(pid_s)) + 1
foot_bands, foot_ntri = {}, {}
for a, c in zip(np.r_[0, cuts], np.r_[cuts, len(pid_s)]):
    foot_bands[int(pid_s[a])] = sorted(set(bnd_s[a:c].tolist()) - {-1})
    foot_ntri[int(pid_s[a])] = int(c - a)
print(f"footprints built for {len(foot_bands)} labelled patches")

dvals = days.values


def event_mask_from_patches(pdf):
    em = np.zeros((T, 10), bool)
    for _, r in pdf.iterrows():
        bb = foot_bands.get(int(r.id), [])
        if not bb:
            continue
        i0 = int(np.searchsorted(dvals, np.datetime64(pd.Timestamp(r.start) - pd.Timedelta(days=30))))
        i1 = int(np.searchsorted(dvals, np.datetime64(pd.Timestamp(r.end) + pd.Timedelta(days=30)), side="right"))
        for b in bb:
            em[max(i0, 0):min(i1, T), b] = True
    return em


def event_mask_tremor():
    em = np.zeros((T, 10), bool)
    p90s = np.full(10, np.nan)
    nrun = np.zeros(10, int)
    for b in range(10):
        xw = tre30_w[:, b]
        p90s[b] = np.nanpercentile(xw, 90)
        m = np.zeros(T, bool)
        m[wi] = xw > p90s[b]
        for a, c in runs_of_true(m):
            if c - a + 1 < 5:
                continue
            nrun[b] += 1
            em[max(a - 30, 0):min(c + 31, T), b] = True
    return em, p90s, nrun


EM = {"A0": event_mask_from_patches(conf), "A1": event_mask_from_patches(subst)}
EM["A2"], p90_tre, n_ets_runs = event_mask_tremor()
print(f"A2 tremor-ETS calendar: band 30-d tremor p90 = {np.round(p90_tre, 0).tolist()}, "
      f"ETS runs per band = {n_ets_runs.tolist()}")

IEV = {d: ~EM[d][wi, :] for d in DEF_NAMES}
EVE = {d: EM[d][wi, :] for d in DEF_NAMES}

groups = {"north/central": [b for b in range(10) if BAND_EDGES[b] >= 44],
          "south": [b for b in range(10) if BAND_EDGES[b] < 44]}
group_of_band = {b: g for g, bs in groups.items() for b in bs}

rows_bands, rows_bursts = [], []
SLIP_PEAKS = {d: {} for d in DEF_NAMES}
OBS = {d: {g: {} for g in groups} for d in DEF_NAMES}

for dname in DEF_NAMES:
    ev, iev = EVE[dname], IEV[dname]
    obs_pack = {g: dict(s_hits=0, s_n=0, t_hits=0, t_n=0) for g in groups}
    for b in range(10):
        lo = BAND_EDGES[b]
        seg = "north" if lo >= 46 else ("central" if lo >= 44 else "south")
        g = group_of_band[b]
        e, i = ev[:, b], iev[:, b]
        tot_pos = pos[:, b].sum()
        med = nanmed(tre30_w[i, b]) if i.any() else np.nan

        smask = i & (bi_w[:, b] > 3.0)
        speaks = [a + int(np.nanargmax(bi_w[a:c + 1, b])) for a, c in runs_of_true(smask) if c - a + 1 >= 5]
        SLIP_PEAKS[dname][b] = speaks
        s_conf = []
        for (a, c), k in zip([r for r in runs_of_true(smask) if r[1] - r[0] + 1 >= 5], speaks):
            tcnt = tre30_w[k, b]
            cf = bool((tcnt >= 3 * med) and (tcnt >= 10))
            s_conf.append(cf)
            rows_bursts.append(dict(definition=dname, band=BAND_LABELS[b], seg=seg, group=g,
                                    kind="slip_burst", start=str(wdays[a].date()), end=str(wdays[c].date()),
                                    peak_day=str(wdays[k].date()), dur_d=c - a + 1,
                                    peak_slip30_mm=float(bi_w[k, b]), tremor30_at_peak=float(tcnt),
                                    band_iev_median_tremor30=float(med),
                                    ratio_to_median=float(tcnt / med) if med > 0 else np.nan, confirmed=cf))
        obs_pack[g]["s_hits"] += int(np.sum(s_conf))
        obs_pack[g]["s_n"] += len(s_conf)

        tmask = i & (tre30_w[:, b] >= 3 * med) & (tre30_w[:, b] >= 30.0)
        t_conf = []
        for a, c in runs_of_true(tmask):
            if c - a + 1 < 5:
                continue
            k = a + int(np.nanargmax(tre30_w[a:c + 1, b]))
            cf = bool(bi_w[k, b] > 3.0)
            t_conf.append(cf)
            rows_bursts.append(dict(definition=dname, band=BAND_LABELS[b], seg=seg, group=g,
                                    kind="tremor_burst", start=str(wdays[a].date()), end=str(wdays[c].date()),
                                    peak_day=str(wdays[k].date()), dur_d=c - a + 1,
                                    peak_slip30_mm=float(bi_w[k, b]), tremor30_at_peak=float(tre30_w[k, b]),
                                    band_iev_median_tremor30=float(med),
                                    ratio_to_median=float(tre30_w[k, b] / med) if med > 0 else np.nan,
                                    confirmed=cf))
        obs_pack[g]["t_hits"] += int(np.sum(t_conf))
        obs_pack[g]["t_n"] += len(t_conf)

        rows_bands.append(dict(
            definition=dname, band=BAND_LABELS[b], seg=seg, group=g, n_tri=int((band_of_tri == b).sum()),
            n_event_days=int(e.sum()), n_iev_days=int(i.sum()), frac_event_days=float(e.mean()),
            frac_pos_slip_interevent=pos[i, b].sum() / tot_pos if tot_pos > 0 else np.nan,
            mean_pos_rate_event_mmd=pos[e, b].mean() if e.any() else np.nan,
            mean_pos_rate_iev_mmd=pos[i, b].mean() if i.any() else np.nan,
            ratio_event_over_iev=(pos[e, b].mean() / pos[i, b].mean())
            if (e.any() and i.any() and pos[i, b].mean() > 0) else np.nan,
            n_tremor=float(cnt_w[:, b].sum()),
            frac_tremor_interevent=cnt_w[i, b].sum() / cnt_w[:, b].sum() if cnt_w[:, b].sum() > 0 else np.nan,
            iev_median_tremor30=med,
            n_slip_bursts=len(s_conf), slip_burst_conf_rate=float(np.mean(s_conf)) if s_conf else np.nan,
            n_tremor_bursts=len(t_conf), tremor_burst_conf_rate=float(np.mean(t_conf)) if t_conf else np.nan))
    for g in groups:
        p = obs_pack[g]
        OBS[dname][g] = dict(slip=(p["s_hits"] / p["s_n"] if p["s_n"] else np.nan, p["s_n"]),
                             tremor=(p["t_hits"] / p["t_n"] if p["t_n"] else np.nan, p["t_n"]))

bands_df = pd.DataFrame(rows_bands)
bursts_df = pd.DataFrame(rows_bursts)

# --- one shift-null loop shared by all definitions --------------------------------------------
null = {d: {g: dict(slip=[], tremor=[]) for g in groups} for d in DEF_NAMES}
for draw in range(N_NULL_A):
    shift = int(RNG_A.integers(365, T - 365))
    t30s = np.column_stack([roll30(np.roll(cnt_day[:, b], shift))[wi] for b in range(10)])
    for dname in DEF_NAMES:
        iev = IEV[dname]
        for gname, gb in groups.items():
            s_h = s_n = t_h = t_n = 0
            for b in gb:
                i = iev[:, b]
                if not i.any():
                    continue
                med = nanmed(t30s[i, b])
                for k in SLIP_PEAKS[dname][b]:
                    s_n += 1
                    if (t30s[k, b] >= 3 * med) and (t30s[k, b] >= 10):
                        s_h += 1
                tmask = i & (t30s[:, b] >= 3 * med) & (t30s[:, b] >= 30.0)
                for a, c in runs_of_true(tmask):
                    if c - a + 1 < 5:
                        continue
                    kk = a + int(np.nanargmax(t30s[a:c + 1, b]))
                    t_n += 1
                    if bi_w[kk, b] > 3.0:
                        t_h += 1
            null[dname][gname]["slip"].append(s_h / s_n if s_n else np.nan)
            null[dname][gname]["tremor"].append(t_h / t_n if t_n else np.nan)

conf_rows = []
for dname in DEF_NAMES:
    for gname in groups:
        for kind in ("slip", "tremor"):
            o, n = OBS[dname][gname][kind]
            nl = np.asarray(null[dname][gname][kind], float)
            nl = nl[np.isfinite(nl)]
            p = (1 + np.sum(nl >= o)) / (1 + len(nl)) if (np.isfinite(o) and len(nl)) else np.nan
            conf_rows.append(dict(definition=dname, group=gname, burst_kind=kind + "_burst", n_bursts=n,
                                  obs_conf_rate=o,
                                  null_mean=np.nanmean(nl) if len(nl) else np.nan,
                                  null_p95=np.nanpercentile(nl, 95) if len(nl) else np.nan, p=p))
conf_df = pd.DataFrame(conf_rows)

for _, r in conf_df.iterrows():
    col = "slip_burst" if r.burst_kind.startswith("slip") else "tremor_burst"
    m = (bands_df.definition == r.definition) & (bands_df.group == r.group)
    bands_df.loc[m, col + "_null_mean"] = r.null_mean
    bands_df.loc[m, col + "_group_p"] = r.p

bands_df.to_csv(os.path.join(RES, "s6_interets_bands.csv"), index=False)
bursts_df.to_csv(os.path.join(RES, "s6_interets_bursts.csv"), index=False)
conf_df.to_csv(os.path.join(RES, "s6_interets_conf.csv"), index=False)

for dname in DEF_NAMES:
    print(f"\n--- definition {dname}: {DEF_DESC[dname]} ---")
    sub = bands_df[bands_df.definition == dname]
    print(sub[["band", "seg", "n_event_days", "n_iev_days", "frac_pos_slip_interevent",
               "mean_pos_rate_event_mmd", "mean_pos_rate_iev_mmd", "ratio_event_over_iev",
               "frac_tremor_interevent", "iev_median_tremor30", "n_slip_bursts", "slip_burst_conf_rate",
               "n_tremor_bursts", "tremor_burst_conf_rate"]].to_string(index=False,
                                                                       float_format=lambda v: f"{v:.3f}"))
    nc = sub[sub.group == "north/central"]
    so = sub[sub.group == "south"]
    print(f"    inter-event share of POSITIVE SLIP: north/central {nc.frac_pos_slip_interevent.mean():.3f}, "
          f"south {so.frac_pos_slip_interevent.mean():.3f}")
    print(f"    inter-event share of TREMOR       : north/central {nc.frac_tremor_interevent.mean():.3f}, "
          f"south {so.frac_tremor_interevent.mean():.3f}")
    print(f"    event/inter-event mean positive rate ratio: north/central {nc.ratio_event_over_iev.mean():.2f}, "
          f"south {so.ratio_event_over_iev.mean():.2f}")

print("\n--- confirmation rates vs rigid tremor-shift null (200 draws) ---")
print(conf_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

print(f"\ninter-event bursts written: {len(bursts_df)} rows")
print(bursts_df.groupby(["definition", "kind"]).agg(n=("confirmed", "size"),
                                                    n_conf=("confirmed", "sum")).to_string())
for dname in DEF_NAMES:
    s = bursts_df[(bursts_df.definition == dname) & (bursts_df.kind == "slip_burst") &
                  (bursts_df.group == "north/central")]
    if len(s):
        print(f"\n{dname}: north/central inter-event SLIP bursts (all {len(s)}):")
        print(s[["band", "start", "end", "peak_day", "dur_d", "peak_slip30_mm", "tremor30_at_peak",
                 "ratio_to_median", "confirmed"]].to_string(index=False, float_format=lambda v: f"{v:.2f}"))

# ==================================================================================================
# PART B - the southern +-5 mm oscillation   (UNCHANGED from the first run)
# ==================================================================================================
hdr("PART B - the southern +-5 mm oscillation: seasonal? a few components? noise?")

# Part B keeps its own RNG stream, reproduced bit-for-bit from the first run (which consumed
# 2 groups x 200 scalar shift draws in Part A before reaching Part B) so the already-documented
# Part B numbers are unchanged by the Part A rewrite.
RNG = np.random.default_rng(20260910)
for _ in range(2 * N_NULL_A):
    RNG.integers(365, T - 365)

south_tri = plat < 44.0
north_tri = (plat >= 46.0) & (plat < 49.5)
print(f"south triangles (lat<44): {south_tri.sum()};  north triangles (46-49.5): {north_tri.sum()}")

# --- (1) sign symmetry per segment per year ---------------------------------------------------
inc_w = inc30[wi, :]
yrs = wdays.year.to_numpy()
sym_rows = []
for segname, sel in (("south(<44)", south_tri), ("north(46-49.5)", north_tri)):
    sub = inc_w[:, sel]
    for y in np.unique(yrs):
        a = sub[yrs == y, :]
        n = np.isfinite(a).sum()
        sym_rows.append(dict(segment=segname, year=int(y), n_triangle_days=int(n),
                             frac_gt5=float(np.sum(a > 5) / n) if n else np.nan,
                             frac_lt_m5=float(np.sum(a < -5) / n) if n else np.nan,
                             frac_abs_gt5=float(np.sum(np.abs(a) > 5) / n) if n else np.nan,
                             mean_mm=float(np.nanmean(a)) if n else np.nan))
sym_df = pd.DataFrame(sym_rows)
sym_df.to_csv(os.path.join(RES, "s6_sign_symmetry_year.csv"), index=False)
print("\n--- (1) sign symmetry: fraction of triangle-days with inc30 > +5 / < -5 mm, per year ---")
print(sym_df.pivot(index="year", columns="segment", values=["frac_gt5", "frac_lt_m5"]).to_string(
    float_format=lambda v: f"{v:.4f}"))
for segname in sym_df.segment.unique():
    s = sym_df[sym_df.segment == segname]
    tot_g = (s.frac_gt5 * s.n_triangle_days).sum() / s.n_triangle_days.sum()
    tot_l = (s.frac_lt_m5 * s.n_triangle_days).sum() / s.n_triangle_days.sum()
    print(f"  {segname:16s} all years: frac>+5 = {tot_g:.4f}   frac<-5 = {tot_l:.4f}   "
          f"ratio pos/neg = {tot_g/tot_l if tot_l>0 else np.nan:.2f}")

# --- (2) spectrum + seasonal cycle -------------------------------------------------------------
bi_win = band_inc30[wi, :]
mo = wdays.month.to_numpy()
spec_rows, welch_store, moy_store = [], {}, {}
for b in range(10):
    x = bi_win[:, b].astype(float)
    x = np.where(np.isfinite(x), x, np.nanmean(x))
    f, P = welch(x, fs=1.0, nperseg=1461, detrend="linear")
    per = np.where(f > 0, 1.0 / np.maximum(f, 1e-12), np.inf)
    welch_store[b] = (per, P)
    sel = (per >= 60) & (per <= 1096)
    ann = (per >= 330) & (per <= 400)
    nb = ((per >= 180) & (per <= 300)) | ((per >= 430) & (per <= 900))
    moy = np.array([x[mo == m].mean() for m in range(1, 13)])
    moy_store[b] = moy
    amp = float(moy.max() - moy.min())
    nulls = np.empty(N_NULL_SEAS)
    for d in range(N_NULL_SEAS):
        xs = np.roll(x, int(RNG.integers(30, 701)))
        m2 = np.array([xs[mo == m].mean() for m in range(1, 13)])
        nulls[d] = m2.max() - m2.min()
    spec_rows.append(dict(band=BAND_LABELS[b],
                          seg="north" if BAND_EDGES[b] >= 46 else ("central" if BAND_EDGES[b] >= 44 else "south"),
                          peak_period_d=float(per[sel][np.argmax(P[sel])]), peak_power=float(P[sel].max()),
                          ann_ratio=float(P[ann].mean() / np.median(P[nb])) if ann.any() and nb.any() else np.nan,
                          seasonal_amp_mm=amp,
                          seasonal_p=(1 + np.sum(nulls >= amp)) / (1 + N_NULL_SEAS),
                          seasonal_null_mean=float(nulls.mean()),
                          moy_peak_month=int(np.argmax(moy) + 1), moy_min_month=int(np.argmin(moy) + 1),
                          std_inc30_mm=float(np.nanstd(x))))
spec_df = pd.DataFrame(spec_rows)
for b in range(10):
    a = inc_w[:, band_of_tri == b]
    fin = np.isfinite(a).sum()
    spec_df.loc[b, "frac_tri_days_gt5"] = float(np.sum(a > 5) / fin) if fin else np.nan
    spec_df.loc[b, "frac_tri_days_lt_m5"] = float(np.sum(a < -5) / fin) if fin else np.nan

# --- (3) coherence ------------------------------------------------------------------------------
C = np.corrcoef(np.nan_to_num(bi_win, nan=0.0).T)
pd.DataFrame(C, index=BAND_LABELS, columns=BAND_LABELS).to_csv(os.path.join(RES, "s6_band_corr.csv"))
sb = [b for b in range(10) if BAND_EDGES[b] < 44]
nb_ = [b for b in range(10) if BAND_EDGES[b] >= 46]


def offdiag_mean(bs):
    v = [C[i, j] for ii, i in enumerate(bs) for j in bs[ii + 1:]]
    return float(np.mean(v)) if v else np.nan


mean_ss, mean_nn = offdiag_mean(sb), offdiag_mean(nb_)
mean_sn = float(np.mean([C[i, j] for i in sb for j in nb_]))

A30 = pd.DataFrame(A_raw).rolling(30, min_periods=30).mean().to_numpy()
sumA30 = np.nansum(A30, axis=1)[wi]
good = np.isfinite(sumA30)
for b in range(10):
    spec_df.loc[b, "corr_sumA30"] = float(np.corrcoef(np.nan_to_num(bi_win[:, b])[good], sumA30[good])[0, 1])
spec_df.to_csv(os.path.join(RES, "s6_south_oscillation_bands.csv"), index=False)

comp_rows = []
for k in range(A_raw.shape[1]):
    ser_cum = -1e3 * np.outer(A_raw[:, k], Ld[:, k])
    inc_k = ser_cum - np.roll(ser_cum, 30, axis=0)
    inc_k[:30, :] = np.nan
    comp_rows.append(dict(component=k,
                          var_cum_south=float(np.nansum(np.nanvar(ser_cum[wi][:, south_tri], axis=0))),
                          var_cum_north=float(np.nansum(np.nanvar(ser_cum[wi][:, north_tri], axis=0))),
                          var_inc30_south=float(np.nansum(np.nanvar(inc_k[wi][:, south_tri], axis=0))),
                          var_inc30_north=float(np.nansum(np.nanvar(inc_k[wi][:, north_tri], axis=0)))))
comp_df = pd.DataFrame(comp_rows)
for c in ["var_cum_south", "var_cum_north", "var_inc30_south", "var_inc30_north"]:
    comp_df[c.replace("var", "frac")] = comp_df[c] / comp_df[c].sum()
comp_df = comp_df.sort_values("var_inc30_south", ascending=False).reset_index(drop=True)
comp_df["rank_south_inc30"] = np.arange(1, len(comp_df) + 1)
comp_df.to_csv(os.path.join(RES, "s6_south_components.csv"), index=False)

epochs = [("2010-2014", 2010, 2014), ("2015-2019", 2015, 2019), ("2020-2026", 2020, 2026)]
stat_rows = []
for name, y0, y1 in epochs:
    m = (yrs >= y0) & (yrs <= y1)
    for segname, sel in (("south(<44)", south_tri), ("north(46-49.5)", north_tri)):
        a = inc_w[m][:, sel]
        fin = np.isfinite(a).sum()
        stat_rows.append(dict(epoch=name, segment=segname,
                              frac_abs_gt5=float(np.sum(np.abs(a) > 5) / fin) if fin else np.nan,
                              frac_gt5=float(np.sum(a > 5) / fin) if fin else np.nan,
                              frac_lt_m5=float(np.sum(a < -5) / fin) if fin else np.nan,
                              std_mm=float(np.nanstd(a))))
stat_df = pd.DataFrame(stat_rows)

print("\n--- (2) spectrum + seasonal cycle per band ---")
print(spec_df[["band", "seg", "std_inc30_mm", "peak_period_d", "ann_ratio", "seasonal_amp_mm",
               "seasonal_null_mean", "seasonal_p", "moy_peak_month", "moy_min_month",
               "frac_tri_days_gt5", "frac_tri_days_lt_m5", "corr_sumA30"]].to_string(
    index=False, float_format=lambda v: f"{v:.3f}"))
print("\n--- (3) coherence between the 10 band-mean inc30 series ---")
print(pd.DataFrame(C, index=BAND_LABELS, columns=BAND_LABELS).to_string(float_format=lambda v: f"{v:.2f}"))
print(f"  mean corr among SOUTH bands (<44)   = {mean_ss:.3f}")
print(f"  mean corr among NORTH bands (>=46)  = {mean_nn:.3f}")
print(f"  mean corr SOUTH x NORTH             = {mean_sn:.3f}")
print("\n  vbICA components ranked by southern inc30 variance:")
print(comp_df[["component", "rank_south_inc30", "var_inc30_south", "frac_inc30_south",
               "frac_inc30_north", "frac_cum_south"]].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
top3 = comp_df.frac_inc30_south.iloc[:3].sum()
print(f"  top-3 components carry {100*top3:.1f}% of the southern inc30 variance "
      f"(components {comp_df.component.iloc[:3].tolist()})")
print("\n--- (4) stationarity of the southern oscillation ---")
print(stat_df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# ==================================================================================================
# FIGURE
# ==================================================================================================
DEF_COL = {"A0": "#8172b3", "A1": "#55a868", "A2": "#c44e52"}
fig = plt.figure(figsize=(17.5, 15))
gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.34)
xb = np.arange(10)
wbar = 0.27

# A1a inter-event share of positive slip
ax = fig.add_subplot(gs[0, 0])
for k, dn in enumerate(DEF_NAMES):
    s = bands_df[bands_df.definition == dn]
    ax.bar(xb + (k - 1) * wbar, s.frac_pos_slip_interevent, width=wbar, color=DEF_COL[dn],
           edgecolor="k", lw=0.3, label=dn)
ax.axvline(3.5, color="k", ls=":", lw=1)
ax.set_xticks(xb); ax.set_xticklabels(BAND_LABELS, rotation=60, fontsize=8)
ax.set_ylabel("fraction of POSITIVE slip\non inter-event days"); ax.set_ylim(0, 1.08)
ax.set_title("A1a  inter-ETS share of modeled positive slip\n(dotted line: south | north+central)", fontsize=10)
ax.legend(fontsize=8, ncol=3); ax.grid(alpha=0.3, axis="y")

# A1b inter-event share of tremor
ax = fig.add_subplot(gs[0, 1])
for k, dn in enumerate(DEF_NAMES):
    s = bands_df[bands_df.definition == dn]
    ax.bar(xb + (k - 1) * wbar, s.frac_tremor_interevent, width=wbar, color=DEF_COL[dn],
           edgecolor="k", lw=0.3, label=dn)
ax.axvline(3.5, color="k", ls=":", lw=1)
ax.set_xticks(xb); ax.set_xticklabels(BAND_LABELS, rotation=60, fontsize=8)
ax.set_ylabel("fraction of tremor detections\non inter-event days"); ax.set_ylim(0, 1.08)
ax.set_title("A1b  inter-ETS share of tremor", fontsize=10)
ax.legend(fontsize=8, ncol=3); ax.grid(alpha=0.3, axis="y")

# A1c event/inter-event rate ratio
ax = fig.add_subplot(gs[0, 2])
for k, dn in enumerate(DEF_NAMES):
    s = bands_df[bands_df.definition == dn]
    ax.bar(xb + (k - 1) * wbar, s.ratio_event_over_iev, width=wbar, color=DEF_COL[dn],
           edgecolor="k", lw=0.3, label=dn)
ax.axhline(1.0, color="k", lw=0.8)
ax.set_xticks(xb); ax.set_xticklabels(BAND_LABELS, rotation=60, fontsize=8)
ax.set_ylabel("mean positive rate,\nevent / inter-event")
ax.set_title("A1c  is the model's slip concentrated in ETS?", fontsize=10)
ax.legend(fontsize=8, ncol=3); ax.grid(alpha=0.3, axis="y")

# A2a / A2b confirmation vs null
for pi, kind in enumerate(["slip_burst", "tremor_burst"]):
    ax = fig.add_subplot(gs[1, pi])
    labs, xs = [], []
    for gi, gname in enumerate(["north/central", "south"]):
        for k, dn in enumerate(DEF_NAMES):
            r = conf_df[(conf_df.definition == dn) & (conf_df.group == gname)
                        & (conf_df.burst_kind == kind)].iloc[0]
            x = gi * 4 + k
            xs.append(x)
            labs.append(f"{dn}\nn={int(r.n_bursts)}")
            ax.bar(x, r.obs_conf_rate, color=DEF_COL[dn], edgecolor="k", lw=0.4)
            ax.plot([x - 0.4, x + 0.4], [r.null_mean] * 2, "k-", lw=2)
            ax.plot([x - 0.4, x + 0.4], [r.null_p95] * 2, "r-", lw=2)
            if np.isfinite(r.obs_conf_rate):
                ax.text(x, r.obs_conf_rate + 0.03, f"p={r.p:.3f}", ha="center", fontsize=7, rotation=90)
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=7)
    ax.text(1, -0.22, "north/central", ha="center", transform=ax.get_xaxis_transform(), fontsize=9)
    ax.text(5, -0.22, "south", ha="center", transform=ax.get_xaxis_transform(), fontsize=9)
    ax.set_ylim(0, 1.25); ax.set_ylabel("confirmation rate")
    ttl = ("A2a  inter-event SLIP bursts confirmed by tremor"
           if kind == "slip_burst" else "A2b  inter-event TREMOR bursts with slip > 3 mm")
    ax.set_title(ttl + "\n(black = shift-null mean, red = null p95)", fontsize=10)
    ax.grid(alpha=0.3, axis="y")

# B2a spectra
ax = fig.add_subplot(gs[1, 2])
for b in range(10):
    per, P = welch_store[b]
    sel = (per >= 40) & (per <= 3000)
    isS = BAND_EDGES[b] < 44
    ax.loglog(per[sel], P[sel], color=("#c44e52" if isS else "#4c72b0"), alpha=0.9 if isS else 0.35,
              lw=1.6 if isS else 1.0, label=BAND_LABELS[b] if isS else None)
ax.axvline(365.25, color="k", ls="--", lw=1, label="1 yr")
ax.axvline(182.6, color="gray", ls=":", lw=1, label="6 mo")
ax.set_xlabel("period (days)"); ax.set_ylabel("power (mm$^2$/cpd)")
ax.set_title("B2a  Welch spectra of band-mean inc30\n(red = south <44N, blue = 44-50N)", fontsize=10)
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3, which="both")

# B2b month-of-year
ax = fig.add_subplot(gs[2, 0])
mm = np.arange(1, 13)
for b in range(10):
    isS = BAND_EDGES[b] < 44
    ax.plot(mm, moy_store[b], "-o", ms=3, color=("#c44e52" if isS else "#4c72b0"),
            alpha=0.9 if isS else 0.35, lw=1.8 if isS else 1.0, label=BAND_LABELS[b] if isS else None)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(mm); ax.set_xlabel("month of year"); ax.set_ylabel("mean inc30 (mm/30 d)")
ax.set_title("B2b  month-of-year cycle, south (red) vs north (blue)", fontsize=10)
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)

# B3 correlation matrix
ax = fig.add_subplot(gs[2, 1])
im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(xb); ax.set_xticklabels(BAND_LABELS, rotation=90, fontsize=7)
ax.set_yticks(xb); ax.set_yticklabels(BAND_LABELS, fontsize=7)
for i in range(10):
    for j in range(10):
        ax.text(j, i, f"{C[i,j]:.2f}", ha="center", va="center", fontsize=5.5,
                color="k" if abs(C[i, j]) < 0.6 else "w")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
ax.set_title(f"B3  band-mean inc30 correlation\nS-S {mean_ss:.2f} | N-N {mean_nn:.2f} | S-N {mean_sn:.2f}",
             fontsize=10)

# B4 stationarity
ax = fig.add_subplot(gs[2, 2])
w = 0.35
for k, (segname, c) in enumerate([("south(<44)", "#c44e52"), ("north(46-49.5)", "#4c72b0")]):
    s = stat_df[stat_df.segment == segname]
    ax.bar(np.arange(3) + (k - 0.5) * w, s.frac_abs_gt5, width=w, color=c, edgecolor="k", lw=0.4, label=segname)
ax.set_xticks(np.arange(3)); ax.set_xticklabels([e[0] for e in epochs])
ax.set_ylabel("fraction of triangle-days |inc30| > 5 mm", labelpad=2)
ax.set_title("B4  stationarity of the >5 mm excursions\n"
             f"(top-3 vbICA comps = {100*top3:.0f}% of southern inc30 var)", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3, axis="y")

fig.suptitle("S6  model-only checks: inter-ETS slip under 3 event-day definitions (A) "
             "and the southern +-5 mm oscillation (B)", fontsize=13, y=0.955)
fig.savefig(os.path.join(FIG, "s6_model_checks.png"), dpi=130, bbox_inches="tight")
plt.close(fig)

hdr("OUTPUTS")
for f in ["results/s6_interets_bands.csv", "results/s6_interets_bursts.csv", "results/s6_interets_conf.csv",
          "results/s6_south_oscillation_bands.csv", "results/s6_sign_symmetry_year.csv",
          "results/s6_south_components.csv", "results/s6_band_corr.csv", "figures/s6_model_checks.png"]:
    print(f"  {f:46s} {os.path.getsize(os.path.join(HERE, f)):>9,d} bytes")

hdr("HEADLINES")
for dname in DEF_NAMES:
    s = bands_df[bands_df.definition == dname]
    nc, so = s[s.group == "north/central"], s[s.group == "south"]
    print(f"{dname}  event-day coverage {100*s.frac_event_days.mean():.0f}%  |  "
          f"inter-event positive slip N/C {nc.frac_pos_slip_interevent.mean():.3f} S {so.frac_pos_slip_interevent.mean():.3f}"
          f"  |  inter-event tremor N/C {nc.frac_tremor_interevent.mean():.3f} S {so.frac_tremor_interevent.mean():.3f}"
          f"  |  rate ratio N/C {nc.ratio_event_over_iev.mean():.2f} S {so.ratio_event_over_iev.mean():.2f}")
for _, r in conf_df.iterrows():
    print(f"{r.definition}  {r.group:14s} {r.burst_kind:13s} n={int(r.n_bursts):4d} obs={r.obs_conf_rate:.3f} "
          f"null={r.null_mean:.3f} (p95 {r.null_p95:.3f}) p={r.p:.3f}")
print(f"southern seasonal amplitude (mean over 4 southern bands): "
      f"{spec_df[spec_df.seg=='south'].seasonal_amp_mm.mean():.2f} mm, "
      f"p = {[round(v,3) for v in spec_df[spec_df.seg=='south'].seasonal_p]}")
print(f"southern band-mean peak periods (d): {spec_df[spec_df.seg=='south'].peak_period_d.round(0).tolist()}")
print(f"top-3 components = {100*top3:.1f}% of southern inc30 variance")
print("done.")
