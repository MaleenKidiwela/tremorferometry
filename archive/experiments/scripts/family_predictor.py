#!/usr/bin/env python
"""
family_predictor.py
===================
Can we PREDICT, from cheap pre-densify signatures, which LFE families will be
continuously-repeating / CWI-suitable ("GOOD") vs gappy/bursty ("BAD") -- so we
can pre-screen the eligible pool before the expensive full densify?

We use families already processed (labels in data/family_quality_flags.csv) as
labeled training data. GOOD=1 vs BAD=0; MARGINAL dropped.

WHAT'S ALREADY KNOWN (not re-litigated here, only carried as baseline features):
  detection-COUNT/time-stat features + template SNR were tested. Pooled AUC looks
  predictive but it is a STATION-LEVEL CONFOUND: per-station AUC ~= 0.5 for all of
  them; only template SNR carries weak within-station signal (~0.30, inverted).

THIS SCRIPT tests the UNTESTED, physically-motivated feature families:
  1. TEMPLATE WAVEFORM SHAPE  (data/<sta>_pnsn_families_100km.npz)
     fs=40 Hz, 80-sample 2 s window, bandpass 2-8 Hz, L2-normalized,
     S-peak pinned at the center sample (template-pre=-1.0, post=+1.0).
  2. CATALOG / LOCATION       (catalogs/pnsn_tremor_cascadia_full.csv,
     data/station_slab2_depth.csv + hardcoded station coords)
  3. Combine ALL (ours + the existing cheap ones) -> logistic + HistGB.

CRITICAL METHODOLOGY:
  Cross-validate with GroupKFold grouped by STATION (whole stations held out).
  Report GROUPED CV AUC, never pooled. Compare to a within-station baseline.
  Report per-feature grouped AUC, combined-model grouped AUC, permutation importance.
  Be adversarial: grouped AUC ~= 0.5 => NOT predictable => must measure via probe.

Outputs (NEW files only):
  data/family_predictor_full_features.csv
  figures/family_predictor_auc.png
  notes/FAMILY_PREDICTOR.md   (written by hand from the printed report)
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import glob
import warnings
import numpy as np
import pandas as pd
from scipy.signal import hilbert
from scipy.stats import kurtosis, skew

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GroupKFold
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve

warnings.filterwarnings("ignore")
RNG = np.random.RandomState(0)

ROOT = "/home/jovyan/tremorferometry"
FS = 40.0          # Hz
NS = 80            # template length (samples) = 2.0 s
CENTER = NS // 2   # S-peak pinned here (pre=-1.0 s, post=+1.0 s)
FMIN, FMAX = 2.0, 8.0

# Station coords (from scripts/build_dvv_map.py STA dict).
STA_COORDS = {
 'B927': (49.2188, -124.8113), 'NLLB': (49.2271, -123.9882), 'B928': (48.834, -125.134),
 'PGC': (48.6498, -123.4521), 'B011': (48.65, -123.448), 'B004': (48.202, -124.427),
 'B013': (47.813, -122.9108), 'HDW': (47.6490, -123.0530), 'GNW': (47.5641, -122.8250),
 'B014': (47.5133, -123.8125), 'B941': (46.9868, -122.219), 'B018': (46.9795, -123.0203),
 'B020': (46.3827, -123.8445), 'B201': (46.3033, -122.2648), 'B204': (46.136, -122.169),
 'B023': (46.1112, -123.0787), 'B022': (45.9546, -123.931), 'B026': (45.3094, -123.8231),
 'COLT': (45.17044, -122.438152), 'COR': (44.5855, -123.3046), 'B028': (44.4937, -122.9638),
 'B030': (43.9713, -122.7717), 'B032': (43.668, -123.3923), 'B033': (43.2917, -123.1245),
 'B036': (42.5058, -123.3817), 'B040': (41.8308, -122.4205), 'B039': (41.4667, -122.4847),
 'B935': (40.4787, -123.5732),
}


# ----------------------------------------------------------------------------
# Template-shape features (one waveform -> dict of scalars)
# ----------------------------------------------------------------------------
def template_features(w):
    """Per-template waveform-shape features. w: float array, L2-normalized, len 80,
    S-peak at center sample. Returns dict of scalars (NaN-safe)."""
    w = np.asarray(w, dtype=np.float64)
    n = len(w)
    out = {}
    if n < 8 or not np.all(np.isfinite(w)) or np.allclose(w, 0):
        return {k: np.nan for k in TEMPLATE_FEATS}

    rms = np.sqrt(np.mean(w ** 2))
    peak = np.max(np.abs(w))
    # crest factor: peak / RMS  (impulsiveness; high => spiky, low => emergent)
    out["t_crest"] = peak / rms if rms > 0 else np.nan
    # waveform-amplitude distribution shape
    out["t_kurtosis"] = kurtosis(w, fisher=True, bias=False)
    out["t_skew"] = skew(w, bias=False)

    # ---- envelope (analytic signal) ----
    env = np.abs(hilbert(w))
    env = np.nan_to_num(env)
    epk = env.max()
    ipk = int(np.argmax(env))
    out["t_env_crest"] = epk / np.sqrt(np.mean(env ** 2)) if env.mean() > 0 else np.nan
    # impulsiveness: duration above half-max envelope (seconds)
    if epk > 0:
        above = env >= 0.5 * epk
        out["t_dur_halfmax_s"] = above.sum() / FS
        # rise time: from first sample >=10% of peak (before peak) to the peak
        thr = 0.10 * epk
        pre = np.where(env[:ipk + 1] >= thr)[0]
        out["t_rise_s"] = (ipk - pre[0]) / FS if len(pre) else np.nan
        # peak position relative to pinned center (samples). Misalignment / which
        # phase actually dominates (S-pinned families should sit near 0).
        out["t_peakpos_s"] = (ipk - CENTER) / FS
    else:
        out["t_dur_halfmax_s"] = np.nan
        out["t_rise_s"] = np.nan
        out["t_peakpos_s"] = np.nan

    # envelope decay rate: slope of log-envelope post-peak (1/s); steeper = faster coda decay
    post = env[ipk:]
    post = post[post > 1e-6 * epk] if epk > 0 else post
    if len(post) >= 4:
        t = np.arange(len(post)) / FS
        y = np.log(post)
        # robust-ish linear fit
        A = np.vstack([t, np.ones_like(t)]).T
        slope = np.linalg.lstsq(A, y, rcond=None)[0][0]
        out["t_decay_per_s"] = -slope          # positive = decaying
    else:
        out["t_decay_per_s"] = np.nan

    # energy concentration: fraction of total energy in central +/-0.25 s (10 samples)
    half = int(round(0.25 * FS))
    lo, hi = max(0, CENTER - half), min(n, CENTER + half)
    out["t_energy_conc"] = float(np.sum(w[lo:hi] ** 2) / np.sum(w ** 2))

    # self-SNR: peak envelope / pre-arrival RMS (pre = first 0.75 s, well before center)
    pre_n = int(round(0.75 * FS))
    pre_rms = np.sqrt(np.mean(w[:pre_n] ** 2)) if pre_n > 1 else np.nan
    out["t_self_snr"] = (epk / pre_rms) if (pre_rms and pre_rms > 0) else np.nan

    # zero-crossing rate (per second) -- proxy for frequency content / roughness
    zc = np.sum(np.abs(np.diff(np.sign(w))) > 0)
    out["t_zcr_hz"] = zc / (n / FS)

    # autocorrelation half-width (lag in s where |autocorr| first drops below 0.5)
    ac = np.correlate(w, w, mode="full")[n - 1:]
    if ac[0] != 0:
        ac = ac / ac[0]
        below = np.where(ac < 0.5)[0]
        out["t_ac_halfwidth_s"] = (below[0] / FS) if len(below) else (n / FS)
    else:
        out["t_ac_halfwidth_s"] = np.nan

    # ---- spectrum (rfft of the windowed template) ----
    win = np.hanning(n)
    sp = np.abs(np.fft.rfft(w * win))
    freqs = np.fft.rfftfreq(n, d=1.0 / FS)
    psd = sp ** 2
    band = (freqs >= FMIN) & (freqs <= FMAX)
    if psd[band].sum() > 0:
        fb = freqs[band]
        pb = psd[band]
        pn = pb / pb.sum()
        out["t_dom_freq"] = fb[int(np.argmax(pb))]
        cen = float(np.sum(fb * pn))
        out["t_spec_centroid"] = cen
        out["t_spec_bw"] = float(np.sqrt(np.sum(((fb - cen) ** 2) * pn)))
        # spectral entropy (normalized 0..1 over band bins)
        nz = pn[pn > 0]
        out["t_spec_entropy"] = float(-np.sum(nz * np.log(nz)) / np.log(len(pn)))
        # spectral flatness (geo mean / arith mean) over band
        out["t_spec_flatness"] = float(np.exp(np.mean(np.log(pb + 1e-20))) / np.mean(pb))
    else:
        for k in ("t_dom_freq", "t_spec_centroid", "t_spec_bw",
                  "t_spec_entropy", "t_spec_flatness"):
            out[k] = np.nan
    return out


TEMPLATE_FEATS = [
    "t_crest", "t_kurtosis", "t_skew", "t_env_crest", "t_dur_halfmax_s",
    "t_rise_s", "t_peakpos_s", "t_decay_per_s", "t_energy_conc", "t_self_snr",
    "t_zcr_hz", "t_ac_halfwidth_s", "t_dom_freq", "t_spec_centroid",
    "t_spec_bw", "t_spec_entropy", "t_spec_flatness",
]


# ----------------------------------------------------------------------------
# Catalog / location features
# ----------------------------------------------------------------------------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def azimuth_deg(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dlon = np.radians(lon2 - lon1)
    x = np.sin(dlon) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360.0) % 360.0


def parse_famll(patch):
    a = patch.split("__")[0].split("_")
    return float(a[0]), float(a[1])


def build_catalog_cell_stats(cat, grid=0.05):
    """Per-0.05deg-cell tremor recurrence stats from the full PNSN catalog.
    Returns dict (ilat,ilon) -> dict of stats. Cells keyed by floored grid index."""
    cat = cat.dropna(subset=["lat", "lon", "time"]).copy()
    cat["time"] = pd.to_datetime(cat["time"], errors="coerce")
    cat = cat.dropna(subset=["time"])
    cat["ilat"] = np.round(cat["lat"] / grid).astype(int)
    cat["ilon"] = np.round(cat["lon"] / grid).astype(int)
    # global month index for active-month-fraction
    cat["ym"] = cat["time"].dt.year * 12 + cat["time"].dt.month
    ym_span_total = cat["ym"].max() - cat["ym"].min() + 1

    stats = {}
    for (ila, ilo), g in cat.groupby(["ilat", "ilon"]):
        ts = np.sort(g["time"].values.astype("datetime64[s]").astype(np.int64))
        n = len(ts)
        nmonths = g["ym"].nunique()
        # inter-event gaps (days) regularity -> CV of gaps; episodic if regular ~13-14 mo? use overall
        if n >= 3:
            gaps = np.diff(ts) / 86400.0
            gaps = gaps[gaps > 0]
            cv = (gaps.std() / gaps.mean()) if (len(gaps) and gaps.mean() > 0) else np.nan
        else:
            cv = np.nan
        stats[(ila, ilo)] = dict(
            cat_count=float(n),
            cat_log_count=float(np.log10(n + 1)),
            cat_active_mo_frac=float(nmonths / ym_span_total),
            cat_inter_cv=float(cv),
        )
    CAT_DEFAULT = dict(cat_count=0.0, cat_log_count=0.0,
                       cat_active_mo_frac=0.0, cat_inter_cv=np.nan)
    return stats, grid, CAT_DEFAULT


CATALOG_FEATS = ["cat_count", "cat_log_count", "cat_active_mo_frac", "cat_inter_cv"]
LOCATION_FEATS = ["loc_dist_km", "loc_az_deg", "loc_sta_slab_km",
                  "loc_fam_lat", "loc_fam_lon"]

# existing cheap features (carried as baseline; from data/family_predictor_features.csv)
CHEAP_FEATS = ["n_det", "span_yr", "n_yr", "n_mo", "frac_active_yr",
               "monthly_cv", "det_per_activemo", "snr"]


# ----------------------------------------------------------------------------
# Main feature-table build
# ----------------------------------------------------------------------------
def build_features():
    q = pd.read_csv(f"{ROOT}/data/family_quality_flags.csv")
    q = q[q["quality"].isin(["GOOD", "BAD"])].copy()
    q["good"] = (q["quality"] == "GOOD").astype(int)

    # ---- station slab2 depth ----
    sd = pd.read_csv(f"{ROOT}/data/station_slab2_depth.csv")
    sd_map = dict(zip(sd["station"].str.upper(), sd["slab2_depth_km"]))

    # ---- catalog cell stats ----
    print("[cat] loading full PNSN tremor catalog ...", flush=True)
    cat = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv",
                      usecols=["time", "lat", "lon"])
    cat_stats, grid, CAT_DEFAULT = build_catalog_cell_stats(cat)
    print(f"[cat] {len(cat_stats)} populated 0.05deg cells", flush=True)

    # ---- existing cheap features ----
    cheap = pd.read_csv(f"{ROOT}/data/family_predictor_features.csv")
    cheap_key = cheap.set_index(["station", "family_id"])

    rows = []
    for sta, grp in q.groupby("station"):
        npz_path = f"{ROOT}/data/{sta.lower()}_pnsn_families_100km.npz"
        templates = {}
        if os.path.exists(npz_path):
            d = np.load(npz_path, allow_pickle=True)
            keys = set(d.files)
        else:
            d, keys = None, set()
            print(f"[tmpl] {sta}: NO NPZ -> template features = NaN "
                  f"({len(grp)} families)", flush=True)

        slat, slon = STA_COORDS[sta]
        sta_slab = sd_map.get(sta.upper(), np.nan)

        nt_found = 0
        for _, r in grp.iterrows():
            patch = r["patch"]
            row = dict(station=sta, family_id=patch, quality=r["quality"],
                       good=int(r["good"]))

            # template features
            if d is not None and patch in keys:
                row.update(template_features(d[patch]))
                nt_found += 1
            else:
                row.update({k: np.nan for k in TEMPLATE_FEATS})

            # location features
            flat, flon = parse_famll(patch)
            row["loc_fam_lat"] = flat
            row["loc_fam_lon"] = flon
            row["loc_dist_km"] = haversine_km(slat, slon, flat, flon)
            row["loc_az_deg"] = azimuth_deg(slat, slon, flat, flon)
            row["loc_sta_slab_km"] = sta_slab

            # catalog-cell features
            cs = cat_stats.get((int(round(flat / grid)), int(round(flon / grid))),
                               CAT_DEFAULT)
            row.update(cs)

            # existing cheap features
            try:
                cr = cheap_key.loc[(sta, patch)]
                for f in CHEAP_FEATS:
                    row[f] = float(cr[f])
            except KeyError:
                for f in CHEAP_FEATS:
                    row[f] = np.nan

            rows.append(row)
        if d is not None:
            print(f"[tmpl] {sta}: {nt_found}/{len(grp)} templates", flush=True)

    df = pd.DataFrame(rows)
    return df


# ----------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------
def grouped_cv_auc_single(x, y, groups, n_splits):
    """Grouped-CV AUC for a SINGLE feature (out-of-fold pooled AUC, sign-agnostic).
    Returns max(auc, 1-auc) so a perfectly-inverted predictor reads as informative,
    and also the raw (signed) auc."""
    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    xx = x.astype(float).copy()
    for tr, te in gkf.split(xx, y, groups):
        # impute test NaN with TRAIN median (no leakage)
        med = np.nanmedian(xx[tr])
        col = xx[te].copy()
        col[~np.isfinite(col)] = med if np.isfinite(med) else 0.0
        oof[te] = col
    m = np.isfinite(oof)
    if m.sum() < 10 or len(np.unique(y[m])) < 2:
        return np.nan, np.nan
    auc_raw = roc_auc_score(y[m], oof[m])
    return max(auc_raw, 1 - auc_raw), auc_raw


def within_station_auc_single(x, y, groups):
    """Average within-station AUC for a single feature (sign-agnostic), weighted by
    number of usable families per station. Only stations with both classes count."""
    aucs, ws = [], []
    for g in np.unique(groups):
        m = (groups == g) & np.isfinite(x)
        if m.sum() < 8:
            continue
        yy = y[m]
        if len(np.unique(yy)) < 2:
            continue
        a = roc_auc_score(yy, x[m])
        aucs.append(max(a, 1 - a))
        ws.append(m.sum())
    if not aucs:
        return np.nan, 0
    return float(np.average(aucs, weights=ws)), len(aucs)


def grouped_cv_model(X, y, groups, model_fn, n_splits):
    """Out-of-fold predicted probabilities under GroupKFold; returns oof, pooled AUC."""
    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    for tr, te in gkf.split(X, y, groups):
        mdl = model_fn()
        mdl.fit(X[tr], y[tr])
        oof[te] = mdl.predict_proba(X[te])[:, 1]
    m = np.isfinite(oof)
    return oof, roc_auc_score(y[m], oof[m])


def main():
    feats_csv = f"{ROOT}/data/family_predictor_full_features.csv"
    df = build_features()
    df.to_csv(feats_csv, index=False)
    print(f"\n[out] wrote {feats_csv}  ({len(df)} rows)\n", flush=True)

    ALL_FEATS = (TEMPLATE_FEATS + CATALOG_FEATS + LOCATION_FEATS + CHEAP_FEATS)
    y = df["good"].values.astype(int)
    groups = df["station"].values
    n_stations = len(np.unique(groups))
    n_splits = min(n_stations, 10)

    print(f"=== DATA ===  n={len(df)}  GOOD={y.sum()}  BAD={(y==0).sum()}  "
          f"stations={n_stations}  base_rate_GOOD={y.mean():.3f}", flush=True)
    print(f"GroupKFold splits (by station) = {n_splits}\n", flush=True)

    # ---- per-feature grouped + within-station AUC ----
    print("=== PER-FEATURE AUC ===")
    print(f"{'feature':22s} {'grpCV':>7s} {'(raw)':>7s} {'within':>7s} "
          f"{'#sta':>5s} {'nfin':>6s}")
    perfeat = []
    for f in ALL_FEATS:
        x = df[f].values.astype(float)
        nfin = int(np.isfinite(x).sum())
        g_abs, g_raw = grouped_cv_auc_single(x, y, groups, n_splits)
        w_abs, nsta = within_station_auc_single(x, y, groups)
        perfeat.append(dict(feature=f, grpcv_auc=g_abs, grpcv_raw=g_raw,
                            within_auc=w_abs, n_sta=nsta, n_finite=nfin))
        print(f"{f:22s} {g_abs:7.3f} {g_raw:7.3f} {w_abs:7.3f} "
              f"{nsta:5d} {nfin:6d}")
    perfeat = pd.DataFrame(perfeat).sort_values("grpcv_auc", ascending=False)

    # ---- combined models, grouped CV ----
    print("\n=== COMBINED MODELS (GroupKFold by station) ===")

    def logit_fn():
        return Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(max_iter=2000, class_weight="balanced",
                                      C=1.0)),
        ])

    def hgb_fn():
        return HistGradientBoostingClassifier(
            max_depth=3, max_iter=300, learning_rate=0.05,
            l2_regularization=1.0, min_samples_leaf=30,
            class_weight="balanced", random_state=0)

    feature_sets = {
        "template_only": TEMPLATE_FEATS,
        "catalog+location": CATALOG_FEATS + LOCATION_FEATS,
        "cheap_only(baseline)": CHEAP_FEATS,
        "ours(template+cat+loc)": TEMPLATE_FEATS + CATALOG_FEATS + LOCATION_FEATS,
        "ALL": ALL_FEATS,
    }
    model_results = {}
    oof_store = {}
    for name, fs in feature_sets.items():
        X = df[fs].values.astype(float)
        for mname, mfn in [("logit", logit_fn), ("HGB", hgb_fn)]:
            oof, auc = grouped_cv_model(X, y, groups, mfn, n_splits)
            model_results[(name, mname)] = auc
            oof_store[(name, mname)] = oof
            print(f"  {name:24s} {mname:6s}  grpCV AUC = {auc:.3f}")

    # within-station baseline: train+test the TEMPLATE model INSIDE each station.
    # This is the decisive confound test -- if the template signal is real and
    # family-intrinsic (not a station-cluster artifact), it must discriminate
    # GOOD vs BAD among families recorded at the SAME station / same instrument.
    # Use the logistic pipeline (stable on small per-station n; HGB over-regularizes).
    from sklearn.model_selection import StratifiedKFold
    print("\n=== WITHIN-STATION baseline (template_only, logit, k-fold per station) ===")
    print("    (decisive confound test: does template shape split GOOD/BAD inside a station?)")
    within_aucs = []
    for g in np.unique(groups):
        m = groups == g
        yy = y[m]
        if len(np.unique(yy)) < 2 or m.sum() < 20:
            continue
        Xs = df.loc[m, TEMPLATE_FEATS].values.astype(float)
        if not np.isfinite(Xs).any():        # e.g. NLLB (no templates)
            continue
        n_sp = int(min(5, yy.sum(), (yy == 0).sum()))
        if n_sp < 2:
            print(f"  {g:6s} skipped (too few minority: GOOD={yy.sum()}, "
                  f"BAD={(yy==0).sum()})")
            continue
        skf = StratifiedKFold(n_splits=n_sp, shuffle=True, random_state=0)
        oof = np.full(len(yy), np.nan)
        try:
            for tr, te in skf.split(Xs, yy):
                mdl = logit_fn()
                mdl.fit(Xs[tr], yy[tr])
                oof[te] = mdl.predict_proba(Xs[te])[:, 1]
            a = roc_auc_score(yy, oof)
            within_aucs.append((g, a, m.sum()))
            print(f"  {g:6s} within AUC = {a:.3f}  (n={m.sum()}, GOOD={yy.sum()}, "
                  f"BAD={(yy==0).sum()})")
        except Exception as e:
            print(f"  {g:6s} skipped ({e})")
    if within_aucs:
        ws = np.array([w for _, _, w in within_aucs])
        wa = np.array([a for _, a, _ in within_aucs])
        print(f"  -> weighted-mean within-station AUC (template/logit) = "
              f"{np.average(wa, weights=ws):.3f} over {len(within_aucs)} stations")

    # ---- NEGATIVE CONTROL: shuffle labels WITHIN each station, re-run grouped CV.
    # Within-station shuffle destroys any real family-intrinsic signal but PRESERVES
    # the per-station base-rate structure. If grouped CV stays high under shuffle,
    # the "signal" is just the station-confound leaking; it must collapse to ~0.5.
    print("\n=== NEGATIVE CONTROL: within-station label shuffle (template_only/HGB) ===")
    Xtmpl = df[TEMPLATE_FEATS].values.astype(float)
    shuf_aucs = []
    for s in range(5):
        rng = np.random.RandomState(100 + s)
        y_sh = y.copy()
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            y_sh[idx] = rng.permutation(y[idx])
        _, a = grouped_cv_model(Xtmpl, y_sh, groups, hgb_fn, n_splits)
        shuf_aucs.append(a)
    print(f"  shuffled grpCV AUC = {np.mean(shuf_aucs):.3f} +/- {np.std(shuf_aucs):.3f} "
          f"(5 shuffles)  [real = {model_results[('template_only','HGB')]:.3f}]")

    # ---- permutation importance for best combined model (ALL feats, HGB) ----
    # Fit on all data with a final model and permute, using grouped-aware single
    # holdout is messy; we report permutation importance on a fitted full model
    # (descriptive, not for AUC claims).
    print("\n=== PERMUTATION IMPORTANCE (ALL feats, HGB fit on full data) ===")
    Xall = df[ALL_FEATS].values.astype(float)
    full = hgb_fn().fit(Xall, y)
    pi = permutation_importance(full, Xall, y, scoring="roc_auc",
                                n_repeats=20, random_state=0, n_jobs=4)
    imp = pd.DataFrame(dict(feature=ALL_FEATS,
                            imp_mean=pi.importances_mean,
                            imp_std=pi.importances_std)
                       ).sort_values("imp_mean", ascending=False)
    for _, r in imp.head(15).iterrows():
        print(f"  {r['feature']:22s} {r['imp_mean']:+.4f} +/- {r['imp_std']:.4f}")

    # ---- screening-gain at best grouped model ----
    best_set, best_mdl = max(model_results, key=model_results.get)
    best_auc = model_results[(best_set, best_mdl)]
    best_oof = oof_store[(best_set, best_mdl)]
    print(f"\n=== BEST GROUPED MODEL: {best_set}/{best_mdl}  AUC={best_auc:.3f} ===")
    # screening: flag BAD = low prob-good. Want to skip BADs without dropping GOODs.
    print("  (screening to flag/skip BAD families at various thresholds)")
    print(f"  {'thr':>6s} {'skip%':>7s} {'BADskip':>8s} {'GOODlost':>9s} "
          f"{'precBAD':>8s} {'recBAD':>8s}")
    n_bad = (y == 0).sum()
    n_good = (y == 1).sum()
    rows_screen = []
    for thr in [0.3, 0.4, 0.5, 0.6, 0.7]:
        flag_bad = best_oof < thr          # predicted BAD -> would skip
        skipped = flag_bad.sum()
        bad_skipped = ((y == 0) & flag_bad).sum()
        good_lost = ((y == 1) & flag_bad).sum()
        prec = bad_skipped / skipped if skipped else np.nan
        rec = bad_skipped / n_bad if n_bad else np.nan
        rows_screen.append((thr, 100 * skipped / len(y), bad_skipped, good_lost,
                            prec, rec))
        print(f"  {thr:6.2f} {100*skipped/len(y):7.1f} {bad_skipped:8d} "
              f"{good_lost:9d} {prec:8.3f} {rec:8.3f}")

    # ----------------------------------------------------------------- figure
    make_figure(perfeat, model_results, oof_store, y, imp, best_set, best_mdl)

    # final structured summary for the report
    print("\n=== SUMMARY (for note) ===")
    top_feat = perfeat.iloc[0]
    print(f"best single feature (grpCV): {top_feat['feature']} = "
          f"{top_feat['grpcv_auc']:.3f} (raw {top_feat['grpcv_raw']:.3f}, "
          f"within {top_feat['within_auc']:.3f})")
    print(f"best combined model (grpCV): {best_set}/{best_mdl} = {best_auc:.3f}")
    print(f"template_only HGB grpCV = {model_results[('template_only','HGB')]:.3f}")
    print(f"cheap_only HGB grpCV    = {model_results[('cheap_only(baseline)','HGB')]:.3f}")
    print(f"ALL HGB grpCV           = {model_results[('ALL','HGB')]:.3f}")
    perfeat.to_csv(f"{ROOT}/data/family_predictor_perfeature_auc.csv", index=False)
    imp.to_csv(f"{ROOT}/data/family_predictor_importance.csv", index=False)
    print(f"[out] wrote per-feature AUC + importance CSVs")


def make_figure(perfeat, model_results, oof_store, y, imp, best_set, best_mdl):
    fig = plt.figure(figsize=(16, 10))

    # (1) per-feature grouped vs within-station AUC
    ax1 = fig.add_subplot(2, 2, 1)
    pf = perfeat.sort_values("grpcv_auc")
    ypos = np.arange(len(pf))
    ax1.barh(ypos, pf["grpcv_auc"], color="#3b6fb0", label="grouped-CV (by station)")
    ax1.scatter(pf["within_auc"], ypos, color="#d1495b", s=22, zorder=3,
                label="within-station")
    ax1.axvline(0.5, color="k", ls="--", lw=1)
    ax1.set_yticks(ypos)
    ax1.set_yticklabels(pf["feature"], fontsize=7)
    ax1.set_xlabel("AUC (sign-agnostic)")
    ax1.set_xlim(0.40, max(0.75, pf["grpcv_auc"].max() + 0.05))
    ax1.set_title("Per-feature predictiveness\n(grouped CV vs within-station)")
    ax1.legend(fontsize=7, loc="lower right")

    # (2) combined-model grouped AUC bars
    ax2 = fig.add_subplot(2, 2, 2)
    keys = list(model_results.keys())
    labels = [f"{a}\n{b}" for a, b in keys]
    vals = [model_results[k] for k in keys]
    cols = ["#2a9d8f" if "HGB" in k[1] else "#e9c46a" for k in keys]
    ax2.bar(range(len(vals)), vals, color=cols)
    ax2.axhline(0.5, color="k", ls="--", lw=1)
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, rotation=60, ha="right", fontsize=6)
    ax2.set_ylabel("grouped-CV AUC")
    ax2.set_ylim(0.4, max(0.75, max(vals) + 0.05))
    ax2.set_title("Combined-model grouped-CV AUC\n(green=HGB, yellow=logit)")
    for i, v in enumerate(vals):
        ax2.text(i, v + 0.005, f"{v:.2f}", ha="center", fontsize=6)

    # (3) ROC curves for key models
    ax3 = fig.add_subplot(2, 2, 3)
    for name in ["template_only", "cheap_only(baseline)", "ALL"]:
        for mdl in ["HGB"]:
            oof = oof_store[(name, mdl)]
            m = np.isfinite(oof)
            fpr, tpr, _ = roc_curve(y[m], oof[m])
            auc = roc_auc_score(y[m], oof[m])
            ax3.plot(fpr, tpr, lw=1.6, label=f"{name}/{mdl} (AUC={auc:.2f})")
    ax3.plot([0, 1], [0, 1], "k--", lw=1)
    ax3.set_xlabel("False positive rate (GOOD flagged BAD)")
    ax3.set_ylabel("True positive rate")
    ax3.set_title("ROC (grouped-CV OOF)")
    ax3.legend(fontsize=7, loc="lower right")

    # (4) permutation importance
    ax4 = fig.add_subplot(2, 2, 4)
    ii = imp.sort_values("imp_mean").tail(15)
    ax4.barh(range(len(ii)), ii["imp_mean"], xerr=ii["imp_std"],
             color="#6a4c93")
    ax4.set_yticks(range(len(ii)))
    ax4.set_yticklabels(ii["feature"], fontsize=7)
    ax4.set_xlabel("permutation importance (AUC drop)")
    ax4.set_title("Permutation importance (ALL/HGB, full-fit)\n"
                  "(descriptive; can reflect station confound)")

    fig.suptitle("Family CWI-suitability predictability from pre-densify signatures",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{ROOT}/figures/family_predictor_auc.png"
    fig.savefig(out, dpi=130)
    print(f"\n[out] wrote {out}")


if __name__ == "__main__":
    main()
