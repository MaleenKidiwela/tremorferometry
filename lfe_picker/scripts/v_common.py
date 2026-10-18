#!/usr/bin/env python
"""P4-0b -- the TRUE Arm B kill test on the PQ-era north.  Shared machinery.  Prefix `v_`.

Merlin's P4 review (results/merlin_review_p4_2026-09-20.md, F1/F2 and the P4-0b prescription)
established that P3D's "seed-and-grow / ARM B" was only a seed-restricted PAIRWISE graph and that
the construction pre-registered in results/p3d/arm_b_preregistration.md steps 2-5 (matched-filter
growth from a seed -> iterated stacks -> acceptance) has NEVER been run.  This module implements
that construction, on the north, where the 130 Bostock families minus the 30 hidden ones are truth.

Nothing here edits an existing script or a frozen json.  Everything is read-only against
`s_common` / `p_common` / `g_common` / `o_common` and the frozen 20 Hz cache.

THE MEASUREMENT (Merlin, binding)
  template      per-station window at the seed's moveout (20 s from predicted S-6 s in pass 0;
                12 s from S-2 s after each iteration), each channel normalised to unit L2
  cc            frozen-runner recipe: full-day normalised correlation per channel by rFFT with a
                FLOAT64 prefix-sum denominator (o_run_month.py); station series = mean over 3 ch
  per-station z lowhalf MAD: m = median(cc); MAD = 1.4826*median(m - cc[cc<=m]); z = (cc-m)/MAD
                (Merlin p1p2 flaw 3: the two-sided day MAD collapses on loud days)
  per-site      mean of the member stations' z, then CLIPPED at 8 (anti-dominance)
  network       EQUAL-WEIGHT SUM over the sites observable at that sample, then converted to
                net_mad by the frozen runner's per-(template, day) median/MAD normalisation.
                The second normalisation is what makes ONE fixed per-template threshold comparable
                across days whose site count differs (per-day availability is NOT intersected --
                intersecting is what starved NET-B).  It is the frozen runner's own `net_mad`.
  hit           net_mad >= T_seed AND >= 2 sites with per-site z >= 3 at the hit sample AND
                >= 3 sites observable; declustered at 30 s; INDEPENDENT hits are >= 60 s apart
  T_seed        per template, frozen FA-only rule (lfe_extension/scripts/c_tiers_by_epoch.py):
                smallest T on a 0.25 grid at or above the 8.0 floor with <= 0.1 false hits per
                template-day on the quiet days

SITES.  PGCq (PGC+B009+B010+B011), VGZ, SYMB, BUCKS, LHLYB, ALBHB, DAOB.  The three held-out
stations of results/holdouts.json (B926, NLLB, LRIV) and the four P3D withheld sites
(B926q incl. CLRS, GOBB, LOPEZ, TURTL) never enter any sum.  Hidden families are never seeded,
never scored and never matched.  2011/2013 are sealed by s_common.cache_index.
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
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
OUT = f"{RES}/p4/p4_0b"
FIG = f"{PICK}/figures"
sys.path.insert(0, HERE)

import s_common as S                                                        # noqa: E402

FS = S.FS
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
EQCAT = "/home/jovyan/tremorferometry/lfe_extension/results/trackC/eq_catalog.csv"

EPISODES = ["S18", "S21", "S22"]

# ------------------------------------------------------------------ sites / stations ------------
SITE_MEMBERS = {
    "PGCq": ["PGC", "B009", "B010", "B011"],
    "VGZ": ["VGZ"], "SYMB": ["SYMB"], "BUCKS": ["BUCKS"],
    "LHLYB": ["LHLYB"], "ALBHB": ["ALBHB"], "DAOB": ["DAOB"],
}
SITE_ORDER = ["PGCq", "VGZ", "SYMB", "BUCKS", "LHLYB", "ALBHB", "DAOB"]
DET_STAS = [s for q in SITE_ORDER for s in SITE_MEMBERS[q]]
SITE_OF = {s: q for q, mem in SITE_MEMBERS.items() for s in mem}
# never in ANY network sum: the three held-out stations + the four P3D withheld sites' members
WITHHELD_STAS = ["B926", "NLLB", "LRIV", "CLRS", "GOBB", "LOPEZ", "TURTL"]
MIN_SITES = 3

# ------------------------------------------------------------------ frozen parameters -----------
TPL0_LEN_S = 20.0          # pass-0 template length (Merlin: longer than 12 s; the true cell may be
TPL0_PRE_S = 6.0           # ~25 km from the centre) -- cut from predicted S - 6 s
TPL1_LEN_S = 12.0          # after each iteration: 12 s from S - 2 s
TPL1_PRE_S = 2.0
ALIGN_ITERS, ALIGN_LAG = 2, 30      # frozen align_stack recipe (+-1.5 s at 20 Hz)
DECLUSTER_S = 30.0
INDEP_S = 60.0
T_FLOOR, T_GRID, T_MAX = 8.0, 0.25, 60.0
FA_PER_TPL_DAY = 0.1
SITE_Z_COINC = 3.0         # >= 2 sites at this per-site z at the hit sample
SITE_Z_CLIP = 8.0
SURV_MIN_INDEP = 20
SURV_MIN_DAYS = 3
SURV_MIN_EPISODES = 2
SURV_MAX_NEAR_SEED = 0.20  # <= 20 % of hits within +-10 min of the seed
NEAR_SEED_S = 600.0
NMOVEOUT_RATIO = 3.0
NMOVEOUT_MIN_DT = 5.0      # wrong cell: predicted S differs >= 5 s at >= half the stations
DEDUPE_OVERLAP = 0.50
DEDUPE_TOL_S = 2.0
SCREEN_DB = 3.0
MATCH_TOL_S = 2.0
# arm (i) relative clauses (q_rediscover.py)
REL_RECALL, REL_PURITY = 0.5, 0.8
# arm (ii) absolute clauses
FREE_MIN_FRAC, FREE_MIN_MATCH = 0.30, 20
# template-free seed rule (Merlin P4-2 step 2, corrected)
SEED_WIN_S, SEED_HOP_S = 12.0, 6.0
SEED_MIN_EPI_HOUR = 5
SEED_EPI_KM = 50.0
SEED_SITE_Z = 3.0
SEED_MIN_SITES = 3
SEED_ZMED_LO, SEED_ZMED_HI = 4.0, 25.0
SEED_EQ_VETO_S = 60.0
SEED_HF_MAX = 3.5
SEED_MAX_PER_HOUR = 3
SEED_MIN_SEP_S = 120.0
SEED_BUDGET = 1500
GIVEN_MIN_HP = 50
GIVEN_SEEDS_PER_FAM = 3
QUIET_TREMOR_KM = 150.0
N_QUIET_DAYS = 20


# ------------------------------------------------------------------ geometry --------------------
def station_table():
    p = f"{S.P3D}/stations_p3d_v2.csv"
    if os.path.exists(p):
        return pd.read_csv(p).set_index("sta")
    return S.station_table()


def region_centre():
    """Box centre from p_common.GRID_LAT/GRID_LON, at the DEPTH OF THE FAMILY FIELD.

    DOCUMENTED DEVIATION.  Slab2 at (48.50, -123.75) returns 21.9 km, ~13 km shallower than the
    Bostock family field this kill test has to find (median 35 km).  Using it would mis-cut every
    template by ~3 s in the common term and ~1 s in the relative moveout.  The lat/lon are Merlin's
    (p_common.GRID_LAT/GRID_LON); the depth is the median of the 130 family depths.  Both the
    slab2 value and the adopted value are written into the freeze."""
    import p_common as P
    lat = 0.5 * (P.GRID_LAT[0] + P.GRID_LAT[1])
    lon = 0.5 * (P.GRID_LON[0] + P.GRID_LON[1])
    F = S.families()
    dep = float(np.median(F.dep.values))
    dep_slab2 = float(np.atleast_1d(P.slab2_interp()(lat, lon))[0])
    return lat, lon, dep, dep_slab2


_TTI = [None]


def _tt_interp():
    """Bilinear Crosson-1D S travel-time table on (epicentral distance, source depth)."""
    if _TTI[0] is None:
        import p_vel1d as V
        from scipy.interpolate import RegularGridInterpolator
        dg = np.arange(0.0, 220.0, 1.0)
        zg = np.arange(10.0, 60.001, 1.0)
        T = np.empty((dg.size, zg.size))
        for i, d in enumerate(dg):
            for j, z in enumerate(zg):
                T[i, j] = V.tt(float(d), float(z), "S")
        _TTI[0] = RegularGridInterpolator((dg, zg), T, bounds_error=False, fill_value=None)
    return _TTI[0]


def cell_tt(lat, lon, dep, stas=None, ST=None):
    """{sta: predicted S delay after origin} = Crosson 1-D + the a5 station correction."""
    stas = DET_STAS if stas is None else stas
    ST = station_table() if ST is None else ST
    F = _tt_interp()
    CR = S.station_corr(stas)
    out = {}
    for s in stas:
        d = float(S.epi_km(lat, lon, ST.lat[s], ST.lon[s]))
        out[s] = float(F([[min(d, 219.0), float(np.clip(dep, 10, 60))]])[0]) + CR[s]
    return out


def family_tt(stas=None, ST=None):
    """{(bkey, sta): predicted S delay} at the 130 Bostock hypocentres."""
    stas = DET_STAS if stas is None else stas
    ST = station_table() if ST is None else ST
    F = S.families()
    G = _tt_interp()
    CR = S.station_corr(stas)
    out = {}
    for b in F.index:
        for s in stas:
            d = float(S.epi_km(F.lat[b], F.lon[b], ST.lat[s], ST.lon[s]))
            out[(int(b), s)] = float(G([[min(d, 219.0), float(np.clip(F.dep[b], 10, 60))]])[0]) + CR[s]
    return out


# ------------------------------------------------------------------ days / availability ---------
def cache_idx(stas=None):
    return S.cache_index(years=list(range(2019, 2027)), stas=DET_STAS if stas is None else stas)


def availability(idx=None, days=None):
    """PER-DAY availability -- NOT intersected across stations (intersecting starved NET-B).

    Returns a DataFrame date, stas, sites, n_sites, usable (n_sites >= MIN_SITES)."""
    idx = cache_idx() if idx is None else idx
    have = {}
    for (s, y, d) in idx:
        t = pd.Timestamp(year=y, month=1, day=1) + pd.Timedelta(days=d - 1)
        have.setdefault(t, set()).add(s)
    rows = []
    for t in sorted(have if days is None else [d for d in days if d in have]):
        ss = sorted(have[t])
        qq = sorted({SITE_OF[s] for s in ss if s in SITE_OF})
        rows.append(dict(date=t, stas=",".join(ss), sites=",".join(qq), n_sites=len(qq),
                         usable=len(qq) >= MIN_SITES))
    return pd.DataFrame(rows)


def episode_days():
    """{date: episode id} for the pooled S18+S21+S22 padded windows."""
    W = S.sse_windows()
    out = {}
    for ep in EPISODES:
        r = W[W.id == ep].iloc[0]
        for t in pd.date_range(r.start, r.end, freq="D"):
            out[t] = ep
    return out


def quiet_days(n=N_QUIET_DAYS, idx=None):
    """Cached days with >= MIN_SITES sites and NO PNSN tremor within QUIET_TREMOR_KM of the region
    centre.  Ranked by site count, then spread evenly over the PQ era (2024-07 onward, the era the
    kill test runs in), so the false-alarm rate is measured on the same network."""
    A = availability(idx)
    A = A[A.usable]
    tr = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"])
    lat, lon, _, _ = region_centre()
    tr = tr[S.epi_km(tr.lat.values, tr.lon.values, lat, lon) <= QUIET_TREMOR_KM]
    bad = set(tr.time.dt.normalize())
    Q = A[[t not in bad for t in A.date]].copy()
    era = Q[Q.date >= pd.Timestamp("2024-07-01")]
    Q = era if len(era) >= n else Q
    Q = Q.sort_values(["n_sites", "date"], ascending=[False, True])
    top = Q[Q.n_sites == Q.n_sites.max()]
    pool = top if len(top) >= n else Q
    pool = pool.sort_values("date").reset_index(drop=True)
    if len(pool) <= n:
        sel = pool
    else:
        sel = pool.iloc[np.linspace(0, len(pool) - 1, n).round().astype(int)]
    return sel.reset_index(drop=True)


def eligible_families():
    """Kept (non-hidden) Bostock families with >= GIVEN_MIN_HP hp events (p_real >= 0.8) inside the
    pooled episodes.  Hidden families are removed here and never re-enter."""
    hid, _, _ = S.holdouts()
    D = S.hp_events(min_p=0.8, episodes=EPISODES)
    D = D[~D.bkey.isin(hid)]
    c = D.groupby("bkey").size()
    keep = sorted(int(b) for b in c[c >= GIVEN_MIN_HP].index)
    return keep, D[D.bkey.isin(keep)].reset_index(drop=True), D


# ------------------------------------------------------------------ robust scales ----------------
def mad_lowhalf(x, axis=-1):
    """(median, 1.4826 x median absolute deviation of the LOWER HALF).  Loud-day proof."""
    import torch
    if isinstance(x, np.ndarray):
        m = np.median(x, axis=axis, keepdims=True)
        d = np.where(x <= m, m - x, np.nan)
        s = 1.4826 * np.nanmedian(d, axis=axis, keepdims=True)
        return m, np.maximum(s, 1e-9)
    m = x.median(dim=axis, keepdim=True).values
    d = torch.where(x <= m, m - x, torch.full_like(x, float("inf")))
    ds, _ = d.sort(dim=axis)
    nlow = (d < float("inf")).sum(dim=axis, keepdim=True).clamp_min(1)
    k = (nlow // 2).clamp_min(0)
    s = 1.4826 * ds.gather(axis, k)
    return m, s.clamp_min(1e-9)


def decluster(vals, thr, width_samp):
    """Indices of declustered peaks of `vals` (1-D torch tensor) above `thr`.

    Same rule as o_run_month.decluster_idx: a sample is kept if it is the maximum of the
    +-width/2 window around it and above the threshold."""
    import torch
    v = vals[None, None]
    w = int(width_samp)
    pooled = torch.nn.functional.max_pool1d(v, kernel_size=2 * w + 1, stride=1, padding=w)[0, 0]
    m = (vals >= thr) & (vals >= pooled - 1e-6)
    return torch.nonzero(m, as_tuple=False).squeeze(1)


def independent(times, sep_s=INDEP_S):
    """Greedy >= sep_s thinning of a sorted time array; returns the kept boolean mask."""
    t = np.asarray(times, float)
    o = np.argsort(t)
    keep = np.zeros(len(t), bool)
    last = -1e18
    for j in o:
        if t[j] - last >= sep_s:
            keep[j] = True
            last = t[j]
    return keep


def one2one(a, b, tol=MATCH_TOL_S):
    """Greedy one-to-one matching of sorted arrays a -> b within tol.  Returns a boolean mask on a
    and the matched index in b (or -1).  Same recipe as g_rediscover.one2one."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    oa, ob = np.argsort(a), np.argsort(b)
    sa, sb = a[oa], b[ob]
    used = np.zeros(len(sb), bool)
    ma = np.zeros(len(sa), bool)
    idx = np.full(len(sa), -1)
    j = 0
    for i, t in enumerate(sa):
        while j < len(sb) and sb[j] < t - tol:
            j += 1
        k, best = j, None
        while k < len(sb) and sb[k] <= t + tol:
            if not used[k] and (best is None or abs(sb[k] - t) < abs(sb[best] - t)):
                best = k
            k += 1
        if best is not None:
            used[best] = True
            ma[i] = True
            idx[i] = ob[best]
    out_m = np.zeros(len(a), bool)
    out_i = np.full(len(a), -1)
    out_m[oa] = ma
    out_i[oa] = idx
    return out_m, out_i


# ------------------------------------------------------------------ template banks ---------------
def save_bank(path, tids, stas, W, off, meta):
    """W (n_tpl, n_sta, 3, n) float32 (nan where absent); off (n_tpl, n_sta) seconds from origin."""
    np.savez_compressed(path, tids=np.asarray(tids), stas=np.asarray(stas, dtype=object),
                        W=np.asarray(W, np.float32), off=np.asarray(off, np.float64),
                        meta=json.dumps(meta, default=str))


def load_bank(path):
    z = np.load(path, allow_pickle=True)
    return (list(z["tids"]), [str(s) for s in z["stas"]], z["W"], z["off"],
            json.loads(str(z["meta"])))


def norm_channels(w):
    """Each channel to unit L2 (the frozen runner's cc = <w_c, x_c>/||x_c|| needs ||w_c|| = 1)."""
    w = np.asarray(w, np.float64)
    r = np.sqrt((w ** 2).sum(-1, keepdims=True))
    out = np.where(r > 0, w / np.maximum(r, 1e-30), 0.0)
    return out.astype(np.float32)


# ------------------------------------------------------------------ freeze helpers ---------------
def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def v_scripts():
    return sorted(glob.glob(f"{HERE}/v_*.py")) + sorted(glob.glob(f"{HERE}/v_*.sh"))


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%Y-%m-%d %H:%M:%S')}]", *a, flush=True)
