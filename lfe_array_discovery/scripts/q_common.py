#!/usr/bin/env python
"""FINAL SPECIFICATION (PLAN.md sec.12 / results/merlin_ruling_step1cprime_2026-09-20.md) --
shared helpers for the q_ stages.  Prefix `q_`.

NOTHING here retunes an engine constant.  The q_ stages import `z_` and `w_` and add only:

  * a cache index that can reach BOTH bed years (2005 in this study's cache20hz, 2012 in
    lfe_catalog_repair/cache20hz), needed by the quantisation A/B of Step 1;
  * the fixed cut/stack window of Step 1 and the +-1.0 s cc-lag re-alignment;
  * the three recall denominators of Merlin's Q4 (all events / 60-s-independent outside the seed
    window / the LOUD subset: network-median 3-s S-energy excess >= 2.0 at the family's four
    nearest admitted sites);
  * the truth-seed labelling of Step 3.

Every accept / abandon constant of the ruling lives here, ONE definition, read by the stage
scripts and hashed into results/freeze_step1c_prime_v3.json.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import z_common as Z                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

FS = W.FS
OUT = f"{W.RES}/final"
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ Step 1 (quantisation A/B) ---
S1_N_FAMILIES = 5             # 5 hidden families ...
S1_MIN_EVENTS_2005 = 300      # ... with >= 300 catalogued events in the 2005 bed ...
S1_MIN_EVENTS_2012 = 100      # ... and >= 100 catalogued 2012 events for the reference stack
S1_SHARED = ("PGC", "B003")   # the shared permanents present in BOTH caches
QWIN_PRE_S = 2.0              # stack window: predicted S - 2 s ...
QWIN_LEN_S = 8.0              # ... 8 s long (S + early coda, inside the 10-s engine window)
CC_LAG_S = 1.0                # cc-peak search for the corrected 2005 time
S1_MIN_2012_FOR_CORRECTION = 30      # a 2012 reference stack needs >= 30 events ...
S1_MIN_REF_COHERENCE = 1.5           # ... and coherence >= 1.5 at that station to be used
S1_ACC_COH_2012 = 2.0         # accept: 2012 coherence >= 2 at >= 1 shared station ...
S1_ACC_FAM_2012 = 3           # ... for >= 3 of the 5 families
S1_ACC_COH_2005_QUANT = 1.2   # accept: median 2005 QUANTISED coherence <= 1.2
S1_ACC_LAG_SPREAD_S = 0.5     # accept: corrected-lag spread (p90 - p10) >= 0.5 s
S1_ACC_COH_POLARIS = 2.0      # accept: corrected-time POLARIS coherence >= 2 ...
S1_ACC_POLARIS_SITES = 2      # ... at >= 2 distinct POLARIS sites
S1_ABANDON_COH_2012 = 1.5     # abandon: the BEST 2012 coherence anywhere is below this

# ------------------------------------------------------------------ Step 2 (corrected oracle) ---
S2_ACC_C30 = 3                # accept C_30 >= 3; below that -> STOP ("this chain")
LOUD_EXCESS = 2.0             # denominator (iii): network-median 3-s S excess >= 2.0 ...
LOUD_N_SITES = 4              # ... at the family's four nearest ADMITTED sites
EX_WIN_S = 3.0                # (the D8 statistic, z_d8.py lines 45-47, verbatim)
EX_NEIGH_S = 60.0
EX_CUT_PRE_S = 65.0
EX_CUT_LEN_S = 130.0

# ------------------------------------------------------------------ Step 3 (bootstrap test) -----
S3_MATCH_TOL_S = 2.5          # a node matches a catalogued event within +- 2.5 s ...
S3_MIN_NODES = 2              # ... and a component needs >= 2 matched nodes of >= 2 distinct events
S3_CHANCE_SHIFT_S = 300.0     # chance count from a +300 s catalogue shift
S3_ABANDON_CHANCE_FRAC = 0.30 # abandon if chance-matched >= 30 % of matched
S3_CONT_MIN_INDEP = 5         # continuation filter: >= 5 independent hits ...
S3_CONT_MIN_DAYS = 2          # ... on >= 2 ETS days ...
S3_CONT_SCRAMBLE_RATIO = 2.0  # ... AND >= 2 x the template-scramble twin
S3_ITER_MIN_INDEP = 5         # v_iterate use_all with min_indep 5
S3_ACC_MIN = 3                # accept: rediscovered >= max(3, 0.3 K)
S3_ACC_FRAC = 0.30
S3_ACC_TWIN_FRAC = 0.10       # accept-B also needs twins <= 10 % of true-seed survivors
S3_ABANDON_TWIN_FRAC = 0.50   # abandon: twins >= 50 % of true-seed survivors under B

# ------------------------------------------------------------------ Step 5 (the blind run) ------
S5_NULL_REALISATIONS = 3


def log(*a):
    W.log(*a)


def marker(path, text):
    with open(path, "w") as f:
        f.write(text if text.endswith("\n") else text + "\n")


# ================================================================== catalogue ====================
_B = [None]


def bostock_all():
    """The full Bostock export on the frozen time base (p_common.bostock), with date/hour."""
    if _B[0] is None:
        import p_common as P
        B = P.bostock().copy()
        B["date"] = B.time.dt.strftime("%Y-%m-%d")
        B["hh"] = B.time.dt.hour.astype(int)
        _B[0] = B
    return _B[0]


def hidden_sets():
    hid, iso, _ = S.holdouts()
    return sorted(int(x) for x in hid), sorted(int(x) for x in iso)


def ets_dates():
    """The 24 ETS-window dates of the 2005 bed (quiet days excluded) as 'YYYY-MM-DD'."""
    return sorted(d.strftime("%Y-%m-%d") for d in W.ets_days())


def bed_dates():
    return sorted(d.strftime("%Y-%m-%d") for d in W.bed_days())


def events_2005(bkeys=None):
    """{bkey: sorted t_origin} for the 2005 ETS-WINDOW days of the bed."""
    B = bostock_all()
    days = set(ets_dates())
    g = B[(B.year == 2005) & B.date.isin(days)]
    if bkeys is not None:
        g = g[g.bkey.isin(set(int(x) for x in bkeys))]
    return {int(b): np.sort(x.t_origin.to_numpy(float)) for b, x in g.groupby("bkey")}


def events_2012(bkeys=None, day0="2012-08-15", day1="2012-10-05"):
    """{bkey: sorted t_origin} for the 2012 bed window (sub-second catalogue times)."""
    B = bostock_all()
    g = B[(B.year == 2012) & (B.date >= day0) & (B.date <= day1)]
    if bkeys is not None:
        g = g[g.bkey.isin(set(int(x) for x in bkeys))]
    return {int(b): np.sort(x.t_origin.to_numpy(float)) for b, x in g.groupby("bkey")}


# ================================================================== two-year cache ===============
_CI = {}


def cache_index(years, stas):
    """{(sta, year, doy): path} over BOTH bed roots (2005 study cache + 2012 repair cache)."""
    key = (tuple(sorted(int(y) for y in years)), tuple(sorted(stas)))
    if key in _CI:
        return _CI[key]
    want_y, want_s = set(int(y) for y in years), set(stas)
    idx = {}
    for r in (W.BED, W.STUDY_CACHE):
        for p in glob.glob(f"{r}/*/*/*.npz"):
            parts = p.split(os.sep)
            sta = parts[-3].split(".")[1]
            yr = int(parts[-2])
            if sta not in want_s or yr not in want_y:
                continue
            idx.setdefault((sta, yr, int(parts[-1][:3])), p)
    _CI[key] = idx
    return idx


class TwoYearCutter:
    """Absolute-time 3-C ZNE cuts out of either cache year (LRU over w_common.load_day)."""

    def __init__(self, years, stas, maxn=48):
        from collections import OrderedDict
        self.idx = cache_index(years, stas)
        self.d = OrderedDict()
        self.maxn = int(maxn)

    def day(self, sta, year, doy):
        k = (sta, int(year), int(doy))
        if k in self.d:
            self.d.move_to_end(k)
            return self.d[k]
        v = W.load_day(sta, int(year), int(doy), self.idx)
        self.d[k] = v
        if len(self.d) > self.maxn:
            self.d.popitem(last=False)
        return v

    def cut(self, sta, t_abs, len_s):
        ts = pd.Timestamp(float(t_abs), unit="s")
        v = self.day(sta, int(ts.year), int(ts.dayofyear))
        if v is None:
            return None
        t0, X = v
        n = int(round(len_s * FS))
        i0 = int(round((float(t_abs) - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        seg = X[:, i0:i0 + n].astype(np.float64)
        if not np.isfinite(seg).all() or (seg ** 2).sum() <= 0:
            return None
        return seg


def stack_at(cut, sta, times, tt, pre_s=QWIN_PRE_S, len_s=QWIN_LEN_S, max_events=0):
    """(stack, N, coherence) of a PLAIN LINEAR rms-normalised 3-C stack at t + tt - pre_s."""
    t = np.asarray(times, float)
    if max_events and len(t) > max_events:
        t = t[np.linspace(0, len(t) - 1, max_events).round().astype(int)]
    segs = []
    for tv in t:
        w = cut.cut(sta, float(tv) + float(tt) - pre_s, len_s)
        if w is not None:
            segs.append(w)
    if len(segs) < 5:
        return None, len(segs), np.nan
    Sg = Z.norm_rows(np.stack(segs))
    if len(Sg) < 5:
        return None, len(Sg), np.nan
    return Sg.mean(0), int(len(Sg)), Z.stack_coherence(np.stack(segs))


def cc_lag(ref, seg, max_lag):
    """(best lag in samples, peak normalised cc) of a 3-C `seg` against `ref` within +- max_lag.

    `seg` must be `max_lag` samples longer than `ref` on BOTH sides."""
    n = ref.shape[1]
    r = ref - ref.mean(1, keepdims=True)
    rn = np.sqrt((r ** 2).sum())
    if rn <= 0:
        return 0, np.nan
    best, bcc = 0, -np.inf
    for k in range(-max_lag, max_lag + 1):
        a = max_lag + k
        w = seg[:, a:a + n]
        w = w - w.mean(1, keepdims=True)
        wn = np.sqrt((w ** 2).sum())
        if wn <= 0:
            continue
        c = float((r * w).sum() / (rn * wn))
        if c > bcc:
            bcc, best = c, k
    return best, bcc


# ================================================================== geometry / travel times ======
_FT = [None]


def fam_tt():
    if _FT[0] is None:
        mem, of, order, det = W.sites()
        _FT[0] = Z.fam_tt_map(det, W.station_table())
    return _FT[0]


def fam_shift(bk, det=None):
    """The pair-seed node axis: median over detection stations of fam_tt, minus 1 s."""
    det = W.sites()[3] if det is None else det
    FT = fam_tt()
    tt = [FT[(int(bk), s)] for s in det if (int(bk), s) in FT]
    return float(np.median(tt)) - 1.0 if tt else 0.0


def family_ll():
    import p_common as P
    F = P.families()
    return {int(b): (float(F.lat[b]), float(F.lon[b])) for b in F.index}


def sites_by_distance(bk, exclude=()):
    """Admitted SITES ordered by epicentral distance of their representative member to family bk."""
    mem, of, order, det = W.sites()
    ST = W.station_table()
    fam = family_ll()
    la, lo = fam[int(bk)]
    rows = []
    for q in order:
        ms = [s for s in mem[q] if s in ST.index and s not in exclude]
        if not ms:
            continue
        d = float(np.median([S.epi_km(la, lo, ST.lat[s], ST.lon[s]) for s in ms]))
        rows.append((d, q, ms[0]))
    return [(q, s, d) for d, q, s in sorted(rows)]


def polaris_sites(bk, n=2):
    """The `n` nearest admitted sites that are NOT the shared permanents (i.e. the POLARIS/C8/PO
    temporary deployment) -- the sites the corrected-time restack is judged on."""
    out = []
    for q, s, d in sites_by_distance(bk):
        if s in S1_SHARED:
            continue
        out.append((q, s, d))
        if len(out) >= n:
            break
    return out


# ================================================================== denominators =================
def independent_outside(times, t0, t1, sep_s=W.INDEP_S, near_s=W.NEAR_SEED_S):
    """Denominator (ii): 60-s-independent catalogued events, excluding those the near-seed clause
    would exclude (within +- NEAR_SEED_S of EITHER end of the seed window)."""
    t = np.sort(np.asarray(times, float))
    if not len(t):
        return t
    keep = W.independent(t, sep_s=sep_s)
    t = t[keep]
    if np.isfinite(t0) and np.isfinite(t1):
        m = ~((np.abs(t - float(t0)) <= near_s) | (np.abs(t - float(t1)) <= near_s))
        t = t[m]
    return t


def event_s_excess(cut, bk, times, n_sites=LOUD_N_SITES):
    """Denominator (iii): the per-event network-median 3-s S-energy excess at the family's
    `n_sites` nearest admitted sites (z_d8's statistic, verbatim).  Returns a DataFrame
    bkey, t_origin, n_sites_used, excess_median, plus one column per site."""
    FT = fam_tt()
    near = sites_by_distance(bk)[:n_sites]
    nw = int(round(EX_WIN_S * FS))
    iS = int(round(EX_CUT_PRE_S * FS))
    lo0 = max(iS - int(EX_NEIGH_S * FS), 0)
    rows = []
    for tv in np.asarray(times, float):
        per = {}
        for q, s, _d in near:
            tt = FT.get((int(bk), s))
            if tt is None:
                continue
            seg = cut.cut(s, tv + tt - EX_CUT_PRE_S, EX_CUT_LEN_S)
            if seg is None:
                continue
            e = (seg ** 2).sum(0)
            cs = np.concatenate([[0.0], np.cumsum(e)])
            win = cs[nw:] - cs[:-nw]
            hi = min(iS + int(EX_NEIGH_S * FS), len(win))
            if hi <= lo0 or iS >= len(win):
                continue
            med = float(np.median(win[lo0:hi]))
            if not np.isfinite(med) or med <= 0:
                continue
            per[q] = float(win[iS] / med)
        r = dict(bkey=int(bk), t_origin=float(tv), n_sites_used=len(per),
                 excess_median=float(np.median(list(per.values()))) if per else np.nan)
        for q, v in per.items():
            r[f"ex_{q}"] = v
        rows.append(r)
    return pd.DataFrame(rows)


# ================================================================== truth seeds (Step 3) =========
def label_true_seeds(members, seeds, hidden, shift_s=0.0, tol=S3_MATCH_TOL_S,
                     min_nodes=S3_MIN_NODES):
    """Label each seed component a TRUE SEED of hidden family k.

    A component is a true seed of k when >= `min_nodes` of its node times fall within +- `tol` of
    `t_origin + median fam_tt - 1 s` of >= `min_nodes` DISTINCT catalogued events of k
    (D3's matching, y_step1diag.py ~550).  `shift_s` shifts the CATALOGUE (the chance control).

    Returns a DataFrame cid, bkey, n_nodes_matched, n_events_matched."""
    ev = events_2005(hidden)
    node = {c: np.sort(g.t.to_numpy(float)) for c, g in members.groupby("cid")}
    rows = []
    for bk in sorted(hidden):
        e = ev.get(int(bk))
        if e is None or not len(e):
            continue
        axis = np.sort(e + fam_shift(bk) + float(shift_s))
        for c in seeds.cid:
            t = node.get(c)
            if t is None or len(t) < min_nodes:
                continue
            j = np.searchsorted(axis, t)
            hit_node, hit_ev = [], set()
            for i, tv in enumerate(t):
                for k in (j[i] - 1, j[i]):
                    if 0 <= k < len(axis) and abs(axis[k] - tv) <= tol:
                        hit_node.append(i)
                        hit_ev.add(int(k))
                        break
            if len(set(hit_node)) >= min_nodes and len(hit_ev) >= min_nodes:
                rows.append(dict(cid=c, bkey=int(bk), n_nodes_matched=len(set(hit_node)),
                                 n_events_matched=len(hit_ev), n_nodes=int(len(t))))
    return pd.DataFrame(rows, columns=["cid", "bkey", "n_nodes_matched", "n_events_matched",
                                       "n_nodes"])


# ================================================================== misc =========================
def jdump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=str)
    return obj


def jload(path, default=None):
    try:
        return json.load(open(path))
    except Exception:                                                       # noqa: BLE001
        return default


def n_surv(path):
    T = W.read_csv_safe(path, W.SURV_COLS)
    return int(T.survives.astype(bool).sum()) if len(T) else 0


def q_scripts():
    return sorted(glob.glob(f"{HERE}/q_*.py")) + sorted(glob.glob(f"{HERE}/q_*.sh"))
