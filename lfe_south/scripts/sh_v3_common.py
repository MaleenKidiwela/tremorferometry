#!/usr/bin/env python
"""SH_V3 (2026-09-29) -- shared helpers for the Shelly (2025) Text S1 reproduction (v3).

Final-catalogue rules, Text S1 "Integrated Catalog":
  * summed cc >= 9 x daily MAD              (our net_mad = (sum - median)/MAD; the daily median of
                                              the sum is ~0 (|med| < 3e-3 vs MAD ~0.38), so this is
                                              his ccsum / MAD to < 1 %)
  * mean cc over available channels >= 0.12 (our meancc = sum / live channels at the peak)
  * absolute summed cc >= 2.75               (S1 sec. 1 "additional minimum threshold")
  * declustering: detections within a 4 s window -> keep only the highest summed cc.
    His published catalogue has min inter-event dt = 4.00 s over ALL families (0 pairs < 4 s), so
    the suppression is applied across all families (greedy by ccsum, |dt| < 4 s); a cross-family-
    only variant is available (same_family=False).
Scoring: t8_score.one_to_one/score and the t8_score_aa.py truth-day filter, reproduced verbatim.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import *
from t8_score import one_to_one, score

V3 = os.environ.get("SH_V3_DIR", f"{RES}/sh_v3")   # v3b (2026-09-29): overridable, default unchanged
TEST = f"{RES}/aa_shelly/stage2_days_test.json"
BUILD = f"{RES}/aa_shelly/stage2_days_build.json"
SRC = "Shelly2025"
MEANCC_MIN = 0.12; CCSUM_MIN = 2.75; DECL_S = 4.0


def his_catalogue():
    S = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
    S.columns = [c.strip('%') for c in S.columns]
    S["t"] = (pd.to_datetime(dict(year=S.year, month=S.month, day=S.day, hour=S.hour, minute=S.minute,
              second=S.second.astype(int))).astype("datetime64[ns]") + pd.to_timedelta(S.second % 1, unit="s"))
    S["fam"] = S.idnum.astype(str)
    return S


def family_days(path, source=SRC):
    dd = json.load(open(path))
    return {k.split("|")[1]: {pd.Timestamp(x) for x in v} for k, v in dd.items() if k.split("|")[0] == source}


def decluster(D, win=DECL_S, same_family=True):
    """greedy by ccsum (desc): a detection is dropped when an already-kept detection lies within
    |dt| < win (any family if same_family else only other families).  Returns bool keep mask."""
    import numba
    if not len(D):
        return np.zeros(0, bool)
    t = D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    fam = pd.factorize(D.fam.astype(str))[0].astype(np.int64)
    cs = D.ccsum.values.astype(np.float64)
    o = np.argsort(t, kind="stable"); ts = t[o]; fs = fam[o]; cso = cs[o]
    lo = np.searchsorted(ts, ts - win, side="right"); hi = np.searchsorted(ts, ts + win, side="left")
    order = np.argsort(-cso, kind="stable")
    keep = _greedy(ts, fs, lo, hi, order, win, same_family)
    out = np.zeros(len(D), bool); out[o] = keep
    return out


_GREEDY = []
def _greedy(*a):
    if not _GREEDY:
        import numba
        @numba.njit
        def g(ts, fs, lo, hi, order, win, same):
            n = ts.size; keep = np.zeros(n, np.bool_)
            for q in range(n):
                i = order[q]; ok = True
                for j in range(lo[i], hi[i]):
                    if j != i and keep[j] and abs(ts[j] - ts[i]) < win and (same or fs[j] != fs[i]):
                        ok = False; break
                keep[i] = ok
            return keep
        _GREEDY.append(g)
    return _GREEDY[0](*a)


def apply_rules(D, k, meancc=True, floor=True, decl=True, same_family=True):
    m = D.net_mad.values >= k
    if meancc: m &= D.meancc.values >= MEANCC_MIN
    if floor: m &= D.ccsum.values >= CCSUM_MIN
    E = D[m]
    if decl:
        E = E[decluster(E, same_family=same_family)]
    return E


def truth_for(T, mode, days_file=TEST):
    """t8_score_aa.py truth-day filter: 'source' = union of the source's test days for every family,
    'fam' = each family's own test days."""
    DF = family_days(days_file)
    if mode == "source":
        U = set().union(*DF.values()); DF = {f: U for f in DF}
    day = T.t.dt.floor("D")
    return T[[(f in DF and d in DF[f]) for f, d in zip(T.fam, day)]].copy(), DF


def ours_for(O, mode, DF):
    if mode == "source":
        return O
    day = O.t.dt.floor("D")
    return O[[(f in DF and d in DF[f]) for f, d in zip(O.fam.astype(str), day)]]


def score_all(O, T, tol=2.0):
    """family-specific (t8_score per family) + family-agnostic (pooled one-to-one) scores."""
    rows = []
    for f in sorted(set(T.fam) | set(O.fam.astype(str))):
        oo = O[O.fam.astype(str) == f].t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
        tt = T[T.fam == f]
        r = score(oo, tt.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9, tt.cc.values, tol)
        r.update(fam=f); rows.append(r)
    F = pd.DataFrame(rows)
    ot = O.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    tt = T.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    pairs = one_to_one_pairs(ot, tt, tol)
    fo = O.fam.astype(str).values; ft = T.fam.astype(str).values
    same = int(sum(fo[i] == ft[k] for i, k in pairs))
    return F, dict(his=len(tt), ours=len(ot), match_fam=int(F.n_match.sum()), match_any=len(pairs), same_fam=same)


def one_to_one_pairs(ours, theirs, tol):
    """t8_score.one_to_one (greedy by |dt|, one-to-one) but returning (i_ours, k_theirs) index pairs
    into the UNSORTED inputs."""
    if len(ours) == 0 or len(theirs) == 0: return []
    oi = np.argsort(ours, kind="stable"); ti = np.argsort(theirs, kind="stable")
    o = ours[oi]; t = theirs[ti]
    j = np.searchsorted(t, o)
    cand = []
    for i, x in enumerate(o):
        for k in (j[i] - 1, j[i]):
            if 0 <= k < len(t) and abs(t[k] - x) <= tol: cand.append((abs(t[k] - x), i, k))
    cand.sort(); uo = set(); ut = set(); P = []
    for dt, i, k in cand:
        if i in uo or k in ut: continue
        uo.add(i); ut.add(k); P.append((oi[i], ti[k]))
    return P
