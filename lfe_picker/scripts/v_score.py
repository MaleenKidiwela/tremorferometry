#!/usr/bin/env python
"""P4-0b scoring -- ONE pass, everything read from the freeze.

ARM (i) GIVEN seeds.  The RELATIVE clauses of `q_rediscover.py`: a family is ASSEMBLED if one of
its grown templates reaches recall >= 0.5 x recall_ref and purity >= 0.8 x purity_ref, where the
reference is the family's own truth-stack template grown on the SAME days with the SAME sites
(withheld stations are never in any sum).  Families whose reference itself has recall 0 are
reported as NO-REFERENCE and counted neither way (q_rediscover's rule).

ARM (ii) TEMPLATE-FREE seeds.  One-to-one +-2 s matching of a grown template's independent hits to
a kept family's hp events.  A template-free seed's origin axis is defined by the REGION-centre
moveout, so a true family's hits arrive on that axis displaced by ONE constant b = t_seed - tau_0
(the same for every event -- the derivation is in the freeze).  b is therefore estimated as a
single nuisance parameter per (template, family) by a grid search over +-30 s, and the search is
CHARGED to a null in which the template's hits are redrawn uniformly inside their own days.  The
strict b = 0 numbers are reported alongside.  Rediscovered = matches >= 20 AND >= 30 % of the
family's hp events; ALIAS = the same grown template meeting the clause for a second family.

ACCEPTANCE (pre-registered, from the freeze): (i) >= 2/3 of eligible kept families assembled;
(ii) >= 3 kept families rediscovered template-free with <= 1 alias each.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

OUT = V.OUT
B_GRID = np.arange(-30.0, 30.0 + 1e-9, 0.1)
N_NULL = 20


def _indep_table(prefix):
    p = f"{prefix}_indep.parquet"
    return pd.read_parquet(p) if os.path.exists(p) else pd.DataFrame(columns=["tid", "t", "date",
                                                                              "net_mad"])


def _hp_by_family():
    keep, D, _ = V.eligible_families()
    A = V.availability()
    usable = {t for t, u in zip(A.date, A.usable) if u}
    ED = V.episode_days()
    D = D[[pd.Timestamp(d) in usable and pd.Timestamp(d) in ED for d in D.dt]]
    return keep, {int(b): np.sort(g.t.to_numpy(float)) for b, g in D.groupby("bkey")}


def score_given(res, keep, EV):
    """res: dict of paths.  Returns the per-family table and the arm-(i) summary."""
    G = _indep_table(res["given"])
    R = _indep_table(res["ref"])
    dd = pd.read_csv(res["given_dedupe"]) if os.path.exists(res["given_dedupe"]) else None
    if dd is not None and len(dd):
        G = G[G.tid.isin(set(dd[dd.kept].tid))]
    rows = []
    for b in keep:
        ev = EV.get(int(b), np.zeros(0))
        rt = R[R.tid == f"REF{b}"].t.to_numpy(float)
        mr, _ = V.one2one(ev, rt, V.MATCH_TOL_S)
        mp, _ = V.one2one(rt, ev, V.MATCH_TOL_S)
        rec_ref = float(mr.mean()) if len(ev) else np.nan
        pur_ref = float(mp.mean()) if len(rt) else np.nan
        no_ref = (not np.isfinite(rec_ref)) or rec_ref <= 1e-9
        best = None
        for tid in sorted(set(G.tid)):
            if not str(tid).startswith(f"G{b}_"):
                continue
            ht = G[G.tid == tid].t.to_numpy(float)
            if len(ht) == 0:
                continue
            m1, _ = V.one2one(ev, ht, V.MATCH_TOL_S)
            m2, _ = V.one2one(ht, ev, V.MATCH_TOL_S)
            rec, pur = float(m1.mean()) if len(ev) else 0.0, float(m2.mean())
            ok = (not no_ref) and rec >= V.REL_RECALL * rec_ref and pur >= V.REL_PURITY * pur_ref
            r = dict(bkey=int(b), tid=tid, n_hits=len(ht), recall=rec, purity=pur, assembled=ok)
            if best is None or (ok, rec) > (best["assembled"], best["recall"]):
                best = r
        row = dict(bkey=int(b), n_hp=len(ev), n_ref_hits=len(rt), recall_ref=rec_ref,
                   purity_ref=pur_ref, no_reference=bool(no_ref),
                   need_recall=V.REL_RECALL * rec_ref if np.isfinite(rec_ref) else np.nan,
                   need_purity=V.REL_PURITY * pur_ref if np.isfinite(pur_ref) else np.nan)
        row.update(best or dict(tid=None, n_hits=0, recall=0.0, purity=np.nan, assembled=False))
        rows.append(row)
    T = pd.DataFrame(rows)
    testable = T[~T.no_reference]
    summ = dict(n_eligible=len(keep), n_no_reference=int(T.no_reference.sum()),
                n_testable=int(len(testable)), n_assembled=int(T.assembled.sum()),
                frac_of_testable=float(testable.assembled.mean()) if len(testable) else np.nan,
                bar="(i) >= 2/3 of eligible kept families assembled",
                pass_i=bool(T.assembled.sum() >= np.ceil(2 * len(keep) / 3)))
    return T, summ


def _match_with_offset(ht, ev, span=30.0, tol=V.MATCH_TOL_S, ncand=5):
    """Best ONE constant offset b for the template's hits against a family's events.

    The candidate b's come from the (event - hit) difference histogram (0.1 s bins, box-filtered
    over +-tol); the exact greedy one-to-one match is then evaluated only at the top `ncand` peaks.
    Exhaustive over the grid in effect, but ~100x cheaper."""
    ht = np.asarray(ht, float); ev = np.asarray(ev, float)
    if len(ht) == 0 or len(ev) == 0:
        return 0, 0.0
    d = (ev[None, :] - ht[:, None]).ravel()
    d = d[np.abs(d) <= span + tol]
    if d.size == 0:
        return 0, 0.0
    bins = np.arange(-span - tol, span + tol + 1e-9, 0.1)
    h, _ = np.histogram(d, bins=bins)
    k = int(round(2 * tol / 0.1))
    box = np.convolve(h, np.ones(k + 1), mode="same")
    order = np.argsort(-box)[:max(ncand * 8, 40)]
    cands, seen = [], []
    for i in order:
        b = float(0.5 * (bins[i] + bins[i + 1]))
        if any(abs(b - u) < tol for u in seen):
            continue
        seen.append(b); cands.append(b)
        if len(cands) >= ncand:
            break
    best = (0, 0.0)
    for b in cands:
        if abs(b) > span:
            continue
        m, _ = V.one2one(ev, ht + b, tol)
        kk = int(m.sum())
        if kk > best[0]:
            best = (kk, b)
    return best


def score_free(res, keep, EV, rng_seed=20260920):
    F = _indep_table(res["free"])
    dd = pd.read_csv(res["free_dedupe"]) if os.path.exists(res["free_dedupe"]) else None
    if dd is not None and len(dd):
        F = F[F.tid.isin(set(dd[dd.kept].tid))]
    rng = np.random.default_rng(rng_seed)
    rows = []
    for tid in sorted(set(F.tid)):
        g = F[F.tid == tid]
        ht = np.sort(g.t.to_numpy(float))
        if len(ht) < V.FREE_MIN_MATCH:
            continue
        d0 = pd.to_datetime(g.date).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
        for b in keep:
            ev = EV.get(int(b), np.zeros(0))
            if len(ev) == 0:
                continue
            k, boff = _match_with_offset(ht, ev)
            k0, _ = V.one2one(ev, ht, V.MATCH_TOL_S)
            nulls = []
            for _ in range(N_NULL):
                hn = np.sort(d0 + rng.uniform(0, 86400, len(ht)))
                nulls.append(_match_with_offset(hn, ev)[0])
            rows.append(dict(tid=tid, bkey=int(b), n_hits=len(ht), n_hp=len(ev),
                             matches=k, offset_s=boff, frac_hp=k / len(ev),
                             matches_b0=int(k0.sum()), frac_hp_b0=float(k0.mean()),
                             null_mean=float(np.mean(nulls)), null_max=int(np.max(nulls)),
                             rediscovered=bool(k >= V.FREE_MIN_MATCH
                                               and k / len(ev) >= V.FREE_MIN_FRAC
                                               and k > np.max(nulls))))
    M = pd.DataFrame(rows)
    alias = {}
    if len(M):
        for tid, g in M[M.rediscovered].groupby("tid"):
            alias[tid] = int(len(g) - 1)
    fam = {}
    if len(M):
        for b, g in M[M.rediscovered].groupby("bkey"):
            g = g.assign(alias=[alias.get(t, 0) for t in g.tid])
            ok = g[g.alias <= 1]
            fam[int(b)] = dict(n_templates=int(len(g)), min_alias=int(g.alias.min()),
                               counts=bool(len(ok) > 0),
                               best_tid=str(g.sort_values("matches", ascending=False).tid.iloc[0]),
                               best_matches=int(g.matches.max()))
    n_ok = sum(1 for v in fam.values() if v["counts"])
    summ = dict(n_templates_scored=int(M.tid.nunique()) if len(M) else 0,
                n_rediscovered_families=n_ok, families=fam,
                alias_counts=alias,
                bar="(ii) >= 3 kept families rediscovered template-free with <= 1 alias each",
                pass_ii=bool(n_ok >= 3))
    return M, summ


def burst_convergence(res):
    p = f"{res['given']}_survival.csv"
    if not os.path.exists(p):
        return {}
    Sv = pd.read_csv(p)
    Sv["bkey"] = [int(str(t)[1:].split("_")[0]) for t in Sv.tid]
    by = Sv.groupby("bkey").frac_near_seed.median()
    return dict(n_families=int(len(by)),
                frac_families_gt50pct_near_seed=float((by > 0.5).mean()) if len(by) else np.nan,
                median_frac_near_seed=float(by.median()) if len(by) else np.nan,
                per_family={int(k): float(v) for k, v in by.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default=OUT)
    a = ap.parse_args()
    s = a.stem
    res = dict(given=f"{s}/run_given_it2", ref=f"{s}/run_ref",
               free=f"{s}/run_free_it2",
               given_dedupe=f"{s}/dedupe_given.csv", free_dedupe=f"{s}/dedupe_free.csv")
    keep, EV = _hp_by_family()
    out = dict(eligible_families=keep, n_eligible=len(keep))
    T, si = score_given(res, keep, EV)
    T.to_csv(f"{s}/score_given.csv", index=False)
    out["arm_i_given"] = si
    try:
        M, sii = score_free(res, keep, EV)
        M.to_csv(f"{s}/score_free.csv", index=False)
    except Exception as e:
        sii = dict(error=f"{type(e).__name__}: {e}", pass_ii=False)
    out["arm_ii_free"] = sii
    out["burst_convergence"] = burst_convergence(res)
    for k, f in (("hp_seed_z_band", "hp_seed_z.json"), ("free_seed_stats", "seeds_free.json"),
                 ("given_seed_stats", "seeds_given.json")):
        p = f"{s}/{f}"
        if os.path.exists(p):
            out[k] = json.load(open(p))
    pi, pii = si.get("pass_i"), sii.get("pass_ii")
    out["verdict"] = ("BOTH PASS -> P4-2 prior ~30 %" if (pi and pii) else
                      "GIVEN PASS / FREE FAIL -> seeds are the bottleneck; P4-2 as an "
                      "injection-calibrated bound, prior <=10 %" if pi else
                      "GIVEN FAIL -> STOP P4-2 (do not run Puget discovery)")
    json.dump(out, open(f"{s}/P4_0B_SCORE.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "eligible_families"},
                     indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
