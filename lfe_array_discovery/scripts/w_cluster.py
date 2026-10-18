#!/usr/bin/env python
"""Stage 2 -- FA-ONLY THRESHOLD and CLUSTERING of the autocorrelation pairs.  Prefix `w_`.

T (prerequisite P1: NO 8.0 floor).  On the quiet hours of the bed, T is the smallest value on a
0.25 grid at or above the ENGINEERING floor 5.0 such that

    (i)  false pairs <= 1 per quiet hour,  AND
    (ii) zero quiet-hour components of size >= 5,

and the report says whether the engineering floor BINDS (T == 5.0 because the FA rule was already
satisfied there) or whether the FA rule binds.

Clustering (Merlin Q1/(a)).  Nodes = window times merged at 2 s; edges = hits; connected
components over the whole DAY (so the "hour +- 1 h" clause is natural).  A component is a SEED
CLUSTER when it has >= 5 members inside one hour (or >= 8 across three consecutive hours) AND its
members span >= 2 distinct 10-minute bins (burst guard).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402


def _components(P, merge_s=W.NODE_MERGE_S):
    """Union-find over the pair list; returns (node_times, labels)."""
    t = np.unique(np.concatenate([P.t_i.to_numpy(float), P.t_j.to_numpy(float)]))
    t.sort()
    # merge node times within merge_s (single linkage along the sorted axis)
    grp = np.zeros(len(t), int)
    g = 0
    for k in range(1, len(t)):
        if t[k] - t[k - 1] > merge_s:
            g += 1
        grp[k] = g
    node_t = np.array([t[grp == k].mean() for k in range(g + 1)])
    idx_of = {}
    for k in range(len(t)):
        idx_of[round(t[k], 3)] = grp[k]
    parent = list(range(g + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(P.t_i.to_numpy(float), P.t_j.to_numpy(float)):
        ia, ib = idx_of[round(a, 3)], idx_of[round(b, 3)]
        if ia != ib:
            parent[find(ia)] = find(ib)
    lab = np.array([find(k) for k in range(g + 1)])
    return node_t, lab


def components_of(P, merge_s=W.NODE_MERGE_S):
    """[(times array)] per connected component, computed day by day."""
    out = []
    if not len(P):
        return out
    for d, g in P.groupby("date"):
        nt, lab = _components(g, merge_s)
        for k in np.unique(lab):
            out.append((d, np.sort(nt[lab == k])))
    return out


def components_detail(P, merge_s=W.NODE_MERGE_S):
    """[(date, node times, max pair net_mad, n_pairs)] per connected component, day by day.

    Same nodes/edges/components as `components_of`; it additionally carries each component's
    STRONGEST pair, which is the ranking key of the 1c' pair-seed rule (Merlin F3)."""
    out = []
    if not len(P):
        return out
    for d, g in P.groupby("date"):
        nt, lab = _components(g, merge_s)
        ti = g.t_i.to_numpy(float)
        tj = g.t_j.to_numpy(float)
        t = np.unique(np.concatenate([ti, tj]))
        t.sort()
        grp = np.zeros(len(t), int)
        gg = 0
        for k in range(1, len(t)):
            if t[k] - t[k - 1] > merge_s:
                gg += 1
            grp[k] = gg
        idx_of = {round(float(x), 3): int(grp[k]) for k, x in enumerate(t)}
        comp_pair = np.array([lab[idx_of[round(float(x), 3)]] for x in ti])
        nm = g.net_mad.to_numpy(float)
        for k in np.unique(lab):
            m = comp_pair == k
            out.append((d, np.sort(nt[lab == k]),
                        float(nm[m].max()) if m.any() else float("nan"), int(m.sum())))
    return out


def qualifies(times, bin_s=W.SEED_BIN_S):
    """Merlin's seed-cluster clause + burst guard.  Returns (bool, n_hour, n_3h, n_bins)."""
    t = np.asarray(times, float)
    if len(t) < 2:
        return False, len(t), len(t), 1
    hh = np.floor(t / 3600.0).astype(int)
    n_hour = max(int((hh == h).sum()) for h in np.unique(hh))
    n_3h = max(int(((hh >= h - 1) & (hh <= h + 1)).sum()) for h in np.unique(hh))
    nb = len(np.unique(np.floor(t / bin_s).astype(int)))
    ok = ((n_hour >= W.SEED_MIN_MEMBERS) or (n_3h >= W.SEED_MIN_MEMBERS_PM1H)) \
        and nb >= W.SEED_MIN_BINS
    return bool(ok), n_hour, n_3h, nb


def threshold(quiet_pairs, n_quiet_hours, out=None):
    """The FA-only rule.  Returns (T, table, binds)."""
    P = quiet_pairs
    grid = np.arange(W.T_ENG_FLOOR, W.T_MAX + 1e-9, W.T_GRID)
    rows = []
    for T in grid:
        Q = P[P.net_mad >= T]
        fa = len(Q) / max(n_quiet_hours, 1)
        ncomp = 0
        if len(Q):
            comps = components_of(Q)
            ncomp = int(sum(1 for _, t in comps if len(t) >= W.QUIET_MAX_COMPONENT))
        rows.append(dict(T=float(T), n_pairs=int(len(Q)), fa_per_hour=float(fa),
                         n_components_ge5=ncomp))
        if fa <= W.QUIET_FA_PAIRS_PER_HOUR and ncomp == 0:
            break
    Tt = pd.DataFrame(rows)
    ok = Tt[(Tt.fa_per_hour <= W.QUIET_FA_PAIRS_PER_HOUR) & (Tt.n_components_ge5 == 0)]
    T = float(ok["T"].iloc[0]) if len(ok) else float(W.T_MAX)
    binds = "engineering_floor" if (len(ok) and T <= W.T_ENG_FLOOR + 1e-9) else "fa_rule"
    if out:
        Tt.to_csv(out + "_grid.csv", index=False)
        json.dump(dict(T=T, binds=binds, n_quiet_hours=int(n_quiet_hours),
                       engineering_floor=W.T_ENG_FLOOR,
                       fa_allowed_per_quiet_hour=W.QUIET_FA_PAIRS_PER_HOUR,
                       quiet_max_component=W.QUIET_MAX_COMPONENT,
                       saturated=bool(not len(ok))), open(out + ".json", "w"), indent=1)
    W.log(f"[T] T = {T:.2f} ({binds}); {n_quiet_hours} quiet hours")
    return T, Tt, binds


def cluster(pairs, T, out):
    """Pairs -> seed clusters.  Writes <out>_clusters.parquet (cid, t) and <out>_seeds.csv.

    SEED RULE (W_SEED_RULE).
      "frozen"  the 1c rule: >= 5 members inside one hour (or >= 8 across +- 1 h) AND >= 2 distinct
                ten-minute bins; ranked by size; top W_SEED_BUDGET.
      "pairs2"  Merlin's 1c' repair (F3): EVERY component with >= 2 nodes is a seed, ranked by its
                maximum pair net_mad then by size; top W_SEED_BUDGET.  The >= 5-member and >= 2-bin
                clauses are removed for that run -- the window-transient veto replaces the burst
                guard and the growth/survival stage carries the false-alarm control.
    """
    P = pairs[pairs.net_mad >= T].copy()
    pair_seed = (W.SEED_RULE == "pairs2")
    comps = components_detail(P) if pair_seed else [(d, t, np.nan, 0)
                                                    for d, t in components_of(P)]
    rows, mem = [], []
    cid = 0
    for d, t, mx, npr in comps:
        if pair_seed:
            if len(t) < 2:
                continue
            ok, n_hour, n_3h, nb = True, *qualifies(t)[1:]
        else:
            ok, n_hour, n_3h, nb = qualifies(t)
        if not ok:
            continue
        cid += 1
        c = f"C{cid:05d}"
        rows.append(dict(cid=c, date=d, n_members=int(len(t)), n_in_hour=n_hour, n_in_3h=n_3h,
                         n_bins=nb, t0=float(t.min()), t1=float(t.max()),
                         t_seed=float(np.median(t)), max_net_mad=float(mx), n_pairs=int(npr)))
        mem.append(pd.DataFrame(dict(cid=c, t=t)))
    SEED_COLS = ["cid", "date", "n_members", "n_in_hour", "n_in_3h", "n_bins", "t0", "t1",
                 "t_seed", "max_net_mad", "n_pairs"]
    Cl = pd.DataFrame(rows, columns=SEED_COLS) if rows else pd.DataFrame(columns=SEED_COLS)
    if len(Cl):
        Cl = (Cl.sort_values(["max_net_mad", "n_members"], ascending=False) if pair_seed
              else Cl.sort_values("n_members", ascending=False)).head(W.SEED_BUDGET
                                                                     ).reset_index(drop=True)
        keep = set(Cl.cid)
        M = pd.concat([m for m in mem if m.cid.iloc[0] in keep], ignore_index=True)
    else:
        M = pd.DataFrame(columns=["cid", "t"])
    Cl.to_csv(out + "_seeds.csv", index=False)
    M.to_parquet(out + "_clusters.parquet", index=False)
    W.log(f"[cluster] {len(P)} pairs >= T -> {len(comps)} components -> {len(Cl)} seed clusters "
          f"(rule {W.SEED_RULE}, budget {W.SEED_BUDGET})")
    return Cl, M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["threshold", "cluster"])
    ap.add_argument("--quiet-pairs", default="")
    ap.add_argument("--pairs", default="")
    ap.add_argument("--n-quiet-hours", type=int, default=0)
    ap.add_argument("--thr", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if a.cmd == "threshold":
        P = pd.read_parquet(a.quiet_pairs)
        nq = a.n_quiet_hours or json.load(open(a.quiet_pairs.replace(".parquet", "_meta.json")))["n_hours"]
        threshold(P, nq, out=a.out)
    else:
        T = json.load(open(a.thr))["T"]
        cluster(pd.read_parquet(a.pairs), T, a.out)


if __name__ == "__main__":
    main()
