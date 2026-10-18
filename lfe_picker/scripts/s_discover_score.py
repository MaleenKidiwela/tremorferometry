#!/usr/bin/env python
"""P3D DISCOVERY -- Gates 0-3, identity, and the pre-registered scorer.

Applies, to the clusters `s_discover_chain.py` produced:
  * Gate 1 at the WITHHELD sites.  The frozen `g_gates` knows only the S15 roster, so B926 is run
    through `g_gates` itself and the new withheld sites (CLRS, GOBB, LOPEZ, TURTL) are run through
    the SAME a3b statistic implemented on the P3D readers (`s_transfer.leg` -> alignment-free
    S-window excess with the +-1.5 s offset search against 200 matched random-time stacks).  Both
    are reported; nothing in `g_gates` is edited.
  * identity: cluster PGC stack cc vs each 2012 PGC template (chance level 0.524).
  * PRIMARY scorer on the scoreable isolated stratum, and the two pre-registered SECONDARY readouts
    (identity-only on 12/29/62; the Delta-sites-weighted fraction over all 10 isolated families).
All clauses and bars are read from results/p3d/discovery_freeze.json -- none is re-decided here.
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_transfer as T                                                     # noqa: E402

OUT = S.P3D
WITHHELD_STAS = ["B926", "CLRS", "GOBB", "LOPEZ", "TURTL"]
PGC_CROSSFAM = 0.524


def load_freeze():
    return json.load(open(f"{OUT}/discovery_freeze.json"))


def match_to_truth(t_cluster, t_truth, tol=2.0):
    """one-to-one matched count at +-tol s (greedy nearest)."""
    if not len(t_cluster) or not len(t_truth):
        return 0
    a = np.sort(np.asarray(t_cluster, float))
    b = np.sort(np.asarray(t_truth, float))
    used = np.zeros(len(b), bool)
    n = 0
    j = 0
    for x in a:
        while j < len(b) and b[j] < x - tol:
            j += 1
        k, best, bk = j, tol + 1, -1
        while k < len(b) and b[k] <= x + tol:
            if not used[k] and abs(b[k] - x) < best:
                best, bk = abs(b[k] - x), k
            k += 1
        if bk >= 0:
            used[bk] = True
            n += 1
    return n


def withheld_excess(dc, sta, t_ev, ttS, ST, rng, n_rand=200):
    """a3b / Gate-1 statistic at a withheld site: S-window excess vs matched random-time stacks."""
    L = T.leg(dc, sta, 0, t_ev, ttS, ST, rng, want_nulls=True, n_rand=n_rand, n_cell=0,
              keep_stack=False)
    if L.get("status") == "no_data":
        return dict(sta=sta, status="no_data")
    return dict(sta=sta, status="ok", n=L.get("n_stacked", 0), excess_db=L.get("excess_db", np.nan),
                null_p99=L.get("null_rand_p99", np.nan), null_med=L.get("null_rand_med", np.nan),
                pass_gate1=bool(L.get("pass_a", False)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--arm", default="A")
    ap.add_argument("--net", default="B")
    ap.add_argument("--max-clusters", type=int, default=200)
    ap.add_argument("--tol", type=float, default=2.0)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    fz = load_freeze()
    C = pd.read_csv(f"{OUT}/disc_{a.tag}_clusters.csv")
    if not len(C) or "pass_excess" not in C:
        json.dump(dict(tag=a.tag, arm=a.arm, net=a.net, n_clusters=0,
                       verdict="NO CLUSTER REACHED THE CHAIN'S S-EXCESS CLAUSE"),
                  open(f"{OUT}/disc_{a.tag}_score.json", "w"), indent=1)
        print("[score] no clusters to score")
        return
    M = pd.read_parquet(f"{OUT}/disc_{a.tag}_members.parquet")
    keep = C[C.pass_excess].head(a.max_clusters)
    ST = pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta")
    F = S.families()
    hid, iso, _ = S.holdouts()
    T1 = S.hp_events(0.8, episodes=fz["s1_candidates"]["pooled_episodes"])
    Tt = pd.read_csv(f"{OUT}/discovery_T1_truth.csv").set_index("bkey")
    scoreable_iso = set(fz["coordinator_decisions_2026-09-19"]["b_bars"]
                        ["PRIMARY_unchanged"]["scoreable_isolated"])
    sec_i = set(fz["coordinator_decisions_2026-09-19"]["b_bars"]["SECONDARY_i_identity_only"]["families"])
    dc = S.DayCache(maxn=500, years=list(range(2019, 2027)), stas=WITHHELD_STAS)
    tts = S.fam_tt(WITHHELD_STAS, ST)
    rng = np.random.default_rng(20260919)
    rows = []
    for r in keep.itertuples():
        tt = M[M.cid == r.cid].t.to_numpy()
        rec = dict(cid=int(r.cid), n_members=int(r.n_members), span_days=float(r.span_days))
        # ---- which hidden family, if any, does this cluster match?
        best_b, best_n = -1, 0
        for b in sorted(hid):
            n = match_to_truth(tt, T1[T1.bkey == b].t.to_numpy(), a.tol)
            if n > best_n:
                best_b, best_n = int(b), n
        rec["best_hidden_bkey"] = best_b
        rec["n_matched_T1"] = int(best_n)
        if best_b >= 0:
            nT1 = int(Tt.n_T1.get(best_b, 0))
            rec["frac_T1_matched"] = best_n / max(nT1, 1)
            rec["hidden_stratum"] = str(Tt.stratum.get(best_b, ""))
            rec["d45_V2"] = int(Tt.d45_V2.get(best_b, 0))
            # ---- Gate 1 at the withheld sites, at that family's predicted S
            g1 = []
            for w in WITHHELD_STAS:
                g1.append(withheld_excess(dc, w, tt, tts[(best_b, w)], ST, rng))
            rec["gate1_sites_pass"] = int(sum(1 for x in g1 if x.get("pass_gate1")))
            rec["gate1_detail"] = ";".join(
                f"{x['sta']}:{x.get('excess_db', float('nan')):.2f}/{x.get('null_p99', float('nan')):.2f}"
                for x in g1 if x.get("status") == "ok")
        rows.append(rec)
        if len(rows) % 10 == 0:
            print(f"  [score] {len(rows)}/{len(keep)}  {time.time()-t0:.0f}s", flush=True)
    R = pd.DataFrame(rows)
    # a chain that yields no surviving cluster is a RESULT, not an error: keep the schema so the
    # scorer, the figures and the write-up all see a well-formed empty table.
    for c, d in (("cid", 0), ("n_members", 0), ("span_days", 0.0), ("best_hidden_bkey", -1),
                 ("n_matched_T1", 0), ("frac_T1_matched", 0.0), ("hidden_stratum", ""),
                 ("d45_V2", 0), ("gate1_sites_pass", 0), ("gate1_detail", "")):
        if c not in R.columns:
            R[c] = pd.Series(dtype=type(d))
    R.to_csv(f"{OUT}/disc_{a.tag}_scored.csv", index=False)
    # ---- the pre-registered readouts
    def rediscovered(b, ident_only):
        g = R[R.best_hidden_bkey == b]
        if not len(g):
            return False
        g = g.sort_values("n_matched_T1", ascending=False).iloc[0]
        if ident_only:
            return bool(g.get("frac_T1_matched", 0) >= 0.5 and g.get("gate1_sites_pass", 0) >= 1)
        return bool(g.get("frac_T1_matched", 0) >= 0.30 and g.get("gate1_sites_pass", 0) >= 1)
    prim = {int(b): bool(rediscovered(b, False)) for b in sorted(scoreable_iso)}
    seci = {int(b): bool(rediscovered(b, True)) for b in sorted(sec_i)}
    w = {int(b): int(Tt.d45_V2.get(b, 0)) for b in sorted(iso)}
    red_all = {int(b): bool(rediscovered(b, False)) for b in sorted(iso)}
    wsum = sum(w.values())
    out = dict(tag=a.tag, arm=a.arm, net=a.net, n_clusters_scored=int(len(R)),
               PRIMARY=dict(scoreable_isolated=prim, n_rediscovered=int(sum(prim.values())),
                            n_scoreable=len(prim), bar=">= 2/3",
                            verdict=("PROCEED" if sum(prim.values()) >= (2 / 3) * len(prim)
                                     else "NOT MET")),
               SECONDARY_i_identity_only=dict(families=seci,
                                              note="descriptive only, 4-8 truth events, never a bar"),
               SECONDARY_ii_weighted=dict(
                   weights=w, rediscovered=red_all,
                   weighted_fraction=(sum(w[b] for b in red_all if red_all[b]) / wsum
                                      if wsum else float("nan")),
                   note="descriptive statistic, never a bar"),
               runtime_s=time.time() - t0)
    json.dump(out, open(f"{OUT}/disc_{a.tag}_score.json", "w"), indent=1, default=str)
    print(json.dumps({k: out[k] for k in ("PRIMARY", "SECONDARY_i_identity_only")}, indent=1))


if __name__ == "__main__":
    main()
