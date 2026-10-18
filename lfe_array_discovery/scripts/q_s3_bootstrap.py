#!/usr/bin/env python
"""FINAL SPEC, STEP 3 (GPU) -- THE TRUTH-SEEDED BOOTSTRAP-CONVERGENCE TEST.

The study's actual question (Merlin, PLAN sec.12 step 3): CAN A TRUE SAME-FAMILY PAIR SEED
BOOTSTRAP TO A FAMILY HERE?  Answered with truth-SELECTED seeds and no blind result read.

  1. Take the stage-7 seed set exactly as frozen in 1c' (every >= 2-node component of pairs
     >= T_pair = 8.0, top 2,000 by max net_mad).
  2. Label a component a TRUE SEED of hidden family k when >= 2 of its nodes fall within +- 2.5 s
     of `t_origin + median fam_tt - 1 s` of >= 2 DISTINCT catalogued events of k (D3's matching,
     y_step1diag.py ~550).  The CHANCE count comes from a +300 s catalogue shift.
     ABANDON if chance-matched >= 30 % of matched (the truth labels are not meaningful).
  3. K = hidden families with >= 1 true seed.  Run the true seeds through
       (A) the FROZEN order and (B) the STAGED order (q_chain), scramble twins through B
           identically.
  ACCEPT-A  rediscovered >= max(3, 0.3 K) under the frozen order -> SELECT the frozen order.
  ACCEPT-B  A fails AND B rediscovers >= max(3, 0.3 K) INCLUDING >= 1 isolated family AND the
            scramble twins yield <= 10 % as many survivors -> SELECT the staged order.
  ABANDON   neither, or twins >= 50 % of true-seed survivors under B (an FA amplifier).

Scoring uses the SAME relative-clause scorer as the blind run (`q_rediscover.score`) against the
Step-2 reference `oracle_reference_2005_v3.csv`.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import q_common as Q                                                        # noqa: E402
import q_chain as QC                                                        # noqa: E402

ABANDON_RC = 20


# ------------------------------------------------------------------ scoring --------------------
def score_run(run_dir, oracle_csv, out_prefix):
    """The blind scorer's rediscovery count, applied to one S3 chain (no null clause here)."""
    import q_rediscover as QR
    import s_common as S

    fam = Q.family_ll()
    hidden, iso = Q.hidden_sets()
    B = Q.bostock_all()
    B = B[B.year == 2005]
    events = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    windows = [(pd.Timestamp(W.STEP1_DAY0).timestamp(),
                pd.Timestamp(W.STEP1_DAY1).timestamp() + 86400)]
    ORA = pd.read_csv(oracle_csv).set_index("bkey")
    ref = ORA[["recall_ref", "purity_ref"]].copy()
    C = (pd.read_parquet(f"{run_dir}/lg_clusters.parquet")
         if os.path.exists(f"{run_dir}/lg_clusters.parquet")
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    clusters = {c: dict(times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                        lon=float(g.lon.iloc[0])) for c, g in C.groupby("cid")} if len(C) else {}
    R, SP, summ = QR.score(clusters, events, fam, sorted(hidden), ref, windows,
                           tol=W.SCORE_TOL_S, r_km=W.SCORE_R_KM,
                           recall_frac=W.SCORE_RECALL_FRAC, purity_frac=W.SCORE_PURITY_FRAC,
                           max_explained=W.SCORE_MAX_EXPLAINED, split_frac=W.SCORE_SPLIT_FRAC)
    if len(R):
        R.to_csv(out_prefix + "_rediscovery.csv", index=False)
    if len(SP):
        SP.to_csv(out_prefix + "_splits.csv", index=False)
    n_all = int(R.rediscovered.sum()) if len(R) else 0
    RI = R[R.hidden.isin(sorted(iso))] if len(R) else R
    n_iso = int(RI.rediscovered.sum()) if len(RI) else 0
    fams = sorted(int(x) for x in R[R.rediscovered].hidden) if len(R) else []
    return dict(n_clusters=len(clusters), n_rediscovered=n_all, n_isolated_rediscovered=n_iso,
                rediscovered=fams, spurious_split_rate=summ.get("spurious_split_rate"))


# ------------------------------------------------------------------ the stage ------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=f"{W.RES}/step1c_prime")
    ap.add_argument("--dir", default=f"{Q.OUT}/s3")
    ap.add_argument("--oracle", default=f"{Q.OUT}/oracle_reference_2005_v3.csv")
    ap.add_argument("--label-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    hidden, iso = Q.hidden_sets()

    seeds = W.read_csv_safe(f"{a.run}/cl_seeds.csv")
    members = pd.read_parquet(f"{a.run}/cl_clusters.parquet")
    if not len(seeds):
        raise SystemExit(f"{a.run}/cl_seeds.csv is empty -- build the stage-7 seed set first")
    Q.log(f"[s3] stage-7 seed set: {len(seeds)} components, {len(members)} nodes "
          f"(rule {W.SEED_RULE}, budget {W.SEED_BUDGET}, T_pair from {a.run}/T.json)")

    lab_p, chn_p = f"{a.dir}/true_seeds.csv", f"{a.dir}/chance_seeds.csv"
    if os.path.exists(lab_p) and os.path.exists(chn_p):
        LAB, CHA = pd.read_csv(lab_p), pd.read_csv(chn_p)
    else:
        LAB = Q.label_true_seeds(members, seeds, hidden, shift_s=0.0)
        CHA = Q.label_true_seeds(members, seeds, hidden, shift_s=Q.S3_CHANCE_SHIFT_S)
        LAB.to_csv(lab_p, index=False)
        CHA.to_csv(chn_p, index=False)
    n_matched = int(LAB.cid.nunique()) if len(LAB) else 0
    n_chance = int(CHA.cid.nunique()) if len(CHA) else 0
    K_fams = sorted(int(x) for x in LAB.bkey.unique()) if len(LAB) else []
    K = len(K_fams)
    K_iso = sorted(set(K_fams) & set(iso))
    chance_frac = (n_chance / n_matched) if n_matched else np.inf
    Q.log(f"[s3] true seeds {n_matched} components over K = {K} hidden families "
          f"({len(K_iso)} isolated); chance (+{Q.S3_CHANCE_SHIFT_S:.0f} s) {n_chance} "
          f"= {chance_frac:.3f} of matched")
    lab = dict(n_seeds=int(len(seeds)), n_matched=n_matched, n_chance=n_chance,
               chance_fraction=float(chance_frac) if np.isfinite(chance_frac) else None,
               K=K, K_families=K_fams, K_isolated=K_iso,
               bar=int(max(Q.S3_ACC_MIN, int(np.ceil(Q.S3_ACC_FRAC * K)))),
               rule=dict(tol_s=Q.S3_MATCH_TOL_S, min_nodes=Q.S3_MIN_NODES,
                         axis="t_origin + median fam_tt - 1 s",
                         chance_shift_s=Q.S3_CHANCE_SHIFT_S,
                         abandon_chance_frac=Q.S3_ABANDON_CHANCE_FRAC))
    Q.jdump(lab, f"{a.dir}/labelling.json")
    if a.label_only:
        return 0
    if n_matched == 0 or chance_frac >= Q.S3_ABANDON_CHANCE_FRAC:
        Q.marker(f"{Q.OUT}/S3_ABANDON.txt",
                 f"S3: truth labels not meaningful -- {n_chance} chance-matched components from a "
                 f"+{Q.S3_CHANCE_SHIFT_S:.0f} s catalogue shift against {n_matched} matched "
                 f"({chance_frac:.3f} >= {Q.S3_ABANDON_CHANCE_FRAC}).")
        Q.jdump(dict(verdict="ABANDON", reason="truth labels not meaningful", **lab),
                f"{a.dir}/s3_verdict.json")
        return ABANDON_RC

    # ---- the true-seed subset, written as its own stage-7 hand-off -------------------------
    keep = sorted(set(LAB.cid))
    Sd = seeds[seeds.cid.isin(keep)].copy()
    Sd = Sd.merge(LAB.groupby("cid").bkey.min().rename("true_bkey"), on="cid", how="left")
    M = members[members.cid.isin(keep)].copy()
    Sd.to_csv(f"{a.dir}/cl_seeds.csv", index=False)
    M.to_parquet(f"{a.dir}/cl_clusters.parquet", index=False)

    res = {}
    for order in ("frozen", "staged"):
        d = f"{a.dir}/{order}"
        r = QC.run(d, f"{a.dir}/cl_seeds.csv", f"{a.dir}/cl_clusters.parquet", order=order)
        sc = score_run(d, a.oracle, f"{a.dir}/{order}_score")
        r.update(sc)
        res[order] = r
        Q.jdump(r, f"{a.dir}/{order}_result.json")
        Q.log(f"[s3/{order}] survivors {r['n_survivors_final']}, located {r['n_located']}, "
              f"rediscovered {r['n_rediscovered']} ({r['n_isolated_rediscovered']} isolated): "
              f"{r['rediscovered']}")

    bar = int(max(Q.S3_ACC_MIN, int(np.ceil(Q.S3_ACC_FRAC * K))))
    A, Bo = res["frozen"], res["staged"]
    twin_n = int(Bo.get("n_twin_final") or 0)
    real_n = int(Bo.get("n_survivors_final") or 0)
    twin_frac = (twin_n / real_n) if real_n > 0 else (np.inf if twin_n else 0.0)
    acc_A = A["n_rediscovered"] >= bar
    acc_B = (Bo["n_rediscovered"] >= bar and Bo["n_isolated_rediscovered"] >= 1
             and twin_frac <= Q.S3_ACC_TWIN_FRAC)
    fa_amp = twin_frac >= Q.S3_ABANDON_TWIN_FRAC

    if acc_A:
        verdict, selected = "ACCEPT_A", "frozen"
    elif fa_amp:
        verdict, selected = "ABANDON_FA_AMPLIFIER", None
    elif acc_B:
        verdict, selected = "ACCEPT_B", "staged"
    else:
        verdict, selected = "ABANDON_NO_CONVERGENCE", None

    out = dict(verdict=verdict, selected_order=selected, bar=bar, **lab,
               frozen=A, staged=Bo,
               twin_survivors=twin_n, real_survivors_staged=real_n,
               twin_fraction=float(twin_frac) if np.isfinite(twin_frac) else None,
               accept=dict(rediscovered_ge=bar, twin_le_frac=Q.S3_ACC_TWIN_FRAC,
                           abandon_twin_ge_frac=Q.S3_ABANDON_TWIN_FRAC,
                           b_needs_isolated=1),
               note="a CONTINUATION FILTER is a compute filter, not a claim; every threshold, "
                    "clause and bar of the frozen chain is unchanged in order B -- only the ORDER "
                    "of frozen operations differs, and the scramble twins go through it "
                    "identically")
    Q.jdump(out, f"{a.dir}/s3_verdict.json")
    Q.log(f"[s3] VERDICT {verdict} (bar {bar} of K = {K}); frozen {A['n_rediscovered']}, "
          f"staged {Bo['n_rediscovered']} ({Bo['n_isolated_rediscovered']} isolated), "
          f"twins {twin_n} vs real {real_n} ({twin_frac:.3f})")
    if selected is None:
        Q.marker(f"{Q.OUT}/S3_ABANDON.txt",
                 ("S3: FA amplifier -- the staged order's scramble twins yield "
                  f"{twin_n} survivors against {real_n} real ({twin_frac:.2f} >= "
                  f"{Q.S3_ABANDON_TWIN_FRAC}); the staged order must be forbidden."
                  if fa_amp else
                  f"S3: bootstrap does not converge -- with truth-SELECTED same-family pair seeds "
                  f"(K = {K} hidden families, bar {bar}) the frozen order rediscovers "
                  f"{A['n_rediscovered']} and the staged order {Bo['n_rediscovered']} "
                  f"({Bo['n_isolated_rediscovered']} isolated, twins {twin_n}/{real_n})."))
        return ABANDON_RC
    Q.marker(f"{Q.OUT}/SELECTED_ORDER.txt", selected)
    return 0


if __name__ == "__main__":
    sys.exit(main())
