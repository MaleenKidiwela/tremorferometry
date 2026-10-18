#!/usr/bin/env python
"""FINAL SPEC, STEP 4 -- results/freeze_step1c_prime_v3.json, a NEW file.

Written BEFORE the one blind run (Step 5) and never mutated (Merlin F8).  It records the order
Step 3 SELECTED, the ETS-day-only counting (F5), the continuation constants, the Step-2 reference,
the unchanged bars, and the sha256 of every script and of every lfe_picker module the engine
imports.  `freeze_step1c_prime.json` and `_v2.json` stay on disk unchanged.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import q_common as Q                                                        # noqa: E402
import z_score as ZS                                                        # noqa: E402

PICK = "/home/jovyan/tremorferometry/lfe_picker/scripts"
IMPORTED = ["v_common.py", "v_engine.py", "v_bank.py", "v_iterate.py", "v_dedupe.py",
            "g_common.py", "g_gates.py", "g_locate.py", "g_rediscover.py", "q_rediscover.py",
            "p_common.py", "s_common.py", "p_vel1d.py"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/freeze_step1c_prime_v3.json")
    ap.add_argument("--run", default=f"{W.RES}/step1c_prime")
    ap.add_argument("--final", default=Q.OUT)
    a = ap.parse_args()
    import s_common as S
    hid, iso, _ = S.holdouts()
    assert sorted(iso) == W.ISOLATED_HIDDEN
    for p in (f"{W.RES}/freeze_step1c.json", f"{W.RES}/freeze_step1c_prime.json",
              f"{W.RES}/freeze_step1c_prime_v2.json"):
        assert os.path.exists(p), f"the earlier freeze {p} must stay on disk (F8)"

    s1 = Q.jload(f"{a.final}/s1_quantisation.json", {})
    s2 = Q.jload(f"{a.final}/oracle_reference_2005_v3.json", {})
    s3 = Q.jload(f"{a.final}/s3/s3_verdict.json", {})
    selected = (s3 or {}).get("selected_order")
    assert selected in ("frozen", "staged"), f"Step 3 selected no order ({selected})"

    git = ""
    try:
        git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(HERE),
                             capture_output=True, text=True).stdout.strip()
    except Exception:                                                       # noqa: BLE001
        pass
    consts = {k: getattr(W, k) for k in dir(W)
              if k.isupper() and not k.startswith("_")
              and isinstance(getattr(W, k), (int, float, str, tuple, list, bool))}
    qconsts = {k: getattr(Q, k) for k in dir(Q)
               if k.isupper() and not k.startswith("_")
               and isinstance(getattr(Q, k), (int, float, str, tuple, list, bool))}
    scripts = (sorted(glob.glob(f"{HERE}/q_*.py")) + sorted(glob.glob(f"{HERE}/q_*.sh"))
               + sorted(glob.glob(f"{HERE}/z_*.py")) + W.w_scripts())
    seeds = W.read_csv_safe(f"{a.run}/cl_seeds.csv")
    T = Q.jload(f"{a.run}/T.json", {})

    F = dict(
        created_utc=pd.Timestamp.now("UTC").isoformat(),
        study="lfe_array_discovery / FINAL SPECIFICATION of the north positive control "
              "(PLAN.md sec.12; results/merlin_ruling_step1cprime_2026-09-20.md)",
        supersedes="results/freeze_step1c_prime.json and _v2.json -- NOT edited; the 1c' stages "
                   "1-5 they cover are reused verbatim (station-hour veto, quiet hours, "
                   "ac.parquet, T_pair = 8.0)",
        status="LAST CHANGE OF ANY KIND (Merlin's ruling): if anything is found wanting after "
               "Step 5, STOP.",
        git_head=git, bed_mode=W.BED_MODE, bed_config=W.bed_cfg(),
        reused_from_1c_prime=dict(
            station_hour_veto=f"{a.run}/station_hour_veto.csv",
            station_hour_veto_sha256=W.sha256_file(f"{a.run}/station_hour_veto.csv"),
            quiet_hours=f"{a.run}/quiet_hours.csv",
            quiet_hours_sha256=W.sha256_file(f"{a.run}/quiet_hours.csv"),
            autocorrelation=f"{a.run}/ac.parquet",
            autocorrelation_sha256=W.sha256_file(f"{a.run}/ac.parquet"),
            T_pair=T, flags_v1=Q.jload(f"{W.RES}/freeze_step1c_prime.json", {}).get("repair", {})
            .get("flags")),
        corrections=dict(
            F5="SURVIVAL COUNTS ETS-WINDOW DAYS ONLY (W_SURV_ETS_ONLY).  w_common.bed_days() "
               "returns all 44 cached days, so the 20 FA-calibration quiet days could supply two "
               "of the '>= 3 distinct days'.  n_indep / n_days / frac_near_seed / the scramble "
               "ratio are now computed on the ETS subset; quiet-day independent hits are REPORTED "
               "as the empirical false-alarm count (n_indep_quiet).  _indep.parquet is unchanged.",
            oracle="the oracle is re-seeded from STEP-1-CORRECTED times spread over the whole "
                   "24-day window (engine cap MAX_STACK = %d), so the near-seed burst guard "
                   "cannot bind; the 1c' v2 'burst guard reported not applied' fix is RETIRED "
                   "(unnecessary once the seed is not a burst)." % W.MAX_STACK,
            quantisation="the 2004-2005 rows of the Bostock export carry whole-second times; the "
                         "2005 event times are re-aligned by the cc peak of the 2012 sub-second "
                         "stack within +- %.1f s at the shared permanents %s"
                         % (Q.CC_LAG_S, list(Q.S1_SHARED)),
            scramble_missing_zero="W_SCRAMBLE_MISSING_ZERO: a template whose twin was dropped "
                                  "upstream has n_scramble = 0, not NaN.  Used ONLY by the staged "
                                  "order's final survival, where the continuation filter "
                                  "legitimately kills twins."),
        selected_order=selected,
        order_definition=dict(
            frozen="handoff -> quiet-day T_seed -> growth -> template-scramble twin -> SURVIVAL "
                   "GATE -> iterate 1 -> quiet -> grow -> survive -> iterate 2 -> quiet -> grow "
                   "-> survive -> dedupe -> locate/Gates 0-2",
            staged="handoff -> quiet-day T_seed -> pass-0 growth -> CONTINUATION FILTER (>= %d "
                   "independent hits on >= %d ETS days AND >= %g x the template-scramble twin; a "
                   "COMPUTE filter, not a claim) -> v_iterate(use_all, min_indep %d) -> growth -> "
                   "iterate -> growth -> THE FROZEN SURVIVAL GATE with every clause, ETS days "
                   "only -> dedupe -> locate/Gates 0-2; the scramble twins go through it "
                   "identically (their continuation filter is the count clause only)."
                   % (Q.S3_CONT_MIN_INDEP, Q.S3_CONT_MIN_DAYS, Q.S3_CONT_SCRAMBLE_RATIO,
                      Q.S3_ITER_MIN_INDEP),
            known_constraint="v_iterate re-stacks only templates with >= 10 independent hits "
                             "(lfe_picker/scripts/v_iterate.py, frozen and NOT edited), so "
                             "min_indep 5 binds only in the range >= 10."),
        step1_quantisation=s1,
        step2_reference=dict(path=f"{a.final}/oracle_reference_2005_v3.csv",
                             sha256=W.sha256_file(f"{a.final}/oracle_reference_2005_v3.csv")
                             if os.path.exists(f"{a.final}/oracle_reference_2005_v3.csv") else None,
                             **{k: s2.get(k) for k in
                                ("C_iso", "C_30", "fallback_active", "pass_bar", "denominators",
                                 "median_recall_all", "median_recall_indep", "median_recall_loud",
                                 "median_purity")}),
        step3_bootstrap={k: (s3 or {}).get(k) for k in
                         ("verdict", "selected_order", "bar", "K", "K_families", "K_isolated",
                          "n_matched", "n_chance", "chance_fraction", "twin_survivors",
                          "real_survivors_staged", "twin_fraction")},
        seeds=dict(rule=W.SEED_RULE, budget=W.SEED_BUDGET, n=int(len(seeds)),
                   path=f"{a.run}/cl_seeds.csv",
                   sha256=W.sha256_file(f"{a.run}/cl_seeds.csv")
                   if os.path.exists(f"{a.run}/cl_seeds.csv") else None,
                   min_nodes=2, rank="max pair net_mad, then component size",
                   T_pair=(T or {}).get("T")),
        bars=dict(pass_isolated_ge=W.ACCEPT_REDISCOVERED,
                  centroid_km=W.SCORE_R_KM, recall_frac=W.SCORE_RECALL_FRAC,
                  purity_frac=W.SCORE_PURITY_FRAC, max_explained=W.SCORE_MAX_EXPLAINED,
                  fallback=(s2 or {}).get("pass_bar"),
                  split_le=W.ACCEPT_SPLIT_RATE, null_gated_le=W.ACCEPT_NULL_GATED,
                  null_survivors_le_frac=ZS.NULL_SURV_PASS_FRAC,
                  fa_broken_ge_frac=ZS.NULL_SURV_BROKEN_FRAC,
                  fail_below_frac_of_ceiling=ZS.FAIL_FRAC_OF_CEILING,
                  stop_rules="Step 5 FAIL (< 2/10 or < 0.4 x C_iso with clean nulls) or FA_BROKEN "
                             "(null survivors >= 50 % of real) -> STOP Arm A north AND WA; the "
                             "study closes.  PASS -> the WA_GO conditions of the first audit "
                             "apply with v3 in place of 1c'."),
        null=dict(realisations=Q.S5_NULL_REALISATIONS, kind="whole-chain site-day shift, members "
                  "shift together", counts="gated families AND survivors per realisation",
                  seed=W.SITESHIFT_SEED, min_days=W.SITESHIFT_MIN_DAYS),
        isolated_hidden=W.ISOLATED_HIDDEN,
        isolated_hidden_sha256=W.sha256_of(W.ISOLATED_HIDDEN),
        hidden_sha256=W.sha256_of(sorted(int(x) for x in hid)),
        constants=dict(engine=consts, final_spec=qconsts),
        flags=W.repair_flags(),
        sha256=dict(scripts={os.path.basename(p): W.sha256_file(p) for p in scripts},
                    imported_lfe_picker={f: W.sha256_file(f"{PICK}/{f}") for f in IMPORTED
                                         if os.path.exists(f"{PICK}/{f}")}),
        step6="ONLY after a PASS: the 2004-07 re-fire confirmation (~900 station-days). A "
              "STATISTIC, never a pass condition. NOT run by this queue.",
        wa="results/WA_GO is NEVER created by this queue, whatever the verdict.")
    json.dump(F, open(a.out, "w"), indent=1, default=str)
    Q.log(f"[freeze/v3] {a.out}: order {selected}; C_30 {s2.get('C_30')} / C_iso "
          f"{s2.get('C_iso')}; bar {(s2.get('pass_bar') or {}).get('value')}; "
          f"seeds {len(seeds)}; flags {W.repair_flags()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
