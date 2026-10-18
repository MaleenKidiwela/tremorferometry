#!/usr/bin/env python
"""Step 0 -- FREEZE.  Writes results/freeze_step1.json BEFORE the queue is launched.

Every constant of the engine, the bed's station/day availability table, the quiet-hour list, the
candidate-hour list, the data-driven site table, the isolated-hidden-family list hash, and the
sha256 of every w_ script and of every lfe_picker script the engine imports.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

PICK = "/home/jovyan/tremorferometry/lfe_picker/scripts"
IMPORTED = ["v_common.py", "v_engine.py", "v_bank.py", "v_iterate.py", "v_dedupe.py",
            "g_common.py", "g_gates.py", "g_locate.py", "g_rediscover.py", "q_rediscover.py",
            "p_common.py", "s_common.py", "p_vel1d.py"]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/freeze_step1.json")
    ap.add_argument("--tag", default="step1")
    ap.add_argument("--deviation", default="")
    ap.add_argument("--spec", action="store_true",
                    help="freeze the SPEC before the bed exists (no availability / hour lists)")
    a = ap.parse_args()
    import s_common as S
    hid, iso, _ = S.holdouts()
    if a.spec:
        consts = {k: getattr(W, k) for k in dir(W)
                  if k.isupper() and not k.startswith("_")
                  and isinstance(getattr(W, k), (int, float, str, tuple, list))}
        F = dict(created_utc=pd.Timestamp.now("UTC").isoformat(),
                 study=f"lfe_array_discovery / {a.tag} -- SPEC frozen BEFORE the bed is downloaded",
                 bed_mode=W.BED_MODE, deviation=a.deviation,
                 bed_config=W.bed_cfg(),
                 isolated_hidden=W.ISOLATED_HIDDEN,
                 isolated_hidden_sha256=W.sha256_of(W.ISOLATED_HIDDEN),
                 hidden_sha256=W.sha256_of(sorted(int(x) for x in hid)),
                 accept=dict(rediscovered_ge=W.ACCEPT_REDISCOVERED,
                             spurious_split_le=W.ACCEPT_SPLIT_RATE,
                             null_gated_le=W.ACCEPT_NULL_GATED),
                 constants=consts,
                 note="availability, quiet-hour and candidate-hour lists are written by the queue "
                      "AFTER the download, into <tag>_bed_availability.csv / _quiet_hours.csv / "
                      "_candidate_hours.csv and freeze_<tag>_full.json; no engine constant changes",
                 sha256=dict(
                     w_scripts={os.path.basename(p): W.sha256_file(p) for p in W.w_scripts()},
                     imported_lfe_picker={f: W.sha256_file(f"{PICK}/{f}") for f in IMPORTED
                                          if os.path.exists(f"{PICK}/{f}")}))
        json.dump(F, open(a.out, "w"), indent=1, default=str)
        print(f"[freeze/spec] {a.out}")
        print(f"  bed {W.BED_MODE}: sites {W.bed_cfg().get('sites')} withheld "
              f"{W.bed_cfg().get('withheld')}")
        return
    mem, of, order, det = W.sites()
    A = W.availability()
    QD = W.quiet_days()
    QH, n_quiet_days = W.quiet_hours()
    CH = W.candidate_hours()
    assert sorted(iso) == W.ISOLATED_HIDDEN
    ST = W.station_table()
    site_tbl = [dict(site=q, members=mem[q], lat=float(ST.lat[mem[q][0]]),
                     lon=float(ST.lon[mem[q][0]])) for q in order]
    consts = {k: getattr(W, k) for k in dir(W)
              if k.isupper() and not k.startswith("_")
              and isinstance(getattr(W, k), (int, float, str, tuple, list))}
    git = ""
    try:
        git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(HERE),
                             capture_output=True, text=True).stdout.strip()
    except Exception:                                                       # noqa: BLE001
        pass
    F = dict(
        created_utc=pd.Timestamp.utcnow().isoformat(),
        study=f"lfe_array_discovery / {a.tag} -- north positive control of the Arm-A engine",
        bed_mode=W.BED_MODE,
        deviation=a.deviation,
        plan="PLAN.md sec.10 (supersedes sec.1-9); results/merlin_review_plan_2026-09-20.md "
             "Q1/(a), F3, Step 0, Step 1, Q5/(b), Honest priors",
        git_head=git,
        bed=dict(root=W.BED, year=W.STEP1_YEAR, day0=W.STEP1_DAY0, day1=W.STEP1_DAY1,
                 n_days=int(len(A)), n_stations_on_bed=len(W.bed_stations()),
                 detection_stations=det, withheld=list(W.WITHHELD),
                 n_station_days=int(A.n_sta.sum()),
                 availability=A.to_dict("records")),
        sites=dict(merge_km=W.SITE_MERGE_KM, rule="1-km single linkage on the station table",
                   table=site_tbl, min_sites=W.MIN_SITES,
                   note="per-day availability is NEVER intersected"),
        prerequisites=dict(
            P1_no_T8_floor=dict(engineering_floor=W.T_ENG_FLOOR, grid=W.T_GRID,
                                autocorr_rule="<= 1 false pair per quiet hour AND zero quiet-hour "
                                              "components of size >= 5",
                                growth_rule=f"<= {W.FA_PER_TPL_DAY} false hits per template-day",
                                binding_reported_in=["step1/T.json",
                                                     "step1/p*_thresholds.csv:floor_binds"]),
            P2_data_driven_sites=True,
            P3_station_z_valid_samples_only=True,
            P4_min_sites=W.MIN_SITES),
        autocorrelation=dict(
            window_s=W.AC_WIN_S, hop_s=W.AC_HOP_S, block_s=W.AC_BLOCK_S,
            same_absolute_time_no_moveout_model=True, lag_resolution_samples=1,
            min_abs_lag_s=W.AC_MIN_LAG_S, decluster_s=W.AC_DECLUSTER_S,
            per_station="normalised 3-C mean cc of window i against the whole hour -> lower-half "
                        "MAD z over the hour (valid samples only)",
            per_site=f"mean of member z, clipped at {W.AC_SITE_Z_CLIP}",
            network="equal-weight site sum, re-normalised over lag by its own median / lower-half "
                    "MAD (net_mad)",
            hit=f"net_mad >= T AND >= {W.AC_MIN_SITE_Z_SITES} sites with site-z >= "
                f"{W.AC_SITE_Z_COINC} AND >= {W.MIN_SITES} sites observable",
            robust_scale_subsample=W.AC_Z_SUB,
            gpu="torch.cuda on the H100 (DOCUMENTED DEVIATION from Merlin's 'cupy': the same "
                "device and the same fp32 matmul; torch is already the engine's runtime in "
                "v_engine, and mixing two GPU runtimes in one process was judged the larger risk)",
            candidate_hour=f">= {W.AC_HOUR_MIN_EPI} PNSN epicentres within {W.AC_HOUR_EPI_KM} km "
                           f"of any active cell",
            active_cells=dict(grid_deg=W.ACTIVE_CELL_DEG, min_epi=W.ACTIVE_CELL_MIN_EPI,
                              km=W.ACTIVE_CELL_KM, dep_km=list(W.ACTIVE_CELL_DEP),
                              source="PNSN 2010-2013, Slab2; NO family coordinate is used"),
            n_candidate_hours=int(len(CH)), hour_budget=W.AC_HOUR_BUDGET,
            candidate_hours=[f"{r.date} {int(r.hour):02d}" for r in CH.itertuples()]),
        quiet=dict(rule=f"no PNSN tremor within {W.QUIET_TREMOR_KM} km and >= "
                        f"{W.QUIET_MIN_DAYS_FROM_TREMOR} days from any tremor day",
                   n_quiet_days_found=int(n_quiet_days),
                   merlin_asked_for=">= 20 quiet days x 2 hours = >= 40 quiet hours",
                   deviation="the 2012 bed has only 52 day-files, so >= 20 quiet DAYS is "
                             "impossible; the binding clause (>= 40 quiet hours on the same "
                             "sites) is met with more hours per quiet day",
                   n_quiet_hours=int(len(QH)),
                   quiet_days=QD[QD.quiet].date.tolist(),
                   quiet_hours=[f"{r.date} {int(r.hour):02d}" for r in QH.itertuples()]),
        clustering=dict(node_merge_s=W.NODE_MERGE_S, seed_min_members=W.SEED_MIN_MEMBERS,
                        seed_min_members_pm1h=W.SEED_MIN_MEMBERS_PM1H,
                        burst_guard_bins=W.SEED_MIN_BINS, bin_s=W.SEED_BIN_S,
                        seed_budget=W.SEED_BUDGET),
        handoff=dict(align_iters=W.ALIGN_ITERS, align_lag_samples=W.ALIGN_LAG,
                     cut_pre_s=W.HANDOFF_CUT_PRE_S, cut_len_s=W.HANDOFF_CUT_LEN_S,
                     onset_search_s=list(W.HANDOFF_ONSET_SEARCH),
                     tpl0_len_s=W.TPL0_LEN_S, tpl0_pre_s=W.TPL0_PRE_S,
                     tpl1_len_s=W.TPL1_LEN_S, tpl1_pre_s=W.TPL1_PRE_S,
                     n_iterations=W.N_ITERATIONS),
        survival=dict(min_independent=W.SURV_MIN_INDEP, independence_s=W.INDEP_S,
                      min_days=W.SURV_MIN_DAYS, max_frac_near_seed=W.SURV_MAX_NEAR_SEED,
                      near_seed_s=W.NEAR_SEED_S, scramble_ratio=W.SCRAMBLE_RATIO,
                      decluster_s=W.GROW_DECLUSTER_S),
        dedupe=dict(jaccard=W.DEDUPE_JACCARD, cc_min_sites=W.DEDUPE_CC_MIN_SITES,
                    tol_s=W.DEDUPE_TOL_S, lag_s=W.DEDUPE_LAG_S,
                    rule="inter-template cc above the cross-family p95 at >= 3 sites OR hit-set "
                         "Jaccard >= 0.5; keep the larger"),
        location=dict(grid_lat0=W.GRID_LAT0, grid_lon0=W.GRID_LON0,
                      halfwidth_km=W.GRID_HALFWIDTH, pick_search_s=W.LOC_PICK_SEARCH_S,
                      min_pick_snr=W.LOC_MIN_PICK_SNR, n_pick_null=W.LOC_N_PICK_NULL,
                      n_synth=W.LOC_N_SYNTH, synth_bar_km=W.LOC_SYNTH_BAR_KM,
                      loso_bar_km=W.LOC_LOSO_BAR_KM, depth_bar_km=W.LOC_DEPTH_BAR_KM,
                      null_ratio=W.LOC_NULL_RATIO, min_sites=W.LOC_MIN_SITES,
                      max_gap_deg=W.LOC_MAX_GAP_DEG),
        gates=dict(gates=list(W.GATES), offset_search_s=W.GATES_OFFSET_SEARCH_S,
                   withheld=list(W.WITHHELD), control_table="north (g_gate1_controls.csv)",
                   gate3="NOT RUN -- untestable inside one ETS window (Merlin F6)"),
        nulls=dict(n_realisations=W.N_SITESHIFT_REALISATIONS, min_shift_days=W.SITESHIFT_MIN_DAYS,
                   seed=W.SITESHIFT_SEED,
                   rule="each SITE's day series circularly shifted by a site-specific whole-day "
                        "offset; members shift together"),
        scoring=dict(scorer="q_rediscover (RELATIVE clauses)",
                     reference="lfe_picker/results/p3/truth_template_reference.csv "
                               "(frozen template, 2012, SNB+B926 withheld)",
                     reference_deviation="the reference run withheld SNB+B926; this run also "
                                         "withholds VGZ, so the reference is if anything a "
                                         "HARDER bar than the engine's own network",
                     windows=W.SCORE_WINDOWS, tol_s=W.SCORE_TOL_S, r_km=W.SCORE_R_KM,
                     recall_frac=W.SCORE_RECALL_FRAC, purity_frac=W.SCORE_PURITY_FRAC,
                     cluster_truth_km=W.CLUSTER_TRUTH_KM,
                     isolated_hidden=W.ISOLATED_HIDDEN,
                     isolated_hidden_sha256=W.sha256_of(W.ISOLATED_HIDDEN),
                     n_hidden_total=len(hid),
                     hidden_sha256=W.sha256_of(sorted(int(x) for x in hid)),
                     ceiling="frozen template detector: 9/30 absolute, 5/10 isolated",
                     accept=dict(rediscovered_ge=W.ACCEPT_REDISCOVERED,
                                 spurious_split_le=W.ACCEPT_SPLIT_RATE,
                                 null_gated_le=W.ACCEPT_NULL_GATED),
                     failure=dict(engine_inadequate="< 2/10 with null <= 1",
                                  fa_rule_broken="null >= 3 per realisation")),
        constants=consts,
        sha256=dict(
            w_scripts={os.path.basename(p): W.sha256_file(p) for p in W.w_scripts()},
            imported_lfe_picker={f: W.sha256_file(f"{PICK}/{f}") for f in IMPORTED
                                 if os.path.exists(f"{PICK}/{f}")},
            truth_template_reference=W.sha256_file(
                "/home/jovyan/tremorferometry/lfe_picker/results/p3/truth_template_reference.csv"),
            holdouts_json=W.sha256_file(
                "/home/jovyan/tremorferometry/lfe_picker/results/holdouts.json")),
    )
    p = a.out
    json.dump(F, open(p, "w"), indent=1, default=str)
    A.to_csv(f"{W.RES}/{a.tag}_bed_availability.csv", index=False)
    QH.to_csv(f"{W.RES}/{a.tag}_quiet_hours.csv", index=False)
    CH.to_csv(f"{W.RES}/{a.tag}_candidate_hours.csv", index=False)
    print(f"[freeze] {p}")
    print(f"  bed {len(det)} detection stations / {len(order)} sites / {len(A)} days "
          f"({int(A.n_sta.sum())} station-days); withheld {list(W.WITHHELD)}")
    print(f"  quiet {n_quiet_days} days -> {len(QH)} quiet hours; candidate {len(CH)} tremor hours")
    print(f"  sites: " + "; ".join(f"{q}={'+'.join(mem[q])}" for q in order))


if __name__ == "__main__":
    main()
