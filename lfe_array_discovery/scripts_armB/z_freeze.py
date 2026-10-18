#!/usr/bin/env python
"""STEP 1c' -- THE FREEZE.  Writes results/freeze_step1c_prime.json BEFORE stage 4 runs.

A NEW file: results/freeze_step1c.json is never edited (Merlin F8 -- version freezes, never mutate).
It records the hashes of every z_ and w_ script and of every lfe_picker module the engine imports,
the station-hour veto table, the hour-matched earthquake-screened quiet-hour list, the D8 verdict,
the oracle-reference PROTOCOL (its results are written later, by stage 6, and hashed into
oracle_reference_2005.json), the pair-seed rule, every repair-flag VALUE, and the pre-registered
bars of stage 10.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import z_score as ZS                                                        # noqa: E402
import z_oracle as ZO                                                       # noqa: E402
import z_quiet as ZQ                                                        # noqa: E402
import z_health as ZH                                                       # noqa: E402
import z_d8 as ZD                                                           # noqa: E402
import z_tpair as ZT                                                        # noqa: E402

PICK = "/home/jovyan/tremorferometry/lfe_picker/scripts"
IMPORTED = ["v_common.py", "v_engine.py", "v_bank.py", "v_iterate.py", "v_dedupe.py",
            "g_common.py", "g_gates.py", "g_locate.py", "g_rediscover.py", "q_rediscover.py",
            "p_common.py", "s_common.py", "p_vel1d.py"]


def jload(p):
    try:
        return json.load(open(p))
    except Exception:                                                       # noqa: BLE001
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/freeze_step1c_prime.json")
    ap.add_argument("--run", default=f"{W.RES}/step1c_prime")
    a = ap.parse_args()
    import s_common as S
    hid, iso, _ = S.holdouts()
    assert sorted(iso) == W.ISOLATED_HIDDEN
    assert os.path.exists(f"{W.RES}/freeze_step1c.json"), "the 1c freeze must stay on disk"

    mem, of, order, det = W.sites()
    A = W.availability()
    QH, nqd = W.quiet_hours()
    CH = W.candidate_hours()
    veto_csv = f"{a.run}/station_hour_veto.csv"
    quiet_csv = f"{a.run}/quiet_hours.csv"
    V = pd.read_csv(veto_csv) if os.path.exists(veto_csv) else pd.DataFrame()
    git = ""
    try:
        git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(HERE),
                             capture_output=True, text=True).stdout.strip()
    except Exception:                                                       # noqa: BLE001
        pass
    consts = {k: getattr(W, k) for k in dir(W)
              if k.isupper() and not k.startswith("_")
              and isinstance(getattr(W, k), (int, float, str, tuple, list, bool))}
    scripts = sorted(__import__("glob").glob(f"{HERE}/z_*.py")) + W.w_scripts()

    F = dict(
        created_utc=pd.Timestamp.now("UTC").isoformat(),
        study="lfe_array_discovery / STEP 1c' -- THE ONE PERMITTED REPAIR of the north positive "
              "control (pair-seeded growth + four defect fixes)",
        supersedes="results/freeze_step1c.json (NOT edited; 1c remains reproducible from it)",
        plan="PLAN.md sec.11 (binding) and sec.10 addenda; prescription steps 1-10 of "
             "results/merlin_review_step1_2026-09-20.md",
        git_head=git, bed_mode=W.BED_MODE, bed_config=W.bed_cfg(),
        repair=dict(
            defects_fixed=dict(
                F1="station-hour health veto (flat / H-V rms ratio outside [1/30,30] / 1-8 Hz "
                   "Welch PSD peak-over-median > 1e4), applied in the autocorrelation AND at "
                   "template stacking",
                F2="hour-matched quiet set (every clock hour >= 1, remainder by the candidate-hour "
                   "histogram) with a ComCat/NRCan earthquake screen (local M>=2.5 within 300 km, "
                   "global M>=5.5 with predicted P or 3.5 km/s surface arrival in the hour)",
                F3="PAIR-SEEDED GROWTH: every >= 2-node component is a seed, ranked by max pair "
                   "net_mad then size, budget 2000; the >= 5-member and >= 2-bin clauses removed; "
                   "handoff per-station minimum 3 -> 2 segments; near-seed clause against EITHER "
                   "end of the seed window",
                F4="scoring reference = the 2005 ORACLE (stage 6), never the transferred 2012 "
                   "truth-template table",
                F6="a station with zero hour-energy median is DROPPED from the hour and no NaN "
                   "can reach the site sum",
                F7="y_step1diag d6native: obspy UTCDateTime.timestamp is a property, not a call",
                F8="this file is a NEW freeze; freeze_step1c.json is not mutated"),
            window_transient_veto=dict(
                factor=W.WIN_TRANSIENT_FACTOR,
                rule="a 10-s window whose network-median energy (each site normalised by its own "
                     "hour median) exceeds the factor is excluded as row i AND as lag j, "
                     "identically in quiet and tremor hours"),
            flags=W.repair_flags()),
        d8=jload(f"{a.run}/d8.json"),
        station_hour_veto=dict(path=veto_csv,
                               sha256=W.sha256_file(veto_csv) if os.path.exists(veto_csv) else None,
                               n_rows=int(len(V)), n_vetoed=int(V.veto.sum()) if len(V) else 0,
                               census=jload(veto_csv.replace(".csv", ".json"))),
        quiet_hours=dict(path=quiet_csv,
                         sha256=W.sha256_file(quiet_csv) if os.path.exists(quiet_csv) else None,
                         n_hours=int(len(QH)), n_days=int(nqd),
                         list=QH.to_dict("records"),
                         screen=jload(quiet_csv.replace(".csv", ".json"))),
        candidate_hours=dict(n=int(len(CH)),
                             by_clock_hour={int(k): int(v)
                                            for k, v in CH.groupby("hour").size().items()}),
        availability=dict(n_days=int(len(A)), n_station_days=int(A.n_sta.sum()),
                          detection_stations=det, withheld=list(W.WITHHELD),
                          sites=list(order), records=A.to_dict("records")),
        T_pair_rule=dict(grid=W.T_GRID, floor=W.T_ENG_FLOOR,
                         fa_pairs_per_quiet_hour=W.QUIET_FA_PAIRS_PER_HOUR,
                         quiet_max_component=W.QUIET_MAX_COMPONENT,
                         accept_range=[ZT.T_LO, ZT.T_HI], loo_range_max=ZT.LOO_RANGE_MAX,
                         mandatory_leave_one_hour_out=True),
        seed_rule=dict(rule=W.SEED_RULE, budget=W.SEED_BUDGET, node_merge_s=W.NODE_MERGE_S,
                       min_nodes=2, rank="max pair net_mad, then component size",
                       removed_clauses=[">= 5 members in one hour", ">= 2 ten-minute bins"],
                       handoff_min_segments=W.HANDOFF_MIN_SEGMENTS,
                       unchanged=["align_stack", "onset search", "20-s pass-0 template",
                                  f">= {W.MIN_SITES} sites"],
                       accept="at least 500 seeds handed off with >= 4 sites"),
        oracle_protocol=dict(
            max_seed_events=ZO.MAX_SEED_EVENTS,
            seed="each of the 30 hidden families, its own catalogued events from its busiest hour",
            axis="t_origin + network-median predicted S delay - 1 s",
            align="align_stack +- 1.5 s (2005 catalog times are 1-s quantised)",
            chain="handoff -> quiet-day T_seed -> growth -> scramble -> survival -> "
                  f"{W.N_ITERATIONS} iterations -> dedupe -> locate/gates, all frozen",
            counts=f"centroid <= {W.SCORE_R_KM} km AND recall > 0 AND kept-explained <= "
                   f"{W.SCORE_MAX_EXPLAINED}",
            output=f"{a.run}/oracle_reference_2005.csv",
            fallback="C_iso < 4 -> the pass bar becomes >= max(ceil(0.5 x C_30), 5) of the 30 "
                     "hidden families",
            note="this is the ONLY place outside the scorer where a hidden family is used"),
        bars=dict(pass_isolated_ge=W.ACCEPT_REDISCOVERED,
                  split_le=W.ACCEPT_SPLIT_RATE, null_gated_le=W.ACCEPT_NULL_GATED,
                  null_survivors_le_frac=ZS.NULL_SURV_PASS_FRAC,
                  fa_broken_ge_frac=ZS.NULL_SURV_BROKEN_FRAC,
                  fail_below_frac_of_ceiling=ZS.FAIL_FRAC_OF_CEILING,
                  stop_rules="1c' < 2/10 (or < 0.4 x C_iso) with clean nulls -> STOP Arm A north "
                             "AND WA; null survivors >= 50 % of real -> FA broken -> STOP Arm A; "
                             "no third repair"),
        d8_accept=dict(peak_within_s=ZD.ACC_PEAK_S, peak_sites_ge=ZD.ACC_PEAK_SITES,
                       control_range=list(ZD.ACC_CONTROL), chain_ratio_le=ZD.ACC_CHAIN,
                       abandon_chain_gt=ZD.ABANDON_CHAIN, abandon_peak_gt_s=ZD.ABANDON_PEAK_S),
        health_accept=dict(hv_range=[ZH.HV_LO, ZH.HV_HI], tone_ratio=ZH.TONE_RATIO,
                           flat_zero_frac=ZH.FLAT_ZERO_FRAC,
                           dense_max_veto_frac=ZH.DENSE_MAX_VETO_FRAC,
                           expect_flagged=ZH.EXPECT_FLAGGED, expect_clean=ZH.EXPECT_CLEAN),
        quiet_accept=dict(local_mag=ZQ.LOCAL_MAG, local_km=ZQ.LOCAL_KM,
                          global_mag=ZQ.GLOBAL_MAG, pre_pad_s=ZQ.PRE_PAD_S,
                          surface_v_km_s=ZQ.SURF_V, n_hours=W.QUIET_MIN_HOURS),
        isolated_hidden=W.ISOLATED_HIDDEN,
        isolated_hidden_sha256=W.sha256_of(W.ISOLATED_HIDDEN),
        hidden_sha256=W.sha256_of(sorted(int(x) for x in hid)),
        constants=consts,
        sha256=dict(scripts={os.path.basename(p): W.sha256_file(p) for p in scripts},
                    imported_lfe_picker={f: W.sha256_file(f"{PICK}/{f}") for f in IMPORTED
                                         if os.path.exists(f"{PICK}/{f}")}),
        expected_cost="~8-10 GPU-h (3 site-shift null realisations dominate)",
    )
    json.dump(F, open(a.out, "w"), indent=1, default=str)
    print(f"[freeze/1c'] {a.out}")
    print(f"  flags {W.repair_flags()}")
    print(f"  quiet hours {len(QH)} on {nqd} days; veto rows {len(V)} "
          f"({int(V.veto.sum()) if len(V) else 0} vetoed); candidate hours {len(CH)}")
    print(f"  D8 verdict {(F['d8'] or {}).get('verdict')}")


if __name__ == "__main__":
    main()
