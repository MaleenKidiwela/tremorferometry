#!/usr/bin/env python
"""STEP 0 -- write `results/freeze_wa_v1.json`: the Washington pre-registration.

Nothing in Step 5 may be scored unless the scripts on disk hash to what is written here.  The
freeze carries the sub-region and window, the admitted and withheld sets, the measured moveout
spread, the engine constants (verified BYTE-IDENTICAL to `freeze_step1c.json` apart from three
documented region constants), the threshold and survival rules, the location protocol, the Gate
settings including the measured HF cut, the whole-chain null, the injection design, the claim
ladder verbatim, the Chestler descriptive-prediction protocol, and the sha256 of every `x_` and
`w_` script and of every input table.

    python scripts/x_freeze.py
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

# The ONLY engine constants that differ from Step 1c, each with its written reason.
ALLOWED_DIFFS = {
    "GRID_LAT0": "Merlin Q4: the Washington grid centre is 48.00 (north bed 48.50)",
    "GRID_LON0": "Merlin Q4: the Washington grid centre is -123.20 (north bed -123.75)",
    "ACTIVE_CELL_DEP": "Merlin F2/Step 0: WA active cells are Slab2 28-45 km (north bed 25-50)",
    "LOC_MIN_SITES": ("Merlin Q4: 'located' needs >= 5 picked sites <= 60 km in WA; w_common's "
                      "own comment records the north bed used 4 because it has only 6 sites"),
}
BED_KEYS = {"BED", "BED_MODE", "STEP1_YEAR", "STEP1_DAY0", "STEP1_DAY1", "WITHHELD",
            "STUDY_CACHE", "STEP1B_SITE_KEYS", "STEP1C_SITE_KEYS", "STEP1C_QUIET_DAYS",
            "STEP1C_BED", "STEP1C_ENERGY_HOURS", "STEP1B_STATIONS_CSV", "ISOLATED_HIDDEN",
            "WITHHELD_A", "WITHHELD_B", "SCORE_WINDOWS", "WORK", "ACTIVE_CELL_MIN_EPI",
            "ACTIVE_CELL_KM", "ACTIVE_CELL_DEG"}


def current_constants(W):
    out = {}
    for k in dir(W):
        if not k.isupper():
            continue
        v = getattr(W, k)
        if isinstance(v, (int, float, str, bool)) or (
                isinstance(v, (list, tuple)) and all(isinstance(x, (int, float, str, bool))
                                                     for x in v)):
            out[k] = list(v) if isinstance(v, tuple) else v
    return out


def check_constants(W):
    """Every engine constant must equal Step 1c's, except the documented region/bed ones."""
    F = json.load(open(f"{X.RES}/freeze_step1c.json"))["constants"]
    cur = current_constants(W)
    diffs, undocumented = {}, {}
    for k, v in F.items():
        if k in BED_KEYS:
            continue
        if k not in cur:
            continue
        if cur[k] != v:
            diffs[k] = dict(step1c=v, wa=cur[k])
            if k not in ALLOWED_DIFFS:
                undocumented[k] = diffs[k]
    return cur, diffs, undocumented


def main():
    W = X.activate()
    C = X.cfg()
    cur, diffs, undoc = check_constants(W)
    if undoc:
        raise SystemExit(f"UNDOCUMENTED engine-constant changes vs freeze_step1c.json: {undoc}")
    X.log(f"[freeze] engine constants vs Step 1c: {len(diffs)} documented differences "
          f"{sorted(diffs)}; every other constant byte-identical")

    hf = json.load(open(f"{X.OUT}/hf_cut.json"))
    corr = pd.read_csv(X.CORR_CSV)
    cells = pd.read_csv(f"{X.OUT}/active_cells.csv")
    tremor_days = X.tremor_days()

    F = dict(
        created_utc=pd.Timestamp.now("UTC").isoformat(),
        study="lfe_array_discovery / STEP 5 -- N-Olympic 2010 discovery, Washington",
        role=("Merlin Step 0: the Washington pre-registration.  Written BEFORE the Step-5 queue "
              "is launched and before any Washington byte is scored."),
        supersedes=None,
        sub_region=dict(name=X.SUBREGION, lat=[X.LAT0, X.LAT1], lon=[X.LON0, X.LON1],
                        centre=C["centre"],
                        active_cell_rule=C["active_cell_rule"], n_active_cells=int(len(cells)),
                        grid=C["grid"]),
        window=dict(day0=X.DAY0, day1=X.DAY1, year=X.YEAR,
                    n_bed_days=int(len(W.availability())),
                    tremor_days=tremor_days, n_tremor_days=len(tremor_days),
                    tremor_day_rule=("PNSN-corrected: a day carrying at least one CANDIDATE HOUR, "
                                     "i.e. >= 5 PNSN epicentres within 50 km of an active cell"),
                    quiet_days=C["quiet_days"], quiet_rule=C["quiet_rule"]),
        sites=dict(admitted_rule=C["admitted_rule"],
                   n_admitted_stations=C["n_admitted_stations"],
                   site_merge_km=C["site_merge_km"], site_rule=C["site_rule"],
                   site_order=C["site_order"], site_members=C["site_members"],
                   site_table=C["site_table"],
                   sites_dropped_by_spread_rule=C.get("sites_dropped_by_spread_rule", []),
                   n_sites=C["n_sites"], n_detection_stations=C["n_detection_stations"],
                   excluded=dict(CL="XG.CL array: 0/16 admitted by the probe, median contrast "
                                    "1.9 dB",
                                 PFB="CN.PFB: local-far -10.6 dB; the S-WA 'far' days are not far "
                                     "from Vancouver Island"),
                   coverage=C["coverage"]),
        moveout=dict(ac_win_s=10.0,
                     spread_s_at_centre_all_admitted=C["moveout_spread_s_centre_all_admitted"],
                     spread_s_at_centre=C["moveout_spread_s_centre"],
                     spread_s_median_over_cells=C["moveout_spread_s_median_over_cells"],
                     rule=("the 10-s window is kept if the admitted-site spread at the sub-region "
                           "centre is <= 10 s; it was 21.12 s, so the SAME written rule used for "
                           "Step 1b/1c was applied -- sites added by increasing median epicentral "
                           "distance while the median spread over the active cells stays <= 10 s")),
        withheld=dict(stations=C["withheld_stations"], sites=C["withheld_sites"],
                      rule=C["withheld_rule"],
                      cell_visibility_frac=C["withheld_cell_visibility_frac"],
                      scope=("never in any autocorrelation sum, growth sum or stack, for the "
                             "whole run; Gate 2 DOES use them for its own independent location, "
                             "as the north design did, and that is recorded here")),
        engine=dict(source="scripts/w_*.py, unchanged; the WA bed is bound in process by "
                           "scripts/x_common.py",
                    constants=cur,
                    identical_to="results/freeze_step1c.json",
                    documented_differences={k: dict(v, reason=ALLOWED_DIFFS.get(k, ""))
                                            for k, v in diffs.items()},
                    min_sites=W.MIN_SITES,
                    borehole_channel_convention=("PB boreholes are cached (EHZ, EH1, EH2) and are "
                                                 "NOT rotated; every statistic in this chain is "
                                                 "rotation-invariant about the vertical, so the "
                                                 "reader requires only (vertical, horizontal, "
                                                 "horizontal)")),
        threshold=dict(rule=("FA-only, no 8.0 floor: T is the smallest value on a 0.25 grid at or "
                             "above the ENGINEERING floor 5.0 with <= 1 false pair per quiet hour "
                             "AND zero quiet-hour components of size >= 5; where it binds is "
                             "reported"),
                       engineering_floor=W.T_ENG_FLOOR, grid=W.T_GRID,
                       quiet_days=len(C["quiet_days"]),
                       quiet_hours_target=W.QUIET_MIN_HOURS,
                       quiet_hours_per_day=W.QUIET_HOURS_PER_DAY,
                       growth_fa_per_template_day=W.FA_PER_TPL_DAY,
                       stop=("if the quiet hours contain >= 5-member components at every T <= 20, "
                             "write STEP5_FAILED with 'site-coherent noise' and stop")),
        seed_cluster=dict(node_merge_s=W.NODE_MERGE_S, min_members=W.SEED_MIN_MEMBERS,
                          min_members_pm1h=W.SEED_MIN_MEMBERS_PM1H, bin_s=W.SEED_BIN_S,
                          min_bins=W.SEED_MIN_BINS, budget=W.SEED_BUDGET),
        survival=dict(min_indep=W.SURV_MIN_INDEP, indep_s=W.INDEP_S, min_days=W.SURV_MIN_DAYS,
                      max_frac_near_seed=W.SURV_MAX_NEAR_SEED, near_seed_s=W.NEAR_SEED_S,
                      scramble_ratio=W.SCRAMBLE_RATIO,
                      scramble_null=("each station's template cut from a DIFFERENT random tremor "
                                     "window and grown identically; real independent hits must be "
                                     ">= 3 x scramble")),
        iterations=dict(n=W.N_ITERATIONS, tpl0_len_s=W.TPL0_LEN_S, tpl1_len_s=W.TPL1_LEN_S,
                        tpl1_pre_s=W.TPL1_PRE_S,
                        third_pass_rule=("a third iteration is run only if the survivors' median "
                                         "S-SNR (w_grow handoff diag `s_snr_db`) rose by >= 0.5 dB "
                                         "from iteration 1 to iteration 2; the rule is applied by "
                                         "the queue and its inputs are written to "
                                         "results/step5/iter_snr.json")),
        dedupe=dict(jaccard=W.DEDUPE_JACCARD, cc_min_sites=W.DEDUPE_CC_MIN_SITES,
                    tol_s=W.DEDUPE_TOL_S, lag_s=W.DEDUPE_LAG_S),
        location=dict(
            model="Crosson 1-D (Puget), g_common.tt_table / o_vel1d",
            corrections=dict(
                source="results/intraslab_eq_2006_2012.csv + traces_eq/ (+-3 min event windows)",
                recipe=("a5_eq_corrections: native rate, 2-15 Hz Butterworth-4 zero-phase, P on Z "
                        "within +-3 s and S on sqrt(H1^2+H2^2) within +-4 s of the Crosson "
                        "prediction, envelope onset + AIC, accepted at peak/noise-rms >= 4"),
                usable_rule="n >= 5 accepted S picks AND residual MAD <= 0.5 s, else 0.0 and FLAGGED",
                n_stations=int(len(corr)), n_usable=int(corr.usable.sum()),
                file="results/x_station_corrections.csv",
                observed=("the MAD clause binds, not the pick count: 35 of 86 stations have "
                          ">= 5 S picks but MAD 0.5-0.9 s and therefore take 0.0")),
            grid=dict(lat0=X.GRID_LAT0, lon0=X.GRID_LON0, halfwidth_km=X.GRID_HALFWIDTH,
                      mode="2-D Slab2 interface + free-depth 3-D sanity run"),
            picks=dict(search_s=W.LOC_PICK_SEARCH_S, min_snr=W.LOC_MIN_PICK_SNR,
                       array_rule="members picked individually, site median"),
            located_clause=dict(min_picked_sites_60km=W.LOC_MIN_SITES,
                                max_gap_deg=W.LOC_MAX_GAP_DEG,
                                shuffled_pick_ratio_lt=W.LOC_NULL_RATIO,
                                n_pick_null=W.LOC_N_PICK_NULL,
                                synthetic_floor_le_km=W.LOC_SYNTH_BAR_KM, n_synth=W.LOC_N_SYNTH,
                                loso_le_km=W.LOC_LOSO_BAR_KM,
                                free_depth_within_km=W.LOC_DEPTH_BAR_KM,
                                not_at_bound=True),
            accuracy=("ESTIMATED from the Step-5b injection at a known cell, never certified from "
                      "the controls (Merlin F7)"),
            known_error=("g_common.epi_km and the Gate-1 moveout null bind cos(48.6) as a default "
                         "argument, a ~1 % east-west distance error at 48.0 N; it is applied "
                         "identically to the data and to every null and is left unchanged so the "
                         "engine stays byte-identical")),
        gates=dict(
            run=[0, 1, 2],
            gate3="UNTESTABLE inside one ETS window (Merlin F6); evaluated on the extension record",
            offset_search_s=W.GATES_OFFSET_SEARCH_S,
            gate0=dict(
                hf_station=f"PB.{C['hf_station']}",
                hf_veto=C["hf_veto"],
                hf_rule=hf["rule"],
                hf_quiet_distribution=dict(n=hf["wa_n"], median=hf["wa_median"],
                                           p25=hf["wa_p25"], p75=hf["wa_p75"],
                                           p95=hf["wa_p95"], p99=hf["wa_p99"],
                                           days=hf["wa_days"]),
                hf_north_reference=dict(station=hf["north_station"], n=hf["north_n"],
                                        median=hf["north_median"], cut=hf["north_cut"]),
                hf_cut_over_quiet_median=hf["cut_over_quiet_median"],
                hf_reported_intraslab_eq_median=hf["eq_median"],
                hf_reported_note=hf["reported_eq_note"],
                distinctness=("the ADOPTED registering-site rule (g_gates.gate0): among the sites "
                              "where the proposal's own stack REGISTERS an arrival (S-window "
                              "excess >= 2 dB), its waveform must be below the cross-template cc "
                              "p95 at most of them (frac_reg_above_chance < 0.5), OR the proposal "
                              "must be >= 5 km from every accepted family.  Merlin's literal "
                              "'>= 2 sites below chance' rule is reported but NOT used (it "
                              "false-passes the alias control).  The accepted list here is the WA "
                              "candidate list of this run itself, leave-one-out by drop_bkey, "
                              "because no WA family has been accepted yet; the cross-template p95 "
                              "is this bank's own (w_dedupe)")),
            gate1=dict(withheld_stations=C["withheld_stations"],
                       min_members=20, n_time_null=200, n_moveout_null=100,
                       min_excess_db=3.0, moveout_min_dt_s=5.0,
                       pass_rule=("random-time null: observed > p99 AND excess >= +3 dB; "
                                  "moveout-shuffled null: excess > the p95 of the shuffled excess"),
                       non_pass=("UNTESTABLE unless a Step-5b injected repeater at matched N and "
                                 "comparable distance PASSED at the same withheld station "
                                 "(Merlin F4); both the raw and the re-labelled verdict are "
                                 "written")),
            gate2=dict(min_sites_60km=W.LOC_MIN_SITES, max_gap_deg=W.LOC_MAX_GAP_DEG,
                       min_pick_snr=W.LOC_MIN_PICK_SNR, n_pick_null=W.LOC_N_PICK_NULL,
                       n_synth=W.LOC_N_SYNTH, synth_bar_km=W.LOC_SYNTH_BAR_KM,
                       loso_bar_km=W.LOC_LOSO_BAR_KM, depth_bar_km=W.LOC_DEPTH_BAR_KM,
                       null_ratio=W.LOC_NULL_RATIO)),
        whole_chain_null=dict(
            n_realisations=W.N_SITESHIFT_REALISATIONS, seed=W.SITESHIFT_SEED,
            min_shift_days=W.SITESHIFT_MIN_DAYS,
            design=("the ENTIRE pipeline re-run with each SITE's day series circularly shifted by "
                    "a site-specific whole-day offset, members shifting together, so site-coherent "
                    "noise is preserved exactly and only true inter-site coherence is destroyed"),
            accept="real gated count >= 5 x the null mean AND null mean <= 1 per realisation",
            stop="null >= 3 gated families per realisation -> the FA rule is broken, stop"),
        injection=dict(
            role=("MECHANICS + LOCATION-ERROR DISTRIBUTION ONLY.  NEVER a sensitivity ceiling: E3 "
                  "step 3 measured injected stacks recovering 4-11x off real at matched s/n"),
            n_cells=4, cell_choice="random active cells of the sub-region, seed 20260920",
            n_events_per_family=40, amplitudes=[0.5, 1.0, 2.0],
            amplitude_definition=("target rms = amp x the donor event's OWN signal-to-noise at its "
                                  "own northern station x the recipient station-day's local noise "
                                  "rms in the 12 s before the insertion point"),
            donors=("northern hp events, lfe_extension/results/trackE/realness_v1_1.parquet with "
                    "hp == True, drawn from distinct northern families; the 30 hidden families of "
                    "lfe_picker/results/holdouts.json are NEVER donors"),
            donor_class_rule=("borehole -> borehole, surface -> surface; within the class the "
                              "donor station is the one whose own epicentral distance is closest "
                              "to the recipient's distance from the planted cell"),
            insertion=("additively into ADMITTED and WITHHELD station-days at the planted cell's "
                       "moveout (Crosson + the frozen station correction), on the 20 Hz 1-8 Hz "
                       "cache, into a SEPARATE cache root whose uninjected days are symlinks"),
            scoring=("the SAME chain is run on the injected root; a planted repeater is recovered "
                     "when a located proposal matches >= 30 % of its planted times at +-2 s; the "
                     "location error and the per-withheld-station Gate-1 verdicts at the planted "
                     "N become the WA Gate-1 control table")),
        claim_ladder=dict(
            candidate=("N >= 20, Gate 0 PASS (incl. the HF veto at a WA 100-Hz borehole; distinct "
                       "from every already-accepted WA family by the adopted registering-site "
                       "rule), Gate 1 PASS at >= 1 withheld site"),
            gated_family=("Gate 1 PASS at >= 2 withheld sites (both nulls: 200 random-time p99 "
                          ">= +3 dB and 100 moveout-shuffled p95); located; whole-chain site-shift "
                          "null x3: real >= 5x null mean, null mean <= 1"),
            confirmed_recurrent=("re-fires in the sealed 2011-08 window or in >= 2 later ETS with "
                                 "Gate 3 PASS -- NOT evaluable in Step 5"),
            source="Merlin review 2026-09-20, answer Q5/(b), verbatim"),
        chestler=dict(
            file=X.CHESTLER,
            protocol=("after the gated list is frozen, report the distribution of nearest-Chestler "
                      "distances for the gated families against 1,000 draws of the same number of "
                      "random active cells; NEVER used to select, merge or locate"),
            n_random_draws=1000),
        sealed_2011=dict(window=["2011-07-20", "2011-09-10"],
                         state="sealed and unread; results/SEALED_2011_DO_NOT_READ",
                         manifest="results/sealed_2011_manifest.sha256"),
        accept=dict(step5=(">= 5 gated families with null mean <= 1 and real >= 5 x null; injected "
                           "repeater recovered, gated and located within 5 km at the middle "
                           "amplitude"),
                    stop=("0-2 gated with a clean null -> report, do not open other windows "
                          "without a diagnosis; null >= 50 % of the real count -> FA broken, stop; "
                          "injected repeater not recovered at any amplitude -> engine/moveout bug, "
                          "stop")),
        compute_expectation=dict(autocorrelation_gpu_h=5.0,
                                 measured_s_per_tremor_hour=48.4,
                                 n_candidate_hours=int(len(W.candidate_hours())),
                                 growth_and_nulls_gpu_days="1-2"),
    )
    sha = dict(x_scripts={os.path.basename(p): X.sha256_file(p) for p in X.x_scripts()},
               w_scripts={os.path.basename(p): X.sha256_file(p) for p in X.w_scripts()},
               inputs={os.path.basename(p): (X.sha256_file(p) if os.path.exists(p) else None)
                       for p in X.input_tables()
                       + [X.BED_JSON, X.CORR_CSV, f"{X.OUT}/hf_cut.json",
                          f"{X.OUT}/active_cells.csv", f"{X.OUT}/admitted_stations.csv",
                          X.CHESTLER]})
    F["sha256"] = sha
    json.dump(F, open(X.FREEZE, "w"), indent=1, default=str)
    X.log(f"[freeze] wrote {X.FREEZE}: {len(sha['x_scripts'])} x_ scripts, "
          f"{len(sha['w_scripts'])} w_ scripts, {len(sha['inputs'])} input tables hashed")


if __name__ == "__main__":
    main()
