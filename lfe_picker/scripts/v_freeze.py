#!/usr/bin/env python
"""P4-0b FREEZE -- written and hashed BEFORE the queue is launched.

Everything the run may use is fixed here: the sha256 of every v_ script, the seed rules, the
threshold rule, the survival / iteration / dedupe / screen / acceptance clauses, the site list and
the per-day availability table, the eligible-family list with hp counts, and the quiet-day list.
Hidden families (results/holdouts.json) are never seeded, scored or matched; the withheld stations
never enter a sum; 2011 and 2013 stay sealed (s_common.cache_index refuses them).
"""
from __future__ import annotations

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
FRZ = f"{V.RES}/p4/p4_0b_freeze.json"


def main():
    os.makedirs(OUT, exist_ok=True)
    hid, iso, ho = S.holdouts()
    keep, D, _ = V.eligible_families()
    lat, lon, dep, dep_s2 = V.region_centre()
    A = V.availability()
    ED = V.episode_days()
    Aep = A[[t in ED for t in A.date]].copy()
    Aep["episode"] = [ED[t] for t in Aep.date]
    Aep.to_csv(f"{OUT}/availability_episode_days.csv", index=False)
    Q = V.quiet_days()
    Q.to_csv(f"{OUT}/quiet_days.csv", index=False)
    cst = V.cell_tt(lat, lon, dep)
    FT = V.family_tt()

    frz = dict(
        created_utc=pd.Timestamp.utcnow().isoformat(),
        spec="results/merlin_review_p4_2026-09-20.md, sections F1/F2/F3, 'P4-0b' and P4-2 steps 2-6",
        preregistration="results/p3d/arm_b_preregistration.md steps 2-5 (never previously executed)",
        what="the TRUE Arm B construction (matched-filter growth from a seed -> iterated stacks -> "
             "acceptance) run on the PQ-era north, where the 130 Bostock families minus the 30 "
             "hidden ones are truth; its result gates whether Puget discovery (P4-2) happens",
        scripts={os.path.basename(p): V.sha256_file(p) for p in V.v_scripts()},
        governance=dict(
            hidden_families=sorted(hid),
            hidden_rule="never seeded, never scored, never matched; removed before any seed table",
            withheld_stations=V.WITHHELD_STAS,
            withheld_rule="never in any network sum in any pass (holdouts.json held_out_stations "
                          "B926/NLLB/LRIV plus the four P3D withheld sites' members)",
            sealed_years=[2011, 2013],
            puget_untouched="no Puget data is read or written; results/p4/ probe outputs of the u_ "
                            "scripts are not modified and no FDSN request is made",
            holdouts_sha256=ho["sha256_of_hidden_family_list"]),
        setting=dict(
            episodes=V.EPISODES,
            n_episode_days=int(len(Aep)), n_usable_episode_days=int(Aep.usable.sum()),
            sites=V.SITE_ORDER, site_members=V.SITE_MEMBERS, stations=V.DET_STAS,
            min_sites=V.MIN_SITES,
            per_day_availability="results/p4/p4_0b/availability_episode_days.csv -- per-day, NOT "
                                 "intersected across stations (intersecting starved NET-B)",
            region_centre=dict(lat=lat, lon=lon, dep_used_km=dep, dep_slab2_km=dep_s2,
                               deviation="lat/lon are p_common.GRID_LAT/GRID_LON; the DEPTH is the "
                                         "median of the 130 Bostock family depths because Slab2 at "
                                         "the box centre returns 21.9 km, ~13 km above the family "
                                         "field, which would mis-cut every template by ~3 s"),
            region_station_constants={s: float(cst[s]) for s in V.DET_STAS},
            quiet_days=[str(t.date()) for t in Q.date],
            n_quiet_days=int(len(Q)),
            quiet_rule=f"cached, >= {V.MIN_SITES} sites, no PNSN tremor within "
                       f"{V.QUIET_TREMOR_KM:.0f} km of the region centre; drawn from the PQ era"),
        eligible_families=dict(
            rule=f"kept (non-hidden) Bostock families with >= {V.GIVEN_MIN_HP} hp events "
                 "(realness_v1_1 p_real >= 0.8) inside the pooled episodes on usable days",
            bkeys=keep, n=len(keep),
            hp_counts={int(b): int(n) for b, n in D.groupby("bkey").size().items()},
            acceptance_i=f">= 2/3 of the {len(keep)} eligible families assembled "
                         f"(>= {int(np.ceil(2*len(keep)/3))})"),
        seed_rule_given=dict(
            n_per_family=V.GIVEN_SEEDS_PER_FAM,
            rule="3 hp events per eligible family, distinct days, distinct episodes preferred, "
                 "ranked by the frozen catalogue's own net_mad (Bostock starts from a LOUD window)",
            moveout="the family's own Bostock cell (Crosson 1-D + a5 station correction)"),
        seed_rule_free=dict(
            windows=f"{V.SEED_WIN_S} s at {V.SEED_HOP_S} s hop on PNSN tremor hours with "
                    f">= {V.SEED_MIN_EPI_HOUR} epicentres within {V.SEED_EPI_KM} km of the centre, "
                    "episode days only",
            energy="per-site 3-C energy z at the region moveout, MAD-normalised per station-day on "
                   "the LOWER HALF of that day's window-energy distribution",
            band=dict(min_sites_z3=V.SEED_MIN_SITES, site_z=V.SEED_SITE_Z,
                      zmed_lo=V.SEED_ZMED_LO, zmed_hi=V.SEED_ZMED_HI),
            eq_veto_s=V.SEED_EQ_VETO_S, eq_catalog=V.EQCAT,
            hf_rule=f"E(10-19)/E(2-8) < {V.SEED_HF_MAX} at the two nearest stations, NATIVE rate; "
                    "an unmeasurable ratio is never a veto",
            cap=dict(max_per_hour=V.SEED_MAX_PER_HOUR, min_sep_s=V.SEED_MIN_SEP_S,
                     budget=V.SEED_BUDGET, order="z_med descending within the band"),
            band_caveat="the frozen 20 Hz cache is 1-8 Hz, so the '2-8 Hz' seed energy is measured "
                        "in the cache's 1-8 Hz band -- FLAGGED, not corrected"),
        growth=dict(
            template=f"{V.TPL0_LEN_S} s from predicted S - {V.TPL0_PRE_S} s at the seed's moveout, "
                     "each channel normalised to unit L2",
            cc="frozen runner recipe (o_run_month.py): full-day normalised correlation per channel "
               "by rFFT with a FLOAT64 prefix-sum denominator; station series = mean over 3 channels",
            per_station_z="lowhalf MAD: m = median(cc), MAD = 1.4826 x (m - q25(cc)), z = (cc-m)/MAD",
            per_site="mean of the member stations' z, then CLIPPED at "
                     f"{V.SITE_Z_CLIP} (anti-dominance)",
            network=f"EQUAL-WEIGHT SUM over the sites observable at that sample, then converted to "
                    f"net_mad by the frozen runner's per-(template, day) median/lowhalf-MAD "
                    f"normalisation -- which is what makes ONE fixed per-template threshold "
                    f"comparable across days whose site count differs",
            hit=f"net_mad >= T_seed AND >= 2 sites with per-site z >= {V.SITE_Z_COINC} at the hit "
                f"sample AND >= {V.MIN_SITES} sites observable",
            decluster_s=V.DECLUSTER_S, independent_s=V.INDEP_S,
            T_rule=f"per template: smallest T on a {V.T_GRID} grid at or above the {V.T_FLOOR} "
                   f"floor with <= {V.FA_PER_TPL_DAY} false hits per template-day on the quiet "
                   f"days (the frozen FA-only rule, c_tiers_by_epoch.py docstring); T is capped at "
                   f"{V.T_MAX} and a template that never reaches the rate is flagged `saturated`",
            note="the per-site clip at 8 bounds the network sum at 8 x n_sites, so net_mad "
                 "saturates near 20 for a 6-site day; this is a property of Merlin's clip"),
        survival=dict(min_independent=V.SURV_MIN_INDEP, min_days=V.SURV_MIN_DAYS,
                      min_episodes=V.SURV_MIN_EPISODES,
                      max_frac_within_10min_of_seed=V.SURV_MAX_NEAR_SEED,
                      n_moveout=f"per template: the SAME seed re-cut with a WRONG cell's timing "
                                f"(predicted S differing >= {V.NMOVEOUT_MIN_DT} s at >= half the "
                                f"stations), grown identically with its own T; require true "
                                f"independent hits >= {V.NMOVEOUT_RATIO}x the scrambled count",
                      applied="ONCE, at pass 0; iterations 1 and 2 re-stack and re-grow"),
        iteration=dict(n=2, align_stack=f"{V.ALIGN_ITERS} iterations, +-{V.ALIGN_LAG/V.FS} s "
                                        "(s_common.align_stack, the frozen recipe)",
                       repick="S re-picked on the aligned stack at the HIGHEST-S-SNR site; that "
                              "single global shift is applied to every site",
                       recut=f"{V.TPL1_LEN_S} s from S - {V.TPL1_PRE_S} s",
                       rule="exactly two; no third (convergence is a forking path)"),
        dedupe=dict(overlap=V.DEDUPE_OVERLAP, tol_s=V.DEDUPE_TOL_S,
                    stack_rule="per-station cross-TEMPLATE |cc| p95 recomputed from the grown bank "
                               "(as g_crossfam.py did), exceeded at >= half the sites where both "
                               "register",
                    common_axis="origin times converted to predicted S at PGC so templates cut at "
                                "different cells are comparable",
                    keep="the larger (most independent hits)"),
        screen=dict(rule=f"alignment-free mean stack of the independent hits at each detection "
                         f"site; >= +{V.SCREEN_DB} dB S excess at >= 3 sites OR >= 50 % of "
                         f"registering sites",
                    status="CHEAP SCREEN, REPORTED, NEVER USED AS ACCEPTANCE (Merlin F3)"),
        scoring=dict(
            arm_i=dict(rule="relative clauses of q_rediscover.py: recall >= "
                            f"{V.REL_RECALL} x recall_ref and purity >= {V.REL_PURITY} x "
                            "purity_ref of the family's own reference template on the same days, "
                            "withheld stations excluded from detection",
                       reference="the family's own TRUTH-STACK 12-s template (align_stack over all "
                                 "its hp events at each site, cut from S-2 s), grown through the "
                                 "identical threshold/growth machinery",
                       reference_deviation="Merlin's wording is 'the family's own transferred "
                                           "template'.  The P3D transfer bank on disk covers only "
                                           "26-184 (family, station) pairs and NONE at "
                                           "PGC/VGZ/ALBHB, so it cannot supply a 7-site reference "
                                           "for the 18 eligible families.  The truth-stack "
                                           "template uses ALL of a family's hp events and is at "
                                           "least as strong as a transferred 2012 template, so the "
                                           "relative bar is CONSERVATIVE.  Absolute recall/purity "
                                           "are reported alongside.",
                       no_reference="a family whose reference recall is 0 is reported as "
                                    "NO-REFERENCE and counted neither way (q_rediscover's rule)"),
            arm_ii=dict(rule=f"one-to-one +-{V.MATCH_TOL_S} s matches of a grown template's "
                             f"independent hits to a kept family's hp events; rediscovered if "
                             f"matches >= {V.FREE_MIN_FRAC:.0%} of that family's hp events AND "
                             f">= {V.FREE_MIN_MATCH} matches",
                        offset="a template-free seed's origin axis is the REGION-centre moveout, so "
                               "a true family's hits arrive displaced by ONE constant "
                               "b = t_seed - tau_0 (identical for every event; derivation: the "
                               "template matches when t + cst_s - pre = tau + tt_s - delta_s, and "
                               "delta_s is fixed by the seed, giving t = tau + (t_seed - tau_0) "
                               "for every station).  b is estimated as a single nuisance parameter "
                               "per (template, family) over +-30 s and the search is CHARGED to a "
                               "null in which the template's hits are redrawn uniformly inside "
                               "their own days (20 realisations); the strict b = 0 numbers are "
                               "reported alongside.",
                        alias="a grown template meeting the clause for a second kept family; a "
                              "rediscovery counts only with <= 1 alias"),
            also_reported=["the seed-z band TRUE hp events occupy (median-over-sites z of their "
                           "12-s region-moveout windows on >= 3 sites) -- this fixes the Puget "
                           "seed rule",
                           "the burst-convergence signature: the fraction of families whose "
                           "given-seed growth puts > 50 % of its hits within +-10 min of the seed",
                           "the free-seed time histogram: max seeds/hour and the top-10-hours "
                           "fraction (> 20 % means the cap is not working)"]),
        acceptance=dict(
            i=f">= 2/3 of the {len(keep)} eligible kept families assembled",
            ii=">= 3 kept families rediscovered template-free with <= 1 alias each",
            consequences=dict(
                both="P4-2 prior ~30 %",
                i_only="seeds are the bottleneck; P4-2 runs as an injection-calibrated bound, "
                       "prior <= 10 %",
                i_fails="STOP -- do not run P4-2 until the construction is repaired and re-passes"),
            report_immediately="given-seed growth returning > 50 % of hits within +-10 min of the "
                               "seed for the majority of families (burst convergence) means the "
                               "decluster/independence guards are not working"),
        outputs=dict(dir=OUT, done_marker=f"{V.RES}/p4/P4_0B_DONE",
                     failed_marker=f"{V.RES}/p4/P4_0B_FAILED",
                     score=f"{OUT}/P4_0B_SCORE.json", log="logs/p4_0b_queue.log"),
    )
    frz["self_sha256_excluded"] = "this file's own hash is not included; the v_ script hashes are"
    json.dump(frz, open(FRZ, "w"), indent=1, default=str)
    print(f"[freeze] -> {FRZ}")
    print(f"  {len(frz['scripts'])} v_ scripts hashed; {len(keep)} eligible families; "
          f"{int(Aep.usable.sum())}/{len(Aep)} usable episode days; {len(Q)} quiet days")
    for k, v in frz["scripts"].items():
        print(f"  {k:24s} {v[:16]}")


if __name__ == "__main__":
    main()
