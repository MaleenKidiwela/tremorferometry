#!/usr/bin/env python
"""ARM B -- the pre-registration.  Writes `results/freeze_armB.json` BEFORE stage 2 runs.

Everything the run may be judged on is fixed here: the weights, the data, the thresholds, the
association and clustering rules, the frozen growth chain it hands off to, the null, and the
ACCEPT BARS.  The engine scripts are SNAPSHOTTED (scripts -> scripts_armB) by the queue before
this runs, and the snapshot's sha256 are recorded, so a concurrent edit of `scripts/w_*.py` by
the Step-1c' work cannot change what Arm B executed.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("W_BED", "c")
import w_common as W                                                        # noqa: E402
import b_armB_pick as BP                                                    # noqa: E402
import b_armB_assoc as BA                                                   # noqa: E402
import b_armB_cluster as BC                                                 # noqa: E402

OUT = f"{W.RES}/freeze_armB.json"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    if os.path.exists(OUT):
        W.log(f"[freeze] {OUT} exists -- NOT overwritten (a freeze is never mutated)")
        return
    cfg = W.bed_cfg()
    mem, of, order, det = W.sites()
    cand = W.candidate_hours()
    qh = W.quiet_hours()
    qh = qh[0] if isinstance(qh, tuple) else qh
    card = None
    try:
        m, dev = BP.load_model("cpu")
        card = BP.model_card(m)
    except Exception as e:                                                  # noqa: BLE001
        card = dict(error=f"{type(e).__name__}: {e}")
    F = dict(
        created_utc=pd.Timestamp.now("UTC").isoformat(),
        study="lfe_array_discovery / ARM B -- north control with LFEDetect(jms) on the 2005 bed",
        role=("Merlin Step-1 audit Q3 and 'Experiments I could not run': the only Arm-B north "
              "control that is not a training-set test.  Written BEFORE any pick is computed."),
        authorisation=("results/merlin_review_step1_2026-09-20.md Q3: 'Arm B's north control on "
                       "the 2005 bed is legitimate only with jms (Nankai + Guerrero + San "
                       "Andreas, no Cascadia), and WA must then run the same weights'"),
        weights=dict(model="LFEDetect", version="jms", version_str="1",
                     training_regions=["Guerrero, Mexico (Frank et al. 2014)",
                                       "San Andreas fault (Shelly 2017)",
                                       "Nankai, Japan (JMA catalog)"],
                     cascadia_in_training=False,
                     excluded_weights=dict(
                         names=["cascadia", "cms", "jcms", "jcs"],
                         reason="trained on Bostock et al. 2015 Cascadia stacks -> a north "
                                "control with them is a TRAINING-SET test"),
                     model_card=card),
        data=dict(root="traces/bed2005_polaris", rate_hz=100, components="3-C raw miniSEED",
                  why_not_cache="the model card is 20 Hz + 1-8 Hz zerophase + peak norm and "
                                "seisbench.annotate applies that chain itself at the native "
                                "rate; the frozen 20 Hz 1-8 Hz cache would band-pass twice",
                  boreholes="PB (EHZ, EH1, EH2) unrotated; seisbench "
                            "flexible_horizontal_components=True accepts Z12 as ZNE",
                  window=f"{cfg['day0']}..{cfg['day1']}", quiet_days=cfg["quiet_days"],
                  n_bed_days=len(W.bed_days())),
        bed=dict(sites=list(order), site_members={k: list(v) for k, v in mem.items()},
                 detection_stations=list(det), withheld=list(W.WITHHELD),
                 withheld_excluded_from_association=True,
                 n_candidate_hours=int(len(cand)), n_quiet_hours=int(len(qh))),
        picks=dict(default_threshold=BP.THR_DEFAULT, low_threshold=BP.THR_LOW,
                   sweep=list(BP.SWEEP), overlap=600, batch_size=512,
                   nesting="picks_from_annotations triggers on prob > t, off at t/2, so the set "
                           "at threshold t is exactly {thr_first <= t <= prob}"),
        association=dict(engine="pyocto.OctoAssociator.from_area",
                         velocity_model="Crosson (1976) Puget layers, Vp/Vs 1.78, verbatim from "
                                        "lfe_catalog_repair/scripts/e2_associate.py",
                         zlim_km=list(BA.ZLIM),
                         zlim_rule="Slab2 interface band prior (PLAN sec.10 ACTIVE_CELL_DEP)",
                         box="detection-station bbox + 0.3 deg (NEVER the family locations)",
                         n_picks=3, n_s_picks=3, n_p_picks=0, min_stations=3,
                         min_stations_rule="e2_associate's documented Vancouver-Island rule "
                                           "(>= 3 stations inside VI, >= 4 elsewhere); a 2-day "
                                           "pilot measured 0 events at 4 stations / thr 0.30 vs "
                                           "15 at 3 stations, and 144 at 3 stations / thr 0.15 "
                                           "(n_stations 113x3, 28x4, 3x5).  The 4-station count "
                                           "is run as a DIAGNOSTIC and reported, never as a gate",
                         pick_match_tolerance=2.5, edt_pick_std=1.5, min_interevent_time=3.0,
                         levels=[BP.THR_DEFAULT, BP.THR_LOW]),
        clustering=dict(algorithm="DBSCAN on epicentres (local km projection)",
                        eps_km=BC.EPS_KM, min_samples=BC.MIN_EVENTS,
                        min_cluster_events=BC.MIN_EVENTS, time_radius="whole window",
                        seed_budget=W.SEED_BUDGET,
                        carry_rule="the DEFAULT-threshold (%.2f) level is the one carried into "
                                   "growth; the lower level (%.2f) is associated and clustered as "
                                   "a diagnostic and is carried ONLY if the default level yields "
                                   "zero clusters" % (BP.THR_DEFAULT, BP.THR_LOW)),
        growth=dict(chain="w_grow handoff -> quiet T_seed -> grow -> scramble -> survive -> "
                          "w_dedupe -> w_locate_gates (Gates 0-2)",
                    passes=1,
                    template="frozen pass-0: %.0f s from onset - %.0f s (TPL0_LEN_S/TPL0_PRE_S). "
                             "The 12-s template (TPL1_LEN_S) is produced only by w_iterate, "
                             "which is the SECOND growth pass and is not run here."
                             % (W.TPL0_LEN_S, W.TPL0_PRE_S),
                    T_seed="FA-only per template on the 2005 quiet days, <= %.2f false hits per "
                           "template-day, engineering floor %.1f" % (W.FA_PER_TPL_DAY,
                                                                     W.T_ENG_FLOOR),
                    survival=dict(min_indep=W.SURV_MIN_INDEP, min_days=W.SURV_MIN_DAYS,
                                  max_near_seed=W.SURV_MAX_NEAR_SEED,
                                  scramble_ratio=W.SCRAMBLE_RATIO),
                    dedupe=dict(jaccard=W.DEDUPE_JACCARD, cc_min_sites=W.DEDUPE_CC_MIN_SITES),
                    location="w_locate_gates: Slab2 2-D grid + 4 controls, Gates 0-2 with "
                             "withheld " + ",".join(W.WITHHELD)),
        null=dict(realisations=1,
                  rule="ONE whole-chain site-day-shift realisation of stages 2-5: each SITE's "
                       "day-series is circularly shifted by a site-specific whole-day offset "
                       "(members together); the picks are moved by the same map "
                       "(w_common.site_shift(1)), so this is bit-equivalent to re-running "
                       "LFEDetect on the shifted day-series",
                  seed=W.SITESHIFT_SEED),
        accept=dict(n_isolated=3, of=10, centroid_km=5.0, recall_abs=0.30,
                    purity_abs_cluster3km=0.50, null_clusters_le=1,
                    statement="PASS = >= 3 of the 10 isolated hidden families rediscovered with "
                              "centroid <= 5 km, one-to-one recall >= 0.30 at +- 2 s and purity "
                              ">= 0.50 under the 3-km cluster truth, AND <= 1 cluster in the "
                              "site-shift null realisation"),
        stop_rule=("Merlin: Arm B(jms) north control < 2/10 -> STOP the study and write it up as "
                   "'array-era discovery not reproducible from the public archive at this SNR "
                   "with these methods'"),
        scoring=dict(scorer="w_score_step1.py (frozen), isolated stratum "
                            f"{W.ISOLATED_HIDDEN}",
                     absolute="bars above, computed in b_armB_score.py",
                     relative="results/step1c_prime/oracle_reference_2005.csv if it exists by "
                              "scoring time; otherwise marked 'reference pending' (Merlin F4a "
                              "rejects the transferred 2012 reference)"),
        deviations=[
            "Templates are the frozen 20-s pass-0 form, not 12 s: the 12-s template is defined by "
            "w_iterate, which is a second growth pass, and the brief authorises ONE growth pass.",
            "Picks are computed on the raw 100 Hz traces with seisbench's own preprocessing "
            "(model card), not on the frozen 20 Hz cache.",
            "Quiet-day picks are computed as well (not used for association) so the tremor-vs-"
            "quiet pick rate is on the record.",
            "The engine is a SNAPSHOT (scripts_armB/) of scripts/ taken at queue start, because "
            "scripts/w_*.py are being edited concurrently for Step 1c'.",
        ],
        sha256=dict(
            b_scripts={os.path.basename(p): sha(p)
                       for p in sorted(glob.glob(f"{HERE}/b_arm*.py"))},
            engine={os.path.basename(p): sha(p)
                    for p in sorted(glob.glob(f"{HERE}/w_*.py"))},
            snapshot_dir=HERE),
    )
    json.dump(F, open(OUT, "w"), indent=1, default=str)
    W.log(f"[freeze] wrote {OUT}")


if __name__ == "__main__":
    main()
