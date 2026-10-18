# PLAN — WA permanent-station TRIAL extension (step 4 of the extension plan)

Step 4 of "dedupe → permanent carriers → permanent templates → trial episode → scan to 2026" for the
Washington families (Royer NW + SW unions, Olympic 2010 Arm A + Arm B). Written 2026-09-29.

## Families
All 329 dup_group representatives of `data/catalogue/families_v1_1.csv` in regions Royer_NW_WA (132),
Royer_SW_WA (79), Olympic_2010 (75), Olympic_2010_armB (43). Detections = each row's `det_source`
(final growth iteration of its own run), matched by `family_id` = `grp`.

## Stations (permanent only)
42 stations: the anchors of `data/extension_stations/{NW_WA,SW_WA,OLY2010}_anchors.csv` (37) plus the
permanent bridge-chain carriers / clear stations of the catalogue not already anchors (UW BLN, SMW,
SQM, OPC; US NLWA). Excluded as temporary: XU XG ZK Z3 7D ZU YW XD ZZ, and also TA/XA/YB/3F/X4/PO/C8.
Sites = stations < 3 km apart (catalogue rule); this reproduces exactly SVI's merge-B clusters
{PGC,B009,B010,B011} and {LRIV,B005,B006,B007}, plus {B013,B943}.

## Windows (`scripts/wa_common.py WINDOWS`)
| label | days | role |
|---|---|---|
| V2010 | 2010-07-14..18 quiet · 2010-08-07..09-04 episode · 2010-09-14..18 quiet | validation, NW + Olympic A/B (their array detections overlap here) |
| V2008 | 2008-05-01..26 episode · 2008-06-16..18 post | validation, SW (SW array detections exist only 2007-01 and 2008-05) + a second NW check |
| T2012 | 2012-08-20..24 quiet · 2012-09-05..10-13 episode · 2012-10-18..22 quiet | trial: first strong WA ETS after the arrays (XU ended 2012-05-31, XG 2011-09) |

Quiet days chosen tremor-free in the regions' PNSN tremor-day tables (2010, 2012). 2008 has no catalogue.

## Recipe (SVI frozen point, `lfe_extension/results/frozen_operating_point.json`, copied)
- 20 Hz, 1–8 Hz Butterworth-4 zero-phase, Z / N|1 / E|2 resolved not rotated, 12-s templates from
  S − 0.5 s, per-channel L2 norm; station cc = mean of channel cc; w = n_channels × SNR weight
  (e5r rule: 10^(S-SNR dB/20), cap 3× median, / median); site merge B; per family-day network MAD;
  record ≥ 8.0, **T_ext 8.25 / T_core 8.75** MAD; ≥ 2 sites per family-day; one-to-one ±2 s own-family
  scoring; site-rolled scrambled control (one shift per site per day in [300, 3000] s); frozen
  cross-family rule (≥ 50 % station sharing, higher raw cc within ±2 s); c_run_month robust energy floor.
- Templates per station EPOCH (bridge node): each station-day is matched only with the template of its
  own node.

## WA equivalents of SVI-specific pieces
- **Site merge B**: SVI named two clusters by hand; here the catalogue's 3-km single-linkage rule, which
  gives the same two clusters + {B013, B943}.
- **HF veto station** (SVI: CN.PGC for every family): no WA station sees all 329 families, so each
  family's HF station is its highest-weight clear template station live that day (all ≥ 40 Hz native),
  E(10–19)/E(2–8) at the native rate (filters designed with fs = native), cut per HF station node =
  99th percentile of build-half own-matched detections (≥ 30), SVI's per-epoch cut rule. Uncalibrated
  node or undefined ratio → never vetoed. In the trial, HF stations with a calibrated node are preferred.
- **Templates**: SVI used Bostock's SeisBench stacks. Here: stack the family's own BUILD-half (2-h
  blocks, 60-s guard; e4 `--half`) array-era detections inside V2008/V2010 on each permanent node
  (≥ 30 events, e5r rule), each event unit-max normalised (Royer restack rule), all-valid windows only.
  Keep only clear-arrival templates (catalogue rule: imp > 5, dur < 3 s, 2–8 Hz designed at fs = 20,
  window [t−12, t+32]). Validation is scored on the TEST half.
- **Trial epochs**: a station on the same node as in validation keeps its template. A new node (sensor /
  band change: e.g. GNW BH50→BH40 2010-12, WISH, SP2, SQM, CPW, NLWA, RER) is templated by the bridge
  rule: pass 1 without it → stack it at pass-1 `cat_core` detections (≥ 30 per family) → clear-arrival
  rule → pass 2 with it. Both passes are reported.

## Outputs
`results/WA_TRIAL_RESULTS.md` (validation table + trial table), `results/score_*/…csv`,
`figures/wa_trial_summary.png`. No gate on tremor windows anywhere; tremor days are a label.

## Added during the build (2026-09-29 17:05Z): per-network-epoch tiers
The first validation days showed ~10 T_core detections per family-day on QUIET days (data ≈ scrambled
control): with 2–3 permanent sites the frozen 8.75 MAD (set on SVI's 9-site network) is not a
0.1-per-template-day false-alarm level. SVI's own recipe re-derives the tiers per network epoch by the
FA-only rule (`lfe_extension/scripts/c_tiers_by_epoch.py`): T_core = smallest T ≥ 8.0 on a 0.25 grid
with quiet-day FA ≤ 0.1 per template-day, T_ext = FA ≤ 1.0. That rule is applied here per window
(V2008 / V2010 / T2012 = one network epoch each, quiet = non-episode observable family-days; V2008's
are post-episode days with no tremor catalogue). Both the frozen pair (`*_core`, `*_ext`) and the
epoch pair (`*_core_ep`, `*_ext_ep`) are reported; no threshold is moved by recall.
Scrambled-control chance matching (`prec_own_ctl`, `recall_ctl`) is reported beside precision/recall:
on ETS days the site-rolled control keeps single-site bursts (B013+B943 is one site carrying two
boreholes), so its COUNT ≈ data on 2-site families, while its matches to the array detections stay ~1 %.

## Step 4b — ADD STATIONS (2026-09-29, user decision after the trial)
The trial showed families separate from the scrambled control only with ≥ 3–4 permanent sites. Step 4b adds
stations for the NW + Olympic families (SW later; 250 of the 329 families, same `fid`s), everything else frozen.
- **Candidates** (`scripts/step4b_setup.py` → `results/step4b/candidates.csv`, every candidate + exclusion reason):
  all stations of `data/extension_stations/{NW_WA,OLY2010}_stations.csv` with an HH/BH/EH ≥ 40-sps epoch in a
  step-4b window (2010–2014), a verified 10-min sample in that year, ≤ 75 km from a family, not in the 42.
  Excluded: discovery-era arrays XU/XG/ZK/Z3/YW/ZU, PH5-only ZZ/YB/3F, GS/PO (no data), stations < 3 km from a
  trial site (merge, no new site: C8 PA array at PGC, TA *A at the anchors, GM at SP2, …). Tier 1 = permanent +
  served later arrays (UW/CN/PB/CC/C8/XD): ALL eligible taken (the planned 95 %-saturation stop ended after a few
  sites). Tier 2 = TA/7D, only for families tier 1 leaves below P(≥3 sites) ≥ 0.5 (≥ 1 expected clear template for
  them, live in a scored window): TA D03D D04D D04E A04D B05D, 7D FN01C FN05C FN07C FN14C FN18C.
  → **52 stations added** (UW 31, CN 6, TA 5, 7D 5, C8 2, PB 1, CC 1, XD 1).
- **Projection** (planning number, not a gate): P(clear | distance, 1c/3c) = the 42-station build clear rate
  (`results/step4b/p_clear_table.csv`), Poisson-binomial P(≥ 3 sites). Families with ≥ 3 clear sites: 25 before →
  28 projected (NW 24 → 27, Olympic 1 → 1, Arm B 0 → 0). The projection says most families CANNOT reach 3 sites
  with the stations that exist 2010–2014 (clear rate ~0.1 beyond 20 km); the run measures it.
- **Windows** (`wa_common.WINDOWS`): V2010 (validation, test half), A2011 = 2011-07-30..08-27 (templating only: NW
  array detections), T2012 (as before), T2013 = 2013-09-04..10-13 + 5 tremor-free days either side, T2014 =
  2014-11-18..12-14 + 5 tremor-free days either side.
- **Bridge templating** (clear-arrival rule only, 2–8 Hz at fs = 20, dur in seconds; no new gate):
  (a) array era: the families' own V2010 BUILD-half + all A2011 detections stacked on every epoch of the added
  stations and on trial-station epochs with no build stack (GNW BH40, SQM EH100, WISH BH40, NLWA 2011, SP2 …) →
  tag `step4b_a`; (b) per episode E: source families = those that separate (median ≥ 3 sites and data/control ≥ 2
  at frozen T_core, tier counts per family-day) in the trial validation, the step-4b V2010 validation, or any
  step-4b trial pass so far; their `cat_core` detections on tremor days of E (accumulated over episodes) stacked
  on epochs without a build stack → tag `step4b_<e>`. First tag listed wins (earlier templates kept).
- **Runs**: V2010 before/after (validation mode, own HF-cut calibration); each episode before (42-station
  templates) / p1 / p2; trial HF cuts = frozen trial cuts + new-node cuts from the step-4b V2010 validation.
  Frozen T_ext 8.25 / T_core 8.75 AND per-epoch tiers; site-rolled scrambled control plus a SECOND independent
  scrambled control (`--n-control 2`) so the separation statistic is applied to a control too.
- **Per family** (`results/step4b/step4b_by_family.csv`, every run, every column): sites (observable, clear
  templates) before/after, recall/precision (V2010), data / control / control2 rates on tremor days (label),
  separation = (data+1)/(control+1) and the same on the control = (control2+1)/(control+1). Ranked, nothing
  dropped; OLY:938 flagged (runaway, B013/B943 pair). Report `results/STEP4B_RESULTS.md`.

## Stack activity — follow each family's ACTIVITY by stacking (2026-09-29, user idea, approved)
Single detections of most WA families cannot be told from the scrambled control on permanent stations (only
families with ≥ 3–4 sites separate). Instead: stack many low-threshold candidate times (windows up to one month)
and ask whether the stack looks like the station template. No new gates; every quantity a column; the scrambled
control gets the same table. Scripts `stack_activity.py` (engine), `stack_activity_score.py`, `stack_activity_report.py`,
`stack_activity_select.py`, queue `stack_activity_queue.sh`; outputs `results/stack_activity/`.
- **Templates**: ALL (family, station-epoch) stacks of ≥ 30 array-era events (clear flag = a column; A_clear uses clear only).
- **Candidates** (fixed before scoring): the family's frozen 1–8 Hz network matched filter, LEAVE-ONE-SITE-OUT (a station's
  candidates come from the other sites only — so its stack is not built from its own noise peaks and random times are a
  valid reference), local maxima declustered ±4 s (4-s block maxima ≥ neighbours), **top-K per day, K = 2 per scored hour**
  (48 per day; 23 on a validation day's test half). Deviation from wa_mf: valid-weighted network mean (sample kept when the
  valid other-site weight ≥ 50 %) instead of requiring every station valid.
- **Stack**: 2–8 Hz (Butterworth-4 zero-phase, |H|² in the frequency domain, designed at fs = 20 Hz = the cache rate), 12-s
  template window, each candidate window unit-max normalised, unweighted. ρ = stack–template correlation (mean over live
  channels). Daily stacks stored → W ∈ {1, 5, 10, 30} d windows = sums of days (calendar blocks from each window's first day).
- **References**: (a) SITE-ROLLED scrambled control ×3 (candidates picked on the rolled network, stack at the rolled data
  positions); (b) 6 × K random times in the same valid scored samples. z = (ρ − random mean)/pooled random SD;
  family **A = Σ z / √n** (Stouffer over stations; A_w network-weighted, A_clear); p_ctl = empirical p vs the 3 controls × all
  windows of the family at that W (each control scored vs the other two).
- **Validation** (known answer): valA = V2008 + V2010 test half (42 stations, build templates); valB = A2011 (independent
  episode, build templates, all halves); valC = V2010 + step-4b stations (build + step4b_a). Active window = **≥ 5** own
  array detections (N = 1, 20 sensitivity); only days the family's array run scanned. AUC + Spearman by W, per family.
  Nothing tuned (the choice of the continuous-series families is a compute budget, not a tuning of the statistic).
- **Apply**: T2012 / T2013 / T2014 (+ quiet days) with build + step4b_a templates (array-era only; the episode-derived
  t2012new / step4b_<e> templates are excluded as circular); tremor days = label only. Continuous 2012-06-01 → 2014-12-31 for
  the top-20 families by W = 1 validation AUC on the top-12 stations (≤ 2 per site) by active-day z gain; valS / valSB
  validate that exact configuration (valSB = build templates only: step4b_a is built from A2011 events).

## Merlin's two deciding tests (2026-09-30, user-requested; queue `scripts/merlin_tests_queue.sh`)
Two flaws Merlin found: wa_mf kept only the 386 clear templates of 11,358 ≥ 30-event stacks (violates rank-don't-discard;
SW's weak carriers were cut), and 'own family' scoring where only ≤ 3-km cluster identity is supported.
- **T1 all templates** (`wa_mf.py --all-templates`, `wa_hf.py --all-templates`, `merlin_t1_score.py`): the network MF with
  EVERY ≥ 30-event build template weighted by its existing SNR weight w (clear flag a column), windows V2008 + V2010 (test
  half) + A2011 (all; truth from step-4b setup), 42 stations; before = the same run clear-only (rerun, identical rules) +
  the existing `mf_val` scored by the same scorer. Frozen T_ext 8.25 / T_core 8.75 + per-epoch tiers (A2011 has no quiet
  days → V2010's), 2 site-rolled scrambled controls. Levels: own (one-to-one ±2 s), ≤ 3-km neighbours (non-transitive,
  lat_use/lon_use), any family. Added rule for BOTH runs: a station-day with > 5 % guarded samples is absent
  (`--drop-dead 0.05`; with every template the all-stations-valid rule let one half-dead station blank half a day for all
  274 families). HF station: clear templates first, then any (by w).
- **T2 withheld-station template recovery** (`merlin_t2.py`): XU/XG stations (V2008 XU, V2010 XU + XG; SH band accepted for
  XG), pairs = (family, window, station with a clear BUILD-half truth template, ≥ 30 events), ≤ 4 per family-window at
  random. Round 0 candidates = stack_activity's permanent-network top-K (verified identical), test-half origins; stack the
  withheld station, cut a template (template rule); rounds 1–2 re-pick with the station in the network (own site, weight
  n_ch × w). Scrambled control = same on the site-rolled network; ceiling = stack at the own TEST-half truth. cc at 2–8 Hz
  (fs = 20) vs the true template (±2 s), rank among all families' truth templates at the station.
- Outputs `results/merlin_tests/` + `results/MERLIN_TESTS_RESULTS.md`. Oregon (XA) not used: no Oregon families here.

## STEP 2 — group detection + member attribution (Merlin idea 2, 2026-09-30; queue `scripts/step2_groups_queue.sh`)
Only ≤ 3-km cluster identity is supported on permanent stations, so detect GROUPS, then ask which member fired as a probability.
- **Groups** (`step2_groups.py --stage groups`): the 329 T1 representatives, complete-linkage on epicentral distance (location or
  approximate position) cut at 3 km (every pair ≤ 3 km, no chaining); no position → own group. 215 groups: 150 single, 49 pairs,
  7×3, 4×4, 3×5, 1×6, 1×23 (21 unlocated Olympic A families at the B013/B943 approximate position + 2 NW).
- **Alignment** (deviation from "own template offset": own S picks are noise-dominated): member origin shift δ by stack
  cross-correlation with the reference member (±10 s; restricted to the co-detection mode ±1 s when that mode is well defined),
  per-station residual τ ±1 s. Group truth = union of member array detections at t + δ, merged ±2 s.
- **Templates** (`--stage tpl`): all members' build-half events (co-detections once) cut at t + τ; wa_templates rule; weights on the
  T1 build scale; every ≥ 30-event template kept. Singleton groups reproduce T1 templates bit-for-bit (checked).
- **MF**: wa_mf `--all-templates --drop-dead 0.05 --n-control 2`, same seed as T1 (identical site shifts, checked) + seed 20261001
  (controls 3–4). Scoring `step2_score.py`: frozen + per-epoch tiers, one-to-one ±2 s vs the union truth, vs T1 detections of the
  members moved to the group frame and merged ±2 s ("T1 summed over members"); separation curve data vs control.
- **Attribution** (`step2_attrib.py`): member templates from their own build stacks at the group window, stations of the 3 nearest
  sites, cc at the aligned position ± 0.5 s (sub-sample lag), P = softmax(β Σ cc) with β fitted on the build half; known answer on
  the test half / A2011 (top-1 vs chance and vs 'most active member', calibration). Stop signal = attribution ≈ uniform.
- **P(real)**: step-1 lfdr refitted on group candidates (`step1_features.py` / `step1_fit.py` with STEP1_* env overrides);
  family-level expected count Σ P(real)·P(member). Report `results/STEP2_GROUPS_RESULTS.md`.

## PNET as independent per-repeat evidence (Merlin step 6, moved ahead of step 3 by the user, 2026-09-30; queue `scripts/pnet_queue.sh`)
- **Model**: `lfe_picker/results/p1/models/pnet_w05.pt` (PhaseNet-style U-Net, 60 s at 20 Hz, 3-C ZNE, trained on Bostock SVI LFEs).
  Inference = `lfe_picker/scripts/p_scan.cf_day`. WA cache = same 20 Hz 1–8 Hz chain; horizontals ROTATED to N/E (PB EH1/EH2 and
  NLWA by StationXML azimuths; HDW single horizontal → Z only; 1-C EHZ stations as (Z, 0, 0)); invalid samples zeroed.
- **A (off-the-shelf)**: CF for the 42 T1 stations on V2008 / V2010 / A2011; for every step-1 candidate (data + 4 scrambled controls)
  the max CF within ±0.5 s of origin + the family's template offset at each station (control read at the site's rolled position);
  site = max over its stations; features `pnet_max`, `pnet_2nd`, `n_sites_pnet03/05` (+ `pnet_mean`, `pnet_rel_max` columns).
  AUC data vs control and truth-matched vs unmatched, by evaluation set × region.
- **B (fine-tune)**: labels = build-half array detections at (family, station node) pairs whose build stack has imp ≥ 3 and dur < 3 s
  (weight 1 if clear, 0.5 otherwise), S at t + template offset, Gaussian σ 0.5 s, loss within ±3 s only; negatives = quiet-day windows
  and the stations' scrambled-control candidate windows (build half), masked ±3 s around every known arrival. 1-C augmentation
  (horizontals zeroed p = 0.3) + polarity flip. Start = SVI weights; 20k steps, Adam 5e-5; checkpoint by internal validation
  (every 10th build-half date). Held out: stations UW.HDW, UW.OSD, PB.B018, PB.B017, US.NLWA, UW.GNW; days V2008/V2010 test half + A2011.
  Evaluation per station: labelled arrival vs quiet-day random time and vs episode-day random time (≥ 3 s from known arrivals); quiet-day
  CF peaks > 0.3 / 0.5 per hour. Pre-set rule "works" = held-out median ΔAUC > 0.02 vs quiet AND > 0 vs episode.
- **C**: `step1_fit.py` unchanged except the feature list (+ `pnet_max`, `pnet_2nd`, `n_sites_pnet03`, `n_sites_pnet05`), once per CF.
- **D**: PNET beam over the family's 3 nearest permanent sites (and nearest 1) on the origin axis, top-K (2 per scored hour, test half)
  → T2 withheld-station stacks, scored by `merlin_t2.py --stage score` on the same pairs/truth/ceiling; run with the off-the-shelf CF
  always (reference) and with the fine-tuned CF only if B "works".
