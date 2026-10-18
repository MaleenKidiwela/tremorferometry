# lfe_picker — an array-level LFE detector for a margin-wide catalog of LOCATED families

Drafted 2026-09-17 (orchestrator), before Merlin review. Status: PLAN ONLY, nothing trained.

## Why
The template extension (`../lfe_extension`, catalog v1 delivered 2026-09-17) follows Bostock's 130 families on the
permanent network at 12–18 % of his completeness, family-dependent, north only, with inter-SSE detections ≤50 %
real at the loose tier. The limit is measured, not a bug: one LFE is ~2σ at a single permanent station and 4–6 MAD in
the network sum; a template threshold that keeps false alarms tolerable in a continuous scan sits at 8–9 MAD. What is
left on the table is the coherent-array gain (~√N_sites) between the single-station statistic and the array noise
floor, plus everything outside the 130 patches. Lin et al. 2024 reached ~90 % of Bostock with a learned single-station
picker + association on this station class; his weights are unavailable and borrowed pretrained pickers gave ~0 here.

## Goal
Reliable, calibrated LFE detections that yield LOCATED repeating families margin-wide (S. Vancouver Island →
Mendocino), with per-family per-day completeness, judged by (i) a held-out Bostock year in the north, (ii) Ducellier /
Shelly families in the far south, (iii) leave-families-out rediscovery, (iv) the three-gate rule for any new family
(withheld-station stack, interface location with controls, PNSN tremor timing). Nothing is called a family that fails
a gate.

## Design (four ideas, ranked by confidence)
1. **Detect-and-locate over a source grid (learned matched-field beamformer).** Inputs: all stations' 3-C windows
   (60 s, 20 Hz, 1–8 Hz, frozen chain), station positions, and the predicted moveout from a 1-D/3-D velocity model to
   every grid cell (Slab2 interface ±10 km, ~3–5 km cells) as physics tokens. Output: P(source at cell c, time t) plus
   per-station P/S pick heads. Coherent evidence across sites is combined before the decision; location is intrinsic.
2. **Self-supervised pretraining by cross-station prediction.** Predict station j's window from the other stations'
   windows given the moveouts (masked-station reconstruction / contrastive). Only source wavefield is predictable this
   way; noise is not. Runs on all 13 years of continuous data (`../lfe_extension/traces`, cache20hz) with no labels.
3. **Family identity as a learned embedding.** A head that maps each detection to a source-signature vector trained
   with metric learning on Bostock family labels, with instrument-response / gain / azimuth augmentation, so families =
   clusters robust to sensor swaps and station changes. Replaces cc-to-a-fixed-template as the family assignment.
4. **Injection curriculum = completeness model.** Our 116 family stacks (and, later, discovered families' stacks)
   planted in real noise at annealed amplitudes with the true moveout; the recovery curve vs amplitude per
   station-day IS the detection threshold / Mc. Train loud→quiet; never let the injection anchor be mistaken for
   natural-event completeness (E4/E5 lesson).
Supporting: positive-unlabelled loss (catalogs are incomplete; plain negatives are poisoned); random station dropout
during training (the network changes over time and by region); calibrated probabilities (reliability on held-out year;
the Track E realness-score framework); station-subset ablations as the completeness axis.

## Labels and data (mostly on disk)
- Positives north: Bostock 2003–2013 events at our stations with P/S from family travel times (SeisBench stack picks +
  Crosson), 269k events × ~9–15 stations; `../lfe_extension` catalog v1 events with p_real ≥ 0.8 (2014–2026, all
  instrument epochs) once Track E delivers. Positives south: Ducellier 2004–2011 (66 templates / 37 located families)
  and Shelly 2018–2024 (27 families) on roster stations (B045, B935, RVIT, PETY, …).
- Negatives / PU background: PNSN-quiet windows, rolled/shifted controls, catalogued earthquakes and blasts (hard
  negatives, lfe_features R9), the dv/v "families" as known non-sources.
- Continuous data: north 2014–2026 (1.4 TB, on disk) + referee windows 2010–2013; south and Puget/Oregon roster
  stations from 2019 (dense network era), SSE windows first (~0.7 TB), continuous later (~7 TB); disk 72 TB free.
- Velocity: Crosson 1-D (validated in relocation), casc1.6 3-D as a sensitivity run; Slab2 for the grid.

## Stages (each gated; stop cheaply)
- **P0 Data/labels/infra**: label tables, window generator, injection library, station-epoch responses, grid + moveout
  tables per region; the frozen 20 Hz / 1–8 Hz chain reused from `../lfe_extension/scripts/c_*.py`.
- **P1 Honest Lin replication (baseline, ~1 day)**: his architecture (`../lfe_catalog_repair/results/external/lin_cas_lfe_code`)
  trained from scratch on OUR stations with Bostock labels; PyOcto association; scored on a held-out Bostock year
  (2011 or 2013; 2010/2012 are development years) at matched false-alarm rate vs the template catalog's 15–20 %.
- **P2 Array model**: pretrain (idea 2) → fine-tune (ideas 1, 4, PU loss, station dropout) → same held-out test; ablate
  pretraining, PU loss, station dropout; report reliability curves.
- **P3 Family discovery on truth**: hide 30 Bostock families from all labels; do they re-emerge as embedding clusters
  (idea 3), stack, and locate within 5 km (median) with the relocation controls? Rediscovery rate = the bar.
- **P4 Transfer**: train north → test far south (Ducellier/Shelly) and reverse; leave-region-out; injections into
  southern noise to measure sensitivity before trusting detections.
- **P5 Margin-wide discovery (Puget, Olympic, Oregon, 2019→)**: proposals → three gates → families → located →
  per-family completeness → catalog v2 (north: template v1 ∪ picker; elsewhere: picker only).
- **P6 Production + product**: monthly runner like Track C; catalog with p_real, location + uncertainty, family
  embedding id, Mc/λ; Thomas harness per SSE per region.

## Pre-registered bars (report met / not met; never tune to them)
- P1/P2: own-event recall of the held-out Bostock year at the false-alarm rate of the template catalog's T_core; P2
  must beat P1 by a margin that survives the reliability check (calibration within ±0.1 per bin on the held-out year).
- P3: ≥ 2/3 of hidden families rediscovered as clusters with ≥ 50 members and located within 5 km (median).
- P4: far-south recall vs Ducellier/Shelly reported; leave-region-out degradation ≤ 30 % relative.
- P5: every accepted family passes all three gates; the fraction rejected is reported.
- Everywhere: per-family per-day completeness from injections validated against Bostock on the held-out year within
  ±0.05 per magnitude bin (as in `../lfe_extension`).

## Risks / known traps (from E2–E6, Codex audit, lfe_extension)
Component order & rotation of boreholes (StationXML azimuths); many-to-one matching; clocks fitted on the evaluation
period; injection anchors ≠ natural completeness; DeepDenoiser-type models delete LFEs; fixed-MAD comparisons across
network sizes; day-shuffle nulls invariant under day-resolution statistics; single-station family-days; the dv/v
families are NOT sources. Compute: H100 95 GB, 32 CPUs; pretraining on ~13 station-years × 20 Hz × 3-C is days, not
weeks. Governance: Merlin reviews this plan; Opus agents build; Fable orchestrates only; everything documented here
(README/PLAN/NOTES/scripts/results/figures).

## Amendments after Merlin review (2026-09-17; full text results/merlin_review_plan_2026-09-17.md)
- RE-SCOPE: the learned model is a source-agnostic SEED/DISCOVERY engine, not a completeness engine. Two premises fell:
  (1) the √N coherent gain is ALREADY realised by the template network sum; (2) Lin's "90 %" is window accuracy — his
  on-disk catalog recovers only ~10–15 % of Bostock's events chance-free (±10 s), less than our template catalog.
  Honest completeness ceiling on the north permanent network ≈ 0.30–0.40 of Bostock (71/130 families have ≤1 site
  within 30 km). Production detector for every accepted family = the frozen template machinery (Track C/D).
- ARCHITECTURE: per-station learned characteristic function + EXPLICIT physics back-projection on the Slab2 interface
  (2-D 3 km cells, depth fixed on Slab2, Crosson + a5 station corrections). NO learned location head supervised at the
  130 family positions (learned E3 pathology). Learned beamformer = ablation only, after the core passes. SSL
  pretraining optional and gated by a cheap linear-probe test (AUC ≥ 0.85; targets ≥ 15 km from inputs; EQ windows out).
- GATE THE PHASE ON A CPU ORACLE (P0-O): (1) seal 2013 as the final referee, 2011 one-shot per stage, 2010/2012 dev,
  frozen 2012 offsets for all; order 2008/2009 ETS windows; (2) run the frozen template on 2011/2013 = P1a baseline;
  (3) per-station headroom oracle (incoherent S-energy excess vs template cc, single-station ROCs, network-combined at
  matched network-wide quiet-day FA): <1.25× template ⇒ P2 discovery-only; ≥1.5× ⇒ sensitivity claim allowed;
  (4) Shelly-family stack test on roster stations before any southern recall claim.
- LABELS: Bostock 2010–2013 (+2008/09) on-disk ≈ 80k events; S = DB time + frozen offset + a5 stack pick where present
  else Crosson + correction, σ 0.5 s, provenance flagged (far boreholes are PREDICTED labels); P auxiliary only. Certain
  negatives = PNSN-quiet station-days away from Bostock events; hard negatives = EQs/blasts, HF-veto transients, T_ext
  quiet-day detections. ETS-day non-event windows = UNLABELLED (excluded or nnPU), never negatives. Excluded from
  training P1–P3: catalog v1 events, dv/v families, time-shifted controls (nulls only). Boreholes rotated to ZNE for the
  learned model; peak-normalised windows; station dropout incl. whole-site drops; holdouts hashed now: 30 families
  (10 isolated ≥5 km, 10 at 2.5–5 km, 10 <2.5 km) and 3 stations (B926, NLLB, LRIV).
- BARS: all comparisons at matched NETWORK-WIDE quiet-day false events/day (T_core ≈ 12/day, T_ext ≈ 116/day) plus a
  moveout-shuffle null ratio. P1b: single-station S recall ≥ 0.30 at ≤ 100 peaks/station-day at PGC. P2: ≥ 1.25× P1a
  only if the oracle licensed it; else discovery-only. P3: ≥ 2/3 of the 10 ISOLATED hidden families rediscovered with
  one-to-one member purity ≥ 50 % and alias test vs kept families. Completeness calibrated on natural Bostock events
  (injection extends the mapping only, frozen on the north). Gates 0–3 quantified in the review (Gate 0 identity/class
  added; Gate 3 needs local-vs-distant tremor contrast).

## Amendments after Merlin's P1/P2 review (2026-09-19; results/merlin_review_p1p2_2026-09-19.md)
- P2 learned CF + back-projection CLOSED as locator and as event-level discovery engine (no location information; shuffled
  moveout keeps 60 % of recall; per-event error 25 km). The CF survives only as a trigger-ablation arm and as a candidate
  Track E realness feature, each with a kill test.
- P3 = template-free MULTI-STATION REPEATER SEARCH judged by Gates 0–3, opened by a POSITIVE CONTROL (Bostock's own 2012
  event times as candidates): pair-linkability kill test (Step 1a) → full chain (Step 1b) → trigger ablation incl. uniform
  tremor-hour windows (Step 2) → pre-registered discovery nulls (Step 3) → freeze → one 2011 look (Step 5).
- Rediscovery clauses made RELATIVE to the truth template's per-family recall/purity (≥0.5× / ≥0.8×).
- Withheld for discovery and Gate 1: SNB + B926. 2013 stays SEALED as the final referee for discovered products
  (orchestrator decision, user to confirm).
- Corrections to earlier claims: P1b "4.1×" is not chance-corrected (≤2.5×, likely <2×); far-borehole headroom closed; P2's
  per-cell MAD normalisation collapsed on the loudest days (lowhalf must be used if the CF is ever reused).

## P3D — dense-north discovery test, 2019–2026 (user go 2026-09-19)
Premise (from P3): proximity and the number of sites at ~40 km hypocentral distance are the limit, not the method. Ten
3-C permanent stations near the families were never used because they post-date Bostock: SYMB (1.4 km, 2016+), BUCKS
(4 km, 2023+), CLRS (7 km, 2014+), PTRF (8 km, 2017+), RNWY (8 km, 2021+), GOBB (11 km), TURTL (14 km, 2023+), LOPEZ /
OSQM (18–19 km, 2021+), ODUC (21 km, 2023+); LRIV still records (dropped for template failure). Z-only helpers: STW, UWFH.
1. Download continuous 2019-01-01 → 2026-09-15 for the ten (+ LRIV already on disk), SSE windows first, then fill
   (~560 GB); verify actual delivery (metadata "open" ≠ data).
2. Template TRANSFER to the new stations (Bostock step 2 / our A-2 protocol): stack each new station at the catalog's
   high-purity detection times per family; accept per (family, station) by registration vs random-time and
   moveout-shuffled nulls; new network epochs (2019-05.., 2021-04.., 2023-..) with tiers re-derived by the FA-only rule;
   production re-run 2019–2026 on the dense network → recall change for known families (reference: the 13-station
   catalog; no Bostock truth post-2013 — Merlin to prescribe the honest metric).
3. Discovery test on the dense network: pair-linkability per SSE episode 2019–2026 and pooled, candidates TEMPLATE-FREE
   (uniform tremor-hour windows; the learned CF only as an optional trigger arm), truth = the template catalog's
   high-purity events of 30 HIDDEN families (hidden from the transferred templates too), judged by Gates 0–3 and the
   relative-clause rediscovery scorer; report the episode where linkable pairs first cross the bar.
Merlin reviews 2–3 before they run; the download (1) starts immediately.

## Amendments after Merlin's P3D review (2026-09-19; results/merlin_review_p3d_2026-09-19.md)
- PREMISE CORRECTED: the ten stations give the median family ONE new site within 45 km (SYMB), not ~12 sites; CLRS = B926
  site, OSQM = B001 site (merge into sites, withhold by site); GOBB/LOPEZ/TURTL/RNWY are far (65–78 km) and dilute.
- P3D = transfer + recall study + MEASUREMENT of the discovery arithmetic. Discovery (Step 6) runs only if the pre-registered
  prediction from Step 3 crosses (predicted pair-shift δ ≥ 1.5σ, i.e. ρ_SYMB ≥ ~3–4× ρ_PGC). Prior ≤10 %; arm B (seed-and-grow,
  Bostock's construction) prior ~25 %, pre-registered separately.
- CF re-test arm (Step 4) required; 4a timing-residual kill test first; reinstatement needs shuffled-recall ≤0.2 AND location
  error ≤10 km AND ≥3/10 isolated hidden within 5 km.
- Recall metric: prediction (S-SNR vs distance model, hashed first) + measurement (certified marginal detections at withheld
  sites); bar ≥ +25 % certified gain concentrated where geometry says. Injection ratios are NOT a bar.
- Time split: pre-2023 stations calibrate on S01–S13, test on S14–S22; BUCKS/TURTL/ODUC build on S18, test S21/S22.

## P3D RESULTS (agent P3D, 2026-09-19; full text results/P3D_RESULTS.md, ledger rows P3D-0..5)
- Premise correction CONFIRMED: median family gains ONE genuinely new site ≤45 km (SYMB); CLRS/OSQM add no
  geometry and are merged into the B926 / B001 sites.
- Step 1 co-located machinery control PASSES on the powered leg (CLRS vs B926 44/52 = 0.85); OSQM vs B001
  is UNTESTABLE at n=5.
- Step 2: 54 (family, SYMB) templates accepted (gate ≥30). **ρ_SYMB/ρ_PGC = 0.78** — the best-placed new
  site registers single LFEs WORSE than PGC. Every new station is below PGC. The hashed class+distance
  prediction missed SYMB's S-SNR by 3.1 dB, the size of the site-term spread (6.5 dB).
- Step 3: ratio δ_dense/δ_old7 = 1.50 → predicted δ_full 0.57σ → **STEP 6 NOT LICENSED**. PGCq alone is the
  best pair configuration (2.44σ); every added site dilutes it.
- Step 4a: the CF arm **SURVIVES** (SYMB 0.95 s, PGCq 1.05 s, VGZ 1.20 s median |residual|; far sites at the
  2.50 s uniform-chance level). 4b licensed, not run.
- Step 5: at matched quiet-day FA the dense network's certified recall gain is **−5.1 %** against a hashed
  prediction of +27.9 %; Spearman(gain, Δsites) ≈ 0. **Bar NOT MET.**
- PQ amendment (V1/V2) hashed before any PQ stack existed. JRDRB does not deliver; LHLYB has 19/59 of its
  build episode. The S18 V2 readout is UNINTERPRETABLE — 5 null realisations give ratios 1.72–5.36 across
  the bar of 4. **Protocol correction: the δ-ratio rule needs a null-realisation average and a minimum
  populated-quantile count.**
- Step 7 Arm B pre-registered and hashed (results/p3d/arm_b_preregistration.md), not run.
- ADDENDUM 2026-09-19b (PQ dip): all four PQ stations declare HHZ dip +90 (positive DOWN). Only the V2 S18
  pair run had read PQ data, and the flip is a NO-OP there (it cancels in a same-station inner product) —
  verified end-to-end. PQ registration ρ built for the first time: **LHLYB 1.086 (n=8), DAOB 0.826 (n=33)**,
  both above SYMB's 0.783 — the first new stations to match PGC. CF arm still SURVIVES on the V2 site set.
  Licence statistic with the protocol fix (≥3 populated quantiles, null-realisation average, ≥3 admissible
  realisations): **V2 on S18 is UNADJUDICABLE** (1 of 4 realisations admissible); Step 6 remains NOT licensed
  on the calibration run, which satisfies every clause over 3 realisations — ratio 1.52 (LB 1.50) → 0.579σ
  against a 1.5σ bar, a factor 2.6 short. δ_old7 varies ±2.6 % across null realisations on the calibration
  set and ±49 % on S18: the denominator is estimable only at the calibration set's size.

## P4 — Puget / Olympic discovery (user go 2026-09-20)
Region 47.0–48.3 N, tremor-active Slab2 cells 28–45 km; all networks (inventory results/p3d_puget_geometry/): median 4–5
sites ≤45 km per active cell from 2023 (47.0–47.5 N best: median 5, p90 10, azimuthal gap 68°), vs 3 in the north with PQ.
No published families → Gates 0–3 are the only judge; pre-registered nulls are primary evidence.
Lessons applied from P3/P3D: (1) SITE NOISE before geometry — rank every candidate station by 2–8 Hz noise vs PGC and by
tremor-hour vs quiet-hour energy contrast on a few days BEFORE bulk download; (2) seed-and-grow is the construction with
demonstrated linkage power (14.5× chance); pairwise kNN percolates on few-site sums — cluster from fixed seeds via matched
filter + iterated stacks instead; (3) the S-excess clause must be "≥3 of N sites at +3 dB", never unanimity of exactly 3;
(4) withheld sites by whole site; (5) Gate 3 is testable (multi-year span); (6) episodes: PNSN Puget SSEs 2020–2026 pooled;
(7) agents launch shell queues and EXIT; the orchestrator watches; one resume for scoring/write-up.
Stages: P4-0 station shortlist + site-noise inventory (small probes); P4-1 SSE-window download for the shortlist; P4-2 seed-
and-grow + gates with pre-registered nulls; P4-3 anything accepted → stacks → location → templates → backfill on the long-
running Puget stations (PB boreholes B0xx since ~2007). Merlin reviews P4-2 design before it runs.

### P4 — Amendments after Merlin review (2026-09-20; full text `results/merlin_review_p4_2026-09-20.md`)
- **Premise correction (F1/F2):** P3D "ARM B / seed-and-grow" was a seed-restricted *pairwise* graph, never Bostock's construction (template growth → iterated stacks → acceptance); its seeds (top-1 % median-z) were earthquakes/bursts (seed_z up to 8.8e5, 456 seeds in one hour). The 14.5× figure is pair-linkage only. `arm_b_preregistration.md` steps 2–5 have never been executed.
- **New gate for Puget discovery = P4-0b, the true Arm B kill test on the north** (data on disk, kept families, PQ-era S18+S21+S22, per-day availability not intersected): given-seed arm (3 hp events per kept family with ≥50 hp events) and template-free arm (corrected seed rule). Accept (i) ≥2/3 eligible families assembled, (ii) ≥3 rediscovered with ≤1 alias. (i) fail → STOP P4-2. Also reports the seed-z band real hp events occupy.
- **P4-0a station proxy must be calibrated** on the 14 north stations with measured single-event ρ (Spearman ≥0.6, LOO error ≤0.25) before it ranks Puget; otherwise select by geometry + tremor contrast and say so. Sub-regions R1 47.0–47.5 N, R2 48.0–48.3 N, R3 47.5–48.0 N (secondary; median 2 sites in 2025). Proceed rule: ≥4 admitted + ≥2 withheld sites per sub-region.
- **Corrected seed rule:** lowhalf MAD; ≥3 sites z≥3; 4 ≤ z_med ≤ 25; EQ veto ±60 s; HF <3.5; ≤3 seeds/hour, ≥120 s apart, ≤1,500 total.
- **Growth statistic:** 20-s template from S−6; per-site z clipped at 8; hit = sum ≥ T_seed (FA-only rule, ≤0.1/template-day on quiet days) AND ≥2 sites z≥3; independent hits ≥60 s apart; survival ≥20 independent, ≥3 days, ≥2 episodes, ≤20 % within ±10 min of seed, N-moveout ratio ≥3; exactly two iterations; dedupe ≥50 % hit overlap or stack corr > cross-family p95.
- **S-excess at detection sites is a screen, not evidence (F3)** — evidence = Gate 1 at withheld sites, N-shift (≥3 realisations), N-quiet, N-moveout. Gate 1 needs an injection-built Puget control table or it can never FAIL (F6). Engine is hard-wired north (F4): re-point by override tables only.
- **Data split (hashed before download):** BUILD 2022–2024, TEST 2025–2026-05 (sealed), RESERVE 2020–2021 (backfill only). Claim rule requires re-fire on TEST.
- **Prior:** discovery 25–35 % if P4-0b passes both arms; ≤10 % if only given-seed; candidate ~50 %; injection-calibrated bound ~90 %.
- **P4-0b OUTCOME (2026-09-20; full text `results/P4_0B_RESULTS.md`): BOTH ARMS FAIL** — given seeds assemble **0/18** eligible kept families (bar ≥12/18; pass-0 survival 54→5 on the ≥20-independent-hit clause alone, median 6 indep hits over 203 days; near-passes families 44 and 75 clear recall and fail PURITY, 0.298 vs 0.341 and 0.122 vs 0.207) and the template-free arm rediscovers **1** family (44, via F187, 34 matches vs null max 2, 0 aliases) against a bar of ≥3, so the pre-registered consequence **STOP P4-2 — do not run Puget discovery** applies; the pre-registered sensitivity ablation (floors ≥8 and ≥5, `results/p4/p4_0b_abl_prereg.json`, labelled not-the-frozen-result) moves the score only to **1/18** because purity and hit count are anti-correlated (Spearman −0.68) and purity falls monotonically across both iterations, so the floor was not the cause; and the measured seed-z band of true hp events puts **64.6 % below z_med 4**, i.e. Merlin's 4–25 band would discard two-thirds of real LFEs on this network.

**P4-0a outcome (2026-09-20): STOP under the pre-registered sub-region rule, all three sub-regions**, at both the 4-day and Merlin's 6-day clause and under both far-tremor readings (`results/P4_0_RESULTS.md` + addendum). Proxy calibrated on the north (LOO Spearman 0.69, LOO error 0.17; distance alone 0.73). Admitted+withheld sites reach 6+2 (R1) / 8+2 (R3) / 5+2 (R2) under the intent far-bar, but cells with ≥4 admitted sites ≤45 km are 2–8 % vs the 50 % bar (all 118 candidate sites would give 34–60 %). Limit = number of sites that are both near the interface and register tremor (~30 needed per sub-region, 6–10 qualify). Open question for the next Merlin consult (bundled with the P4-0b result to save a Fable call): whether a discovery attempt restricted to the few cells that do have ≥4 admitted sites ≤45 km is legitimate, or whether STOP is final. P4-1 download is NOT launched.

### P4 — CLOSED (2026-09-20, Merlin close-out `results/merlin_review_p4_closeout_2026-09-20.md`)
STOP for Puget discovery is over-determined: (1) Puget has 6–10 quality sites per sub-region at 50–65 km from the tremor field, covering ≤22 % of active cells at ≥3 sites (the 50 % clause itself was unattainable even with all 118 sites — Merlin's own spec error); (2) on the north, single-event seeds starve at pass 0 (13/18 families never reach 20 hits in 203 days) at an untested 8.0 growth floor; (3) the template-free seed rule selects catalogued LFEs 2–7 % of the time regardless of z-band — Puget, with no given seeds, has no entry point. RETRACTED: the "purity collapse is a property of the construction" reading (P4_0B_RESULTS §3/§6) — 97.5 % of hp events are co-detected by neighbouring family templates 0.1–1.6 km away, so hp-only single-label purity is not a noise measure; re-scored with a 3-km cluster truth (post-hoc, for the record). Credited: F187 → family 44 from a non-catalogued seed. No further P4 experiment changes the decision (multi-year north given-seed = optional positive control for lfe_extension only, prior 45 %; F1184 Gate 1 descriptive, prior 10 %). Standing conclusion for the discovery line (P3D + P4): new LFE families cannot be found on this network without a pre-existing catalogue with the tools we have; the only route that has worked is growth from known events/stacks. Parked items: front estimator v2 (must de-duplicate co-detected events / cluster families ≤3 km before fitting); far-south existence test (Shelly's 27 families on B045/B935/RVIT/PETY, ~50–60 GB download, needs user authorisation; Merlin prior ~40 % that ≥10/27 register at ≥2 stations).
