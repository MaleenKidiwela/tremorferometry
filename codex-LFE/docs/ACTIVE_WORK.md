# Active objective

Continue until a finalized full-margin LFE catalog. Intermediate pilots and failures are checkpoints, not completion. Never label unverified candidates a final scientific catalog.

Scope: all usable archive history through the present; initial cutoff 2026-09-14T18:03:52Z. Exact roster: 183 stations in results/target_network_from_map/stations.csv. Metadata is not waveform availability.

## Storage instruction

User explicitly requested monthly download/process/delete processing. Do not download the whole archive. See results/full_history_work_plan_v1/streaming_policy.json and monthly_processing_queue.csv. Process newest to oldest, starting at the present cutoff. Verify actual station availability and revalidate templates across instrument/network changes. One active calendar month, with network-coherent daily/subdaily units and a provisional 100 GB temporary waveform ceiling. Preserve boundary padding, snippets, provenance, coverage, controls, family state and verified outputs. Cleanup is limited to newly staged production-owned raw files after verification; existing shared caches are not disposable. Policy and queue exist; production runner integration remains unfinished.

## Immediate priority: picker validation

User explicitly reaffirmed picker testing before production. Published labeled-stack benchmark failed on a nonbucket trace identifier and is not restarted pending reuse audit; prototype lives in codex_lfe.stack_picker_benchmark (results/picker_labeled_stacks_jms_v1.log). Uses jms, whose documented training excludes Cascadia, nominal observed P/S labels (not predicted travel times), original split labels, fixed threshold/tolerance grid, actual prediction coverage, CPU/GPU parity and horizontal permutation sensitivity. Does not establish continuous false-pick rate or full-margin validity. Existing 24-window review packet remains unlabeled; never score pending windows as negatives. Reuse existing experiments described in docs/CLAUDE_TEST_AUDIT.md.

## Current state

- Northern expanded 2012 training/calibration running: results/north_full_family_bank_2012_v1.log; ten instruments, seven physical sites, PB.B926 reserved for independent checking. Reuses existing cached waveforms.
- Southern repaired January bank has two four-site families. February 1–14 calibration failed for equal-channel, vector, expanded bank and three-site variants. February 15–29 remains sealed. Failure requires diagnosis, not a lower unvalidated threshold.
- Southern continuous acquisition totals 70 station-days, approximately 902 MB. Full archive downloading has not started.
- Historical metadata and 26,017 station-epoch monthly work blocks exist. All remain awaiting validation; no claimed full-history processing.
- Northern experimental checkpoint: 61 candidates on 17 days, not final events. Southern 61,441 published reference records are external benchmarks, not our detections.
- Latest tests: 98 passed. GPU correlation numerical validation passed; production GPU integration remains unfinished.

Remaining: validated regional detection and discovery, independent event checks, timing/location solutions, monthly processing integration, all-history scans, sensitivity/coverage assessment and versioned release. Do not contact authors without explicit messaging authorization.

Latest user instruction: check prior tests to avoid duplication. See results/prior_picker_test_audit_v2/reuse_ledger.json. E6 was completed and needs to be included alongside E2/E4 and our corrected interleaved comparison. Completed: 56 E2 S-reference rows reproduced and rescored one-to-one without inference or downloads. See results/prior_picker_test_audit_v2/README.md. At score .3 / ±3 s, cascadia full-stack reference recovery changes 6.11% to 4.19%. Phase truth and coverage remain separate gaps.

Workspace boundary reaffirmed by user: all new code, downloaded data, outputs, logs and notes must remain inside codex-LFE. Other projects are read-only sources. Do not run legacy entry points that write back into their original results directories; copy/adapt read-only calculations into codex-LFE instead.

Resumed after user objected to stopping: remain active through intermediate checkpoints. Northern expanded calibration finished (47 usable source groups; 82 candidates at 9 MAD; 71 same-family matches among 2,222 eligible references; worst control ratio 1.22%). No full-history production. Specific missing cross-region-weight P/S labeled-stack diagnostic is running as picker_labeled_stacks_jms_v2, with per-trace checkpoints and nonbucket HDF5 support; prior failure preserved. All modifications remain within codex-LFE.

## Latest completed picker work

- Published-stack jms timing: results/picker_labeled_stacks_jms_v2 and picker_labeled_stacks_tail_v1. Tail check covers every labeled P/S arrival and reuses previously covered results; production inference code unchanged. Score .3 / ±.5 s: P 569/1067 (53.3%), S 663/1453 (45.6%). Family bootstrap and figure in tail directory. These are stack-label diagnostics, not continuous precision. Systematic early residuals remain; no correction deployed.
- Compared installed VariableLengthPhaseNet forward/skip logic with upstream v0.7.0: inspected methods unchanged. Version diff saved in prior_picker_test_audit_v2. No software-version timing correction inferred.
- Southern transfer: results/southern_picker_transfer_v1, jms on five stations × Feb 1,7,14 2024, all using preexisting raw files. Response metadata only was newly fetched into southern_picker_response_2024_v1. All 15 station-days complete with explicit coverage and probability arrays. BK.RVIT Feb14 has only 1411 valid inference windows vs typical2877; preserve this exclusion.
- Reference packet coincidence uses January-derived stack peaks (NOT actual S labels): at score .3 / ±2 s, zero matches over the small two-family cohort across all five station tables. Do not infer an operational false-pick rate or physical inability from this underpowered proxy test.
- Building results/southern_picker_review_v1: 30 pending windows in six paired blocks across five southern surface/borehole stations. Predictions reused from transfer run, no new inference. Initial review is blinded to model outputs and sampling stratum; labels must not be fabricated.
- Asked asynchronously whether additional reviewed phase labels/results exist outside the three audited projects; no reply yet. Continue independent useful work while awaiting any paths.

Completed southern review v2: 29 pending windows, six paired blocks, all29 have full P/S prediction coverage; one BK.RVIT packet excluded for missing waveform coverage. See results/southern_picker_review_v2/review.html and verification.json. First all-five-station sampling attempt failed before outputs and is preserved under v1.

Component check on review windows: BK.RVIT median max/min RMS ratio6.53, median maximum cross-component CC0.951; other sites median ratios1.07–2.02. Descriptive only, no station automatically excluded. Original alternate-event southern template QC halves share17–20% of eight-second waveform interval unions for the two four-site families. This is potential dependence, not proof of false templates. A new training-only UTC-date-parity consistency check is running in results/southern_template_date_check_v1.log; same saved January packets, no downloads, no production template replacement.

## Continuous southern benchmark resumed

The roster-matched Ducellier family 080424.05.060 now has a working one-hour continuous development scan, using published 35-second templates, response-corrected velocity, verified component transforms and 1.5–9 Hz filtering. Code: codex_lfe/ducellier_continuous_check.py; reference comparison: codex_lfe/ducellier_continuous_report.py. Results: results/ducellier_continuous_development_v1 and results/ducellier_continuous_next_day_v1. Seven stations requested, five usable; NC.KRMB returned no data and NC.KHBB lacked valid full-triplet coverage. These are not assertions of archive-wide absence.

2010-02-17 18:00 UTC: shared observed/control exposure 2,415 seconds. At 8 MAD, seven candidates all match published times within 0.5 seconds, among eight eligible published records; eight shifts yield one candidate total. At 10 MAD, three candidates all match, and controls yield zero. This is selected development data, not unbiased precision/recall. Control boundary exclusion leaves eight of sixteen hourly references eligible. 2010-02-18 18:00 UTC: zero observed candidates at all four thresholds; published catalog also has no records in that hour. Absence from the reference does not prove a negative truth label. Existing matcher and evaluation tests: 14 passed.

Picker validation has NOT passed for production. The new result establishes a useful matched-filter benchmark on the requested array, while neural picker phase validation and full-margin discovery remain unresolved. Do not substitute published records for our detections. Continue bounded validation and monthly runner implementation; no full-history scan has started.

The April 25 04:00 UTC hour completed: at 8 MAD, ten candidates, nine matches among ten eligible published records (0.5-second tolerance). A full November 22, 2010 day is now running under the unchanged grid; all seven requests succeeded, totaling 144,243,712 bytes. Each hour test used about 5 MB. These remain development benchmarks because the three dates were selected from the published catalog.

Monthly storage lifecycle implemented in codex_lfe/monthly_storage.py, not yet integrated into production acquisition/detection. Enforces one active month, descending completed-month order, staging budget, retained artifact categories, atomic output commitment and checksum verification before deleting only ledger-owned staged raw files. Five tests pass, including changed-output/raw preservation, symlink protection, restart, budget and unfinished-month checks. Scientific validation remains separate; no existing raw archives were deleted. Single-writer integration is required.

## Full-day edge validity and source grouping correction

The one-template November 22 scan completed before conservative response-taper masking: 8 MAD gives nine candidates and five matches among seven published records, with one total candidate in eight shifted controls. These preliminary counts must not be used as final validation: ObsPy's default response-removal taper spans 5% of each input trace, so complete template windows touching that taper need exclusion. `codex_lfe/response_validity.py` now defines the metadata-only mask, with a passing edge-support test. The multi-template rescorer will apply it to saved correlations without redownloading or repeating response removal/correlation computation.

The comparison bank has seven template identifiers but FIVE published Plourde families: three templates belong to family 026. Source mapping and checksum are saved in `comparison_bank_source_families.csv` and `comparison_bank_grouping_note.json`. Do not sum template-level counts into unique LFE counts. Cross-template physical-clock alignment and duplicate resolution remain required.

An independently specified next-day panel is running on November 23, 2010. Its fixed design was written before reading that day's reference events: 8 MAD, seven template groups, PB.B039 withheld from detection, minimum three other supporting sites, and eight shifts. The held-out station provides descriptive independent waveform scores, not an automatic phase-truth label. See `results/ducellier_next_day_evaluation_design_v1.json`. New waveform acquisition was limited to this second day in the same November month; no full archive scan has started.

## Completed southern next-day panel

Corrected results are in `results/ducellier_multifamily_november_edges_v2` and `results/ducellier_multifamily_november23_holdout_edges_v2`; summary/figure in `results/ducellier_validation_summary_v1`. Complete-template response-taper exclusion leaves approximately 22.50 hours of shared exposure per template. At the preselected 8 MAD on November 23, template 080424.05.060 yields 11 candidate rows and all four eligible published references; 080421.14.048 yields seven rows and all three references; 080422.12.039 yields four rows and two of three references. Other groups have few or no references. These per-template counts are not unique event totals.

Training-template correlation establishes provisional offsets within published family 026: add +5.35 seconds to template 080421.15.050 and -3.95 seconds to 080422.13.003 to align with 080421.14.048. Network lag cycle closes; some individual stations disagree. Offset signs have a passing synthetic test. Grouping is diagnostic only, not physical-origin determination. After within-family grouping, the next-day panel has 23 observed candidate rows across four published families, with cross-family duplicates unresolved.

Withheld PB.B039 confirmation is mixed: median CC is negative for family 026 and modestly positive for 027/029. This does not certify individual candidates or prove unsupported rows false. Three nearby pairs match weak template self-correlation sidelobes within 0.15 seconds; they are flagged for review, not deleted. All 23 candidate rows now have waveform/template overlay figures and snippets in `results/ducellier_candidate_waveform_review_v1/index.html`; all review labels remain pending. An inspected strong 029 example shows clear BK.YBH waveform similarity but weak borehole signal. This is a template audit, not blinded phase truth.

## Northern evaluation resumed without duplicate training

A provenance bug was found: the pilot sealed `input_manifest.json` in the training artifact, then extended that same file during calibration. Reconstructing its training-date subset reproduces the original sealed hash EXACTLY (40 station-days). Recovered artifact: `results/north_full_family_training_recovered_v1`; neither training nor calibration was rerun. The pilot now writes a separate immutable `training_input_manifest.json`; legacy artifacts remain supported and changed inputs still fail verification.

The frozen 9-MAD evaluation is running on the bank's already declared September 18–20, 2012 dates, reusing cached raw recordings and the 47-group training bank. Results/log: `results/north_full_family_frozen_evaluation_v1`. These dates are evaluation for this bank, not claimed untouched by all other projects. PB.B926 remains reserved for independent-site checking. No production archive scan has started.

## Northern evaluation and independent-site results

The frozen northern evaluation completed: 30 candidates, 28 any-source temporal matches, 18 same-source matches among 681 eligible references, and ten ambiguous family assignments. Mean shifted-control ratio is 1.25%, but the maximum is 6.67% (two control candidates versus 30 observed), exceeding the frozen 5% criterion. Do not claim that gate passed or retune the threshold on these dates.

PB.B926 date-separated training produced coherent templates for 15 of the 47 source groups. Only two selected evaluation candidates belong to those supported groups; their signed CCs are 0.125 and 0.445, with positive excess over 100 shifted checks. The 61 eligible published-reference checks have median CC 0.088. This supplies limited independent support, not validation of the remaining candidates. Results: `results/north_full_family_B926_validation_v1`. The immutable training-provenance fix has three passing artifact/end-to-end tests, including preservation under runtime-manifest updates.

## First current central-roster panel

`results/central_recent_picker_panel_v1` contains one hour from September 13, 2026, selected around the completed availability probe. Six roster stations returned approximately 6 MB in total and have full core P/S prediction coverage. UO.UCHR has no current triplet in the saved metadata (last epoch ends March 26, 2026); PB.B027 returned an empty waveform response. Neither is interpreted as absence of LFEs. Initial acquisition handling failed on the missing UCHR probe, and that attempt's log is preserved; completed downloads were reused on restart.

The panel uses jms with response-corrected velocity, verified ZNE, a fixed five-second taper and the existing 1–8 Hz inference path. No phase-accuracy labels exist for this hour. At score 0.3 there are 17 individual P/S peaks across six stations. Controlled homogeneous-velocity association yields no four-site seeds at score 0.1 or 0.3, and no three-site seeds at 0.3. At 0.1/three sites there are 11 observed seeds versus 7–13 per shifted control, so this configuration does not distinguish the observed hour from chance alignment. Do not promote these seeds to templates or catalog events. The velocity model is explicitly diagnostic, not a calibrated Oregon location solution. Independent tremor context is being checked before choosing any further activity-specific panel.


## Central activity-matched panel completed

PNSN public tremor metadata selected September 2, 2026 09:00–10:00 UTC and eight nearest successfully probed roster stations before waveform inspection. Matching September 13 quiet data used the same panel. At threshold 0.1, three-site associations were 69 observed versus 59–73 shifted controls in the active hour; four-site associations were 16 versus 13–20. At 0.3, one three-site candidate versus 0–2 controls. No discovery validation passed. Comparison and waveform review: `results/central_activity_comparison_v1`. The diagnostic association location is not a validated event hypocenter.

A fixed 15-second window-grid stability diagnostic now uses the same saved preprocessed waveforms, with no new downloads or threshold optimization: `results/central_grid_stability_v1`. This checks a visible stitching discontinuity and cannot provide missing phase truth.


## Grid-stability and daily-plan follow-up

Fixed 15-second inference-grid perturbation retained 18/41 active-panel S picks at score >=0.3 within 0.5 seconds; installed SeisBench overlap averaging retained 4/10. These are stability measurements, not accuracy scores. With averaging, active associations were 34 versus 25–33 controls (three sites) and 7 versus 3–7 (four sites); threshold 0.3 yielded no observed associations. Neither configuration is approved.

The September daily plan now has 14 descending daily units, 2,394 station-epoch request records and 168 missing-metadata station intervals. No production downloads were made. Boundary/gap tests pass (3); current output is `results/september_daily_work_plan_v2`. The earlier v1 plan is superseded by the general partial-day gap handling in v2.

Published Cascadia stack station geography supports three same-code roster matches within 100 m: CN.PGC, CN.VGZ, CN.SNB, despite the source PO network field. Historical instrument/response and horizontal orientation equivalence remain unverified. See `results/published_stack_roster_geography_v1`.


## Stronger 2026 central panel and specific reference-data dependency

The cached PNSN count metadata selected March 2, 2026 (79 reported observations), then its highest-energy hour, 04:00–05:00 UTC (22 observations), before waveform inspection. Eight roster instruments were acquired, but B031/B032 are within 1 km: the corrected analyses select B031 by station ID and count seven physical sites. The initial cross-station run was stopped and marked invalid; no partial v1 scores are accepted.

Corrected March neural association gives 38 three-site candidates versus 31–41 shifted controls at score 0.1, and 7 four-site candidates versus 1–7 controls. At score 0.3 it gives zero observed candidates. Cross-station coherence at 0.5 gives 10 observed versus 11–15 controls. Its seven higher-CC seeds cluster at 04:06:56–04:08:42; the strongest waveform is consistent with the P arrival and subsequent signal from USGS event us7000s1cq, M6.0 near Japan. IASP91 predicts P arrivals at the Oregon stations around 04:07:07–04:07:10. This is a contamination flag, not a blanket exclusion of possible LFEs during that interval. No seed was promoted to a template bank.

Review artifacts: `results/central_validation_summary_v1/index.html`, `results/central_strong_waveform_review_v1/index.html`, `results/central_strong_conventional_event_audit_v1`. Seven new tests pass across the daily-unit planner (3), physical-site selector (1), cross-station normalization/lag sign (1), and MAT-reference importer (2).

A more specific Washington reference lead is now documented: `LFEsAll.mat`, with `LFEs[].name` and `LFEs[].peakTimesMo`, is used in the public ArianeDucellier/longrange scripts. Its README explicitly says underlying data are available on request. The file was not found in the workspace; the user was asked asynchronously for a local path if already available. A derived coordinate-only pickle contains 69 unidentified rows and is not treated as the paper's final 43-family event catalog. The importer is implemented and tested, preserving unresolved time convention and family selection. No authors have been contacted; concrete drafts are in `docs/REFERENCE_DATA_REQUEST_DRAFT.md`.

## Northern source identity audit, 2026-09-14

Reused cached Bostock table, SeisBench metadata, and prior published-stack clock mappings; no waveform downloads or picker inference. `results/north_source_identity_audit_v1` preserves strict rounding and coordinate-interval evidence separately. Strict decimal rounding uniquely matches 128/130 reference coordinate groups. Allowing the printed-coordinate half-unit interval matches 129/130, including all 47 usable training families. The additional match is published family 47 at longitude -123.655, a decimal rounding boundary. The remaining group (48.51, -123.92, 32.94 km; 3,446 source rows) is unresolved. These are inferred identity links, not authoritative family assignments or independent phase truth; frozen evaluation and production status remain unchanged. Earlier clock-transfer work is recorded and was not repeated.

## Continued work after the intermediate-stop correction

The expanded northern calibration loss audit (`results/north_full_bank_calibration_loss_v1`) reused the existing diagnostic implementation on this bank for the first time: 2,134/2,222 eligible references were below the fixed 9-MAD network threshold, 84 passed joint network/site criteria, and four passed network threshold but failed site support. These are per-reference windows, not counts of unique final events. The frozen evaluation was not rescored or retuned.

The training family audit (`results/north_family_distinguishability_v1`) evaluated 1,081 pairs. Two conspicuous families share 126 matched training anchor times out of 140/139 records, with mean fixed-anchor template CC 0.741 across three physical sites. No pair met the declared simultaneous >=80% bidirectional temporal overlap and >=0.8 CC on three-site flag. No family was merged.

Storage integration now supports explicit empty units and months without epochs. `catalog_products.py` verifies retained experimental products, full declared-roster coverage, support-window coverage, family references, controls and clocks before `commit_verified_unit` enables owned staging cleanup. Details and limitations are in `docs/RETAINED_UNIT_PRODUCTS.md`. The tested lifecycle preserves original waveforms and retained products. This does not validate a scientific production method.

## Northern phase transfer and current-era tests

Built 47 separate-date three-station training-stack panels from cached September 2012 data (`north_published_phase_training_stacks_v1`). Published waveform comparison under unique rounding-interval identities gave 165 P/S transfers; 35 initially had sufficient whole-date support to measure stability, and 15 passed all criteria, supporting five model-dependent family clocks. The follow-up `north_date_stability_completion_v1` completes 25 previously unmeasurable checks using count-balanced whole dates. It reproduces each original full stack to <=1e-10, reuses all original full-stack scores and previously measurable checks, and gives 32 passing transfers and nine model-dependent family clocks. No family has three passing phase-transfer sites. Individual-event origin times remain unvalidated. Two PGC peak anchors clearly align with P rather than S.

`north_bank_epoch_compatibility_v1` shows that the 2012 bank cannot run unchanged in 2026: PGC changed BH->HH, SNB selected metadata ends March 2021, zero complete original triplet sets remain, and only 19/47 families retain >=3 station-identity sites under an unvalidated station-subset change.

Selected an August 2026 northern panel using only public tremor metadata: latest 2026 northern day with >=20 reported windows is August 18 (43); highest-energy hour 10-11 UTC has 11 windows. The comparison is August 20, same hour, the next later date with no reported windows in the bounds. Eight distinct physical sites: CN.PGC, CN.VGZ, UW.RNWY, UW.LOPEZ, UW.OLGA, UW.GUEM, PB.B001, PB.B005. Processed comparison first, then activity; total raw acquisition 18,012,160 bytes. jms three-site association at .1 gives 57 activity vs 42-50 controls and 39 comparison vs 27-34 controls; four-site gives 12 vs 6-14 and 6 vs 3-7. At .3 both hours have one three-site candidate and zero four-site candidates. Both waveform panels inspected; labels remain pending. `north_recent_review_v1` contains tables and figures. Conventional-event context is retained separately; neither strict candidate coincides with the listed nearby earthquake arrival times, which does not certify LFEs. No candidates promoted.

Archived B011 neural training was found and reused: refined-label experiment reported 0.024 pick precision/0.817 recall under its proxy labels, with random segment splits. The PGC RF report and the prior E1 Lin-resource audit were also read. Lin checkpoint absence was already documented; no reimplementation or archived training rerun was attempted.

A new bounded jms P/S adaptation experiment (`north_picker_adaptation_v1`) uses cached published stacks, spatially grouped source splits (2 km 3D connected groups), and phase-randomized noise spectra from disjoint portions of the August 20 panel. Eight epochs, no original weights overwritten. Unknown source phases are masked. Synthetic source-held-out test F1: P 0.2592->0.3185; S 0.4051->0.5884. Synthetic noise false P peaks increase 3->7; S decrease 3->1. These are development results, not field precision. Natural-waveform replay and a three-source held-out phase-proxy check are running; no production method has passed.

## Natural transfer failure and continuous network-score check

The synthetic adapter failed its real-waveform source-held-out proxy check (`north_adapted_natural_reference_v1`): zero P or S matches at score 0.3; original jms gives zero P and 1/432 S within 0.5 s (2/432 within 1 s). These are transferred phase proxies, not manual truth. Descriptive pooled proxy-versus-shift AUC at 1 s is P 0.629 baseline/0.584 adapted and S 0.678/0.661. Windows are correlated; these are not independent accuracy estimates. The synthetic adapter is excluded from production.

The adapted 2026 northern active panel has seven three-site associations at 0.1 versus 1-4 shifted controls, but zero at 0.3; the comparison has five versus 4-6. No candidates promoted. Original-model grid-preserving 30 s shift controls reuse observed associations exactly: at 0.3 the active candidate count is one versus 0-2 controls; the comparison is one versus zero in eight controls. Small control counts do not establish reliability.

The continuous original-model S-score scan (`north_network_probability_scan_v1`) uses 663 fixed source nodes, eight physical sites, per-site log-score normalization and a four-site support gate. At 4 MAD the active hour gives six proposals versus 3-21 per shift; comparison gives eleven versus 6-17. At 6/8 MAD there are no observed proposals in either panel while some shifts produce proposals. The scan does not discriminate the active panel. Its travel-alignment and nonwrapping time-shift sign test passed. Homogeneous-model origin times/coordinates remain unvalidated. No raw bytes downloaded.

A new fixed experiment (`north_natural_picker_adaptation_v1`) uses cached real waveforms and all 32 date-stable phase transfers. Training dates Sep 7/10/13/16, selection dates Sep 8/11/14, evaluation dates Sep 18/19/20, 2012. Source families may recur across dates; no source-holdout claim. Unknown-source phase windows are masked, with known-proxy local supervision and explicitly reference-free background contexts. Eight epochs, original jms initialization, no tuning on test results. Three new timing/masking tests passed. This is developmental evidence only and does not certify full-margin production.

The real-waveform adaptation finished: 3,557 training, 2,892 validation and 1,357 test windows. Epoch 4 selected by validation. Test P: baseline 0/62 matches with one unmatched scored peak; selected model 1/62 with 94 unmatched peaks. Test S: baseline 0/208 with three unmatched peaks; selected 0/208 with 88 unmatched peaks. The incomplete proxies do not define true precision, but there is no useful demonstrated recovery. Model not deployed and no redundant 2026 replay performed.

The real-data monthly retention replay completed (`monthly_panel_replay_v1`): August 20 then August 18, eight stations per one-hour unit, 178 individually hashed retained files per unit. Verified original acquisition receipts, full waveform and P/S core coverage, event support windows, controls and product schemas, restarted from ledger, then deleted 18,012,160 bytes of newly staged copies. All original evaluation raw files remain byte-identical. This isolated two-unit plan is not full-month or full-roster production. Monthly planner/storage/product tests: 21 passed.

The southern cross-family template audit completed: 18 previously untested cross-family pairs, no three-site high-CC alias flags, no families merged (`southern_cross_family_audit_v1`). All seven station identities in this panel are geographically distinct in the fixed roster; no nearby-site double count was introduced. Existing within-family clock curves were not rerun. Weak similarity is not proof that every candidate row is a distinct physical event.

The withheld PB.B039 cohort-stack check completed (`southern_withheld_stack_check_v1`). Primary nonoverlapping cohorts: family 026 five windows, CC 0.015 versus shifted maximum 0.186; 027 four windows, CC 0.046 versus 0.105; 029 nine windows, CC 0.102 versus 0.117; 031 one window, CC -0.029 versus 0.113. Family 029 is near the upper shifted tail (2/101 descriptive rank fraction) but modest in absolute waveform correlation; other groups are not strongly supported. Controls are correlated, four families/two cohorts were examined, and no significance or individual-event certification is claimed. Family 026 and 029 waveform figures inspected. All-rows stacks are also retained; no held-out station timing, polarity or detector threshold was optimized. No new raw bytes downloaded.

## Follow-up questions withdrawn, 2026-09-15

The user requested deletion of the two pending file-path questions (Washington catalog and existing reviewed labels). Cleared pending_input in ACTIVE_STATUS.json. Do not wait for answers or resend those questions. This withdraws the requests; it does not supply missing reference data or authorize author contact. No tool is available to delete previously displayed question cards from the application queue.

## Extended northern validation resumed, 2026-09-15

Running `north_extended_family_validation_v1`: all ten consecutive days September 21-30, 2012, submitted newest first, three worker processes. Original 47-group bank and fixed 9-MAD threshold; all 100 station-day files already cached, no downloads. Family-specific anchor clocks are now explicit in this separate validation output, and source hypotheses are not merged across incompatible clocks. Observed and 16 shifted realizations share each family's exact common valid exposure. This extends validation onto a different date block; it does not erase the previous failure or claim these dates were never examined by another project. Two tests passed for cross-family clock isolation, one-to-one matching and half-open coverage. Full-margin production has not started. The two withdrawn user questions remain withdrawn.

Complementary run `north_interleaved_family_validation_v1` includes all four unused calendar dates within the original September 7-20 design: September 9/12/15/17. Fixed 9 MAD, unchanged bank, two additional workers; no activity-based selection. `north_extended_B926_validation_v1` consumes each completed day from both runs, reusing the frozen 15-group independent-site bank and recording exact-lag signed correlations with 100 date-level shifts. No training repeats and no downloads. Missing reference rows are not treated as reviewed negative events.

While the validation scans run, retained-product verification now accepts explicitly bounded, evidenced waveform padding outside the core. Padding is hashed, cannot override core gaps, and does not change half-open event ownership. This closes the documented boundary-support limitation. Product/storage checks: 25 passed. Existing retained replay files were not rewritten.


## 2026-09-15: user-defined temporal scope and closed reference questions

User prefers January 2011–August 2026, with January 2019 fallback depending on network capability. Completed exact monthly metadata exposure and an updated descending 188-month / 21,490 station-epoch plan in `results/network_temporal_scope_v1`; 2019 fallback has 92 months. Physical-site metadata counts rise from 54 (January 2011) to 86.57 monthly mean (January 2019) to 159 (August 2026). Metadata is not actual waveform availability or regional sensitivity. Chart: `network_history.png`. Three exact exposure/coverage tests passed.

User confirms no additional reviewed P/S labels or validation results and no copy of LFEsAll.mat. Both previous input requests remain withdrawn; no user response is pending. Small 30-second waveform snapshots are now running via `python -m codex_lfe.temporal_waveform_probes`, in August 2026, January 2019, January 2011 order, two concurrent requests within only one active snapshot month. Maximum 2 MB per station request, all receipts and small validation snippets retained within codex-LFE. This is not full-month production.

Additional northern four-date validation completed: 80 family hypotheses / 75 same-family matches / 1,911 eligible references; worst shifted count 1 of 80, descriptive check passed, recall 3.92%, production false. Extended ten-day validation and frozen B926 checks continue independently.
