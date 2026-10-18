# Independent CPU results — updated 2026-09-14

Southern work has advanced to a [four-station training bank](docs/SOUTHERN_PREPARATION.md) for one named family. Published references are imported separately; no southern detections yet. Latest software suite: 94 passed.

Latest catalog increment: [northern candidate catalog v1](docs/NORTH_CATALOG_V1.md), **61 candidates over 17 scanned days**. New active-date validation yields 12 candidates, 11/121 same-group reference matches and zero shifted detections across 16 controls. The five-day quiet block yields none, with five eligible references. Eight-second templates retained after the shorter-window experiment lost recovery. Scientific release remains unvalidated.

Latest follow-up: [calibration waveform/template review](docs/PACKET_REVIEW.md) covers all 336 eligible references and includes 19 review panels. Missed-reference median fixed-lag correlations are 0.083–0.184 across the four sites, above shifted-window medians near zero. Weak signals versus template mismatch remain unresolved; no detector change is adopted. Existing software checks: 90 passed.

Latest method experiment: [training-derived site weighting](docs/SITE_WEIGHTING_EXPERIMENT.md)
completed in 103 seconds on four CPU workers. Both arms select 9 MAD; weighted
matching recovers 29/336 calibration references versus 31/336 for equal weighting,
with identical reference eligibility and control exposure. **Retain equal
weighting.** No evaluation scan was performed. Current software tests: 90 passed.

Latest follow-up: [miss diagnosis completed](docs/MISS_DIAGNOSIS.md). At the
current neural threshold, only one calibration reference and no evaluation
references have four-site S-pick support, even with ±5-second arrival windows.
The one supported calibration reference was associated. For template matching,
298/336 eligible calibration references fall below the network threshold, seven
more fail simultaneous site support, and the remaining 31 match the saved
catalog. The next priority is waveform/template sensitivity, not wider association
tolerances or GPU scaling. Current software tests: 85 passed.

Latest: [the revised interleaved comparison](docs/INTERLEAVED_METHOD_RESULTS.md)
is complete using six CPU processes. The matcher passed calibration at 9 MAD and
recovered 13/342 eligible evaluation references. Corrected geographic-component
neural inference covered 50 station-days in 27.2 seconds; 28 observed/control
association runs followed. On a smaller common set with qualified model clocks
and exact prediction coverage, the matcher recovers 3/51 evaluation references
and neural association recovers none. These development results support retaining
the matcher as the known-repeat baseline, but do not establish a full-margin
method winner. Current software tests: 82 passed.

The [first shared-date comparison](docs/NORTH_2010_COMPARISON_SETUP.md) and
[legacy-input association](docs/REUSED_ASSOCIATION_RESULTS.md) remain historical
diagnostics. The latter reused borehole picks without verified horizontal
orientation; the corrected experiment above supersedes it for method assessment.
No bulk waveform download or GPU inference was needed. Response metadata was
downloaded to verify instrument orientation and gains.

**A margin-wide catalog is a plausible research objective. A reliable, uniformly complete catalog has not been demonstrated.** Published northern, central Oregon, Olympic Peninsula, and southern results provide regional starting points; see the primary sources in [the benchmark registry](docs/REGIONAL_BENCHMARKS.md). Reliability must be established by region, station era, source family, and detection strength. GPU access helps scale searches, but does not establish those properties.

## What was built and run

This folder contains an independent masked waveform loader, training-only template construction, site-balanced matched filtering, repeated shifted-site controls, final-event one-to-one reference evaluation, an unused-station diagnostic, and a published neural-picker inference baseline. The original experiments below use raw waveforms and public source tables, without importing a previous detector or corrected event clock. The later, explicitly labeled reuse experiment additionally consumes Claude's saved neural picks and travel-time table; its event association is rerun after roster restriction. Inputs are read-only, with provenance recorded.

The first three-group northern pilot failed its calibration gate and produced zero test candidates. Only eight reference events were eligible on its test dates. That is a failed, underpowered pilot—not evidence that full-margin detection is impossible.

The broader pilot selected twelve geographically separated coordinate groups using training data. Templates used September 7, 10, 13, and 16, 2012; calibration used September 8, 11, and 14; fresh test dates were September 18–20. Four detector sites were used. The threshold was frozen before the test scan. Completed output: [north_cpu_v2_parallel](results/north_cpu_v2_parallel/report.json).

| Frozen threshold 9 | Calibration | Test |
|---|---:|---:|
| Final candidates | 31 | 19 |
| Correct-source-group reference matches | 28 | 17 |
| Eligible reference events | 705 | 284 |
| Reference events excluded by eligibility | 255 | 103 |
| Shifted-control detections, summed over 16 scans | 1 | 3 |
| Largest exposure-adjusted control/observed ratio | 3.22% | 5.26% |

Test reference recovery is **17/284 = 6.0%**. Reference agreement is **17/19 = 89.5%**, which is not a precision estimate: the reference is incomplete, and coincident matches alone do not certify events. Seventeen of the nineteen candidates belong to a single source group. This demonstrates limited detection capability, not broad family coverage.

The calibration gate passed. On test data, the largest control ratio slightly exceeds the 5% calibration criterion: one control event versus nineteen observed events. Averaging controls gives roughly 1%, but neither number is a calibrated false-discovery rate. Sparse counts and the exposure normalization limit interpretation.

The LRIV station was excluded from detection. Only three of twelve groups passed its training template QC; none contained the nineteen test candidates. The held-out check therefore has **zero usable candidate observations and is inconclusive**, with only one usable reference positive-control observation. See [the full diagnostic](results/north_cpu_v2_heldout_lriv/report.json). This provides no independent confirmation or rejection of the candidates.

No candidate falls within ±120 seconds of an earthquake in the downloaded regional USGS table. This is only a descriptive proximity screen, comparing an arrival anchor with earthquake origin times; it cannot exclude conventional earthquakes or substitute for phase association.

## Calibration miss diagnosis

The raw-waveform diagnostic covers all **705 eligible calibration references**, using the existing ±1-second matching tolerance and frozen templates. Only **28** have a network score reaching threshold 9; the same **28** have a sample that also passes the three-site support rule. The final catalog recovers 28 correct-group references. Thus, at this operating point, neither relaxing site support nor changing final grouping promises additional reference recovery: the immediate bottleneck is subthreshold network correlation.

Across all references, 314 have no supported sample anywhere in that window, so site support does matter at weaker scores. The result does not prove that weaker detections are reliable or identify whether template averaging, station choice, waveform variability, or reference timing causes the low scores. Next work should compare training-derived template representations and station networks on calibration data, preserving fresh dates for evaluation. See [per-family measurements](results/north_cpu_v2_calibration_diagnostic/family_summary.csv) and [verification](results/north_cpu_v2_calibration_diagnostic/verification.json).

## Separated-training-date audit

An audit compares September 7/10 with September 13/16 at the original fixed template windows, without fitting time lags. Only three of twelve source groups have at least ten valid windows in each half: 12 station/group comparisons are measurable, while 34 lack sufficient temporal coverage. These missing comparisons are inconclusive, not failed coherence tests.

For all three measurable groups, PGC, SNB, and NLLB pass the existing 0.25 correlation bar; VGZ fails (separated-date correlations 0.220, 0.100, and 0.208). This motivates a calibration experiment with training-selected station membership. It does not yet justify removing VGZ from every family or claim improved recovery. The windows themselves were previously selected from all training dates, so this is a post-fit audit rather than fully independent validation. [Detailed results](results/north_cpu_v2_temporal_qc/temporal_qc.csv).

Template extraction now has an explicit sinc interpolation option. A synthetic half-sample-shifted 8 Hz signal is preserved within 0.0003 absolute error; legacy linear extraction reduces its RMS to under 40% in this example. This identifies an avoidable signal-processing loss, but improved field detection remains to be tested. Saved pilots retain legacy extraction and their original results.

## Station-removal experiment: no improvement

The separated-date QC hypothesis was tested by removing VGZ only from the three measured failing source groups. All twelve groups, original templates and timing, calibration dates, thresholds, and 16 control realizations were retained. The scan completed in 149 seconds; no test dates were processed.

All 705 reference events remain eligible in both variants. At threshold 9 both yield 31 candidates and 28 correct-group matches, with the same maximum control ratio (3.22%). At threshold 8, matches decrease from 58 to 56 and the maximum control ratio increases from 7.15% to 7.36%. No tested threshold improves shared-reference recovery. **Do not adopt this station-removal rule on this evidence.** Temporal template coherence alone was insufficient to predict a beneficial removal.

[Full comparison](results/north_cpu_v2_station_ablation/comparison.csv), [shared-reference comparison](results/north_cpu_v2_station_ablation/shared_reference_comparison.csv), and [recorded decision](results/north_cpu_v2_station_ablation/decision.json). Next: test band-limited template extraction with original station membership, rather than combining it with an unsupported station change.

## Verification and reproducibility

`python -m pytest tests -q`: **57 passed**. Tests include known multistation injections with gaps and a strong single-site artifact, duplicate-safe matching, time-basis checks, masked FFT correlation, torch CPU parity, metadata epochs/orientation, and fractional-sample alignment. They verify software behavior, not field precision or completeness.

The broader scan completed in approximately 335 seconds with three CPU day workers and reusable preprocessing caches. This is not a full-margin throughput forecast. The neural smoke test ran one station-hour and produced one unassociated pick; this demonstrates executable inference only.

Saved pilots used linear fractional-sample alignment. The current generic loader uses tested band-limited sinc alignment; existing pilot configs retain legacy linear alignment for reproducibility. New experiments should explicitly choose `fractional_alignment: sinc`. CUDA has not been exercised on this pod.

## What must happen before a reliable margin-wide release

1. Diagnose missed reference events and improve training templates/site selection, then freeze and test on new dates. Alternating-event template QC is not an independent temporal holdout because nearby waveform windows can overlap.
2. Obtain machine-readable central Oregon and Olympic benchmark events/templates. Oregon channel metadata has been acquired for the historical 2009 and a modern 2024 window; waveform availability and regional detection performance remain unverified. No new waveforms were downloaded.
3. Implement and validate new-family discovery, multistation association, and location. The neural picker currently emits unassociated picks. Known-family matching alone cannot establish coverage of unknown sources.
4. Establish useful unused-station confirmation, blinded event review, earthquake discrimination, and regional completeness/contamination estimates. Preserve failed or inconclusive diagnostics.
5. Replace the pilot's block-edge exclusions with validated padded scanning; the present guards discard about 21% of nominal time. Improve merged-event exposure accounting beyond summed family-seconds when station/family availability varies.
6. Recover original family identities and validated origin times, then calibrate amplitudes/magnitudes where data permit. Current event times are PGC training-stack arrival anchors, and family IDs are exact published-coordinate groups.

The next GPU session can reproduce the matcher and run neural discovery experiments using [the handoff](docs/GPU_HANDOFF.md). It should not simply scale these nineteen experimental candidates into a claimed complete catalog.

## Template interpolation calibration

Sinc extraction alone was tested with original station membership and legacy daily waveform alignment. Calibration selected threshold 9 again: 30 candidates and 27 correct-group matches, versus 31 and 28 for linear extraction. Threshold 8 gives 71 candidates and 59 matches versus 70 and 58; its maximum control ratio remains above 5% (7.05%). All 705 references are shared. This does not demonstrate a useful sensitivity gain at the selected operating point. No test dates were scanned. The signal-preservation unit test remains valid; field detection improvement is not established.

The user has clarified that the intended detection network is their existing 198 long-term stations. Exact roster identification takes priority over expanding these pilots to external or temporary stations.

The pilot now accepts an explicit station-roster CSV and rejects outside stations before waveform loading. Five additional tests cover full network identity, outside stations, duplicates, missing networks, wildcards, and empty rosters. The exact 198-station input remains pending; enforcement capability is implemented, but network membership of prior pilots is not yet verified.

## Network roster resolved from the user-designated map

The current `figures/borehole_dvv_map.html` has 183 station markers, mapped to 183 unique network/station identities: BK 21, CC 11, CN 18, NC 13, PB 30, UO 32, US 1, UW 57. The user clarified that their original 198 estimate may be high. This exact map is now the authoritative detection roster, superseding the earlier pending-roster status. PGC's network was resolved from independently acquired StationXML; other codes have unique exact matches in project station tables. Physical epoch verification remains separate.

PGC, VGZ, and SNB are in the roster; NLLB and LRIV are not. Earlier pilot results used NLLB and therefore do not establish performance of the target network. A three-station northern calibration config is prepared and membership-validated, but has not yet been scanned. See [roster](results/target_network_from_map/stations.csv) and [provenance](results/target_network_from_map/manifest.json). No prior map certification labels were reused to validate LFEs.

## First calibration constrained to the 183-station roster

`map_network_north_calibration_v1` uses only CN.PGC, CN.VGZ, and CN.SNB. Its manifest records the authoritative roster hash and membership. Eleven of twelve training-selected source groups pass template QC. Calibration selects threshold 9: **27 candidates, 23 correct-group reference matches, 648 eligible reference events, and zero detections across 16 shifted controls**. Recovery is 3.55% on eligible references; 312 of the 960 selected-group references are excluded (including the group without adequate templates). Zero controls do not establish zero false alarms. No test dates were scanned. At threshold 8, 52 candidates include 40 correct-group matches, but the worst control ratio is 9.60%, above the predeclared criterion.

This establishes a runnable baseline within the actual target network. It does not validate full-margin precision or completeness. [Calibration scores](results/map_network_north_calibration_v1/calibration_scores.json).

## Independent channel-metadata audit of all target stations

New provider metadata was acquired for exactly the map roster over 2009-01-01 through 2026-09-11. All **183** stations have at least one qualifying three-component epoch with independent orientations and sample rates at least 20 Hz. In 2012, **63** stations have qualifying epochs; in 2024, **177** do, with **172** covering at least 90% of that year's time in metadata. Fractions use unions of overlapping instrument alternatives, not their sum. This is metadata coverage, not waveform uptime or LFE sensitivity.

BK/NC initially rejected date-only timestamps. A retry using complete timestamps succeeded; the failure and retry records are preserved. No waveforms were downloaded. Existing raw day files in `lfe_catalog_repair/traces` cover 16 roster stations, mostly in the north, but their full contents have not all been verified.

See [regional metadata coverage](results/target_network_from_map/regional_metadata_coverage.csv), [station-year fractions](results/target_network_from_map/station_year_metadata_coverage.csv), and [audit provenance](results/target_network_from_map/metadata_audit.json). Station-location bins are planning aids, not source detectability maps. A modern-period regional pilot should be evaluated before assuming this network provides comparable historical sensitivity.

## Picker benchmark tooling and unreviewed packet

Built a development review interface and strict P/S scorer, with 70 total tests passing. The 24-window northern packet has no human labels. Both jms and jcms inference complete and cover all core windows, but scoring the pending review correctly fails. This is engineering progress only; no picker accuracy or winner has been established. See [validation protocol](docs/PICKER_VALIDATION.md) and [verification](results/picker_review_development_v1/engineering_verification.json).
