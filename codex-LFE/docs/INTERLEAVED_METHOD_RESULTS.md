# Interleaved northern method comparison — September 12, 2026

Follow-up: [the miss diagnosis is now complete](MISS_DIAGNOSIS.md). It traces the
losses through saved neural scores/picks and recomputed template correlations;
this report preserves the preceding comparison results.

The revised CPU comparison is complete. Keep matched filtering as the known-repeat
baseline; neither arm has demonstrated enough recovery or independent event
validation to select a full-margin production method. The corrected neural arm
remains a discovery candidate, with its missed-reference and location behavior
requiring diagnosis. These are development results on previously examined dates.

Six date processes completed the matched-filter scan in 254.5 seconds and corrected
neural inference for 50 station-days in 27.2 seconds. Matching used two FFT threads
per process; inference used two Torch threads. No GPU or new raw-waveform download
was needed. Small response-level StationXML files were downloaded to establish
instrument gains and orientation. Timings exclude training and metadata acquisition
and are not full-margin throughput estimates.

## Fixed network and dates

The detector uses CN.PGC, CN.SNB, CN.VGZ, PB.B011, and PB.B926, all inside the
183-station map roster. PGC and B011 share one physical-site vote and control
shift, leaving four site groups. This northern geometry does not represent the
whole margin. SNB uses its actual 2010 BH channels.

The [configuration](../configs/method_comparison_north_2010_interleaved.yaml) uses:

- Training: August 19, 23, 27, 31 and September 4, 2010.
- Calibration: August 21, 25, 29 and September 2.
- Evaluation: August 22, 26, 30 and September 1, 3, 11.

Twenty-four source-coordinate groups were selected from training counts alone.
Eleven passed the existing four-site template QC. Selected-source evaluation
references shrink from 1,411 before template QC to 434 in usable source groups
and 342 after matcher masks. Coordinate groups are not authoritative published
family identities. No references on September 11 does not establish a noise-only day.

## Matched-filter result

Threshold 9 MAD passed the calibration gate and was frozen before evaluation.
The gate requires at least 20 calibration candidates and a maximum shifted-control
ratio no greater than 5%; this inherited engineering rule is not a demonstrated
catalog precision requirement.

| Frozen threshold 9 | Calibration | Evaluation |
| --- | ---: | ---: |
| Eligible source-matched references | 336 | 342 |
| Candidates | 34 | 15 |
| One-to-one same-source matches | 31 | 13 |
| Reference recovery | 9.23% | 3.80% |
| Control detections summed across 16 realizations | 1 | 1 |

Evaluation reference agreement is 13/15, **not precision**. An unmatched candidate
is not automatically false; a temporal reference match does not independently
certify identity. Lowering to 8 MAD failed the calibration gate: 56 candidates,
49 reference matches, and up to three detections per control (5.36%). Its
evaluation result is exploratory, not a replacement frozen operating point.

[Matcher report](../results/method_interleaved_mf_v1/report.json) and
[threshold selection](../results/method_interleaved_mf_v1/selection_before_test.json).

## Why the neural inputs were rerun

The legacy E2 path passed native streams to SeisBench `annotate` without inventory
rotation. Flexible component mapping aliases channel 1 to N and 2 to E. Actual
StationXML azimuths are 11°/101° at B011 and 299°/29° at B926. The latter is a
61° discrepancy from nominal north for channel 1. Reusing those inputs unchanged
would not establish performance with verified geographic components.

The new CPU path verifies component response compatibility, equalizes gains,
rotates to geographic ZNE, preserves gap masks, and saves prediction coverage and
probability arrays. It also changes filtering/stitching and extracts local maxima
instead of the legacy threshold-interval peaks. **This is not an isolated rotation
ablation, so differences cannot be attributed to orientation alone.** The model
weights are the existing published Cascadia LFEDetect weights; there was no training.

[Inference manifest](../results/method_interleaved_neural_corrected_v1/manifest.json),
[CPU driver](../codex_lfe/parallel_neural_panel.py), and response inventories
[CN](../cache/metadata/cn_north_2010_response.xml) /
[PB](../cache/metadata/pb_b011_b926_2010_response.xml).

Association was rerun for all observed and shifted configurations with four
physical sites supplying S support. The existing Crosson model and 20–50 km
depth bounds remain assumptions. Matching tolerance is 1.2 seconds;
`edt_pick_std=0.5` is inactive with `exponential_edt=False`. Observed and control
events exclude the first/last 3,120 seconds of each date, covering shifts up to
3,000 seconds plus association overlap without depending on unprocessed neighbours.

| Pick threshold / station IDs | Observed | Constant-site controls, 3 seeds | Hourly-site controls, 3 seeds | Within 1 km of depth bounds |
| --- | ---: | --- | --- | ---: |
| 0.3 / 4 | 71 | 9, 6, 10 | 1, 6, 5 | 24 |
| 0.3 / 5 | 35 | 1, 1, 0 | 0, 2, 1 | 10 |
| 0.5 / 4 | 15 | 0, 0, 1 | 0, 1, 1 | 7 |
| 0.5 / 5 | 4 | 0, 0, 0 | 0, 1, 0 | 0 |

For 0.3/five stations, mean control counts are 1.90% and 2.86% of observed.
Sparse controls do not prove that 97–98% of events are real. Boundary locations
need waveform and model checks; these associations are candidates.
[All 28 runs](../results/method_interleaved_association_corrected_v1/runs.csv).

## Conditional comparison on the same references

Published reference timestamps and matcher arrival anchors do not directly equal
neural-association origin times. Single-event training-clock fitting accepted no
multisite clocks. We therefore transferred manually observed P/S labels from the
local published-stack archive to our training-only stacks. Predicted archive
labels were not substituted for observed labels.

Source joins require reciprocal nearest coordinates within 1 km and a 0.25 km
ambiguity margin. Phase transfer requires correlation at least 0.6 and a peak
advantage of at least 0.1. P uses the vertical component; S uses three components.
Within-site P/S estimates must agree within 1 second. Family clocks require at
least two independent sites, with one vote per site and spread at most 1 second.
Travel-time subtraction uses the Crosson model. Six source groups qualify.

These are **model-dependent development clocks**, not validated absolute origin
times. There are no measurable early/late stability comparisons among the accepted
phase transfers: the required source activity does not span both training halves.
Published stacks may include these years, and pretrained weights may include
Cascadia sources. This is not an independent phase-picking benchmark.
[Clock provenance](../results/method_interleaved_published_clock_ps_v2/manifest.json).

Shared scoring requires usable matcher intervals, an extra 60-second context
margin, the same daily 3,120-second guard, and exact neural S-prediction coverage
at all five stations around modeled arrival times ±5 seconds. This leaves 188
calibration references across five groups and 51 evaluation references across
three groups. Fifty of the 51 evaluation references come from two groups.

| Arm / setting | Calibration matches / 188 | Evaluation matches / 51 |
| --- | ---: | ---: |
| Matched filter, frozen 9 MAD | 20 | 3 |
| Neural + association, 0.3 / 4 stations | 1 | 0 |
| Neural + association, 0.3 / 5 stations | 0 | 0 |
| Neural + association, 0.5 / 4 stations | 0 | 0 |
| Neural + association, 0.5 / 5 stations | 0 | 0 |

These are one-to-one temporal matches within 2 seconds. At 1 second, matcher
evaluation matches fall to two; at 5 seconds they remain three. Neural evaluation
matches remain zero at all three tolerances. The one neural calibration match is
2.62 km from the reference coordinates. Matcher locations inherit template
coordinates, so zero separation is not independent location validation.

The candidate inventories and controls differ between arms. This table measures
conditional reference recovery; it does not compare precision at equal contamination,
validate discovery outside the template bank, or establish a universal winner.
The wider matcher result (13/342) and this subset (3/51) have different denominators.

[Full comparison](../results/method_interleaved_shared_corrected_v1/reference_agreement.csv),
[matched pairs](../results/method_interleaved_shared_corrected_v1/matched_pairs.csv), and
[coverage/clock manifest](../results/method_interleaved_shared_corrected_v1/manifest.json).

## Decision and next experiment

Keep the current template matcher as the reproducible known-repeat baseline.
Do not scale the current picker/association settings into a claimed reliable
catalog on this evidence. The concrete next comparison is to diagnose where
reference events disappear: single-station phase evidence, four-site support,
association tolerance/travel-time error, or template mismatch. Start with these
saved probabilities, picks, and existing waveforms; another inference run is not
required for that diagnosis.

Review a balanced set of recovered references, missed references, neural-only
candidates, and shifted controls with method identity hidden during labeling.
Obtain waveform-supported phase labels and retain ambiguous cases. Use the result
to choose one targeted change, such as training-derived station weighting or
template extension, then freeze it for new dates and a second region. Template
matching alone cannot discover families absent from its bank; a full-margin
method must also demonstrate a validated way to build and extend that bank.

## Verification

`python -m pytest tests -q`: **82 passed**. An additional artifact audit checked all
28 association runs and 948 assignments against original pick identities, phases,
scores, and shifted times, including roster membership, four physical-site S
support, station count requirements, and daily guards.
[Verification record](../results/method_interleaved_shared_corrected_v1/verification.json).
These checks validate software/provenance, not scientific event identity.
