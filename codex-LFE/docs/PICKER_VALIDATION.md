# Picker validation protocol — development version 1

The current priority is choosing the method through independent validation, before fleet-scale detection. A working model, a good-looking probability trace, or a high match fraction against a non-exhaustive event catalog does not settle picker quality.

See the [Claude test audit](CLAUDE_TEST_AUDIT.md) before scheduling new experiments. Continuous published-model inference, corrected injections, and association controls already exist. The review packet below supplements that evidence; it should not delay reusing existing picks for roster-restricted association and corrected reference scoring.

## Current executable deliverable

- `results/picker_review_development_v1/review.html`: offline review interface with 24 initially unlabeled windows from CN.PGC and CN.SNB. Six reference-enriched and six random times are shown at both stations, in randomized order. Each scored interval is 60 seconds with 30 seconds of context on either side. All inputs belong to the user's map roster.
- `review_template.json`: pending review records, no inferred P/S labels. `sampling_key.json` contains the sampling origin; keep it out of the initial review. The HTML initially hides source type, station identity, UTC date, and model predictions. This is display blinding, not cryptographic concealment.
- `waveforms.npz` and `manifest.json`: original waveform inputs, preprocessing, source hashes, and sampling design. Windows share paired time-block IDs across stations.
- `picker_review_jms_cpu_v1` and `picker_review_jcms_cpu_v1`: CPU inference on the same padded waveforms. Full score arrays, low-threshold picks, prediction-coverage intervals, weight metadata, response/orientation checks, and provenance are saved separately from the review display. Both models cover all 24 scored windows. No accuracy metrics have been calculated.

This packet uses September 8, 2012, already used for northern matcher calibration. It is a development/workflow test, not held-out performance evidence. It contains no boreholes or central/southern data. The 12 random station-minutes are far too little exposure to establish rare false-pick rates. Reference-enriched examples were sampled from published records, not selected by either competing neural picker. Random windows are not presumed negative.

## Review procedure

1. A reviewer initially sees waveforms with no picker scores, source timestamps, or sampling labels. Inspect the scored interval and use the surrounding context.
2. Mark identifiable P or S arrivals, nominal UTC time (entered through relative-time clicks), timing uncertainty, and event class (`lfe` or `ordinary_earthquake`). Record ambiguous identity or phase as an ambiguous window instead of inventing a label. An ordinary-earthquake phase can be correctly timed yet still be a false LFE proposal.
3. Set each window to `complete`, `ambiguous`, or `unusable`; untouched windows remain `pending`. Complete means the full scored interval was reviewed. A complete window can have no identifiable LFE arrival. Reviewer identity is required. Editing labels returns the window to pending.
4. Obtain a second independent review for the decision benchmark, then adjudicate disagreements using multistation waveform context while still concealing competing model identities/scores. Preserve both original reviews and the adjudication provenance. The initial single-station interface does not itself establish LFE identity or perform adjudication.
5. Export JSON; keep labels versioned. Do not overwrite uncertain judgments with catalog-derived arrival predictions. Catalog origins and PGC template anchors are not P/S truth.

## Scoring contract

`codex_lfe.picker_evaluation` accepts predicted picks, explicit model-coverage intervals, and review JSON. It rejects overlapping station review windows, invalid clocks, duplicate review IDs, missing reviewer identities for complete reviews, and labels outside their windows. It scores only completed windows with full P and S prediction coverage. Missing predictions caused by missing inference coverage are not treated as negatives. All models must ultimately be compared on the same usable window intersection; coverage exclusions remain a separate metric.

Matching is one-to-one within the same station, phase, and window. Duplicate predictions cannot increase recall. Label uncertainty wider than the chosen tolerance excludes the whole window, preventing the removal of an uncertain label from silently creating a false negative example. Report exclusions alongside scores; do not interpret the remaining identifiable-arrival cohort as completeness for all LFEs.

P and S are reported separately: recall, review-panel precision, unmatched picks per reviewed station-hour, signed timing bias, median absolute timing error, and 90th-percentile absolute timing error. Random, reference-enriched, and picker-enriched panels remain separate. Unmatched picks are false relative to the reviewed identifiable-LFE label set; unreviewed weak signals can still exist. The first implementation does not yet estimate phase confusion or clustering-aware confidence intervals; wrong phases produce an unmatched pick and a missed arrival. Ambiguous maximum-cardinality timing assignments are another stated limitation.

Development timing comparisons: report all three tolerances, **0.25, 0.5, and 1.0 seconds**, rather than choosing the most favorable after seeing results. Review uncertainty can change eligible windows between tolerances; report that change explicitly and use identical cohorts between pickers at each tolerance. Final timing requirements must support association/location accuracy and be frozen before test evaluation.

Inference retains picks down to model score 0.05. Do not request lower-threshold metrics from these exported pick tables; regenerate from retained probabilities if needed. For a first curve, use the predeclared score grid 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9. Equal model score is not an equal operating point: choose comparable false-pick rates on development data, then freeze each model's threshold before the final test.

Example after independent review is exported:

```bash
python -m codex_lfe.picker_evaluation \
  --picks results/picker_review_jms_cpu_v1/picks.parquet \
  --coverage results/picker_review_jms_cpu_v1/prediction_coverage.parquet \
  --review path/to/adjudicated_review.json \
  --threshold 0.3 --tolerance-s 0.5 --split development \
  --output results/review_jms_score_new
```

The pending template intentionally fails scoring; it is not a set of negative labels. The same runner works for jcms with the corresponding inference files. This example threshold is not a selected operating point.

## Expand to a decision benchmark before choosing weights

- Select north, central, and south panels within the 183-station roster, including borehole/surface instruments where available, more than one noise/activity condition, ordinary earthquakes, and recording artifacts.
- Include randomly sampled continuous blocks for operational false-pick estimates and candidate-enriched panels for weak-arrival inspection. Keep their sampling weights/strata explicit. Add examples missed by both pickers through independent random review.
- Freeze development/test dates and hold out whole event/time blocks and source families; colocated instruments are not independent validation sites. Audit pretrained-weight provenance separately. `jms` metadata excludes Cascadia, whereas `jcms` includes it; a clean local split does not erase pretrained exposure.
- Choose annotation/exposure sample sizes from the target uncertainty and tolerated false-pick rates before final scoring. The 24-window engineering packet is not that final sample-size decision. Use station/day or event-block uncertainty estimates rather than treating all picks as independent.
- For each picker, run the same association/location baseline at development-selected operating points. Check predicted arrivals at withheld roster sites, location stability, conventional-earthquake confusion, and identical shifted controls through the complete associator. Small residuals alone do not validate associations.
- Add independent-waveform injections into real noise to measure controlled recovery; retain natural-event validation as the deciding evidence.

Remaining work: central/southern and borehole review packets, expert reviews/adjudication, shared-cohort operating-point selection, source-family split audit, phase-confusion diagnostics, clustering-aware uncertainty, and physical association/location validation. No picker winner is declared.
