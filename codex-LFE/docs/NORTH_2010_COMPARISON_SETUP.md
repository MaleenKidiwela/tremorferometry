# Shared northern development panel — September 12, 2026

Historical first-panel report. The revised panel has since been scanned, and
corrected neural inputs and shared-reference scoring are complete; see
[the current results](INTERLEAVED_METHOD_RESULTS.md). Statements below about an
unscanned revised configuration describe the earlier state.

This experiment advances method selection. It does not authorize bulk picking
or select the neural, matched-filter, or hybrid arm in advance.

## Completed result and revised next panel

The first matcher run completed on CPU in 487.5 seconds. Its calibration gate
failed at every operating point. At 9 MAD it recovered two of 79 eligible
known-source references, with no detections in the 16 shifted controls; two
candidates are below the inherited minimum-count gate. At lower thresholds,
control detections rose sharply. This is not a production operating point.

The first date split was unsuitable for ranking methods: of 1,119 reference
events on evaluation dates, only seven belong to the twelve selected groups,
four belong to usable template groups, and three survive waveform/guard masks.
Consequently its zero known-source recovery on evaluation cannot establish a
general sensitivity limit. I should have checked this source overlap before
running the scan.

On the same dates, the neural threshold-0.3/four-station run has 22 calibration
and six evaluation candidates; five of the six evaluation candidates lie near
the imposed depth bounds. These unreviewed counts do not rank it above the
matcher. [All diagnostic tables and limitations](../results/method_comparison_north_2010_summary_v1/report.json)
are saved, including the full reference/source-inventory counts.

A revised, **unscanned** development configuration is ready:
[method_comparison_north_2010_interleaved.yaml](../configs/method_comparison_north_2010_interleaved.yaml).
It distributes training dates across the episode and selects 24 groups using
training counts only. This is retrospective catalog construction: templates may
use dates before and after an evaluation date. The reference-only preflight
finds **939 calibration and 1,411 evaluation events in those selected groups**,
before template QC and waveform exclusions. The 3,564 evaluation events outside
that inventory form a separate discovery question; they are not automatically
new sources relative to published or pretrained catalogs. This revised design
uses already examined development data, not a fresh final test.

The immediate next work is to inspect training-stack phase timing, verify usable
template/source overlap for the revised panel before scanning, and then compare
known-repeat recovery separately from discovery outside the template inventory.
The [nine multistation training-stack panels](../results/method_comparison_north_2010_stack_review_v1/index.html)
and numerical stacks are ready for phase review. Two panels were visually
inspected: some station packets repeat clearly across the date halves, while
other channels are much noisier. This is not an adjudicated P/S-label set; no
phase labels or accepted origin clocks were created.

Verification: 78 software tests passed. The report collector also completed
against actual nonempty and empty association outputs. No GPU inference or raw
waveform download was required.

## Fixed inputs and splits

Configuration: [method_comparison_north_2010.yaml](../configs/method_comparison_north_2010.yaml).
The date split was written before the new matched-filter scan:

- Template building: August 19, 21, 23, 25, 2010.
- Calibration: August 26, 28, 30.
- Evaluation: September 1, 3, 5, 11.

All inputs come from CN.PGC, CN.SNB, CN.VGZ, PB.B011, PB.B926. PGC/B011 share
one site weight and one control shift. Templates use each instrument's native
three channels, including the actual 2010 BH channels at SNB. The 2012 HH
channel selection was not copied onto 2010 data. Four distinct sites are required.

Reference counts were inspected to ensure available training/evaluation examples;
families are selected from training counts only, with the existing 5 km separation
rule and maximum of twelve source groups. Nine groups pass training-template QC
at four sites. Three exclusions remain in the reported coverage limitation.

These dates have already been examined in Claude's experiments and in the
35-day neural-association run. They are separated for fitting this matcher, but
are **development data**, not a clean final test of method generalization.
September 11 has no entries in this reference table; that does not make it a
known noise-only day.

## What is being measured

The matched-filter arm uses the existing 8 s template construction, explicit
gap masks, four-site support, equal site weights, thresholds 5–12 MAD, and 16
site-shift realizations per scan. Threshold selection uses calibration alone.
The existing 5% maximum-control-ratio/20-detection gate is an engineering rule
inherited from the pilot, not an agreed scientific precision requirement.

The neural comparison reuses the completed grouped-site association output on
these dates. Its thresholds are 0.3/0.5 with four/five station IDs and four
physical-site groups. No new neural inference is needed for this stage.

Shared dates and stations remove major confounders, but do **not** by themselves
produce equal-exposure, equal-contamination operating points. The matcher masks
gaps and 190 s at each half-hour edge and uses 60–180 s control shifts; neural
association currently has incomplete inference-coverage records and 300–3,000 s
station/hour shifts. Report these differences explicitly. Do not rank methods
by comparing unadjusted total candidate counts or control percentages.

## Training-only clock test

The published Bostock timestamps have an unresolved physical time convention.
The matcher uses a training-stack arrival peak; neural association reports a
model-dependent origin time. An empirical correction fitted on all 35 days
would leak evaluation information into a comparison.

`codex_lfe.training_clock` instead restricts calibration to the four training
dates and S picks. It estimates per-source/per-station residual-histogram peaks
relative to the existing Crosson S travel-time table, independently for the
early and late training halves. It requires at least five counts in the peak
bin, contrast at least three, and agreement within 1 s. A family origin clock
additionally needs two independent site groups agreeing within 1 s.

The experiment found **four passing single-station/source pairs, zero passing
multisite family clocks**. No event-origin recovery ranking will be computed
from these unaccepted corrections. This may reflect limited/weak training picks
or model/phase inconsistency; it does not establish a fundamental limit.
The QC cutoffs are explicit development choices, not calibrated significance
tests. Source hashes, exact reference IDs, station diagnostics, and rejected
family clocks are retained in
[the clock output](../results/method_comparison_north_2010_clock_v1/manifest.json).

Next timing evidence should come from reviewed phase arrivals on coherent
training stacks and multistation moveouts, with the uncertainty preserved.
Do not relax this check after seeing the evaluation recovery solely to obtain
a favorable comparison.

## Commands

Run from `codex-LFE`:

```bash
python -m codex_lfe.pilot \
  --config configs/method_comparison_north_2010.yaml \
  --output results/method_comparison_north_2010_mf_v1 --workers 3
python -m codex_lfe.training_clock \
  --config configs/method_comparison_north_2010.yaml \
  --output results/method_comparison_north_2010_clock_v1
```

Neither command downloads waveforms or requires a GPU. Both require fresh output
directories to avoid silently replacing an earlier experiment.
