# Where the northern development references are lost

Follow-up: [the proposed training-weight experiment is complete](SITE_WEIGHTING_EXPERIMENT.md).
It did not improve calibration recovery; equal weighting remains the baseline.

Completed September 12, 2026 using saved neural probabilities/picks and cached
waveforms. The diagnostic took 75.9 seconds, including four parallel CPU date
workers for template correlations. No GPU inference, waveform download, or
detector retuning was performed.

The evidence points to insufficient detection scores and multistation pick
support under the current settings. It does not establish a fundamental signal
limit or show that a different picker cannot work. Association is not discarding
a large pool of correctly timed, four-site-supported references in this subset.

## Neural picking and association

We traced the same 188 calibration and 51 evaluation references used in the
[conditional comparison](INTERLEAVED_METHOD_RESULTS.md). The reference clocks
remain model-dependent and the dates remain development data. Each station's
modeled P/S arrival was checked against its saved model scores and exported
peaks, with full finite prediction coverage including the peak guard. Both P
and S diagnostics are saved; the table emphasizes S because the current
association settings require four physical sites to supply S support.

At score 0.3 and within ±2 seconds of the expected S arrival:

| Station | Calibration reference windows with an S peak / 188 | Evaluation reference windows with an S peak / 51 |
| --- | ---: | ---: |
| CN.PGC | 13 | 3 |
| CN.SNB | 29 | 15 |
| CN.VGZ | 12 | 3 |
| PB.B011 | 27 | 12 |
| PB.B926 | 22 | 10 |

Only **one calibration reference and zero evaluation references** have S peaks
at all four physical sites. PGC and B011 count as one site. The sole supported
calibration reference is `bostock_row_196360`; it is also the one reference
successfully matched by the existing four-station association result. None has
all five station IDs supplying S peaks at this threshold. The current association
requires at least four/five S picks as well as the subsequent station/site checks.

Widening the arrival window to ±5 seconds still gives one calibration and zero
evaluation references with four-site support at score 0.3. Raw probability
crossings give the same four-site counts. At individual stations, crossings and
exported peaks differ by at most one reference window in each split. Peak
extraction alone therefore does not explain the low support in these windows.

A diagnostic score of 0.1 raises four-site support within ±2 seconds to **5/188
and 2/51**. That is still limited. No new association run at 0.1 was performed
and this is not a selected operating threshold. Reference-centered windows may
share peaks during bursts, so these numbers are window-evidence counts, not
one-to-one single-station pick recall or verified events.

All references have full S-prediction coverage for these checks. Shifting the
expected arrivals by −300, −120, +120, and +300 seconds retains that coverage.
At score 0.3, none of the four shifted realizations has four-site support in
either split. At score 0.1 and ±2 seconds, calibration control counts are 0, 0,
1, 1; evaluation controls are all zero. These small shifted-window counts
describe background coincidence. They are not association controls, known
noise labels, or false-discovery estimates.

As a limited timing cross-check, direct S arrivals transferred from published
labeled stacks are available at PGC for all 239 references. They yield 13/188
and 4/51 windows with score-0.3 S peaks, compared with 13/188 and 3/51 for modeled
arrivals. Thus replacing the PGC model arrival with the transferred-stack arrival
does not remove the large loss there. Transfer coverage at the other CN stations
is sparse, and no equivalent borehole label transfer was established. Absolute
timing and velocity-model errors remain possible contributors.

[Per-station counts](../results/method_interleaved_miss_diagnosis_v1/phase_summary.csv),
[four-site support counts](../results/method_interleaved_miss_diagnosis_v1/event_support_summary.csv),
and [per-reference support](../results/method_interleaved_miss_diagnosis_v1/event_support.parquet).

## Template matching

We recomputed the existing training-template correlations around all **336
eligible calibration references**, using the original half-hour normalization,
gap/edge masks, four-site support, ±1-second reference tolerance, and frozen
threshold of 9 MAD. No thresholds or templates were changed.

| Stage | Reference count |
| --- | ---: |
| Network score never reaches 9 MAD in the reference window | 298 |
| Network score reaches 9, but simultaneous four-site support fails | 7 |
| Network score and four-site support pass together | 31 |

The 31 passing reference IDs exactly equal the 31 references matched by the
saved final calibration catalog. Peak suppression and final source assignment
cause no additional reference loss in this particular calibration result.
The dominant measured loss is the template-match score: 298/336, or 88.7%,
do not reach the network threshold. Removing the site-support gate alone could
only admit seven more reference windows here, before considering contamination.

This does **not** distinguish weak coherent signals from poor template shape,
phase-window choice, site weighting, or source variability. That distinction
requires waveform/template review and a controlled sensitivity experiment.
It also does not explain references excluded earlier by template QC or coverage.

[Stage counts](../results/method_interleaved_miss_diagnosis_v1/matcher_loss_summary.csv)
and [per-reference scores](../results/method_interleaved_miss_diagnosis_v1/matcher_calibration_reference_scores.parquet).

## Consequence for the next experiment

Prioritize waveform/template sensitivity over widening association tolerances.
Review recovered and missed reference waveforms against their training templates,
including VGZ and the borehole stations, to identify coherent evidence that the
current scores miss. Keep the existing waveform-review labels explicitly pending.

The next bounded method experiment should test one training-derived sensitivity
change, with the same calibration references and shifted-site controls. Claude's
station-weighting work is a concrete candidate to reuse after checking its
definition and transfer to these five stations. Do not lower thresholds or drop
a required site solely to improve reference counts; contamination must be measured
alongside recovery. Fresh dates and another region remain necessary before a
full-margin method choice. No GPU is needed for this next bounded stage.

## Reproduction and verification

Run from `codex-LFE`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
python -m codex_lfe.miss_diagnosis \
  --matcher results/method_interleaved_mf_v1 \
  --neural results/method_interleaved_neural_corrected_v1 \
  --shared results/method_interleaved_shared_corrected_v1 \
  --clocks results/method_interleaved_published_clock_ps_v2 \
  --association results/method_interleaved_association_corrected_v1 \
  --workers 4 --output results/method_interleaved_miss_diagnosis_reproduction
```

The [manifest](../results/method_interleaved_miss_diagnosis_v1/manifest.json)
records input/code hashes and waveform provenance. **85 software tests pass**,
including gap/day-boundary rejection, nearby-site deduplication, and simultaneous
score/site-support requirements. An additional artifact check verifies the exact
31 reference IDs and 45,345 unique phase-window records.
[Verification](../results/method_interleaved_miss_diagnosis_v1/verification.json).
These are software and consistency checks, not human event validation.
