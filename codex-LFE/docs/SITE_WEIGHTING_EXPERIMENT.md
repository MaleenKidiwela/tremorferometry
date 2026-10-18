# Training-derived site weighting — September 12, 2026

**Completed: retain equal weighting.** The weighted arm recovers 29/336
calibration references at its selected threshold, versus 31/336 for the baseline.
Both select 9 MAD under the same engineering gate. The scan took 103.1 seconds
using four date workers and two FFT threads per worker. No evaluation dates were
scanned. This weight rule does not demonstrate a sensitivity gain on this panel.

This calibration-only experiment tests whether stronger training signals should
receive more weight in the network correlation. It follows the
[miss diagnosis](MISS_DIAGNOSIS.md), which found that 298/336 eligible calibration
references did not reach the existing network threshold.

## What is reused, and what changes

Claude's E5/E5r code uses a training-stack amplitude ratio as a weight, with a
cap at three times the median and median normalization. The formula is explicit
in `lfe_catalog_repair/scripts/e5r_templates.py`: energy in a four-second packet
window divided by energy in a window 15 seconds earlier, converted to an amplitude
ratio. We reuse that idea, not the old fitted weights or performance claims.

Our packet center is the frozen training-template energy peak. It has not been
independently labeled S, so the measurement is called **stack contrast**, not
calibrated single-event SNR. The preceding window can contain signal. Weighting
by this contrast is a hypothesis tested on calibration data, not an assumed optimum.

For each of the 11 usable source groups and five instruments, we reconstruct the
full training stack from saved review stacks. Those files omit a date half when
it contains only one event; the missing singleton windows are restored from the
cached training waveforms. This restores the exact original training population.
All 55 reconstructed template windows agree with the frozen templates to a
maximum absolute error of 4.44e-16.

Station contrasts are averaged within each physical site. The 44 resulting
source/site weights are capped at three times their global median and normalized
by that median. PGC/B011 retain their existing mean correlation, one site weight,
and one shared shift. No station is removed. The largest site contribution is
36.0% of a source group's total weight; none receives a majority.

The sole detector change is a weighted mean across physical sites instead of
the equal mean. The resulting network statistic is normalized by its own block
MAD. All site-support checks still use the existing unweighted per-site series.
Observed and all 16 shifted-site controls receive the same frozen weights.

Templates, source groups, arrival anchors, four-site support, channel processing,
half-hour blocks, gaps/guards, thresholds, reference matching, and control seeds
are inherited from the baseline. Only the four calibration dates are scanned.
No evaluation dates, new inference, or waveform downloads are involved.

[Frozen design](../results/method_interleaved_site_weighting_v1/design_before_scan.json),
[station contrasts](../results/method_interleaved_site_weighting_v1/training_station_contrast.csv),
and [physical-site weights](../results/method_interleaved_site_weighting_v1/training_site_weights.csv).

## Comparison rule

Compare one-to-one same-source reference recovery under the same existing
calibration gate: at least 20 detections and a maximum exposure-adjusted
control/observed ratio no greater than 5%. Also report both arms at the original
9-MAD operating point. These empirical ratios are not true false-discovery rates
or catalog precision. A gain here would warrant further validation, not automatic
production adoption. These dates have already been examined.

## Reproduction

Run from `codex-LFE` with a fresh output directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
python -m codex_lfe.site_weighting \
  --run results/method_interleaved_mf_v1 \
  --stacks results/method_interleaved_stack_review_v1 \
  --workers 4 --output results/method_interleaved_site_weighting_reproduction
```

The software suite has **90 passing tests**, including exact equal-weight
regression, weight-scale invariance, invalid-weight rejection, and rejection of
a large single-site artifact under four-site support. These verify software
behavior, not scientific event identity.

## Completed result

| Threshold | Equal-weight candidates / matches | Weighted candidates / matches | Equal maximum control ratio | Weighted maximum control ratio |
| --- | ---: | ---: | ---: | ---: |
| 8 MAD | 56 / 49 | 53 / 46 | 5.36% | 7.55% |
| **9 MAD, selected in both arms** | **34 / 31** | **32 / 29** | **2.94%** | **3.12%** |
| 10 MAD | 20 / 18 | 21 / 19 | 5.00% | 4.76% |

All 336 eligible references and every observed/control exposure row are identical
between arms. At 9 MAD, both have the same single control detection, in the same
one of 16 realizations. The weighted arm gains two matched reference IDs and loses
four; 27 are shared. The result therefore is not a uniform rescaling of candidate
scores. At 10 MAD there is a one-match gain, but choosing that threshold does not
improve on the baseline's 31 matches at its accepted operating point.

[Full seven-threshold comparison](../results/method_interleaved_site_weighting_v1/comparison.csv),
[decision](../results/method_interleaved_site_weighting_v1/decision.json), and
[changed reference IDs](../results/method_interleaved_site_weighting_v1/changed_reference_ids.json).

No weighting change is adopted in the baseline configs. The optional weighted
implementation is retained for reproducibility. The next diagnostic should examine
the actual missed waveform packets against their templates: this experiment rejects
one contrast-derived weighting rule, not template improvement or all possible
weighting schemes. It also does not justify loosening the contamination criterion.

The [artifact verification](../results/method_interleaved_site_weighting_v1/verification.json)
recomputed matched-reference and control counts at all seven thresholds, confirmed
equal observed/control exposure, and reproduced all station-to-site weights from
the frozen formula. The full 90-test suite passed after the implementation change.
