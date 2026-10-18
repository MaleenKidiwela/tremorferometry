# Independent audit of Claude's tests

Reviewed 2026-09-11. This is an artifact and source-code audit, with recomputed association/control counts and station membership. It is not a rerun of GPU inference or an independent validation of every detected event. Existing files in `lfe_catalog_repair` were left unchanged. Paths below are relative to the workspace root.

## Assessment

Claude performed substantial, useful tests. The later corrected results justify continuing with published LFE neural picks, physical association, and template extension as a candidate architecture. They do not establish a reliable full-margin catalog, prove a fundamental detection ceiling, or validate performance on the user's station roster. My earlier review was incomplete: I should have reconciled these tests before proposing additional picker experiments.

## What the artifacts establish

| Experiment | Checked evidence | Assessment |
| --- | --- | --- |
| E2 continuous LFEDetect evaluation | `lfe_catalog_repair/scripts/e2_validate_pgc2010.py`; `results/e2_validation_pgc2010.csv` under that directory | Real continuous inference and useful reference-coincidence measurements exist. The scoring is not an independent phase-accuracy benchmark; specific issues below. |
| E3 loose association controls | `lfe_catalog_repair/results/e3_nullhr_rates.csv` | At threshold 0.3 and at least four stations in vi2010, hourly shifted controls produce 36.53% as many associations as observed data. At threshold 0.2 this is 61.10%. This is strong evidence that loose association needs improvement, not a calibrated fraction of false events. |
| E3 tremor contrast | Same CSV, vi2017, threshold 0.3, three stations | Observed ETS/quiet count contrast is 5.444; shifted controls yield 5.395. Activity contrast alone does not validate individual LFEs. |
| E4 corrected injections | `lfe_catalog_repair/results/e4_inject_recall.csv`; `scripts/e4_inject.py`, `scripts/e4_chanorder.py` | Corrected injections already exist. At threshold 0.1, published-stack injections have excess recovery 0.133, 0.900, 0.983 at injected RMS ratios 0.2, 0.5, 1.0 (60 trials per cell). These are controlled stack-in-noise sensitivity measurements, not natural-event completeness. |
| E4 tighter association controls | `lfe_catalog_repair/results/e4_assoc_null_table.csv`, underlying observed and shifted event CSVs, `scripts/e4_assoc_null_score.py` | Recounted underlying files; summary ratios reproduce. A useful candidate operating point exists, with limited reference recovery and incomplete independent validation. See table below. |
| E3 real-pick relocation diagnostic | `lfe_catalog_repair/results/e3_realpick_summary_t0.1.csv` | 123 cases have three stations; only one has four. This does not robustly establish location performance for four-or-more-station events or isolate association as the sole cause of location errors. |
| Older feature/picker work | `lfe_features/RESULTS.md`, R4–R9 | Reports window-classifier AUC 0.975 and per-sample AUC 0.848 with pick precision 0.03. These are different tasks and labels. I reviewed these reported results, not independently reproduced their training/splits. They do not prove a universal single-station limit or that clustering automatically yields clean LFE families. |

### Recomputed E4 association counts

The stored tight configuration uses a 1.2 s matching tolerance, 0.5 s pick standard deviation, and depth search restricted to 20–50 km. These choices changed together; their individual effects have not been isolated. Reduced search volume itself affects chance-association rates.

| Pick threshold | Minimum stations | Observed candidates | Constant-shift mean / observed | Hourly-shift mean / observed |
| --- | ---: | ---: | ---: | ---: |
| 0.3 | 4 | 2,430 | 17.54% (5 seeds) | 19.96% (2 seeds) |
| 0.3 | 5 | 646 | 2.60% (5 seeds) | 3.79% (2 seeds) |
| 0.5 | 4 | 384 | 4.58% (5 seeds) | 5.21% (2 seeds) |
| 0.5 | 5 | 82 | 1.22% (2 seeds) | 2.44% (1 seed) |

For threshold 0.3 / five stations, the underlying constant-shift counts are 16, 19, 15, 13, 21 and hourly counts are 28, 21. The stored Bostock coincidence recall is 1.78% within ±2 s (15,931 reference events), using family timing corrections and nonexclusive temporal matching. This is a low-recovery candidate-seed result. It is not evidence of 96–97% true precision. Very few control realizations were completed at threshold 0.5 / five stations.

The E4 script divides counts by 34 days, while the vi2010 configured interval, August 10 through September 13 inclusive, spans 35 calendar days. The saved picks cover that interval. Use counts and ratios here; establish usable waveform exposure before quoting operational rates. The common divisor cancels in the displayed ratios. The accepted Bostock reference events lie inside the observation interval; I checked this rather than assuming a reference-date mismatch.

## Specific picker-scoring problems

1. **Matching is not one-to-one.** `match_frac` independently finds nearest neighbors in each direction. Running the existing function with references `[0, 1]`, one pick `[0.5]`, and tolerance 1 s returns recall 1.0 and agreement 1.0. One-to-one recall can only be 0.5. This demonstrates a scoring limitation; its aggregate impact on real data has not yet been quantified.
2. **P rows use S references.** Every entry in `gts` is a predicted S-arrival series, including the concatenated `either` series. The subsequent loop scores both P and S against that same dictionary. Consequently the P rows do not validate P timing. This does not invalidate S rows by itself.
3. **Timing is fitted on the evaluation period.** The script estimates a residual-histogram peak for each reference series using Cascadia-model S picks from the full evaluation interval, then scores that interval after applying the offset. The initial Bostock-to-Lin clock correction is separate; the additional pick-derived adjustment still makes this a development comparison rather than an untouched timing test.
4. **Reference agreement is labeled precision.** Published event lists and modeled arrivals are incomplete, uncertain phase truth. An unmatched pick is not automatically a false LFE. Conversely a temporal match, especially within 3–5 s during bursts, does not prove event identity.
5. **Exposure and independence need explicit accounting.** The scorer divides by 184 calendar days and does not impose a per-reference inference-coverage mask. Cascadia-inclusive pretrained weights also require a training-provenance audit before calling a Cascadia stack benchmark independent.

## Corrections already made by Claude that should be retained

- The E3 injection component-order error was identified and E4 reran corrected injections. Do not use the original E3 injection curve as current evidence.
- E4 briefly argued that real-event signal amplitudes were much larger than the injection anchor. E5 withdrew that comparison because incoherent window energy and coherent template amplitude are different quantities (`lfe_catalog_repair/RESULTS.md`, E4 withdrawal and E5 step 0). That withdrawal is appropriate. It does not establish that fine-tuning can never help.
- E4/E5 contain borehole-transfer and training-derived station-weight experiments. Those deserve reuse and targeted validation; starting again with generic injection or model-loading demonstrations would duplicate completed work.

## Transfer to the user's 183-station roster

I compared the `net`/`sta` pairs in `e2_picks_vi2010_cascadia.parquet` with `codex-LFE/results/target_network_from_map/stations.csv`.

- Included: **CN.PGC, CN.SNB, CN.VGZ, PB.B011, PB.B926**.
- Outside the roster: CN.LZB, CN.PFB, CN.YOUB, CN.NLLB, UW.LRIV.

The 646 threshold-0.3/five-station events have at least three, four, or five allowed assigned stations in **483, 219, and 26** cases respectively. Counts come from `e4_assoc_vi2010_cascadia_t0.3_n5_assign.parquet`. This is only an assignment-membership diagnostic: the original association used outside stations, so none of these counts is a valid roster-only catalog result. Both observed and shifted association must be rerun after restricting inputs. Five-station northern performance cannot be extrapolated across the margin.

## Revised next work

September 12 follow-up: [the corrected interleaved experiment](INTERLEAVED_METHOD_RESULTS.md)
has now completed. Source inspection found that E2 native-stream `annotate` calls
used flexible 1→N / 2→E component aliases without geographic rotation. Response
metadata records B011 horizontal azimuths 11°/101° and B926 azimuths 299°/29°.
This is an input-convention mismatch. The new CPU inference path verifies gains
and rotates to ZNE, while also changing filtering/stitching and peak extraction;
it does not isolate the effect of orientation. Original E2/E4 artifacts remain
intact. The low recovery of either legacy or corrected associations is not proof
of a fundamental neural-picker sensitivity limit.

The tasks below describe the original follow-up plan; completed steps and the
remaining event-identity/missed-reference work are recorded in the linked report.

1. Reuse existing picks to reproduce the northern association baseline on the five allowed stations, including both shift modes. Check waveform coverage and pick preprocessing first. Use CPU association before requesting GPU inference. Keep location search assumptions explicit and inspect boundary solutions.
2. Re-score existing reference comparisons with phase-correct arrivals, one-to-one matching, independent clock calibration, and usable-exposure masks. Call these reference-agreement measurements; use independent multistation review for event identity and precise phase accuracy.
3. Inspect candidate waveforms and withheld allowed-site evidence, where geometry permits, including conventional earthquakes and artifacts. Preserve uncertainty and ambiguous events rather than assigning automatic true/false labels.
4. Only then run the missing cross-region, borehole, and held-out-family comparisons. Compare neural association alone with template extension at comparable measured contamination. The existing 24-window review packet is a development tool, not the entire validation requirement or a prerequisite to reusing completed experiments.

No model winner, full-margin completeness, fundamental sensitivity ceiling, or validated final catalog is claimed by this audit.
