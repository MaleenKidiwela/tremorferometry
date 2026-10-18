# Northern candidate catalog v1 — September 14, 2026

The first consolidated catalog increment contains **61 experimental candidates across 17 scanned days**, spanning August 20–September 11, 2010, with detections in nine of eleven usable source coordinate groups. It is not a continuous-period or margin-wide catalog. Training dates and unscanned dates are not represented as zero-event days.

[Browse the catalog](../results/north_candidate_catalog_v1/index.html) · [CSV](../results/north_candidate_catalog_v1/candidate_catalog.csv) · [Parquet](../results/north_candidate_catalog_v1/candidate_catalog.parquet)

## Frozen method and held-out evidence

The four-second packet experiment recovered 24/336 calibration references, versus 31/336 for the eight-second baseline. Both selected 9 MAD. The baseline was retained, then frozen before the new scans.

| Date stratum | Candidates | Same-group matches / eligible references | Maximum control/observed ratio |
| --- | ---: | ---: | ---: |
| Development calibration, 4 days | 34 | 31 / 336 | 2.94% |
| Previously examined evaluation, 6 days | 15 | 13 / 342 | 6.67% approximately |
| New quiet block, September 6–10 | 0 | 0 / 5 | Undefined |
| New active-episode days, August 20 and 24 | 12 | 11 / 121 | 0 observed across 16 controls |

The quiet block had substantial usable waveform exposure but only five eligible references. After seeing that lack of reference support, the two remaining unused cached days within the original active episode were scanned with the same frozen method. This date-selection sequence is recorded; these are project-held-out dates, not a claim that no other project examined them.

The new active block recovers 9.1% of eligible references. The small counts and zero shifted controls do not establish event precision or a zero false-alarm rate. The older evaluation stratum exceeds the 5% engineering criterion. Neither calibration nor previously examined evaluation rows are presented as independent validation of the final method.

## Contents and verification

Every candidate carries its arrival-anchor time, source-group assignment, detector score, support, development/held-out stratum, and one-to-one matched reference when available. Source-group coordinates are inherited published group locations, not newly estimated individual hypocentres. Times are PGC training-stack arrival anchors, not validated origin times. No magnitude is fabricated.

Coverage intervals, all shifted-control candidates, per-realization exposure, daily counts, source groups, and per-family/reference metrics accompany the table. The file records 55 one-to-one reference matches in total; this pooled number is bookkeeping, not independent precision validation.

Artifact checks verify unique candidate IDs, exactly the frozen threshold, observation-only rows, UTC round trips, one-to-one reference matching, and that every candidate lies within its source group's usable coverage. Software suite: **90 passed**.

[Verification](../results/north_candidate_catalog_v1/verification.json) · [Validation by stratum](../results/north_candidate_catalog_v1/validation.json) · [Provenance](../results/north_candidate_catalog_v1/manifest.json)

## Remaining catalog work

Independent event confirmation, origin-time solutions and wider source discovery remain required. No additional roster station has cached waveforms covering all ten northern training/validation days, so a new unused-site check requires waveform acquisition. Northern sensitivity is low and cannot be represented as completeness.

Regional cache inventory finds 13 northern, two Washington, one Oregon and zero southern-edge roster stations with existing day files (planning latitude bins, not sensitivity regions). The southern published release is now normalized separately: 61,441 external reference records, 27 original named families, with the documented template-start clock retained. These references are not added to our candidate catalog.

A first southern station/date/volume plan is saved in `results/catalog_regional_preflight_v1`. It requires provider availability and physical-site checks before acquisition. Central/Olympic reference acquisition and regional discovery remain outstanding. See [the build plan](CATALOG_BUILD_PLAN.md).
