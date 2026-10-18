# Roster-restricted neural association: development results

Historical legacy-input experiment. The [September 12 corrected comparison](INTERLEAVED_METHOD_RESULTS.md)
supersedes this result for method assessment: the saved borehole neural inputs
had not applied their recorded horizontal orientations. These outputs remain
useful for reproducing the original association and control diagnostics.

Completed 2026-09-11 using existing Cascadia-model picks and CPU association.
No new waveform downloads or GPU inference. This is evidence for method
selection, not a validated catalog or a selection of the final method.

## Result after accounting for nearby instruments

Five allowed stations represent four proximity groups: CN.PGC and PB.B011 are
290 m apart; CN.SNB, CN.VGZ, and PB.B926 form the other groups. Controls shift
PGC/B011 together. Every retained event requires S picks at all four groups.
The period is August 10–September 13, 2010, with the outer 3,120 seconds excluded
from reported events in observed and control runs. This is existing development
data, not an independent test period.

| Pick threshold | Minimum station IDs | Candidates | Constant-site control counts | Hourly-site control counts |
| --- | ---: | ---: | --- | --- |
| 0.3 | 4 | **109** | 5, 7, 7 | 6, 8, 6 |
| 0.3 | 5 | **31** | 1, 0, 1 | 0, 1, 2 |
| 0.5 | 4 | **14** | 1, 0, 0 | 0, 0, 0 |
| 0.5 | 5 | **2** | 0, 0, 0 | 0, 0, 0 |

The mean control/observed ratios at threshold 0.3 are **5.81–6.12%** for four
station IDs and **2.15–3.23%** for five. These are empirical control ratios, not
measured false-discovery rates. Zero controls in three realizations do not
demonstrate perfect precision.

At threshold 0.3 / four station IDs, **45 of 109** candidate depths are within
1 km of the imposed 20–50 km depth bounds. At five station IDs it is **12 of 31**.
This flags location/model/geometry problems; neither low residuals nor a
depth-constrained map verifies LFE identity. The unreviewed candidates are
worth investigating but cannot yet establish a reliable catalog.

## Why the first run is superseded for decision purposes

The archived v1 run required distinct station IDs and shifted every instrument
independently. At threshold 0.3 / four stations it yielded 400 candidates and
control means 2.25–2.83% of observed counts. Requiring four physical groups leaves
109 candidates; shifting nearby sensors together gives the more appropriate
5.81–6.12% control ratios above. The two changes are documented rather than
presented as an improvement in the detector. Five station IDs still represent
four physical groups.

The grouped control addresses one concrete noise-dependence problem. A 1 km
grouping cutoff does not prove more distant sites are independent. Legacy picks
lack a complete inference-coverage mask, and the input header audit found gaps
and overlaps. Operational exposure and false-pick rates remain unresolved.

## Implication for the method decision

**Continue evaluating neural association as a source of candidate seeds. Do not
select it or the hybrid on this result alone.** Increasing thresholds or requiring
all five instruments reduces both control detections and candidate yield. We
have not yet established true-event recovery or location quality at those settings.

Our existing matched-filter results use different dates and subsets. Comparing
their headline recall to these counts would be misleading. The next decision
experiment is the shared-data comparison in [METHOD_COMPARISON.md](METHOD_COMPARISON.md),
with timing reconciliation, candidate waveform evidence, and location sensitivity
addressed before claiming superiority. A hybrid must demonstrate additional
verified recovery over its neural seeds to justify its extra complexity.

## Artifacts and checks

- [Grouped run and control counts](../results/map_network_association_sites_development_v2/comparison.json).
- [Run manifest](../results/map_network_association_sites_development_v2/manifest.json),
  including source hashes, parameters, environment, and preliminary status.
- [Verification](../results/map_network_association_sites_development_v2/verification.json):
  28 completed runs; source-pick identity, final counts, and four-group S support checked.
- [Local waveform header inventory](../results/map_network_association_input_audit/waveform_headers.csv):
  171 files present, four missing station-days, 89 gap/overlap records. Header
  availability is not inference coverage or signal quality.
- Software suite: **75 tests passed**. This checks software behavior, not event identity.

The runner and invocation are recorded in the comparison specification. Original
Claude outputs were read without modification. All new outputs are in `codex-LFE`.
