# Prior picker work: reuse audit

All new artifacts are confined to codex-LFE. Older projects were read without modification. This is a source/artifact audit and saved-pick re-score, not a new inference run.

## Work already completed

- E2 continuous LFEDetect inference and reference comparisons: reuse saved picks; do not repeat inference to correct scoring.
- E4 corrected noise injections: retain these instead of the superseded component-swapped E3 curves.
- E4 association controls and borehole/template transfer: reuse with original station, geometry and clock caveats.
- E6 physical-site controls, site merging, withheld VGZ/SNB checks and template variants: completed GPU outputs exist even though the older summary omitted them.
- codex-LFE corrected interleaved neural/matched-filter comparison: already performed, not a missing first test.
- lfe_features R4–R11: existing classifier and custom-picker tests; different tasks and approximate arrival labels.

## New calculation using saved predictions

Reproduced all 56 selected legacy S-reference recall rows exactly, then applied existing maximum-cardinality one-to-one matching. Retained the original fitted time offsets to isolate the matching effect. No P predictions were compared with S references.

| Model | Reference subset | Old recovery | One-to-one recovery |
|---|---|---:|---:|
| cascadia | bostock_stack | 6.11% | 4.19% |
| cascadia | bostock_near40 | 8.45% | 5.44% |
| jcms | bostock_stack | 5.36% | 3.69% |
| jcms | bostock_near40 | 7.07% | 4.54% |

Table uses model score 0.3 and ±3 seconds. These are reference-coincidence fractions with legacy evaluation-fitted clocks, not independent phase recall or true precision. Original inputs remain unchanged. Across all 56 audited settings, the largest absolute recall inflation is 10.90 percentage points.

## Work still missing

A full-margin picker decision needs phase-specific independent labels, usable inference exposure, representative station/instrument coverage and independent event checks. The existing 24-window review packet has zero completed human labels. Pending windows cannot be negatives. Existing corrected injections constrain controlled sensitivity but do not establish natural-event completeness.

The attempted jms published-stack diagnostic failed on a nonbucket HDF5 identifier before writing scientific scores; it has not been restarted. Its design and failure record are retained. No duplicate generic injection, denoiser, model-loading or station-merging experiment is authorized by this audit alone; identify a concrete unresolved question before scheduling new inference.

Artifacts: reuse_ledger.json, matching_audit.json, saved_S_picker_matching_comparison.csv, workspace_boundary.json.
