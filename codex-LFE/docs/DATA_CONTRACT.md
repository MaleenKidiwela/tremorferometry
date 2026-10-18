# Data and interpretation contract

## Raw and preprocessed waveforms

Full network.station.location.channel IDs; exact UTC start; sampling rate; separate sample mask; source-byte hash; processing-code hash; filter settings; instrument/site identity. Missing or flat-lined samples and conflicting overlaps are invalid. Filtering operates on contiguous valid runs and excludes edges. Data with a within-channel rate or sample-phase change are rejected for epoch splitting.

Native channel matching is permitted only against templates made from the same channel identities and a validated instrument epoch. Native `1` and `2` are never silently interpreted as geographic `N` and `E`. Neural inference requires StationXML orientation and compatible component responses. The present gain-equalization path does not claim full response removal.

Fractional start times are resampled onto a common UTC grid, with a tested interpolation step. The first two detector pilots used linear interpolation; `pilot.py` preserves that setting for configurations without a `fractional_alignment` key. The general loader and new unused-site/neural runs now use a windowed-sinc interpolator, tested across 1–8 Hz and fractional offsets. This avoids the strong high-frequency attenuation that linear interpolation can introduce. A future detector experiment should explicitly set `fractional_alignment: sinc`. Changing alignment changes the cache key and requires a new experiment record.

## Source and time identity

`family_id` in the first pilot is a hash of an exact coordinate group in the original public Bostock table. That table lacks an original family ID, so nearby/distinct groups are not automatically merged or declared to be the original families.

`time_basis = arrival_anchor:CN.PGC:training_stack_peak` means precisely that. These times are not origin times or guaranteed S-onset picks. The peak offset is learned from training waveforms only. Published reference times retain their original unresolved convention until explicitly transformed with the recorded training offset.

Every comparison or temporal grouping checks time-basis compatibility. A later location solution must provide a new, documented transformation to origin UTC with timing uncertainty.

## Candidate versus event

Neural outputs are unassociated phase candidates. Matcher outputs are experimental event candidates with retained competing families. A close score tie has `family_id = null` and `assignment_status = ambiguous`. The current grouping rule is conservative and provisional; it can merge closely spaced real events and requires validation.

A catalog release requires a scientifically validated status, source mapping, reproducible final selection, and uncertainty. Initial files explicitly say `experimental_unverified`. No physical magnitude or individual hypocentre is assigned by this implementation.

The first template-QC implementation alternates training events between two stacks. Nearby events can share waveform context, so these QC halves are not guaranteed statistically independent. The QC threshold is only a training consistency screen. Entire held-out days and unused-site checks carry the generalization evidence. A future stricter QC split should use separated temporal blocks and report families with inadequate independent training coverage.

## Evaluation

The input to evaluation is the final selected set. Each detection and reference can match at most once. Correct-family and time-only matches are distinct. Maximum-cardinality matching is used; matching inside dense ambiguous clusters is not a claim of a unique physical correspondence.

Eligibility is based on valid waveform exposure, not whether a reference event was recovered. Reports include all selected reference groups, eligible counts, excluded counts, and family-level masks. A source with no usable template must be counted as a coverage failure, not quietly advertised as complete.

`control_to_observed` compares identically selected observations and shifted-site controls. In the initial multi-family summary, exposure is summed family-seconds; interpreting the final merged-catalog ratio is approximate if availability mixes differ. Inspect the per-family/block exposure tables. A release workflow must either balance these masks or stratify and combine the estimates explicitly.

## Reproducibility

Each run has immutable parameters, input hashes, package versions, and source hashes. `north_cpu_v1/source_snapshot` additionally preserves the exact original source. Never overwrite a nonempty output directory. Changing the detector after examining test results creates a new development experiment, not a second independent test of the original hypothesis.

Template extraction now supports `template_alignment: sinc` in pilot configuration. This is separate from `fractional_alignment`, which controls daily waveform alignment. The sinc template option uses a 32-sample padded kernel and rejects windows whose padding touches invalid data. Existing configurations retain linear extraction so saved experiments remain reproducible. A synthetic 8 Hz, half-sample-offset example verifies amplitude preservation; improved field recovery has not yet been demonstrated.
