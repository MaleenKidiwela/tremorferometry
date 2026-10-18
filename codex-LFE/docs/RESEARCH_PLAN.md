# Research plan and acceptance decisions

## Objective

Produce a catalog of individual LFEs and their source-family assignments across Cascadia, with quantified reliability, coverage, and sensitivity. Geographic coverage, event reliability, source-family discovery, and temporal completeness are separate deliverables. A full-margin search cannot be reported as full-margin detection sensitivity merely because every station was processed.

This is a new implementation. No detector limit, station-count requirement, magnitude threshold, or negative result in the parent repository is accepted as a physical constraint without an independent test.

## Evidence motivating the approach

* [Lin et al. (2024)](https://seismica.library.mcgill.ca/article/view/1134) demonstrate a long-duration northern Cascadia LFE catalog using neural detection and subsequent association/location.
* [Ducellier and Creager (2022)](https://doi.org/10.1029/2021JB022986) extend southern Cascadia source families over eight years. [Author manuscript](https://arxiv.org/abs/2212.14581).
* [Münchmeyer et al. (2024)](https://seismica.library.mcgill.ca/article/view/1185) demonstrate source discovery and cross-regional transfer, providing a starting point where templates are absent.
* [Shelly et al. (2025)](https://doi.org/10.1029/2025GL116116) demonstrate iterative template improvement at the southern edge of Cascadia.

These motivate a hybrid discovery/template workflow. None proves uniform performance in untested central Cascadia or with our available station epochs.

Additional primary-source checking found published Olympic Peninsula and central Oregon catalogs. See [regional benchmarks](REGIONAL_BENCHMARKS.md). Their machine-readable releases still need acquisition; existence of these benchmarks is already established and changes the regional validation strategy.

## Method-selection prerequisite

The user has prioritized settling the method before scaling. [METHOD_DECISION.md](METHOD_DECISION.md) proposes neural discovery + physical association/location + matched-filter extension, with a bounded comparison against both component approaches. This prerequisite supersedes moving directly from the current small matcher experiments to fleet scans. Architecture recommendation and demonstrated performance remain separate.

## Ordered milestones

1. **Auditable data and scoring.** UTC alignment, native channel identity, masks, orientation metadata, instrument epochs, signed correlation, independent sites, final-row one-to-one matching, and saved provenance. Unit tests include deliberate counterexamples, not just successful detections.
2. **Independent northern demonstration.** Training-only templates from original reference records; entire held-out days; fixed threshold selection; repeated site-preserving nulls; an unused recording site for a positive waveform check. Record failures and exclusions. Existing pretrained Cascadia stacks are not independent evidence for this pilot.
3. **Source discovery pilot.** Run at least two discovery baselines, including a model not trained on Cascadia. Associate physically consistent phases with a documented velocity model, then verify candidates by held-out-site waveforms and repeated sources. A phase pick is not an event. A tremor-correlated event count is not proof of correct associations.
4. **North / central / south regional replication.** Choose fixed short periods covering active and quieter conditions from records independent of candidate performance. Include at least one genuinely unfamiliar central region. Verify data and metadata availability first. Run calibration and held-out evaluation independently in each region.
5. **Improve local templates.** Iterate only using training events that pass independent-site checks. Keep a sealed evaluation period and maintain a separate family-disjoint evaluation for source discovery. Template stacking itself can amplify correlated noise; a prettier stack is not evidence unless unseen data reproduce it.
6. **Historical expansion.** Process blocks with fixed station/channel epochs. Rebuild or validate templates across instrument changes. Compare overlapping station subsets to quantify changes in sensitivity. Keep acquisition gaps distinct from nondetection.
7. **Catalog release.** Publish verified event rows, ambiguity flags, family provenance, coverage/exposure, empirical control results, uncertainty, evaluation exclusions, and rejected candidates. Release boundaries follow validated evidence, not a desired event count.

## Validation rules

* Score the final catalog after temporal grouping and family resolution. Use one-to-one matching. Report any-family temporal recovery and correct-family recovery separately.
* Published catalogs are incomplete. Their match fraction is reference agreement, not event precision; unmatched events require additional evidence.
* Keep station groups that share a site together in every null. Use repeated positive/negative shifts within short blocks to preserve local rate variation. Remove wraparound edges. Apply identical event selection to observations and controls. Report exposure per family/block/realization.
* The initial null experiment estimates a control-count ratio, not an event probability. Cluster dependence, shared LFEs in shifted controls, and differing family availability can bias interpretations. A zero control count is not proof of zero false alarms.
* Use held-out sites, blinded waveform inspection stratified by score/region/noise state, and known conventional-earthquake contamination checks. Reviewers should not know which candidates are shifted controls.
* Injection tests must include independent templates, real noise, range of amplitudes, missing stations, and correlated site artifacts. Injection completeness must not be equated with natural-event completeness without reference validation.
* Train/calibration/test dates never overlap. Thresholds and family-selection rules freeze before test processing. If a test informs a redesign, that test becomes development data; reserve a new final evaluation period.
* Do not force source locations onto the interface to obtain an attractive map. Distinguish unconstrained locations from results using a stated geological prior.

## Provisional quantitative targets

These are engineering decisions to be tested, not values inferred from the old work:

* At least three independently located supporting sites in the first matcher; a larger count alone does not certify an event.
* Repeated shifted-control ratio no greater than 5% at a chosen operating point, with enough detections and control exposure to make the estimate informative. The first pilot conservatively screens on the maximum ratio across 16 realizations; it does not call this a confidence bound.
* Aim for at least 95% verified event precision on a blinded, stratified independent assessment, with an uncertainty interval and region-specific reporting. Do not claim this from catalog matching alone.
* Measure recall by family, amplitude proxy, geometry, and noise state rather than selecting a universal acceptable recall in advance. Assess whether that sensitivity supports the intended science before pooling regional rates.

## CPU work versus GPU work

Data validation, scoring, controls, modest template scans, model smoke tests, and pilot design are possible on CPUs. GPU availability enables broader discovery and template-bank scanning. A GPU does not fix provenance, association ambiguity, or invalid validation. There is no planned bulk waveform download until a concrete regional roster and byte estimate have been written; modest metadata requests are part of ordinary setup.
