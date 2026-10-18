# Method decision before regional expansion

Status: recommendation for discussion and a bounded method comparison; not a validated production method. The user has explicitly prioritized settling the method before scaling catalog production. Further fleet scans and incremental tuning of the current matched-filter pilot are deferred until this decision is supported.

The [Claude test audit](CLAUDE_TEST_AUDIT.md) reconciles existing picker, injection, and association experiments. Reuse the saved northern picks and corrected tests first. The immediate missing comparison is association restricted to the five available roster stations, with matching controls and corrected scoring; additional generic picker demonstrations would duplicate existing work.

The [bounded comparison specification](METHOD_COMPARISON.md) defines the current development experiment and evidence needed to rank the methods. The user has reiterated that method selection remains the objective; bulk downloads and margin-wide inference are deferred.

## Recommended architecture

Use **LFE-specific neural discovery, physical multistation association/location, and template matching to extend verified repeating sources**. All waveform evidence must come from the 183 stations in the designated map. Use overlapping regional subsets within that roster; requiring every station across the margin to observe an event is inappropriate.

1. Verify station/channel epochs and preprocess consistently, retaining gaps and response/orientation provenance.
2. Run a published LFE-specific phase detector to propose arrivals. Compare a model trained without Cascadia against a Cascadia-inclusive model. Do not interpret output scores as event probabilities without calibration.
3. Associate arrivals using predicted P/S travel times in a documented regional velocity model. Search source position and origin time; retain uncertainty and competing associations. Test sensitivity to velocity-model errors and noise. Do not force sources onto the slab to make their locations look plausible.
4. Verify associated candidates with waveform evidence at withheld roster stations where geometry permits, conventional-earthquake checks, and repeated sources. Waveform similarity and compatible travel times define families; map certification labels and rounded coordinates do not.
5. Build station-specific templates from verified training events. Search continuous records using multistation matched filtering to recover weaker repeats. Freeze template updates before independent evaluation. Preserve validated events without a reliable family rather than discarding all nonrepeating events.
6. Resolve duplicate events across regional overlaps and competing templates. Keep origin time distinct from station arrival anchors. Record family provenance, location uncertainty, quality evidence, and station/epoch exposure.

## Why this architecture

Pure known-family matching is constrained by its template inventory. LFE-specific neural detection can propose previously unknown sources, but association and validation are still required. The proposed hybrid assigns discovery and weak-repeat recovery distinct roles rather than treating either detector score as proof.

Published support: [Münchmeyer et al. 2024](https://seismica.library.mcgill.ca/article/view/1185) demonstrate transferable LFE phase detection, association, and location, including previously uncataloged sources. [Lin et al. 2024](https://seismica.library.mcgill.ca/article/view/1134) demonstrate neural detection with regional 3D velocity-model location in southern Vancouver Island. [Ducellier and Creager](https://arxiv.org/abs/2212.14581) provide a southern Cascadia template-based catalog and an independently inspectable implementation. These support components; they do not establish this combined workflow's performance on our permanent-station roster.

The present three-station calibration is a baseline, not evidence that matched filtering is intrinsically limited to its 3.55% recovery. Our template construction, family selection, timing, and station subset are limited. Neither negative station-removal results nor interpolation changes settle the architecture.

## Comparison needed to settle the method

Use the same modest waveform panel, roster station subsets, time splits, and final-event scoring for three arms:

- Published-family matched filtering with appropriately verified templates.
- Neural phase detection plus physical association/location.
- The same neural discovery followed by training-only family construction and matched-filter extension.

Include active and quieter intervals, a known-source case and a withheld-source-family case, and more than one region to test transfer. A known-family arm cannot be credited with discovering a withheld family merely because it reports a nearby time. Neither published references nor low shifted-control counts alone measure precision.

Choose the operating points on development data, then compare final-event precision assessed with independent waveform review, known-source recall, validated new-family recovery, duplicate/association errors, location uncertainty, and compute/data cost on untouched dates. Repeat shifts through the full association or matching chain. Apply the same independent review standard to all arms and include ambiguous/unreviewable events in reporting. Add real-noise injections as a sensitivity diagnostic, not a substitute for natural-event validation.

Prefer the simplest arm meeting the scientific requirements. Adopt the hybrid only if its extension adds verified recovery at comparable contamination; neural association alone remains a legitimate outcome if matching adds no useful gain. If neural discovery misses coherent sources, a bounded cross-station-correlation discovery comparison is the next hypothesis, not an automatic bulk fallback.

## What is not settled yet

The specific neural weights, association algorithm and velocity model, family clustering rule, template-building procedure, operating thresholds, and evidence-based station requirements must be selected by this comparison. A model smoke test and software tests do not settle them. Current code includes inference and matching, but does not yet implement a validated full association/location and family-discovery workflow.

Next deliverable: a concrete method-comparison specification and a verified association/location baseline. No production or full-margin completeness claim follows from selecting an architecture on paper.
