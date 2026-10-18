# Full-margin catalog execution plan

The requested endpoint is a finalized Cascadia LFE event-and-family catalog preferably from January 1, 2011 through August 31, 2026 (exclusive cutoff September 1, 2026 UTC), with January 1, 2019 as the fallback start if regional network capability requires it, on the exact 183-station roster. Every new artifact stays in `codex-LFE`. Intermediate diagnostic results do not complete the task.

## 1. Reuse existing evidence before new tests

Completed: prior picker/matcher audits, one-to-one rescoring of saved neural picks, published-stack coverage repair, southern transfer diagnostics and source-array checks. Existing northern training/calibration and previous archive caches are being reused. See `PICKER_STATUS_20260914.md` and `results/prior_picker_test_audit_v2`.

## 2. Validate regional detection and phase timing

In progress. Neural phase picking is not approved for full-margin production. The southern matched-filter panel now has a completed development day and fixed next-day test with an unused station, using seven published templates representing five families. Controls and published-reference recovery show useful detections, but independent-site confirmation is uneven. There are 23 next-day candidate rows for review, not 23 certified unique LFEs. Complete waveform overlays and preserved snippets are available in `results/ducellier_candidate_waveform_review_v1`.

The expanded northern bank has 47 usable training groups and a calibration-selected 9-MAD threshold. Its declared September 18–20, 2012 evaluation is complete: 30 candidates, 18 same-family matches, and a worst shifted-control ratio of 6.67%, above the frozen 5% criterion. PB.B926's date-separated independent template check supports 15 of 47 groups, containing only two evaluation candidates. This does not validate the whole bank. Central Oregon/Olympic LFE phase-reference acquisition remains incomplete. A PNSN-activity-selected September 2, 2026 central panel and September 13 quiet panel have now been processed with the same eight roster stations. Their associations do not exceed station-shift controls, so these seeds are not promoted to templates. Public sources and failed access attempts are documented separately.

Proceed through independent waveform checks, conventional-earthquake discrimination, family ambiguity and defensible timestamps. Preserve failures and uncertainty. Do not lower thresholds or redefine truth to make a test pass.

## 3. Establish source coverage across the margin

Known northern and southern template banks do not establish coverage of central Cascadia. Validate discovery on actual roster stations in those gaps, reusing public examples where obtainable. Distinguish a lack of detections from low sensitivity, missing recordings and an untested source region. Family locations borrowed from publications must be labeled as such; they are not new individual-event hypocenters.

## 4. Integrate bounded monthly processing

The current scoped descending queue contains 188 months and 21,490 station-epoch work blocks, from August 2026 back through January 2011. The January 2019 fallback contains 92 months. See `results/network_temporal_scope_v1/descending_month_queue.csv` and `station_month_work_plan.csv`. These supersede the temporal scope of the preserved older all-history plan. They are metadata plans, not proof of data availability. Start with August 2026, then move backward. Keep one active calendar month, subdivided into network-coherent daily/subdaily units, with required boundary padding. The provisional temporary waveform ceiling is 100 GB.

`codex_lfe.monthly_storage` implements staging ownership, a budget, one active month, backward ordering, atomic retained-output commitment and checksum-verified cleanup. It is not yet connected to a production detector. Integrate it only with validated configurations and explicit coverage/failure records. Retain event/family tables, review snippets, templates, controls/sensitivity, response epochs, request receipts, hashes and restart state. Delete only the production job's temporary staged raw files after output verification. Existing shared/training/evaluation caches remain preserved.

The real-data integration replay is complete in `results/monthly_panel_replay_v1`: two existing August hours, processed in descending order, each with eight stations. It verified 356 retained files, restarted from the ledger, and deleted 18,012,160 bytes of newly staged copies while preserving every original waveform. The ledger's completed month refers only to that isolated two-unit diagnostic plan. It does not mark August complete in the production queue. Planner, storage and product-verification checks pass (21 tests).

## 5. Run the full usable history and verify the release

Carry source-family state across month boundaries. Revalidate transfer when stations or instruments change. Assign events to half-open core intervals to avoid duplicate boundary counts. Keep unavailable, failed, unsupported and processed intervals distinct. Resume from verified checkpoints rather than starting over.

Finalize unique event/family identifiers, timing conventions and uncertainties, geographic/source coverage, detection sensitivity, unresolved cases, code/configuration/model provenance and a versioned release. Keep the requested end fixed at August 31, 2026; later validation examples are not automatically part of the catalog. A northern subset, external reference catalog, successful pilot or running job is not the requested finished product.

Current machine-readable checkpoint: `results/ACTIVE_STATUS.json` (a timestamped snapshot, not a live monitor). Detailed findings and active jobs: `docs/ACTIVE_WORK.md` and the corresponding result logs.

## Current scientific decision

Both new neural adaptation experiments are complete and excluded from production. Synthetic stack/noise gains did not transfer to natural reference windows. Subsequent training on real September 2012 waveforms, with separate training/selection/test dates, matched 1/62 P proxies and 0/208 S proxies on the evaluation dates, with 94/88 unmatched scored peaks. These incomplete proxies are not true precision labels, but no useful field recovery was demonstrated. Original weights and frozen matched-filter evaluations remain unchanged.

The northern August continuous network-score experiment also failed to separate the active hour from time-shift controls. Further threshold tuning on these same evaluated panels would not provide independent validation.

The southern cross-family template audit found no strong aliases among 18 cross-family pairs. A subsequent PB.B039 withheld-station stack check gives modest cohort-level support for family 029 and weak support for the other tested families. It does not certify the 23 individual candidate rows. Results and control distributions are retained without selective promotion.

The user confirms there are no additional reviewed P/S labels or validation results beyond the located tests, and no local copy of `LFEsAll.mat`. Both former input requests are closed. Import support and unsent reference requests remain available; no author contact is authorized or has occurred. Continue validation and network assessment from existing evidence. A missing reference is an unresolved scientific limitation, not a pending user task. No finalized full-margin catalog exists yet.

## Historical network assessment

The completed metadata assessment in `results/network_temporal_scope_v1` integrates simultaneous selected three-component epochs month by month, grouping colocated roster stations within 1 km. January 2011 has 54 physical site groups (60 station codes); January 2019 averages 86.57 sites (95.57 station codes); August 2026 has 159 sites (171 codes). Much of the growth occurred after 2019. These are potential metadata coverage counts, not waveform uptime or validated source sensitivity. Regional bins describe station locations, not source detection coverage.

Preserve 2011 as the preferred start while testing regional/epoch capability; do not discard early northern or southern coverage solely to impose a common date. A final full-margin interpretation requires measured regional sensitivity and explicit unsupported intervals. Existing 2012 northern templates also encounter changed channels and station losses later, so denser later networks still require template/response transfer validation. Bounded 30-second waveform snapshots at August 2026, January 2019, and January 2011 check archival accessibility, one month at a time.

Four additional northern test dates completed with 80 family hypotheses, 75 same-family reference matches among 1,911 eligible references (3.92% recall), and 0–1 hypotheses per shifted-control realization. This passes the descriptive control comparison but does not establish useful completeness or certify unique physical events. The original frozen evaluation failure remains part of the record.
