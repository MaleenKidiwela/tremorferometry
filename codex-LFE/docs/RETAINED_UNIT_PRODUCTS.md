# Retained unit products

`codex_lfe.catalog_products.verify_unit` checks experimental output integrity before owned waveform cleanup. `MonthlyStorage.commit_verified_unit` connects that check to atomic ledger commitment. This is an integration path, not a validated scientific production method.

Each unit retains five inputs plus the generated verification JSON:

| Product | Format and required content |
|---|---|
| events | Parquet: unique event_id, time_s, time_basis, nullable family_id, quality_status, supporting_stations, support_windows |
| families | Parquet: unique non-null family_id; additional source coordinates and provenance may be retained |
| coverage | Parquet: station, start_s, end_s, status, evidence_id |
| controls | Parquet: unique positive integer realization, nonnegative integer count, nonnegative exposure_s |
| provenance | JSON: time_basis and coverage_evidence keyed by evidence_id, each with description and source_sha256 |
| verification | Generated JSON with all product hashes; explicitly marks scientific validation as not performed |

Event times must lie inside the half-open unit core. Every station in the explicitly supplied roster must have a complete partition of coverage intervals, without overlaps or silent gaps. The production adapter must pass the exact 183-station roster; smaller scopes are permissible only for explicitly labeled integration tests. Coverage statuses are processed, no_metadata, provider_no_data, waveform_gap, or insufficient_context. Request and processing failures cannot pass verification.

Supporting stations are comma-separated full network.station identifiers. `support_windows` is a JSON list of objects with station, start_s, and end_s. These are the actual waveform windows used by the detector, including moveout; the verifier does not substitute an event's anchor time for its station arrival windows. Every support window must be fully covered by processed intervals.

Optional `support_coverage.parquet` permits actual retained waveform support before or after the event core. It has the same columns as core coverage, uses only `processed` intervals, and requires explicit `support_coverage_bounds` (`start_s`, `end_s`) in provenance. Every padding interval must lie outside the core and inside those bounds, have evidence, and use a declared roster station. Overlapping padding intervals are rejected. Padding cannot override core gaps or failed coverage. Event timestamps remain in the half-open core; only their waveform support may extend beyond it. The extra Parquet is included in verification and ledger hashes before cleanup. Without this optional product, the original core-only rule applies.

Accepted quality statuses are experimental_unverified and reviewed_candidate. Neither the integrity verifier nor the monthly ledger certifies a final LFE catalog. Coverage evidence checksums preserve traceability but the verifier does not independently rerun waveform acquisition, infer missing metadata, or authenticate a provider receipt.

A unit with no waveforms is registered through `begin_unit`, then retains empty event/control tables with their schemas, a family table, and explicit coverage/provenance. It still passes through verification, commitment, and cleanup before the month can finish. No dummy waveform is needed. Only ledger-owned staged files can be deleted; retained products and original cached data remain intact.

Validation covers missing roster entries, overlaps, gaps, failed requests, absent evidence, unknown families, unit-end event boundaries, gaps inside support windows, and an integrated commit/restart/cleanup cycle that preserves the original waveform file.

Padding tests cover continuous support across both core boundaries, missing padding support, attempts to override core coverage, undeclared stations, overlaps, absent evidence/bounds and out-of-bounds intervals. These checks validate retained support structure; they do not independently authenticate waveform data or certify detections.
