# Picker and detection validation status

Full-history production has not started. No configuration is validated for the full margin. The target remains all usable history through the present on the exact 183-station roster, processed backward with one active month.

## Existing work reused

The audit reused saved continuous picks, corrected injections, association controls and completed E6 station-withholding experiments. Rescoring reproduced 56 legacy S-reference rows before enforcing one-to-one matching. At score 0.3 / ±3 seconds, Cascadia-model full-stack reference recovery changed from 6.11% to 4.19%. This comparison retains the legacy fitted clock and does not measure true phase precision. See `results/prior_picker_test_audit_v2/README.md`.

The jms stack diagnostic now covers all 1,067 nominal P labels and 1,453 nominal S labels. At score 0.3 / ±0.5 seconds, recovery is 53.3% for P and 45.6% for S; at ±1 second it is 65.3% and 59.0%. Row splits share families, all P labels are aligned at sample 800, and stacked traces are not individual-event recordings. A systematic early residual remains. The tail-only repair reused covered results and evaluated only previously uncovered labels. Production inference was not changed. See `results/picker_labeled_stacks_tail_v1`.

## Southern benchmark correction

The Shelly 2025 detection array has no station overlap with the requested roster. Our BK.RVIT and PB.B045/B935 test stations were used in later published stacking/location stages. Consequently, that experiment tested transfer to a different array, not a replication of the published detector. See `results/prior_picker_test_audit_v2/shelly2025_station_method_audit.json` and the [published supporting methods](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029/2025GL116116&file=2025GL116116-sup-0001-Supporting+Information+SI-S01.pdf).

The January date-separated template consistency check is complete. Neither four-site family survives the original correlation threshold under this stricter split; one family retains three sites. This weakens the earlier training-QC evidence. It does not prove the individual templates are false. No production templates were replaced. BK.RVIT component imbalance/correlation remains an unresolved descriptive finding.

The alternative Ducellier–Creager benchmark uses actual roster stations and author-provided templates built from 2007–2009 waveforms. Small 2010 reference/control packets show coherent template matches but few neural picks at score 0.3. Full response removal did not rescue neural sensitivity in that small panel. Template timestamps are matching-window start times, not independently labeled phase arrivals or asserted physical origins. Source data and checksums are in `results/southern_benchmark_inventory_v1` and `results/ducellier_roster_packet_check_v1`. Primary source: [published code and data](https://github.com/ArianeDucellier/catalog), pinned commit `94717ad64bf808808ecedaeec4d92bb5fdfcf044`.

## Continuous matched-filter development

Published 35-second templates, response-corrected velocity, verified component orientations, 1.5–9 Hz filtering, signed equal-channel correlations and independent-site support are now tested on continuous recordings. Eight nonwrapping time-shift controls share the observed exposure mask. Thresholds 6, 8, 10 and 12 MAD were specified before scanning.

| Interval (UTC) | Usable sites | Observed candidates at 8 MAD | Published matches within 0.5 s | Eligible published records |
|---|---:|---:|---:|---:|
| 2010-02-17 18:00–19:00 | 5 | 7 | 7 | 8 |
| 2010-02-18 18:00–19:00 | 5 | 0 | 0 | 0 |
| 2010-04-25 04:00–05:00 | 5 | 10 | 9 | 10 |

These are development diagnostics, not measured precision. The February active hour has 16 total published records, but control boundary exclusions leave eight in common exposure. Each interval provides 2,415 seconds of shared exposure. Missing/invalid NC.KRMB and NC.KHBB packets are recorded separately from absence of events. The following-day hour was scanned with unchanged settings; its empty reference catalog is not a human-reviewed negative label.

The full November 22, 2010 day and seven-template extension are complete, including response-taper validity correction and a fixed November 23 test. Seven stations supplied 144 MB; the family extension reuses those raw recordings. The extra published templates/reference tables total 2.3 MB. See the corresponding `ducellier_continuous_*` and `ducellier_multifamily_*` result directories and logs for completion state. Family-specific template clocks must not be merged as physical origins.

## Remaining gates

Independent phase/event review is still missing. The southern blinded packet has 29 windows, complete prediction coverage, and no submitted labels; pending windows must not be scored as negatives. Additional existing reviewed-label paths were requested asynchronously.

Needed before a finalized catalog: regional detection/discovery validation, independent event checks, defensible timing/location and family linkage, actual historical availability and sensitivity assessment, monthly pipeline integration, full-history scans and a verified versioned release. The storage lifecycle has five passing ownership/budget/restart/cleanup tests, but it is not yet a production runner.

All new artifacts stay inside `codex-LFE`. Prior project caches remain read-only inputs and are not disposable production staging.

## Subsequent completed tests

Two bounded jms adaptations are excluded from production. The synthetic source-held-out test improved stack/noise S F1 from 0.405 to 0.588, but its natural source-held-out check matched zero P or S proxies at score 0.3. A separate real-waveform adaptation used 3,557 training, 2,892 selection and 1,357 evaluation windows on disjoint dates. The selected model matched 1/62 P and 0/208 S proxies, with 94/88 unmatched peaks in scored intervals. Reference coverage is incomplete, so these are proxy-agreement statistics, not true precision/recall. See `results/north_natural_picker_adaptation_v1`.

The continuous network S-score scan produced six active-hour proposals versus 3-21 per shift at 4 MAD, and none at 6/8 MAD. Its timing alignment test passed; it did not demonstrate discrimination. The original and adapted pickers remain unvalidated for full-margin production.

The storage path now has a completed real-data replay of two eight-station August panels. It verified 356 retained files and deleted only 18 MB of staged copies, preserving original evaluation waveforms. This is infrastructure validation, not a production archive scan. The current full execution plan is `CATALOG_BUILD_PLAN.md`.


## Latest central and northern decisions

The expanded northern bank's frozen September 18–20, 2012 evaluation produced 30 candidates and 18 same-family matches. The worst shifted-control ratio was 6.67%, above its fixed 5% criterion. B926 independently supports 15 training groups containing only two evaluation candidates. That is limited independent evidence, not approval for the full bank.

Central September and stronger March 2026 panels do not establish useful discrimination from controls. A fixed 15-second inference-grid shift retained only 18/41 active September S peaks within 0.5 seconds at score 0.3; library overlap averaging retained 4/10. The strongest March cross-station waveform is consistent with a distant earthquake P arrival. Neither phase-only association nor the tested cross-station statistic is approved for production. See `results/central_validation_summary_v1/index.html` and `results/central_strong_conventional_event_audit_v1`.

The next reference-data dependency is concrete: Washington `LFEsAll.mat` has been identified in published analysis scripts, but its underlying data are distributed on request and no local copy was found in the workspace. The MAT importer is ready. Reviewed phase/event references and a validated central method remain unresolved; the full-history scan is not running.
