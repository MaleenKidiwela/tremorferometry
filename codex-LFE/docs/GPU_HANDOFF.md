# Required network constraint

Use `results/target_network_from_map/stations.csv`: the user's designated current map contains 183 stations. The original estimate was 198; exact map membership is authoritative. The next northern configuration is `configs/map_network_north_calibration.yaml`, which enforces membership and excludes NLLB/LRIV. Use `--calibration-only` for development. Channel epochs and waveform availability still need verification. Previous four-station pilot scores are historical diagnostics, not validation of this three-station subset.

# GPU pod handoff

Working directory: `/home/jovyan/tremorferometry/codex-LFE`.

## Before launching work

1. Read `RESULTS.md`, `README.md`, `docs/RESEARCH_PLAN.md`, and the latest CPU experiment's `report.json` and `selection_before_test.json`.
2. Run `python -m codex_lfe.doctor --config configs/cpu_pilot_broader.yaml` and `python -m pytest tests -q`.
3. Confirm CUDA using `python -c 'import torch; print(torch.cuda.is_available()); print(torch.cuda.get_device_name(0))'`. Do not silently substitute CPU for a requested CUDA run.
4. Run the CUDA parity check: `python -m codex_lfe.doctor --cuda-parity`. It compares a known injected signal, masks, and correlation values against the CPU reference.

## Available GPU entry points

* Matcher: `python -m codex_lfe.pilot --config configs/cpu_pilot_broader.yaml --backend torch --device cuda --output results/north_gpu_parity`. This is a reproducibility/performance check, not fresh held-out scientific evidence after the CPU results were viewed.
* Neural picker: the README command with `--device cuda` and an empty new output directory. First verify a CPU/GPU probability comparison on a short identical window. The picker emits unassociated picks; do not count them as LFEs.

The existing matcher deliberately uses high-precision normalization. It is a correctness reference and has not been throughput-optimized. Profile before choosing fleet concurrency or predicting wall time.

## Next scientific work

1. The LRIV unused-site check is complete but inconclusive: none of the candidate groups passed LRIV template QC. Find a suitable independent station using training data, and investigate northern misses without revising the frozen test scores. Reserve new dates for any tuned detector.
2. Prepare explicit central and southern station-epoch manifests. Choose regional evaluation dates before inspecting candidate outputs. Existing local raw data can be reused if byte provenance and timing pass the independent loader.
3. Compare `jms` (no Cascadia in weight metadata) with `jcms` for discovery. Validate association using site-shift controls through the same associator and locator, plus unused-site signal checks.
4. Develop candidate family templates from the discovery training period. Freeze them and their operating thresholds before the regional test periods.
5. Only then expand acquisition and the full-record scan.

## Unfinished capabilities

The project does not yet contain validated template-free association/location, fleet download orchestration, source-family discovery, absolute magnitude calibration, or a release-quality full-margin catalog. These are explicit remaining tasks, not assumed features of a GPU run.

The current suite has 57 passing CPU tests, including torch/NumPy parity and fractional-sample resampling. Existing pilot configurations default explicitly to legacy linear alignment to reproduce the saved experiments; the generic waveform loader now defaults to band-limited sinc alignment. Set `fractional_alignment: sinc` explicitly for a new pilot and treat it as a new experiment. CUDA remains untested. The interrupted `north_cpu_v2` directory is an audit artifact; use the completed `north_cpu_v2_parallel`.

For a new band-limited pilot, explicitly set both `fractional_alignment: sinc` and `template_alignment: sinc`. These change preprocessing and template construction and require a new calibration; they are not a transparent replacement for the frozen v2 run.

## Current target-network baseline and acquisition evidence

The roster-constrained northern calibration is complete at `results/map_network_north_calibration_v1`: threshold 9, 27 candidates, 23 correct-group reference matches out of 648 eligible references, zero detections in 16 controls. No new test evaluation has occurred. Treat all used calibration dates as development data.

`results/target_network_from_map/channel_triplet_epochs.csv` combines independently downloaded metadata for the 183 exact map stations. The station-year audit shows 63 stations with qualifying epochs in 2012 versus 177 in 2024 (172 with >=90% metadata coverage). Prefer designing a modern-period north/central/south pilot, with dates selected from independent activity records and fresh waveform verification, before historical expansion. Do not assume metadata presence establishes waveform availability, usable signal, or complete source coverage.

Metadata requests are reproducible with `python -m codex_lfe.fleet_metadata --roster results/target_network_from_map/stations.csv --start 2009-01-01 --end 2026-09-11 --output results/metadata_new`. It queries only target stations and never requests waveforms. Review existing successful metadata before repeating downloads.

## Picker validation comes before scaling

Read `docs/PICKER_VALIDATION.md`. Two frozen model variants have CPU inference results for the same 24-window northern development review packet. No independent labels or picker performance estimates exist yet. The scorer requires both completed review and valid prediction coverage. The current suite has 70 passing tests. Use `codex_lfe.review_inference` with `--device cuda` to reproduce inference and check CPU/GPU parity after the method-comparison design is frozen; do not treat inference counts as accuracy. Central/southern and borehole panels and event association validation remain required.
