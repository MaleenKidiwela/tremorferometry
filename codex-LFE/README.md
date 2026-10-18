# codex-LFE

An independent effort to build a reliable event-and-family LFE catalog across the Cascadia margin.

**Current status: experimental workflow and regional pilots, not a validated catalog.** This folder does not import the previous project's detector, corrected reference times, family scores, or scientific conclusions. Existing original waveform files and external catalog/model downloads are reused as read-only inputs. Source hashes, parameters, date splits, and software versions accompany each experiment.

**First catalog increment:** [61 northern candidates across 17 scanned days](docs/NORTH_CATALOG_V1.md), with CSV/Parquet exports, coverage and separate development/held-out validation. [Active build plan](docs/CATALOG_BUILD_PLAN.md).

Measured results and limitations: **[CPU results](RESULTS.md)**. The completed broader pilot is `results/north_cpu_v2_parallel`.

**Current priority:** settle and validate the catalog method before regional expansion. See [the method recommendation and comparison](docs/METHOD_DECISION.md). The existing runs are baselines, not a settled production method.

## Picker validation benchmark

The next method-selection step is independent picker validation. [Protocol and commands](docs/PICKER_VALIDATION.md); [offline development review packet](results/picker_review_development_v1/review.html). The packet contains 24 unlabeled northern windows. Both jms and jcms have been run on identical waveforms, but no performance score is reported until review is supplied. Pending/ambiguous windows cannot silently become negative labels. This packet is a workflow test, not a full-margin benchmark.

## Target network

The user-designated `figures/borehole_dvv_map.html` contains **183 stations**. The exact normalized roster is [stations.csv](results/target_network_from_map/stations.csv), with source hashes and network mapping in the accompanying manifest. Earlier estimates of 198 stations are superseded by this map membership. No map family-certification labels are imported as evidence of LFE quality.

The roster-constrained northern calibration configuration is `configs/map_network_north_calibration.yaml`, containing PGC, SNB, and VGZ and enforcing that roster. NLLB and LRIV from earlier diagnostics are outside the map. Its first run is saved in `results/map_network_north_calibration_v1`; it gives 27 candidates and 23 correct-group calibration matches at threshold 9. Earlier four-station results do not validate this constrained network. Run development configurations with `--calibration-only` until a fresh evaluation period is selected.

## Run the checks

From this folder, in the existing pod environment:

```bash
python -m pytest tests -q
python -m codex_lfe.doctor --config configs/cpu_pilot.yaml
```

For a separate environment, `pip install -e '.[test]'` installs the basic dependencies. Neural inference additionally needs `torch` and `seisbench`. The matcher can run entirely on a CPU.

## Independent northern pilot

The broader follow-up uses twelve source groups and dates spread across the episode:

```bash
python -u -m codex_lfe.pilot --config configs/cpu_pilot_broader.yaml --workers 3 --output results/north_broader_new
```

Its design is recorded in `results/north_broader_design_before_run.json`. The original three-group run below is preserved as a failed coverage/power experiment; its test dates contained only eight eligible reference events for the selected groups. The broader pilot uses new evaluation dates, September 18–20, 2012. `--workers` parallelizes independent days within each split and never launches test dates before calibration completes.

```bash
python -u -m codex_lfe.pilot --config configs/cpu_pilot.yaml --output results/north_cpu_new
```

This builds templates from raw miniSEED on September 7–9, 2012, selects a threshold on September 10–11, and evaluates September 12–14. Three source groups are selected using training counts and geographic separation only. It uses four distinct recording sites. The output directory must be empty. Content-addressed preprocessing caches are reusable.

The pilot intentionally uses **arrival time at the PGC training-stack peak**, not an unverified origin-time correction. The reference table has coordinates but no original family identifiers; our IDs identify exact published coordinate groups. A scientifically reviewed original-family mapping and origin-time solution are still required.

Every threshold uses the same final event grouping and family ambiguity rule on observations and all shifted controls. The score counts each detection and reference event at most once. Unassigned events cannot earn correct-family recall. Missing-data exposure and reference exclusions are explicit. No event magnitude is fabricated from correlation or raw counts.

The threshold is frozen before scanning test dates. If calibration fails its predeclared gate, the test run uses the declared diagnostic threshold and remains a **failed-calibration experiment**, not a released catalog. Passing this gate alone would also be insufficient for release.

On a GPU pod, the same matcher supports:

```bash
python -u -m codex_lfe.pilot --config configs/cpu_pilot.yaml --backend torch --device cuda --output results/north_gpu_new
```

The torch FFT path has CPU parity tests. CUDA correctness and throughput must be measured on the GPU pod; it is not yet an optimized fleet-scale implementation.

## Neural source-discovery preparation

`codex_lfe.neural` runs local published LFEDetect weights on CPU or CUDA. It verifies channel orientation and compatible component responses against StationXML before rotating. It rejects incomplete windows and discards prediction-window edges. It emits **unassociated phase picks**, never event rows masquerading as a catalog.

The initial smoke test uses `jms` weights: the weight metadata states that training used Japan, Mexico, and San Andreas, excluding Cascadia. This is an out-of-region discovery baseline; it is not a claim that these weights perform best in Cascadia.

```bash
python -u -m codex_lfe.neural \
  --raw ../lfe_catalog_repair/traces/CN.PGC/2012/257.mseed \
  --inventory cache/metadata/cn_pilot_2012.xml \
  --channels CN.PGC..BHZ,CN.PGC..BHN,CN.PGC..BHE \
  --site PGC --date 2012-09-13 --duration-s 3600 \
  --output results/neural_cpu_new
```

Add `--device cuda` and omit `--duration-s` for a full-day GPU run. Metadata and model hashes are recorded. Our preprocessing and central-window stitching differ from the authors' full workflow and need empirical comparison.

## What remains

See [the research plan](docs/RESEARCH_PLAN.md), [data conventions](docs/DATA_CONTRACT.md), and [GPU handoff](docs/GPU_HANDOFF.md). The central and southern pilots, new-family discovery/association, validated origin times, calibrated amplitudes, and catalog release validation remain required. No existing result here supports a uniform full-margin completeness claim.

The unused-station check runs after the broader experiment:

```bash
python -m codex_lfe.heldout --run results/north_cpu_v2_parallel --config configs/heldout_station.yaml --output results/north_broader_heldout_new
```

It learns a separate LRIV template from training days and tests exact predicted lags at the frozen candidate times. It reports per-group and pooled effects against repeated shifted times, including a reference-event positive control. Low event counts remain explicitly inconclusive.

## Diagnose missed references

```bash
python -m codex_lfe.diagnose --run results/north_cpu_v2_parallel --output results/calibration_diagnostic_new
```

This inspects raw calibration waveforms with the frozen template bank. It records maximum network scores and scores satisfying the existing site-support rule within the original ±1-second reference tolerance. These reference-centered measurements are diagnostic upper bounds before peak suppression and family assignment; they are not additional detections. It never reads test waveforms or changes the frozen catalog.

Separated-training-date audit: `python -m codex_lfe.temporal_qc --run results/north_cpu_v2_parallel --output results/temporal_qc_new`. It compares early and late training stacks at fixed windows and labels insufficient coverage explicitly. Results remain post-fit diagnostics. New pilot configs can opt into `template_alignment: sinc` separately from daily `fractional_alignment: sinc`; recalibration is required.

Training-derived station removal can be tested on calibration dates alone:

```bash
python -m codex_lfe.station_ablation \
  --run results/north_cpu_v2_parallel \
  --qc results/north_cpu_v2_temporal_qc/temporal_qc.csv \
  --output results/station_ablation_new
python scripts/compare_calibration.py \
  --baseline results/north_cpu_v2_parallel \
  --variant results/station_ablation_new \
  --output results/station_ablation_new/shared_reference_comparison.csv
```

The experiment removes only measured temporal-coherence failures, preserves the template bank and clock, and reruns all original thresholds and controls. It refuses removals that eliminate the anchor or leave fewer than three sites. The comparison uses shared eligible reference IDs and verifies identical reference times and identities before matching.

Calibration-only template interpolation experiment:

```bash
python -m codex_lfe.pilot --config configs/template_sinc_calibration.yaml \
  --calibration-only --workers 3 --output results/template_sinc_new
python scripts/compare_calibration.py --baseline results/north_cpu_v2_parallel \
  --variant results/template_sinc_new --separate-training-anchors \
  --output results/template_sinc_new/shared_reference_comparison.csv
```

`--calibration-only` prevents test scans. This config changes template extraction alone; daily waveform interpolation remains linear to isolate the effect. The comparison verifies that original published reference times agree after subtracting each run's training-derived anchor, then scores each catalog on its own clock.

## Enforce the target station roster

Use the extracted CSV with explicit `network,station` columns and pass `--station-roster path/to/roster.csv` to `codex_lfe.pilot`. It validates membership before waveform loading and records the roster hash and stations used. Duplicate entries, wildcards, missing network codes, and outside stations are rejected. The authoritative map supplies 183 entries. Membership does not establish continuous operation or usable channel epochs. Earlier northern diagnostics included NLLB and used LRIV for a separate check; those stations are outside the user roster.

Channel metadata is now available for all 183 target stations. See [regional coverage](results/target_network_from_map/regional_metadata_coverage.csv); these metadata counts do not establish waveform uptime or detection sensitivity.
