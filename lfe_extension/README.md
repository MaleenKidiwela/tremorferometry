# lfe_extension

Extension of the Bostock et al. (2015) southern-Vancouver-Island LFE family catalog (130 families,
2003–2013) forward to 2026 on the permanent network, by network template matching with iterated
templates, plus independent family locations. See `PLAN.md` (design, bars, tracks), `NOTES.md` (dated
ledger), `scripts/`, `results/`, `figures/`, `logs/`. Raw waveforms in `traces/<NET.STA>/<YYYY>/<DOY>.mseed`.

Inputs (read-only): `../lfe_catalog_repair/` (E4–E6 code, 2010/2012/2017 traces, Bostock + SeisBench
references), `../codex-LFE/` (station metadata, audit), `../catalogs/pnsn_tremor_cascadia_full.csv`.

## Track B — raw waveform download (running)

| file | what |
|---|---|
| `results/station_epochs.csv` | one seismometer band per station epoch (HH > BH > EH), from IRIS StationXML |
| `results/epoch_crosscheck.csv` | vs `../codex-LFE` triplet epochs (0 disagreements) |
| `results/download_queue.csv` | ordered station-days: P0 2010 ETS (PB) → P1 2011 ETS → P2 2013 ETS → P3 continuous 2014-01-01..2026-09-15 |
| `results/download_manifest.csv` | one row per finished station-day — **this is the file Track C follows** |
| `results/download_progress.json` | refreshed every 5 min (per station-year counts, GB, MB/s, ETA) |
| `results/download_qc_<YYYY-MM>.csv` | per-month, per-station QC, written automatically when a month completes |
| `figures/download_coverage.png` | station × month usable fraction |

```bash
bash scripts/download_run.sh                 # start/restart (resumes; never re-fetches a present day)
tail -f logs/download_run.log                # live log
cat results/download_progress.json           # ETA, MB/s, current month
kill $(cat logs/download_run.pid)            # stop, safely resumable
```

`status` in the manifest: `data` | `nodata` (verified empty, `.mseed.nd` marker) | `error` (3 tries +
one retry pass) | `present_elsewhere` (already in `../lfe_catalog_repair/traces`, path in `note`) |
`no_epoch` (station dead). `coverage_frac` = min over channels of seconds/86400; `n_frag` = max
fragments over channels.

## Track A — method (steps 0–3 + 3b done, 2026-09-16)

Full write-up: **`results/TRACK_A_RESULTS.md`**; dated ledger in `NOTES.md`; the spec executed is
`results/merlin_review_plan_2026-09-16.md`. Every script used was copied from `../lfe_catalog_repair/scripts`
into `scripts/` and edited there; the parent folders were never written to.

**Frozen operating point — `results/frozen_operating_point.json`** (what Track C must use):
site-merge **B**, **12-s** templates, **T_ext = 8.25 MAD**, **T_core = 8.75 MAD**, **HF-ratio veto ON**
(PGC E(10–19 Hz)/E(2–8 Hz) over the scored window, cut 3.6384), one-to-one ±2 s matching, ≥ 2-site rule,
measurement chain unchanged (20 Hz, 1–8 Hz Butterworth-4 zero-phase). Own recall of Bostock's 2012 test half
0.184 at T_ext / 0.146 at T_core; ETS-day own precision 0.49 / 0.60.

| file | what |
|---|---|
| `results/e6_matched_fa_roc.csv` | step 0: own recall at matched quiet-day FA for all E6 variants (reproduces Merlin's table) |
| `results/e6_rescore/` | step 0: all nine E6 passes re-scored with the corrected `--e6` key labelling, nearest-neighbour and one-to-one |
| `results/a1_hf_ratio_*_decision.csv` | step 1: HF-veto FA reduction vs recall cost |
| `results/a1_frozen_fa_vs_T_*.csv` | step 1: false alarms vs threshold (quiet and post-SSE) |
| `results/a2_transfer_table.csv`, `a2_reverse_*.csv` | step 2: frozen 2012 → 2010 transfer and the reverse |
| `results/a3_reject_tests*.csv` | step 3: the six template-iteration rejection tests (+ the frozen baseline) |
| `results/a3b_quietday_stack_test.csv` | step 3b: withheld-station test of the quiet-day detections |
| `figures/a0_…, a1_…, a2_…, a3_…, a3b_…` | ROC, FA vs T, HF-ratio distributions, 2010-vs-2012 recall, iteration, quiet-day stack test |

## Track C — production catalog (built 2026-09-16, NOT yet started on the full record)

Full schema and restart instructions: **`results/trackC/README.md`**. The runner reproduces the frozen
Track A configuration **bit for bit** (2012 quiet 253/253 and ETS 3426/3426 detections, Δ = 0).

| script | what |
|---|---|
| `scripts/c_month_ready.py` | is a month runnable from `results/download_manifest.csv`? which stations contribute? |
| `scripts/c_cache_month.py` | 20 Hz cache for a month (10–12 workers) + `results/trackC/daily_noise_<YYYY-MM>.parquet` |
| `scripts/c_run_month.py` | the network matched filter for a month → `detections_`, `family_days_`, `stacks_`, `run_` |
| `scripts/c_qc_month.py` | monthly QC → `qc_<YYYY-MM>.csv`, `qc_summary.csv`, `figures/trackC_monthly_qc.png` |
| `scripts/c_registration.py` | per-station registration of the T_core stacks against the templates (epoch check) |
| `scripts/c_eq_catalog.py` | one-off regional earthquake catalog for the EQ flag |
| `scripts/c_watch.py` + `scripts/c_run_all.sh` | resumable single-writer watcher: cache → run → QC → `months_done.csv` |

```bash
nohup bash scripts/c_run_all.sh > /dev/null 2>&1 &   # start the production watcher
tail -f logs/trackC_watch.log ; cat results/trackC/months_done.csv
kill $(cat logs/trackC_watch.pid)                    # stop — resumable at any point
```

## Track D — the product (completeness model, catalog, science QC)

Full write-up: **`results/trackD/RESULTS.md`**; the spec executed is steps 6–8 of
`results/merlin_review_plan_2026-09-16.md`. Everything is `scripts/d*`, `results/trackD/`, `figures/d*`;
Track C's runner and outputs are read, never written. `d3_thomas_harness.py`, `d3_thomas_base.py` and
`d3_thinned_reference.py` are COPIES of `../lfe_catalog_repair/scripts/e5_thomas_harness.py`,
`e4_thomas_harness.py` and `e6_thinned_reference.py`.

```bash
bash scripts/d1_run_ablation.sh              # D1(a): the 4 station-set GPU passes (referee years)
bash scripts/d1d_run_withheld.sh <YYYY-MM>…  # D1(d): the withheld-station passes for the seasonal test
bash scripts/d_run_all.sh                    # everything else, from whatever Track C months are on disk
python scripts/d4_results.py                 # regenerate results/trackD/RESULTS.md from the outputs
python scripts/d_sweep_qc.py --check-ready 8 # read-only anomaly scan of the Track C sweep
python scripts/d3b2_front_v2.py --stage referee      # §3b-2: design + freeze on the 2012 referee
python scripts/d3b2_front_v2.py --stage production --alpha 0.05 \
       --catalog results/trackD/lfe_extension_catalog_v1_1.parquet \
       --out results/trackD/thomas_front_v2_v1_1.csv        # then every SSE window
python scripts/d2_assemble.py --tag v1_1                     # a new catalog version, v1 untouched
python scripts/d3_run_sse.py --catalog …_v1_1.parquet --suffix _v1_1
python scripts/d5_family_cumulative.py                       # per-family cumulative + hp curve
python scripts/d4_results.py --tag v1_1 --prev-tag v1        # RESULTS.md with the change log on top
python scripts/d3b2_figures.py                       # the front-speed summary figure
python scripts/d_fix_ledger.py [--apply]     # reconcile months_done.csv with the catalog
```

**Epochs and thresholds come from Track C.** `scripts/d_tiers.py` is the single reader of
`results/trackC/tiers_by_epoch.csv` (network epochs + the per-epoch FA-only T_core / T_ext) and maps
each epoch onto the STATION SET the ablation measures — S15 (15 stations), S14 (LRIV dropped),
S13 (LRIV + SNB gone). Every D1–D4 step reads epochs and thresholds through it, so the completeness
model, the floor λ and the catalog can never disagree with the `cat_core` / `cat_ext` flags.

### The product

**Current version: v1.1** (rebuilt 2026-09-17 after Track C's cc-floor fix and the per-epoch tier
re-derivation — see §0.0 of RESULTS.md). v1 is preserved unchanged under its own filenames, and its
completeness / floor / Thomas outputs are in `results/trackD/v1_snapshot/`.

| file | what |
|---|---|
| `results/trackD/lfe_extension_catalog_v1_1.parquet` (+ `.csv`) | **the catalog** — one row per event, both tiers, with Mc(f,d) and λ(f,T,epoch,season) attached |
| `results/trackD/family_table_v1_1.csv` | per family: Bostock hypocentre + SeisBench id, our relocation and its verdict, template stations, per-epoch recall, template-less flag |
| `results/trackD/family_day_table_v1_1.parquet` | per (family, day): observable, n_sites, Mc, λ, day terms |
| `results/trackD/catalog_v1_1_meta.json`, `catalog_month_crosscheck_v1_1.csv` | counts, time basis, per-month agreement with `months_done.csv` |
| `results/trackD/family_cumulative_daily_v1_1.parquet`, `figures/family_cumulative_v1_1/`, `figures/d5_family_cumulative_grid_v1_1.png` | per-family cumulative counts 2011–2026 with Track E's re-scored high-purity curve (`scripts/d5_family_cumulative.py`, a parameterised copy of Track E's) |
| `results/trackD/lfe_extension_catalog_v1.parquet`, `family_table_v1.csv`, `family_day_table_v1.parquet`, `v1_snapshot/` | the previous version, untouched |

### The completeness model (Merlin Step 6)

| file | what |
|---|---|
| `results/trackD/completeness_family_epoch.csv` | (a) per family per Track C epoch, scored at that epoch's own thresholds — the table the catalog and family table merge on |
| `results/trackD/completeness_family_stationset.csv` | the ablation itself: per family per station set at one common threshold, plus the matched-FA recall |
| `results/trackD/completeness_{epoch,by_epoch}_summary.csv`, `_crossyear.csv`, `_matched_fa.csv`, `tiers_used.json` | pooled recall per set and per epoch, the matched-FA thresholds, the held-out-2010 cross-year check, and the epoch/threshold table actually used |
| `results/trackD/mc_family_day.parquet`, `mc_model_params.csv`, `mc_model_validation.csv` | (b) per-family-day Mc at both tiers, the fitted probit, and the pre-registered ±0.05 out-of-sample validation |
| `results/trackD/floor_family_epoch_season.csv`, `floor_epoch_season_summary.csv` | (c) the false-alarm floor λ per family / tier / epoch / season, corrected by the Track A step-3b real-repeat fraction |
| `results/trackD/daily_noise_ref.parquet` | referee-year daily noise, same estimator as Track C |
| `results/trackD/floor_seasonal_ratio.csv`, `season_diagnostic_*.csv` | (d) the per-epoch winter/summer ratio and what drives it: quiet-day tail shape, ambient noise vs MAD scale, within-season noise slope |
| `results/trackD/season_stack_test*.csv` | (d) the withheld-station test per season — Track A step 3b's test with a pooled incoherent statistic and a new moveout-shuffled null |
| `figures/d1_completeness_model.png`, `figures/d1d_season_diagnostic.png` | validation, per-epoch recall, family dependence, Mc, day terms, seasonal floor and its drivers |

### Science QC (Merlin Step 8)

| file | what |
|---|---|
| `results/trackD/sse_windows.csv`, `sse_windows_sub.csv` | SSE windows from PNSN tremor in the family box (≥ 50/day, gaps ≤ 20 d / ≤ 5 d sensitivity, span ≥ 10 d, ±3 d pad) |
| `results/trackD/thomas_per_sse.csv` | the Thomas harness per SSE per tier (`cat_core`, `cat_ext`, and a `cat_core_hp` ≥ 5-site purity leg) + the 2012 referee control |
| `results/trackD/thomas_bars.csv`, `thomas_estimator_yield.csv` | the two pre-registered bars on the 2012 referee, and how often each estimator works at this completeness |
| `results/trackD/thomas_front_v2*.csv/json`, `figures/d3b2_front_v2_summary.png` | §3b-2 (post-registration, approved 2026-09-17): the purity-matched main-front estimator — its 2012-referee design, the frozen α, and the per-SSE speed / CI / direction / window / nulls |
| `results/trackD/catalog_ledger_disagreements.csv`, `sweep_qc_by_month.csv` | the month-by-month cross-check against `months_done.csv` and the sweep anomaly scan |
| `results/trackD/thinned_reference_*.csv` | Bostock 2012 thinned per family at OUR recall + OUR quiet-day floor — the fair reference for every decimation-sensitive metric |
| `figures/d3_sse_<id>_<tier>.png`, `figures/d3_thomas_summary_<tier>.png` | space-time plot per SSE, and the across-SSE summary |
