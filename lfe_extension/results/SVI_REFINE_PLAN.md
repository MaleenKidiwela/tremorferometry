# SVI refine plan: Bostock (2015) families on permanent stations, 2014-01 → 2026-09

Written 2026-10-01 (read-only planning; nothing was run, downloaded, deleted or committed). For the reviewer,
then the build agent. All paths are relative to `/home/jovyan/tremorferometry/` unless they start with `results/`
(= `lfe_extension/results/`). "Scrambled control" = the same detector run with each station site's data shifted
to a different time of the same day, so nothing lines up across sites except by chance.

## 0. In short

- The extension exists and is self-consistent: catalogue **v1.1**, 115 families, 640,355 candidates ≥ 8×MAD,
  74,330 `cat_core` / 244,375 `cat_ext` events, 2014-01-01 → 2026-09-15 (`results/trackD/RESULTS.md` §0.0).
- It predates the southern standard. Missing: a per-repeat confidence fitted on non-template days with a
  scrambled control and a held-out year, PNET with/without, a station-health audit 2014–2026, rates from ΣP with
  step checks, a > 3 % data-error guard, burst muting and line notching, and cut windows.
- Proposed order: cheap audits first (they decide whether the matched filter needs a re-run), then confidence +
  PNET, then rates, then cut windows, then a delete list for the user. GPU ≈ 10–20 h (≈ 14–25 h if the matched
  filter is re-run), wall ≈ 1–2 days. Nothing can be deleted until the Washington (WA) extension decides what it
  needs: 18 of the 19 SVI stations are in its roster.

## 1. Current state (from the files)

**Detector** (`results/frozen_operating_point.json`, `results/trackC/README.md`, `scripts/c_run_month.py`)
- Network matched filter. 1,144 (family, station) templates from Bostock's 2012 events on the 2012 build half
  (alternate 2-h blocks; `lfe_catalog_repair/results/e5r_templates_2012.npz`), 116 families, 12-s templates.
- 20 Hz, 1–8 Hz Butterworth-4. Site merge "B": PGC+B009+B010+B011 = one site, LRIV+B005+B006+B007 = one site.
- Every peak ≥ 8.0 network MAD (per family-day) is stored. Declustered within a family at the template length.
  Cross-family assignment, ≥ 2-site rule, PNSN tremor flag, earthquake flag.
- HF veto at PGC: E(10–19)/E(2–8 Hz) with cut 3.6384 on BH and 3.5129 on HH (RULE-M; `results/a4_hf_veto_cuts.csv`).
- Template iteration was tested and rejected: held-out 2010 recall −0.012 … −0.017 (`results/TRACK_A_RESULTS.md` step 3).
  The 2012 templates are kept in every epoch. The per-epoch merged sets in `results/a4_templates_merged_*.npz`
  are on disk but not used.

**Stations and epochs** (`results/station_epochs.csv`, `results/trackC/tiers_by_epoch.csv`, `results/trackC/templates_by_epoch.csv`)
- 15 template stations: CN PGC, VGZ, NLLB, SNB; UW LRIV; PB B001, B003–B007, B009–B011, B926.
- Raw data without templates: CN LZB (ends 2014-06-17), YOUB (ends 2014-09-03), PFB (ends 2018-12-20), C8 SOKB
  (no data at either provider).
- Sensor changes: NLLB BH→HH 2017-02-23, VGZ sensor swap 2017-03-03, SNB sensor swap 2017-06-08, PGC BH→HH
  2017-08-10 (all one CNSN upgrade to Trillium 120 PH), and a PGC response sub-change 2015-01-12.
- LRIV was dropped from 2019-05-30: its stack cc was −0.22 against +0.54 for a matched control (A-2 step 4).
- SNB ends 2021-03-28.

| epoch | months | stations | T_core | T_ext | T_core recall (vs Bostock) |
|---|---|---|---|---|---|
| E0 | 2014-01..2017-02 | 15 | 8.75 | 8.25 | 0.146 |
| E1 | 2017-03..2018-12 | 15, new CN sensors | 9.00 | 8.00 | 0.129 |
| E2 | 2019-01..2021-03 | 14 (no LRIV) | 9.00 | 8.00 | 0.129 |
| E3 | 2021-04.. | 13 (no LRIV, SNB) | 9.25 | 8.00 | 0.122 |

Thresholds come from false alarms only: quiet-day rate ≤ 0.1 per template-day for T_core; ≤ 1.0, and ≤ 1.5 after a
slow-slip event, for T_ext (`results/trackD/RESULTS.md` §0).

**Validation against Bostock's catalogue.** His catalogue covers ETS windows only, and 2014+ has no truth at all.
- 2012 test half: recall 0.184 / 0.146 (T_ext / T_core).
- 2010 held out: 0.212 / 0.172 (ratio to 2012 1.07 / 1.09); the reverse direction gives 0.96 / 0.945.
- 2011 and 2013 (lfe_picker `results/P0O_RESULTS.md` §O2): 0.200 / 0.162 (18,663 events) and 0.199 / 0.165 (14,873 events).
- Completeness model on held-out 2010: within ±0.05 (worst case 0.024).
- Thomas main-front speed: met (12 %). Re-rupture interval: not met against his raw catalogue, but inside his
  catalogue thinned to our recall.
- Quiet-day candidates: ≥ 86 % noise at T_ext, ~22 % real repeats at T_core (Track A step 3b).
- Seasonal floor: unstable in all 8 epoch-tier cells (Track D §1c).

**Per-event score: Track E `p_real`** (`results/trackE/RESULTS.md`, `results/trackE/realness_v1_1.parquet`)
- Gradient boosting with 18 features (per-site z, confirmation at low-weight sites, burst counts, physics).
- Trained on 2012 ETS days doy 243–269. These include the template build-half blocks, so training was in-sample.
- 2010 AUC 0.855, but the calibration bar failed in level (deviation 0.261; 0.060 after a −0.93 log-odds anchor).
  It is therefore delivered as a bracket, `p_real` / `p_real_cons`.
- No scrambled control was a training class. The site-rolled control lacked rel_amp/hfr/wsh_pgc, so it could
  only be used as a reference.
- No PNET, and no model for the 13–14-station networks.

**Known weaknesses**
1. Recall is 12–21 %. The PGC site carries > 40 % of most families' weight (A step 3).
2. 14 families have no template, and 15 of 130 have no events.
3. The ETS-only truth plus the 2013 → 2014 gap leaves no overlap year to check forward calibration.
   (SH_CONF had a join year; SVI does not.)
4. Track D flagged 193 site-months above 2λ (VGZ 30, PGC 26, B001/LRIV/B926/NLLB 21) and 54 months with an
   HF-ratio shift (`results/trackC/qc_summary.csv`). None of them has been reviewed.
5. The winter excess on confident events (E0 DJF/JJA 6.70) is unresolved.
6. The watcher `scripts/c_watch.py` (PID 629057) has been idle for 14 days and is still running.

**Southern fixes: present / missing in the SVI chain**

| fix (source) | SVI status |
|---|---|
| gap / zero-energy guard (Track C cc-floor fix) | present: per-channel floor 1e-3 × station-day median, dilated by the template length (`c_run_month.py` `--energy-floor`) |
| flat or repeated-sample mask, dead-channel reference, ≤ 1 % flat per window (`lfe_south/scripts/sh_v2_common.py`) | **missing**. A stuck channel that is not exactly constant (like KMR HHN) passes |
| burst muting (K = 15, quiet neighbour sites) | **missing** |
| line notching (2–9 Hz combs) | **missing** |
| stage fails on > 3 % station-day errors (`lfe_south/results/duc_fix/DIAGNOSIS.md`) | **missing**. The cache stage logs errors and carries on ("never kill the month"). Only 1 error is logged in the whole record, but a whole station-day masked by the energy floor is never counted |
| station-health / sensor-swap audit (`SH_KMRFIX`, `SH_CHECKS`) | **missing**. Only download QC per month, and registration for 2012 and 2018-06 |
| rates from ΣP, step checks at every station change (`lfe_south/results/duc_fin/D2011.md`) | **missing**. Counts are `cat_*` at per-epoch thresholds, so they step by construction at epoch boundaries |
| cut windows −15…+30 s, verify 20/20, stack check (`lfe_south/scripts/cut_windows.py`) | **missing**. Only 16-s 20 Hz daily T_core stacks exist (`results/trackC/stacks_*.npz`, 2.0 GB) |

## 2. What "refine" should change, ranked by value per cost

| # | step | value | cost | needs raw? |
|---|---|---|---|---|
| R0 | stop the idle watcher (`touch results/trackC/STOP`); freeze v1.1 as the baseline | removes a silent writer | 0 | no |
| R1 | data-completeness audit, 153 months | decides whether a re-run is needed | CPU < 15 min | no |
| R2 | station-health audit 2014–2026 (KMR method) | catches dead/noisy/swapped sensors | CPU 2–4 h | **yes** |
| R3 | A/B test of the southern cleaning (flat mask, burst mute, line notch) | decides R4 | CPU ~1 h + GPU ~1 h | cache |
| R4 | (only if R3 says so) full matched-filter re-run v2 | fewer false alarms | GPU 3–5 h | cache |
| R5 | per-repeat confidence to the SH_CONF / DUC_FIN standard (replaces Track E) | the core product: P per repeat | GPU 6–12 h | cache + raw (HF) |
| R6 | PNET evidence, reported with and without | user request | GPU ~4 h | cache |
| R7 | firing rate per family from ΣP, with step checks | goal (0): rates | CPU < 1 h | no |
| R8 | cut windows for every candidate and every station + verify + stack check | the dv/v product, and lets raw go | CPU 2–4 h | **yes** |
| R9 | delete list for the user (nothing deleted) | storage | CPU minutes | no |

**R1, completeness** (pattern: `duc_fix/DIAGNOSIS.md`)
- Per station-day, join `results/download_manifest.csv` (status `data`), the `cache20hz/` npz, `run_<ym>.json`
  `n_windows_masked_by_energy_floor` and `family_days_<ym>` live stations.
- Report per station-month the fraction lost at each stage. A stage-month fails at > 3 %.
- Output `results/refine/r1_completeness.csv`. If the losses are > 3 %, rebuild those days before anything else.

**R2, station health** (`lfe_south/results/SH_KMRFIX_RESULTS.md`, `SH_CHECKS_RESULTS.md`)
- From raw counts, per channel per day:
  - DC level and drift; longest constant run; flat fraction; std;
  - 2–8 Hz and 0.1–0.5 Hz levels in nm/s (StationXML gain);
  - line-comb prominence in 2–9 Hz; sample rate, channel and location codes.
- Flag days against the station's own running median.
- Per station-month from the saved `stacks_*.npz`:
  - stack cc against the template and SNR, against a random-offset scrambled control;
  - best lag within ±2 s (timing), sign (polarity), channel permutation (swap), best horizontal rotation.
- Review first: the 193 site-months and 54 HF-shift months above; LRIV before and after 2019-05-30 (is the
  negative cc a polarity or component swap that could be fixed?); B004 2020 (4 GB on disk against ~9 GB in
  other years); B011 absent 2018-06-10..25.
- Output: `results/refine/r2_health_daily.parquet`, `r2_health_monthly.csv` and a `bad_for_dvv` window list in
  the KMR style. This step must finish before any raw file is listed for deletion: the KMR check could not be
  completed after its raw data were gone.

**R3/R4, cleaning**
- Port `sh_v2_common.prep_day`: flat mask, dead references, burst mute, line notch.
- Run it with the frozen templates and thresholds on the truth years (2010–2013) plus 6 production months
  (one per epoch, one winter, one summer).
- Adopt v2 for the full record only if, at a matched scrambled-control rate, it keeps ≥ 85 % of the Bostock
  matches (veto rule, memory `justify-criteria-on-known-families`) and lowers the control rate.
- Otherwise keep v1.1 and store the cleaning flags as columns (rank, do not discard).
- If v2 is adopted: re-derive the tiers with `c_tiers_by_epoch.py --rederive`, and cut windows for the union of
  the v1.1 and v2 candidates.

**R5, confidence** (copy `lfe_south/scripts/sh_conf_*.py`, adapt)
- Candidates: every ≥ 8×MAD peak (640,355 in v1.1), HF-vetoed and cross-family losers included. Nothing is cut.
- Labels: Bostock event of the same family (`m_same`) or any family (`m_any`) within ±2 s; one model per label.
- Features:
  - network: net_mad, net_cc, rel_amp, n_sta, n_sites, live channels;
  - per site and per station: cc and z against that site's own random origins;
  - hfr_rel; density within ±1 h and ±10 min, and other families within ±4 s;
  - the 2-site low-weight confirmation (from Track E);
  - PNET (R6); family id.
  - Never tremor or season.
- Scrambled control: one per day (`c1`), same detector, same features, every month and every truth year.
- Day roles, fixed before fitting:
  - **in-sample, never fitted or scored:** the 2012 build-half 2-h blocks.
  - **fit, out-of-fold by day:** 2010, 2011 and the 2012 test blocks.
  - **held out:** 2013 (decision D2).
  - **PNET fine-tune:** every 10th truth day, never fitted or scored.
- **The SVI-specific step: station-set strength.** Re-run the truth years with the S14 and S13 station sets
  (Track D already did 2010/2012 in `scripts/d1_run_ablation.sh`; add 2011/2013), so the model sees each
  production network on days with truth. Bars per station set: AUC; calibration slope 0.8–1.2;
  matched/ΣP 0.8–1.2; control P > 0.5 ≤ 1 %.
- The 2017 sensor swaps cannot be simulated on truth days. Check them through R2 and the R7 step tests.
- Model: HistGradientBoosting + isotonic, GroupKFold by day, ablations FULL / NO_PNET / NO_FAM / MF_ONLY.
- Output: `results/refine/conf/P_by_candidate.parquet`, keyed `month` + `det_row` = Track E's `det_uid`.
- Reuse of Track E: keep its per-site pass (`e_persite.py`, self-check ≤ 9e-8) as the feature source, but replace
  its model. It was trained in-sample on 2012, had no proper scrambled control and has no station-set handling.
- Caveat: P = probability of matching Bostock's ETS-only catalogue. Between ETS windows it is a lower bound on
  "real", not a measure of it.

**R6, PNET**
- Weights: `lfe_picker/results/p1/models/pnet_w05.pt`. They are already SVI-trained, on 2008 + 2009 + 2010 + 2012
  at S15 minus B926/NLLB/LRIV, without the 30 hidden families (`lfe_picker/results/P1_RESULTS.md` l.41).
- So PNET values on 2010/2012 at the training stations are in-sample. Score the with-PNET model on
  2011/2013, the 3 untrained stations and the 30 hidden families only.
- Report FULL against NO_PNET, off-the-shelf against fine-tuned (fine-tune on 2011 tenth-days at all 18 stations),
  and permutation importance, as in SH_CONF §4.
- Input: 20 Hz 1–8 Hz cache with PB EH1/EH2 rotated to N/E from StationXML azimuths. The frozen chain does not
  rotate.

**R7, rates**
- Per family per day: ΣP(same) and ΣP(any) with Bernoulli + fold spread, the scrambled-control ΣP, and the old
  `cat_core`/`cat_ext` columns.
- Firing rate ≈ ΣP(same) × (his/ΣP for that family and station set, measured on truth days).
- Step checks: ±21 d and ±45 d ratios at every change, against the same calendar windows in other years
  (`t4_step_checks.csv` method). Changes: 2014-06-17, 2014-09-03, 2015-01-12, 2017-02-23, 03-03, 06-08, 08-10,
  2018-12-20, 2019-05-30, 2021-03-28, plus every R2 health flag.
- Do not over-read the steps: a step next to an SSE is not separable from the activity (SH_KMRFIX §4).

**R8, cut windows**
- Profile `svi_ext` in a copy of `cut_windows.py`. Window −15…+30 s around the expected S at each station:
  t + ttP + (S−P) + moveout (`results/trackC/README.md`). Offset 0 for stations without a template (PFB, LZB, YOUB).
- Raw counts at native rate, all channels, HDF5 per month + parquet index.
- Include the truth years 2010–2013 around Bostock's catalogued events, which ties the dv/v series back to 2010.
- Verify 20 random windows sample-identical against the miniSEED, and run the stack check per family.

## 3. Data on disk, what the refine needs, delete draft

| store | size | content | needed by refine | other readers |
|---|---|---|---|---|
| `lfe_extension/traces/` 2014–2026 | 1,820 GB | raw miniSEED, 19 stations (largest: VGZ 159, PGC 148, B003 129, B926 129, B010 119, LRIV 117, NLLB 115, B009 115) | R2, R5 HF features, R8 | **WA**: `lfe_extension_wa/scripts/wa_common.py` TRACE_ROOTS; 18 of 19 stations (all but NLLB) are in `lfe_extension_wa/results/step4b/setup/stations.csv`. **lfe_array_discovery**: 38,185 symlinks (T7 "extension" tranche, B001, B003–B007, LRIV, VGZ, PFB, 2007–2026 = 991 GB of this). **lfe_picker**: 194 symlinks (2023–2025) |
| `lfe_extension/traces/` 2010, 2011, 2013 | 76 GB | ETS truth windows | R5 (fit / held out), R8 | lfe_picker (2011/2013 cache source) |
| `lfe_extension/cache20hz/` | 1,694 GB | 20 Hz 1–8 Hz float32, 2014–2026 (+ B009/B010 2010) | R3–R6 | `lfe_picker/scripts/q_pairlink2.py` (3rd fallback cache; no lfe_picker job running) |
| `lfe_catalog_repair/traces` + `cache20hz` | 99 + 31 GB | 2010 / 2012 / 2017 truth windows | R5 | read-only input of several projects. **Not ours to delete** |
| `lfe_picker/cache20hz` 2011, 2013 | part of 416 GB | same chain (`o_cache_month.py`) | R5 | lfe_picker |
| `lfe_extension/results/` | 3.7 GB | catalogue, Track C detections / family_days / stacks (2.0 GB), Track E | all | data/catalogue/families_v1_2.csv rows SVI:0..129 |

fault_tomography, `scripts/` (the dv/v fleet) and `data/` contain no reader of these paths
(`data/extension_stations/on_disk_traces.csv` is an inventory and will go stale).

**Delete draft** (to be produced by a `deletelist` stage with a needs-check, as in `cut_windows.py`; nothing goes
before R2 + R8 are verified and the user approves the exact list)
1. `lfe_extension/cache20hz/` (1,694 GB): regenerable from raw (`scripts/c_cache_month.py`). Delete after R3–R6.
2. `lfe_extension/traces/CN.NLLB/2014..2026` (115 GB): SVI is its only reader.
3. Hold: the other 18 stations, 2014–2026 (1,705 GB). Their only remaining readers are WA's planned scan to
   2026 and the array_discovery / picker symlinks. Delete them only once WA has cut its own windows, or the
   user accepts a re-download for WA (decision D6).

Keep: the 2010/2011/2013 truth raw (76 GB), all of `results/`, the cut store, the logs.

## 4. Compute, storage, queue

| stage | resource | estimate | basis |
|---|---|---|---|
| R1 | CPU | < 15 min | file joins |
| R2 raw pass | CPU, 24 workers | 2–4 h | ~69,800 station-days, 1.82 TB read |
| R2 stacks | CPU | < 30 min | `stacks_*.npz` 2.0 GB |
| R3 | CPU ~1 h + GPU ~1 h | | truth years + 6 months |
| R4 (optional) | GPU | 3–5 h | Track C 60–95 s per month uncontended × 153 |
| R5 features + scrambled control + station-set reruns | GPU | 6–12 h | Track E per-site pass 9.7 GPU-h; SH_CONF 71 min per 21k station-days (SVI ≈ 3.4× that) |
| R6 PNET | GPU | ~4 h + ~5 min fine-tune | SH_CONF 73 min per 21k station-days |
| R7 + fit / eval | CPU | < 1 h | SH_CONF fit 4.5 min |
| R8 cut + verify + stack check | CPU, 6–12 workers | 2–4 h | Shelly 521 GB in 32 min; Ducellier 708 GB in 86 min |

- Totals: GPU ≈ 10–20 h (≈ 14–25 h with R4); wall ≈ 1–2 days with CPU stages overlapped.
- Storage during the run: +0.1–0.3 TB (prep masks 50–150 GB if R3/R4; features 5–15 GB; cut store ~30–130 GB:
  Shelly 1.5 % of raw, Ducellier 7 %). There are 65 TB free. A pilot month gives the real cut size first.
- Queue: launch-and-exit shell queue `scripts/svi_refine_queue.sh` with markers `results/refine/SVI_REFINE_{PID,HEARTBEAT,DONE,FAILED}`
  and a `flock` lock.
  - GPU stages wait until `lfe_south/results/DUC_FIN_DONE` exists. It appeared 2026-10-01 21:51:56Z, and the GPU
    was idle at 22:02Z.
  - They also wait for no other GPU queue (e.g. a Ducellier follow-up) to hold a lock. At most 3 concurrent
    matched-filter passes on the H100 (Track A-2 OOM).
  - Every data-reading stage applies the > 3 % station-day error guard.
  - CPU stages R1/R2 can start now.

## 5. Open decisions for the user (each with a recommendation)

- **D1. Re-run the matched filter with the southern cleaning?** Recommend: decide from R3's A/B on the truth years.
  Re-run only if it passes the 85 % known-family veto and lowers the scrambled-control rate. Otherwise keep v1.1
  and add the flags as columns.
- **D2. Year roles.** Recommend: fit on 2010 + 2011 + the 2012 test blocks (out-of-fold by day) and hold out all of
  2013. 2013 was lfe_picker's sealed referee; using it here is a separate study, but it needs your OK.
  Alternative: hold out 2011 and add 2013 to the fit.
- **D3. Model the network-strength effect with S14/S13 reruns on the truth years.** Recommend yes. Without it, P
  for E2/E3 (60 % of the record) is an extrapolation, and Track E failed in exactly this way (log-odds shift −0.93).
- **D4. PNET.** Recommend: the off-the-shelf SVI model plus a 2011 fine-tune, both reported with and without PNET,
  and scored only where PNET never trained (2011/2013, B926/NLLB/LRIV, the 30 hidden families).
- **D5. LRIV after 2019-05-30.** Recommend: keep it dropped unless R2 finds a fixable cause (polarity, component
  swap, orientation). In that case re-template it from 2019–2020 high-P repeats and test against B926 as the
  no-change control.
- **D6. Raw deletion scope.** Recommend: after R8 delete the cache (1.69 TB) and the NLLB raw (115 GB). Hold the
  other 18 stations (1.71 TB) until the WA scan to 2026 has cut its windows; re-downloading would take ~1 day at
  the 63 MB/s seen in Track B. Alternative: delete all now and let WA re-fetch.
- **D7. Firing-rate definition.** Recommend: rate = ΣP(same) × his/ΣP per family and station set, with the
  scrambled-control ΣP subtracted and shown. Keep `cat_core`/`cat_ext` counts and Track D's λ and Mc as columns.
  Between ETS windows, report the rate as "Bostock-like repeats" (a lower bound).
- **D8. Out of scope here:** template iteration (rejected on held-out data), the 14 template-less families,
  relocation, and the far-station-only iteration (proposed in Track A step 3, never run; park it unless the
  confidence results ask for it).
