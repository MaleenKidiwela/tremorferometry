# Merlin review of PLAN.md — 2026-09-16 (verbatim; Fable 5.1 as Merlin)

**Verdict**: The approach (network matched filter seeded by Bostock's 130 families on the permanent stations, with iterated templates) is viable as what the on-disk evidence says it is — a PGC-centric detector with ~15–20 % family-dependent recall of Bostock-class events at ~45 % ETS-day precision — and NOT as a Bostock-equivalent continuous catalog. The plan repeats the E6 error of comparing variants at fixed MAD instead of at matched false-alarm rate, carries a chance-control bar that is unfalsifiable by construction, and, most seriously, pre-registers a false-alarm budget (≤1/template/day) that for a 12-year continuous run yields a catalog that is >85 % false alarms. The single decision that matters most: define the product as a two-tier catalog with a per-event false-alarm rate attached, choose thresholds by FA alone, and re-scope the science to what ~20 % completeness supports before any GPU pass is called "Track A".

## Flaws (severity-ordered)

**F1 — FATAL as pre-registered: the FA budget is an ETS-window bar applied to a continuous catalog.** PLAN.md (≤1/template/day); `e6_fa_table_e6_mergeB.csv` (quiet FA 0.97/template/day at T=8, own recall 0.189). 116 templates × 4,600 days × ~1/day ≈ 530k false detections 2014–2026, against ~5–6k real detections per major SSE × 12 SSE windows in the S-VI box (PNSN tremor: 2015-12, 2017-03, 2018-06, 2019-08, 2019-11, 2020-11, 2021-01, 2021-10, 2022-03, 2023-04, 2024-08, 2025-09) ≈ 60–70k. Inter-SSE "healing", tidal and re-rupture analyses would run on a catalog ~90 % noise outside SSE windows. Bostock's own catalog has no events outside ETS windows, so "Bostock's way" was never continuous.

**F2 — FATAL as a bar: "rolled control ≤5 % of ETS-day detections with the SITE-rolled control" is unfalsifiable.** `e4_template_test.py` line 665 rolls each whole site by one shift, so within-site coherence is preserved; a real LFE detected mainly by the PGC site survives the roll. `e6_fa_table_e6_sites.csv` T=10: ETS-day site-rolled 1.12 vs data 2.17 (51 %) while any-family precision is 0.85. The site-rolled control measures whether BETWEEN-site coherence is necessary, not chance. On quiet days it is the right noise-model check (site-rolled/data 0.93–1.01 for every variant). Every variant fails this bar forever.

**F3 — MAJOR: site-merge conclusions were fixed-T comparisons; at matched quiet-day FA the ranking changes.** Recomputed from `e6_fa_table_*.csv` (log-interpolated own recall of Bostock test-half events, 2012, 15 stations):

| variant | recall @ FA 1.0 | @ FA 0.5 | @ FA 0.1 | own precision (ETS) @ FA 1.0 |
|---|---|---|---|---|
| sites (no merge) | 0.165 | 0.140 | 0.092 | 0.42 |
| mergeA (site w = max) | 0.171 | 0.144 | 0.105 | 0.47 |
| **mergeB (site w = √Σw²)** | **0.190** | **0.163** | **0.117** | 0.45 |
| mergeB + cap PGC 25 % | 0.167 | 0.145 | 0.099 | 0.47 |
| mergeB + 2-site coincidence | 0.181 | 0.156 | 0.113 | 0.49 |
| lev_a (wider far-borehole pick) on B | 0.187 | 0.161 | 0.116 | 0.45 |
| lev_b (≥5 dB stations only) on B | 0.199 | 0.169 | 0.125 | 0.44 |
| lev_c (2010+2012 templates) on B | 0.199 | 0.172 | 0.125 | 0.46 |

mergeA gains nothing (re-parameterisation of T); mergeB is a real +15 % relative ROC shift; capping and coincidence gating lose; levers b/c add ~+5 %. Best on-disk configuration: 0.20 at 1 FA/template/day, 0.125 at 0.1. No waveform-level merge is needed or tested.

**F4 — MAJOR: the E6 score JSONs on disk are mislabeled.** `e5r_score.py` line 162 builds the key as `L{L}_T{int(T)}` unless `--e6` is passed; keys `L8_T7, L8_T8, L8_T9, L8_T10` for T_ops {7.5, 8, 8.5, 9, 10}, so "T7" holds 7.5 and "T8" holds 8.5 (8.0 overwritten). Use the `e6_fa_table_*.csv`; re-score with `--e6` before anything is written up.

**F5 — MAJOR: the ≥30 % recall bar is unsupported; A2 is misframed as a sensitivity lever.** Best matched-FA recall 0.20 on 15 stations in 2012. Recall is flat in template event count above ~50; PGC stack S-SNR saturates ~12 dB by N≈100; templates already have N_med 115–160. Post-2021 network is weaker than 2012: 8 sites, median nearest site 22 km, 20/130 families with azimuthal gap >180° within 80 km, 28 with <4 sites within 60 km. Iteration's real job is TRANSFER (new instrument epochs; 14 template-less families), not recall. Expect ≤+0.02.

**F6 — MAJOR: Track C preprocessing differs from everything validated.** Validated chain (`e4_cache20.py`): 20 Hz, 1–8 Hz Butterworth-4 zero-phase, component order Z/N|1/E|2 resolved but NOT rotated. Freeze on the validated chain, or treat 2–8 Hz / 40 Hz as a pre-registered variant at matched FA on 2012, then freeze. Rotation to ZNE is harmless for same-station matched filtering while azimuths are constant within an epoch; do it for the relocation stacks.

**F7 — MAJOR: "operating point chosen on 2012 quiet days + Bostock recall" tunes to the referee.** T by FA criteria alone; recall reported at the resulting T.

**F8 — MAJOR: single-station family-days are scored** (`e4_template_test.py` 577–579). Needs a ≥2-site observability rule and an "unobservable" flag.

**F9 — MAJOR: completeness model lacks the day-level sensitivity term and the FA floor.** Referee windows are Aug–Oct; 6 of 12 extension SSEs are Nov–Apr. MAD threshold is noise-adaptive → completeness varies by season/epoch. Relative amplitude is in units of each station's quiet-day rms → shifts at every response change.

**F10 — MAJOR: relocation as planned inherits the E3 pathology.** (i) PGC/B009/B010/B011 = ONE site (one travel time). (ii) Far boreholes/LRIV have systematic S residuals vs Crosson +1.7 to +3.25 s (B001 +2.35, B003 +2.98, B004 +3.25, B005 +1.70, LRIV +1.70); corrections from Bostock's hypocentres would be circular. (iii) E3's 3.1 km was on the 21-station dense-array stacks; no synthetic floor for our geometry. (iv) 5 km bar conflates horizontal and depth.

**F11 — MINOR: matching is nearest-neighbour, not one-to-one** (`e4_template_test.py` 771–773). Implement greedy-by-time one-to-one; report the delta once.

**F12 — MINOR: event time base.** Bostock DB time + per-family offset (+2.5 to +20.5 s). Extension needs its own origin convention once A3 passes.

**F13 — MINOR: roster arithmetic.** LZB (→2014-06-17), YOUB (→2014-09-03), SOKB (→2015-01-18) cover NO extension-period SSE (first is 2015-12-19); referee years only. PFB (→2018-12-20) covers 3 SSEs, SNB (→2021-03-28) 8. Epochs (`codex-LFE/results/target_network_from_map/channel_triplet_epochs.csv`): PGC BH to 2017-08-10 then HH (+ service gap 2015-01-12); VGZ break 2017-03-03 (inside cached 2017 window); SNB 2011-03-17 and 2017-06-08.

**F14 — MINOR: EQ veto by catalog only.** Quiet-day detections 06–11 UTC concentrated and day-clustered — consistent with uncatalogued micro-earthquakes or site transients. SPECULATIVE: an HF-ratio discriminator (10–20 Hz over 2–8 Hz energy at the PGC site) is a cheap physics-based veto to test.

## Prescription (ordered so the cheapest step can kill the approach)

**Step 0 — Freeze the measurement chain (CPU, hours).** Keep 20 Hz / 1–8 Hz / existing component handling for all of Track A; one-to-one matching; re-score all eleven E6 passes with `--e6` (fixes F4); write the matched-FA ROC table into `lfe_extension/results/` as the baseline record; add a `≥2 sites with usable templates` rule to the network loop (F8), record `n_sites` per family-day. Accept: re-scored table reproduces the numbers above within ±0.005. Failure → audit `e5r_score.py` before proceeding.

**Step 1 — A1 design decision and operating points (one GPU pass ~10 min for 12-s).** Adopt mergeB (within-site weighted-mean cc for {PGC,B009,B010,B011} and {LRIV,B005,B006,B007}; site weight √Σw²; weights w = 10^(S-SNR/20) capped 3× median — `e5r_weights_snr_2012.csv`); no PGC cap, no coincidence gate. Run mergeB with 12-s templates (must be shown at matched FA). Test the HF-ratio veto on existing mergeB detections: cut at the 99th percentile of build-half own-matched detections; evaluate on quiet days and test-half recall. Define two thresholds by FA ONLY: T_ext = smallest T (0.25 steps) with quiet-day FA ≤1.0 AND post-SSE (doy 270–279) FA ≤1.5 per template-day; T_core = smallest T with quiet-day FA ≤0.1. Report recall at both; never adjust. Accept: mergeB(+12 s) own recall at FA 1.0 in [0.18, 0.24]; quiet-day site-rolled/data in [0.77, 1.3]. HF veto adopted only if quiet FA falls ≥30 % at ≤2 % test-half recall cost. Failure: recall <0.16 → does not reproduce E6; diff against `e6_run_gpu1.sh`.

**Step 2 — Cross-year transfer with FROZEN templates and thresholds (2010 data on disk; make-or-break).** First run 2012 with only the 7 persistent stations (PGC-site, VGZ, SNB, NLLB, B011, B926, LRIV) → ablation prediction. Then run 2012-built templates for those stations on 2010 (doy 222–256), frozen T_ext/T_core in MAD units, 2010 ground truth with the 90 accepted 2010 offsets (`e3_bostock_family_offsets.csv`). Also the reverse (2010 templates → 2012) using existing E4/E5 2010 stacks. Accept: 2010 own recall at frozen T within 0.7–1.3× the 7-station 2012 ablation value; 2010 quiet-day FA (doy 222–226, 250–256) within 0.5–2× the 2012 quiet-day FA. Failure: recall <0.5× prediction → check time base/moveout, then placement; if it persists the method cannot be frozen across years → per-SSE re-stacking with per-SSE calibration.

**Step 3 — A2 template iteration protocol.** Selection: detections at T_core on tremor days (PNSN ≥30/day in box; Bostock ETS days in referee years), gated ≥2 sites each ≥2.5 MAD of their own site series (`--coinc-gate` path), within-family declustered at template length, cross-family assigned. For any site carrying >40 % of a family's weight (PGC), select with that site EXCLUDED from the sum (leave-one-site-out). Stacking: reuse `e5r_templates.py` cut_windows/align_stack (rms-normalised windows, linear stack, cc alignment ±1.5 s, 2 iterations, one shift per station for all three components); optional parabolic sub-sample refinement; never per-component shifts, never cc-weighting, never cc-thresholding of events. ≤2 re-stacks; stop when cc(iter k, k−1) ≥0.98 over S−1..S+7 at every station. Reject if ANY: (a) q-control cc to SeisBench stack falls >0.05 or below 0.6; (b) held-out-year own recall at matched FA drops >0.01 (adopt only if ≥+0.01; test 2012→2010, 2010→2012, then 2011/2013; Bostock also covers 2008/2009); (c) S pick moves >0.5 s or station-differential moveout changes >0.3 s; (d) cross-family collapse: max cc to any OTHER family's template rises >0.1 or exceeds 0.7 (116×116 per station); (e) spectral narrowing: >50 % of 1–8 Hz energy in one 1-Hz bin; (f) stack S-SNR falls >1 dB. Expected ≤+0.02 recall; do not chase more.

**Step 4 — A4 epoch transfer and identity test (start on disk with 2017 VGZ).** For each (station, epoch boundary): select events in the new epoch's first SSE with the network EXCLUDING that station at T_core; stack the station; compare with previous-epoch template over S−1..S+7. PGC 2017-08-10 has a natural control (B009/B010/B011 unchanged). VGZ 2017-03-03 sits inside the cached 2017 window (SSE 2017-03-18..04-03) — test now. Lost stations: no re-stack; completeness via Step 6 ablation. Accept: |cc| ≥0.6 (cc ≤ −0.6 → polarity flip, invert); S pick within ±0.5 s; registration ≥0.8× previous epoch. Failure: drop station for that family-epoch; if >50 % of families fail at a station → audit StationXML. Track B: re-order queue so one SSE window per epoch arrives first — 2015-12-05..2016-01-15, 2018-06-10..07-25, 2021-10-01..11-05, plus one winter quiet month (2016-12) — before oldest-first fill.

**Step 5 — A3 relocation with controls.** PGC site = one station. Station corrections for LRIV and B001–B007 from PNSN intraslab earthquakes (M≥2, 25–50 km depth, inside the family box) — not from Bostock. Crosson 1-D (`d3_stage4_locate.py`), S on horizontals + P on Z; Run A free depth, Run B Slab2 prior ±5 km. Controls: (i) synthetic floor (Bostock hypocentres forward-modelled, 0.3 s S / 0.2 s P noise, ±0.5 s correction errors); (ii) shuffled-pick null; (iii) spread vs Bostock's; (iv) fraction at depth bounds. Accept (families with ≥4 sites within 60 km, gap ≤180°): horizontal median ≤5 km AND observed median <0.5× shuffled-null; spread ≥0.7× Bostock's; ≤5 % at bounds. Depth: report median |Δz| (expect 5–10 km). Failure: within 1.3× of null or spread <0.5× → centroid collapse; keep Bostock's.

**Step 6 — A5 completeness model (three components).** (1) Epoch/station-set ablation of the 2012 test half (and 2010) → recall per family per epoch. (2) Daily sensitivity: store each station's daily median 8-s window energy; calibrate per-station excess vs Bostock magnitude (E4 slope 0.38); per-family-day Mc at T_ext and T_core. (3) Floor λ_f(T) per family per epoch per season from quiet days (no PNSN tremor in box, ≥3 days from a tremor day). Every event carries net_mad, n_sites≥2.5, rel_amp, Mc(f,d), λ_f(T,d). Accept: predicted vs measured own recall on held-out 2010 within ±0.05 per magnitude bin; winter quiet-day FA at T_ext within 2× summer.

**Step 7 — Track C production rules.** ≥2 sites; per family-day MAD; two tiers (T_core everywhere; T_ext only inside tremor windows or with λ attached); frozen cross-family assignment (share ≥50 % stations → raw weighted cc, else max MAD); EQ flag; HF veto only if adopted; per-station stacks at T_core detections saved monthly; monthly QC: site quiet-day rate >2× its λ → flag.

**Step 8 — Track D re-scope.** Bars: main-front first-burst speed within 20 % of Bostock on 2012 (met); re-rupture interval median within 5 %. Reported, never bars: secondary-front count and re-rupture family fraction vs per-family-thinned Bostock (`e6_thinned_reference.py`); tides — no claim without an interface tidal-stress model.

## Bars: verdicts
- ETS-day site-rolled ≤5 % — WRONG/unfalsifiable; replace with quiet-day site-rolled/data ∈ [0.77, 1.3] + reported one-to-one precision.
- Quiet-day FA ≤1/template/day — keep as T_ext; ADD T_core ≤0.1 and post-SSE window.
- Recall ≥30 % — unsupported; restate as within 0.7–1.3× of the station-set ablation prediction; honest expectation 0.18–0.24 at FA 1.0 on 15 stations, lower post-2021.
- Location ≤5 km — split; controls mandatory. EQ veto ≤0.5 % — report.
- MISSING: cross-year frozen-transfer bar; epoch identity bar; completeness-model validation bar; seasonal FA stability bar.

## Experiments Merlin could not run
1. mergeB + 12 s at matched FA. 2. HF-ratio veto. 3. Frozen 2012→2010 transfer (decisive). 4. Leave-one-site-out selection counts for PGC templates (<50/family/SSE → fallback). 5. Synthetic relocation floor on our geometry. 6. Winter quiet-day FA.

## What would change Merlin's mind
Step 2 recall ≥0.9× prediction AND Step 3 adding ≥+0.05; winter FA ≤1.5× summer and HF veto cutting FA ≥50 %; stack test showing quiet-day detections register at independent sites.

Key files: `lfe_catalog_repair/scripts/e4_template_test.py` (577–579, 660–670, 771–773); `e5r_score.py` (162); `results/e6_fa_table_e6_*.csv`; `results/e5r_template_qc_2012.csv`; `codex-LFE/results/target_network_from_map/channel_triplet_epochs.csv`.
