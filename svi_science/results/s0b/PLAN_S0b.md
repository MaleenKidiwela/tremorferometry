# S0b — episode consistency with cross-family de-duplication: PLAN (frozen before any de-duplicated count)

Written 2026-10-07 ~8:55 am PT, before any de-duplicated count or consistency number was computed.
Parent result: `lfe_extension_wa/results/episode_rates/RESULTS.md` and `robustness/ROBUSTNESS.md` (first run, untouched).
Script: `svi_science/scripts/s0_s0b_dedup.py`. It imports `lfe_extension_wa/scripts/episode_rates.py` as a module with
`EPR_OUT` pointed at a folder under `svi_science/results/s0b/`; `episode_rates.py` itself is not edited. Only the per
(family, date) count files are rebuilt; episodes, exposure, windows, tests, neighbourhoods, bootstrap (2,000) and
permutations (1,000), seeds and bins are the first run's, unchanged.

## Why
The first run counted every candidate at P ≥ cut, including candidates that a neighbouring family also detected (the same
event seen by two or more families 0.1–1.6 km apart). Such double counts make neighbours' residuals move together by
construction, which could produce part of "residuals of families < 3 km apart correlate 0.73 (SVI)".

## De-duplication rules
- **SVI**: keep a candidate only if `keep_xfam8 == 1` in `lfe_extension/results/refine/conf/P_by_candidate.parquet`.
  Definition (to be verified in code before use, step 1): per run, day and series, a candidate is dropped iff a candidate of
  a station-sharing family (|A ∩ B| / min(|A|, |B|) ≥ 0.5 over the 2012 usable-template station sets) lies within the
  frozen window and has a higher raw network cc (ties: earlier wins). Applied to data and scrambled control the same way
  (check). The flag is applied before the P cut (it is a property of the candidate), so it acts on every cut alike.
- **Washington-NW**: the production WA candidates (`wa_prod/conf/P/prod/P_*.parquet`, net = core) carry no such flag, and NW
  families have no position-group merge (only Olympic A families are merged at ±2 s; `pos_group` is empty for all 132 NW
  families). Rule, frozen here, the same as SVI's: per day and series, among all families present in the core-net rows
  (any region), drop a candidate iff a candidate of a station-sharing family (same 0.5 overlap rule over `tpl_library`
  station sets ∩ that era's core set) lies within ±2 s with higher raw `net_cc` (ties: earlier wins). Scrambled (c1) and
  time-reversed (rev) controls get the same rule within their own series.
  - Duplicate rate reported at ±2, ±4 and ±8 s, for P_own_WAVE_PNET ≥ 0.75 data candidates of NW families, era E3 and all
    eras; also with the pool restricted to NW families only (column).
  - **Materiality bar**: if the ±2 s rule removes > 5 % of NW E3 data repeats at P ≥ 0.75, the NW E3 consistency is re-run
    with de-duplicated counts. Otherwise NW is reported as not materially affected (no re-run).

## Re-run
- Same frozen episode list (`episodes_svi_final.csv`, `episodes_nw_final.csv`), same windows, same rules.
- Counts: SVI = `count_svi` logic on `keep_xfam8 == 1` rows; NW (if material) = `count_wa` logic on kept rows.
  Acceptance of the re-implementation: with de-duplication switched off, the rebuilt count files must equal the first run's
  `counts_svi_fid_date.parquet` / `counts_wa_nw_fid_date_half.parquet` exactly.
- Consistency run with `EPR_ROBUST=1`, so the headline (P ≥ 0.75, variant `main`) and the 0.5 column (variant `cut_050`) come
  from the same run. Primary block: SVI fixed-S13 (as before); NW E3 primary for Washington.
- Before / after table: tested, over-dispersed (BH q < 0.05), consistent, median χ²/df, residual correlation by distance bin
  (< 3, 3–6, 6–10, 10–20, 20–40 km), distance and carrier coefficients with permutation p. "Before" = the robustness run V0
  (`robustness/`, identical to the first run) and V1 (P ≥ 0.5).

## Pre-registered criterion
"Neighbours move together" **survives de-duplication** in a block if, at P ≥ 0.75 de-duplicated:
(i) the distance coefficient keeps its sign (negative: correlation falls with distance) with permutation p < 0.01, AND
(ii) the mean residual correlation stays positive in both the 3–6 km and 6–10 km bins (not only < 3 km).
Evaluated for SVI fixed-S13 (headline) and, if re-run, NW E3. The 0.5 column is reported, not part of the criterion.
If (i) holds but (ii) fails, the verdict is "only co-located (< 3 km) coupling survives". If (i) fails: "does not survive".
Reported alongside, not part of the criterion: change in the < 3 km bin (the expected direct effect of de-duplication),
over-dispersed counts, the scrambled-control duplicate rate (data vs control).

---------------------------------------------------------------------------------------------------------------------------
# Results (added after the runs, 2026-10-07 ~9:05 am PT)

Files: `before_after.csv`, `verdict.json`, `common_pairs_corr_by_bin.csv`, `keep_xfam8_check.{csv,json}`,
`nw/nw_duplicates.json`, `nw/nw_duplicate_rates.csv`, `coincidence_diagnostic.{csv,json}`, `prep_checks.json`, full
episode_rates outputs in `run_dedup/`. "Before" = `lfe_extension_wa/results/episode_rates/robustness/` (V0 / V1).

## Step 1 — keep_xfam8 definition (verified in code and data)
- Code: `lfe_extension/scripts/svi_ref_mf.py:474-477` → `c_common.local_max_keep_cc(t, net_cc, fam, share_ok_matrix(), w=2.0)`,
  per run × day × series. The window is **±2 s, not 8 s** ("8" is Merlin's E5r-audit finding F8). Families must share ≥ 50 %
  of their 2012 template stations; the higher raw network cc wins. Computed on all candidates before any P cut.
- Applied to the scrambled control the same way (per series). Data: every one of the 156,050 (prod) / 105,811 (prod_S13)
  keep = 0 rows has a station-sharing neighbour candidate within ±2 s; 0 without. Duplicate fraction at P ≥ 0.75: data
  23 % (3,088 / 13,466), scrambled control 22 % (18 / 83).
- The re-implemented raw counts equal the first run's count files exactly (SVI and NW, 0 differing rows).

## Step 3 — Washington-NW has the same double counting (material)
- NW families carry no cross-family flag and no merge (pos_group empty; only Olympic A is merged). Frozen ±2 s rule applied:
  **31 % of NW E3 data repeats at P ≥ 0.75 are duplicates** (5,723 / 18,525; ±4 s 38 %, ±8 s 57 %; NW-only pool 29 %).
  E1 31 %, E2 31 %. Scrambled control 24 %. Material (> 5 %) → NW E3 re-run with de-duplicated counts.
- Winners: 92 % NW, 8 % Olympic A families; winner distance median 19.5 km (only 963 of 5,723 within 3 km).
- Diagnostic (not pre-registered), time-shifted neighbours: share of P ≥ 0.75 repeats with a station-sharing other-family
  candidate within ±2 s — SVI 61 % true vs 5–8 % with the pool shifted ±30–120 s (genuine simultaneous detections);
  NW E3 76 % true vs 22–27 % shifted (about a third of NW "duplicates" are random coincidences on busy days). The scrambled
  control shows the same true-vs-shifted gap (57 % vs 1–10 % SVI): co-detection is a detector property, in data and control.
  SVI duplicate neighbours are 0–3 km (2,765), 3–6 (1,879), 6–10 (1,549), 10–20 (1,021), > 20 km (927); median 5.9 km.

## Step 2/4 — before vs after (P ≥ 0.75 headline; 0.5 in `before_after.csv`)
| block | tested | over-disp. | consistent | median χ²/df | < 3 | 3–6 | 6–10 | 10–20 | 20–40 km | distance coef (p) | carrier coef (p) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SVI fixed-S13 before | 75 | 27 | 48 | 1.81 | 0.73 | 0.49 | 0.39 | 0.19 | −0.05 | −0.517 (0.001) | −0.035 (0.48) |
| SVI fixed-S13 after | 68 | 24 | 44 | 1.70 | 0.54 | 0.47 | 0.37 | 0.19 | −0.07 | −0.463 (0.001) | −0.027 (0.50) |
| NW E3 before | 89 | 58 | 31 | 4.03 | 0.31 | 0.11 | 0.16 | 0.06 | −0.03 | −0.111 (0.001) | +0.230 (0.001) |
| NW E3 after | 73 | 50 | 23 | 3.67 | 0.31 | 0.13 | 0.17 | 0.06 | −0.04 | −0.116 (0.001) | +0.269 (0.001) |
At P ≥ 0.5: SVI 54 → 44 of 93 → 88 over-dispersed, median χ²/df 3.05 → 2.42, distance −0.497 → −0.451 (p 0.001);
NW E3 93 → 77 of 112 → 104, 6.43 → 4.86, −0.140 → −0.135 (p 0.001).
Same family pairs only (`common_pairs_corr_by_bin.csv`), SVI 0.75: < 3 km 0.73 → 0.54, 3–6 0.53 → 0.47, 6–10 0.41 → 0.37,
10–20 0.19 → 0.19.

## Verdict (pre-registered criterion): SURVIVES in both SVI fixed-S13 and NW E3
Distance coefficient negative with p = 0.001 and both 3–6 and 6–10 km bins positive after de-duplication. De-duplication
removes about a quarter of the < 3 km SVI correlation (0.73 → 0.54) and little beyond 3 km, so co-detection inflated the
closest pairs but does not make the "neighbours move together within ~20 km" pattern. Fewer families reach sum N ≥ 20
(SVI 75 → 68, NW 89 → 73). The NW carrier term stays significant (+0.27).
Deviation: the ±2 s rule also removes some random coincidences (strongly in NW), i.e. it thins busy days slightly more
than quiet ones; this works against over-dispersion, so the surviving result is conservative on that point.
