# Washington GO/NO-GO on the 21 core stations (Merlin's review, steps 1–2)

Written 2026-10-02 (~20:45Z). Everything here is under `lfe_extension_wa/results/wa_gonogo/`. Scripts: `scripts/gonogo_truth_frac.py`
(Step A), `scripts/gonogo_run.py` (setup, completeness checks, matched filter and HF stage through the unchanged frozen
`wa_mf.py` / `wa_hf.py`), `scripts/gonogo_score.py` (the T1 scorer `merlin_t1_score.py` with the changes listed in its
docstring), `scripts/gonogo_report.py` (tables). Logs: `logs/gonogo_*.log`. No downloads, no deletions, no git, nothing
existing overwritten. "Scrambled control" = the same detector with each station site's data shifted to a different
time of the same day, so nothing lines up between sites except by chance.

## In short

- **Verdict: GO by the fixed rule, with thin margins.** Counted per family-episode window, as T1 counted (instruction:
  "family-episode windows as in T1"), NW has **118** windows ≥ 2× above BOTH scrambled controls (V2010 81 + A2011 37;
  rule ≥ 100) and NW own-family precision **0.226** (rule ≥ 0.2) on the 21 core stations
  (`score_core21/summary.csv`, `comparison.csv`).
- How thin: (1) counted as **distinct families it is 81** (every A2011 passer also passes in V2010), which would fall
  in the 50–100 "borderline" band; (2) the scrambled control tested against the other scrambled control passes 13
  NW family-windows, so ~105 of the 118 are above what chance gives; (3) A2011 alone has own precision **0.197**; (4)
  part of the truth was made with the same permanent stations (Step A), see "Caveats".
- Same numbers on the original 42-station T1 run (scored here the same way, V2010 test half + A2011): 162 family-windows
  (115 distinct families), own precision 0.237. So the 21 core keep **73 % of the passing family-windows (118/162) and
  70 % of the passing families (81/115)**, at nearly the same precision. All 81 core passers also pass on 42 stations.
- **Which regions survive on the 21 core:** NW 81 families (list in `per_family.csv`, `pass_any_core21`); Arm B 37 of
  41 pass the ratio test but precision fell 0.246 → 0.118 and recall is 0.006 (weak, as before); Olympic Arm A 4 of 28
  (unchanged; data ≈ control: 241.9 vs 240.8 per family-day); **SW has no truth in V2010 or A2011** (SW truth is
  2007–2008 only), so it is not in the acceptance test; a V2008 supplement gives SW 73/73 families ≥ 2× both controls on
  the core, precision 0.837 → 0.605, recall 0.010 → 0.003, and **0 SW families have 2 clear core sites**.
- Step A reproduces Merlin's truth-independence counts (within one station at the median; exact for 10 of 11 months on
  his family set). 2011-04 truth has **0 array stations** (made only by permanent stations); 2010-08 (V2010) is ~60 %
  array, 2011-08 (A2011) ~50 %.

## Step A — how independent is the array truth of the permanent stations?

Per family-month: of the family's growth-template stations (`tpl_it<N>_stations.csv` next to its `det_it<N>.parquet`;
`NW_lag15/grow/tpl_it14_stations.csv`, `SW_lag15/grow/tpl_it7_stations.csv` and the other source runs of the WA
representatives), how many are LIVE that month (≥ 10 non-empty raw day files on disk), split into array (temporary
networks, `wa_common.TEMP_NETS`) and permanent. `truth_array_frac` = array / live. File: `truth_array_frac.csv`
(1,583 family-months with ≥ 1 own truth detection; every station code resolved; 0 % missing), month summary `truth_array_frac_months.csv`.
`lfe_array_discovery/traces_sealed` (that project's sealed holdout) is excluded; including it changed no median.

| month | WA NW families (median array / permanent, array frac) | WA SW families | Merlin's set: all 280 grown NW_lag15 families | Merlin |
|---|---|---|---|---|
| 2007-01 | 5 / 3 (0.63) | 8 / 1 (0.89) | 5 / 3 | 6 / 3 |
| 2007-10 | 4 / 3 (0.57) | – | 4 / 3 | 4 / 3 |
| 2008-05 | 6 / 4 (0.60) | 9 / 1 (0.90) | 6 / 3 | 6 / 3 |
| 2009-02 | 3 / 4 (0.43) | – | 3 / 3 | 3 / 3 |
| 2009-04 | 3 / 4 (0.43) | – | 3 / 3 | 3 / 3 |
| 2009-05 | 5 / 4 (0.56) | – | 4 / 3 | 4 / 3 |
| 2010-03 | 4 / 4 (0.50) | – | 3.5 / 4 | 3 / 4 |
| 2010-08 | 6 / 4 (0.60) | – | 5 / 4 | 5 / 4 |
| 2011-04 | **0 / 4 (0.00)** | – | 0 / 4 | 0 / 4 |
| 2011-08 | 4 / 4 (0.50) | – | 4 / 4 | 4 / 4 |
| 2011-12 | 4 / 4 (0.50) | – | 3.5 / 4 | 3 / 4 |

- Why the WA-representative medians differ by one: Merlin's numbers are reproduced on the full grown family set of
  NW_lag15 it14 (all 280 cids), not on the 132 WA representatives (which come from 5 NW source runs with their own
  template-station lists). The only remaining miss is 2007-01 array 5 vs 6: XU.W030 has 5 raw days that
  month (below the ≥ 10 rule); counting it gives 6. The .5 values are ties between 3 and 4. Not material: the pattern
  (array-majority 2007–2008, about half 2009–2011, none in 2011-04) is the same.
- The permanent stations that carry the NW truth are PB.B013, B943, B001, UW.GNW, SQM, PB.B007, UW.DOSE, LRIV; all
  except SQM are in the 21 core. So precision against truth is partly self-referential, most of all in A2011.

## Step B — T1 scorer on the 21 core stations

Matched filter: unchanged frozen chain (`wa_mf.py --all-templates --drop-dead 0.05 --n-control 2`, build-half
templates `results/templates/tpl_build.npz`, every ≥ 30-event template weighted by its SNR weight) on the existing
20 Hz trial cache, network = the 21 core stations. **PB.B926 is not in the 42-station trial roster and has no
build-half template, so 20 stations / 13 sites actually contribute** (`setup/core.json`). Days: V2010 (39 days, test
half scored) + A2011 (29 days, all scored), GPU 4.6 min. HF cuts calibrated on the V2010 build half (V2008 is not in
this test), the same for both runs, so the 42-station numbers here differ slightly from MERLIN_TESTS_RESULTS.md where
V2008 was in the calibration (V2010 NW own precision 0.372 and A2011 0.207 reproduce it). Per-epoch tiers (not used
in the headline): core 23.0 MAD (21.25 on 42 stations), `score_*/tiers_by_epoch.csv`. Headline tier = frozen core
8.75 MAD, catalogue flags, one-to-one ±2 s. "≥ 2× both controls" = (data + 1) / (max(control 1, control 2) + 1) ≥ 2 on
the family's label days (PNSN tremor days of its region, observable), raw tier counts, as T1 (T1 used control 1 only;
that count is in `summary.csv` as `sep2_core_data`: NW 133 on core, 168 on 42). "Control-vs-control" = control 2
against control 1, the false-pass level.

Data completeness: on-disk core station-days 1,244 of 1,428 on the 68 days, 0 missing or unreadable in the cache
(`guard_station_days.csv`); after the run, 1,176 cached usable station-days, 1,156 loaded, 20 dropped by the T1
rule (> 5 % of the day guarded; the same 20 as the 42-station run), 0 lost for any other reason
(`mf_core21/postguard.csv`); HF stage 527/527 station-day jobs returned ratios. (The first completeness check counted
116 extra station-days that exist only in `lfe_array_discovery/traces_sealed`, the sealed holdout; they were excluded,
`guard_station_days_with_sealed_root.csv`.)

| test | region | window | network | obs fam | fam-windows ≥ 2× both controls (distinct families) | control-vs-control passes | own prec | ≤3-km prec | own prec, scrambled control | data / control per fam-day | own recall |
|---|---|---|---|---|---|---|---|---|---|---|---|
| main | NW | V2010 | 42-station T1 (original) | 132 | 115 (115) | 10 | 0.372 | 0.374 | 0.009 | 8.35 / 3.91 | 0.104 |
| main | NW | A2011 | 42-station T1 (original) | 132 | 47 (47) | 10 | 0.207 | 0.208 | 0.003 | 8.97 / 7.84 | 0.131 |
| main | NW | both | 42-station T1 (original) | 132 | 162 (115) | 20 | 0.237 | 0.238 | 0.004 | 8.72 / 6.24 | 0.122 |
| main | SW | V2010 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| main | SW | A2011 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| main | SW | both | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| main | Olympic Arm A | V2010 | 42-station T1 (original) | 28 | 4 (4) | 0 | 0.053 | 0.090 | 0.001 | 209.45 / 205.63 | 0.245 |
| main | Olympic Arm A | A2011 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| main | Olympic Arm A | both | 42-station T1 (original) | 28 | 4 (4) | 0 | 0.053 | 0.090 | 0.001 | 209.45 / 205.63 | 0.245 |
| main | Arm B | V2010 | 42-station T1 (original) | 41 | 37 (37) | 3 | 0.246 | 0.246 | 0.002 | 3.13 / 0.92 | 0.003 |
| main | Arm B | A2011 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| main | Arm B | both | 42-station T1 (original) | 41 | 37 (37) | 3 | 0.246 | 0.246 | 0.002 | 3.13 / 0.92 | 0.003 |
| main | NW | V2010 | 21 core (20 usable) | 132 | 81 (81) | 4 | 0.318 | 0.319 | 0.009 | 10.13 / 6.56 | 0.124 |
| main | NW | A2011 | 21 core (20 usable) | 132 | 37 (37) | 9 | 0.197 | 0.198 | 0.003 | 9.40 / 8.39 | 0.128 |
| main | NW | both | 21 core (20 usable) | 132 | 118 (81) | 13 | 0.226 | 0.227 | 0.004 | 9.70 / 7.64 | 0.126 |
| main | SW | V2010 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| main | SW | A2011 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| main | SW | both | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| main | Olympic Arm A | V2010 | 21 core (20 usable) | 28 | 4 (4) | 0 | 0.048 | 0.080 | 0.001 | 241.87 / 240.83 | 0.249 |
| main | Olympic Arm A | A2011 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| main | Olympic Arm A | both | 21 core (20 usable) | 28 | 4 (4) | 0 | 0.048 | 0.080 | 0.001 | 241.87 / 240.83 | 0.249 |
| main | Arm B | V2010 | 21 core (20 usable) | 41 | 37 (37) | 7 | 0.118 | 0.121 | 0.004 | 5.09 / 3.08 | 0.006 |
| main | Arm B | A2011 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| main | Arm B | both | 21 core (20 usable) | 41 | 37 (37) | 7 | 0.118 | 0.121 | 0.004 | 5.09 / 3.08 | 0.006 |
| V2008 supplement | NW | V2008 | 42-station T1 (original) | 132 | 112 (112) | 11 | 0.339 | 0.340 | 0.005 | 8.08 / 4.35 | 0.109 |
| V2008 supplement | SW | V2008 | 42-station T1 (original) | 73 | 73 (73) | 12 | 0.837 | 0.848 | 0.000 | 3.33 / 0.01 | 0.010 |
| V2008 supplement | Olympic Arm A | V2008 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| V2008 supplement | Arm B | V2008 | 42-station T1 (original) | – | no family-window (no truth in this window) | | | | | | |
| V2008 supplement | NW | V2008 | 21 core (20 usable) | 132 | 86 (86) | 8 | 0.287 | 0.288 | 0.004 | 8.96 / 6.38 | 0.124 |
| V2008 supplement | SW | V2008 | 21 core (20 usable) | 73 | 73 (73) | 16 | 0.605 | 0.605 | 0.000 | 2.10 / 0.03 | 0.003 |
| V2008 supplement | Olympic Arm A | V2008 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |
| V2008 supplement | Arm B | V2008 | 21 core (20 usable) | – | no family-window (no truth in this window) | | | | | | |

Per region, template coverage on the 21 core (all 329 families, build-half templates `results/templates/tpl_qc_build.csv`; sites merged at 3 km):

| region | families | any ≥ 30-event template on core | ≥ 2 core sites with a ≥ 30-event template | ≥ 1 clear core site | ≥ 2 clear core sites | pass (≥ 2× both controls, V2010 or A2011) on core | of these with ≥ 2 clear sites | pass on 42-station T1 | pass on both |
|---|---|---|---|---|---|---|---|---|---|
| NW | 132 | 132 | 132 | 108 | 42 | 81 | 25 | 115 | 81 |
| SW | 79 | 73 | 73 | 14 | 0 | 0 | 0 | 0 | 0 |
| Olympic Arm A | 75 | 28 | 28 | 21 | 1 | 4 | 1 | 4 | 4 |
| Arm B | 43 | 41 | 41 | 7 | 1 | 37 | 1 | 37 | 37 |

NW own precision by the family-month's truth_array_frac (Step A; V2010 -> 2010-08, A2011 -> 2011-08):

| window | truth_array_frac | family-windows | own prec 42-station | own prec core | pass 42-station | pass core |
|---|---|---|---|---|---|---|
| A2011 | 0.45-0.55 | 66 | 0.218 | 0.210 | 24 | 16 |
| A2011 | 0.55-0.65 | 14 | 0.072 | 0.071 | 8 | 7 |
| A2011 | <= 0.45 | 48 | 0.411 | 0.385 | 14 | 12 |
| A2011 | > 0.65 | 4 | 0.400 | 0.400 | 1 | 2 |
| V2010 | 0.45-0.55 | 18 | 0.722 | 0.543 | 17 | 11 |
| V2010 | 0.55-0.65 | 74 | 0.362 | 0.337 | 63 | 43 |
| V2010 | <= 0.45 | 6 | 0.716 | 0.493 | 6 | 5 |
| V2010 | > 0.65 | 34 | 0.307 | 0.227 | 29 | 22 |

Columns: obs fam = families with ≥ 1 observable scanned episode day; own / ≤3-km precision = catalogue detections at
8.75 MAD matched to own / any ≤ 3-km family's array detection; the scrambled control column is its chance level;
data / control per fam-day at 8.75 MAD on observable episode days. Per family-window, all columns:
`score_core21/by_family.csv`, `score_orig42/by_family.csv`; per family side by side (core templates, pass flags,
label-day counts for data and both controls, precision): `per_family.csv`.

V2008 supplement (not part of the acceptance; only way to see SW): `score_core21_v2008/`, `score_orig42_v2008/`. Its
completeness check fails the literal 3 % rule (27 of 539 = 5.0 %), but 26 of the 27 are PB.B926, which has no
template and contributes nothing; on the 20 usable stations 1 of 513 (0.19 %, PB.B018 2008-05)
(`SUPPLEMENT_V2008_GUARD_FAILED`). Treat SW numbers as indicative.

## Caveats (read before acting)

- **Self-referential truth.** NW own precision rises where the truth leaned on permanent stations (table above:
  A2011 truth_array_frac ≤ 0.45: 0.385 on core; V2010 > 0.65: 0.227). This pattern fits partial circularity but is
  confounded (families near permanent stations are also easier to detect). The most array-made subset (V2010, > 0.65,
  34 families) still has 0.227 own precision and 22 passers on the core.
- **What the core loses**: 34 of the 115 NW families that pass on 42 stations. The 42-station run had, beyond the
  core, UW.BLN, SQM, SMW, STW, US.NLWA on all 68 days and GMW, OSD, SVOH, SP2, WISH and others in V2010
  (`results/mf_mt1_all/run_*.json`); BLN carries 65 NW ≥ 30-event templates, SQM 53, GMW 33
  (`results/wa_prod_plan/inventory_by_station.csv`, NW30). These are the tier-B download stations.
- **Recall** is low in both networks (NW own 0.126 on core, 0.122 on 42): a pass here means a family's repeats can be
  told from the scrambled control on episode days, not that most repeats are found.
- **SW cannot be judged** in this test (no truth in V2010/A2011); the supplement suggests few, clean single repeats
  but 0 families with 2 clear core sites.
- **Olympic Arm A** stays at data ≈ control (241.9 vs 240.8 per family-day): still not separable by single repeats.
- What would change the verdict: counting distinct families (81) instead of family-windows makes it borderline; a
  stricter control rule (subtracting the 13 control-vs-control passes) gives ~105, still ≥ 100; A2011 alone is just
  below the precision bar (0.197).
