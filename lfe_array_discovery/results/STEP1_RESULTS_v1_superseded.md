# Step 1a — north PERMANENT bed (6 sites, 46–80 km): the autocorrelation stage has no SNR margin

> **SUPERSEDED NUMBERS (added 2026-09-20 15:xxZ by the orchestrator).** The run described below used 9 of 12 stations / 5 of 6 sites because of the hour-cut t0 bug (cache midnight offsets +0.0084 s PB / +0.0150 s LRIV; fixed 2026-09-20). The corrected re-run (`w_step1a_refresh.sh`, summary `results/STEP1_DONE`, 15:06Z) on the full 12 stations / 6 sites gives the same verdict: T = 7.25 (FA rule), 563 tremor hours → 199,579 pairs ≥5.0 → 3 seed clusters → 0 pass-0 survivors → 0/10 isolated hidden families (ceiling 5/10); D2 true-pair z 2.13 vs null 1.98. Diagnostics D1–D3 below stand. Full Step 1a/1b/1c write-up follows when 1c completes.

### *uninformative about the engine; informative about distance*

Run 2026-09-20, agent E. Frozen before launch in `results/freeze_step1.json` (sha256 of 12 `w_`
scripts + 13 imported lfe_picker scripts; nothing in lfe_picker edited). Engine per PLAN §10 /
Merlin Q1(a). **No constant was retuned during or after the run.**

## 1. What was run

`lfe_catalog_repair/cache20hz`, 2012-08-15..10-05, 52 day-files, **618 station-days**,
12 detection stations merged by the data-driven 1-km rule into **6 sites**; VGZ/SNB/B926 withheld
(Gates suite). Chain: network autocorrelation → clusters → align_stack handoff → matched-filter
growth → survival + template-scramble null → 2 iterations → dedupe → locate → Gates 0–2.
**The three site-shift null realisations were NOT run** (coordinator decision: the outcome is
already determined at the autocorrelation stage, and a null on a chain that produces nothing
measures nothing).

| site | members | median hypocentral distance to the 130 Bostock families |
|---|---|---|
| PGCq | B009+B010+B011+PGC | 52.6 km (40–76) |
| LRIVq | B005+B006+B007+LRIV | 64.9 km (33–117) |
| B003 | B003 | 70.1 km (40–108) |
| B004 | B004 | 74.2 km (46–111) |
| B001 | B001 | 77.4 km (44–133) |
| NLLB | NLLB | 91.4 km (55–137) |

## 2. Engine prerequisites (Merlin F3) — all satisfied, and the FA rule BINDS

| prerequisite | state |
|---|---|
| P1 no `T = 8.0` floor | **T = 7.75 derived by the FA rule**, engineering floor 5.0 never binding; growth-stage `T_seed` 6.25 / 7.75 with `floor_binds = False` for both templates |
| P2 data-driven `SITE_MEMBERS` | 1-km single linkage reproduces the hand-written merge exactly (6 sites) |
| P3 `station_z` over valid samples only | implemented in the `w_` subclass |
| P4 `MIN_SITES` 4, availability never intersected | 6/6 sites usable on all 52 days |

**T derivation** (45 quiet hours on 9 quiet days — the bed has only 9 days ≥3 days from any PNSN
tremor day within 100 km, so Merlin's "≥20 quiet days × 2 h" was met as "≥40 quiet hours" with 5
hours per quiet day; deviation recorded in the freeze):

| T | pairs | false pairs / quiet hour | quiet components ≥5 |
|---|---|---|---|
| 6.50 | 853 | 18.96 | 6 |
| 6.75 | 436 | 9.69 | 2 |
| 7.00 | 216 | 4.80 | 0 |
| 7.50 | 61 | 1.36 | 0 |
| **7.75** | **31** | **0.69** | **0** |

## 3. Result

| stage | count |
|---|---|
| candidate tremor hours (≥5 PNSN epicentres ≤50 km of an active cell) | 563 over 32 days |
| autocorrelation, median cost | 0.81 s / station-hour (H100, 3.8 GB), 11 min total |
| pairs recorded at net_mad ≥ 5.0 | 248 397 |
| pairs at T = 7.75 | 3 740 → **6.64 / tremor hour vs 0.69 / quiet hour (9.6× excess)** |
| connected components → **seed clusters** (≥5 in an hour, ≥2 ten-minute bins) | **2** |
| pass-0 templates with ≥4 sites | 2 |
| growth hits over 52 days | 16 (12 + 4), 9 + 4 independent |
| **survivors** (≥20 independent, ≥3 days, ≤20 % near seed, scramble ≥3×) | **0** |
| iterations / dedupe / locate / Gates 0–2 | 0 templates throughout |
| **isolated hidden families rediscovered** | **0 / 10** (ceiling: the frozen template detector itself 5/10) |
| spurious-split rate | 0.00 |

Verdict from `w_score_step1.py`: `ENGINE_INADEQUATE` **on this bed** — i.e. Merlin's
"<2/10 with null ≤1" failure signature. Accept bars were ≥3/10, split ≤0.1, null ≤1.

There *is* a real tremor-related excess of coherent pairs (9.6× more at T than in quiet hours), so
the statistic is not inert; those pairs simply do not organise into repeating families.

## 4. Why — three diagnostics (`results/step1_diag.json`, `w_diag.py`, KEPT families only)

**D1 — dynamic range.** A window correlated with itself (a perfect repeat) scores
**z_self = 14.4**, because the lower-half MAD of an hour's cc distribution is ≈0.07. So a pair of
occurrences needs **cc ≳ 0.35** to reach z = 5 at one station. That is the ceiling of the whole
statistic on this bed.

**D2 — true-pair z.** For the five largest kept families, at the TRUE lag of catalogued same-family
event pairs, cut moveout-aware at the family's own predicted S (the most favourable cut that
exists): **median z = 2.13** against a same-width jitter null of **1.98**; only **16 %** of true
pairs exceed their own null maximum. Measured at window lengths 3, 4, 6, 8, 10, 15 and 20 s the
site-median z is flat at 1.6–2.4 — window length is not the free parameter that rescues this.
Real LFE pairs here sit at cc ≈ 0.10–0.15, a factor 2.5–3.5 below what z = 5 requires.

**D3 — pair invariant.** Of the 2 000 strongest pairs the engine reported, **0** have a lag equal to
an inter-event interval of a kept family at a common offset; the lag-jitter null gives 0, 2, 0. The
pairs that pass T are not LFE pairs.

**Geometry.** The bed's S-moveout spread across the 6 sites is **median 13.8 s** (9.4–24.3 s;
93 % of families > 10 s). Merlin sized the 10-s common window for the **45-km-aperture WA array**
(5–6 s spread). On a 46–137 km permanent network a single common window cannot contain S at all six
sites for 93 % of the family field, so four of six sites contribute noise to every sum — an
additional √3 loss on top of the cc deficit.

## 5. What this does and does not establish

- It does **not** establish that the Arm-A engine cannot discover families. Every stage ran, the FA
  rule bound above the engineering floor, the handoff/growth/iterate/dedupe/locate/gates path was
  exercised end to end (smoke tests and this run), and the failure is located precisely at the
  first stage, for a measurable physical reason.
- It **does** establish that **single-pair network autocorrelation needs stations at LFE distances
  of order 10–40 km**, not 46–137 km. This is the P4-0b starvation lesson in a new place: the north
  PERMANENT network is the wrong bed for a positive control, because Bostock's families were
  discovered on temporary stations that this bed does not contain.
- Consequence: Step 1 is re-run as **Step 1b** on the same 2012 window with the POLARIS/CN dense
  stations added (`results/freeze_step1b.json`). Constants are carried over unchanged — the window,
  hop and T are NOT retuned before that result is seen.
- **For Merlin to audit**: (i) is the 10-s common window still the right choice once the bed has
  10–40 km stations (the WA arrays' own aperture is 40 km, spread ≈5–6 s, so probably yes); (ii)
  should the statistic fall back to a moveout-aware, cell-scanned pairing on wide-aperture beds;
  (iii) does the 9.6× tremor excess of coherent pairs carry usable information that the
  components/seed-cluster rule discards.

## 6. Defect found in the inherited engine

`v_engine.DayData` stores a single `self.NF` (the last station's rFFT length) and uses it for every
station in `station_z`. It is harmless only while every station-day has an identical sample count;
a short or partial day raises `RuntimeError: The size of tensor a ... must match ...`. Fixed in the
`w_` subclass (per-station `NFs`); **`v_engine.py` on disk still carries it** — for the P4-0b
close-out.

## 7. Files

`results/freeze_step1.json`, `results/step1_score.json`, `results/step1_diag.json`,
`results/step1_diag_D1_dynamic_range.csv`, `results/step1_diag_D2_true_pair_z.csv`,
`results/step1_geometry.json`, `results/step1/` (T.json, T_grid.csv, ac*.parquet, cl_seeds.csv,
bank*.npz, p*_survival.csv, dedupe.csv, lg_*), `logs/step1_queue.log`, `logs/step1a_tail.log`.
