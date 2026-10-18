# NOTES — lfe_extension_wa (dated ledger)

## 2026-09-29 — build, smoke test, launch
- Families: 329 representatives (NW 132, SW 79, Olympic A 75, Arm B 43). Stations: 42 permanent
  (37 anchors + BLN, SMW, SQM, OPC, NLWA), 35 sites (3-km rule = SVI merge-B clusters + {B013,B943}).
- Windows: V2010 (NW + Olympic validation), V2008 (SW validation — deviation: SW array detections exist
  only in 2007-01/2008-05, so no single validation episode covers all three regions), T2012 trial.
- Station-days: 4,914 planned, 4,667 metadata-live; 2,438 already on disk (reused), 2,229 fetched
  (~12 MB each, ~27 GB). CN from NRCan (IRIS returns no CN data for 2008).
- 50-sps stations (GNW, SQM, WISH BH50) are resampled 50→20 with resample_poly(2, 5) after the
  1–8 Hz filter at 50 Hz (the frozen e4_cache20 only handled integer factors).
- Smoke (6 days of 2010-08, build-half events of those days only): 8,714 (family, node) stacks,
  132 clear, 15 families with ≥ 2 clear sites (all NW; B013/B943 carry most). MF 3 s/day on 15
  families; matched detections sit at Δt = 0.00 s from the array detections (timing chain verified);
  one day: fam 20 45 T_core detections, 31 match own array detections; site-rolled control 12, 0 match.
  Smoke numbers are NOT the result (templates and test on overlapping days).
- Expectation to watch: the clear-arrival rule on permanent stations is strict — few families will
  reach 2 clear sites (SW and Olympic especially). That is a finding, not a bug; perm_weak stations are
  not used (rule = clear arrival only, per the design).
- Queue launched: `scripts/wa_trial_queue.sh` (pid in results/WA_TRIAL_PID).

## 2026-09-29 17:06Z — first stages done (fetch, cache, templates, validation MF/HF/score)
- Data: 2,134 station-days fetched (26.3 GB) + 2,438 reused; 4,575 cached (1 empty). 91 nodata (CPW 49 in 2012, DOSE 16, SP2 8).
- Templates (build half, V2008+V2010): 13,830 (family, node) stacks; **386 clear** on 20 stations (B943 104, B013 102, BLN 45, HDW 39,
  B001 35, B018 13, NLWA 11, GNW 7, PGC 6, VGZ 6 …); 180 families with ≥1 clear, **61 with ≥2 clear sites** (NW 55, SW 2, Oly A 2, Arm B 2).
- HF cuts: 8 HF nodes, 4 calibrated (≥30 own-matched build detections).
- Validation (test half, observable scanned days; 116 family-window rows observable, 16,917 array detections):
  frozen T_core 8.75: recall 0.158, precision (own) 0.250 [NW V2008 0.157/0.343, NW V2010 0.154/0.315]; scrambled control
  chance-match precision 0.007. BUT quiet-day FA at 8.75 = 8.5 (V2008 post days) / 12.5 (V2010) per family-day and the control COUNT ≈
  data count (12.2 vs 13.9 per family-day): the 2–3-site permanent networks are dominated by single-site bursts. SVI's per-epoch FA rule
  gives T_core_ep 20.75 / 17.75 → recall 0.003, precision 0.39. So the permanent-only matched filter FINDS the known events (25–35 %
  precision vs 0.7 % chance) but at the frozen threshold most detections are not the families' events. Provisional — see WA_TRIAL_RESULTS.md.

## 2026-09-29 19:30Z — step 4b (add stations, NW + Olympic) built, smoke-tested, launched
- Candidates: 705 stations scanned, 77 eligible; 52 selected (tier 1: 42 = UW 31, CN 6, C8 2, PB 1 (B926), CC 1 (OBSR),
  XD 1 (MB05); tier 2 TA/7D only where the only option: TA D03D D04D D04E A04D B05D, 7D FN01C FN05C FN07C FN14C
  FN18C). Older permanent non-anchor stations (OFR, OBC, ALKI, MEAN, TOLT, ELW, SNB, LZB, …) are included: the
  request lists "newer permanent stations" but also "every station with verified data"; they are the ones live in
  the array era. C8 PA01–12 (at PGC), TA *A anchors-twins, GM (SP2) merge into existing sites → excluded.
- Projection (empirical clear rate by distance): ≥ 3 clear sites 25 → 28 families (NW 24 → 27; Olympic 1 → 1;
  Arm B 0 → 0). Most families cannot reach 3 sites with 2010–2014 stations; the run measures it.
- Code changes: `wa_common` WA_SETUP env + windows A2011/T2013/T2014; `wa_cache` WA_EXTRA_NODES; `wa_fetch`
  --days, C8 → NRCan fallback, NRCan 404 = no data; `wa_mf` first-tag-wins + `--n-control 2`; `wa_templates`
  keeps all QC columns when every stack is < 30 events (smoke crash).
- Smoke (2010-08-20/21, 2011-08-15, 2013-08-31, 2013-09-20; NOT the result): fetch 220 station-days 2 GB in 55 s;
  cache nodes correct (SNB BH40→HH100 2011, LEBA BH→HH 2013, OBC/OCP 2013 re-epoch); 3-day array templates 9,697
  stacks / 7 clear (SQM EH100@2011 4, RMW 2, RCS 1); every stage ran through to the report.
- Queue `scripts/step4b_queue.sh` launched 19:29Z: 9,143 station-days to fetch.

## 2026-09-29 20:30Z — stack activity (family activity by stacking) built, smoke-tested, launched
- Engine `stack_activity.py`: all 11,358 ≥30-event build templates (386 clear), 274 families, leave-one-site-out top-K
  (2/scored h) candidates, 2–8 Hz stacks for data / 3 site-rolled controls / 6 random-times draws, stored per day.
- Bugs fixed in smoke: NpzFile re-decompressed the template array on every access (10 min stall); candidate concat
  index; per-station loop → vectorised gather; multinomial → uniform draws with rejection; ±4 s declustering by 4-s
  block maxima. Speed 42 → ~34 s/day on the shared GPU (42 stations).
- Smoke (2010-08-20/21, NOT a result): full-network candidates hit own array detections 4.6 % (control 0.0 %; up to 62 %
  for Olympic families); W = 1 AUC(active ≥ 5) A 0.67 vs control 0.53 (NW 0.74 vs 0.54).
- Queue `scripts/stack_activity_queue.sh` launched 20:29Z (pid in results/STACK_ACT_PID): valA now; waits for step 4b;
  then valB/valC, tr2012-14, continuous-series selection + fetch (smoke estimate ~2.5k station-days ≈ 30 GB, ~80 % of
  2012–2014 already on disk), valS/valSB, cont, report.

## 2026-09-30 04:00Z — Merlin's two deciding tests built, smoke-tested, launched
- T1 (all templates): `wa_mf.py --all-templates --drop-dead 0.05`, `wa_hf.py --all-templates`, scorer `merlin_t1_score.py`
  (3 matching levels, separation on label days, 2 controls), merged setup `results/merlin_tests/t1_setup` (A2011 truth from
  step 4b). Scorer reproduces WA_TRIAL on `mf_val` exactly (recall 0.158, prec own 0.250, any 0.439, control 0.007; epoch
  tiers 20.75 / 17.75). Smoke (2010-08-20, 2011-08-15; not a result): 11,358 templates, 274 families observable (was ≤ 61),
  14–20 s/day; first attempt found SP2 half-dead → valid fraction 0.49 for every family under the all-stations-valid rule
  → dead-station rule (> 5 % guarded = absent) for both before and after. With every template the frozen cross-family rule
  becomes a competition among nearly all families (Arm B smoke: all catalogue detections removed) → RAW-tier tables too.
- T2 (withheld station): `merlin_t2.py` setup / round / score. XG raw traces are SH 50 sps → SH band accepted (wa_cache
  BANDS only had HH/BH/EH: first smoke cached 0 of 209 XG days). Round-0 candidates identical to stack_activity valA
  (24/24 on every family checked). Smoke (3 days, ≥ 10-event truth; NOT a result): recovered-vs-true median cc Olympic 0.89
  (control 0.14, ceiling 0.91), NW 0.29 (ceiling 0.42 with ~12 events), SW 0.14 (4 pairs). Round 1 changes little: the
  withheld station is ~3–10 % of the network weight.
- Queue `scripts/merlin_tests_queue.sh` (pid in results/MERLIN_T_PID), T1 and T2 as independent chains; ETA ~07:30Z.

## 2026-09-30 06:18Z — STEP 1: per-repeat probability by local FDR (Merlin) built, smoke-tested, launched
- Scripts `step1_features.py` (per-site cc/z at each candidate's aligned or rolled position, HF ratio, CPU 12 workers),
  `step1_fit.py` (densities, calibrated GBM data-vs-control, π0 from quiet days, LOCO refits, checks a–e, rate product),
  `step1_report.py`; queue `step1_prob_queue.sh` (stages mf_c34 → feat → fit → report). Extra scrambled controls = `mf_step1c34`
  (seed 20261001; its data series identical to mt1_all's, checked per day).
- Smoke (4 days; not a result): recomputed network cc equals the MF's to ≤ 2e-6 for data and all 4 controls; pipeline runs to the report.

## 2026-09-30 08:07Z — STEP 2 (Merlin): ≤ 3-km group detection + member attribution, built, smoke-tested, launched
- Groups: complete linkage ≤ 3 km on lat_use/lon_use → 215 groups (150×1, 49×2, 7×3, 4×4, 3×5, 1×6, 1×23); the 23-group = 21
  unlocated Olympic A families sitting at the B013/B943 approximate position + 2 NW (a position-rule group; flagged).
- Found: families' detection times have their own origin conventions (offsets up to ±9 s) and the own S picks of the ≥ 30-event
  stacks are noise-dominated (within-group offset spread ~9 s) → members aligned by stack cross-correlation (δ, per-station τ ±1 s),
  with the co-detection mode as a prior; where the stack peak is strong δ = co-detection offset to ≤ 0.15 s. Reliable alignment
  for 41/114 non-reference members; the rest are noise-like on permanent stations.
- Group templates: 7,587 ≥ 30-event (181 groups), median ×1.54 the events of the best member in multi-member groups; singleton
  templates identical to T1 (4,913 compared, max diff 0); smoke MF: singleton detections + all controls identical to mt1_all.
- Bug caught: NpzFile re-decompression again (singleton check loop) → arrays read once.
- Smoke (4 days, NOT a result): attribution top-1 A2011 0.95 vs chance 0.49 vs most-active 0.80; V2010 test 0.66 vs 0.32 / 0.70.
  P(member) is as peaked on the scrambled control as on data (softmax always picks someone) → peakedness is not evidence.
- Queue `scripts/step2_groups_queue.sh` (pid results/STEP2_PID), ETA ~08:45Z; report results/STEP2_GROUPS_RESULTS.md.

## 2026-09-30 13:05Z — PNET (Merlin step 6) built, smoke-tested, launched
- Input check: WA cache = PNET's training chain (20 Hz, 1–8 Hz); PB EH1/EH2 must be rotated (B013 31.5/301.5, B943 301.5/31.5, …);
  rotation gains ~0.01 AUC at B943/B001 over raw EH1/EH2 on one day; NLWA 4/94; HDW one horizontal → Z only.
- Smoke (1 day, 2010-08-20, NOT a result): the SVI model returns ~0 on 1-C (Z, 0, 0) input (median 0.002–0.007 at BLN / HDW / Z-only B943).
  At the labelled S (t + template offset) vs random times of the same episode day: B943 AUC 0.42–0.54, B001 0.49–0.58, NLWA 0.43–0.70,
  PGC 0.53–0.55 — near chance off the shelf. CF profile around B943 labels peaks ~1 s before the envelope peak (onset) → window centred
  on the template start (= envelope peak − 0.5 s), ±0.5 s.
- Smoke step-1 features (3 days): pnet_max AUC data vs control NW 0.50–0.55, SW 0.63; matched vs unmatched NW 0.62–0.70 (step-1 P 0.75–0.91).
- Fine-tune smoke (2 days, random-window validation, leaky, NOT a result): internal AUC 0.73 → 0.92 in 400 steps (50 steps/s on the H100).
- Queue `scripts/pnet_queue.sh` pid in results/PNET_PID; ETA ~15:15Z.
