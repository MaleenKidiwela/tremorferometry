# NOTES — lfe_array_discovery
## 2026-09-20
- Workspace created after the user's observation that the array waveforms behind every published family set are public (FDSN verified). PLAN.md v0 drafted; Merlin review requested.
- Merlin review received and saved verbatim (`results/merlin_review_plan_2026-09-20.md`); PLAN.md §10 amendments written (sites not stations; N-Olympic 2010-08 first, Puget 2008-05; Step 1 north positive control of the autocorrelation engine on disk before any WA byte is scored; engine prerequisites from P4-0b; claim ladder; withheld sets; download order). Awaiting user go (download ≈150–300 GB for the first window; Step 1 needs no download).
- 2026-09-20T13:07Z user GO. Two Opus agents launched: D (download queue + metadata/site table, order per PLAN §10) and E (autocorrelation engine w_ + Step 1 north positive control on lfe_catalog_repair/cache20hz). Station list snapshot saved results/fdsn_stations_40-50N_2003-2014.txt.

### 2026-09-20 — agent D: download + metadata queue BUILT and LAUNCHED (pid in logs)
**Scripts (all `scripts/d_*`, hashed in `results/freeze_windows.json`):** `d_common.py` (box/network
constants, 1-km site merge, Slab2 active cells, per-day band resolution), `d_meta.py` (tranche 1
channel-level StationXML + flags + site table), `d_queues.py` (all queues + freeze + quiet days),
`d_fetch.py` (REUSES `lfe_picker/scripts/r_download_day.fetch_day` verbatim behind an 8-thread pool,
symlinks days already present in ../lfe_extension|lfe_picker|lfe_catalog_repair traces),
`d_cache.py` (REUSES `o_cache_month.cache_one` + the rational twin `u_cache.cache_one_rational`),
`d_eq.py` (ComCat intraslab catalog + ±3 min windows), `d_seal.py` (sha256 + DO-NOT-READ marker),
`d_download_queue.sh` (the single nohup queue, per-tranche markers `results/markers/DL_*_DONE|FAILED`,
final `results/DL_ALL_DONE`).

**Tranche 1 (done, `results/markers/DL_META_DONE`):** 21,757 channel epochs → 665 stations in the box
with a 2006–2012 epoch, 281 sites at 1 km (`results/{channels_all.csv,channels_2006_2012.csv,
station_flags.csv,sites_array_era.csv,skipped_ehz_only.csv}`). Key results, all of which change the
plan's arithmetic:
- **XG is only half 3-C**: 102 of 188 stations are 3-C (EH?/SH? at **50 Hz → rational chain**),
  **86 are EHZ-only** and are excluded from every sum. The 102 3-C stations fall in **7 sites** at
  1 km (GC+SL chain together) / **8 named arrays** (BH CL DR GC LC PA SL TB).
- **UW is overwhelmingly EHZ-only in the array era**: 187 of 207 box stations. **HDW is EHZ+EH2
  (L-4C vertical) until 2022-09** → it can never enter a 3-C sum in 2007–2011; it is queued anyway
  (explicitly named, `three_c=False`) so the withheld design can still use Z-only coherence.
  GNW BH 50 Hz until 2010-09-10 then 40 Hz from 2010-12-01 (a gap in between), LRIV BH 40 Hz from
  2009-08-27, DOSE BH 40 Hz from **2010-07-21** (one day into the 2010 window).
- ZK 97/97 and Z3 47/47 are 3-C EH at 100 Hz (frozen chain); XU 68/68 3-C but 57 at 50 Hz
  (rational); PB 19/19 3-C EH 100 Hz, sensor HS-1-LT, **all dip −90 → zsign +1** (no PQ +90 case in
  this box); CN VGZ/PFB HH 100 Hz; TA/YW BH 40 Hz.
- **Trap found and fixed in the smoke test:** many XG/XU stations declare BOTH EH? and SH? for the
  same epoch while only one has data (XG.CL13 2010-08-20 → nodata on EH, data on SH). Every queue row
  now carries an `alts` list and `d_fetch.py` falls back through the declared bands (deleting the
  premature `.nd`) before giving up. Without this ~40 % of XG would have been recorded as nodata.

**Queues** (`results/queues/q_*.csv`, `results/queue_summary.csv`), 141,473 station-days total:
tranche     group              queue  n_stations  n_stations_3c  n_days  n_station_days  n_rational       outroot
    T2A      thin     q_T2A_thin.csv          33             32      53            1748         477        traces
    T2B     quiet    q_T2B_quiet.csv          32             31      20             640         180        traces
    T2C   xg_rest  q_T2C_xg_rest.csv          91             91      73            5651        5651        traces
    T2D        xu       q_T2D_xu.csv          25             25      73            1606        1446        traces
     T3    sealed    q_T3_sealed.csv         122            121      53            5062        3631 traces_sealed
     T4 puget2008 q_T4_puget2008.csv         183            182      77            8218        4224        traces
     T5  cafe2007  q_T5_cafe2007.csv          68             67      74            5019        4366        traces
     T7 extension q_T7_extension.csv          17             17    7038          113529        1349        traces

**Quiet days.** 2010 = catalog-quiet, 20/20 found in 2010-06/07 (no PNSN epicentre within 100 km of
any of the 403 tremor-active Slab2 cells on the day or ±3 d) — no extension into 05/10 needed.
2008 (20) and 2007 (15) are **literature-quiet only** (PNSN starts 2010-01-06): ≥30 d from the
published ETS windows 2008-05-01..06-05 and 2007-01-10..02-20. All recorded in
`results/freeze_windows.json`.

**Seal.** The 2011-07-20..09-10 set goes to `traces_sealed/` and those days are **removed from the
extension queue** so the only copy on disk is the sealed one; `d_seal.py` then writes
`results/sealed_2011_manifest.sha256` + `results/SEALED_2011_DO_NOT_READ`. Nothing caches it.

**Tranche 6 catalog already on disk:** 131 intraslab events M≥2.0, 25–50 km, in the box,
2006-06-01..2012-06-30 (`results/intraslab_eq_2006_2012.csv`).

**Smoke test** (3 XG + 2 PB station-days, then the cache on all of them): fetch OK, sample rates seen
**50 Hz (XG SH?/EH? → rational twin) and 100 Hz (PB EH? → frozen chain)**; 20-Hz npz are (3, 1728000)
float32, order Z/N/E, `chain` recorded on every row. PB B001/B003 2010-232 came from an existing
`lfe_extension/traces` day by **symlink**, not re-fetched.

**How to check:** `tail -f logs/download_queue.log`; `ls results/markers/`;
`wc -l results/download_manifest.csv`; `python -c "import pandas as pd;
M=pd.read_csv('results/download_manifest.csv'); print(M.groupby(['tranche','status']).size())"`.
Final marker `results/DL_ALL_DONE`.
- Agent D report (queue pid 1080878, log logs/download_queue.log, markers results/markers/DL_*): T1 metadata done (665 stations, 281 sites). Tranches: T2A thin 33 stns/1,748 sta-days; T2B 20 catalog-quiet 2010 days; T2C XG members 91/5,651; T2D XU 25/1,606; T3 sealed 2011 122/5,062 → traces_sealed/; T4 2008 183/8,218; T5 2007 68/5,019; T6 131 intraslab EQs; T7 extension 17 permanents 113,529 sta-days (tail). Rate 4.15 sta-day/s → discovery tranches in ~2.5 h, extension tail 8–10 h.
- **Findings that amend the freeze:** (1) XG is 102/188 3-C (86 EHZ-only excluded) → **7 sites @1 km, all 50 Hz = rational chain** → the rational-vs-frozen equivalence check (P4-0 §5) is owed BEFORE any XG-derived template product; (2) **UW HDW is EHZ+EH2 until 2022-09** → cannot be a 3-C withheld site in the array era → N-Olympic withheld set must be re-chosen (B013/B943 cluster + GNW + one 3-C PB borehole, e.g. B018 or B017 — compute coverage at Step 0); (3) many XG/XU declare both EH? and SH? with data on only one → fetch uses an `alts` fallback list.

## 2026-09-20 — Agent E: array-discovery engine built (w_) + STEP 1 launched
**Engine.** New scripts under `scripts/`, prefix `w_`, importing (never editing) the lfe_picker
`v_`/`g_`/`q_`/`p_`/`s_` code: `w_common.py` (constants, bed loaders, data-driven site merge,
quiet/candidate-hour selection, in-process patching of the `v_` modules), `w_autocorr.py` (GPU
network autocorrelation), `w_cluster.py` (FA-only T + components/seed clusters), `w_grow.py`
(align_stack handoff, growth, template-scramble null, survival), `w_iterate.py`, `w_dedupe.py`,
`w_locate_gates.py`, `w_null_siteshift.py`, `w_score_step1.py`, `w_diag.py`, `w_freeze.py`,
`w_step1_queue.sh`.

Merlin's four prerequisites are installed by `w_common.patch_v()` in process, so nothing in
lfe_picker changes: **P1** no T=8.0 floor (FA-only, engineering floor 5.0, and every threshold file
carries a `floor_binds` column); **P2** data-driven 1-km single-linkage SITE_MEMBERS — on this bed
it reproduces the hand-written merge exactly (6 sites: B001, B003, B004, LRIVq=B005+B006+B007+LRIV,
NLLB, PGCq=B009+B010+B011+PGC); **P3** station_z median/MAD over VALID samples only; **P4**
MIN_SITES 4, per-day availability never intersected.

Two defects found in the inherited engine while wiring it up:
- `v_engine.DayData` keeps a SINGLE `self.NF` (the last station's rFFT length) but uses it for every
  station in `station_z`. Harmless only while every station-day has identical sample counts; a
  short/partial day raises `size of tensor a must match tensor b`. Fixed in the `w_` subclass
  (per-station `NFs`); **`v_engine.py` on disk still has it** — flagged for the P4 close-out.
- the gap bug (C-8) is the same code path and is fixed in the same subclass.

**Bed / freeze.** `results/freeze_step1.json` (all constants + sha256 of the 12 `w_` scripts and of
the 13 imported lfe_picker scripts + availability table + quiet/candidate hour lists + hidden-list
hashes) written BEFORE launch; hash check re-run as the queue's first step.
Bed: 12 detection stations / 6 sites / 52 days (2012-08-15..10-05), 618 station-days, withheld
VGZ/SNB/B926. **Quiet days: only 9 exist on the bed** (no PNSN tremor within 100 km and >= 3 days
from any tremor day) — Merlin's ">= 20 quiet days x 2 h" is impossible here, so the binding clause
(>= 40 quiet hours on the same sites) is met with 5 hours per quiet day = 45 quiet hours; the
deviation is written into the freeze. Candidate hours: 563 over 32 days.

**Smoke test (all stages, end to end).** 1 tremor hour autocorrelation 0.8–8.5 s (warm/cold) at
3.8 GB GPU, 12.4 k pairs at T_record 4.0; clustering 3 components -> 1–3 seed clusters; handoff 2/2
clusters with >= 4 sites in 3.4 s; quiet-day FA rule gave T_seed 6.00–7.75 with the **engineering
floor binding for 0/2 templates** (the FA rule binds — the P4-0b complaint is fixed); growth 3.3 s
per station-day for a 2-template bank; iterate 2 min for 2 templates / 323 hits; dedupe 2 s; locate
+ gates 12 s. Empty-input paths exercised and safe.

**⚠ PROVISIONAL FINDING from the smoke tests (to be settled by the queue, not by me).** Three
independent probes say single-pair network autocorrelation has almost no SNR margin on THIS bed:
(i) a quiet hour produced MORE pairs (3 439 at z >= 4) than the busiest tremor hour (1 056);
(ii) the per-station z at the TRUE lag of catalogued kept-family event pairs, cut moveout-aware at
the family's own predicted S, is 1.5–2.5 at EVERY window length 3–20 s, i.e. at the max-of-41-
samples noise level; (iii) of the 300 strongest pairs of a tremor hour, 23 have a lag matching a
real inter-event interval at a common offset, vs 9–20 for a lag-jitter null. The reason is
mechanical: a self-correlation (perfect repeat) scores z = 14 here, because the lower-half MAD of
the hour's cc distribution is ~0.07 — so a pair needs cc >= 0.35 to reach z = 5, and real LFE pairs
at 20–90 km hypocentral inside tremor sit at cc ~ 0.1–0.15. The bed's S-moveout spread is 12.3 s,
twice the 5–6 s Merlin's 10-s common window was sized for (that number was for the 45-km WA array).
`w_diag.py` measures all three as a queue stage so the write-up rests on the frozen run, not on
these probes. Do NOT read this as the Step-1 verdict yet.

**Launched** `nohup bash scripts/w_step1_queue.sh > logs/step1_queue.log 2>&1 &` (pid 1085937);
order: quiet-hour T -> autocorrelation on 563 tremor hours -> clusters -> handoff/growth/scramble/
survival -> 2 iterations -> dedupe -> locate + Gates 0-2 -> 3 whole-chain site-shift null
realisations -> scoring (q_rediscover relative clauses vs the 10 isolated hidden families) +
engine diagnostics. Markers `results/STEP1_DONE` / `results/STEP1_FAILED`.

## 2026-09-20 (later) — Agent E: Step 1a closed, empty-path crash fixed, Step 1b launched
**Crash fix.** The queue died at `grow1` because a bank with zero templates wrote a header-less
`p1_thresholds.csv` and `pd.read_csv` raised `EmptyDataError`. Every stage now survives "nothing
survived": `w_common.read_csv_safe()` + explicit column sets (`THR_COLS`, `SURV_COLS`, `HIT_COLS`,
`DEDUPE_COLS`), early empty-bank returns in `w_grow.quiet/grow/survive`, and empty-safe reads in
`w_iterate`, `w_dedupe`, `w_locate_gates`, `w_score_step1`. Re-ran the Step-1a tail (no nulls, per
the coordinator): it now runs clean through iterate → dedupe → locate → gates → score → diag.

**Step 1a result** → `results/STEP1_RESULTS.md`. 563 tremor hours → 248 397 pairs at ≥5.0 →
**T = 7.75 by the FA rule** (not the 5.0 engineering floor; 45 quiet hours, 0.69 false pairs/hour,
0 quiet components ≥5) → 3 740 pairs at T (**6.64/tremor hour vs 0.69/quiet hour = 9.6× excess, so
the statistic is not inert**) → **2 seed clusters** → 16 growth hits → **0 survivors** → 0/10
isolated hidden families. Frozen diagnostics (`results/step1_diag.json`): **z_self = 14.4** (the
perfect-repeat ceiling, because the hour's cc lower-half MAD is ≈0.07, so a pair needs cc ≥ 0.35
for z = 5); **true-pair z = 2.13 vs jitter null 1.98**, only 16 % above their own null max;
**0/2000 of the strongest pairs** match a real inter-event interval (null 0, 2, 0). Site geometry:
median hypocentral 52.6–91.4 km, **median S-moveout spread 13.8 s** against the 10-s common window
Merlin sized for the 45-km WA array. Verdict: uninformative about the engine, informative about
distance — P4-0b starvation again.

**Step 1b bed — what the FDSN actually has.** The premise that Bostock-era dense stations are alive
in 2012 is **false**: the POLARIS BC array (GLBC, SHVB, TWKB, MGCB, TSJB, JRBC, KELB, KHVB, TWBB,
TWGB, SILB, SSIB, LCBC, GOWB) ran **2002–2007** and was gone by 2012 (`results/bostock21_epochs.csv`);
network PO in 2012 has only VI55/VI57 at 230 km. Of the SeisBench Bostock 21-station list
(`results/bostock_seisbench_stations.csv`), only LZB, PFB, PGC, SNB, SOKB, VGZ, YOUB exist in 2012 —
and a delivery probe (metadata ≠ delivery) shows **C8 SOKB 0/8 and CN PFB 0/8 at BOTH IRIS and
NRCan**, i.e. the single best station (24.8 km, 67/130 families) cannot be had.
Deliverable and downloaded: **CN LZB (28.0 km, 54 families ≤25 km — the only new close site)**,
CN YOUB (59.8 km, withheld), and C8 PA02/PA03/PA04/PA05/PA12 + BPCB from NRCan (partial days).
Step-1b sites by a **written rule** (maximal deliverable subset, added by increasing epicentral
distance, with median S-moveout spread ≤ the frozen 10-s window): **LZB + PGCq (now up to 9
members) + BPCB + LRIVq + B003 = 9.7 s**, against 13.6 s for Step 1a's six. Withheld:
VGZ/SNB/B926 **+ YOUB** (coverage rule: dropping YOUB costs the fewest families with ≥5 close sites
within 60 km, and it is the least-covering of the two candidates). **No engine constant retuned**
— window, hop, thresholds, survival, dedupe, gates are byte-identical to `freeze_step1.json`;
only the bed changes, and the deviation is written into `results/freeze_step1b.json` for Merlin.

**⚠ RECOMMENDATION for the coordinator / Merlin — Step 1c.** The bed the positive control was meant
to be exists, but in **2004–2005, not 2012**: Bostock's catalog has 59 659 (2004) and 57 504 (2005)
events, and a delivery probe on 2004-07-15 / 2005-09-05 returns waveforms for **11 of the 21 dense
stations** — GLBC (25.7 km), SHVB (26.3), MGCB (30.8), TSJB (30.8), JRBC (31.7), KHVB (34.1),
TWBB (35.6), TWGB, SILB, SSIB, LCBC — from IRIS under `C8` and NRCan under `PO`. That is a genuine
25–36 km dense bed with a moveout spread the 10-s window was designed for. I did **not** launch it:
it changes the window and needs a pre-2010 candidate-hour and quiet-day rule (PNSN starts 2010-01;
Merlin's Q3 sanctions the array hourly-energy flag for *hour selection only*), which is a method
decision for Merlin, not for me.

**Launched** `nohup bash scripts/w_step1b_queue.sh > logs/step1b_queue.log 2>&1 &`
(download ≤4 concurrent so agent D keeps its 8 → frozen-chain cache → freeze_step1b → quiet T →
autocorrelation → clusters → growth → survival + scramble → 2 iterations → dedupe → locate/gates →
3 site-shift nulls → scoring → diagnostics). Markers `results/STEP1B_DONE|FAILED`. Same acceptance
bars: ≥3/10 isolated, spurious-split ≤0.1, null ≤1 per realisation. The number to read first is
**D2 true-pair z at LZB (28 km)** against Step 1a's 2.13 — it is what predicts whether any 25–35 km
dense bed can work, and so whether Step 1c is worth the download.

## 2026-09-20 (Step 1c) — the 2005-09 POLARIS bed, built and chained behind 1b
**Authorised by the coordinator** as a recorded deviation for Merlin's audit: the 2012 window has no
dense array (POLARIS BC ran 2002–2007), so the 2005-09 ETS is the only north data with array-class
geometry.

**Tremor-hour source: the ARRAY-ENERGY FLAG, not WTD.** The Slow Earthquake Database is reachable
(`http://www-solid.eps.s.u-tokyo.ac.jp/~sloweq/`, `https://eqdb.sakura.ne.jp/`) and its index sheet
names the catalogue `WTD(Cascadia)-Tremor-Cascadia` (2005-01-01..2014-12-31), but its only public
interface is a POST to `php/sqlIntrTrem.php` whose parameter shape I could not reproduce — it
returns `[]` even for `PNSN-Tremor-Cascadia` inside its own date range — and `/download/` is 403.
So hour selection uses Merlin's Q3 fallback (`w_energy_hours.py`): per station, the hourly median
3-C band energy of the frozen 1–8 Hz cache divided by that station's own quiet-day median; per site
the median of members; network = median over sites; an hour is a candidate at or above the p95 of
the same statistic on the quiet days. Hour selection only — never an event seed. Window from the
literature (2005-09-03..09-16 ± 5 d = **2005-08-29..09-21**, the 14-month sequence 2003-02 /
2004-07 / 2005-09 / 2007-01); quiet days **2005-06-25..07-14** (20 days, ≥50 d before onset),
recorded as "literature-quiet". **No Bostock event time enters window, hour or day selection** —
only the scorer.

**Delivery probe** (`w_probe2005.py`, 4 days, 5 networks × 3 bands × 2 providers;
`results/step1c_probe.csv`). Delivering, with median epicentral distance to the 130 families and
families within 25 km: **GLBC 25.7 km/62, SHVB 26.3/58, LZB 28.0/54, TWKB 29.5/42, MGCB 30.8/45,
TSJB 30.9/48, JRBC 31.7/44, TWBB 35.6/41, PGC 37.8/24, VGZ 40.8/36, SILB 46.3/17, SSIB 46.5/11,
LCBC 47.4/20, B005+B006+B007 54.7/11, SNB 57.2/1, B003 60.3, B004 66.2, B001 69.4** — C8 and PO at
IRIS and NRCan in roughly equal measure. Not delivering in 2005: SOKB, PFB, KELB, KHVB, TWGB, GOWB,
BPCB, YOUB, GOBB, OBC, B009/B010/B011 (which only start mid-September 2005). Excluded by frozen
rules: **UW STW** (Z-only EHZ — never summed into a site statistic) and **UW SQM** (BH 50 Hz;
50/20 is not an integer, `o_cache_month` raises, and the P4-0 rational twin is not part of the
frozen chain) — both written to `results/step1c_excluded.csv`.

**Bed (`results/step1c_bed.json`, frozen before launch).** 1-km merge, then the SAME written spread
rule as 1b: **12 detection sites** — GLBC, SHVB, TWKB, MGCB, TSJB, TWBB, PGC, SILB, SSIB, LCBC,
B005q(B005+B006+B007), B003 — **median S-moveout spread 9.81 s ≤ the frozen 10-s window**, against
13.8 s on the Step-1a permanent bed. **Six sites at ≤ 36 km, where Step 1a had none below 52.6 km**;
**130/130 families have ≥5 detection sites within 60 km epicentral**. Withheld: **VGZ + SNB**
(B926 does not exist in 2005) **+ the dense pair LZB and JRBC**, chosen by the coverage rule —
removing them still leaves 130/130. B001 and B004 are downloaded but enter no sum (outside the
spread rule); recorded. **880 station-days** (24 ETS days + 20 quiet days × 20 stations).

**Scoring deviation (frozen).** The relative clauses need the frozen template detector's per-family
recall/purity, which exists only for 2012 (`truth_template_reference.csv`); no 2005 reference can
exist, because the frozen templates live on stations absent in 2005. The 2012 bar is therefore
TRANSFERRED unchanged and fixed before scoring, and the ABSOLUTE recall/purity are reported
alongside in `step1c_score_isolated_detail.csv` so the bar can be re-read. Acceptance unchanged:
≥3/10 isolated, spurious-split ≤0.1, site-shift null ≤1 per realisation.

**Launched** `nohup bash scripts/w_step1c_queue.sh` — Phase A (no GPU, running now): download ≤4
concurrent + frozen-chain cache → array-energy hours → full freeze; **Phase B: a shell wait on
`results/STEP1B_DONE|STEP1B_FAILED`, polled every 300 s**; Phase C (GPU): quiet T →
autocorrelation → clusters → growth → survival + scramble → 2 iterations → dedupe → locate/gates →
3 site-shift nulls → scoring → diagnostics. Markers `results/STEP1C_DONE|FAILED`, log
`logs/step1c_queue.log`. **No engine constant retuned** — only the bed changes.
The number to read first is again **D2 true-pair z**, now at 25–36 km, against Step 1a's 2.13.

### Correction (same day, after the Step-1c launch) — LZB has no 2005 data
A targeted 6-day delivery retry for the key permanents finished after the bed was frozen
(`results/step1c_delivery_retry.csv`): **PGC NRCan/CN.BH 2/6, VGZ 1/6, SNB 2/6 — but LZB 0/6**,
and PFB / SOKB / B009 / B011 all 0/6. When I assembled `results/step1c_stations.csv` I had added
PGC, VGZ and LZB by hand with `n_days` 2/1/1 *before* that retry returned; the LZB value was my
guess, not a measurement (the first probe had already given LZB 0/4). **My error — an unmeasured
number went into a file that fed the bed setup.**

Consequence: LZB is one of the two dense stations the coverage rule picked as WITHHELD, and it has
no data, so it contributes no Gate-1 evidence. The run is unaffected computationally — LZB is in no
detection sum either way, and a station with no data is inert wherever it sits — so I have **not**
re-frozen the bed: re-freezing to move a dataless station between two groups would change nothing
and would break the discipline the freeze exists for. What it does change is the honest reading of
**Gate-1 power**, which rests on **JRBC** (dense, 31.7 km, 44 families ≤ 25 km, 3/4 days) plus the
partial VGZ (1/6) and SNB (2/6) — three withheld stations, one of them strong, not four.
`results/step1c_stations.csv` is left exactly as it stood when the bed was set up; the retry is
recorded beside it rather than edited into it.
- 2026-09-20T14:40:59Z Step 1b and 1c queues were SIGTERM-killed (-15) at ~14:40Z when agent E's session ended (children of its background shell); relaunched detached with setsid by the orchestrator; 1b re-runs from the start (no per-step guards), 1c resumes its fetch then waits for STEP1B_*.

## 2026-09-20 (late) — the "silent deaths" of the Step-1b queue: cause found, and it was mine
**Cause — not OOM, not the code, not the cgroup.** While fixing two unrelated bugs I ran
`pkill -f w_step1b_queue.sh` twice (14:42 and 14:45 UTC). `pkill -f` matches the whole command
line, so it also killed the coordinator's *detached* relaunches — and the queue had **no signal
trap**, so SIGTERM killed bash silently: no `=== STEP` line for the stage it was entering, no
`STEP1B_FAILED` marker, bash and python both gone ~60 s after start. Every observed symptom
follows from that. The evidence: `/sys/fs/cgroup/memory.events` shows **`oom_kill 0`** (and no
nested cgroups have any), `dmesg` is not permitted but no kill is needed to explain it, and the
**foreground reproduction the coordinator asked for passes cleanly** — `w_freeze.py` exit 0, then
`w_autocorr.py --hours quiet` exit 0, 45 hours, **41 294 pairs** (it was 0 before the t0 fix
below), `oom_kill` still 0. My error, twice over: the `pkill`, and shipping queues with no trap.

**Also fixed in the same pass — an absolute-time bug that had corrupted Step 1a as well.**
`DayLoader.hour` indexed each hour from the *file's own* `t0`, and `autocorr_hour` then required
every station's start to match the first station's to 1e-6 s. The frozen cache does not store a
common `t0`: it is midnight +0.0000 s for the CN/C8 broadbands, **+0.0084 s** for the PB boreholes
and **+0.0150 s** for UW.LRIV — all far below one 20-Hz sample (50 ms), but not equal. The
equality test therefore dropped every station whose offset differed from the first one's:
**Step 1a ran on 5 sites / 9 stations, not 6 / 12** (NLLB, PGC and LRIV dropped from every hour),
and Step 1b on 3 sites, below MIN_SITES — which is why it returned **0 pairs over 563 tremor
hours**. The cut is now indexed from the ABSOLUTE target time per station and rounded to the
nearest sample; the residual (≤ half a sample, 25 ms, constant per station and therefore exactly
cancelling inside each station's own autocorrelation) is recorded per hour as
`max_grid_residual_s`. Step 1a is being re-run (`scripts/w_step1a_refresh.sh`, nulls excluded per
the coordinator); **`results/STEP1_RESULTS.md` is superseded until it finishes** — its
diagnostics D1/D2/D3 are unaffected (w_diag never used the equality test), but its T, pair and
cluster counts were computed on the reduced network. Invalid outputs archived under
`results/discarded/`.

**Two more bed bugs found and fixed before anything was scored.**
(i) *Site ids are not stable.* The id is a data-driven tie-break (member with the most cached
files), so the {B005,B006,B007,LRIV} cluster is `LRIVq` on the Step-1a roots and `B006q` once the
study cache joins them. My hard-coded `STEP1B_SITES` matched on the id, so that site was silently
dropped: 5 sites became 4. Sites are now selected by **MEMBER station** (`STEP1B_SITE_KEYS`, and
`site_keys` in `step1c_bed.json`), which is stable; a missing member logs a warning.
(ii) *Roster leak.* `w_fetch_bed2012b.roster()` read `W.extra_stations()`, which is the **union**
of the 1b and 1c rosters (it must be, so the station table knows every bed's coordinates), so the
1b fetch tried all 27 stations × 52 days of 2012 — including stations that exist only in 2005. It
now reads its own csv; the stray `C8.GLBC` / `C8.SOKB` 2012 trees (no-data markers only, no
waveforms) were removed.

**Hardening delivered.** `w_step1b_queue.sh` rewritten with (1) **single-instance `flock`** on
`locks/step1b.lock` — a second launch exits 0 instead of racing, which matters because
`w_step1a_refresh.sh` execs this queue at its end; (2) **signal traps** on TERM/INT/HUP plus an
EXIT trap, all writing `results/STEP1B_FAILED` with the signal name, the stage that was running
and the last 120 log lines; (3) **resume guards** — every stage is skipped when its own output
file exists, so the relaunch skipped fetch, freeze and ac_quiet and went straight to the
threshold. `scripts/w_gpu_gate_watchdog.sh` implements the release rule the marker cannot cover
(a SIGKILL leaves no trap): if for 600 s there is no `w_step1b_queue.sh` process, no marker, and
the GPU is idle (<2000 MiB, <5 %), it writes `STEP1B_FAILED` itself so Step 1c can never be
stranded. It is a separate script on purpose — `w_step1c_queue.sh` was already running, and
editing a bash script that bash is mid-way through reading corrupts its execution.

**Note on the station counts the coordinator queried**: "16 detection stations / 5 sites /
740 station-days" is the CORRECT Step-1b bed after the member-key fix (B003, B005–B007, LRIV,
BPCB, LZB, B009–B011, PA02–PA05, PA12, PGC); the "1404 station-days" fetch line was the roster
leak, now 8 stations × 52 = 416 downloaded, 330 with data.

**Live**: 1a-refresh (ac_tremor 401/563) → 1b **pid 1101224** (T = 7.25 by the FA rule on the
5-site bed, ac_tremor 102/563, verified alive across two foreground 110-s waits) → 1c **pid
1097670** (Phase A fetch 500/880, Phase B gate armed) → watchdog.

### 2026-09-20 — agent D: delivery report, Step 2 site table, Step 4 site probe, Step 0 coverage
Full write-up `results/P_STEP2_SITES.md`. New scripts `scripts/d_{report,sites,probe,coverage,equiv}.py`
(the probe REUSES `lfe_picker/scripts/u_fartremor.hf_one` verbatim with `u_common.TRACES` repointed
at this study's traces/, so no lfe_picker file is touched and nothing is re-downloaded).

**Delivery.** T2A–T2D 9,367/9,645 station-days (93–99 %), T3 sealed 4,453/5,062 (88 %), T4 7,380/8,218
(90 %), T5 4,530/5,019 (90 %), T6 9,143 event windows, T7 still streaming. 357 GiB on disk so far;
zero hard errors, everything missing is `nodata` at the archive. 57 station-tranche entries below 70 %
delivery (`results/delivery_low_stations.csv`).

**Five findings that change the plan.**
1. **ZK is vertical-only in the DATA** — StationXML declares EH1/EH2/EHZ but 1,604 of 1,642 delivered
   ZK station-days carry one channel (EHZ). Merlin's Step-2 failure signature fires for ZK, not XG.
   The 2008 window's 3-C power is XU + PB and nothing else. Every 3-C count in the site table is now
   DELIVERED 3-C (`n_channels >= 3` in the manifest), never declared 3-C.
2. **ZK also has no quiet-day baseline**: it delivered only 2008-03-05, 03-11 and 2008-05-05..05-22 —
   the ETS itself. The 2008 probe needed an in-window quiet fallback (6 quietest days by the array's
   own flag) for ZK/Z3/YW; each station records which baseline it used.
3. **XG = 7 sites at 1 km / 50 Hz / rational — confirmed**, and XG *is* genuinely 3-C in the data.
   But the **SL array died 2010-06-29** (before the window) and **XG.CL is admitted 0/16** by the
   probe (1.9 dB contrast; it is the easternmost array). Six XG arrays clear the ≥5-admitted-member
   bar: BH 9, DR 11, GC 13, LC 13, PA 7, TB 12 — exactly Merlin's Step-4 acceptance threshold.
4. **Merlin's withheld rule is infeasible as written** (3 permanent sites within 60 km of ≥80 % of a
   sub-region's cells): only 2 such sites exist in N-Olympic, 1 elsewhere — a 60 km radius cannot
   cover a 150 km-long sub-region. Two costed options are in §5 of the write-up. Option A (keep every
   extension permanent in the sum) costs no coverage but leaves only ONE withheld site that can
   actually see the cells (B013+B943, 41.5 %; GNW 8.5 %; B018 0 %) → Gate 1 "PASS at ≥2 withheld
   sites" is unreachable for N-Olympic. Option B (withhold B001 + B005-cluster + VGZ from the
   DISCOVERY sum only) buys two real tests for 94.6 → 83.1 % coverage. **Recommended: B.**
   Puget 2008 has essentially no withheld capacity — B018 (85.6 %) is the only admitted permanent
   site that sees the cells, and no UW station is admitted there (HDW is EHZ/EH2 until 2022-09).
5. **The 2008 hourly-energy flag works.** With no tremor catalog before 2010-01-06, the Merlin Q3
   fallback (per-station hourly median 6 s power over its quiet-day median → network median → hot
   above the quiet p95 = 2.91 → tremor day at ≥6 hot hours) reconstructs the published 2008-05 ETS
   from the waveforms alone: it fires 2008-05-05..05-23, 05-29, 06-10 and nowhere else, with 0–2 hot
   hours on every literature-quiet day except 2008-03-11 (7 — probably real tremor on a "quiet" day;
   excluded from the baseline). `results/hour_energy_flag_2008.csv`.

**Probe results.** 2010: 97/131 stations admitted (≥4/6 days ≥+3 dB and local−far ≥+3 dB) — XG 65/85,
XU 19/22, PB 10/19, UW 2/3, CN 1/2. CN.PFB rejected on local−far (−10.6 dB: the "far" S-WA days are
louder at PFB, so that day set is not far for Vancouver Island). UW.DOSE UNTESTABLE (installed
2010-07-21, after every T2B quiet day). 2008: 52/82 admitted on the contrast clause alone (no far days
exist without a catalog) — XU 40/57, PB 10/16, YW 2/6, ZK 0/2, UW GNW fails at 1.7 dB.

**Coverage (admitted sites only).** N-Olympic 2010: 94.6 % of its 130 active cells have ≥5 admitted
sites within 60 km hypocentral (median 10), 50.8 % within 30 km epicentral — Merlin's 100 %/75 %
becomes 94.6 %/50.8 % after the probe and the delivered-3-C audit. Puget/S-WA 2008: 84.2 % / 44.5 %
on the SAME 2010–11 cells, because no 2008 tremor catalog exists (stated in the write-up).

**Rational equivalence: PASS.** XG cc ≥ 0.99985 on every component, PB bit-identical. Caveat for the
freeze: peak sample differences on XG reach 0.38–0.57 × RMS and the RMS ratio is 0.996, so XG products
must carry `chain=rational` and XG-vs-PB amplitudes are not trustworthy below ~1 %.

**Sealed: OK.** 5,062/5,062 files hashed, 0 unhashed, 0 missing; 0 cache20hz entries and 0 unsealed
manifest rows inside 2011-07-20..09-10. Only `os.walk` + sha256 and a day-count inventory have touched
`traces_sealed/`; no waveform sample has been read.

New results files: `P_STEP2_SITES.md`, `delivery_by_{tranche,network}.csv`, `delivery_low_stations.csv`,
`sites_array_era_final.csv`, `sites_by_window.csv`, `site_probe_{2010,2008}.csv`,
`site_probe_by_array.csv`, `probe_station_days_{2010,2008}.csv`, `probe_days_{2010,2008}.json`,
`hour_energy_flag_2008.csv`, `coverage_by_subregion.csv`, `withheld_candidates.csv`,
`withheld_proposal.json`, `rational_equivalence{,_verdict}.csv`, `sealed_check.txt`.

### 2026-09-20 — agent F: WA Step 0 FREEZE written and the Step-5 discovery queue BUILT + parked
**Freeze `results/freeze_wa_v1.json`** (written before launch; hashes 12 `x_` scripts, 21 `w_`
scripts and 16 input tables). Sub-region N-Olympic 47.8–48.3 N / −123.7..−122.7, window
2010-07-20..09-10, **145 active cells** (Slab2 28–45 km, ≥30 PNSN 2010–11 epicentres within 8 km),
29 PNSN-corrected tremor days, 20 catalog-quiet days → 40 quiet hours, **454 candidate hours**.

**Engine.** No `w_` file and nothing in lfe_picker was edited. `scripts/x_common.py` binds the WA
bed onto `w_common` IN PROCESS (the same trick `w_common.patch_v` uses on `v_common`), and
`x_stage.py` runs each frozen `w_` stage unchanged. Reason for the adapter rather than a fifth bed
mode: `freeze_step1{,b,c}.json` hash `w_common.py` and the 1b/1c queues re-import it at every
stage. `x_freeze.check_constants()` verifies every engine constant against `freeze_step1c.json`
byte-for-byte; **4 documented differences only** — GRID_LAT0/LON0 48.00/−123.20 (Merlin Q4),
ACTIVE_CELL_DEP 28–45 km (Merlin F2), LOC_MIN_SITES 5 (Merlin Q4 "located" bar for WA).

**Bed.** 91 admitted stations (probe-admitted, minus XG.CL 0/16, minus CN.PFB, minus withheld) →
17 sites at 1 km. **The admitted-site moveout spread at the sub-region centre is 21.12 s, not
≤10 s**, so the Step-1b/1c spread rule was applied verbatim (add sites by increasing median
epicentral distance while the median spread over the active cells stays ≤ AC_WIN_S): **11 sites /
80 detection stations, spread 9.92 s**; dropped GNW(68 km), PL01(69), B014(73), W060(83),
B017(118), B018(121). Coverage stays **97.9 % of the 145 cells with ≥5 sites ≤60 km** (median 8).
Withheld (Option B, never in any sum or stack): **B001, the B005 cluster (B005/B006/B007/LRIV),
VGZ** — cell visibility 0.99 / 0.84 / 0.47.

**Two Step-0 measurements.** (i) **HF veto = PB.B003** (B004 fails the probe at −5.4 dB).
Its quiet-day HF-ratio distribution (400 windows, 20 days) has median 0.731 against PGC's 2.315,
so the frozen north cut 3.5 transfers to **1.106** — the same 1.51× multiple of each site's own
quiet-noise median. Reported, not used: the catalogued intraslab events are HF-*depleted* at B003
(median 0.042), so this veto must never be described as an earthquake rejector on this bed.
(ii) **Station corrections** (a5 recipe on `traces_eq`, 131 events × 86 stations, 2,485 picks):
**only 20/86 stations are usable — the MAD ≤0.5 s clause binds, not the ≥5-pick count** (35
stations have 40–110 S picks but MAD 0.5–0.9 s). The other 66 take 0.0 and are flagged. A
site-pooled table is written as REPORTED-NOT-USED evidence; pooling does not rescue it (site MAD
0.32–0.85 s).

**Queue `scripts/x_step5_queue.sh`** (pid 1122190, launched `nohup setsid -f`, parked at stage 1):
flock `locks/step5.lock`; TERM/INT/HUP/EXIT traps that always write `results/STEP5_FAILED` with
the signal, the stage and the last 120 log lines; per-stage output guards so a relaunch resumes.
**It touches no GPU until BOTH `results/WA_GO` exists AND no `w_step1` process is running**, and
it never calls `pkill`. Stages: hash check → gate → quiet-hour T (with the Merlin Step-3 stop
rule: ≥5-member components at every T ≤ 20 → STEP5_FAILED "site-coherent noise") → autocorrelation
on the candidate hours → clusters → handoff/growth/scramble/survival → 2 iterations (+1 only if
the survivors' median `s_snr_db` rose ≥0.5 dB, written to `iter_snr.json`) → dedupe → locate +
Gates 0–2 → site-shift null ×3 → Step-5b injections → `results/step5/CANDIDATES.csv`.

**Gates on the WA bed.** `g_gates` is region-configured in process (`x_common.gate_context`):
WA station table, WA corrections, this study's cache/traces, the WA site map, the HF veto above,
`accepted` = this run's own candidate list (leave-one-out, plus a far sentinel so a one-candidate
run still evaluates), `TPL` = this run's iterated templates re-cut to S−1 s, `ccp95` = the bank's
own cross-template p95 from `w_dedupe`. Gate 1 non-PASS is **UNTESTABLE** until the Step-5b
injections supply a WA control, then re-labelled FAIL only where an injected repeater at matched
N and distance passed at that same station. Gate 3 not run (Merlin F6).

**Step 5b injection** (`x_inject.py`): 4 random active cells × 3 amplitudes (0.5/1/2 × the donor
event's own SNR) × ~40 distinct northern hp events (trackE `realness_v1_1`, hp; the 30 hidden
families never donors), donor class matched (borehole→borehole, surface→surface) and distance
matched, planted into admitted AND withheld station-days at the cell's moveout in a **separate
cache root whose uninjected days are symlinks** — the real cache is never modified. Mechanics +
location error only, never a sensitivity ceiling.

**Smoke test** (ran with `X_SMOKE=1` + `results/SMOKE_OK`, gate skipped, SMOKE_OK then deleted;
`results/step5_smoke/`, marker `STEP5_SMOKE_DONE`): 1 tremor hour of autocorrelation **48.4 s →
4,444 pairs**; 4 quiet hours 8–36 s each → T = 6.50 (FA rule); clustering 0.5 s → 5 components →
2 seed clusters; handoff 9 s (2/2 clusters ≥4 sites); quiet-T + 2 growth days 14 s + 9 s (0 hits,
as expected on quiet days); iterate/dedupe/locate on the empty bank all survive "nothing survived";
fake-proposal gate pass 23–93 s → Gate 0 PASS, Gate 1 all six withheld stations FAIL (random
times), Gate 2 UNTESTABLE on geometry — i.e. the whole gate wiring runs; injection builder planted
8 snippets into 1 station-day and the planted stack at B003 reads **sdb 4.15 dB vs 0.90 dB** in the
unmodified cache. One real bug was caught and fixed by the smoke (empty `iterN_diag.csv` crashed
the third-iteration rule). **Cost projection: 454 × 48.4 s ≈ 6.1 GPU-h for the autocorrelation
pass, ~18 GPU-h more for the three whole-chain nulls** — consistent with Merlin's 5 GPU-h + 1–2
GPU-days.

New scripts: `scripts/x_{common,setup,hf_cut,corrections,stage,locate_gates,inject,null_siteshift,
candidates,freeze,smoke}.py`, `scripts/x_step5_queue.sh`. New results: `freeze_wa_v1.json`,
`step5_bed.json`, `x_station_corrections.csv`, `step5/{active_cells,admitted_stations,eq_picks,
hf_cut.json,eq_corrections_by_site_REPORTED_NOT_USED}.csv`.
- 2026-09-20T16:27:16Z Step 1c FAIL (0 seed clusters at T=12.5; nulls 0/0/0). Agent G launched for diagnostics D0–D7 + STEP1_RESULTS.md v2. WA_GO withheld.

### 2026-09-20 — agent G: Step-1 diagnostics D0–D7 and the v2 Step-1 write-up

**D0 (the crash).** `w_diag.py` died at the end of both the 1b and the 1c queue on
`S.fam_tt(det)` with `KeyError: 'BPCB'` / `'GLBC'`: `s_common.fam_tt` defaults to the
lfe_picker P3D roster table, which has none of the Step-1b/1c bed stations.  Fixed by
passing `w_common.station_table()` (the roster extended with their FDSN coordinates).
Two more defects fixed in the same file: D1 cut its self-correlation at the family's
median event over the whole bed (outside the loaded hour on a dense bed → D1 empty), and
the chosen hour could fall outside the bed's own days.  `w_diag.py` is a diagnostic and
is called by no discovery stage; its sha256 in `results/freeze_wa_v1.json` was amended
(with a dated `amendments` block) so the parked WA queue does not abort on a hash check
for a file it never runs.

**D4 — the T = 12.50 cause.** One quiet hour, **2005-07-09 23**,
supplies 100 % of the quiet-hour pairs above net_mad 8, all inside a single ~70 s burst
(23:48:03–23:49:10) seen on 9 of 9 sites at |lag| 13–56 s — an earthquake, not tremor
and not site noise (top site carries 27 % of the clipped sum).  Leave-one-hour-out:
dropping it gives T = 7.50; every other hour leaves
T at 12.50.  The 2005 quiet days were literature-quiet only (no
2005 tremor catalogue) and carry no earthquake screen.

**D4b — an engine defect found on the way.** `w_autocorr` normalises by
`E = ew.clamp_min(1e-10*ew.median()).sqrt()`; a station-hour that is more than half
flat has median 0, so E = 0, cc = inf, station z = NaN, and because `sbuf += z*vm`
makes NaN*0 = NaN the valid mask cannot rescue it — the site mean and then the whole
network sum are NaN at every lag and the ENTIRE hour records zero pairs.  On the 2005
bed 26 of
299 tremor hours
(8.7%) were voided this
way, and the correspondence with the zero-pair hours is exact (offenders: B006q in all
26, SHVB 8, B003 2).  Real bug, too small to explain 0/10.

**D5 — recoverability.** 934 catalogued same-family same-hour
pairs from 99 kept families, measured exactly as `w_autocorr`
does: median net_mad 2.11 (±0.5 s) against a matched
null 1.78; fraction reaching T = 12.5
0.000, 7.75 0.001, 5.0
0.006.  An ORACLE moveout-aware cut (each station at its
own predicted S, sites summed on the inter-event-lag axis) does not rescue it either.
No distance trend: the 21–36 km POLARIS sites read the same z as PGC at 42 km and the
boreholes at 46–92 km.  Self-lag control 25.5 —
the dynamic range is there; the pairs are not.

**D3 (the one encouraging number).** On the 2005 bed 544 of the 2,000 strongest pairs
the engine reported carry a kept family's inter-event interval at a common offset,
against a lag-jitter null of 91 — a 5.7x enrichment, where the 2012 permanent bed gave
0 against 0.  The pairing stage works on array geometry; the FA threshold (set by that
one earthquake) is what discards the pairs.

**D7 — hour flag.** 80.1% of the kept-family
2005 events fall in the 299 selected hours, all of which ran.
Not starved.

**Write-up.** `results/STEP1_RESULTS.md` rewritten as v2 (frozen verdicts for 1a/1b/1c,
beds and deviations, D0–D7, what is and is not shown, files); the old text is kept as
Appendix A and at `results/STEP1_RESULTS_v1_superseded.md`.  **No frozen verdict or
constant changed; `results/WA_GO` was not created and the WA queue stays parked.**  The
one change the diagnostics name is a quiet-hour robustness rule (hour-level outlier
screen or a per-hour FA quantile instead of the pooled mean) — stated as a finding, not
an authorised retune, and necessary-but-not-sufficient: at T = 7.50 the 1c run would
have had 3,833 tremor pairs instead of 25, but D5 says they would still not be the
catalogued LFE pairs.

- 2026-09-20T17:01:01Z STEP1DIAG_DONE; STEP1_RESULTS.md v2 written (v1 kept as appendix). Merlin consult launched on the Step 1 outcome and the repair/redesign question.
- 2026-09-20T17:26:31Z Merlin Step-1 audit saved; PLAN §11 written. Launching agent H (Step 1c′ repair queue) and agent I (Arm B jms north control + WA station-health census + hour-matched quiet set).

## 2026-09-20 (agent I) — Arm B `jms` north control launched; WA health census + hour-matched quiet sets done

**(A) Arm B — the only north control that is not a training-set test.** Merlin (§11, Q3): the
installed LFEDetect weights `cascadia`/`cms`/`jcms`/`jcs` are trained on Bostock 2015 Cascadia
stacks; only `jms` is clean. Model card printed from the installed package and recorded in
`results/armB/jms_model_card.json` and `results/freeze_armB.json`: *"trained on LFE stacks from
Guerrero, Mexico (Frank et al. 2014), the San Andreas fault (Shelly 2017) and LFEs from Nankai,
Japan (JMA catalog)"* — **no Cascadia**. Pre-registration `results/freeze_armB.json` was written
BEFORE the first pick (bars: ≥3/10 isolated at centroid ≤5 km, absolute recall ≥0.30 at ±2 s,
purity ≥0.50 under the 3-km cluster truth, null clusters ≤1; STOP if <2/10).

Queue `scripts/b_armB_queue.sh` (flock, TERM/INT/HUP/EXIT traps → `results/ARMB_FAILED`, per-stage
guards, markers `ARMB_DONE|FAILED`, log `logs/armB_queue.log`), launched detached. The engine is a
**snapshot** (`scripts_armB/`) of `scripts/` taken at queue start, because agent H is editing
`scripts/w_*.py` for 1c′ concurrently; every `W_*` 1c′ repair flag is unset so the handoff is the
frozen 1c chain. GPU shared with the 1c′ queue; no FDSN downloads.

Chain: LFEDetect(jms) on the **raw 100 Hz** `traces/bed2005_polaris` (seisbench does the model's own
1–8 Hz/20 Hz/peak-norm chain; using the frozen cache would band-pass twice) → PyOcto (Crosson 1-D
from `e2_associate`, Slab2 28–45 km prior, search box = station bbox + 0.3°, never the family
locations) → DBSCAN 5 km / whole window, ≥20 events → frozen `w_grow` handoff/quiet-T_seed/ONE
growth pass/scramble/survive → `w_dedupe` → `w_locate_gates` → `b_armB_score` → one whole-chain
site-shift null of stages 2–5 (picks moved by the same day map `DayLoader` uses, which is
bit-equivalent to re-running the picker on shifted days).

Smoke (1 station-day, GLBC doy 259): **683 picks in 2.7 s**, 125 of them ≥0.30. Borehole B003
(EH1/EH2, unrotated; seisbench `flexible_horizontal_components` accepts Z12) 453 picks. Pilot on
2 ETS days × 14 stations: **SHVB returns 0 picks** — an independent confirmation of Merlin's F1
(the 4.76 Hz tone dominates the peak normalisation). Association cliff measured on the pilot:
thr 0.30 → 0 events at the 4-station rule, 15 at 3 stations; thr 0.15 → 144 events at 3 stations
(n_stations 113×3, 28×4, 3×5). The 3-station minimum is `e2_associate`'s documented VI rule and is
what E2's own vi2010 run used (median n_stations = 3 over 9,515 events), so it is pre-registered as
the gate with the 4-station count reported as a diagnostic. Deviations recorded in the freeze:
templates are the frozen **20-s pass-0** form (the 12-s template is defined by `w_iterate`, i.e. a
second growth pass, which "ONE growth pass" excludes); quiet-day picks are computed too (rate
contrast), and the carried threshold level is the default 0.30 unless it yields zero clusters.

**(B) WA preparation (CPU, no WA_GO, `freeze_wa_v1.json` and `x_step5_queue.sh` untouched).**
`scripts/b_wa_health.py` censused **222,936** station-hours (2010, 9,289 station-days) and
**136,560** (2008, 5,690) against Merlin's flat / H-V-outside-[1/30,30] / PSD-peak-median>1e4 rule.
**No array in either year exceeds 20 % vetoed hours** (worst 2010 UW 4.4 %, XG.DR 2.8 %; worst 2008
UW=GNW 12.2 %, YW 7.0 %), so his "the site-probe admission was blind and the coverage numbers are
overstated" conditional does NOT fire — the WA caches carry no SHVB-class station. Two stations
exceed 30 %: XG.DR09 (30.1 %, H/V) and YW.FACB (35.3 %). Coverage recomputed without them
(`b_wa_coverage.py`) moves exactly one published number: **S-Olympic 2008 94.0 % → 90.5 %**;
N-Olympic 2010 stays 94.6 %, Puget-SWA 2008 stays 84.2 %, and the `freeze_wa_v1` 11-site/145-cell
figure stays 97.2 % with 11/11 sites surviving.

`scripts/b_wa_quiet.py` replaced the {0,23} quiet hours (F2): **40 hours on 20 quiet days, all 24
clock hours present, matched to the candidate-hour histogram, earthquake-screened** (ComCat M≥2.5
within 300 km of the sub-region centre for the hour ±10 min; M≥5.5 anywhere for the hour + 2 h) —
2010 from 454 candidate hours (3 local, 43 global events → 113 slots excluded), 2008 from 113
array-energy hot hours (44 local, 254 global → 728 slots). **Finding for v2: the 2008 candidate-hour
set is diurnal** (~86 % in 15–22 UTC = 07–15 local), which is a cultural-noise signature, not a
tremor one — the hour-matched quiet set is therefore essential there and the 2008 energy-hour flag
itself should be audited before any Puget pass. All of it is written up for whoever writes v2 in
`results/WA_V2_INPUTS.md`.

## 2026-09-20 ~17:55Z — Step 1c′ built and LAUNCHED (agent H): the one permitted repair

Merlin's Step-1 audit (`results/merlin_review_step1_2026-09-20.md`, PLAN.md §11) rules that the 1c
FAIL is a failure of the frozen **seed construction**, not of array discovery or the data chain, and
prescribes exactly one repair — **pair-seeded growth** plus four defect fixes. Built and launched as
`scripts/z_step1cp_queue.sh` (pid 1165933, `nohup setsid -f`, flock `locks/step1cp.lock`,
TERM/INT/HUP/EXIT traps → `results/STEP1CP_DONE|FAILED`, every stage guarded by its own output file).

**Engine defect fixes** — all behind FLAGS whose values are recorded in the new freeze, so
`freeze_step1c.json` still reproduces 1c byte-for-byte (verified: with the flags off, the hour
2005-08-30 16h reproduces the 1c pairs to `max |Δnet_mad| = 0.0`, and the voided hour 2005-09-12 10h
still returns 0 pairs):

| flag | defect | 1c | 1c′ |
|---|---|---|---|
| `W_HEALTH_VETO` | F1 station health | none | per-station-hour veto table, honoured in the autocorrelation **and** at template stacking |
| `W_QUIET_HOURS_CSV` | F2 quiet set | clock hours {0,23}, no EQ screen | hour-matched 0–23 + ComCat/NRCan screen |
| `W_WIN_TRANSIENT=10` | burst guard | none | 10-s window with network-median energy >10× excluded as row i **and** lag j |
| `W_DROP_ZERO_MEDIAN=1` | F6 NaN void | 26/299 hours voided | flat station dropped from the hour, NaN never reaches the site sum |
| `W_SEED_RULE=pairs2`, `W_SEED_BUDGET=2000`, `W_HANDOFF_MIN_SEG=2`, `W_NEAR_SEED_MODE=either` | F3 seed rule | ≥5 members + ≥2 bins, budget 600, ≥3 segments | every ≥2-node component, ranked by max pair net_mad, budget 2000, ≥2 segments, near-seed against either end |

F7 (`y_step1diag.py` D6-native: `t0.timestamp()` on an obspy `UTCDateTime` — `timestamp` is a
property) fixed. F8 honoured: `freeze_step1c_prime.json` is a NEW file; `freeze_step1c.json` is
untouched.

**New scripts** (`z_` prefix): `z_common` (cut/stack helpers), `z_d8` (time-base / physics / chain
diagnostic), `z_health` (station-hour veto census), `z_quiet` (hour-matched + earthquake-screened
quiet set), `z_tpair` (T + mandatory leave-one-hour-out), `z_oracle` (oracle ceiling and the 2005
scoring reference), `z_null` (site-shift null ×3 counting gated families **and survivors**),
`z_score` (verdict against the oracle reference), `z_freeze`, `z_checks` (the accept/abandon gates),
`z_step1cp_queue.sh`.

**Smoke tests before launch (all passed).**
- D8 on 2 and 4 busiest kept family-hours: envelope peaks within ±3 s at 5 (2 h) and 9 (4 h) sites;
  +30 s control excess **0.983** (accept 0.9–1.1); chain ratio **1.033** median / 1.42 max — the
  frozen 20 Hz 1–8 Hz cache and the raw 100 Hz 1–8 Hz agree to ~0.001 in stack coherence, exactly
  Merlin's own measurement. Per-site median S-energy excess 0.98–1.82, p90 1.0–3.5 — the physics
  reading, not the chain. Two D8 bugs found and fixed in the smoke (the +30 s control fell outside
  the cut; the chain variant tag carried the native rate, which is 100 Hz for C8/PO and 40 Hz for
  CN and made the comparison pivot empty).
- Health census on three ETS days: **SHVB vetoed 71/72 hours** (H/V rms ratio median 890, 1–8 Hz
  PSD peak/median 7.5e7 — the 4.76 Hz tone), **GLBC 0/72, TWKB 2/48** → Merlin's F1 confirmed
  independently, and his accept clause is met.
- Quiet set (live ComCat + NRCan query, no waveform download): **40 hours, all 24 clock hours, 2 per
  day, 47 hours earthquake-vetoed — and 2005-07-09 23h is excluded**, the hour that set T = 12.50 on
  its own (M 5.8 Tateyama at 23:37 + M 5.9 Palu at 23:59; predicted P ≈ 23:48, exactly the burst).
- `z_tpair` replayed on the 1c quiet pairs reproduces T = 12.50 with LOO min 7.50 and names
  2005-07-09 23h as the dominating hour → it would have ABANDONED with that reason.
- Autocorrelation: the 1c-voided hour 2005-09-12 10h now yields **702 pairs on 12 sites** (B006
  dropped as flat), 167/71,801 windows transient-vetoed.
- Seeds/growth on 3 seeds: 3/3 handed off with ≥4 sites including a **2-node** seed (impossible in
  1c); over all 44 bed days → 141 hits, 1 survivor, the 2-node seed (`frac_near_seed` 0.20).
- Oracle seed construction: **30/30 hidden families** seeded from their busiest hour (1,092 events),
  all 10 isolated among them; the empty-chain scorer path writes a valid reference.

**Cost projection** (measured: growth ≈ 1.06 s/day + 0.0115 s per template-day; autocorrelation
≈ 2–10 s/hour): real pass ≈ 1 h, three site-shift nulls ≈ 3–4 h, autocorrelation ≈ 1.7 h over the
four passes, D8/health/quiet/oracle ≈ 1 h → **≈ 6–9 GPU-h**, matching Merlin's 8–10.

**Pre-registered gates the queue enforces** (each writes its reason into `STEP1CP_FAILED` and stops):
D8 chain ratio > ×1.5 → re-cache; envelope peaks displaced > 5 s → time base; a dense site losing
> 30 % of its hours → the bed shrinks; < 24 clock hours fillable → quiet set; an hour with zero pairs
for a non-veto reason, or the burst hour surviving the screen with pairs ≥ 8 → the fixes did not
take; T outside 6.5–8.5 or LOO range > 1.0 → the dominating hour; C_30 < 3 → the control bed is
unusable; < 500 seeds handed off with ≥ 4 sites → the pair-seeded construction has no material.
Scoring: PASS ≥ 3/10 isolated (or the C_iso < 4 fallback ≥ max(⌈0.5·C_30⌉, 5) of 30) AND split
≤ 0.10 AND null gated ≤ 1 AND null survivors ≤ 10 % of real; FAIL < 2/10 (or < 0.4·C_iso) with clean
nulls → STOP Arm A north and WA; null survivors ≥ 50 % of real → FA broken. **No `results/WA_GO` is
created by this queue and no process is killed by it** (the T7 download queue, the parked Step-5
queue and any Arm-B GPU job keep running and share the H100).
- 2026-09-20T18:00:05Z ARMB_DONE: 0/10, ARMB_INADEQUATE (5 clusters vs 3 in the null; nothing located). Arm B closed. Waiting on STEP1CP.
- 2026-09-20T18:11:14Z STEP1CP_FAILED at oracle_check: C_30 = 0 (frozen chain recovers none of 30 hidden families even when seeded with ≤40 of their own events). Orchestrator did NOT accept this as a bed verdict: it contradicts Merlin's envelope stacks (4–15 MAD at predicted S) and P4-0b truth-stack recall 0.45; whole 30-template bank gave 0–19 hits/ETS day (~1 % recall) → suspected defect (2005 time base offsets of −2…−6 s vs ±1.5 s alignment search; channel/cut). Agent H resumed for diagnosis a–c before the abandon stands.

## 2026-09-20 ~18:30Z — Step 1c′ stopped at the oracle gate; the "bed unusable" reading is WRONG (agent H)

`STEP1CP_FAILED` fired at `oracle_check`: C_30 = 0. The coordinator refused to let that stand as a
bed verdict. Diagnosed a–c inline (`z_diag_oracle.py`, `z_diag_onset.py`, `z_diag_netmad.py`); full
write-up `results/STEP1CP_RESULTS.md`, tables `results/step1c_prime/evidence/`.

- **(a) Time base is off but is not the cause.** 77.5 % of 40 families (30 hidden + 10 kept
  controls) have envelope-peak offsets > 1.5 s at ≥ half their sites; median-of-median 2.79 s, p90
  6.04 s — Merlin's family-45 observation generalises. **Widening `align_stack` does nothing**
  (stack S-SNR 0.855 dB at ± 1.5 s → 0.952 dB at ± 6 s; peak height 3.39 → 3.46 MAD), because
  align_stack aligns events to each other while the offset is a per-family shift against the
  *prediction*. The alignment window was therefore NOT widened.
- **(b) Templates are correct.** Four cut rules (frozen blind argmax, ± 4 s, exact predicted S,
  ± 10 s) give the same mean recall 0.010–0.014 — the onset picker is not the bottleneck. Channels,
  rotation, `t0` and template start are sound; strong templates detect 38–68 % of their own seed
  events at purity 0.70–0.96.
- **(c) Recall is ~1 % over three days but 24 % on the family's own day.** Sampling the detector's
  full `net_mad` series at every catalogued event: median 2.46–3.12 vs null 1.88–2.05, and
  **11–43 % (median 24 %) of events pass ALL detection clauses**. The `n_site_z3 ≥ 2` coincidence
  clause costs nothing (`frac_ALL` = `frac_nm5`).
- **(d) A burst guard applied to a seed chosen as a burst.** 21/28 oracle templates were vetoed by
  `frac_near_seed > 0.20` while carrying 20–46 hits on 2–8 days at scramble ratios up to 17 —
  Merlin's protocol seeds from the family's *busiest hour*. Fixed as a documented defect
  (`z_oracle.oracle_survival`, `results/freeze_step1c_prime_v2.json`): for the oracle only the guard
  is reported not applied, and hits/days are counted **outside** the seed window — stricter than the
  frozen rule, and it cannot manufacture a PASS (a higher ceiling only raises the bar). Result is
  still 0/28: the binding clause moves to `n_indep ≥ 20` (1–13 outside the seed window).
- **(e) The decisive control.** The same templates rebuilt from **400** events spread over the window
  give **28–125 independent hits on 9–20 distinct days for 4 of 6 families** (recall 1.6–8.1 %,
  purity 0.28–0.60; 2 over-trigger at purity 0.01). **The bed supports survival — a ≤ 40-event
  pass-0 stack does not.**

**Conclusion.** The limit is the survival gate sitting at pass 0, *before* the iterations that would
enrich the stack — not the data, and not the seed rule Merlin's repair changed. Since the 1c′ pass-0
templates come from 2-node pair seeds (thinner than the 40-event oracle), stages 7–10 as frozen
would return ~0 survivors for a reason unrelated to the repair. Moving or staging the survival gate
is a chain change = a second repair, which Merlin's ruling forbids, so the queue was **left stopped**
rather than relaunched. Stages 1–5 are complete and guarded: a decision needs only stage 6 onward
re-run (~6–8 GPU-h). Question put to Merlin in the write-up: is the correct ceiling protocol a richer
oracle stack, and is the pass-0 survival gate the defect F3 was pointing at?
- 2026-09-20T18:30:03Z Agent H diagnosis: bed sound, survival gate starved on 24 days (400-event stacks survive, ≤40-event stacks do not). Merlin consulted for the go/no-go.
- 2026-09-20T18:52:01Z Merlin ruling saved; PLAN §12 (final spec). Launching agent J (fresh) to run Steps 1–5 as one chained queue with mechanical abandon/STOP markers.
- 2026-09-20T19:06:32Z MANDATE CHANGE (user): reproduce Bostock on his stations, debug until it works; hidden 30 sealed for one final blind check; Merlin only at the end. Agent K launched (REPRO_BOSTOCK.md ledger).
- 2026-09-20T19:23:15Z User: Bostock is definitely reproducible; we likely imposed requirements/bars he never used. Agent K instructed to implement Brown 2008 / Bostock 2012 / Royer & Bostock 2014 LITERALLY (fetch the method sections; parameter table), strip every invented clause (FA-only quiet-day T, ≥20-indep/≥3-day survival, scramble ratio, burst guard, ≥4-site, per-hour renormalisation, 12-s S−2 cut) from the development loop, and measure recovery of the 100 kept families under his recipe before changing anything.
- 2026-09-20T19:31:58Z Agent K: literal Bostock recipe gives 35 % true pairs (top-100 78 %) in 2 h; our added clauses were the blocker; 2005 sub-second times recovered. k_repro_queue pid 1180966 running.
- 2026-09-20T19:58:12Z Agent K: 88/96 kept families in the pair stage; 14/19 seeded families rediscovered in the growth control; 4 defects fixed (gap windows, cc vs nmad ranking, 1/sqrt(0), 0.5-s row-grid jitter → Bungum–Husebye realignment). Full grow run pid 1188403.
- 2026-09-20T20:44:03Z K_REPRO2 DONE: pair stage 88/96 families; pure groups (≥0.5) 75; ≥0.8: 58; grown 1,010 → 89 after >10 % shared-detection merge; **9 distinct kept families rediscovered** at ≥30 %; 67 lost between grouping and growth, 15 at grouping. Median recall_all of grown templates 0.02. Agent K resumed to diagnose the grouping→growth loss (merge collapse 1,010→89; group size; iterations; realignment; threshold).
- 2026-09-20T20:55:48Z Agent K: 50/96 kept families rediscovered pre-merge; merge was the artefact (no waveform merge; location later); 15-station run + 6×MAD first pass + ceiling running pid 1196829.
- 2026-09-20T21:56:26Z DL_ALL_DONE: agent D's download queue complete (T1–T7 incl. the 2007–2026 extension record for 17 Olympic permanents).
- 2026-09-20T22:36:45Z K_REPRO3 DONE (15 stations): pair stage 92/96; pure groups 77; ≥0.8 60; **58/96 kept families rediscovered**; lost at growth 23, at grouping 15.
- 2026-09-20T22:43:28Z K ledger §8: ceiling 89/96 (recall_all 0.31); 58 of 89 = 65 %; two-stage threshold refuted; locator: 29 unmatched candidates on the interface with ≥3 days. K_REPRO4 pid 1205708 running (permanent-bed recipe, blind hidden check, classification).
- 2026-09-20T23:34:42Z K_BLIND_DONE (provisional, under audit): hidden 0/30 (log says 'of 11'), isolated 0/10 vs kept 58/96 — implausible for a blind recipe → scorer audit ordered (corrected times exist only for kept families; family-key join; raw whole-second re-score). Permanent-bed lines (33/77 and 3) awaiting clarification.
- 2026-09-20T23:39:56Z BLIND 22/30 hidden, 6/10 isolated (scorer env-leak fixed; first true blind measurement). Permanent 2012 bed: discovery 3/77 vs ceiling 33/77; pair purity 1.7× chance vs 11× on arrays. K_REPRO5 pid 1213730 (classification vs correct hidden score; 2010 bed).
- 2026-09-20T23:41:13Z Killed the obsolete parked x_step5_queue. Agent L launched: N-Olympic 2010 with the reproduced Bostock recipe + location-based verification.
- 2026-09-20T23:58:13Z North final: hidden 22/30 (ceiling 29); 18 candidate-new → 7 conservative → 1 spatially distinct (grp 1215). NORTH_DISCOVERY_SUMMARY.md written. K_REPRO6 (2010 permanent bed, 8 stations) running.
- 2026-09-21T00:01:09Z Agent L: WA N-Olympic 2010 queue pid 1222985 (~6 h): 96 3-C stations (XG 64, XU 19, PB 10, UW 2, CN 1; DR09 excluded), 454 PNSN hours / 29 tremor days, 40 quiet hours; 8×MAD = 4.5–5.3 sd; smoke hour 14,277 detections in 3 s GPU; PB borehole azimuths added for rotation; station terms usable for 19/96; pass-0 templates only; site-shift null through the identical chain; first real hour 21,754 detections at 4.9 sd.
- 2026-09-21T00:13:44Z K_REPRO6: 2010 permanent bed discovery 5/60 vs ceiling 50/60 (added to NORTH summary + ledger §11.6).
- 2026-09-21T03:51:40Z Resume file written: notes/RESUME_2026-09-21_lfe_array_discovery.md (user approaching weekly limit; may continue from another account).

---

## 2026-09-20 — agent L: Bostock's recipe on the Washington N-Olympic 2010 array bed

One detached queue (`scripts/l_wa2010_queue.sh`).  Bed: 96 probe-admitted 3-C stations (CN 1, PB 10, UW 2, XG 64, XU 19), 2010-07-20…2010-09-10, 454 PNSN tremor hours, 40 hour-matched quiet hours.

* pair stage: 2,582,086 detections above 8 × MAD (4.90 sd), tremor/quiet pair-count ratio **4.18×**;
* grouping: 1,705 components with ≥3 members; grown 1,500;
* **50 candidate families** against **17** from the whole-chain site-shift null;
* no catalogue truth exists here — nothing is a rediscovery count.  Full ledger `results/WA_2010_RESULTS.md`.

Engine reused unedited (`k_auto`/`k_group`/`k_grow`/`k_locate`); the WA bed is supplied by `l_common.py` (station table from `channels_2006_2012.csv`, PB borehole azimuths added to `s_common.AZ12` so the boreholes are rotated instead of dropped, net-qualified cache index because C8.PA02… shadows XG.PA02…, station-hour health mask, day-shift null).
- 2026-09-21T07:32:40Z L_WA2010_DONE: 50 vs 17 raw; grid-edge filter → interior 4 vs 3; strict 2(+1) vs 0 (grp 1656, 1193, 1152). Launching agent M: sealed-2011 re-fire + wider-grid relocation.

## 2026-09-21 — agent M: the 2011 SEAL IS OPEN; confirmation re-fire launched

**Seal opened 2026-09-21T07:48:19Z** (`results/SEALED_2011_OPENED_2026-09-21T074819Z`, reason
recorded, note appended to `results/SEALED_2011_DO_NOT_READ`).  Integrity: **5,062/5,062 manifest
entries sha256-clean** in 67 s (`scripts/m_seal_open.py`); 317 of them are symlinks into
`lfe_extension/traces` and were followed — `find -type f` misses those, which is why a naive
count says 4,745.  `traces_sealed/` is read-only from here on; it is cached into
`cache20hz/<NET.STA>/2011/<DOY>.npz` and never modified.

Launched as ONE detached queue `scripts/m_refire2011_queue.sh` (flock `locks/m_refire2011.lock`,
TERM/INT/HUP/EXIT traps → `results/M_REFIRE2011_FAILED`, per-stage output guards, `env VAR=…
python` per stage, markers `results/M_REFIRE2011_DONE|FAILED`, log `logs/m_refire2011_queue.log`).
Stages: seal → cache (4,398 station-days, T3, frozen chain + rational twin for the 50 Hz XG/XU) →
`b_wa_health.py --year 2011` → bed → re-fire → wide-grid re-location → report.

**Nothing is tuned.**  `k_grow.matched_filter` is called verbatim on the template waveforms read
straight out of `results/l_wa2010/grow/tpl_it0.npz` (and the null's); a 2010 template station that
does not exist in 2011 is simply dropped and the survival fraction is recorded per template.

2011 bed: **90 stations** = the 96 of 2010 minus XG.DR02, XG.LC11, XG.LC13, XG.LC14, XG.TB14,
XU.W030 (absent in 2011), before the 2011 station-level health bar.  **26 PNSN tremor days**
(429 candidate hours, 2011-08-07…09-01) — the identical code reproduces the 2010 run's 454 hours /
29 days exactly, so the hour rule is unchanged.  **Only 7 strictly-quiet days** exist in the sealed
window (no PNSN tremor within 100 km of the region centre, ≥3 days away): 07-22…07-24 and
09-07…09-10.  The brief asked for 20; no more data may be downloaded, so this is a **documented
bed-forced deviation** — the false-alarm reference is run over all 27 non-tremor days and the
pre-stated ≤0.1/template-day quiet clause is evaluated on the 7 strict ones (with 7 days that
clause means *zero* quiet-day detections).

CONFIRMATION RULE, fixed before the seal was opened and applied mechanically to the 50 real and
the 17 null candidate templates alike (the null count = the false-confirmation rate):
≥3 distinct 2011 tremor days with detections **and** quiet-day rate ≤0.1/template-day **and** a
2011 location that is interior to the wider grid and within 10 km (3-D) of the 2010 location.

Wider relocation grid **46.8–49.0 N / −125.2…−121.8 / 15–60 km** at the 2010 spacing (0.03°,
0.04°, 2 km).  Both years are re-stacked from their OWN detections with `k_grow.build_templates`
(unchanged) before locating — locating a template file twice gives the same answer by
construction, so only a stack built at the 2011 detection times can answer "does it locate in the
same place in 2011?".  The 2010 null chain is re-stacked with `L_NULL_SHIFT=1` on its own shifted
network; both template sets are re-fired on the SAME real 2011 data.

Smoke test before launch (foreground): 121 station-days cache in 8 s; grp 1656/1193/1152 over six
2011 days → 30 detections, all three templates fire, ~1.7 det/template-day; the wide-grid
re-stack+relocate path runs end-to-end.

New scripts: `m_seal_open.py`, `m_cache2011.py`, `m_common.py` (2011 bed adapter over
`l_common`), `m_bed2011.py`, `m_refire.py`, `m_relocate.py`, `m_report.py`,
`m_refire2011_queue.sh`.  Result → `results/WA_2011_REFIRE_RESULTS.md`, tables
`results/m_2011/confirm_{real,null}.csv`.
- 2026-09-21T07:55:10Z Agent M: seal opened 2026-09-21T07:48:19Z (manifest 5,062/5,062 clean); 4,398 station-days cached; 90/96 stations exist in 2011; 26 tremor days (429 h); only 7 catalog-quiet days inside the window (FA also over all 27 non-tremor days); smoke: 1656/1193/1152 all fire in 2011 (~1.7 det/template-day over 6 days). Queue pid 1254314, ~2.5 h, markers M_REFIRE2011_DONE|FAILED.

## 2026-09-21 — agent M: sealed-2011 re-fire (confirmation of the WA 2010 candidates)
Seal opened and recorded (SEALED_2011_OPENED_2026-09-21T074819Z); 3052 station-days cached; 89 stations (of the 96 of 2010); 26 tremor days, 7 strictly-quiet days (bed-forced deviation from the requested 20 — documented in `results/WA_2011_REFIRE_RESULTS.md` §1).  Unchanged templates, unchanged 8×MAD matched filter, wider relocation grid 46.8–49.0 / −125.2…−121.8 / 15–60 km.
**Result: REAL 0/50 confirmed, NULL 0/17 confirmed (false-confirmation rate 0.00).**  Interior 2010 candidates: grp 429 not confirmed, grp 1152 not confirmed, grp 1193 not confirmed, grp 1656 not confirmed.
Full write-up `results/WA_2011_REFIRE_RESULTS.md`; tables `results/m_2011/confirm_{real,null}.csv`.
- 2026-09-21T09:10:15Z M_REFIRE2011_DONE: CONFIRMED 0/50 real, 0/17 null; real templates fire 3× more on non-tremor days than tremor days → site/cultural repeaters, not LFEs. WA 2010 N-Olympic = no family. Next: Puget/S-Olympic 2008 CAFE bed with a per-station tremor-hour audit and template-level tremor/quiet screens; 2007-01 as its confirmation bed.
- 2026-09-21T09:39:33Z Agent N: Puget/S-Olympic 2008 queue pid 1271745 (~2–3 h): 51 stations (PB 10, XU 40, YW 1); hour audit → 115 tremor hours on 20 days (05-05..05-26, 0.9 % daytime; old flag 86 %), 0/456 quiet hours fire; tremor/quiet pair ratio 8.2× (2010: 4.2×); wide grid + edge flag from the start; 5-clause screen real vs null; 2007-01 confirmation chained.
- 2026-09-21T09:53:46Z User restated the goal: reproduce Royer & Bostock 2014 / Chestler & Creager 2017 WA families. Agent O launched: 2010 Olympic re-run with per-array stacking (arrays as arrays), Chestler check, 2011 re-fire.

## 2026-09-21 — agent N: Puget / S-Olympic 2008 (CAFE XU + PB), the last full-array WA bed
**Tremor-hour audit replacing the 86 %-daytime energy flag** (per-station diurnal-normalised
2–8 Hz excess; ≥ 60 % of the near-field array simultaneously above its OWN quiet-day p95; plus a
spatial-coherence test — the excess must fit an amplitude decay from one active cell better than a
uniform field, p ≤ 0.05 against 200 station-label shuffles): **115 hours on
20 days**, only 0.009 of them in 15–22 UTC
(the old flag: 0.86), weekday fraction 0.774 (5/7 = 0.714 by calendar);
**0 of 456 literature-quiet-day hours fire**.  The selected days
reconstruct the published May-2008 Washington ETS from the waveforms alone.  The array's diurnal
noise floor is ~10× higher at 19 UTC than at 09 UTC, so the night-heavy clock distribution is
detection-limited, not a claim that tremor is nocturnal (quiet-day false-positive rate by UTC hour
is flat).  Bed: 51 3-C stations {'PB': 10, 'XU': 40, 'YW': 1}, 19
literature-quiet days.
Pair stage → grouping → growth → the five-clause screen (per-template quiet-day FA; tremor
co-timing; wide-grid INTERIOR location + shuffled-station control + Slab2; ≥ 3 days; second-half
re-stack stability), real against a whole-chain site-shift null through byte-identical code:
**CANDIDATES REAL 0 vs NULL
0**.
2007-01 confirmation re-fire (tranche T5, read by no discovery step of this project):
**REAL 0/0 confirmed, NULL
0/0** (2007 audit: 16 hours on
9 days).
Full write-up `results/WA_2008_RESULTS.md`; new scripts `n_common.py`, `n_hours.py`, `n_auto.py`,
`n_pairstats.py`, `n_run.py`, `n_chunk.py`, `n_locate.py`, `n_screen.py`, `n_relocate.py`,
`n_refire.py`, `n_confirm.py`, `n_report.py`, `n_notes.py`, `n_puget2008_queue.sh`.
- 2026-09-21T10:15:26Z N_PUGET2008_DONE.
- 2026-09-21T10:16:45Z N_PUGET2008_DONE: 0/0 candidates; pair ratio 0.89 over the full run; grown templates fire equally on quiet and tremor days (noise repeaters). Launching agent P: per-station noise-repeater census + quiet-subset 2008 bed.
- 2026-09-21T10:29:32Z Agent O: 2010 Olympic arrays-as-arrays queue pid 1297832 (~7–8 h; markers O_OLY2010_DONE|FAILED; results/WA_2010_ARRAY_RESULTS.md): 96 stations → 22 sites (10 beams incl. XU BS/PL mini-arrays and PB clusters + 12 singles); Chestler's method fetched from her thesis (MBBP beams 5–9 Hz vertical, 18×MAD peaks, 3 stations/array 15-s 2–8 Hz, 2nd-root stack 45×MAD, 150-window templates ×3, R_com>0.5 dedupe) and run literally as Arm B (non-fatal); slowness-grid beam rejected by its quiet-hour control (max-over-grid bias) → zero-slowness stack; array coherence at LFE incidence 1.9 (1–2 Hz) → 1.06 (6–8 Hz); queue picks member/beam/siteavg by tremor/quiet ratio on a common station set; then screen real vs null, Chestler check both directions, 2011 re-fire.
- 2026-09-21T10:30:12Z Agent P: Puget 2008 quiet-subset queue pid 1298829 (~1.5–2 h; markers P_PUGET2008Q_DONE|FAILED; results/WA_2008Q_RESULTS.md). Census smoke: CAFE surface sites ~14 dB louder than PB boreholes on quiet days; worst XU.BD03 ratio 0.13. Provisional 16-station quiet subset: pair tremor/quiet ratio 10.8× (vs 0.89× on 51 stations; north 11×) — caveat: fewer stations in quiet blocks. PB-only covers 7.2 % of cells at ≥5 stations → screen counts are a lower bound.
- 2026-09-21T10:33:17Z User: quiet-hour logic was not in the papers — keep a criterion only if justified. Clarified: the recipe has no quiet-hour criterion (8×MAD of the block's own distribution); quiet hours are (1) a diagnostic ratio and (2) screen clause (a). Agent Q launched: apply every screen clause to the 58 rediscovered Bostock families (must pass ≥85 %) and to the north noise population; drop unjustified clauses, reformulate borderline ones data-driven; re-screen WA 2008/2010 with the justified set.

## 2026-09-21 — agent P: Puget / S-Olympic 2008, the QUIET-SUBSET bed (station quality as the lever)
Agent N's 51-station CAFE+PB run gave a pair-stage tremor/quiet ratio of
0.89 and grown templates firing equally on quiet and tremor days —
the pair stage was dominated by repeating cultural/site noise.  Royer & Bostock (2014) built
COMPOSITE ARRAYS out of chosen stations, so this run tests station QUALITY, not count.
**Per-station noise-repeater census** (`scripts/p_census.py`): each station's own 3-C mean cc
recomputed at exactly agent N's stored (i, j) pairs and compared with its own median + 8 × MAD
over 4,000 random pairs of the same block (validated — the mean over stations reproduces the
stored network cc to 7e-8).  8 of
51 stations are
noise repeaters (tremor/quiet contribution ratio < 2 AND a top-3 share on quiet hours above
twice chance).
**Quiet subset** (rule fixed before any subset result: the PB boreholes + CAFE/YW with ratio ≥ 3
and low top-3 share, floor 12 / cap 20): **20 stations**
(10 PB + 10 CAFE/YW);
coverage 53.0 % of the
sub-region's active cells with ≥ 5 subset stations within 60 km hypocentral
(51-station bed: 90.5 %).
**Pair stage, unchanged recipe:** 51 stations 0.89× →
**quiet subset 0.86×** → boreholes only
0.45× (north POLARIS 2005: 11×; N-Olympic 2010: 4.2×).
Caveat recorded in the write-up: the quiet blocks run on fewer stations than the tremor blocks
(CAFE caching on calendar-remote days), which raises the absolute 8 × MAD bar on quiet hours.
**Five-clause screen (agent N's, unchanged) real vs whole-chain site-shift null: CANDIDATES
REAL 0 vs NULL
0**; network false-alarm rate
11.228 per template-day on quiet days against
12.526 on tremor days.
**2007-01 confirmation:** REAL 0/0,
NULL 0/0.
Bed as run: 20 stations {'PB': 10, 'XU': 10},
115 audited tremor hours, 19
literature-quiet days.  Write-up `results/WA_2008Q_RESULTS.md`; scripts `p_census.py`,
`p_subset.py`, `p_report.py`, `p_notes.py`, `p_puget2008q_queue.sh`; the only change to an
existing file is the optional `N_STATION_SUBSET` whitelist in `n_common.admitted()` (unset →
agent N's numbers reproduce byte-for-byte).
- 2026-09-21T10:52:00Z P_PUGET2008Q_DONE.
- 2026-09-21T10:53:04Z P_PUGET2008Q_DONE: quiet subset ratio 0.86, boreholes-only 0.45 → station quality is not the lever; 0 candidates. Fetching Royer's UBC thesis for the WA family facts.
- 2026-09-21T10:54:22Z Literature check for Royer & Bostock 2014 WA families: abstract only (4 composite arrays on VI + WA; EarthScope TA/FA/PBO + POLARIS + permanents; depths 29–46 km); ScienceDirect/Wiley/ResearchGate/UBC cIRcle/Semantic Scholar all blocked or abstract-less; Bostock et al. 2019 review confirms WA LFEs (Royer & Bostock 2014) but gives no family count/stations. WA family count, region and station set remain unknown → we cannot say how many families the 2007/2008 CAFE bed 'should' yield.
- 2026-09-21T11:05Z **Agent Q — SCREEN VALIDATION ON THE NORTH BED (`results/SCREEN_VALIDATION_NORTH.md`, workspace `results/q_screenval/`, queue `scripts/q_sv_queue.sh` → `Q_SV_DONE`).** The five-clause WA screen (`n_screen.py`) was applied, exactly as coded, to the grown templates of the **58 Bostock families the same recipe rediscovered blind** on the 2005 POLARIS bed, to the 64 templates of the 18 located sealed families, the 18 candidate-new, and the 1,201 unlocatable/unmatched (noise). Two missing measurements were made first: the unchanged K_REPRO3 templates re-fired at 8×MAD over the **20 literature-quiet 2005 days** (43,224 detections, `det_quiet15.parquet`) and **first/second-half re-stacks + re-locations** of every template (`reloc_h1/h2.csv`).
  **Result: the screen as coded keeps 2 of the 58 (3 %); without clause (e), 7 (12 %).** Per-clause pass fraction among the 58 true families / noise rejection: (c1) located ≥5 stations **86 % / 67 % → the only JUSTIFIED clause**; (d) ≥3 days 90 % but rejects only 20 % of noise; (a) quiet FA 84 % / 45 %; (b) co-timing 57 % (81 % with a base-rate-matched hour set) ; (c2) shuffled control 53 %; (c4) |Δz|≤8 km 48 %; (c3) interior 38 %; **(e) 2nd-half re-stack 22 %** — a true family's re-stacked location moves a median **35.5 km** (h1-vs-h2 33.9 km), so (e) is not a property Bostock's own families have.
  **Also: the quiet-day network false-alarm rate that condemned the 2008 bed is 1.480/tpl-day on the NORTH bed vs 1.854 on tremor days (ratio 1.25), against 1.525 vs 1.541 (1.01) at Puget 2008** — the diagnosis does not separate the bed where the recipe works from the one where it does not (only the 2011 re-fire's 0.35 does).
  Data-driven reformulations measured: |Δz| width **19 km** (85th pct of known families; still rejects 77 % of noise); quiet-FA cap **0.25/tpl-day**; **the day-label-shuffle replacement for (a) proposed in the brief is much WORSE (22 % of true families) — not adopted**; (c3) and (e) dropped (grid-adequacy diagnostic and split-half location uncertainty instead); (c2) and (b′) demoted to tiers.
  **Recommended screen** tier A = (c1)+(|Δz|≤19 km)+(≥3 days)+(FA≤0.25/tpl-day, ≥3×); tier B = +control; tier C = +base-rate-matched co-timing. North yield 33/21/13 of 58 families vs 11.0/1.3/0.8 % of the noise population.
  **Re-screened WA: 2008 (agent N) 0/0 as coded → tier A 17 real vs 2 null, tier B/C 4 real vs 0 null. 2010 (agent L) 50 vs 17 (2.9×) as published → tier B 234 vs 59 (4.0×), tier C 127 vs 27 (4.7×)**; on a 96-station bed the location gates are free, so the thresholds must be re-derived per bed, and **clause (a) has never been run on the 2010 bed — the 2010 templates should be re-fired on the 40 hour-matched quiet hours.**
  Literature (Bostock 2012 itself is genuinely closed — Unpaywall/OpenAlex/S2 all `closed`, UBC thesis bot-blocked; quotes are from later papers of the same group): the ONE family-acceptance rule positively attested for Bostock 2012 is a **dedupe** — Royer, Thomas & Bostock (2015): the 2012 SVI catalogue "was culled to **93 geographically independent LFE families**" (130 templates → 93). Bostock 2015's only stated culling is the 4-s declustering; its cc>0.4 / >70 %-energy / cc>0.25 rules are per-CHANNEL, not per family. Plourde 2015: dedupe at >10 % shared detections or cc>0.3, plus a de facto "sufficiently clear to be picked" gate (59 templates → 34 families). **Thomas & Bostock 2015 (Tectonophysics): ≥100 detections in a 3-day window, >10 %-shared-times dedupe, HYPOINVERSE location — but they KEPT 5 of 18 families they could not locate, and their own quiet-period estimate is "2–3 false detections per day" for HIGH-QUALITY templates, i.e. 20–30× our 0.1/day cap, treated as acceptable.** Chestler 2017 §2.1–2.2: 18×MAD beam peaks → 2nd-root stack at 45×MAD → 150-window template ×3 iterations, the one qualitative gate being "If the template has clear P- and S-wave arrivals…", 3-D grid search + double-difference location, R_com>0.5 dedupe (57→43). **No quiet-day false-alarm gate, no minimum-day rule, no Slab2 agreement, no interior-to-grid clause and no split-half stability test appears anywhere.** Six of our eight coded conditions are project inventions; the two most destructive, (c3) and (e), have no published counterpart — while the one rule that IS attested (the dedupe) is the one we do not apply.
- 2026-09-21T11:03:40Z Agent Q: screen validation — as coded it keeps 2/58 real Bostock families; justified tier A/B/C defined; WA 2008 → 17 vs 2 (tier A), 4 vs 0 (B/C); WA 2010 → 234 vs 59 (B), 127 vs 27 (C). Confirmation rules also unjustified. Launching agent R: calibrate confirmation on the north, re-confirm 2008 on 2007 and 2010 on 2011 with the justified tiers, Chestler check.
- 2026-09-21T11:15:58Z BOOKKEEPING at 98 % weekly limit: resume file §10 finalised (running jobs, established numbers, next actions); memory snapshot refreshed. No further agents to be launched this session.
- 2026-09-21T11:16:17Z Agent R stopped before launch (user: no agents; new agent reads results). Only o_oly2010 (pid 1297832) runs, unattended.

## 2026-09-21 -- agent O: the 2010 N-Olympic bed run with the arrays AS ARRAYS
Goal (PLAN.md sec.17): reproduce Chestler & Creager 2017 / Sweet et al.  Sites defined by geometry:
96 stations -> 22 sites (10 arrays
covering 84 stations + 12 singles).
Beam = zero-slowness delay-and-sum; the Ghosh slowness grid was measured to give the SAME semblance
gain on quiet hours as on tremor hours (max-over-grid bias), and intraslab-earthquake S waves show the
1-km arrays are only ~25 % coherent at 1-2 Hz and incoherent above 4 Hz.
Site treatment chosen by measurement: siteavg ([{'mode': 'member', 'n_blocks_tremor': 12, 'det_per_hour_tremor': 4349.75, 'traces_tremor': 96.0, 'thr_sd_tremor': 4.9169816075540425, 'n_blocks_quiet': 12, 'det_per_hour_quiet': 3139.75, 'traces_quiet': 78.0, 'thr_sd_quiet': 5.220667853711181, 'tremor_quiet_ratio': 1.3853810016721078}, {'mode': 'beam', 'n_blocks_tremor': 12, 'det_per_hour_tremor': 130.16666666666666, 'traces_tremor': 21.0, 'thr_sd_tremor': 5.167747842044173, 'n_blocks_quiet': 12, 'det_per_hour_quiet': 131.41666666666666, 'traces_quiet': 21.0, 'thr_sd_quiet': 5.248150965867303, 'tremor_quiet_ratio': 0.9904882688649335}, {'mode': 'siteavg', 'n_blocks_tremor': 12, 'det_per_hour_tremor': 8633.25, 'traces_tremor': 96.0, 'thr_sd_tremor': 5.1128297081545035, 'n_blocks_quiet': 12, 'det_per_hour_quiet': 545.9166666666666, 'traces_quiet': 78.0, 'thr_sd_quiet': 5.217725969313501, 'tremor_quiet_ratio': 15.81422683559762}]).
Pair stage: tremor/quiet ratio 8.224534056252118 (2010 member-wise run 4.2x, north 11x).
Screen: real candidates 0 vs null 0.
Chestler check: {"n_screened_located": 0, "median_km_3d": null, "min_km_3d": null, "random_median_km": null, "random_p05_km": null, "frac_random_at_least_as_close": null, "n_within_5km": 0, "n_within_10km": 0, "n_within_20km": 0} / {"n_chestler_coords": 69, "n_interior_located_real": 113, "n_interior_located_null": 67, "n_coords_hit_10km_real": 18, "n_coords_hit_5km_real": 5, "n_coords_hit_10km_null": 14, "n_coords_hit_5km_null": 0, "median_nearest_km_real": 13.571941731615428, "median_nearest_km_null": 14.914688279822334}.
2011 confirmation: real 0 / null 0.
Arm B (Chestler's own chain): {"n_seeds": 756, "n_keep": 739}.
Full ledger results/WA_2010_ARRAY_RESULTS.md; workspace results/o_oly2010/.
- 2026-09-23T00:21:01Z New session/account (Sonnet 5 orchestrating). O_OLY2010 had finished 2026-09-21T18:23Z, unread since. Key unprocessed finding: Arm B (Chestler's literal chain) produced 739 dedup'd families but was never located/screened/null-tested. Old screen's Chestler check (before tier fix) found 5/69 published coordinates with an interior-located template within 5 km vs 0/69 for the null -- promising but on the wrong (unjustified) population. Agent R's north confirmation-calibration was never run. Launching agent S to do all of: locate+screen Arm B output, calibrate confirmation on north (agent R's brief), correctly re-screen+re-confirm WA 2010 (tier A/B/C) and WA 2008 (tier A) via 2011/2007, redo Chestler check on the properly-screened candidate sets with a matched null.
- 2026-09-23T00:46:50Z Agent S launched s_final_queue.sh pid 1400449 (~8-10 h; markers S_FINAL_DONE|FAILED, interim S_FINAL_STAGE12_DONE). Stage-1 calibration running (354 templates: 177 true from 58 families + 177 matched noise; only 128/354 have >=3 H1 detections - recorded as a property of any half-window protocol). Dedupe smoke on 2008 real: 160 grown -> 129 after literature dedupe (70 merged on >10% shared detections, 0 on stack cc; max pairwise stack cc 0.156 -> Plourde's cc>0.3 separates nothing here, matching the north) -> 125 located -> 41 accepted at 15-58 km (16 edge-railing, reported not rejected). Arm-B locate crash fixed; Arm-B null patch = 22 site offsets across 96 members; 2011 Arm-B re-fire uses 45xMAD/2nd-root, not 8xMAD. Tier A/B/C on agent O's array chain folded in as a comparison row.
- 2026-09-23T01:00:38Z USER CORRECTION (binding): the 85%-of-known-families test is a VETO on proposed cuts, NOT a licence — calibrating cuts to known families cannot protect undiscovered ones that are smaller/deeper/noisier. NEW RULE: rank, do not discard. Every grown family keeps its row; every criterion becomes a column; the null gets the identical table; evidence = real-vs-null separation curve + ranking, not a pre-chosen threshold. Only the papers' own dedupe and location may reduce a list, and both are reported before/after. Memory: rank-dont-discard.md. Pushed to both running agents.
- 2026-09-23 ~01:20Z **RANK-DON'T-DISCARD REBUILD (agent S, queue pid 1408435).** Queue stopped mid-stage-2 (kill -TERM on its own pid) and redesigned: no filtering anywhere. Every grown template of every chain keeps a permanent row in results/sf_final/table_<bed>_<chain>.csv (sf_table.py); every criterion is a COLUMN; the site-shift null gets the byte-identical table; evidence = sf_separation.py real-vs-null counts swept across thresholds per column, AUC +- SE, lift, ranked by strongest discriminator. No thresholds, no candidate count. Only dedupe and location reduce a list and both retain flagged rows.
  **Depth window measured and REJECTED as a cut:** literal 15-58 km is vacuous on the north (locator grid is 18-56 km); its transferable content ("reject a depth-rail solution") keeps only 33/50 = 66% of located true families, below the 85% veto bar; and 100% of its rejections across all four WA chains are the bottom grid rail (dep=59 km), 0 are genuinely interior -- a grid-rail clause with a physical label. The >=5-station bar passes the veto (50/58 = 86%) and stays a reported flag.
  **Stage 1 calibration (information only now):** rule = >=1 tremor day + quiet FA <=1.0/tpl-day; the 10-km location-move clause used by EVERY previous WA confirmation cannot be calibrated at all. Only 18/58 families are testable on a half-window; the calibrated rule confirms 16/18 true families but ALSO 60/84 noise templates (71%) -> sealed-window re-fire is a weak discriminator, now columns not a gate.
  **Dedupe collapse = the strongest real-vs-null contrast found so far** (from the stage-2 tables already on disk): Olympic 2010 real 1485 grown -> 780 kept from 11,668 overlapping pairs; null 1486 grown -> 1247 kept from 349 pairs (33x more overlap on the real network); max template-pair cc real 0.644 vs null 0.475. Puget 2008 real 160 -> 129 (70 pairs, max cc 0.156); null 36 -> 35 (1 pair). CONTROL STILL OWED: whether the real bed simply has more total detections, which would inflate sharing on its own -- normalise the overlap rate by detection count before trusting it.
  **Smoke separation:** 2008 real vs null separates on detection count (AUC 0.72 +- 0.07, lift 3.9x at >=8 detections) but REVERSES on median cc (AUC 0.27 -- null templates more self-similar). Reverse separations reported explicitly.

## 2026-09-23 — agent U: the LOCATION step rebuilt and calibrated on Bostock's published coordinates
Queue `scripts/u_locate_queue.sh` (pid in `logs/u_locate_queue.log`; markers `results/U_LOCATE_BASELINE_DONE`
then `U_LOCATE_DONE|FAILED`); code `scripts/uu_{common,stacks,picks,locate,resid,dd,eval,adopt,wa,report}.py`;
output `results/u_locate/` and `results/WA_LOCATION_REBUILD.md`.  Detection and growth untouched.

**Baseline, measured for the first time** (58 blind-rediscovered Bostock families, representative template
per family, truth = the published family coordinates in `lfe_catalog_repair/results/e5r_families_2012.csv`,
verified to be Bostock 2015's own catalogue triples): the current locator (`k_locate.py` on `tpl_it0`)
locates 44/58, **median epicentral error 27.6 km, p90 87.9 km, median |depth| error 15.9 km (bias -11.7 km),
median 3-D 35.2 km, 75 % of solutions railing on a grid face**; the chance level (same solutions, family
labels permuted) is 39.8 km 3-D.  Only 11 of 58 land within 10 km.

**Diagnosis.** With the published hypocentre fixed the S pick residual is BIMODAL: at stack envelope peaks
>= 8 MAD the observed time matches the Crosson 1-D prediction to 0.15 s (90 % within 1 s) -- the travel times
are right -- but over all admitted stations 22 % of picks are > 4 s out, and an unweighted least-squares fit
lets those drag the solution onto a grid face.  Station medians span only -0.8..+1.9 s.

**Measured on the 58 (development ladder, pilot on 176 kept-matching templates):** wide grid alone does
nothing; robust L1 misfit 27.6 -> 10.5 km; SNR^2 weights -> 9.1; horizontal envelope + LOO station terms +
nested refinement -> 8.5; IRLS -> 7.8; soft Slab2 depth prior -> 7.3 (rail 0.04); **station-pair envelope
cross-correlation DIFFERENTIAL times (Chestler's own step) -> 4.0 km epicentral, p90 19 km, 3-D 6.4 km,
rail 0.09**.  DISCARDED because measured worse: AIC onset picks (16.9 km), P picks / S-P (19.6 / 12.6 km),
envelope back-projection (15.6 km), Vp/Vs 1.72.  Also measured: Bostock's own published family depths sit a
median **9.7 km BELOW** this project's Slab2 field, so the `|dz| <= 8 km of Slab2` screen clause is wrong.
**Arm B (739 Chestler-chain families) had NEVER been located** -- `n_locate` crashed with KeyError 'B005'
(logs/o_armb_locate.log): it used the 22 beam SITES as the station table while Arm B correlates 41 member
stations.  Fixed here with a `wa_mem` bed adapter.
- 2026-09-23 ~02:30Z **AGENT U - LOCATOR REBUILT AND CALIBRATED ON BOSTOCK'S PUBLISHED LOCATIONS** (queue u_locate_queue.sh pid 1428881, markers U_LOCATE_BASELINE_DONE done -> U_LOCATE_DONE|FAILED, report results/WA_LOCATION_REBUILD.md).
  **BASELINE, never measured until now:** the current k_locate, run on the 58 blind-rediscovered Bostock families and compared with their PUBLISHED coordinates (e5r_families_2012.csv = Bostock 2015's own triples): locates 44/58, median epicentral error **27.6 km**, p90 87.9, median |depth| 15.9 (bias -11.7), median 3-D 35.2, **75% railing**, only 11/58 within 10 km. Chance (permuted labels) = 39.8 km. So our locator was barely better than chance and that is why every location-based WA test was ambiguous.
  **DIAGNOSIS: travel times are fine, the PICK was wrong.** At stack envelope peaks >=8 MAD the observed S matches the Crosson 1-D prediction to 0.15 s (90% within 1 s; station medians -0.8..+1.9 s). But 22% of admitted picks are >4 s out and unweighted least squares let them rail the solution.
  **LADDER (cumulative, 176 kept-matching templates):** wide grid 0 -> robust L1 10.5 km -> SNR^2 weights 9.1 -> horizontal envelope + LOFO station terms + nested refinement 8.5 -> IRLS 7.8 -> soft Slab2 depth prior 7.3 (rail 0.04) -> **station-pair envelope cross-correlation DIFFERENTIAL times (Chestler's own step) 4.0 km epicentral, p90 19 km, 3-D 6.4 km, rail 0.09, 35/58 within 10 km.** ~7x better than baseline.
  **DISCARDED because measured worse:** AIC onset picks (16.9 km - the LFE stack onset is emergent), P picks and S-P (19.6/12.6), envelope back-projection (15.6), Vp/Vs 1.72/1.85, hard Slab2 depth.
  **ADOPTED:** L9_xcdiff (min median 3-D error at >= baseline coverage) + L12_comb as max-coverage companion (56/58 located, 39/58 within 10 km); per-family bootstrap radius + shuffled-station control carried as COLUMNS, nothing filtered.
  **BY-CATCH 1:** Bostock's own published family depths sit a median **9.7 km BELOW** this project's Slab2 field -> the old "|dz| <= 8 km of Slab2" screen clause was wrong on its face, independently of the rail pathology.
  **BY-CATCH 2:** the 739 Arm-B (Chestler-chain) families had NEVER been located - n_locate crashed with KeyError 'B005' because it used the 22 beam sites as the station table while Arm B correlates 41 member stations. Fixed with a wa_mem bed adapter, so the Chestler-coordinate check will see Arm B for the first time.

## 2026-09-23 -- agent S: Washington re-scored under RANK-DON'T-DISCARD (no screen, no candidate count)
User rule, 2026-09-23: criteria calibrated on the 58 known Bostock families can only catch cuts that
are too strict for THAT population and protect no undiscovered family, so nothing may be filtered.
Every grown template of agent O's array-based 2010 chain, agent N's 2008 chain and Chestler's own
Arm B now keeps its row in `results/sf_final/table_<bed>_<chain>.csv`, with every criterion as a
column and a byte-identical table for the whole-chain site-shift null; the only list-reducing steps
are the published dedupe and the location, both recorded as flags with the rows retained.
Evidence = real-vs-null separation curves and AUC per column (`sep_<bed>_*`), not thresholds.
The largest separation IN FAVOUR OF THE REAL CHAIN anywhere is on **2010 Arm B, Chestler's own chain (agent O)**, column `rms`: AUC 0.815 (0.5 = indistinguishable), i.e. |AUC-0.5| = 0.315, best lift 5.70x (favours real).
New work behind it: Arm B was located for the first time (agent O's locate stage had crashed on a
member-station KeyError) and given the site-shift null it never had; every template of every chain
was re-fired on the sealed 2011 / 2007 windows, not just a screened subset; the confirmation rule
was calibrated on Bostock's own families (north half-split) and shown to be a weak discriminator
(16/18 true vs 60/84 noise), so it is reported as columns, not a gate;
the depth window was measured and rejected as a cut (100 % of its Washington rejections are the
bottom grid rail; its transferable form keeps only 66 % of located true north families).
Full ledger `results/WA_FINAL_ASSESSMENT.md`.
- 2026-09-23 ~04:50Z **S_FINAL_DONE — AND ITS VERDICT IS WRONG; CORRECTED HERE BY THE ORCHESTRATOR.** results/WA_FINAL_ASSESSMENT.md §6 concludes "real and null are not distinguishable by any statistic this project can compute". Its own table contradicts that: the `dedupe_n_partners` column was given direction "lo" in sf_separation.py, so a strong REAL-favouring signal was printed as a reverse separation (AUC 0.242).
  **Recomputed in the correct direction from results/sf_final/table_oly2010_{real,null}.csv (orchestrator, verified):** mean partners real 15.71 vs null 0.47 (33x); median 2 vs 0; frac>0 0.703 vs 0.301. **AUC (hi favours real) = 0.758 — the strongest discriminator on the 2010 Olympic bed, beating rms 0.642 and comparable to Arm B rms 0.815.** Lift by threshold: >=1 partner 2.3x, >=2 5.2x, >=5 45.3x, >=10 122x, **>=20 partners: 280 real templates and ZERO null templates; >=50: 161 vs 0; >=100: 60 vs 0 (complete separation).**
  **Detection-richness control (the one I flagged as owed) PASSES:** the real bed is 10x detection-richer overall (mean n_det 172 vs 17), but at MATCHED detection count the signal survives — median partners real vs null by n_det bin: 0-20 1 vs 0; 20-50 2 vs 0; 50-100 5 vs 0; 100-300 19 vs 3. So the sharing is not an artefact of detection count.
  **Physical reading:** many templates sharing detections = the same source independently re-discovered from different seed windows, which is what an LFE family looks like to this chain. The Olympic 2010 array bed is NOT a null result.
  **Two corrections owed before the assessment can stand:** (1) fix the direction flag for dedupe_n_partners (and audit every other column's direction) and rewrite §6; (2) every location-derived column in the tables (lat/lon/dep, edge, dep_rail, near_chestler_km) was computed with the OLD locator, which agent U has now measured at 27.6 km median error against Bostock's published coordinates vs 4.0 km for the rebuilt one — so the Chestler comparison in §3 is not yet meaningful. Recompute after U_LOCATE_DONE.

### 2026-09-23 16:02 UTC — u_locate_queue crash at `north_eval2`, fixed and resumed (agent W)

**True cause.** Not a short frame, a missing family or an empty variant: `uu_dd.py` wrote its
differential-time observation table as `results/u_locate/north/loc_ZZdd_obs.csv`
(`a.out.replace(".csv","_obs.csv")`), and `uu_eval.py` globs `loc_*.csv` in that directory and
treats every match as a location table. The obs table has columns `i,j,sta,ddt,cc` — no `grp` —
so `R.merge(L, on="grp")` raised `KeyError: 'grp'` and killed the queue 4 s after the
double-difference pass succeeded. The locator ladder was never at fault; the crash was purely a
filename collision in the evaluator's glob.

**Fix.** (1) `uu_dd.py` now writes `dd_ZZdd_obs.csv`, which cannot be globbed as a variant (the
stale file was renamed). (2) `uu_eval.py` is hardened independently of that: `why_not_a_location_
table()` checks `grp,lat,lon,dep` and emptiness before scoring, every variant is wrapped in
try/except, and a variant that is unreadable, empty, not a location table or that locates nothing
is now **reported as a ladder row with a `status` column** instead of aborting the run —
rank-don't-discard applies to variants too. `score()` also tolerates a missing `located` column
and NaN coordinates. Re-running the evaluator reproduced all 36 previously measured rows to the
last digit, so the hardening changed no number.

**Double-difference verdict (new, measured).** `ZZdd` scores 3-D median 7.40 km vs `L9_xcdiff`
6.44 km, epicentral 4.73 vs 4.00 km, 35/58 within 10 km either way → **NOT ADOPTED**, reported and
not used. Only 1,123 differential times over 179 linked pairs survive cc>=0.4 / >=4 shared
stations out of 3,363 candidate pairs, so most families get no differential constraint at all
(median |move| 0.00 km).

**Priority question wired in.** `uu_wa.py` now joins `results/sf_final/table_oly2010_{real,null}
.csv` onto every located Washington row (`dedupe_n_partners`, `sf_top` and the old location as
extra COLUMNS — no family is filtered) and adds an explicit sub-population block: the 280 real
templates with >=20 dedupe partners, the null's equal-size best-ranked head as the matched
comparator, the empty `null[>=20]` row kept so the table is identical, and both whole populations
for scale — with bootstrap radii, distances to Chestler's 69 published coordinates and counts
within 5/10/20 km in both directions. `uu_report.py` renders it as §9b.

**Resume.** Sentinel `results/u_locate/north/.ladder_done` added so a relaunch never re-walks the
ladder. The relaunch skipped stacks, baseline, picks, residuals, all 34 ladder variants, the
2.1 h bootstrap and the dd pass in 30 s and went straight to the Washington stages.
- 2026-09-23 ~16:15Z **WA locator queue FIXED and RESUMED (agent W), pid 1505597, ~8-12 h.** Cause was a FILENAME COLLISION, not a data problem: uu_dd.py wrote its differential-time observations as loc_ZZdd_obs.csv, and uu_eval.py globs loc_*.csv and treats every match as a location table; that file has columns i,j,sta,ddt,cc and no grp, so the merge raised KeyError 'grp' 4 s AFTER the dd pass had succeeded. Fix: dd observations renamed dd_*_obs.csv, and uu_eval hardened (why_not_a_location_table() checks grp/lat/lon/dep and emptiness; every variant wrapped; unreadable/empty/zero-solution variants become a ladder row with a status column instead of killing the queue). Re-running the evaluator reproduced all 36 previously measured rows to the last digit - the hardening changed no number.
  **New measured verdict: double-difference NOT ADOPTED.** ZZdd 3-D median 7.40 km vs L9_xcdiff 6.44; epi 4.73 vs 4.00; 35/58 within 10 km either way; only 1,123 differential times over 179 of 3,363 candidate pairs survive, median |move| 0.00 km. Reported, not used.
  **Resume worked as designed:** sentinel .ladder_done; the relaunch skipped stacks/baseline/picks/residuals/all 34 ladder variants/the 2.1 h bootstrap/the dd pass in 30 s. Remaining: WA real -> null -> Arm-B (wa_mem adapter verified, 123 member-station stacks, no KeyError B005), then the Chestler comparison for both L9_xcdiff and L12_comb.
  **Priority question wired in as sec.9b:** uu_wa.py joins sf_final/table_oly2010_{real,null}.csv onto every WA row (dedupe_n_partners, sf_top, old location as COLUMNS, nothing filtered) and reports the 280 real >=20-partner templates against the null's equal-size best-ranked head, with bootstrap radii, distance to Chestler's 69 coordinates and counts within 5/10/20 km in both directions. Report refreshes after each WA population, so it is readable before the queue ends.
- 2026-09-24T16:13:51Z PLAN written: notes/MARGIN_CATALOG_PLAN_2026-09-24.md — phases 1-4, decision points D1-D4, and an explicit rule against starting a new region while a running one has an unresolved decision point (the Oregon mistake).
- 2026-09-24T16:15:22Z FRAMING CORRECTION (user): a published POSITIVE result is itself testable without the authors' catalogue — if they found families in a region and we find none there, WE failed. Retracts 'Royer & Bostock cannot be tested'. The two Puget zeros are bug reports, not evidence about Washington. Scoring rule for catalogue-less regions: existence of located families above the null inside the published region and depth band; their numbers are reference lines. Written into notes/MARGIN_CATALOG_PLAN_2026-09-24.md.
- 2026-09-24 **ROYER & BOSTOCK 2014 FULL TEXT OBTAINED (user supplied notes/1-s2.0-S0012821X13004780-main.pdf + supplement notes/royer2014_supp.pdf). THEY LOCATED 176 WASHINGTON FAMILIES. Every previous "we cannot test Royer & Bostock" statement is RETRACTED, and the y_wa_finish queue was stopped because it was running the wrong method.**
  **What they report:** 4 subarrays NW / SW / NVI / SVI processed independently. **NW = 122 final templates, SW = 54** (initial clusters 224 NW, 70 SW; initial autocorrelation detections 4,915 NW, 3,306 SW). Each template has 100s-1000s of independent detections. Located with **Hyp2000**, relocated with **hypoDD** (Table S1: IDAT=3 catalog+cross-correlation, 40 iterations, DAMP 40-50). **Depths: NW 31-46 km, SW 32-49 km** (SVI 29-40, NVI 33-39). Locations appear in Fig 4, Fig S1 (hypoDD) and Fig S2 (depth profiles) - NO coordinate table, so their positions are only digitisable from figures.
  **Their exact recipe, which is NOT what we ran:**
   - NW autocorrelation: **7 three-component stations B001, SQM, BS11, W020, W040, GNW, PL11**; episodes 2007-01-17..30, 2008-05-01..23, 2010-08-07..24, 2011-08-04..23.
   - SW autocorrelation: **7 three-component stations N050, S030, W070, N060, PL11, S040, N070**; episodes 2007-01-14..02-01, 2008-05-01..25.
   - **15-second windows lagged by 0.5 s** (7,170 windows/hour), 1-8 Hz, 40 sps, 24-h segments, hour-by-hour, network autocorrelation stacked over all stations and components, **8xMAD**, then cull overlapping detections keeping high SNR.
   - Grouping: **hierarchical cluster analysis on a similarity matrix of "super traces"** (Rowe et al. 2002), not connected components.
   - Then network cross-correlation using **subsets of 7-10 high-SNR 3-C stations**, 8xMAD over a **24-hour** period, restack, iterate to convergence.
  **Our deviations, all three material:** we used **51 stations** (they used 7 - the Ducellier dilution lesson again), **6-second windows** (they used 15 - we took 6 s from Bostock 2015/Plourde, the wrong paper for this bed), and **connected-component grouping** (they used hierarchical clustering). Band, lag, rate and threshold matched.
  **Consequence: the two Puget zeros (WA_2008_RESULTS.md, WA_2008Q_RESULTS.md) and the stopped y_wa_finish run tested a method nobody published. Washington has never been attempted with the method that produced the 176 families.**
- 2026-09-24 **ALL SIX PAPERS READ AT LAST (agent Z). notes/PAPER_METHODS_TABLE.md + notes/CATALOGUE_SUMMARY.md.** Every paper was already on disk in notes/ while the orchestrator asserted they were unavailable and reconstructed parameters from secondary sources.
  **Family counts and whether the authors LOCATED them:** Royer & Bostock 2014 122 NW / 54 SW / 47 NVI (initial 224/70/54), ALL located (Hyp2000 then hypoDD, retained 77%/72%/91%/92%), depths 31-46 NW, 32-49 SW, 33-39 NVI. Bostock 2012: 140 templates from 5,775 detections, ALL located (hypoinverse, ~16 P + 22 S per template), 27-37 km, deliberately NO dedupe. Chestler & Creager 2017: 57 -> 43 families, 32,387 LFEs, ALL located (3-D Preston grid search VP/VS=sqrt3 + Waldhauser-Ellsworth double-difference). Thomas & Bostock 2015: 18 families, >15,000 detections, only 13 LOCATED and just 4 with free depths (31-34 km), HYPOINVERSE. Bostock 2015: 130 templates, 269,586 detections, did NOT locate (inherited). Plourde 2015: 59 templates -> 34 families, 9,900 LFEs, all located 28-47 km. Ducellier 2022: used Plourde's unreduced 66 -> 60 kept, did NOT locate (adopted Plourde's).
  **TOP 5 MISMATCHES between our code and the papers:** (1) Washington window: ours WIN_S=6.0 (k_auto.py:42) at 20 sps (k_common.py:28) vs Royer's 15 s at 40 sps - and our own docstring cites the 6.0 to BOSTOCK 2012, a different region. (2) Washington stations: our 52 probe-admitted (n_common.py:171-208) vs their 7 NAMED; SQM, GNW, W020 are not even downloaded for NW, N070 missing for SW. (3) Grouping: our connected components (k_group.py:51-65,96-99) vs their hierarchical clustering on a super-trace similarity matrix - components CHAIN, clustering does not. (4) Chestler front end entirely absent: no 18xMAD beam peaks, no boxcar cross-array stack, no 2nd-root stack, no 45xMAD, no Rcom>0.5 dedupe. (5) Ducellier MAD block: our BLOCK_H=12 (t7_duc_mf.py:26) vs his 1 HOUR (paper, reviewer reply, and his find_all_LFEs_parallel.py:220,227,321) - a one-line fix that plausibly explains part of the 0.836 recall shortfall.
  **Two corrections to our own claims:** our 4-s declustering is Bostock 2015's MAGNITUDE-catalogue cull, not any detector's published rule; and the Olympic "61 of 69" is coverage against derived_coordinate_rows_not_a_catalog.csv, NOT 61 of Chestler's 43 families - the comparison must be restated.

## 2026-09-25 — agent AB: Royer & Bostock (2014) Washington — code written, bug audit done, NOTHING RUN

Brief: reproduce their **122 NW + 54 SW** located Washington LFE families with their own method, on
their own 14 named stations, after auditing the Washington chain for the zero-fill / no-validity-mask
bug class that invalidated the southern Shelly run.  Mid-task instruction change: **write the code and
do not run it.**

**Disclosure, and a lesson.**  I had launched the stage-0 quantifier `ab_gapscan.py` at 03:33Z before
that instruction arrived and issued `kill -TERM` against the pid the shell handed back.  That pid was
the backgrounding subshell, not python; the python process survived as an orphan and **ran to
completion at 03:41Z**.  `ps -p <pid>` returned nothing, so I reported it stopped when it was not.
It is read-only -- it opens cached `.npz` files and writes two CSVs, no download, no GPU, no chain
stage -- and its output (`results/ab_royer/gapscan_{cache20,hours}.csv`) is used in the audit because
it CORRECTS my header-based reading.  Lesson for every future kill: `nohup cmd &` inside the harness
gives you the subshell's pid; kill the process group (`kill -TERM -<pgid>`) or match on the script and
verify with `pgrep -af <script>`, not `ps -p`.  Nothing else was executed: no download, no GPU, no
detection stage, no 40 Hz cache.

### Stage 0 — bug audit verdict (`results/WA_BUG_AUDIT.md`, read-only)
**The southern bug is NOT in the Washington chain.** The south *floors* the cc denominator
(`sd = sqrt(max(M*S2-S1*S1, 1e-30))`, `lfe_south/scripts/t11_mf.py:27`) so a constant run divides an
FFT-round-off numerator by 1e-15 and explodes; Washington sets a rejected window's denominator to
`+inf` (`k_auto.py:76`, `k_grow.py:167-169`), which forces cc = 0.  Every Washington cc is
Cauchy-Schwarz bounded by 1.  **Three weaker forms of the same class ARE present:**
- **D1** the guard is RELATIVE (`good = ew > 1e-4 * emed`) and degenerates to `ew > 0` whenever a
  component is more than half dead (`emed == 0` exactly).  `k_auto.py:69-70`, `k_grow.py:166-169`,
  `o_auto.py:66-76`.  Gap-edge windows carry the band-pass impulse response and pair at high cc.
- **D2** masked pairs enter the median/MAD pool as exact zeros (`k_auto.py:92-100`,
  `k_grow.py:202-204`), deflating `med + 8*MAD` → false-alarm inflation.
- **D3** the network divisor is a constant channel count (`k_auto.py:90`) or a per-station scalar
  (`k_grow.py:199-200`), not the per-sample count of contributing channels → dilution where a station
  is in a gap, and a non-stationary statistic inside the block whose MAD is taken.
- `o_beam.py:126` has the same relative-guard defect (`x/x.std()` after only an `std <= 0` test), so a
  mostly-gap member enters the delay-and-sum beam at equal weight.  `uu_*.py` is clean for the adopted
  locator (it reads stacks, and the one floor-in-a-log defect, `uu_picks.py:68-70`, is in the AIC onset
  variant that was measured and discarded).
- **D4** `k_common.load_day_native` zero-fills with no mask and no guard, but only the `--chain native`
  sweeps use it; no number on record comes from it.

**Two further defects, both found on the way and both worth as much as the bug.**  (a) **The exact
zero runs DO survive the band-pass and the decimation** -- `cache20hz/UW.GNW/2010/224.npz` is 41.93 %
non-zero on all three components with a **50,174 s** all-zero run, exactly matching its
`coverage_frac` of 0.418.  So the southern pathology's raw ingredient is in this cache too; only the
guard differs.  (My first inference from the manifest, that filtering smeared the zeros away, was
wrong and is retracted.)  (b) **The manifest column built to detect this is simply wrong**:
`cache_T*.parquet` records `frac_nonzero = 1.000000` for that same file (`o_cache_month.py:53,92`), and
cache-wide it claims only 10 of 23,744 station-days below 1.0 and none below 0.5.  It cannot be used as
a gap indicator anywhere.  `coverage_frac` and `e8_Z/N/E` in the same manifest ARE right and both
flagged the day; nothing downstream reads either.

**Measured exposure of D1, two resolutions, and the finer one is the one that matters** because
`k_auto` takes its MAD over the HOUR.  *Per station-day* (`e8_Z/N/E == 0` in `cache_T*.parquet`, which
IS `emed`): 81 station-days cache-wide; **24 of 3,867 (0.62 %)** in the 2008 Puget admitted roster
(YW.FACB 19, PB.B018, XU.N050, XU.BS11); **3 of 6,817 (0.04 %)** in the 2010 Olympic roster (UW.GNW,
UW.LRIV, XG.PA02); 1 of 756 Royer episode station-days.  *Per station-HOUR on Royer's own bed*
(`gapscan_hours.csv`): **58 of 18,138 (0.32 %)**, i.e. **5x** the day-level exposure, and on **four
stations the day view calls clean** -- PL11 (11 h in NW08, 11 h in SW08), S040 (11 h), N060 (9 h), all
on 2008-05-09, plus GNW (14 h on 2010-08-12, 1 h each in NW07/NW08).  The PL11/S040/N060 hours have
`frac_zero = 0` and no zero-energy 15 s window at all, yet **one or two individual components are dead
for more than half the hour**, so that component's `emed` is exactly 0, its guard becomes `ew > 0`, and
the 19-55 % of its windows the guard should have rejected (`frac_win_tiny`) are all admitted.  A single
dead component inside a live three-component station is invisible to every all-component zero test this
project has applied.  2008-05-09 is exactly the day with `coverage_frac` 0.75 (PL11, S040) / 0.84
(N060): the gap is real, per-component, and only a per-component test sees it.  **Per-hour exposure on
the 2008 Puget and 2010 Olympic rosters was NOT measured** -- extending `ab_gapscan.py` to them is a
one-argument change and is the first thing to run if anyone wants sec.3's verdict tightened.

**What this invalidates: nothing on record.**  The 2010 Olympic 61-of-69 coordinate match and the 280
dedupe-rich templates sit on 0.04 % exposure and their site-shift null received the identical
treatment (null 3-of-69 and a null richness maximum of 16).  The 2008 zeros are **not** explained by
this bug — D1/D2 push counts UP and that result was zero families; its named cause stands (6 s window,
20 Hz, 52 substituted stations: `PAPER_METHODS_TABLE.md` M1/M2).  The locator calibration is untouched.
What is owed is not carrying the defect forward, and the `ab_` chain does not.

### The data gap is much smaller than the brief assumed
GNW, W020 and N070 are all on disk — the earlier "missing" claim was about the probe-admitted roster,
not the archive.  Of **833** station-days Royer needs (7 NW stations x 75 episode days + 7 SW x 44,
PL11 shared), **638 raw are on disk** and only **195 must be fetched**: **UW.SQM** entirely (75 days;
absent from raw and from the 20 Hz cache) and the **2011-08-04..23 NW episode** for the other six
(120 days — cached at 20 Hz but the raw was freed, and a 40 sps cache can only be built from raw).
At the measured 19.7 MB/station-day that is **~3.9 GB** to download; the 40 sps cache is **~37 GB**.
69 TB free.  Native rates: XU 50 Hz, PB.B001 100 Hz, UW.GNW 40/50 Hz — all need the rational chain.
Open risk recorded, not assumed away: the paper calls SQM three-component; if FDSN returns only a
vertical, NW runs with 6 of 7 and that is a reported column.

### The parameter mismatches that matter (full table: `results/AB_ROYER_PLAN.md` sec.1)
40 sps not 20; **15 s** windows not 6 s (the 6 s came from Bostock 2012, a different paper and region);
their **7 named stations** per subarray, not 51-52 probe-admitted ones; all **four** NW episodes
(2007, 2008, 2010, 2011) and **both** SW episodes, not just 2008-05 and 2010-08; **hierarchical
clustering on a super-trace similarity matrix**, not union-find connected components (which chain, and
which is why we produced giant components instead of 224 clusters); restack normalised to **unit
maximum amplitude across three components**, not rms; **all** newly registered detections in the stack,
not the best 200; **iterate to convergence**, not a fixed 3; 8 x MAD over the **24-hour period** at the
growth stage.  Everything the paper is SILENT on — linkage, cut, minimum cluster size, iteration count,
the numeric cull rule, the template length — is a **swept column**, and the cluster cut is chosen by
the mean silhouette, a data-internal criterion that never looks at 224/70 (the cut nearest the
published count is reported beside it, never substituted).

### Code written, syntax-checked, unrun
`scripts/ab_common.py` (paper parameters, the masked 40 sps cache API, the whole-chain site-shift null
remap), `ab_gapscan.py` (stage 0 quantification -- the one script that did run, read-only), `ab_fetch.py` (195 station-days, <= 8 FDSN
connections), `ab_cache40.py` (40 sps 1-8 Hz **with a per-sample validity mask** built from the trace
spans before the fill, widened by the filter transient), `ab_auto.py` (sec.3.1, both readings of the
lag axis), `ab_cluster.py` (sec.3.2, super traces -> lagged all-pairs cc -> similarity matrix ->
hierarchical clustering sweep -> shift-and-stack), `ab_grow.py` (sec.3.3), `ab_locate.py` (L9_xcdiff +
L12_comb via `uu_*` with a `royer` bed installed by monkeypatching `uu_common.bed`, input decimated to
the locator's calibrated 20 Hz so a calibrated instrument is not fed an uncalibrated rate),
`ab_null.py`, `ab_score.py`, `ab_royer_queue.sh` (flock, TERM/INT/HUP/EXIT traps writing
`results/AB_ROYER_FAILED`, per-stage output guards, `env AB_STAGE=n python` per stage, 120 s heartbeat,
markers `AB_ROYER_AUDIT_DONE` / `AB_ROYER_AUTOCORR_DONE` / `AB_ROYER_DONE` / `AB_ROYER_FAILED`, pid in
`AB_ROYER_PID`).  `bash -n` and `python -m py_compile` pass on all of it; imports and array shapes
checked on synthetic noise only (window_stats shapes and mask semantics, null_remap determinism and
episode containment, `autocorr_hour` on 4 synthetic stations -> 12 channels, 371 rows, 0 detections at
8 x MAD on white noise, `unit_max3`, the cluster sweep and silhouette on a synthetic block matrix).

### Cost, and the decision left open
To the end of stage 6 for both subarrays: **~3-5 days** (stage 5 growth dominates, 20-40 h for NW);
the three null realisations add **4-8 days** and are queued one at a time with scoring re-run after
each.  Peak GPU 8-20 GB, so it coexists with the two southern queues.  **The one decision that wants
review before launch:** the growth station pool.  The default is the subarray's own 7 (the literal
lower bound of "subsets of 7 to 10 high SNR, three-component stations"); widening to the 25-station
regional pool of sec.2.1 is the change most likely to lift the final count toward 122 and costs
~4,760 extra station-days (~95 GB raw, ~200 GB cache).  It is a `--pool` argument, no code change.

Launch, when approved:
`cd /home/jovyan/tremorferometry/lfe_array_discovery && nohup setsid -f bash scripts/ab_royer_queue.sh >> logs/ab_royer_queue.log 2>&1 < /dev/null`
- 2026-09-25 ~04:15Z **AGENT AB: Washington code WRITTEN, NOT RUN (user instruction). results/AB_ROYER_PLAN.md (parameters + downloads + launch command at sec 8), results/WA_BUG_AUDIT.md.**
  **Bug audit verdict: the southern unbounded-cc bug is NOT in the Washington chain** - k_auto.py:76 and k_grow.py:167-169 set a rejected window's denominator to +inf, not to a 1e-15 floor, so every cc stays Cauchy-Schwarz bounded. **It invalidates nothing on record**, and it could not explain the 2008 zeros anyway (the bug inflates counts; that result was zero). The named cause of the Puget zeros stands: 6-s windows, 20 sps, wrong stations, wrong grouping, missing second stage.
  **Three weaker forms of the class ARE present and should be fixed before the run:** D1 the dead-channel guard is relative (ew > 1e-4*emed) and degenerates to ew > 0 when a component is >50% dead (also o_beam.py:126); D2 masked pairs still enter the median/MAD pool, deflating the 8xMAD bar; D3 the network divisor is a constant channel count, not the per-sample valid count. uu_* is clean for the adopted locator.
  **Two harder findings.** (i) Exact zero runs survive filtering: GNW 2010-08-12 is 41.9% non-zero with a 50,174-s constant run, and the manifest column built to catch exactly that (frac_nonzero) reads 1.000000 on that file - unusable as a gap indicator. (ii) Per-HOUR exposure on Royer's bed is 58/18,138 hours (0.32%), 5x the day-level view, and on four stations the day test calls clean (PL11, S040, N060, GNW) because one dead component inside a live 3-C station passes every all-component zero test. Day-level: 24/3,867 (2008 Puget), 3/6,817 (2010 Olympic); per-hour unmeasured there.
  **Downloads much smaller than briefed:** GNW, W020 and N070 are already on disk. Of 833 needed station-days, 638 raw are present; **195 to fetch ~3.9 GB** (UW.SQM entirely, 75 d; plus the 2011-08-04..23 NW episode for six stations, 120 d - cached at 20 Hz but raw was freed and 40 sps needs raw). The 40-sps cache will be ~37 GB.
  **Parameter mismatches to fix in the run:** 40 sps not 20; 15-s windows not 6; their 7 named stations not 51; all four NW episodes and both SW; hierarchical clustering not connected components; unit-max-amplitude restack not rms; all detections not best-200; iterate to convergence not a fixed 3.
  **Disclosure:** ab_gapscan.py (read-only) had been launched before the no-run instruction arrived and survived kill -TERM as an orphan (the pid killed was the subshell's), finishing 03:41Z; its two CSVs are the source of the per-hour numbers above and the queue's stage-0 guard will skip it. Nothing else executed.
