# FLEET BROADBAND/LAND dv/v PIPELINE — authoritative (restructured 2026-07-18 per Merlin)

## §0 SCOPE
Governs the **178-station BROADBAND fleet** (BH/HH land stations, `data/broadband_fleet_order.csv`, region-
stratified). Boreholes are governed by `notes/FINAL_PIPELINE.md` (different selection/band/provider rules — do
NOT apply it here). Sections marked **[FROZEN]** change only via a merlin consult + a §7 changelog entry.
Goal: reliable LFE-coda dv/v from each station → more receivers → higher-resolution 4-D δβ/β on the interface.

## §1 RESUME PROTOCOL  ← READ THIS FIRST AFTER A CRASH
**STEP ZERO — before launching ANYTHING, list what's still running:**
`ps -eo pid,etime,cmd | grep -E '[f]leet_broadband_all|[f]leet_station|[f]leet_download_ahead|[f]leet_trace_janitor|[d]ensify_gnw|[d]iscover_gpu'`

- **If `fleet_broadband_all` is ALIVE → TOUCH NOTHING.** It self-resumes each station. Just monitor:
  `tail -f logs/fleet_broadband.log`. Do NOT launch a second batch (two instances race the same station + GPU).
- **If DEAD (pod restart):** environment checklist first —
  1. `/home/jovyan/envs/tremorferometry/bin/python` exists? (pod restart can wipe /opt/conda — see mem: pod-env-ephemeral)
  2. cupy CUDA-12 libs present? `ls /opt/conda/lib/python3.13/site-packages/nvidia/*/lib` (the recurring libcufft fix)
  3. GPU visible? `nvidia-smi`
  Then relaunch nohup'd, in this order (each is idempotent; batch self-resumes via `STATION DONE` sentinels):
  ```
  cd /home/jovyan/tremorferometry
  nohup bash scripts/fleet_broadband_all.sh   > logs/fleet_broadband_top.out 2>&1 &
  nohup bash scripts/fleet_download_ahead.sh  > /dev/null 2>&1 &
  setsid bash scripts/fleet_trace_janitor.sh  >/dev/null 2>&1 < /dev/null &   # (survives session end)
  ```
- **Mid-station truncated checkpoint:** if a resumed station fails INSTANTLY at score/cluster with a
  parquet/CSV read error, the crash truncated that station's checkpoint. Delete ONLY that station's stale
  intermediate (`data/<s>_pnsn_candidates.parquet` / `_cand_baseline.parquet` / `<s>_disc.summary.csv`) and let
  it re-run. densify/stacks/dvv have no guard and simply re-run from scratch — expected, NOT a bug.
- **RED LINES:** never a 2nd batch instance · never clean `logs/` (it IS the resume state) · never delete
  `data/long_window_daily_*` / `daily_dvv_*` / `*_causality_cert.csv` (the science products).

## §2 LIVE STATE  (mutable — any number older than the as-of date is a HINT; recompute from disk)
- **COMPLETE 2026-07-21:** 178/178 processed · **INCLUDE 124 / FLAG 55** (incl tier-2). Network total **12,703
  certified families / 187 stations** (30 borehole + 124 broadband + 3 anchor). Map `figures/borehole_dvv_map.html`
  (164 stns). END-OF-RUN METHODOLOGY + deep-resolution study: **`notes/METHODOLOGY_END_OF_RUN_2026-07-21.md`**
  (deep dv/v = large-scale ~200 km regional INDEX/bound, not a fine 4-D map; 0.77 checkerboard retired; captured
  fraction ~5%; scripts in `fault_tomography/inversion/`). Top INCLUDEs: OTR 155, YELM 155, WISH 153, DOORS 146.
  Re-run queue (6, `data/rerun_queue.txt`) incl JCC (still 0 candidates via NCEDC — needs a look).
- **KSXB re-run DONE (2026-07-18):** HHZ → 78 certified / 300 (26%) → INCLUDE (was 2, invalid accelerometer run).
  Recovered via bug-6 fix. Re-run queue now 6 (FISH, JCC, GOBB, SYMB, KMR, KHMB — KMR/KHMB likely recover too).
- **Recompute from disk:** progress = `grep -c 'RESULT:' logs/fleet_broadband.log` · per-station verdicts =
  `GATE:` lines in `logs/fleet_<s>.log` · current station = the running `fleet_station.sh` proc · scoreboard =
  scan `logs/fleet_broadband.log` for `RESULT:` (INCLUDE/FLAG + certified/densified/survival).
- **Queues:** `data/rerun_queue.txt` (bug-flag re-runs) · `data/tier2_survival_borderline.txt` (survival-borderline).
- **Pending end-of-run decisions:** survival-bar policy · tier-2 inversion inclusion (needs robustness check) ·
  CN/NRCan sparse recovery · short-record stations · final coverage map → 4-D inversion inputs.
- **UPDATE this section on:** any bug fix / policy change / queue change / every ~10 stations. Do NOT edit
  [FROZEN] sections without a §7 changelog entry.

## §3 RECIPE [FROZEN] — per station, all 40 Hz. Driver: `scripts/fleet_station.sh <NET> <STA> [START] [CHAN=auto]`
1. **BAND** `pick_band.py <NET> <STA>` → ONE band. Prefers broadband BH/HH (+100) over short-period EH/SH;
   requires **high-gain SEISMOMETER (2nd char 'H'), NEVER HN=accelerometer / L=low-gain**; aggregates span
   PER BAND across epochs (picks the longest era). Provider: **BK/NC → NCEDC**, else IRIS/EarthScope.
2. **DOWNLOAD** `download_broadband.py --net N --sta S --chan <band> --start .. --end 2026-08-01` — ALWAYS runs
   (resumable; skips existing days + `.nd` nodata sentinels). Never skip a partial fragment.
3. **DETECT@40** `discover_nllb_pnsn_driven.py --station S --bbox <lat±0.9 lon±1.35> --pnsn <catalog> --fs 40
   --candidates-only`.  Guard: **<100 candidates → FLAG** (no tremor/data).
4. **SCORE@40** copy `picker_broadband_pgc.joblib`→`tremor_picker_<s>.joblib`; `score_candidates.py --net N
   --sta S --target-fs 40` (**--target-fs 40 MANDATORY** — picker trained at 40 Hz).
5. **SELECT** `adaptive_cand_threshold.py S 30000` (rank-based top-30k; P floor **0.05** = garbage bound only).
6. **CLUSTER@40** `discover_gpu.py --station S --fs 40 --cc-threshold 0.80 --min-family-members 3 **--min-years 1**`
   (min-years 1 MANDATORY; default 3 = 4×-yield trap).  Guard: **<50 families → FLAG** (skip densify).
7. **COVERAGE-SELECT** `select_families_coverage.py <disc> - 300 snr` — rank by **SNR** (validated AUC 0.68), cap
   300. **NEVER rank by picker p_lfe (anti-predictive AUC 0.40).** (appB `score_family_stacks_picker` is DROPPED.)
8. **DENSIFY (fwd only)** `densify_gnw_gpu.py --templates-npz <disc>.npz --summary-csv <disc>_sel300.summary.csv
   --network N --station S --fs 40 --min-snr 0 --top-n 0 --despike-mad 8 --out-prefix mf_<s>p90f40_`.
9. **STACKS** `build_long_window_3comp.py --mf-csv-glob "data/mf_<s>p90f40_*.csv" --fs 40 --cc-min 0.80
   --min-det 20 --despike-mad 8 --out-prefix data/long_window_daily_<TAG>`  (TAG = STA uppercase).
10. **dv/v** `dvv_roll30cal.py --station S --npz ..._Z.npz --window 2 4 --origin-anchor` → `daily_dvv_<TAG>_Z_2to4.csv`
    (30-cal-day trailing rolling stack, SVD-Wiener, coda stretch 2–4 s vs per-era all-time ref).
11. **CAUSALITY** `finalize_causality.py <TAG> <TAG>` → `<s>_causality_cert.csv` + `<s>_3comp_summary.json`.
12. **GATE + rolling buffer** (see §4). Map: `build_borehole_dvv_map.py` auto-discovers INCLUDEs (≥20 cert).
- **Exit codes** (fleet_station.sh): 1 no-coords · 2 download · 3 detect · 4 score · 5 adaptive · 6 cluster ·
  7 select · 8 densify/stacks · 9 dvv · 10 finalize.  FLAG-and-continue (exit 0): <100 cand, <50 families.
- **Stage logs:** `logs/{dl,cand,score,adaptive,cluster,select,densify,stack,dvv,finalize}_<s>.log`.
- **Batch:** `fleet_broadband_all.sh` (sequential, GPU-serialized via wait_gpu_free, skips STATION-DONE).
  `fleet_download_ahead.sh` pre-fetches ≤3 ahead. `fleet_trace_janitor.sh` frees DONE stations' raw traces.

## §4 GATES & TAXONOMY [FROZEN]
- **Inclusion gate (frozen 2026-07-13, rule-before-result):** INCLUDE iff **≥20 causality-certified families
  AND survival ≥15%**, where **survival = certified ÷ DENSIFIED (NSEL)**, not clustered. Else FLAG.
- **DIAGNOSIS PRECEDENCE — before calling any FLAG a "genuine weak site", check IN ORDER:**
  1. **Band** — accelerometer (HN)? wrong/stub era? (KSXB, FISH failure modes) → re-run with correct band.
  2. **Provider** — BK/NC must use NCEDC (clear `.nd` sentinels on provider change) (JCC failure mode).
  3. **Download completeness** — fragment / IRIS-sparse (HEBO, MYRA failure modes).
  4. **Only then** genuine-weak. (Most of the 22 FLAGs were mis-diagnosable as "weak" first.)
- **CLASSES:**
  - **INCLUDE** — ≥20 cert & ≥15% survival. On the map.
  - **GENUINE-WEAK** — complete data, correct band, few coherent LFEs / low survival (TKEY, CAVE, KCPB, KUHN,
    PERY 0, DOSE). Real; keep flagged.
    - **Coda-less is PER-STATION, not per-array** (corrected: I predicted the whole CC "PR" array would be
      coda-less after PR01/PR03 — WRONG. PR02 came in 131 certified @44% INCLUDE). So PR01 & PR03 are coda-less
      but PR02 is strong; don't assume siblings share it — each PR0x must be checked on its own causality tail.
      Non-PR CC stations also fine (CARB 88, RUSH 31, VOIT 69 INCLUDE).
    - **CODA-LESS sub-type** (LRIV 0/128, PR01 0/300, PR03 0/146 — CONFIRMED 3×): families cluster fine (100k dv/v rows, millions of matched det,
      cc 0.84) but the **causality distribution has NO high tail — max ratio < 1.4 across ALL families**
      (LRIV max 1.35, median 1.00). One-number tell: `max(causality) < 1.4` → coda-less site (emergent/
      microseism-dominated, no impulsive onset; LRIV = N. Olympic coast/Strait). NOT a bug — the gate is
      right. Distinct from short-tail-but-real (CINE max 3.16, 12 certified) and data-sparse (<50 families).
      DIAGNOSTIC, don't re-guess from detection count: KSXB(INCLUDE 78) has 15-23M det/yr — 10× LRIV's 1-2.7M.
      Det count and family count are NOT discriminators; the causality-tail shape is.
  - **DATA-SPARSE** — IRIS/archive holds only a fraction of the record. This is a **per-station** condition,
    NOT a network rule: CN is MOSTLY well-served at IRIS (TOFB 89, CBB 84, SNB 115, MGRB 117, TXDB 127 cert),
    only individual CN stations are sparse — **MYRA (161 days) and SHPB (48 day-files)** so far (2 of 15 CN
    processed). Don't wire NRCan mid-batch for 2 stations; bank them in `data/cn_sparse_nrcan.txt` for an
    end-of-run NRCan retry. Escalate to a real NRCan provider only if the sparse list grows large.
  - **SHORT-RECORD** — recent EHZ→HHZ upgrade, only post-2020 broadband (BHW, FISH2, PERY). Real limit.
  - **TIER-2** — ≥20 certified but <15% survival. **COMPLETED & banked** (dv/v exists & is reliable), NOT a
    failure — BUT **not yet an inversion input**: decide inclusion later (robustness check: ratio distribution
    + spatial clustering + shift-null).
    List: `data/tier2_survival_borderline.txt` (regenerate: `scripts/tier2_list.sh`).
  - **BUG-FLAG** — invalid result from a since-fixed bug; re-run. List: `data/rerun_queue.txt`.
- **Rule:** every coincidence/identity claim carries a **+120 s shift-null** (raw coincidence measures the
  tremor catalog, not LFE identity). Cap 300 = densify BUDGET; dv/v uses ALL certified families.

## §5 THE 8 BUGS [FROZEN] — do NOT reintroduce (symptom → fix → regression tell)
1. **Family selection by picker p_lfe** (anti-predictive, AUC 0.40) → rank by **SNR** (AUC 0.68), `select_..._coverage.py 'snr'`. *Tell: any selection step consulting p_lfe for ranking.*
2. **Survival = clustered denominator** → wrongly FLAGged BBO(116). Use **certified ÷ densified**. *Tell: survival % uses cluster count.*
3. **Download skipped partial fragments** (>100 files) → HEBO failed. **Always resume** + `.nd` nodata sentinel. *Tell: a conditional around the download step.*
4. **BK/NC not at IRIS** → 0 candidates. **NCEDC provider** in download_broadband + pick_band. *Tell: a BK/NC station FLAGs with 0 downloaded files.*
5. **pick_band scored per-EPOCH** → short stub era won (FISH ran 3.5 yr BHZ vs 8 yr HHZ). **Aggregate span per BAND**. *Tell: a 2-band station picks the shorter era.*
6. **pick_band treated HN (accelerometer) as broadband** → KSXB/KMR/KHMB ran strong-motion. **Require 2nd char 'H'**. *Tell: `auto-band -> HN?/BN?/EN?`.*
7. **GPU not serialized** across the sequential batch → contention. **wait_gpu_free** before GPU stages. *Tell: two densify/discover_gpu running.*
8. **appB family-stack picker gate** (`pred=='LFE'` P≥0.5) is useless/anti-predictive → **dropped; causality is the sole verdict.** *Tell: score_family_stacks_picker back in the fleet driver.*

## §6 FILE MAP & STORAGE
- **DELETABLE (janitor's job):** raw traces `data/waveforms/<NET>.<STA>/` of STATION-DONE stations only.
- **NEVER DELETE:** `logs/` (resume state) · `data/long_window_daily_*` (~152 GB daily stacks — the durable
  product) · `data/daily_dvv_*_Z_2to4.csv` · `data/<s>_causality_cert.csv` · `data/<s>_3comp_summary.json` ·
  `data/mf_<s>p90f40_*.csv` · `figures/borehole_dvv_map.html` · the two queue files. (mem: dont-delete-without-ok)
- **Required inputs:** `lfe_features/models/picker_broadband_pgc.joblib` · `catalogs/pnsn_tremor_cascadia_full.csv`
  · `data/broadband_fleet_order.csv` · `src/tremorferometry` (PYTHONPATH=src). GPU + ~32 CPU + FDSN needed.

---
# §7 — CHANGELOG (append-only; newest at bottom; historical — may contain superseded intermediate states)

## MULTI-BAND STATIONS — policy (2026-07-15)
Some stations record several vertical bands. Pick ONE band at download (§3 step 1). Prefer the longest broadband
era (native-40 BHZ ideal). If neither band spans the record and a real sensor swap splits it → PGC-style per-era
references (two segments, one figure). SHB example: BHZ40 2010-2019 (era-1) + HHZ100 2020-2026 (era-2).

## SHAKEDOWN: SHB (station 1, 2026-07-15) — DRIVER VALIDATED
CN SHB ran end-to-end autonomously → 56/331 certified (17%) INCLUDE (era-1 BHZ 2010-2019). SHB two-era like PGC;
HHZ era-2 added later (57 certified). Driver validated → fleet launched.

## 2026-07-15 — EHZ episode (SUPERSEDED by the broadband-only decision)
Batch hit short-period EHZ (SMW): adaptive floor 0.2 starved it (39 candidates); broadband picker scores EHZ low.
Chain of investigation (floor 0.2→0.05, then a wrongly-reinstated appB P≥0.5 gate) → Merlin corrected: the
family-stack gate was the bug, not picker transfer. NET RESULT: user DROPPED short-period; fleet = 178 broadband
(BH/HH) only; ~61 EHZ-only stations OUT. select_families_coverage rank-not-label fix stands (see bug 1/8 in §5).

## 2026-07-15 (late) — SELECTION FIX (bug 1): SNR rank, not p_lfe
HEBO FLAGged 10/507 (2%) — CLRS calibration showed family-stack p_lfe anti-predictive (AUC 0.404; certified 0.201
vs non-cert 0.207). SNR AUC 0.681, SNR-top-300 retains 56/56 SHB certified. Fixed → 'snr' mode; appB dropped.

## 2026-07-15 (night) — NCEDC provider (bug 4)
JCC (BK) 0 candidates from IRIS. BK+NC live at NCEDC (BK.JCC returns BHZ/HHZ from NCEDC, 204 from IRIS). Fixed
download_broadband._client(net) + pick_band; cleared false .nd sentinels. Recovered DANT(100), PETY(87).

## 2026-07-16 — BAND-PICK epoch bug (bug 5)
FISH's HHZ split into epochs → short BHZ won → ran 3.5 yr. Fixed: aggregate span per band. Survival denominator
also fixed to certified/densified (bug 2; had wrongly FLAGged BBO 116). <60 guard → <50 (borderline GOBB 56/SYMB 57).

## 2026-07-17 — ACCELEROMETER mis-pick (bug 6)
NC.KSXB picked HN? (accelerometer). pick_band checked band code not instrument code. Fixed: require 2nd char 'H'.
KSXB/KMR/KHMB (NC) ran on accel → INVALID FLAGs → re-run queue (7 total: FISH,JCC,GOBB,SYMB,KSXB,KMR,KHMB).

## 2026-07-17 — TIER-2 policy (survival-borderline)
KRP (38 certified, 13% survival) & BHW (21, 7%) pass count, miss survival. User: keep FLAGged-but-COMPLETED,
decide inversion inclusion later on evidence (their families are causality-passed = reliable; low survival =
noisy pool, not bad families). Boreholes never triggered survival-only exclusion (their flags were all count).

## 2026-07-18 — restructured this doc (Merlin); added MEMORY.md ★★★ resume pointer; materialized
`scripts/fleet_trace_janitor.sh` (was in-memory only); scope-tagged FINAL_PIPELINE as borehole-only.

## 2026-07-18→19 — batch 83→104/178 (INCLUDE 69, FLAG 34). New FLAG diagnostics settled during monitoring:
  • CODA-LESS sub-type CONFIRMED 2× (LRIV 0/128, PR01 0/300): full pipeline OK, families cluster, dv/v rows
    exist, but causality distribution has NO tail — `max(causality) < 1.4` across ALL families (LRIV 1.35,
    PR01 1.32) vs real stations 3-13. Emergent/microseism-dominated sites (LRIV=N.Olympic coast, PR01=CC
    volcano near Rainier). NOT a bug. Diagnostic is the tail SHAPE, not det/family count (KSXB 15-23M det/yr,
    INCLUDE 78; LRIV 1-2.7M, FLAG 0). → §4 GENUINE-WEAK/CODA-LESS.
  • CN is NOT a systematic NRCan case (earlier note overstated it): 10+ CN stations are strong INCLUDE from
    IRIS (TOFB 89, CBB 84, SNB 115, MGRB 117, TXDB 127); only per-station sparse ones (MYRA 161 days, SHPB 48
    day-files) need NRCan. Banked in `data/cn_sparse_nrcan.txt` for end-of-run retry; don't wire NRCan mid-batch.
  • Count-short/survival-OK near-miss = mirror of tier-2 (count-OK/survival-short). <20 cert but survival ≥15%;
    families ARE reliable, just family-supply-limited. FORMALIZED (2): TRIPT 17/84 @20%, CRBN 14/60 @23% →
    `data/count_short_survival_ok.txt`. End-of-run inversion-inclusion candidates (same logic as tier-2).
  • Top INCLUDEs this stretch: OTR 155, WISH 153, MARQ 141, JEDS 137, RADR 119, LBR 119, TOLE 117.
