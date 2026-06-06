# Nisqually / GNW dv/v — crash-recovery & resume note

**Point a fresh session here after a crash.** Last updated 2026-05-29 (session: densify speed/crash overhaul).

## Goal
Produce a **GNW coda 1–3 s dv/v plot in the same style as the existing NLLB one**
(`figures/smoke_dvv_nllb58_coda_1to3.png`), via LFE coda interferometry at station GNW
(Cascadia). Coda window 1–3 s. Standard per-patch dv/v — NOT the Nisqually dual-ref variant.
The blocker has always been the **matched-filter densification** of the 53 GNW seeds.

---

## ✅✅ PIPELINE COMPLETE (2026-05-29 ~08:27) ✅✅
**Deliverable done:** `figures/smoke_dvv_GNW_coda_1to3.png` + `data/daily_dvv_GNW_coda_1to3.csv`
(371,330 dv/v measurements, 53 patches, mean cc 0.986, span 1995-2026). Built end-to-end on the
GPU: densify all 32 years in ~17 min, stack → `data/long_window_daily_GNW.npz` (382,646 stacks),
dv/v via the NEW parallel `scripts/dvv_coda_parallel.py` (~2.5 min vs 40+ min single-threaded).
- Coverage continuous 2000-2026 (~11-16k meas/yr); 1995-1999 sparse (QC-dropped noisy early data).
- ~~Around 2001 Nisqually: dv/v fell +0.073%→+0.030%, suggestive precursor~~ **← SUPERSEDED.**
  That apparent dip was **largely an obspy-FFT-resample artifact** (see the resampler gotcha in the
  GPU-solution section). On `resample_poly` stacks the 2000–2006 dv/v is essentially FLAT — no clear
  Nisqually signal. **Any precursor test must use resample_poly stacks** (`long_window_daily_GNW_resp.npz`),
  e.g. `dvv_coda_dual_ref.py --eq-date 2001-02-28`, and should not rely on the obspy-resampled originals.
- Scale TODO (100s of stations): parallelize stack like dvv_coda_parallel; use resample_poly in
  the stacker; move the min-gap dedup off the serial main loop (GPU or vectorized).
- **⚠️ INSTRUMENT ARTIFACT at 2019 (NOT an earthquake):** GNW switched **40 Hz → 100 Hz in early
  2019** (2018 files=40 Hz, 2019=100 Hz). dv/v steps ~0.15% DOWN across that boundary (c334:
  +0.06% 2019-01 → −0.08% 2019-07) — a resampling/response artifact vs the native-40 Hz all-time-
  mean reference, not a real velocity change. **dv/v is not comparable across the 2019 boundary.**
  Fix options: per-era reference, or exclude the rate-transition, or response-correct.
  **Exact sample-rate eras (verified): 50 Hz 1995–2010, 40 Hz 2011–2018, 100 Hz 2019–2026.**
  **2000–2001 was stably 50 Hz** so the pre-Nisqually decline is NOT this artifact (still verify
  no response/metadata change near 2001-02).
  **FIX IMPLEMENTED:** `scripts/dvv_coda_perera.py` (per-era reference, parallel) — each day
  referenced to its own era's mean → 2019 step removed. Outputs `*_perera.{csv,png}`. Each era is
  centered on its own ref (no absolute cross-era level), within-era variation preserved.
  Also: `scripts/dvv_coda_perera.py --era-bounds 2011,2019` is the knob.
- Family waveform / "correlation function" viz: ad-hoc script made
  `figures/smoke_dvv_GNW_c334_corrfunc.png` (ref stack + daily-stack time-lapse image for the
  cleanest family 47.750_-123.050__c334, mean cc 0.993, 8,370 daily stacks — coherent over 30 yr).
- Window comparison done: dv/v at `--window 1.0 3.0` (mean cc 0.99, 371k meas) AND `2.0 4.0`
  (mean cc 0.88, ~117k meas — later coda much noisier, ~2/3 drop below cc>=0.8). Products:
  `data/daily_dvv_GNW_coda_2to4{,_perera}.csv` + `figures/smoke_dvv_GNW_coda_2to4{,_perera}.png`.
  **1–3 s is the higher-SNR window at GNW; 2–4 s is a sensitivity cross-check.**
- **30-day TRAILING (backward) stacks** — `scripts/build_trailing_stacks.py` (daily npz → trailing
  npz; each date = count-weighted combo of the prior N days, causal, re-normalized). Made
  `data/long_window_trail30_GNW.npz` (median 2798 det/stack vs 100 daily); per-era dv/v
  `..._trail30_perera.{csv,png}` → mean cc 0.997, much smoother, step-free.
- **INSTRUMENT RESPONSE REMOVAL** — inventory saved `data/UW.GNW.response.xml` (fetched IRIS FDSN;
  confirms channel swap **BHZ→HHZ 2019-05-07** + earlier response epochs = physical cause of the
  2019 step). `scripts/build_long_window_resp.py` = stacker with `remove_response(VEL, pre_filt
  [0.5,1,15,18])` + spawn (fixes the 135 GB fork blowup) + resample-FIRST-then-deconv (~2.75×
  faster: 138 s/yr; full record ~1 h → `data/long_window_daily_GNW_resp.npz`). After it builds,
  run dv/v (ALL-TIME ref, no per-era) to test if response removal flattens 2019 on its own.
  Re-fetch inventory if missing: `Client("IRIS").get_stations(network="UW",station="GNW",
  channel="*HZ",level="response").write("data/UW.GNW.response.xml","STATIONXML")`.
- **EPOCH-AWARE fix (no response removal needed for the step):** the "massive peak at the 2019
  instrument swap" in the trailing dv/v was the **30-d trailing window STRADDLING the BHZ→HHZ
  change** (mixing two instruments into one stack) + the per-era boundary being on calendar-year
  2019 instead of the real **2019-05-07**. Fixed both: `build_trailing_stacks.py --epoch-bounds
  "2010-09-10,2019-05-07"` (window never pools across a change) and `dvv_coda_perera.py
  --era-bounds "2010-09-10,2019-05-07"` (date-based, exact). Peak GONE without response removal →
  it was a straddle artifact, not signal. Products: `data/long_window_trail30ea_GNW.npz`,
  `figures/smoke_dvv_GNW_trail30ea_1to4_perera.png` (1–4 s, mean cc 0.995). Exact transition
  dates from inventory: 50Hz<2010-09-10, 40Hz 2010-12-01..2019-05-07, 100Hz HHZ >=2019-05-07.
- **RESULT — response removal does NOT flatten the instrument steps (tested both ways):**
  all-time-ref dv/v on response-removed stacks shows sharp era offsets (50Hz +0.04%, 40Hz +0.26%,
  100Hz −0.39%); the 2019 BHZ→HHZ step is **~−0.65%** with BOTH the fast resample-first deconv
  AND the correct native-rate deconv (`--deconv-native`, tested on 2017–2021: 40Hz +0.333%, 100Hz
  −0.317%, step −0.650%). So my resample-first speedup was NOT the cause — deconvolution genuinely
  can't remove it (residual per-instrument PHASE / response-metadata imperfections → uniform coda
  stretch). ALL 53 patches step in unison → systematic/instrumental, not geology. **CONCLUSION:
  instrument epochs are not cross-comparable in absolute dv/v; PER-ERA referencing (compare each
  epoch to its own mean) is the necessary fix — response removal alone does not stitch eras.**
  Implication for 100s-of-stations: do dv/v per-instrument-epoch. Per-era response-removed product:
  `dvv_coda_perera.py --npz long_window_daily_GNW_resp.npz --era-bounds 2010-09-10,2019-05-07`.
- **★★ BIGGEST GOTCHA — obspy FFT resample injects a SPURIOUS dv/v drift (use scipy resample_poly) ★★**
  Decisive test: within the SINGLE stable response epoch 2000-06-21..2006-11-01 (constant response,
  self-referenced — so nothing instrumental can vary in time), c334 dv/v differs ONLY by resampler:
  **obspy `tr.resample` swings ~0.16% (+0.095% in 2000 → −0.067% in 2003); scipy `resample_poly`
  is FLAT (±0.02%).** Response is constant there, so deconvolution is a common filter (no relative
  effect) → the swing is PURELY the resampler. A real velocity change would appear in BOTH; only
  obspy drifts → **obspy's FFT resample is the artifact** (ringing on gappy/zero-filled data).
  CONSEQUENCES:
  - The apparent "2001 Nisqually dip" + 2002–2005 decline were **largely obspy-resample artifacts**.
    With resample_poly the 2000–2006 dv/v is **flat → no clear Nisqually dv/v signal** once fixed.
  - `build_long_window_daily_all51.py` uses obspy `tr.resample` (line ~63) → the ORIGINAL stacks
    `long_window_daily_GNW.npz` and everything derived (original deliverable `daily_dvv_GNW_coda_1to3`,
    `trail30ea`, all no-resp plots) are CONTAMINATED by this drift. **`build_long_window_resp.py`
    uses `resample_poly` → its `long_window_daily_GNW_resp.npz` is clean.**
  - **ACTION: always use resample_poly for dv/v stacks.** Trust the resample_poly products. The
    canonical clean base = `long_window_daily_GNW_resp.npz` (resample_poly + response removal); a
    `--no-deconv` resample_poly build would isolate resampler vs response if ever needed.
  - This supersedes the earlier "suggestive 2001 precursor" framing in this doc — that dip was
    mostly the obspy artifact. Re-evaluate the precursor ONLY on resample_poly stacks.
  Demonstrated in `figures/smoke_dvv_GNW_c334_stableepoch.png`.
- **★★ NISQUALLY VELOCITY-DROP — FINAL VERDICT (after peeling away every artifact) ★★**
  Question: is there a coseismic dv/v drop at the 2001-02-28 Nisqually M6.8 at GNW? Answer after
  full analysis: **a real, REGION-WIDE ~0.04–0.05% coseismic-shaped drop exists (min at the EQ,
  recovery), in clean resample_poly+per-era+deseasonalized data — present both in the NW set AND in
  the 22 near-epicenter families. Consistent with a coseismic origin; significance vs the noisy
  early-era background is the one thing still UNTESTED (needs the null/false-positive test).**
  Chain of reasoning (each step removed an artifact that had inflated/created an apparent signal):
  1. Original all-time-ref dv/v showed a big 2000→2001 decline → looked like a precursor.
  2. INSTRUMENT STEPS (50/40/100 Hz) — all-time ref blends eras; per-era referencing needed.
  3. OBSPY RESAMPLE ARTIFACT — the original stacks' big 2002–2005 decline was an obspy-FFT-resample
     artifact (resample_poly is flat in the same stable-response window). Use resample_poly.
  4. ALL-TIME-REF INFLATION — referencing 50 Hz data to a later-instrument-dominated mean inflated
     the apparent dip; per-era referencing shrinks it.
  5. After all that, on clean (resample_poly + response-removed + per-era) data: cross-patch median
     drops +0.022%→−0.020% at the EQ (Feb–Mar 2001), recovers — coseismic SHAPE, ~0.04%.
  6. NOT SEASONAL: GNW seasonal cycle is only ~0.006% peak-to-peak (2011–2026 climatology);
     deseasonalized, the dip survives (+0.027%→−0.017%).
  7. BUT NOT ANOMALOUS vs noise: deseasonalized monthly cross-patch median std = 0.094% (noisy
     early/50 Hz era); the 2001 dip ranks only 20th-lowest of 313 months; comparable wiggles recur
     2002–2004 with no EQ.
  8. SPATIAL TEST: per-family deseasonalized drop, mapped with station→patch paths. FIRST PASS
     (SNR≥15, only ~6 near-epicenter families) looked NW-clustered/away-from-rupture → I wrongly
     concluded "not coseismic." **RETRACTED:** that was UNDER-SAMPLING the south. Re-ran with the
     **22 SNR≥12 families within 40 km of the epicenter** (`scripts/run_gnw_circle.sh` →
     `data/daily_dvv_GNWcircle_perera.csv`, `figures/smoke_dvv_GNWcircle_perera.png`): the
     near-epicenter families ALSO show the ~0.05% 2001 dip (deseas: +0.028% Sep2000 → −0.019% Feb2001
     → recovery). So the drop is **REGION-WIDE** (NW + near-epicenter), consistent with a coseismic
     response (M6.8 shakes the whole region), NOT NW-only. Updated map: `smoke_dvv_GNW_pathmap_nisqually_v2.png`.
     LESSON: GNW LFE families are 5:1 N:S (deep ETS band under the Olympics); at SNR≥15 the south is
     barely sampled — drop to SNR≥12 (150 families) to cover the Nisqually source region.
  Caveats: weak per-family significance & spatial corr (r≈0.2); deep plate-interface LFEs (paths go
  to ~35–40 km sources, map-view ≠ shallow epicentral sampling); n=1 event (no false-positive rate).
  METHODS/figures: deseasonalize = subtract each family's 2011–2026 monthly climatology; per-family
  drop = deseas median(2001-03..08) − median(2000-09..2001-02). Figures: `smoke_dvv_GNW_deseasonalized.png`
  (cross-patch median, full+zoom), `smoke_dvv_GNW_perfamily_deseason.png` (all families deseas),
  `smoke_dvv_GNW_nisqually_dropmap.png` (per-family drop map), `smoke_dvv_GNW_pathmap_nisqually.png`
  (station→patch ray paths colored by drop; uses cartopy). Data: `daily_dvv_GNW_resp_perera_1to4.csv`.
  Remaining rigorous test if revisited: NULL/false-positive test (dv/v change across 100s of random
  non-EQ dates — what fraction give ≥0.04% drop?) + dual-ref with error bars on resample_poly stacks.

## ★★ GPU SOLUTION (2026-05-29 late session) — USE THIS ★★
**There is an idle NVIDIA L40S (46 GB) on this pod.** GPU matched filter is **bit-identical**
to the CPU path (`max|Δcc|=0.0`, same detections) and **~180× faster** (44 ms/day vs 7.9 s/day).
This is the answer (also for the planned 100s-of-stations scale).

- **cupy installed**: `pip install --user "cupy-cuda12x<14"` → cupy **13.6.0** (cupy 14 drags
  numpy 2.x and breaks scipy/numba — see [[pod-env-ephemeral]]; 13.x is numpy-1.26-safe).
  CUDA 12.4 driver, cuFFT wheels already present.
- **Scripts (built this session):**
  - `scripts/densify_gnw_gpu.py` — CPU worker pool (read+resample_poly+bandpass) feeds the GPU
    (one rfft/day + batched template multiply + batched irfft + GPU sliding-norm + cc>1 fix +
    threshold + sparse-candidate extraction; only detections cross PCIe). Per-year chunks +
    resume (skip years whose `data/mf_gnw_YYYY.csv` exists; atomic .tmp+rename). Per-day QC
    BEFORE top-N cap. **GPU only imported in main(); pool workers never touch cupy.**
  - `scripts/run_gnw_pipeline.sh` — densify(retry-resume) → concat year csvs → stack → dv/v plot.
- **Validated**: year 2015 = 365 days in **0.53 min (11.6 days/s)**, 1.41M rows, cc∈[0.8,1.0],
  1 QC-dropped. Full 10,662-day record projects to **~15–20 min**.
- **RESUME COMMAND after a crash** (skips finished years automatically):
  ```bash
  cd /home/jovyan/tremorferometry
  export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  bash scripts/run_gnw_pipeline.sh    # or just: python scripts/densify_gnw_gpu.py --workers 12
  ```
  Then check `logs/gnw_gpu_pipeline.log`. Output: per-year `data/mf_gnw_<YEAR>.csv`, then
  `data/mf_gnw_all.csv`, `data/long_window_daily_GNW.npz`, `figures/smoke_dvv_GNW_coda_1to3.png`.

The CPU-only fixes below are the FALLBACK if the GPU is ever unavailable.

### RUN STATUS (2026-05-29 ~07:15, autonomous run while user asleep)
- Pipeline running via `run_gnw_pipeline.sh` (background). Per-year files in `data/mf_gnw_YYYY.csv`.
  Done so far: 1995, 1996, 1997, 1998, 2015. Resumes by skipping existing year files.
- **Speed gotcha (TODO for 100s-of-stations scale):** the GPU correlation is 44 ms/day, but the
  MAIN process is the bottleneck — it runs `_dedup_min_gap` (pure-Python min-gap peak-pick) +
  QC + cap serially per day, so gappy 50 Hz years run only ~1.5 days/s (12→24 workers barely
  helped — workers are NOT the limit). Fast 40 Hz years (2010-2019) should hit ~10/s.
  To scale: vectorize the dedup or do peak-pick on the GPU; overlap GPU with CPU dedup.
- **QC dropping ~76-83% of 1995-1998 days** (station's noisy first years; clean 2015 dropped 1/365).
  WATCH whether 2000-2001 (Nisqually precursor baseline) also drop heavily — if so the pre-EQ
  baseline is sparse and QC (`--qc-median-cc 0.96`, `--qc-median-count 2000`) may need relaxing.
- dvv_coda_51.py now takes `--station` (title fix); run final plot with `--station GNW`.
- **DENSIFY DONE: all 32 years (1995-2026) in ~17 min**, `data/mf_gnw_all.csv` = 32,993,627 rows
  (2.5 GB). 2000-2026 are CLEAN (2001 dropped only 15/347 days) — Nisqually baseline is solid.
- **STACK step memory gotcha:** `build_long_window_daily_all51.py` loads the full 33M-row CSV
  THEN forks 30 workers → COW + Python refcounting balloons each worker to ~8 GB → ~135 GB
  total (72% of the 187 GB cap). It held stable but it's close. **If the stack ever OOMs, re-run
  it ALONE with fewer workers (densify is saved):**
  ```bash
  python scripts/build_long_window_daily_all51.py --mf-csv data/mf_gnw_all.csv \
    --network UW --station GNW --out data/long_window_daily_GNW.npz --workers 10
  python scripts/dvv_coda_51.py --npz data/long_window_daily_GNW.npz --window 1.0 3.0 \
    --station GNW --out-csv data/daily_dvv_GNW_coda_1to3.csv \
    --out-fig figures/smoke_dvv_GNW_coda_1to3.png
  ```
  (For scale: the stacker should use spawn + load detections lazily per-worker, and resample_poly.)

---

## ★ SPEED + CRASH FIXES (CPU fallback) — the heart of this session ★

User goal set verbatim: *"densify as fast as possible without crashing the terminal.
densify should pay attention to removing anomalous days, fix the cc>1 bug, cap the detections."*

### Why it was slow (ROOT CAUSE = pod CPU cap, not the algorithm)
- **The pod is capped at 32 CPUs.** `cat /sys/fs/cgroup/cpu.max` → `3200000 100000` →
  quota/period = **32**. BUT `os.cpu_count()` / `nproc` / `sched_getaffinity` all report
  **176** (the *node's* cores) — they LIE about the pod's real budget.
- Past runs used **`--workers 48`** → oversubscribed the 32-CPU quota → CFS throttling +
  context-switch thrash → only ~2.4× effective speedup instead of ~30×. That is why a
  ~10.7 s/day job projected to **~13 h**. **Use `--workers 30` (≤32). This is the #1 fix.**
- Also set thread env so each worker uses ONE thread (else 30 workers each spawn FFT/BLAS
  threads → re-oversubscribe): `export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
  OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` before python starts.

### Per-day cost breakdown (profiled on a 100 Hz GNW day, single thread)
- read 0.10 s | **resample 100→40 = 2.77 s** | merge+detrend+bandpass 0.69 s |
  sliding_norm 0.12 s | **correlation (53 templates) = 6.84 s** | total ≈ 10.5 s.

### Algorithmic fixes (validate in `/tmp/test_fast.py`, then port to the module)
1. **`resample_poly` instead of obspy `tr.resample`** — obspy's FFT resample = 2.74 s/day;
   scipy `resample_poly(x, up, down)` = **0.11 s (25× faster)**, keeps fs=40 (NLLB-comparable).
   Ratios: 100→40 = up2/down5, 50→40 = up4/down5, 40→40 = skip. Apply per-trace BEFORE merge;
   set `tr.stats.sampling_rate = 40` after. Fall back to obspy for weird/jittery rates.
2. **Reuse the day FFT.** `scan_day_multi` (matched_filter_fast.py:117) calls
   `fftconvolve(data, t0[::-1])` **inside the per-template loop** → recomputes the ~3.5M-pt
   day FFT 53× per day. Fix: `D = rfft(data, N)` ONCE (N = `scipy.fft.next_fast_len(n+m-1)`),
   precompute template FFTs `Tk` once (cache by N — full days share one N), then per template
   `irfft(D*Tk[i], N)[m-1:m-1+(n-m+1)]`. ~1.5–3× on correlation. (Naive reuse with a non-fast N
   was SLOWER — must use next_fast_len.) Looped irfft (not a batched (53,N) array) keeps RAM low.
3. **cc>1 fix** — already in code (matched_filter_fast.py:126-127): `cc[~isfinite]=0;
   cc[|cc|>1]=0`. KEEP. (Zero-filled gaps made win_std→1e-12 floor → cc exploded; GNW was
   99.9% cc>1 artifacts before this.)

### Expected result: ~30 workers under the cap + resample_poly → ~100–300 days/min →
full 10,662-day GNW span in **~1 h (vs ~13 h)**, memory flat.

### The LOCKED densify design (what to build / rebuild `densify_nllb_seeds.py` + `scan_many_days_multi` to do)
- **`--workers 30`**, thread env pinned to 1 (above).
- **resample_poly + day-FFT reuse** (above).
- **forward-chronological** (NOT `--reverse`) so the 2001 Nisqually era is reached early and
  durable. (`--reverse` + delete-on-start = never reaches 2001; that was the old trap.)
- **YEARLY chunks + resume**: write `data/mf_gnw_YYYY.csv` per year, accumulate that year in RAM
  (fine now — see cap), one `to_csv`, free. **Resume = skip years whose file already exists**
  (glob, no parsing). A year (~27 min modern / faster old) ≈ inside a crash window; lose ≤1 year
  per crash. Stacker reads all year files (or concat at end).
- **top-N = 100 cap per (family, day)**: keep the 100 highest-cc detections. Stacks only need
  ~20 (min-det), √N saturates ~100; the cap drops the noisy marginal cc≈0.80 ones FIRST so a
  capped stack is often *cleaner*. **Quiet/pre-2001 family-days (<100 det) are untouched →
  coverage preserved.** Shrinks output ~4× (2.2M→0.58M on the 149-day partial).
- **Per-day QC (anomalous days) — KEEP, and run it on the FULL day BEFORE the cap.** Gates:
  median per-template count > 2000 OR median cc > 0.96 → drop the whole day (telemetry/instrument
  glitch saturates all 53 families near-equally at cc~0.98, e.g. GNW 1995-07-17). **Order matters:
  if you cap first, every family shows 100 det and the count-gate signature vanishes.** Capping
  does NOT clean artifacts (it'd keep 100×53 of cc~0.98 garbage that pollutes the stack) — only
  dropping the day removes it.
- **NO streaming, NO bounded-submission.** Those were band-aids for the millions-of-detections
  firehose; the cap removes the root cause. With the cap GNW is ~10–30M rows total (few GB RAM).
- **cc threshold stays 0.8 at densify** (do NOT raise to 0.85+). Stacker groups per (date,template)
  with **min-det=20 PER FAMILY** (`build_long_window_daily_all51.py:131`). At cc≥0.8, 86% of
  family-days have ≥20 det; cc≥0.85 → 52%; cc≥0.9 → 3% (collapses), worst in quiet/pre-2001. So
  thin via the top-N cap, not by raising cc. (Adaptive ">2000/family→0.85" ≈ no-op: max ever seen
  was 2085/family/day.) If you later want stricter, set the *stacker's* `--cc-min 0.85` — reversible,
  no MF re-run; the 0.80–0.85 band on disk is insurance for sparse pre-2001 days.

### This crash (2026-05-29) = pod EVICTION, not OOM
Died at **77 GB** (cap `memory.max` ≈ 187 GB) with a full pod restart (new jupyterhub 04:55) —
a cgroup OOM would kill only python and leave the pod up. The *recurring historical* crashes WERE
OOMs (the firehose), now fixed by the cap. Eviction cadence is ONE data point (~20 min) — not a
proven rate; the pod may be stable for hours. Yearly chunks cover us either way.

### Why NLLB ran much faster than GNW
NLLB (CN.NLLB) was **40 Hz native for 2005–~2020** → skipped resampling for most of its record.
GNW (UW.GNW): **50 Hz pre-~2010 (gappy!), 40 Hz ~2010–2019, 100 Hz 2020+** → pays resampling
almost everywhere. After resample everything is 40 Hz so the correlation cost is identical; the
only per-day difference is the (now `resample_poly`-cheap) resample.

---

## Current state on disk
- `data/gnw_pnsn_candidates.parquet` + `data/gnw_pnsn_families.npz` + `.summary.csv` (has `snr`)
  + `.members.parquet` — 20,612 families/seeds, built by `scripts/discover_nllb_pnsn_driven.py
  --station GNW`. **SNR≥15 → 53 seeds.** ✓
- `data/waveforms/UW.GNW` — 10,662 day-files, **1995-06-29 → 2026-05-28**, network **UW**. ✓
- `data/mf_gnw_seeds_top.csv` — **STALE/PARTIAL**: only the newest **149 days**
  (2025-12-30…2026-05-28), `--reverse`, pre-cap. cc≤0.999 (post-fix, clean) but incomplete and
  superseded by the new yearly-chunk design. Delete/ignore; new run writes `mf_gnw_YYYY.csv`.
- `figures/smoke_gnw_family_map.png` — 53-seed map, done. ✓
- NLLB reference (done, but mildly cc>1-contaminated — ran pre-fix): `data/mf_nllb_seeds_top.csv`
  (12 GB, cc>0.7), `mf_nllb_seeds_cc08.csv`, `long_window_daily_nllb58.npz`,
  `daily_dvv_nllb58_coda_1to3.csv`, `figures/smoke_dvv_nllb58_coda_1to3.png`.
- `/tmp/test_fast.py` — validation harness for the optimized day-scan (correctness vs current +
  single-day speed + parallel throughput ≤32 workers). RUN THIS before the final densify.

## DATE SPAN — RESOLVED
User directive "densify every year" → did the FULL record 1995–2026 (covers the 2001 precursor).
The earlier 2010–2026-vs-1999 question is moot. DONE.

## Re-run the pipeline (scripts generic; all built/validated this session)
The whole chain is `scripts/run_gnw_pipeline.sh` (densify→concat→stack→dv/v). To re-run pieces:
1. Densify (GPU): `python scripts/densify_gnw_gpu.py --workers 24` → per-year `data/mf_gnw_YYYY.csv`
   (cc≥0.8, top-100/family-day, QC). Resumes by skipping existing year files.
2. Concat + Stack: concat `mf_gnw_[12]*.csv` → `data/mf_gnw_all.csv`, then
   `python scripts/build_long_window_daily_all51.py --mf-csv data/mf_gnw_all.csv --network UW
   --station GNW --out data/long_window_daily_GNW.npz --workers 10` (use ≤10 — see STACK memory gotcha).
3. dv/v (PARALLEL, the fast one): `python scripts/dvv_coda_parallel.py --station GNW --window 1.0 3.0
   --out-csv data/daily_dvv_GNW_coda_1to3.csv --out-fig figures/smoke_dvv_GNW_coda_1to3.png`.
   (The serial `scripts/dvv_coda_51.py` does the same math but took 40+ min — use the parallel one.)
4. (Nisqually precursor extra — RECOMMENDED NEXT) `python scripts/dvv_coda_dual_ref.py
   --npz data/long_window_daily_GNW.npz --eq-date 2001-02-28` (pre-EQ-referenced dv/v).

## Env gotcha
`/opt/conda` is EPHEMERAL overlay — a crash wipes pip-installed pkgs. Install with
`pip install --user <pkg>` (lands in `~/.local`, persists). Do NOT `--ignore-installed`.
Working pin: numpy 1.26.4 (conda), scipy 1.14.0, cartopy 0.25.0, obspy 1.5.0 (~/.local).
This session env was INTACT (conda not wiped). Work is done in a REPL — no `.ipynb` — so
in-flight commands don't survive a crash; only on-disk artifacts + this note do.

## HDW multi-station run — GPU discovery + station-centered box (this session)
Goal: run the full pipeline on UW.HDW (Olympic Pen., 47.649, -123.053), coverage-select
families, densify, dv/v. New tooling built + validated this session:

- **GPU family discovery (the big win).** `scripts/discover_gpu.py` replaces the slow/OOM-y
  CPU stage-2. Three speedups, all validated equivalent:
  1. LOAD-EACH-DAY-ONCE (group candidate cuts by calendar day; ~11k loads not ~100k).
  2. `resample_poly` not obspy FFT-resample (0.12 s vs 1.1 s/100Hz-day; far less RAM).
  3. GPU all-pairs max-shifted CC (`src/tremorferometry/repeater_gpu.py`): 11.6 s → 0.06 s
     per 2000-window bin, max|Δ|=1.8e-7 vs CPU, ZERO cluster-flips at cc≥0.8.
  Single GPU process; CPU workers only cut waveforms (KB each) → no N²-per-worker blowup.
  HDW stage-2: 29,721 families in **89 s at anon 2 G** (CPU run was ~100 min / 54 G / 50-cap).
  **It writes the `snr` column** (template env-peak / pre-pulse-RMS, pre = T[:34]) that the
  committed `discover_nllb_pnsn_driven.py` silently STOPPED writing — without it,
  select_coverage_families.py / plot_family_map.py break at the next step.
- **Stage-1 also fixed for memory.** `discover_nllb_pnsn_driven.py` now has `_load_day_poly`
  (resample_poly) used by `_candidates_one_day`, + `--candidates-only` flag (run stage-1,
  save parquet, exit; GPU does stage-2). 24-worker obspy stage-1 hit anon 99 G (near-OOM,
  killed); resample_poly + 16 workers → anon **4 G**, 1543 days in ~95 s.
- **PNSN catalog was DOWNLOAD-CUT at 47.5°N** (the `--bbox 47.5 ...` in 00b_fetch). HDW sits
  ~16 km north of that false floor, so the south looked empty. The MASTER catalog already on
  disk is `catalogs/pnsn_tremor_cascadia_full.csv` (748k rows, lat 39.5–50.9, 2010–2026) —
  USE THIS, not the 47.5-cut `pnsn_tremor_2014-2026.csv`.
- **Station-CENTERED box** (per user): center on the station, choose ±km. HDW used ±100 km
  N/S+E/W → bbox `46.748 48.550 -124.390 -121.716`. 43% of within-60km tremor is SOUTH of HDW.
  Plot the box first: `scripts/plot_station_box.py`.

HDW pipeline state (DONE up to selection):
  data/hdw_pnsn_candidates_100km.parquet (1.72M, 34% south)
  data/hdw_pnsn_families_100km.npz/.summary.csv (29,721 fams; SNR≥10: 96, 23 south)
  data/hdw_coverage_selection_100km.summary.csv (37 fams, 10/12 az, 4/4 rings, med SNR 19.8)
  figures: smoke_hdw_family_map_100km.png, smoke_hdw_coverage_selection_100km.png
NEXT: densify the 37 selected at HDW (GPU matched filter, full record) → stack → coda dv/v,
same chain as GNW (densify_gnw_gpu.py is the template; point it at the 37-family templates).

Generic multi-station recipe (for the 100s-of-stations goal):
  1. download_hdw.py-style waveform pull.  2. discover_nllb_pnsn_driven.py --candidates-only
  --pnsn cascadia_full --bbox <station-centered>.  3. discover_gpu.py (GPU stage-2 + snr).
  4. select_coverage_families.py.  5. plot_family_map.py.  6. densify (GPU).  7. dv/v.
