# MARGIN_WORKFLOW.md — authoritative playbook for whole-Cascadia coda dv/v

**Purpose:** run the LFE-coda dv/v workflow on EVERY useful station along the Cascadia margin.
This doc is self-contained: if the session resets, READ THIS FIRST and you can resume blind.
History/lessons live in `notes/2026-05-29_Notes.md`, `notes/2026-05-30_Notes.md`,
`notes/NISQUALLY_GNW_RESUME.md`. Memory index: `~/.claude/.../memory/MEMORY.md`.

Goal: for each station, produce a per-era, response/instrument-QC'd coda dv/v time series and a
coverage-balanced family map. **Done (2026-06-06): 16 stations** — GNW/HDW/CPW/PGC/NLLB + 11 PB boreholes
(B018/B023/B941/B013/B011/B928/B004/B026/B014/B204/B028). See [[margin-station-status]] memory.

---
## 0a. REFINEMENTS SINCE THIS DOC WAS FIRST WRITTEN (2026-06 — apply on top of the steps below)
- **ENV:** use `/home/jovyan/envs/tremorferometry/bin/python` (py3.11, full stack), NOT base conda (py3.13, no obspy).
  For any GPU step (discover_gpu, densify): `export CUDA_PATH=/opt/conda/targets/x86_64-linux` or cupy's JIT fails
  (`cannot open cuda_fp16.h`). See [[pod-env-ephemeral]].
- **DESPIKE:** pass **`--despike-mad 8`** to BOTH `densify_gnw_gpu.py` AND `build_long_window_resp.py` (winsorize
  beyond 8·MAD/day). Suppresses glitch-contaminated station-years that otherwise QC-drop ~90% of days (PGC 2013-14).
  Default off so other stations unaffected; on for all new runs. See [[glitch-days-false-dvv-drop]].
- **FAMILY SELECTION:** coverage-balanced (top-2/az×dist cell) **+ add back the top-10% highest-SNR families** of the
  eligible pool (keep strongest repeaters). Do it UPFRONT before densify. See [[family-selection-rule]].
- **min-det STATION-DEPENDENT, decide UPFRONT from density** (while traces still on disk → no re-download): check
  `mf.groupby([template,day]).size().median()`; set min-det ≈ that — ~20 rich central, ~8 moderate, **~5 sparse**
  (southern OR/CA, distant eastern). Match threshold to per-day firing rate or dv/v comes out gappy.
- **GAPS:** diagnose before "filling" — **multi-month blocks = station OUTAGES (unfixable**, e.g. B202 dropped);
  scattered single days = min-det starvation (fixable by lower min-det). Check gap structure first.
- **dv/v plots:** auto-scale y-limits to each station's amplitude (robust ~1.7× p98 of cross-patch median, clamp
  [0.025,0.12]); don't hard-code ±0.05. Most PB boreholes are single EHZ era → all-time ref (`dvv_coda_parallel`).
- **PATH MAP:** `scripts/plot_pb_path_map.py` — station→coverage-families map, auto-discovers done PB stations.
- **SLAB DEPTHS:** `data/station_slab2_depth.csv` (Slab2 interface depth per station; ScienceBase grid is firewalled,
  interpolated from input DB `data/cas_slab2_input_04-18.csv`).
- Background jobs: ALWAYS launch tracked (not bare `&`) so completion notifies; remember CUDA_PATH on GPU jobs.

---
## 0. ENVIRONMENT & HARD CONSTRAINTS (read every session)
- **Pod is capped at 32 CPUs** (`/sys/fs/cgroup/cpu.max`=`3200000 100000`); `nproc`/`os.cpu_count()` lie (176).
  Keep total workers ≤30 across all running jobs or pools thrash. Pin threads: `export OMP_NUM_THREADS=1
  MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.
- **GPU:** 1× NVIDIA L40S (46 GB), cupy-cuda12x **<14** (cupy 14 pulls numpy 2 → breaks scipy/numba). Only ONE GPU job at a time.
- **Memory watchdog:** watch `anon` (real, OOM-risk), NOT `memory.current` (incl. reclaimable page cache).
  Kill any job if `awk '/^anon /{print $2}' /sys/fs/cgroup/memory.stat` > ~120 GB.
- **Env is ephemeral:** a crash wipes /opt/conda; `pip install --user` to persist. Work is in a REPL — only
  on-disk artifacts survive a crash.
- **`PYTHONPATH=src`** for all scripts (they `import tremorferometry...`).
- **Always `resample_poly`, never obspy FFT resample** (obspy injects a spurious dv/v drift + OOMs). All
  pipeline scripts already use resample_poly EXCEPT the OLD `build_long_window_daily_all51.py` (obspy →
  contaminated; its `long_window_daily_GNW.npz` is tainted — don't merge it with resample_poly stacks).
- **Monitor loops must use the bracket trick** (`pgrep -f "[d]ensify..."`) or they self-match and hang forever.
  Never put a `pkill -f X` in the SAME shell command as a launch that spells out `X` (it kills its own launcher).

## 1. CATALOG (one-time, already done)
Master tremor catalog = **`catalogs/pnsn_tremor_cascadia_full.csv`** (748k rows, lat 39.5–50.9, 2010–2026,
full Cascadia). USE THIS. The older `pnsn_tremor_2014-2026.csv` is **clipped at lat 47.5N** — do not use it.
Re-fetch wider/newer if needed: `scripts/00b_fetch_pnsn_tremor.py --start .. --end .. --bbox <wide> --out ..`.
(Tremor catalog starts ~2010 → families are 2010+ LFEs; still matched-filtered across the full waveform record.)

---
## 2. PER-STATION PIPELINE (the core loop)
For station `<STA>` on network `<NET>` at `(LAT, LON)`. Replace placeholders. Run steps in order.

### 2a. SCOUT first (mandatory — never skip)
```
PYTHONPATH=src python -c "
from obspy.clients.fdsn import Client
c=Client('EARTHSCOPE')
inv=c.get_stations(network='<NET>',station='<STA>',level='channel',channel='?HZ',starttime='1980-01-01',endtime='2027-01-01')
for net in inv:
 for s in net:
  print(s.code, s.latitude, s.longitude)
  for ch in s:
   if ch.code.endswith('Z'): print(' ',ch.code,round(ch.sample_rate),'Hz',ch.start_date,ch.end_date)"
```
Record: **coords**, **vertical channel** (prefer EHZ/HHZ/BHZ; borehole PB stations = EHZ), **sample rate(s)**,
and **data span incl. END date** (many UW short-period sites are decommissioned; PB borehole B-codes are the
active replacements). NOTE every sensor/rate change date — these become the per-era boundaries (step 2g).
"to current" is FALSE if end_date is closed.

### 2b. Station-centered box (±100 km default)
`dlat = 100/111 = 0.901`; `dlon = 100/(111*cos(LAT_rad))`. bbox = `LAT-dlat LAT+dlat LON-dlon LON+dlon`.
Preview: `python scripts/plot_station_box.py --station <STA> --station-lat LAT --station-lon LON --ns-km 100 --ew-km 100 --out figures/smoke_<sta>_box.png`

### 2c. Download waveforms
```
PYTHONPATH=src nohup python scripts/download_station.py --network <NET> --station <STA> \
  --start <YYYY-MM-DD> --end <YYYY-MM-DD> --workers 8 --client EARTHSCOPE >> logs/<sta>_download.log 2>&1 &
```
Resumable; picks one vertical channel/day (EHZ>HHZ>BHZ>SHZ); fast (~12–40 day/s). Layout: `data/waveforms/<NET>.<STA>/<year>/<jday>.mseed`.

### 2d. Stage-1 candidate detection (after download)
```
PYTHONPATH=src python scripts/discover_nllb_pnsn_driven.py --station <STA> \
  --pnsn catalogs/pnsn_tremor_cascadia_full.csv --bbox <BBOX> \
  --candidates-out data/<sta>_pnsn_candidates_100km.parquet --candidates-only --workers 16
```
(`--candidates-only` exits after saving; uses resample_poly loader — 16 workers, anon ~4 G. 24+ workers OOMs.)

### 2e. GPU stage-2 clustering
```
PYTHONPATH=src python scripts/discover_gpu.py --station <STA> \
  --candidates data/<sta>_pnsn_candidates_100km.parquet \
  --out data/<sta>_pnsn_families_100km.npz --max-bin-candidates 2000 --workers 24
```
Writes `.npz` + `.summary.csv` (WITH `snr` col = env-peak/pre-pulse-RMS) + `.members.parquet`. ~1–2 min, anon ~2 G.

### 2f. Coverage selection + candidate maps → PAUSE for user confirm
**NEW (2026-06-10) — TEMPLATE-SHAPE PRE-SCREEN before selection.** Drop noise-like / contaminant families
(they come out gappy/low-cc = wasted densify) using the discovery-stage templates in `<sta>_pnsn_families_100km.npz`.
Verified rule (grpCV AUC 0.92, see [[family-cwi-predictor]] / `scripts/family_predictor.py`): flag families with
**spectral centroid > 4.3 Hz AND kurtosis > 4** as BAD and exclude from the eligible pool BEFORE coverage+top-10%
selection (97% precision, removes ~47% of BADs, loses 0.4% of GOODs). GOOD LFEs are low-frequency (centroid ~3.4 Hz)
and emergent (low kurtosis); contaminants are high-freq and spiky. Compute centroid/kurtosis per template from the
npz, intersect the kept set with the SNR-floored pool, then proceed:
```
PYTHONPATH=src python scripts/select_coverage_families.py --summary data/<sta>_pnsn_families_100km.summary.csv \
  --station-lat LAT --station-lon LON --min-snr 10 --az-sectors 12 --dist-rings "0,20,40,65,100" --k 2 \
  --out data/<sta>_coverage_selection.summary.csv --out-fig figures/smoke_<sta>_coverage_selection.png
PYTHONPATH=src python scripts/plot_family_map.py --station <STA> --station-lat LAT --station-lon LON \
  --summary data/<sta>_pnsn_families_100km.summary.csv --min-snr 10 --out figures/smoke_<sta>_family_map.png
```
After the top-10% add-back merge (0a), draw the SELECTION map (greyed pool + highlighted final selection, same
cartopy basemap as the family map — use this, not select_coverage_families.py's bare scatter, which predates the
top-10% merge and lacks coastlines):
```
PYTHONPATH=src python scripts/plot_selection_map.py --station <STA> --station-lat LAT --station-lon LON \
  --summary data/<sta>_pnsn_families_100km.summary.csv --selected data/<sta>_coverage_selection.summary.csv \
  --min-snr <floor> --out figures/smoke_<sta>_coverage_selection.png
```
**SHOW the user the candidate maps and WAIT for confirmation before densifying** (the convention this project uses).
Optional tweaks if asked: add a southern tail, lower floor, k=3 per cell.
⚠ **The SNR floor is SENSOR-DEPENDENT — do NOT hard-code 10.** Template SNR = peak/pre-pulse-RMS; borehole (PB)
sensors record less-impulsive LFE templates → compressed SNR (B018: median 3.2, max 11) vs surface (CPW: 90th-pct 19,
max 39). Check the station's own SNR distribution first (`df.snr.describe(percentiles=[.9,.99])`) and pick a floor at
roughly its top ~1% / a comparable pool size (~300). B018 needed SNR≥6, not 10. Families are still real repeaters; low
template SNR does NOT preclude dv/v (the dv/v uses daily-stack coda, not template SNR).

### 2g. Densify (GPU matched filter, full record)
```
PYTHONPATH=src nohup python scripts/densify_gnw_gpu.py --templates-npz data/<sta>_pnsn_families_100km.npz \
  --summary-csv data/<sta>_coverage_selection.summary.csv --min-snr 0 \
  --network <NET> --station <STA> --out-prefix mf_<sta>_ --workers 20 --top-n 100 --max-raw-det 3000000 --despike-mad 8 \
  >> logs/<sta>_densify.log 2>&1 &
```
Per-year chunks, resumes by skipping existing `mf_<sta>_YYYY.csv`. QC drops anomalous days (median count>2000 or
median cc>0.96); top-100 cap/family-day; `--max-raw-det 3000000` glitch-day guard (speeds glitchy analog-era years).
Then concat: `python -c "import pandas as pd,glob; d=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('data/mf_<sta>_[12]*.csv'))],ignore_index=True); d.to_csv('data/mf_<sta>_all.csv',index=False)"`.

### 2h. Daily long-window stacks
```
PYTHONPATH=src python scripts/build_long_window_resp.py --mf-csv data/mf_<sta>_all.csv --network <NET> --station <STA> \
  --no-deconv --workers 16 --despike-mad 8 --out data/long_window_daily_<STA>.npz   # min-det per density (5 sparse / 8 mod / 20 rich)
```
`--no-deconv` = resample_poly + bandpass (clean). 13 s window (−3..+10 s, 520 samp @ 40 Hz), per (family,day). spawn (no COW).
For response removal: drop `--no-deconv`, add `--inv data/<NET>.<STA>.response.xml` (default order = resample→deconv).
REUSE: merge stacks of the SAME resampler with `scripts/merge_stacks.py` instead of re-stacking.
⚠ **`--min-det` (default 20) is STATION-DEPENDENT** — like the SNR floor. A station with SELECTIVE families (few
detections/family-day) gets a sparse/unusable dv/v at 20. After densify, check the distribution
(`mf.groupby(['template',day]).size().median()`); if it's small (e.g. PGC median **6**), drop min-det to ~8 (PGC: 9,010
→ 94,382 meas, 10×; mean cc 0.978→0.897, acceptable). Prolific-family stations (GNW, hundreds/day) keep min-det=20.

### 2i. Per-era coda dv/v + metadata QC (ALWAYS)
```
PYTHONPATH=src python scripts/dvv_coda_perera.py --npz data/long_window_daily_<STA>.npz --window 1.0 4.0 \
  --station <STA> --era-bounds "<rate/sensor-change dates from 2a>" \
  --out-csv data/daily_dvv_<STA>_perera.csv --out-fig figures/smoke_dvv_<STA>_perera.png --workers 24
PYTHONPATH=src python scripts/plot_dvv_metadata.py --dvv-csv data/daily_dvv_<STA>_perera.csv \
  --network <NET> --station <STA> --out figures/smoke_dvv_<STA>_metadata.png
```
**Read the metadata overlay** to decide which eras are usable. Per-era removes STEPS at changes; it CANNOT fix a
sensor that is internally unstable for a whole epoch (e.g. CPW ROCK1 2011–2019 = discard). See [[metadata-overlay-qc]].

### 2j. Cleanup + log
Delete raw waveforms (`rm -rf data/waveforms/<NET>.<STA>`, re-downloadable); keep detections/stacks/dv/v/families.
**Append to `notes/<today>_Notes.md`:** station, coords, sensor, span, family counts, dv/v verdict, any lesson.

---
## 3. MARGIN-SCALE NOTES
- **Station selection:** walk the tremor belt (~lat 40–50N). Networks: UW/CC (PNSN WA/OR), CN (Canada/VI),
  PB (PBO borehole — active replacements), NC/BK (N California). Center each box on the station; the master
  catalog covers the whole margin. Prefer stations near/over the tremor band; pair decommissioned surface sites
  with co-located active borehole (e.g. CPW↔B018) for long records.
- **Parallel multi-station runs:** to process MANY stations at once, use the **Workflow tool** — fan out one
  pipeline per station (`parallel`/`pipeline`), each station independent. Mind: only ONE GPU job at a time (serialize
  stage-2/densify on GPU), and the 32-CPU cap (don't run >~2 CPU-heavy stages concurrently). Downloads + stage-1
  (CPU) can overlap; GPU steps must queue.
- **State / what's done:** a station is DONE if `data/daily_dvv_<STA>_perera.csv` + `figures/smoke_dvv_<STA>_metadata.png`
  exist. Quick audit: `ls data/daily_dvv_*_perera.csv`. Family sets: `data/<sta>_pnsn_families_100km.summary.csv`.

## 4. CRASH RECOVERY
1. Read this doc + `notes/<latest>_Notes.md` + `MEMORY.md`.
2. `ls data/daily_dvv_*_perera.csv` (done stations), `ls data/mf_*_all.csv` (densified), `ls data/waveforms/`
   (downloaded), `ls logs/` (in-flight jobs). `pgrep -af "discover_gpu|densify_gnw_gpu|build_long_window|download_station"`.
3. Every step is resumable (densify skips existing year files; download skips existing days). Re-run the step that
   was in flight. Catalog + families on disk survive; only the REPL/in-flight commands are lost.

## 5. SCRIPT INVENTORY
- `download_station.py` — generic FDSN waveform pull (any NET/STA).
- `discover_nllb_pnsn_driven.py --candidates-only` — stage-1 candidate detection.
- `discover_gpu.py` — GPU stage-2 clustering (writes snr).  `src/tremorferometry/repeater_gpu.py` — GPU all-pairs CC.
- `select_coverage_families.py` — az×distance coverage selection.  `plot_family_map.py` / `plot_candidate_map.py` / `plot_station_box.py` — maps.
- `densify_gnw_gpu.py` — GPU matched-filter densify (generic via --station/--network/--out-prefix).
- `build_long_window_resp.py` — daily stacks (resample_poly; --no-deconv or --deconv).  `merge_stacks.py` — merge stacks w/o re-stacking.
- `dvv_coda_perera.py` (per-era) / `dvv_coda_parallel.py` (all-time) — coda dv/v.  `plot_dvv_metadata.py` — metadata-overlay QC.
- `plot_dvv_perfamily_deseason.py` — per-family overlay (--no-deseason, --zoom-end-year, cc/err QC, gap-break).
