---
name: glitch-days-false-dvv-drop
description: Glitchy/fragmented station-years cause FALSE dv/v drops; fix with --despike-mad in densify + stacking
metadata: 
  node_type: memory
  type: reference
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**Some station-years have impulsive-glitch + fragmented-mseed contamination (telemetry trouble) that produces a FALSE dv/v drop.** Found 2026-06-05 at **PGC 2013–2014** (the "sudden drop 2013–2015" the user flagged).

**Mechanism:** impulsive glitches → matched filter false-triggers thousands of times/day → median count/fam-day exceeds the densify QC (`--qc-median-count 2000`) → ~90% of days dropped → only ~40 glitchy survivors/yr → their biased stacks read a deep low (PGC: −0.074%). It LOOKS like a sudden velocity drop but is a data artifact. It survives all reference schemes (per-era, all-time) because they share the corrupted data — so "survives both references" does NOT prove a feature is real.

**How to detect (add to QC):**
- traces-per-day: `len(obspy.read(dayfile))` — clean=1, fragmented=50+ (PGC 2013/2014 ≈ 55).
- days-with-detections per year in `mf_<sta>_YYYY.csv`: healthy ≈ 350; if ~40 with high median/fam-day (58–63 vs normal 6–8), suspect glitch-QC-drop.
- outlier-sample fraction (>8×MAD after bandpass): PGC 2013 had 5–100× more than clean years.
- NOTE: data is usually NOT actually missing — only ~0.4% real gaps; it's glitch-contaminated, not absent.

**Fix (implemented):** both `scripts/densify_gnw_gpu.py` and `scripts/build_long_window_resp.py` now take **`--despike-mad K`** (winsorize/clip bandpassed samples beyond K·MAD/day; 0=off, default off so other stations unaffected). Use **`--despike-mad 8`**, and **set it on BOTH densify and stacking** so detections and coda stacks agree. Result at PGC: 2013 42→365 days, 2014 38→364 days, 0 QC-dropped; false low −0.074% → mild −0.020% continuous with neighbors (`smoke_dvv_PGC_alltime_1to3_despike.png`). Backups in `data/backup_pgc_glitchy/`.

**Scan of all done stations (2026-06-05):** glitch problem was ISOLATED to PGC 2013-2014 (fixed). GNW/HDW/B018/NLLB clean. GNW 2026 flagged but FALSE POSITIVE (partial current year + GNW's normal high density ~72/fam + 1 trace/day). CPW 2012-2013 flagged but = the already-known **ROCK1 bad-sensor era 2011-2019** (broad elevated density 44-83/fam across 2012-2017, not a sharp 2-yr glitch) → degraded sensor, not despike-fixable, already marked discard; CPW waveforms deleted so re-download needed to test. Scan rule (post-2010 yr, <150 days-with-det AND >30 det/fam-day) has false positives → always confirm vs the station's OWN density baseline + traces/day before concluding glitch.

**Still TODO:** regen PGC canonical per-era product with despike; auto-flag fragmentation in the pipeline QC. See [[margin-station-status]], [[metadata-overlay-qc]].
