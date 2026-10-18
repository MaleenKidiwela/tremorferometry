# lfe_south — dated ledger

## 2026-09-23 — agent T, set-up and launch
* Verified the southern truth on disk (results/SOUTH_TRUTH_INVENTORY.md). Confirmed again that
  `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` is a 1,205-byte saved 404 page; the working
  Shelly copy is `lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv`.
* Ducellier's published parameters read from his own code, not from the paper:
  `src/get_data.py` (detrend → hann taper 5 s → remove_response VEL, pre_filt (0.2,0.5,10,15),
  water_level 80 → bandpass 1.5–9 zerophase → interpolate(lanczos a=10) → decimate to 20 Hz) and
  `src/find_all_LFEs_parallel.py` __main__ (TDUR 10 s, 8×MAD, freq0 1.0 → 1-s decluster,
  12-hour blocks, per-channel normalised cc averaged over channels with a complete gapless block).
* Our `ncc` was verified against his own `correlate.simple` on random data: max |difference|
  6.0e-16, argmax at the injected lag. GPU (cupy 14.1.1, H100) available.
* **Z-only census (the fidelity trap):** 37 of the 57 rows of his `stations_permanent.txt` (65 %)
  are single-component verticals. His detector correlates PER CHANNEL, so they are full members of
  his network. Two consequences, both recorded as fidelity edits, not method changes:
  (a) `t3_cache.py` / `t6_duc_cache.py` cache 1-C and 2-C station-days and use `resample_poly`
      with the exact rational fs/20 (`chain=rational` flagged) where integer decimation does not
      apply; (b) the ported discovery recipe `ks_common.load_day` accepts <3 components
      (the north copy rejected them after a 1-C station-day crashed the 2010 autocorrelation).
* Station geometry derived from FDSN, not assumed: 30,501 channel epochs (EarthScope + NCEDC,
  NC/BK/PB/UW/CC/US/UO, 38.6–42.4 N). 204 stations alive 2026-06-01; **59 alive in both 2009 and
  2026** (50 of them 3-C today) — so the southern bed can be reproduced AND extended on the same
  instruments, unlike POLARIS in the north.
* Stage-1 smoke test (Shelly family 23, 4 stations, 3 days, 121 events): the stack is real — the
  envelope peaks at tau = −5.65 s with peak/median 2.28 and the two nearest 3-C stations agree on
  tau (RVIT −7.0 s, KMR −6.7 s). The dB excess statistic is weak at 121 events (max z 1.97), so
  the per-family detection target was raised from 150 to 400 events (download 4,456 → 12,606
  station-days) and the S window set to Merlin's exact (−1, +5) s against the (−12, −6) s noise
  window. Nothing was filtered on the smoke result.
* **Rank-don't-discard** (user, mid-task): stage 1 reports measured columns (excess dB, percentile
  against the 200 nulls, z, best tau, stations tested and their distances) and a ranked
  distribution; the "≥2 stations above null p99" and "≥10 of 27" bars are printed as reference
  lines. Carry-forward = any station above its own null median.
* **Reproduction reframe** (user, mid-task): the primary scoring axis is now event-by-event
  reproduction of Ducellier and Shelly, before any extension claim. Stages 2 and 4 added.
* Download volumes (recorded as they complete in results/download_manifest.csv,
  duc_download_manifest.csv and ext_download_manifest.csv): stage 1 plan 12,606 station-days
  (~440 GB estimated); Ducellier reproduction 2,800 station-days over 2008-04-01..2008-05-10
  (~70 GB); extension month-by-month thereafter.

## 2026-09-23 01:23 UTC — queue launched
`nohup setsid -f bash scripts/t_south_queue.sh` — **pid 1412422**, detached (PPID 1), flock
`locks/t_south.lock`, log `logs/t_south_queue.log`, markers
`results/T_SOUTH_REPRO_A_DONE | T_SOUTH_STAGE1_DONE | T_SOUTH_REPRO_DONE | T_SOUTH_DONE | T_SOUTH_FAILED`.
Stage order: Ducellier reproduction (bed A then bed B) → existence test → templates → Shelly
reproduction → blind discovery → extension → report. Each stage guarded by its own output file.

**End-to-end smoke of the reproduction chain, before launch** (3 of his 45 stations — BK.YBH,
NC.LRB SHZ vertical-only, PB.B039 — on 2008-04-25 alone, 2 blocks): 184 detections, **59 matched
one-to-one at ±2 s against his own catalogue — pooled recall 0.295, precision 0.321**, 9 of the
36 families active that day with ≥1 match, one family (080429.15.005) at recall 0.974 on 39 of his
events. With 3/45 stations and 1/40 days. The chain reproduces his catalogue.

Channel-selection fidelity fix found by that smoke: BK.YBH carries BH (40 Hz) and HH (80 Hz)
simultaneously, and a single decimation ratio for the whole stream is wrong. `t6_duc_cache.py` now
(a) selects exactly the channels HIS `stations_permanent.txt` names for each station (his templates
are keyed `<STA>_<CHA>`), and (b) applies his `Fraction(int(fs), 20)` interpolate/decimate per
native rate, flagging `chain=frozen|rational` per station-day (NC.LRB SHZ at 50 Hz → rational 5/2,
exactly his formula).
- 2026-09-23 ~01:45Z ORCHESTRATOR: agent T launched t_south_queue.sh pid 1412422. Southern truth verified on disk: Ducellier 66 families 39.45-41.72N, 15-46.8 km, 477,715 permanent-network detections 2004-2011 + 107,924 FAME, WITH templates/station lists/all parameters (1.5-9 Hz, 20 Hz, response-removed VEL, per-channel cc, 8xMAD, 1-s decluster, 12-h blocks, read from his code not the paper). Shelly 2025: 27 families in a 7x15 km zone 40.03-40.09N, 22-29 km, 61,441 detections 2018-2024, times+locations only. catalogs/mtj_ducellier is a 404 page, not data.
  Geometry: 59 stations alive in BOTH 2009 and 2026 (50 now 3-C) -> the southern bed can be reproduced AND extended on the same instruments, unlike the north where POLARIS is gone.
  Z-ONLY TRAP CONFIRMED AND FIXED: 37 of his 57 station rows (65%) are single-component verticals; caching, matched filter and the ported discovery loader now accept 1-C/2-C by construction (fidelity edit, logged).
  SMOKE: our normalised cc matches his correlate.simple to 6e-16; his templates on 3 of 45 stations for 1 of 40 days gave 184 detections, 59 matched one-to-one at +-2 s (recall 0.295, precision 0.321), one family at recall 0.974 -> the reproduction chain works.
  Downloads: 2,800 station-days (~70 GB) for the reproduction window 2008-04-01..05-10 (64/66 families with >=20 events) in flight; existence test 12,606 station-days (~440 GB) next.
- 2026-09-23 ~16:00Z **DUCELLIER REPRODUCED — the strongest validation in the project.** results/duc_repro_bedA_score.csv, 66 families, HIS exact station set, 2008-04-01..05-10, his templates through our frozen chain, scored one-to-one at +-2 s against his published detections: **mean recall_all 0.836 (median 0.878), recall_indep 0.812, recall_loud 0.988 (median 1.000), precision 0.708, and 63 of 66 families pass Bostock's >=30% criterion.** Bed B (his stations PLUS extra modern ones) is markedly WORSE: recall 0.621, precision 0.400, 43/66 passing — adding stations HURT. That inversion needs its own paragraph in SOUTH_RESULTS.md (likely uncorrelated noise added to the network sum shifting the 8xMAD bar) and it is a caution for every "more stations is better" assumption in this project (cf. the north, where 11->15 stations helped).
  Consequence: our chain reproduces an INDEPENDENT catalogue in a DIFFERENT region built by DIFFERENT authors. Washington conclusions drawn with the same chain can now be trusted on their own terms.
- 2026-09-23 ~16:00Z **BOTH QUEUES DIED OVERNIGHT; ~11 h and ~9 h of dead time.** south t_south_queue failed 04:46Z at stage tpl_build (ValueError broadcast (6,240) vs (3,240) — 1-C/2-C stations meeting 3-C template code; 37 of Ducellier's 57 rows are vertical-only); WA u_locate_queue failed 06:18Z at stage north_eval2 after writing loc_ZZdd.csv. Agents V and W launched to fix, resume (NOT redo completed stages) and relaunch. LESSON: a detached queue that dies leaves no notification — the orchestrator must keep a live marker/pid watch armed at all times, not only when a report is expected.

## 2026-09-23 16:11 UTC — agent V: tpl_build crash diagnosed, fixed, queue resumed (pid 1509294)
* **True cause of the 04:46Z crash** (`ValueError: operands could not be broadcast together with
  shapes (6,240) (3,240)` in `t9_templates.py:60`): NOT the 1-C/2-C stations. It is a station whose
  CHANNEL COUNT CHANGES BETWEEN DAYS. `t3_cache.py` keeps every channel of the preferred band, so
  **BK.YBH** — which carries two co-located sensors in HH (HH1/HH2/HH3 plus HHE/HHN/HHZ) — is cached
  6-C on 405 days and 3-C on 13. `t9_templates.py` accumulated per-day stacks without checking the
  count, so the first 3-C day hit the 6-C stack. A header-only scan of all 122 cached station dirs
  found exactly **two** such stations: BK.YBH (6/3) and **NC.LTC (1-C on 104 days, 2-C on 14)**,
  which would have crashed the same stage a few stations later. Stage-1 (`t4_stack_test.py`) never
  crashed because it already takes the station's MODAL count and skips off-count days.
* **Fix** (`scripts/t9_templates.py`): the per-station channel count is now taken from the count
  stage 1 actually stacked tau on (`ncomp` in `stage1_existence.csv`; constant per station there —
  30 stations 1-C, 22 3-C, 1 6-C), with a modal-count fallback; off-count days are skipped and
  counted; and a real guard raises a message naming the station, the day, the two counts and the
  cached channel names instead of a broadcast error. Meta now carries `ncomp`, `ncomp_station`,
  `n_days`, `n_days_offcount` per template. `t11_mf.py` no longer silently drops a station whose
  cached count differs from its template — it counts and logs it.
* **Second trap found before it fired**: `results/stage2_days_test.json` (the held-out days the
  Shelly reproduction is scored on) was never written — `t9` only dumped the half it was building,
  and the `--half all` pass dumps nothing. `t9` now writes BOTH halves whenever it splits; the
  build half is byte-identical to the one already on disk (md5 32a6c953…), test half disjoint,
  1,461 build / 1,418 test family-days.
* **Smoke before relaunch**: templates built for BK.YBH (6-C, 21 templates, 8 off-count days
  skipped), NC.KHBB (3-C, 36) and NC.LGP (1-C, 48) — shapes (6,240)/(3,240)/(1,240); and
  `t11_mf.py` ran on the 1-C templates on a held-out day with GPU=True.
* **Queue relaunched** `nohup setsid -f bash scripts/t_south_queue.sh` — **pid 1509294**, flock
  held, TERM/INT/HUP + EXIT traps (an unexpected exit now writes `results/T_SOUTH_FAILED`), and a
  new `results/T_SOUTH_HEARTBEAT` (stage + pid + UTC) rewritten at every stage so a dead queue is
  visible from disk. It skipped the completed download, bed-A/bed-B reproduction and all stage-1
  work (no station-day was re-fetched) and restarted at `tpl_build`. Remaining: tpl_build, tpl_all,
  Shelly reproduction on held-out days, blind discovery, extension to 2026, report.
* **Bed A vs bed B, verified from `duc_repro_bed{A,B}_score.csv` and now the HEADLINE section of
  `results/SOUTH_RESULTS.md`**: bed A mean recall_all 0.836 (median 0.878), recall_indep 0.812,
  recall_loud 0.988 (median 1.000), precision 0.708, 63/66 at Bostock's 30 % line; pooled recall
  0.873, precision 0.667. Bed B 0.621 / 0.609 / 0.862 / 0.400, 43/66; pooled 0.513 / 0.305.
* **Why bed B is worse — measured, not asserted** (`t13_report.py` recomputes it from the detection
  parquets, so it cannot go stale): channels averaged per family-block median 10 → 21, network
  mean-cc at a detection 0.154 → 0.105, block MAD 0.0150 → 0.0113, detections 10,841 → 13,923.
  On the 4,206 detections both beds share the signal falls ×0.672 while the MAD falls ×0.628
  (independent-noise prediction ×0.615) — so the added channels' noise is independent and they do
  carry LFE energy, about half a his-station channel each. The surplus bed-B detections are NOT
  mis-timed copies of his events (within ±2 s of one of his: 0.679 bed A vs 0.307 bed B; widening
  to ±10 s changes nothing). Reading: averaging ~2× more channels of half the per-channel cc drops
  the averaged SIGNAL faster than the averaged NOISE (~1/sqrt(n)), so the quiet events fall under
  the 8×MAD bar while network-scale noise excursions cross it. The cure to TEST (nothing was
  discarded on this result) is weighting each channel by its own measured cc, not deleting stations.

## 2026-09-25 — agent AA: the Shelly stage was broken in OUR code, and the paper is now read

Shelly et al. (2025) and Sweet, Creager & Houston (2014) are now folded into
`notes/PAPER_METHODS_TABLE.md` (rows **S25** / **SW14** in every subsection, mismatches **M4b**,
**M11**, **M12**) and `notes/CATALOGUE_SUMMARY.md` (Sweet is the new **§4**; Shelly's **§10** is
rewritten from the paper). His actual parameters: **2–9 Hz**, **6 s** candidate window, network
**correlation sum ≥ 9 × the DAILY MAD**, keep only the strongest detection **in each fixed 4 s
window** (that is also his only cross-family dedupe), templates from **2–3 iterations** of
detect → stack the highest-correlation-sum detections → re-pick P and S, locations by grid search in
the Yoon & Shelly (2024) layered P model then **hypoDD**. Station codes, sample rate, response
removal, stack size and the channel-selection rule for the detection sum are **SILENT** (his full
methods are in Supporting Information S1, which is not on disk).

**The 12,848,734 detections were three faults of ours, not a statement about Shelly.**
1. `t11_mf.py`'s `ncc()` normalises by `sd = sqrt(max(M*S2 - S1*S1, 1e-30))` built from cumulative
   sums, while the 20 Hz cache **zero-fills gaps** and stores no validity mask
   (`t3_cache.py:31-32`). In a constant run the true window variance is 0, the cumsum difference
   goes negative from round-off, the clamp gives `sd = 1e-15`, and the "normalised" cc reaches
   **1e10–1e15**. Measured on family 23, 2022-05-11: BK.RVIT 33.3 % constant samples (|cc|max
   2.4e15), NC.KIP 50.1 % (2.8e10); the channel-summed series then has **MAD 25.7** where a
   20-channel normalised sum cannot exceed 20, **16.6 %** of samples clear 8×MAD, and that ONE day
   yields **128,922** detections for 27 families against Shelly's **71**.
2. The network stack was summed **at a common lag** even though each station's template is cut at
   its own offset `tt_s + tau_s` (`t9_templates.py:93`) — spread **1.49–5.85 s** inside family 1, so
   true events never stacked and our times sat on no published axis. After shifting each channel by
   `tt+tau-pre` our peak lands at his published times with **median dt 0.00 s**.
3. Scoring: `t8_score.py` used **all 61,441** of his events as the denominator while the filter ran
   only on the **385 held-out days**, where only **5,021** of them live — the recall **cap was
   0.082**, so "pooled recall 0.026" was 0.32 of the achievable maximum.
Ducellier's path is immune to 1 and 2 (his cache carries `valid` and he requires a gapless block,
`t7_duc_mf.py:113`; his templates are all cut at the same absolute catalogue time on every station),
which is exactly why the same chain reproduced him at recall 0.836 and failed here.

**Fixed, as COPIES** (`t11_mf.py` is in use by the running `t_south_queue.sh` extension stage, pid
1509294, and was not touched): `scripts/t11_mf_aa.py` (live-window mask, moveout shift, `--tmin 9`,
fixed-4-s-window dedupe, `--hp2` for 2–9 Hz, per-family-day diagnostics + threshold sweep + the
statistic at every published time), `scripts/t8_score_aa.py` (`--days-mode source|fam|all`),
`scripts/t9_templates_aa.py` (6 s templates, writes its day-split JSON under `results/aa_shelly/`),
`scripts/aa_restack.py` (Shelly's template iteration on BUILD days only), `scripts/aa_report.py`.
Smoke on 2022-05-11: **128,922 → 21** detections, **11** of them his (precision 0.52, recall 0.155
that day). Full re-run: `scripts/aa_shelly_queue.sh`, outputs in `results/aa_shelly/`, report
`results/aa_shelly/AA_SHELLY_REPORT.md`, markers `results/AA_SHELLY_DONE|FAILED`.

⚠ **Still owed:** `t11_mf.py` itself carries fault 1 and 2, so the 995,201 core / 1,302,991 ext
extension detections in `SOUTH_RESULTS.md` §5 are inflated by the same bug and must be re-run once
the running queue finishes. Nothing was deleted.
- 2026-09-25 ~03:30Z **ORCHESTRATOR: the 12.8 M Shelly detections were THREE bugs in our code, all now measured and fixed in copies (agent AA; queue aa_shelly_queue.sh pid 1800120).**
  (1) **Zero-fill gap bug:** t11_mf.py ncc() divides by sd=sqrt(max(M*S2-S1*S1,1e-30)) while t3_cache.py:31-32 zero-fills gaps with NO valid mask -> constant runs give sd~1e-15 and "normalised" cc of 2.4e15. Family 23 on 2022-05-11: RVIT 33% / KIP 50% constant, summed-cc MAD 25.7 where 20 channels cap at 20, 16.6% of samples over 8xMAD, **128,922 detections that day against Shelly's 71**.
  (2) **Common-lag summation:** cc series were summed at a common lag although per-station template offsets tt+tau span 1.49-5.85 s within one family -> incoherent stack on the wrong time axis.
  (3) **Scoring denominator:** all 61,441 of his events were scored against 385 days holding only 5,021, so the recall CAP was 0.082 - the reported 0.026 was 0.32 of the achievable maximum, not 0.026 of his catalogue.
  **Ducellier's path is immune to (1) and (2)** - valid mask, gapless blocks, and his released templates share one absolute cut time. That contrast is the whole explanation for 0.836 vs 0.026.
  **Fix verified on the same day: 128,922 -> 21 detections, 11 of them his (precision 0.52), his times landing at dt median 0.00 s.** Ablation with normalisation fixed but moveout not: 477 detections, precision 0.023 - BOTH fixes are required.
  **OWED: t11_mf.py still carries both bugs, so the extension counts in SOUTH_RESULTS.md section 5 are inflated and must be re-run after pid 1509294 finishes.**

## 2026-09-30 — DUC_FIX: Ducellier extension re-run with the FIXED detector (queue launched)
* Queue `scripts/duc_ext_fixed_queue.sh` (pid in `results/DUC_FIX_PID`, flock `locks/duc_fix.lock`, markers
  `results/DUC_FIX_{DONE,FAILED,HEARTBEAT}`, log `logs/duc_ext_fixed_queue.log` + `logs/duc_fix_{fetch,cache,prep}.log`).
  New files only: `duc_fix_common.py`, `duc_epochs.py`, `duc_fix_data.py`, `duc_fix_clean.py` (copy of sh_v2_common for
  cache_duc: valid mask kept, no 2 Hz high-pass, line band 1.5-9 Hz), `duc_fix_prep.py`, `duc_fix_mf.py`,
  `duc_fix_stage1.py`, `duc_fix_retemplate.py`, `duc_fix_report.py`. Outputs `results/duc_fix/`, report
  `results/DUC_EXT_FIXED_RESULTS.md`.
* His released code read again: the MAD is HOURLY (`find_all_LFEs_parallel.py` loops hours; his 12-h is job chunking);
  channel enters an hour only when complete; his catalogue rule `cc*nchannel >= threshold_perm` (threshold_cc.txt,
  his analysis scripts); threshold 0 = his six eliminated families. His loc '0' = '00'. Both now implemented.
* Station epochs (FDSN StationXML, key = loc/channels/sr/sensor/azimuths): 63 'his' epochs, 43 'new'. The big one:
  2011-10-27..12-08 every NC analogue SHZ@50 Hz -> digital EHZ@100 Hz (same L-4C), then loc-code moves 2013-2024,
  sensor swaps KRP/KHBB -> Trillium Horizon (2021/22), YBH/GASB/HOPS -> STS-7 (2022-24), GBB -> MBB-2, LTC -> Trillium
  Compact. Sensor swaps INSIDE 2008-2011 (KHBB/KHMB/KRP/KSXB/KRMB 40T->3ESPC/STS-2) keep his template (he ran them).
  Azimuth changes on his epochs (KRMB 2008-11..2009-09, KRP N=180 2010-11..2015-10) rotated to the template
  orientation as his rotate_data does. New epochs re-templated (his recipe: top-150 by cc, RMS linstack) from pass-1
  detections made WITHOUT the station, kept per family-channel only if the stack beats 20 scrambled controls
  (station data shifted by whole days); else dropped for that epoch. Months touched are re-scanned (pass 2).
* Smoke (4 days 2008-04-24..27, his stations): burst muting HURTS on this sparse network (mean recall 0.789 / pooled
  precision 0.552 with muting vs 0.851 / 0.936 without; old unfixed run 0.750 / 0.558) - a near-station LFE on one
  station with quiet neighbours 20-50 km away reads as a burst. Stage 1 therefore scores both and carries the variant
  that reaches >= 0.836 with the higher pooled F1 under his cc x nchannel rule (`results/duc_fix/detector_choice.json`; raw pooled numbers are dominated by his eliminated family 080328.09.029, 10,031 raw detections in 40 days).
* Speed-ups that do not change a sample (verified max diff 0 / 1e-9 relative): lanczos evaluated directly on the 20 Hz
  grid (identity for integer ratios), evalresp memoised and evaluated only < 15 Hz (his pre_filt is exactly 0 above),
  remove_response FFTs on the GPU.
* Old extension outputs marked SUPERSEDED (`results/ext/SUPERSEDED_Ducellier2022.md`, SOUTH_RESULTS.md §5 note); not deleted.
