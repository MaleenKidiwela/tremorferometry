# lfe_oregon — work log

## 2026-09-28 — setup, papers, cache, cross-station detector, chain

### Papers
- Savard & Bostock is **2015** BSSA 105(4):2128 (T&B's "2014" citation is the AGU abstract). The paper is closed access; the user supplied the PDF → `notes/oregon_sources/`.
- Found Savard's original MATLAB on GitHub → cloned read-only to `notes/oregon_sources/CS-quake-detection/`. It has an Oregon branch (`region='OR'`, `Dsta_OR.mat`: 42 stations including ALP0 and BWO0 but **not ALS0**). Distances check out: ALP0–BWO0 is 32.99 km in Dsta_OR vs 32.96 km from our coordinates (UTM 10).
- Rubin & Armbruster 2013: not obtained (Wiley 403 to curl). The recipe doesn't need it.
- Recipe and every silent point: `METHOD_FROM_PAPER.md`.

### Cache (`scripts/or_cache40.py`, 12 workers, OMP 1)
- 1–10 Hz, 40 sps, validity mask (copied from ab_cache40 logic). 1299 station-days built, 81 noraw (FACM 46, FACK 29, HSO 4, B027 1, FACO 1). The median frac_valid is 0.99997.
- Forced deviations:
  - IU.COR is 20 sps BH, so its band is 1–9 Hz before upsampling.
  - PB boreholes are rotated with StationXML azimuths (`meta/pb_azimuths.json`); B031 Z (dip +90) is sign-flipped.
  - EHZ-only UW stations (BRO, HSO, MPO, RNO, SSO) store N/E as zeros.
- FAILURE → fix: XA.BWO0 2009-09-09 and 09-10 interleave **40 Hz packets with the 100 Hz HH stream** under one channel code ("mixed sampling rates"). Fix: keep the rate with the most samples; the dropped packets become mask gaps. Rebuilt OK.
- Low-coverage station-days (min frac_valid): B027 0.002, FACO 0.013, MAP0 0.17, YAC1 0.21, EUO 0.21, FACK 0.21 (per-day minima; see results/cache40_manifest.csv).

### Cross-station detector (`scripts/xs_detect.py`)
- A literal port of main_scc/cleanwf/lfedetect/mccc2/ccsearch/gettwin/shift.
- Unit test: three shifted wavelets are aligned exactly by mccc2 + fshift (sign convention OK, σ≈0).
- Lag limits (s): ALP0–ALS0 ±2.53, ALS0–BWO0 ±2.85, ALP0–BWO0 **±13.7** (D = 33 km > 25 km, so 1.5·D/3.6).
- One-hour test (09-04 00 h): 11 s per hour.
- **ALS0 passes cleanwf in only 35–67 % of windows.** Cause: the f-weighted energy above 6 Hz is > 70 % (median 73 % on N). ALS0 has 54 % of its 1–10 Hz power above 6 Hz. The 1–10 Hz band (T&B) combined with a 6 Hz edge (CODE, tuned for 1–8 Hz) removes it. This is a known sensitivity: D9/S4 column. ALP0 and BWO0 pass > 94 %.
- **Station floor:** the code's ≥4 stations is impossible with 3 stations → `--min-sta 3` (S1).

### Cross-station detections, 3 days 09/03–09/05, ALP0/ALS0/BWO0 (`results/xs/`)
| day | windows | pass stage 1 (Cenv) | double-count suppressed | ccsearch fail | **detections** | Q2-lost |
|---|---|---|---|---|---|---|
| 09-03 | 172,800 | 537 | 50 | 458 | **9** | 11 |
| 09-04 | 172,800 | 832 | 130 | 642 | **13** | 34 |
| 09-05 | 172,800 | 2,294 | 142 | 2,048 | **26** | 55 |
| total | | 3,663 | | | **48** | |

- Most failures happen at ccsearch. Quirk Q1 (the full channel set is never accepted) matters a lot with 3 stations: a window whose envelope stage kept exactly 3 channels (one per station) can **never** detect (403 such windows over 3 days).
- Ablation `--allow-full-set` (the paper's "largest subset S3"), `results/xs_fullset/`: **37 / 43 / 109 = 189** detections. This is a column; the literal default stays 48.
- n_ch in the final set (literal run): 3 channels in 19 of 48 detections. With 3 channels, σ has 1 degree of freedom, so the σmax test is weak for these (sigr_max 0.00 happens).
- Sanity checks:
  - stage-1 envelope cc stays within [0, 0.91] (lower bound 0 because MATLAB takes the max of mask·ccf).
  - rmean ≤ 0.50; sigr_max median 0.12–0.19 s (< 0.3).
  - **0 mask violations** (every detecting station's 18 s window is fully valid).
  - At most 4 detections per hour; detections fall in 6 / 11 / 17 distinct hours per day, rising into 09-05.
- Figures: `figures/xs/xs_per_hour.png` and 6 record sections `figures/xs/xs_*.png` (all 29 stations, N/E/Z, S picks). Single-event windows are noisy; ALS0 shows the clearest S burst.
- Runtime: about 3.5 min per day on CPU.

### Chain (written; smoke-tested only)
`scripts/xs_chain_queue.sh <TAG>` runs `xs_grow.py` → `xs_scan46.py` → `xs_locate.py` → `xs_review.py --mode stacks`, with markers `results/<TAG>/CHAIN_DONE|CHAIN_FAILED`. GPU is used through the imported `ab_grow.matched_filter`. The lfe_array_discovery code is imported; its paths are patched in-process to our cache, and none of its files are modified.

Smoke test, 1 day (09-05), 1 iteration, the 26 candidates from that day (`results/smoke_grow`):
- 149 network detections at 8×MAD. The threshold is about 0.15 (mean channel cc, 9 channels). The median is 4 detections per template (range 1–20).
- **All 26 templates re-detect themselves** (cc = 1.0 at |dt| < 0.5 s). The median non-self cc is 0.16.
- FAILURE → fix: the first admission null (99th percentile) admitted 18 stations after one iteration. That is too loose: 29 stations × 15 iterations of tests, and the member half-cc median is only 0.18. Changed to the **99.9th percentile of 2000 random-time half-stack draws** (for incoherent noise the null does not depend on stack size). Result: 1 station admitted (UO.EUO, half-cc 0.20 vs null 0.16). Per-station null99.9 ranges 0.13–0.71; UW.TOLO is highest because of coherent tonal noise. The rule is still SILENT-S11; every station's statistic is a column.
- Downstream plumbing, 2 days, templates with ≥ 5 detections (12): 1165 detections. The 10 %/5 s elimination keeps **1 of 12**: all share 29–49 % of their times with the largest. These 1-iteration, 3-station templates all catch the same tremor bursts. This is an outcome to watch in the full run, not something to tune.
- The Hyp2000 free and fixed-to-McCrory-2012 runs work. With rmin 1.5 the smoke template had no usable picks (that output, results/smoke_hypo, was my own minutes-old smoke output and I removed it; r1 kept). With rmin 1.0 (plumbing test only) it produced a junk location (rms 7 s, erh 10 km), as expected for 31 stacked events.
- The review sheet renders (`figures/smoke_review/fam_000005.png`): ALS0 and BWO0 show a coherent S at +2–3 s.

### Launching the full chain needs
1. Coordinator OK for GPU time. It is light: 48 templates × 3 days × ≤ 15 iterations, then the passing templates × 46 days. Roughly tens of minutes.
2. `nohup bash scripts/xs_chain_queue.sh chain_v1 &` (literal candidates). Optionally `chain_fullset` with `--xs results/xs_fullset`.
3. For location as T&B did it: the **PNSN O0** and **Williams et al. (2011)** 1-D layer tables. Until then GSC and Crosson are placeholders.

## 2026-09-28 22:00Z — Savard-style location of INDIVIDUAL detections (scripts/xs_event_locate.py → results/chain_paper/event_locs/)
Method: every one of the 33,382 det46 detections of the 26 families; per station an 18 s N/E window centred on the template's expected S
(clear template pick rS≥1.5 where there is one, else the 1-D USGS_OR straight-ray S from a family centre, origin anchored to the clear picks;
centre = interface grid fit to ≥3 clear picks, else the seed centroid; pass 2 re-centres on the pass-1 median), band 1–8 Hz (S&B's band;
cache is 1–10), then xs_detect's literal cleanwf → envelope cc > 0.43 → envelope MCCC → ±4 s mask → ccsearch (±1 s, σ<0.3 s, Q1 kept) →
t_S = tabs − tdel − tdel0, weight fix(σ/0.1). Hyp2000 S-only + Savard's bogus P (weight 3, WET 1 .5 .2 0), depth fixed (ZTR z T) at the
McCrory 2012 depth under the family centre; USGS_OR Vp/Vs 1.71 headline, GSC 1.73 column; other Hyp2000 settings default (ERR 0.15, default
DIS/RMS weighting). Null = same code with every detection moved +1800 s (dropping shifted times within 20 s of any real detection; 18,452).
Runtime ~2 min per pass on 12 CPUs; anon memory ≤ 18 GB.

| | real (33,382) | null +1800 s (18,452) |
|---|---|---|
| ≥3 stations after ccsearch | 1,936 (5.8 %) | 634 (3.4 %) |
| exactly 3 | 1,585 | 550 |
| ≥4 stations (paper minimum) | 351 (1.05 %) | 84 (0.46 %) |
| Hyp2000 solution | 145 | 44 |
| **erh < 3 km** | **58 (0.17 %)** | 19 (0.10 %) |
| erh<3 and ≥4 distinct sites (B027/B028, B031/B032 co-located) | 32 (0.096 %) | 8 (0.043 %) |
| erh<3 and ≥2 seed stations (ALP0/ALS0/BWO0) | 30 (0.090 %) | 3 (0.016 %) |
- 351 → 145: Hyp2000 abandons events whose 4th S is zero-weighted by default distance/residual weighting ("ONLY 3 READINGS").
- Scaled to equal size, the null predicts ~34 of the 58 erh<3 real solutions (~60 % noise-consistent); the null cloud sits at the same place
  (median 44.42 N −123.58, real 44.53 N −123.60) → the event epicentres are set by the network geometry around the seed stations.
  The number of seed stations is the column that separates real from null (30 vs 3).
- Per family (erh<3, USGS_OR): 102:6, 3:5, 116:5, 62:4, 117:4, 13:4, 37:3, 27:3, 120:3, eight families with 2, five with 1, four with 0 (68, 50, 58, 52).
  72 % of the located events are co-detections (within 5 s) of another family, i.e. the same physical event located once per family.
- Separation: 9 families have ≥3 events; medians 2–39 km apart, r50 spread 3–20 km each, SE of medians ~5–10 km; **no pair is significantly
  separated (all z_sep < 3; 1 group)**, clouds overlap (d/spread<1: 2 groups). GSC gives the same picture (family_summary_with_gsc.csv).
- T&B pattern (13 located, front north Aug 30, 5 km/day south): not testable, since there are too few locations and all of them fall in 44.45–44.65 N.
- Reading: the individual MF detections are far below single-event visibility. An envelope stack of the top 300 detections of family 68 only has
  peak/median 1.10–1.15 at the seed stations, and the linear stack's SNR (5–8) grows as √N. So Savard's per-event S timing, which works on events
  selected for being visible across stations, gives mostly noise-level locations here. **This test cannot confirm or refute the 26.** It supports
  neither more than 1 resolved source region nor the 8–26 range. Files: family_summary.csv, family_pairs.csv, hyp_pass{1,2}.csv, meas_*.parquet,
  stimes_*.parquet, map.png; null in event_locs_null/.
