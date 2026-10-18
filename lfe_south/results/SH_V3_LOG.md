# SH_V3 — Shelly (2025) Mendocino reproduction following Text S1 literally (agent, 2026-09-29)

Inputs: `results/SH_V2_LOG.md`, `notes/shelly_sources/Shelly2025_GRL_SI.txt` (Text S1), notes 2026-09-25
entry "09-29 SHELLY SOURCES" (a–g). v2 files untouched; new scripts `scripts/sh_v3_{common,rescore,restack,mix,report}.py`,
`scripts/sh_v3_queue.sh`; outputs `results/sh_v3/`.

## Step A — v2 detections re-scored with the Text S1 final-catalogue rules (`sh_v3_rescore.py`)

The rules ("Integrated Catalog"): ≥ 9×daily MAD of the summed cc; mean cc over available channels ≥ 0.12;
summed cc ≥ 2.75 (the floor S1 gives for detection); detections within 4 s → keep only the highest summed cc.
- The v2 parquets already hold `ccsum`, `meancc` (= sum / live channels at the peak) and `nlive`, so every rule
  is applied directly. Our net_mad = (sum − median)/MAD. Per family-day fits give a median of the sum of
  |med| < 3e-3 against MAD ≈ 0.38, so this equals his ccsum/MAD to within 1 %.
- His catalogue: the minimum inter-event dt is **4.00 s across all families**, with 0 pairs under 4 s, and
  cc_div_thr ≥ 1.13 = 9/8. So `thresh` is his 8×MAD building threshold and the catalogue is ≥ 9×. The main
  declustering variant `S1` therefore suppresses any detection within |dt| < 4 s of a stronger kept one, in
  any family (greedy by ccsum). `S1x` applies the suppression only between different families.
- `v2` rows reproduce the SH_V2 tables exactly (it0 9× source: 32208 ours / 9082 matched), which checks this
  re-implementation.
- family-specific = our family f matched one-to-one to his family f (the v2 metric). any-fam = all our
  detections matched one-to-one (±2 s) to all his events, pooled over families. same-fam share = the fraction of
  any-fam matches whose families agree. Declustering runs on the full file; the fam-mode day filter comes after it.

**Effect of each rule:** mean cc ≥ 0.12 removes **nothing** at ≥ 8× (with ~18 live channels, 8× MAD already
means mean cc > 0.12). The 2.75 floor removes 6 % at 9× (21 % at 8×). Declustering carries almost all of the change.

it0 (v2 templates from his times), held-out test days:

| variant | xMAD | mode | ours | his | recall fam | precision fam | recall any-fam | precision any-fam | same-fam share |
|---|---|---|---|---|---|---|---|---|---|
| v2 | 9 | source | 32208 | 26895 | 0.338 | 0.282 | 0.428 | 0.358 | 0.74 |
| +floor | 9 | source | 30171 | 26895 | 0.329 | 0.293 | 0.416 | 0.371 | 0.75 |
| **S1** | 9 | source | 14744 | 26895 | **0.267** | **0.487** | **0.406** | **0.740** | 0.66 |
| S1x | 9 | source | 15694 | 26895 | 0.269 | 0.461 | 0.410 | 0.702 | 0.66 |
| v2 | 9 | fam | 4427 | 5021 | 0.323 | 0.366 | 0.345 | 0.392 | 0.92 |
| **S1** | 9 | fam | 2237 | 5021 | **0.248** | **0.557** | 0.275 | 0.616 | 0.90 |
| v2 | 8 | source | 74583 | 26895 | 0.456 | 0.165 | 0.587 | 0.212 | 0.71 |
| S1 | 8 | source | 27367 | 26895 | 0.326 | 0.320 | 0.531 | 0.522 | 0.61 |
| v2 | 10 | source | 17237 | 26895 | 0.244 | 0.381 | 0.303 | 0.474 | 0.77 |
| S1 | 10 | source | 9147 | 26895 | 0.203 | 0.598 | 0.295 | 0.867 | 0.69 |
| S1 | 12 | source | 4107 | 26895 | 0.109 | 0.713 | 0.147 | 0.964 | 0.74 |

it1 (v2 iteration 1 templates):

| variant | xMAD | mode | ours | his | recall fam | precision fam | recall any-fam | precision any-fam | same-fam share |
|---|---|---|---|---|---|---|---|---|---|
| v2 | 9 | source | 43664 | 26895 | 0.360 | 0.222 | 0.466 | 0.287 | 0.71 |
| **S1** | 9 | source | 18180 | 26895 | **0.268** | **0.397** | **0.440** | **0.651** | 0.61 |
| v2 | 9 | fam | 5795 | 5021 | 0.359 | 0.311 | 0.388 | 0.336 | 0.91 |
| **S1** | 9 | fam | 2672 | 5021 | 0.260 | 0.488 | 0.289 | 0.544 | 0.90 |
| S1 | 10 | source | 10999 | 26895 | 0.214 | 0.522 | 0.333 | 0.815 | 0.64 |

What the tables show:
- In v2, **half of our 9× detections were the same event detected by another family**. Declustering removes
  54 % of them (32208 → 14744).
- Family-agnostic precision at 9× rises from 0.36 to **0.74**: three quarters of our catalogue entries are
  events that he also has.
- Family-specific recall falls from 0.34 to 0.27. His family label is the one family whose template won the
  4 s window. Ours wins with a different family in ~1/3 of the matched events (same-fam share 0.66 in source
  mode, 0.90 in fam mode).
- Family-agnostic recall of 0.41 means we find about 41 % of his events at all.
- Rate: 38 per test day, against his 70. The shortfall is recall, not the threshold. At 8× we detect as many
  events as he does, but the precision is only 0.52.
- Per-family and full-sweep tables: `results/sh_v3/rescore_v2it{0,1}.csv` and `rescore_v2it{0,1}_fam.csv`.

## Step B — Text S1 templates (queue `scripts/sh_v3_queue.sh`, launched 2026-09-29 10:18Z, pid 3597880)

**What the S offsets were.** The stage-1 offsets tt_s + best_tau are *not* S picks. In `t4_stack_test.py`,
tau is the lag that maximises the energy ratio E(−1..+5 s)/E(−12..−6 s). It is searched over ±30 s
**independently for every family × station** on top of iasp91 tt_S. Within a family, tau should be one constant
(his catalogue-time convention) if every offset were S. Instead, 91 of 270 v2 template offsets deviate > 1.5 s
from their family median and ~25 deviate > 6 s (NC.KIP, KTJ, KCR, KLR, BK.HALS, PB.B934 up to ±35 s). Those are
noise lags, and some are P-like. The median-centred IQR is ±1 s.

**v3 S arrivals.**
- The predicted S is iasp91 tt_S at his Table S1 hypocentre plus the family constant (median stage-1 tau).
- The pick is an automatic onset on the stack of the top-100 events: the horizontal-channel energy envelope,
  peak within ±3 s, onset at 25 % of peak above noise, with the pick_q column.
- The template is 6 s beginning 1 s before that pick.
- Smoke test (60 build days): pick − stage-1 offset median |1.07| s, and B934's −8.4/+29 s stage-1 lags come back
  to a consistent S.

**Recipe as run:**
- **it0:** his published events on each family's own build days, ranked by his ccsum. 100 per channel with valid
  data, drawn from up to 500 candidates. Build days give only 117–~200 candidates per family, so "500" is never
  reached.
- Per-channel shifts ≤ 0.1 s (±2 samples at 20 sps), re-estimated twice. The final stack keeps the 50 of 100
  with the best leave-one-out cc.
- **it1:** MF on the build days with it0, then the same stacking from our detections at ≥ 8×MAD and summed
  cc ≥ 2.75, ranked by summed cc. The it0 S picks are kept.
- **mix:** the 6 families with the lowest build-day stack quality (median cc50) take it1, the other 21 take it0.
  This is S1's "second iteration for six low-SNR families".
- Detector: t11_mf_v2.py unchanged, with v2 cleaning (dead floors, burst muting, line notching, gap guard).
- Catalogue: 9× plus the step-A rules, via `sh_v3_rescore.py`. Scored on the held-out test days exactly as in
  v2, in source and fam modes, with family-specific and family-agnostic matching.

**Stations (S1 f): 13 of 17 in cache.** BK.BLAS BRIC HALS JEWT PRDS SPRL RVIT, NC.KMR KPR KBU, PB.B049 B934 B933.
- **Missing:** BK.GHOP (FDSN from 2022-05), NC.KMPB (HH from 2003), PB.B932 (EH from 2008), BK.HUNT (from 2022-06).
  All four are at NCEDC/IRIS but were never downloaded.
- Several cached stations cover only part of the build/test days: BLAS 96/406 build days, HALS 116, JEWT 138,
  SPRL 116, KPR 151, PRDS 208, RVIT 318.
- v2's non-S1 stations are dropped: USAL, KBN, KBSB, KCR, KIP, KLR, KSM, KTJ, B047, B935.

**Deviations from S1 that remain:**
1. The S picks are automatic, not manual.
2. The cache is 20 sps, with a 1–8 Hz cache filter plus a 2 Hz high-pass, so the band is effectively 2–8 Hz
   (his is 100 Hz, 2–9 Hz).
3. Four S1 stations are missing, and there are station-day gaps.
4. The initial templates are stacks of his final catalogue on build days, not his manual windows.
5. Scoring is on held-out test days only, not 2018–2024-11 continuously.
6. The 6 low-SNR families are chosen by our stack quality; his six are not named.
7. The it1 stack reuses the it0 S picks; he re-picked on each stack.

Queue ETA ~3 h: tpl0 ~30 min CPU, then MF ~7 s/day × (385 test + 406 build + 385 test), restack it1, rescore
and report. Results are appended below automatically.

<!-- SH_V3_RESULTS_BEGIN -->
## SH_V3 results (auto, sh_v3_report.py, 2026-09-29 12:36Z)

### templates it0

pairs 351; ok 351 over 27 families / 13 stations; status {'ok': 351}; median n_stack 100; median cc50 0.416
S pick: pick_dt (pick - predicted) median 0.55 s, IQR 0.00..1.10; pick_q median 29.2 (q<3: 7); pick - stage-1 offset: median 0.55 s, |.|>1.5 s 67/218, |.|>5 s 6
per station: BK.BLAS n=27 cc50=0.32 q=9.6; BK.BRIC n=27 cc50=0.49 q=54.9; BK.HALS n=27 cc50=0.33 q=9.4; BK.JEWT n=27 cc50=0.46 q=27.9; BK.PRDS n=27 cc50=0.43 q=39.9; BK.RVIT n=27 cc50=0.39 q=16.0; BK.SPRL n=27 cc50=0.51 q=28.9; NC.KBU n=27 cc50=0.53 q=85.1; NC.KMR n=27 cc50=0.43 q=34.5; NC.KPR n=27 cc50=0.41 q=32.4; PB.B049 n=27 cc50=0.45 q=83.2; PB.B933 n=27 cc50=0.32 q=20.4; PB.B934 n=27 cc50=0.34 q=14.7

### templates it1

pairs 351; ok 351 over 27 families / 13 stations; status {'ok': 351}; median n_stack 100; median cc50 0.452
per station: BK.BLAS n=27 cc50=0.38; BK.BRIC n=27 cc50=0.52; BK.HALS n=27 cc50=0.37; BK.JEWT n=27 cc50=0.53; BK.PRDS n=27 cc50=0.50; BK.RVIT n=27 cc50=0.42; BK.SPRL n=27 cc50=0.55; NC.KBU n=27 cc50=0.55; NC.KMR n=27 cc50=0.47; NC.KPR n=27 cc50=0.50; PB.B049 n=27 cc50=0.46; PB.B933 n=27 cc50=0.35; PB.B934 n=27 cc50=0.36

mix: it1 for the 6 lowest-cc50 families [11, 27, 9, 14, 23, 22] (build-day stack quality), it0 for the rest

### it0: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 103468 | 26895 | 0.572 | 0.149 | 0.684 | 0.178 | 0.74 | 4 |
| v2 | 8 | fam | 14400 | 5021 | 0.596 | 0.208 | 0.625 | 0.218 | 0.92 | 10 |
| v2 | 9 | source | 50745 | 26895 | 0.440 | 0.233 | 0.526 | 0.279 | 0.77 | 9 |
| v2 | 9 | fam | 7688 | 5021 | 0.464 | 0.303 | 0.486 | 0.317 | 0.93 | 14 |
| v2 | 10 | source | 29453 | 26895 | 0.330 | 0.302 | 0.391 | 0.357 | 0.79 | 12 |
| v2 | 10 | fam | 4579 | 5021 | 0.349 | 0.383 | 0.362 | 0.397 | 0.95 | 20 |
| v2 | 12 | source | 12301 | 26895 | 0.182 | 0.397 | 0.210 | 0.460 | 0.82 | 17 |
| v2 | 12 | fam | 2000 | 5021 | 0.192 | 0.481 | 0.198 | 0.497 | 0.95 | 23 |
| +floor | 8 | source | 94102 | 26895 | 0.559 | 0.160 | 0.668 | 0.191 | 0.74 | 5 |
| +floor | 8 | fam | 12918 | 5021 | 0.588 | 0.229 | 0.617 | 0.240 | 0.92 | 10 |
| +floor | 9 | source | 49613 | 26895 | 0.436 | 0.236 | 0.521 | 0.282 | 0.77 | 9 |
| +floor | 9 | fam | 7471 | 5021 | 0.461 | 0.310 | 0.483 | 0.325 | 0.93 | 14 |
| +floor | 10 | source | 29359 | 26895 | 0.330 | 0.302 | 0.390 | 0.357 | 0.79 | 12 |
| +floor | 10 | fam | 4570 | 5021 | 0.349 | 0.383 | 0.362 | 0.398 | 0.95 | 20 |
| +floor | 12 | source | 12301 | 26895 | 0.182 | 0.397 | 0.210 | 0.460 | 0.82 | 17 |
| +floor | 12 | fam | 2000 | 5021 | 0.192 | 0.481 | 0.198 | 0.497 | 0.95 | 23 |
| S1 | 8 | source | 33357 | 26895 | 0.410 | 0.331 | 0.642 | 0.518 | 0.64 | 13 |
| S1 | 8 | fam | 4596 | 5021 | 0.404 | 0.441 | 0.439 | 0.479 | 0.92 | 22 |
| S1 | 9 | source | 18058 | 26895 | 0.336 | 0.501 | 0.506 | 0.753 | 0.66 | 21 |
| S1 | 9 | fam | 2860 | 5021 | 0.334 | 0.586 | 0.359 | 0.630 | 0.93 | 25 |
| S1 | 10 | source | 11643 | 26895 | 0.263 | 0.608 | 0.381 | 0.881 | 0.69 | 25 |
| S1 | 10 | fam | 1936 | 5021 | 0.263 | 0.683 | 0.281 | 0.728 | 0.94 | 25 |
| S1 | 12 | source | 5786 | 26895 | 0.151 | 0.703 | 0.207 | 0.964 | 0.73 | 27 |
| S1 | 12 | fam | 1008 | 5021 | 0.153 | 0.763 | 0.162 | 0.807 | 0.95 | 27 |

separation at his 26895 test-day times (best within +-2 s): p50 8.56; fraction >= 8/9/10: 0.576/0.443/0.331

### it1: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 122987 | 26895 | 0.587 | 0.128 | 0.703 | 0.154 | 0.71 | 4 |
| v2 | 8 | fam | 16147 | 5021 | 0.615 | 0.191 | 0.649 | 0.202 | 0.91 | 7 |
| v2 | 9 | source | 63875 | 26895 | 0.457 | 0.193 | 0.552 | 0.232 | 0.73 | 7 |
| v2 | 9 | fam | 9218 | 5021 | 0.481 | 0.262 | 0.506 | 0.276 | 0.92 | 11 |
| v2 | 10 | source | 38515 | 26895 | 0.346 | 0.242 | 0.419 | 0.293 | 0.75 | 10 |
| v2 | 10 | fam | 5767 | 5021 | 0.368 | 0.321 | 0.384 | 0.334 | 0.94 | 18 |
| v2 | 12 | source | 17050 | 26895 | 0.197 | 0.310 | 0.237 | 0.374 | 0.77 | 14 |
| v2 | 12 | fam | 2683 | 5021 | 0.213 | 0.399 | 0.221 | 0.414 | 0.95 | 21 |
| +floor | 8 | source | 114578 | 26895 | 0.574 | 0.135 | 0.688 | 0.162 | 0.71 | 4 |
| +floor | 8 | fam | 15402 | 5021 | 0.608 | 0.198 | 0.642 | 0.209 | 0.91 | 7 |
| +floor | 9 | source | 63062 | 26895 | 0.454 | 0.193 | 0.548 | 0.234 | 0.73 | 7 |
| +floor | 9 | fam | 9131 | 5021 | 0.478 | 0.263 | 0.504 | 0.277 | 0.92 | 11 |
| +floor | 10 | source | 38440 | 26895 | 0.346 | 0.242 | 0.419 | 0.293 | 0.74 | 10 |
| +floor | 10 | fam | 5760 | 5021 | 0.368 | 0.321 | 0.384 | 0.334 | 0.94 | 17 |
| +floor | 12 | source | 17050 | 26895 | 0.197 | 0.310 | 0.237 | 0.374 | 0.77 | 14 |
| +floor | 12 | fam | 2683 | 5021 | 0.213 | 0.399 | 0.221 | 0.414 | 0.95 | 21 |
| S1 | 8 | source | 34859 | 26895 | 0.392 | 0.302 | 0.660 | 0.509 | 0.59 | 10 |
| S1 | 8 | fam | 4775 | 5021 | 0.386 | 0.406 | 0.426 | 0.448 | 0.91 | 20 |
| S1 | 9 | source | 19449 | 26895 | 0.324 | 0.448 | 0.531 | 0.734 | 0.61 | 18 |
| S1 | 9 | fam | 2992 | 5021 | 0.320 | 0.537 | 0.352 | 0.590 | 0.91 | 24 |
| S1 | 10 | source | 12781 | 26895 | 0.255 | 0.537 | 0.408 | 0.859 | 0.62 | 23 |
| S1 | 10 | fam | 2032 | 5021 | 0.254 | 0.627 | 0.274 | 0.678 | 0.93 | 26 |
| S1 | 12 | source | 6610 | 26895 | 0.150 | 0.612 | 0.233 | 0.950 | 0.64 | 24 |
| S1 | 12 | fam | 1110 | 5021 | 0.155 | 0.700 | 0.166 | 0.750 | 0.93 | 26 |

separation at his 26895 test-day times (best within +-2 s): p50 8.69; fraction >= 8/9/10: 0.592/0.460/0.348

### mix: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 105613 | 26895 | 0.576 | 0.147 | 0.689 | 0.175 | 0.74 | 4 |
| v2 | 8 | fam | 13974 | 5021 | 0.596 | 0.214 | 0.627 | 0.225 | 0.92 | 9 |
| v2 | 9 | source | 52733 | 26895 | 0.443 | 0.226 | 0.531 | 0.271 | 0.76 | 8 |
| v2 | 9 | fam | 7758 | 5021 | 0.462 | 0.299 | 0.484 | 0.313 | 0.93 | 14 |
| v2 | 10 | source | 30982 | 26895 | 0.334 | 0.290 | 0.398 | 0.345 | 0.78 | 11 |
| v2 | 10 | fam | 4736 | 5021 | 0.353 | 0.374 | 0.367 | 0.389 | 0.95 | 20 |
| v2 | 12 | source | 13028 | 26895 | 0.183 | 0.378 | 0.214 | 0.442 | 0.82 | 17 |
| v2 | 12 | fam | 2091 | 5021 | 0.194 | 0.465 | 0.201 | 0.483 | 0.95 | 23 |
| +floor | 8 | source | 97826 | 26895 | 0.563 | 0.155 | 0.673 | 0.185 | 0.74 | 5 |
| +floor | 8 | fam | 13263 | 5021 | 0.589 | 0.223 | 0.619 | 0.234 | 0.92 | 9 |
| +floor | 9 | source | 51912 | 26895 | 0.439 | 0.228 | 0.527 | 0.273 | 0.76 | 8 |
| +floor | 9 | fam | 7657 | 5021 | 0.459 | 0.301 | 0.481 | 0.315 | 0.93 | 14 |
| +floor | 10 | source | 30913 | 26895 | 0.333 | 0.290 | 0.397 | 0.346 | 0.78 | 11 |
| +floor | 10 | fam | 4731 | 5021 | 0.353 | 0.375 | 0.367 | 0.389 | 0.95 | 20 |
| +floor | 12 | source | 13028 | 26895 | 0.183 | 0.378 | 0.214 | 0.442 | 0.82 | 17 |
| +floor | 12 | fam | 2091 | 5021 | 0.194 | 0.465 | 0.201 | 0.483 | 0.95 | 23 |
| S1 | 8 | source | 33484 | 26895 | 0.405 | 0.325 | 0.647 | 0.520 | 0.63 | 13 |
| S1 | 8 | fam | 4558 | 5021 | 0.392 | 0.432 | 0.431 | 0.475 | 0.91 | 22 |
| S1 | 9 | source | 18285 | 26895 | 0.333 | 0.489 | 0.511 | 0.752 | 0.65 | 21 |
| S1 | 9 | fam | 2812 | 5021 | 0.323 | 0.577 | 0.351 | 0.626 | 0.92 | 24 |
| S1 | 10 | source | 11881 | 26895 | 0.262 | 0.594 | 0.389 | 0.880 | 0.68 | 25 |
| S1 | 10 | fam | 1942 | 5021 | 0.260 | 0.671 | 0.278 | 0.719 | 0.93 | 25 |
| S1 | 12 | source | 5904 | 26895 | 0.150 | 0.682 | 0.211 | 0.962 | 0.71 | 26 |
| S1 | 12 | fam | 1014 | 5021 | 0.150 | 0.744 | 0.161 | 0.795 | 0.94 | 26 |
<!-- SH_V3_RESULTS_END -->
