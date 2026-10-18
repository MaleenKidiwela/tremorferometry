# SH_V3B — v3 rerun with the full 17-station Text S1 set (agent, 2026-09-29)

**What changed from v3 (results/SH_V3_LOG.md).** The four S1 stations v3 lacked are added: BK.GHOP (HH, from
2022-05-05), NC.KMPB (HH), PB.B932 (EH) and BK.HUNT (HH, from 2022-06-28). Everything else is the v3 recipe and code.
- **Days:** every Shelly2025 build + test day (stage2_days_{build,test}.json, 572 days in 2018-01..2024-11) inside each
  station's FDSN epoch. That is 1844 station-days planned: B932 572, KMPB 572, GHOP 361, HUNT 339.
- **Download:** unchanged `t2_fetch.py`. It takes native-rate BH/HH/EH/SH velocity channels, with NCEDC for BK/NC and IRIS
  for PB, as in the existing cache. A probe on 4 days found IRIS has no GHOP/KMPB waveforms, so NCEDC is the only
  source. Some days are empty at NCEDC for both stations (e.g. 2023-03-01) and are marked `.nd`.
- **Cache:** unchanged `t3_cache.py`: 1–8 Hz Butterworth, resample_poly to 20 sps, HH > BH > EH, same
  `cache20hz/<NET.STA>/<yyyy>/<doy>.npz` layout. The existing cache has **no response removal** (resp/ is used only by
  the Ducellier cache_duc), so none is applied here. Coordinates come from the existing `results/fdsn_channels_south.csv`,
  which already lists all four stations.
- **Prep** (`scripts/sh_v3b_data.py prep` → `results/sh_v3b/prep/`):
  - v2 `prep_day` (dead floors, burst muting, line notching, refs) is run on the v2 prep station set (23) plus the
    4 new stations, into `prep_full/`.
  - The merged day file keeps the v2 prep entries of the old stations **byte-identical**, verified on 2 smoke days,
    and adds the new stations' entries. So v3 → v3b changes only the station set.
  - The one side effect is small: the old stations' burst-mute masks do not see the new stations as neighbours.
- **Code changes, defaults unchanged:**
  - `sh_v2_common.PREP` honours `SH_PREP_DIR`.
  - `sh_v3_common.V3` honours `SH_V3_DIR`.
  - `sh_v3_report.py` honours `SH_V3_LOG` / `SH_V3_NAME`. With `SH_V3_BASE` set, it adds a v3 vs v3b table (S1 rules,
    8/9/10×, source + fam, family-specific + any-family) and per-family recall deltas (S1 9×, mix and it0).
  - `sh_v3_restack.py` already listed all 17 S1 stations and picks up any that are in the cache.
- **Smoke** (2 days, 2022-06-30 / 07-02):
  - fetch: 8/8 station-days, 0.23 GB, 4 s. cache: 8 built (HH/HH/HH/EH, 3 channels each). prep: new stations
    muted 0.3–8 %, 0 dead.
  - restack his-mode: 17/17 S1 stations found and templates built on all 4 new stations. MF ran on the merged prep
    (|cc|>1 = 0). rescore and report comparison block ran.
- **Queue:** `scripts/sh_v3b_queue.sh`, outputs in `results/sh_v3b/`. log `logs/sh_v3b_queue.log`, markers
  `results/SH_V3B_{DONE,FAILED,HEARTBEAT}`, pid `results/SH_V3B_PID`.
- **Deviations:** the v3 list 1–7 still applies, except item 3: now 17/17 stations, though the station-day gaps remain.
  GHOP and HUNT cover only 2022-05+ (≈ 60 % of the days). The mix family choice is recomputed from v3b it0 cc50, so it
  may differ from v3.

### data obtained (auto, 2026-09-29 13:32Z)

| station | planned days | miniSEED days | no-data days | fetch errors | cached (built/cached) | prepped | station-dead days | band |
|---|---|---|---|---|---|---|---|---|
| BK.GHOP | 361 | 332 | 29 | 0 | 332 | 332 | 0 | HH |
| BK.HUNT | 339 | 253 | 86 | 0 | 253 | 253 | 0 | HH |
| NC.KMPB | 572 | 499 | 73 | 0 | 499 | 499 | 1 | HH |
| PB.B932 | 572 | 563 | 9 | 0 | 563 | 563 | 0 | EH |

<!-- SH_V3_RESULTS_BEGIN -->
## SH_V3B results (auto, sh_v3_report.py, 2026-09-29 15:43Z)

### templates it0

pairs 459; ok 459 over 27 families / 17 stations; status {'ok': 459}; median n_stack 100; median cc50 0.409
S pick: pick_dt (pick - predicted) median 0.70 s, IQR 0.05..1.25; pick_q median 28.2 (q<3: 10); pick - stage-1 offset: median 0.55 s, |.|>1.5 s 67/218, |.|>5 s 6
per station: BK.BLAS n=27 cc50=0.32 q=9.6; BK.BRIC n=27 cc50=0.49 q=54.9; BK.GHOP n=27 cc50=0.55 q=86.7; BK.HALS n=27 cc50=0.33 q=9.4; BK.HUNT n=27 cc50=0.38 q=14.5; BK.JEWT n=27 cc50=0.46 q=27.9; BK.PRDS n=27 cc50=0.43 q=39.9; BK.RVIT n=27 cc50=0.39 q=16.0; BK.SPRL n=27 cc50=0.51 q=28.9; NC.KBU n=27 cc50=0.53 q=85.1; NC.KMPB n=27 cc50=0.38 q=29.0; NC.KMR n=27 cc50=0.43 q=34.5; NC.KPR n=27 cc50=0.41 q=32.4; PB.B049 n=27 cc50=0.45 q=83.2; PB.B932 n=27 cc50=0.36 q=14.3; PB.B933 n=27 cc50=0.32 q=20.4; PB.B934 n=27 cc50=0.34 q=14.7

### templates it1

pairs 459; ok 459 over 27 families / 17 stations; status {'ok': 459}; median n_stack 100; median cc50 0.454
per station: BK.BLAS n=27 cc50=0.37; BK.BRIC n=27 cc50=0.51; BK.GHOP n=27 cc50=0.60; BK.HALS n=27 cc50=0.36; BK.HUNT n=27 cc50=0.43; BK.JEWT n=27 cc50=0.52; BK.PRDS n=27 cc50=0.51; BK.RVIT n=27 cc50=0.43; BK.SPRL n=27 cc50=0.55; NC.KBU n=27 cc50=0.57; NC.KMPB n=27 cc50=0.42; NC.KMR n=27 cc50=0.47; NC.KPR n=27 cc50=0.51; PB.B049 n=27 cc50=0.47; PB.B932 n=27 cc50=0.41; PB.B933 n=27 cc50=0.35; PB.B934 n=27 cc50=0.37

mix: it1 for the 6 lowest-cc50 families [3, 11, 27, 9, 14, 22] (build-day stack quality), it0 for the rest

### it0: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 134579 | 26895 | 0.767 | 0.153 | 0.846 | 0.169 | 0.80 | 7 |
| v2 | 8 | fam | 18444 | 5021 | 0.762 | 0.207 | 0.783 | 0.213 | 0.94 | 8 |
| v2 | 9 | source | 77159 | 26895 | 0.625 | 0.218 | 0.704 | 0.245 | 0.81 | 8 |
| v2 | 9 | fam | 11181 | 5021 | 0.624 | 0.280 | 0.644 | 0.289 | 0.95 | 13 |
| v2 | 10 | source | 49057 | 26895 | 0.489 | 0.268 | 0.553 | 0.303 | 0.82 | 11 |
| v2 | 10 | fam | 7332 | 5021 | 0.493 | 0.337 | 0.507 | 0.347 | 0.95 | 16 |
| v2 | 12 | source | 22549 | 26895 | 0.288 | 0.343 | 0.323 | 0.385 | 0.84 | 15 |
| v2 | 12 | fam | 3514 | 5021 | 0.297 | 0.424 | 0.303 | 0.433 | 0.96 | 22 |
| +floor | 8 | source | 125261 | 26895 | 0.755 | 0.162 | 0.836 | 0.180 | 0.80 | 7 |
| +floor | 8 | fam | 17110 | 5021 | 0.748 | 0.220 | 0.769 | 0.226 | 0.94 | 9 |
| +floor | 9 | source | 76855 | 26895 | 0.625 | 0.219 | 0.704 | 0.246 | 0.81 | 8 |
| +floor | 9 | fam | 11141 | 5021 | 0.624 | 0.281 | 0.643 | 0.290 | 0.95 | 13 |
| +floor | 10 | source | 49057 | 26895 | 0.489 | 0.268 | 0.553 | 0.303 | 0.82 | 11 |
| +floor | 10 | fam | 7332 | 5021 | 0.493 | 0.337 | 0.507 | 0.347 | 0.95 | 16 |
| +floor | 12 | source | 22549 | 26895 | 0.288 | 0.343 | 0.323 | 0.385 | 0.84 | 15 |
| +floor | 12 | fam | 3514 | 5021 | 0.297 | 0.424 | 0.303 | 0.433 | 0.96 | 22 |
| S1 | 8 | source | 37973 | 26895 | 0.562 | 0.398 | 0.809 | 0.573 | 0.70 | 15 |
| S1 | 8 | fam | 5361 | 5021 | 0.520 | 0.487 | 0.558 | 0.523 | 0.93 | 23 |
| S1 | 9 | source | 23534 | 26895 | 0.485 | 0.555 | 0.685 | 0.782 | 0.71 | 25 |
| S1 | 9 | fam | 3651 | 5021 | 0.453 | 0.623 | 0.483 | 0.665 | 0.94 | 27 |
| S1 | 10 | source | 16307 | 26895 | 0.392 | 0.647 | 0.541 | 0.892 | 0.73 | 27 |
| S1 | 10 | fam | 2613 | 5021 | 0.369 | 0.709 | 0.392 | 0.753 | 0.94 | 27 |
| S1 | 12 | source | 8879 | 26895 | 0.238 | 0.720 | 0.317 | 0.961 | 0.75 | 27 |
| S1 | 12 | fam | 1490 | 5021 | 0.230 | 0.777 | 0.243 | 0.819 | 0.95 | 27 |

separation at his 26895 test-day times (best within +-2 s): p50 9.92; fraction >= 8/9/10: 0.772/0.628/0.490

### it1: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 160078 | 26895 | 0.781 | 0.131 | 0.858 | 0.144 | 0.78 | 6 |
| v2 | 8 | fam | 21488 | 5021 | 0.786 | 0.184 | 0.810 | 0.189 | 0.93 | 8 |
| v2 | 9 | source | 96396 | 26895 | 0.647 | 0.180 | 0.728 | 0.203 | 0.78 | 8 |
| v2 | 9 | fam | 13703 | 5021 | 0.651 | 0.239 | 0.673 | 0.247 | 0.93 | 9 |
| v2 | 10 | source | 63390 | 26895 | 0.512 | 0.217 | 0.588 | 0.249 | 0.78 | 8 |
| v2 | 10 | fam | 9242 | 5021 | 0.523 | 0.284 | 0.538 | 0.292 | 0.94 | 13 |
| v2 | 12 | source | 30671 | 26895 | 0.312 | 0.274 | 0.359 | 0.315 | 0.79 | 12 |
| v2 | 12 | fam | 4611 | 5021 | 0.329 | 0.358 | 0.337 | 0.367 | 0.95 | 19 |
| +floor | 8 | source | 149369 | 26895 | 0.770 | 0.139 | 0.849 | 0.153 | 0.78 | 6 |
| +floor | 8 | fam | 19944 | 5021 | 0.773 | 0.194 | 0.798 | 0.201 | 0.93 | 8 |
| +floor | 9 | source | 95974 | 26895 | 0.646 | 0.181 | 0.728 | 0.204 | 0.78 | 8 |
| +floor | 9 | fam | 13644 | 5021 | 0.651 | 0.240 | 0.672 | 0.247 | 0.93 | 9 |
| +floor | 10 | source | 63390 | 26895 | 0.512 | 0.217 | 0.588 | 0.249 | 0.78 | 8 |
| +floor | 10 | fam | 9242 | 5021 | 0.523 | 0.284 | 0.538 | 0.292 | 0.94 | 13 |
| +floor | 12 | source | 30671 | 26895 | 0.312 | 0.274 | 0.359 | 0.315 | 0.79 | 12 |
| +floor | 12 | fam | 4611 | 5021 | 0.329 | 0.358 | 0.337 | 0.367 | 0.95 | 19 |
| S1 | 8 | source | 39615 | 26895 | 0.539 | 0.366 | 0.820 | 0.557 | 0.66 | 13 |
| S1 | 8 | fam | 5380 | 5021 | 0.495 | 0.462 | 0.535 | 0.499 | 0.93 | 23 |
| S1 | 9 | source | 25242 | 26895 | 0.469 | 0.499 | 0.707 | 0.753 | 0.66 | 21 |
| S1 | 9 | fam | 3774 | 5021 | 0.435 | 0.578 | 0.468 | 0.623 | 0.93 | 26 |
| S1 | 10 | source | 17897 | 26895 | 0.383 | 0.575 | 0.574 | 0.862 | 0.67 | 22 |
| S1 | 10 | fam | 2765 | 5021 | 0.362 | 0.658 | 0.387 | 0.703 | 0.94 | 27 |
| S1 | 12 | source | 10048 | 26895 | 0.238 | 0.636 | 0.352 | 0.942 | 0.68 | 24 |
| S1 | 12 | fam | 1631 | 5021 | 0.233 | 0.719 | 0.248 | 0.762 | 0.94 | 26 |

separation at his 26895 test-day times (best within +-2 s): p50 10.11; fraction >= 8/9/10: 0.786/0.650/0.514

### mix: Text S1 rules, held-out test days (one-to-one +-2 s)

| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 8 | source | 139115 | 26895 | 0.768 | 0.148 | 0.850 | 0.164 | 0.80 | 7 |
| v2 | 8 | fam | 18939 | 5021 | 0.764 | 0.203 | 0.788 | 0.209 | 0.94 | 8 |
| v2 | 9 | source | 80403 | 26895 | 0.627 | 0.210 | 0.709 | 0.237 | 0.80 | 8 |
| v2 | 9 | fam | 11573 | 5021 | 0.628 | 0.272 | 0.649 | 0.282 | 0.94 | 11 |
| v2 | 10 | source | 51391 | 26895 | 0.491 | 0.257 | 0.560 | 0.293 | 0.81 | 9 |
| v2 | 10 | fam | 7613 | 5021 | 0.496 | 0.327 | 0.510 | 0.336 | 0.95 | 15 |
| v2 | 12 | source | 23751 | 26895 | 0.291 | 0.330 | 0.329 | 0.372 | 0.83 | 15 |
| v2 | 12 | fam | 3663 | 5021 | 0.301 | 0.413 | 0.309 | 0.423 | 0.96 | 22 |
| +floor | 8 | source | 129628 | 26895 | 0.756 | 0.157 | 0.840 | 0.174 | 0.80 | 7 |
| +floor | 8 | fam | 17568 | 5021 | 0.751 | 0.215 | 0.775 | 0.221 | 0.94 | 9 |
| +floor | 9 | source | 80091 | 26895 | 0.627 | 0.210 | 0.709 | 0.238 | 0.80 | 8 |
| +floor | 9 | fam | 11534 | 5021 | 0.628 | 0.273 | 0.648 | 0.282 | 0.95 | 11 |
| +floor | 10 | source | 51391 | 26895 | 0.491 | 0.257 | 0.560 | 0.293 | 0.81 | 9 |
| +floor | 10 | fam | 7613 | 5021 | 0.496 | 0.327 | 0.510 | 0.336 | 0.95 | 15 |
| +floor | 12 | source | 23751 | 26895 | 0.291 | 0.330 | 0.329 | 0.372 | 0.83 | 15 |
| +floor | 12 | fam | 3663 | 5021 | 0.301 | 0.413 | 0.309 | 0.423 | 0.96 | 22 |
| S1 | 8 | source | 38245 | 26895 | 0.553 | 0.389 | 0.812 | 0.571 | 0.68 | 14 |
| S1 | 8 | fam | 5286 | 5021 | 0.507 | 0.481 | 0.548 | 0.520 | 0.93 | 23 |
| S1 | 9 | source | 23894 | 26895 | 0.478 | 0.538 | 0.690 | 0.777 | 0.69 | 24 |
| S1 | 9 | fam | 3676 | 5021 | 0.443 | 0.605 | 0.477 | 0.651 | 0.93 | 26 |
| S1 | 10 | source | 16571 | 26895 | 0.388 | 0.630 | 0.546 | 0.887 | 0.71 | 26 |
| S1 | 10 | fam | 2606 | 5021 | 0.360 | 0.693 | 0.385 | 0.742 | 0.93 | 27 |
| S1 | 12 | source | 9075 | 26895 | 0.236 | 0.701 | 0.323 | 0.958 | 0.73 | 26 |
| S1 | 12 | fam | 1498 | 5021 | 0.227 | 0.761 | 0.241 | 0.808 | 0.94 | 27 |

### SH_V3B vs baseline (sh_v3): Text S1 rules, held-out test days

| iter | xMAD | mode | ours base | ours new | recall fam base | recall fam new | prec fam base | prec fam new | recall any base | recall any new | prec any base | prec any new |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| it0 | 8 | source | 33357 | 37973 | 0.410 | 0.562 | 0.331 | 0.398 | 0.642 | 0.809 | 0.518 | 0.573 |
| it0 | 8 | fam | 4596 | 5361 | 0.404 | 0.520 | 0.441 | 0.487 | 0.439 | 0.558 | 0.479 | 0.523 |
| it0 | 9 | source | 18058 | 23534 | 0.336 | 0.485 | 0.501 | 0.555 | 0.506 | 0.685 | 0.753 | 0.782 |
| it0 | 9 | fam | 2860 | 3651 | 0.334 | 0.453 | 0.586 | 0.623 | 0.359 | 0.483 | 0.630 | 0.665 |
| it0 | 10 | source | 11643 | 16307 | 0.263 | 0.392 | 0.608 | 0.647 | 0.381 | 0.541 | 0.881 | 0.892 |
| it0 | 10 | fam | 1936 | 2613 | 0.263 | 0.369 | 0.683 | 0.709 | 0.281 | 0.392 | 0.728 | 0.753 |
| it1 | 8 | source | 34859 | 39615 | 0.392 | 0.539 | 0.302 | 0.366 | 0.660 | 0.820 | 0.509 | 0.557 |
| it1 | 8 | fam | 4775 | 5380 | 0.386 | 0.495 | 0.406 | 0.462 | 0.426 | 0.535 | 0.448 | 0.499 |
| it1 | 9 | source | 19449 | 25242 | 0.324 | 0.469 | 0.448 | 0.499 | 0.531 | 0.707 | 0.734 | 0.753 |
| it1 | 9 | fam | 2992 | 3774 | 0.320 | 0.435 | 0.537 | 0.578 | 0.352 | 0.468 | 0.590 | 0.623 |
| it1 | 10 | source | 12781 | 17897 | 0.255 | 0.383 | 0.537 | 0.575 | 0.408 | 0.574 | 0.859 | 0.862 |
| it1 | 10 | fam | 2032 | 2765 | 0.254 | 0.362 | 0.627 | 0.658 | 0.274 | 0.387 | 0.678 | 0.703 |
| mix | 8 | source | 33484 | 38245 | 0.405 | 0.553 | 0.325 | 0.389 | 0.647 | 0.812 | 0.520 | 0.571 |
| mix | 8 | fam | 4558 | 5286 | 0.392 | 0.507 | 0.432 | 0.481 | 0.431 | 0.548 | 0.475 | 0.520 |
| mix | 9 | source | 18285 | 23894 | 0.333 | 0.478 | 0.489 | 0.538 | 0.511 | 0.690 | 0.752 | 0.777 |
| mix | 9 | fam | 2812 | 3676 | 0.323 | 0.443 | 0.577 | 0.605 | 0.351 | 0.477 | 0.626 | 0.651 |
| mix | 10 | source | 11881 | 16571 | 0.262 | 0.388 | 0.594 | 0.630 | 0.389 | 0.546 | 0.880 | 0.887 |
| mix | 10 | fam | 1942 | 2606 | 0.260 | 0.360 | 0.671 | 0.693 | 0.278 | 0.385 | 0.719 | 0.742 |

per-family recall change, mix, S1 9x (family-specific one-to-one; his test-day events per family)

| fam | his (source) | recall base | recall new | delta | prec base | prec new | his (fam days) | recall base (fam) | recall new (fam) | delta (fam) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 954 | 0.265 | 0.395 | +0.130 | 0.444 | 0.493 | 197 | 0.234 | 0.376 | +0.142 |
| 2 | 1157 | 0.174 | 0.242 | +0.068 | 0.369 | 0.412 | 203 | 0.177 | 0.251 | +0.074 |
| 3 | 591 | 0.237 | 0.284 | +0.047 | 0.373 | 0.444 | 199 | 0.281 | 0.357 | +0.075 |
| 4 | 400 | 0.177 | 0.258 | +0.080 | 0.206 | 0.320 | 130 | 0.108 | 0.146 | +0.038 |
| 5 | 586 | 0.143 | 0.299 | +0.155 | 0.275 | 0.438 | 167 | 0.168 | 0.293 | +0.126 |
| 6 | 1358 | 0.416 | 0.591 | +0.175 | 0.661 | 0.741 | 194 | 0.397 | 0.603 | +0.206 |
| 7 | 1072 | 0.299 | 0.411 | +0.113 | 0.505 | 0.527 | 195 | 0.292 | 0.415 | +0.123 |
| 8 | 1691 | 0.512 | 0.710 | +0.198 | 0.858 | 0.847 | 195 | 0.451 | 0.692 | +0.241 |
| 9 | 319 | 0.229 | 0.285 | +0.056 | 0.104 | 0.141 | 108 | 0.231 | 0.241 | +0.009 |
| 10 | 798 | 0.327 | 0.393 | +0.066 | 0.400 | 0.414 | 194 | 0.345 | 0.397 | +0.052 |
| 11 | 298 | 0.326 | 0.305 | -0.020 | 0.180 | 0.260 | 143 | 0.301 | 0.301 | +0.000 |
| 12 | 484 | 0.430 | 0.496 | +0.066 | 0.289 | 0.299 | 162 | 0.444 | 0.500 | +0.056 |
| 13 | 876 | 0.301 | 0.424 | +0.122 | 0.446 | 0.498 | 200 | 0.295 | 0.420 | +0.125 |
| 14 | 1415 | 0.424 | 0.462 | +0.038 | 0.590 | 0.686 | 200 | 0.395 | 0.475 | +0.080 |
| 15 | 862 | 0.275 | 0.423 | +0.148 | 0.342 | 0.386 | 207 | 0.275 | 0.411 | +0.135 |
| 16 | 601 | 0.295 | 0.429 | +0.135 | 0.277 | 0.321 | 176 | 0.256 | 0.375 | +0.119 |
| 17 | 621 | 0.258 | 0.380 | +0.122 | 0.394 | 0.428 | 186 | 0.204 | 0.328 | +0.124 |
| 18 | 657 | 0.402 | 0.472 | +0.070 | 0.318 | 0.308 | 186 | 0.457 | 0.484 | +0.027 |
| 19 | 1811 | 0.512 | 0.644 | +0.133 | 0.779 | 0.780 | 189 | 0.630 | 0.693 | +0.063 |
| 20 | 955 | 0.284 | 0.374 | +0.090 | 0.342 | 0.331 | 198 | 0.247 | 0.308 | +0.061 |
| 21 | 647 | 0.457 | 0.620 | +0.162 | 0.648 | 0.681 | 196 | 0.485 | 0.587 | +0.102 |
| 22 | 1244 | 0.474 | 0.609 | +0.135 | 0.539 | 0.520 | 181 | 0.536 | 0.657 | +0.122 |
| 23 | 2852 | 0.223 | 0.485 | +0.262 | 0.731 | 0.777 | 210 | 0.252 | 0.519 | +0.267 |
| 24 | 1095 | 0.342 | 0.422 | +0.080 | 0.463 | 0.435 | 192 | 0.349 | 0.432 | +0.083 |
| 25 | 1129 | 0.366 | 0.439 | +0.074 | 0.525 | 0.569 | 195 | 0.354 | 0.426 | +0.072 |
| 26 | 553 | 0.378 | 0.514 | +0.136 | 0.666 | 0.633 | 207 | 0.338 | 0.517 | +0.179 |
| 27 | 1869 | 0.210 | 0.577 | +0.368 | 0.708 | 0.647 | 211 | 0.152 | 0.531 | +0.379 |

mix: families with recall up / down / same (source, S1 9x): 26 / 1 / 0; median delta +0.122

per-family recall change, it0, S1 9x (family-specific one-to-one; his test-day events per family)

| fam | his (source) | recall base | recall new | delta | prec base | prec new | his (fam days) | recall base (fam) | recall new (fam) | delta (fam) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 954 | 0.264 | 0.397 | +0.133 | 0.445 | 0.493 | 197 | 0.234 | 0.376 | +0.142 |
| 2 | 1157 | 0.175 | 0.249 | +0.074 | 0.372 | 0.409 | 203 | 0.182 | 0.261 | +0.079 |
| 3 | 591 | 0.274 | 0.321 | +0.047 | 0.386 | 0.516 | 199 | 0.327 | 0.397 | +0.070 |
| 4 | 400 | 0.182 | 0.273 | +0.090 | 0.196 | 0.311 | 130 | 0.115 | 0.169 | +0.054 |
| 5 | 586 | 0.147 | 0.300 | +0.154 | 0.281 | 0.440 | 167 | 0.174 | 0.299 | +0.126 |
| 6 | 1358 | 0.416 | 0.591 | +0.175 | 0.655 | 0.741 | 194 | 0.397 | 0.603 | +0.206 |
| 7 | 1072 | 0.299 | 0.421 | +0.121 | 0.498 | 0.521 | 195 | 0.282 | 0.415 | +0.133 |
| 8 | 1691 | 0.512 | 0.710 | +0.198 | 0.857 | 0.847 | 195 | 0.451 | 0.692 | +0.241 |
| 9 | 319 | 0.194 | 0.273 | +0.078 | 0.120 | 0.299 | 108 | 0.167 | 0.213 | +0.046 |
| 10 | 798 | 0.338 | 0.405 | +0.066 | 0.398 | 0.417 | 194 | 0.351 | 0.433 | +0.082 |
| 11 | 298 | 0.352 | 0.376 | +0.023 | 0.244 | 0.381 | 143 | 0.329 | 0.350 | +0.021 |
| 12 | 484 | 0.432 | 0.500 | +0.068 | 0.289 | 0.297 | 162 | 0.438 | 0.506 | +0.068 |
| 13 | 876 | 0.396 | 0.539 | +0.143 | 0.473 | 0.515 | 200 | 0.390 | 0.515 | +0.125 |
| 14 | 1415 | 0.428 | 0.461 | +0.034 | 0.602 | 0.683 | 200 | 0.430 | 0.470 | +0.040 |
| 15 | 862 | 0.284 | 0.436 | +0.152 | 0.343 | 0.388 | 207 | 0.285 | 0.420 | +0.135 |
| 16 | 601 | 0.296 | 0.431 | +0.135 | 0.276 | 0.321 | 176 | 0.261 | 0.375 | +0.114 |
| 17 | 621 | 0.320 | 0.467 | +0.147 | 0.403 | 0.433 | 186 | 0.269 | 0.409 | +0.140 |
| 18 | 657 | 0.403 | 0.472 | +0.068 | 0.317 | 0.304 | 186 | 0.452 | 0.489 | +0.038 |
| 19 | 1811 | 0.513 | 0.646 | +0.133 | 0.776 | 0.780 | 189 | 0.630 | 0.698 | +0.069 |
| 20 | 955 | 0.287 | 0.375 | +0.088 | 0.342 | 0.332 | 198 | 0.253 | 0.308 | +0.056 |
| 21 | 647 | 0.454 | 0.623 | +0.168 | 0.634 | 0.681 | 196 | 0.480 | 0.587 | +0.107 |
| 22 | 1244 | 0.490 | 0.621 | +0.131 | 0.608 | 0.598 | 181 | 0.575 | 0.685 | +0.110 |
| 23 | 2852 | 0.194 | 0.484 | +0.290 | 0.810 | 0.776 | 210 | 0.243 | 0.519 | +0.276 |
| 24 | 1095 | 0.343 | 0.421 | +0.078 | 0.465 | 0.434 | 192 | 0.349 | 0.432 | +0.083 |
| 25 | 1129 | 0.365 | 0.444 | +0.079 | 0.511 | 0.565 | 195 | 0.354 | 0.415 | +0.062 |
| 26 | 553 | 0.378 | 0.514 | +0.136 | 0.661 | 0.633 | 207 | 0.338 | 0.517 | +0.179 |
| 27 | 1869 | 0.202 | 0.535 | +0.333 | 0.786 | 0.702 | 211 | 0.152 | 0.455 | +0.303 |

it0: families with recall up / down / same (source, S1 9x): 27 / 0 / 0; median delta +0.131
<!-- SH_V3_RESULTS_END -->
