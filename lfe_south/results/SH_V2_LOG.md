# SH_V2 — Shelly (2025) Mendocino reproduction, detector v2 (agent, 2026-09-28)

Inputs read: `results/aa_shelly/AA_SHELLY_REPORT.md`, `notes/PAPER_METHODS_TABLE.md` (S25 rows),
`notes/CATALOGUE_SUMMARY.md` §10, `notes/2026-09-25_Notes.md` §Q/§V. Nothing lfe_south was running
(ps 22:5xZ: only `ab_grow_fast.py --sub NW` NW_lag15 on the GPU, 38/96 GB GPU mem, anon 15.8 GB).

## Step 1 — diagnosis of the ~10,000-detection family-days (row C = 6 s, 2–9 Hz, `mfC_diag.csv`)

**Cause: dead channels carrying float residue plus a multi-million-count glitch, in the
day-long cumsum energy.** The AA guard (`std > 1e-4 x median std`, `>=1 sample change`) does not see
them: the channel is not constant (residue ~1e-13 changes every sample) and the cumsum S2 after a
~4e6-count glitch has absolute round-off ~1e-3, so residue windows read a garbage non-zero std and
pass. Result: per-channel "normalised" |cc| up to 50–211.

| measurement | value |
|---|---|
| family-days with stat_max > 1 (mean of normalised cc — impossible if cc is bounded) | **664 / 10,395 (6.4 %)** |
| detections at 9xMAD from those family-days | **1,859,401 of 1,886,890 (98.5 %)** — mean 2,800 / family-day |
| distinct days involved | 27 (all 26–27 families fire on each; e.g. 2023-09-21: 10,128–10,131 det per family at every threshold 6–12) |
| channel responsible, every one of the 27 days (scan of all 23 stations, all channels) | **BK.BRIC HHE** (2023-01..2024-05, 19 days), **BK.RVIT HHE** (2022-01/02, 6 days), **NC.KCR EHZ + NC.KSM EHZ** (2019-09..11, 4 days; both dead together) — always cache channel index 0 |
| those channels: median |x| / max |x| | 1e-17 … 3e-13 (dead residue) / 6e2 … 7.5e6 (glitch; BRIC 2023-09-21: 13,312 samples > 1e3, peak 3.9e6) |
| per-channel |cc|max / samples with |cc|>1 | 31–211 / 1.0e5–1.36e6 per day (BRIC 2023-09-21: 203.6, 749,142) |
| exact (np.corrcoef) cc at the same argmax | 0.005 (the true value) |
| clean-day channels (5 quiet + 5 active days) with |cc|>1 | **0** (BRIC HHE is also dead on 2 of them but without a glitch -> harmless there; still counted as a "live" channel) |
| burst blocks (60 s rms > 10x day median), mean fraction | quiet days 0.0079, active 0.0052 — not the cause |

**Why quiet ≈ active.** It is entirely the 664 broken family-days (broken quiet 2,672 vs broken
active 2,881 det/family-day at 9x). On the 9,731 clean family-days the detector already separates:

| xMAD | clean quiet (n=3,428) | clean active (n=6,303) | Shelly active |
|---|---|---|---|
| 6 | 101.3 | 125.4 | |
| 7 | 13.0 | 23.6 | |
| 8 | 2.62 | 7.85 | |
| **9** | **0.93** | **3.86** | **4.01** |
| 10 | 0.49 | 2.25 | |
| 12 | 0.22 | 0.95 | |

(6xMAD ≈ 4σ — ~100/day is the Gaussian tail over 1.7 M samples, not a pathology.)

**Row D (template iteration) collapse explained.** `aa_restack.py` took the top-200 net_mad
detections per family from `mfA_build.parquet`: **100 % of them came from broken family-days**
(541/10,962 build family-days carry 1,733,362 of 1,787,349 build detections; seed "meancc" median
6.8). The re-stacked templates were stacks of glitch windows -> separation p50 fell 7.4 -> 2.9.
Also, the stacker has no guard: dead-channel windows are rms-normalised to unit amplitude and stacked.

Diagnostic scripts/outputs (scratch): `diag1.py`, `diag2.py`, `scan.csv` in the session scratchpad.

## Steps 2–3 — detector v2 and the literal recipe (queue launched 2026-09-28 23:17Z)

New files only (nothing edited): `scripts/sh_v2_common.py` (cleaning), `sh_v2_prep.py`, `t11_mf_v2.py`,
`sh_v2_restack.py`, `sh_v2_report.py`, `sh_v2_queue.sh`. Outputs `results/sh_v2/`.

**Cleaning (per station-day, identical in detector and every template stack):**
1. flat runs >= 0.5 s on the raw cache invalid, padded 2 s (filter transient); 2 Hz high-pass -> 2–9 Hz.
2. **dead floors first** — a channel whose plain daily-median 2 s energy < 1e-4 x its station's loudest
   channel is dead for the day (ref = 0, excluded from the burst ratio); station floor 1e-6 x the day's
   median station level. (Needed beyond the WA `channel_refs`: there, BRIC HHE's 11-min glitch alone
   passes the 10-min live rule and becomes the channel's reference.)
3. burst muting, ab_mute.py port with the NW_mute values (K=15, 2 s Q<2, 60 s Q<1.3, pad 3 s, 5
   nearest other sites, coordinates from `fdsn_channels_south.csv`); line notching (prominence >= 4,
   2–9 Hz, 0.06 Hz notches).
4. exact block-local window energies (no day-long cumsum), relative guard e > 1e-3 ref_c, window valid
   only if every sample valid; per-sample live-channel count; any residual |cc|>1 counted and zeroed.

**Smoke test (5 days: 3 pathology, 1 quiet, 1 active; row-C templates hp-filtered):** detections at
9x per day (all 27 families) new vs row C vs Shelly: 2023-09-21 **51** vs 263,341 vs 68;
2019-09-03 184 vs 22,619 vs 85; 2022-02-19 44 vs 22,409 vs 28; quiet 2023-02-16 **2** vs 22 vs 7;
active 2022-06-16 224 vs 165 vs 219. |cc|>1 samples 0; family-days with mean cc>1: 0. Separation at
his 407 times: p50 7.58 (row C 6.78), fraction >= 9: 0.334 (0.251). BK.BRIC HHE, BK.RVIT HHE, NC.KCR,
NC.KSM all flagged dead on their pathology days. Muted fraction 3–7 % per day (B933 58 % on 2019-09-03)
— the queue's report checks his S-arrival windows' muted fraction vs the station-day fraction.

**Recipe (literal where the main text states it; SILENT items flagged):** 2–9 Hz; 6 s templates, 2 s
pre (row C cut, stage-1 offsets tt_s+best_tau, not re-measured); network **correlation SUM** over live
channels with station moveout; threshold x DAILY MAD, sweep 6–12 as columns (9 = his); fixed 4 s
bins. Iteration 0 = his published times on each family's OWN build days; iterations 1–3 = the 200
highest correlation sums >= 9x MAD of our own build-day detections, family's own build days only
(N=200 and the P/S window rule are SILENT in the main text — S1 not on disk; one 6 s window on all
channels kept). Scored on held-out test days: `--days-mode source` (as AA) and `--days-mode fam`
(each family only on its own test days — strict; in source mode a family is also scored on days that
were other families' test days but its OWN build days, i.e. days its templates saw).

**Queue:** `scripts/sh_v2_queue.sh`, pid 3328101 (setsid, lock `locks/sh_v2.lock`, heartbeat child
closes fd 9, killed in EXIT trap), log `logs/sh_v2_queue.log`, markers `results/SH_V2_{DONE,FAILED,HEARTBEAT}`.
Order: prep (572 days) -> tpl_it0 -> [mf test, score source+fam, report, mf build, restack] x it0..it3.
Estimated ~7 h (MF ~8 s/day on GPU, restack 4 workers). Results are appended below automatically.

Step 4 (locations vs his) not started.

<!-- SH_V2_RESULTS_BEGIN -->
## SH_V2 results (auto, sh_v2_report.py, 2026-09-29 06:59Z)

### prep (mute / lines / dead floors)

station-days 6805; muted fraction mean 0.0513 (median 0.0304); station-days with a dead channel 210; dead stations 22; lines notched 5718 (station-days with lines 414)
most-muted stations: PB.B933 0.121, BK.BRIC 0.118, PB.B934 0.083, PB.B935 0.060, NC.KIP 0.053, BK.HALS 0.042, BK.RVIT 0.039, NC.KSM 0.037
dead-channel station-days by station: BK.BRIC 163, BK.RVIT 25, NC.KCR 10, NC.KSM 7, NC.KMR 2, BK.HALS 1, BK.JEWT 1, NC.KIP 1
his S-arrival windows muted: 0.0369 (n=5905) vs station-day muted fraction 0.0664 (60 random catalogue days)

### iteration 0 (templates `sh_v2/tpl_it0`), test days

family-days 10395; |cc|>1 samples 0; family-days with mean cc > 1: 0; median live channels 18; max det/family-day at 9x 42

| xMAD | quiet det/fam-day (n=3684) | active det/fam-day (n=6711) | Shelly active | max det/fam-day |
|---|---|---|---|---|
| 6 | 124.15 | 151.65 | 4.01 | 963 |
| 7 | 16.65 | 29.48 | 4.01 | 255 |
| 8 | 3.08 | 9.42 | 4.01 | 98 |
| 9 | 0.84 | 4.34 | 4.01 | 42 |
| 10 | 0.32 | 2.39 | 4.01 | 33 |
| 11 | 0.14 | 1.47 | 4.01 | 20 |
| 12 | 0.08 | 0.94 | 4.01 | 15 |

score, days-mode source (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 26895 | 1475095 | 19728 | 0.734 | 0.013 | 27/27 | 0 |
| 7 | 26895 | 259174 | 15894 | 0.591 | 0.061 | 27/27 | 0 |
| 8 | 26895 | 74583 | 12272 | 0.456 | 0.165 | 27/27 | 3 |
| 9 | 26895 | 32208 | 9082 | 0.338 | 0.282 | 27/27 | 8 |
| 10 | 26895 | 17237 | 6569 | 0.244 | 0.381 | 27/27 | 16 |
| 11 | 26895 | 10378 | 4755 | 0.177 | 0.458 | 27/27 | 19 |
| 12 | 26895 | 6566 | 3374 | 0.125 | 0.514 | 27/27 | 21 |

per family at 9x (source): 1: R 0.34 P 0.22 (324/954/1459); 2: R 0.28 P 0.22 (326/1157/1507); 3: R 0.26 P 0.20 (156/591/768); 4: R 0.16 P 0.09 (64/400/689); 5: R 0.11 P 0.10 (65/586/649); 6: R 0.30 P 0.48 (405/1358/852); 7: R 0.33 P 0.26 (351/1072/1344); 8: R 0.54 P 0.73 (910/1691/1245); 9: R 0.06 P 0.11 (20/319/186); 10: R 0.37 P 0.19 (297/798/1539); 11: R 0.37 P 0.12 (110/298/906); 12: R 0.32 P 0.23 (153/484/680); 13: R 0.44 P 0.26 (389/876/1482); 14: R 0.32 P 0.40 (456/1415/1133); 15: R 0.21 P 0.20 (183/862/898); 16: R 0.35 P 0.09 (209/601/2394); 17: R 0.38 P 0.19 (233/621/1219); 18: R 0.32 P 0.21 (210/657/1012); 19: R 0.61 P 0.45 (1101/1811/2437); 20: R 0.32 P 0.19 (308/955/1608); 21: R 0.44 P 0.28 (284/647/1029); 22: R 0.46 P 0.36 (568/1244/1571); 23: R 0.20 P 0.70 (582/2852/828); 24: R 0.55 P 0.25 (603/1095/2442); 25: R 0.36 P 0.27 (409/1129/1509); 26: R 0.36 P 0.35 (201/553/571); 27: R 0.09 P 0.66 (165/1869/251)

score, days-mode fam (own test days only) (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 5021 | 114459 | 3713 | 0.739 | 0.032 | 27/27 | 0 |
| 7 | 5021 | 25202 | 2969 | 0.591 | 0.118 | 27/27 | 2 |
| 8 | 5021 | 9054 | 2262 | 0.451 | 0.250 | 27/27 | 8 |
| 9 | 5021 | 4427 | 1620 | 0.323 | 0.366 | 27/27 | 19 |
| 10 | 5021 | 2512 | 1157 | 0.230 | 0.461 | 27/27 | 21 |
| 11 | 5021 | 1578 | 848 | 0.169 | 0.537 | 27/27 | 23 |
| 12 | 5021 | 1052 | 624 | 0.124 | 0.593 | 27/27 | 25 |

per family at 9x (fam (own test days only)): 1: R 0.32 P 0.32 (64/197/198); 2: R 0.29 P 0.25 (58/203/235); 3: R 0.21 P 0.33 (41/199/126); 4: R 0.15 P 0.23 (20/130/87); 5: R 0.11 P 0.14 (18/167/130); 6: R 0.31 P 0.66 (61/194/92); 7: R 0.34 P 0.39 (66/195/168); 8: R 0.48 P 0.67 (93/195/138); 9: R 0.06 P 0.09 (6/108/65); 10: R 0.38 P 0.33 (73/194/221); 11: R 0.35 P 0.19 (50/143/261); 12: R 0.33 P 0.46 (53/162/115); 13: R 0.41 P 0.34 (82/200/241); 14: R 0.27 P 0.54 (54/200/100); 15: R 0.17 P 0.33 (35/207/105); 16: R 0.31 P 0.19 (55/176/284); 17: R 0.30 P 0.28 (56/186/197); 18: R 0.30 P 0.30 (55/186/186); 19: R 0.74 P 0.62 (140/189/224); 20: R 0.33 P 0.30 (65/198/214); 21: R 0.43 P 0.49 (84/196/172); 22: R 0.48 P 0.52 (86/181/164); 23: R 0.25 P 0.79 (53/210/67); 24: R 0.56 P 0.37 (107/192/292); 25: R 0.30 P 0.33 (58/195/174); 26: R 0.34 P 0.47 (71/207/150); 27: R 0.08 P 0.76 (16/211/21)

separation: statistic at his 26895 test-day times (best within +-2 s): p10/p50/p90 4.75/7.68/12.66; fraction >= 6..12: 6:0.742, 7:0.595, 8:0.458, 9:0.339, 10:0.245, 11:0.177, 12:0.126; dt median 0.00 s; noise p99.9 median 4.77

### iteration 1 (templates `sh_v2/tpl_it1`), test days

family-days 10395; |cc|>1 samples 0; family-days with mean cc > 1: 0; median live channels 18; max det/family-day at 9x 82

| xMAD | quiet det/fam-day (n=3684) | active det/fam-day (n=6711) | Shelly active | max det/fam-day |
|---|---|---|---|---|
| 6 | 149.54 | 183.46 | 4.01 | 3195 |
| 7 | 24.07 | 39.62 | 4.01 | 1046 |
| 8 | 4.99 | 12.69 | 4.01 | 295 |
| 9 | 1.40 | 5.74 | 4.01 | 82 |
| 10 | 0.52 | 3.21 | 4.01 | 34 |
| 11 | 0.25 | 2.01 | 4.01 | 30 |
| 12 | 0.14 | 1.34 | 4.01 | 22 |

score, days-mode source (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 26895 | 1782087 | 19871 | 0.739 | 0.011 | 27/27 | 0 |
| 7 | 26895 | 354575 | 16212 | 0.603 | 0.046 | 27/27 | 0 |
| 8 | 26895 | 103558 | 12685 | 0.472 | 0.122 | 27/27 | 3 |
| 9 | 26895 | 43664 | 9675 | 0.360 | 0.222 | 27/27 | 6 |
| 10 | 26895 | 23452 | 7265 | 0.270 | 0.310 | 27/27 | 8 |
| 11 | 26895 | 14409 | 5439 | 0.202 | 0.377 | 27/27 | 15 |
| 12 | 26895 | 9479 | 4026 | 0.150 | 0.425 | 27/27 | 17 |

per family at 9x (source): 1: R 0.35 P 0.19 (337/954/1787); 2: R 0.28 P 0.17 (328/1157/1880); 3: R 0.28 P 0.19 (167/591/870); 4: R 0.20 P 0.08 (82/400/1082); 5: R 0.13 P 0.03 (74/586/2187); 6: R 0.34 P 0.36 (462/1358/1282); 7: R 0.38 P 0.22 (407/1072/1820); 8: R 0.54 P 0.71 (920/1691/1291); 9: R 0.10 P 0.08 (33/319/392); 10: R 0.36 P 0.15 (290/798/1893); 11: R 0.41 P 0.11 (123/298/1117); 12: R 0.35 P 0.16 (171/484/1099); 13: R 0.47 P 0.22 (412/876/1842); 14: R 0.36 P 0.34 (510/1415/1513); 15: R 0.26 P 0.15 (227/862/1512); 16: R 0.39 P 0.08 (233/601/3002); 17: R 0.40 P 0.15 (248/621/1671); 18: R 0.37 P 0.15 (240/657/1613); 19: R 0.62 P 0.41 (1118/1811/2758); 20: R 0.38 P 0.18 (363/955/2072); 21: R 0.45 P 0.18 (292/647/1624); 22: R 0.47 P 0.26 (589/1244/2277); 23: R 0.23 P 0.68 (670/2852/988); 24: R 0.56 P 0.20 (617/1095/3038); 25: R 0.37 P 0.22 (421/1129/1907); 26: R 0.36 P 0.21 (198/553/936); 27: R 0.08 P 0.68 (143/1869/211)

score, days-mode fam (own test days only) (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 5021 | 139621 | 3766 | 0.750 | 0.027 | 27/27 | 0 |
| 7 | 5021 | 33369 | 3095 | 0.616 | 0.093 | 27/27 | 2 |
| 8 | 5021 | 12085 | 2429 | 0.484 | 0.201 | 27/27 | 7 |
| 9 | 5021 | 5795 | 1803 | 0.359 | 0.311 | 27/27 | 15 |
| 10 | 5021 | 3315 | 1323 | 0.263 | 0.399 | 27/27 | 20 |
| 11 | 5021 | 2053 | 979 | 0.195 | 0.477 | 27/27 | 24 |
| 12 | 5021 | 1380 | 713 | 0.142 | 0.517 | 27/27 | 22 |

per family at 9x (fam (own test days only)): 1: R 0.35 P 0.28 (68/197/239); 2: R 0.30 P 0.22 (60/203/271); 3: R 0.25 P 0.33 (49/199/148); 4: R 0.19 P 0.17 (25/130/145); 5: R 0.12 P 0.08 (20/167/252); 6: R 0.37 P 0.55 (71/194/129); 7: R 0.42 P 0.34 (82/195/239); 8: R 0.48 P 0.66 (93/195/140); 9: R 0.10 P 0.10 (11/108/114); 10: R 0.39 P 0.28 (75/194/264); 11: R 0.43 P 0.18 (61/143/347); 12: R 0.35 P 0.30 (57/162/187); 13: R 0.46 P 0.32 (92/200/285); 14: R 0.31 P 0.50 (62/200/124); 15: R 0.24 P 0.31 (49/207/159); 16: R 0.36 P 0.18 (63/176/355); 17: R 0.33 P 0.24 (62/186/257); 18: R 0.37 P 0.24 (68/186/286); 19: R 0.75 P 0.59 (142/189/242); 20: R 0.40 P 0.28 (79/198/285); 21: R 0.45 P 0.39 (88/196/225); 22: R 0.50 P 0.45 (91/181/204); 23: R 0.31 P 0.79 (66/210/84); 24: R 0.58 P 0.32 (112/192/351); 25: R 0.35 P 0.30 (69/195/232); 26: R 0.35 P 0.34 (73/207/213); 27: R 0.07 P 0.83 (15/211/18)

separation: statistic at his 26895 test-day times (best within +-2 s): p10/p50/p90 4.77/7.80/13.25; fraction >= 6..12: 6:0.749, 7:0.608, 8:0.474, 9:0.361, 10:0.271, 11:0.202, 12:0.150; dt median 0.00 s; noise p99.9 median 4.80

### iteration 2 (templates `sh_v2/tpl_it2`), test days

family-days 10395; |cc|>1 samples 0; family-days with mean cc > 1: 0; median live channels 18; max det/family-day at 9x 1218

| xMAD | quiet det/fam-day (n=3684) | active det/fam-day (n=6711) | Shelly active | max det/fam-day |
|---|---|---|---|---|
| 6 | 176.15 | 216.98 | 4.01 | 6927 |
| 7 | 36.61 | 56.87 | 4.01 | 3862 |
| 8 | 10.18 | 20.11 | 4.01 | 1824 |
| 9 | 3.33 | 8.42 | 4.01 | 1218 |
| 10 | 1.22 | 4.09 | 4.01 | 648 |
| 11 | 0.47 | 2.30 | 4.01 | 263 |
| 12 | 0.19 | 1.43 | 4.01 | 80 |

score, days-mode source (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 26895 | 2105090 | 19771 | 0.735 | 0.009 | 27/27 | 0 |
| 7 | 26895 | 516512 | 16064 | 0.597 | 0.031 | 27/27 | 0 |
| 8 | 26895 | 172460 | 12559 | 0.467 | 0.073 | 27/27 | 3 |
| 9 | 26895 | 68749 | 9555 | 0.355 | 0.139 | 27/27 | 5 |
| 10 | 26895 | 31927 | 7232 | 0.269 | 0.227 | 27/27 | 8 |
| 11 | 26895 | 17156 | 5365 | 0.199 | 0.313 | 27/27 | 14 |
| 12 | 26895 | 10331 | 3983 | 0.148 | 0.386 | 27/27 | 15 |

per family at 9x (source): 1: R 0.35 P 0.18 (330/954/1869); 2: R 0.28 P 0.17 (322/1157/1943); 3: R 0.28 P 0.18 (163/591/886); 4: R 0.22 P 0.07 (87/400/1184); 5: R 0.11 P 0.00 (63/586/25009); 6: R 0.34 P 0.29 (467/1358/1630); 7: R 0.38 P 0.22 (404/1072/1847); 8: R 0.54 P 0.71 (914/1691/1282); 9: R 0.10 P 0.07 (31/319/463); 10: R 0.36 P 0.16 (290/798/1839); 11: R 0.40 P 0.11 (118/298/1121); 12: R 0.36 P 0.12 (176/484/1447); 13: R 0.46 P 0.22 (407/876/1874); 14: R 0.34 P 0.30 (488/1415/1615); 15: R 0.26 P 0.14 (220/862/1623); 16: R 0.39 P 0.08 (236/601/3046); 17: R 0.40 P 0.15 (248/621/1700); 18: R 0.37 P 0.13 (243/657/1859); 19: R 0.61 P 0.40 (1101/1811/2720); 20: R 0.38 P 0.18 (361/955/2045); 21: R 0.45 P 0.16 (288/647/1781); 22: R 0.47 P 0.24 (583/1244/2460); 23: R 0.24 P 0.66 (685/2852/1043); 24: R 0.56 P 0.20 (611/1095/3098); 25: R 0.34 P 0.20 (382/1129/1924); 26: R 0.32 P 0.15 (176/553/1209); 27: R 0.09 P 0.69 (161/1869/232)

score, days-mode fam (own test days only) (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 5021 | 171940 | 3738 | 0.744 | 0.022 | 27/27 | 0 |
| 7 | 5021 | 49715 | 3092 | 0.616 | 0.062 | 27/27 | 2 |
| 8 | 5021 | 19218 | 2423 | 0.483 | 0.126 | 27/27 | 6 |
| 9 | 5021 | 8534 | 1798 | 0.358 | 0.211 | 27/27 | 12 |
| 10 | 5021 | 4204 | 1326 | 0.264 | 0.315 | 27/27 | 20 |
| 11 | 5021 | 2299 | 965 | 0.192 | 0.420 | 27/27 | 22 |
| 12 | 5021 | 1441 | 706 | 0.141 | 0.490 | 27/27 | 22 |

per family at 9x (fam (own test days only)): 1: R 0.35 P 0.27 (69/197/251); 2: R 0.30 P 0.22 (61/203/275); 3: R 0.22 P 0.31 (43/199/138); 4: R 0.23 P 0.16 (30/130/182); 5: R 0.11 P 0.01 (18/167/2695); 6: R 0.37 P 0.47 (71/194/151); 7: R 0.42 P 0.32 (82/195/253); 8: R 0.49 P 0.67 (95/195/142); 9: R 0.11 P 0.09 (12/108/132); 10: R 0.39 P 0.28 (75/194/268); 11: R 0.43 P 0.18 (62/143/342); 12: R 0.36 P 0.24 (59/162/246); 13: R 0.47 P 0.32 (94/200/293); 14: R 0.29 P 0.48 (59/200/124); 15: R 0.24 P 0.29 (49/207/169); 16: R 0.36 P 0.18 (64/176/356); 17: R 0.33 P 0.23 (62/186/270); 18: R 0.37 P 0.20 (68/186/337); 19: R 0.74 P 0.58 (140/189/241); 20: R 0.41 P 0.29 (81/198/283); 21: R 0.46 P 0.38 (90/196/237); 22: R 0.51 P 0.44 (93/181/209); 23: R 0.30 P 0.81 (63/210/78); 24: R 0.58 P 0.31 (112/192/359); 25: R 0.32 P 0.27 (63/195/231); 26: R 0.32 P 0.26 (66/207/251); 27: R 0.08 P 0.81 (17/211/21)

separation: statistic at his 26895 test-day times (best within +-2 s): p10/p50/p90 4.76/7.77/13.20; fraction >= 6..12: 6:0.746, 7:0.603, 8:0.470, 9:0.357, 10:0.269, 11:0.200, 12:0.148; dt median 0.00 s; noise p99.9 median 4.81

### iteration 3 (templates `sh_v2/tpl_it3`), test days

family-days 10395; |cc|>1 samples 0; family-days with mean cc > 1: 0; median live channels 18; max det/family-day at 9x 1397

| xMAD | quiet det/fam-day (n=3684) | active det/fam-day (n=6711) | Shelly active | max det/fam-day |
|---|---|---|---|---|
| 6 | 190.56 | 232.14 | 4.01 | 6974 |
| 7 | 41.27 | 62.49 | 4.01 | 3840 |
| 8 | 11.64 | 22.31 | 4.01 | 1977 |
| 9 | 3.91 | 9.37 | 4.01 | 1397 |
| 10 | 1.45 | 4.48 | 4.01 | 839 |
| 11 | 0.57 | 2.43 | 4.01 | 391 |
| 12 | 0.24 | 1.47 | 4.01 | 138 |

score, days-mode source (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 26895 | 2259923 | 19653 | 0.731 | 0.009 | 27/27 | 0 |
| 7 | 26895 | 571459 | 15887 | 0.591 | 0.028 | 27/27 | 0 |
| 8 | 26895 | 192570 | 12432 | 0.462 | 0.065 | 27/27 | 3 |
| 9 | 26895 | 77306 | 9410 | 0.350 | 0.122 | 27/27 | 4 |
| 10 | 26895 | 35398 | 7100 | 0.264 | 0.201 | 27/27 | 7 |
| 11 | 26895 | 18382 | 5301 | 0.197 | 0.288 | 27/27 | 12 |
| 12 | 26895 | 10770 | 3947 | 0.147 | 0.366 | 27/27 | 14 |

per family at 9x (source): 1: R 0.35 P 0.17 (330/954/1902); 2: R 0.27 P 0.16 (315/1157/1977); 3: R 0.27 P 0.18 (162/591/898); 4: R 0.22 P 0.07 (87/400/1327); 5: R 0.08 P 0.00 (47/586/30144); 6: R 0.34 P 0.22 (466/1358/2115); 7: R 0.37 P 0.21 (398/1072/1858); 8: R 0.54 P 0.72 (916/1691/1278); 9: R 0.10 P 0.06 (31/319/500); 10: R 0.36 P 0.16 (288/798/1814); 11: R 0.40 P 0.10 (119/298/1158); 12: R 0.36 P 0.05 (173/484/3675); 13: R 0.45 P 0.21 (398/876/1872); 14: R 0.34 P 0.29 (475/1415/1618); 15: R 0.25 P 0.13 (214/862/1661); 16: R 0.39 P 0.08 (233/601/3002); 17: R 0.39 P 0.14 (244/621/1704); 18: R 0.36 P 0.12 (234/657/1973); 19: R 0.60 P 0.40 (1091/1811/2712); 20: R 0.38 P 0.18 (362/955/2045); 21: R 0.43 P 0.15 (280/647/1857); 22: R 0.46 P 0.23 (575/1244/2471); 23: R 0.24 P 0.64 (683/2852/1065); 24: R 0.56 P 0.20 (609/1095/3084); 25: R 0.32 P 0.18 (359/1129/2007); 26: R 0.29 P 0.12 (159/553/1347); 27: R 0.09 P 0.67 (162/1869/242)

score, days-mode fam (own test days only) (one-to-one +-2 s):

| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |
|---|---|---|---|---|---|---|---|
| 6 | 5021 | 188014 | 3715 | 0.740 | 0.020 | 27/27 | 0 |
| 7 | 5021 | 56521 | 3039 | 0.605 | 0.054 | 27/27 | 2 |
| 8 | 5021 | 22194 | 2382 | 0.474 | 0.107 | 27/27 | 6 |
| 9 | 5021 | 9885 | 1789 | 0.356 | 0.181 | 27/27 | 12 |
| 10 | 5021 | 4684 | 1309 | 0.261 | 0.279 | 27/27 | 16 |
| 11 | 5021 | 2452 | 946 | 0.188 | 0.386 | 27/27 | 22 |
| 12 | 5021 | 1473 | 697 | 0.139 | 0.473 | 27/27 | 21 |

per family at 9x (fam (own test days only)): 1: R 0.34 P 0.27 (67/197/251); 2: R 0.30 P 0.21 (61/203/289); 3: R 0.22 P 0.32 (43/199/133); 4: R 0.20 P 0.14 (26/130/187); 5: R 0.11 P 0.01 (18/167/3536); 6: R 0.37 P 0.44 (71/194/163); 7: R 0.41 P 0.32 (79/195/248); 8: R 0.50 P 0.69 (98/195/143); 9: R 0.12 P 0.09 (13/108/138); 10: R 0.39 P 0.28 (75/194/269); 11: R 0.45 P 0.18 (65/143/359); 12: R 0.36 P 0.09 (59/162/646); 13: R 0.46 P 0.31 (92/200/293); 14: R 0.30 P 0.49 (60/200/122); 15: R 0.23 P 0.27 (47/207/175); 16: R 0.36 P 0.18 (64/176/359); 17: R 0.34 P 0.23 (63/186/271); 18: R 0.35 P 0.19 (66/186/350); 19: R 0.74 P 0.58 (140/189/242); 20: R 0.40 P 0.28 (80/198/282); 21: R 0.45 P 0.37 (89/196/241); 22: R 0.51 P 0.44 (92/181/208); 23: R 0.30 P 0.80 (63/210/79); 24: R 0.59 P 0.31 (113/192/366); 25: R 0.32 P 0.27 (63/195/237); 26: R 0.30 P 0.23 (63/207/274); 27: R 0.09 P 0.79 (19/211/24)

separation: statistic at his 26895 test-day times (best within +-2 s): p10/p50/p90 4.72/7.73/13.20; fraction >= 6..12: 6:0.741, 7:0.596, 8:0.465, 9:0.351, 10:0.265, 11:0.197, 12:0.147; dt median 0.00 s; noise p99.9 median 4.81
<!-- SH_V2_RESULTS_END -->
