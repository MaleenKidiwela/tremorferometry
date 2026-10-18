# Merlin's two deciding tests — T1 all templates, T2 withheld-station template recovery

Written 2026-09-30 04:53Z by `scripts/merlin_tests_report.py` (numbers only; interpretation to be added after review). Queue `scripts/merlin_tests_queue.sh`; outputs `results/merlin_tests/`; log `NOTES.md`.

## Summary (plain language)

T1 asks: if the network matched filter keeps EVERY station template that stacked ≥ 30 array-era events (weighted by its signal-to-noise weight) instead of only the 386 clear ones, do we find more of the known array detections while the scrambled control (stations shifted to different times) stays flat? T2 asks: if we pretend an array station is new, can we rebuild its template for a family only from times the permanent network picks, i.e. without the station?

- T1 NW: observable families 55 → 132; own-family recall over ALL scanned truth (core) 0.074 → 0.119 (≤3 km: 0.078 → 0.123); own precision 0.290 → 0.254 (scrambled control's chance level 0.004); detections per family-day data / control 11.5 / 10.1 → 8.5 / 5.6; families ≥ 2× above the control 33 (control 2) → 280 (control 31).
- T1 SW: observable families 2 → 73; own-family recall over ALL scanned truth (core) 0.003 → 0.010 (≤3 km: 0.003 → 0.010); own precision 0.714 → 0.837 (scrambled control's chance level 0.000); detections per family-day data / control 3.3 / 0.5 → 3.3 / 0.0; families ≥ 2× above the control 2 (control 0) → 73 (control 12).
- T1 Olympic: observable families 2 → 28; own-family recall over ALL scanned truth (core) 0.021 → 0.226 (≤3 km: 0.021 → 0.432); own precision 0.053 → 0.052 (scrambled control's chance level 0.001); detections per family-day data / control 120.4 / 119.2 → 209.4 / 205.6; families ≥ 2× above the control 1 (control 0) → 4 (control 0).
- T1 Arm B: observable families 2 → 41; own-family recall over ALL scanned truth (core) 0.002 → 0.002 (≤3 km: 0.002 → 0.003); own precision 0.571 → 0.224 (scrambled control's chance level 0.002); detections per family-day data / control 0.9 / 0.4 → 3.1 / 0.9; families ≥ 2× above the control 2 (control 0) → 37 (control 3).
- T2 NW (627 family–station pairs, 117 families): recovered-vs-true template median cc 0.27 (round 0) → 0.28 (round 2); cc > 0.7 in 0.11 of pairs; own family best-matching in 0.21 (top 3: 0.38; best is a ≤3-km neighbour: 0.23). Scrambled control: median cc 0.14, own best 0.01. Ceiling (stack at the true test-half detections): median cc 0.67.
- T2 SW (208 family–station pairs, 66 families): recovered-vs-true template median cc 0.15 (round 0) → 0.15 (round 2); cc > 0.7 in 0.00 of pairs; own family best-matching in 0.02 (top 3: 0.08; best is a ≤3-km neighbour: 0.04). Scrambled control: median cc 0.14, own best 0.03. Ceiling (stack at the true test-half detections): median cc 0.64.
- T2 Olympic (45 family–station pairs, 26 families): recovered-vs-true template median cc 0.86 (round 0) → 0.86 (round 2); cc > 0.7 in 0.76 of pairs; own family best-matching in 0.58 (top 3: 0.76; best is a ≤3-km neighbour: 0.76). Scrambled control: median cc 0.12, own best 0.04. Ceiling (stack at the true test-half detections): median cc 0.87.
- T2 Arm B (123 family–station pairs, 35 families): recovered-vs-true template median cc 0.15 (round 0) → 0.14 (round 2); cc > 0.7 in 0.00 of pairs; own family best-matching in 0.00 (top 3: 0.03; best is a ≤3-km neighbour: 0.01). Scrambled control: median cc 0.14, own best 0.01. Ceiling (stack at the true test-half detections): median cc 0.57.
- T2 ALL (1003 family–station pairs, 244 families): recovered-vs-true template median cc 0.21 (round 0) → 0.22 (round 2); cc > 0.7 in 0.10 of pairs; own family best-matching in 0.16 (top 3: 0.29; best is a ≤3-km neighbour: 0.19). Scrambled control: median cc 0.14, own best 0.02. Ceiling (stack at the true test-half detections): median cc 0.66.

## T1 — all templates vs clear only (V2008 + V2010 test half, A2011 all)

Same chain as WA_TRIAL: frozen 1–8 Hz network matched filter, site merge B, ≥ 2 sites per family-day, record ≥ 8.0 MAD, T_ext 8.25 / T_core 8.75, plus per-epoch tiers (SVI FA-only rule on quiet days; A2011 has none → V2010's), HF veto calibrated on V2008/V2010 build-half own matches, frozen cross-family rule on the stations of the templates in use, one-to-one ±2 s own-family scoring. Two site-rolled scrambled controls (control 1 is the chance level; control 2 applies the separation rule to a control). One rule added for BOTH before and after: a station-day with > 5 % of its samples guarded is treated as absent (under the all-stations-valid network rule one half-dead station, e.g. SP2 on 2010-08-20, blanked half of the day for all 274 families once every template is used). The existing trial run (`val`) is scored by the same scorer for reference (1 control → no control separation column).

**Tier `core`, cat** (catalogue = tier & HF veto & cross-family rule; the scrambled control = tier only). `obs fam` = families with ≥ 1 observable scanned episode day; `truth obs / all` = own array detections on observable days / on all scanned days (scored half); recall = on observable days, `recall all` = over ALL scanned truth (so a run that makes more families observable is not penalised). Levels own / ≤3 km / any.

| region | run | obs fam | truth obs / all | recall own · 3km · any | recall all (own) | prec own · 3km · any | control prec own · 3km · any | data / control per fam-day | sep ≥1.5 · 2 · 3 · 5 (data / control) |
|---|---|---|---|---|---|---|---|---|---|
| NW | existing trial run (clear only, V2008+V2010) | 55 | 15901 / 33211 | 0.155 · 0.162 · 0.306 | 0.074 | 0.328 · 0.331 · 0.490 | 0.009 · 0.011 · 0.124 | 11.93 / 10.19 | 32/– · 22/– · 17/– · 8/– |
| NW | before: clear only (rerun) | 55 | 30751 / 66517 | 0.161 · 0.168 · 0.305 | 0.074 | 0.290 · 0.291 · 0.428 | 0.006 · 0.007 · 0.086 | 11.51 / 10.09 | 45/4 · 33/2 · 24/1 · 9/0 |
| NW | after: ALL templates | 132 | 66517 / 66517 | 0.119 · 0.123 · 0.335 | 0.119 | 0.254 · 0.255 · 0.332 | 0.004 · 0.005 · 0.058 | 8.51 / 5.63 | 317/44 · 280/31 · 245/8 · 194/2 |
| SW | existing trial run (clear only, V2008+V2010) | 2 | 364 / 7950 | 0.069 · 0.069 · 0.080 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 · 0.000 · 0.000 | 2.17 / 0.35 | 2/– · 2/– · 2/– · 2/– |
| SW | before: clear only (rerun) | 2 | 361 / 7950 | 0.069 · 0.069 · 0.080 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 · 0.000 · 0.000 | 3.35 / 0.54 | 2/0 · 2/0 · 2/0 · 2/0 |
| SW | after: ALL templates | 73 | 7833 / 7950 | 0.010 · 0.010 · 0.079 | 0.010 | 0.837 · 0.848 · 0.891 | 0.000 · 0.000 · 0.000 | 3.33 / 0.01 | 73/12 · 73/12 · 73/2 · 73/1 |
| Olympic | existing trial run (clear only, V2008+V2010) | 2 | 384 / 7698 | 0.432 · 0.432 · 0.440 | 0.022 | 0.053 · 0.054 · 0.312 | 0.001 · 0.001 · 0.065 | 120.45 / 119.16 | 1/– · 1/– · 1/– · 1/– |
| Olympic | before: clear only (rerun) | 2 | 384 / 7698 | 0.430 · 0.430 · 0.438 | 0.021 | 0.053 · 0.054 · 0.312 | 0.001 · 0.001 · 0.065 | 120.45 / 119.16 | 1/0 · 1/0 · 1/0 · 1/0 |
| Olympic | after: ALL templates | 28 | 7056 / 7698 | 0.247 · 0.471 · 0.788 | 0.226 | 0.052 · 0.089 · 0.170 | 0.001 · 0.017 · 0.070 | 209.45 / 205.63 | 4/1 · 4/0 · 4/0 · 4/0 |
| Arm B | existing trial run (clear only, V2008+V2010) | 2 | 268 / 6844 | 0.045 · 0.045 · 0.075 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 · 0.000 · 0.333 | 0.88 / 0.39 | 2/– · 2/– · 1/– · 0/– |
| Arm B | before: clear only (rerun) | 2 | 268 / 6844 | 0.045 · 0.045 · 0.075 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 · 0.000 · 0.333 | 0.88 / 0.39 | 2/0 · 2/0 · 1/0 · 0/0 |
| Arm B | after: ALL templates | 41 | 6783 / 6844 | 0.003 · 0.003 · 0.134 | 0.002 | 0.224 · 0.224 · 0.355 | 0.002 · 0.008 · 0.075 | 3.13 / 0.92 | 39/3 · 37/3 · 37/0 · 37/0 |
| ALL | existing trial run (clear only, V2008+V2010) | 61 | 16917 / 55703 | 0.158 · 0.165 · 0.301 | 0.048 | 0.250 · 0.252 · 0.439 | 0.007 · 0.009 · 0.111 | 13.95 / 12.24 | 37/– · 27/– · 21/– · 11/– |
| ALL | before: clear only (rerun) | 61 | 31764 / 89009 | 0.162 · 0.169 · 0.302 | 0.058 | 0.254 · 0.256 · 0.411 | 0.006 · 0.007 · 0.084 | 12.80 / 11.39 | 50/4 · 38/2 · 28/1 · 12/0 |
| ALL | after: ALL templates | 274 | 88189 / 89009 | 0.110 · 0.132 · 0.333 | 0.109 | 0.151 · 0.170 · 0.250 | 0.002 · 0.013 · 0.066 | 19.13 / 16.19 | 433/60 · 394/46 · 359/10 · 308/3 |

**Tier `core`, raw** (RAW = tier only, no HF veto / cross-family rule (with every template nearly all families share >= 50 % of their stations, so the frozen cross-family rule becomes a competition among all families); the scrambled control = tier only). `obs fam` = families with ≥ 1 observable scanned episode day; `truth obs / all` = own array detections on observable days / on all scanned days (scored half); recall = on observable days, `recall all` = over ALL scanned truth (so a run that makes more families observable is not penalised). Levels own / ≤3 km / any.

| region | run | obs fam | truth obs / all | recall own · 3km · any | recall all (own) | prec own · 3km · any | control prec own · 3km · any | data / control per fam-day | sep ≥1.5 · 2 · 3 · 5 (data / control) |
|---|---|---|---|---|---|---|---|---|---|
| NW | existing trial run (clear only, V2008+V2010) | 55 | 15901 / 33211 | 0.249 · 0.252 · 0.325 | 0.119 | 0.295 · 0.300 · 0.511 | 0.009 · 0.011 · 0.124 | 11.93 / 10.19 | 32/– · 22/– · 17/– · 8/– |
| NW | before: clear only (rerun) | 55 | 30751 / 66517 | 0.251 · 0.254 · 0.322 | 0.116 | 0.253 · 0.256 · 0.439 | 0.006 · 0.007 · 0.086 | 11.51 / 10.09 | 45/4 · 33/2 · 24/1 · 9/0 |
| NW | after: ALL templates | 132 | 66517 / 66517 | 0.184 · 0.187 · 0.362 | 0.184 | 0.246 · 0.248 · 0.365 | 0.004 · 0.005 · 0.058 | 8.51 / 5.63 | 317/44 · 280/31 · 245/8 · 194/2 |
| SW | existing trial run (clear only, V2008+V2010) | 2 | 364 / 7950 | 0.069 · 0.069 · 0.082 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 · 0.000 · 0.000 | 2.17 / 0.35 | 2/– · 2/– · 2/– · 2/– |
| SW | before: clear only (rerun) | 2 | 361 / 7950 | 0.069 · 0.069 · 0.083 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 · 0.000 · 0.000 | 3.35 / 0.54 | 2/0 · 2/0 · 2/0 · 2/0 |
| SW | after: ALL templates | 73 | 7833 / 7950 | 0.012 · 0.012 · 0.089 | 0.012 | 0.836 · 0.845 · 0.900 | 0.000 · 0.000 · 0.000 | 3.33 / 0.01 | 73/12 · 73/12 · 73/2 · 73/1 |
| Olympic | existing trial run (clear only, V2008+V2010) | 2 | 384 / 7698 | 0.458 · 0.458 · 0.469 | 0.023 | 0.055 · 0.056 · 0.312 | 0.001 · 0.001 · 0.065 | 120.45 / 119.16 | 1/– · 1/– · 1/– · 1/– |
| Olympic | before: clear only (rerun) | 2 | 384 / 7698 | 0.458 · 0.458 · 0.469 | 0.023 | 0.055 · 0.056 · 0.312 | 0.001 · 0.001 · 0.065 | 120.45 / 119.16 | 1/0 · 1/0 · 1/0 · 1/0 |
| Olympic | after: ALL templates | 28 | 7056 / 7698 | 0.550 · 0.805 · 0.907 | 0.504 | 0.049 · 0.145 · 0.271 | 0.001 · 0.017 · 0.070 | 209.45 / 205.63 | 4/1 · 4/0 · 4/0 · 4/0 |
| Arm B | existing trial run (clear only, V2008+V2010) | 2 | 268 / 6844 | 0.045 · 0.045 · 0.078 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 · 0.000 · 0.333 | 0.88 / 0.39 | 2/– · 2/– · 1/– · 0/– |
| Arm B | before: clear only (rerun) | 2 | 268 / 6844 | 0.045 · 0.045 · 0.078 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 · 0.000 · 0.333 | 0.88 / 0.39 | 2/0 · 2/0 · 1/0 · 0/0 |
| Arm B | after: ALL templates | 41 | 6783 / 6844 | 0.006 · 0.008 · 0.161 | 0.006 | 0.084 · 0.094 · 0.429 | 0.002 · 0.008 · 0.075 | 3.13 / 0.92 | 39/3 · 37/3 · 37/0 · 37/0 |
| ALL | existing trial run (clear only, V2008+V2010) | 61 | 16917 / 55703 | 0.246 · 0.249 · 0.319 | 0.075 | 0.250 · 0.254 · 0.474 | 0.007 · 0.009 · 0.111 | 13.95 / 12.24 | 37/– · 27/– · 21/– · 11/– |
| ALL | before: clear only (rerun) | 61 | 31764 / 89009 | 0.250 · 0.253 · 0.319 | 0.089 | 0.235 · 0.238 · 0.428 | 0.006 · 0.007 · 0.084 | 12.80 / 11.39 | 50/4 · 38/2 · 28/1 · 12/0 |
| ALL | after: ALL templates | 274 | 88189 / 89009 | 0.184 · 0.207 · 0.366 | 0.183 | 0.126 · 0.185 · 0.308 | 0.002 · 0.013 · 0.066 | 19.13 / 16.19 | 433/60 · 394/46 · 359/10 · 308/3 |

**Tier `core_ep`, cat** (catalogue = tier & HF veto & cross-family rule; the scrambled control = tier only). `obs fam` = families with ≥ 1 observable scanned episode day; `truth obs / all` = own array detections on observable days / on all scanned days (scored half); recall = on observable days, `recall all` = over ALL scanned truth (so a run that makes more families observable is not penalised). Levels own / ≤3 km / any.

| region | run | obs fam | truth obs / all | recall own · 3km · any | recall all (own) | prec own · 3km · any | control prec own · 3km · any | data / control per fam-day | sep ≥1.5 · 2 · 3 · 5 (data / control) |
|---|---|---|---|---|---|---|---|---|---|
| NW | existing trial run (clear only, V2008+V2010) | 55 | 15901 / 33211 | 0.002 · 0.002 · 0.005 | 0.001 | 0.917 · 0.917 · 0.917 | 0.000 · 0.000 · 0.250 | 0.04 / 0.01 | 30/– · 24/– · 8/– · 2/– |
| NW | before: clear only (rerun) | 55 | 30751 / 66517 | 0.002 · 0.003 · 0.005 | 0.001 | 0.922 · 0.922 · 0.935 | 0.000 · 0.000 · 0.077 | 0.04 / 0.02 | 42/8 · 32/1 · 10/0 · 3/0 |
| NW | after: ALL templates | 132 | 66517 / 66517 | 0.000 · 0.000 · 0.002 | 0.000 | 0.952 · 0.952 · 0.952 | – · – · – | 0.01 / 0.00 | 56/0 · 56/0 · 16/0 · 3/0 |
| SW | existing trial run (clear only, V2008+V2010) | 2 | 364 / 7950 | 0.000 · 0.000 · 0.000 | 0.000 | – · – · – | – · – · – | 0.00 / 0.00 | 0/– · 0/– · 0/– · 0/– |
| SW | before: clear only (rerun) | 2 | 361 / 7950 | 0.000 · 0.000 · 0.000 | 0.000 | – · – · – | – · – · – | 0.00 / 0.00 | 0/0 · 0/0 · 0/0 · 0/0 |
| SW | after: ALL templates | 73 | 7833 / 7950 | 0.000 · 0.000 · 0.000 | 0.000 | – · – · – | – · – · – | 0.00 / 0.00 | 5/0 · 5/0 · 0/0 · 0/0 |
| Olympic | existing trial run (clear only, V2008+V2010) | 2 | 384 / 7698 | 0.029 · 0.029 · 0.029 | 0.001 | 0.145 · 0.145 · 0.816 | 0.000 · 0.000 · 0.039 | 2.95 / 2.89 | 0/– · 0/– · 0/– · 0/– |
| Olympic | before: clear only (rerun) | 2 | 384 / 7698 | 0.029 · 0.029 · 0.029 | 0.001 | 0.145 · 0.145 · 0.816 | 0.000 · 0.000 · 0.039 | 2.95 / 2.89 | 0/0 · 0/0 · 0/0 · 0/0 |
| Olympic | after: ALL templates | 28 | 7056 / 7698 | 0.020 · 0.035 · 0.049 | 0.018 | 0.409 · 0.540 · 0.638 | 0.003 · 0.029 · 0.072 | 1.29 / 0.98 | 16/2 · 14/1 · 7/1 · 2/0 |
| Arm B | existing trial run (clear only, V2008+V2010) | 2 | 268 / 6844 | 0.000 · 0.000 · 0.004 | 0.000 | 0.000 · 0.000 · 1.000 | – · – · – | 0.02 / 0.00 | 1/– · 1/– · 0/– · 0/– |
| Arm B | before: clear only (rerun) | 2 | 268 / 6844 | 0.000 · 0.000 · 0.004 | 0.000 | 0.000 · 0.000 · 1.000 | – · – · – | 0.02 / 0.00 | 1/0 · 1/0 · 0/0 · 0/0 |
| Arm B | after: ALL templates | 41 | 6783 / 6844 | 0.000 · 0.000 · 0.002 | 0.000 | – · – · – | – · – · – | 0.02 / 0.00 | 9/0 · 9/0 · 5/0 · 1/0 |
| ALL | existing trial run (clear only, V2008+V2010) | 61 | 16917 / 55703 | 0.003 · 0.003 · 0.005 | 0.001 | 0.389 · 0.389 · 0.850 | 0.000 · 0.000 · 0.083 | 0.10 / 0.08 | 31/– · 25/– · 8/– · 2/– |
| ALL | before: clear only (rerun) | 61 | 31764 / 89009 | 0.003 · 0.003 · 0.006 | 0.001 | 0.532 · 0.532 · 0.877 | 0.000 · 0.000 · 0.057 | 0.08 / 0.06 | 43/8 · 33/1 · 10/0 · 3/0 |
| ALL | after: ALL templates | 274 | 88189 / 89009 | 0.002 · 0.003 · 0.006 | 0.002 | 0.441 · 0.564 · 0.656 | 0.003 · 0.029 · 0.072 | 0.08 / 0.06 | 86/2 · 84/1 · 28/1 · 6/0 |

**Tier `ext`, cat** (catalogue = tier & HF veto & cross-family rule; the scrambled control = tier only). `obs fam` = families with ≥ 1 observable scanned episode day; `truth obs / all` = own array detections on observable days / on all scanned days (scored half); recall = on observable days, `recall all` = over ALL scanned truth (so a run that makes more families observable is not penalised). Levels own / ≤3 km / any.

| region | run | obs fam | truth obs / all | recall own · 3km · any | recall all (own) | prec own · 3km · any | control prec own · 3km · any | data / control per fam-day | sep ≥1.5 · 2 · 3 · 5 (data / control) |
|---|---|---|---|---|---|---|---|---|---|
| NW | existing trial run (clear only, V2008+V2010) | 55 | 15901 / 33211 | 0.186 · 0.195 · 0.369 | 0.089 | 0.247 · 0.250 · 0.405 | 0.007 · 0.009 · 0.116 | 19.11 / 17.00 | 27/– · 20/– · 13/– · 1/– |
| NW | before: clear only (rerun) | 55 | 30751 / 66517 | 0.191 · 0.200 · 0.366 | 0.088 | 0.214 · 0.216 · 0.346 | 0.005 · 0.006 · 0.082 | 18.73 / 16.92 | 40/2 · 29/1 · 18/1 · 2/0 |
| NW | after: ALL templates | 132 | 66517 / 66517 | 0.141 · 0.146 · 0.393 | 0.141 | 0.220 · 0.221 · 0.299 | 0.004 · 0.004 · 0.059 | 11.64 / 8.31 | 288/41 · 255/23 · 206/7 · 137/1 |
| SW | existing trial run (clear only, V2008+V2010) | 2 | 364 / 7950 | 0.080 · 0.080 · 0.099 | 0.004 | 0.659 · 0.659 · 0.682 | 0.000 · 0.000 · 0.000 | 2.75 / 0.75 | 2/– · 2/– · 2/– · 0/– |
| SW | before: clear only (rerun) | 2 | 361 / 7950 | 0.080 · 0.080 · 0.100 | 0.004 | 0.659 · 0.659 · 0.682 | 0.000 · 0.000 · 0.000 | 4.23 / 1.15 | 2/0 · 2/0 · 2/0 · 0/0 |
| SW | after: ALL templates | 73 | 7833 / 7950 | 0.013 · 0.014 · 0.102 | 0.013 | 0.667 · 0.679 · 0.724 | 0.000 · 0.000 · 0.000 | 3.60 / 0.07 | 73/23 · 73/21 · 73/5 · 73/1 |
| Olympic | existing trial run (clear only, V2008+V2010) | 2 | 384 / 7698 | 0.443 · 0.443 · 0.456 | 0.022 | 0.047 · 0.048 · 0.285 | 0.001 · 0.001 · 0.065 | 140.46 / 138.64 | 1/– · 1/– · 1/– · 1/– |
| Olympic | before: clear only (rerun) | 2 | 384 / 7698 | 0.440 · 0.440 · 0.453 | 0.022 | 0.047 · 0.048 · 0.285 | 0.001 · 0.001 · 0.065 | 140.46 / 138.64 | 1/0 · 1/0 · 1/0 · 1/0 |
| Olympic | after: ALL templates | 28 | 7056 / 7698 | 0.249 · 0.477 · 0.800 | 0.228 | 0.046 · 0.079 · 0.156 | 0.001 · 0.017 · 0.069 | 246.84 / 242.82 | 4/0 · 4/0 · 4/0 · 3/0 |
| Arm B | existing trial run (clear only, V2008+V2010) | 2 | 268 / 6844 | 0.056 · 0.056 · 0.086 | 0.002 | 0.395 · 0.395 · 0.553 | 0.000 · 0.000 · 0.143 | 1.45 / 0.98 | 1/– · 0/– · 0/– · 0/– |
| Arm B | before: clear only (rerun) | 2 | 268 / 6844 | 0.056 · 0.056 · 0.086 | 0.002 | 0.395 · 0.395 · 0.553 | 0.000 · 0.000 · 0.143 | 1.45 / 0.98 | 1/0 · 0/0 · 0/0 · 0/0 |
| Arm B | after: ALL templates | 41 | 6783 / 6844 | 0.004 · 0.005 · 0.155 | 0.004 | 0.139 · 0.139 · 0.243 | 0.003 · 0.007 · 0.074 | 4.14 / 1.74 | 38/6 · 37/5 · 37/3 · 36/0 |
| ALL | existing trial run (clear only, V2008+V2010) | 61 | 16917 / 55703 | 0.188 · 0.195 · 0.360 | 0.057 | 0.202 · 0.205 · 0.378 | 0.006 · 0.007 · 0.107 | 21.15 / 19.09 | 31/– · 23/– · 16/– · 2/– |
| ALL | before: clear only (rerun) | 61 | 31764 / 89009 | 0.192 · 0.200 · 0.362 | 0.068 | 0.195 · 0.197 · 0.339 | 0.005 · 0.006 · 0.080 | 20.07 / 18.26 | 44/2 · 32/1 · 21/1 · 3/0 |
| ALL | after: ALL templates | 274 | 88189 / 89009 | 0.128 · 0.150 · 0.382 | 0.127 | 0.139 · 0.155 · 0.233 | 0.002 · 0.012 · 0.065 | 23.68 / 20.36 | 403/70 · 369/49 · 320/15 · 249/2 |

Per window (tier core, catalogue):

| window | region | run | obs fam | recall own (obs) | recall all (own) | prec own · 3km · any | control prec own | data / ctl per fam-day | sep≥2 data / ctl |
|---|---|---|---|---|---|---|---|---|---|
| V2008 | NW | val | 55 | 0.157 | 0.075 | 0.343 · 0.345 · 0.502 | 0.007 | 11.08 / 9.59 | 10 / – |
| V2008 | NW | mt1_before | 55 | 0.156 | 0.075 | 0.343 · 0.345 · 0.502 | 0.007 | 11.10 / 9.61 | 10 / 1 |
| V2008 | NW | mt1_all | 132 | 0.108 | 0.108 | 0.342 · 0.343 · 0.420 | 0.005 | 8.08 / 4.35 | 112 / 11 |
| V2008 | SW | val | 2 | 0.069 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 | 2.17 / 0.35 | 2 / – |
| V2008 | SW | mt1_before | 2 | 0.069 | 0.003 | 0.714 · 0.714 · 0.743 | 0.000 | 3.35 / 0.54 | 2 / 0 |
| V2008 | SW | mt1_all | 73 | 0.010 | 0.010 | 0.837 · 0.848 · 0.891 | 0.000 | 3.33 / 0.01 | 73 / 12 |
| V2008 | ALL | val | 57 | 0.153 | 0.051 | 0.347 · 0.349 · 0.505 | 0.007 | 10.81 / 9.30 | 12 / – |
| V2008 | ALL | mt1_before | 57 | 0.153 | 0.051 | 0.346 · 0.348 · 0.504 | 0.006 | 10.94 / 9.43 | 12 / 1 |
| V2008 | ALL | mt1_all | 205 | 0.076 | 0.075 | 0.351 · 0.352 · 0.428 | 0.005 | 6.25 / 2.68 | 185 / 23 |
| V2010 | NW | val | 55 | 0.154 | 0.074 | 0.315 · 0.318 · 0.479 | 0.012 | 12.92 / 10.88 | 12 / – |
| V2010 | NW | mt1_before | 55 | 0.157 | 0.076 | 0.320 · 0.322 · 0.481 | 0.012 | 13.05 / 10.94 | 12 / 1 |
| V2010 | NW | mt1_all | 132 | 0.104 | 0.104 | 0.372 · 0.373 · 0.485 | 0.009 | 8.35 / 3.91 | 117 / 10 |
| V2010 | Olympic | val | 2 | 0.432 | 0.022 | 0.053 · 0.054 · 0.312 | 0.001 | 120.45 / 119.16 | 1 / – |
| V2010 | Olympic | mt1_before | 2 | 0.430 | 0.021 | 0.053 · 0.054 · 0.312 | 0.001 | 120.45 / 119.16 | 1 / 0 |
| V2010 | Olympic | mt1_all | 28 | 0.247 | 0.226 | 0.052 · 0.089 · 0.170 | 0.001 | 209.45 / 205.63 | 4 / 0 |
| V2010 | Arm B | val | 2 | 0.045 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 | 0.88 / 0.39 | 2 / – |
| V2010 | Arm B | mt1_before | 2 | 0.045 | 0.002 | 0.571 · 0.571 · 0.762 | 0.000 | 0.88 / 0.39 | 2 / 0 |
| V2010 | Arm B | mt1_all | 41 | 0.003 | 0.002 | 0.224 · 0.224 · 0.355 | 0.002 | 3.13 / 0.92 | 37 / 3 |
| V2010 | ALL | val | 59 | 0.162 | 0.045 | 0.202 · 0.204 · 0.407 | 0.008 | 17.33 / 15.40 | 15 / – |
| V2010 | ALL | mt1_before | 59 | 0.166 | 0.046 | 0.205 · 0.207 · 0.408 | 0.008 | 17.45 / 15.45 | 15 / 1 |
| V2010 | ALL | mt1_all | 201 | 0.114 | 0.112 | 0.093 · 0.125 · 0.210 | 0.002 | 41.53 / 37.75 | 158 / 13 |
| A2011 | NW | mt1_before | 55 | 0.165 | 0.074 | 0.257 · 0.259 · 0.377 | 0.004 | 10.77 / 9.89 | 11 / 0 |
| A2011 | NW | mt1_all | 132 | 0.131 | 0.131 | 0.206 · 0.207 · 0.277 | 0.003 | 8.97 / 7.84 | 51 / 10 |
| A2011 | ALL | mt1_before | 55 | 0.165 | 0.074 | 0.257 · 0.259 · 0.377 | 0.004 | 10.77 / 9.89 | 11 / 0 |
| A2011 | ALL | mt1_all | 132 | 0.131 | 0.131 | 0.206 · 0.207 · 0.277 | 0.003 | 8.97 / 7.84 | 51 / 10 |

Per-epoch tiers `val`: V2008 T_ext_ep 16.00 / T_core_ep 20.75 (FA at 8.75 = 8.51); V2010 T_ext_ep 13.50 / T_core_ep 17.75 (FA at 8.75 = 12.45)

Per-epoch tiers `mt1_before`: A2011 T_ext_ep 13.50 / T_core_ep 17.75 (FA at 8.75 = –); V2008 T_ext_ep 16.00 / T_core_ep 20.75 (FA at 8.75 = 9.16); V2010 T_ext_ep 13.50 / T_core_ep 17.75 (FA at 8.75 = 12.45)

Per-epoch tiers `mt1_all`: A2011 T_ext_ep 16.75 / T_core_ep 21.25 (FA at 8.75 = –); V2008 T_ext_ep 18.50 / T_core_ep 22.25 (FA at 8.75 = 23.93); V2010 T_ext_ep 16.75 / T_core_ep 21.25 (FA at 8.75 = 26.95)

Per family: `results/merlin_tests/t1/score_<run>/by_family.csv` (every family-window, every column).

## T2 — withheld-station template recovery (array era: V2008 XU, V2010 XU + XG)

Pairs = (family, window, XU/XG station with a clear truth template from the family's BUILD-half array detections, ≥ 30 events), ≤ 4 stations per family-window drawn at random. Round 0: candidates = top-K (2 per scored hour, ±4 s declustered) peaks of the family's PERMANENT-network matched filter (all ≥ 30-event templates; identical to stack_activity's candidates), TEST-half origin times only; the withheld station is stacked there and a template is cut by the template rule. Rounds 1–2: the station joins the network with its recovered template (its own site, weight n_ch × SNR weight) and candidates are re-picked. cc = 2–8 Hz (designed at fs = 20 Hz) normalised cc between the recovered stack and the TRUE template (build half), best lag within ±2 s; rank among all families' truth templates at that station (any lag). Scrambled control = the same procedure on the site-rolled network. Ceiling = the station stacked at the family's own TEST-half array detections. `perm_clear_sites` = median number of permanent sites with a CLEAR template for the family (the carrying network).

| group | value | series | round | pairs | families | median cc (±2 s) | cc > 0.7 | cc > 0.5 | own best | top 3 | best ≤3 km | median rank / competitors | candidate purity (own) | windows stacked |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| region | ALL | ceiling | -1 | 1003 | 244 | 0.659 | 0.38 | 0.82 | 0.68 | 0.92 | 0.72 | 1.0 / 200 | – | 119 |
| region | ALL | control | 0 | 1003 | 244 | 0.139 | 0.00 | 0.00 | 0.01 | 0.02 | 0.03 | 85.0 / 200 | 0.000 | 771 |
| region | ALL | control | 1 | 1003 | 244 | 0.139 | 0.00 | 0.00 | 0.01 | 0.02 | 0.03 | 78.0 / 200 | 0.000 | 771 |
| region | ALL | control | 2 | 1003 | 244 | 0.139 | 0.00 | 0.00 | 0.02 | 0.04 | 0.03 | 80.0 / 200 | 0.000 | 771 |
| region | ALL | data | 0 | 1003 | 244 | 0.207 | 0.09 | 0.17 | 0.20 | 0.32 | 0.22 | 21.0 / 200 | 0.032 | 773 |
| region | ALL | data | 1 | 1003 | 244 | 0.218 | 0.10 | 0.17 | 0.17 | 0.30 | 0.19 | 18.0 / 200 | 0.032 | 773 |
| region | ALL | data | 2 | 1003 | 244 | 0.216 | 0.10 | 0.17 | 0.16 | 0.29 | 0.19 | 20.0 / 200 | 0.031 | 773 |
| region | Olympic_2010 | ceiling | -1 | 45 | 26 | 0.873 | 0.73 | 0.91 | 0.89 | 0.96 | 0.96 | 1.0 / 200 | – | 100 |
| region | Olympic_2010 | control | 0 | 45 | 26 | 0.130 | 0.00 | 0.00 | 0.02 | 0.09 | 0.16 | 93.0 / 200 | 0.000 | 816 |
| region | Olympic_2010 | control | 1 | 45 | 26 | 0.126 | 0.00 | 0.00 | 0.04 | 0.09 | 0.22 | 103.0 / 200 | 0.000 | 816 |
| region | Olympic_2010 | control | 2 | 45 | 26 | 0.122 | 0.00 | 0.00 | 0.04 | 0.07 | 0.27 | 107.0 / 200 | 0.000 | 816 |
| region | Olympic_2010 | data | 0 | 45 | 26 | 0.856 | 0.76 | 0.80 | 0.69 | 0.80 | 0.80 | 1.0 / 200 | 0.055 | 816 |
| region | Olympic_2010 | data | 1 | 45 | 26 | 0.862 | 0.76 | 0.80 | 0.62 | 0.78 | 0.78 | 1.0 / 200 | 0.055 | 816 |
| region | Olympic_2010 | data | 2 | 45 | 26 | 0.864 | 0.76 | 0.80 | 0.58 | 0.76 | 0.76 | 1.0 / 200 | 0.054 | 816 |
| region | Olympic_2010_armB | ceiling | -1 | 123 | 35 | 0.568 | 0.16 | 0.67 | 0.75 | 0.89 | 0.81 | 1.0 / 200 | – | 78 |
| region | Olympic_2010_armB | control | 0 | 123 | 35 | 0.138 | 0.00 | 0.00 | 0.01 | 0.02 | 0.07 | 73.0 / 200 | 0.000 | 816 |
| region | Olympic_2010_armB | control | 1 | 123 | 35 | 0.137 | 0.00 | 0.00 | 0.02 | 0.03 | 0.05 | 66.0 / 200 | 0.000 | 816 |
| region | Olympic_2010_armB | control | 2 | 123 | 35 | 0.138 | 0.00 | 0.00 | 0.01 | 0.05 | 0.02 | 59.0 / 200 | 0.000 | 816 |
| region | Olympic_2010_armB | data | 0 | 123 | 35 | 0.147 | 0.00 | 0.00 | 0.01 | 0.02 | 0.03 | 64.0 / 200 | 0.004 | 816 |
| region | Olympic_2010_armB | data | 1 | 123 | 35 | 0.151 | 0.00 | 0.00 | 0.01 | 0.02 | 0.02 | 63.0 / 200 | 0.004 | 816 |
| region | Olympic_2010_armB | data | 2 | 123 | 35 | 0.142 | 0.00 | 0.00 | 0.00 | 0.03 | 0.01 | 64.0 / 200 | 0.002 | 816 |
| region | Royer_NW_WA | ceiling | -1 | 627 | 117 | 0.675 | 0.42 | 0.83 | 0.61 | 0.90 | 0.63 | 1.0 / 200 | – | 129 |
| region | Royer_NW_WA | control | 0 | 627 | 117 | 0.138 | 0.00 | 0.00 | 0.01 | 0.02 | 0.01 | 85.0 / 200 | 0.000 | 816 |
| region | Royer_NW_WA | control | 1 | 627 | 117 | 0.139 | 0.00 | 0.00 | 0.01 | 0.02 | 0.01 | 79.0 / 200 | 0.000 | 816 |
| region | Royer_NW_WA | control | 2 | 627 | 117 | 0.140 | 0.00 | 0.00 | 0.01 | 0.03 | 0.02 | 82.0 / 200 | 0.000 | 816 |
| region | Royer_NW_WA | data | 0 | 627 | 117 | 0.268 | 0.10 | 0.21 | 0.25 | 0.43 | 0.27 | 6.0 / 200 | 0.048 | 816 |
| region | Royer_NW_WA | data | 1 | 627 | 117 | 0.286 | 0.11 | 0.21 | 0.22 | 0.40 | 0.24 | 6.0 / 200 | 0.048 | 816 |
| region | Royer_NW_WA | data | 2 | 627 | 117 | 0.277 | 0.11 | 0.22 | 0.21 | 0.38 | 0.23 | 7.0 / 200 | 0.045 | 816 |
| region | Royer_SW_WA | ceiling | -1 | 208 | 66 | 0.644 | 0.31 | 0.84 | 0.82 | 0.97 | 0.88 | 1.0 / 201 | – | 109 |
| region | Royer_SW_WA | control | 0 | 208 | 66 | 0.141 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 89.5 / 201 | 0.000 | 686 |
| region | Royer_SW_WA | control | 1 | 208 | 66 | 0.144 | 0.00 | 0.00 | 0.02 | 0.02 | 0.02 | 83.5 / 201 | 0.000 | 685 |
| region | Royer_SW_WA | control | 2 | 208 | 66 | 0.138 | 0.00 | 0.00 | 0.03 | 0.04 | 0.04 | 88.0 / 201 | 0.000 | 686 |
| region | Royer_SW_WA | data | 0 | 208 | 66 | 0.148 | 0.00 | 0.00 | 0.03 | 0.06 | 0.05 | 90.5 / 201 | 0.011 | 685 |
| region | Royer_SW_WA | data | 1 | 208 | 66 | 0.157 | 0.00 | 0.00 | 0.02 | 0.06 | 0.04 | 73.5 / 201 | 0.010 | 684 |
| region | Royer_SW_WA | data | 2 | 208 | 66 | 0.150 | 0.00 | 0.00 | 0.02 | 0.08 | 0.04 | 74.5 / 201 | 0.011 | 685 |

By number of permanent sites carrying the family (clear templates):

| group | value | series | round | pairs | families | median cc (±2 s) | cc > 0.7 | cc > 0.5 | own best | top 3 | best ≤3 km | median rank / competitors | candidate purity (own) | windows stacked |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| perm_clear_sites | 0 | ceiling | -1 | 198 | 72 | 0.599 | 0.17 | 0.69 | 0.80 | 0.95 | 0.82 | 1.0 / 200 | – | 93 |
| perm_clear_sites | 0 | control | 0 | 198 | 72 | 0.143 | 0.00 | 0.00 | 0.01 | 0.01 | 0.01 | 87.0 / 200 | 0.000 | 696 |
| perm_clear_sites | 0 | control | 1 | 198 | 72 | 0.137 | 0.00 | 0.00 | 0.01 | 0.02 | 0.01 | 78.0 / 200 | 0.000 | 696 |
| perm_clear_sites | 0 | control | 2 | 198 | 72 | 0.134 | 0.00 | 0.00 | 0.02 | 0.02 | 0.02 | 77.0 / 200 | 0.000 | 696 |
| perm_clear_sites | 0 | data | 0 | 198 | 72 | 0.147 | 0.02 | 0.03 | 0.03 | 0.05 | 0.04 | 80.0 / 200 | 0.006 | 696 |
| perm_clear_sites | 0 | data | 1 | 198 | 72 | 0.156 | 0.02 | 0.03 | 0.02 | 0.06 | 0.03 | 74.0 / 200 | 0.006 | 696 |
| perm_clear_sites | 0 | data | 2 | 198 | 72 | 0.147 | 0.02 | 0.03 | 0.03 | 0.07 | 0.04 | 74.0 / 200 | 0.005 | 696 |
| perm_clear_sites | 1 | ceiling | -1 | 435 | 112 | 0.663 | 0.43 | 0.81 | 0.69 | 0.91 | 0.74 | 1.0 / 200 | – | 114 |
| perm_clear_sites | 1 | control | 0 | 435 | 112 | 0.140 | 0.00 | 0.00 | 0.01 | 0.03 | 0.04 | 82.0 / 200 | 0.000 | 816 |
| perm_clear_sites | 1 | control | 1 | 435 | 112 | 0.142 | 0.00 | 0.00 | 0.02 | 0.03 | 0.05 | 75.0 / 200 | 0.000 | 816 |
| perm_clear_sites | 1 | control | 2 | 435 | 112 | 0.144 | 0.00 | 0.00 | 0.02 | 0.05 | 0.05 | 81.0 / 200 | 0.000 | 816 |
| perm_clear_sites | 1 | data | 0 | 435 | 112 | 0.203 | 0.13 | 0.21 | 0.22 | 0.33 | 0.25 | 24.0 / 200 | 0.028 | 816 |
| perm_clear_sites | 1 | data | 1 | 435 | 112 | 0.209 | 0.14 | 0.22 | 0.20 | 0.31 | 0.23 | 21.0 / 200 | 0.029 | 816 |
| perm_clear_sites | 1 | data | 2 | 435 | 112 | 0.207 | 0.15 | 0.21 | 0.19 | 0.31 | 0.22 | 23.0 / 200 | 0.027 | 816 |
| perm_clear_sites | 2 | ceiling | -1 | 196 | 35 | 0.683 | 0.41 | 0.86 | 0.56 | 0.87 | 0.59 | 1.0 / 200 | – | 148 |
| perm_clear_sites | 2 | control | 0 | 196 | 35 | 0.135 | 0.00 | 0.00 | 0.01 | 0.01 | 0.01 | 85.5 / 200 | 0.000 | 816 |
| perm_clear_sites | 2 | control | 1 | 196 | 35 | 0.138 | 0.00 | 0.00 | 0.01 | 0.02 | 0.02 | 84.0 / 200 | 0.000 | 816 |
| perm_clear_sites | 2 | control | 2 | 196 | 35 | 0.138 | 0.00 | 0.00 | 0.01 | 0.03 | 0.03 | 82.0 / 200 | 0.000 | 816 |
| perm_clear_sites | 2 | data | 0 | 196 | 35 | 0.243 | 0.13 | 0.22 | 0.16 | 0.32 | 0.20 | 10.0 / 200 | 0.049 | 816 |
| perm_clear_sites | 2 | data | 1 | 196 | 35 | 0.261 | 0.14 | 0.22 | 0.14 | 0.32 | 0.17 | 8.0 / 200 | 0.051 | 815 |
| perm_clear_sites | 2 | data | 2 | 196 | 35 | 0.259 | 0.14 | 0.22 | 0.13 | 0.31 | 0.16 | 7.5 / 200 | 0.051 | 816 |
| perm_clear_sites | >=3 | ceiling | -1 | 174 | 25 | 0.698 | 0.48 | 0.93 | 0.68 | 0.96 | 0.68 | 1.0 / 200 | – | 137 |
| perm_clear_sites | >=3 | control | 0 | 174 | 25 | 0.137 | 0.00 | 0.00 | 0.01 | 0.03 | 0.02 | 85.0 / 200 | 0.000 | 816 |
| perm_clear_sites | >=3 | control | 1 | 174 | 25 | 0.138 | 0.00 | 0.00 | 0.00 | 0.03 | 0.01 | 81.5 / 200 | 0.000 | 816 |
| perm_clear_sites | >=3 | control | 2 | 174 | 25 | 0.145 | 0.00 | 0.00 | 0.01 | 0.05 | 0.01 | 80.5 / 200 | 0.000 | 816 |
| perm_clear_sites | >=3 | data | 0 | 174 | 25 | 0.327 | 0.04 | 0.15 | 0.37 | 0.59 | 0.39 | 2.5 / 200 | 0.067 | 816 |
| perm_clear_sites | >=3 | data | 1 | 174 | 25 | 0.340 | 0.05 | 0.16 | 0.29 | 0.52 | 0.32 | 3.0 / 200 | 0.069 | 816 |
| perm_clear_sites | >=3 | data | 2 | 174 | 25 | 0.329 | 0.05 | 0.18 | 0.28 | 0.47 | 0.31 | 4.0 / 200 | 0.067 | 816 |

Per pair and round: `results/merlin_tests/t2/pair_results.csv`; candidates `round<r>_cand.parquet`; recovered stacks/templates `round<r>.npz` + `round<r>_qc.csv`; truth `truth_qc.csv`, `truth_stacks.npz`.
