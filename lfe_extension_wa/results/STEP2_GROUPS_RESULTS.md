# STEP 2 — detect ≤ 3-km GROUPS of families, then attribute which member fired

Written 2026-09-30 08:33Z by `scripts/step2_report.py` (numbers; interpretation after review). Queue `scripts/step2_groups_queue.sh`; outputs `results/step2_groups/`; figure `figures/step2_groups.png`.

## Summary (plain language)

Families closer than 3 km cannot be told apart reliably by location, and T1 showed that their permanent-station detections overlap. Step 2 treats each such cluster as one GROUP: one template per station built from all members' array-era events, one matched filter per group, scored against the union of the members' array detections. Then, for every group detection, it asks which member fired, using only the stations of the 2–3 nearest sites, as a probability. The scrambled control (each site shifted to a different time of the same day) is the same one T1 used.

- **Groups** (complete linkage, every pair ≤ 3 km): 215 groups from 329 families; sizes 1: 150, 2: 49, 3: 7, 4: 4, 5: 3, 6: 1, 23: 1. 65 multi-member groups hold 179 families. The 23-member group is 21 UNLOCATED Olympic Arm A families whose approximate position is the B013/B943 borehole pair, plus 2 NW families there — a group by position rule, not by evidence.
- **Group templates**: 7587 ≥ 30-event templates for 181 groups; in multi-member groups the median template stacks 335 events vs 207 for the best single member (median gain ×1.54); clear templates 253. Member alignment reliable for 41 of 114 non-reference members (co-detection prior or stack-cc peak ≥ 1.2× the next peak); the others' origin shift is uncertain (their stacks are noise-like on permanent stations).
- **Detection ALL** (core tier, raw): observable groups 180 (T1 summed) → 181 (group MF); recall over all scanned truth 0.181 → 0.138; precision 0.093 → 0.150 (scrambled control's chance 0.003 / 0.003); per group-day data/control 26.47/22.67 → 13.37/10.73; group-windows ≥ 3× above the control 256 (control 11) → 255 (control 7).
- **Detection Royer_NW_WA** (core tier, raw): observable groups 104 (T1 summed) → 104 (group MF); recall over all scanned truth 0.182 → 0.152; precision 0.239 → 0.220 (scrambled control's chance 0.005 / 0.004); per group-day data/control 10.49/6.99 → 9.20/6.53; group-windows ≥ 3× above the control 188 (control 9) → 187 (control 6).
- **Detection Royer_SW_WA** (core tier, raw): observable groups 52 (T1 summed) → 52 (group MF); recall over all scanned truth 0.012 → 0.007; precision 0.837 → 0.810 (scrambled control's chance 0.000 / 0.000); per group-day data/control 4.14/0.01 → 2.89/0.01; group-windows ≥ 3× above the control 52 (control 2) → 52 (control 1).
- **Detection Olympic_2010** (core tier, raw): observable groups 8 (T1 summed) → 9 (group MF); recall over all scanned truth 0.523 → 0.276; precision 0.030 → 0.055 (scrambled control's chance 0.003 / 0.002); per group-day data/control 538.52/531.18 → 163.84/161.11; group-windows ≥ 3× above the control 2 (control 0) → 3 (control 0).
- **Detection Olympic_2010_armB** (core tier, raw): observable groups 16 (T1 summed) → 16 (group MF); recall over all scanned truth 0.007 → 0.002; precision 0.084 → 0.200 (scrambled control's chance 0.008 / 0.000); per group-day data/control 7.45/2.34 → 1.76/0.33; group-windows ≥ 3× above the control 14 (control 0) → 13 (control 0).
- **Attribution V2010_test** (1231 group detections matched to one member's array detection): the top member is the true one 0.767 [0.74–0.79] vs chance 0.371 and vs 'always the most active member' 0.768; mean P(true) 0.737 (uniform 0.371); log-loss skill vs uniform 0.439; with the activity prior added top-1 0.807, skill 0.488.
- **Attribution V2008_test** (488 group detections matched to one member's array detection): the top member is the true one 0.887 [0.86–0.91] vs chance 0.490 and vs 'always the most active member' 0.867; mean P(true) 0.850 (uniform 0.490); log-loss skill vs uniform 0.600; with the activity prior added top-1 0.928, skill 0.693.
- **Attribution A2011** (1319 group detections matched to one member's array detection): the top member is the true one 0.931 [0.92–0.94] vs chance 0.492 and vs 'always the most active member' 0.888; mean P(true) 0.900 (uniform 0.492); log-loss skill vs uniform 0.721; with the activity prior added top-1 0.941, skill 0.745.
- Stop-signal check (pooled test sets, n = 3038): attribution is NOT uniform — top-1 0.857 vs chance 0.443 and vs most-active 0.836, log-loss skill 0.587 → member identity is carried by the nearest permanent sites for these (array-era) events.
- **Family-level expected count** (all families, 13334 episode family-days): Spearman with the true array count Σ P(real)·P(member) 0.45 (uniform split 0.45; step-1 family ΣP 0.51); ΣE 16954 vs truth 88189 (step-1 ΣP 16991).
- **Family-level expected count** (members of multi-member groups, 7120 episode family-days): Spearman with the true array count Σ P(real)·P(member) 0.41 (uniform split 0.41; step-1 family ΣP 0.51); ΣE 4884 vs truth 51493 (step-1 ΣP 9454).

Caveats: the array catalogue is the truth but incomplete; merged co-detections (the same LFE in two members' catalogues) are counted once; attribution temperature is fitted on the BUILD half, whose events are inside the member templates (in-sample for the templates) — the TEST half and A2011 are the honest numbers.

## Method

- **Group rule**: complete-linkage clustering of the epicentral distances (lat_use/lon_use: location, else the catalogue's approximate position), cut at 3.0 km → every pair in a group ≤ 3 km, no chaining. Family without a position (OLY:526) = own group. Single-linkage (for reference) would chain into groups of 20, 24 and 33.
- **Alignment** (deviation, data-driven): each family's detection time is on its own origin convention (offsets between co-located families up to ±9 s). A member's shift δ to the reference member (most build events) = the lag maximising the summed 3-C normalised cross-correlation of their build-half long stacks over all shared station-nodes (±10 s), searched within ±1 s of the co-detection mode when ≥ 20 of the member's detections have a reference detection within 8 s and the mode bin holds ≥ max(5, 3× uniform); per-station residual τ within δ ± 1 s. The own S picks ('template offsets') were not used: on the ≥ 30-event stacks they are noise-dominated at most stations. Where δ is strong it equals the co-detection offset to ≤ 0.15 s.
- **Truth**: union of member array detections in the group frame (t + δ), merged ±2 s; a group-day is scanned when every member whose windows include it was scanned. Templates: all members' build-half events (V2008 + V2010 half 0, co-detections once) cut at t + τ, wa_templates rule, weights on the T1 build scale. Every ≥ 30-event template kept (clear = a column).
- **Group MF**: `wa_mf.py --all-templates --drop-dead 0.05 --n-control 2` (seed as T1 → identical site shifts) + a second run with seed 20261001 (2 more scrambled controls, as step 1). Scoring `scripts/step2_score.py` (see its docstring).
- **Attribution** `scripts/step2_attrib.py`: member templates cut from each member's own build stack at the group window (≥ 30 events at that node); nearest 3 sites with ≥ 2 member templates live that day; eligible members = template at every chosen station; cc per station at the aligned position and ±0.5 s (sub-sample lag by parabola); P = softmax(β Σ cc0), β by maximum likelihood on build-half candidates matched to one member.
- **P(real)**: step-1 lfdr (`step1_features.py` + `step1_fit.py`, same features, 4 scrambled controls, fit on build half) refitted on the group candidates → `results/step2_groups/step1/`.

Checks: singleton group templates vs T1 build templates: 4913 compared, max |diff| 0.0, offset mismatches 0. Singleton group detections vs T1: 169052 vs 169052, identical = True.

### Group sizes

|   n_members |   groups |   families |   max_pair_km |   with_unlocated |
|------------:|---------:|-----------:|--------------:|-----------------:|
|       1.000 |  150.000 |    150.000 |         0.000 |            6.000 |
|       2.000 |   49.000 |     98.000 |         2.981 |           11.000 |
|       3.000 |    7.000 |     21.000 |         2.766 |            2.000 |
|       4.000 |    4.000 |     16.000 |         2.438 |            4.000 |
|       5.000 |    3.000 |     15.000 |         2.668 |            2.000 |
|       6.000 |    1.000 |      6.000 |         2.452 |            0.000 |
|      23.000 |    1.000 |     23.000 |         2.120 |            1.000 |

| region            |   groups |   multi |   families |
|:------------------|---------:|--------:|-----------:|
| Olympic_2010      |       38 |       7 |         68 |
| Olympic_2010_armB |       16 |      11 |         40 |
| Royer_NW_WA       |      104 |      36 |        144 |
| Royer_SW_WA       |       57 |      11 |         77 |

Member alignment (non-reference members):

| region            |   members |   codet_prior |   reliable |   unaligned |
|:------------------|----------:|--------------:|-----------:|------------:|
| Olympic_2010      |        38 |            14 |         17 |           7 |
| Olympic_2010_armB |        26 |             4 |          7 |           0 |
| Royer_NW_WA       |        28 |             6 |         10 |           0 |
| Royer_SW_WA       |        22 |             4 |          7 |           0 |

## Group detection vs T1 summed over members

`grp` = group matched filter; `t1union` = T1 family detections of the members in the group frame, merged ±2 s. Recall = union truth matched one-to-one ±2 s on observable scanned days (`recall all` = over ALL scanned truth); `recall al` = truth of well-aligned members only; precision own = matched / ours, any = within ±2 s of any group's truth; control = the scrambled control's chance precision. sep = (data+1)/(control+1) on label days; control column = the same on control 2 vs control 1.

**Tier `core`, raw** (tier only)

| region            | run     |   n_groups_obs |   n_families_obs |   n_truth_test |   n_truth_scored_all |   recall |   recall_all |   recall_alignok |   prec |   prec_any |   prec_ctl |   data_per_gday |   ctl_per_gday |   sep2_data |   sep2_ctl |   sep3_data |   sep3_ctl |   sep5_data |   sep5_ctl |
|:------------------|:--------|---------------:|-----------------:|---------------:|---------------------:|---------:|-------------:|-----------------:|-------:|-----------:|-----------:|----------------:|---------------:|------------:|-----------:|------------:|-----------:|------------:|-----------:|
| ALL               | grp     |            181 |              293 |          82540 |                83015 |    0.139 |        0.138 |            0.159 |  0.150 |      0.272 |      0.003 |          13.370 |         10.731 |         278 |         33 |         255 |          7 |         210 |          2 |
| ALL               | t1union |            180 |              291 |          82507 |                83015 |    0.182 |        0.181 |            0.189 |  0.093 |      0.192 |      0.003 |          26.467 |         22.672 |         285 |         38 |         256 |         11 |         212 |          3 |
| Olympic_2010      | grp     |              9 |               37 |           6125 |                 6514 |    0.294 |        0.276 |            0.445 |  0.055 |      0.186 |      0.002 |         163.843 |        161.112 |           3 |          0 |           3 |          0 |           3 |          0 |
| Olympic_2010      | t1union |              8 |               35 |           6092 |                 6514 |    0.559 |        0.523 |            0.734 |  0.030 |      0.122 |      0.003 |         538.523 |        531.184 |           2 |          0 |           2 |          0 |           2 |          0 |
| Olympic_2010_armB | grp     |             16 |               40 |           6098 |                 6098 |    0.002 |        0.002 |            0.003 |  0.200 |      0.477 |      0.000 |           1.759 |          0.335 |          13 |          2 |          13 |          0 |          12 |          0 |
| Olympic_2010_armB | t1union |             16 |               40 |           6098 |                 6098 |    0.007 |        0.007 |            0.008 |  0.084 |      0.403 |      0.008 |           7.451 |          2.342 |          14 |          3 |          14 |          0 |          13 |          0 |
| Royer_NW_WA       | grp     |            104 |              144 |          63131 |                63131 |    0.152 |        0.152 |            0.169 |  0.220 |      0.336 |      0.004 |           9.197 |          6.533 |         210 |         24 |         187 |          6 |         143 |          2 |
| Royer_NW_WA       | t1union |            104 |              144 |          63131 |                63131 |    0.182 |        0.182 |            0.186 |  0.239 |      0.356 |      0.005 |          10.495 |          6.987 |         216 |         26 |         188 |          9 |         145 |          3 |
| Royer_SW_WA       | grp     |             52 |               72 |           7186 |                 7272 |    0.007 |        0.007 |            0.008 |  0.810 |      0.873 |      0.000 |           2.893 |          0.011 |          52 |          7 |          52 |          1 |          52 |          0 |
| Royer_SW_WA       | t1union |             52 |               72 |           7186 |                 7272 |    0.012 |        0.012 |            0.013 |  0.837 |      0.875 |      0.000 |           4.141 |          0.014 |          53 |          9 |          52 |          2 |          52 |          0 |

**Tier `core`, cat** (tier & HF veto & cross-group rule)

| region            | run     |   n_groups_obs |   n_families_obs |   n_truth_test |   n_truth_scored_all |   recall |   recall_all |   recall_alignok |   prec |   prec_any |   prec_ctl |   data_per_gday |   ctl_per_gday |   sep2_data |   sep2_ctl |   sep3_data |   sep3_ctl |   sep5_data |   sep5_ctl |
|:------------------|:--------|---------------:|-----------------:|---------------:|---------------------:|---------:|-------------:|-----------------:|-------:|-----------:|-----------:|----------------:|---------------:|------------:|-----------:|------------:|-----------:|------------:|-----------:|
| ALL               | grp     |            181 |              293 |          82540 |                83015 |    0.099 |        0.098 |            0.114 |  0.167 |      0.261 |      0.003 |          13.370 |         10.731 |         278 |         33 |         255 |          7 |         210 |          2 |
| ALL               | t1union |            180 |              291 |          82507 |                83015 |    0.120 |        0.119 |            0.122 |  0.100 |      0.172 |      0.003 |          26.467 |         22.672 |         285 |         38 |         256 |         11 |         212 |          3 |
| Olympic_2010      | grp     |              9 |               37 |           6125 |                 6514 |    0.182 |        0.171 |            0.276 |  0.063 |      0.179 |      0.002 |         163.843 |        161.112 |           3 |          0 |           3 |          0 |           3 |          0 |
| Olympic_2010      | t1union |              8 |               35 |           6092 |                 6514 |    0.355 |        0.332 |            0.458 |  0.032 |      0.102 |      0.003 |         538.523 |        531.184 |           2 |          0 |           2 |          0 |           2 |          0 |
| Olympic_2010_armB | grp     |             16 |               40 |           6098 |                 6098 |    0.001 |        0.001 |            0.001 |  0.192 |      0.346 |      0.000 |           1.759 |          0.335 |          13 |          2 |          13 |          0 |          12 |          0 |
| Olympic_2010_armB | t1union |             16 |               40 |           6098 |                 6098 |    0.002 |        0.002 |            0.003 |  0.203 |      0.311 |      0.008 |           7.451 |          2.342 |          14 |          3 |          14 |          0 |          13 |          0 |
| Royer_NW_WA       | grp     |            104 |              144 |          63131 |                63131 |    0.111 |        0.111 |            0.123 |  0.223 |      0.305 |      0.004 |           9.197 |          6.533 |         210 |         24 |         187 |          6 |         143 |          2 |
| Royer_NW_WA       | t1union |            104 |              144 |          63131 |                63131 |    0.121 |        0.121 |            0.121 |  0.250 |      0.328 |      0.005 |          10.495 |          6.987 |         216 |         26 |         188 |          9 |         145 |          3 |
| Royer_SW_WA       | grp     |             52 |               72 |           7186 |                 7272 |    0.007 |        0.007 |            0.008 |  0.833 |      0.883 |      0.000 |           2.893 |          0.011 |          52 |          7 |          52 |          1 |          52 |          0 |
| Royer_SW_WA       | t1union |             52 |               72 |           7186 |                 7272 |    0.010 |        0.010 |            0.011 |  0.841 |      0.864 |      0.000 |           4.141 |          0.014 |          53 |          9 |          52 |          2 |          52 |          0 |

**Tier `core_ep`, raw** (tier only)

| region            | run     |   n_groups_obs |   n_families_obs |   n_truth_test |   n_truth_scored_all |   recall |   recall_all |   recall_alignok |    prec |   prec_any |   prec_ctl |   data_per_gday |   ctl_per_gday |   sep2_data |   sep2_ctl |   sep3_data |   sep3_ctl |   sep5_data |   sep5_ctl |
|:------------------|:--------|---------------:|-----------------:|---------------:|---------------------:|---------:|-------------:|-----------------:|--------:|-----------:|-----------:|----------------:|---------------:|------------:|-----------:|------------:|-----------:|------------:|-----------:|
| ALL               | grp     |            181 |              293 |          82540 |                83015 |    0.002 |        0.002 |            0.002 |   0.348 |      0.468 |      0.000 |           0.107 |          0.065 |         106 |          0 |          58 |          0 |          17 |          0 |
| ALL               | t1union |            180 |              291 |          82507 |                83015 |    0.002 |        0.002 |            0.002 |   0.131 |      0.211 |      0.005 |           0.185 |          0.155 |          69 |          0 |          22 |          0 |           6 |          0 |
| Olympic_2010      | grp     |              9 |               37 |           6125 |                 6514 |    0.014 |        0.013 |            0.022 |   0.268 |      0.406 |      0.000 |           2.216 |          1.707 |           2 |          0 |           2 |          0 |           2 |          0 |
| Olympic_2010      | t1union |              8 |               35 |           6092 |                 6514 |    0.022 |        0.020 |            0.025 |   0.113 |      0.195 |      0.005 |           5.368 |          4.813 |           4 |          0 |           2 |          0 |           0 |          0 |
| Olympic_2010_armB | grp     |             16 |               40 |           6098 |                 6098 |    0.000 |        0.000 |            0.000 | nan     |    nan     |    nan     |           0.051 |          0.000 |           5 |          0 |           4 |          0 |           2 |          0 |
| Olympic_2010_armB | t1union |             16 |               40 |           6098 |                 6098 |    0.000 |        0.000 |            0.000 |   1.000 |      1.000 |    nan     |           0.049 |          0.000 |           6 |          0 |           3 |          0 |           1 |          0 |
| Royer_NW_WA       | grp     |            104 |              144 |          63131 |                63131 |    0.001 |        0.001 |            0.001 |   0.612 |      0.673 |      0.000 |           0.032 |          0.006 |          90 |          0 |          48 |          0 |          11 |          0 |
| Royer_NW_WA       | t1union |            104 |              144 |          63131 |                63131 |    0.000 |        0.000 |            0.000 |   0.958 |      0.958 |    nan     |           0.012 |          0.000 |          55 |          0 |          17 |          0 |           5 |          0 |
| Royer_SW_WA       | grp     |             52 |               72 |           7186 |                 7272 |    0.000 |        0.000 |            0.000 | nan     |    nan     |    nan     |           0.016 |          0.000 |           9 |          0 |           4 |          0 |           2 |          0 |
| Royer_SW_WA       | t1union |             52 |               72 |           7186 |                 7272 |    0.000 |        0.000 |            0.000 | nan     |    nan     |    nan     |           0.003 |          0.000 |           4 |          0 |           0 |          0 |           0 |          0 |

**Tier `ext`, raw** (tier only)

| region            | run     |   n_groups_obs |   n_families_obs |   n_truth_test |   n_truth_scored_all |   recall |   recall_all |   recall_alignok |   prec |   prec_any |   prec_ctl |   data_per_gday |   ctl_per_gday |   sep2_data |   sep2_ctl |   sep3_data |   sep3_ctl |   sep5_data |   sep5_ctl |
|:------------------|:--------|---------------:|-----------------:|---------------:|---------------------:|---------:|-------------:|-----------------:|-------:|-----------:|-----------:|----------------:|---------------:|------------:|-----------:|------------:|-----------:|------------:|-----------:|
| ALL               | grp     |            181 |              293 |          82540 |                83015 |    0.165 |        0.164 |            0.189 |  0.136 |      0.258 |      0.003 |          16.935 |         13.903 |         263 |         37 |         223 |         12 |         171 |          4 |
| ALL               | t1union |            180 |              291 |          82507 |                83015 |    0.215 |        0.213 |            0.222 |  0.087 |      0.187 |      0.003 |          32.437 |         28.119 |         264 |         35 |         224 |         11 |         167 |          2 |
| Olympic_2010      | grp     |              9 |               37 |           6125 |                 6514 |    0.303 |        0.285 |            0.458 |  0.049 |      0.172 |      0.002 |         189.778 |        187.240 |           3 |          0 |           3 |          0 |           2 |          0 |
| Olympic_2010      | t1union |              8 |               35 |           6092 |                 6514 |    0.581 |        0.544 |            0.746 |  0.027 |      0.114 |      0.003 |         627.113 |        618.268 |           2 |          0 |           2 |          0 |           1 |          0 |
| Olympic_2010_armB | grp     |             16 |               40 |           6098 |                 6098 |    0.005 |        0.005 |            0.006 |  0.139 |      0.393 |      0.015 |           2.621 |          0.942 |          13 |          0 |          13 |          0 |          11 |          0 |
| Olympic_2010_armB | t1union |             16 |               40 |           6098 |                 6098 |    0.012 |        0.012 |            0.014 |  0.071 |      0.333 |      0.007 |          10.004 |          4.433 |          14 |          3 |          13 |          2 |          13 |          0 |
| Royer_NW_WA       | grp     |            104 |              144 |          63131 |                63131 |    0.185 |        0.185 |            0.205 |  0.188 |      0.309 |      0.004 |          12.518 |          9.400 |         194 |         22 |         155 |          6 |         106 |          3 |
| Royer_NW_WA       | t1union |            104 |              144 |          63131 |                63131 |    0.221 |        0.221 |            0.225 |  0.202 |      0.323 |      0.004 |          14.285 |         10.242 |         196 |         18 |         157 |          4 |         101 |          1 |
| Royer_SW_WA       | grp     |             52 |               72 |           7186 |                 7272 |    0.010 |        0.010 |            0.012 |  0.610 |      0.695 |      0.000 |           3.179 |          0.052 |          53 |         15 |          52 |          6 |          52 |          1 |
| Royer_SW_WA       | t1union |             52 |               72 |           7186 |                 7272 |    0.016 |        0.016 |            0.018 |  0.654 |      0.715 |      0.000 |           4.466 |          0.086 |          52 |         14 |          52 |          5 |          52 |          1 |

By group size (core, raw):

| size_class   | run     |   n_groups_obs |   n_truth_test |   recall_raw_core |   recall_all_raw_core |   prec_raw_core |   prec_ctl_core |   data_per_gday_core |   ctl_per_gday_core |   sep3_core_data |   sep3_core_ctl |
|:-------------|:--------|---------------:|---------------:|------------------:|----------------------:|----------------:|----------------:|---------------------:|--------------------:|-----------------:|----------------:|
| 1            | grp     |            118 |          36696 |             0.162 |                 0.160 |           0.169 |           0.003 |               10.452 |               7.604 |              177 |               6 |
| 1            | t1union |            118 |          36696 |             0.162 |                 0.160 |           0.169 |           0.003 |               10.452 |               7.604 |              177 |               6 |
| 2            | grp     |             47 |          33310 |             0.136 |                 0.135 |           0.144 |           0.004 |               17.815 |              15.464 |               62 |               1 |
| 2            | t1union |             46 |          33277 |             0.185 |                 0.184 |           0.169 |           0.004 |               22.057 |              17.432 |               60 |               5 |
| 23           | grp     |              1 |           3732 |             0.178 |                 0.178 |           0.072 |           0.003 |              185.239 |             182.296 |                0 |               0 |
| 23           | t1union |              1 |           3732 |             0.595 |                 0.595 |           0.025 |           0.003 |             1729.662 |            1709.099 |                0 |               0 |
| 3-6          | grp     |             15 |           8802 |             0.041 |                 0.041 |           0.598 |           0.004 |                2.481 |               0.713 |               16 |               0 |
| 3-6          | t1union |             15 |           8802 |             0.084 |                 0.084 |           0.459 |           0.004 |               10.577 |               2.654 |               19 |               0 |

By window and region (core, raw):

| window   | region            | run     |   n_groups_obs |   n_truth_test |   recall_raw_core |   recall_all_raw_core |   prec_raw_core |   prec_ctl_core |   data_per_gday_core |   ctl_per_gday_core |   sep3_core_data |   sep3_core_ctl |
|:---------|:------------------|:--------|---------------:|---------------:|------------------:|----------------------:|----------------:|----------------:|---------------------:|--------------------:|-----------------:|----------------:|
| A2011    | Olympic_2010      | grp     |              2 |           1023 |             0.008 |                 0.008 |           0.001 |           0.000 |              203.052 |             201.569 |                0 |               0 |
| A2011    | Olympic_2010      | t1union |              2 |           1023 |             0.247 |                 0.247 |           0.004 |           0.000 |             1113.310 |            1107.862 |                0 |               0 |
| A2011    | Royer_NW_WA       | grp     |            104 |          31345 |             0.171 |                 0.171 |           0.177 |           0.004 |               10.055 |               8.982 |               28 |               3 |
| A2011    | Royer_NW_WA       | t1union |            104 |          31345 |             0.203 |                 0.203 |           0.192 |           0.004 |               11.005 |               9.673 |               27 |               5 |
| V2008    | Olympic_2010      | grp     |              2 |            440 |             0.009 |                 0.009 |           0.001 |           0.001 |              145.913 |             142.717 |                0 |               0 |
| V2008    | Olympic_2010      | t1union |              2 |            440 |             0.152 |                 0.152 |           0.004 |           0.001 |              744.630 |             732.717 |                0 |               0 |
| V2008    | Royer_NW_WA       | grp     |            104 |          15388 |             0.139 |                 0.139 |           0.299 |           0.005 |                8.523 |               5.148 |               81 |               2 |
| V2008    | Royer_NW_WA       | t1union |            104 |          15388 |             0.169 |                 0.169 |           0.329 |           0.005 |               10.021 |               5.445 |               81 |               3 |
| V2008    | Royer_SW_WA       | grp     |             52 |           7162 |             0.007 |                 0.007 |           0.810 |           0.000 |                2.953 |               0.011 |               52 |               1 |
| V2008    | Royer_SW_WA       | t1union |             52 |           7162 |             0.012 |                 0.012 |           0.837 |           0.000 |                4.226 |               0.015 |               52 |               2 |
| V2010    | Olympic_2010      | grp     |              9 |           4662 |             0.384 |                 0.354 |           0.102 |           0.003 |              157.650 |             154.701 |                3 |               0 |
| V2010    | Olympic_2010      | t1union |              8 |           4629 |             0.666 |                 0.611 |           0.095 |           0.009 |              330.665 |             323.816 |                2 |               0 |
| V2010    | Olympic_2010_armB | grp     |             16 |           6098 |             0.002 |                 0.002 |           0.200 |           0.000 |                1.759 |               0.335 |               13 |               0 |
| V2010    | Olympic_2010_armB | t1union |             16 |           6098 |             0.007 |                 0.007 |           0.084 |           0.008 |                7.451 |               2.342 |               14 |               0 |
| V2010    | Royer_NW_WA       | grp     |            104 |          16398 |             0.130 |                 0.130 |           0.337 |           0.009 |                8.725 |               4.564 |               78 |               1 |
| V2010    | Royer_NW_WA       | t1union |            104 |          16398 |             0.156 |                 0.156 |           0.365 |           0.010 |               10.298 |               4.853 |               80 |               1 |
| V2010    | Royer_SW_WA       | grp     |              1 |             24 |             0.000 |                 0.000 |         nan     |         nan     |                0.000 |               0.000 |                0 |               0 |
| V2010    | Royer_SW_WA       | t1union |              1 |             24 |             0.000 |                 0.000 |         nan     |         nan     |                0.000 |               0.000 |                0 |               0 |

Group-run per-epoch tiers:

| window   |   n_quiet_famdays | quiet_kinds                  |   T_ext_ep |   T_core_ep | derived   |   fa_at_8p25 |   fa_at_8p75 |
|:---------|------------------:|:-----------------------------|-----------:|------------:|:----------|-------------:|-------------:|
| A2011    |                 0 | none: V2010 epoch tiers used |     14.750 |      18.500 | False     |      nan     |      nan     |
| V2008    |               543 | post                         |     16.500 |      20.250 | True      |       13.884 |       11.449 |
| V2010    |              1810 | quiet                        |     14.750 |      18.500 | True      |       15.581 |       12.517 |

Separation curve (core, raw; number of group-windows with label days whose data/control ratio ≥ x, and the same on the control):

|   sep ≥ |   n_ctl grp |   n_ctl t1union |   n_data grp |   n_data t1union |
|--------:|------------:|----------------:|-------------:|-----------------:|
|   1.000 |     267.000 |         261.000 |      387.000 |          387.000 |
|   1.250 |      65.000 |          64.000 |      340.000 |          343.000 |
|   1.500 |      44.000 |          45.000 |      313.000 |          317.000 |
|   2.000 |      33.000 |          38.000 |      278.000 |          285.000 |
|   3.000 |       7.000 |          11.000 |      255.000 |          256.000 |
|   5.000 |       2.000 |           3.000 |      210.000 |          212.000 |
|  10.000 |       0.000 |           0.000 |      164.000 |          164.000 |

## Attribution — which member fired

Fitted on 5682 build-half candidates: cc_k3 β = 8.52; cc_k2 β = 9.63; cc_k1 β = 10.04; ccmax_k3 β = 10.00; cc_lag_k3 β = 8.30, -1.01; cc_k3_prior β = 6.86, 1.81; prior_only β = 1.00.

All groups (candidates matched to exactly one eligible member's array detection). `top1` = the highest-P member is the true one; `chance` = mean 1/n eligible; `most_active` = always pick the member with the most build events; skill = 1 − log-loss / uniform log-loss.

| eval_set   | model       |   n_cand |   top1 |   top1_lo |   top1_hi |   chance |   most_active_top1 |   mean_P_true |   mean_P_unif |   logloss_skill_vs_unif |   mean_n_elig |
|:-----------|:------------|---------:|-------:|----------:|----------:|---------:|-------------------:|--------------:|--------------:|------------------------:|--------------:|
| fit        | cc_k3       |     5682 |  0.900 |     0.892 |     0.908 |    0.457 |              0.904 |         0.841 |         0.457 |                   0.650 |         2.961 |
| fit        | cc_k2       |     5682 |  0.879 |     0.870 |     0.887 |    0.457 |              0.904 |         0.818 |         0.457 |                   0.612 |         2.961 |
| fit        | cc_k1       |     5682 |  0.797 |     0.787 |     0.808 |    0.457 |              0.904 |         0.723 |         0.457 |                   0.457 |         2.961 |
| fit        | ccmax_k3    |     5682 |  0.841 |     0.831 |     0.850 |    0.457 |              0.904 |         0.774 |         0.457 |                   0.538 |         2.961 |
| fit        | cc_lag_k3   |     5682 |  0.899 |     0.891 |     0.907 |    0.457 |              0.904 |         0.842 |         0.457 |                   0.650 |         2.961 |
| fit        | cc_k3_prior |     5682 |  0.931 |     0.924 |     0.937 |    0.457 |              0.904 |         0.879 |         0.457 |                   0.724 |         2.961 |
| fit        | prior_only  |     5682 |  0.904 |     0.896 |     0.911 |    0.457 |              0.904 |         0.613 |         0.457 |                   0.354 |         2.961 |
| V2010_test | cc_k3       |     1231 |  0.767 |     0.742 |     0.790 |    0.371 |              0.768 |         0.737 |         0.371 |                   0.439 |         5.621 |
| V2010_test | cc_k2       |     1231 |  0.774 |     0.750 |     0.797 |    0.371 |              0.768 |         0.739 |         0.371 |                   0.478 |         5.621 |
| V2010_test | cc_k1       |     1231 |  0.617 |     0.589 |     0.643 |    0.371 |              0.768 |         0.567 |         0.371 |                   0.283 |         5.621 |
| V2010_test | ccmax_k3    |     1231 |  0.764 |     0.740 |     0.787 |    0.371 |              0.768 |         0.736 |         0.371 |                   0.483 |         5.621 |
| V2010_test | cc_lag_k3   |     1231 |  0.765 |     0.741 |     0.788 |    0.371 |              0.768 |         0.737 |         0.371 |                   0.441 |         5.621 |
| V2010_test | cc_k3_prior |     1231 |  0.807 |     0.785 |     0.829 |    0.371 |              0.768 |         0.780 |         0.371 |                   0.488 |         5.621 |
| V2010_test | prior_only  |     1231 |  0.768 |     0.744 |     0.791 |    0.371 |              0.768 |         0.564 |         0.371 |                   0.372 |         5.621 |
| V2008_test | cc_k3       |      488 |  0.887 |     0.856 |     0.912 |    0.490 |              0.867 |         0.850 |         0.490 |                   0.600 |         2.098 |
| V2008_test | cc_k2       |      488 |  0.883 |     0.852 |     0.909 |    0.490 |              0.867 |         0.842 |         0.490 |                   0.548 |         2.098 |
| V2008_test | cc_k1       |      488 |  0.816 |     0.779 |     0.847 |    0.490 |              0.867 |         0.766 |         0.490 |                   0.436 |         2.098 |
| V2008_test | ccmax_k3    |      488 |  0.877 |     0.845 |     0.903 |    0.490 |              0.867 |         0.835 |         0.490 |                   0.555 |         2.098 |
| V2008_test | cc_lag_k3   |      488 |  0.875 |     0.843 |     0.901 |    0.490 |              0.867 |         0.849 |         0.490 |                   0.601 |         2.098 |
| V2008_test | cc_k3_prior |      488 |  0.928 |     0.902 |     0.948 |    0.490 |              0.867 |         0.891 |         0.490 |                   0.693 |         2.098 |
| V2008_test | prior_only  |      488 |  0.867 |     0.834 |     0.894 |    0.490 |              0.867 |         0.639 |         0.490 |                   0.338 |         2.098 |
| A2011      | cc_k3       |     1319 |  0.931 |     0.916 |     0.943 |    0.492 |              0.888 |         0.900 |         0.492 |                   0.721 |         2.083 |
| A2011      | cc_k2       |     1319 |  0.901 |     0.884 |     0.916 |    0.492 |              0.888 |         0.881 |         0.492 |                   0.672 |         2.083 |
| A2011      | cc_k1       |     1319 |  0.817 |     0.795 |     0.836 |    0.492 |              0.888 |         0.772 |         0.492 |                   0.433 |         2.083 |
| A2011      | ccmax_k3    |     1319 |  0.917 |     0.900 |     0.930 |    0.492 |              0.888 |         0.886 |         0.492 |                   0.696 |         2.083 |
| A2011      | cc_lag_k3   |     1319 |  0.933 |     0.919 |     0.946 |    0.492 |              0.888 |         0.899 |         0.492 |                   0.720 |         2.083 |
| A2011      | cc_k3_prior |     1319 |  0.941 |     0.927 |     0.952 |    0.492 |              0.888 |         0.908 |         0.492 |                   0.745 |         2.083 |
| A2011      | prior_only  |     1319 |  0.888 |     0.870 |     0.904 |    0.492 |              0.888 |         0.638 |         0.492 |                   0.331 |         2.083 |

Primary model (cc, 3 nearest sites) by region / group size / tier:

| eval_set   | region            | size_class   | tier   |   n_cand |   top1 |   chance |   most_active_top1 |   mean_P_true |   logloss_skill_vs_unif |
|:-----------|:------------------|:-------------|:-------|---------:|-------:|---------:|-------------------:|--------------:|------------------------:|
| fit        | Olympic_2010      | ALL          | all    |      898 |  0.697 |    0.311 |              0.716 |         0.675 |                   0.398 |
| fit        | Olympic_2010_armB | ALL          | all    |      350 |  0.937 |    0.400 |              0.969 |         0.828 |                   0.736 |
| fit        | Royer_NW_WA       | ALL          | all    |     3920 |  0.946 |    0.496 |              0.937 |         0.890 |                   0.779 |
| fit        | Royer_SW_WA       | ALL          | all    |      514 |  0.879 |    0.455 |              0.936 |         0.768 |                   0.570 |
| fit        | ALL               | 2            | all    |     4567 |  0.947 |    0.500 |              0.950 |         0.893 |                   0.778 |
| fit        | ALL               | 23           | all    |      395 |  0.316 |    0.070 |              0.359 |         0.270 |                   0.211 |
| fit        | ALL               | 3-6          | all    |      720 |  0.925 |    0.395 |              0.912 |         0.828 |                   0.730 |
| fit        | ALL               | ALL          | core   |     4392 |  0.902 |    0.455 |              0.908 |         0.850 |                   0.661 |
| V2010_test | Olympic_2010      | ALL          | all    |      746 |  0.676 |    0.294 |              0.665 |         0.645 |                   0.367 |
| V2010_test | Olympic_2010_armB | ALL          | all    |       21 |  0.619 |    0.343 |              0.762 |         0.562 |                   0.124 |
| V2010_test | Royer_NW_WA       | ALL          | all    |      464 |  0.920 |    0.497 |              0.935 |         0.892 |                   0.735 |
| V2010_test | ALL               | 2            | all    |      800 |  0.951 |    0.500 |              0.969 |         0.938 |                   0.854 |
| V2010_test | ALL               | 23           | all    |      358 |  0.324 |    0.070 |              0.302 |         0.262 |                   0.189 |
| V2010_test | ALL               | 3-6          | all    |       73 |  0.918 |    0.437 |              0.863 |         0.865 |                   0.579 |
| V2010_test | ALL               | ALL          | core   |     1034 |  0.765 |    0.361 |              0.762 |         0.736 |                   0.451 |
| V2008_test | Olympic_2010      | ALL          | all    |        4 |  0.000 |    0.286 |              0.000 |         0.001 |                  -4.526 |
| V2008_test | Royer_NW_WA       | ALL          | all    |      480 |  0.898 |    0.492 |              0.873 |         0.859 |                   0.704 |
| V2008_test | Royer_SW_WA       | ALL          | all    |        4 |  0.500 |    0.500 |              1.000 |         0.627 |                   0.103 |
| V2008_test | ALL               | 2            | all    |      424 |  0.899 |    0.500 |              0.882 |         0.861 |                   0.666 |
| V2008_test | ALL               | 23           | all    |        2 |  0.000 |    0.071 |              0.000 |         0.000 |                  -3.415 |
| V2008_test | ALL               | 3-6          | all    |       62 |  0.839 |    0.435 |              0.790 |         0.800 |                   0.632 |
| V2008_test | ALL               | ALL          | core   |      321 |  0.885 |    0.490 |              0.885 |         0.854 |                   0.556 |
| A2011      | Olympic_2010      | ALL          | all    |        9 |  0.000 |    0.307 |              0.000 |         0.005 |                  -3.871 |
| A2011      | Royer_NW_WA       | ALL          | all    |     1310 |  0.937 |    0.493 |              0.894 |         0.906 |                   0.792 |
| A2011      | ALL               | 2            | all    |     1167 |  0.937 |    0.500 |              0.915 |         0.904 |                   0.753 |
| A2011      | ALL               | 23           | all    |        4 |  0.000 |    0.065 |              0.000 |         0.000 |                  -2.076 |
| A2011      | ALL               | 3-6          | all    |      148 |  0.912 |    0.438 |              0.696 |         0.890 |                   0.755 |
| A2011      | ALL               | ALL          | core   |      905 |  0.944 |    0.492 |              0.907 |         0.919 |                   0.745 |

Co-detected events (the same LFE in ≥ 2 members' catalogues): top member in the true set vs chance |set|/n:

| eval_set   |   n_cand |   top1_in_set |   chance |
|:-----------|---------:|--------------:|---------:|
| A2011      |      375 |         0.960 |    0.897 |
| V2008_test |      172 |         0.988 |    0.836 |
| V2010_test |      463 |         0.687 |    0.389 |
| fit        |     1196 |         0.875 |    0.678 |
| guard      |       23 |         0.826 |    0.514 |

How peaked P(member) is on data vs the scrambled control (all record-tier candidates; the control has no true member):

| eval_set   | series   |   n_cand |   mean_Pmax |   mean_Pmax_unif |   mean_norm_entropy |
|:-----------|:---------|---------:|------------:|-----------------:|--------------------:|
| V2010_test | data     |    12527 |       0.772 |            0.319 |               0.305 |
| V2010_test | control  |    11888 |       0.761 |            0.310 |               0.310 |
| V2008_test | data     |     8156 |       0.759 |            0.312 |               0.365 |
| V2008_test | control  |     7600 |       0.748 |            0.300 |               0.365 |
| A2011      | data     |    32008 |       0.819 |            0.353 |               0.305 |
| A2011      | control  |    30509 |       0.817 |            0.347 |               0.302 |

Reliability of P(member) (every candidate × member pair):

| P bin   |   frac_true A2011 |   frac_true V2008_test |   frac_true V2010_test |   mean_P A2011 |   mean_P V2008_test |   mean_P V2010_test |   n A2011 |   n V2008_test |   n V2010_test |
|:--------|------------------:|-----------------------:|-----------------------:|---------------:|--------------------:|--------------------:|----------:|---------------:|---------------:|
| 0.0-0.1 |             0.018 |                  0.019 |                  0.026 |          0.012 |               0.017 |               0.012 |  1140.000 |        360.000 |       4985.000 |
| 0.1-0.2 |             0.123 |                  0.090 |                  0.185 |          0.146 |               0.135 |               0.142 |   106.000 |         67.000 |        405.000 |
| 0.2-0.3 |             0.242 |                  0.341 |                  0.234 |          0.252 |               0.256 |               0.248 |    66.000 |         41.000 |        222.000 |
| 0.3-0.4 |             0.345 |                  0.385 |                  0.279 |          0.357 |               0.357 |               0.346 |    55.000 |         39.000 |        140.000 |
| 0.4-0.5 |             0.358 |                  0.457 |                  0.417 |          0.450 |               0.453 |               0.448 |    67.000 |         35.000 |        108.000 |
| 0.5-0.6 |             0.594 |                  0.485 |                  0.351 |          0.549 |               0.548 |               0.553 |    64.000 |         33.000 |         77.000 |
| 0.6-0.7 |             0.655 |                  0.615 |                  0.422 |          0.644 |               0.646 |               0.649 |    58.000 |         39.000 |         90.000 |
| 0.7-0.8 |             0.787 |                  0.667 |                  0.524 |          0.748 |               0.746 |               0.749 |    61.000 |         30.000 |         63.000 |
| 0.8-0.9 |             0.835 |                  0.912 |                  0.583 |          0.850 |               0.865 |               0.852 |    85.000 |         57.000 |         60.000 |
| 0.9-1.0 |             0.988 |                  0.985 |                  0.983 |          0.988 |               0.983 |               0.992 |  1045.000 |        323.000 |        769.000 |

Per group (test sets pooled, ≥ 20 candidates):

|   gid | region       |   n_members |   n_cand |   top1 |   chance |   mean_P_true | true_member_counts                                                                                   |
|------:|:-------------|------------:|---------:|-------:|---------:|--------------:|:-----------------------------------------------------------------------------------------------------|
|  1031 | Olympic_2010 |           2 |      395 |  0.982 |    0.500 |         0.981 | 280:388;91:7                                                                                         |
|  1000 | Olympic_2010 |          23 |      364 |  0.319 |    0.070 |         0.257 | 259:108;251:99;277:87;285:12;271:12;217:8;108:7;244:7;243:6;274:6;249:5;228:2;279:2;66:1;281:1;213:1 |
|  1032 | Royer_NW_WA  |           2 |      326 |  0.982 |    0.500 |         0.940 | 45:313;110:13                                                                                        |
|  1025 | Royer_NW_WA  |           2 |      232 |  0.884 |    0.500 |         0.856 | 5:231;103:1                                                                                          |
|  1037 | Royer_NW_WA  |           2 |      185 |  1.000 |    0.500 |         0.995 | 35:185                                                                                               |
|  1012 | Royer_NW_WA  |           3 |      150 |  1.000 |    0.500 |         0.997 | 117:150                                                                                              |
|  1038 | Royer_NW_WA  |           2 |      143 |  0.979 |    0.500 |         0.939 | 26:110;61:33                                                                                         |
|  1041 | Royer_NW_WA  |           2 |      137 |  0.993 |    0.500 |         0.974 | 47:136;93:1                                                                                          |
|  1039 | Royer_NW_WA  |           2 |      129 |  1.000 |    0.500 |         0.986 | 63:129                                                                                               |
|  1046 | Royer_NW_WA  |           2 |      126 |  1.000 |    0.500 |         0.991 | 87:126                                                                                               |
|  1059 | Royer_NW_WA  |           2 |       85 |  0.988 |    0.500 |         0.969 | 129:85                                                                                               |
|  1029 | Royer_NW_WA  |           2 |       81 |  0.593 |    0.500 |         0.554 | 106:49;38:32                                                                                         |
|  1034 | Royer_NW_WA  |           2 |       77 |  0.805 |    0.500 |         0.706 | 104:50;123:27                                                                                        |
|  1040 | Royer_NW_WA  |           2 |       61 |  1.000 |    0.500 |         0.974 | 28:61                                                                                                |
|  1062 | Royer_NW_WA  |           2 |       60 |  0.967 |    0.500 |         0.965 | 11:60                                                                                                |
|  1054 | Royer_NW_WA  |           2 |       56 |  0.857 |    0.500 |         0.781 | 21:56                                                                                                |
|  1058 | Royer_NW_WA  |           2 |       54 |  0.611 |    0.500 |         0.545 | 112:31;127:23                                                                                        |
|  1048 | Royer_NW_WA  |           2 |       52 |  1.000 |    0.500 |         0.993 | 23:52                                                                                                |
|  1013 | Royer_NW_WA  |           3 |       50 |  0.880 |    0.333 |         0.771 | 20:49;99:1                                                                                           |
|  1050 | Royer_NW_WA  |           2 |       45 |  0.911 |    0.500 |         0.827 | 128:45                                                                                               |
|  1015 | Royer_NW_WA  |           3 |       37 |  0.459 |    0.333 |         0.431 | 73:21;52:16                                                                                          |
|  1011 | Royer_NW_WA  |           3 |       31 |  1.000 |    0.500 |         0.993 | 4:31                                                                                                 |
|  1023 | Royer_NW_WA  |           2 |       27 |  0.667 |    0.500 |         0.635 | 125:27                                                                                               |
|  1063 | Royer_NW_WA  |           2 |       25 |  1.000 |    0.500 |         0.973 | 100:25                                                                                               |
|  1055 | Royer_NW_WA  |           2 |       22 |  0.773 |    0.500 |         0.631 | 69:11;70:11                                                                                          |

## P(real) on group candidates (step-1 lfdr refitted)

| eval_set   | level           |      n |     sum_P |   slope_own |   matched_over_sumP |
|:-----------|:----------------|-------:|----------:|------------:|--------------------:|
| V2010_test | group           |  33310 |  3279.780 |       0.604 |               1.475 |
| V2010_test | family (step 1) | 115301 |  6031.465 |       0.650 |               1.287 |
| V2008_test | group           |  15848 |  2192.286 |       0.621 |               1.347 |
| V2008_test | family (step 1) |  13943 |  2494.646 |       0.556 |               1.546 |
| A2011      | group           |  67835 | 14735.434 |       0.275 |               0.473 |
| A2011      | family (step 1) |  61999 |  8464.658 |       0.544 |               1.047 |

Scrambled control vs data:

| eval_set   | level           |   data_n |   data_P_gt05 |   data_sumP |   ctl_P_gt05 |   ctl_P_gt05_never_seen |   ctl_sumP_per_real |
|:-----------|:----------------|---------:|--------------:|------------:|-------------:|------------------------:|--------------------:|
| V2010_test | group           |    56562 |         0.045 |    4598.409 |        0.009 |                   0.010 |            2528.560 |
| V2010_test | family (step 1) |   172854 |         0.021 |    8795.923 |        0.004 |                   0.005 |            5952.251 |
| V2008_test | group           |    33336 |         0.056 |    2950.702 |        0.008 |                   0.009 |            1377.723 |
| V2008_test | family (step 1) |    93754 |         0.029 |    5110.635 |        0.007 |                   0.007 |            3437.870 |
| A2011      | group           |   114826 |         0.128 |   21219.657 |        0.094 |                   0.058 |           16898.866 |
| A2011      | family (step 1) |   294824 |         0.028 |   24821.464 |        0.014 |                   0.036 |           20152.425 |

Group-level expected count (ΣP per group-day vs the union truth):

| eval_set   | region            |   group_days |      sumP |   sumP_ctl |   truth_own |   ratio_P_over_own |   spearman_own |   spearman_ncand_own |
|:-----------|:------------------|-------------:|----------:|-----------:|------------:|-------------------:|---------------:|---------------------:|
| V2010_test | ALL               |         2778 |  3279.780 |   1404.196 |       27182 |              0.121 |          0.473 |                0.454 |
| V2010_test | Royer_NW_WA       |         2068 |  2222.934 |    437.395 |       16398 |              0.136 |          0.476 |                0.453 |
| V2010_test | Royer_SW_WA       |           28 |     0.119 |      0.000 |          24 |              0.005 |         -0.130 |               -0.130 |
| V2010_test | Olympic_2010      |          234 |   999.863 |    932.927 |        4662 |              0.214 |          0.557 |                0.487 |
| V2010_test | Olympic_2010_armB |          448 |    56.864 |     33.874 |        6098 |              0.009 |          0.219 |                0.250 |
| V2008_test | ALL               |         3790 |  2192.286 |    572.415 |       22990 |              0.095 |          0.428 |                0.405 |
| V2008_test | Royer_NW_WA       |         2392 |  1941.989 |    416.414 |       15388 |              0.126 |          0.503 |                0.473 |
| V2008_test | Royer_SW_WA       |         1352 |   112.043 |     12.611 |        7162 |              0.016 |          0.203 |                0.175 |
| V2008_test | Olympic_2010      |           46 |   138.254 |    143.391 |         440 |              0.314 |         -0.016 |               -0.103 |
| A2011      | ALL               |         3074 | 14735.434 |  10545.535 |       32368 |              0.455 |          0.457 |                0.395 |
| A2011      | Royer_NW_WA       |         3016 | 12695.799 |   8682.846 |       31345 |              0.405 |          0.456 |                0.394 |
| A2011      | Olympic_2010      |           58 |  2039.635 |   1862.689 |        1023 |              1.994 |         -0.191 |               -0.310 |

### Family-level expected count Σ P(real)·P(member)

Per member family-day (episode days the family's array run scanned; scored half): `E` = attributed, `E_unif` = P(real) split equally over the group, `P_T1` = step 1's family-level ΣP; Spearman against the family's own array count.

| window   | region            | families                       |   n_famdays |   truth |     sum_E |   sum_E_unif |   sum_P_T1 |   spearman_E |   spearman_E_unif |   spearman_P_T1 |
|:---------|:------------------|:-------------------------------|------------:|--------:|----------:|-------------:|-----------:|-------------:|------------------:|----------------:|
| ALL      | ALL               | all families                   |       13334 |   88189 | 16954.289 |    17251.468 |  16990.769 |        0.451 |             0.449 |           0.508 |
| ALL      | ALL               | members of multi-member groups |        7120 |   51493 |  4884.350 |     5181.529 |   9453.906 |        0.410 |             0.408 |           0.513 |
| A2011    | ALL               | all families                   |        3828 |   33306 | 11695.057 |    12040.729 |   8464.658 |        0.399 |             0.419 |           0.476 |
| A2011    | ALL               | members of multi-member groups |        1856 |   17214 |  3138.128 |     3483.800 |   3784.436 |        0.352 |             0.402 |           0.512 |
| V2008    | ALL               | all families                   |        4934 |   23962 |  1957.561 |     1978.397 |   2494.646 |        0.393 |             0.392 |           0.452 |
| V2008    | ALL               | members of multi-member groups |        2304 |   12587 |   572.669 |      593.505 |   1352.317 |        0.353 |             0.350 |           0.469 |
| V2010    | ALL               | all families                   |        4572 |   30921 |  3301.671 |     3232.342 |   6031.465 |        0.440 |             0.430 |           0.491 |
| V2010    | ALL               | members of multi-member groups |        2960 |   21692 |  1173.553 |     1104.224 |   4317.153 |        0.396 |             0.381 |           0.481 |
| ALL      | Olympic_2010      | all families                   |         784 |    7056 |  1152.548 |     1044.868 |   3481.801 |        0.518 |             0.427 |           0.373 |
| ALL      | Olympic_2010      | members of multi-member groups |         672 |    6252 |   523.939 |      416.259 |   3020.002 |        0.535 |             0.420 |           0.347 |
| ALL      | Olympic_2010_armB | all families                   |        1148 |    6783 |    57.388 |       63.452 |     78.523 |        0.133 |             0.132 |           0.180 |
| ALL      | Olympic_2010_armB | members of multi-member groups |        1008 |    6440 |    48.895 |       54.959 |     74.383 |        0.104 |             0.103 |           0.164 |
| ALL      | Royer_NW_WA       | all families                   |        9504 |   66517 | 15631.831 |    16021.718 |  13269.206 |        0.463 |             0.463 |           0.525 |
| ALL      | Royer_NW_WA       | members of multi-member groups |        4608 |   35021 |  4297.758 |     4687.645 |   6284.821 |        0.416 |             0.419 |           0.542 |
| ALL      | Royer_SW_WA       | all families                   |        1898 |    7833 |   112.522 |      121.430 |    161.238 |        0.200 |             0.201 |           0.278 |
| ALL      | Royer_SW_WA       | members of multi-member groups |         832 |    3780 |    13.759 |       22.666 |     74.700 |        0.181 |             0.182 |           0.350 |

Files: `results/step2_groups/` — groups.csv, members.csv (δ, alignment quality), align_station.parquet (τ), setup/ (group families / truth / scan days), score/ (by_group, summary, sep_curve, tiers, HF cuts, singleton check), attrib/ (per-day cc), attrib_pairs.parquet (every candidate × member with P), attrib_*.csv, attrib_expected_*.parquet, step1/ (group P(real)). Templates `results/templates/*_step2grp.*`; MF `results/mf_step2grp`, `results/mf_step2grp_c34`.

![step2](../figures/step2_groups.png)
