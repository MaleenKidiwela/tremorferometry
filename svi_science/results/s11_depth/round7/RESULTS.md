# S11 round 7: front timing, per-episode zones with E4, joint arrival-time + coda, the WA large-episode hint

Written Tue 2026-10-13, ~9:50 am PT. Rules: `../PLAN.md` §11 (written 9:26 am PT before any round-7 number) and
Amendment 1 (post-code-test, labelled: ballistic kernels rescaled to the observed S−P; slope p not read; the WA zone
numbers were seen in a code test). Queue 9:41–9:42 am PT, every stage rc 0 (an earlier 9:39 launch was stopped by PID
after the far-grid error; its ZONES_SVI_DONE output was kept). Decision log `../NOTES.md` (decisions 32–33 + result lines).
Tables `{svi,wa_nw}/tables/r7_*`; figures `figures/fig_r7_{front,zones,joint,hint}.png`. dv/v in %. No new download, no
new event-level measurement. Everything is provisional.

**Machinery checks**
- The generic E4 model reproduces round-5 E4 ΔD and σ to 1e-13 (both regions). The hint fit reproduces the round-6 headline
  ΔD (≤ 2e-13) and the round-6 WA duration contrast (+0.072 ± 0.026, p 0.013).
- Day-pooled values agree with stacking the same days (corr 0.81 SVI, 0.70 WA; median offset ≤ 0.004 %).
- S3 and own-detection front dates agree within 1 day for 95–97 % of the (family, major) pairs.

## In plain words

**1. Front timing: no front-locked change, and front timing does not bias ΔD.**
- Inside an episode, dv/v does not change between the first days after the front reaches a family (lag 0–1) and the later
  days (lag ≥ 3):

  | region | late − early | front-scrambled p | slope | 2σ bound |
  |---|---|---|---|---|
  | SVI | +0.003 ± 0.007 % | 0.80 | −0.001 ± 0.001 %/day | < 0.016 % (< 0.004 %/day) |
  | WA | +0.006 ± 0.012 % | 0.39 | +0.003 ± 0.002 %/day | < 0.030 % (< 0.007 %/day) |

- The bound is about 5–10 × tighter than Merlin's forecast (~0.03 % per lag-day).
- Same answer with S3-only fronts, own-detection fronts, near or far stations, and high or low deep sensitivity.
- **Bias on ΔD:** even at the 2σ upper slope, the different sampling lags of the majors shift ΔD by ≤ 0.09σ (SVI) and
  ≤ 0.39σ (WA) per major, ≤ 0.16σ on the composite. Rule (0.5σ): **no bias**.
- One column to note: against the between-ETS baseline, WA reads −0.04 ± 0.02 % (1.8σ) late minus early. SVI reads
  +0.01 ± 0.01. Not significant.

**2. Zones with E4, one major at a time: no zone differs, pooled or in any single episode. Slipped minus not slipped
cannot be measured well yet.**
- **Pooled contrasts:**
  - SVI (2 zones): +0.00 ± 0.25 %, p 0.998.
  - WA (3 zones): +0.05 / −0.04 / −0.10 ± 0.08–0.10 %, p 0.47.
  - 2σ bounds: < 0.49 % (SVI), < 0.25 % (WA).
- **Per episode:**
  - SVI: smallest p 0.39.
  - WA: five majors have p 0.001–0.004 (2008-04, 2011-07, 2024-01, 2024-08, 2025-09). None passes the Bonferroni line
    (0.0004).
- **Why WA is not read as a zone signal.** The WA per-major χ² sum is high (109 on 48 dof), but the WA injection rule
  fails: a one-episode shallow pulse leaks into the zone contrasts at 0.75 (rule ≤ 0.35). So the WA per-episode zone
  scatter is most likely unmodelled shallow change.
- **SVI gates:** recovery 1.05, cross-zone leak 0.004, shallow leak 0.35 — all pass.
- **Controls:** pass in both regions (93–97 % of control contrasts within 2σ). Leaving out one site keeps 93–96 % within 2σ.
- **Slipped − not slipped:** S3 coverage says nearly every zone slips in every major. Episodes with both classes: SVI 2
  (S3) / 3 (family activity); WA 5 / 8.

  | region | slipped − not | permutation p | 2σ bound |
  |---|---|---|---|
  | SVI | −0.19 ± 0.38 % | 0.51 | < 0.96 % |
  | WA | +0.09 ± 0.36 % | 0.75 | < 0.82 % |

  Not slip-locked. Too few partial episodes to be useful.
- Caveat: the SVI per-zone composites are both +0.13 % (σ 0.20 / 0.08), while one D gives +0.001 ± 0.023. Splitting D by
  zone loosens the neighbour (burst) episodes.

**3. Joint arrival-time + coda: the ballistic rows separate the depth layers that the coda alone cannot. Nothing is
detected at any depth.**
- σ per major (shallow 0–5 / mid-crust 5–25 / deep 25–45 km):

  | region | coda only (S, M, D free) | ballistic only | joint | joint + far-pair dS |
  |---|---|---|---|---|
  | SVI | 0.14 / 0.34 / 0.085 | 2.7 / 1.1 / 0.48 | 0.11 / 0.072 / 0.055 | 0.10 / 0.048 / 0.034 |
  | WA | 0.14 / 0.20 / 0.091 | 5.7 / 2.3 / 0.56 | 0.11 / 0.076 / 0.059 | 0.10 / 0.058 / 0.035 |

  - Ballistic adds most to the mid crust (σ ÷ 3–7) and adds to the deep (÷ 1.5; ÷ 2.5 with dS).
  - The coda carries the shallow term.
  - Injections: deep pulse 0.99 (SVI) / 0.91 (WA); leak of deep into the shallow and mid terms ≤ 0.08; shallow leak into
    deep ≤ 0.03.
- **Composites of the majors (joint):**

  | region | shallow | mid | deep |
  |---|---|---|---|
  | SVI | +0.03 ± 0.04 | +0.06 ± 0.03 | +0.000 ± 0.025 |
  | WA | −0.029 ± 0.011 (2.7σ; coda and ballistic agree in sign) | −0.01 ± 0.01 | +0.021 ± 0.011 (1.9σ) |

  - No group reaches 3σ, so nothing is detected.
  - The WA shallow value is the closest. It is −0.014 ± 0.009 with dS.
  - Deep 2σ: < 0.05 % (SVI), < 0.04 % (WA). With dS: < 0.05 / < 0.015 %.
- **Decade (≥ 2015 minus ≤ 2014), joint:**

  | region | deep | mid | shallow | deep 2σ bound |
  |---|---|---|---|---|
  | SVI | +0.015 ± 0.014 | −0.010 ± 0.030 | +0.03 ± 0.04 | < 0.044 % |
  | WA | −0.013 ± 0.026 | −0.06 ± 0.04 | +0.06 ± 0.05 | < 0.066 % |

  - Ballistic adds to the deep and mid decade terms (joint σ ≤ 0.8 × coda), not to the shallow.
- Source-locus nuisance column: the ballistic part loses its depth resolution, but the joint deep σ is still 0.07–0.08.
  Unscaled-kernel column: within 5 % of the headline.

**4. The WA large-episode hint holds inside WA but does not replicate in SVI. Verdict: WA-ONLY (pre-registered).**
- **WA, top third by duration** (full set +0.072 ± 0.026 %):

  | split | half 1 | half 2 | halves differ by |
  |---|---|---|---|
  | family | odd +0.04 ± 0.06 | even +0.05 ± 0.03 | 0.2σ |
  | station | A +0.06 ± 0.03 | B +0.08 ± 0.04 | 0.4σ |

  - All four halves are positive and agree. The tremor top third gives the same.
- **Not a stack-length effect.** Episode duration is barely correlated with events per cell (r 0.14). The duration slope
  keeps 2.4σ once log events is included.
- **SVI replication fails:**
  - duration top third −0.006 ± 0.028 % (p 0.80);
  - tremor +0.030 ± 0.029 (p 0.23; that one was already seen in round 6).
- **Reading:** a coherent WA-wide +0.07 % excess during long or tremor-rich WA episodes, seen by independent families and
  stations. It is absent in SVI. Still a regional hint, not a finding.

## Figures
- `figures/fig_r7_front.png`: dv/v against days since the front, shape and absolute, both regions.
  - Red: data, points coloured by weight.
  - Grey: front-scrambled draws.
  - Blue dashed: time-of-day scrambled control.
- `figures/fig_r7_zones.png`: zone ΔD per major and slipped − not, coloured by weight.
- `figures/fig_r7_joint.png`: σ per major and decade σ for shallow / mid / deep (coda, ballistic, joint).
- `figures/fig_r7_hint.png`: top-third contrasts, WA full, halves, SVI, coloured by weight; scrambled 95 % lines.

## What failed and next
1. **WA zone injection rule** (shallow leak 0.75 into contrasts). Next: add the E4 + smooth regional shallow term, or a
   per-episode shallow term constrained by the ballistic rows (task 3 shows the joint fit holds the shallow leak at 0.03),
   to the zone model.
2. **Slipped vs not slipped:** too few partial episodes. Next: use bursts (inter-ETS partial slip) as the slipped/not
   contrast, or zone labels from tremor.
3. **WA hint:** next check whether it is WA-wide shallow. The joint shallow composite is −0.03 % (2.7σ, opposite sign),
   so a joint per-episode split of the hint by depth group is the cheapest test.

## Files
- Scripts `svi_science/scripts/s11_r7_{common,front,zones,joint,hint,figures}.py`, queue `s11_r7_queue.sh`, resume
  block 12 in `scripts/resume_after_reboot.sh`.
- Verdicts:
  - `{svi,wa_nw}/tables/r7_front_verdict.json`, `r7_zones_verdict.json`, `r7_joint_verdict.json`;
  - `hint_verdict.json`.
- Code tests (not results): `{svi,wa_nw}/tables/test/`. The SVI front test there ran on the full SVI set
  (`--test 2000`).
