# Step 1c′ — why the oracle ceiling returned C_30 = 0, and what it does and does not mean

Agent H, 2026-09-20 ~18:30Z. Queue `z_step1cp_queue.sh` (pid 1165933) stopped at stage 6 with
`results/STEP1CP_FAILED`: *"ABANDON at oracle ceiling: C_30 = 0 (< 3) — the frozen chain cannot
recover even perfectly seeded families on this bed."* The coordinator ruled that this must not
stand as a bed verdict without evidence, because it contradicts two prior measurements. It does
not stand. **The bed, the cache, the time base and the templates are sound. What C_30 = 0 measures
is that a ≤ 40-event pass-0 stack cannot clear the frozen survival bar — and the real chain's
pass-0 stacks are thinner still.**

Everything below is measured on the frozen 1c′ inputs (the same station-hour veto table and the
same cache); nothing was re-downloaded and no frozen constant was retuned. Raw tables are in
`results/step1c_prime/evidence/`.

---

## a. Time base — real, several seconds, and NOT the cause

`scripts/z_diag_oracle.py` part a, on the 30 hidden families **plus 10 kept controls**: per station,
the family's events are cut at `t_origin + fam_tt`, rms-normalised and stacked; the envelope peak is
located within ± 10 s of the prediction.

| quantity | value |
|---|---|
| families whose peak offset exceeds 1.5 s at half or more of their sites | **77.5 %** (31/40) |
| median over families of the median \|offset\| | **2.79 s** |
| p90 of that quantity | **6.04 s** |
| median stack S-SNR, no alignment | 0.855 dB |
| median stack S-SNR, align ± 1.5 s (frozen) | 0.855 dB |
| median stack S-SNR, align ± 6 s | 0.952 dB |
| median envelope-peak height (MAD), no / ± 1.5 s / ± 6 s | 3.43 / 3.39 / 3.46 |

So the 2005 travel-time prediction is systematically off by ~3 s (p90 6 s) for many families —
Merlin's family-45 observation generalises. **But widening the alignment search buys nothing**
(0.855 → 0.952 dB, peak height unchanged): `align_stack` aligns *events to each other at one
station*, and a per-family offset against the *prediction* is common to all of a family's events, so
there is nothing for it to fix. The frozen ± 1.5 s window is therefore **not** the defect, and it was
not widened.

The full D8 stage (50 busiest kept family-hours) agrees: per-site median peak offsets 0.35–3.1 s,
**11 of 12 sites within ± 3 s**, +30 s control excess **1.0009**, chain ratio **1.028** (frozen 20 Hz
1–8 Hz vs raw 100 Hz). D8 verdict ACCEPT.

## b. Template sanity — the templates are correct; the onset picker is not the bottleneck

`scripts/z_diag_onset.py` rebuilt the oracle templates under four cut rules and grew them over the
seed days (8 isolated families, T = 5):

| cut rule | mean self-detection | mean recall | median purity |
|---|---|---|---|
| `blind` — the frozen argmax in −4…+16 s (what 1c′ ran) | 0.054 | 0.0096 | 0.70 |
| `narrow` — argmax within ± 4 s of the predicted S | 0.142 | 0.0143 | 0.00 |
| `ttcut` — no search, cut exactly at the predicted S | 0.135 | 0.0139 | 0.34 |
| `peakcut` — argmax within ± 10 s of the predicted S | 0.129 | 0.0128 | 0.27 |

An oracle cut that uses the family's own travel times performs the same as the blind argmax. The
component mapping, rotation, per-file `t0` and template start are sound, and the strong templates
detect the events they were built from (`self_detect_frac` 0.38–0.68 for H0008/H0029) at
**purity 0.70–0.96** against the catalogue — far above the ~8 % chance rate in those dense hours.

## c. Recall — ~1 % over three days, but 24 % on the family's own day

`scripts/z_diag_netmad.py` computes the detector's **full `net_mad` series** for the oracle template
over the seed day and samples it at every catalogued event of that family (± 1 s), against 4,000
random times:

| family | events that day | `net_mad` median | null median | ≥ 5 | passes **all** clauses¹ | > null p99 |
|---|---|---|---|---|---|---|
| 8 | 358 | 2.46 | 2.05 | 11.5 % | 11.2 % | 14.8 % |
| 12 | 83 | 3.12 | 2.05 | 43.4 % | 43.4 % | 45.8 % |
| 29 | 192 | 2.24 | 1.88 | 19.8 % | 19.8 % | 19.8 % |
| 62 | 129 | 2.57 | 2.04 | 27.9 % | 27.9 % | 31.0 % |

¹ `net_mad ≥ 5` **and** ≥ 2 sites at site-z ≥ 3 **and** ≥ 4 sites observable — the `n_site_z3 ≥ 2`
coincidence clause costs nothing here (`frac_ALL` = `frac_nm5` to three decimals).

**So a correct template detects 11–43 % (median 24 %) of its family's catalogued events on that
family's own day.** Over three ETS days the recall falls to ~1 % because the other days are far
poorer. Stack coherence versus N (5 → 160 events) is flat (peak MAD 4.5 → 3.9, S-SNR 2.6 → 1.3 dB):
the stack does not improve with more events from the *same* hour, consistent with tremor being made
of the same family's other events.

## d. Why nothing survived — a burst guard applied to a seed chosen as a burst (defect, fixed)

The oracle's own survival table (`results/step1c_prime/oracle/p0_survival_frozenrule.csv`) shows 28
templates with **20–46 hits on 2–8 distinct days and scramble ratios up to 17** — and 0 survivors:

| failing clause | templates |
|---|---|
| `frac_near_seed > 0.20` (burst guard) | **21 / 28** |
| `n_indep < 20` | 17 / 28 |
| `n_days < 3` | 8 / 28 |
| `scramble_ratio < 3` | 10 / 28 |

Merlin's oracle protocol seeds each family from its **busiest hour**, i.e. deliberately from its
biggest burst, and the frozen survival rule then vetoes a template for detecting events near its own
seed window. That is a protocol self-contradiction, and it was fixed as a **defect**
(`z_oracle.oracle_survival`, documented in `results/freeze_step1c_prime_v2.json`): for the oracle
only, the burst guard is *reported, not applied*, and independent hits and distinct days are counted
**outside** the seed window — strictly stricter than the frozen rule. The fix cannot manufacture a
PASS: the ceiling enters the bars as 0.5 × C_30 and 0.4 × C_iso, so a higher ceiling only makes the
1c′ bar harder.

**Under the fix the count is still 0 / 28** — the binding clause simply moves to `n_indep ≥ 20`:
outside the seed window the templates yield only **1–13** independent hits over the 44-day bed.
So C_30 = 0 is *not* an artefact of the burst guard alone.

## e. The decisive control — the bed DOES support survival, with a richer stack

The same oracle templates rebuilt from **400** catalogued events spread over the whole window
(instead of ≤ 40 from one hour), grown over all 44 days
(`results/step1c_prime/evidence/e_allevent_stacks.log`):

| family | stacked | hits | independent | distinct days | recall | purity |
|---|---|---|---|---|---|---|
| 8 | 400 | 135 | **125** | **20** | 0.040 | 0.28 |
| 10 | 400 | 97 | **90** | **18** | 0.054 | 0.43 |
| 12 | 400 | 78 | **65** | **14** | 0.081 | 0.60 |
| 29 | 400 | 30 | **28** | **9** | 0.039 | 0.57 |
| 15 | 400 | 731 | 550 | 25 | 0.016 | 0.01 (over-triggering) |
| 62 | 400 | 959 | 687 | 24 | 0.015 | 0.01 (over-triggering) |

**Four of six families clear the frozen survival bar (≥ 20 independent hits on ≥ 3 days) outright**
once the pass-0 stack is rich enough. Two over-trigger and would be caught by the scramble null.

---

## Conclusion

1. **"The 2005 POLARIS bed is unusable" is NOT supported.** The chain, the cache, the time base, the
   channel mapping and the templates are all sound, and with a rich stack the bed yields 28–125
   independent detections on 9–20 days per family.
2. **What C_30 = 0 actually measures** is that a pass-0 stack of ≤ 40 events from one hour reaches
   only 1–13 independent hits outside its seed window — below the frozen `≥ 20 independent hits on
   ≥ 3 days` survival bar. The limit is the **survival gate being applied at pass 0, before the
   iterations that would enrich the stack**, not the data.
3. **Consequence for the 1c′ discovery run.** Its pass-0 templates come from 2-node pair seeds —
   thinner than the 40-event oracle — so the frozen chain cannot do better than the ceiling it just
   failed. Running stages 7–10 as frozen would return ~0 survivors for a reason that has nothing to
   do with the seed rule Merlin's repair changed.
4. **Not acted on unilaterally.** Moving or staging the survival gate (e.g. survive after iteration
   1, or seed pass-0 growth without the ≥ 20-hit gate and apply it post-iteration) is a change to
   the chain, i.e. a second repair, which Merlin's ruling forbids without his say-so. The queue is
   left stopped at `STEP1CP_FAILED`; stages 1–5 are complete and guarded, so a decision only needs
   stage 6 onwards re-run (~6–8 GPU-h).

**Question for Merlin:** given (e), is the correct ceiling protocol a richer oracle stack, and is
the pass-0 survival gate — not the seed rule — the defect his F3 diagnosis was pointing at?

### Files
`results/step1c_prime/evidence/` (a_offsets_by_site, a_per_family_offsets, b_onset_rules,
bc_templates, c_netmad_at_catalogued_events, c_clause_budget, c_stack_vs_N, e_allevent_stacks.log);
`results/step1c_prime/oracle/p0_survival_frozenrule.csv` (the frozen-rule survival table);
`results/step1c_prime/oracle_reference_2005.{csv,json}`; `results/freeze_step1c_prime_v2.json`;
scripts `z_diag_oracle.py`, `z_diag_onset.py`, `z_diag_netmad.py`.
