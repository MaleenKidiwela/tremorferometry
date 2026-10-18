# GATES — a reusable, tested family-acceptance pipeline (Gates 0-3) + the P3 rediscovery scorer

2026-09-17, lfe_picker, agent GATES.  Spec: `results/merlin_review_plan_2026-09-17.md` § "The gates,
quantified" and § P3, plus `PLAN.md` (Amendments).  Built and validated on the 2010/2012 truth data
(`../lfe_catalog_repair/cache20hz`) with `../lfe_extension`'s frozen Track A run as the proposal source.
Everything new lives in `lfe_picker/` with the prefix `g_`; `../lfe_extension`, `../lfe_catalog_repair`
and `../codex-LFE` were read ONLY and every borrowed routine was COPIED (`a5_env.py`, `align_stack`,
`g_locate.py` from `a5_locate.py`/`d3_stage4_locate.py`, Crosson from `o_vel1d.py` = `d3_vel1d.py`).
Nothing was deleted; nothing with a `p_*` / `o_*` / `p1` / `p2` prefix was touched.

---

## 0. Headline

| | |
|---|---|
| **Gate 1 power** (own-event proposals, withheld VGZ/SNB/B926) | **0.77** at N=20, **0.91** at N=50, **0.92** at N=100, **0.92** at N=200 |
| **Gate 1 false-pass** (3 negative controls, 150 proposals) | **1 / 150 = 0.7 %** (one moveout-shuffled proposal; 0/50 time-shifted, 0/50 dv/v) |
| **Gate 0 alias false-pass** | **6 / 50 = 12 %** with the adopted rule — **37 / 50 = 74 %** with Merlin's rule as literally written |
| **Gate 0 class (HF) veto** | fires on **6/6** HF-rich proposals (median ratio 5.8-7.7); **0/734** spurious firings |
| **Gate 2** | **untestable below N ≈ 100**; testable 21/62 at N=200, of which 10 pass; median location error 5.6 km |
| **Gate 3 power / false-pass** | **54/60** on catalog-v1 families; **0/60** random-time; **0/60** activity-covarying (Merlin F13 closed) |
| **Runtime** | **3.7 s median**, 9.5 s p90, 26 s max per proposal for Gates 0+1+2 on a warm cache (Gate 3 ≈ 0.1 s) |
| **Rediscovery scorer** | calibrated: the frozen template detector's own clusters score only **9/30** hidden families |

**The one-line verdict.** Gate 1 and Gate 3 are strong, cheap and well-calibrated and can screen
thousands of proposals. Gate 0 works only after one correction to Merlin's rule (below). Gate 2 is real
but has almost no power below ~150 members on the Vancouver Island network, so a P5 proposal must be
carried to N ≳ 150 before a location can be demanded of it.

---

## 1. What the pipeline is

`scripts/g_gates.py` takes a **Proposal** — a set of detection times on the ORIGIN axis plus a proposed
interface cell (lat, lon, dep; NaN depth → Slab2) — and a region/station configuration, and returns
per-gate statistics, all nulls, and PASS / FAIL / UNTESTABLE per gate and overall. It knows nothing
about templates, families or matched filters: any detector can feed it. Usage: `results/gates/USAGE.md`.

Predicted S at a station = `t_origin + ttS(Crosson 1-D, cell → station) + a5 station correction`
(`a5_station_corrections.csv`, the intraslab-earthquake terms, never Bostock hypocentres). Waveforms are
the frozen chain (1-8 Hz Butterworth-4 zero-phase at the native rate → 20 Hz), except the HF class test,
which must read native-rate traces (the 20 Hz cache stops at 10 Hz).

### Gate definitions, with the exact numbers in force

**Gate 0 — identity / class.**
* `min_dist_km` to every accepted family (configurable list; default = the 130 Bostock families, with a
  leave-one-out `drop_bkey` so a re-proposal of a known family is judged against the OTHERS).
* multi-site stack cc to the NEAREST accepted family's 12-s template, max over ±1 s, per station, against
  that station's **cross-family cc p95** recomputed here for all 15 stations
  (`results/gates/g_crossfam_cc.csv`: PGC 0.524, B010 0.482, B009 0.449, B011 0.420, SNB 0.402, B005 0.388,
  VGZ 0.382, B926 0.348, NLLB 0.342, LRIV 0.307, B003 0.300, B007 0.292, B006 0.259, B004 0.258, B001 0.219).
* HF class ratio E(10-19 Hz)/E(2-8 Hz) at PGC (B011 fallback) over the 12-s window from S-1, **median over
  members**; veto at ≥ 3.5 (Track A's calibrated cut; LFE core median 0.247, ours 0.29).
* **PASS** = (`min_dist ≥ 5 km` OR waveform-distinct) AND the HF veto does not fire.

**Gate 1 — withheld-station stack.** Withheld = `VGZ, SNB, B926` (three distinct sites), excluded from
detection for the whole run. N ≥ 20 usable members. Statistic
`10 log10 E(S-1..S+5) / E(S-12..S-6)` on the 3-C 20 Hz stack at the predicted S — **no alignment
refinement**, deliberately: an alignment fitted on the data would import the coherence the gate tests.
Nulls: **200** matched random-time stacks on the member days (require observed > null p99 AND excess
≥ +3 dB) and **100** moveout-shuffled stacks at interface cells whose predicted S at that station differs
by ≥ 5 s (require the true-moveout excess > the shuffled p95). Positive control:
`results/gates/g_gate1_controls.csv`, 3 stations × 26 Bostock families × N ∈ {20,50,100} = 234 tests;
if no control at matched N and comparable distance passes, the verdict is **UNTESTABLE**, not FAIL.
Gate 1 PASSes when ≥ 1 withheld SITE passes both nulls.

**Gate 2 — location with controls.** Per-station stacks, **cc-aligned** (`align_stack`, the e5r recipe a5
picked on), S on √(H1²+H2²) within ±3 s and P on Z within ±2.5 s of the prediction, envelope-onset + AIC,
accepted at peak/noise-rms ≥ 3; sites collapse by median (PGC cluster and LRIV cluster = one site each).
Requires ≥ 4 picked sites within 60 km and azimuthal gap ≤ 180°. Estimator = the a5/d3 cost function
(σ_S 0.30 s, σ_S-P 0.50 s) on the **2-D Slab2 interface grid** (124 × 124 km at 1 km, depth fixed on Slab2)
plus a **free-depth** 3-D sanity run. Controls: shuffled-pick null (200 permutations of the pick vector
among the selected sites; require observed RMS S-residual < 0.5 × the null median), synthetic floor at
that exact geometry (0.3 s S noise, ±0.5 s per-site statics, 50 realisations; bar ≤ 3 km),
leave-one-site-out spread ≤ 5 km, free depth within ±10 km of Slab2 and not at a bound.

**Gate 3 — tremor timing (soft).** `E_local` = (fraction of members on days with PNSN tremor within 25 km
of the cell) / (fraction of EVALUATED days that are such days); `E_far` the same for the 75-150 km annulus.
Bars: `E_local ≥ 3`, `E_local/E_far ≥ 1.5`, p < 0.01 from 1,000 circular shifts of ≥ 30 d of the member
series on the evaluated calendar. A proposal that passes Gates 0-2 and fails Gate 3 is returned as
`PASS_inter_ETS_candidate`, never as FAIL.

---

## 2. Two places the spec could not be implemented as written

**(a) Gate 0's waveform clause, as literally written, false-passes 74 % of aliases.** Merlin's rule is
"≥ 5 km from every accepted family, **or** multi-site stack cc below the per-station cross-family p95 at
≥ 2 sites". At the far, low-SNR sites a proposal's stack is noise, so its cc is *trivially* below chance;
almost every proposal therefore finds two such sites. Measured on the alias control (a neighbouring
family's own events proposed as a new family 1-2 km away, which **must** fail): the literal rule passes
**37/50**, and it also passes re-proposals of accepted families at distance 0.00 km. The adopted rule
counts only sites where the proposal's **own stack registers an arrival** (S-window excess ≥ 2 dB) and
requires the cc to be below chance at **≥ 2 registering sites and at more than half of them**. That rule
fails 44/50 aliases (false-pass 12 %) and keeps 100 % power on proposals ≥ 5 km from any accepted family.
Both verdicts are stored (`g0_literal` vs `g0`), so the decision is auditable.

**(b) "Stack HF ratio" is measured per member, not on the coherent stack.** At 10-19 Hz the per-member
alignment error is many half-periods, so the coherent stack is depleted in the HF band:
`results/gates/g_hf_stack_vs_member.csv` measures the member-median / stack-ratio bias at **2.6× for LFE
families** and 1.3× for HF-rich transients. The 0.247 median / 3.5 veto were calibrated by Track A on
*per-detection* ratios, so the member median is the statistic those numbers describe. Both are stored.

**(c) Gate 3's circular-shift null does not exist on a single ETS window.** A shift of ≥ 30 d needs a span
of ≥ 61 d; the 2012 development window is 52 d, and every day in it is a local-tremor day, so `E_local`
degenerates to 1. Gate 3 therefore returns UNTESTABLE there **by construction**, and it was validated
instead on the long span the gate is actually for: `lfe_extension` catalog v1, 2014-01..2026-09 (4,641 d).

---

## 3. Validation suite

Proposal sets (`results/gates/proposals_*.parquet`, built by `g_make_proposals.py` /
`g_make_g3_proposals.py`):

| set | what | n |
|---|---|---|
| `pos_bost_N{20,50,100,200}` | the 116 templated Bostock families, members = Bostock's own ETS-2012 event times | 116 / 109 / 102 / 62 |
| `pos_tpl_N{20,50,100}` | the same families, members = the FROZEN template detections of the run with VGZ/SNB/B926 **excluded from the network sum** (`e5r_detections_a3b_noVSB`, 12 s, T_ext 8.25, ETS days, own-matched) | 81 / 50 / 14 |
| `neg_moveout_N50` | (b-i) the same times at a WRONG cell (S offset differs ≥ 5 s at every withheld station) | 50 |
| `neg_timeshift_N50` | (b-ii) the same times shifted by a random ±300-3600 s | 50 |
| `neg_dvv_N50` | (b-iii) 50 causality-certified dv/v "families" (`../data/*_causality_cert.csv`) with their stored detection times (`../data/mf_*.csv`, 2012 ETS days) at their tremor-cell coordinates, Slab2 depth, **±30 s constant-offset search** (the dv/v time is not an origin time; the same search is applied to every null) | 50 |
| `neg_hf_N50` | HF-rich quiet-day detections (PGC HF ratio > 3.5) proposed as families — the positive control for the Gate-0 class clause | 6 |
| `alias_N50` | (c) a neighbouring family's own events proposed as a NEW family 1-2 km away | 50 |
| `g3_pos / g3_neg_random / g3_neg_activity` | catalog-v1 families 2014-2026; random-time controls; **activity-covarying** controls (members drawn only on days with tremor anywhere within 300 km) | 60 each |

## Table 1 — verdicts per proposal set

| set | n | N | g0_PASS | g0_FAIL | g0_UNTE | g1_PASS | g1_FAIL | g1_UNTE | g2_PASS | g2_FAIL | g2_UNTE | overall_PASS | runtime_s_med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pos_bost_N20 | 116 | 20 | 70 | 34 | 12 | 88 | 26 | 2 | 0 | 0 | 116 | 0 | 2.678 |
| pos_bost_N50 | 109 | 50 | 61 | 42 | 6 | 99 | 10 | 0 | 0 | 1 | 108 | 0 | 4.005 |
| pos_bost_N100 | 102 | 100 | 49 | 49 | 4 | 94 | 8 | 0 | 4 | 5 | 93 | 1 | 5.809 |
| pos_bost_N200 | 62 | 200 | 27 | 34 | 1 | 57 | 5 | 0 | 10 | 11 | 41 | 5 | 8.417 |
| pos_tpl_N20 | 81 | 20 | 41 | 36 | 4 | 63 | 18 | 0 | 0 | 1 | 80 | 0 | 2.672 |
| pos_tpl_N50 | 50 | 50 | 24 | 25 | 1 | 46 | 4 | 0 | 0 | 2 | 48 | 0 | 3.657 |
| pos_tpl_N100 | 14 | 100 | 5 | 9 | 0 | 14 | 0 | 0 | 4 | 0 | 10 | 1 | 5.130 |
| neg_moveout_N50 | 50 | 50 | 47 | 0 | 3 | 1 | 49 | 0 | 0 | 0 | 50 | 0 | 3.620 |
| neg_timeshift_N50 | 50 | 50 | 9 | 0 | 41 | 0 | 50 | 0 | 0 | 0 | 50 | 0 | 3.599 |
| neg_dvv_N50 | 50 | 50 | 25 | 0 | 25 | 0 | 50 | 0 | 0 | 0 | 50 | 0 | 11.276 |
| neg_hf_N50 | 6 | 44 | 0 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6.946 |
| alias_N50 | 50 | 50 | 6 | 44 | 0 | 46 | 4 | 0 | 0 | 2 | 48 | 0 | 3.637 |

## Table 2 — Gate 1 (withheld-station stack): power and false-pass rate

| set | N | n | n_testable | untestable | g1_pass | rate_of_testable | rate_of_all | med_excess_db | med_best_dist_km |
|---|---|---|---|---|---|---|---|---|---|
| pos_bost_N20 | 20 | 116 | 114 | 2 | 88 | 0.772 | 0.759 | 3.649 | 43.275 |
| pos_bost_N50 | 50 | 109 | 109 | 0 | 99 | 0.908 | 0.908 | 5.482 | 40.279 |
| pos_bost_N100 | 100 | 102 | 102 | 0 | 94 | 0.922 | 0.922 | 7.508 | 42.687 |
| pos_bost_N200 | 200 | 62 | 62 | 0 | 57 | 0.919 | 0.919 | 9.457 | 44.720 |
| pos_tpl_N20 | 20 | 81 | 81 | 0 | 63 | 0.778 | 0.778 | 4.235 | 48.432 |
| pos_tpl_N50 | 50 | 50 | 50 | 0 | 46 | 0.920 | 0.920 | 6.074 | 45.636 |
| pos_tpl_N100 | 100 | 14 | 14 | 0 | 14 | 1.000 | 1.000 | 7.205 | 36.825 |
| neg_moveout_N50 | 50 | 50 | 50 | 0 | 1 | 0.020 | 0.020 | 0.933 | 54.834 |
| neg_timeshift_N50 | 50 | 50 | 50 | 0 | 0 | 0.000 | 0.000 | 0.471 | 45.989 |
| neg_dvv_N50 | 50 | 50 | 50 | 0 | 0 | 0.000 | 0.000 | 0.343 | 57.096 |
| neg_hf_N50 | 44 | 6 | 0 | 0 | 0 |  | 0.000 |  |  |
| alias_N50 | 50 | 50 | 50 | 0 | 46 | 0.920 | 0.920 | 5.940 | 47.307 |

## Table 3 — Gate 0 (identity/class) by distance to the nearest ACCEPTED family (leave-one-out for the positives)

| set | nn_km | n | g0_pass | g0_pass_rate | literal_rule_pass_rate | med_frac_above | med_hf |
|---|---|---|---|---|---|---|---|
| pos_bost_N50 | 0-2.5 | 59 | 19 | 0.322 | 0.898 | 0.600 | 0.297 |
| pos_bost_N50 | 2.5-5.0 | 28 | 20 | 0.714 | 0.964 | 0.250 | 0.277 |
| pos_bost_N50 | 5.0-inf | 22 | 22 | 1.000 | 1.000 | 0.183 | 0.260 |
| pos_bost_N200 | 0-2.5 | 37 | 7 | 0.189 | 0.865 | 0.667 | 0.292 |
| pos_bost_N200 | 2.5-5.0 | 14 | 9 | 0.643 | 0.929 | 0.400 | 0.291 |
| pos_bost_N200 | 5.0-inf | 11 | 11 | 1.000 | 1.000 | 0.286 | 0.281 |
| pos_tpl_N50 | 0-2.5 | 28 | 6 | 0.214 | 0.893 | 0.750 | 0.230 |
| pos_tpl_N50 | 2.5-5.0 | 13 | 9 | 0.692 | 0.923 | 0.400 | 0.242 |
| pos_tpl_N50 | 5.0-inf | 9 | 9 | 1.000 | 1.000 | 0.200 | 0.268 |
| alias_N50 | 0-2.5 | 50 | 6 | 0.120 | 0.740 | 0.800 | 0.241 |

## Table 4 — Gate 2 (location with controls)

| set | N | n | testable | untestable | g2_pass | med_sites_60km | med_dist_to_proposal_km | med_null_ratio | med_synth_km | med_loso_km |
|---|---|---|---|---|---|---|---|---|---|---|
| pos_bost_N20 | 20 | 116 | 0 | 116 | 0 | 1.000 |  |  |  |  |
| pos_bost_N50 | 50 | 109 | 1 | 108 | 0 | 2.000 | 3.544 | 0.081 | 1.413 | 7.603 |
| pos_bost_N100 | 100 | 102 | 9 | 93 | 4 | 3.000 | 4.011 | 0.121 | 1.000 | 6.082 |
| pos_bost_N200 | 200 | 62 | 21 | 41 | 10 | 3.000 | 5.588 | 0.080 | 1.413 | 4.115 |
| pos_tpl_N20 | 20 | 81 | 1 | 80 | 0 | 2.000 | 4.287 | 0.053 | 1.996 | 5.376 |
| pos_tpl_N50 | 50 | 50 | 2 | 48 | 0 | 2.000 | 4.551 | 0.068 | 0.499 | 7.819 |
| pos_tpl_N100 | 100 | 14 | 4 | 10 | 4 | 3.000 | 5.588 | 0.050 | 0.706 | 1.823 |
| neg_moveout_N50 | 50 | 50 | 0 | 50 | 0 | 0.000 |  |  |  |  |
| neg_timeshift_N50 | 50 | 50 | 0 | 50 | 0 | 0.000 |  |  |  |  |
| neg_dvv_N50 | 50 | 50 | 0 | 50 | 0 | 0.000 |  |  |  |  |
| neg_hf_N50 | 44 | 6 | 0 | 0 | 0 |  |  |  |  |  |
| alias_N50 | 50 | 50 | 2 | 48 | 0 | 2.000 | 4.859 | 0.139 | 1.704 | 13.901 |

## Table 5 — Gate 3 (tremor timing, long-span validation on catalog v1 2014-2026)

| set | n | N | g3_pass | rate | med_E_local | med_E_far | med_ratio | frac_p_lt_0p01 |
|---|---|---|---|---|---|---|---|---|
| g3_pos | 60 | 300 | 54 | 0.900 | 8.126 | 1.424 | 5.480 | 0.950 |
| g3_neg_random | 60 | 300 | 0 | 0.000 | 0.954 | 1.005 | 0.932 | 0.000 |
| g3_neg_activity | 60 | 300 | 0 | 0.000 | 2.510 | 2.503 | 1.001 | 1.000 |

## Table 6 — Gate 1 POSITIVE-CONTROL table (known Bostock families, their own event times)

| sta | N | n | pass_ | rate | med_excess_db | med_dist_km |
|---|---|---|---|---|---|---|
| B926 | 20 | 26.000 | 12.000 | 0.462 | 3.016 | 49.002 |
| B926 | 50 | 26.000 | 19.000 | 0.731 | 4.318 | 49.002 |
| B926 | 100 | 26.000 | 21.000 | 0.808 | 6.030 | 49.002 |
| SNB | 20 | 26.000 | 10.000 | 0.385 | 1.930 | 60.189 |
| SNB | 50 | 26.000 | 15.000 | 0.577 | 4.254 | 60.189 |
| SNB | 100 | 26.000 | 16.000 | 0.615 | 4.995 | 60.189 |
| VGZ | 20 | 26.000 | 7.000 | 0.269 | 2.139 | 38.190 |
| VGZ | 50 | 26.000 | 15.000 | 0.577 | 3.662 | 38.190 |
| VGZ | 100 | 26.000 | 19.000 | 0.731 | 5.027 | 38.190 |

## Runtime

median per proposal (gates 0+1+2): **4.0 s**, p90 9.2 s, max 25.8 s (n = 740 proposals)


---

## 4. Reading the tables

**Gate 1 is the workhorse.** Power rises 0.77 → 0.91 → 0.92 between N = 20, 50 and 100 and the median
excess rises 3.6 → 5.5 → 9.5 dB; the template-detection proposals behave the same as the Bostock-time
proposals (0.78 / 0.92 / 1.00), i.e. the ~60 % per-event purity of a real detection list costs the gate
nothing. Track A 3b's ETS control (76/81 registering at N ≈ 62, +7.6 dB) is reproduced. False passes are
1 in 150 across three independent negative mechanisms, and the one failure is a moveout-shuffled proposal
whose wrong cell happened to leave a real arrival inside the 6-s S window. **The 50 dv/v families fail
0/50 even with a ±30 s offset search** — the strongest form of the "they are not sources" result, because
the search was given every chance to find an arrival and the nulls were charged the same freedom.

**Gate 1 power does not track distance.** Per-station pass rate is 0.77 (< 30 km), 0.52 (30-45),
0.50 (45-60), 0.59 (60-80), 0.87 (> 80 km) — station quality and family radiation dominate geometry over
this aperture (`results/gates/g_gate1_per_station.csv`). Every positive proposal had at least one withheld
station within 60 km, so the ">60 km" blind spot never binds in this region; in a new region it must be
checked with the `all_withheld_far` flag, which is reported per proposal.

**Gate 0's power is a function of how close the nearest accepted family is, as it must be.** 100 % of
proposals ≥ 5 km from any accepted family pass; 64-71 % at 2.5-5 km; 19-32 % below 2.5 km. The last number
is the gate working, not failing: at < 2.5 km on this network a stack is genuinely not distinguishable
from its neighbour's (cross-family cc reaches 0.52 at PGC), and the honest answer is "not demonstrably
new". 56 % of the Bostock families are in that regime, which is the real ceiling on any margin-wide
discovery claim in the north.

**Gate 2 is N-limited, not geometry-limited.** With perfect picks the network geometry would admit
113 of 130 families (≥ 4 sites within 60 km and gap ≤ 180°, `g_gate2_geometry_ceiling.csv`), but the
*pick* yield at 20 Hz needs a big stack: 0/116 testable at N=20, 1/109 at N=50, 9/102 at N=100 and 21/62
at N=200. Of the 21 testable at N=200, 10 pass all four controls; the failing check is almost always
leave-one-site-out (median 4.1 km against a 5 km bar). Where Gate 2 does return a location, the controls
are comfortable (shuffled-pick null ratio 0.08 against a 0.5 bar; synthetic floor 1.4 km against 3 km) but
the **actual** distance from the proposal to the located cell is 5.6 km (median, N=200) — larger than the
a5 relocation's 3.6 km, which used full-N template stacks and only the 14 best-geometry families. The
controls therefore certify that a location is data-driven; they do **not** certify it to 5 km.

**Gate 3 closes Merlin F13.** The activity-covarying control passes the circular-shift test 60/60
(p < 0.01 every time) and passes the E_local ≥ 3 bar in many cases, yet **0/60** survive the gate, because
`E_local/E_far` is 1.00 for them and 5.5 for real families. The local-vs-distant contrast, not the shift
null, is the load-bearing clause.

**Runtime.** 3.7 s median / 9.5 s p90 / 26 s max per proposal for Gates 0+1+2 on a warm day cache
(632 proposals), Gate 3 ≈ 0.1 s. Cost is dominated by (i) native-trace reading for the HF class test and
(ii) the 300 null stacks per withheld station; both are pure NumPy on cached day arrays, so the H100 was
not needed. At this rate P5 can screen ~1,000 proposals/hour per process, and the sets here ran 11-way in
parallel on the 32-CPU cap. Sort proposals by time so the `DayStore` LRU works.

### Example figures
`figures/g_example_PASS_posb200_f1.png` (PASS: Gate 0 PASS, Gate 1 PASS, Gate 2 PASS),
`figures/g_example_FAIL_alias_alias_f10.png` (FAIL: a real source — Gate 1 passes at +7 dB — rejected by
Gate 0 as an alias of the family 2.0 km away),
`figures/g_example_UNTESTABLE_posb20_f0.png` (UNTESTABLE: Gates 0 and 1 pass at N=20, Gate 2 has no picks).

---

## 5. The P3 rediscovery scorer

`scripts/g_rediscover.py`. A hidden family counts as rediscovered when some proposed cluster has a
centroid within 5 km, matches ≥ 30 % of that family's Bostock events one-to-one at ±2 s inside the
evaluated windows, has member purity ≥ 50 % against it, and is ≤ 50 % explained by a KEPT family's events.
Also reported: the spurious-split rate of the kept families (a kept family appearing as ≥ 2 clusters, each
explaining ≥ 20 % of its events at ≥ 50 % purity).

**Calibration (this is the number P3 must be read against).** Feeding the scorer the best cluster set that
exists on disk — the frozen template detector's own per-family detections at T_ext, placed at the exact
Bostock position, with the 30 holdout families of `results/holdouts.json` hidden — gives

| stratum | n | rediscovered | rate | median one-to-one recall | median purity |
|---|---|---|---|---|---|
| isolated (NN ≥ 5 km) | 10 | 5 | 0.50 | 0.34 | 0.55 |
| mid (2.5-5 km) | 10 | 1 | 0.10 | 0.24 | 0.49 |
| close (< 2.5 km) | 10 | 3 | 0.30 | 0.29 | 0.59 |
| **all** | **30** | **9** | **0.30** | | |

spurious-split rate **0.01** (1 of 100 kept families). The binding clause is the ≥ 30 % one-to-one recall,
which is **above the frozen template detector's own recall** (0.18-0.20 per family): a detector that
matched the templates exactly would still score ~0.3, not 1.0. **Merlin's P3 bar of "≥ 2/3 of the 10
isolated hidden families" therefore demands a per-family recall the best detector on disk does not have**,
and P3 should be reported against this 0.50 ceiling as well as against the bar.

---

## 6. Known blind spots (Merlin's list, plus what the validation added)

1. **Intraslab multiplets with interface-consistent moveout.** Untested and untestable here: nothing in
   the suite is a real non-LFE repeating source at an interface-consistent location. Gate 1 would pass one;
   the only defences are the Gate-0 HF class veto (which does fire on HF-rich transients, 6/6) and Gate 2's
   free-depth check (±10 km of Slab2). A real multiplet 15 km below the interface with a similar moveout
   would be caught only if Gate 2 is testable, i.e. only at N ≳ 150.
2. **Aliases.** Reduced from 74 % to 12 % false-pass by the corrected Gate-0 rule, not to zero. The
   residual 6/50 are proposals whose own stack registers at only 2-3 sites, where the majority test is
   weak. Recommended P5 practice: treat a PASS at `n_sites_registering < 4` and `min_dist < 3 km` as a
   candidate, not a family.
3. **Families whose withheld stations are all > 60 km.** Never binding in this region (0 positives), but
   the flag `all_withheld_far` is computed per proposal and must be honoured in Puget / Oregon, where the
   withheld pair has to be chosen per region.
4. **Activity-covarying false alarms.** Gate 3's contrast clause closes the version of this that Merlin
   flagged (0/60 activity-covarying controls pass). It does **not** close a false-alarm process driven by
   *local* tremor specifically — e.g. a detector that fires on tremor energy at the proposed cell. Gate 1's
   moveout-shuffled null is the defence there, and it is the null that 49/50 wrong-cell proposals fail on.
5. **Gate 2 accuracy is not certified by Gate 2's controls.** Median 5.6 km proposal-to-solution distance
   at N=200 while all four controls pass. The controls test that the location is data-driven and stable,
   not that it is accurate to 5 km.
6. **Gate 0's accepted-family list is the science.** Everything Gate 0 says is relative to that list; in a
   region with no catalogued families, Gate 0 degenerates to the HF class test plus "≥ 5 km from nothing",
   i.e. it passes everything. In P5, the accepted list must grow as families are accepted, and proposals
   must be re-screened against it.
7. **The Gate-1 control table is regional.** `g_gate1_controls.csv` is built from Bostock families and the
   three Vancouver Island withheld stations. A new region needs its own control table, or Gate 1 there can
   only return FAIL/UNTESTABLE honestly by distance-matching, which it cannot do without one.
8. **2008/2009 and the far south** are outside this validation entirely (no CN data pre-2010; no roster
   station for the southern families — P0-O § O4).

---

## 7. Files

| | |
|---|---|
| pipeline | `scripts/g_gates.py`, `g_common.py`, `g_locate.py`, `a5_env.py` |
| P3 scorer | `scripts/g_rediscover.py` |
| calibration | `scripts/g_crossfam.py` → `results/gates/g_crossfam_cc.csv`; `scripts/g_controls.py` → `results/gates/g_gate1_controls.csv` |
| validation | `scripts/g_make_proposals.py`, `g_make_g3_proposals.py`, `g_validate.py`, `g_report.py`, `g_hf_check.py`, `g_figures.py` |
| results | `results/gates/g_val_<set>.{csv,json}`, `g_val_all.csv`, `g_table_*.csv`, `g_tables.md`, `g_gate1_per_station.csv`, `g_gate2_geometry_ceiling.csv`, `g_hf_stack_vs_member.csv`, `g_rediscovery_demo*.csv` |
| figures | `figures/g_example_{PASS,FAIL_alias,UNTESTABLE}_*.png` |
| usage | `results/gates/USAGE.md` |
