# S11 round 5: the prior test, a smooth regional shallow term, the WA yearly cycle, 2019-11 per instrument, kernel band

Written Mon 2026-10-12, ~11:15 pm PT. The rules are in `../PLAN.md` §9, written ~10:25–10:30 pm PT before any round-5
number. Amendment 1 (~10:35 pm PT) is **post-code-test**: it was written after a 40-series code test whose numbers were seen.
The queue ran 10:41–10:54 pm PT and every stage had rc = 0. The decision log is `../NOTES.md` (decisions 27–28 and the
result lines).

**Where things are**
- Tables: `{svi,wa_nw}/tables/`.
- Figures: `figures/fig_r5_{forms,dD,q1,pgc,band}.png`.
- Selection: `selected.json`; smooth-Sc adoption: `sc_adoption.json`.

**Conventions**
- dv/v is in %.
- ΔD = the deep-weighted change of one slow-slip episode, minus its neighbours within ±3 yr.
- L-SEL = two3_20, the kernel chosen by the lapse fit. It is the round-5 headline (decision 27).
- NOM = ell40.
- "Scrambled control" = each station's data shifted to another time of the same day.
- Everything here is provisional.

**Checks of the machinery**
- The round-5 model class reproduces round-4 H and HnoSc exactly (max |Δ| 8e-17).
- The re-simulated kernels reproduce the stored 2–5, 2–3.5 and 3.5–6 s coefficients exactly (0.0).
- The re-measured 2–5 s graph reproduces the cell values (max |Δ| 7e-9 %).

## In plain words

**1. The prior test passes for the first time.** The fix was one receiver-type shallow term per station (or per site),
no ridge on any term, and a deep episode term with no slow or yearly part. This is form **E4**.
- **The prior test.** Adding a 1 % or 2 % ridge now moves every major and the composite by 0.034–0.072σ (rule < 0.1σ).
  This holds in both regions, under both kernels.
- **Precision got better, not worse:**

  | | SVI | WA |
  |---|---|---|
  | σ per major, L-SEL | 0.056 % (15/15 ≤ 0.1) | 0.065 % (23/24 ≤ 0.1) |
  | σ per major, round-4 HnoSc | 0.082 % | 0.124 % |

- **Composites (L-SEL):** SVI +0.001 ± 0.023 %, WA −0.004 ± 0.009 %.
- **A yearly receiver cycle common to every station no longer leaks** (T4, failed in every round-4 form). Composite
  leak is 0.04–0.05 and the median major leak is 0.12–0.19 (rule ≤ 0.35).
- **One item fails the task-1 pass rule:** WA under L-SEL recovers a deep pulse at 0.75 (rule 0.8–1.1).

**2. A smooth regional shallow term works without the trade-off.**
- The term is a B-spline with 2-yr knots, made orthogonal to the per-station terms. It meets the goal in all four
  region × kernel combinations.
- Leak 0.04–0.11 (goal ≤ 0.35). σ per major 0.058–0.092 %, below round-4 HnoSc in every combination.

**3. The WA yearly cycle is not deep in any plausible kernel. The verdict is MIXED (near + far, two3_20).**
- The cycle survives a per-station receiver term, but its peak moves by 64 days.
- Its lapse ratio is above 1 (1.40 ± 0.21), where a deep change predicts 0.69.
- The PGC-site stations carry their own cycle (3.3σ) with an undecided lapse ratio.

**4. SVI 2019-11 appears at all four co-located instruments. It is SITE-WIDE, not one instrument, and no metadata
changed.**
- In the new model it is ≤ 1.2σ (L-SEL).
- Its lapse profile falls with lapse time, which a deep change would not do.
- The most likely reading is a local change at the PGC site.

**5. The Q3 conclusion holds across the kernel band. The Q1 call does not, except "not deep".**
- Q3: the composite is within 2σ of zero for every plausible kernel. Size bound 0.08–0.11 % (SVI), 0.05–0.07 % (WA).
- Q1: the Q1/σ ≥ 3 statement and the exact verdict change between two3_20 and two5_40. "Not deep" holds everywhere.

## Task 1: the prior test (PLAN §9.1)

Forms:
- E1 = Ua + D0; E2 = Ua + Dp; E3 = U1 + D0; E4 = U1 + Dp.
  - Ua = a_s × basis per unit; U1 = coefficient 1 per unit; Dp = D(e) orthogonal to the slow + annual basis.
- E2b, E4b fold the > 45 km layer into the deep term (Amendment 1).
- Every form has M and B and no ridge (exact basis).

**Selection** (score = passes summed over 4 region × kernel combinations): E4 31, E2 27, E4b 27, E1 26, E2b 25, E3 25.
**Selected: E4.**

**E4 against the rules:**

| rule | SVI NOM | SVI L-SEL | WA NOM | WA L-SEL |
|---|---|---|---|---|
| σ per major (≤ 0.1 count) | 0.072 (13/15) | **0.056 (15/15)** | 0.081 (21/24) | **0.065 (23/24)** |
| composite | −0.015 ± 0.028 | **+0.001 ± 0.023** | −0.014 ± 0.008 | **−0.004 ± 0.009** |
| G-prior: ridge 1 % / 2 %, max over majors + composite (< 0.1σ) | 0.055 ✔ | 0.072 ✔ | 0.045 ✔ | 0.034 ✔ |
| T1 deep pulse (0.8–1.1) | 0.92 ✔ | 0.90 ✔ | 0.82 ✔ | **0.75 ✘** |
| T2 common shallow leak (≤ 0.35) | 0.27 ✔ | 0.09 ✔ (fixed) | 0.28 ✔ | 0.04 ✔ (fixed) |
| T3 one-site leak (≤ 0.35) | ≤ 0.19 ✔ | ≤ 0.08 ✔ | ≤ 0.12 ✔ | ≤ 0.11 ✔ |
| T4 receiver annual, all stations (comp / major) | 0.04 / 0.19 ✔ | 0.05 / 0.15 ✔ | 0.03 / 0.15 ✔ | 0.04 / 0.12 ✔ |
| T5 shallow annual → composite | 0.003 ✔ | 0.010 ✔ | 0.005 ✔ | 0.008 ✔ |
| G-ctrl (B2 control) | ✔ (+0.21 ± 0.13) | ✔ | ✔ | ✔ |
| B1 reference-path control (column) | ✔ | ✘ 11/14 majors within 2σ (composite +0.03 ± 0.10) | ✔ | ✘ 20/23 (composite −0.01 ± 0.01) |
| G-loso | ✔ | ✔ | ✔ (worst EQUIL 0.92) | ✔ |
| G-fam | ✔ | ✔ | ✔ | ✔ |
| G-site (≤ 50 %) | ✘ PGC 65 % | ✘ PGC 65 % | ✔ B013+B943 31 % | ✔ |
| majors at ≥ 3σ | 0 | 0 | 1 (2023-03, 3.0σ) | 0 |

**What fixed the prior test.** The per-group diagnostic puts the remaining sensitivity in the unit terms (≤ 0.06σ). M and
B are below 0.04σ.

**Other forms and columns:**
- E2 (kernel-shallow units) fails the prior test, at 0.13–0.37σ.
- The round-4 forms fail it at 0.29–0.57σ.
- E2b and E4b, with the > 45 km layer folded into the deep term, collapse the WA deep-pulse recovery to 0.29. They are
  rejected.
- The calibrated-penalty column (Pcal) fails ×2 / ×0.5 at 0.27–0.34σ. The exact basis is the better route.

**Task-1 pass rule: not met, on one item only.** WA L-SEL T1 is 0.75.
- Obstacle: with only 27 WA episodes, the eight constraints that make D orthogonal to the slow + annual basis remove part
  of a one-episode pulse.

**WA 2023-03.** Under NOM it is +0.243 ± 0.081 % (3.0σ), and every gate passes in that column. It is not the headline:
- under L-SEL it is +0.180 ± 0.069 (2.6σ);
- over the plausible kernels it is 2.3–2.6σ;
- with the smooth regional shallow term it is 1.2–2.1σ.

So it is a candidate, not a claim.

## Task 2: a smooth regional shallow term (PLAN §9.2)

The form is E4 + Sc_spl: cubic B-splines with knots every k = 2, 3 or 4 yr, orthogonal to the unit basis (13 / 10 / 8
columns kept in SVI). The goal is met for every k in all four combinations. **Adopted: k = 2 yr** (smallest mean L-SEL σ).

| E4 + Sc 2-yr | SVI NOM | SVI L-SEL | WA NOM | WA L-SEL |
|---|---|---|---|---|
| σ per major (round-4 HnoSc) | 0.086 (0.111) ✔ | 0.058 (0.082) ✔ | 0.092 (0.128) ✔ | 0.064 (0.124) ✔ |
| T2 leak (≤ 0.35) | 0.092 ✔ | 0.057 ✔ | 0.109 ✔ | 0.043 ✔ |
| composite | −0.037 ± 0.032 | −0.002 ± 0.024 | −0.004 ± 0.011 | +0.007 ± 0.009 |
| G-prior | 0.102 (just over 0.1) | 0.079 ✔ | 0.067 ✔ | 0.082 ✔ |
| T1 / T3 / T4 / T5 / G-ctrl / G-loso / G-fam | ✔ all | ✔ all | ✔ all | T1 0.75 ✘, rest ✔ |

- **What it buys.** The smooth term costs almost no precision.
- **What it changes.** It takes up the WA 2023-03 excess (2.6σ → 1.6σ, L-SEL). So that candidate is partly a regional 2-yr
  shallow wiggle.
- **Known logging bug.** `sc_adoption.json` "failing" lists every row, not only failing ones. Its "goal_all" is correct.

## Task 3: the WA yearly cycle (Q1) (PLAN §9.3)

The cycle was refitted with one receiver-type term per unit (coefficient 1; exact basis; no M/B, as in round 4).
Two cell sets were used: (a) near head cells; (b) near + far, which is the headline for WA because only it contains the PGC
site.

**WA (b), two3_20: MIXED.**
- **Pooled:**
  - Q1-R is 0.128 ± 0.041 % (3.1σ, peak day 157). The round-4 form gave 0.091 ± 0.025, peak day 221.
  - Date permutation p 0.024. G1 and G2 pass.
  - **G3 fails:** the scrambled control is 0.202 ± 0.093 (2.2σ).
- **Pooled lapse ratio:** 1.40 ± 0.21. Deep predicts 0.69; a receiver-type change predicts 1.0. Call: shallow/receiver-like.
- **Per group:**

  | group | share | Q1-R | lapse ratio | call |
  |---|---|---|---|---|
  | P = PGC site (B009, B010, B011, PGC; far for WA) | 0.49 | 0.216 ± 0.066 (3.3σ, day 127) | 0.80 ± 0.33 (deep 0.69) | undecided |
  | B = other PB boreholes | 0.44 | 0.107 ± 0.062 (1.7σ) | 1.72 ± 0.39 | shallow-like |
  | O = other stations | 0.07 | 0.04 ± 0.06 | 1.23 ± 0.28 | shallow-like |

- **Per station** (`fig_r5_q1.png`): most WA stations sit on or above the ratio-1 line, not on the deep line.

**WA (a), the round-4 Q1 cells: UNDECIDED.**
- Q1-R is 0.254 ± 0.079 (3.2σ, peak day 183); permutation p 0.020; G1–G3 pass.
- But the pooled lapse ratio is 1.54 ± 0.34 (deep 0.81, shallow 1.11). That is not ≥ 1σ closer to either, so it is
  undecided.
- Only 20 near cells are at the PGC site.

**SVI (a): UNDECIDED.**
- With the station term, Q1-R falls to 0.039 ± 0.027 (1.45σ, day 124). The round-4 form gave 0.052 ± 0.015.
- G2 fails.
- The pooled lapse ratio is 0.55 ± 0.11: deep-like, but below even the deep prediction (0.85).
- **The PGC group alone:** Q1-R 0.099 ± 0.042, lapse ratio 0.52 ± 0.12.

**Reading: the WA cycle is not deep.**
- Its lapse ratio is above 1 wherever it is measured with enough cells.
- Its phase moves when station terms are allowed.
- What is left at ≥ 3σ sits at the PGC-site stations, which round 4 already linked to this cycle.

**Limitation.** For groups that include far cells, the two3_20 shallow prediction is ill-defined: a_s ≈ 0 gives ratios of
−62 to 2.1. The calls then rest on the deep and receiver (1.0) predictions.

## Task 4: SVI 2019-11 per instrument (PLAN §9.4)

**Per instrument** (apparent, 2–5 s, coefficient 1):

| instrument | ΔE(2019-11) | scrambled control |
|---|---|---|
| CN.PGC (vault) | +0.362 ± 0.177 % | +0.47 ± 0.37 |
| PB.B009 | +0.220 ± 0.122 % | −0.54 ± 0.36 |
| PB.B010 | +0.080 ± 0.128 % | +0.26 ± 0.30 |
| PB.B011 | +0.158 ± 0.121 % | −0.33 ± 0.33 |

- Common value +0.184 ± 0.066 % (2.8σ), χ² 1.8 (p 0.62), all four positive. **Call: SITE-WIDE ("the site or real").**
- Pooled scrambled control: −0.08 ± 0.15.
- Vault minus borehole mean: +0.21 ± 0.19 (1.1σ). The vault reads somewhat more, but not significantly.

**Metadata** (FDSN EarthScope / NRCan, s10 epochs):
- No channel epoch starts or ends within ±120 d of 2019-11-15.
- PGC has had the same Trillium 120 PH since 2017-08-10.
- The boreholes have had HS-1-LT / Q330 since 2005-11-02.
- All four instruments have the same 12 data days, 11-09 to 11-25 UTC.
- So it is not an instrument artefact.

**Model (E4)** — every non-PGC-site cell plus one instrument:
- L-SEL: +0.017 to +0.096 (each ≤ 1σ). The non-PGC-site cells alone give +0.042 ± 0.108; everything gives
  +0.111 ± 0.089 (1.2σ).
- NOM: +0.257 ± 0.105 (2.4σ).
- Round 4 HnoSc gave +0.22 ± 0.12.

**Sub-episode parts.**
- 11-22…25, when the carrier families fired: +0.02 to +0.12 at every instrument (≤ 1σ).
- ≤ 11-10 at B009: +0.40 ± 0.20 (one series).
- There were no data days 11-11…21.

**Lapse profile** (pooled site; 2–3 / 3–4 / 4–5 / 5–6 s):
- Values: +0.14 / +0.15 / +0.08 / −0.00 % (σ 0.12–0.08).
- The change falls with lapse time. The best shape is timing-like (χ² 0.95, vs deep 1.78, shallow 1.65, uniform 1.80),
  but it is not decisive (Δχ² < 4).
- Per instrument: PGC, B009 and B011 also prefer the falling (timing) shape; B010 is flat. None is decisive.

**Reading: not instrument-specific.**
- All four co-located sensors see it.
- The far stations do not (round 4: far +0.03 ± 0.20).
- The profile falls with lapse time, against the deep shape.
- The most likely cause is a change local to the PGC site, in the shallow structure or the coupling near the site. A deep
  change and an instrument artefact are both disfavoured.

## Task 5: kernel error as a systematic band (PLAN §9.5)

Plausible kernels: SVI {two3_20, two5_40, ell40}; WA {two3_20, two5_40}.

| quantity | SVI band (plausible) | WA band (plausible) | holds? |
|---|---|---|---|
| Q3 composite, E4 | −0.015 … +0.001 (σ 0.023–0.031) | −0.016 … −0.004 (σ 0.009–0.011) | ✔ within 2σ of 0 for every plausible kernel, and for all of V6 + AVG6 |
| 2σ bound / recovery T1 (0.65–0.90 SVI; 0.54–0.75 WA) | 0.073 % → **0.08–0.11 %** | 0.037 % → **0.05–0.07 %** | ✔ |
| Q3 composite, E4 + Sc 2-yr | −0.037 … −0.002 | −0.004 … +0.007 | ✔ (bound 0.11–0.16 / 0.04–0.06 %) |
| Q1-R, set a | 0.039–0.110 (z 1.1–2.3) | 0.214–0.254 (z 2.0–3.2) | ✘ verdicts UNDECIDED / SHALLOW / UNDECIDED (SVI), UNDECIDED / UNDECIDED (WA) |
| Q1-R, set b | 0.028–0.062 (z 1.5–1.9) | 0.091–0.128 (z 1.7–3.1) | ✘ WA SHALLOW (two5_40) / MIXED (two3_20); no DEEP anywhere |
| SVI 2019-11 ΔD (E4) | +0.11 … +0.26 (1.2–2.4σ) | — | ✔ below 3σ for every kernel |

**What holds across the band:**
- An average slow-slip episode changes the deep-weighted velocity by less than about 0.1 % (SVI) and 0.07 % (WA).
- No yearly cycle is deep.
- No single episode reaches the wording rule in the headline.

## What failed, the measured obstacle, and the next adaptation

1. **WA deep-pulse recovery under L-SEL is 0.75 (rule 0.8).**
   - Obstacle: with 27 episodes, making D orthogonal to eight slow + annual functions removes part of a one-episode pulse.
   - Next: orthogonalise D to the yearly pair only (c1, s1); that is what T4 needs, and the slow part is already
     high-passed by ΔD. Re-run T1/T4. Meanwhile, sizes are quoted divided by the measured recovery (task-5 band).
2. **WA Q1 verdict is not stable across the band (MIXED vs SHALLOW).**
   - Obstacles:
     - the far-cell shallow prediction is undefined under two3_20;
     - the PGC-site group is measured only as far stations for WA.
   - Next:
     - fit the PGC-site stations jointly with their own SVI near cells and the WA far cells, with ONE receiver term per
       station shared across both regions. If the cycle is in the station term, it is receiver-side; if a deep annual
       still remains, it is path;
     - use deep-vs-receiver (1.0) calls for far groups.
3. **WA scrambled control for Q1 fails (0.20 ± 0.09, 2.2σ, set b).**
   - Next: list which stations carry the control's yearly cycle (the same cells' per-station in-phase amplitude).
4. **SVI 2019-11: site-wide at PGC, not deep, not an instrument.**
   - Next, cheapest first:
     - add a PGC-site episode term (one parameter) to E4 and check that ΔD(2019-11) at every other station stays ~0;
     - check PGC-site groundwater / precipitation and local disturbance records for 2019-11 (an environmental site
       change).
5. **WA 2023-03** (3.0σ under NOM with every gate passing; 2.6σ L-SEL; 1.6σ with the smooth regional term).
   - Next: the same dive as 2019-11 (per site, near vs far, per instrument, metadata).
6. **SVI G-site** (PGC 65 % of the weight) stays a flag.
7. **The B1 reference-path control** (a column, not a gate) has 79–87 % of majors within 2σ under L-SEL. Its composites
   are clean.

## Files
- **Scripts:**
  - `svi_science/scripts/s11_r5_{common,model,q1,pgc,figures}.py`;
  - queue `s11_r5_queue.sh`;
  - resume block 10 in `scripts/resume_after_reboot.sh`.
- **Code-test outputs (not results):** `{svi,wa_nw}/tables/test/`.
- **Task 1:** `r5_forms.csv` (per form × setup, all rules and gates, per-group prior diagnostics `diag_ridge_*`),
  `r5_forms_rules.csv`, `r5_forms_onesite.csv`, `forms_<form>_<setup>_{dD,loso,ctrl_*}.csv`, `forms_reproduction_check.csv`.
- **Task 2:** `r5_sc.csv`, `r5_sc_rules.csv`, `sc_*`.
- **Task 5:** `r5_band_q3.csv`, `r5_band_T1_lsel.csv`, `band_dD_<form>_<setup>.csv`, `q1_band.csv`.
- **Task 3:** `q1_verdict.json`, `q1_groups_<set>_<setup>.csv`, `q1_stations_*.csv`, `q1_G1_*.csv`.
- **Task 4** (`svi/tables/`): `pgc_instrument_apparent.csv`, `pgc_instrument_model.csv`, `pgc_split.csv`,
  `pgc_metadata_{channels,events}.csv`, `pgc_data_days.csv`, `pgc_kernel_profile.csv`, `pgc_profile.csv`,
  `pgc_verdict.json`.
- **Figures:**
  - `fig_r5_forms.png`: σ per major and the prior shift per form;
  - `fig_r5_dD.png`: ΔD per episode, E4 + smooth Sc, weight-coloured;
  - `fig_r5_q1.png`: per-station yearly cycle 2–3.5 vs 3.5–6 s and the Q1 band;
  - `fig_r5_pgc.png`: 2019-11 per instrument, episode series and lapse profiles;
  - `fig_r5_band.png`: the Q3 composite over kernels.
