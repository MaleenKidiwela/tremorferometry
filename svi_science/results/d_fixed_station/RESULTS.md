# D. Did NW LFE activity change across eras once the station network is held fixed? — results

Run 2026-10-07, ~3:15–3:20 pm PDT. Rules frozen in `PLAN.md` (3:15 pm PDT) before any number was computed. Script
`svi_science/scripts/d_fixed_station.py` (`once`, `stats`), figures `svi_science/scripts/d_figures.py`.
"Time-reversed control" = the same detector with every template reversed in time (main false-alarm reference).
"Scrambled control" = each station's data shifted to a different time of the same day. Count = P_own_WAVE ≥ 0.75; "net" =
data minus time-reversed control. Fixed network: core9 = the 9 PB boreholes (B001, B003–B007, B013, B017, B018), primary;
all10 = core9 + GNW, secondary. "Strict day" = every fixed station of the family's template set is live that day.

## Plain answer

**No change in NW LFE activity between 2007–11, 2012–14 and 2015–26 can be detected once the network is held fixed.
Most of the production drop 0.48 / 0.35 / 0.28 was the station network.**
- With the 9 fixed boreholes on strict days, net repeats per family-tremor-day are **0.27 / 0.26 / 0.30** (production
  WAVE on the same days: 0.47 / 0.31 / 0.25). The E2→E3 drop disappears; E3 is, if anything, slightly higher.
- Per PNSN tremor epicentre (from 2009-08): fixed 0.0112 / 0.0057 / 0.0082, production 0.0126 / 0.0067 / 0.0066.
- No era pair passes the pre-registered rule (core9 per-epicentre permutation p = 0.041 / 0.070 / 0.040; rule p < 0.01).
- So only bounds can be given. At 99 %, with fixed boreholes on strict days and per tremor epicentre:
  - E2 vs E3: between −67 % and +5 %;
  - E1 (2009-08–2011 only) vs E3: between −53 % and +144 %.
- Per ETS episode, 2008–2026: no trend on the fixed network (+0.9 %/yr, 95 % CI −5.4 to +5.9 %/yr; 60 major episodes).
  A steady change of more than ~6 %/yr would have been seen.
- One provisional item, not pre-registered: inside 2015–26 only, the fixed-network rate per epicentre falls about 9.5 %/yr
  (p 0.007). The production 19-station network, which is also unchanged through 2015–26, does not confirm it
  (−5.2 %/yr, p 0.15). Not a result yet (see §5).

![era rates](../../figures/d_era_rates.png)
![episode series](../../figures/d_episode_series.png)

## 1. Net repeats per family-tremor-day (measure 1)

All checks reproduce the post_v2 summaries exactly (`checks.json`, 15/15: production WAVE view, core9 and all10, all days
and strict days). Clusters (runs of tremor dates with gaps ≤ 5 days, the resampling unit): E1 43, E2 55, E3 220.

| network | days | E1 2007–11 | E2 2012–14 | E3 2015–26 | E1/E3 (95 %) p | E2/E3 (95 %) p |
|---|---|---|---|---|---|---|
| production (WAVE ≥ 0.75) | all | 0.471 | 0.306 | 0.246 | 1.91 (1.08–2.94) 0.056 | 1.24 (0.64–1.80) 0.38 |
| core9 | all | 0.374 | 0.237 | 0.267 | 1.40 (0.88–2.11) 0.20 | 0.89 (0.48–1.35) 0.60 |
| core9 | **strict** | **0.271** | **0.260** | **0.300** | 0.90 (0.41–1.67) 0.73 | 0.87 (0.47–1.27) 0.53 |
| all10 | all | 0.465 | 0.280 | 0.295 | 1.58 (0.99–2.34) 0.10 | 0.95 (0.59–1.38) 0.81 |
| all10 | strict | 0.341 | 0.290 | 0.330 | 1.04 (0.43–2.00) 0.91 | 0.88 (0.55–1.25) 0.56 |

- **Scrambled control instead of time-reversed** (same table, `era_rates.csv` variant `scr075`): core9 strict
  0.259 / 0.254 / 0.283. The fixed-network controls are larger than production's (the WAVE model was refit on 9–10
  stations): core9 time-reversed 2.0 / 3.8 / 2.5 % of data, scrambled 5.6 / 6.0 / 8.0 %. Production: ≤ 0.6 % and ≤ 1.6 %.
- **Counted once per source** (families_v1_3 `nw_source_group`; member families' candidates within ±2 s merged):
  - Every rate drops by about 22 % (core9 strict 0.215 / 0.201 / 0.229).
  - The ratios do not change (E1/E3 0.94, E2/E3 0.88).
  - Before merging, the raw counts equal rates_daily exactly in all three networks.
- **E1 caveat (definitional).** Before 2009-08-06, "tremor day" means every family-day inside a WECC ETS window. That dilutes
  E1: strict 2008-04 gives 0.12 and 2009-04 gives 0.16. Restricted to days from 2009-08-06, E1 is higher in every network:
  - production 0.556, core9 all 0.480, core9 strict 0.404, all10 strict 0.511;
  - E1/E3 for core9 strict = 1.35, p 0.30.
  - This rests on two ETS: 2010-03 (0.97) and 2010-08 (0.70).
  - The 2011-07 ETS has no strict day at all: B017 and B018 were down 2011-07-20..09-10.
  - E1 strict coverage is 11 % (2007) and 33 % (2008) of days. So E1 can only be stated as a bound.

### How much of the production drop was the network (measure 5)
The share of the production ln-ratio that is removed when the network is held fixed (measure 1):

| comparison | production | core9 all | core9 strict | all10 all | all10 strict |
|---|---|---|---|---|---|
| E1/E3, all tremor days | 1.91 | 1.40 (48 % removed) | 0.90 (116 %) | 1.58 (30 %) | 1.04 (95 %) |
| E1/E3, days from 2009-08 | 2.26 | 1.79 (28 %) | 1.35 (63 %) | 1.95 (18 %) | 1.55 (46 %) |
| E2/E3 | 1.24 | 0.89 (155 %) | 0.87 (166 %) | 0.95 (124 %) | 0.88 (160 %) |

- The E2→E3 step is entirely the network. On the fixed network it reverses slightly (n.s.).
- The E1→E3 step is about half to all network, depending on which E1 days are used.
- What remains of it is the 2010 ETS pair, and it is not significant.

### n_sta_day covariate
- **Day-level weighted fit** (rate ~ era + live fraction, at the pooled mean live fraction 0.955): core9 E1/E3 1.47
  (95 % 1.01–2.22) and E2/E3 0.87 (0.49–1.30).
  - The live-fraction slope is +0.021 per +0.1 (−0.066 to +0.100), i.e. not resolved (`wls_covariate.csv`).
- **Strata** (`strata_live.csv`):
  - In E2 and E3, days with a borehole down have much lower rates than strict days (core9 E3: 0.17 vs 0.30).
  - In E1 it is the opposite (0.47 vs 0.27). The E1 non-strict days contain the 2011-07/08 ETS and the 2010-08 build
    days, so this split mixes station loss with which ETS fell on which days.
  - This is why strict days are the clean comparison.

## 2. Net repeats per PNSN tremor epicentre ≤ 30 km (measure 2, from 2009-08-06; detection-independent normaliser)

| network | days | E1 (2009-08–11) | E2 | E3 | E1/E3 (99 %) p | E2/E3 (99 %) p | E1/E2 p |
|---|---|---|---|---|---|---|---|
| production | all | 0.0126 | 0.0067 | 0.0066 | 1.90 (0.67–2.69) 0.011 | 1.00 (0.61–1.32) 0.99 | 0.023 |
| core9 | all | 0.0109 | 0.0052 | 0.0072 | 1.51 (0.70–2.17) **0.041** | 0.72 (0.37–1.09) **0.070** | **0.040** |
| core9 | strict | 0.0112 | 0.0057 | 0.0082 | 1.36 (0.47–2.44) 0.15 | 0.69 (0.33–1.05) 0.047 | 0.067 |
| all10 | all | 0.0131 | 0.0061 | 0.0080 | 1.64 (0.76–2.35) 0.026 | 0.77 (0.44–1.13) 0.11 | 0.027 |
| all10 | strict | 0.0138 | 0.0064 | 0.0091 | 1.52 (0.37–2.60) 0.074 | 0.71 (0.42–1.03) 0.049 | 0.039 |

- PNSN epicentres per family-tremor-day are 58.7 / 54.1 / 46.3 by era (core9, major episodes). Within the E3 major
  episodes they show no trend (ρ +0.12, p 0.43).
- Counted once per source: same picture (core9 strict 0.0087 / 0.0044 / 0.0063).

**Pass rule (frozen):** a "real change" needs (i) core9 per-epicentre p < 0.01, (ii) the same direction in all10, and
(iii) survival on strict days.
- E1–E3: (i) fails (p 0.041).
- E2–E3: (i) fails (p 0.070).
- E1–E2: (i) fails (p 0.040).
- **Verdict: no change detectable for every era pair.**
- Bounds (99 %, core9 strict per epicentre): E1/E3 0.47–2.44, E2/E3 0.33–1.05, E1/E2 0.65–5.0.
- Secondary verdict on measure 1 (strict, all tremor days): also no change. Bounds E1/E3 0.34–1.95, E2/E3 0.37–1.42.
- The direction is the same in core9 and all10 everywhere: E1 high, E2 lowest, E3 in between. Only the size of the effect
  falls short.

## 3. Deep / shallow ratio per era (measure 3)

Deep = down-dip position above the NW median (`deep_shallow.csv`). Ratio = deep rate / shallow rate.

| network | days | E1 | E2 | E3 | permutation p (E1–E3 / E2–E3 / E1–E2) |
|---|---|---|---|---|---|
| production | all | 1.27 | 2.11 | 2.12 | 0.21 / 0.99 / 0.12 |
| core9 | all | 1.46 | 2.82 | 2.24 | 0.24 / 0.56 / 0.072 |
| core9 | strict | 1.49 | 2.87 | 2.08 | 0.48 / 0.40 / 0.25 |
| all10 | strict | 1.27 | 3.01 | 2.00 | 0.37 / 0.30 / 0.18 |

- Deep families are about twice as active per tremor day as shallow ones in every era.
- No era difference reaches p < 0.05.
- With the Hyp2000-depth split the contrast is weaker (core9 strict 1.05 / 1.73 / 1.33) and again shows no era difference.
- The E2 upper CI is wide (to ~15), because there are few shallow events in E2.

## 4. Per-ETS-episode series 2007–2026 (measure 4; `episode_series.csv`, `episode_trend.csv`)

MAJOR episodes with ≥ 100 family-tremor-days. Slope is of ln rate vs year, weighted by family-tremor-days, with an
episode bootstrap 95 % CI. Episodes with net ≤ 0 (2–4 per series) are left out of the log fit and counted in the table.

| network | days | measure | episodes | Spearman ρ (perm p) | slope %/yr (95 %) |
|---|---|---|---|---|---|
| **core9** | **strict** | per family-tremor-day | 60 (2008–2026) | −0.04 (0.77) | **+0.9 (−5.4 to +5.9)** |
| core9 | strict | per epicentre | 58 (2009–2026) | −0.09 (0.53) | −1.3 (−6.2 to +4.2) |
| all10 | strict | per family-tremor-day | 60 | −0.03 (0.84) | +0.1 (−6.0 to +4.9) |
| core9 | all | per family-tremor-day | 66 (2007–2026) | −0.18 (0.16) | −2.5 (−7.0 to +1.5) |
| production | all | per family-tremor-day | 65 | −0.19 (0.14) | −3.9 (−8.2 to +0.3) |
| production | all | per epicentre | 62 | −0.19 (0.15) | −3.4 (−6.8 to +0.8) |

- On strict days, the fixed network shows no trend over 18 years.
- Production and fixed "all days" lean negative. That lean comes from the network-change step and the E1 2010 pair.

## 5. Provisional (not pre-registered): trend inside 2015–26

These rows were added in the script but are not in PLAN.md. Treat them as post hoc.

| network | days | measure | ρ (p) | slope %/yr (95 %) |
|---|---|---|---|---|
| core9 | strict | per epicentre | −0.42 (0.007) | −9.5 (−14.7 to −4.6) |
| core9 | strict | per family-tremor-day | −0.22 (0.16) | −6.8 (−17.5 to +1.8) |
| all10 | strict | per epicentre | −0.36 (0.024) | −8.0 (−14.1 to −1.9) |
| production (E3 = a fixed 19-station set) | all | per epicentre | −0.22 (0.15) | −5.2 (−10.6 to +0.5) |
| production | all | per family-tremor-day | −0.04 (0.77) | −2.0 (−12.3 to +8.9) |

- The decline shows up only per epicentre and only on the borehole network.
- The PNSN epicentre density does not trend.
- The 19-station production network, which is also unchanged through E3, does not confirm it.
- Borehole sensor ageing is one possible cause (not tested).
- Status: a lead to check, not a finding.

## Files
- `PLAN.md`: pre-registration.
- `era_rates.csv`: network × day set × variant × measure × era, with rate, 95 % CI, numerator, denominator and clusters.
  Variants: rev075, scr075, once_rev075, rev050, rev085.
- `era_ratios.csv`: ratio, 95 % and 99 % CI, permutation p.
- `deep_shallow.csv`, `strata_live.csv`, `wls_covariate.csv`, `episode_series.csv`, `episode_trend.csv`, `checks.json`.
- `once_{prod,core9,all10}.parquet`: per family-day once-per-source counts.
- `d.log`.
- Figures: `svi_science/figures/d_era_rates.png`, `d_episode_series.png`.
