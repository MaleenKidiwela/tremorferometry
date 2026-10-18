# S11 lever 4: stations 80–150 km away (both regions)

Written Mon 2026-10-12, ~9:35 pm PT. The rules were fixed in `PLAN.md` (§1–§9 plus Amendment 1) before any far-station
number. The decision log is `../NOTES.md`. "Scrambled control" means station data shifted to another time of the same day.
dv/v is in %. ΔD is the deep-weighted change of one slow-slip episode relative to its neighbours within ±3 years. The B2
section is post hoc; the coordinator asked for it after the run.

## In plain words
- **The far stations see the families.** Most far pairs do not: 104 of 1,896 SVI pairs and 136 of 2,734 WA pairs are usable.
  - Usable pairs sit at 85–110 km: Nanaimo, the Olympic boreholes, the PGC boreholes, SNB and VGZ.
  - Beyond ~115 km the coda is not seen.
- **Far cells added.** They give 1,041 cells in SVI and 1,297 in WA. That is 15 % and 24 % of the inversion weight.
- **The kernels look as forecast.**
  - In SVI, far pairs carry deep weight a_d 0.97 and shallow weight a_s 0.17. The mid-crust weight a_m is −0.39 and the
    below-45 km weight a_b is 0.32. Merlin forecast 0.91, 0.23, −0.34 and 0.20.
  - In WA, far pairs carry a_d 0.59, a_s 0.06 and a_b 0.64. These paths dive below 45 km.
  - The unit sum is exact to 2e-14. a_d is negative for 0 % of SVI pairs and 7 % of WA pairs; the stop rule is 10 %.
- **The lever works as pre-registered, in both regions.** Each comparison is with the same model on near cells only.

  | | SVI | WA |
  |---|---|---|
  | shallow-pulse leak | 0.47 → 0.37 | 0.60 → 0.41 |
  | σ per major episode | 0.102 → 0.076 % | 0.132 → 0.090 % |
  | majors with σ ≤ 0.1 % | 7 → 14 of 15 | 2 → 20 of 32 |

  WA now passes the precision target (5 majors with σ ≤ 0.1 %) for the first time.
- **No episode-averaged deep change, at tighter precision.** The average of all majors is:
  - SVI −0.000 ± 0.034 %;
  - WA −0.001 ± 0.017 %.

  So an average slow-slip episode changes the deep-weighted velocity by less than 0.07 % (SVI) and 0.035 % (WA), at 2σ.
- **The new model separates deep from shallow per episode.** It adds one regional shallow term per episode
  (Amendment 1, "T+Sc"). Far pairs make it workable:

  | | SVI | WA |
  |---|---|---|
  | σ per major | 0.179 → 0.125 % | 0.249 → 0.127 % |
  | leak of a one-episode shallow change into ΔD | 0.025 | 0.021 |
  | the same leak, without the regional term | 0.37 | 0.41 |

  - An injected deep pulse comes back at 1.13 (SVI) and 1.06 (WA).
  - Price: σ about ×1.5. Upside: the leak problem is gone.
  - No major reaches 3σ in SVI under this model.
- **With the round-3 graph estimator (B2) for the near cells plus the far cells, the shallow leak finally passes its
  rule.**

  | | SVI | WA |
  |---|---|---|
  | design-level leak (rule ≤ 0.35) | 0.28 | 0.22 |
  | σ per major | 0.099 % | 0.115 % |
  | majors with σ ≤ 0.1 % | 8 of 15 | 12 of 32 |
  | composite | +0.036 ± 0.049 % | −0.032 ± 0.031 % |

  With the regional shallow term as well, the leak is 0.018 / 0.017, σ 0.131 / 0.137 %.
- **SVI 2019-11 (round-3 candidate +0.53 ± 0.15 %): the far stations see it, but at about half the size.**
  - Far cells alone: +0.28 ± 0.14 %. These are 34 cells at NLLB, B001, B005–7, LRIV, STW and UWFH: no PGC, other paths,
    same sign.
  - Graph near cells plus far cells: +0.32 ± 0.11 % (2.9σ).
  - With the regional shallow term free: +0.18 ± 0.15 %.
  - Leaving out any site, PGC included, keeps 0.26–0.38 %.
  - It remains a candidate, not a finding: below 3σ once far pairs enter, and the prior test fails.
- **WA single episodes.** Two majors reach 3σ with near and far cells agreeing in sign:
  - 2014-08: +0.49 ± 0.14 % (near +0.55, far alone +0.35);
  - 2023-03: −0.25 ± 0.08 % (near −0.15, far alone −0.35).

  But WA majors scatter 1.8× more than their σ (χ² 3.1) while the scrambled control does not. Scaled by that, neither
  passes 3σ.
- **The annual cycle: the two regions now disagree. This is new.**
  - **SVI.** Far cells carry almost none of the near cycle: far/near ratio 0.35 ± 0.26. Shallow change predicts 0.52;
    deep change predicts 1.73 and is rejected at 5.4σ. This is the first kernel-free test that says SVI's annual cycle is
    not deep.
  - **WA.** Far cells carry the full cycle: ratio 0.95 ± 0.14. Deep predicts 0.94; shallow predicts 0.27 and is rejected
    at 4.7σ.
  - But the WA far cells are mostly the SVI stations PGC/B009–B011. There the annual amplitude is about the same whether
    the families are near (SVI, 0.05–0.09 %) or far (WA, 0.05–0.16 %). The cycle follows the station, not the path.
  - Also, WA's scrambled control fails the Q1 gate (0.225 ± 0.110 %, 2.05σ).
  - So the WA "deep-like" ratio is most likely a receiver-side effect the kernel under-predicts at 100 km, not deep change.
- **Q1 and Q2 with far pairs (headline F1).**

  | | SVI | WA |
  |---|---|---|
  | Q1 deep annual | 0.099 ± 0.037 %, day 261 | 0.108 ± 0.043 %, day 227 |
  | Q1 with kernel-averaged coefficients | 0.047 ± 0.054 % | 0.048 ± 0.057 % |
  | Q2 decade contrast | +0.053 ± 0.057 % | −0.011 ± 0.058 % |

  - Neither Q1 is called deep: Q1/σ is below 3, and the kernel-error B leak is 0.31–0.43.
  - Q2 bounds are |deep decade change| < 0.17 % (SVI) and < 0.13 % (WA) at 2σ. σ(Q2) improves ×0.6 versus N1.

## Numbers
Configurations:
- **N0:** near cells, existing model. It reproduces the stored results to 1e-17.
- **N1:** near cells plus free regional mid-crust (M) and below-45 km (B) smooth terms.
- **F1 (headline):** near plus far cells, M and B free.
- **F0:** near plus far cells, no M or B.
- **FO:** far cells only.

All values use the NOM (ell40) kernel unless stated.

| quantity | SVI N0 | SVI N1 | **SVI F1** | WA N0 | WA N1 | **WA F1** |
|---|---|---|---|---|---|---|
| cells (far) | 3,531 | 3,531 | 4,572 (1,041) | 2,466 | 2,466 | 3,763 (1,297) |
| Q1 % (day) | 0.114 ± 0.037 (229) | 0.158 ± 0.055 (250) | **0.099 ± 0.037 (261)** | 0.164 ± 0.071 (195) | 0.138 ± 0.072 (199) | **0.108 ± 0.043 (227)** |
| Q2 % | +0.059 ± 0.074 | +0.083 ± 0.097 | **+0.053 ± 0.057** | −0.003 ± 0.089 | +0.037 ± 0.092 | **−0.011 ± 0.058** |
| Q3 composite % | −0.000 ± 0.041 | +0.009 ± 0.042 | **−0.000 ± 0.034** | +0.012 ± 0.021 | +0.009 ± 0.021 | **−0.001 ± 0.017** |
| σ per major, median % (≤ 0.1) | 0.091 (12/15) | 0.102 (7/15) | **0.076 (14/15)** | 0.131 (3/32) | 0.132 (2/32) | **0.090 (20/32)** |
| M / B annual % | — | 0.14 ± 0.09 / 0.13 ± 0.09 | 0.10 ± 0.07 / 0.04 ± 0.05 | — | 0.30 ± 0.17 / 0.13 ± 0.06 | 0.08 ± 0.10 / 0.07 ± 0.03 |

Waveform-level injections. Each was injected into the day stacks, run through the full raw chain, and the natural values
subtracted. Kernel-error generation used two5_40 and ell80; inversion used ell40.

| injection (rule) | SVI N0 | SVI N1 | **SVI F1** | WA N0 | WA N1 | **WA F1** |
|---|---|---|---|---|---|---|
| deep pulse, median recovery (0.8–1.1) | 1.14 | 1.14 | **1.13** | 1.06 | 1.06 | **1.07** |
| common shallow pulse, median leak (≤ 0.35) | 0.45 | 0.47 | **0.37** | 0.60 | 0.60 | **0.41** |
| shallow pulse in all majors → composite | 0.57 | 0.57 | 0.48 | 0.22 | 0.22 | 0.17 |
| A deep step: Q2 recovery, ell40 / two5_40 / ell80 (0.7–1.3) | 0.97 / 1.19 / 0.90 | 0.95 / 0.88 / 0.76 | **0.90 / 1.00 / 0.74** | 0.83 / 0.66 / 0.71 | 0.84 / 0.60 / 0.70 | **0.90 / 0.69 / 0.70** |
| B shallow annual → Q1 leak, ell40 / two5_40 / ell80 (< 0.25) | 0.04 / 0.12 / 0.16 | 0.03 / 0.45 / 0.35 | **0.03 / 0.38 / 0.32** | 0.10 / 0.41 / 0.27 | 0.11 / 0.47 / 0.30 | **0.07 / 0.43 / 0.31** |
| mid-crust annual → Q1 leak | 0.58 | 0.04 | 0.03 (F0: 0.64) | 0.05 | 0.01 | 0.01 (F0: 0.40) |

Pre-registered verdicts (PLAN §1):
- **"Lever works": yes in both regions.** F1/N1 leak factor 0.78 (SVI) / 0.69 (WA). Composite-σ factor 0.82 / 0.84.
  Per-major σ factor 0.75 / 0.68.
- **"Deep annual ambiguity reduced": no in both regions.** σ(Q1) falls (×0.67 / ×0.60). But the kernel-error B leaks fall
  only ×0.83–1.06; the rule is ≤ 0.8.
- **Freeing M and B costs kernel robustness.** The B leak is 0.12 in N0 and 0.45 in N1 (SVI two5_40). It is needed anyway:
  without M, a mid-crust annual change leaks 0.40–0.64 into Q1 (F0).

Regional-shallow model (T+Sc) and the B2 combination. B2 leaks are design-level: noise-free linear injection through the same
model. In B0 the design leak is 0.03–0.07 below the waveform leak, because of the reference over-read, which B2 does not have.

| model | SVI σ/major (≤ 0.1) | SVI leak | SVI composite | WA σ/major (≤ 0.1) | WA leak | WA composite |
|---|---|---|---|---|---|---|
| N1 T+Sc (waveform) | 0.179 (0) | 0.067 | +0.074 ± 0.093 | 0.249 (0) | 0.095 | +0.048 ± 0.043 |
| **F1 T+Sc (waveform)** | **0.125 (0)** | **0.025** | +0.003 ± 0.067 | **0.127 (3)** | **0.021** | +0.020 ± 0.031 |
| N1 B2 (design) | 0.103 (7) | 0.36 | +0.020 ± 0.037 | 0.161 (0) | 0.47 | −0.011 ± 0.018 |
| **F1 B2 near + B0 far (design)** | **0.099 (8)** | **0.28** | +0.036 ± 0.049 | **0.115 (12)** | **0.22** | −0.032 ± 0.031 |
| F1 B2 + Sc (design) | 0.131 (0) | 0.018 | −0.024 ± 0.072 | 0.137 (4) | 0.017 | −0.017 ± 0.036 |
| FO far only (design) | 0.135 (0) | 0.11 | +0.025 ± 0.068 | 0.150 (2) | 0.06 | −0.032 ± 0.040 |

SVI 2019-11 across configurations (dv/v %):

| N0 | F1 | F1 T+Sc | far only | N1 B2 | F1 B2 | F1 B2 + Sc | F1 B2 without PGC |
|---|---|---|---|---|---|---|---|
| +0.19 ± 0.12 | +0.21 ± 0.08 | +0.21 ± 0.14 | +0.28 ± 0.14 | +0.54 ± 0.15 | +0.32 ± 0.11 | +0.18 ± 0.15 | +0.30 ± 0.13 |

Gates on F1:

| gate (rule) | SVI | WA |
|---|---|---|
| scrambled control Q1 (< 2σ) | **pass** 0.06 ± 0.12 % (359 far control cells; far-only χ² 1.01) | **fail** 0.225 ± 0.110 % (2.05σ) |
| scrambled control Q3 (≥ 90 %, composite < 2σ) | **pass** 15/15, +0.04 ± 0.12 | **pass** 0.91, −0.01 ± 0.04 |
| leave-one-site-out, every site incl. each far site | **pass** (Q1 0.08–0.12; Q3 100 %) | **pass** (Q1 0.09–0.13; Q3 ≥ 94 %) |
| family split Q1 / Q3 | **pass / pass** (Q1 had failed in near-only) | **pass / pass** |
| far/near annual amplitude (deep iff ≥ 1σ closer) | ratio 0.35 ± 0.26; deep 1.73 (z −5.4), shallow 0.52 (z −0.7) → **shallow** | ratio 0.95 ± 0.14; deep 0.94 (z 0.1), shallow 0.27 (z 4.7) → **deep-like** (see caveats) |
| lapse ratio 3.5–6 / 2–3.5 s, far cells | uninformative (2–3.5 s amplitude ≈ 0: −7.9 ± 21.9) | 1.26 ± 0.30 vs deep 0.58; the shallow prediction is undefined (a_s changes sign between windows) → uninformative |
| lapse ratio, near cells | 0.55 ± 0.11 (deep 0.81, shallow 1.15) | 1.24 ± 0.16 (= shallow 1.23) |
| site dominance (≤ 50 %) | **fail** PGC 57 % (was 67 %) | **pass** B013+B943 34 % |
| prior ×2/×0.5 (< 0.1σ) | composite pass (≤ 0.01σ); majors **fail** 0.27/0.48σ; Q1 **fail** 0.15σ | composite 0.09σ pass; majors **fail** 0.34/0.55σ; Q1 **fail** 0.26/0.32σ |

Same-station annual amplitudes are in `station_near_far_annual.csv` (`figures/fig_dist_station.png`). At the PGC-site
boreholes the cycle as far stations for WA families (B009 0.16 ± 0.04, B011 0.15 ± 0.05, B010 0.09 ± 0.05 %) is at least
as large as for their own near SVI families (0.05–0.08 %). The kernel gives those far pairs a_s ≈ 0.01–0.04.

## What failed, the measured obstacle, and the next adaptation
1. **Prior sensitivity: every model fails the prior rule** (majors move 0.3–0.6σ).
   - Obstacle: ΔD is a contrast, but the 1 % prior also fixes the common level of D(e), M, B and Sc, which trade off with
     the smooth shallow terms.
   - Next: replace the prior on D(e) and Sc(e) by an exact sum-to-zero constraint within each ±3-yr window, so that
     only contrasts are estimated. Re-run the prior test.
2. **WA per-major excess scatter** (χ² 3.1 in F1, control normal).
   - Next: an empirical σ per major from a delete-one-site / delete-one-family jackknife, before any WA single episode is
     named. Check 2014-08 and 2023-03 with it.
3. **The annual cycle follows the station, not the path, at the PGC boreholes.**
   - Obstacle: either the kernel under-predicts receiver-side sensitivity at 100 km, or something at the station is
     annual.
   - Next: a joint fit per station of near and far cells with a free receiver term not scaled by a_s, versus the
     kernel-scaled shallow term. Also inspect which stations carry WA's failing scrambled-control annual (0.225 %).
4. **Kernel-error B leak stays 0.3–0.4 once M and B are free.**
   - Next: the round-3 kernel average (AVG6, leave-one-kernel-out) applied to far pairs. Round-3 tables show this passes
     for near cells.
5. **Far pairs beyond ~115 km are unusable** (cc20 ≈ 0).
   - Next: decade-half or year stacks for Q2 only; or stacking co-located families (≤ 3 km) per far station, after testing
     coda coherence between them.
6. **B2 is not yet applied to the far cells.**
   - Next: run the round-3 graph measurement on the far stores (grid −28…+10 s), about 15 min. Then re-run the waveform-level
     pulses for the combined B2 model. The leak above is design-level.
7. **2019-11 needs the PGC timing glitch ruled out.** Far stations without PGC see +0.28 ± 0.14 %.
   - Next: a non-PGC near + far composite with the sum-to-zero constraint.

## Files
- Tables (per region): `<region>/tables/`
  - `far_pairs.csv`, `far_stations_summary.csv`, `inventory_summary.json`;
  - `series.csv`, `episode_raw.parquet`, `kernel_coef.parquet`, `kernel_summary.csv`;
  - `inv_models.csv`, `inv_q3.csv`, `inv_responses.csv`, `inv_summary.json`, `inv_prior.csv`;
  - `gates.csv`;
  - `b2_models.csv`, `b2_q3.csv`, `b2_loso.csv`, `b2_2019_11.csv` (SVI), `far_cells_2019_11.csv` (SVI).
- Figures: `figures/`
  - `fig_dist_inventory.png`: pick-SNR vs distance, coloured by causality;
  - `fig_dist_annual.png`: near vs far cells by day of year, coloured by CC, weighted monthly lines;
  - `fig_dist_q3.png`: ΔD per major for N1, F1, F1 T+Sc, coloured by weight;
  - `fig_dist_b2.png`: graph-estimator combinations and far-only;
  - `fig_dist_lever.png`: leak, σ and B leak for N0 → N1 → F1;
  - `fig_dist_station.png`: same-station near vs far amplitude.

  No error bars on dv/v.
- Scripts: `svi_science/scripts/s11_dist_{common,inventory,windows,measure,kernel,model,invert,gates,b2,figures}.py`, queue
  `s11_dist_queue.sh`.
