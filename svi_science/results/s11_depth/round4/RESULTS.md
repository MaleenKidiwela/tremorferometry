# S11 round 4: combined model, the SVI 2019-11 candidate, does the lapse fit pick the right kernel?

Written Mon 2026-10-12, ~10:30 pm PT. The rules are in `../PLAN.md` §8 (written ~9:30 pm PT before any round-4 number),
plus Amendment 1 (~9:45 pm PT, after a natural-value smoke test, before any injection) and Amendment 2 (~9:55 pm PT,
**post-test**: written after a code test of the dive whose numbers were seen). The queue ran 9:57–10:12 pm PT; every stage
had rc = 0. Decision log: `../NOTES.md` (decisions 24–26 and the result lines). Tables: `{svi,wa_nw}/tables/`. Figures:
`figures/fig_r4_{dD,suite,lapse,2019_11}.png`.

Conventions:
- dv/v is in %. ΔD = the deep-weighted change of one slow-slip episode minus its neighbours within ±3 years.
- "Scrambled control" = station data shifted to another time of the same day.
- NOM = ell40 kernel. L-SEL = two3_20, the kernel the lapse fit chooses.
- H = the pre-registered combined model. HnoSc = H without the regional shallow term (Amendment 1).
- Everything here is provisional.

## In plain words

**1. The combined model works as a leak fix. Its price is precision, and the price comes from one term.**
- H combines:
  - near and far graph cells (B2, ≥ 50-event nodes, cc-weighted);
  - one regional shallow term per episode (Sc);
  - smooth shallow terms per site;
  - a receiver term per station;
  - free mid-crust and > 45 km terms;
  - an exact sum-to-zero on the episode terms instead of the prior.
- **Shallow leak.** H removes the shallow-to-deep leak completely (SVI 0.006, WA 0.006; rule ≤ 0.35). A deep pulse comes back
  at 0.98 (SVI) and 0.91 (WA).
- **Precision.** Once the 1 % prior is gone, the deep term and the regional shallow term trade off. σ per major ETS is
  0.20 % in both regions, and the composite σ is 0.18 % (SVI) and 0.12 % (WA).
- **What this tells us about round 3.** The good precision of the earlier "T+Sc" model came from its 1 % prior on Sc. That
  prior is also why its prior test failed.
- **Without Sc (HnoSc) the precision comes back and the leak stays under the rule:**

  | | SVI | WA-NW |
  |---|---|---|
  | σ per major | 0.111 % (4/15 ≤ 0.1) | 0.128 % (5/24 ≤ 0.1, acceptance met) |
  | shallow leak | 0.29 | 0.26 |
  | composite | −0.007 ± 0.041 % | −0.011 ± 0.017 % |

  So an average slow-slip episode changes the deep-weighted velocity by less than 0.08 % (SVI) and 0.03 % (WA) at 2σ.
- **With the data-chosen kernel (two3_20) the trade-off shrinks.** Even with Sc free, σ per major is 0.099 % (SVI, 9/15 ≤ 0.1)
  and 0.132 % (WA). The leak is 0.008 / 0.023.

**2. The sum-to-zero does NOT fix the prior test.**
- Majors still move 0.17–0.57σ when the remaining prior (on the smooth S, R, M, B terms) is doubled or halved. The rule is
  < 0.1σ.
- What has improved: in H the composite moves only 0.002–0.06σ at SVI.
- The remaining dependence is on the conditioning of the smooth terms, mainly the near-duplicate pair "S per site
  (× a_s)" and "R per station (× 1)".

**3. The per-station receiver term absorbs a yearly cycle that sits at one site. It cannot absorb one shared by every
station.**
- One site at a time, the receiver-type annual cycle leaks:
  - ≤ 0.09 into the composite (SVI) and ≤ 0.02 (WA);
  - ≤ 0.26 into a major (B926 is the worst).
- The same cycle at all stations leaks 0.38 (SVI) and 0.11 (WA) into the composite, and 0.69–0.78 into single majors (rule
  ≤ 0.35). Without R these are 0.43 / 0.11, so R barely helps here.
- Reason: a cycle common to every station is the same shape in the data as an annual pattern in the episode terms.
- This was a design (injection) test. Q1 was not refitted with the receiver term.

**4. SVI 2019-11: undecided by the pre-set rule. It is not a clock or timing artefact. The evidence points to something local
to the PGC site.**

| model form | ΔD(2019-11) |
|---|---|
| H | −0.24 ± 0.26 % (Sc absorbs it) |
| HnoSc | +0.22 ± 0.12 % (1.8σ) |
| round-3 model on the near graph cells | +0.57 ± 0.15 % (round-3 value reproduced) |

- **Timing artefact ruled out.**
  - The per-event realignment already removes the PGC −6 ms clock step before stretching. 2019-11 node stacks sit
    0.0 ± 0.2 ms from their series median.
  - Realigning each node again leaves ΔD unchanged: +0.221 → +0.217.
  - Without the PGC vault: +0.21 ± 0.12.
  - All three timing rules are negative.
  - The stretch bias would be −0.029 % per ms of delay, so a residual shift of more than a few ms would have mattered.
- **Who carries it.** The PGC site (B009, B011, PGC) carries 77 % of HnoSc's value.
  - Without the PGC site: +0.06 ± 0.15.
  - Near stations: +0.70 ± 0.30. Far stations: +0.03 ± 0.20. A deep change of this size would show up at the far stations
    too.
- **Other blocking tests.**
  - Leave-one-family-out: minimum 1.4σ (family 56).
  - Lapse ratio: 0.68 ± 0.45 (deep 0.84, shallow 1.17) → undecided.
  - Scrambled controls are clean (+0.02 ± 0.25 and −0.24 ± 0.36). Look-elsewhere is normal (2 of 53 episodes at ≥ 1.8σ;
    3.4 expected).
- **Context.**
  - No metadata change at the carrying stations within ±120 d.
  - Local earthquakes: only M ≤ 2.6.
  - Tremor near the carrier families was weaker than in a typical major (2,082 vs a median of 3,061 PNSN windows).
  - The carrier families (10, 56, 87, 75, 43, 54, 17) were silent during the early burst (2019-11-10 5:30 am PST, 413
    events in 5 h). They fired only on 11-22…11-25.
  - The ballistic S−P term for 2019-11 (+0.060 ± 0.027 %) equals what a real +0.22 % deep change would predict (+0.060 %).
    Consistent, but weak support.

**5. The lapse fit does pick the right kernel. two3_20 is adopted for Q1 in both regions.**
- **Synthetic test.** Data were generated with each of the six kernels.
  - When the truth is two3_20, the fit picks it 98 % (Gaussian noise) / 86–98 % (residual bootstrap) of the time.
  - When the truth is another kernel, the fit picks two3_20 at most 1–11 % of the time.
  - Confidence that two3_20 is the true kernel: 0.72–0.96 (SVI), 0.86–0.92 (WA).
- **Q1 under two3_20:**

  | | SVI | WA-NW |
  |---|---|---|
  | Q1 deep annual | 0.052 ± 0.015 % (3.4σ), day 219 | 0.151 ± 0.043 % (3.5σ), day 208 |
  | fails | G2 (halves peak 70 d apart), G3 (scrambled control 0.24 ± 0.11) | G4, G5 |
  | passes | — | G1, G2, G3; date permutation p 0.010; in-sample B leak 0.03 |
  | verdict | not called deep | meets the §3/§8.3 wording rule — a first |

  - **But in WA the two supporting tests point to shallow.** The lapse ratio is 1.24 ± 0.16, equal to the shallow prediction
    1.11 (deep 0.81). The station pattern fits neither form. The lever-4 far/near test found that this cycle follows the
    PGC-site stations, and prior ×0.5 moves Q1 by −0.21σ.
  - So WA Q1 is reported as "deep by the pre-set rule, disputed by its own supporting tests". It is not claimed.

## Checks against the PLAN §8 rules

### Task 1: combined model (NOM unless stated; `r4_verdict.json`, `r4_gates.csv`, `r4_rules.csv`, `r4_onesite_receiver.csv`)

| rule | SVI H | SVI HnoSc | WA H | WA HnoSc |
|---|---|---|---|---|
| σ per major (≤ 0.1) | 0.199 (0/15) | 0.111 (4/15) | 0.196 (1/24) | 0.128 (5/24) ✔ acceptance |
| composite | +0.057 ± 0.177 | −0.007 ± 0.041 | +0.312 ± 0.115 | −0.011 ± 0.017 |
| T1 deep pulse 0.8–1.1 | 0.98 ✔ | 0.96 ✔ | 0.91 ✔ | 0.93 ✔ |
| T2 common shallow leak ≤ 0.35 | 0.006 ✔ (fixed) | 0.29 ✔ | 0.006 ✔ (fixed) | 0.26 ✔ |
| T3 one-site leak ≤ 0.35 | ≤ 0.07 ✔ | ≤ 0.20 ✔ | ≤ 0.04 ✔ | ≤ 0.11 ✔ |
| T4 receiver annual, all stations | 0.38 / 0.69 ✘ | 0.25 / 0.71 ✘ | 0.11 / 0.75 ✘ | 0.11 / 0.74 ✘ |
| T5 shallow annual → composite | 0.04 ✔ | 0.01 ✔ | 0.01 ✔ | 0.01 ✔ |
| T6 kernel error | ✘ (T1 0.80–1.16) | ✘ (T2 0.42 at two5_40) | ✘ (T1 0.72–1.40) | ✘ (T1 0.80–1.28; T2 0.53) |
| G-ctrl scrambled control | ✔ | ✔ | ✔ | ✔ |
| G-loso | ✔ | ✔ | ✘ (B001: 58 %) | ✔ |
| G-fam | ✔ | ✔ | ✘ (odd −0.16, even +0.37) | ✔ |
| G-prior < 0.1σ | ✘ 0.35σ | ✘ 0.43σ | ✘ 0.26σ | ✘ 0.45σ |
| G-site ≤ 50 % | ✘ PGC 65 % | ✘ PGC 65 % | ✔ B013+B943 31 % | ✔ |

- No major passes the wording rule in any form.
- WA H has one major at ≥ 3σ (2019-05), but it fails G-loso, G-fam and G-prior. WA H also shows excess scatter of the majors
  (χ² 3.1), which disappears without Sc (χ² 1.3).
- **L-SEL column (two3_20):**

  | | SVI H | SVI HnoSc | WA H | WA HnoSc |
  |---|---|---|---|---|
  | σ per major | 0.099 | 0.082 | 0.132 | 0.124 |
  | composite | +0.014 ± 0.048 | +0.008 ± 0.033 | +0.061 ± 0.028 | −0.038 ± 0.018 |

  - In-sample: deep pulse 0.89–0.97; leak 0.008–0.08 (with Sc ≤ 0.02).
  - Failing: the prior test, T4, and T6. T6 fails because with a weak-scattering truth the deep pulse comes back at only
    0.50–0.77.
- **What each term buys** (natural values, `r4_models.csv`, SVI):

  | form | σ per major | composite |
  |---|---|---|
  | H | 0.199 | ± 0.177 |
  | without Sc | 0.111 | ± 0.041 |
  | 1 % prior instead of the constraint | 0.149 | ± 0.067 |
  | without R | 0.173 | — |
  | without M/B | 0.173 | — |
  | without R and Sc | 0.076 | ± 0.031 |

  Most of the cost is Sc; R adds about 0.03 % per major.

### Task 2: 2019-11 (`svi/tables/dive_*.csv`, `dive_verdict.json`)

| item (pre-set) | H | HnoSc |
|---|---|---|
| (a) carriers | — | PGC site 77 % (B009 37, B011 24, PGC vault 21 %); near 89 %, far 11 %; families 10, 56, 87, 75, 43, 54, 17 |
| (b) subsets | near +1.16 ± 0.81; far −0.19 ± 0.33; boreholes +1.41 ± 0.81 | near +0.70 ± 0.30; far +0.03 ± 0.20; boreholes +0.23 ± 0.23; without PGC site +0.06 ± 0.15 |
| (c) own controls | B2 −0.16 ± 0.47; B1 +0.02 ± 0.76 ✔ | B2 +0.02 ± 0.25; B1 −0.24 ± 0.36 ✔ |
| (d) leave-one-family-out (rule ≥ 2σ) | min z 0.15 ✘ | min z 1.40 (family 56) ✘ |
| (e) lapse (apparent, near graph) | — | 2–3.5 s +0.17 ± 0.08; 3.5–6 s +0.11 ± 0.05; ratio 0.68 ± 0.45 → undecided |
| (g) timing rules i / ii / iii | — | negative / negative / negative |
| (g) metadata ±120 d | — | none (FDSN failed for CN.PFB only) |
| verdict | undecided | undecided |

- The H column's blocking tests are the same four as for HnoSc (below).
- **Supporting detail:**
  - (b) **Borehole-only.**
    - In H (with Sc): +1.41 ± 0.81 near and far, +2.88 ± 1.33 near only. These large values are an Sc trade-off artefact of
      small subsets.
    - In HnoSc: near boreholes +0.55 ± 0.40.
  - (e) **Lapse.**
    - Each window inverted with its own kernel: 2–5 s +0.53 ± 0.18; 2–3.5 s +0.27 ± 0.17; 3.5–6 s +0.17 ± 0.21.
    - B0 path, apparent: +0.09 / +0.06 ± 0.05.
  - (f) **Arrival times.**
    - Clock term in 2019-11: PGC −6.0 ± 2.5 ms; boreholes ≤ 0.6 ms.
    - S lag vs the series median: PGC −6.4 ms; others −1.7…0 ms.
    - Region S−P: +0.060 ± 0.027 %. Expected if real: +0.22 % × 0.27 (the share of S time spent in 25–36 km) = +0.060 %.
  - (g) **Timing details.**
    - Per-event realignment lag in 2019-11: PGC −3.4 ms, B009 −4.1 ms. These offsets are removed before stacking.
    - Node-stack lag after that: |lag| ≤ 0.2 ms at the PGC site.
    - Node realignment: +0.221 → +0.217 (HnoSc) and −0.235 → −0.220 (H).
    - The shift-free estimator (eps2p) on the 5-ms grid is too coarse: +0.65 ± 0.29 vs +0.15 ± 0.13 on the same cells.
  - (h) **Seismicity.**
    - Earthquakes: 5 events with M ≥ 2, largest M2.6, 41–99 km away.
    - PNSN tremor within 30 km of the carriers: 2,082 windows in the 2019-11 unit vs a median of 3,061 for the other majors.
  - (i) **The burst.**
    - Largest early burst: 2019-11-10 5:30 am PST, 413 events in 5 h, 14 families, along-strike x ≈ 24 km.
    - Carrier-family ΣP was 0 until 11-22.
    - Kernel-free split: UTC ≤ 11-10 +0.22 ± 0.14 % (3 series); ≥ 11-11 +0.07 ± 0.06 %; whole episode +0.08 ± 0.05 %.
- **Blocking tests (both forms):** < 3σ; model gates (prior, site, T4); leave-one-family-out < 2σ; near and far do not agree
  (far ≈ 0).

### Task 3: lapse-selection recovery (`{svi,wa_nw}/tables/lapse_confusion.csv`, `lapse_verdict.json`)

| | SVI | WA-NW |
|---|---|---|
| data Δχ² (two3_20 = 0) | two5_40 0.6, ell40 5.4, ell80 9.7, ell10 14.3, ell20 23.6 | two5_40 8.9, ell80 20.0, ell40 22.5, ell20 36.9, ell10 73.9 |
| P(pick two3_20 \| true two3_20), N1 / N2 | 0.98 / 0.86 | 0.98 / 0.98 |
| max P(pick two3_20 \| other truth), N1 / N2 | 0.01 / 0.11 | 0.05 / 0.08 |
| confidence (uniform prior), N1 / N2 | 0.96 / 0.72 | 0.92 / 0.86 |
| noise-only control | no preference (0.12–0.23) | picks two5_40 91 % (structural bias, but never two3_20) |
| adoption rule | ✔ adopt two3_20 | ✔ adopt two3_20 |
| Q1 (adopted) | 0.052 ± 0.015 % day 219 (3.4σ) | 0.151 ± 0.043 % day 208 (3.5σ) |
| G1 / G2 / G3 | ✔ / ✘ / ✘ (control 0.24 ± 0.11) | ✔ / ✔ / ✔ (control 0.18 ± 0.09) |
| date permutation p | 0.14 | 0.010 |
| G4 / G5 | ✘ / ✔ (0.55 ± 0.11, closer to deep) | ✘ / ✘ (1.24 ± 0.16 = shallow) |
| kernel-error, posterior-weighted | Q2 recovery 0.92, B leak 0.02 | Q2 recovery 0.80, B leak 0.04 |
| kernel-error range over all kernels | B leak 0.01–0.30 | B leak 0.03–0.47 |
| Q2 decade | −0.006 ± 0.030 % | +0.007 ± 0.050 % |
| F1 column (near + far) | Q1 0.055 ± 0.017 % | Q1 0.091 ± 0.026 % |

## What failed, the measured obstacle, and the next adaptation
1. **Precision with the regional shallow term.**
   - Obstacle: D(e) and Sc(e) are nearly collinear under ell40, because a_d and a_s co-vary across cells. Under two3_20 they
     are much less so.
   - Next:
     - adopt two3_20 for Q3 too, since its selection is now validated (σ per major 0.099 % with Sc free);
     - or tie Sc to a smooth regional curve with a few knots, so it cannot copy D episode by episode.
2. **Prior test (fails in every form).**
   - Obstacle: the conditioning prior on the near-duplicate S-per-site and R-per-station terms.
   - Next:
     - drop S_site wherever a station's R exists (R absorbs it), or fix their null directions exactly as for D/Sc;
     - report prior-marginal σ.
3. **T4: an annual cycle common to every station.**
   - Obstacle: by construction it cannot be told from an annual pattern in the episode terms.
   - Next: project the annual component out of D(e) (ΔD then has no annual part), and test by injection.
4. **T6: kernel-error size.**
   - The deep-pulse size stays uncertain by about ±30 % with the kernel (0.72–1.40).
   - With the selection now validated, state sizes under two3_20 with the posterior-weighted recovery (0.92 SVI, 0.80 WA).
5. **SVI 2019-11.**
   - Next: per-instrument graph values at B009, B011 and the PGC vault for 2019-11-22…25 (when the carrier families fired),
     and a co-located vault-vs-borehole contrast.
   - A real deep change would also appear at the far stations. A local PGC-site change would not.
6. **WA Q1 (rule met, supporting tests disagree).**
   - Next: refit Q1 with the per-station receiver term in the smooth model. Q1 under that term was not run this round.
     - If the cycle survives with R, it is not tied to the stations.
     - If it vanishes, it follows the stations, as the far/near test suggested.

## Files
- Scripts: `svi_science/scripts/s11_r4_{common,measure,farkernel,lapse,combined,dive,figures}.py`; queue `s11_r4_queue.sh`.
- Code-test outputs (not results): `svi/tables/test/`, `*_test.*`.
- Figures (dv/v points coloured by weight, no error bars):
  - `figures/fig_r4_dD.png`: ΔD per episode, H, both regions;
  - `figures/fig_r4_suite.png`: injection rules per generating kernel;
  - `figures/fig_r4_lapse.png`: selection confusion matrices;
  - `figures/fig_r4_2019_11.png`: the dive. Its carriers panel shows the H form; the HnoSc carriers are in
    `dive_carriers_noSc.csv`.
