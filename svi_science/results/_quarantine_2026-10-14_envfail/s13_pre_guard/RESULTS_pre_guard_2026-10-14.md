# S13 RESULTS — long-period (20–50 s) stacks on our LFE times

## Pilot (2026-10-13, ~10:50 am PT; rules PLAN.md §3, verdict `tables/pilot_verdict.json`)

**Plain answer: yes. Our SVI LFE times carry a coherent 20–50 s pulse with the shape and polarity of a thrust on the
plate interface. It is seen in both ETS tested, with almost the same size.** Provisional (no Merlin review yet).

What the data show (pooled SVI families, P ≥ 0.75, counted once):

| episode | LFEs | VLF moment per LFE (N m) | Mw equiv. | z vs scrambled control | p | VR thrust | VR free MT | ±1 day stacks (N m) |
|---|---|---|---|---|---|---|---|---|
| 2011-08 | 1,445 | 5.3e12 | 2.4 | 5.6 | 0.002 (floor) | 0.25 | 0.26 | −0.7e12, −0.4e12 |
| 2012-08 | 2,864 | 5.1e12 | 2.4 | 7.9 | 0.002 (floor) | 0.39 | 0.40 | −0.6e12, +0.4e12 |

- Raw stacks first: `figures/s13_pilot_2011-08_pooled.png`, `figures/s13_pilot_2012-08_pooled.png`
  (black stack, red thrust template × fitted moment, grey scrambled-control 2.5/97.5 %, blue ±1 day).
  - The pulse is clearest on the verticals of the nearest quiet stations: VGZ, LRIV, PGC (2011), SNB (2012).
    - Example: VGZ Z 2.3 nm peak vs 1.6 nm control 99 %.
    - LRIV Z 1.8 nm vs 0.9 nm (2012).
    - The up–down shape and timing follow the template.
  - Horizontals are 5–40× noisier (tilt), and far stations (NLLB, YOUB, 64–100 km) predict < 0.2 nm. Both behave as
    the physics expects.
- The thrust explains essentially everything a free deviatoric moment tensor explains (VR 0.25 vs 0.26; 0.39 vs 0.40).
  - Strike-slip fits poorly (VR 0.04 / 0.14).
  - A steep dip-slip fits with the opposite sign (z −4.8 / −5.8).
- Odd-day and even-day halves agree: 5.0e12 (z 3.2) and 5.4e12 (z 4.9) in 2011; 5.0e12 (z 5.5) and 5.3e12 (z 5.5) in 2012.
- The P ≥ 0.5 set gives the same moment per LFE (4.6e12 and 4.9e12; z 6.1 and 9.5).
- Size: 5e12 N m per catalogued LFE is ~12× the LFE's own moment (Bostock mag_med 1.7 ≈ 4.5e11 N m).
  - Catalogued neighbours within ±10 s are only 0.16–0.19 per LFE, so the moment is not an overlap of catalogued LFEs.
  - It is the slow moment that moves with the LFEs: Ide-type VLF moment.
  - The absolute scale uses a whole-space Green's function with a free-surface factor 2, good to about ×2.
  - Relative values between episodes and families are the real product.
- Per family (~60–110 LFEs per episode):
  - fam 10, 3, 0, 56 reach z ≥ 3 in at least one episode (fam 10: z 3.5 and 6.7; fam 3: 4.6 and 3.8).
  - fam 43 does not (z 1.5 and 2.6).
  - Family moments per LFE range 6e12–1.7e13. Per-family-per-episode numbers are therefore possible for the busier
    families.
- Cross-check with Hutchison & Ghosh 2018 (`tables/hg2018_crosscheck.csv`): around their four 2011 SVI-area VLFE
  templates, our LFEs within ±10 min number 14, 10, 3, 8. The day-rate expectation is 4.4, 4.4, 1.7, 5.5, so about
  2× more LFEs at their VLFE times (small numbers).

Pass rules:

| rule | 2011-08 | 2012-08 | result |
|---|---|---|---|
| P1 detection | yes | yes | PASS |
| P2 thrust sign and pattern | all parts pass | VR, mechanism and sign pass; sign agreement on the larger template components 0.67 (rule ≥ 0.70) | PASS in 2011-08 |
| P3 halves | yes | yes | PASS |
| P4 families (≥ 2 of 5 at z ≥ 3) | — | — | PASS (4 of 5) |

**Pilot works** under the pre-set definition (P1 + P2 + P3 in ≥ 1 episode).

Caveats and known problems:
- PGC dropped out of 2012-08: its 2012 files have 6–12 s gaps every few minutes, and the pilot prep split on every gap.
  Fixed in Amendment 2 for the scale-up.
- LZB is 30–150× noisier than its neighbours at 20–50 s (response or site). It is weighted out, not removed.
- Control p-values sit at the 1/501 floor.
- The test runs (`tables/test_*`, `STACK_test_*`, `figures/test_*`) used 20 control draws. They are code checks, not
  results.

## Scale-up — launched 2026-10-13 10:57 am PT (Amendments 1–2)
Queue `scripts/s13_scale_queue.sh`, PID in `logs/scale_queue.pid`, log `logs/scale_queue.log`. Done when
`SCALE_QUEUE_DONE` exists. It covers:
- SVI: 39 majors 2003–2026;
- WA-NW: 56 majors 2007–2026;
- between-ETS sets: SVI BE 76 / BE_broad 688; WA BE 220 / BE_broad 934.

Then write this section from `tables/scale_verdict.json`, `tables/scale_{svi,wa}_episodes.csv`, `scale_*_families.csv`,
`scale_between.csv`, and `figures/s13_scale_*.png`.

## Scale-up (before gap fill) — written 2026-10-14 ~6 pm PT

Inputs: `SCALE_QUEUE_DONE` run of 2026-10-14 5:39 pm PT (all SVI stacks redone after the 3:40 pm PGC sample-rate fix).
`tables/scale_verdict.json`, `tables/scale_{svi,wa}_episodes.csv`, per-episode `scale_<region>_<ep>_{units,stacomp}.csv`,
`scale_between.csv`. **Before gap fill:** no `scalef_*` table existed when the summary ran. The fill queue
(Amendment 3, started 5:39 pm PT) will re-stack some episodes tonight, and these numbers will move a little.
Coverage check: 40 SVI stacks (39 episodes + between-ETS) and 57 WA stacks (56 + between-ETS) are all on disk. Every
built episode has a row with a value (13 SVI and 13 WA majors have < 30 LFEs at P ≥ 0.75 and were not built, as
pre-set). There are no NaN rows from missing stacks.

Two problems found while writing this up. They change some verdict numbers, so they come first:
1. **Flat channel → degenerate fits.**
   - UW.LRIV N and Z are flat (control SD exactly 0) in 2019-04 and partly in 2019-05; UW.OPC is flat in WA 2009-04.
   - A zero SD gives that channel infinite weight, so the fitted moment collapses to ~0. Affected: 27 of 1,623 units
     (every unit of SVI 2019-04 and WA 2019-04, WA 2009-04 zones 2–3 and some families, WA 2019-05 zone 1).
   - These units make the as-run A2 values meaningless: `scale_verdict.json` gives a "weighted mean moment per LFE"
     of 5.9 N m (SVI) and 707 N m (WA). Those are the 2019-04 values (5.9 N m × 57 LFEs, 707 N m × 376 LFEs), not
     the ~1e12 N m ones. They also distort the family A4 numbers.
   - All A2/A4 numbers below drop units with control SD < 1e10 N m. The pre-set A1/A3/A6 counts are unaffected
     except WA 2019-04, where the as-run "detection" (z 3.9) comes from this artefact.
   - `s13_scale_summary.py` and `s13_stack2.py` have no guard against this. The fill-queue summary will repeat it
     unless a flat-channel guard is added.
2. **UW.SP2 E is near-dead in 2012-04 → 2019-03.**
   - Its control SD is 0.04–0.35 nm, against 7–20 nm on SP2 N, so it carries 33–70 % of the fit weight in 19 WA
     episodes and in both WA between-ETS sets.
   - A near-zero channel pulls the moment toward 0. To first order it scales the control SD the same way, so z is
     about unchanged.
   - Re-fit without it (stored stacks, same weights otherwise; `scripts/s13_writeup_checks.py` →
     `tables/writeup_bgf_refit.csv`; the all-channel re-fit reproduces every stored moment exactly): those WA
     moments rise 1.2–2.6×. For example, 2012-08 goes from 1.2 to 3.2e12, 2017-02 from 1.4 to 3.1e12, and
     WA_BE_broad from 0.96 to 1.5e12.
   - New z needs the control draws, which are not stored, so these re-fits are a column, not a verdict.

Raw figure first: `figures/s13_scale_both_m0_per_lfe.png` (moment per LFE per episode vs time, both regions, colour = z
vs scrambled control, size ∝ √LFEs, no error bars; `scripts/s13_scale_fig_both.py`). Per-region figures from the
summary are `figures/s13_scale_{svi,wa}_episodes.png` and `_vs_size.png`. They still plot the degenerate 2019-04
point near 0.

### What the data show

**S Vancouver Island: the LFE-locked 20–50 s pulse is in nearly every episode, at about the same size per LFE.**
- **Detections (A1):** 20 of 38 valid episodes reach z ≥ 3, and none reaches z ≤ −3.
  - All 38 have a positive moment (38 of 38 > 0).
  - Of the 19 episodes with ≥ 300 LFEs, 17 are detected. The two misses are 2003-02 (z 2.9) and 2009-05 (z 2.9).
  - Of the 19 small episodes (31–286 LFEs), 3 are detected: 2020-05, 2025-05, 2025-09.
  - Strongest: 2018-06 (z 9.5), 2024-08 (9.1), 2019-08 (7.8), 2012-08 (7.8; re-run with PGC, 4.7e12), 2005-09 (6.6),
    2023-04 (6.1).
  - 2011-08 reproduces the pilot: 5.31e12 vs 5.33e12, z 5.6.
- **Moment per LFE (A2):**
  - Median 5.1e12 N m (IQR 4.1–6.1e12, ≈ Mw 2.4), about 12× the LFEs' own moment.
  - Weighted mean over all 38 episodes is 4.4e12, with χ²/dof 2.7 (p 1e-7). By the pre-set rule (p < 0.01) the
    partition therefore *changes*.
  - The change sits in the three early Bostock-era episodes, which are also the three largest LFE counts:
    2003-02 2.8e12 (4,034 LFEs), 2004-07 1.9e12 (15,692), 2005-09 3.7e12 (10,406).
  - **2008–2026 alone (34 episodes): weighted mean 5.3e12, χ²/dof 1.03, p 0.42. The slow moment per LFE is
    constant across 18 years of ETS.** Dropping only 2004-07 gives χ²/dof 1.4, p 0.05.
  - Likely reading (not tested): the 2003–05 backward catalogue counts more LFEs for the same slow slip, so the
    per-LFE normalisation depends on catalogue density. It is not evidence that the slip changed.
- **Size (A3):**
  - The episode total (moment per LFE × n) tracks GNSS slip: ρ 0.81, p 0.001, n 15 (all 2014+, five of them with
    < 0.5 mm offset). It also tracks tremor (ρ 0.86, p 0.0002, n 20) and LFE count (ρ 0.97, n 15). It does not
    track duration (ρ 0.37, p 0.10, n 20).
  - The moment per LFE does not track any of them (ρ −0.11 to −0.21, p ≥ 0.36).
  - Plainly: bigger episodes have more LFEs, each LFE carries the same slow moment, and the GNSS correlation comes
    through the LFE count. n = 15 is small, and the GNSS series covers only 2014+.
  - Fraction of geodetic moment: for the 10 episodes with Tier-B M0 > 1e18 N m, the LFE-locked 20–50 s total is
    0.15 % of the GNSS moment (median; range 0.04–0.30 %).
- **Families and zones (A4):**
  - 31 of 111 families have a multi-episode weighted moment at z ≥ 3, and none at z ≤ −3. In single
    family-episodes, 61 of 605 reach z ≥ 3 and 1 reaches z ≤ −3 (about 0.8 expected by chance in each tail).
  - Best measured: fam 10 9.5e12 (z 15, 20 episodes, χ²/dof 0.73), fam 7 11.7e12, fam 3 9.9e12, fam 16 12.5e12,
    fam 93 13.2e12, fam 56 9.5e12 per LFE.
  - Families differ: between-family χ²/dof 3.0 (p 1e-23). This still holds after inflating σ by the within-family
    scatter (χ²/dof 1.55, so σ × 1.25): χ²/dof ≈ 1.9. Among the 41 well-measured families (σ < 2.5e12) the IQR is
    4.8–9.4e12.
  - No trend with depth (ρ −0.27, p 0.09) or latitude (ρ −0.14, p 0.37).
  - The family values sit above the pooled 5e12 because each family is fit at its own location rather than the
    centroid. Relative station amplitudes carry the whole-space ×2 uncertainty (Merlin F5), so family-to-family
    differences of < 2× are not yet safe.
  - Zones (along-strike terciles, SE → NW): 5.2e12 (z 23, 19 episodes at z ≥ 3), 6.1e12 (z 18, 13), and 4.9e12
    (z 13, 8). Between zones χ²/dof 3.4 on 2 dof, p 0.03: about the same along strike.
- **Between ETS (A5):**
  - **The 688 between-ETS LFEs (BE_broad, 2014–2026, forward catalogue) carry slow moment: 4.2e12 N m per LFE,
    z 3.7.** The ±1-day stacks give −1.5e12 and −0.3e12, and VR is 0.08.
  - That is 0.79 of the 2008–2026 ETS value (5.3e12): between-ETS LFEs carry about the same slow moment per LFE as
    ETS LFEs.
  - The 76 tremor-free BE events give z −2.0. Too few: the minimum detectable moment is 3σ = 1.3e13 per LFE, so
    this set says nothing either way.
- **vs dv/v (A6):** no relation. Episode total vs round-7 ΔD: ρ 0.08 (n 26). Moment per LFE vs ΔD: ρ −0.02.

**Washington NW: a weaker, less clean signal, and it lives in the SE third of the families.**
- **Detections (A1):** 14 of 56 episodes reach z ≥ 3 as run; 13 of 55 once the flat-channel 2019-04 is dropped.
  4 episodes sit at z ≤ −3 (2014-06, 2014-11, 2018-08b1, 2020-09).
  - Positive moment in 34 of 55. Among the 30 episodes with ≥ 300 LFEs: 12 detected, 21 positive.
  - Detected: 2009-04 (z 5.0), 2017-02 (5.1), 2011-07 (4.3), 2021-09 (4.2), 2012-08, 2013-09, 2014-07, 2018-05,
    2019-03, 2022-10b1, 2024-01, 2024-08, 2025-09 (3.0–4.1).
- **The pooled WA unit mixes signs; zone 1 does not.**
  - Zone 1 (SE tercile, centred 47.74 N, −122.61, about 33 km; families 106, 38, 126 near 47.65 N −122.8 at 43 km)
    gives 1.5e12 per LFE. Its weighted z is 16, with 16 of 47 episodes at z ≥ 3 and none at z ≤ −3.
  - Zones 2 and 3 (toward the Olympics and the Strait) give −0.4e12 each (z −2.8, −2.7). They have no episode at
    z ≥ 3, and zone 3 has 3 at z ≤ −3.
  - Several pooled "negative" episodes are positive in zone 1: 2020-09 (pooled z −3.6, zone 1 +4.2), 2023-03
    (−2.6 / +4.2), 2008-04 (−2.2 / +5.7).
  - The consistently negative NW zones point to a template or location mismatch there (WA family locations,
    origin shift, or GF amplitude pattern) rather than to physics. This is a hypothesis; it is not tested.
- **Moment per LFE (A2):** detected-episode median 1.5e12 (1.9e12 with SP2 E removed). The weighted mean over all
  55 is 0.64e12, with χ²/dof 5.5, so the partition *changes* by the rule. The WA scatter is not trustworthy, though
  (sign mixing, SP2 E, wide controls; see checks). Over the 13 detected episodes: 1.4e12, χ²/dof 1.2, p 0.28.
- **Size (A3):** the episode total tracks tremor, the headline (ρ 0.55, p 0.009, n 22; as run 0.55, p 0.010, n 23).
  It also tracks LFE count (ρ 0.43, p 0.04, n 24). It does not track duration (ρ 0.26, n 24) or the WA-station GNSS
  offset (ρ 0.13, n 14). The moment per LFE tracks nothing (|ρ| ≤ 0.17).
- **Families (A4):** 8 of 72 at z ≥ 3 and 2 at z ≤ −3 (fam 52, fam 37). The best are fam 106 (4.0e12, z 13.5) and
  fam 38 (3.6e12, z 12.4). These two are nearly co-located, so possibly one source. Between families χ²/dof 7.7,
  within families 2.5: the WA controls under-state the scatter.
- **Between ETS (A5):** WA_BE_broad (934 LFEs, 2015–2026) gives 0.96e12 per LFE at z 3.2 (1.5e12 without SP2 E), a
  pass by the pre-set rule. Scaled by the WA ±1-day control width (below), it is only about 2.2–2.5σ, so it counts
  as marginal. WA_BE (220) gives z 0.8.
- **vs dv/v (A6):** total vs ΔD ρ −0.06 (n 24); moment per LFE vs ΔD ρ −0.37 (p 0.07). With SP2 E removed,
  ρ −0.42 (p 0.035), but that is a post-hoc column and not a finding.

**SVI vs WA.** SVI LFEs carry about 3× more 20–50 s moment per LFE than WA LFEs (5.3e12 vs about 1.5–1.9e12; WA zone
1 gives 1.5e12 and its best families 3.6–4.0e12). Part of the gap may be method rather than physics:
- WA sign mixing across zones;
- the SP2 E down-weighting;
- the WA catalogue's different P (P_own_pnet) and family definitions.
SVI is the clean, constant result; WA confirms the signal exists in the SE third of the NW families.

### Checks
- **±1-day stacks** (same LFEs moved to the next or previous day; ratio to the scrambled-control SD):
  - SVI: SD 0.94 and 0.90, mean 0.0, no |z| ≥ 3 in 76 stacks. The scrambled control is calibrated.
  - WA: SD 1.29 and 1.48, two of 110 at |z| ≥ 3 (max 5.3). The WA scrambled control under-states the noise by
    30–50 %, so WA z values should be read about 1.4× smaller. At z ≥ 4.2, WA keeps 4 detections (2009-04, 2011-07,
    2017-02, 2021-09) and 0 negatives.
- **Odd/even-day halves:**
  - SVI: in 16 of 20 detected episodes both halves reach z ≥ 2, and both halves are positive in all 20. Across all
    episodes the odd-half and even-half z agree (ρ 0.69, p 0.0002).
  - WA: 4 of 13 have both halves at z ≥ 2 and 7 of 13 have both positive. Odd vs even ρ is 0.04. The halves do
    not replicate within WA episodes, which fits the zone sign mixing.
- **Mechanism (detected episodes):**
  - SVI: thrust VR median 0.23 vs free moment tensor 0.26. Thrust beats strike-slip in 20 of 20 and steep dip-slip
    in 19 of 20, and the steep dip-slip fits with the opposite sign (median z −4.4), as in the pilot.
  - WA: VR 0.06 vs 0.10; thrust beats strike-slip 13 of 13 and steep dip-slip 9 of 13.
- **PGC in/out:** not available (no such unit in this run).
- **Neighbour column `n_nb10`:** not meaningful in the scale-up. `s13_stack2.py` counts neighbours over the whole row
  file, where each LFE appears again in its half, zone and family units, so it reads 2–4 (SVI) instead of the
  pilot's 0.16–0.19.

### Caveat (Merlin F4) and scope
- **Lower bound, by construction.** Window QC drops any window whose 20–50 s RMS exceeds 3× the median. That
  removes windows containing a discrete VLFE above about Mw 3.6 (VGZ) to 3.8 (PGC). The "slow moment per LFE" here
  is therefore the small, continuous, LFE-locked accompaniment with the large discrete VLFEs cut out: a lower bound
  on the time-locked VLF moment. It bears directly on A3. The 0.15 % of the GNSS moment is a floor, also limited to
  the 20–50 s band and to catalogued P ≥ 0.75 LFEs. These numbers cannot predict how many discrete VLFEs occur
  (that is the between-ETS scan, `MERLIN_VLFE_SCAN.md`).
- **Absolute scale:** whole-space Green's function × 2 free surface, good to about ×2. Relative values between
  episodes are the product. Between families or regions they are less safe, because relative station amplitudes
  carry order-1 GF errors (F5).
- **Scope of "constant":** episode-vs-episode stacks over each whole ETS (2–4 weeks), at the LFE source area
  (30–45 km, SVI families), 2008–2026, at a resolution of about ±1e12 N m per LFE per large episode. It says nothing
  about changes within an episode or on time scales shorter than one ETS.
- **Before gap fill.** The `scalef_*` restacks (fill queue running now) will refresh some episodes. The flat-channel
  guard should go into the summary first.

### Score against the pre-set rules (Amendment 1)

| rule | SVI | WA-NW |
|---|---|---|
| A1 detection, pooled z ≥ 3 (count) | **20 / 38** (0 at z ≤ −3; 38/38 positive; 17/19 of the ≥ 300-LFE episodes) | **13 / 55** (14/56 as run, incl. the flat-channel 2019-04); 4 at z ≤ −3; ~4 at the ±1-day-scaled z ≥ 4.2 |
| A2 partition constant (χ² p < 0.01 → changes) | **changes** (χ²/dof 2.7, p 1e-7; as-run verdict value broken by 2019-04). Driven by 2003–05; **2008–26 constant** (5.3e12, χ²/dof 1.03, p 0.42) | **changes** (χ²/dof 5.5); not interpretable (sign mixing, controls too narrow); detected-only constant (1.4e12, p 0.28) |
| A3 size, p < 0.05 (headline SVI total vs GNSS; WA total vs tremor) | **PASS** ρ 0.81, p 0.001 (n 15), but through the LFE count; per LFE vs GNSS ρ −0.11 | **PASS** ρ 0.55, p 0.009 (n 22) |
| A4 between families (χ²/dof) | families differ (3.0, p 1e-23; ≈ 1.9 after σ inflation); 31/111 families at z ≥ 3; zones about equal (p 0.03) | differ (7.7); 8/72 at z ≥ 3, 2 at z ≤ −3; signal in zone 1 only (z 16); zones 2–3 slightly negative |
| A5 between ETS, z ≥ 3 | **PASS**: BE_broad 4.2e12 per LFE, z 3.7 = 0.79 × ETS; BE (76) uninformative | **PASS by the rule** (BE_broad z 3.2, 0.96e12), marginal after the ±1-day scaling (~2.3σ) |
| A6 vs dv/v (column) | no relation (ρ 0.08 / −0.02, n 26) | no relation (ρ −0.06 / −0.37, p 0.07, n 24) |

**Plain summary.** Southern Vancouver Island LFEs carry a constant ~5e12 N m (about Mw 2.4) of 20–50 s thrust moment
each, ETS after ETS from 2008 to 2026, and also between ETS (0.8× the ETS value). Episode totals scale with GNSS slip
because the LFE count does. In Washington the signal is clear only in the SE third of the NW families, at about a third
of the SVI size per LFE. The NW zones need a template/location check before WA numbers are trusted.
