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

## Scale-up, SVI (after gap fill, channel guard v2) — rewritten 2026-10-14 ~7:40 pm PT

**WA-NW: withdrawn pending the origin-time fix.** Merlin F1: `s13_units.wa_shift` moves every WA origin about 16 s
early (the WA catalogue `time` is already origin-like; `s_off` is a within-template offset). Every WA number from the
scale-up (detections, zones, families, between-ETS, the SVI-vs-WA ratio) is withdrawn and not quoted here. WA is re-run
from scratch after the fix. The old WA text is kept in
`../_quarantine_2026-10-14_envfail/s13_pre_guard/RESULTS_pre_guard_2026-10-14.md`.

**What changed since the 6 pm write-up.**
- The two problems found then are fixed in the code, not patched in the tables:
  1. flat channels (UW.LRIV N/Z 2019-04/05 with control SD 0) that collapsed fits to ~0;
  2. a near-dead channel taking most of the weight (UW.SP2 E; WA only).
- Rule: PLAN.md Amendment 4 v2 (Merlin). In each unit, drop a station-component if any of these holds:
  - n_used < 5;
  - σ ≤ 0;
  - σ < 0.05 × the median σ of the station's other two components;
  - σ < 0.02 × the median σ over all channels in the unit.
- There is no weight cap. The same weights go to the data, the scrambled-control draws and the ±1-day stacks. Every
  drop is logged with its reason. A first cap rule (v1, 5:52 pm) was withdrawn before use; its outputs are in
  quarantine.
- Re-runs (v2):
  - all 29 fill-affected SVI row files as `scalef_*` on the gap-filled data;
  - the 10 other SVI row files with a v2 drop as `refit_*`;
  - svi_2005-09 has no v2 drop and keeps its first-pass table.
  - The summary picks `scalef_` > `refit_` > `scale_`: 28 / 10 / 1 episodes. No unit was excluded.
- Checks on the fix:
  - SVI 2019-04 is now 3.4e12 per LFE (was ~0).
  - Where v2 drops nothing (16 units in the refits), M0 and z equal the first pass exactly.
  - The gap fill is complete: `FILL_QUEUE_DONE`, touched by hand after the v2 SVI `scalef_*` finished.
- New columns: `max_share` (largest single-channel share of (template/σ)²; median 0.55 over episodes) and the
  leave-one-station-out M0 range (median spread 37 % of M0).
- What moved. The as-run broken headline (5.9 N m) is now 4.4e12 (all 39) and 5.4e12 (2008–26, previously 5.3e12).
  Detections went from 20/38 to 21/39. Between-ETS went from z 3.7 to 3.2 (4.4e12 per LFE).

Inputs: `tables/scale_verdict.json` (SVI; WA marked withdrawn), `tables/scale_svi_{episodes,families,zones}.csv`,
`tables/scale_between.csv`, `tables/scale_guard_excluded.csv` (empty), `tables/guard2_drops_scale_svi.csv`,
`tables/guard2_zvariant.csv`. Raw figure: `figures/s13_scale_svi_episodes.png` (moment per LFE vs time, colour = z vs
scrambled control, size ∝ √LFEs, no error bars). `figures/s13_scale_both_m0_per_lfe.png` still shows the pre-guard
values and the withdrawn WA points; do not use it.

### What the data show (S Vancouver Island)

**The LFE-locked 20–50 s pulse is in nearly every episode, at about the same size per LFE.**
- **Detections (A1):**
  - 21 of 39 built episodes reach z ≥ 3; none reaches z ≤ −3; 38 of 39 moments are positive.
  - Large episodes: 18 of the 19 with ≥ 300 LFEs are detected. The one miss is 2009-05 (z 2.9).
  - Small episodes: 3 of the 20 with < 300 LFEs are detected (2020-05, 2025-05, 2025-09). 2019-04 (57 LFEs) is now a
    normal non-detection: 3.4e12, z 1.2.
  - Strongest: 2018-06 (z 9.1), 2024-08 (9.1), 2019-08 (8.0), 2005-09 (6.6), 2012-08 (6.2), 2023-04 (6.1).
  - 2011-08 still reproduces the pilot: 5.31e12, z 5.6.
- **Moment per LFE (A2):**
  - Median 5.2e12 N m (IQR 4.0–6.7e12, ≈ Mw 2.4), about 12× the LFEs' own moment.
  - Weighted mean over all 39 is 4.4e12 (σ 0.18e12), with χ²/dof 2.7 (p 9e-8). By the pre-set rule the partition
    *changes*.
  - The change sits in the three early Bostock-era episodes (the three largest LFE counts): 2003-02 2.9e12,
    2004-07 1.9e12, 2005-09 3.7e12.
  - **2008–2026 alone (35 episodes): weighted mean 5.4e12, χ²/dof 1.09, p 0.34. Constant.** Dropping only 2004-07
    gives χ²/dof 1.39, p 0.06.
  - Likely reading (not tested): the denser 2003–05 backward catalogue lowers the per-LFE normalisation.
- **Size (A3):**
  - The episode total tracks GNSS slip: ρ 0.81, p 0.0006, n 15 (2014+ only). It also tracks tremor (ρ 0.85, n 20)
    and LFE count (ρ 0.97, n 15). It does not track duration (ρ 0.39, p 0.09).
  - The moment per LFE tracks none of these (ρ −0.13 to −0.20, p ≥ 0.40). The GNSS link runs through the LFE count.
  - For the 10 episodes with geodetic M0 > 1e18 N m, the LFE-locked total is 0.15 % of it (median; 0.04–0.33 %;
    a floor, see caveat).
- **Families and zones (A4):**
  - Families: 29 of 111 have a multi-episode weighted moment at z ≥ 3, none at z ≤ −3. In single family-episodes,
    51 of 607 reach z ≥ 3 and 1 reaches z ≤ −3.
  - Best measured (per LFE): fam 10 9.7e12 (z 15, 20 episodes, χ²/dof 0.7), fam 7 11.7e12, fam 3 10.4e12,
    fam 16 12.7e12, fam 93 12.8e12, fam 56 8.9e12.
  - Between families χ²/dof 3.0. Among the 38 well-measured families (σ < 2.5e12) the IQR is 4.9–9.6e12.
  - No clear trend among the 38 well-measured families: depth ρ −0.30 (p 0.07), latitude ρ −0.24 (p 0.14). The GF
    ×2 caveat (F5) still applies to family-to-family differences.
  - Zones (SE → NW): 5.3e12 (z 22, 20 episodes at z ≥ 3), 6.2e12 (z 18, 14), 4.9e12 (z 12, 8). About the same along
    strike.
- **Between ETS (A5):** BE_broad (688 LFEs, 2014–2026) gives **4.4e12 per LFE, z 3.2**, a pass by the rule. That is
  0.81 of the 2008–26 ETS value. The 76 tremor-free BE events give z −1.8 (too few to say anything).
- **vs dv/v (A6):** no relation. Total vs round-7 ΔD ρ 0.11; per LFE ρ 0.03 (n 26).

### Checks (SVI, guarded tables)
- **±1-day stacks** (same LFEs moved one day): SD 1.00 and 0.94 of the scrambled-control SD, mean ≈ 0, no |z| ≥ 3 in
  82 stacks (41 units). The scrambled control is calibrated.
- **Odd/even-day halves:** in 16 of 21 detected episodes both halves reach z ≥ 2, and both are positive in all 21.
  Odd-half vs even-half z across episodes: ρ 0.70 (p 1e-6).
- **Mechanism (21 detected):** thrust VR median 0.20 vs free MT 0.25. Thrust beats strike-slip in 21 of 21 and steep
  dip-slip in 20 of 21; steep dip-slip fits with the opposite sign (median z −3.8).
- **`n_nb10`:** still not meaningful in the scale-up (it counts each LFE's own copies in the half/zone/family units).

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
- **Channel guard (Amendment 4 v2).** Its station-ratio clause also drops healthy verticals (SNB Z, KHVB Z, VGZ Z)
  whose 20–50 s horizontals are tilt-noisy. That costs some sensitivity; it should not bias the moment (Z-added-back
  re-fit, `tables/guard2_zvariant.csv`: median ratio 1.00 on 40 pooled/between units, 33 within 10 %). Units resting on
  one channel are visible in `max_share` and the leave-one-station-out range; e.g. 2003-02 depends on PGC Z
  (max share 0.89; without PGC the moment goes to −1.3e13).

### Score against the pre-set rules (Amendment 1), SVI only (WA withdrawn)

| rule | SVI (guard v2, after gap fill) |
|---|---|
| A1 detection, pooled z ≥ 3 (count) | **21 / 39** (0 at z ≤ −3; 38/39 positive; 18/19 of the ≥ 300-LFE episodes) |
| A2 partition constant (χ² p < 0.01 → changes) | **changes** (4.4e12, χ²/dof 2.7, p 9e-8), driven by 2003–05; **2008–26 constant** (5.4e12, χ²/dof 1.09, p 0.34) |
| A3 size, p < 0.05 (total vs GNSS) | **PASS** ρ 0.81, p 0.0006 (n 15), through the LFE count; per LFE vs GNSS ρ −0.13 |
| A4 between families (χ²/dof) | families differ (3.0); 29/111 at z ≥ 3, 0 at z ≤ −3; zones about equal |
| A5 between ETS, z ≥ 3 | **PASS**: BE_broad 4.4e12 per LFE, z 3.2 = 0.81 × ETS; BE (76) uninformative |
| A6 vs dv/v (column) | no relation (ρ 0.11 / 0.03, n 26) |

**Plain summary.** Southern Vancouver Island LFEs carry a constant ~5e12 N m (about Mw 2.4) of 20–50 s thrust moment
each, ETS after ETS from 2008 to 2026, and about the same between ETS (0.8× the ETS value). Episode totals scale with
GNSS slip because the LFE count does. The guard fix changed the broken headline values and moved the rest only slightly.
Washington is withdrawn until the WA origin times are fixed.
