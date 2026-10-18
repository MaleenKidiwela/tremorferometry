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

**WA-NW: the first-pass WA numbers are withdrawn. The corrected re-run (Amendment 5) is reported in the next
section; it fails its pre-set lag rule, so it is not a headline either (see "Amendment 5 results").** Original note: Merlin F1: `s13_units.wa_shift` moves every WA origin about 16 s
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
- **`n_nb10`:** the first-pass column was not meaningful (it counted each LFE's own copies in the half/zone/family
  units). Fixed in Amendment 5 (step 6): see the next section (median 0.06 per LFE within ±10 s; 0.38 in 2004-07).

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

**Plain summary (superseded in wording by the provisional statement in "Amendment 5 results" below).** Southern
Vancouver Island LFEs carry ~5e12 N m (about Mw 2.4) of 20–50 s thrust moment each. By the pre-set rule the per-LFE value
*changes* between episodes; the change sits in 2003–05, and 2008–2026 alone looks constant (a post-hoc split).
Between-ETS LFEs carry about 0.8× the ETS value. Episode totals scale with GNSS slip because the LFE count does.


## Amendment 5 results (Merlin steps 2–6, Q1c) — written 2026-10-14 ~8:15 pm PT

Rules: PLAN.md Amendment 5 (6:30 pm PT) and its addendum (UW.SP2 E instrument veto 2011-09-05 → 2019-05-22, 7:35 pm PT),
both written before the numbers. Summary: `scripts/s13_v2_summary.py` → `tables/v2_verdict.json`, `tables/v2_*.csv`.
New code only (`s13_units_v2.py`, `s13_stack2_v2.py`; the v2 stacker reproduces `s13_stack2.py` exactly with kmax 1).

### SVI — provisional statement (Merlin step 7)
Raw figure first: `figures/s13_scale_svi_episodes.png` (moment per LFE per episode, colour = z vs scrambled control, no
error bars).
1. **Pre-set verdict (A2): the per-LFE moment *changes* between episodes**: 4.4e12 N m over all 39, χ²/dof 2.7, p 9e-8.
2. **Post-hoc split (labelled as such):** 2008–2026 (35 episodes) gives 5.4e12, χ²/dof 1.09, p 0.34, i.e. no scatter
   beyond the noise. The low values are the three 2003–05 episodes (2003-02 2.9e12, 2004-07 1.9e12, 2005-09 3.7e12).
   They coincide with a catalogue regime change: 336–628 LFEs/day in 2003–05 vs 40–120/day in 2007–13. Step 4 below
   tests this and does not support a simple dilution.
3. **Resolution:** about ±15–20 % per large episode. That holds only for the ~12 episodes whose control SD is ≤ 1e12
   (9 of them in 2008–26: weighted mean 5.3e12, χ²/dof 0.69). The ~20 small episodes (< 300 LFEs) say nothing
   individually.
4. **PGC dependence (Q1c).** PGC (Z in all but one) is the largest single channel in 22 of 39 episodes, holding 23–89 % of the
   template-weighted energy (`max_share`); ≥ 45 % in 16 of them (e.g. 2003-02 89 %, 2009-05 88 %, 2008-05 78 %, 2018-06 67 %).
   - Leaving CN.PGC out (refit of the stored guarded stacks, σ scaled by the measured control/white-noise ratio
     6.3): the 2008–26 weighted mean is 5.65e12 vs 5.39e12 with PGC (+5 %); all 39: 4.86e12 vs 4.44e12 (+9 %). Both
     are below Merlin's 30 % bar.
   - Single episodes that rest on PGC flip when it is removed: 2003-02 −1.3e13 (z −1.5), 2008-05 −5.1e12 (z −1.3).
     PGC-out: 35/39 positive, 13 at z ≥ 3. Table `tables/v2_svi_pgc_out.csv`.
5. **Lower bound (F4).** Window QC removes windows that contain a discrete VLFE above Mw ~3.6–3.8. The number is the
   continuous LFE-locked accompaniment, a lower bound.
6. **What the number means (Q4).** It is the kernel-weighted 20–50 s thrust moment within ~±15 s of each catalogued
   P ≥ 0.75 LFE.
   - (a) The band-limited VLF-to-LFE moment ratio is ≈ 12 (5e12 vs 4.5e11). To first order this is independent of
     catalogue completeness: the fixed neighbour count is 0.06 per LFE within ±10 s (median), and the
     kernel-weighted neighbour sum is slightly negative.
   - (b) As a moment rate, ≈ 5e12 / ~30 s ≈ 2e11 N m/s at LFE times. The episode-mean geodetic rate is
     1.3–1.9e12 N m/s. So the LFE-locked 20–50 s fluctuation is ≈ 10 % of the mean SSE rate: most slip at the LFE
     sites is released smoothly at these periods, in discrete VLFEs cut by F4, or outside the band.
   - (c) Summed over an episode it is 0.1–0.3 % of the GNSS moment: a floor, and ×2 in absolute scale. It must not
     be read as "0.15 % of the slip is VLF".
7. **Between ETS:** say "LFEs in inter-ETS tremor bursts" (the set is LFEs outside episodes ± 2 d with local tremor),
   not "LFEs without slip". Composition and jackknife in step 5.
8. **Channel guard note.** guard2's station-ratio clause (σ < 0.05 × the station's other two components) also drops
   healthy quiet verticals whose 20–50 s horizontals are tilt-noisy (CN.SNB Z in 266 SVI units, C8.KHVB Z 177,
   CN.VGZ Z 101, C8.SHDB Z 33). Adding them back leaves M0 unchanged: median ratio 1.00, 33 of 40 pooled/between
   units within 10 % (`tables/guard2_zvariant.csv`). It costs sensitivity, not bias. It is also why SNB and KHVB are
   missing from the single-channel table below.

### SVI step 4 — 2003–05 P bins and thinning
Raw figure: `figures/s13_v2_svi_pbin.png` (moment per LFE by P bin, colour = z, no error bars). Table
`tables/v2_svi_pbin.csv` (thrust, 200 draws, same stations and data as the first pass).

| episode | P .75–.90 | P .90–.97 | P .97–1 | high/low |
|---|---|---|---|---|
| 2004-07 | 1.2e12 (z 2.1, n 5,715) | 2.0e12 (3.7, 7,130) | 3.0e12 (3.4, 2,847) | 2.4 |
| 2005-09 | 3.4e12 (4.3, 4,572) | 3.8e12 (4.8, 4,363) | 4.4e12 (3.0, 1,471) | 1.3 |
| 2012-08 | 3.1e12 (2.8, 1,504) | 5.0e12 (3.3, 774) | 7.1e12 (4.7, 586) | 2.3 |
| 2013-09 | 5.5e12 (4.1, 771) | 4.3e12 (2.6, 351) | 8.8e12 (4.3, 237) | 1.6 |

2004-07 thinning: top 3,000 by P 3.2e12 (z 3.7); random 3,000 2.4e12 (z 2.8).

- **Pre-set outcome: (c), "the accompaniment scales with LFE detection strength".** The high/low ratio is ≥ 1.5 in a
  2003–05 episode (2004-07: 2.4) and in a 2012–13 episode (2012-08: 2.3).
  - (a) "dilution by weak detections" fails: 2004-07 top-P is 3.0–3.2e12, not ≈ 5e12.
  - (b) "genuine era difference" fails only narrowly: the 2004-07 high bin is 3.03e12, just over the pre-set 3e12
    cap.
- What the bins show plainly:
  - (i) In both eras the high-P third carries about 2× the moment per LFE of the low-P third. Each ratio is only
    1.6–2.2 σ, and the 2013-09 bins are not monotonic.
  - (ii) At the same P, 2003–05 still reads lower than 2012–13: high bin 3.0 / 4.4e12 in 2004-07 / 2005-09 vs
    7.1 / 8.8e12.
  - (iii) 2004-07 is not low-P-heavy (36 % of its LFEs in the lowest bin vs 52 % in 2012-08). So the catalogue mix
    of weak detections does not explain the 2003–05 drop.
- Reading, provisional: the drop is partly an era effect at fixed P, and per-LFE moment rises with P. P is a
  detection probability, a proxy for LFE strength. Testing "scales with strength" properly needs LFE amplitudes, not P.

**Single-channel Z, 2003–2013** (pooled unit, each station's vertical alone at the unit's best lag, 1e12 N m per LFE;
σ = channel control σ × the measured control/white-noise ratio; channels with σ > 3e12 omitted; full table
`tables/v2_svi_single_Z_2003_2013.csv`):

| episode | CN.PGC | CN.VGZ | UW.LRIV | C8.TWBB | C8.SHDB | CN.NLLB |
|---|---|---|---|---|---|---|
| 2003-02 | 2.6 ± 1.0 | — | — | — | — | — |
| 2004-07 | 2.1 ± 0.5 | 4.0 ± 2.6 | — | 5.0 ± 1.9 | 1.3 ± 0.8 | −0.8 ± 1.9 |
| 2005-09 | 3.7 ± 0.7 | — | — | 4.3 ± 2.4 | — | — |
| 2008-05 | 4.9 ± 1.5 | — | — | — | — | — |
| 2009-05 | 3.8 ± 1.7 | — | — | — | — | — |
| 2010-08b1 | 4.2 ± 1.4 | — | 6.2 ± 1.7 | — | — | — |
| 2011-08 | 5.6 ± 2.3 | 5.6 ± 1.8 | 4.6 ± 1.7 | — | — | — |
| 2012-08 | 3.9 ± 1.5 | 4.5 ± 1.3 | 5.7 ± 1.3 | — | — | — |
| 2013-09 | 7.1 ± 2.6 | 2.7 ± 2.2 | 5.9 ± 1.4 | — | — | — |

PGC alone carries the 2003–05 low values (2.1–3.7). The independent 2004-07 channels split: TWBB 5.0 ± 1.9 and VGZ
4.0 ± 2.6 are near the later value, SHDB 1.3 ± 0.8 is low. So the 2003–05 drop rests mainly on PGC Z, which has one
sensor epoch 2000–2015. That is one more reason to keep it provisional. (± values are channel noise only; this table
exists to show which channel carries what. It is not a fit result.)

### SVI step 5 — inter-ETS tremor-burst LFEs (BE_broad, 688 LFEs, 2014–2026)
Table `tables/v2_svi_between.csv` (P = P_same, 500 draws). Burst composition: 36 % of the events fall on 5 days
(2020-09-13 77, 2015-10-08 64, 2015-10-10 42, 2014-03-23 33, 2020-09-14 32; `tables/svi_be_burst_days.csv`). The
2020-09-13 burst is 12 days before the WA 2020-09-25 ETS.

| set | n | M0 per LFE | z | ratio to ETS (2008–26, 5.39e12) |
|---|---|---|---|---|
| all | 688 | 4.3e12 | 3.2 | 0.80 ± 0.25 |
| minus 2020-09-13/14 | 579 | 4.1e12 | 2.9 | 0.75 ± 0.26 |
| minus 2015-10-08/10 | 573 | 4.1e12 | 2.8 | 0.75 ± 0.27 |
| minus top-5 days | 440 | 4.4e12 | 2.7 | 0.81 ± 0.31 |
| tau_next ≥ 30 d | 488 | 2.6e12 | 1.6 | 0.49 ± 0.31 |
| 2014–2019 | 357 | 3.4e12 | 2.2 | 0.63 ± 0.29 |
| 2020–2026 | 331 | 4.0e12 | 2.6 | 0.74 ± 0.28 |

P-weighted, all: 4.5e12.

- **Pre-set pass:** all variants positive; minus-top-5 z 2.7 ≥ 2; both year halves positive. **PASS**: no single
  burst carries the signal.
- The set ≥ 30 days before the next ETS is weaker (0.49 ± 0.31, z 1.6). That is not a rule, and it is not
  distinguishable from the full set at 1σ.
- Ratio wording: inter-ETS burst LFEs carry 0.80 ± 0.25 of the ETS value per LFE (0.3–1.3 at 2σ). That is consistent
  with "the same" but cannot exclude half.

### Step 6 — neighbour counts (fixed)
Pooled rows only, self excluded (`tables/v2_svi_nb_fixed.csv`, from the row files).
- SVI median per LFE: 0.06 within ±10 s, 0.17 within ±20 s, 0.58 within ±60 s. Busiest episodes: 2004-07 0.38,
  2005-09 0.29, 2003-02 0.19.
- The kernel-weighted sum (20–50 s band-pass autocorrelation) is −0.05 (median) and −0.23 in 2004-07. Neighbours
  fall mostly on the kernel's negative side lobe, so they do not inflate the per-LFE moment; if anything they lower
  it slightly in the dense 2003–05 catalogue.

### WA-NW — corrected re-run (origin = time + 0.5 s, deduped ±2 s, Slab2 depths, SP2 E vetoed)
Raw figures:
- `figures/s13_v2_wa_episodes.png`: headline ±8 s grid; moment per LFE, zones and best lag per episode; colour = z.
- `figures/s13_v2_wa_c6.png`: the pre-registered column grid −14…+2 s and the grand-stack lag scans.

Tables: `tables/v2_wa_episodes.csv`, `v2_wa_units_all.csv`, `v2_wa_families.csv`, `v2_wa_between.csv`,
`v2_wa_{zone1,pooled}_grand_lagscan.csv`.

Inputs:
- 28,537 of 37,061 selected NW rows survive the ±2 s cross-family dedupe.
- 32 of 132 families take the Slab2 depth.
- 54 majors are built (25 with ≥ 300 LFEs).
- UW.SP2 E is blanked 2011-09-05 → 2019-05-22.

**Pre-set rule W1 fails.** On the headline ±8 s grid, the pooled best lags split into two groups.
- 30 of 54 sit at −4…−8 s with positive moment. 13 are pinned at +8 s, all with *negative* moment (the kernel's
  opposite-sign side lobe). The mode is +8, and 33 % of units sit at ±8 (rule: mode within ±4, ≤ 10 % at ±8).
- Failure branch as written: the grand stack over all episodes, scanned over ±40 s, peaks at **−7 s for zone 1
  (4.0e12 per LFE)** and **−5 s pooled (2.6e12)**. The opposite lobe is at about +8…+10 s (−3.7e12). The residual
  origin offset is therefore **about −5 to −7 s**: the physical origin is ~5–7 s before the catalogue `time` + 0.5.
  That is what the pre-run near-station check predicted (clear S arrivals → origin 5–8 s before `time`).
- By the rule as written, **WA stops here and its numbers stay withdrawn as a headline.** The ±8 s grid around
  `time + 0.5` reaches the wrong-sign lobe at its +8 s edge, which is the method property Merlin flagged (a ±15 s
  timing error flips the sign in this band).

**Other rules and A1–A6 on the headline grid** (reported, not claimed):
- W2 nb ±2 s: 0.00 (pass by construction). Neighbours within ±10 s remain: 0.81 per LFE.
- W3: ±k-day z SD 1.14 over 540 shifted stacks (pass, at the edge; mean 0.04; 1 at |z| ≥ 3).
- W4 zones: positive in 17/23, 20/24, 18/24 of the big episodes (pass).
- W5 odd/even ρ: 0.45 (pass).
- A1: 19/54 at z ≥ 3, 4 at z ≤ −3, 34/54 positive.
- A2: 1.8e12, χ²/dof 6.9.

**Column grid c6 (lag −14…+2 s = origin ≈ time − 5.5 s ± 8 s, each unit with its own scrambled-control z;
pre-registered as a column before the run).** The WA signal is clean there:
- **Detection:** 21/54 episodes at z ≥ 3 (18/25 of the ≥ 300-LFE episodes), **0 at z ≤ −3, 51/54 positive**.
- **Moment per LFE:** weighted mean **2.7e12 N m per LFE** (median 2.9e12), χ²/dof 1.53.
- **Zones** (all positive, none at z ≤ −3):
  - zone 1 (SE): 4.1e12, positive in 23/23 big episodes, 21 unit-episodes at z ≥ 3;
  - zone 2: 2.0e12, 22/24, 3 at z ≥ 3;
  - zone 3 (NW): 3.3e12, 22/24, 8 at z ≥ 3.
  The first-pass "signal only in zone 1, zones 2–3 negative" was the origin error.
- **Halves:** odd/even ρ 0.57.
- **Size:** the episode total tracks tremor (ρ 0.84, p 0.0002, n 23), LFE count (ρ 0.91) and the WA GNSS offset
  (ρ 0.62, p 0.016, n 15). The moment per LFE tracks none of them (|ρ| ≤ 0.2), except GNSS at ρ 0.49 (n 15, p > 0.05).
- **Families:** 15 of 59 at z ≥ 3, 1 at z ≤ −3. Best: fam 38 7.2e12 (z 14), fam 106 6.7e12 (z 12), fam 126 5.7e12
  (z 8).
- **Between ETS:** WA_BE_broad (918) 2.4e12 per LFE, z 4.2; WA_BE (220) 3.1e12, z 2.0. On the headline grid both
  sit on the negative lobe (z −4.4, −2.8).
- **vs dv/v:** no relation (ρ 0.18).

Caveats:
- The column grid was not the headline, so these are column numbers, not a pass.
- The ±k-day calibration (W3) was computed on the headline grid only.
- Neighbours 2–10 s apart remain after the ±2 s dedupe (0.81 per LFE, kernel-weighted +0.20). The WA per-LFE value
  may therefore be inflated by up to ~20 % relative to SVI (SVI −0.05).
- WA per LFE ≈ 0.5× SVI (2.7e12 vs 5.4e12; zone 1 4.1e12). The region-to-region ratio carries the GF caveat (F5)
  and the different catalogues (P_own_pnet, family definitions).

**Decision needed (not taken here):** adopt the WA origin as `time − 5.5 s` (equivalently the c6 grid). The prior came
from the template arrivals before the run, and the grand-stack scan confirmed it. Merlin should sign off first,
because the choice was confirmed by the data it is applied to. The stacks already exist (`data/stacks/wa2_*.pkl`;
every c6 column has its own controls), so adopting it needs only a summary pass, not a re-run.

| rule | WA headline (±8 s about time + 0.5) | WA column c6 (not a verdict) |
|---|---|---|
| W1 lag mode ±4 s, ≤ 10 % at ±8 | **FAIL** (mode +8, 33 % at ±8; grand-stack peak −5 / −7 s) | mode −5 (by construction near grid centre −6) |
| W2 nb ±2 s < 0.05 | pass (0.00; ±10 s 0.81) | same |
| W3 ±k-day SD 0.85–1.15 | pass (1.14) | not computed |
| W4 zones 2–3 positive in most big episodes | pass (20/24, 18/24) | pass (22/24, 22/24) |
| W5 odd/even ρ > 0.4 | pass (0.45) | 0.57 |
| A1 detections | 19/54, 4 at z ≤ −3 | 21/54, 0 at z ≤ −3 |
| A2 per LFE | 1.8e12, χ²/dof 6.9 | 2.7e12, χ²/dof 1.5 |
| A3 total vs tremor | ρ 0.63, p 0.002 | ρ 0.84, p 0.0002 |
| A5 BE_broad | −2.4e12, z −4.4 (wrong lobe) | 2.4e12, z 4.2 |

### Score (Amendment 5)
- SVI step 4: outcome **(c)**, (b) missed narrowly; 2003–05 stays low at fixed P and rests mainly on PGC Z.
- SVI step 5: **PASS** (ratio 0.80 ± 0.25).
- PGC-out: +5 % (2008–26), +9 % (all), within 30 %.
- n_nb: fixed.
- WA: **W1 FAIL → stays withdrawn as a headline**. On the pre-registered c6 column it is clean and positive in all
  zones (2.7e12 per LFE); adopting that needs a sign-off.
