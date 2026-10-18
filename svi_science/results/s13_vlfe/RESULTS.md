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

## Amendment 6 results (Merlin WA-origin review, steps 1–5) — written 2026-10-14 ~9:35 pm PT

Rules: PLAN.md Amendment 6 (8:50 pm PT) and addendum 6a (9:05 pm PT). Both were written before any row or stack
number. Code: `s13_offsets_v3.py`, `s13_units_v3.py`, `s13_units_v3_svi.py`, `s13_v3_queue.sh` (56 stacks, 8:57–9:20 pm PT,
none failed), `s13_v3_summary.py` → `tables/v3_verdict.json`, `tables/v3_*.csv`.

### WA-NW v3 — raw figures first
- `figures/s13_v3_wa_episodes.png`, four panels:
  1. moment per LFE per major on the ±8 s grid;
  2. the same at the fixed lag +2 s;
  3. the three zones;
  4. best lag.
  Colour = z vs scrambled control; no error bars.
- `figures/s13_v3_wa_grandscan.png`: v2 (dotted) vs v3 (solid) grand-stack lag scans.

Plainly: once each family's catalogue time constant is removed, the WA 20–50 s signal sits where SVI's does.
- **The lag is +2…+3 s after the origin**, the same as SVI.
- **It is positive in 50 of 51 majors.**
- **The size is 0.77× SVI.**
- The v3 grand stack peaks at +3 s pooled (3.9e12) and at +2 s in zone 1 (5.5e12). In v2 the peaks were 2.6e12 and
  4.0e12, at −5 / −7 s. The opposite-sign lobe now sits at +17…+18 s, outside the ±8 s grid.

**Status by the pre-set rules: "WA present, size provisional".**
- W1 passes, but two caps fixed before the stacks keep it from being a headline:
  - step 1: the leave-one-out anchor SD is 3.16 s (rule < 3 s);
  - step 3: the dedupe removed 44 % of rows (rule 25–35 %).
- Both caps concern the catalogue (anchor and duplicates), not the stacks.

| rule (Amendment 6) | v3 |
|---|---|
| **W1** pooled lag mode in +1…+4 s, ≤ 10 % at ±8 | **pass**: mode +3 s, median +2 s, 2 of 51 (4 %) at ±8. Lag counts: +1: 9, +2: 11, +3: 14, +4: 8; 9 outside +1…+4. Under the CVM v1.7 anchor every lag moves +2.0 s (mode +5) |
| W2 n_nb2 < 0.05 | pass (0.00; n_nb10 0.14; kernel-weighted −0.01, v2 +0.20) |
| W3 ±k-day z SD 0.85–1.15 | pass: 1.01 (510 shifted stacks, none at abs(z) ≥ 3); no σ scaling |
| W4 zones 2 and 3 positive in > half of the 20 big episodes | pass: zone 1 18/18, zone 2 17/20, zone 3 19/20 |
| Size: pooled weighted mean ≥ 3.5e12 | pass: **4.1e12 N m per LFE** (χ²/dof 0.96, p 0.56), **0.77 ± 0.04 × SVI** (5.39e12, 2008–26) |
| W5 odd/even ρ > 0.4 | pass: 0.61 (n 48) |
| A1 detections | 25/51 at z ≥ 3 (15/20 big), none at z ≤ −3, 50/51 positive |
| A3 episode total vs tremor / LFE count / WA GNSS | ρ 0.89 (p 0.0002, n 23) / 0.90 / 0.69 (p 0.004, n 15). Per LFE vs the same: ρ 0.12 / −0.02 / 0.50 (p 0.06) |
| A4 families | 48 families; 14 at z ≥ 3, none at z ≤ −3 |
| A5 between ETS | WA_BE (195): 3.9e12, z 3.0, lag +1 (0.95× WA ETS). WA_BE_broad (629): 2.4e12, z 3.4, lag +1 (0.58×) |
| A6 vs dv/v | ρ 0.33 total (p 0.11), 0.01 per LFE: no relation |

- **Fixed lag +2 s column (no scan):**
  - 31/51 at z ≥ 3 (18/20 big), none at z ≤ −3, 48/51 positive.
  - Weighted mean 3.9e12 (0.72 ± 0.04 × SVI), χ²/dof 2.8. The scatter is above the noise: a single fixed lag fits
    some episodes less well than their own best lag.
  - Zones: 5.2 / 3.8 / 3.1e12. Odd/even ρ 0.66.
- **Zones on the ±8 s grid** (all positive, none at z ≤ −3): zone 1 (SE) 5.4e12, zone 2 3.97e12, zone 3 (NW) 3.4e12.
- **Weak episodes:** 2017-03 (lag −4, z 1.7) and 2020-10 (lag −5, z 1.4) are the big episodes outside the W1 window.
  One small episode in 2019 is negative.
- **Caveats:**
  - The absolute moment carries the whole-space ×2 Green's-function uncertainty, so the WA/SVI ratio is the robust
    product.
  - 65 % of WA rows now use a borrowed cluster location (step 2), which changes the Green's functions for them.

### Per-family lag prediction (the independent check)
`tables/v3_family_lag_check.csv`, `figures/s13_v3_family_lags.png`.
- The prediction comes from the catalogue alone: the v2 lag = offset_f + 2 s.
- The observed value is the v2 family grand stack (all v2 episodes, joint fit, wide grid). The v2 stacks used origin
  `time + 0.5` and never saw the offsets.
- Families with ≥ 3 v2 episodes, sorted by grand-stack z (approximate):

| fam | v2 episodes / LFEs | offset_f (s) | predicted v2 lag | observed v2 lag | v3 lag (expect ≈ +2) | note |
|---|---|---|---|---|---|---|
| 38 | 27 / 1,471 | −10.3 | −8.3 | −8 | +2 | cluster best member |
| 106 | 26 / 1,134 | −8.7 | −6.7 | −7 | +2 | GF → fam38 |
| 75 | 27 / 1,289 | −4.6 | −2.6 | −14 | +2 | v2 used its own location 66 km away; v3 uses fam38's |
| 126 | 19 / 789 | −5.6 | −3.6 | −5 | +2 | GF → fam38 |
| 37 | 6 / 186 | −2.9 | −0.9 | 0 | +3 | anchor family |
| 52 | 11 / 332 | −4.6 | −2.6 | −2 | +3 | GF → fam73 |
| 35 | 8 / 249 | −3.9 | −1.9 | −1 | +2 | |
| 124 | 9 / 295 | −9.8 | −7.8 | −8 | 0 | |
| 20 | 15 / 528 | −9.2 | −7.2 | −7 | +2 | |
| 73 | 4 / 102 | −9.7 | −7.7 | −8 | +2 | |
| 58 | 11 / 339 | −8.5 (default) | −6.5 | −11 | −3 | unanchored component; its own clear pick says −14.4 |
| 46 | 9 / 293 | −6.2 | −4.2 | −22 | +29 | weak (z ≈ 3 / 1.5); cluster span 136 km |

- Among the families with z ≥ 3 in v2 (6), 5/6 are predicted within 2 s (median residual −0.04 s, MAD 0.8 s). The
  miss is fam75, explained by its location.
- In v3 the same 6 families sit at median lag +2 s (5/6 within ±2 s of +2).
- The strong families fam38/106/126 are not anchors (they are flagged as duplicates). So their match tests the
  anchor constant independently.
- Where a family has a network offset, the prediction holds to ~1 s even at z 2–3 (35, 124, 20, 73). The misses are
  a default-offset family (58) and a weak family in a 136 km cluster (46).

### Catalogue facts (stated plainly)
- **Family time constants.** The WA catalogue `time` is not a common origin. Each family has its own constant,
  −14.5 … −0.3 s relative to `time + 0.5` (SD 3.5 s, IQR −10.4 / −8.5 / −4.6 s, network families).
  - Duplicate pairs sit at exactly fixed offsets (302 of 325 pairs with IQR ≤ 0.1 s; network residual rms 0.06 s).
  - The offsets are stable in time: fam38–fam75 is 5.65 s in 2007–11, 2012–17 and 2018–24; the 10 biggest pairs
    drift ≤ 0.05 s.
  - 89 % of rows sit on the anchored network; 36 families take the default −8.5 s.
- **Duplicates.**
  - 98 of the 132 NW families sit in 30 duplicate clusters, which hold 93 % of the rows. The catalogue therefore has
    about **64 distinct sources**.
  - On corrected times, 44 % of the selected rows are cross-family duplicates within ±2 s (v2 removed 23 %).
- **Mislocated families.**
  - 27 of the 30 clusters span more than 15 km (up to 136 km) for the same events. So at least one member of each is
    mislocated by tens of km.
  - 66 families (65 % of rows) now take the best member's location; the median move is 23 km, and 44 move > 15 km.
  - fam75 (catalogue depth 16.6 km, 66 km from fam38) and fam121 → fam38's location; fam77 (depth 0.01 km) →
    fam58's location.
  - 32 families are shallower than 20 km or have no depth.
  - The single anchor outlier is fam49 (B quality): its clear pick disagrees with its network offset by 11.8 s.
    Without it, the leave-one-out SD would be ≈ 1.9 s (post hoc, not a rule).

### SVI step 5 — net_mad terciles within P ≥ 0.9
Raw figure: `figures/s13_v3_svi_nmt.png` (moment per LFE by net_mad third; colour = z; no error bars). Table
`tables/v3_svi_nmt.csv`. Thrust model, 200 draws.

| episode | n per third | low third | middle | high third | high/low | (high − low)/σ |
|---|---|---|---|---|---|---|
| 2004-07 | 3,326 | 1.8e12 (z 2.4) | 2.3e12 (3.0) | 3.3e12 (3.9) | 1.9 | 1.3 |
| 2005-09 | 1,945 | 2.6e12 (2.0) | 4.0e12 (3.3) | 5.3e12 (4.2) | 2.0 | 1.5 |
| 2012-08 | 454 | 4.4e12 (2.4) | 3.9e12 (2.2) | 10.1e12 (5.4) | 2.3 | 2.2 |
| 2013-09 | 196 | 3.4e12 (1.5) | 6.4e12 (3.0) | 9.4e12 (3.8) | 2.7 | 1.8 |

- **Pre-set verdict: "per-LFE VLF moment scales with detection amplitude".** The high/low ratio is ≥ 1.5 in 4 of 4
  episodes (rule: ≥ 3 of 4).
- Combined, high − low = 2.6e12 ± 0.9e12 (2.9σ, inverse-variance over the 4 episodes); each episode alone is
  1.3–2.2σ. The middle third is between low and high in 3 of 4 episodes.
- **Composition caveat:**
  - P still rises with net_mad inside P ≥ 0.9 (ρ 0.55–0.81; P medians 0.93 → 0.95 → 0.98 across the thirds).
  - net_mad is a network detection statistic (template-match height), not a calibrated amplitude.
- Wording: "LFEs detected more strongly carry about twice the 20–50 s moment per LFE". It holds in both eras (2003–05
  and 2012–13) and at fixed P ≥ 0.9.
- At the same net_mad third, 2004-07 stays lower than 2012–13 (high third 3.3e12 vs 9–10e12). The era difference of
  Amendment 5 is not removed by the amplitude split.

### Score (Amendment 6)
- **WA:**
  - Rules: W1 pass (lag mode +3 s); W2–W5 pass; size 4.1e12 = 0.77 ± 0.04 × SVI.
  - **Status "present, size provisional"** because of the two pre-set catalogue caps (anchor leave-one-out SD
    3.16 s; dedupe 44 %).
  - Independent check: catalogue-only offsets predict 5/6 strong families' VLF lags within 2 s.
- **SVI:** the per-LFE moment scales with detection amplitude (4/4, combined 2.9σ). It is still provisional, because P
  and net_mad remain correlated.

## Scan pilot (2011) — between-ETS VLFE scan, S Vancouver Island (2026-10-14 ~11:45 pm PT)
Rules: `PLAN_SCAN.md` (Merlin design `MERLIN_VLFE_SCAN.md`, amendments S1–S6 written before the numbers they affect).
Verdict file `tables/scan_verdict.json` (first pass with the duration bug: `scan_verdict_pass1.json`).

**Plain answer: NO-GO. On this network (5–7 broadband stations, 2011), a single-event 20–50 s scan cannot tell
discrete VLFEs from noise.**
- The scrambled controls fire as often as the real data at every threshold up to z = 30.
- So the scan gives no between-ETS VLFE count and no rate bound yet.
- It does show that the detector itself works:
  - injected thrusts of Mw ≥ 3.6 are found and located to ~6 km;
  - the two largest H&G 2011 templates come out at z ≈ 18 with Mw 3.4;
  - the ETS stacks locate within 13–17 km of the LFE centroid.
- What fails is the background: one-station glitches and network-wide long-period waves outnumber real events.

Raw figures: `figures/scan_z_2011.png` (max z over the grid for 2011 and for the 2011-08 ETS; masks shaded; ETS
marked; H&G times as red ticks; detections coloured by network CC), `figures/scan_separation.png` (separation
curves), `figures/scan_injection.png`, `figures/scan_map_2011.png`, `figures/scan_loc6e.png`.

### Go/no-go items (PLAN_SCAN §7)
| item | rule | result | pass |
|---|---|---|---|
| controls at the headline threshold | (a) ≤ 1/yr and (b) ≤ 2/yr | no listed threshold meets it; at z ≥ 7 (highest listed): (a) 714/yr, (b) 462/yr between ETS | NO |
| 6a ETS re-find | ETS detection density ≥ 10× rest | 2.2 vs 1.35 per live day = 1.6× | NO |
| ETS sign check | positives ≥ 10× negatives in ETS | 33 vs 21 | NO |
| 6b mechanical consistency | mean M̂0 at LFE times = 5.3e12 ± 20 % | 3.84e12 (0.73×); 5.18e12 with gains set to 1 | NO (no bug: same within ~1.5 σ) |
| 6c H&G re-found | ≥ 3 of 4 (z ≥ threshold, ≤ 40 km) | 2 of 4 | NO |
| 6d completeness (summer, standard GF) | ≥ 80 % at Mw 3.6, ≥ 50 % at 3.3 | 82 % at 3.6, 24 % at 3.3 | NO (3.3 fails) |
| 6d location | median ≤ 30 km at Mw 3.6 | 5.9 km (95 %: 29 km) | yes |
| 6e ETS-stack location | ≤ 20 km | 17 km (2012-08), 13 km (2018-06) | yes |

### What the data show
- **The background is not Gaussian, and the controls match the real data.** The rates below use the headline mask and
  the headline filter (a cluster = the largest |z| within ±150 s, either sign).
  - z ≥ 4: real positive 4,362/yr, negative 4,408/yr, control (a) 2,703/yr, control (b) 4,346/yr between ETS.
  - z ≥ 7: real positive 515/yr, negative 492/yr, (a) 714/yr, (b) 462/yr.
  - Extended (diagnostic, `tables/scan_separation_tail.csv`): at z ≥ 15, real 85/yr, (a) 226/yr, (b) 78/yr. At z ≥ 30
    the three are 34, 69 and 27 per year.
  - Real between-ETS clusters sit on control (b) everywhere. In the ETS, 2 positives at z ≥ 15 compare with 2.5 per
    polarity-scrambled draw.
- **Where the false alarms come from (diagnosis, post hoc).**
  1. **Single-station transients.** At z ≥ 15 the largest single-channel share of the fit is 0.99 (median).
     - The best node sits under that station (LRIV Z, PFB Z, VGZ Z), so the "excess share" filter (S5) does not catch
       it.
     - Requiring z ≥ 7 also WITHOUT the dominant station (`scripts/scan_loso.py`) cuts control (a) from 714 to
       4/yr (2/yr at z_loso ≥ 10, 0.3/yr at 15).
  2. **Network-wide long-period waves.** After that cut, real and control (b) stay equal (z_loso ≥ 7: 69 vs 66/yr
     between ETS).
     - 60 % of these clusters fall inside Merlin's literal M ≥ 4.5 teleseism mask, which S1 had relaxed.
     - Most of the rest sit at grid-edge nodes (sources outside the network).
  - Even with all three restrictions (z_loso ≥ 7, literal mask, best node inside the LFE footprint): real between-ETS
    7.6/yr vs (b) 12.7/yr (`tables/scan_posthoc_combo.csv`). There is still no excess over the polarity-scrambled
    background.
- **Detector and location work on injections** (`tables/scan_injection_summary.csv`; 1-s and 20-s sources pooled;
  headline filter at z ≥ 7).
  - Summer recovery: Mw 3.4 47 %, 3.6 82 %, 3.8 77 %, 4.0 80 %.
  - Winter recovery: Mw 3.4 43 %, 3.6 74 %, 3.8 83 %, 4.0 95 %.
  - Completeness Mw at 50 % / 80 %: summer 3.42 / 3.59; winter 3.45 / 3.73.
  - Location error: median 6–12 km for Mw 3.4–4.0 (95 % radius 15–35 km). The moment comes back within 10 % (median).
  - The perturbed-GF set (free surface 1, Vs 3.4 km/s, 30-s source, Mw 3.6) is recovered only 5 % of the time. A
    slower or weaker-than-modelled source of Mw 3.6 is invisible.
  - Each filter keeps ≥ 93 % of recovered Mw ≥ 3.4 injections, except Merlin's max_share ≤ 0.5, which keeps 22 %.
    PGC Z carries most of the fit weight, so real events have share ≈ 0.6. That filter was replaced by the
    pre-declared excess-share rule (keeps 99 %).
- **H&G 2011 templates** (`tables/scan_hg_refind.csv`):
  - 2 (Mw 3.7): z 18.0, our Mw 3.42, 13 km away. Re-found, but its cluster is HF-loud: the 2–8 Hz envelope jumps
    15–33× the day median 9–24 s after the origin at all five stations, with no catalogued earthquake.
  - 4 (Mw 3.5): z 18.8, Mw 3.44, 13 km. Re-found. Its opposite-sign side lobe 15 s later is slightly larger (−20.8),
    so the joint cluster is labelled negative.
  - 3 (Mw 3.3): z 4.4 at 44 km. Missed.
  - 5 (Mw 3.4): z 6.6 (< 7), 15 km. Missed at the threshold.
  - Our Mw is 0.1–0.3 below theirs for all four.
- **ETS share of the LFE-locked VLF moment (6a column).** Σ M̂0 of the 33 ETS-time headline positives = 7.2e15 N m,
  i.e. 0.93 × the S13 time-locked total (7.7e15). Most of those 33 are probably false alarms, given the controls, so
  this is not yet a measurement.
- **Sign ambiguity in this band.** The band-passed step template is antisymmetric in time, so every event gives an
  opposite-sign match 13–17 s away at 60–110 % of the main |z|. Without known origin times (unlike the S13 stacks), the
  sign test is weak.

### Caveats
- One year, one region, the 2011 network: SNB was too noisy Jan–Oct, LZB never usable, YOUB has only 58 days, and PGC
  was excluded for 3 months by the noise rule.
- Gains come from S13 stacks, which share the whole-space Green's function.
- The post-hoc diagnostics are not tests. They show where an amended design would have to change.
- The real-run `m0_20s` column uses the old τ = 20 s (S6) and enters no rule.

### Adaptations for a re-design (for Merlin; nothing run)
These follow NO-GO branch (ii) plus the diagnosis above:
1. Require the network statistic to stay significant without the dominant station.
2. Restore the literal M ≥ 4.5 teleseism mask, which leaves 31 % of the time, or detect surface-wave trains from the
   data.
3. Restrict the grid to the LFE footprint plus a buffer.
4. Raise the low corner to 25 s.
5. Add stations (C8/XU eras, UW 2021+).
If (b) still matches the real data, use the fallback: the stack-at-candidate-times path (S13 A5 extended to
P 0.5–0.75 candidates and PNSN between-ETS tremor times), which answers the complementary question at coarse
resolution.

## Scan pilot 2 (2011) — Merlin's second-pilot design (2026-10-15 ~1:50 am PT)
Rules: `PLAN_SCAN.md` Amendment S7 (S7.1–S7.13, written before the numbers they govern; S7.12 after the cleaner failed,
S7.13 when the far-reference download was approved). Verdict file `tables/scan2_verdict.json`.

**Plain answer: still NO-GO.** The detector and location work, the decay is weak enough for far stations to act as a
referee, and the far stations do see the teleseismic body waves. But no pre-declared veto variant passes its validation,
no threshold pair gets both controls under target, and at the pair with the lowest false-alarm rate the scan finds
nothing in 2011.

Raw figures: `figures/scan2_z_2011.png` (cleaned-run z(t), mask v2, ETS, H&G and the classified times),
`figures/scan2_veto.png` (veto plane, points coloured by z), `figures/scan2_separation.png`,
`figures/scan2_injection.png`, `figures/scan2_decay_stacks.png`, `figures/scan2_classify_sections.png` (record sections
of 08-18 and the 08-22 triple, near stations and references).

### Go/no-go items (S7.9)
| item | rule | result | pass |
|---|---|---|---|
| controls at the chosen pair | (a) + (b) ≤ 1/yr between ETS | no pair meets (a) ≤ 0.3 and (b) ≤ 0.7; best: z ≥ 20, z_loso ≥ 10 → (a) 1.4, (b) 3.9/yr | NO |
| completeness | ≥ 80 % at Mw 3.8, location ≤ 15 km | 58 % (48 family injections), median 6.4 km | NO (location yes) |
| H&G 2 and 4 | headline detections | z 18.0 / 18.8 < 20 at that pair (both pass every other rule) | NO |
| ETS excess | ≥ 3 in ETS, Poisson p < 0.01 | 0 in ETS (and 0 between ETS) at that pair | NO |
| veto | validated (S7.3) | no variant passes (3 of 8 body-wave clusters flagged; 7 needed) | NO |

### What the data show
- **References (S7.1, S7.13).** On disk: GNW all year (teleseism amplitude vs PGC 1.00, 27 events), DOSE (column), and the
  ETS-days set (ELW and WISH excluded by the amplitude rule). Far set fetched after approval: IU.COR, US.NEW, US.HAWA
  and the substitutes UW.OMAK (for CN.PNT) and TA.I02D (for CN.BBB) — PNT and BBB have no 1-sps channel; no 1-sps station
  exists north of the network. All five pass QC (0.77–0.93 vs PGC).
- **Decay (S7.2).** The LFE-time VLF stack is not detected at GNW in 2011-08 or 2012-08. The upper bound on the GNW/PGC
  amplitude ratio is 0.35 and 0.17 (rule ≤ 0.5: pass; whole-space prediction 0.19 / 0.13). So a local VLFE should look
  small at 100+ km, and a wave from outside should not.
- **Step cleaner (S7.5): harmful, not used.** It removed 35 % of LRIV's year (10–13 % of PGC, VGZ, PFB). Mw ≥ 3.6
  injection recovery fell from 0.82 to 0.23, because a real nearby event looks like "one loud station, others quiet" at
  that station's best node. Control (a) at z ≥ 7 did not improve (1,061/yr vs 714/yr in pilot 1; the cluster definition
  differs). Pilot 2 therefore uses the uncleaned pilot-1 runs, with z_loso handling single-station events (S7.12).
- **Masks (S7.6).** The corrected masks cover 52 % of 2011 (pilot 1: 34 %). Most of the increase is Ms = mb − 0.3
  (head2 47 %). Live time left: 0.44 yr, of which 0.41 yr between ETS and 11.1 days in the ETS.
- **Veto (S7.3).** The far stations see the teleseismic S waves, but only at 1–7 noise units each. So the median
  ref_excess of the 8 known body-wave clusters is 1.0–4.1, and only 3 of 8 reach the pre-declared 3. Without the far
  set, 0 of 8 are flagged, because GNW alone has n_ref = 1 at most of these times.
  - H&G 2 and 4 are never flagged (ref_excess 0.4 / 0.5, ref_ratio 0.34 / 0.46). Family injections with Mw ≥ 3.6 are
    flagged 0 % (V1, V2) and 2 % (V3) of the time with the far set.
  - Diagnostic only, NOT a rule (post hoc): the amplitude ratio separates the two groups better than the excess does.
    The body-wave clusters sit at ref_ratio 0.67–1.11, injections mostly below 0.5, H&G 2 and 4 at 0.34 and 0.46.
    (`tables/scan2_veto_validation.csv`, `scan2_validation_events.csv`)
- **Known body-wave list (S7.4).** The pre-declared rule finds 8 of the 47 pilot-1 survivors, not Merlin's 11. All 8 are
  M 5.1–5.5 teleseisms, mostly 100–250 s after S (`tables/scan2_bw_known_list.csv`).
- **Threshold grid (S7.8, `tables/scan2_threshold_grid.csv`; the veto is inactive because none validated).**
  - z ≥ 10, z_loso ≥ 10: (a) 1.4/yr, (b) 7.8/yr between ETS. Real data: 1 between-ETS cluster (2.5/yr) and 3 in the
    ETS window (11 live days).
  - z ≥ 20, z_loso ≥ 10: (a) 1.4/yr, (b) 3.9/yr. Real data: 0 clusters.
  - Mw-3.8 recovery is 58–62 % at z_loso ≥ 10. z_loso costs completeness because PGC carries most of the weight.
  - At z ≥ 10, z_loso ≥ 10 the five real headline clusters (raw observation, not a result) are:
    - H&G 2 (Mw 3.42) and H&G 4 (3.44);
    - 2011-08-28 01:08:29 UTC (Aug 27 18:08 PDT): Mw 3.6, in the ETS, HF-loud, ref_ratio 0.36 — reference-consistent
      with a local source; a new candidate;
    - 2011-09-04 15:15 UTC: inside ETS + 2 d;
    - 2011-10-22 16:11 UTC (09:11 PDT), between ETS: ref_excess 4.4 and ratio 1.02, which looks like a wave from
      outside.
- **Polarity-scrambled control (b).** Between ETS, (b) fires more often than the real data at high z (e.g. 3.9/yr vs 0
  at z ≥ 20, z_loso ≥ 10). With ~12 live components, random-sign templates still match parts of real long-period
  wavefields. As a false-alarm referee for this network, (b) is conservative.

### The 2011-08-18 21:19:46 UTC cluster (14:19:46 PDT) and the 08-22 triple (S7.10)
- **08-18: rule class "unresolved"** (`tables/scan2_classify.csv`).
  - Checks: not single-station (z_loso 30.2 without LRIV); no catalogue event in any mask tier, including ComCat M ≥ 3
    within 2,500 km; ref_excess 3.9 ≥ 3, so not "local"; but no validated veto, so not "travelling" by rule.
  - The evidence points to a wave from outside, not a northern VLFE:
    - the far eastern stations see it at 7–12 noise units (NEW 7.5, HAWA 12.5, OMAK 7.2; DOSE 14);
    - the plane-wave fit over 18 stations has VR 0.69, slowness 0.29 s/km (~3.4 km/s, surface-wave speed) and
      back-azimuth ~283° (arriving from the WNW, offshore);
    - its best node is on the north grid edge.
  - The extended-grid test was not run: it applies only to the local class.
- **08-22 08:21 / 08:25 / 08:28 UTC (01:21–01:28 PDT, in the ETS): rule class "local candidate"** (median ref_excess
  0.7–1.3 < 3, n_sta 5). But:
    - the references are as loud as the near network (ref_ratio 1.2–1.5);
    - NEW, HAWA and OMAK see them at 3.4–6.6 noise units;
    - all three best nodes are on grid edges outside the footprint.
  - So they are not footprint VLFEs, and they would fail the headline. The "local" label rests on the unvalidated
    median-excess veto.

### Caveats
- One year (2011), one network, 0.44 live years after masks. The ETS-days references cover only 29 days.
- No 1-sps reference exists to the north. On the 8 body-wave days only GNW and the far set were live.
- The S7.12 fp_interior rule was written after the cleaned runs gave 0.66; the fallback runs give 0.92, so it is used.
- The three veto variants were fixed before any veto number. The ratio-only separation above is a post-hoc observation.

### Next decision (coordinator/Merlin; nothing run)
- (i) Pre-register a ratio-based veto (e.g. ref_ratio ≥ 0.5 with ≥ 3 references live) and validate it on another year
  (2012, which has more year-round references).
- (ii) Or take Merlin's fallback, stacking at candidate times between ETS.
- The 2011-08-28 01:08 UTC ETS cluster is worth a look whichever route is taken.

## Scan pilot 3 (2012, ratio veto) — out-of-sample test of the 2011 ratio veto (2026-10-15 ~4:00 am PT)
Rules: `PLAN_SCAN.md` Amendment S8 (written ~2:00 am PT, before anything ran on 2012). Verdict `tables/scan3_verdict.json`.

**Plain answer: still NO-GO, but the 2012 data show where it can work.**
- The ratio veto fixed on 2011 (R1: median reference/near amplitude ratio > 0.551 with ≥ 2 references) removes the
  network-coherent background almost completely in 2012:
  - control (b) drops to 0 between ETS at every threshold with z_loso ≥ 10;
  - the real between-ETS clusters drop to 0 as well;
  - all 7 known body-wave clusters of 2012 are flagged.
- But as pre-registered it also flags 29.5 % of the injected local Mw ≥ 3.6 events (rule ≤ 15 %), so it fails validation.
- Mw-3.8 completeness is 51 % (rule ≥ 80 %), and the ETS has 1 detection (rule ≥ 3).
- The scan found **no between-ETS VLFE in 2012** at the lowest false-alarm setting. It found one ETS candidate, on
  2012-09-19 19:18:46 PDT.

Raw figures: `figures/scan3_z_2012.png` (max z over the grid, mask v2, ETS, R1-flagged clusters, headline detections
coloured by z), `figures/scan3_veto.png` (ratio vs z for real clusters coloured by z, control (a), injections, body waves,
the fixed threshold), `figures/scan3_separation.png` (rates vs z threshold, with R1 solid / without dotted, z_loso ≥ 5,
7, 10), `figures/scan3_injection.png`.

### Set-up facts (2012)
- References: all 17 pass QC in 2012 (amplitude vs PGC 0.83–1.84; WISH 0.87, excluded in 2011, now passes).
  - GNW, DOSE, WISH, SP2, NLWA, D03D and D04D are live most of the year, so a typical time has ~10–11 live references
    (2011: 6 at most times).
- Noise rule: PGC excluded October–December 2012 (Z 46–82 nm); LZB excluded all year.
- Mask v2 = 46 % of 2012 (2011: 52 %). Live time 0.52 yr; between ETS 0.41 yr; ETS window 14.7 live days.
- Filter re-validation (215 recovered injections, Mw ≥ 3.4): all criteria ≥ 0.90 except max_share (0.34), which is
  replaced by share_excess (0.995). fp_interior keeps 0.89 → used.

### Go/no-go items (S8.8)
| item | rule | 2012 result | pass |
|---|---|---|---|
| R1 validated (S8.5) | ≥ 80 % body waves flagged AND ≥ 85 % family injections Mw ≥ 3.6 kept | 7/7 flagged; 70.5 % kept (n 132) | NO |
| controls at the chosen pair | (a) ≤ 0.3, (b) ≤ 0.7, (a)+(b) ≤ 1 per yr between ETS | no pair qualifies; best: z ≥ 5, z_loso ≥ 10 → (a) 0.76, (b) 0.0, sum 0.76 | NO ((a) > 0.3) |
| completeness | ≥ 80 % at Mw 3.8, location ≤ 15 km | 51 % (n 49), median 6.4 km | NO (location yes) |
| ETS excess | ≥ 3 in ETS and Poisson p < 0.01 | 1 in ETS, 0 between ETS; p = 0.048 | NO |
| H&G 2012 | — | none listed (Table S1 not on disk) | n/a |

### What the data show
- **The veto removes nearly all background (`figures/scan3_veto.png`).** In 2012, R1 flags:
  - 99 % of the real z ≥ 7 clusters that pass every other rule outside the mask (89 of 90);
  - 94 % of the polarity-scrambled ones;
  - 0–50 % of the few control (a) survivors (2–5 per shift).
  - So most of what the 2011 scan called "real = control (b)" is long-period energy that the references see as loudly as
    the near network: waves arriving from outside.
- **Why it failed on the injections.** In 2012, R1 flags 29.5 % of the family injections with Mw ≥ 3.6. The flagged
  fraction falls with size: 31 % at Mw 3.6, 29 % at 3.8, 9 % at 4.0. The median injection ratio is 0.29–0.45 (2011:
  0.14–0.23).
  - The reason is the reference set. In 2011 the median was taken over GNW plus the five far stations (300–480 km). In
    2012 five more regional references at 100–150 km are live all year, and a local event is not small there.
  - So the same threshold means something different when the reference mix changes.
  - The column variant on the same six references as 2011 (GNW + far set; `R1_fixed6`, not a rule) passes the 2012
    validation: 7/7 body waves flagged, 94.7 % of injections kept. Its best pair is z ≥ 20, z_loso ≥ 10: (a) 0.76,
    (b) 0.0/yr, Mw-3.8 recovery 55 %, 1 ETS detection, 0 between. Without the far set (`R1_no_far`), 76 % of injections
    are flagged.
- **Separation (`scan3_separation.png`, tables `scan3_threshold_grid{,_noveto,_r1fixed6,_r1nf}.csv`).**
  - Without a veto at z ≥ 20, z_loso ≥ 10: real between ETS 7.3/yr vs (b) 7.1/yr, i.e. no separation, as in 2011.
  - With R1, at z_loso ≥ 10 and any z threshold: real between 0, (b) 0, (a) 0.76/yr.
  - (a) is the limit. It amounts to about 2 clusters summed over the five shifted controls, so it is at the
    small-number floor.
  - In the ETS window the real rate is 25/yr (1 detection in 14.7 live days).
- **Completeness is limited by z_loso and the cuts, not by detection.**
  - 98 % of the family Mw-3.8 injections reach z ≥ 5. Only 80 % reach z_loso ≥ 10, because PGC carries most of the
    weight and is excluded for three months.
  - fp_interior (0.88), isolation (0.88) and R1 (0.71) cut the rest.
  - Without any veto, the best recovery is 76 % at z ≥ 5, z_loso ≥ 3, where control (b) runs at 423/yr.
- **The one 2012 headline detection (ETS):** 2012-09-20 02:18:46 UTC = **2012-09-19 19:18:46 PDT**.
  - z 23, z_loso 13, Mw 3.45, best node 48.25 N −123.45 (inside the footprint).
  - Ratio 0.28 over 16 references, so it is reference-consistent with a local source.
  - HF-loud (2–8 Hz envelope ≥ 10× the day median at the two nearest stations), like H&G 2 in 2011.
  - Our catalogue has 16 P ≥ 0.75 LFEs within ±10 min (1.7 expected from that day's 125 LFEs) and 1 within ±60 s. A local
    VLFE inside an LFE burst; a candidate, not a confirmed event.
- **Clusters that only pass without a veto** (raw observation):
  - between ETS:
    - 2012-06-26 13:51 UTC (06:51 PDT): z 47, Mw 3.7, ratio 1.05;
    - 2012-12-05 00:50 UTC (Dec 4 16:50 PST): Mw 3.7, ratio 0.70;
    - 2012-12-31 05:16 UTC (Dec 30 21:16 PST): Mw 3.7, ratio 0.87;
  - and 2012-04-22 13:00 UTC (in the 2012-04 major ±2 d, ratio 1.15) and 2012-09-07 17:35 UTC (ETS, ratio 1.11).
  - All have the references as loud as the near network. None has a catalogued LFE within ±10 min.
  - They look like waves from outside, not between-ETS VLFEs.

### Caveats
- One year. Control (a) and the ETS count are small numbers (a few clusters).
- The body-wave validation list has 7 members.
- The R1 threshold came from 10 events in 2011. The 2012 failure is about the reference mix, not the idea. The
  fixed-6-reference version passing in 2012 is a post-hoc observation, not a test.
- Gains are from 2011 calibration stacks that include the 2012-08 ETS.

### Next decision (coordinator/Merlin; nothing run)
- (i) Pre-register the far-set ratio (GNW + ≥ 300 km references only, threshold re-derived from 2011 + 2012 body waves
  vs H&G/injections) and test it on a third year (2013 or 2014 with the 2014-12 H&G templates).
  - It still needs the completeness fixed: z_loso ≥ 10 alone caps Mw-3.8 recovery at ~80 %.
  - Lower z_loso is now affordable because R1-type vetoes kill control (b).
- (ii) The fallback (between-ETS stacks at candidate times, PLAN.md Amendment 7) is running and answers the complementary
  question.

## Between-ETS stacks at candidate times (Merlin's scan fallback; 2026-10-15 ~4:50 am PT)
Rules: PLAN.md Amendment 7 (written ~2:20 am PT, before any row, prep or stack). Summary `scripts/s13_be3_summary.py` →
`tables/be3_summary.csv`, `tables/be3_verdict.json`. Machinery as the scale-up: `s13_stack2_v2.py` (guard2, instrument
vetoes, scrambled control = random time on the same day, ±1-day stacks, ±8 s lag grid). Stations CN.PGC, CN.VGZ, CN.SNB,
CN.NLLB, UW.LRIV. Between ETS = outside every SVI major ±2 d, 2003-01-01 → 2026-09-15.

**Plain answer:**
- Between-ETS catalogued LFEs (P ≥ 0.75) and inter-ETS burst LFEs (BE_broad) carry about the ETS-sized 20–50 s
  accompaniment: 4.0e12 and 3.8e12 N m per LFE (ETS 5.4e12). Each is only z 1.8–2.6, so "about the same, cannot exclude
  half or zero" is the honest reading.
- Lower-confidence candidates (P 0.1–0.5) carry none: at 2σ, at most 0.11–0.22 of the ETS value. That is below what their
  mean P would predict if a fraction P of them were ETS-like LFEs.
- Tremor windows with no catalogued LFE carry little: 1.0e12 N m per window, z 1.8 (≤ 0.39 of the ETS value at 2σ).
- No unit reaches the pre-set "present" bar (z ≥ 3).

**Scope.** These stacks see 20–50 s moment only at times with 2–8 Hz activity: a catalogued candidate, or a tremor window
in the SVI footprint. They cannot see VLF slip with no 2–8 Hz activity. They measure the continuous LFE-locked
accompaniment; window QC drops windows holding a discrete VLFE above Mw ~3.6–3.8 (F4).

**Sampling.** The class units use a random 20 % of 7-day blocks: 224 of 1,118 eligible blocks, 1,488 between-ETS days.
Because the candidates bunch on a few days, the sample holds 8 % of the P ≥ 0.75 class (859 of 10,362) and 14–18 % of
the others. The PNSN units use every between-ETS tremor day in the footprint, 2010+ (403 days). BE_broad uses its own 688
LFEs.

Raw figures first:
- `figures/s13_be3_classes.png`: ρ = moment per candidate / ETS value for every unit, coloured by z, no error bars.
  Right axis = VLF:LFE ratio, ×12.
- `figures/s13_be3_stacks_{C075,C050,C030,C010,PNSN_noLFE,PNSN_noCand,BE_broad}.png`: per station and component; black
  stack, red thrust × M̂0, grey scrambled 2.5/97.5 %, blue ±1 day.
  - In C075 the clearest thrust-shaped pulse is on PGC Z, inside the control band. Most channels are noise-dominated, as
    the z values say.

| unit | n | mean P | M0 per candidate (N m) | z | ±1-day z | ρ = M0/5.39e12 | VLF:LFE (= 12ρ), 2σ range | ratio to BE_broad | wording (pre-set) |
|---|---|---|---|---|---|---|---|---|---|
| C075 (P ≥ 0.75) | 859 | 0.87 | 4.0e12 | 1.8 | −1.0 / 0.7 | 0.74 ± 0.41 | 8.9 (−0.9 … 18.7) | 1.05 | undetermined |
| C050 (0.5–0.75) | 1,706 | 0.61 | 3.1e12 | 1.7 | −1.1 / 2.5 | 0.57 ± 0.33 | 6.8 (−1.1 … 14.8) | 0.81 | "ETS-like" (rule), also consistent with 0 |
| C030 (0.3–0.5) | 3,273 | 0.38 | −0.7e12 | −0.8 | 0.3 / −0.9 | −0.14 ± 0.18 | −1.7 (−5.9 … 2.6) | −0.20 | below ETS |
| C010 (0.1–0.3) | 11,046 | 0.18 | −0.4e12 | −0.7 | −0.7 / −0.4 | −0.07 ± 0.09 | −0.8 (−3.0 … 1.3) | −0.10 | below ETS |
| PNSN_noLFE (tremor, no P ≥ 0.75 within ±5 min) | 7,333 | — | 1.0e12 | 1.8 | −0.9 / 0.9 | 0.19 ± 0.10 | 2.2 (−0.2 … 4.6) | 0.26 | below ETS |
| PNSN_noCand (no P ≥ 0.1 within ±5 min) | 4,758 | — | 0.9e12 | 1.4 | −1.2 / 1.1 | 0.17 ± 0.12 | 2.1 (−0.7 … 4.9) | 0.25 | below ETS |
| BE_broad (anchor, 5 stations) | 688 | 0.86 | 3.8e12 | 2.6 | −1.0 / 0.7 | 0.71 ± 0.28 | 8.5 (1.9 … 15.1) | 1 | ETS-like |

(± = control SD / 5.39e12. The ±1-day stacks never reach |z| 3. The ETS value is the SVI 2008–26 weighted mean,
5.39e12 N m per LFE ≈ 12 × 4.5e11.)

What the numbers show:
- **Catalogued between-ETS LFEs.**
  - C075 and BE_broad agree with each other (ratio 1.05) and sit at 0.7–0.75 of the ETS value.
  - BE_broad with these 5 stations gives 3.8e12 (z 2.6). With 10 stations (Amendment 5) it gave 4.3e12 (z 3.2), so the
    station set costs SNR, not size.
  - The pre-set "present" bar (z ≥ 3) is not reached. The control SD per LFE is ~2× the ETS-stack value: winter noise,
    and only 5 stations.
- **The lower-confidence classes fall off faster than their P.** Read as a mixture (a fraction f of real LFEs with the ETS
  accompaniment), ρ = f:
  - C050: f 0.57 ± 0.33 vs mean P 0.61 — consistent.
  - C030: f ≤ 0.22 (2σ) vs mean P 0.38, 2.9σ below.
  - C010: f ≤ 0.11 vs mean P 0.18, 2.7σ below.
  - Either P overstates the real fraction at low P, or the weak real LFEs carry less VLF moment. The second matches the
    Amendment 5 finding that per-LFE moment rises with P.
- **Tremor days vs quiet days (column splits, 2010+; no rule).** Candidates on PNSN tremor days carry more than those on
  quiet days in every class:
  - C010: 2.4e12 (z 1.7) vs 0.3e12;
  - C030: 1.9e12 vs −1.2e12;
  - C050: 4.1e12 (z 2.3) vs 2.5e12.
  - C075 on tremor days: 4.2e12, z 2.7.
- **Tremor with no catalogued LFE.** The PNSN windows give 1.0e12 N m per window (z 1.8). That is 0.26 of the BE_broad
  per-LFE value and ≤ 0.39 of the ETS value at 2σ. A tremor window with no catalogued LFE carries at most a few tenths of
  one LFE's VLF accompaniment, not the large VLF moment the H&G between-ETS activity would imply.
  - Caveat: the origin is the network 2–8 Hz envelope peak, which falls a median 102 s after the PNSN epoch time (5–95 %:
    −127 … +283 s). A burst lasting minutes has no single origin, so some VLF moment can be smeared outside the ±8 s lag
    scan. The best lag is +2 s (noLFE) and +6 s (noCand).
- **Eras (column):** C075 2003–08 −1.2e12 (z −0.2, n 555), 2009–13 3.7e12, 2014–26 4.9e12 (z 2.5). The 2003–08 between-ETS
  catalogue (dense, few stations) carries the least, as in the ETS (Amendment 5, step 4). C030 and C010 are ≤ 0 in every
  era.

Relation to the Hutchison & Ghosh reproduction (`hg_repro/RESULTS.md`):
- Their between-ETS matched-filter detections do not rise above their own noise-template background, and they carry no
  more of our LFEs than random times.
- This fallback shows the converse at coarse resolution. Where between-ETS 2–8 Hz activity exists, catalogued LFEs carry
  the ETS-like VLF accompaniment, while tremor without catalogued LFEs and low-P candidates carry little.
- Neither route supports a large population of VLF events between ETS in SVI.

Caveats:
- All z < 3: these are sizes with uncertainties, not detections.
- 20 % day sample: the full C075 set (10,362, 12× the sample) would reach z ≈ 1.8 × √12 ≈ 6 if size and noise scale as
  in the sample (less if the extra LFEs are day-clustered). It is the obvious next run if a firmer between-ETS number is
  wanted.
- The 4.5e11 LFE moment is assumed (Bostock mag_med 1.7), not measured.
- Absolute moments carry the whole-space ×2 uncertainty; the ratios between units do not.
- WA was not run (not quick: it needs a WA prep of the drawn days).
