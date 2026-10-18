# newvlfe RESULTS — template-free search for new VLFEs, 2013 + 2018 (written 2026-10-16 ~10:45 am PT)

Rules: `PLAN.md` (pre-registered ~10:05 am PT). Log: `NOTES.md`. "Scrambled control" = the same scan on data scrambled
so no real event can line up. There are two kinds: polarity draws (polB, 4 per year) and station-day shifts of all stations
(shift, 3 per year). "Known" = best node ≤ 15 km from a template location (H&G T1–T5, C0–C5 events, 2018 relay events).
"New" = more than 15 km. ETS = SVI episodes (start − 2 d to end + 3 d). Live time: real 0.418 yr in ETS and 1.004 yr
between; controls 2.93 yr in ETS and 7.03 yr between.

Raw figures first:
- `figures/nv_z_vs_ratio.png`: every gated candidate, z vs far/near ratio; real points coloured by vr_thrust, controls grey,
  B2 rings.
- `figures/nv_map.png`: B0/B2 locations, real vs controls, with 15-km circles around the templates.
- `figures/nv_rank_ets.png`: rank curves, real vs expected from each control type, per year.
- `figures/nv_separation.png`: rates vs z threshold.
- `figures/nv_wave_*.png` / `nv_ctrlwave_*.png`: raw 0.02–0.05 Hz Z record sections (near red, mid blue, far black).

## What the data show
- **Near the known patches, the template-free scan finds the ETS VLFEs on its own.**
  - At tier B1 (catalogue battery + ≥ 3 far references + no single-station dominance + vr_thrust ≥ 0.12), z ≥ 7 gives
    90 real events in the ETS. That is 215 per live year, against 44 per year in the pooled controls (4.9×) and 13 per year
    in the polarity controls (16×).
  - Both years show it. 2013 (clean network): 32 real vs 4.5 expected from polB at z ≥ 5. 2018: 110 vs 16.
  - Between ETS the known patches show no excess: 15.9/yr real vs 32/yr pooled controls (7.0/yr polB).
- **Away from the patches (> 15 km), the excess in the ETS is small and not clean.**
  - B1, z ≥ 5: 45 real vs 32.0 expected from the pooled controls (1.41×, p = 0.017).
  - **Implied new ETS events: 45 − 32.0 ≈ 13, 95 % range about −1 to +27** (Poisson on the real count plus the control
    count). The estimate is +17 with the polB controls alone and +7 with the shift controls alone.
  - **All of it is at the lowest z.** By z band: z 5–7 gives +14.9 (32 real vs 17.1 expected); z 7–10 gives 0.0 (8 vs 8.0);
    z ≥ 10 gives −1.9 (5 vs 6.9).
  - **All of it is in 2018.** 2018: 36 real vs 23.3 expected (+12.7). 2013: 9 vs 8.7 (+0.3), and 2013 real tracks the polB
    expectation at every z (rank curve).
  - **RETRACTED (W2b section below): LRIV was already off in every 2018 scan run, so it cannot explain the 2018 excess.**
    Original text: **2018 is the year UW.LRIV was broken through the ETS.** It shows 3.5–9.8 µm square-wave glitches in the fit window
    (`nv_wave_2018-06-29_041643.png`, `2018-06-28_175116`, `2018-07-09_044903`). LRIV is broken at 92 % of the 2018 new
    real rows and at 79–85 % of the 2018 control rows. With LRIV working, the new-ETS count is 12 real vs 12.8 expected.
  - The HF class does not discriminate: 43 of 44 new ETS B2 real rows are HF-quiet, and 215 of 220 controls.
- **Between ETS: nothing new.** New B1 at z ≥ 5 gives 26.9/yr real vs 46.2/yr pooled controls (23.9/yr polB). At every
  tier and z threshold the real rate is at or below the control rate. The same holds for "other" (NW-only) episodes.

## The new ETS candidates (`tables/nv_new_candidates.csv`, 44 B2 rows in ETS; controls `nv_new_controls.csv`)
- 19 of the 44 form one group: 2018-06-25 – 07-15 PDT, at 48.5–48.9 N, 124.26–124.53 W. That is 20–45 km west of relay
  cluster RC4, near PFB and SVI families 89/98/108 and Royer_NW_WA:2226.
  - 11 have SVI/WA-NW LFEs within ±10 min.
  - 6 have patch-specific LFE coincidence at p ≤ 0.05 (local families ≤ 20 km against random same-day times), e.g.
    2018-06-29 04:16:43 and 04:34:52 PDT (7 patch LFEs each, p 0.03–0.04), 2018-06-30 07:31:21 PDT (5, p 0.03) and
    2018-07-03 09:48:02 PDT (2, p 0.01).
  - Across all 44 new ETS rows: 43 % have LFEs within ±10 min (controls 32 %), and 14 % have patch p < 0.05 (controls 3 %).
    So about 6 rows against about 1.4 expected.
  - In waveforms none shows a pulse above the background. PFB (8–20 km) has 65–90 nm wiggles like those before and after,
    and LRIV is clipped. 8 of the 44 were also seen by a production template (MF 5.6–7.3); 4 of those 8 are production tier
    events.
- The top-z new row, 2018-11-14 01:23:59 PST (z 26, 48.07 N 122.77 W, at the east grid edge, no LFEs), is a coherent
  ~140–210 nm pulse across the near and regional network that reaches WISH first. A polarity-scrambled control fires
  2 s later at the same time. So it does not depend on the thrust template: it is a real wave of unknown origin, not a
  footprint VLFE.
- 2013-09-28 08:04:50 PDT (z 10.9) is a 20-s wavetrain of similar amplitude on near and mid stations, with the far
  stations at 30–80 nm. It is background noise.

## Why the controls run at ~44/yr here and at ≤ 1/yr in the catalogue
1. **The selector is different.** The catalogue's controls must first pass a matched filter on an empirical 90-s waveform at
   ≥ 5.5 MAD with two templates agreeing. That leaves a few hundred control detections a year, and the battery passes
   0.1–0.3 % of them. The template-free scan passes every 20–50 s energy burst that fits a synthetic thrust: 54k real and
   487k control positive candidates at z ≥ 5. The same kind of battery passes 0.2 % (polB) to 0.4 % (shift) of those, which
   still leaves tens per year. It is the volume of the synthetic-GF selector, not a weaker battery.
2. **The shift control is too generous for this battery.** Shifting stations to different days destroys the near-to-far
   coherence of real long-period waves, which is what the far/near test catches. The shift rows that survive the
   single-station test pass far/near 87 % of the time; real and polB rows pass 13 % and 11 %. Between ETS the shift
   controls exceed the real data (known 66 vs 16/yr). So the shift does not let events line up; it makes a different,
   incoherent noise population that the battery cannot reject.
3. **The polarity control is the closer match, but in the ETS it partly fits real signal.** Between ETS polB ≈ real (new
   9.7 vs 9.0/yr at B1 z ≥ 7). In the ETS polB rises to 25/yr for new and 13/yr for known, because scrambled-sign grid
   templates still fit real ETS long-period energy. In 2018 the LRIV glitches add to this.

## Validation of the added tests (PLAN §2; 93 production tier events of 2013 + 2018)
n_far ≥ 3 keeps 1.00, share_excess ≤ 0.25 keeps 0.98, vr_thrust ≥ 0.12 keeps 0.85, and the dense-ring plane-wave test
keeps 1.00, so all are used. z_loso ≥ 5 keeps 0.58 and is a column only. B1 and B2 differ by ≤ 4 events everywhere.

## Scope
- Two years (2013, 2018), scan-v2 network, 20–50 s.
- Events must fit the synthetic thrust grid at z ≥ 5 (detection Mw ~3.4–3.6 in the footprint).
- The "known" class uses 15 km around 36 template locations, and scan location errors are 6–12 km (median), so some
  "new" rows at 15–30 km may be mislocated patch events.
- Point 5 (stacking new events into templates) was not run, by instruction.

## W2b re-run of 2018 and correction (written 2026-10-16 ~11:00 am PT; `tables/nv_w2b_zbins.csv`, `nv_w2b_west.csv`)
Run: scan v2 for 2018 with the catalogue's W2b station-component-month exclusions (`catalogue/tables/qc_exclusions.csv`,
2018 = 33 UW.LRIV N/E/Z months; read only). One GPU pass with 4 polarity draws (variant 0 = real), then the same battery.
Marker `W2B2018_DONE` (10:54 am PT).

- **W2b changes nothing.** The scan's own monthly noise rule had already given LRIV zero weight in every 2018 month (all
  three components, in the real run and in every control). The candidate list is identical (72,678 real), and so is every
  count below.
- **Correction:** LRIV appears in the 2018 waveform figures, but it never entered the 2018 detection or the battery. The
  earlier "LRIV-broken" explanation is retracted. The "LRIV ok: 12 vs 12.8" comparison was really 2013 vs 2018.
- **2018 ETS, new (> 15 km), tier B1, real vs expected from the polarity controls** (live 0.255 yr real, 1.018 yr control):

| z | real | expected | excess |
|---|---|---|---|
| 5–7 | 27 | 13.2 | +13.8 |
| 7–10 | 5 | 3.2 | +1.8 |
| ≥ 10 | 4 | 1.5 | +2.5 |
| all ≥ 5 | 36 | 18.0 | **+18** |

  For comparison: known patches 110 vs 16.2 (+94). With the shift controls the 2018 new expectation is 30.3 (excess +6),
  but that control is too generous for this battery (see above).
- **The western group** (B2, ETS, 48.45–48.95 N, west of 124.2 W):
  - 17 real rows vs 6 in the polarity controls, i.e. 1.5 expected in the real live time.
  - 11 of the 17 have SVI/WA-NW LFEs within ±10 min.
  - **6 have local-family LFE coincidence at p ≤ 0.05**, against 2 in the controls (0.5 expected).
  - The 6 are 2018-06-26 21:39:52, 06-29 04:16:43, 06-29 04:34:52, 06-29 18:34:10, 06-30 07:31:21 and 07-03 09:48:02 PDT,
    at 48.52–48.79 N, 124.26–124.39 W, z 5.9–7.5, with 2–7 local LFEs each (p 0.010–0.039).
  - None is a production tier event. Only one (06-30 07:31) has an MF detection within ±30 s (RC1, 5.6 MAD).
  - They lie 20–30 km west of relay cluster RC4. Scan location errors are 6–12 km median, up to ~30 km, and the network sits to
    the east, so these may be RC4-type events placed too far west rather than a new patch.
  - At z 6–7.5 no single station shows a pulse above its background in the raw waveforms, so a stack is the test.

## Point 5, western group: stacks vs RC4 (written 2026-10-16 ~11:10 am PT; `scripts/nv_stackw.py`, `nv_stackw_coh.py`, `nv_stackw_ref.py`)
Stacks: W6 (the 6 LFE-coincident candidates) and W17 (all 17 western events), built with the catalogue relay recipe (90-s
0.02–0.05 Hz velocity from origin − 9 s, unit per channel, ±5 s alignment, 2 passes) on the 2018 production channels.
Raw figures: `figures/nv_stackw_section.png` (W6, W17, RC4), `nv_stackw_amp.png`, `nv_stackw_loc.png`.
- **Waveform match with the bank** (channel-mean CC at the best lag within ±20 s; `tables/nv_stackw_cc.csv`):
  - The best match for W6 is RC4 at 0.26 (lag +9 s). For W17 it is C2 at 0.23, then RC4 at 0.22 (lag +10 s).
  - Same-source references: one relay member against its own leave-one-out stack gives 0.25–0.49 (median 0.27); a relay
    leave-one-out stack against its primary gives 0.72–0.95.
  - Different relay clusters against each other give 0.06–0.41 (median 0.22).
  - So the stacks match RC4 only at the level of two different clusters, at a lag of 9–10 s.
- **Station amplitude pattern** (20–50 s Z, median over members; `tables/nv_stackw_amp.csv`): the W pattern correlates
  0.67 with RC4 and 0.45 with RC0, while RC0, RC2 and RC4 correlate 0.88–0.95 with each other. W is relatively larger at
  PFB and NLLB. Signal-to-pre-event ratios are 0.9–2.5 for W and 0.7–1.3 for RC4: single stations sit at noise level for
  both.
- **Stacked scan location** (Σ num / √Σ den over members; `tables/nv_stackw_loc.csv`): W6 48.52 N 124.26 W, W17 48.70 N
  124.26 W, RC4 48.52 N 123.99 W (it reproduces its own events), RC0 48.43 N 123.85 W. W6 is about 20 km due west of RC4.
- **Internal coherence** (member vs the stack of the others, best of ±5 s; `tables/nv_stackw_coh.csv`, `nv_stackw_ref.csv`):
  - W17 members: median 0.18 (0.10–0.38).
  - Selection-matched noise, i.e. the 208 between-ETS scan candidates at the same western nodes against the W17 stack:
    median 0.11, 95th percentile 0.22. 5 of 17 members exceed that (0.85 expected).
  - W6 members: 3 of 6 exceed the W6 reference 95th percentile (0.19); 0.3 expected.
  - Relay members RC0/RC1 for comparison: 0.22–0.49, all above every random time.
- Reading: the western group is not RC4 shifted west. Its waveform match to RC4 is at the different-cluster level, its
  amplitude pattern differs, and its stack locates 20 km west. But it is only weakly self-similar: about a third of the
  members carry a shared waveform beyond what the scan selection gives. The relay acceptance test below decides whether
  that is a stable source.
- Relay acceptance on 2018 (report only; W2b QC on, so LRIV is excluded): 4 clusters by complete linkage 15 km (WC0 7
  members, WC1 5 — it holds 4 of the 6 LFE-coincident events — WC2 3, WC3 2) plus leave-one-out stacks. Running:
  `scripts/nv_wbank_queue.sh`, marker `STACKW_DONE`.

## Relay acceptance of the western clusters, 2018 (written 2026-10-16 ~11:20 am PT; report only)
Marker `STACKW_DONE` (11:15 am PT). Tables: `tables/nv_wbank_verdict.json`, `nv_wbank_loo.csv`, `nv_wbank_tier_events.csv`,
`nv_wbank_bank.csv`. Bank: WC0 (7 events), WC1 (5), WC2 (3), WC3 (2), each with leave-one-out stacks; catalogue relay rule;
W2b on; z* = 6.33 MAD (prod2018 threshold at 0.5 noise detections per template-day).
- **Leave-one-out recovery: 2 of 17 (12 %; rule ≥ 85 %).**
  - Re-found: W4 (2018-06-26 21:39:52 PDT, 6.6 MAD) and W14 (2018-07-09 04:49:03 PDT, the z-21 scan event, 8.2 MAD).
  - The other 15 reach only 3.5–5.9 MAD with templates that do not contain them. Of the 6 LFE-coincident events, 1 is re-found.
  - For comparison, the 2018 relay bank re-found 15 of 17 (88 %).
- **Scrambled controls: 0.11 per year between ETS** (1 event in 16 catalogues), under the 1-per-year rule.
- **Detections beyond the original 17** (tier = ≥ 5.5 MAD, a second template agreeing, battery′): 8 events.
  - 5 sit at the known patches and are also detected by production templates: 2018-07-08 23:56 PDT (C0/C1/RC/T agree),
    07-09 00:49, 00:57 and 01:04 PDT (RC0/RC1/RC2/RC4 agree, around the RC4 members of 00:04 and 00:53), and 07-09 02:47
    PDT (C2/RC5/RC6).
  - 1 is in NW Washington: 2018-05-15 15:40 PDT, 47.98 N 123.04 W, with RC1 agreeing.
  - 1 is in a NW-only episode: 2018-08-10 05:14 PDT, with vr_thrust 0.09.
  - 1 is between ETS: 2018-04-13 14:02 PDT, 48.43 N 124.26 W. Only W templates agree, and it fails the Q2 mechanism clause
    (vr_thrust 0.117 < 0.12).
  - No new western event repeats with the W templates alone.
- **Verdict: not accepted.** Taken one at a time, the western events do not find each other. What the W templates detect
  beyond their own members is mostly the known-patch VLFEs of 2018-07-08/09.
- So the 2018 western group is not a stable, repeating source at this sensitivity. It is a set of low-z scan candidates, a
  few sharing a weak common waveform, that mostly coincide with LFE activity near the patches. No 2013 or forward-year run
  (the rule needed a pass).
- **Overall answer to "only repeats, no new VLFEs":**
  - The denser v2 network lets the template-free scan re-find the known-patch ETS VLFEs on its own, but it finds no new
    source that passes the catalogue's acceptance rule.
  - Scope: 2013 and 2018; 20–50 s; events that fit the synthetic thrust grid at z ≥ 5, about Mw 3.4–3.6 in the footprint.

## NV3 — 2011, 2013, 2016, 2018, 2019 with W2b, ladder, emergence, blind eye labels (written 2026-10-16 ~7:30 pm PT)
Rules: PLAN Amendment NV3 (~4:55 pm PT) + §5 amendment (label cap, ~5:15 pm PT). Marker `NV3_DONE`. The labels step first
failed on a pandas bug (groupby-apply dropped `year`); it was fixed and re-run alone, and the earlier steps were kept.
Tables: `tables/nv3_separation.csv` (every level × z × known/new × period, per year and pooled), `nv3_ladder_counts.csv`,
`nv3_labels_with_key.csv`. Raw figures: `figures/nv3_separation.png`, `nv3_rank.png`, `nv3_map_level8.png`,
`figures/nv3_labels/<ID>.png` (238 blind figures), `nv3_wave_2011-05-01_113520.png`, `nv3_wave_2016-06-13_074128.png`.
Live time: real 1.11 yr in ETS and 2.40 yr between; polarity controls 4.45 / 9.61 yr; shift controls 3.33 / 7.21 yr.

### What the data show
- **Near the known patches, the scan alone finds the ETS VLFEs in every year.** Level 8 (all 8 checks), z ≥ 5, ETS, known
  patches: 449 real (404 per yr) against 201 per yr in the polarity controls.
  - Excess over the polarity controls by year: 2011 +54, 2013 +20, 2016 +20, 2018 +90, 2019 +41.
  - Between ETS at the known patches there is nothing beyond the controls: 213 vs 217 per yr.
- **Away from the patches (> 15 km), the ETS excess is small and comes mostly from 2018–2019.** Level 8, z ≥ 5, ETS, new:
  189 real vs 130 per yr polarity (145 shift), +44 vs polarity.
  - By year: 2011 +5.5, 2013 0.0, 2016 +6.8, 2018 +15.8, 2019 +16.3.
  - Without earthquake-wave-flagged rows: 48 vs 29 per yr (1.5×). The excess sits at z 5–7; at z ≥ 7 it is 12 vs 10/yr.
- **Between ETS, nothing new in five years.** Level 8, new, z ≥ 5: 139 vs 131 per yr (p 0.17). With flagged rows removed
  the real rate is at or below the polarity controls at every z: z ≥ 7 gives 5.0 vs 7.9 per yr, z ≥ 10 gives 0.4 vs 3.0.
- **The only pre-set cell that "separates":** ETS, new, level 8, z ≥ 20, flag-excluded, 3 real vs 5 control (0.64/yr).
  It does not hold at z ≥ 15 (3 vs 11; ratio 1.9). The 3 rows are:
  - 2011-05-01 11:35 PDT (z 41): a coherent wavetrain on every station, far stations at ~50 % of near;
  - 2018-07-09 04:49 PDT: the western-group event that failed the relay acceptance;
  - 2018-11-14 01:24 PST: a regional wave that a polarity control also fires on 2 s later.
  None is a new VLFE.
- **The earthquake-wave flag does not separate real from control.** It fires on 80 % of real and 77–100 % of control
  label rows, mostly through the "M ≥ 4.5 within 90 min" clause (a high base rate, disclosed in PLAN). The far ≥ near
  emergence clause alone is the useful half.
- Median near emergence of the labelled rows: real 2.16, control 2.15.

### Blind eye labels (119 real "new" level-8 rows, year-stratified from 600, + 119 control rows; key opened after labelling)
| | candidate | eq-wave | spike | wavefield |
|---|---|---|---|---|
| real (119) | **5** | 4 | 14 | 96 |
| control, all (119) | **1** | 2 | 46 | 70 |
| control, polarity only (59) | 1 | 1 | 11 | 46 |

- Candidates are 5 of 119 real vs 1 of 119 control (Fisher p = 0.21). Against polarity controls only: 5/119 vs 1/59
  (p = 0.67). This is no clear excess.
- The spike difference comes from the shift controls (35 of their 60 rows are spikes): shifting stations to different
  days makes single-station transients line up. The polarity controls look like the real data (mostly wavefield).
- The 5 real candidates:
  - 2011-08-13 08:50 PDT (ETS, 47.98 N 122.90 W);
  - 2011-08-17 00:05 PDT (ETS, 47.98 N 123.04 W; 7 local-family LFEs within ±10 min, p 0.04);
  - 2013-02-07 01:25 PST (between, 48.88 N 123.31 W, storm percentile 0.93, no LFEs);
  - 2018-06-25 17:42 PDT (ETS, western group W0);
  - 2019-04-05 16:50 PDT (between, 47.98 N 123.04 W, no LFEs).
- 3 of the 5 sit at 47.98 N, 122.90–123.04 W, 17–21 km from the H&G T1 location (47.83 N 123.08 W), i.e. at the
  "known/new" edge given 6–12 km location errors. The one control candidate (2016-01-15 09:25 PST, 48.07 N 123.04 W) is
  at the same spot.
- No production template detected any of the 5 (no MF ≥ 5.5 MAD within ±30 s).
- **No new source holds up** (candidates neither concentrated beyond the controls nor repeating with LFEs), so no
  stacking/relay test was run. The 47.98 N spot is the only place where a stack could be tried. Expect it to be the H&G T1
  source or noise, since a control also produces a candidate there.

### The 2016-06-13 07:41 PDT candidate (catalogue RESULTS line 743, "the best candidate")
- **It is in the pool:** real row 2016-06-13 07:41:28 PDT, z 5.5, z_loso 3.9, ratio_far 0.49, **ladder level 8**, near
  emergence 3.7 vs far 2.3.
- **It is not "new":** its best node (48.43 N 123.58 W) is 9 km from a template location, so it falls in the known class
  and outside the label set.
- It carries the eq-wave flag only through an M ≥ 4.5 event within 90 min; the far emergence is lower than near.
- In the raw section (`figures/nv3_wave_2016-06-13_074128.png`) there is a modest pulse at PGC/VGZ/SNB/LRIV in the fit window
  (44–52 nm), not clearly above the surrounding wiggles. It is consistent with a weak event at the known patches, between
  ETS.

### Scope
- Five years (2011, 2013, 2016, 2018, 2019), scan v2 with W2b, 20–50 s.
- Events must fit the synthetic thrust grid at z ≥ 5 (about Mw 3.4–3.6 in the footprint).
- "New" means more than 15 km from 36 template locations, given scan location errors of 6–12 km median.
- Eye labels are on a 119/119 sample, by one labeller, with the scan numbers visible on the figures but kind, time and
  location hidden.
