# S11 ballistic timing, round 3: WA-NW replication, ETS-phase composite, region-wide change, joint P (rules fixed first)

Written Mon 2026-10-12, 7:30 pm PT, before any round-3 number was computed. Rounds 1–2: `../PLAN.md`, `../RESULTS.md`,
`../round2/PLAN.md`, `../round2/RESULTS.md`. Round 2's one success: pooling years and families gives a decade bound
δVs/Vs < 0.014 % (2σ) on the B926 leg; per-episode floors stay 3–5 ms; P3 (3C principal direction) cc ~0.6, σ_P 13–20 ms.
User (10-12 evening): maximum effort and perseverance; coarse resolution is fine if the finding is new.
"Scrambled control" = the c2 windows (each station's data shifted to another time of the same day), as rounds 1–2.

Scripts `svi_science/scripts/s11_ball3_*.py`; outputs here (`tables/`, `figures/`, `logs/`). CPU ≤ 8 workers, RAM < 40 GB,
one region in memory at a time; nothing deleted; other jobs' files not touched. Plots: points coloured by cc, no error
bars, line graphs; times in PT.

## 0. One measurement pass per region (`s11_ball3_measure.py --region SVI|WA_NW`)
Round-1/2 machinery unchanged (S pick, P pick, ±100 ms search = ±4 samples at 40 Hz, parabolic refinement, P-weighted
stacks, 200 bootstrap draws, leave-bin-out reference = station-epoch grand stack minus the bin's rows, n ≥ 20, an epoch
is never crossed). Inputs: SVI = S10-multi tables/windows (76 families, 7 core stations); WA-NW = the S11 WA inventory
and windows (`wa_nw/tables/inv/`, 108 families, every station with windows, ≤ 80 km).

**Estimators per stack:** `S` (headline, round 1: horizontals, [tS − 0.5, tS + 1.5] s); `S3c` (column: 3C,
[tS − 0.5, tS + 2.5] s); `P3` (headline P: 3C principal direction, [tP − 0.5, min(tP + 2.5, tS − 1)] s); `P3s` (column:
same direction, short window [tP − 0.3, tP + 1.2] s); `Pz` (column). Bootstrap σ of each and of S − P3, S − P3s, S3c − P3.

**Bin schemes** (each row gets one label per scheme; stacks per (family, station, epoch, scheme, label)):
- `EP` episode unit (majors and bursts), `YR` calendar year, `Y2` two-year bins — as round 2.
- `PH` (**headline ETS phase**). Family onset t0(f, e) for each major e = the 10 % quantile of the family's data-row times in
  e (the S3 definition, applied to the repeat rows). Labels: `BTW` = rows in bursts (LFE activity between majors);
  `ARR` = major rows with t − t0 ≤ 1 day (front arrival at the family, incl. the first 10 %); `AFT` = major rows with
  t − t0 > 1 day (behind the front).
- `PH5` (column): `BPRE` bursts starting ≤ 60 days before a major's window start; `BPOST` bursts ≤ 60 days after a major's
  window end (nearer major wins); `BMID` other bursts; `EARLY` t − t0 ≤ 0; `ARR0` 0 < t − t0 ≤ 1 d; `AFT1` 1–3 d; `AFT3` > 3 d.
- `PHE` (column, "episode start and end"): `BTW`; major rows in the first / middle / last third of the major's
  [start, end] window = `E1`/`E2`/`E3` (rows outside the window go to E1 / E3).
- `PHp1…PHp200` (**episode-date permutation of PH**, S and P3 only, no bootstrap, no control): BTW unchanged; for each
  draw a random derangement π of the region's majors; fake onset t0'(f, e) = start(e) + [t0(f, π(e)) − start(π(e))],
  or the median of that offset over the families of π(e) if f was not active in π(e). Seed 20261013 + draw.
- Control rows (c2 windows) get PH/PH5/PHE labels from their own time with the data onset of their family-major.

**Pairs on common rows** (exact row match, common bootstrap draws; S, S3c and P3 differences):
- SVI: the 18 round-2 pairs (6 PGC-quartet + 12 far × quartet).
- WA-NW: anchors {B005, B006, B007, LRIV, B013, B943} × every other station with windows (anchor-anchor pairs once).

**Site-joint P and S (task 4)**, groups: SVI PGC site {PGC, B009, B010, B011}; WA G1 {B005, B006, B007, LRIV},
G2 {B013, B943}. Per (family, group, scheme, label): rows common to all members that have the label (members with fewest
common rows dropped until ≥ 20 common rows and ≥ 2 members remain); per member the stack and its leave-bin-out reference;
joint lag = peak of the equal-weight mean of the members' cc(τ) curves, parabolic refinement; S, P3 and S − P3 with
common bootstrap draws. Members' own lags on the common rows are kept.

## 1. WA-NW replication of the long-term pooled bound (`s11_ball3_clock.py`, `s11_ball3_longterm.py --region WA_NW`)
- **Clock terms** as round 2 §1, generalised to groups: per (family, bin) LS of the within-group pair differences, gauge =
  median of members (mean for 2), leave-family-out average over families (≥ 3 other families). G1 = B005/B006/B007/LRIV,
  G2 = B013/B943 (SVI: the PGC quartet, reproducing round 2).
- **References:** `G1BH` = mean of the clock-corrected X − A over A ∈ {B005, B006, B007} (single epoch from 2005; headline);
  `G2` = mean over {B013, B943} (headline, independent site); columns `LRIV`, `G1all` (incl. LRIV).
- **Items, YR bins headline, Y2 column:** dS for every far station X against G1BH and G2; dP3 the same; co-located dS
  inside G1 and G2 (method control); G1BH − G2 (two borehole sites ~60 km apart); S − P3 per station (task 3).
- Same statistics as round 2 §3: demean per series, ≥ 3 years per series; year decomposition with within-series
  year-label permutation (10,000 draws); **decade contrast = years ≥ 2015 minus ≤ 2014** (≥ 2 years each side; WA data
  start 2007); split contrast and linear trend as columns; scrambled control.
- **Detection (per far station):** |Δ| ≥ 3σ AND permutation p < 0.01 AND scrambled-control p ≥ 0.01 AND the same sign with
  p < 0.05 with BOTH references G1BH and G2 (two independent sites). Otherwise a 2σ bound (formal; conservative =
  larger of formal and permutation SD).
- **Conversions** as round 2: T_S = 1.75/0.75 · SP_obs, T_P = SP_obs/0.75; far leg only: −Δ/T_S,X; uniform on both legs:
  −Δ/(T_S,X − T_S,ref), reported only when |T_S,X − T_S,ref| ≥ 2 s.
- **Pass rule 1 (replication):** (a) method control: co-located decade contrasts inside G1 and G2 within 3σ of zero; and
  (b) at least one WA far leg with ≥ 20 series gives a 2σ leg bound ≤ 0.03 % (2 × SVI's 0.014 %) with both references,
  or a detection. Every pair is reported (`tables/wa_pairs_decade.csv`).

## 2. ETS-phase composite, both regions (`s11_ball3_phase.py`)
- **Items:** S − P3 per station; inter-station dS per pair (clock-corrected, references as §1 and round 2: SVI far pairs
  with the borehole mean B009/B010/B011, WA far stations with G1BH and G2); absolute station S lag (column; carries the
  detector frame).
- Per item: demean each series over its phase labels (≥ 2 labels); family-common phase curve = 1/σ² mean per label;
  one-way decomposition with within-series label permutation (10,000).
- **Pre-registered contrasts** (per series, then 1/σ² pooled, σ inflated by √(χ²/dof) if > 1):
  `C1 = ARR − BTW` (front arrival vs between-ETS), `C2 = AFT − ARR` (behind the front vs arrival),
  `C3 = AFT − BTW`.
- **Region-wide pooling:** per-series contrasts converted to δ units (S − P: constant Vp/Vs −C/SP, Vs only −C/T_S,
  Vp only +C/T_P; dS: far-leg −C/T_S,X), then pooled over all series of the region; σ = the larger of the formal error
  and the random-effects error from the scatter of per-station means.
- **Controls:** (i) within-series label permutation; (ii) episode-date permutation (PHp, 200 draws): the pooled contrast
  recomputed with fake onsets, p = fraction of draws with |C| ≥ |C_obs|; (iii) scrambled control (c2) through the same
  path.
- **Detection:** |C| ≥ 3σ AND label-permutation p < 0.01 AND date-permutation p < 0.01 AND scrambled-control p ≥ 0.01 AND
  the same sign in the other region (p < 0.05). Region-only (passes all but the cross-region clause) = "candidate".
  Otherwise a 2σ bound in ms and %. PH5 and PHE are columns (no detection rule on them).

## 3. Region-wide (common-mode) change (`s11_ball3_common.py`)
- **S − P time series per station** (YR, P3 headline; P3s and site-joint columns), 2003–2026 SVI, 2007–2026 WA; family-common
  year curve and decade contrast per station (rules as §1).
- **Pooled over stations:** per-series decade contrasts in δ units, pooled region-wide; σ = larger of formal and the
  random-effects error across stations (so instrument-level P effects such as B010 vs B011 enter the error).
- **Catalogue timing drift vs real change:** absolute station-level S and P3 lags (YR) carry the detector frame O(b);
  O_S(b) and O_P(b) = family/station-common year terms of the absolute S and P lags. Catalogue drift moves S and P alike
  (O_S − O_P = 0) and leaves S − P flat; a real uniform change moves S − P. Reported: decade contrast of O_S, O_P,
  O_S − O_P and of pooled S − P; "catalogue drift" is called if O_S (or O_P) decade contrast ≥ 3σ while S − P is within 2σ.
- **Joint uniform (a = δVs/Vs, b = δVp/Vp):** weighted LS on per-series decade contrasts with design S − P: −a·T_S + b·T_P;
  dS(X − ref): −a·(T_S,X − T_S,ref); dP3(X − ref): −b·(T_P,X − T_P,ref); σ inflated by √(χ²/dof); within-series
  year-permutation p for a and b (2,000 draws). Also b − a (Vp/Vs change).
- **Per ETS:** EP bins (majors), region-wide δ(S − P)/SP per major (pooled as above), σ per major; χ² of majors about zero
  against the scrambled control's; 32/15-major composite = the ETS mean (relative to bursts it is §2 C1).
- **Success rule 3:** a region-wide 2σ bound on uniform δVs/Vs (Vs only) < 0.03 % or constant-Vp/Vs < 0.05 % over the
  decade in either region, or a detection (|Δ| ≥ 3σ, permutation p < 0.01, control p ≥ 0.01, same sign in both regions).

## 4. Joint multi-station P alignment (`s11_ball3_jointp.py`)
- Compare site-joint P3 with the best single member on the same common rows: median bootstrap σ_P, median P cc, rail rate,
  S − P decade σ and year-scale robust scatter.
- **Pass rule 4:** median σ_P(site-joint) ≤ 0.7 × median σ_P(best single member) on YR or EP bins at ≥ 2 of the 3 sites
  (PGC, G1, G2), and the control's joint P cc stays < 0.3. If it passes, the site-joint S − P becomes a headline column in
  §3 (the pre-registered headline stays P3 single-instrument).

## 5. Order and conduct
Queue (`s11_ball3_queue.sh`, nohup, DONE markers, idempotent; restart line in `scripts/resume_after_reboot.sh`):
measure SVI → measure WA-NW → clock (both) → longterm WA → phase (both) → common (both) → jointp → figures.
A failed rule leads to the next adaptation, not to a stop; each failure is written up with the next thing to try.
Amendments are added below, time-stamped, before the number they affect.

## Amendment 1 (Mon 2026-10-12, 7:40 pm PT — after single-family test runs `*_test.*`; no pooled number seen)
Tests: SVI family 10 and WA family 77 (`svi/tables/*_test.*`, `wa_nw/tables/*_test.*`). Seen: the round-3 engine reproduces
round 2 exactly on data stacks (family 10: 509 matching EP/YR/Y2 stacks, max |ΔlagS| = |ΔlagP3| = 0; control stacks differ
only by the random subsample); family-10 PH lags per station; site-joint vs best-member σ_P medians (SVI PGC site and WA
G1/G2: joint within 0.9–1.0 × best member, i.e. the bootstrap σ looks dominated by row sampling shared by co-located
instruments); WA G1 P3 cc 0.84–0.88, σ_P 6–8 ms on family 77. No rule changed. Timing: 95 s (SVI fam 10), 186 s (WA fam 77)
on one worker. Queue launched 7:39 pm PT (PID in `logs/queue.pid`).

## Amendment 2 (Mon 2026-10-12, 7:45 pm PT — while the SVI measurement runs; no pooled number seen)
§2 phase composite, dS items: the episode-date permutation draws (PHp) carry no clock terms (their labels are fake), so
the **headline phase dS is the uncorrected reference mean** (real and permuted labels treated identically); the
clock-corrected value is a column (`dScorr`). SVI phase dS reference = the borehole mean (B009/B010/B011) for X = B926,
SNB, VGZ (PGC vs that mean is the co-located control item); WA references G1BH and G2 (LRIV vs G1BH and B013 vs B943 are
co-located controls, excluded from the region-wide pool). Date-permutation p uses the data series' bootstrap variances as
weights and drops a series from a draw if its fake ARR/AFT stack has < 20 rows. Label-permutation draws: 10,000 for
region-wide items and per-station S − P, 2,000 for other per-station items.

## Amendment 3 (Mon 2026-10-12, 7:50 pm PT — after the SVI clock/longterm/phase/common stages ran as a pipeline check; before any WA-NW number)
Seen for SVI: the long-term stage reproduces round 2 (B926 leg −0.13 ± 1.29 ms); phase contrasts, common-mode decade numbers
and the per-major region-wide S − P χ² (data 2.8 vs control 0.8). Because per-major excess scatter cannot by itself tell a
common (regional) change from station-level P/S timing effects, one test is added (no rule above is changed):
**split-network coherence per major** — per-major region-wide S − P (constant Vp/Vs, demeaned over majors per series) is
computed separately from two disjoint station sets (SVI: PGC site {PGC, B009, B010, B011} vs far {B926, VGZ, SNB}; WA:
borehole groups G1+G2 {B005, B006, B007, LRIV, B013, B943} vs all other stations). Pearson r across majors; permutation p
from 10,000 shuffles of the major labels of one set. A common per-ETS S − P change is called only if r > 0 with p < 0.01 in
both regions; the same on the scrambled control is reported. The SVI result of this test is labelled "computed after the
SVI χ² was seen".

## Amendment 4 (Mon 2026-10-12, 7:55 pm PT — after the SVI split-network numbers; before any WA-NW number)
SVI station-split coherence came out r = 0.48 (p 0.005), but the large per-major values sit in small majors with 1–4
active families, which both station halves see. A source-side effect of those few families (e.g. a depth or locus shift
of the active patch) is common to all stations, so the station split cannot separate it from a regional velocity change.
Added: **family-split coherence** — the same per-major region-wide S − P from two disjoint family halves (even vs odd
family id), Pearson r and 10,000-shuffle p. A regional per-ETS change needs r > 0 with p < 0.01 for BOTH the station split
and the family split, in both regions. Columns: Spearman r, and both splits restricted to majors with σ ≤ 0.05 % in both
halves. Diagnostic column (no rule): per-major station-level δ(S − P) in ms regressed on SP (velocity-type ∝ SP vs a
constant offset, depth-type).

## Amendment 5 (Mon 2026-10-12, 8:06 pm PT — after rule 1 was judged on the WA year bins: FAIL; written before any number below)
Rule 1 failed in WA-NW because WA families have few repeats per year: only ≤ 25 family series per far station reach
≥ 3 year stacks with n ≥ 20 (SVI: 43–88), and the co-located controls show decade drifts (B013−B943 +1.4 ± 0.3 ms,
B013 S − P −7.7 ± 1.7 ms against B943 −1.4 ± 1.8 ms). **Next adaptation (both regions, columns to the pre-registered
year-bin headline):** stack each family's rows directly in two decade halves, `DH` = {E: ≤ 2014, L: ≥ 2015} (one stack per
half, n ≥ 20 each), and in 4-year bins `Q4` (≤ 2006, 2007–10, 2011–14, 2015–18, 2019–22, 2023–26) for the time course.
Same estimators, pairs, references, site-joint and scrambled control (`s11_ball3_measure.py --schemes DH,Q4 --tag _dh`).
Analysis (`s11_ball3_dh.py`): per series contrast C = L − E (σ² = σ_E² + σ_L²); group clock terms estimated on the
contrasts of within-group pairs (median gauge, leave-family-out mean over ≥ 3 other families); reference means as §1;
1/σ² pooling with χ² inflation; p from 10,000 random sign flips of the per-series contrasts; scrambled control the same.
Rule 1 is re-judged on DH as "passes on the DH adaptation" with the same thresholds (co-located within 3σ; a far leg with
≥ 20 series and 2σ leg bound ≤ 0.03 % with both G1BH and G2, or a detection with the §1 clauses). Region-wide S − P and the
joint (a, b) fit are repeated on DH contrasts (columns to §3). Interpretation rule fixed now: an S-time shift seen
against BOTH other references (G1BH and B001-type far stations) at a site whose S − P is unchanged at one co-located
instrument is a site timing (clock) term, not a path change.
Amendment 5a (8:08 pm PT, before any DH number): with only two bins the leave-bin-out reference of L is the E stack and
vice versa, so lag(L) − lag(E) = 2 × the shift; the DH contrast is therefore C = [lag(L) − lag(E)] / 2, σ² = σ_E² + σ_L².
