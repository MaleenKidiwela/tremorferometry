# G2 — Local earthquakes and deep LFEs under SVI and Washington, both directions. PLAN (pre-registered)

Written Sun 2026-10-11 ~6:40 am PT, **before any G2 statistic was computed**. Looked at while writing: the download
row counts, the merged catalogue's class counts and its list of M ≥ 4 events (time, place, depth, distance; no LFE
count, no rate, no completeness estimate). Changes after this point go to "Revisions" with a PT time and the reason.
Science only: no catalogue file is changed. Scripts `svi_science/scripts/g2_local_*.py`; outputs here (`tables/`,
`figures/`). CPU ≤ 8 cores, RAM < 60 GB. "Scrambled control" = the same detector on the same day with each station
site's data shifted to another time of that day. Times shown to the reader are PT; computations use UTC.

Reused from study G (`results/g_trigger/`): its LFE / candidate / observable-day / episode tables, the random-date
machinery, the ETS local-background test, the named-earthquake rule and the eye-check code. Main lesson of G: earthquake
shaking makes false matched-filter hits; P removes most of them, the scrambled control does not see them, and P ≥ 0.75
still leaks a few rows inside regional aftershock sequences → record sections are a veto.

## 0. LFE data and usage rules (as G, PLAN §0 there)

- **SVI 2003-01-01 … 2026-09-15**: `g_trigger/tables/svi_events.parquet` = `r23/tables/events_2003_2026.parquet`
  (backward 2003–2013 under RESULTS_BACKWARD rules: primary detector, era rows P ≥ 0.90, frozen rows P ≥ 0.75,
  counted once; forward 2014+ keep_xfam8 and P ≥ 0.75; data and scrambled-control rows through the same steps).
  Network segment `seg` (CN4f, ERA, B2f, B3f, S15f, FWD). Multi-site `n_z2`. Observable family-days per segment.
  Episodes 2003–2026: `r23/tables/episodes_2003_2026.csv` (145; pre-2010 = 6 Bostock windows + 39 envelope-tremor
  minors = the pre-2010 envelope flag, `validated == False`). Between-episode rates are a control-subtracted bracket;
  nothing is trended or compared across network segments.
- **WA-NW 2007-01 … 2026-09**: `g_trigger/tables/wa_events.parquet` (core net, NW, P_own ≥ 0.75, counted once
  = dedup D False), eras E1/E2/E3, n_z2, observable family-days, episodes `episodes_nw_final.csv` (118). Before
  2009-08-06 only the 3 WECC windows are listed → 2007-01 … 2009-08-05 is a separate column, not headline.
- **Tremor**: PNSN epicentres 2009-08 … 2026-05-27 (`catalogs/pnsn_tremor_cascadia_full.csv`), as columns.
- **P-weighting**: every LFE count is reported two ways: Σ P over the rows passing the usage rules (**headline**)
  and the plain count (column). Low-P rows are not added (they carry the earthquake false hits, lesson of G).
- **Artefact gauge**: raw matched-filter candidates of any P (data and scrambled), as in G.

## 1. Local earthquake catalogue (`g2_local_build.py` → `tables/eq_local.parquet`)

1. ComCat (USGS FDSN, PNSN `uw` + `us`), all event types, M ≥ 1.0, box 45.9–50.5 N, 126.8–120.4 W, 2003-01-01 …
   2026-10-11; NRCan FDSN event service (QuakeML, for the event type), same box, by month with retries (0 failed
   months). Raw 33,540 + 12,655 rows.
2. **Duplicates**: ComCat–NRCan pairs within 5 s and 25 km (4,988 pairs). The row kept follows the border (US side
   ComCat, Canadian side NRCan, crude polygon along 49 N / Haro Strait / mid Juan de Fuca Strait); its magnitude is
   used; a fixed depth (NRCan "operator assigned", ComCat no depth error) is replaced by the other row's free depth.
3. **Removed**: every non-earthquake type (explosion, quarry blast, mining explosion, induced) in either catalogue;
   events within 15 km of Mt Baker, Glacier Peak, Mt Rainier, Mt St Helens, Mt Adams (volcanic swarms). Neither
   catalogue holds LFE / tremor entries by type (ComCat types: earthquake / explosion / quarry blast / mining
   explosion; NRCan: earthquake / mining explosion / induced). Remaining LFE-like risk = mislocated tremor logged as an
   earthquake during ETS; checked in the eye check and by the column "near-interface events within 20 km of a family
   during ETS".
4. **Classes** (Slab2 `cas_slab2_dep_02.24.18`, the repo standard; McCrory 2012 drafted surface = column),
   dz = hypocentre depth − slab-top depth at the epicentre: **intraslab** dz > +5 km; **near-interface** |dz| ≤ 5 km;
   **crustal** dz < −5 km and depth < 20 km; **overriding-deep** dz < −5 km and depth ≥ 20 km (deep crust / forearc
   mantle; reported, not a headline class); offshore-no-slab; unknown-depth (fixed depth; only in "all").
5. **Distance to the families**: minimum hypocentral distance to the region's families (SVI 115; WA-NW quality A/B
   primaries, 100). Region set R = events within **150 km** of R's families. An event can be in both region sets.
6. **Stress estimates at the nearest family** (order of magnitude, ranking variables only): M0 = 10^(1.5 M + 9.1);
   static Δσ_s = M0 / (4π r³); dynamic Δσ_d = μ PGV / β with the far-field Brune S pulse
   PGV = 0.6 · M0 · ω_c² / (4π ρ β³ r) · exp(−π f_c r / (Q β)), f_c = 0.37 β (16 Δσ / 7 M0)^(1/3), Δσ = 3 MPa,
   β 3.7 km/s, ρ 2,900, μ 30 GPa, Q 300. Magnitudes are treated as Mw (ML / Md mixing noted).
7. **Declustering**: Gardner & Knopoff (1974) windows, largest first, on all kept events within 150 km of either
   region (M ≥ 1.0). Mainshocks are the trial / rate unit; aftershocks are not counted separately.
8. **Completeness Mc_R** (estimated after this plan): maximum curvature + 0.2 (0.1-M bins) on all kept events of
   the region set, in 3-year blocks 2003–2026; **Mc_R = the largest block value**, never below 1.5. Same estimate on
   ETS days vs other days (detection check, §3).

## 2. Q1 — do local earthquakes trigger LFEs? (per region SVI, WA-NW)

**Trials**: GK mainshocks with M ≥ Mc_R in region set R, origin inside the region's catalogue span with an episode
list (SVI 2003-01-01 … 2026-05-27; WA-NW 2009-08-06 … 2026-05-27; WA 2007-01 … 2009-08-05 column).
**Windows** (t0 = origin; UTC day 0 = day of t0): hours **H6 = [t0, t0 + 6 h] (headline)**, H1 = [t0, t0 + 1 h] and
H24 = [t0, t0 + 24 h] (columns); "before" = same length ending at t0; days **D3 = days +1…+3 (headline delayed)**,
D10 = +1…+10 (column). Status (G rule): between (every day touched by [t0 − 6 h, t0 + 6 h] outside each listed
episode ± 2 d and observable), ets (a touched day in an episode), pad, not observable.

**Hours test** (between trials): per trial Σ P in the window for data, scrambled control, data with n_z2 ≥ 2, plus
counts, raw candidates (gauge), PNSN tremor ≤ 100 km of the region reference. Stacked over trials:
S = Σ P_data, X = Σ (P_data − P_c1), S_mz = Σ P_data,n_z2≥2, BA = Σ [(P_data − P_c1)(after) − (…)(before)].
**Random-date control**: per trial 2,000 windows at the same UTC time of day on random days of the same segment / era,
between episodes ± 2 d, observable, not within ±3 d of the trial day and not within ±1 d of any region-set
earthquake with M ≥ 3.0; 10,000 random stacks (one random window per trial) → p_R(S), p_S(X), p_mz(S_mz), p_BA.
One-sided, p = (1 + #{random ≥ observed}) / 10,001. **Clustering check** (column): segment-wise circular shift of all
trial times by one common random whole-day offset, 2,000 draws, statistic mean X per eligible trial.
**ETS column**: ETS-status trials vs their own ±24 h local tiles (G method) for H1 and H6.

**Delayed test** (G test 2, unchanged): rate = Σ(data − c1) counted-once LFEs per observable family-day on days
+1…+3 (+1…+10) vs −3…−1 (−10…−1), between-episode days only, ≥ 2 (≥ 5 for D10) such days each side; random-date
control on all eligible same-segment days; multi-site leg n_z2 ≥ 2; scrambled control alone must not pass.

**PASS (per stack)**: p_R < 0.01 AND p_S < 0.01 AND p_mz < 0.01 (delayed: G's rule) AND the eye check (§5) does not
show the excess to be earthquake energy. **Headline stacks** (per region): H6 and D3, between episodes, for
**all classes, crustal, intraslab, near-interface** → 16 headline tests (≈ 0.16 chance passes expected at 0.01 per
leg if the legs were independent; they are not, so fewer). **Columns**: magnitude (Mc–2.5, 2.5–3.5, ≥ 3.5),
distance (< 30, 30–75, 75–150 km), dynamic-stress tier (≥ 1 kPa, 0.1–1, < 0.1), static ≥ 1 Pa, overriding-deep,
McCrory classes, per segment / era, H1 / H24 / D10, ETS column, Q3 bins.
**Dose-response**: G test 3 on Δσ_d (Poisson N ~ μ(1 + βσ), 99 % profile bound; Spearman ρ with permutation p);
"scales" only if ρ > 0 with p < 0.01 and the β lower bound > 0. Same on Δσ_s as a column.
**Named earthquakes**: G rule (a)–(e) on H6 (between: random windows; ETS: local tiles).
**Upper bounds** (if no pass): 99 % bound on extra Σ P per earthquake and relative to the random mean, per stack.

## 3. Q2 — does deep slow activity trigger local earthquakes? (per region R, per class)

**Events**: GK mainshocks, M ≥ Mc_R, region set R; classes crustal, intraslab, near-interface (headline),
overriding-deep and all (columns). Day span as in §2 (SVI 2003-01-01 … 2026-05-27; WA-NW 2009-08-06 … 2026-05-27).

**(a) ETS episodes.** E = days inside R's listed episodes (all episodes; majors-only column), O = days outside every
episode ± 2 d; pad days are in neither. RR = (N_E / |E|) / (N_O / |O|).
Controls (10,000 each): (i) **shuffled episode dates by circular shift** of the whole episode list by one random
whole-day offset modulo the span (keeps spacing and durations); (ii) **random placement**: each episode's start
drawn uniformly, overlaps rejected. p one-sided (RR ≥ observed); the low tail (rate lower in ETS) is reported as a
column. Columns: 10 days after each episode end (vs O), 10 days before each start; PNSN-tremor days (≥ 10
epicentres ≤ 100 km of the region reference, 2010-01-06+) vs tremor-free days; split pre-2010 / post-2010 (SVI).
**Multi-site LFE leg**: LFE-active days (≥ 5 counted-once data LFEs with n_z2 ≥ 2 in R that day) vs LFE-silent
observable days, same two controls applied to the active-day list.
**PASS**: p < 0.01 against both controls AND the LFE leg p < 0.01 against both, AND the detection check holds.

**(b) Between-episode deep LFE bursts.** Burst = counted-once data LFEs outside every episode ± 2 d, chained by
gaps < 24 h, ≥ 3 LFEs with ≥ 2 having n_z2 ≥ 2 (the Q5 definition plus a multi-site requirement). Burst window =
UTC days from first to last LFE plus the next day. RR_B = rate in burst windows / rate on the other between-episode
days. Controls: (i) circular shift of the burst list along the region's between-episode observable days, (ii)
uniform random placement on between-episode days of the same segment / era. Scrambled-control bursts (same
definition on c1) reported. PASS: p < 0.01 against both. Column: windows extended to +3 d.

**Detection check (veto)**: Mc on E days vs O days (same estimator, bootstrap 1,000 → 90 % range); masking index =
(N with Mc_R − 0.5 ≤ M < Mc_R) / (N with M ≥ Mc_R) in E vs O. If Mc(E) > Mc_R, the ETS test is re-run at
M ≥ Mc(E) and that run is the headline.

**Upper bounds** (if no pass): 99 % bound on RR from the control distribution (observed RR / 1st percentile of the
shuffled RR, and Poisson-exact on the counts), per region × class.

## 4. Q3 — where?

Position of each earthquake relative to R's LFE zone (`g2_local_build.py`): Slab2 depth at the epicentre vs the
families' slab-depth 5–95 % range (SVI 34.6–43.3 km; WA 33.1–49.9 km) → **up-dip / LFE-zone / down-dip**, only within
the families' along-strike span ± 15 km (the rest = "outside along strike"); **beneath** = intraslab with epicentre
≤ 15 km from a family; **near** = any class ≤ 30 km hypocentral from a family. Q1 (H6, D3) and Q2 (a, b) are re-run
in each bin with the same controls. These are columns: a bin result counts as "concentrated" only if it passes the
§2 / §3 rule by itself and the complementary bins do not.

## 5. Eye check (veto only)

2–8 Hz envelope record sections (20-Hz caches, G code), PGC broadband, LFE rows (data / scrambled / n_z2 ≥ 2),
PNSN tremor, for: the top 5 between trials per region by (Σ P_data − random mean) in H6; the top 3 ETS trials by
excess over local; every named earthquake; the largest local events (M ≥ 4.3 within 100 km of either region:
2003-04-25 M4.8 intraslab Olympic, 2009-01-30 M4.7 intraslab Kitsap, 2015-12-29 M4.7 intraslab Sidney/Victoria,
2019-07-12 M4.6 Monroe, 2020-01-24 M4.5 off-VI intraslab, 2023-10-08 M4.3 intraslab, 2025-03-03 M4.5 crustal).
Verdicts: "LFE-like bursts with moveout", "earthquake / aftershock energy", "ongoing ETS", "nothing".

## 6. Literature (from memory, labelled, not checked here)

Cascadia: slow-slip periods and local seismicity — reports that intraslab/crustal seismicity rates do not obviously
change during ETS in northern Cascadia (Vidale et al. 2011-type stacks; Wech & Creager-era descriptions) vs a few
reports of small rate changes; Mexico (Guerrero): slow slip events followed by / accompanied by intraslab and
interface seismicity changes (e.g. the 2014 Papanoa M7.3 during an SSE; Radiguet et al. 2016); New Zealand Hikurangi:
SSEs accompanied by swarms of small earthquakes (Delahaye et al. 2009; Bartlow et al. 2014-type studies); triggering of
tremor by local/regional earthquakes is rarer than by teleseismic surface waves. All will be labelled "from memory".

## Revisions
- **R1 (Sun 10-11 ~6:50 am PT, after the tests ran; implementation detail, no rule chosen by result).** The delayed test
  uses G's `delayed()` code unchanged: it needs ≥ 3 (D3) / ≥ 5 (D10) between-episode days after and ≥ 5 before (the
  plan said ≥ 2 for D3), and its random-day pool excludes ±3 d around region-set earthquakes with **M ≥ 3.0** (G
  excluded ±3 d around every listed earthquake, which with ~4,000 local trials would leave no days).
- **R2 (same).** Burst "random placement" control (§3b) draws from all between-episode observable days of the region
  span, not per network segment: the earthquake catalogue does not depend on the LFE network segment.
- **R3 (same).** The §3 columns (10 days after / before, pre/post 2010, tremor days) were run with the circular-shift
  control only (placement not run); they are columns, not headline.
- **R4 (~6:55 am PT, POST-HOC, labelled so in RESULTS).** `g2_local_posthoc.py`: (1) the rate at which the
  named-earthquake rule fires on random pseudo-earthquakes of the same status (46 trials met rule a–d, almost all
  during ETS, where the ±24 h local background has only 8 tiles); (2) burst starts within 6 h after a local
  earthquake vs random; (3) M ≥ 3.5 intraslab / near-interface earthquakes ≤ 15 km from a family that fall inside a
  between-episode burst window.
- **R5 (~6:57 am PT).** Eye check: figures made for the §5 list (top-5 between, top-3 ETS, the top named by Poisson p,
  the largest local events) plus SVI 2006-01-15 M3.6 (found by R4). The remaining ~40 named trials were not drawn:
  R4 shows the named rule fires on them at the random-pseudo-earthquake rate.
