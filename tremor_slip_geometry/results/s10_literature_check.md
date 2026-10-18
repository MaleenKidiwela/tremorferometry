# S10 — Literature check on the candidate results (2026-09-10)

Method: WebSearch/WebFetch against the published literature (publisher sites, PMC, Semantic Scholar API,
author/institutional PDF copies). Where a paywall blocked the full text I say so explicitly and mark the
statement as ABSTRACT-ONLY. Nothing below is inferred from memory; every number is either quoted from a
source I actually read or flagged as unverified.

Scope note: candidates are the ones handed to me (R9, R11, R8, R3/R5/R7, and "what did Gualandi do").
R10 (slip band brackets tremor band) is not on the list but a near-exact prior publication turned up while
searching, so it is recorded at the end as an unsolicited but important finding.

Timing note: RESULTS.md was updated with the S11 and S12 results while this check was running. The
verdicts for Candidates 1 and 2 and the final ranking below have been revised against the **new**
RESULTS.md (R9 with the energy exponent and the reversed smoothing caveat; R11 with the resolution test
passing). The prior-work tables are unaffected.

---

## Candidate 1 — Tremor rate is a SUB-LINEAR function of modelled slip rate (R9)

Our claim (as of the S11 update to RESULTS.md): daily PNSN tremor count ~ (modelled potency rate)^0.40
[0.31, 0.45], 27/27 events below 1; radiated tremor energy ~ rate^0.64 (2017+, 18 events, 14/18 below 1,
range 0.3-0.8); identical in north and central and in every epoch; robust to window saturation, network
growth, station count and smoothing (S11 simulation shows the smooth-rate fit is unbiased and reverses the
earlier "upper bound" caveat). Physical content: tremor productivity per unit slip falls as slip
accelerates.

### Prior work

| Paper | Data | Result |
|---|---|---|
| Aguiar, Melbourne & Scrivner (2009), *JGR* 114, B00A05, [doi:10.1029/2008JB005909](https://doi.org/10.1029/2008JB005909) | 23 Cascadia ETS events, 11 yr, GPS + tremor duration | "A highly **linear** relationship between moment release, as estimated from GPS, and total duration of the tremor events." Event-integrated totals, not daily rates. (Full text 403; result quoted from the USGS/publisher record.) |
| Bartlow, Miyazaki, Bradley & Segall (2011), *GRL* 38, L18309, [doi:10.1029/2011GL048714](https://doi.org/10.1029/2011GL048714) | Aug 2009 central Cascadia ETS, NIF inversion of GPS + PNSN tremor | "We find an excellent correlation between the instantaneous region of high slip rate and tremor epicenters for all times during the ETS event." Purely qualitative — **no correlation coefficient, no exponent**. Also: "tremor and slip started within one day of one another." (Read in full.) |
| Hawthorne & Rubin (2013), *JGR* 118, [doi:10.1002/jgrb.50103](https://doi.org/10.1002/jgrb.50103) | 5 Cascadia SSEs 2007–2011, borehole strain rate vs tremor envelope amplitude | Moment rate roughly **proportional** to tremor envelope amplitude on 30 min – 16 h timescales, but only over a factor 1.5–2 of amplitude variation. Amplitude, not count; sub-daily, not daily. (ABSTRACT/summary-level; full text 403.) |
| Wech & Bartlow (2014), *GRL* 41, 392–398, [doi:10.1002/2013GL058607](https://doi.org/10.1002/2013GL058607) | Cascadia, seismic + geodetic + strainmeter | Tremor "falls silent at reduced slip speeds when the migrating slip front pauses"; implies a **slip-rate threshold** for tremor genesis and that "the same section of fault can slip both with and without detectable tremor, limiting tremor's use as a proxy for slip." Threshold, not exponent. (ABSTRACT/summary-level; escholarship copy was empty, AGU 403.) |
| Frank & Brodsky (2019), *Sci. Adv.* 5, eaaw9386, [doi:10.1126/sciadv.aaw9386](https://doi.org/10.1126/sciadv.aaw9386) | Guerrero, Mexico; 2.5 yr, millions of LFEs, 8 geodetic SSEs | Calibration **Ṁ₀ᵍᵉᵒ = α Aᵝ**, A = median LFE displacement amplitude; "the best-fit values of α and β are 10^13.87 ± 0.19 and **3.17 ± 0.29**". Inverted, that is A ∝ Ṁ₀^0.32 — i.e. a **sub-linear exponent ≈ 0.32 for LFE AMPLITUDE vs slip rate**. They use LFE *count* per day only as a proxy for slow-slip *duration*, and never test count vs slip rate. (Read via PMC6774729.) |
| Bletery & Nocquet (2020), *Nat. Commun.* 11, 2159, [doi:10.1038/s41467-020-15494-4](https://doi.org/10.1038/s41467-020-15494-4) | Sep–Oct 2013 Cascadia SSE, 62 GPS, daily unsmoothed inversion | "The daily number of tremors correlates fairly well with the moment release rate (**coefficient of correlation: 0.63**) and even better with **daily maximum slip rate (0.76)**." Linear correlation on one 47-d event; **no exponent fitted**. (Read via PMC7195424.) |
| Passarelli, Selvadurai, Rivalta & Jónsson (2021), *Sci. Adv.* 7, eabg9718, [doi:10.1126/sciadv.abg9718](https://doi.org/10.1126/sciadv.abg9718) | Global compilation of tremorgenic + swarmgenic SSEs (Cascadia, Guerrero, Nankai, …) | "log(M₀ᵍᵉᵒᵈ) = 15.3 (14.9, 15.9) + **0.20 (0.17, 0.24)** log(M₀ˢᵉⁱˢ), with R² = 0.4"; deep-only subset gives 0.29 (0.27, 0.32), R² = 0.7. This is **between event-integrated moments across events**, not within an event, and the regressand is the geodetic moment. |
| Michel, Gualandi & Avouac (2019a), *PAGEOPH* 176, 3867–3891, [doi:10.1007/s00024-018-1991-x](https://doi.org/10.1007/s00024-018-1991-x) | 352 cGPS, 2007.0–2017.6, vbICA | 64 SSEs detected that "spatially and temporally match with the tremors activity". Match, no scaling. |
| Gualandi (2025), *GJI* 242(2) ggaf198, [doi:10.1093/gji/ggaf198](https://doi.org/10.1093/gji/ggaf198) | Cascadia vbICA, near-real-time | "The correlation with the maximum slip potency is typically of the order of **0.7** over the all Cascadia subduction." No exponent. (Read in full, see Candidate 5.) |
| Frank, Rousset, Lasserre & Campillo (2018), *Sci. Adv.* 4, eaat0661 (preprint [arXiv:1711.11528](https://ar5iv.arxiv.org/html/1711.11528)) | Guerrero 2006 SSE | Uses a **threshold** on daily LFE amplitude sum to label slipping days; states "the surface displacement rate is proportional to the slip rate … the evolution of the slip rate during slow slip mirrors the observed complex time-history of the low-frequency seismicity." **No calibration constant, no exponent.** |
| Huang & Hawthorne (2022), *Nat. Commun.* 13, 5697, [doi:10.1038/s41467-022-33158-3](https://doi.org/10.1038/s41467-022-33158-3) | Parkfield LFEs | LFE moment–duration is **linear**, matching SSE moment–duration; about source scaling, not rate-vs-rate. |

Nankai/Mexico: I searched specifically for an Obara or Frank paper reporting an *exponent* between tremor
rate and slip rate and **found none**. The recurring statement in the literature is the *linear* one
("daily counts of tremor increase proportionally to the size of transient aseismic slip", as quoted by
Passarelli et al. 2021), used as a working assumption rather than a measurement.

### Verdict: **B — quantitative extension of something known.**

That tremor rate tracks slip rate is thoroughly established (Bartlow 2011; Hawthorne & Rubin 2013;
Bletery & Nocquet 2020; Gualandi 2025). What is **not** in the literature is a fitted exponent for
**daily tremor count or daily radiated tremor energy vs modelled slip rate within events**. The field's
default is proportionality (exponent 1), stated as an assumption rather than a measurement — Passarelli
et al. (2021) quote it as "daily counts of tremor increase proportionally to the size of transient
aseismic slip". The two published sub-linear numbers are for different observables and different
aggregations: Frank & Brodsky's beta = 3.17 +- 0.29 for per-LFE median amplitude in Guerrero, and
Passarelli's 0.20-0.29 between event-integrated moments across a global compilation. Nobody has measured
the within-event, daily-resolution exponent in Cascadia.

**The Frank & Brodsky comparison is the single most valuable thing in this section — and it needs care.**
Their fit is Mdot_geo = alpha * A^3.17, i.e. **A ∝ Rate^0.315**. Two readings of our S11 energy exponent
are possible and they disagree:

- *Summed-amplitude reading*: daily **summed** energy ∝ R^0.64, and energy ~ amplitude^2, so summed
  amplitude ∝ R^0.32 — **numerically identical to 1/3.17 = 0.315**, across two subduction zones and two
  independent methods. If defensible, this is a headline-grade coincidence.
- *Per-detection reading*: Frank & Brodsky's A is explicitly the *median displacement amplitude of all the
  LFEs during the SSE*, i.e. a **per-event** quantity. Our per-detection energy exponent is
  0.64 - 0.40 = **0.24**, so per-detection amplitude ∝ R^0.12 — a factor ~2.6 shallower than theirs, not
  a match.

Decide which of the two is like-for-like, state it explicitly, and do not let the 0.32 coincidence into
the abstract until the aggregation is matched. (Guerrero LFE templates are per-family detections; PNSN
epicentres are 5-min windows that may contain several sources — that is the crux.)

### To be claimable we need

S11 has already retired most of the blockers I would otherwise have listed here (smoothing bias shown by
simulation to be absent, network growth / station count / saturation controlled, an amplitude-based
energy version added for 2017+). What remains:

1. **Model independence.** Everything so far uses the 13 fixed vbICA patterns. Repeat the exponent against
   at least one independent slip model — Bartlow (2020) or, better, the Costantino et al. (2026) denoised
   daily inversion, whose whole point is that "scaling laws are strongly influenced by the chosen
   detection threshold and the signal-to-noise ratio". A reviewer will ask whether 0.40 is a property of
   the basis. This is now the main open threat.
2. **Match the Frank & Brodsky aggregation** (above) so the numbers can be placed side by side; and
   report the per-detection energy exponent (0.24) explicitly, since that is the quantity a rate-and-state
   or asperity model actually predicts.
3. **Per-unit-area normalisation.** Normalise counts and energy by the concurrent active footprint so the
   exponent is a rate density, not a total. Wech & Bartlow's (2014) threshold picture predicts that most
   of the count increase at high slip rate is *area* recruitment, and that would be a physically distinct
   explanation for sub-linearity from "each patch saturates".
4. **Frame against the linear default explicitly.** The paper's contribution is "the assumption of
   proportionality that underlies tremor-as-slip-proxy (Wech & Creager 2011; Frank et al. 2018;
   Rousset et al.) is wrong by a factor that matters" — quantify what a 0.40 exponent does to a
   tremor-derived slip estimate over the observed dynamic range.
5. **Cite Wech & Bartlow (2014) as the mechanism candidate.** Sub-linearity plus their silent-slip
   threshold is a coherent story (a threshold-plus-saturation nonlinearity); say so rather than leaving
   the exponent unexplained.

## Candidate 2 — Persistent tremor WITHOUT resolvable slip in central Oregon, 44.3–45.4 N (R11)

### Prior work

| Paper | Data | Result |
|---|---|---|
| Brudzinski & Allen (2007), *Geology* 35, 907–910, [doi:10.1130/G23740A.1](https://doi.org/10.1130/G23740A.1) | Single-station tremor envelopes + tanh fits to PANGA GPS, whole margin | ETS — tremor **and GPS-detected slow slip** — is present all along Cascadia including central Oregon. Three recurrence zones: Wrangellia 14 ± 2 mo, **Siletzia 19 ± 4 mo**, Klamath 10 ± 2 mo; "in the Siletzia zone … station COR in Corvallis, Oregon … 19 ± 3 months". Their Fig. 2E/2F are North and South Oregon GPS *slow-slip episodes*. (Read in full.) |
| Schmidt & Gao (2010), *JGR* 115, B00A18, [doi:10.1029/2008JB006045](https://doi.org/10.1029/2008JB006045) | 1998–2008 GPS, 16 slip transients | "**Limited station coverage south of Portland (45.5° N latitude) restricted the analysis to events on the northern half of the subduction zone.**" The GNSS-coverage caveat for exactly your latitude band is already in print. (Summary-level; AGU full text 403.) |
| Wells, Blakely, Wech, McCrory & Michael (2017), *Geology* 45, 515–518, [doi:10.1130/G38835.1](https://doi.org/10.1130/G38835.1) | PNSN tremor 2009–Jul 2015, tremor density vs forearc faults | "intensity is **lowest in the central Oregon forearc**"; transverse low-tremor-density zones named Gales Creek (~45.5 N), Tillamook-Yamhill, **Corvallis (~44.5 N)**, Wildlife Safari; "low tremor occurrence beneath the central Oregon segment related to **Siletzia**, the accreted oceanic basalt terrane". K-S test: densities within 5 km of major forearc faults have "<<10⁻⁶ chance of being a random sample". Also: "it is possible that **rare slow slip may occur in the absence of detectable tremor** (Wech and Bartlow, 2014)" — the converse of your claim. (Read in full.) |
| Wech (2021), *JGR* 126, e2021JB022523, [doi:10.1029/2021JB022523](https://doi.org/10.1029/2021JB022523) | 180,000 epicentres 2017–2021 with amplitudes | "Southern Cascadia is most active, but the highest tremor energy rates occur in northern Cascadia. **Tremor in central Cascadia is systematically weaker and less frequent.**" |
| Bartlow (2020), *GRL* 47, e2019GL085303, [doi:10.1029/2019GL085303](https://doi.org/10.1029/2019GL085303) | GNSS decomposed into ETS / inter-ETS, whole margin | Time-averaged **ETS slip rate resolved along the whole margin including Oregon**; "slip correlates well with tremor locations"; "the highest resolution image of the ETS zone derived from geodetic data to date". This is the strongest published counterexample to "no resolvable slip in central Oregon". (Summary-level; AGU 403.) |
| Nuyen & Schmidt (2021), *G-cubed* 22, e2020GC009477, [doi:10.1029/2020GC009477](https://doi.org/10.1029/2020GC009477) | GNSS + tremor, southern Cascadia | Low-amplitude long-term SSE in southern Cascadia; the paper is about the *southern* gap, not 44–45.4 N. |
| Michel, Gualandi & Avouac (2019a), *PAGEOPH* 176, 3867–3891 | vbICA 2007–2017, 352 stations | Checkerboard shows "10–20-km-wide haloes of intermediate magnitude slip … between the zero and nonzero patches". A resolution test exists in the vbICA family — yours does not yet. (Summary-level; Caltech PDF had a TLS error, Springer paywalled.) |
| Wech & Bartlow (2014), *GRL* 41, 392–398 | — | Slip without tremor observed; the reverse asymmetry to yours. |
| Gualandi (2025), *GJI* 242 ggaf198 | — | **No** discussion of Oregon, **no** resolution or checkerboard test anywhere in the paper (verified by full-text search). |

### Verdict: **B — a quantitative and now resolution-controlled extension; NOT C, but two published results still have to be answered.**

The S12 update changes this materially. The version of R11 in RESULTS.md now reports 213 GNSS stations
recovered, nearest-station distance and misfit comparable to the north, a synthetic 30-km patch recovered
at 98.6 % potency with a detection threshold 0.91x the north, and a fixed-basis handicap of only ~20 %.
That is exactly the test Schmidt & Gao's (2010) coverage caveat demanded, and it passes. So the claim is
no longer "a statement about the model".

Three separate sub-claims still get different verdicts:

- *"Central Oregon is where slip-tremor co-location is weakest on the margin"* — **C in kind, B in
  degree.** The central Oregon tremor deficit and its association with Siletzia and forearc
  block-bounding faults is published (Brudzinski & Allen 2007; Wells et al. 2017; Wech 2021: "Tremor in
  central Cascadia is systematically weaker and less frequent"). Your rank-rho 0.40-0.46 vs 0.70-0.77 is
  a new measurement of a known phenomenon.
- *"Persistent tremor without resolvable slip at 44.3-45.4 N, 27-38 km, in all epochs"* — **A**, subject
  to the two reconciliations below. I found no prior report of tremor-without-geodetic-slip as a
  *persistent, localised* feature anywhere in Cascadia. The published asymmetry runs the other way:
  slip-without-tremor (Wech & Bartlow 2014; Hall et al. 2018; and Wells et al. 2017 explicitly write
  "it is possible that rare slow slip may occur in the absence of detectable tremor").
- *"a tremor-productive, slip-poor fault segment"* — **wrong wording, drop "tremor-productive".** Wells
  (2017) and Wech (2021) establish central Oregon as the *least* tremor-productive part of the margin.
  The correct framing is tremor-**present**, slip-**absent**: even the margin's weakest tremor has no
  modelled slip under it.

### To be claimable we need

1. **Answer Bartlow (2020) directly.** This is now the decisive external test, and it is not yet done.
   Bartlow decomposes GNSS into ETS and inter-ETS components and inverts for a *time-averaged* ETS slip
   rate over the whole margin including Oregon — which beats the per-event SNR problem by construction —
   and reports that "slip correlates well with tremor locations". If her time-averaged map has slip at
   44.3-45.4 N and yours has none over 16 yr, that is a discrepancy the paper must explain, not omit. If
   her map also shows a deficit there, that is corroboration and should be the headline citation.
2. **Answer Brudzinski & Allen (2007).** They report GPS slow-slip *episodes* at central-Oregon stations
   (their Fig. 2E/2F, station COR, Siletzia zone, 19 +- 3 month recurrence) using single-station tanh
   fits. Single-station detection is not a slip inversion, so this is answerable — but it has to be
   answered, because on its face it is a published detection of Oregon slow slip.
3. **Show the null is a null, not a gap in the basis' time functions.** S12b (a transient with its own
   time history) is the right remaining test; until it closes, the claim is provisional.
4. **State the bound in Mw.** "No resolvable slip" should be "no slip above Mw X in any 30-d window,
   2010-2026, at 95 % confidence", using the S12 recovery curve. That is the number a reader can use.
5. **Connect to the published mechanism, don't invent one.** Wells et al.'s (2017) fluid-escape model —
   forearc faults draining the overpressured megathrust and locally reducing tremor — makes a prediction
   for slip too. If slip is also absent there, that is either a stronger version of their model or
   evidence against it. Frame the result as a test of Wells et al. (2017) and it becomes interesting
   rather than anomalous.

## Candidate 3 — Within-event tremor DEEPENING (+2.3 km first→last fifth); shallower when most intense (R8)

### Prior work

| Paper | Data | Result |
|---|---|---|
| **Wech & Creager (2011)**, *Nat. Geosci.* 4, 624–628, [doi:10.1038/ngeo1215](https://doi.org/10.1038/ngeo1215) | ~95,000 WECC epicentres, 2006–mid-2011, northern Washington wedge | **The direct conflict.** Their Fig. 2b plots "median daily downdip distance" through four ETS episodes — the same statistic you compute — and they conclude: "individual slow-slip events are **initiated deep on the plate interface and migrate upwards**"; "Both ETS and inter-ETS tremor swarms initiate about **80 km** away from the updip reference line and slowly migrate updip, but ETS swarms propagate further, reaching shallower depths." Also: recurrence decreases exponentially with depth, "from the **15-month** ETS periodicity on the updip side to **~2 months** on the downdip side"; minor swarms <20 h concentrate downdip. (Read in full.) |
| Luo & Liu (2021), *Nat. Commun.* 12, 1959, [doi:10.1038/s41467-021-22232-x](https://doi.org/10.1038/s41467-021-22232-x) | Cascadia ETS model + observations | "Major ETS typically nucleate at mid-deep portion of ETS zone (>38 km depth) then **gradually expand up-dip** and fill in the entire ETS zone." (Read via PMC8010077.) |
| Peng, Rubin, Bostock & Armbruster (2015), *JGR* 120, 4317–4332 | Cross-station tremor, S. Vancouver Island | For the 2011 ETS: "tremor and slip initiated at the downdip edge of the ETS zone in northern Oregon, and epicenters **migrated updip over the next 5 days**." (Summary-level.) |
| Hall, Schmidt & Houston (2019), *G-cubed* 20, [doi:10.1029/2019GC008510](https://doi.org/10.1029/2019GC008510) | 2010 and 2012 Mw 6.8 ETS, ~72 GPS, extended NIF | "tremor occurs on the **leading edge** of propagating slipping regions, well ahead of the highest slip rates, independent of the along-strike propagation direction"; "We also observe **ongoing slip continuing several days after tremor has passed**." The nearest published statement about the late stage of an event — and it is about *slip* outlasting tremor, not tremor deepening. (Abstract read via Semantic Scholar; AGU full text 403.) |
| Rubin & Armbruster (2013), *G-cubed* 14, 5371–5392, [doi:10.1002/2013GC005031](https://doi.org/10.1002/2013GC005031) | POLARIS 2003–2006, sub-km relative locations | Secondary fronts "start at or within about 1 km of the main tremor front and propagate back along strike at rates of 10–20 km/h". |
| Peng & Rubin (2017), *EPSL* 464, [doi:10.1016/j.epsl.2017.01.046](https://doi.org/10.1016/j.epsl.2017.01.046) | 383 SSFs in LFE catalogues; 693 SSFs in the Cascadia tremor catalogue Aug 2009–Dec 2015 | Secondary slip fronts occur behind the main front and "propagate both along dip and along strike", "less commonly propagate up or downdip at or within 1–2 km of the main front". A late down-dip *tail* is therefore not unprecedented in kind. |
| Houston, Delbridge, Wech & Creager (2011), *Nat. Geosci.* 4, 404–409, [doi:10.1038/ngeo1157](https://doi.org/10.1038/ngeo1157) | 5 largest ETS 2004–2009, northern Washington | RTRs migrate back through fault that "has just ruptured in the past three days"; main propagation 7–12 km/day. Retrograde along strike, not down-dip. |
| Sweet, Creager, Houston & Chestler (2019), *G-cubed* 20, [doi:10.1029/2018GC007998](https://doi.org/10.1029/2018GC007998) | 4 LFE families spanning the transition zone, W. Washington | "LFE swarm duration, recurrence interval, and event size **decrease systematically with increasing depth**"; along-dip contrasts in tidal response. Deeper = smaller and shorter, which argues against a deep late tail being a big-amplitude effect. |
| Bombardier, Cassidy, Dosso & Kao (2024), *JGR* 129, [doi:10.1029/2024JB029159](https://doi.org/10.1029/2024JB029159) | 896 tremor episodes, 17 yr, northern Cascadia | Major episodes up-dip and minor episodes down-dip of the forearc mantle corner (bimodal in S. Washington and Vancouver Island). |
| **Wech (2021)**, *JGR* 126, e2021JB022523 | 2017–2021 amplitudes | In northern Cascadia there is "a clear shift to **higher magnitudes at shallower** portions of the megathrust as a result of a higher proportion of MeL > 1.25 events updip." **This is an alternative, purely spatial explanation for "tremor is shallower when most intense".** |
| Zhang & Houston (2025), *AGU Advances* 6, [doi:10.1029/2025AV001746](https://doi.org/10.1029/2025AV001746) | ~330,000 tremors, ~13,700 extracted migrations | Migration controlled by stress diffusion, diffusivity 10³–10⁵ m²/s, ~30-km-wide propagation zone; barriers at 48.5 N and 42.5 N. Nothing on systematic within-event deepening. |

### Verdict: **A for the specific statistic, but it currently CONFLICTS with the standard picture, and half of it (the intensity–depth part) is probably C.**

- The exact quantity — median tremor depth, first vs last fifth, over 27 slip-anchored ETS with a
  shifted-window null — is not published. Nobody has done the 16-yr, 27-event, phase-normalised version.
- But Wech & Creager (2011) plotted essentially the same statistic (median daily down-dip distance through
  four ETS) and reported the **opposite sign** for the ETS body: initiation deep, migration up-dip. Luo &
  Liu (2021) and Peng et al. (2015) say the same. You cannot publish +2.3 km of deepening without
  reconciling it head-on with Fig. 2b of Wech & Creager (2011).
- "Tremor is shallower when most intense" (rho −0.11, 19/27, p 0.05) is very likely a re-expression of
  Wech (2021)'s **spatial** depth–amplitude gradient (more MeL > 1.25 up-dip) plus the fact that the
  intense days are the up-dip main-front days. As stated it is **C**.

### To be claimable we need

1. **Reproduce Wech & Creager's metric on your 27 events** — down-dip distance from a fixed up-dip
   reference line, northern Washington wedge, median per day — and show where and why the sign flips.
   Most plausible reconciliation: their result is dominated by the *initiation* ramp (first days) in a
   restricted wedge; yours integrates the whole window, margin-wide, and is dominated by the *tail*.
   Split your phase axis into "onset ramp" and "tail" and report both. If the deepening is entirely in
   the last two fifths and the first two fifths shallow, the two results are compatible and that is the
   paper.
2. **Kill the detectability explanation.** Deep tremor is weaker and more continuous (Wech & Creager 2011;
   Sweet et al. 2019; Wech 2021). A growing/improving network preferentially adds deep, weak detections,
   and the *late* part of an event is the low-amplitude part — so an improving detector deepens the late
   median for free. Required: (a) fixed station subset, (b) an amplitude-thresholded repeat using
   Wech (2021) MeL for 2017+, (c) show the effect in single years, not only pooled.
3. **Kill the geometry aliasing.** Slab2 depth at a fixed down-dip position varies along strike; an event
   that propagates into a steeper segment deepens its median depth with no along-dip motion at all.
   Recompute with depth referenced to the local band's own concurrent p50 (you already have the
   machinery from O1) and with down-dip distance instead of depth.
4. **Justify 2.3 km against location error.** PNSN epicentre errors are 5–10 km horizontal (Wech &
   Creager 2008, cited in Bartlow et al. 2011) and depths are worse; you are quoting Slab2 depth at the
   epicentre. The argument has to be that the *median of thousands* is unbiased between the first and last
   fifth — which is exactly what a time-varying detector breaks. Show the epicentre-count and
   station-count per fifth.
5. Drop or re-frame the "shallower when most intense" half unless it survives a null that shuffles daily
   depth within the same day's spatial footprint (i.e. removes the Wech-2021 spatial gradient).

---

## Candidate 4 — Sixteen-year stationarity of the tremor–slip relationship (R3, R5, R7)

### Prior work — is there a published claim of secular evolution this would test?

| Paper | Data | Result |
|---|---|---|
| **Meng, Liu, Adhikari & Luo (2026)**, *GRL* 53, e2025GL120690, [doi:10.1029/2025GL120690](https://doi.org/10.1029/2025GL120690) | Cascadia tremor + SSE, stress perturbations | **Yes — this is the paper your null tests.** From the abstract: "we systematically examine major stress perturbations … and identify two distinct patterns: an incongruent tremor seasonal variation that cannot be explained by hydrological forcing and **a long-term change in the locations where ETS events nucleate and their propagation directions**. … the long-term ETS evolution likely reflects spatiotemporal stress redistribution on subduction fault interface through stress feedback of past major ETS events." **ABSTRACT-ONLY** — AGU returned 403 for both the HTML and the OA PDF, so I could not read the magnitude, latitudes or years of the claimed change. This must be obtained before writing. |
| Brudzinski & Allen (2007), *Geology* 35, 907–910 | 1989–2006 | Recurrence treated as a **fixed** terrane-controlled property: 14 ± 2 / 19 ± 4 / 10 ± 2 months. Implicitly a stationarity claim. |
| Wech & Creager (2011), *Nat. Geosci.* 4, 624–628 | 2006–2011 | Recurrence as a fixed function of depth (15 mo → ~2 mo). |
| Bartlow (2020), *GRL* 47, e2019GL085303 | Whole GNSS record | Long-term **average** ETS slip rate — a stationarity assumption by construction; segmentation is strong and along-strike variable. |
| Michel, Gualandi & Avouac (2019b), *Nature* 574, 522–526, [doi:10.1038/s41586-019-1673-6](https://doi.org/10.1038/s41586-019-1673-6) | 2007–2017 Cascadia SSE catalogue | M₀ ∝ T³ scaling. A *scaling* claim, not a secular one — your R5 (potency per detection stable over 16 yr) is adjacent but not a test of it. |
| Costantino, Radiguet, El Yousfi & Socquet (2026), *GRL* 53, e2025GL117446, [doi:10.1029/2025GL117446](https://doi.org/10.1029/2025GL117446) | 15 yr (2007–2023) denoised-GNSS daily slip inversion | "scaling laws are **strongly influenced by the chosen detection threshold and the signal-to-noise ratio**" — the direct warning that any apparent secular change in a threshold-defined catalogue is suspect. Reinforces your X5. |
| Wech (2021) + PNSN blogs | 2009– | Catalogue is non-stationary by construction: algorithm updated Feb 2018, new method gives +55 % epicentres with less scatter. |

I found **no** published claim of a secular drift in the *tremor–slip relative geometry* (across-strike
offset, depth offset, band edges, lead/lag, potency per detection, recurrence phase). Meng et al. (2026) is
the only secular-evolution claim for the Cascadia ETS zone that I could verify, and it is about **nucleation
location and propagation direction**, which your R3/R7 do **not** currently test.

### Verdict: **A for the specific bounds — but as a stand-alone null it is only publishable if it is aimed at Meng et al. (2026).**

The bounds themselves (depth offset −0.04, [−0.16, +0.09] km/yr; tremor depth percentiles ≤0.02 km/yr;
lead/lag ±0.6 d/yr; recurrence phase 0.00 d/yr, p 0.61) are new numbers. The strongest asset is the
back-slip control that *does* drift at +2.9 km/yr — that is the honest statement of what a spurious drift
looks like in this data set, and it is a genuinely useful contribution to a field that has just published
a secular-change claim.

### To be claimable we need

1. **Read Meng et al. (2026) and test its actual quantities.** Nucleation latitude per event vs year, and
   propagation direction (N/S) vs year, with the same shifted-window and back-slip controls. A null on
   quantities nobody claimed to change is not publishable; a null on *their* quantities is.
2. **Power analysis.** Inject synthetic drifts (0.1, 0.25, 0.5 km/yr in depth offset; 1 d/yr in lag) and
   report the smallest slope detected at 95 %. Without this the bounds are uninterpretable.
3. **Express the bounds physically.** −0.04 ± 0.1 km/yr over 16 yr is <2 km total, which is well inside
   one 28-km mesh cell and inside the PNSN depth error — say so, so a reader knows the null is a real
   constraint and not just a wide interval.
4. **Tie the network-growth control to the documented catalogue changes** (Feb 2018 algorithm change,
   +55 % epicentres) rather than to a generic "network growth"; and repeat the key trends on a fixed
   station subset. This also rescues R5's "×1.3 tremor-days growth" caveat.

---

## Candidate 5 — What did Gualandi (2025) itself compare against tremor? And Costantino et al. (2026)?

### Gualandi (2025), *GJI* 242(2), ggaf198, [doi:10.1093/gji/ggaf198](https://doi.org/10.1093/gji/ggaf198)

Read in full (Gold OA PDF). Single-author. Findings, verbatim where quoted:

- **Tremor is used only as a 1-D visual/temporal comparison.** The comparison lives entirely in Fig. 7:
  latitude binned into 100 bins, maximum slip potency rate per bin per epoch, PNSN tremor plotted as black
  dots on top, plus a normalised daily tremor-rate curve in the top panel.
- **The only quantitative tremor numbers in the paper**: *"The correlation with the maximum slip potency
  is typically of the order of 0.7 over the all Cascadia subduction. The time lags observed are typically
  between 0 and a few days, with positive values indicating that the slip potency rate precedes the
  tremors."*
- **Co-location offsets: none.** No distances in km, no depth comparison, no per-event geometry, no
  along-strike or along-dip offsets anywhere in the paper.
- **Scaling: none.** No exponent, no rate-vs-rate fit beyond the single correlation coefficient.
- **Caveats the paper itself raises, and which bear directly on our R4 and R9**: *"The time lags observed
  suggest that the aseismic deformation may slightly anticipate the seismic signature, but further
  analysis on the effect of the adopted filter are needed."* And on the low-pass (1/14 and 1/28 d⁻¹
  tested): *"more aggressive filters tend to increase the duration of the SSEs since they retain only
  lower frequencies, and may anticipate the rising of the slip rate for past events. It is thus important
  to be aware of these effects and conduct a detailed analysis for forecasting purposes and for scaling
  law analyses."*
- **No resolution test.** Full-text search for "resolution" outside the reference list and for
  "checkerboard" returns nothing. **No Oregon-specific discussion** either.
- Method: vbICA with missing data (Gualandi et al. 2016), automatic IC selection by cumulative PSD
  excluding 1 and 2 yr⁻¹ neighbourhoods, leave-one-out smoothing parameter, Tarantola linear inversion
  on triangular patches with a Radiguet-type exponential prior covariance; 18–37 ICs tested daily
  (29 in the Jan 2025 example; 36 tested / 17 retained on 27 Jan 2025); 232 GNSS stations retained;
  near-real-time from Jul/Aug 2024, but the latitude–time figures span the full record (the text discusses
  "large SSEs involving both the central (44°N–47°N) and northern (47°N–50°N) sections … e.g. in 2011,
  2016 and 2019").
- He explicitly disclaims novelty on the tremor correlation: *"This was expected, since this correlation
  has already been reported by several previous studies on Cascadia (Wech et al. 2009; Bartlow et al.
  2011; Michel et al. 2019a; Bartlow 2020)."*

**Bottom line: Gualandi (2025) leaves essentially everything our study measures unmeasured.** One margin-wide
correlation (~0.7) and a qualitative lag statement he himself flags as possibly a filter artefact. That is
the gap statement for the paper — and it also means our R4 (|lag| ≤ 5 d in 63 % of events, slip not lagging)
is the *first* controlled version of the lag he reports, and it must address the filter caveat he raises.

### Costantino, Radiguet, El Yousfi & Socquet (2026), *GRL* 53, e2025GL117446, [doi:10.1029/2025GL117446](https://doi.org/10.1029/2025GL117446)

Abstract read (Semantic Scholar); **full text 403 at AGU, so the tremor question is not fully settled.**

- 15 yr (2007–2023) of Cascadia slip from deep-learning-denoised daily GNSS, inverted for slip on the
  interface; SSEs extracted by ST-DBSCAN at four slip-rate thresholds → four catalogues.
- Headline results: "an upper bound for slow-slip moment rates exists"; "scaling laws are strongly
  influenced by the chosen detection threshold and the signal-to-noise ratio"; "slow slip nucleates as a
  two-dimensional expanding crack, propagating laterally when encountering the along-dip limits of the
  transition zone".
- **No tremor comparison appears in the abstract.** Their earlier work (Costantino et al. 2023,
  *Commun. Earth Environ.* 4, 435; preprint [arXiv:2305.19720](https://arxiv.org/abs/2305.19720)) is
  titled "…validated against tremors in Cascadia", so the group does use tremor as validation — but I
  could not verify whether the 2026 paper contains a tremor figure. **Get the PDF before claiming
  novelty against it.**
- Their "2-D expanding crack propagating laterally when it hits the along-dip limits of the transition
  zone" is the closest published competitor to our R10 and R8: it is a statement about along-dip
  confinement of slip within an event, derived independently.

---

## Unsolicited but important: R10 is largely already published

R10 states the slip up-dip edge sits ~11 km shallower than the concurrent tremor up-dip edge in 27/27 events.

**Hall, Houston & Schmidt (2018)**, *G-cubed* 19, 2706–2718,
[doi:10.1029/2018GC007694](https://doi.org/10.1029/2018GC007694): 80 GPS stations, **six** major ETS
2007–2016, northern Cascadia: *"slow slip of up to 2 cm extends **updip of tremor by about 15 km**
beneath the Olympic Peninsula … **Slip updip of tremor is a persistent feature of all six ETS events**"*,
with along-strike variation that "correspond[s] to changes in lithology of the overlying crust".

Peng & Rubin (2017) restate it as "slow slip may extend 10–20 km updip of the updip extent of tremor and
the along-dip dimension of the slow slip zone is about 60 km".

**Verdict on R10: C.** Our 11 km over 27 events / 16 yr is an extension of Hall et al.'s 15 km over 6
events / 9 yr — worth a sentence and a citation, not a headline. The genuinely new part of R10 is the
**no-drift bound on the edges** (0.03–0.6 km/yr), which belongs with Candidate 4.

---

## Ranked publishability

(Revised after the S11 / S12 update to RESULTS.md, which landed while this check was running.)

| Rank | Candidate | Verdict | Best venue | Why |
|---|---|---|---|---|
| **1** | **C1 — sub-linear tremor-vs-slip-rate scaling (R9)** | B | **GRL** stand-alone, or lead result of a **JGR-SE** paper | It is a *number* where the field has an assumption, and S11 has already killed the smoothing, saturation and network-growth objections and added an independent energy-based exponent (0.64) that is sub-linear too. Two observables agreeing is what makes it publishable. The Frank & Brodsky (2019) comparison gives it instant context — get the aggregation right first. Only real blocker left: repeat it on a slip model that is not the 13 fixed vbICA patterns. |
| **2** | **C2 — persistent tremor without resolvable slip, central Oregon (R11)** | B (was C before S12) | **GRL** if Bartlow (2020) is answered; otherwise **G-cubed**/**JGR-SE** section | S12 promoted this from a caveat to a result: the resolution test passes (98.6 % recovery, threshold 0.91x the north). The published asymmetry runs the *other* way (slip without tremor), so tremor-without-slip as a persistent localised feature would be new. But it is not defensible until Bartlow (2020)'s time-averaged Oregon ETS slip and Brudzinski & Allen (2007)'s central-Oregon GPS episodes are explicitly addressed, and S12b closes. Drop "tremor-productive" — Wells (2017) and Wech (2021) say the opposite. |
| **3** | **C4 — 16-yr stationarity (R3/R5/R7)** | A (null) | **GRL** *only if* framed as a direct test of Meng et al. (2026); otherwise a **JGR-SE** section | Nulls sell when they refute a named claim, and Meng et al. (2026) have just published a secular-change claim for this exact fault zone (nucleation locations and propagation directions). That is the target. The +2.9 km/yr back-slip artefact is a strong, quotable honesty result. Needs the power analysis and needs Meng's actual quantities tested. |
| **4** | **C3 — within-event deepening (R8)** | A for the statistic; C for the intensity half | **GRL** if it survives the controls; otherwise cut | Highest interest, highest risk, and the *least* controlled of the four now that S11/S12 have hardened the others. It contradicts Wech & Creager (2011) *Nat. Geosci.* — who plotted essentially the same statistic and found up-dip migration — so it is either a real refinement (onset ramp vs late tail) or a detectability artefact. The "shallower when most intense" half is probably Wech (2021)'s spatial depth-amplitude gradient and should be dropped. The S11 machinery (fixed station subset, MeL energies 2017+) can be pointed straight at this; do that before writing. |
| **5** | **C5 — what Gualandi did** | not a result | — | The **gap statement** for whichever paper goes first: Gualandi (2025) reports one margin-wide correlation (~0.7) and a lag "between 0 and a few days" that he himself flags as possibly a filter artefact, with no co-location offsets, no scaling, and no resolution test anywhere in the paper. Quote it in the introduction; it is the strongest single justification for the study. |

**Recommended package.** One **JGR Solid Earth** paper — "Relative geometry and scaling of tectonic tremor
and a 16-year vbICA slow-slip model in Cascadia" — with C1 as the headline, C2 as the second result now
that S12 supports it, C4 as the stationarity section aimed at Meng et al. (2026), R1/R2/R4 as the
co-registration foundation, R10 cited to Hall et al. (2018) as confirmation rather than discovery, and C3
held back until the detectability and geometry controls are run. If a fast GRL is wanted, **C1 alone** is
the cleanest — it is closest to done, and the sub-linear exponent measured two independent ways (count and
energy) is a self-contained, quotable result.

**Things I could not verify online** (must be obtained before writing):
- Meng et al. (2026) full text — magnitude, latitudes and years of the claimed long-term ETS change. AGU 403.
- Costantino et al. (2026) full text — whether it contains a tremor comparison. AGU 403.
- Bartlow (2020) full text and slip map — **now the highest-priority retrieval**, because it is the decisive external test for C2. AGU 403.
- Hall et al. (2018, 2019) full texts — only abstracts read; AGU 403.
- Wech & Bartlow (2014) full text — only the publisher/USGS summary read; no slip-rate threshold value obtained.
- Aguiar et al. (2009) full text — only the summary read; the regression slope and R^2 of the moment-vs-tremor-duration fit not obtained.
- Michel et al. (2019a) full text — Caltech PDF TLS error, Springer paywalled; the checkerboard "10-20 km halo" statement is second-hand.

---

# 2026-09-10 S10b

Two follow-ups requested by the coordinator: (1) Candidate 6, recurrence predictability; (2) download the
Bartlow (2020) supporting data and extract the central-Oregon statements from Brudzinski & Allen (2007)
and Schmidt & Gao (2010). Same rules as above: quoted text is text I actually read; anything paywalled is
marked ABSTRACT-ONLY.

## Candidate 6 — episode size vs the FOLLOWING recurrence interval

Our claim: episode tremor count correlates with the time to the NEXT episode (Spearman +0.33, p 0.0005,
195 pairs, 12 half-degree bands 44-50 N) more than with the time since the previous one (+0.16); the
model-slip size proxy shows nothing.

### FIRST — the label in the candidate statement is inverted

Under Shimazaki & Nakata (1980), spelled out verbatim by Gao et al. (2012, p. 357):

> "For the time predictable model, the time to the next event is equal to the static stress drop released
> by the most recent earthquake divided by the stressing rate, which implies that **a larger event requires
> a longer recovery time until the next one**. For the slip predictable model, the fault slip released for
> an event is proportional to the time interval since the last event. Thus, **a longer interevent period
> produces a larger event**."

Size -> the FOLLOWING interval is the **TIME-predictable** signature. Size <- the PRECEDING interval is
**slip-predictable**. Our result (+0.33 with the following interval, +0.16 with the preceding) is therefore
**weakly TIME-predictable**, which is exactly the phrase Gao et al. already used for northern Cascadia.
Publishing this as "slip-predictable" would be caught immediately, and Gao et al. (2012) is the paper a
reviewer would cite to do it.

### Prior work

| Paper | Data | Result |
|---|---|---|
| **Gao, Schmidt & Weldon II (2012)**, *BSSA* 102(1), 352-360, [doi:10.1785/0120110096](https://doi.org/10.1785/0120110096) | **11** Cascadia SSEs centred near the Olympic Peninsula, from the Schmidt & Gao (2010) catalogue; recurrence intervals from the starting days of successive transients at GPS station **ALBH**; size = geodetic seismic moment | **THE prior publication — same region, same sign, same magnitude.** p. 357: "we find that the seismic moment is **anticorrelated to the preevent recurrence interval** for each event with a **R2 of 0.48**, but has a **weak positive dependence on the postevent recurrence interval with an R2 of 0.16**. This suggests that SSEs in northern Cascadia are **moderately antislip predictable and weakly time predictable**." Abstract (p. 352): "The recurrence statistics of slow slip events on the northern Cascadia subduction zone are weakly time predictable and moderately antislip predictable (that is, the event size and preevent recurrence interval are anticorrelated), which may indicate that **healing between events strengthens the fault with time**." Fig. 6a/6b are the two scatter plots. Interpretation (p. 357): "the fault is gradually healed and strengthened with time between two sequential SSEs. The stronger the fault is at the end of the inter-SSE period, the smaller the strain release for the next event, which would explain the antislip predictable behavior we see." They also note "The weak time predictability for Cascadia SSEs is consistent with the lack of evidence for time predictability for seismogenic faults (Murray and Segall, 2002; Weldon et al., 2004)" and close with "we hope to explore this more as the catalog of events grows." (Read in full.) |
| Shimazaki & Nakata (1980), *GRL* 7, 279-282 | — | Defines the two models. |
| Gualandi, Avouac, Michel & Faranda (2020), *Sci. Adv.* 6, eaaz5548, [doi:10.1126/sciadv.aaz5548](https://doi.org/10.1126/sciadv.aaz5548) | Cascadia slip history 2007-2017 | Slow earthquakes are a **low-dimensional (<5) nonlinear chaotic** system, not stochastic; predictability horizon "of the order of days to weeks". System-level predictability — **not** a size-vs-interval test. |
| Michel, Gualandi & Avouac (2019b), *Nature* 574, 522-526 | 2007-2017 SSE catalogue | M0 ~ T^3. Scaling, not recurrence predictability. |
| Michel, Gualandi & Avouac (2019a), *PAGEOPH* 176, 3867-3891 | 64 SSEs, vbICA | Recurrence for Mw~6 "ranges between a few months to about two years" (as quoted by Gualandi 2025). No predictability test found. |
| **Wech & Creager (2011)**, *Nat. Geosci.* 4, 624-628 | ~95,000 epicentres 2006-2011, N. Washington | Recurrence "**decreasing exponentially with depth**", "from the 15-month ETS periodicity on the updip side to ~2 months on the downdip side", within an explicit **stress-threshold / stress-transfer** model in which each depth band fails when the stress transferred from below exceeds its own threshold. Not a size-vs-interval test — but this is precisely the physical model that *predicts* time-predictability, so it is the mechanism to cite. |
| Brudzinski & Allen (2007), *Geology* 35, 907-910 | Whole margin, 1989-2006 | Recurrence as a fixed terrane property (14+-2 / 19+-4 / 10+-2 months). No predictability test. |
| Bartlow (2020), *GRL* 47, e2019GL085303 | Whole margin | Time-averaged ETS slip rate. No predictability test. |
| Nuyen & Schmidt (2021), *G-cubed* 22, e2020GC009477 | 2006-2019 GNSS | One low-amplitude long-term SSE in southern Cascadia (Mw 6.4, 15-35 km); "the total fault slip amplitude is an order-of-magnitude smaller in Cascadia when compared to large events elsewhere, making long-term slow slip detection challenging". Not a predictability test. |
| Hawthorne & Bartlow (2018), *JGR* 123, 4243-4265, [doi:10.1029/2017JB015124](https://doi.org/10.1029/2017JB015124) | Cascadia moment-rate power spectrum from GPS inversions, borehole strain and tremor beamforming | Spectrum flattens at 1-10 d then decays as f^-nm, nm = 1.1-1.4 (tremor) / 0.4-1.5 (strain). Spectral structure, **not** recurrence predictability. (ABSTRACT-ONLY.) |
| Beeler, Thomas, Burgmann & Shelly (2013), *JGR* 118, 5976-5990, [doi:10.1002/2013JB010118](https://doi.org/10.1002/2013JB010118) | San Andreas LFEs + tidal shear stress | Uses LFE **seismicity rate as a proxy for strain rate** to fit fault rheologies producing creep at 400-600 C. Rheology from rate; no recurrence-predictability test. Relevant as the physical machinery that would convert a recurrence relation into a friction statement. (ABSTRACT-ONLY.) |
| Obara & Kato (2016), *Science* 353, 253-257 | Review | Slow-earthquake recurrence and its relation to megathrusts. No Cascadia size-vs-interval test found. (ABSTRACT-ONLY.) |
| Boso peninsula: Ozawa (2014) *GRL* [doi:10.1002/2014GL060072](https://doi.org/10.1002/2014GL060072); Fukuda (2018) *JGR*; Gardonio et al. (2018) *JGR*; Kobayashi (2019) *EPS* [doi:10.1186/s40623-019-1058-y](https://doi.org/10.1186/s40623-019-1058-y) | 1996-2018 GPS + repeating earthquakes | Recurrence 4-6 yr pre-Tohoku, then 0.6, 2.2, 4.4 yr; the cumulative-moment vs recurrence-time relation is reported as favouring a **slip-predictable** model — the *opposite* pole to Cascadia. Different region, different depth class. (ABSTRACT/summary-level.) |
| Japan renewal-process study (2026), *EPS*, [doi:10.1186/s40623-026-02377-z](https://doi.org/10.1186/s40623-026-02377-z) | Kii, Yaeyama, Boso, Bungo | Mean recurrence 0.96 / 0.51 / 3.97 / 3.16 yr; CoV 0.25-0.54, "SSEs are quasi-periodic, being neither completely deterministic nor highly random". Periodicity, not size-vs-interval. |
| Truttmann, Poulet, Wallace, Herwegh & Veveakis (2024), *GRL* 51, e2023GL107741 | Hikurangi GNSS | Deterministic chaos with short-term predictability. Not size-vs-interval. |
| Nankai/Japan tremor forecasting, *PEPS* (2022), [doi:10.1186/s40645-022-00523-1](https://doi.org/10.1186/s40645-022-00523-1) | Japanese tremor episodes | Renewal-process forecast; 95.3 % / 70.7 % of events inside the 95 % / 68 % intervals. Timing forecast, not a size-vs-interval test. |

**Did anyone report the same sign (size predicts the following interval)?** Yes — Gao, Schmidt & Weldon
(2012), for northern Cascadia, with R2 = 0.16, and they named it "weakly time predictable". Our Spearman
+0.33 is the same sign and the same order of strength.

### Verdict: **C — already published for northern Cascadia, same sign, same magnitude, same name.**

The headline relation is not new. What is new, and what the write-up would have to be built on:

1. **Sample size and spatial extent.** 195 pairs across 12 half-degree bands, 44-50 N, 2010-2026 versus
   **11 events at one GPS station (ALBH), Olympic Peninsula only, 1998-2008**. Gao et al. explicitly
   invited this: "we hope to explore this more as the catalog of events grows."
2. **A tremor-count size proxy** instead of geodetic moment, and a formal p-value where they gave only R2.
3. **A genuine DISAGREEMENT on the other side — this is the interesting part.** Gao et al. found the size
   vs PRE-event interval relation **moderately ANTI-correlated** (R2 = 0.48, r ~ -0.69) and built their
   whole fault-healing interpretation on it. We find it **weakly POSITIVE (+0.16)**. Those are opposite
   signs on the relation they considered the *stronger* of the two. If that reversal survives, the result
   is "the anti-slip-predictable behaviour that motivated the fault-healing interpretation of Gao et al.
   (2012) does not survive an order-of-magnitude larger sample" — which is publishable, and is a much
   better paper than restating the time-predictable half.
4. **The size-proxy discrepancy inside our own study.** Gao et al. used geodetic moment and *did* recover
   the weak positive time-predictable signal. We get +0.33 from tremor count and **nothing** from the
   model-slip proxy. Until that is explained the result is a statement about the tremor catalogue, not
   about the fault.

### To be claimable we need

1. **Relabel to time-predictable** everywhere, and cite Gao et al. (2012) as the prior.
2. **Reproduce Gao et al. (2012) exactly** on the overlapping region and epoch (Olympic Peninsula,
   ALBH-defined intervals, geodetic moment as size) before generalising. Then show what changes with 195
   pairs. The pre-event sign reversal is the deliverable.
3. **Fix the independence problem — this is the biggest threat to p = 0.0005.** 195 pairs from 12
   half-degree bands are not independent: our own R10/S1b concludes "the north is ONE tremor band, not a
   set of local patches", so the 46-50 N bands are near-duplicates of each other. The effective n is far
   below 195. Required: a block bootstrap over episodes (not band-pairs), or a per-band statistic combined
   with Fisher/Stouffer, and a re-quoted p. Expect the significance to fall substantially.
4. **Rule out the episode-definition artefact.** A positive size-vs-next-interval correlation appears
   automatically if the episode detector merges and splits episodes: a merged (hence large) episode
   swallows the following one and lengthens the apparent next interval. Test with a fixed detection
   threshold, with declustering, and with the episode definition varied over a range.
5. **Explain the tremor-count vs model-slip disagreement**, using Wech (2021) MeL energies 2017+ as an
   independent third size proxy.
6. **Control the detection-rate drift**, as for Candidate 1 (PNSN Feb 2018 algorithm change, +55 %
   epicentres): a secular growth in counts imposes a spurious correlation between "size" and calendar
   time, and hence with any interval that also drifts.
7. **Offer the mechanism, don't just report the correlation.** Wech & Creager's (2011) exponential
   recurrence-vs-depth law under a fixed failure threshold *predicts* time-predictability; if the band
   exponents order with depth as they predict, that is a physical result rather than a statistic.

**Venue:** as stated, this is a **section**, not a paper — the headline is Gao et al. (2012). If item 3
holds up and item 2 delivers the pre-event sign reversal, it becomes a short **GRL** or a **BSSA** note
(BSSA being where Gao et al. published).

---

## Bartlow (2020) supporting data — DOWNLOADED

**Success.** `tremor_slip_geometry/results/external/results_for_paper.mat`, 169,711 bytes,
sha256 `093e48cd51bbfc9db16ab8a17b5f98b4ce09259fcde58c4b3a45a42f199294d4` — matches the Mendeley manifest
exactly. Provenance recorded in `results/external/PROVENANCE.md`.

- Dataset: Bartlow, N. (2020), "Time series from PBO and UNR used in Bartlow (2019)", Mendeley Data V1,
  [doi:10.17632/mc49zmg7n7.1](https://doi.org/10.17632/mc49zmg7n7.1), CC-BY-4.0.
- Download URL (via the public API file manifest):
  `https://data.mendeley.com/public-files/datasets/mc49zmg7n7/files/66392424-e44d-4d66-b949-44eaa5b1799a/file_downloaded`
- The other two files were not downloaded (not requested, and large): `PBO_stations.mat` 787 MB,
  `UNR_stations_new.mat` 1.51 GB. Their manifest entries are in PROVENANCE.md if they are ever needed.

### Contents

    Lats, Lons, ID_G                   (443,)      GPS station lat, lon, 4-char code
    Ev_SSE, Nv_SSE, Uv_SSE             (443,)      long-term-average ETS velocity field
    sig_Ev_SSE, sig_Nv_SSE, sig_Uv_SSE (443,)      1-sigma
    el                                 (3607, 3)   mesh triangle connectivity
    slip_rate                          (3607,)     LTA-ETS slip rate; min 0.000, max 0.042, mean 0.004
    strike, rake, ds_fact, ss_fact     (3607,)     per-triangle geometry
    Pole                               (2,)        JdF-NA Euler pole

### The decisive test is still blocked, and here is exactly why

**The mesh node coordinates are not in the deposit.** `slip_rate` is indexed by triangle, and the node
coordinates (`nd`, `nd_ll` in her published script) live in the Green's-function files
`homog_GFs_Nov2019.mat` / `Heter_GFs_May2019`. The repository README states the code "requires Green's
Functions, **not computed specifically for this study**", and those files are neither on Mendeley nor in
the GitHub repo (https://github.com/noelbartlow/Long-Term-Cascadia-ETS-velocities holds only
`Cascadia_secular_fits.m`, `simple_inversion_1comp_2019.m`, `remove_outliers.m`, `README.md`). So
**Bartlow's Oregon slip rate cannot be read off the public data**. Three routes remain: digitise her
published figure; email the author for the mesh/GF file; or rebuild the mesh from her source.

### What the deposit DOES settle — her station distribution, and it helps Candidate 2

443 stations, 38.83-53.91 N, -129.96 to -117.00; **330 west of -120** (the `ind_GPS` cut in her script).
Per-degree counts west of -120:

| band | 39-40 | 40-41 | 41-42 | 42-43 | 43-44 | 44-45 | 45-46 | 46-47 | 47-48 | 48-49 | 49-50 | 50-51 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n | 37 | 40 | 21 | 21 | 21 | 26 | 31 | 38 | 35 | 29 | 15 | 4 |

**Central Oregon 44-46 N = 57 stations; north 46-48 N = 73; ratio 0.78.** A mild deficit, not a hole — and
consistent with our S12 result that the central-Oregon detection threshold is 0.91x the north. Bartlow's
own network does not support a severe modern Oregon coverage gap.

---

## Central-Oregon statements in the two source papers (exact quotes)

### Schmidt & Gao (2010), *JGR* 115, B00A18, [doi:10.1029/2008JB006045](https://doi.org/10.1029/2008JB006045)

PDF obtained from `http://faculty.washington.edu/dasc/PAPERS/2008JB006045.pdf` (13 pp.) and read in full.

- **Abstract, p. 1 (of 13):** "The 16 slip transients have sufficient station coverage to solve for the
  slip distribution on the plate interface. GPS time series are inverted for fault slip using the Extended
  Network Inversion Filter. **Limited station coverage south of Portland (45.5 deg N latitude) restricts our
  analysis to events on the northern half of the subduction zone.**"

- **Paragraph [15], PDF p. 5 — the key passage, and it is TIME-BOUNDED:** "Addition of GPS stations on
  Vancouver Island helps to improve resolution in 2004. Fault plane resolution greatly improves in
  Washington beginning in 2005 with the densification of GPS benchmarks in Puget Sound (47.5 deg N latitude).
  **In Oregon, resolution is poor prior to 2007 from inadequate station coverage.** Therefore, we focus our
  analysis on events that occur north of Portland (>45.5 deg N) because the GPS stations are too widely
  distributed in the south for proper inversions. Likewise, we also limit our analysis to regions south of
  49.5 deg N on Vancouver Island."

  **This materially changes what I wrote in the first pass.** The Oregon caveat is explicitly limited to
  **"prior to 2007"** — before the PBO build-out. Our study window is **2010-2026**, entirely after it. So
  Schmidt & Gao's caveat does *not* bind our epoch, and citing it as a general Oregon-resolution objection
  to R11 would be a misreading of the source. Correct the Candidate 2 section accordingly.

- **PDF p. 7 (2003 event), on how their low resolution actually fails:** "We attributed the patchiness of
  the slip distribution to the fault resolution because **slip tended to be a maximum in the low-resolution
  regions between stations**. Therefore the slip was likely more continuous than visualized."

  Their poor-resolution failure mode produces spurious slip **maxima** between stations, not slip **holes**.
  A persistent slip *hole* is the harder artefact to manufacture — this cuts in favour of R11.

- **No formal resolution test.** The strings "resolution test", "checkerboard", "poorly resolved" and
  "not resolved" do not occur in the paper.

- **PDF p. 11, relevant to Candidate 3's geometry control:** "the dip of the plate interface is shallower
  (~7 deg) [in northern Washington] compared to the dip in **central Oregon (~15 deg)** and Vancouver Island
  (~20 deg)."

### Brudzinski & Allen (2007), *Geology* 35, 907-910, [doi:10.1130/G23740A.1](https://doi.org/10.1130/G23740A.1)

Read in full (4 pp.).

- **p. 908, col. 1:** "Previous observations from Oregon are sparse; a few seismic stations in southern
  Oregon showed a period of nonvolcanic tremor in 2005 (McCausland et al., 2005), and **data obtained using
  one GPS station in central Oregon showed slow slip episodes with a longer recurrence interval** (Szeliga
  et al., 2004)."
- **p. 908, col. 1:** "nonvolcanic tremors that correlate with slow slip episodes are apparent in several
  new locations along the subduction zone, **particularly along central Cascadia**, despite more limited
  observatories (Figs. 2D-2F)."
- **p. 908, col. 1:** "Corresponding seismic and GPS data availability ranges from 1 to 8 yr, with **30
  stations reporting slow slip episodes and 55 stations reporting nonvolcanic tremor** (Fig. 1), of the
  >300 stations that have been investigated with our automated techniques."
- **p. 908, col. 2:** "in the **Siletzia zone** in the central part of Cascadia, data from station **COR in
  Corvallis, Oregon**, extend the history of nonvolcanic tremor back by a factor of two to 1989 (**Fig. 2E**),
  giving a robust estimate of **19 +- 3 months** for the recurrence interval at this station. The average
  interval across the Siletzia zone (**19 +- 4 months**) is longer than those observed on Vancouver Island to
  the north (14 +- 2 months), and is nearly twice as long as that from California to the south (10 +- 2 months)."
- **p. 909, col. 2 / p. 910:** "The central Siletzia zone, with an ~18 month recurrence interval,
  corresponds to the relatively low-lying and young Coastal Range block of central and northern Oregon and
  southern Washington (mostly thick Siletzia terrane) (Trehu et al., 1994)."
- **Fig. 2 panels E and F** are labelled "North Oregon 19 +- 3 mo" and "South Oregon 19 +- 3 mo". The GPS
  station codes listed in the figure are `albh, nano, cpxf, newp, ddsn, ybhb, ntka`. **`newp` is almost
  certainly Newport, Oregon (~44.58 N, -124.06 W)** — inside the 44.3-45.4 N band of R11, and `ddsn` is
  also an Oregon code. *Caveat:* I could not verify the panel-to-station mapping from the PDF text layer
  (the label order in the text stream is not the panel order); the definitive check is **Fig. DR4, the
  station map in GSA Data Repository item 2007224**, which I did not retrieve.
- **Method caveat that limits how much this weighs against R11 (p. 908, col. 1):** their GPS detection is a
  **single-station hyperbolic-tangent fit** with an f-test at 99 % confidence, and they say so explicitly:
  "These methods are **not replacements for network solutions**, but are simply used as surrogates to
  perform a uniform investigation of ETS over the entire subduction zone while network density is still
  heterogeneous." A single-station tanh detection is not a slip inversion and cannot localise slip to
  44.3-45.4 N at 27-38 km depth. So Brudzinski & Allen is weaker evidence against R11 than Bartlow (2020)
  would be — it establishes that *something* moves the Oregon coast on ETS timescales, not that there is
  resolvable slip in our cells.
- No slip-amplitude estimate for Oregon is given anywhere in the paper.

### Net effect on Candidate 2: **strengthened again**

Four things now point the same way: (i) Schmidt & Gao's Oregon caveat is explicitly **pre-2007** and does
not bind 2010-2026; (ii) their poor-resolution failure mode makes spurious slip **maxima**, not holes;
(iii) Bartlow's own station distribution gives central Oregon 0.78x the northern density, matching our S12
threshold ratio of 0.91x; (iv) Brudzinski & Allen's Oregon "slow slip episodes" are single-station tanh
fits that the authors themselves decline to treat as network solutions. The one decisive external test
still outstanding is **Bartlow's slip map itself**, which the public deposit cannot georeference.

---

## Revised ranked publishability (including Candidate 6)

| Rank | Candidate | Verdict | Best venue | Why |
|---|---|---|---|---|
| **1** | C1 — sub-linear tremor-vs-slip-rate scaling (R9) | B | GRL stand-alone, or JGR-SE headline | Unchanged. A number where the field has an assumption, measured two independent ways (count 0.40, energy 0.64), with the smoothing/saturation/network objections already retired by S11. One blocker: a non-vbICA slip model. |
| **2** | C2 — persistent tremor without resolvable slip, central Oregon (R11) | B | GRL, once Bartlow's slip map is obtained | Strengthened twice today. S12 passes the resolution test; the Schmidt & Gao caveat turns out to be pre-2007 and does not apply to our window; low resolution in that class of inversion makes slip maxima, not holes. Remaining: Bartlow's Oregon slip values (blocked by the missing mesh — digitise the figure or email the author), S12b, and drop "tremor-productive". |
| **3** | C4 — 16-yr stationarity (R3/R5/R7) | A (null) | GRL only if aimed at Meng et al. (2026); else a JGR-SE section | Unchanged. Needs Meng's actual quantities and a power analysis. |
| **4** | C3 — within-event deepening (R8) | A for the statistic; C for the intensity half | GRL if it survives; else cut | Unchanged, and now the least controlled of the set. Conflicts with Wech & Creager (2011). |
| **5** | **C6 — size vs following interval** | **C** | a section, not a paper; **BSSA** note or short **GRL** only if the pre-event sign reversal holds | **Gao, Schmidt & Weldon (2012) already published this for northern Cascadia, same sign, R2 = 0.16, and named it "weakly time predictable".** Our contribution is 195 pairs vs their 11, a tremor-count proxy, and — the actually interesting part — a **sign reversal on the pre-event (slip-predictable) relation** they measured at R2 = 0.48 and built their fault-healing interpretation on. Two hard blockers: the label is inverted, and the 195 pairs are not independent (our own R10 says the north is one band), so p = 0.0005 will not survive a block bootstrap intact. |
| **6** | C5 — what Gualandi did | not a result | — | The gap statement for the introduction. |

**Unchanged recommendation:** one JGR-SE paper led by C1, with C2 second, C4 as the stationarity section,
C6 as a short recurrence subsection framed against Gao et al. (2012), and C3 held back pending controls.

**Still could not verify online:** Meng et al. (2026) full text; Costantino et al. (2026) full text;
**Bartlow (2020) full text and the georeferenced slip map** (data downloaded but the mesh is not public);
Hall et al. (2018, 2019); Wech & Bartlow (2014); Aguiar et al. (2009); Michel et al. (2019a);
Hawthorne & Bartlow (2018); Beeler et al. (2013); Obara & Kato (2016); Brudzinski & Allen GSA Data
Repository item 2007224 (Fig. DR4 station map).
