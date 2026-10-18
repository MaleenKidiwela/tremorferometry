# The published Cascadia LFE catalogues — what each one is, and whether we are reproducing it

Written 2026-09-24 by agent Z, from the full-text papers on disk; **§4 (Sweet et al. 2014) and §10 (Shelly et al. 2025) added 2026-09-25 by agent AA** from the two PDFs that arrived after agent Z's pass — sections 4–10 are therefore renumbered from the original 4–9 (not from abstracts or secondary
summaries) plus our own results files. Companion file: `notes/PAPER_METHODS_TABLE.md`, which holds
the verbatim parameter-by-parameter table and the mismatch list.

**Gloss of terms used throughout.**
*LFE* = low-frequency earthquake, a very small earthquake (magnitude ~1–2) that repeats in the same
spot thousands of times and is depleted in high frequencies. *Family* = one such spot; all the
events in a family have the same waveform. *Template* = the average waveform of a family, built by
stacking many detections. *Network autocorrelation* = sliding every short window of a day's data
against every other window and summing across stations, to find pairs that look identical — this is
how families are *discovered* with no prior knowledge. *Matched filter* / *network
cross-correlation* = sliding a known template through the data to find more copies — this is how a
family is *extended* once you have it. *MAD* = median absolute deviation, a robust measure of
scatter; "8×MAD" means a detection is logged when the correlation exceeds 8 of those units above
the median. *Declustering* = keeping only the strongest detection in a short time interval so one
event is not logged twice.

Ordered geographically, north to south.

---

## 1. Royer & Bostock 2014 — northern Vancouver Island (NVI)

**What they found.** 54 initial templates, reduced to a final suite of **47** ("Our final suites of
templates number 122, 54 and 47 for the NW, SW and NVI datasets", §3.3), from **4267 initial
detections** (§3.1). Data span two short tremor episodes only: **5–11 Sept 2006** and **13–18 June
2007**, "each episode lasting under 1 week" (§2.2). Depths: "hypocentres fall between **33 and 39
km**" (§5). Each template has "100's to 1000's of independent detections" (§3.3).
**They located them — yes, all of them**, first with **Hyp2000**, then relocated with **hypoDD**;
"**91 %** … of the original events remain in hypoDD solutions for the NVI … dataset" (§5). Hyp2000
nominal errors 1.8 km horizontal, 2.16 km vertical.

**Method in one paragraph.** A 27-station portable POLARIS deployment on two perpendicular arms
(§2.2); network autocorrelation used "**different combinations of 7 of the 10 three-component
stations VI10, VW03, VI11, VW02, VI52, VI53, VI05, VI08, VI06, VI04**" (§3.1). All data cut into
24-hour files, **1–8 Hz, 40 samples per second** (§3); inside each hour, **15-second windows lagged
by 0.5 s** (7170 windows per hour) correlated against every other window on each station-component,
summed across all stations and components; detection when the sum exceeds **8 × MAD** of that hour
(§3.1); then the detections were "culled to exclude overlapping detections and retain only those
with high SNR". Detections were turned into families by joining every channel into one "**super**
trace", cross-correlating all such traces, and running a **hierarchical cluster analysis** on the
resulting similarity matrix (§3.2). Templates were then improved by iterating network
cross-correlation with "**subsets of 7 to 10 high SNR, three-component stations**", now thresholded
at **8 × MAD over a 24-hour period** (note: the MAD window changes from 1 hour to 24 hours between
the two stages), restacking after each pass, for "several iterations" (§3.3).

**What they released.** Waveforms via IRIS/GSC. **No catalogue.** The supplementary material is
one table of hypoDD settings and four figures — no family list, no detection times. SILENT.

**Are we using their method? NO.** Our chain uses **6-second windows** (`lfe_array_discovery/
scripts/k_auto.py:42`) where they used 15; **20 samples per second** (`k_common.py:28`) where they
used 40; and **connected-component graph grouping** (`k_group.py:51-65,96-99`) where they used
hierarchical clustering. We have also never run the NVI stations.

**Status.** **Not attempted.** No NVI run exists.

**To reproduce from here.** Download the 10 named POLARIS NVI stations for the two 2006/2007 tremor
weeks; run at 15 s / 40 sps / 1–8 Hz; score against 47 final templates at 33–39 km depth.

---

## 2. Bostock et al. 2012 — southern Vancouver Island (SVI) — *the catalogue we reproduced*

**What they found.** "An initial suite of **5775 LFEs** detected in 2004 and 2005 at a select set of
7 stations is grouped into **140 families**" (Abstract; §2.1 ¶6). 46 tremor days analysed. "Most
templates comprise between **1000–5000 independent detections** made over the 8 ETS episodes"
(§2.2 ¶7). "The LFEs define a relatively tight surface **between 27 and 37 km depth**" (§3 ¶9).
**They located them — yes, all 140**, with **hypoinverse** [Klein 2002], from "averages of
approximately **16 P-times and 22 S-times** per template" (§3 ¶8). They deliberately did *not*
deduplicate: "Although some of the 140 composite waveform templates exhibit overlap in contributing
detections, we will consider each template as representing a distinct LFE location" (§3 ¶8).

**Method in one paragraph — note how it differs from §1 above, same authors, one year earlier.**
Discovery used the 2004 (July 4–26) and 2005 (Sept 3–25) ETS episodes and the POLARIS/GSC "anchor"
line "**TWBB, TSJB, TWKB, MGCB, KELB**, and permanent stations **PFB, LZB and PGC**" with "5–10 km
station spacing" (§2 ¶4) — the paper lists eight codes but calls them "the 7 three-component anchor
stations" (§2.1 ¶5), an inconsistency in the paper itself. Band **1–10 Hz** (not 1–8), **40 sps**,
24-hour files. Windows were **6 seconds with 0.5 s overlap**, 7188 per hour, each correlated against
every other window in the same hour; summed over stations and components; detection at **8 × MAD**
(§2.1 ¶5). Grouping was "a **hierarchical cluster analysis** [Rowe et al. 2002] based on correlation
coefficients of composite wave trains" (§2.1 ¶6). Extension used "subset templates comprising **6 or
7 stations**", again at 8 × MAD, alternating detection and restacking; "**After 5 such iterations**
neither template SNR nor numbers of detections continued to show significant improvement" (§2.2 ¶7).
No declustering rule is published.

**What they released.** Waveform access acknowledgement only. **No catalogue in the paper** —
SILENT. (The 130-family successor set from Bostock et al. 2015 is what our project works from; that
paper is also SILENT on release.)

**Are we using their method? PARTLY.**
- window 6 s vs their 6 s — **matches** (`k_auto.py:42`)
- threshold 8 × MAD vs their 8 × MAD — **matches** (`k_auto.py:45`)
- band **1–8 Hz vs their 1–10 Hz** (`d_equiv.py:5-8`)
- sample rate **20 sps vs their 40 sps** (`k_common.py:28`)
- stations: our 11-station `CLEAN11` = GLBC, TWKB, MGCB, TSJB, TWBB, PGC, SILB, SSIB, LCBC, B003,
  B005 (`k_common.py:49`), and the queues run 15 (`k_repro3_queue.sh:13`) — **only 5 of their 8
  anchors are in it** (TWKB, MGCB, TSJB, TWBB, PGC; KELB, PFB and LZB are missing from CLEAN11,
  LZB is added in the 15-station set)
- grouping **connected components vs their hierarchical clustering** (`k_group.py:51-65`)
- iterations **3 vs their 5** (`k_grow.py:234`)
- declustering **4 s, which we imported from Bostock et al. 2015** (`k_auto.py:46`) where this paper
  publishes no rule at all

**Status. REPRODUCED — the project's anchor result.** Blind on the 2005 array bed with 15 stations:
**58 of 96** open ("kept") families and **22 of 30** sealed hold-out families rediscovered, against
ceilings of 89 and 29 respectively (`lfe_array_discovery/results/REPRO_BOSTOCK.md`,
`NORTH_DISCOVERY_SUMMARY.md`). The same chain on permanent-station beds discovers only 3/77 (2012)
and 5/60 (2010) — dense arrays find families, permanent networks only run them.

**To reproduce more exactly from here.** Move the cache to 40 sps and 1–10 Hz, add KELB/PFB, swap
connected components for hierarchical clustering, and raise iterations to 5.

---

## 3. Chestler & Creager 2017 — Olympic Peninsula

**What they found.** "Initially, we found **57 LFE families** beneath the Olympic Peninsula …, each
with **250–1700 detections** … There are **34,479 LFEs** in total" (§3.1); after removing duplicate
families, "our final catalog includes … **32,387 unique LFEs sorted into 43 families**" (§3.1).
Moments span "M0 ≈ 1.4 × 10^10 to 1.9 × 10^12 N m (**Mw = 0.7–2.1**)" (§3.2). Time span: the Array
of Arrays deployment, June 2009 – Sept 2010 and Aug 2011.
**They located them — yes, all 43**: "we perform a **grid search using the 3D velocity model from
Preston et al. [2003]** to locate each LFE family", then "improve relative locations between nearby
families using the **double difference method of Waldhauser and Ellsworth [2000]**" (§2.1). The
paper does not print a depth range in the text (cross-section only) — SILENT.

**Method in one paragraph — a completely different recipe from the two above.** They did not
autocorrelate raw station records. They started from Ghosh et al.'s array **beams**: eight 1-km
arrays of 10–20 three-component stations, **vertical component only, 5–9 Hz**, 1-minute windows
with 50 % overlap, delay-and-sum beamforming over a slowness grid (§2.1). On each beam they took
peaks "**larger than 18 times the median absolute deviation**", spread each peak into a 2-second
boxcar, stacked the boxcars across arrays, and used the peaks of that as candidate LFEs. Only then
did they autocorrelate: "**15 s windows** of three-component data, band-passed between **2 and 8
Hz**", using **three stations per array**, stacked "using an **nth-root stack where n = 2**", with
a detection bar of "**45 times the MAD**" (§2.1) — a far higher bar on a different statistic than
the 8 × MAD used by Bostock and Royer. Templates were built from "**up to 150** of the time windows
that correlate best", refined over "**three iterations**", and family membership used a correlation
threshold that slides with how many stations were live: **0.05 / 0.1 / 0.14 / 0.15 / 0.2** for <6,
7–9, 10–12, 13–15, >15 stations. Duplicate families were removed when they shared more than half
their detections ("When two LFE families have **Rcom > 0.5**, the family with the lower median
correlation coefficient is removed", §2.2).

**What they released.** "All the waveforms from the Array of Arrays and CAFE experiments can be
found at the IRIS Data Management Center" (Acknowledgments). **Catalogue: SILENT** — only Table S1
in the supporting information (characteristic moment per family). No detection times, no
coordinates, in the paper.

**Are we using their method? NO.** We do build array beams (`lfe_array_discovery/scripts/
o_beam.py:49-52`, 1-minute windows, 50 % overlap, slowness grid — faithful), but every step after
that is the Bostock recipe, not Chestler's:
- **no 18 × MAD beam-peak stage and no boxcar cross-array stack at all**
- autocorrelation window **6 s vs their 15 s** (`o_auto.py` imports `k_auto.WIN_S`, `k_auto.py:42`)
- band **1–8 Hz vs their 2–8 Hz** for the autocorrelation, and **5–8 vs their 5–9 Hz** for the beam
  (a documented cache limitation, `o_beam.py:18-20`)
- **linear stack at 8 × MAD vs their 2nd-root stack at 45 × MAD** (`k_auto.py:45`)
- template stack **200 vs their 150** (`k_grow.py:39`); iterations 3 vs their 3 — matches
- **no station-count-dependent family cc ladder, and no Rcom > 0.5 dedupe** in the running chain
- location: **Crosson 1-D grid search with VP/VS = 1.78** (`uu_common.py:177,220`) plus our own
  damped least-squares double-difference (`uu_dd.py:30-39`) vs their **Preston 3-D model with
  VP/VS = √3 ≈ 1.732** and Waldhauser–Ellsworth

**Status. A SPATIAL MATCH ONLY, and the phrasing matters.** Our 2010 Olympic arrays-as-arrays run
produced 1,479 candidate families of which 788 located; **61 of 69** published Olympic coordinates
have one of ours within 10 km, against **3 of 69** for the whole-chain site-shift null
(`lfe_array_discovery/results/WA_LOCATION_REBUILD.md` §9). Two honest caveats: (a) those 69
coordinates are *not* Chestler's 43 — they are a 69-row coordinate file recovered from a public
repository whose provenance we could not confirm, and our own file name says so
(`codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv`);
(b) "61 of 69 coordinates have a neighbour within 10 km among our 788" is a *coverage* statement,
not 61 identified families. An earlier, different comparison on the same 2010 bed was outright
negative: our candidates sat *farther* from those coordinates than random draws
(`WA_2010_RESULTS.md` §5d, "fraction of random draws at least as close: **1.000**").

**To reproduce from here.** Implement the missing front end — 18 × MAD beam peaks, boxcar
cross-array stack, 15-s 2–8 Hz autocorrelation, 2nd-root stack, 45 × MAD — then the Rcom > 0.5
dedupe, and score against 57→43 families and 32,387 LFEs.

---

## 4. Sweet, Creager & Houston 2014 — ONE family at the downdip edge of tremor, western Washington

**What they found.** A single family, worked in unprecedented detail: "We analyze an **isolated** LFE
family located at the **downdip edge of the main episodic tremor and slip (ETS) zone beneath western
Washington State**. The **9000** individual LFEs from this repeating family cluster into **198
swarms** that recur roughly every week" (Abstract). After three named culls the catalogue is
**8942 events, October 2006 – September 2011** (§2). It is "the **deepest (most downdip) and most
frequently active** family we have yet found", on "the far eastern edge of the region of PNSN tremor
detections (2009–2013) … and at the downdip edge of 122 LFE families located in northern Washington
[Royer and Bostock, 2013]. … **Their nearest LFE family is 10 km south of ours, and their next
nearest is 20 km west**" (§1) — and "we find **no other detectable LFEs within 10 km** of this LFE
family" (§6), which is what makes it a clean laboratory.
**They located it — yes**, 700 of its events individually: "we employed a **double-difference
relocation method** [Waldhauser and Ellsworth, 2000] to find precise relative locations for the
**700 most highly correlating events** … nearly **4 million** autocorrelation-obtained differential
times (**1.4 million for P waves, 2.6 million for S waves**) … a **variance reduction of 73 %**"
(§4), with "formal errors of **100 m**". The cloud is "approximately **2 km long and 500 m wide**,
elongated parallel to the relative plate convergence direction", "**<300 m** vertical spread", with
a core "roughly circular with a radius of **300 m**" holding 46 % of the events. Absolute depth in
kilometres is **SILENT** — depth is reported only relative to interface models: "LFE depths are
**within 1 km of two different plate interface models** [Preston et al., 2003; McCrory et al., 2012]
and exhibit a similar dip angle to the east", while "The plate interface model of Audet et al. [2010]
is about **5 km shallower** than the LFEs" (§4). Latitude/longitude are likewise given only in
Figure 1 (a map with BS, W020, B013, FACB).

**Method in one paragraph — and this is the fuller statement of Chestler & Creager's front end.**
Same group, same CAFE data, three years earlier. "We first identified this LFE through **visual
inspection** of the data. We chose **one of the highest signal-to-noise instances** of this LFE and
used that reference event to build up a template. A **15 s window containing the P and S waves** for
the reference event was **autocorrelated with 12 h of data for all available channels from the 1 km
aperture Big Skidder (BS) array of six 3-component stations** … The resulting **fully normalized
autocorrelation functions—one for each channel—were then combined using a third-root stack** to
identify time windows that correlated best **across all 18 channels**. **Templates were built by
linearly stacking the 80 best correlating windows, forcing identical time offsets for each
station/channel. This process was repeated three times.** A complete catalog was then produced using
a stack of fully normalized **running autocorrelations** between each stacked template and its
**continuous 5 year long seismogram**. **For each day, we calculate the median of the absolute
deviation (MAD)** … a detection is a time when the stacked autocorrelation function **exceeded
10\*MAD**, producing a typical threshold in which the **mean correlation typically exceeded
0.2–0.3**. Because the **S minus P time is about 6 s**, and there are often spurious detections when
the P wave of the template aligns with the S wave of an LFE …, we **require the time between adjacent
LFE detections to be greater than 6 s**" (§2). Then three culls, each with its count: −60 within 30 s
of a PNSN catalogue earthquake, −(30 logging days in 2007), −180 isolated events visually judged
non-LFE → **9600 → 8942**. Amplitudes are measured on horizontals with the **response deconvolved to
displacement, 2–8 Hz** (§5).

**What they released.** Supporting information is "**Readme, Figures S1–S10**" only. There is **no
data availability section and no catalogue release — SILENT**. The BS array is CAFE, whose waveforms
are at the IRIS DMC.

**Are we using their method? NO — and this is Olympic mismatch M4b.** Our Olympic/array chain runs
**6 s** windows where they use **15 s**, a **linear** stack where they use a **third-root** stack
(Chestler uses a 2nd-root — the two papers differ from each other), **8×MAD** where they use
**10×MAD**, **200** stacked windows where they use **80** with forced identical offsets, and a **4 s**
declustering radius where their **6 s** minimum separation is *derived from the S−P time* to kill the
P-on-S alias (`lfe_array_discovery/scripts/o_auto.py:237`, `k_auto.py:42,46`, `k_grow.py:39`).

**Status. Not attempted as a reproduction.** This family sits at the downdip edge east of the
Olympics, inside the region our Chestler run covers spatially, but we have never targeted it. It is
the cheapest single-family positive control on the margin: one family, one 1-km array, a five-year
continuous record, and four falsifiable behaviours that need no catalogue file at all.

**To reproduce from here.** Download the six BS stations (CAFE) for 2006-10 → 2011-09 plus W020,
B013, FACB; run 15 s / third-root / 80-window / 3-iteration templates at 10×MAD daily with a 6 s
minimum separation; then check the four published behaviours: **198 swarms**, **8-day** recurrence,
swarm-moment-vs-time-to-next-swarm **r = 0.52 ± 0.03** (p 10⁻⁹–10⁻¹³, against 0.16 ± 0.03 for time
*since*), in-swarm migration at **1 km/h** (5 updip, 2 downdip, 2 none of 8 swarms), amplitudes
**0.2–18 nm** → **Mw 0.3–1.5** with an **exponential** distribution (A₀ ≈ 4.5 nm), total moment
**2.6 × 10¹⁴ N m (Mw 3.6)**, LFE locking efficiency **≤ 20 %** concentrated in two 50-m-radius spots,
stress drops **1–20 kPa**, and LFE activity **1–3 days after** vigorous updip tremor in each of the
five 2007–2011 ETS episodes.

---

## 5. Royer & Bostock 2014 — northwest Washington (NW), eastern flank of the Olympics

**What they found.** "**4915** … initial detections" → "**224** … initial templates" → a final suite
of "**122**" templates (§3.1–3.3). Data from major tremor episodes in 2007, 2008, 2009, 2010 and
2011, autocorrelated on four of them: **17–30 Jan 2007, 1–23 May 2008, 7–24 Aug 2010, 4–23 Aug
2011** (§3.1). Depths: "those below Washington state fall between **31 and 46 km** depth (NW)" (§5).
**They located them — yes.** Hyp2000 for all, then hypoDD, which retained "**77 %** of the original
events" for NW (§5). Hyp2000 nominal errors **2.0 km horizontal, 1.3 km vertical**.

**Method in one paragraph.** The station pool was "**25 stations** skirting the eastern flanks of
the Olympic Peninsula to the west of Puget Sound" (§2.1), mostly from the CAFE Flexible Array.
Network autocorrelation used exactly **seven named three-component stations: B001, SQM, BS11, W020,
W040, GNW and PL11** (§3.1). Everything else is as in §1 above: 24-hour files, **1–8 Hz, 40 sps**,
**15-second windows lagged 0.5 s**, **8 × MAD over the hour**, cull overlapping and low-SNR
detections, hierarchical clustering of super-trace similarity into 224 initial templates, then
iterated network cross-correlation with 7–10 high-SNR 3-C stations at **8 × MAD over 24 hours**,
converging to 122.

**What they released.** **Nothing.** No catalogue, no coordinates, no detection times — SILENT.
Their *positive claim* (families exist here, on the interface, at 31–46 km) is nonetheless testable
on our side by existence.

**Are we using their method? NO — and this is the failure the project already hit.**
- window **6 s vs their 15 s** (`k_auto.py:42`)
- sample rate **20 sps vs their 40 sps** (`k_common.py:28`)
- stations: **52 probe-admitted stations** chosen by our own quiet-day contrast probe
  (`lfe_array_discovery/results/site_probe_2008.csv`, selection at `n_common.py:171-208`) **vs
  their 7 named**. Worse, **three of their seven are not in our set at all** — SQM, GNW and W020
  are absent (GNW was explicitly rejected by our probe, `n_report.py:200`); only B001, BS11, PL11
  and W040 survive.
- grouping **connected components vs their hierarchical clustering** (`k_group.py:51-65`)
- location **Crosson 1-D + our own LSQR double-difference vs their Hyp2000 + hypoDD**
  (`uu_locate.py:37-41`, `uu_dd.py:30-39`); their hypoDD settings are published in Table S1 (IDAT 3,
  five blocks of 8 iterations = 40, DAMP 50/50/40/40/40) and we have never used them.

**Status. FAILED — and the failure is now attributed to the method, not to Washington.** The
51-station Puget/S-Olympic 2008 run found nothing usable: the pair stage fired as often on a quiet
day as inside the ETS (tremor/quiet ratio **0.89**, `WA_2008Q_RESULTS.md` intro) against 11× on the
Vancouver Island bed where the same code reproduces Bostock. The project plan already records this
as "BUG REPORTS against our code, not evidence about Washington"
(`notes/MARGIN_CATALOG_PLAN_2026-09-24.md`). The present document identifies the bug concretely:
**the run used Bostock's 6-second Vancouver Island window and our own 52-station probe set, not
Royer & Bostock's 15-second window and their seven named stations.**

**To reproduce from here.** Re-cache the four NW tremor episodes at 40 sps, run at 15 s / 1–8 Hz on
exactly B001, SQM, BS11, W020, W040, GNW, PL11 (downloading SQM, GNW, W020, which we do not have),
group by hierarchical clustering, then extend with 7–10-station subsets at 24-hour MAD. Reference
lines: 4915 detections, 224 initial templates, 122 final, depths 31–46 km.

---

## 6. Royer & Bostock 2014 — southwest Washington (SW), southwest of Puget Sound

**What they found.** "**3306** … initial detections" → "**70** … initial templates" → "**54**"
final (§3.1–3.3). Only two tremor episodes were available: "**14 Jan – 1 Feb 2007**" and "**1–25
May 2008**" (§2.1, §3.1). Depths: "**32 and 49 km** depth (SW)" (§5).
**They located them — yes**, Hyp2000 then hypoDD, which retained "**72 %**" (§5). Hyp2000 errors
here are the worst of the four subarrays: **3.7 km horizontal, 3.5 km vertical**.

**Method in one paragraph.** Identical to §4 except the station pool ("**25 stations** running
approximately east–west to the southwest of Puget Sound", §2.1) and the autocorrelation set: "the
**7 three-component stations N050, S030, W070, N060, PL11, S040 and N070**" (§3.1). Note PL11 is the
only station shared with the NW set — the two Washington subarrays were processed independently.

**What they released.** Nothing — SILENT.

**Are we using their method? NO** — same five mismatches as §4. Of their seven SW stations, **six
are in our admitted set** (N050, S030, W070, N060, PL11, S040) and **N070 is not**. So SW is closer
to recoverable than NW on stations, but still runs at 6 s / 20 sps / connected components.

**Status. Not attempted separately.** The 2008 Puget run overlaps this footprint and found nothing
(see §4).

**To reproduce from here.** Same as §4 with the SW seven and the two 2007/2008 episodes; reference
lines 3306 detections, 70 initial templates, 54 final, depths 32–49 km.

---

## 7. Thomas & Bostock 2015 — central Oregon (~44°N)

**What they found.** "a total of **18 families** and over **15,000** individual LFE detections"
(§3; the abstract says "nearly 16 thousand events"). Scan span: "**46 days of continuous data from
August 3 through September 17, 2009**" (§2).
**They located them — only partly.** "Of these 18, we were unable to locate **five** families
because the templates were not of sufficient quality to accurately pick arrivals" (§3), leaving
**13 locatable** (Fig. 3 caption). Of those 13, "most locations have poor depth control so we **fix
their depths to the plate interface**"; only **four** have free depths, "ranging from **31 to 34
km**" with the Williams et al. (2011) model, "**increase[d] by an average of 4 km**" with the PNSN
O0 model (§3). Software: **HYPOINVERSE**.
Falsifiable extra: the slip front reached the northernmost families "on **August 30th** and
propagates south at a rate of **5 km/day**" (§3).

**Method in one paragraph — and a published negative result we must respect.** They first tried
exactly the Bostock/Royer approach and it failed outright: "We applied network autocorrelation to 3
days of continuous data (09/03–09/05) … the network autocorrelation method **failed to identify a
single LFE template**. This is likely due to the station spacing (approximately 30 km on average)"
(§2). So they switched to a **cross-station** method (Savard & Bostock 2014) applied to just
**three** stations, "**ALP0, ALS0, and BWO0** … due to their generally high signal-to-noise ratios
and proximity to the 2009 SSE" (§2). Data were **1–10 Hz, 40 sps, 24-hour files**; network
cross-correlation at **8 × MAD**; **phase-weighted stacking on the intermediate iterations** (a
stacking scheme that up-weights parts of the stack that agree in phase — neither Bostock 2012 nor
Royer 2014 uses it), iterating "until no additional detections are registered (generally **10–15
iterations**)", with stations added as coherent energy appeared. A template needed "a **minimum of
100 detections in the 3-day window**". Families were deduplicated so that "**no family shares more
than 10 % of detection times**", matched within **5 s**. The autocorrelation/cross-correlation
window length is **SILENT** — they say only that they "experimented with window durations".

**What they released.** IRIS waveform acknowledgement only. **Catalogue: SILENT.**

**Are we using their method? NO.** We have no cross-station detector at all; our chain is the
network-autocorrelation one that this paper reports **does not work in Oregon**. Specific
differences: **6 s window vs SILENT** (`k_auto.py:42`); **20 sps vs their 40** (`k_common.py:28`);
**1–8 Hz vs their 1–10 Hz** (`d_equiv.py:5-8`); **linear stacking only, no phase-weighted stacking**
(`k_grow.py:50-136`); **3 iterations vs their 10–15** (`k_grow.py:234`); **no 10 %/5 s dedupe** in
the running chain; **no ≥100-detection template gate**.

**Status. Queued, not started** (`notes/MARGIN_CATALOG_PLAN_2026-09-24.md` row 6).

**To reproduce from here.** Implement the Savard & Bostock cross-station detector, seed from ALP0,
ALS0, BWO0 over 2009-08-03→09-17, add phase-weighted stacking and the 10–15-iteration convergence
rule. Reference lines: 18 families, >15,000 detections, 13 locatable, 4 free depths at 31–34 km,
and the 5 km/day southward front.

---

## 8. Plourde et al. 2015 — southern Cascadia (40–41.8°N) — *the seed set behind Ducellier*

**What they found.** "A total of **59 initial templates** from Cascadia produced waveforms
sufficiently clear to be picked. These 59 templates were found to represent **34 distinct LFE
families**. These families comprised **9900 LFEs**, or a mean of **290 per family**, with **870** as
the largest number of detections in a single family" (§3). Time span: **21 March – 30 April 2008**.
"Cascadia LFEs extend from **40°N to 41.8°N**" and "Hypocentral depths range from **28 km to 47
km**, a larger distribution than in northern Cascadia" (§3). Three further families sit far to the
south on strike-slip faults: two on the Maacama Fault at **15 km** and **21 km** depth and one on
the northern Bucknell Creek Fault at **16 km**, with 869, 302 and 325 repeats (§4).
Falsifiable extras: "The main front of tremor epicenters migrates along strike at **9 km d⁻¹**"
and one rapid tremor reversal "at **10–20 km h⁻¹**" (Abstract).
**They located them — yes, all 34** (every family is labelled with its depth in Fig. 1), using "a
**one-dimensional velocity model modified from Verdonck and Zandt [1994]**" (§2). **The location
software is SILENT.**

**Method in one paragraph.** **25 stations from the EarthScope Mendocino Experiment plus 7
permanent Northern California Seismic Network stations** — the individual codes are SILENT (§2).
Band **1.5–9 Hz**, sample interval **Δt = 0.025 s** (i.e. 40 sps). Crucially there is **no
autocorrelation stage to copy**: candidates came from "a combination of **autodetection methods**
[Savard and Bostock, 2015] **and visual identification**", with "windows of **60 s** around these
LFEs … used as initial templates" (§2). Templates were then improved by iterated network
cross-correlation and stacking, "**phase-weighted stacking … until the final iteration**" and linear
after that; per-station taper windows let them work with "a **short window of 12 s**" (§2).
**The detection threshold is SILENT — this paper never states a MAD multiple.** Uniqueness used two
independent tests: templates are the same family if they share "**>10 % of detections**" (allowing
25 s of lag) **or** if their cross-correlation exceeds **0.3** (§2).

**What they released.** Waveform DOI 10.7914/SN/XQ_2007 only; **catalogue SILENT in the paper** —
but Plourde privately shared it with Ducellier, and that shared copy (66 template families) is on
our disk inside `catalogs/ducellier_catalog/data/Plourde_2015/templates_list.txt`, with latitudes
39.45–41.72 N, depths 15.0–46.8 km and published uncertainties (`lfe_south/results/
SOUTH_TRUTH_INVENTORY.md` §2).

**Are we using their method? NO** — we have no cross-station/visual candidate stage, and this paper
supplies no threshold to copy. We do use their **12 s template length** (`k_grow.py:33`, cited to
them) and, in the south, their **1.5–9 Hz band** (`lfe_south/scripts/t6_duc_cache.py:25`).

**Status.** Reproduced **indirectly**, through Ducellier — see §9. The 34/66 Plourde families are
the objects our southern reproduction recovers.

**To reproduce directly from here.** Would need the Savard & Bostock cross-station detector plus a
threshold we would have to choose ourselves (the paper gives none), so a *direct* reproduction is
not well posed. Scoring against 59→34 families, 9900 LFEs and 28–47 km depths is.

---

## 9. Ducellier & Creager 2022 — southern Cascadia, eight years — *our best external validation*

**What they found.** They extended Plourde's April-2008 snapshot to **2004–2011**. They worked from
Plourde's **unreduced 66** families, not his published 34: "For each of **66 LFE families** they
provided hypocentral locations, lists of stations and channels used to detect LFEs, and the timing
of all LFE detections. [Plourde et al.] have later reduced the number of LFE families to 37 …, but
we chose to use the initial detections". Six were later dropped ("We eliminate **six families** for
which too few permanent stations have good templates"), leaving **60** carried to the end. On disk
the released catalogues hold **477,715 detections (2004–2011, permanent networks)** and **107,924
(2007–2009, FAME array)** (`lfe_south/results/SOUTH_TRUTH_INVENTORY.md` §2).
**They did not locate anything new** — "We use the same locations as P[l]ourde et al. (2015). We
have not tried to relocate the families using the new templates, except for family A" (reply to
reviewers). Family A alone was retried and may sit "**30 to 50 kilometers east**" at "**20-25
kilometers**" depth, shallower than the plate boundary, possibly on a crustal fault.

**Method in one paragraph — a pure matched-filter paper, no discovery at all.** Stations: FAME
(temporary) plus three permanent networks (BK, NC, PB), one- *and* three-component, channels
BHZ/EHZ/HHZ/SHZ, "**less than 100 kilometers away from the epicenter** of an LFE family". Waveforms
were detrended, Hann-tapered, **instrument response removed** (the only one of these papers to say
so for detection), **band-passed 1.5–9 Hz**, and **resampled to 20 Hz**. Templates: an 80-s cut
trimmed to **60 s**, RMS-normalised and linearly stacked; later, for the permanent networks, rebuilt
from "the **150** LFE detection times with the best cross-correlation value" and shortened to
"**25 to 40 seconds** … average duration equal to **30 seconds**". Detection: per-channel normalised
cross-correlation, **averaged over channels**, thresholded at "**eight times the median absolute
deviation (MAD) of the stacked cross-correlation**", computed **per one-hour block** ("The threshold
is thus different for each hour of data", reply to reviewers). Declustering keeps "only LFEs that
are separated by at least **one second**". A second, family-specific filter then requires
**cross-correlation × number of live channels** to exceed a per-family bar, calibrated so half the
survivors are confirmed by the other catalogue; this removed "about **84 %**" of the FAME detections
and "about **63 %**" of the network ones.

**What they released — everything.** Waveform DOIs, and "The Python scripts used to analyze the
data and make the figures can be found on the first author's Github account, accessible through
Zenodo". That repository is on our disk at `catalogs/ducellier_catalog/`: his detector source, his
station tables, his **66 families of templates**, and both detection catalogues. It is the only one
of these nine with released code *and* catalogue *and* templates.

**Are we using their method? YES, almost exactly — one parameter differs.**
- band 1.5–9 Hz, 20 Hz, response removed — **matches** (`lfe_south/scripts/t6_duc_cache.py:25,81-82`)
- his own templates, his own station list — **matches** (`t7_duc_mf.py`)
- per-channel normalised cc averaged over channels with a complete gapless record — **matches**
  (`t7_duc_mf.py:113-125`)
- 8 × MAD — **matches** (`t7_duc_mf.py:26`)
- 1-second declustering, his `clean_LFEs` run-chaining — **matches** (`t7_duc_mf.py:44-52`)
- **MAD block: ours is 12 hours, his is 1 hour** (`t7_duc_mf.py:26`, `BLOCK_H = 12`). His paper says
  "we download **one hour** of seismic data" and his code loops hour by hour
  (`catalogs/ducellier_catalog/src/find_all_LFEs_parallel.py:220,227` and the MAD at line 321).
  The 12-hour figure in his script is only how he chunked the job, not the statistic's window. A
  12-hour MAD sets a different bar than an hourly one, and our stricter "gapless for the whole
  block" channel rule drops a channel for 12 hours where his loses only the affected hour.
- his second, per-family `cc × n_channels` filter is **not implemented**

**Status. REPRODUCED.** Scored one-to-one at ±2 s against his published detections over
2008-04-01…05-10, on his stations with his templates: **mean recall 0.836** (median 0.878),
precision 0.708, pooled recall 0.873; **63 of 66 families** reach Bostock's ≥30 % reference line;
the loud decile is recovered essentially completely (mean recall **0.988**)
(`lfe_south/results/SOUTH_RESULTS.md` §0–1). Adding 25 modern stations he did not use made it
**worse** (recall 0.621, precision 0.400) because averaging more channels dilutes the network mean
faster than it lowers the 8 × MAD bar — documented with the arithmetic in the same file.

**To close the last gap from here.** Set `BLOCK_H = 1` in `lfe_south/scripts/t7_duc_mf.py:26`,
re-derive the gapless-channel rule per hour, and re-score; then add his per-family
`cc × n_channels` filter to compare precision on equal terms.

---

## 10. Shelly et al. 2025 — Mendocino Triple Junction (~40.0°N)

**Read 2026-09-25 (agent AA).** ⚠ One caveat travels with every parameter below: the paper says
"The methods are **briefly summarized here; full methods are available in Supporting Information
S1**" (§2), and **S1 is not on disk**. Everything here comes from the main text, the figure captions
or the acknowledgements; anything else is SILENT and must not be filled in from a neighbouring paper.

**What they found.** "we examine its constituent low-frequency earthquakes (LFEs) by developing **27
stacked LFE waveform templates** that we use to detect events from **2018 to 2024**" (Abstract). The
target is not the main Cascadia band at all: it is "an anomalous zone in southern Cascadia, **50–100
km west of the main tremor band**, near the expected **southern edge of the subducting Gorda slab at
the Mendocino triple junction**" (Abstract), the zone Wech (2021) found. The family count is honest
about being a choice: "**Although largely arbitrary, this number was chosen as sufficient to
collectively detect activity during the 2018–2024 episodes** of tremor identified in this area by
PNSN monitoring" (§2).
**They located them — all 27, twice.** "we determined initial locations for each family using a
**grid search for minimum mean residual using the layered P-wave velocity model from Yoon and Shelly
(2024)**, following the approach of Shelly and Hardebeck (2010). We refined these locations with
**hypoDD** (Waldhauser & Ellsworth, 2000), using **differential times derived both from phase pick
times and from waveform cross-correlation**. We relocated the LFE families **together with nearby
seismicity from 01 January 2015 to 11 October 2024**, using phase pick times obtained from the
NCEDC" (§2). A 3-D model was considered and declined: the differences "would likely be minor,
particularly for relative locations" (§3).
**The geometry — the paper's headline claim.** "Double-difference-relocated hypocenters reveal that
the 27 LFE family hypocenters form an **arcuate alignment at 22–29 km depth, extending nearly 20 km
along a strike of ∼290° and dipping toward the north-northeast at approximately 45°**" (§3), with
"the majority of the identified LFE families locat[ing] near the **downdip edge and in the middle
along strike**". It sits "**∼15 km south of the southern edge of Gorda slab seismicity**", separated
from it by "a ∼15 km-wide **aseismic zone of pronounced low-P-velocity**" (§3), and along "a
pronounced velocity contrast" in the Furlong et al. (2024) model (Fig. 3). The dip is the point:
"the **northeast-dipping geometry** of the Mendocino LFE sources **is inconsistent with occurrence on
the subduction megathrust**" (§3) — hence the hypothesis that these LFEs "may **demark the southern
edge of Cascadia subduction**".
**Falsifiable temporal behaviour.** "small bursts of activity every few days in each family, with a
**typical recurrence time of only about 2 days**", "Only minor variation … among families" (§3),
restated as "typically every **2–3 days**" (§4); "**Notably absent from our data set are any families
that exhibit large, infrequent bursts of activity associated with episodic slow slip**, as observed
for much of the main Cascadia zone" (§3) — so the ETS signature the rest of the margin shows is
absent here by observation, not by our failure to see it; and "frequent episodes of **migration**
among LFE sources … with **migration velocities of ∼40–80 km/hr**" (§3). A geodetic stack was tried
and failed: "this initial attempt was **unsuccessful** in resolving any deformation" (§2).

**Method in one paragraph.** Stations are named only by network — "operated by the U.S. Geological
Survey (**NC** network), UC Berkeley (**BK** network), and the Earthscope Consortium (**PB**
network)" (Acknowledgements); **codes, count and component mix are SILENT**, and one rule tells us
his detection set is smaller than his stacking set: "Channels added in the final stacking (**only
used in stacking and not for event detection**)" (Fig. 2b). Band **2–9 Hz** (Fig. 2); candidate
window **6 s** (Fig. 2a); **sample rate and response removal SILENT**. Templates come from "a
**multi-step, iterative procedure of cross-correlation of candidate templates, stacking the
detections with the highest correlation sum across the network, and picking phase arrivals (P and S)
on the stacked waveforms**. The stacked waveforms are then used to construct new templates with
higher signal-to-noise ratios, with the process repeating", after Shelly and Hardebeck (2010), with
"**Iterative detection and stacking (2–3 iterations)**" (§2); the seeds were **manual**: "initial
waveform templates by **manual selection** of waveforms during times of PNSN-detected tremor …
multiple nearby time windows were tried as templates; we **retained the one that produced the
strongest detections** at other times" (§2). Detection is a network **correlation sum** thresholded
at "**9 times daily median absolute deviation**" (Fig. 4a) — the MAD is **daily**. The only
declustering, and the only dedupe between his densely-spaced families, is "we **retain only the
detection with highest correlation sum within each 4 s window**" (§3).

**What they released.** "**The full catalog of LFE detections is available in Shelly (2025)**" — a
USGS data release (doi:10.5066/P1TCKK7G, CC0), on our disk as
`lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv`: **27 families,
61,441 detections, 2018-01-01 → 2024-11-30**, locations 40.027–40.088 N, −123.713 to −123.532,
depth **22.2–28.6 km** (`lfe_south/results/SOUTH_TRUTH_INVENTORY.md` §1). **Templates are not
released; code is SILENT.** Waveforms via NCEDC and SAGE.

**Are we using their method? NO — five parameters differ, and two of ours were outright bugs.**
- band **1–8 Hz** (`lfe_south/scripts/t3_cache.py:17`) vs his **2–9 Hz**
- template **12 s**, 2 s pre-arrival (`t9_templates.py:29-30`) vs his **6 s**
- threshold **8×MAD** (`t11_mf.py:39`) vs his **9×MAD** — both daily, which does match
- declustering a **4 s greedy exclusion radius** (`t11_mf.py:29-35`) vs his **fixed 4 s window** max
- templates from **one** stacking pass at his published times, with no P/S picks and no re-detection,
  vs his **2–3 iterations** — and every (family, station) template enters the average with **unit
  weight**, including ones whose own stack SNR is 0.45 (measured p10 0.84 / median 1.50 / p90 3.52
  across our 270)
- **BUG 1**: `t11_mf.py`'s `ncc()` divides by `sd = sqrt(max(M·S2 − S1·S1, 1e-30))` from cumulative
  sums, and the 20 Hz cache **zero-fills gaps** with no validity mask (`t3_cache.py:31-32`), so a gap
  is an exactly-constant run where the true variance is 0, the clamp yields `sd = 1e-15` and the
  "normalised" cc reaches **2.4 × 10¹⁵**. That is the whole of the 12.8 M detections.
- **BUG 2**: each station's cc was summed **at a common lag** although each template is cut at its own
  offset `tt_s + tau_s` (spread **1.49–5.85 s** within family 1), so the network stack was incoherent
  and its time axis was not Shelly's.
- Ducellier's path was immune to both — his cache carries a `valid` mask and he requires a gapless
  block (`t7_duc_mf.py:113`), and his templates are all cut at the *same* absolute catalogue time on
  every station, so no per-station shift is needed. That is why the same chain reproduced him at
  recall 0.836 and failed here.

**Status. The 12.8 M-detection result is RETRACTED as a bug report against our code, not a statement
about Shelly.** The published failure (`SOUTH_RESULTS.md` §2: 12,848,734 detections against his
61,441, pooled recall 0.026, precision 0.000) is now explained in three parts: the unbounded cc
(BUG 1, worth ~600× too many detections — 128,922 on 2022-05-11 alone against his 71), the missing
moveout (BUG 2, which prevented true events from stacking), and a **scoring denominator fault** —
`t8_score.py` used all **61,441** of his events while the filter only ran on the **385 held-out
days**, where just **5,021** of them live, so the recall **cap was 0.082** and "0.026" was 0.32 of the
achievable maximum. With the live-window mask, the moveout shift and his own 9×MAD, the same day
gives **21** detections of which **11 match his 71** (precision 0.52, recall 0.155 on that day). The
full re-run over all 385 held-out days is in flight (`lfe_south/scripts/aa_shelly_queue.sh`, outputs
under `lfe_south/results/aa_shelly/`, report `AA_SHELLY_REPORT.md`) and will report a threshold
sweep, the statistic's value at every one of his published times, and the quiet-day false-alarm rate.
The existence result is untouched by all of this: **all 27 families** show a station above its own
200-draw null 99th percentile, median max-station z **37.1** (`SOUTH_RESULTS.md` §4).

**To reproduce properly from here.** Re-cache at **2–9 Hz**, build **6 s** templates, run his
**2–3 template iterations** (`aa_restack.py` is the first one), threshold at **9×MAD daily** with the
fixed-4-s-window dedupe, and score only on days actually processed. His station list, sample rate,
channel-selection rule and stack size are SILENT, so each must be swept as a column, not guessed.

---

## Summary table

| Catalogue | Region | Families | Located? | Released | Our method matches? | Status |
|---|---|---|---|---|---|---|
| Royer & Bostock 2014 NVI | N. Vancouver Is. | 47 final (54 initial) | Yes — Hyp2000 + hypoDD (91 % kept) | Nothing | **Being made literal** (ab_royer v2: 15 s, 40 sps, hierarchical clustering, Hyp2000+hypoDD) | **Added to ab_royer 2026-09-25** (was omitted); see notes/2026-09-25_Notes.md |
| Bostock et al. 2012 SVI | S. Vancouver Is. | **140** | Yes — hypoinverse, all 140 | Nothing | **Partly** (window + 8×MAD match; band, rate, stations, grouping, iterations differ) | **Reproduced: 58/96 blind, 22/30 sealed** |
| Chestler & Creager 2017 | Olympic Peninsula | **57 → 43** | Yes — 3-D Preston grid search + DD, all 43 | Waveforms only; Table S1 moments | **No** (no 18×/45×MAD stages, no 2nd-root, 6 s vs 15 s, no Rcom dedupe) | Spatial coverage only: 61/69 coords within 10 km vs 3/69 null |
| Sweet, Creager & Houston 2014 | W. Washington, downdip edge | **1** (8942 events) | Yes — hypoDD, 700 events, 100 m errors | Nothing (no data-availability section) | **No** (6 s vs 15 s; linear vs third-root; 8× vs 10×MAD; 200 vs 80 windows; 4 s vs 6 s separation) | Not attempted — the cheapest single-family positive control on the margin |
| Royer & Bostock 2014 NW | NW Washington | **122** final (224 initial) | Yes — Hyp2000 + hypoDD (77 % kept) | Nothing | **No** (6 s vs 15 s; **52 stations vs their 7 named**, 3 of which we lack) | **Literal rebuild in progress (ab_royer v2)**: autocorr 3258 vs 4915 (−34 %, open bug); clustering sweep; Hyp2000+hypoDD built — notes/2026-09-25_Notes.md |
| Royer & Bostock 2014 SW | SW Washington | **54** final (70 initial) | Yes — Hyp2000 + hypoDD (72 % kept) | Nothing | **No** (same as NW; 6 of their 7 stations available) | **Literal rebuild in progress (ab_royer v2)**: autocorr 2936 vs 3306 (−11 %) after cull fix (v1 was 7981) |
| Thomas & Bostock 2015 | Central Oregon | **18** | **Partly — 13 of 18**, only 4 with free depth | Nothing | **No** (they used cross-station, not autocorrelation; no PWS; 3 vs 10–15 iterations) | Queued |
| Plourde et al. 2015 | S. Cascadia 40–41.8N | **59 → 34** (66 unreduced) | Yes — all 34, 1-D Verdonck & Zandt; software SILENT | Nothing in paper; 66 templates reached us via Ducellier | **No** (no cross-station stage; their threshold is SILENT) | Reproduced indirectly, via Ducellier |
| Ducellier & Creager 2022 | S. Cascadia, 2004–2011 | **66 used → 60 kept** | No — adopted Plourde's (family A excepted) | **Everything**: code, templates, both catalogues | **Yes, almost** — only the MAD block (12 h vs his 1 h) and his second filter differ | **Reproduced: recall 0.836, 63/66** |
| Shelly et al. 2025 | Mendocino ~40.0N | **27** | Yes — all 27, grid search + hypoDD | Times + locations, CC0; no templates, no code | **No** (1–8 Hz vs 2–9; 12 s vs 6 s; 8× vs 9×MAD; 1 vs 2–3 iterations; greedy vs fixed-window dedupe) | **Re-run FINISHED 2026-09-25, still FAILS**: recall 0.48 @8×MAD but precision 0.006, quiet-day rate ≈ active-day rate (≈180 vs his 4 det/fam-day); residual pathology days (~10k det); existence 27/27 stands |

---

## Where the margin-wide catalogue actually stands — honestly

Two of the nine catalogues are genuinely reproduced: **Bostock's southern Vancouver Island set**
(58 of 96 open families and 22 of 30 sealed hold-outs, rediscovered blind) and **Ducellier's
southern Cascadia set** (mean recall 0.836, 63 of 66 families past the reference line, using his own
templates and stations). Those two anchor opposite ends of the margin and they are real results.
Everything between them is weaker than the project ledger has been implying. The Olympic result is
a *coverage* statistic against a 69-row coordinate file whose provenance our own filename disclaims
— not 61 of Chestler's 43 families identified — and an earlier comparison on the same data was
negative. Washington is a **failed run, not a null result**: this audit shows it was executed with
Bostock's 6-second Vancouver Island window and a 52-station probe set instead of Royer & Bostock's
15-second window and their seven named stations, three of which we have never even downloaded, so
it says nothing yet about whether their 122 NW and 54 SW families are recoverable. Central Oregon
has not been started, and the one paper covering it reports that the method we would bring
**fails there by design**. Mendocino has an existence signal on all 27 families, and the detector that
appeared to produce 12.8 million false alarms has since been traced to two numerical defects in our own matched
filter (an unbounded cross-correlation wherever the cache zero-fills a gap, and a network stack
summed without the per-station moveout) plus a recall denominator twelve times larger than the
days we actually ran — on one held-out day the fixed chain gives 21 detections, 11 of them his.
And the two papers added on 2026-09-25 move one region from 'no recipe' to 'recipe in hand': Shelly's
Mendocino parameters are now on the table (2–9 Hz, 6 s, 9×MAD daily, fixed-4-s dedupe, 2–3 template
iterations), and Sweet et al. 2014 supplies the full version of the array front end that Chestler
compresses into a paragraph (15 s windows, third-root stack, 80-window templates, 10×MAD daily, a 6 s
minimum separation derived from S−P), together with a one-family positive control that needs no
catalogue file to score. The single most
consequential finding of this audit is that these nine catalogues are **not one recipe** — window
lengths run 6 s to 60 s, thresholds run 8× to 45×MAD on three different statistics, and two of the
regions were built with a cross-station detector we have not implemented at all — and that our
chain has been applying one region's parameters to all of them.
