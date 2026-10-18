# Sources that fill gaps in Royer & Bostock (2014, EPSL 402:247-256)

Compiled 2026-09-25 (web search + on-disk PDFs). No code was run on data. `METHOD_FROM_PAPER.md` was
read but not edited.

**Bottom line: no public catalogue, detection list, template waveform set or template-location table
exists for the NW/SW (or NVI) templates.** Every follow-up paper and supplement that was checked
distributes figures only. The only table in the thesis appendices is the hypoDD parameters (same as
Table S1) plus a per-family tidal table for the **SVI** catalogue (Bostock 2012 families, not
Royer's).

## Files in this folder
| file | what |
|---|---|
| `royer2014_phd_thesis.pdf` | A.A. Royer, *Studies of seismic deconvolution and low-frequency earthquakes*, PhD (Geophysics), UBC, 2014 (graduation 2014-09), 134 pp. hdl:2429/49978, DOI 10.14288/1.0165914, CC BY-NC-ND 2.5 CA. MD5 c6d3b8997685a373dc137cb80848bd20 (matches DSpace) |
| `royer2014_phd_thesis.txt` | pymupdf text dump with `=====PDFPAGE n=====` markers (PDF page = printed page + 13) |

Download route (the open.library.ubc.ca web UI is bot-blocked): DSpace REST
`https://circle.library.ubc.ca/rest/handle/2429/49978?expand=bitstreams` then
`https://circle.library.ubc.ca/rest/bitstreams/203388/retrieve`. Full text is also available as JSON
from `https://oc-index.library.ubc.ca/collections/24/items/1.0165914` (field `FullText`).
Royer has **no UBC MSc**. Her BSc (2007) and MSc (2009) are both from Université Joseph Fourier
(thesis title page), so the PhD thesis is the only UBC thesis.

## Per source

### S1. Royer PhD thesis (UBC 2014)
* **Ch. 3** (pp. 35-51) is the EPSL paper almost verbatim. Its §3.3.1-3.3.3 (PDF pp. 52-53) give the
  same text as the paper: 15 s windows lagged by 0.5 s, 7170/h, 8×MAD, "culled to exclude overlapping
  detections and retain only those with high SNR", 4915/3306/4267, 224/70/54, "subsets of 7 to 10
  high SNR, three-component stations", "several iterations", 122/54/47. The paper's silences stay
  silent: no linkage, cut, minimum size, SNR definition, template length or iteration count.
  NVI "different combinations of 7 of the 10 three-component stations" is also not expanded.
* **Ch. 3 §3.5** (PDF p. 57): "one using a standard linearized inversion [Hyp2000, Klein, 2002] and
  one using the double difference location algorithm [hypoDD, Waldhauser, 2001]". hypoDD uses
  "a combination of ordinary phase picks from our LFE catalog and differential travel times from
  phase correlation of P- and S-waves". hypoDD retention: NVI 91 %, SVI 92 %, NW 77 %, SW 72 %.
  §3.7 (PDF p. 63): "the 1-D, **fixed V_P/V_S** model used for LFE location". The model itself is
  **not named**.
* **Ch. 1** (PDF p. 16-17), on the growth loop: "The template can be used as an updated Green's
  function estimate and the same process is then iterated **until there is no further improvement
  in SNR**." It also describes the template as "a **small windowed section** of a Green's function
  estimate … at a number of channels".
* **Ch. 4 §4.2** (PDF p. 73), the tidal study that became Royer, Thomas & Bostock 2015:
  > "The northern Washington (NW) catalogue was **regrouped into 100 geographically independent LFE
  > families** and includes a total of **60,527 detections** from ETS episodes between **2007-2011**
  > … The southern Washington (SW) catalogue, **regrouped into 49** geographically independent LFE
  > families, includes a total of **22,868 detections** from ETS episodes during the period
  > **2007-2008**."
  SVI: "culled to 93 geographically independent LFE families with a total of 148,847 detections
  (2003-2013)".
* **Fig. 4.1** (PDF p. 74) maps all family locations with a colour-coded family ID ordered by latitude
  N to S. Implied numbering: SVI 1-93, NW 94-193, SW 194-242. The examples are NW family 123 and SW
  family 242. Figures 3.4/3.5/A1/A2 map the Hyp2000/hypoDD epicentres and depth profiles. These are
  **figures only** and could be digitised. No coordinates are tabulated.
* **Fig. 4.2**: the 2007 NW/SW tremor-episode time axis runs as "days from 16/01/2007".
* **Appendix A**: Table A1 = hypoDD parameters (identical to EPSL Table S1). Figs A1-A4 = hypoDD
  epicentres, Hyp2000 depth profiles, condition numbers, strike/dip/rake histograms.
  **Appendix B**: Table B1 = per-family UDSS/FNS N_ex and detection counts for **SVI only**.

### S2. Royer, Thomas & Bostock (2015), JGR 120:384-405, doi:10.1002/2014JB011430
Wiley returned 403 to every fetch route and OpenAlex lists no open-access copy. Its content is thesis
Ch. 4 (submitted July 2014, per the thesis preface). It gives the 100/60,527 NW and 49/22,868 SW
numbers quoted above, as the web search summary confirms. No data repository was found. Whether the
published SI adds anything beyond thesis Appendix B was **not verified**.

### S3. Bostock, Royer, Hearn & Peacock (2012), G-cubed, doi:10.1029/2012GC004391 (on disk, notes/)
This is the paper Royer defers to.
* §2: data "filtered between **1 and 10 Hz**", 40 sps. Royer uses 1-8 Hz.
* §2.1 [5]: "**6-second** long windows with 0.5 s overlap … 7188 windows per hour". Royer uses
  15 s / 7170. Records are "**summed** to produce a network autocorrelation record", and detections
  are taken "when the network autocorrelation exceeds a threshold of 8 times the median absolute
  deviation". The paper also notes: "we are not able to detect events that occur less than twice per
  hour."
* [6]: "5775 initial LFE detections … grouped into 140 families using a hierarchical cluster analysis
  [Rowe et al., 2002] based on correlation coefficients of composite wave trains constructed by
  lining up all station/component records corresponding to an individual detection window". **No
  linkage, cut or minimum size is given.**
* §2.2 [7]: "Using subset templates comprising **6 or 7 stations** … Selection of stations for the
  subset template is based on the **SNR level of the template waveforms and … the availability of
  station data for a given tremor episode**." Stations operating outside the seed episodes "were
  incorporated within the LFE template … using detections within the appropriate recording period".
  Iterations: "**After 5 such iterations** neither template SNR nor numbers of detections continued
  to show significant improvement." Templates hold 1000-5000 detections.
* §3 [8]: "Although some of the 140 composite waveform templates **exhibit overlap in contributing
  detections**, we will consider each template as representing a distinct LFE location." P and S
  "first arrival times were picked from high SNR waveform templates", giving about 16 P and 22 S
  picks per template, inverted with **hypoinverse** (Klein 2002). **The velocity model is not
  given.**

### S4. Bostock, Thomas, Savard, Chuang & Rubin (2015), JGR 120:6329-6350, doi:10.1002/2015JB012195 (on disk)
* §2: the SVI catalogue "comprises those templates from Bostock et al. [2012] and **Royer and Bostock
  [2013]**, as well as new additions for a total of 130 distinct source regions". Declustering:
  "culled such that **only the detection with the strongest correlation in any given 4 s
  interval** (approximate maximum duration of the S wave coda) is retained". This is the group's
  declustering rule.
* Step 1: template windows are "**6 s** in length", 1-8 Hz, "aligned using predicted times from the
  1-D model used for routine earthquake location by the **Geological Survey of Canada**".
* **Table 1** (p. 6332), the velocity model: layer tops at 0 / 1 / 6 / 30 / 45 km, with
  V_P = 5.0 / 6.0 / 6.7 / 7.1 / 7.75 km/s; "**A Poisson solid is assumed**", which gives
  V_P/V_S = √3 ≈ 1.732. This is the only 1-D fixed-V_P/V_S model the Bostock group has published for
  this work. It is a strong **candidate** for Royer's "1-D, fixed V_P/V_S model", but that link is an
  inference, not a statement.

### S5. Thomas & Bostock (2015), Tectonophysics 659:111-123 (on disk: `1-s2.0-S0040195115003807-main.pdf`)
* Paraphrasing Royer: "Detections are then culled to eliminate duplicates and retain only
  high-quality detections (Royer and Bostock, 2014)."
* The group's family-dedup rule (used for central Oregon): "no family shares more than **10 % of
  detection times**, defined as a difference between detection times of **<5 s**, with any other
  family."
* HYPOINVERSE with the **PNSN O0 model** and the Williams et al. (2011) local model. This shows the
  group used PNSN's model when working in US territory, which is weak evidence for a PNSN model in
  Washington.

### S6. Sweet, Creager & Houston (2014), G-cubed, doi:10.1002/2014GC005238 (on disk)
p. 2: the family "lies … at the downdip edge of **122 LFE families** located in northern Washington
[Royer and Bostock, 2013]. Four of their LFE families were at similar depths … Their nearest LFE
family is **10 km south** of ours, and their next nearest is **20 km west**." This is a geometric
anchor only. Sweet's own template was built from a **15 s** window. No Royer data are included.

### S7. Chestler & Creager (2017), JGR, doi:10.1002/2016JB013717 (on disk)
Its catalogue has 43 families and 32,387 LFEs from Array of Arrays data (2009-2011). The paper says
its locations are "similar to LFEs in other catalogs beneath the Olympic Peninsula from Royer and
Bostock [2014]" and that Royer "did not use data from the AofA". No Royer data are included, and
only waveforms are declared public (IRIS).

### S8. Savard & Bostock (2015), BSSA 105(4), electronic supplement
URL: https://www.seismosoc.org/Publications/BSSA_html/bssa_105-4/2014301-esupp/index.html
Figs S4/S5 show cross-station LFE epicentres for **northern Washington 2010 and 2011**. The web
search summary gives 2176 and 1568 detections. The supplement has PNG figures only: no catalogue
and no mention of Royer's templates.

### S9. Rowe et al. (2002): the correct citation
Both Royer & Bostock 2014 (reference list) and Bostock 2012 cite **Rowe, C.A., Aster, R.C., Borchers,
B., Young, C.J., 2002. An automatic, adaptive algorithm for refining phase picks in large seismic data
sets. BSSA 92, 1660-1674, doi:10.1785/0120010224**. They do **not** cite the JVGR Soufrière Hills
paper. The Soultz/Soufrière "Using automated, high-precision repicking…" papers are separate Rowe
et al. works (PAGEOPH 2002 for Soultz).
The BSSA paper could **not be obtained**: it is paywalled (GSW 403), with no OA copy on OpenAlex,
Semantic Scholar or CiteSeerX. The nearest same-group source is Young, Merchant & Aster (2001),
"Comparison of cluster analysis methods for identifying regional seismic events" (23rd SRR,
https://www.ldeo.columbia.edu/res/pi/Monitoring/Doc/Srr_2001/Screen/08-06.pdf):
* similarity = "absolute value of the normalized waveform correlation (i.e. best lagged
  cross-correlation divided by …)", with distance = 1 - similarity;
* Lance-Williams dendrograms. Of the six linkages tested, the **"flexible" method** (weights of
  Ludwig & Reynolds 1988, i.e. β = -0.25) "consistently produced the clearest dendrograms … we use
  it as our default";
* "Deciding where to cut the stems of a dendrogram is a **subjective process** … set the level at
  the highest possible similarity value which still yields distinct groups."

This is context only. It does not establish what Rowe 2002 or Royer actually used.

### S10. Pre-2009-08 western-Washington tremor catalogues
* Wech & Creager's map (`tunk.ess.washington.edu/map_display`, archived 2011-2015) served epicentres
  through a CGI (`/cgi-bin/get_tremor.py`, not archived). Its calendar `min: 20090801`
  (`cal_stuff.js`, Wayback 2012-02-09) shows the online database **starts 2009-08-01**. The
  "Download" form e-mailed requests and was never automatic.
* Wech, Creager & Melbourne (2009), JGR 114:B10316, doi:10.1029/2008JB006090: about 15,000 epicentres
  for the Jul 2004, Sep 2005, **Jan 2007, May 2008** ETS episodes and the **Feb 2007-Apr 2008
  inter-ETS** window (per the search summary). Wiley was blocked and no machine-readable SI was
  found. Wech & Creager (2008), GRL 35:L20302: method paper, same status.
* The local `catalogs/pnsn_tremor_cascadia_full.csv` starts 2010-01-06.
* **Result: no public day list or epicentre file for 2007-2008 was found. Nothing was downloaded.**

## Gap table

| gap in Royer & Bostock 2014 | answer found (source) | status |
|---|---|---|
| released catalogue / detections / template waveforms | none anywhere (S1-S8) | **open**: not public |
| template locations (coordinates) | figures only: thesis Figs 3.4, 3.5, A1, A2, 4.1 (S1). Sweet's family has Royer families 10 km S and 20 km W (S6) | open (could be digitised) |
| final detection counts | NW 60,527 in 100 regrouped families (2007-2011); SW 22,868 in 49 families (2007-2008) (S1 Ch. 4, S2) | **found** (a count target, after regrouping) |
| family regrouping / dedup after growth | "regrouped into geographically independent families" 122→100 NW, 54→49 SW (S1). Group rule elsewhere: families sharing >10 % detections within 5 s are duplicates (S5). SVI templates allowed to overlap in detections (S3) | partly found (rule inferred from group practice) |
| clustering linkage / threshold / min size | Rowe 2002 = **BSSA** paper (S9), not obtained. Group-adjacent practice: flexible linkage, subjective cut (S9) | **open** |
| cull rule (overlap) | Bostock-group declustering: keep strongest detection per 4 s (S4, for growth detections) | partly found |
| SNR definition / threshold | none (S1, S3) | open |
| autocorrelation window | 15 s (Royer). Bostock 2012 used 6 s; template windows 6 s in S4 | Royer's value stands; template length open |
| growth template length | "small windowed section" (S1 Ch. 1). 6 s in S4 (SVI magnitudes) | open (6 s = group precedent) |
| number of growth iterations | Bostock 2012: 5, stopping when SNR and detection count stop improving (S3). Thesis: "until there is no further improvement in SNR" (S1) | **found** (stop rule), count ≈5 by precedent |
| growth station subset rule | rank by template SNR **and** per-episode data availability; stations enter via detections in their recording period (S3) | found (precedent) |
| NVI "combinations of 7 of 10" | no expansion (S1) | open |
| tremor dates for growth, 2007-2009 | the NW catalogue spans ETS 2007-2011 and SW 2007-2008 (S1). Pre-2009-08 tremor catalogue not public (S10) | open for day bounds |
| velocity model / Vp/Vs for Hyp2000 | "1-D, fixed V_P/V_S model" (S1, EPSL §7). Candidate: GSC 1-D model, Poisson (Vp/Vs 1.732), S4 Table 1. PNSN O0 used by the group in Oregon (S5) | partly found (model named nowhere) |
| pick method | manual picks of "unambiguous zero-phase" P/S on Hilbert-rotated templates (S1). Bostock 2012: first-arrival picks on high-SNR templates, ~16 P + 22 S per template | found (manual, qualitative) |
| hypoDD parameters | thesis Table A1 = EPSL Table S1 | confirmed |

## Contradictions and tensions with `lfe_array_discovery/results/ab_royer/METHOD_FROM_PAPER.md`
1. **B3, Vp/Vs sweep with headline 1.78.** The paper and thesis say "1-D, **fixed** V_P/V_S". The only
   Bostock-group 1-D model published (S4 Table 1) is a Poisson solid, **1.732**. 1.73 is the better
   supported headline, and the GSC 5-layer model should be added as a swept model alongside
   Crosson (1976).
2. **Step 5c, declustering "SILENT".** The group's published rule is to keep the strongest
   detection per 4 s (S4). This is precedent, not Royer's own statement.
3. **Step 5e, iterations "SILENT".** The precedent is 5 iterations, stopping when template SNR and
   detection count stop improving (S3, S1). This is consistent with `--max-iters 8 --tol 0.01`, but
   the stop criterion in the sources is **SNR plus count**.
4. **Final counts.** METHOD scores against 122 / 54 templates. The downstream catalogue is 100 / 49
   *regrouped* families with 60,527 / 22,868 detections (S1). If our pipeline merges co-located
   templates, 100 / 49 is the matching target, and the detection totals are a new check.
5. **B5, the Rowe reference.** METHOD cites "Rowe et al. 2002" without venue. It is **BSSA 92:1660**
   (S9), not the JVGR Soufrière Hills paper named in the task brief.
6. **Station subset (A4/B1).** Bostock 2012 also ranks by **per-episode availability**, not only by
   SNR. Stations are added per recording period. METHOD's "rank by SNR, top ≤10" matches in spirit,
   but availability per episode should drive the subset.
7. **B2, pre-2009-08 inter-ETS days.** METHOD says there is no public machine-readable catalogue.
   This is confirmed: the Wech database starts 2009-08-01 (S10).
