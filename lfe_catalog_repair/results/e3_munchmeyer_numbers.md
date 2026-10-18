# E3 step 0 — Münchmeyer et al. (2024) LFEDetect paper: the numbers we needed
Source: **arXiv:2311.13971v1**, "Deep learning detects uncataloged low-frequency earthquakes across
regions", Münchmeyer, Giffard-Roisin, Malfante, Frank, Poli, Marsan, Socquet (23 Nov 2023).
PDF fetched 2026-09-11 to `results/external/papers/arxiv_2311.13971.pdf` (3.1 MB, 15 pages).
**The supplement is NOT in the arXiv package** (`/src/2311.13971/anc` → HTTP 404; the v1 PDF is the
main text + references only), and every methodological detail we most wanted — the PyOcto
parameters, the pick thresholds, the event filters — is deferred to it ("we provide an overview of
the workflow with further details in the supplement"). What follows is what the main text states,
plus what we could compute ourselves from the catalogues on disk.

## 1. Their Cascadia DL catalogue vs the Bostock reference, same period
Their Cascadia study periods (§3, verbatim): **2003-02-26 to 2003-03-10, 2004-07-02 to 2004-07-27,
2005-09-03 to 2005-09-25** = 12 + 25 + 22 = **59 days**.

| quantity | value | source |
|---|---|---|
| DL Cascadia catalogue | **8,369 events** = 141.8 /day | paper §3 |
| Bostock et al. (2015) in the SAME 59 days | **148,650 events** = 2,519.5 /day | computed here from `results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv` (269,423 events, 2003-02-26..2013-10-13) |
| ratio DL / reference | **0.056 (~18x smaller)** | — |
| per-window (Bostock) | 31,578 / 59,606 / 57,466 | 2,632 / 2,384 / 2,612 per day |

So Merlin's recollection is confirmed and quantified: **their Cascadia DL catalogue is 5.6 % of the
template-matching reference over the identical window**, and the paper frames the deficit as a
strength ("uncataloged sources"), not as recall. Their own conclusion says it plainly: template
matching "leads to a more sensitive model, as evidenced by higher event counts... it allows
template matching to identify LFEs with fewer stations than deep learning."
**This is the single most important calibration for E2/E3: our vi2010 cascadia N>=4 thr-0.3 run
produced 100 events/day against 2,663 Bostock events/day in that window — i.e. we are at the same
place on the same curve as the authors of the picker, not below it.**

## 2. Their recall / precision definitions — and their chance floor
Verbatim (§4): *"In Cascadia, for **60 %** of the LFEs from our workflow, we find an LFE in the
reference catalog **within 10 s**; for Guerrero for 70 % of the events. For Nankai, only 6 % of our
LFEs are in the reference catalog, however our catalog also substantially exceeds the original
catalog in the total number of events. **Conversely, we recover 32 % of the cataloged events.**
Note that a loose threshold for matching is justified due to uncertainties in the origin times due
to inaccurate locations."*

* So their headline numbers are **one-way match fractions at a ±10 s tolerance**, with NO chance
  correction and NO shifted-catalogue null. The 60 % is a **precision-like** quantity (fraction of
  DL events with a reference event nearby); the 32 % **recall** sentence follows the Nankai
  discussion and reads as Nankai's, not Cascadia's — and it cannot be Cascadia's, because
  32 % of 148,650 = 47,600 > their entire 8,369-event catalogue.
* Cascadia recall implied by their own numbers is therefore at most
  8,369 x 0.60 / 148,650 = **3.4 %**.
* **We measured the chance floor of their criterion.** Drawing 2x10^5 uniform random times inside
  each of their three windows and asking for a Bostock event within ±10 s gives
  **0.310 / 0.249 / 0.289, mean 0.283** (the Poisson estimate 0.44 is too high because Bostock's
  events are burst-clustered, so random times mostly land between bursts). Their 60 % is thus
  **2.1x chance**, and chance-corrected only **(0.60-0.283)/0.60 = 53 %** of their "matches" are
  informative. At **±2 s** — the tolerance E2 used — the same chance floor is **0.088**.
  => E2's precision-vs-Bostock of 0.15-0.22 at ±2 s (chance-corrected excess 0.13-0.22) is
  **not obviously worse than the published Cascadia result once both are put on the same footing**;
  it is measured against a chance floor we computed and they did not.
* Neither Lin's 90/95 % nor Münchmeyer's 60 % is a continuous-record precision. Münchmeyer's pick
  residuals: "average residuals are about 0.2 s" with Cascadia the lowest of the four regions,
  and at 0 dB SNR the residual standard deviations are **0.3-1.3 s** (§2, Fig. S3) — an explicit
  statement that pick precision degrades with SNR, consistent with F2's SNR-limited reading.

## 3. Their workflow and PyOcto parameters
* Workflow (§3): SeisBench LFEDetect P+S picks on continuous data → **PyOcto** association →
  **NonLinLoc with a 1-D velocity model** for absolute location → *"To avoid false detections, we
  filter the events based on the number of phase picks and the location residuals."*
  A PhaseNet/INSTANCE earthquake catalogue is run on the same waveforms and events detected as
  both are removed.
* **`n_picks`, `pick_match_tolerance`, velocity tolerance, the pick threshold and the
  pick-count / residual filter values are NOT in the main text — they are in the unavailable
  supplement.** We therefore cannot check E2's `n_picks=3 / n_s_picks=3 / tol 2.5 s / edt std
  1.5 s` against theirs. Recorded as an open item; it needs the supplement (or the published JGR
  version), not another guess. Our own operating point must be justified by our own null (step 1).
* Their event-level depths are explicitly not trusted: *"In all catalogs, the event depth is poorly
  constrained... We suggest that this is caused by the high timing uncertainty of the picks. In
  particular, a high P pick uncertainty will cause poor depth constraints as the P to S time is
  indicative of depth. Template matching catalogs alleviate this problem by locating LFE families
  instead of single events."* — the same conclusion E2 reached from the -24 km depth bias, and an
  independent endorsement of the family-location product Merlin's step 7 re-targets us to.
* San Andreas was **excluded** from their analysis because a small-aperture borehole array plus
  high pick uncertainty gave locations "broadly distributed in space instead of along the fault" —
  i.e. the published workflow also fails on unfavourable geometry.

## 4. Do they report a chance-association control? YES — and this is the direct precedent for step 1
Verbatim (§4): *"First, we rule out spurious detection. To this end, we **scramble the picked
arrival times of each station by applying small random shifts. We choose constant shifts for each
station for every hour.** This destroys the exact times, while keeping the pick distribution, the P
to S times within a station, and the higher number of picks within tremor bursts intact. The
scrambled "catalogs" only contain about **5 to 10 %** of the number of events contained in the
original catalogs and show no spatial coherence (Figure S8). **Even these numbers are still likely
an overestimation of the false positive rate**, as events recorded at many stations are likely to
be unperturbed by our scrambling procedure. The more than 30 % additional detections can therefore
not be attributed to spurious associations."*

Differences from our step-1 null, which matter when comparing the numbers:
1. **Their shift is per station PER HOUR; ours is per station over the whole record.** Theirs
   re-randomises 24x/day, so it preserves the hourly rate envelope (and hence burst structure)
   more tightly than ours, which preserves only the daily envelope. Neither is uniformly more
   conservative: hourly re-shuffling can leave sub-hour bursts partially co-aligned.
2. **"Small random shifts"** — the magnitude is not stated in the main text (supplement again).
   Ours is ±[300, 3000] s, which is >> the ~30-60 s association aperture, so real coincidence is
   certainly destroyed.
3. Their 5-10 % is quoted for the catalogue **after** their pick-count + residual filters, on their
   own (denser, 2003-2005 PO/CN) networks — not at a stated N>=3/4/5 station cut. Our comparison
   must therefore be made at matched N, which step 1 does.
4. They read a **5-10 % null as clean**, and their own caveat says it over-estimates the false-alarm
   rate. That is close to Merlin's `<=15 % -> PASS` branch, so the published precedent and the
   audit's acceptance bar agree on the scale.

## 5. Items still open because the supplement is missing
* PyOcto configuration (n_picks / n_s_picks / tolerances / velocity model tolerance).
* Pick probability threshold(s) and the "number of phase picks and location residuals" filter.
* The exact magnitude of the scrambling shifts, and whether the 5-10 % is per region or pooled.
* Figure S6 (match statistics), S8 (scrambled catalogues), S10 (uncataloged detections vs family
  stacks) are all referenced but not in the arXiv PDF. **S10's stated result is directly relevant
  to E3 step 4**: *"For Cascadia, the distribution of correlation values for these uncataloged
  detections are identical to the noise distribution, i.e., the new detections do not match known
  sources."* — the authors already ran the template-as-classifier test on their own uncataloged
  Cascadia detections and got a NULL result. Step 4 must therefore be read with the expectation
  that most associated events will NOT match a Bostock family template; the informative question is
  the excess over our own shifted-pick null, not the absolute fraction.
