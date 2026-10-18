# S11 — why we did what we did (plain-language reasoning, kept up to date)
Last updated 2026-10-12 2:25 pm PT. Detailed decision log: NOTES.md. Designs/reviews: MERLIN_DESIGN.md, MERLIN_REVIEW_1.md, MERLIN_FORECAST.md.

## The question
An LFE family's coda (2–5 s after S) samples velocity along all scattered paths between a deep source (~35 km) and the station.
A dv/v value mixes deep (near-source, 25–45 km) and shallow (near-station, top few km) change. The user wants them separated, to see
whether the deep fault zone changes velocity, e.g. episode by episode.

## Step-by-step reasoning
1. **Measurement recipe.** The family-10/PGC tests picked stretching anchored at the picked S, 3 components, 2–5 s. 3C plus the longer window gave the steadiest values.
2. **Denoised vs raw stacks.** The user's recipe used SVD-Wiener denoised day stacks. Our own injection test had shown that denoising shrinks real
   changes, so we ran both paths and let a pre-set test choose. Denoised recovered 26 % of an injected change, raw 114 %. → raw is the headline.
   The denoised path looked steadier only because it flattens real change.
3. **Kernel choice.** The measurement is anchored at S, so the sensitivity must be the coda's extra travel time after S. We computed
   S-pinned kernels by simulating scattered waves through the regional velocity model. Check: the pieces add up to 1 (passed exactly).
   - Result: about 55 % of the sensitivity is in 25–45 km and 33 % in the top 5 km.
   - The deep share rises with family depth and distance; that is the lever.
4. **"Deep" = a 25–45 km layer, not the fault surface itself.** A thin layer at the interface has sensitivity that flips sign for a third of the
   pairs, so it can't be inverted honestly.
5. **First check, no kernels needed.** PGC's surface vault and the 200 m boreholes show the same yearly cycle. That rules out an instrument
   artefact. Later correction: it says nothing about depth, because the kernels predict they would agree whatever the depth.
6. **Pilot result.** No long-term deep change (+0.04 ± 0.12 %). Per-episode deep could not be resolved: it traded off against a shallow change
   shared by all sites.
7. **Why, and the fix (Merlin forecast).** Every family–station series has its own unknown zero. With the shallow term free every episode,
   those zeros swallow any per-episode deep signal; no amount of data fixes that. What *is* answerable: is an episode faster or slower at depth
   than its neighbours (ΔD), assuming shallow change is smooth in time (yearly cycle plus slow trend)?
   - Rationale: rain and groundwater don't jump in step with slow-slip episodes.
   - Price: about half of any shallow change that happened exactly during an episode would show up as "deep". We declared this up front.
8. **More data, chosen by expected gain.**
   - Keep smaller stacks (n ≥ 20) weighted by their measured precision: a noise-derived cut, not tuned.
   - Fix a bug that dropped 32 borehole cells.
   - Add the 39 missing families.
   - Distant stations held back: speculative gain, plus a new leak from the middle crust.
9. **Result (SVI).**
   - Precision about ±0.09 % per big episode; no episode stands out; the average of 15 big episodes is 0.000 ± 0.041 %.
   - Scrambled control and leave-one-site-out (including leave-PGC-out) pass.
   - Injection checks show some mixing: deep pulse returned at 114 %, shallow leak 45 %. So the reading is: no episode-locked deep change above
     about 0.2 % per big ETS, and below about 0.08 % on average. The mixing would matter only if something had shown up.
10. **Yearly cycle.** The station pattern says it is shallow. Stations that see the deep layer most show the smallest cycle. Its deep share is
    reported only as an upper limit (< 0.28 %).
11. **Next: Washington NW** (user said "proceed as you think best"). It is an independent region with its own deep zone, so it can replicate the
    null with the identical pipeline and pre-set checks. Distant stations stay deferred.

## Principles we kept to
- Pass rules are written before each run; any rule changed after seeing data is logged as a revision with its evidence.
- Every claim needs a control (scrambled control, injection, leave-one-out) and a stated limit.
- Separation is measured by injection through the whole chain, never assumed from kernel numbers.
- Plots follow the user's style: no error bars on dv/v; colour by CC or weight.

## Round 3 (2026-10-12 evening): three ways to keep going
12. **Is the deep term just a kernel choice?** In Washington, a wrong scattering model turned part of a shallow yearly cycle
    into "deep". Averaging four kernel models fixed that, but it was tested with kernels inside its own average. Round 3 tests
    it honestly: leave one kernel model out of the average, make the fake signal with that left-out model, and check whether
    the average still recovers it. Two new models lie outside the old range (much stronger scattering), so this is a real
    out-of-sample test. A second route lets the data choose the kernel: the late and early parts of the coda weight deep and
    shallow layers differently, so the kernel that best explains both at once is preferred.
13. **Can REDPy-style stacking fix the over-read?** The episode-vs-reference measurement reads a pulse 14 % too large, because
    the pulse also leaks into the neighbours' references with the opposite sign. An all-pairs graph of episode stacks has no
    shared reference, so it should read the pulse at its true size (a test on two families gave 0.30 for 0.30). Weighting each
    LFE by how well it matches the family's waveform should cut the noise by ~10 %. The shallow leak (a shallow change exactly
    in an episode partly counted as deep) is a property of the model, so it is not expected to go away; it is measured again.
14. **Does one part of the fault behave differently?** Coarse zones along strike (3 per region) can be resolved only if the
    noise allows it, so the expected precision is computed first from the cell errors and kernels alone, without looking at
    any dv/v value. To split the deep sensitivity along strike properly, the kernel simulation was re-run to record where
    along the fault each scattered path spends its time at 25–45 km depth.

## Lever 4 (2026-10-12 evening): stations 80–150 km away
15. **Why far stations.** A station far from the family sees coda that spent more of its time deep and less near the
    receiver, so its deep and shallow weights differ from those of the near stations. Two different mixes of the same two
    unknowns can be told apart; one mix cannot. The price: far signals are weak, and the mid-crust (5–25 km) and the mantle
    below 45 km now matter, so both get their own smooth terms.
16. **First check: can far stations see the families at all?** Stack the repeats on a few big days per year at every far
    station and pick S. About 1 pair in 18 shows a clear S and coda above the scrambled control, almost all at 85–110 km
    (Nanaimo, the PGC boreholes for Washington families, Olympic boreholes for Vancouver Island families). Beyond ~115 km the
    coda is not seen. Enough for > 100 cells per region, so the study continues.
17. **A sharper question the far stations can answer.** If one regional shallow change per episode is allowed (rain, say,
    hitting every station in the same episode), the near network alone can just barely separate it from deep change
    (±0.18 % per episode). Far stations should sharpen that, because their deep/shallow mix differs most. This is the
    cleanest test of the lever and is reported beside the pre-set ones.
15. **What round 3 showed (9:30 pm PT).** The deep term cannot be made independent of the kernel model. But what changes with
    the kernel is mainly its size (about ±30 %), not whether something is there. The data themselves favour a kernel with
    strong scattering near the surface, in both regions. The graph estimator fixed the 14 % over-read. The shallow leak
    remains, because it comes from the model, not the measurement. Zones along strike show no part of the fault behaving
    differently, and their differences are almost immune to the shallow leak. One SVI episode (2019-11) stands out in the
    graph estimator at 3.6σ. It fails the pre-set checks, so it is a candidate to test, not a finding.
18. **What the far stations bought.** Adding them cut the shallow-to-deep leak by about a quarter to a third and the
    per-episode error by a quarter to a third, in both regions. Still, no slow-slip episode stands out on average:
    < 0.07 % (SVI) and < 0.035 % (WA). Allowing one regional shallow change per episode now costs only ×1.5 in error and
    removes the leak almost entirely. Paired with the graph estimator, the leak passes its rule for the first time.
19. **The yearly cycle, tested without kernels.** If the cycle were deep, far stations (more deep weight) should show
    more of it. In SVI they show less, so SVI's cycle is not deep. In WA they show the same amount. But those far stations
    are the PGC boreholes, where the cycle looks the same whatever the source, so it follows the station rather than the
    path. The two regions disagree, and the station-side reading explains both.
20. **2019-11.** Stations far from Vancouver Island, none of them PGC, also see a speed-up in that episode, about half
    as large (+0.28 ± 0.14 %). Combined it is +0.32 ± 0.11 %, just under 3σ. It is still a candidate.
