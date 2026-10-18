# RESULTS — PNSN tremor vs Gualandi (2025) slip, relative geometry 2010-2026 (2026-09-10; S1-S16 complete)

Scope: north (>=46N) and central (44-46N) Cascadia, 27 tremor-confirmed slip patches Mw>=6.0 (unbridged catalog, thrust slip
within the trailing 30 d > 5 mm); the south is excluded (model at its noise floor: tremor-confirmed 23 % vs 4 % chance; S6-B
shows the southern signal is a sign-symmetric two-component common mode at the model's own timescale, not seasonal, stationary).
All numbers: see NOTES.md ledger and results/*.csv. Methods: PLAN.md (pre-registered), scripts/s*.py (docstrings).

## What holds (survives the mirror / time-shift / completeness controls)
R1  Tremor and modeled slip are CO-LOCATED at the model's resolution. Per event, tremor sits a median 4.7 km down-dip and
    0.9 km deeper than the slip-weighted centroid (81 % down-dip), IQR 9 km / 1.7 km — inside one 28-km mesh cell. Back-slip
    "events" give 18.6 km / 5.3 km with 3x the scatter, so this is not mesh geometry. (S1)
R2  Inter-ETS background tremor around these patches is DEEPER than the tremor during the events (group z3 = -0.92, 24/27
    events; per-event not significant). Consistent with down-dip inter-ETS tremor in the literature; here from a slip-anchored
    reference. (S1 null diagnostics)
R3  NO secular drift 2010-2026 of any relative-geometry quantity, with POWER stated (80 %): depth offset -0.04 km/yr [-0.16,
    +0.09], smallest detectable 0.2 km/yr = 3.2 km over 16 yr (inside the PNSN location error); across-strike offset (MDD 15 km);
    lead/lag (MDD 14 d); tremor depth percentiles (north |slope| <= 0.02 km/yr); slip edges; recurrence phase. Nucleation latitude
    and propagation direction (model-free tremor-event catalog, 191 events) also show no change, but the power there is poor
    (MDD ~110 km, 0.25 in direction fraction), so this does NOT test Meng et al. 2026's claim; it does fix the tremor-slip relative
    geometry to within one location error. The back-slip control DOES drift (+2.9 km/yr) from network growth — the artefact any
    drift claim must beat. (S1, S2, S4, S13)
R4  Slip and tremor are SYNCHRONOUS at the daily scale: the vbICA slip-rate track and the tremor track agree in timing to
    |lag| <= 5 d in 63 % of events (18.5 % under the null, p 4e-7), and in propagation direction in 78 %. The model has no
    latency relative to tremor; the ~35-d "tremor first" seen in catalog start dates is the trailing-window threshold. (S4)
R5  Potency per tremor detection in confirmed events is stable over 16 years (log10 Q5 = -4.1..-4.7, a factor ~3, no trend),
    while the back-slip control scatters over 3 decades; along-strike, slip-per-tremor has a persistent pattern (peaks
    42.8-44.3N and 46.3N, trough 47.5-49N) that is the same in all three epochs. Caveat: tremor-days per unit time grew x1.3
    in the north over the record, so "stable" means to within that factor. (S1, S3)
R7  Recurrence phase is stable: the lag between tremor-episode onset and slip-patch start shows no drift 2010-2026 (0.00 d/yr,
    p 0.61); band-level association is real but modest against chance (79 % of episodes matched vs 61 % null; 3 of 20 bands
    survive Holm, all north 46-47.5N). (S5)
R6  The model places 35-65 % of its positive slip outside event windows, but that slip is featureless: no inter-event slip
    burst north of 44N survives once real ETS days are excluded, and none is tremor-confirmed above chance. The geodesy resolves
    no inter-ETS creep that tremor can confirm at the 1-deg-band scale. (S6-A, definition A1)

## Round 2 — mechanism-oriented results (S7, S9; S8 pending)
R8  Within an event, the sparse waning TAIL of tremor is deeper than the main phase: last-tenth depth +4.0 km across-strike
    (p < 0.005), +2.4 km after latitude-detrending (p 0.005); the early up-dip migration of Wech & Creager 2011 is neither
    confirmed nor contradicted (-1.6 km, p 0.22). Tremor is slightly shallower when most intense (19/27 events; survives a spatial
    control). Detection geometry in the sparse tail cannot be excluded -> supporting observation, not a headline. (S9, S9b)
R9  Tremor is a SUB-LINEAR function of slip rate: daily tremor count ~ (modeled potency rate)^0.40 [0.31, 0.45] (27/27 events
    below 1), radiated tremor energy ~ rate^0.64 (2017+, 18 events, 14/18 below 1; range 0.3-0.8). Identical in north and central
    and in every epoch. Robust to window saturation, network growth, station count, and smoothing (the smooth-rate fit is
    unbiased; S9's "upper bound" caveat was reversed by simulation). Energy per detection rises with slip rate, so counts
    under-represent strong tremor, yet energy remains sub-linear. MODEL-FREE confirmation: with raw GNSS displacement rates as
    the regressor (no slip patterns), b(count) = 0.28 [0.10, 0.51], b(energy) = 0.31 [0.04, 0.71], statistically identical to the
    model-based values at matched bandwidth (paired p 0.92), b = 1 excluded. Physical content: tremor productivity per unit slip
    falls as slip accelerates (10x slip rate -> ~2.5x detections, ~4.5x tremor energy). No published exponent exists for daily
    tremor vs slip rate within events; the field assumes proportionality (S10). CLAIMABLE. (S9, S11, S11b)
R10 The modeled slipping band brackets the tremor band: slip up-dip edge ~11 km shallower than the concurrent tremor up-dip edge
    in 27/27 events, down-dip edge ~5 km deeper in 22/27 (upper bounds, broad patterns); no edge drifts (bounds 0.03-0.6 km/yr).
    In the north, a patch's tremor is indistinguishable from the whole northern band on the same days: the north is ONE tremor
    band, not a set of local patches. (S9, S1b)
R11 Central Oregon (44.3-45.4N, 27-38 km): tremor persists over triangles whose modeled slip is bottom-quartile in every epoch,
    and slip-tremor co-location is weakest there (rank rho 0.40-0.46 vs 0.70-0.77). Not a coverage/resolution artefact (S12) and
    not a temporal-basis artefact (S12b): the nearest GNSS stations see a tremor-synchronous transient (10 +- 3 mm composite) of
    which the model captures 63 % (north 77 %), residual 1.2 sigma; unmodelled slip bound 6 mm/30 d on a 30-km patch. What the
    model withholds is AMPLITUDE: potency rate per km^2 is half that of the adjacent slip+tremor patches, a rank statement at the
    resolution limit. Read as "tremor with ~half the slip per area", not "tremor without slip"; central Oregon is also the
    least tremor-productive segment (Wells 2017). Persistent slip-without-tremor: none north of 42N. Modeled slip extends ~5 km
    below the tremor band in the north (S7; cf. Hall et al. 2018 for the up-dip side). (S7, S12, S12b)
R12 No persistent along-strike barrier survives correction (best candidate 48.88N, p_Bonf 0.09). Modeled slip extents are wider
    than tremor extents (89 % of events) and run ~70 km further south (model smoothing), so slip does not stop where tremor
    stops in this model. Termination sites do not depend on propagation direction. (S7)

## Round 5 — aggregation match and paper deliverables (S15, S16)
R13 The Frank & Brodsky (2019) amplitude-slip relation REPRODUCES in Cascadia in sign and sub-linearity, and its VALUE
    reproduces once the aggregation granularity and the estimator are matched. Their statistic is the MEDIAN per-detection
    amplitude in a time bin (beta = 3.17 +- 0.29, i.e. amplitude ~ rate^0.315). Computed the same way here on the 18
    amplitude-era events (2017+), with our windows, the R7 regressor, Theil-Sen fits and the >= 365-d shifted-tremor null:
    A1 = sqrt(median daily energy) ~ R7^0.121 [0.032, 0.151] at DAILY resolution (18/18 events below 1, 14/18 above their own
    null, sign p 0.015), rising to 0.166 (3-d bins), 0.227 [0.177, 0.436] (7-d bins — CI contains 0.315) and 0.296 with one
    point per event, the granularity closest to their 8 SSEs. The estimator matters by a factor ~9 because the per-event
    r^2 is only 0.11: the direction-symmetric RMA slope is 0.271 (beta_RMA 2.54, bracketing 3.17) while F&B's own regression
    direction gives beta 1.01 (all 18 events below 2.88). Read for the paper as "consistent with Frank & Brodsky at matched
    aggregation, with the estimator stated", never as a single number. Two things must NOT be said: (i) that our energy
    exponent 0.64 "equals" their 0.32 via 0.64/2 — sqrt(daily summed energy) IS b(E_sum)/2 = 0.344 identically, an
    arithmetic identity, not an amplitude measurement (the LITERAL summed amplitude gives 0.561 = count x amplitude);
    (ii) that their 0.32 agrees with our COUNT exponent 0.40 — on the same 18 events the like-for-like amplitude exponent is
    0.26 BELOW the count exponent (paired -0.262 [-0.292, -0.155], p 1.5e-05), and the three exponents are rank-correlated
    0.88-0.92, i.e. three readings of one relation. Decomposition: count R^0.37 + per-detection amplitude R^0.12 -> summed
    amplitude R^0.56; count R^0.37 + per-detection energy R^0.23 -> summed energy R^0.69. The whole daily energy
    distribution shifts with slip rate (p10 0.151 ... p90 0.257), and the >= 8-station control LOWERS the amplitude exponent
    (-0.031, p 0.007) while it RAISED the count exponent, so it is not a detection-threshold artefact. Caveat: PNSN
    `energy`/`magnitude` are ONE observable (log10 E = 1.781 M, residual 0.052 dex) describing a 150/300-s window, i.e. a
    windowed envelope amplitude, not F&B's per-LFE displacement amplitude; a strict comparison still needs a per-LFE
    amplitude catalog. (S16; results/s16_fb_summary.md, s16_fb_aggregation.csv, s16_fb_events.csv)
R14 Paper deliverables exist: Figs 1-4 (figures/s16_paper_fig1..4.png), Table 1 (the 27 events, results/s16_table1.csv),
    Table 2 (every exponent with CI, null and provenance, results/s16_table2.csv) and caption drafts
    (results/s16_captions.md), all built from the saved S1-S16 CSVs. The only quantity that had to be regenerated is the
    per-event DAILY (R7, count, energy) tracks, which no step had saved (scripts/s16_daily_tracks.py ->
    results/s16_daily_tracks.csv); it reproduces every stored per-event exponent to 8e-17
    (results/s16_track_validation.csv). The S15A snapshot check and the count exponent are unchanged (S15).

## What does NOT hold (and why)
X1  "Tremor efficiency changed, especially in central Oregon (x0.64)." The margin-wide change flips sign with the network
    normalisation; in tremor-days the decomposition shows the modeled slip rate is flat in every segment and 100 % of the change
    is tremor-day growth (central x1.92 vs north x1.30 per unit time, threshold-limited). It measures PNSN growth in Oregon
    relative to Washington. (S3, S3b)
X2  "Southern tremor shallowed by 1.7 km." Tracks completeness 0.63 -> 1.5. (S2)
X3  "Tremor-confirmed inter-ETS slip (p 0.005)." Artefact of defining event days from the strict confirmed set; vanishes under
    A1/A2. (S6-A first pass)
X5  "ETS recurrence shortened by ~30 d/yr in bands 43.0 and 44.5N." Sits in the fastest detection-growth bands and vanishes under
    the rate-normalisation control (p 0.5-0.7). Slip patches also "recur" 2-4x faster than tremor episodes because the 5-mm
    catalog fragments events (S5); not physics.
X4  Along-strike offsets from union footprints (Q2) and the slip-centroid propagation SPEED (v_s ~ v_t/6): mixture and chaining
    artefacts of the 13-fixed-pattern model; not physical. (S1, S4)

## Open items
O1  RESOLVED (S1b): the 2022 outliers (tremor 20-29 km down-dip / 3.5-5.4 km deeper) are window selection (short threshold windows
    that caught the deepest northern tremor days; relative to concurrent north-wide tremor they are normal, -0.5..+0.1 km), an
    unusual vbICA component mix (non-k1/k2 events: Q1 19.6 vs 4.1 km) and small size. Not a geometry change. Lesson: reference
    depth offsets to the concurrent north-wide tremor depth.
O5  Central Oregon tremor-without-slip: needs GNSS-coverage check (station density vs the model's resolution there) before any
    physical reading; natural target for the dv/v stage.
O3  Any efficiency-change question needs a fixed station subset or an amplitude-thresholded tremor measure.
O4  dv/v bridge deferred by the user.

## Publishability memo (orchestrator, 2026-09-10; S9b and S13 still pending)
LEAD RESULT (GRL-scale): "Tremor is a sub-linear function of slow-slip rate in Cascadia, and the relation has not changed in
16 years." Daily tremor count ~ slip-rate^0.3-0.4 and tremor energy ~ slip-rate^0.6-0.65 across 27 ETS events (north + central),
identical in both segments and in 2010-14 / 2015-19 / 2020-26, obtained with a 13-component vbICA slip model AND with raw GNSS
displacement rates; b = 1 excluded; saturation, network growth, station count and smoothing excluded (S9, S11, S11b). Literature
assumes proportionality; no exponent has been published for daily tremor vs slip rate within events (S10). Supporting results
for the same paper: co-location within one mesh cell and day-scale synchrony (R1, R4; cf. Hall et al. 2019 on sign), and the
stationarity bounds (R3, R5) with the S13 power table — frame as "relative geometry fixed to within one location error"; do NOT
frame as a refutation of Meng et al. 2026 (nucleation power is ~110 km).
Interpretation to develop (carefully): fewer detections/less energy per unit slip at high slip rate means the tremor-generating
population saturates or that fast slip is accommodated increasingly aseismically — compare with Frank & Brodsky 2019 (per-LFE
amplitude ~ rate^0.32 in Guerrero) after matching aggregation (their per-LFE median vs our per-day sum), and with rate-and-state
expectations for a tremor asperity population.
TO DO before drafting — STATUS 2026-09-10 after S15/S16: (b) DONE (R13); (c) DONE (R14: Figs 1-4, Tables 1-2, captions); (d) DONE, null result (S15A: the shipped A_raw IS the 2026-08-18 solution, bit-identical in float32 — the S12b 'stale snapshot' caveat is WITHDRAWN); (e) DONE (the stationarity rows of Table 2 plus the S13 power numbers in R3); (a) STILL OPEN. Original list: (a) a third tremor measure from a PROPER LFE catalog (S15: the picker-ranked candidate counts are not independent of PNSN, rho 0.97, and reproduce the energy exponent 0.68 — useful as a consistency check only) (S14 showed the project's densify detections are
false-alarm-dominated stacking material, ETS/non-ETS ratio 1.2, and cannot give an exponent — a published Cascadia LFE catalog or a
MAD-thresholded, declustered re-detection is needed); (b) match Frank & Brodsky's aggregation; (c) per-event figure set: 27 event panels of log tremor vs log slip rate; the exponent-vs-epoch/segment panel; the
model-vs-GNSS agreement panel; (d) re-derive A_raw from the 2026-08-18 solution (snapshot mismatch flagged by S12b) and confirm
the 27-event set and exponents are unchanged; (e) power/bounds table for the stationarity claims.
SECOND RESULT (section, not paper): central Oregon slips with its tremor but with about half the slip per unit area of its
neighbours, at the resolution limit; not a coverage artefact (S7, S12, S12b, S10b).
DROPPED: efficiency change (network growth), predictability (published, mislabeled), within-event deepening unless S9b reconciles
it with Wech & Creager 2011.
