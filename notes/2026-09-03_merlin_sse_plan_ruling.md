# 2026-09-03 — MERLIN RULING on the six-stage SSE-model vs db/b comparison plan (Fable 5.1 confirmed)
GOAL (user /goal): find a pattern in how SSE model variables may relate to db/b.
VERDICT: viable as a DIFFERENTIAL, episode-relative comparison, with two gates: (i) admit episodes only when
TREMOR-SUPPORTED; (ii) test stress/strain only as FOOTPRINT MEANS, never per cell. Reduced version = three
pre-registered numbers + one max-statistic exploratory scan; defensible in one day.

FLAWS (his computations): F1 MAJOR model episode catalog in 42.5-45N is NOT slow slip (tremor coincidence 29% = chance;
 winter phase-locked R=0.44; sign-symmetric rate) -> hydrological/seasonal leakage through vbICA; N 78%, C 80%, S 55%
 supported. F2 MAJOR per-cell stress at 10 km above the interface = mesh-roughness noise (corr(shear change, local slip)
 POSITIVE +0.45 inside a slipping patch; normal stress sign-less); footprint means keep the sign (shear median -0.73 kPa
 at 10 km, 73% negative, Sep 2024). F3 MAJOR earthquake/glitch offsets inside amplitudes (comp 5: M7.0 Cape Mendocino
 2024-12-05..07 jumps; others 2015-08-29, 2017-08-06, 2024-09-15, ...). F4 MAJOR lag range aliases recurrence in the
 south (median 235-242 d S; 362-374 d N/C); lags <15 d unidentifiable; db/b persistent component (acf 0.38 at 90-180 d).
 F5 MAJOR per-episode noise 2-10x expected signal (30-d change SD deseasoned deep cells 0.5/0.8/1.3%; spatial corr 0.17
 at 50-100 km, ~0 beyond; footprint of 300 km ~4-8 independent cells -> footprint-mean 30-d noise ~0.3%, 20-episode
 composite ~0.07%). F6 minor stored trailing-frame stress contrasts ~7x weaker than delta-A_raw (must reconcile).
 F7 minor 89 cells >15 km from nearest patch (exclude). F8 minor rectified level invariants (eps_eq, max_shear, von
 Mises, max_shear_stress) are GAUGE artifacts (V-shape at the 2019-24 mean crossing) — the -0.18 corr in the 47.6N cell
 figure is that class; use increment-tensor invariants or drop. F9 spec: 13-rank truncation -> model propagation is
 superposition arithmetic; take propagation from tremor. NOT a flaw: strain sign (slip_enu rotates negative dip-slip
 up-dip; thrust gives negative resolved shear); receiver rake 90 convention; mask/null/controls carry over.

RULINGS: Q1 threshold 3x negative-side RMS on footprint-mean raw increment (2x/4x sensitivity); ADMISSION GATE: trailing
 30-d tremor count in the episode's latitude band within +-15 d of the model peak in the band's top decile; segment
 with admitted fraction <0.40 gets no primary test (S-C fails; S borderline). Min footprint 5 cells (patch_dist<=15 km);
 min duration 10 consecutive days of raw 7-d rate above threshold; overlaps: space-time linkage 100 km / 45 d, one
 episode per component, segment assignment >60% moment. Onset = LOCAL onset from the RAW 7-d rate (recover 13-vector
 per cell by regressing stored slip_thrust_cum on A_trail, apply to A_raw; onset = first day rate > 25% of episode
 peak), +-5 d; local-onset alignment PRIMARY, episode-onset SECONDARY.
 Q2 identifiable lags -90..+150 d N/C, -60..+120 d S, with EXCLUSIVE masking (drop lag if within 45 d of another
 admitted episode's onset-end in the same segment, applied identically to shifted catalogs); disqualify features at
 |L| > 0.5x median recurrence; the scan is exploratory.
 Q3 PRIMARY (p<0.01 on BOTH nulls): P1 footprint response = mean ON [+15,+60] - PRE [-90,-30] (shifted +15 d for the
 trailing lag) of episode-masked-deseasoned footprint-mean db/b, local-onset aligned, pooled over admitted episodes
 (sign count + mean); P2 locality = P1 minus the same on the RING; P3 stress functional = pooled fit of footprint-mean
 30-d db/b change on the footprint-mean 6-component strain increment (one shared 6-vector) + footprint-mean slip
 increment, variance explained. Null A = group-rigid circular shift; Null B (NEW, mandatory) = integer-year shift
 (1..14) + uniform +-15 d jitter (preserves seasonal phase); 4000 draws each. Everything else EXPLORATORY with a
 max-statistic FWER (record max |stat| over the whole scan per null draw). Variables 24 -> 12: 6 strain, dilatation,
 local-orientation shear/normal/Coulomb, slip increment, cumulative slip.
 Q4 PCA doable ONLY on footprint-mean responses in 30-d lag bins (-90..+150 = 8 bins), each episode demeaned by PRE
 and divided by its own PRE-window SD (removes roster heteroscedasticity); tests: leading eigenvalue vs Null-A
 eigenvalues (p<0.01) AND split-half (odd/even episodes within segment) mode-shape corr above the null's 99th pct.
 Drop episode clustering (per-episode SNR<=1) and cell-resolution PCA/ICA.
 Q5 per-cell lasso MEANINGLESS (all level variables in the 6-dim span of the cell's strain; everything in the 13-dim
 span of A_trail). Correct: ONE shared linear functional of the strain tensor per segment (6 coefs) + local slip-
 increment coef + footprint-mean slip-increment coef + per-cell intercepts, OLS on 30-d increments, H0 vs Null A;
 interpret the 6-vector by its ANGLE to the named functionals (Coulomb, normal, shear, dilatation, max-shear). Fit on
 footprint means or 100-km-smoothed cells, never raw cells.
 Q6 ring = cells >=100 km from every footprint cell, <=300 km, patch_dist<=15 km, not within +-45 d of an admitted
 episode in the ring's own segment; calibrate with Null A (pre-window footprint/ring corr |r|<0.15 across draws, else
 widen to 150 km); require >=10 ring cells spanning >=100 km else exclude from P2.
 Q7 carry: Null A, four controls (cc, rows/day, SITES footprint stations, day-of-year R), patch-scale aggregation,
 2011-04..2026-04. Deseason mandatory but EPISODE-MASKED (climatology from days outside [-30,+90] of any admitted
 episode in the segment; month with <40% days retained -> segment seasonally inseparable). Family-drift gate NOT
 required for PRE-baselined contrasts; run P1 once on splits/gateD0.3 as sensitivity. Glitch list: days |dA_k| > 6x
 robust SD + fixed windows 2024-12-05 +120 d, 2014-03-10, 2021-12-20, 2022-12-20, 2012-10-28 each +60 d; exclude
 episodes whose ON touches a glitch day of a component carrying >30% of footprint slip. Receiver: LOCAL patch strike/dip
 (median 333/11.4 N; changes shear medians 2x). State fields: stored trailing frame + 30-d differences are already
 matched; raw-increment rebuild only for RATE variables (slip) — do not "fix" the state fields.
 Q8 DROP: model-derived propagation; 24-var list; per-cell lasso; episode clustering; response-vs-propagation/depth/
 roster regressions (keep moment + segment, exploratory); size/depth composites as primaries; S-C from any primary.
 ONE-DAY VERSION: admission + glitch gates on a segment-level catalog (~55 admitted N/C/S), local-onset alignment,
 episode-masked deseason, P1/P2/P3 with Nulls A+B + four controls, one max-statistic exploratory table.

PRESCRIPTION (order): (1) admission + glitch gates -> accept N,C retain >=70%, >=40 episodes total, day-of-year R<0.3
 per admitted segment; <30 episodes -> insufficient, stop, do not pad with S-C. (2) reconcile F6 BEFORE any stress test:
 recover 13-vectors per field, compare Sep-2024 ON-PRE from delta-A_trail vs delta-A_raw vs direct Green's re-eval;
 delta-A_trail must equal stored to float32; ratio to delta-A_raw must match the slip-series ratio within 20%, else
 build error -> report. (3) footprint-mean predictors with local orientation at 0/10/20 km; accept: inside-footprint
 shear-change medians negative at all offsets for >=80% of admitted N/C episodes; <60% -> stress untestable, slip-rate
 hypotheses only. (4) episode-masked deseason + response cube with exclusive lag masking; per-lag N >=15 N+C pooled.
 (5) P1/P2/P3 + Nulls A,B + controls; DETECTION = p<0.01 both nulls, same sign N and C, every leave-one-segment-out
 p<0.05, all controls p>=0.05; locality rule (local-onset peak >30% above episode-onset peak AND P2 detected -> local;
 P1 not P2 -> regional; P1 fails -> 95% bound from Null A). STOP-ALL: day-of-year control R>0.3 or SITES detects.
 (6) exploratory scan + PCA only after (5), FWER + eigenvalue/split-half nulls; sensitivity map at 100 km only.
COULD-NOT-RUN: S segment one vs two populations (peak >=8 mm vs <8 mm coincidence); S-C pseudo-slip vs precipitation/
 snow index; point stress vs kernel-weighted volumetric strain (average strain over G.npz kernel weights).
