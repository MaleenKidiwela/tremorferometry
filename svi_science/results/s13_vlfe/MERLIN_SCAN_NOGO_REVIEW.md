# MERLIN — scan pilot NO-GO review and second-pilot design (2026-10-15 ~0:20 am PT; saved by main session, condensed)
Verdict: salvageable cheaply. Detection works (inj Mw 3.6 median z 11, 3.8 z 35, location 6 km; H&G 2,4 z 18); the background is
network-coherent travelling long-period waves that a 5–7-station 100-km network cannot tell from a near thrust pulse.
Route: re-score the SAME 2011 run with (1) far-reference travelling-wave veto (stations 100–700 km, same band), (2) best node
inside the LFE footprint, (3) single-station step cleaner instead of hard LOSO, (4) corrected earthquake masks as columns.
Not empirical templates (template noise never the problem); not vertical-only (a Rayleigh wave is a uniform uplift).
Flaws: "60 % in m45 mask" = base rate (69 %); only 11/47 survivors credibly M 5–5.5 body waves (S1 masked only surface-wave
windows); 36 unknown origin, 31/47 on grid edges, 3–8-min pairs → waves from outside (M3–4.5 regional, unmasked); hard LOSO
halves completeness (Mw 3.6 0.78→0.33); LRIV dominant in 26/47 (gain 1.46, glitchy); go/no-go items H&G ≥ 3/4, sign ≥ 10×,
ETS density ≥ 10× were unachievable/meaningless (H&G 3,5 below completeness; antisymmetric template; 2–4 ETS events) — replace;
HF-quiet must be a class, never a cut (H&G 2 is HF-loud with LFEs).
Second pilot (PLAN_SCAN Amendment S7 before any number):
1 references: on disk GNW (full 2011), ETS-days set (OBSR, LEBA, RATT, ELW, D04D, WISH, SP2, B05D, NLWA, D03D); ASK USER to download
  1-sps LHZ 2011 for ~5 stations at 300–700 km (IU.COR, US.NEW, US.HAWA, CN.PNT, CN.BBB; ~30 MB each). GNW quality check.
2 measure VLFE far-field decay with the S13 stacker at GNW/DOSE (2011-08) and 2012 full-year set (ratio ≤ 0.5 at ≥ 100 km).
3 veto columns: ref_excess, ref_ratio, plane-wave VR/slowness, footprint-interior; validation: H&G 2 & 4 survive, ≥ 9/11 known
  body-wave clusters flagged.
4 single-station step cleaner (|z_s| ≥ 8 with others < 2.5 → NaN ±600 s), re-run detector (~40 min cupy); (a) ≤ 20/yr at z ≥ 7,
  injections Mw ≥ 3.6 keep ≥ 85 %.
5 mask corrections as columns (body-wave windows; Ms = mb − 0.3; ComCat M ≥ 3.0 within 2,500 km 'reg3').
6 threshold pair (z × z_loso) maximising Mw-3.8 recovery s.t. (a) ≤ 0.3/yr, (b) ≤ 0.7/yr.
7 GO: (a)+(b) ≤ 1/yr between ETS (strong ≤ 0.5); ≥ 80 % at Mw 3.8, location ≤ 15 km; H&G 2 AND 4 survive; ETS ≥ 3 detections and
  Poisson p < 0.01 vs between-ETS rate. Classify the 2011-08-18 21:19 UTC cluster (z 42.9, Mw 3.86, N grid edge 49.24 N, no
  catalogue/H&G match) — if local, a candidate northern VLFE; extend grid to 49.6 N.
8 cost ~1 day, ~1–1.5 h compute.
If still failing: stack at candidate times (P 0.5–0.75 / 0.3–0.5 / 0.1–0.3 between-ETS candidates; PNSN tremor windows with no
catalogued LFE; BE_broad anchor) → between-ETS VLF:LFE ratio per class vs ETS 12 (cannot see VLF with no 2–8 Hz activity).
