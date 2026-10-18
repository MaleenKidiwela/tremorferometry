# 2026-09-02 — MERLIN RULING: 3-5 s window "full inversion" request (model: Claude Fable 5.1, confirmed)

VERDICT: 3-5 s is buildable from the stacks on disk and worth a PILOT, but not for "multiple scattering" — both windows
sit in the same mixed-scattering regime (lapse/mean-free-time 1.4-3.1 at ell=40 km; 1.6 s apart on a 15-35 s lapse).
What 3-5 s buys: (i) first window free of the 0-2 s template imprint + S-tail (2-4 s window energy is 70-96% in its
first 0.5 s), (ii) a DECISIVE LAPSE-PROPORTIONALITY TEST on the 47.2N decline. The decline (~2%) = ~1 sample (25 ms)
of alignment drift at tau_bar 1.3 s -> first hypothesis to confront (Flaw 4).

Q1 measurable: YES marginally. 30-d-stack RMS(3-5)/RMS(8-10): HOQUI 1.92, YELM 1.93 (88-92% of certified fams >=1.41),
  PR02 7.4. Gate = per-family median amplitude ratio >=1.41 AND odd/even split-half corr >=0.5 (>=2 yr) AND injection
  recovery >=0.7. Do NOT gate with the causality/mirror ratio (mirror window is coherent pre-S energy 6-13x floor).
Q2 kernel: SAME kernel form for both windows, built in the S-PINNED DIFFERENTIAL frame from path RESIDENCE
  (K_eff = [t*K_coda - t_S*K_ray]/(t-t_S)), tau_bar 1.3 s (2-4) / 2.9 s (3-5), ell in {20,40,80}; MC-RTE random walks
  cheapest honest implementation (all orders). Never "single-scatter for 2-4, RTE for 3-5".
Q3 anchor: production IS S-pinned (fleet_station.sh --origin-anchor; t_s=+1.0 s). June _calT 3-5 files are t=0-anchored
  -> NEVER reuse. Predicted eps35/eps24: homogeneous 1.0; shallow ~0.95; constant time shift ~0.45; slab-confined ~undefined.
Q4 product: replica RESDIR res_catalog_3to5 OK but LAST; the decisive product is a DATA-LEVEL paired test (no inversion).
  Joint multi-window inversion: NO for now (3-4 s shared data, window-dependent site terms, shared frame error).
Q5 gain: a diagnostic ("decline is/isn't lapse-proportional"), NOT a deeper map, NOT independent replication.
Q6 cost: dv/v 1-3 min/station certified-only @28 workers -> 69-stn core roster 2-3 h; census <1 h; kernel diag ~half day.
ORDER: (0) differential-kernel diagnostic -> (1) split-network test (queue #1 unchanged) -> (2) 3-5 s census + pilot
  on the 47.2N core roster -> (3) paired ratio / shift+stretch test with injection calibration -> (4) fleet chain only if
  proportional AND corrected G exists.

FLAWS (his, severity-ordered; magnitudes from HIS coarse MC = provisional, verify in Step 0 — don't-over-conclude):
 F1 FATAL-for-interface-amplitude/sign: G not in the S-pinned differential frame. Direct ray and coda share the slab
    fraction (~0.06); residual amplified by t_S/tau_bar ~8-18 -> S-pinned on-interface sensitivity ~+0.1 (34-38% of pairs
    NEGATIVE) vs naive 1.14/0.47; top-5-km sensitivity ~1.1 never negative. Checkerboards/T/injections can't see it
    (same G). "Deep stable" bound ~10x looser than stated. Condition that would collapse it: stacks NOT aligned on S —
    CHECKED 2026-09-02: build_long_window_3comp.py aligns on the Z matched-filter detection time -> S pinned -> stands.
 F2 MAJOR: coded "single-scatter ellipsoid" = scatterer (kink) density, not travel-time RESIDENCE (theory 01 eq 7-8 vs 15);
    source cell gets zero weight; captured fraction f=0.047 is a kink fraction.
 F3 MAJOR: window bookkeeping 1 s off (build_G W1,W2=2,4 past DIRECT arrival; data 2-4 s stack time = 1-3 s past S);
    energy centroid 1.3 s past S (not midpoint). G.npz is geometrically the matching shell for 3-5 s, not 2-4 s.
 F4 MAJOR: 2% decline = ~1 sample; trigger-floor-dominated stacks (~200 noise triggers/day/family); slow noise-trigger
    fraction drift changes stack shape vs all-time ref -> reads as stretch.
 F5 minor: mirror window is coherent (certification stays conservative; not a 3-5 s gate).  F6: June 3-5 evidence unusable.
 F7 traps: build_G_captured reads pairs.src_lat/src_lon/src_dep — pairs.csv on disk LACKS them (use G.npz cxy/depth);
    assembler hardcodes _Z_2to4 (needs WINDOW param); cycle-skip at 3-5 s: eps_max 0.03 safe (0.05 skips); SVD-Wiener
    shrinks low-SNR rows toward ref -> uncalibrated ratio <1 uninterpretable.

PRESCRIPTION (accept / failure signatures):
 S0 differential residence-kernel diagnostic: new script under fault_tomography/kernels/ (don't modify kernel.py); per
    production pair on the build_G volume grid (0.2 deg x 1 km, IFACE +-1.5 km): residence of singly-scattered paths
    (shell scatterers, both legs, image source), direct-ray residence, K_eff for tau_bar 1.3/2.9 s, ell {20,40,80};
    project to 872 cells -> res_catalog_kernels/G_diff_w24.npz, G_diff_w35.npz. Report per row on-interface net, sign,
    source-cell weight, top-5-km net. ACCEPT: reproduces his MC within 2x (interface net ~0.1, 30-40% negative; z<5 km
    ~1.1). FAIL (report, don't proceed to S5): >50% rows negative or median |net|<0.05 both windows -> reframe as
    volumetric/shallow+regional inversion.
 S2 census: scripts/census_window_snr.py -> data/window_snr_census.csv (tag,fam,snr24,snr35,snr45,nrows) + grand-stack
    RMS-by-0.5-s-bin profile. FAIL: fleet median snr35<1.2 with <30% fams >=1.41 -> coda-rich stations only.
 S3 3-5 s dv/v on 69-stn core roster (patterns[16]>0.5): dvv_roll30cal.py --window 3 5 --origin-anchor --eps-max 0.03
    --n-eps 301 --cert-csv ... --workers 28 -> data/daily_dvv_<TAG>_Z_3to5.csv; <=2 stations concurrently. FAIL: rail
    fraction >20% at a station -> exclude it.
 S4 per-window injection calibration (~10 stns x 5 fams, S-pinned stretch +0.01 over 2 yr, differential protocol) ->
    r24, r35 vs snr35. ACCEPT r24>=0.9, r35>=0.7 at snr35>=2, monotone. FAIL r35<0.5 -> filter-dominated, stop+report.
 S5 paired lapse test: per family D_w = mean(eps_w 2021-22) - mean(eps_w 2014-16), corrected by r_w; Theil-Sen slope of
    D35 vs D24 + 30-d block bootstrap; also cell-aggregated; optional 2-param (shift dt0, stretch) grid on 2-5 s.
    PREDICTIONS: proportional slope ~1.0 (0.9-1.05); constant time-shift ~0.45; scatter -> composition change.
    ACCEPT CI width <0.2 + same verdict both levels. FAIL CI spans 0.45 and 1.0 -> widen to 102-stn support set.
    OUTCOMES: ~0.45 -> decline = ~1-sample alignment/composition drift, existing 2-4 s slow-field results contaminated
    (becomes the priority finding). ~1 -> real lapse-proportional change; given F1 shallow-vs-slab 1.1 vs 0.1, suspects =
    shallow regional / instrument response / processing, NOT the megathrust; only then S6 with G_diff.
 S6 fleet 3-5 s chain replica (only if proportional): RESDIR res_catalog_3to5, G_diff_w35 as G.npz, assembler WINDOW env,
    re-oracle lambda. FAIL: deep transmission <0.1 -> publish only the S5 data-level result.
NEVER: overwrite res_catalog/G.npz, any _Z_2to4*.csv, family_days.parquet; reuse June _3to5_cal*; mixed kernels.
