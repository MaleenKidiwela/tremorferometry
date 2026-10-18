# 2026-09-02 — MERLIN RULING: split-network coherence test on the 47.2N deep decline (Fable 5.1 confirmed)

VERDICT: design (A) station-split is UNDERPOWERED and tests the wrong thing; do NOT run first. Merlin read the raw data
(no re-solve): the q16 (deep 47.2) level change is carried by families whose station-common-removed decade-scale level
changes are 0.3-3% in eps (10-100x what a deep-interface signal can produce at 5% capture), concentrated at B017
(post-2020 collapse), B023, B927, KSXB, STOR, LEBA. With that population robustly handled the q16 level change is 0 to
-0.9 patch-units at 1-3 sigma (station-cluster), and its shape is a 2020->2021 STEP, not a 2017-22 decline.
THE decision: gate family-differential slow drift BEFORE any slow-field claim.

FLAWS: F1 FATAL — slow-claim noise model wrong by ~10x: for 3,022 family-eras with >=60 rows in both W1=2011-16 and
 W2=2022-26 (50 stns), station-daily-family-mean removed, D=mean(W2)-mean(W1): median|D| 0.23% (bh 0.19, bb 0.28), p90
 0.74%, 1207/3022 >0.3%. A -2% q16 field predicts |d eps| ~0.08% (median fcap 0.049). Regression on the 20-patch basis
 w/ station FE: q16 -2.45 +-0.57 (model SE) but station-JACKKNIFE SE 1.61; clipped 0.3%: -0.86+-0.58; excluding |D|>0.3%:
 +0.24+-0.78; orig-file stns only: -0.72+-1.91; within-station label permutation p=0.16. sig_agg (0.036%) implies
 ~18 sigma; empirical station-cluster noise makes it ~1.5 sigma. => sig_agg/sig_frame are NOT the noise for slow
 features (incl. deep-index trend -0.037%/yr + "claimable tile" roster).
 F2 MAJOR — B017 pathological post-2020, no era boundary/veto: station-mean eps -0.5 (2013-18) -> +0.05 (2020) -> +1.62
 (2021), +1.33...; rows/yr 18k -> 3-4.6k; cc 0.92-0.99 -> 0.71-0.93; per-family W2-W1 spans -0.89..+3.07% (4% spread in
 ONE station). B023 similar (drift 0 -> +0.5 -> -1.15%, family spread 3.2%).
 F3 MAJOR — fixed-roster (bridge-only) replicate shows +0.4 (2011-17) -> +0.26 (2018) -> +0.17 (2020) -> -0.61+-0.18
 (2021) -> -0.3: a 2020/21 STEP; the product's 2016-20 ramp is roster-shaped (WISH/RADR/SP2 ending 2019-05, B017
 pre-collapse). Per-station sign support for the step: B017/B023/B201 (+), marginal B943/B005, ~0 B013/B941/STOR/SWF2/
 LEBA/NLWA, RAIN/B011 flipping.
 F4 MAJOR — the 08-31 era whitelist (visually confirmed STATION-common step) is blind to the harmful FAMILY-DIFFERENTIAL
 step (freq-dependent response change x family spectra). RAIN boundary 2022-06-09 (not whitelisted): support families
 +0.09/+0.13 -> -0.04..-0.14 while non-support go the opposite way.
 F5 MAJOR — (A) simulated (300 stratified halves): half estimates -0.50+-0.37, both-negative 87%, corr(A,B) across
 splits = -0.47 (complementary halves anti-correlated by construction) -> corr-between-halves is an INVALID statistic.
 F6 minor — YACT boundary 2020-10-17 whitelisted -> identifies nothing across the transition; B201 sparse 2024-26.
 F7 minor — BCD not converged on slow modes (dc/sig 0.147 > 0.1 at sweep 30); subset reruns must use same sweep count.
 F8 minor — deconvolved A16 is a high-gain contrast (1.89 y16 - 0.77 y6 - 0.80 y17 ...); run coherence tests on RAW
 projections y16 (and y6+y16 column sum), not A16. Daily T for q16 flat 0.79-0.82 across the transition (not a mechanism).

RULINGS: Q1 order: (0) family-drift census + gate (no rerun) -> (1) gated joint rerun -> (2) design (E) family-split, one
 split = empirical slow-claim noise -> (3) station leave-one-out (top 5) -> (4) design (A) only if a signal survives.
 Decisive = (1)+(3). Designs (B),(C),(D) dropped (subsumed). Q2 (if A): unit station; strata bridge/file-type/borehole;
 each half >=8 bridge stns; N=6 splits. Q3 statistic: L = mean(2022-26) - mean(2011-16) of each half's raw T-normalized
 y16 (NOT a slope); noise = SD over splits of (L_A-L_B)/2; bar |z|>=3; label permutation (>=500) at data level p<0.01;
 NO circular shift; NO per-split corr(A,B). Q4 compare raw y16/T and y6+y16; own deconv only SECONDARY. Q5 no re-gauge
 needed for L. Q6 lambda FIXED 0.464; failure: subset T16 <0.5 or 100-km checker <0.45 -> re-oracle + note.
 Q7 OUTCOMES: CANDIDATE = gated L <= -1.0 AND same sign in every single-station leave-out with |dL|<50% AND family-split
 halves same sign with |L_full| >= 3 SD AND (if A) |z|>=3 in >=5/6 splits -> next: end-to-end amplitude calibration,
 then lapse-proportionality. ARTIFACT = gated L within +-0.5 of 0, OR any leave-out flips sign / |dL|>50%, OR gate
 sensitivity (0.2/0.3/0.5%) changes L >2x or flips sign -> write-up: 47.2N deep interface velocity-STABLE within the
 empirical bound (2x station-cluster SE of gated L). AMBIGUOUS = else (e.g. L -0.7..-1.0 at 1.5-2.5 sigma) -> run (A)
 N=6; if |z|<3 report as bound with the 2020/21 step as unresolved.
 Q8 confounds split-by-station CANNOT touch (state all six): family-differential measurement drift as population
 property (30-40% of support rows from unif05/eps05 files); shared LFE sources (source-side coda change coherent across
 halves); shared processing chain (all-time ref, cc>0.7 gate -> biased survivors when cc degrades); kernel/model error
 (95% of coda sensitivity off-interface); non-whitelisted instrument changes (F4); regional environmental forcing
 (2015-21 drought, 2021 heat dome). Bookkeeping drift weak but present: Spearman |D| vs |dlog n_stack| 0.14 (p<1e-3).

PRESCRIPTION: (1) family-drift census (~1 h, no rerun): per family-era, subtract station daily family-mean; D where
 >=60 rows in both windows; 365-d rolling-mean range for all; flag |D|>0.3% or range>0.4%; table by station/src.
 FAIL: flagged fraction <10% -> reconcile before anything else. (2) station review of top-flagged + every
 non-whitelisted era_table boundary (support vs non-support +-1 yr; RAIN 2022-06-09 first); PROPOSE vetoes to user
 (B017 date veto from ~2020-07-01 or boundary; B023) — do NOT apply without user's word. (3) gated joint rerun
 (separate RESDIR, --lam 0.464, 30 sweeps, T + deconv), gates 0.3 (primary) + 0.2 + 0.5; report y16/T and y6+y16
 annual means + L; if day-gating removes >5% frames relax MIN_SITES to 8 and record. (4) design (E): odd/even
 family-eras within every station (sorted fam id), gated, fixed lambda, own T -> SD((L_A-L_B)/2) = empirical
 slow-claim noise; annual-difference variance vs averaging window (= the "split-half tau" the checkerboard ledger
 needs). FAIL: halves differ in sign -> bound write-up. (5) station leave-one-out, top-5 support-row bridge
 contributors after gating (B018, STOR, NLWA, B943, LEBA/B017): all same sign, |dL|<50%, else ARTIFACT.
 (6) (A) N=6 only if candidate/ambiguous. (7) write-up under Q7 with Q8's six confounds + explicit statement that
 sig_agg/sig_frame are not the noise for slow features.
Merlin's predictions for (3): L ~0+-0.8 -> artifact; L ~ -0.9+-0.6 -> ambiguous; L <= -1.5 stable across gates -> candidate.
Merlin scratch scripts: scratchpad/levelchange_check{,2,3,4}.py, b017_check.py (session scratchpad).

## ADDENDUM (03:45, reply to census): D-ONLY absolute gate 0.3% CONFIRMED primary; 365-d RANGE criterion WITHDRAWN (his
calibration error: expected noise range of 15 yearly means ~0.46%; D noise-expected median ~0.06% vs observed 0.226 -> D
discriminates at multi-year scale). Optional statistical variant (low priority): SE_D = sig_f*sqrt(1/N1+1/N2), N_k =
max(1,min(n_k, span_k/75d)), flag |D|/SE_D > 4. CAVEAT for pre-registration: |D| exclusion truncates signal+drift
symmetrically -> gated |L| is a LOWER bound; if |L| grows monotonically with gate value (0.2 -> 0.5), truncation is
eating signal -> no "artifact" call from the 0.2 run alone; |L| ~0 at all three -> artifact. (b) non-bridge family-eras:
gate only spans >= 3 yr via span-halves criterion |dH| > 0.3% (first-half minus second-half mean of dw), reported vs span
as calibration; ONE extra chain only if primary lands candidate/ambiguous. (c) 0.2/0.5 runs D-only, same lambda/sweeps;
per gate write the retained bridge family-eras per station inside q16 support -> step-5 leave-one-out roster from the
GATED network.
