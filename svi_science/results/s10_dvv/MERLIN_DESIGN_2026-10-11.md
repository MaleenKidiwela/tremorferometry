# S10 dv/v from LFE repeat stacks — Merlin design (2026-10-11), adopted for the pilot

Verdict: viable only as per-big-episode, within-response-epoch, pairwise dv/v on single-epoch stations: PB boreholes (B011 2005-11→, B926 2007-11→, B001–B010) plus PGC BHZ 2000-06→2015-01-12. Daily dv/v is not viable. Between episodes there are almost no real repeats (282 candidates at P ≥ 0.75 on zero-tremor days in 12 yr), and below P 0.75 the quiet-day candidates are 36–71 % scrambled control.

Rules:
- Never reference across a response-epoch boundary. Each (station, epoch) is its own series. Cross-era continuity comes only from never-changed boreholes. The all-time reference in scripts/dvv_roll30cal.py l.58 leaked sensor steps before.
- P-weighting: clean on tremor days down to P ≥ 0.5 (≤ 5 % control weight). Not clean on quiet days (≥ 36 %). Never use "all candidates" (25–43 % even on tremor days).

Pilot (go/no-go):
1. Inventory N(P ≥ 0.75 counted-once; era P ≥ 0.90) and ΣP(P ≥ 0.5) per (family, station, big episode). Stations: B011, B009, B010, B926, PGC, VGZ, SNB, NLLB, each with its response epoch. Families: 10 + one of 43/56/87. Accept ≥ 8 episodes with N ≥ 100 at ≥ 2 single-epoch stations; otherwise S10 is a bound.
2. Windows:
   - 40 Hz, 2–8 Hz, despike MAD 8, per-window RMS normalisation.
   - cc-realign on S−1..S+2 s (±0.5 s) to the family-station grand stack. Re-pick S once per (family, station, epoch).
   - Build identical windows for the scrambled-control rows.
   - Sources: lfe_extension/results/cut/svi_ext/ for 2010-08→2026; lfe_extension/traces_backward for 2003–2010-07.
3. N_eff: disjoint random halves n = 25/50/100/200 in the two biggest episodes. N_eff = smallest n with half-split stretch SD ≤ 0.05 % (2–4 s and 3–5 s). Accept N_eff ≤ 150 on Z at ≥ 2 stations.
4. Unit: P-weighted stack per (family, station, big episode), headline weights P ≥ 0.5. Also compute P ≥ 0.75 equal-weight and P ≥ 0.3 weighted. Causal 5–7 d sub-stacks only where ΣP_excess ≥ N_eff. No 30-d boxcar, no SVD-Wiener, no all-time reference.
5. Estimator:
   - Pairwise doublets between episode stacks within an epoch, as in dvv_pairwise.py: S-anchored stretch, Huber-weighted graph, zero mean per epoch.
   - Windows 2–4 s and 3–5 s. Mirror window S−15..S−(S−P+1) s, i.e. before P.
   - Run the scrambled-control stack through the identical path.
   - eps = −dv/v.
6. Gates:
   - (a) coda/mirror RMS > 1.5.
   - (b) control-stack stretch within ±2× floor; |mirror| < |coda| in ≥ 70 % of pairs.
   - (c) P ≥ 0.75-only vs P-weighted: |Δ| < floor.
   - (d) matched-N subsampling, 50×.
   - (e) lapse slope 3–5 s on 2–4 s ≥ 0.6; sign agreement ≥ 70 %.
   - (f) differential step injection of ±0.1 % and ±0.3 %: recovery 0.8–1.2; leakage < 20 %.
7. GO if 1, 3 and (a)–(f) pass for both pilot families at ≥ 2 stations. Then scale up, and run WA-NW on B001–B007.

Reach: same-source 2005–2026 history on boreholes (~16 episodes); PGC 2003–2015 (9 episodes); front vs tail within an episode; healing as episode-to-episode differences; SVI vs WA-NW on identical instruments; whether the ETS-locked systematic vanishes with a P-clean population. Expected precision ~0.1–0.2 % per family-station-episode, and ~0.03–0.05 % per episode with ≥ 10 family-stations.
