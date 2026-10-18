# MERLIN — WA origin sign-off review (2026-10-14 ~8:45 pm PT; saved by main session, condensed verbatim points)
Verdict: SVI reportable tonight as provisional. WA is real, but "origin = time − 5.5 s" is the wrong model: the WA catalogue
`time` carries a family-specific constant (SD 3.3 s, range ≈ −7…+4 s relative to fam38); fam75 (6 % of WA rows) is mislocated
by ~60 km with inverted template polarity (also fam121, fam77 suspect); ~30 % of rows remain cross-family duplicates after
correction. Do not promote the c6 column; derive per-family offsets from the catalogue alone, re-dedupe, re-stack, W1 on the
original rule.
Evidence: cross-family neighbour pairs at fixed offsets (IQR ≤ 0.1 s): fam38–fam75 879 pairs −5.65 s; fam26–fam117 434 +9.2 s;
fam58–fam77 397 −4.2 s. Pair network: 73 families, 21,924 rows, rms 1.5 s. Smearing attenuation 0.78 at 30 s. Catalogue-only
prediction matches observed VLF lags per family (fam37 −1.9 vs −2…+1; fam58 −11.1 vs −10.5; fam106 −7.3 vs −7/−8; fam126 −4.5 vs
−5). Absolute anchor from 321 clear template S picks: VLF lag − δ = +1…+4 s for strong families (SVI mode +2 s).
(2) WA ≈ 0.5× SVI not interpretable: "same order, lower bound 2.7e12"; expected 0.7–1.0× after correction.
(3) "Scales with LFE strength" fair as provisional: composition 0.97–1.03; control rows ≤ 1.4 % per P bin; P tracks net_mad
(Spearman 0.76); hi−lo ≈ 2.8σ over 4 episodes. Wording: "LFEs detected more strongly carry about twice the 20–50 s moment".
Clean test: net_mad terciles within P ≥ 0.9.
Prescription: (1) Amendment 6, per-family offsets (pair medians ≥ 8 co-detections ±12 s, WLS √n, anchor to clear-pick δ on
well-located families; acceptance rms < 2 s, network–arrival SD < 3 s; failure SD > 4 s → stop). (2) a-priori location-suspect
flag (≥ 30 % co-detection with a family > 15 km away; depth < 20 km) → duplicate cluster takes best-located member's location.
(3) rebuild rows: to_f = time + 0.5 + offset_f; dedupe ±2 s on corrected times (n_nb2 = 0, n_nb10 < 0.15, 25–35 % removed).
(4) re-stack 54 + between, s13_stack2_v2.py, kmax 5, ±8 s grid; W1 (lag mode +1…+4, ≤ 10 % at ±8), zones 2–3 positive in > half
of ≥ 300-LFE episodes, pooled ≥ 3.5e12, W3 0.85–1.15, A1–A6; fixed-lag +2 s column. (5) SVI net_mad terciles within P ≥ 0.9 for
2004-07, 2005-09, 2012-08, 2013-09 (top/bottom ≥ 1.5 in ≥ 3/4 → scales with detection amplitude).
Change mind: corrected W1 still failing; fam38–fam75 offset drifting > 1.5 s 2008→2024 (→ per family-epoch); terciles ratio ≈ 1.
