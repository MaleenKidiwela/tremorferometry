#!/usr/bin/env python
"""Analyse the fleet-wide direct stack lapse test (stack_lapse_fleet.csv): classify stations, null for the classification
counts, map, relation to noise-field proxy (pre-arrival RMS ratio), coda/direct change, modeled station strain, latitude."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
SD = "fault_tomography/inversion/res_catalog_sse"
import glob
files = [f for f in sorted(glob.glob(f"{SD}/stack_lapse_fleet_s*.csv")) if "families" not in f]
T = pd.concat([pd.read_csv(f) for f in files]) if files else pd.read_csv(f"{SD}/stack_lapse_fleet.csv"); T = T[T.n_fam >= 10].drop_duplicates("tag").copy()
T["sig24"] = T.e24.abs() > 3 * T.e24_se
T["lapse_consistent"] = (np.sign(T.e24) == np.sign(T.e35)) & (T.ratio35_24.between(0.5, 1.5)) & (T.sign_agree >= 0.7)
T["flip"] = (np.sign(T.e24) != np.sign(T.e35)) & T.sig24
S = pd.read_csv(f"{SD}/station_model_strain_ONPRE.csv"); S = S[S.depth_km == 1.0].set_index("tag")
T["dil_model"] = T.tag.map(S.dil_change_1e8); T["ezz_model"] = T.tag.map(S.ezz_change_1e8)
print(f"stations {len(T)} | |e24| > 3 SE: {T.sig24.sum()} | of those lapse-consistent (same sign, ratio 0.5-1.5, >=70% families agree): {(T.sig24 & T.lapse_consistent).sum()}; sign-flip: {T.flip.sum()}")
sig = T[T.sig24]
print(f"significant stations: e24 mean {sig.e24.mean():+.3f}, positive {100*(sig.e24>0).mean():.0f}%; median |e24| {sig.e24.abs().median():.3f}%")
print(f"time-half replication among significant (early vs late e24 same sign): {100*np.mean(np.sign(sig.e24_early)==np.sign(sig.e24_late)):.0f}% of {int(np.isfinite(sig.e24_early*sig.e24_late).sum())} with both")
print(f"Spearman(e24, lat) {spearmanr(T.e24, T.lat).correlation:+.2f} | (e24, noise ratio ON/PRE) {spearmanr(T.e24, T.noise_ratio).correlation:+.2f} | (|e24|, |1-noise ratio|) {spearmanr(T.e24.abs(), (T.noise_ratio-1).abs()).correlation:+.2f} | (e24, coda/direct change) {spearmanr(T.e24, T.cdr).correlation:+.2f} | (e24, modeled dilatation) {spearmanr(T.e24, T.dil_model, nan_policy='omit').correlation:+.2f} | (e24, modeled ezz) {spearmanr(T.e24, T.ezz_model, nan_policy='omit').correlation:+.2f}")
lc = T[T.sig24 & T.lapse_consistent]
print(f"\nLAPSE-CONSISTENT significant stations ({len(lc)}): borehole {lc.borehole.sum()}, broadband {(~lc.borehole).sum()}; positive {(lc.e24>0).sum()}, negative {(lc.e24<0).sum()}")
print(lc.sort_values("e24")[["tag", "lat", "borehole", "n_fam", "n_ep", "e24", "e24_se", "e35", "ratio35_24", "sign_agree", "noise_ratio", "cdr", "e24_early", "e24_late", "dil_model"]].round(3).to_string(index=False))
fl = T[T.flip]
print(f"\nSIGN-FLIP stations ({len(fl)}):"); print(fl.sort_values("e24")[["tag", "lat", "borehole", "n_fam", "e24", "e35", "noise_ratio", "cdr"]].round(3).to_string(index=False))
T.to_csv(f"{SD}/stack_lapse_fleet_classified.csv", index=False)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(17, 6))
a = ax[0]; a.scatter(T.lon, T.lat, c="#e5e7eb", s=10)
sc = a.scatter(sig.lon, sig.lat, c=sig.e24, cmap="RdBu_r", vmin=-0.4, vmax=0.4, s=[60 if x else 25 for x in sig.lapse_consistent], edgecolor=["k" if x else "none" for x in sig.lapse_consistent], linewidth=.6)
a.set_aspect(1 / np.cos(np.radians(47))); plt.colorbar(sc, ax=a, shrink=.7, label="ETS-locked stack stretch, 2-4 s (%)"); a.set_title("significant stations (|e24|>3 SE); black edge = lapse-consistent (2-4 & 3-5 s agree)", fontsize=8)
a = ax[1]; a.scatter(T.e24, T.e35, c=["#dc2626" if b else "#1d4ed8" for b in T.borehole], s=12); a.plot([-1, 1], [-1, 1], "k--", lw=.5); a.axhline(0, color="k", lw=.4); a.axvline(0, color="k", lw=.4)
a.set_xlim(-0.7, 0.7); a.set_ylim(-0.7, 0.7); a.set_xlabel("2-4 s stack stretch (%)"); a.set_ylabel("3-5 s stack stretch (%)"); a.set_title("lapse test per station (red = borehole); velocity change lies near the diagonal", fontsize=8)
a = ax[2]; a.scatter(T.noise_ratio, T.e24, c=["#dc2626" if b else "#1d4ed8" for b in T.borehole], s=12); a.axhline(0, color="k", lw=.4); a.axvline(1, color="k", lw=.4); a.set_xlabel("pre-arrival RMS ratio ON/PRE (noise-field proxy)"); a.set_ylabel("2-4 s stack stretch (%)"); a.set_title("response vs noise-field change", fontsize=8)
fig.suptitle("Fleet-wide direct stack test of the ETS-locked change (ON [+15,+60] vs PRE [-90,-30] mean stacks, admitted slip episodes within 150 km)", fontsize=10)
fig.savefig("figures/sse_stack_lapse_fleet.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_stack_lapse_fleet.png, stack_lapse_fleet_classified.csv")
