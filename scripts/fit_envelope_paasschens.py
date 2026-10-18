#!/usr/bin/env python
"""Estimate the S-wave scattering mean free path ell from the stacked LFE coda envelope, per station and fleet-wide,
by fitting the Paasschens (1997) 3-D isotropic radiative-transfer energy solution (direct delta + scattered term)
plus intrinsic absorption and a noise floor to the RMS-by-0.5-s-bin profiles in data/window_snr_profile.csv
(grand stacks, S at stack +1.0 s, RMS relative to the 8-10 s floor).

Model: E(t) = A * exp(-beta*t/ell) * exp(-2*pi*f*t/Qi) * [ direct: exp(0) delta(t - r/beta)/(4 pi r^2 beta) -> bin
       energy; scattered: (1 - r^2/(beta t)^2)^(1/8) / (4 pi ell beta t / 3)^(3/2) * G((beta t/ell)(1 - r^2/(beta t)^2)^(3/4)) ]
       + N,  G(x) = exp(x) sqrt(1 + 2.026/x);  fitted in log space over the S bin (+1.0..1.5 s stack) and 2.5..9.5 s
       (the 1.5-2.5 s S-tail is excluded: forward-scattered/pulse energy that isotropic RTE does not model).
r = hypocentral distance (per-station median over its pairs, from G.npz offset + cell depth). Grid search over
ell (2-300 km), Qi (30-2000 at f = 5 Hz), N. Output: data/envelope_ell_fit.csv + printed fleet summary.
"""
import numpy as np, pandas as pd

BETA, F0 = 3.5, 5.0
prof = pd.read_csv("data/window_snr_profile.csv")
gz = np.load("fault_tomography/inversion/res_catalog/G.npz", allow_pickle=True)
dep = dict(zip(gz["cell"], gz["depth_km"]))
pr = pd.DataFrame(dict(tag=gz["pair_tag"], off=gz["offset"], dep=[dep[c] for c in gz["pair_cell"]]))
pr["r"] = np.sqrt(pr.off ** 2 + pr.dep ** 2)
rmed = pr.groupby("tag").r.median()

ELL = np.logspace(np.log10(2), np.log10(300), 60)
QI = np.array([30, 50, 100, 200, 400, 800, 2000, 1e9])
NREL = np.logspace(-8, -1, 36)                                   # noise energy relative to direct-bin energy


def paasschens(r, t, ell):
    """scattered energy density (isotropic RTE, Paasschens 1997 eq. 36 without absorption), t > r/beta"""
    a = 1 - (r / (BETA * t)) ** 2
    a = np.clip(a, 1e-12, None)
    x = (BETA * t / ell) * a ** 0.75
    G = np.exp(x) * np.sqrt(1 + 2.026 / np.maximum(x, 1e-9))
    return a ** 0.125 / (4 * np.pi * ell * BETA * t / 3) ** 1.5 * np.exp(-BETA * t / ell) * G


def model_bins(r, ell, qi, bins):
    """energy per 0.5-s stack bin; bins = stack-time bin starts; S at stack 1.0 => lapse t = r/beta + (stack - 1)"""
    tS = r / BETA
    out = np.zeros(len(bins))
    for k, b0 in enumerate(bins):
        tt = tS + (np.linspace(b0, b0 + 0.5, 6)[:-1] + 0.05) - 1.0
        tt = tt[tt > tS + 1e-3]
        e = paasschens(r, tt, ell) * np.exp(-2 * np.pi * F0 * tt / qi) if len(tt) else np.zeros(0)
        out[k] = e.mean() * 0.5 if len(tt) else 0.0
        if b0 <= 1.0 < b0 + 0.5:                                   # direct arrival falls in this bin
            out[k] += np.exp(-r / ell) * np.exp(-2 * np.pi * F0 * tS / qi) / (4 * np.pi * r ** 2 * BETA)
    return out


def fit(profile_rel, r):
    """profile_rel: Series bin_t0 -> RMS rel to floor. Returns best (ell, qi, nrel, rms_logres)."""
    bins = np.array(sorted(profile_rel.index))
    use = (bins >= 1.0) & (bins < 1.5) | (bins >= 2.5) & (bins < 10.0)
    eobs = profile_rel.values ** 2
    floor = eobs[(bins >= 8.0) & (bins < 10.0)].mean()
    best = None
    for ell in ELL:
        for qi in QI:
            m = model_bins(r, ell, qi, bins)
            if not np.all(np.isfinite(m)) or m.max() <= 0:
                continue
            md = m / m[(bins >= 1.0) & (bins < 1.5)][0]           # relative to direct bin
            for nrel in NREL:
                mm = md + nrel
                mm = mm / mm[(bins >= 8.0) & (bins < 10.0)].mean() * floor
                res = np.log(mm[use]) - np.log(eobs[use])
                s = float(np.sqrt(np.mean(res ** 2)))
                if best is None or s < best[3]:
                    best = (ell, qi, nrel, s)
    return best


rows = []
for tag, g in prof.groupby("tag"):
    if tag not in rmed.index:
        continue
    s = g.set_index("bin_t0").rms_rel
    ell, qi, nrel, rms = fit(s, float(rmed[tag]))
    rows.append(dict(tag=tag, r_km=float(rmed[tag]), ell_km=ell, Qi=qi, noise_rel=nrel, rms_logres=rms,
                     snr35=float(s.get(3.0, np.nan))))
df = pd.DataFrame(rows).sort_values("ell_km")
df.to_csv("data/envelope_ell_fit.csv", index=False)
fl = prof.groupby("bin_t0").rms_rel.median()
ell, qi, nrel, rms = fit(fl, float(pr.r.median()))
print(f"FLEET-median profile (r = {pr.r.median():.0f} km): ell = {ell:.1f} km, Qi = {qi:g}, noise/direct = {nrel:.1e}, "
      f"rms log-residual {rms:.2f}")
print(f"  mean free time {ell / BETA:.1f} s; lapse/mft at 3-5 s window (t = r/beta + 3 s): "
      f"{(pr.r.median() / BETA + 3) / (ell / BETA):.1f}")
print(f"PER-STATION ({len(df)}): ell median {df.ell_km.median():.1f} km, IQR {df.ell_km.quantile(.25):.1f}-"
      f"{df.ell_km.quantile(.75):.1f}, Qi median {df.Qi.median():g}, fit rms median {df.rms_logres.median():.2f}")
print("  ell distribution (km):", np.percentile(df.ell_km, [10, 25, 50, 75, 90]).round(1))
print("wrote data/envelope_ell_fit.csv")
