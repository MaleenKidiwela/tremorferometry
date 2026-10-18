"""B4 steps 2-5 of Merlin's ROUND-2 prescription (2026-08-22): deconvolved, normalized patch product.

Consumes transmission_T.npz (v2: single-patch daily T + T_qp at 2009/2010/2016/2023) and
inversion_4d_daily.npz. Factorization (Merlin ruling d): T_qp = C · diag(dia) with C = column-
diagonal-normalized shape matrix (epoch-mean); deconvolved coefficients
    a(t) = diag(T_day(t))^-1 · C^-1 · y(t),   y_q(t) = <pat_q, MODEL(t)> / <pat_q, pat_q>
UNDAMPED 20x20 inverse (ruling a). Steps:
  2  2009 shape acceptance: 2009 column-normalized shape within the 2010-2023 spread -> include
     2009-10, else mask 2009-10 in the deconvolved deep product.
  3  Noise propagation at the 4 rep days: FULL per-day posterior covariance (day_solve want_cov=
     "full", iid-weight posterior per Merlin spec; station-common not included — noted caveat),
     projected onto patch functionals, propagated through C^-1 -> per-patch post-deconvolution
     sigma_frame and aggregated sigma (Neff = record_span/75 d). Mask = aggregated sigma > 0.05%.
  4  Product: deconvolved daily patch series, 4-epoch shape systematic (spread of a(t) recomputed
     with each epoch's own C), re-gauged to zero-mean over R=2019-2024 (item-6b), sigma table.
  5  Aliasing probe: blob midway between two deep patch centers on 2023 rep day; deconvolved
     coefficients must split ~evenly and preserve summed amplitude within ~20%.
Outputs: res_catalog/transmission_deconv.npz + transmission_deconv.png. Nothing overwritten.
NOTE: amplitude claims STILL await the end-to-end BCD ramp injection (original step 4), which must
run THROUGH this deconvolution pipeline.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invert_dvv_4d_daily as inv

RD = inv.OUT
R_WIN = inv.R_WIN
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
Tday = np.asarray(zt["T"], float)                  # (T, P) single-patch daily diagonals
Tqp = np.asarray(zt["Tqp"], float)                 # (4, P, P) epochs 2009/2010/2016/2023
patterns = np.asarray(zt["patterns"], float)
band = zt["patch_band"]; P = len(band)
selfdot = (patterns ** 2).sum(1)
days = pd.to_datetime(zp["days"]); T = len(days)
MODEL = np.asarray(zp["MODEL"], float)
LAMF, LAMS = float(zp["lam_f"]), float(zp["lam_s"])
rep_dates = [str(x) for x in zt["rep_days"]]
rep = [int(days.get_indexer([pd.Timestamp(x)], method="nearest")[0]) for x in rep_dates]
deep = band == "deep"
inR = np.asarray((days >= R_WIN[0]) & (days <= R_WIN[1]))

# ---- shape matrices: column-diagonal-normalized ----
with np.errstate(divide="ignore", invalid="ignore"):
    Cs = np.array([M / np.diag(M)[None, :] for M in Tqp])       # (4,P,P); nan cols = dead epochs
k09 = rep_dates.index([d for d in rep_dates if d.startswith("2009")][0])
kmod = [k for k in range(len(rep_dates)) if k != k09]
spread_mod = Cs[kmod].max(0) - Cs[kmod].min(0)                  # 2010-2023 entrywise spread
dev09 = np.abs(Cs[k09] - Cs[kmod].mean(0))
thr = np.maximum(spread_mod, 0.05)
exceed = np.argwhere(np.nan_to_num(dev09) > thr)                 # nan-safe (Merlin E-i)
dia09 = np.diag(Tqp[k09])
# headline statistic restricted to usable-in-2009 columns (dia > 0.2) — Merlin E-ii
usable = dia09 > 0.2
dev09_restricted = float(np.nanmax(dev09[:, usable])) if usable.any() else np.nan
restricted_exceeds = bool((np.nan_to_num(dev09[:, usable]) > thr[:, usable]).any())
mask09_needed = bool(len(exceed))                                 # conservative: mask either way
mask09_reason = ("shape instability" if restricted_exceeds
                 else "missing illumination (weak/dead 2009 columns only)")
exceed_tab = np.array([(int(r), int(c), float(dev09[r, c]), float(dia09[c]))
                       for r, c in exceed if np.isfinite(dev09[r, c])])
print(f"STEP 2 — 2009 shape: restricted max dev (cols dia09>0.2) {dev09_restricted:.3f}, "
      f"2010-23 max spread {np.nanmax(spread_mod):.3f}; {len(exceed)} exceeding entries -> "
      f"{'MASK 2009-10 deep' if mask09_needed else 'INCLUDE 2009-10'} [{mask09_reason}]")
C_mean = Cs[kmod].mean(0)                                        # primary shape (mature epochs)
Cinv = np.linalg.inv(C_mean)

# ---- step 3: noise propagation at rep days ----
PROJ = patterns / selfdot[:, None]                               # (P, Mf) patch functionals
span_d = (days.max() - days.min()).days
neff = max(span_d / 75.0, 1.0)
sig_frame = np.zeros((len(rep), P))
SDEC = np.zeros((len(rep), P, P))                                # FULL coeff covariance (Merlin r3-A)
for k, t in enumerate(rep):
    _, _, _, Cov = inv.day_solve(t, np.zeros(inv.P), LAMF, LAMS, want_cov="full")
    Spatch = PROJ @ Cov @ PROJ.T                                 # patch-coeff covariance
    dscale = np.clip(Tday[t], 0.05, None)                        # daily diagonal division
    SDEC[k] = (Cinv @ Spatch @ Cinv.T) / np.outer(dscale, dscale)
    sig_frame[k] = np.sqrt(np.clip(np.diag(SDEC[k]), 0, None))
    amp = np.sqrt(np.diag(Cinv @ Spatch @ Cinv.T) / np.clip(np.diag(Spatch), 1e-12, None))
    print(f"STEP 3 {rep_dates[k]}: deconv sigma_frame deep med {np.median(sig_frame[k][deep]):.3f}% "
          f"(max {sig_frame[k][deep].max():.3f}) | noise amplification C^-1 med {np.median(amp):.2f} "
          f"max {amp.max():.2f}")
sig_frame_rep = sig_frame[kmod].mean(0)                          # mature-epoch mean per patch
sig_agg = sig_frame_rep / np.sqrt(neff)
mask = sig_agg > 0.05
print(f"STEP 3 mask (agg sigma > 0.05%): {mask.sum()}/{P} masked "
      f"({(mask & deep).sum()}/{deep.sum()} deep) | Neff {neff:.0f} | "
      f"deep agg sigma med {np.median(sig_agg[deep]):.4f}% [PROVISIONAL-v1, pre-mirror-v2]")

# ---- r3-A: merged adjacent-deep tiles from the coefficient covariance ----
Smat = SDEC[kmod].mean(0)                                        # mature-epoch coeff covariance
sc_all = np.asarray(zt["patch_sc"], float)
dq = np.where(deep)[0][np.argsort(sc_all[np.where(deep)[0]])]
merged = []
for q1, q2 in zip(dq[:-1], dq[1:]):
    var_m = 0.25 * (Smat[q1, q1] + Smat[q2, q2] + 2 * Smat[q1, q2])
    corr = Smat[q1, q2] / np.sqrt(max(Smat[q1, q1] * Smat[q2, q2], 1e-18))
    s_m = np.sqrt(max(var_m, 0)) / np.sqrt(neff)
    merged.append((int(q1), int(q2), float(corr), float(s_m), bool(s_m <= 0.05)))
    print(f"  merged deep q{q1}+q{q2} (lat {zt['patch_lat'][q1]:.1f}/{zt['patch_lat'][q2]:.1f}): "
          f"err-corr {corr:+.2f}, merged sig_agg {s_m:.4f}% -> {'PASS' if s_m <= 0.05 else 'fail'}")
merged = np.array(merged)

# ---- r3-D: precision-weighted deep-margin aggregate (ALL 10 deep tiles, frozen weights) ----
Sdd = Smat[np.ix_(np.where(deep)[0], np.where(deep)[0])]
ones = np.ones(deep.sum())
wi = np.linalg.solve(Sdd, ones); w_deep = wi / wi.sum()          # min-variance weights
sig_agg_deep = float(1.0 / np.sqrt(wi.sum()) / np.sqrt(neff))
print(f"AGGREGATE: precision-weighted deep-margin mean, sigma_agg {sig_agg_deep:.4f}%; "
      f"weights {np.round(w_deep, 2)} (WA-dominated as expected)")

# ---- step 4: deconvolved product ----
Y = np.array([PROJ @ np.nan_to_num(MODEL[t]) for t in range(T)])          # raw patch coeffs
dscale = np.clip(Tday, 0.05, None)
A = (Cinv @ Y.T).T / dscale                                                # primary series
A_ep = np.stack([(np.linalg.inv(Cs[k]) @ Y.T).T / dscale for k in kmod])   # per-epoch variants
sys_shape = A_ep.max(0) - A_ep.min(0)                                      # shape systematic
mu = A[inR].mean(0)                                                        # re-gauge over R
A = A - mu[None, :]
A_ep = A_ep - A_ep[:, inR].mean(1)[:, None, :]
if mask09_needed:
    pre11 = days < pd.Timestamp("2011-01-01")
    A[np.ix_(pre11, deep)] = np.nan

# ---- step 5: aliasing probe (blob midway between two deep centers, 2023 day) ----
sc = np.asarray(zt["patch_sc"], float)
dq = np.where(deep)[0]; o = np.argsort(sc[dq])
q1, q2 = dq[o[len(o) // 2]], dq[o[len(o) // 2 + 1]]
s_mid = 0.5 * (sc[q1] + sc[q2])
a = inv.cxy - inv.cxy.mean(0)
u = np.linalg.eigh(a.T @ a)[1][:, -1]; s = a @ u
dep_cells = np.asarray(zp["depth"], float) >= 30.0
blob = np.zeros(inv.Mf)
blob[dep_cells] = 1.0 * np.exp(-0.5 * ((s[dep_cells] - s_mid) / (125.0 / 2.355)) ** 2)
t23 = rep[rep_dates.index([d for d in rep_dates if d.startswith("2023")][0])]
idx = inv.BY_DAY[t23]
mrec, _, _, _ = inv.day_solve(t23, np.zeros(inv.P), LAMF, LAMS,
                              d_override=inv.Gc[inv.ROW[idx]] @ blob)
arec = (Cinv @ (PROJ @ mrec)) / np.clip(Tday[t23], 0.05, None)
atru = Cinv @ (PROJ @ blob)                                       # truth in the same basis
tot_r, tot_t = arec[[q1, q2]].sum(), atru[[q1, q2]].sum()
print(f"STEP 5 aliasing probe (deep, between q{q1}/q{q2}): recovered split "
      f"{arec[q1]:.2f}/{arec[q2]:.2f} vs truth {atru[q1]:.2f}/{atru[q2]:.2f}; "
      f"summed amp ratio {tot_r / tot_t:.2f} (bar 0.8-1.2)")

# aggregate series (frozen w over ALL deep tiles; NaN pre-2011 via the deep mask09)
di = np.where(deep)[0]
a_agg = A[:, di] @ w_deep
a_agg_area = A[:, di].mean(1)                                    # secondary, naive estimand
np.savez_compressed(f"{RD}/transmission_deconv.npz",
                    A=A.astype(np.float32), days=days.astype(str).values,
                    sys_shape=sys_shape.astype(np.float32), sig_frame=sig_frame_rep,
                    sig_agg=sig_agg, mask=mask, C_mean=C_mean, Cinv=Cinv, Cs=Cs,
                    SDEC=SDEC.astype(np.float32), merged_pairs=merged,
                    w_deep=w_deep, agg=a_agg.astype(np.float32),
                    agg_area=a_agg_area.astype(np.float32), sig_agg_deepmean=sig_agg_deep,
                    dev09_max=dev09_restricted, dev09_exceed=exceed_tab,
                    mask09=mask09_needed, mask09_reason=mask09_reason,
                    provisional="v1-pre-mirror-v2",
                    probe_scope=("amplitudes for structure unaligned with the tiling carry up to "
                                 "~20% ONE-SIDED LOW bias (can hide signal, never fake it)"),
                    patch_band=band, patch_lat=zt["patch_lat"], patch_sc=sc,
                    probe=np.array([arec[q1], arec[q2], atru[q1], atru[q2], tot_r / tot_t]),
                    neff=neff)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
for i, bd in enumerate(("shallow", "deep")):
    qi = np.where(band == bd)[0]
    latv = np.asarray(zt["patch_lat"], float)
    cm = plt.cm.viridis((latv[qi] - latv.min()) / (latv.max() - latv.min() + 1e-9))
    for q, c in zip(qi, cm):
        ser = pd.Series(A[:, q], index=days).rolling("180D", min_periods=30).mean()
        ax[i].plot(ser.index, ser.values, lw=0.9, color=c,
                   ls="--" if mask[q] else "-")
    ax[i].set_ylabel(f"{bd} deconv δβ/β (%)"); ax[i].axhline(0, color="#94a3b8", lw=0.5)
ax[0].set_title("Deconvolved + T-normalized patch series (180-d mean; dashed = noise-masked; "
                "amplitudes await end-to-end BCD calibration)", fontsize=10)
fig.savefig(f"{RD}/transmission_deconv.png", dpi=110, bbox_inches="tight")
print(f"wrote {RD}/transmission_deconv.npz + transmission_deconv.png")
