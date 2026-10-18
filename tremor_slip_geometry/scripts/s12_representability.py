#!/usr/bin/env python
"""
S12 — can the 13-pattern vbICA slip model represent slip under the central-Oregon
tremor cluster?  (PLAN.md Round 3, S12)

WHY
---
S7 found (a) a persistent cluster of nine triangles in central Oregon, 44.33-45.42N,
Slab2 depth 27-38 km, that are TREMOR-ONLY (high tremor-day rank, low modeled-slip rank) in
all three epochs and are nowhere near a mesh boundary, and (b) that the rank co-location of
modeled slip with tremor occupancy is far weaker there (Spearman 0.40-0.46) than in the north
(0.71-0.77) or the south (0.70-0.73).  Three explanations compete:

  (1) MODEL-BASIS limitation.  Gualandi (2025) describes the entire 2007-2026 slip history as
      13 FIXED spatial patterns times 13 scalar time functions.  Slip therefore lives in a
      13-dimensional subspace of the 742-dimensional triangle space.  If a compact slip patch
      at 44.3-45.4N is close to orthogonal to that subspace, the model could NEVER have placed
      slip there no matter what the GNSS data said, and the "decoupling" is an artefact of the
      basis.
  (2) COVERAGE limitation.  If the Oregon Coast Range is poorly instrumented, or if the
      stations that do exist fit badly, the inversion has no information there.
  (3) A PHYSICAL result: a tremor-productive, slip-poor (or aseismic-slip-poor) stretch of the
      interface.

This script separates (1) and (2) from (3).  It does NOT test (3) directly; it can only
establish whether (1) and (2) are large enough to explain what S7 saw.

--------------------------------------------------------------------------------------
INPUTS (all read-only, all outside this folder)
--------------------------------------------------------------------------------------
  /home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt     9-column triangular mesh
                                                              (lon lat h) x 3, 742 patches
  /home/jovyan/cascadia-sse-strain/data/2026-08-18/           one dated Gualandi solution:
      X.mat            Xd: station llh (639 rows = 213 stations x [e,n,u]), name, timeline,
                       unit ('mm'), ts, var_ts
      ICA.mat          ICA_essential: U (639 x 26 spatial responses), var_U (639 x 26 data
                       variances), S (26), V (7169 x 26 time functions)
      ind_comps.mat    the 13 retained component indices (1-based: 5,10,11,15,16,17,20,21,
                       22,23,24,25,26)
      misfit_comps.mat (9 sigma) x (13 component) cross-validation misfit grid; its per-column
                       argmin recovers the per-component smoothing hyperparameter that the
                       public download otherwise omits
      options.mat      inversion settings: lambda_strike = lambda_dip = lambda0 = 21.09 km,
                       exponential prior, sigma grid 1e-4..1e-1, m0 = 0, fault.nu
      fitresult.mat    trajectory-model fit (NOT used here)
  tremor_slip_geometry/results/s7_occupancy_triangles.csv     per-triangle S7 table: tri (0-based,
                                                              same order as the mesh), lat, lon,
                                                              depth_km, area_km2, seg,
                                                              class_persistent, S_tot_mm,
                                                              D3nb_tot_days, N_det_tot,
                                                              on_mesh_edge
  fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz   patch_depth_km, used only to
                                                              confirm depth = -centroid_llh[:,2]

Green's functions are rebuilt with cutde (Nikkhoo & Walter 2015) through
sse_strain.Solution.displacement_gf, and the per-component linear inversion is redone with
Solution.invert (Tarantola-Valette, exponential prior, per-component sigma0 from the
misfit_comps argmin).  Thrust-positive slip patterns are P = -Ld with Ld the dip-slip block of
L converted to metres (data_unit_to_m = 1e-3 because X['unit'] == 'mm').

--------------------------------------------------------------------------------------
DEFINITIONS
--------------------------------------------------------------------------------------
Segments: north lat >= 46, central 44 <= lat < 46, south lat < 44 (project convention).
Distances between triangles are 3-D distances between centroids in the solution's own local
frame (fault.centroid_xyz_km), i.e. ON-FAULT distances -- the same metric the inversion's prior
covariance uses.

A1 REPRESENTABILITY.  For each triangle j a unit thrust "test patch"
        u_j[i] = exp(-0.5 * d_ij^2 / sigma^2),   sigma = 30 km and 50 km,
peak 1.0 at j, thrust-only (pure dip direction, thrust sense).  P = -Ld (742 x 13, rank 13).
Least-squares projection onto the pattern span, u_hat = P (P^T P)^-1 P^T u_j (computed from the
thin QR of P):
        r_j        = 1 - ||u_j - u_hat||^2 / ||u_j||^2       (fraction of the patch's L2 energy
                     the 13 patterns can express; a random 742-vector scores 13/742 = 0.018)
        pot_frac_j = sum_i area_i u_hat_i / sum_i area_i u_j_i   (fraction of the patch's POTENCY
                     the basis reproduces -- the quantity that matters for "would this patch have
                     shown up as slip at all")
        leak_j     = fraction of sum_i area_i |u_hat_i| lying farther than 2 sigma from j
                     (how much of what the basis does reproduce is smeared elsewhere)
A low r_j / pot_frac_j means the fixed basis could not have placed a patch of that size there.
CAVEAT: this is generous to the model.  In the real model a component's amplitude drives its
strike-slip pattern too, so expressing a thrust-only patch is at least as hard as this says.

A2 PATTERN ENERGY.  Per component k, w_k[i] = area_i |P[i,k]|; the central-Oregon share is
sum over 44 <= lat < 46 of w_k divided by its total.  The aggregate "slip capacity" weights each
component by the standard deviation of its amplitude time function (S[k] * V[:,k]), so a
component that never fires contributes little.  Compared against each segment's share of the
mesh AREA (the null share) and against its share of tremor-days D3nb and of raw detections.

A3 GNSS COVERAGE.  Station coordinates ARE recoverable (X['llh'], 639 rows decimated by 3 ->
213 stations; X['name'] gives the 4-character names).  Reported: nearest-station horizontal
distance from every triangle centroid; number of stations within 50 and 100 km; and, per 1-degree
latitude band, the station count and the count within 150 km of the deep ETS band (triangles with
Slab2-equivalent mesh depth 30-45 km).

A4 MISFIT.  misfit_comps.mat is a cross-validation grid, not a per-station misfit, so the
per-station misfit is computed directly from the reproduced inversion:
        misfit_s = sqrt( sum_k sum_{c in e,n,u} (G L[:,k] - U[:,idx_k])_sc^2
                         / sum_k sum_c U[:,idx_k]_sc^2 )
i.e. the fraction of each station's 3-component ICA spatial response, summed over the 13 retained
components, that the slip model fails to explain.  A chi-form normalised by var_U is also given.

A5 SYNTHETIC RECOVERY.  Runs the actual inversion machinery on the test patches, so it measures
what the GNSS NETWORK + PRIOR can do WITHOUT the 13-pattern restriction.  For each component k
(its own sigma0 and Cd = diag(var_U[:,idx_k])) build Cm0 (exponential, lambda 21.09 km),
GC = G Cm0, A = GC G^T + Cd, and recover all 742 test patches at once:
        M_rec = GC^T A^-1 (G M_true),     M_true = [0 ; -U_sigma]  (strike 0, thrust = u_j)
Reported per triangle, median over the 13 components:
        res_pot   = recovered potency / true potency
        res_peak  = recovered thrust at j / 1.0  (peak retention = how much the smoothing blurs)
        det_rel   = 2 sqrt(Cm_post[n+j, n+j]) / res_peak_j, the 2-sigma peak slip a patch at j
                    must have to be detected, normalised per component by its median over the
                    north SLIP+TREMOR triangles (unit-free relative detectability; > 1 = harder
                    to see than a typical northern ETS patch)
The ICA-level noise enters through Cd = var_U, which is the data covariance the published
inversion itself uses; because the estimator is linear the signal recovery is deterministic and
the noise contribution is fully described by the posterior covariance, so no Monte-Carlo is
needed.

A6 CLASS ROBUSTNESS.  If low representability were the cause of S7's TREMOR-ONLY class, dividing
each triangle's accumulated slip by its potency fraction (i.e. undoing the basis's inability to
put slip there) should move those triangles out of the class.  The within-segment percentile rank
of S_tot_mm / pot_frac30 is recomputed and compared with the raw rank.

--------------------------------------------------------------------------------------
OUTPUTS (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------
  results/s12_representability_triangles.csv   742 rows, one per triangle
  results/s12_summary.csv                      medians per segment and per S7 persistent class
  results/s12_pattern_energy.csv               per component + aggregate rows
  results/s12_stations.csv                     213 rows: name, lon, lat, misfit, distances
  results/s12_misfit_by_latitude.csv           1-degree bands: station counts, coverage, misfit
  figures/s12_representability.png             9 panels

Run:  OMP_NUM_THREADS=4 python scripts/s12_representability.py     (~1 min, deterministic)
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution          # noqa: E402
from sse_strain.geodesy import llh2localxy                      # noqa: E402

import matplotlib                                               # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                 # noqa: E402
from matplotlib.lines import Line2D                             # noqa: E402

HERE = Path(__file__).resolve().parent.parent          # tremor_slip_geometry/
ROOT = HERE.parent                                     # tremorferometry/
MESH = "/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt"
SOLDIR = "/home/jovyan/cascadia-sse-strain/data/2026-08-18"
S7CSV = HERE / "results" / "s7_occupancy_triangles.csv"
NPZ = ROOT / "fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.npz"

SIGMAS_KM = (30.0, 50.0)
ETS_DEPTH = (25.0, 45.0)          # the band both tremor and modeled ETS slip occupy
DEEP_BAND = (30.0, 45.0)          # "deep ETS band" for the station-coverage question
COVER_KM = 150.0

np.set_printoptions(suppress=True)


def seg_of(lat):
    return np.where(lat >= 46.0, "north", np.where(lat >= 44.0, "central", "south"))


def rank_pct(x):
    """Percentile rank (average ties), 0-100, matching S7's convention."""
    s = pd.Series(x)
    return (s.rank(method="average") / len(s) * 100.0).values


# ======================================================================== load
print("[1] loading mesh + solution")
fault = load_fault_ascii(MESH)
sol = load_solution(SOLDIR, fault)
NTRI = fault.n_patch
area = fault.area
lon = fault.centroid_llh[:, 0]
lat = fault.centroid_llh[:, 1]
depth = -fault.centroid_llh[:, 2]          # height is negative below sea level

zz = np.load(NPZ, allow_pickle=True)
assert np.allclose(depth, zz["patch_depth_km"], atol=1e-6), "depth sign/ordering mismatch"
assert np.allclose(lat, zz["patch_lat"], atol=1e-9)
print(f"    mesh {NTRI} triangles; depth = -centroid_llh[:,2] confirmed against the npz")

G = sol.displacement_gf()
L = sol.invert(G=G, verbose=False)
idx = sol.component_index
K = idx.size
frac_misfit_comp = sol.misfit_check(G=G)
print(f"    G {G.shape}, L {L.shape}, K={K}, sigma0 source: {sol.sigma0_source}")
print(f"    per-component fractional misfit: "
      f"{np.round(frac_misfit_comp, 3)}")

Ls, Ld = sol.component_slip_patterns()      # metres per unit amplitude
P = -Ld                                     # thrust-positive, 742 x 13
print(f"    thrust pattern matrix P {P.shape}, rank {np.linalg.matrix_rank(P)}")

s7 = pd.read_csv(S7CSV)
assert len(s7) == NTRI and np.allclose(s7["lat"].values, lat)
seg = seg_of(lat)
assert (seg == s7["seg"].values).all()

xyz = fault.centroid_xyz_km
DTRI = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))   # 742 x 742, km

out = pd.DataFrame({
    "tri": np.arange(NTRI),
    "lat": lat, "lon": lon, "depth_km": depth, "area_km2": area,
    "seg": seg,
    "class_s7": s7["class_persistent"].values,
    "D3nb_tot_days": s7["D3nb_tot_days"].values,
    "N_det_tot": s7["N_det_tot"].values,
    "S_tot_mm": s7["S_tot_mm"].values,
    "on_mesh_edge": s7["on_mesh_edge"].values,
})

# the two reference sets
TO_C = ((out["class_s7"] == "TREMOR-ONLY") & (out["seg"] == "central")).values
ST_N = ((out["class_s7"] == "SLIP+TREMOR") & (out["seg"] == "north")).values
out["s7_cluster"] = np.where(TO_C, "central-TREMOR-ONLY",
                             np.where(ST_N, "north-SLIP+TREMOR", ""))
print(f"    S7 reference sets: central TREMOR-ONLY {TO_C.sum()} tri "
      f"(tri {list(np.where(TO_C)[0])}), north SLIP+TREMOR {ST_N.sum()} tri")

# =============================================================== A1 representability
print("[2] A1 representability of unit thrust patches in the 13-pattern span")
Q, _ = np.linalg.qr(P)                       # 742 x 13 orthonormal basis of span(P)
PATCH = {}
for sg in SIGMAS_KM:
    U_ = np.exp(-0.5 * (DTRI / sg) ** 2)     # column j = test patch centred on j
    C_ = Q.T @ U_
    Uh = Q @ C_
    r = (C_ ** 2).sum(0) / (U_ ** 2).sum(0)
    pot = (area @ Uh) / (area @ U_)
    far = DTRI > 2.0 * sg
    aUh = np.abs(Uh) * area[:, None]
    leak = (aUh * far).sum(0) / aUh.sum(0)
    tag = int(sg)
    out[f"r{tag}"] = r
    out[f"pot_frac{tag}"] = pot
    out[f"leak{tag}"] = leak
    PATCH[sg] = U_
    print(f"    sigma {sg:.0f} km: r  north {np.median(r[seg=='north']):.3f}  "
          f"central {np.median(r[seg=='central']):.3f}  south {np.median(r[seg=='south']):.3f}"
          f"   (random-vector baseline {K/NTRI:.4f})")
    print(f"                 pot_frac  north {np.median(pot[seg=='north']):.3f}  "
          f"central {np.median(pot[seg=='central']):.3f}  south {np.median(pot[seg=='south']):.3f}")
    print(f"      central TREMOR-ONLY  r {np.median(r[TO_C]):.3f}  pot {np.median(pot[TO_C]):.3f}"
          f"   |  north SLIP+TREMOR  r {np.median(r[ST_N]):.3f}  pot {np.median(pot[ST_N]):.3f}")

# ================================================================ A2 pattern energy
print("[3] A2 where the model's slip capacity is")
Amp = sol.component_amplitudes()             # 7169 x 13
amp_sd = Amp.std(0)
W = area[:, None] * np.abs(P)                # 742 x 13, area-weighted |thrust|
cap = W * amp_sd[None, :]                    # amplitude-weighted "typical potency" capacity
masks = {s: (seg == s) for s in ("north", "central", "south")}
band = (depth >= ETS_DEPTH[0]) & (depth <= ETS_DEPTH[1])

rows = []
for k in range(K):
    row = {"component": int(idx[k] + 1), "amp_sd": amp_sd[k],
           "sigma0_used": np.nan, "frac_misfit": frac_misfit_comp[k],
           "peak_thrust_m_per_unit_amp": np.abs(P[:, k]).max(),
           "lat_of_peak": lat[np.argmax(np.abs(P[:, k]))]}
    for s, m in masks.items():
        row[f"energy_frac_{s}"] = W[m, k].sum() / W[:, k].sum()
    rows.append(row)
pe = pd.DataFrame(rows)

agg = {"component": "AGGREGATE_amp_weighted", "amp_sd": np.nan, "sigma0_used": np.nan,
       "frac_misfit": np.nan, "peak_thrust_m_per_unit_amp": np.nan, "lat_of_peak": np.nan}
for s, m in masks.items():
    agg[f"energy_frac_{s}"] = cap[m].sum() / cap.sum()
pe = pd.concat([pe, pd.DataFrame([agg])], ignore_index=True)

ref_rows = []
for name, val in (("MESH_AREA", area),
                  ("TREMOR_DAYS_D3nb", out["D3nb_tot_days"].values.astype(float)),
                  ("TREMOR_DETECTIONS", out["N_det_tot"].values.astype(float))):
    rr = {"component": name, "amp_sd": np.nan, "sigma0_used": np.nan, "frac_misfit": np.nan,
          "peak_thrust_m_per_unit_amp": np.nan, "lat_of_peak": np.nan}
    for s, m in masks.items():
        rr[f"energy_frac_{s}"] = val[m].sum() / val.sum()
    ref_rows.append(rr)
pe = pd.concat([pe, pd.DataFrame(ref_rows)], ignore_index=True)
print(pe[["component", "energy_frac_north", "energy_frac_central",
          "energy_frac_south"]].round(4).to_string(index=False))

cap_density = cap.sum(1) / area              # per-km^2 capacity, for the map
out["capacity_density"] = cap_density
out["capacity_density_rel"] = cap_density / np.median(cap_density[ST_N])

# =================================================================== A3 GNSS coverage
print("[4] A3 GNSS coverage")
st_llh = sol.station_llh
names = [str(x)[:-1] for x in np.asarray(sol.X["name"], dtype=object).ravel()[::3]]
assert len(names) == st_llh.shape[0]
sx, sy = llh2localxy(st_llh[:, 1], st_llh[:, 0], fault.origin)
sx, sy = np.ravel(sx), np.ravel(sy)
DH = np.hypot(xyz[:, None, 0] - sx[None, :], xyz[:, None, 1] - sy[None, :])   # 742 x 213

out["nearest_station_km"] = DH.min(1)
out["nearest_station"] = [names[i] for i in DH.argmin(1)]
out["n_stn_50km"] = (DH < 50).sum(1)
out["n_stn_100km"] = (DH < 100).sum(1)

deep = (depth >= DEEP_BAND[0]) & (depth <= DEEP_BAND[1])
near_deep = DH[deep, :].min(0) < COVER_KM

# ========================================================================== A4 misfit
print("[5] A4 per-station misfit")
res2 = np.zeros(G.shape[0])
dat2 = np.zeros(G.shape[0])
chi2 = np.zeros(G.shape[0])
for k in range(K):
    d = sol.U[:, idx[k]]
    pred = G @ L[:, k]
    res2 += (pred - d) ** 2
    dat2 += d ** 2
    chi2 += (pred - d) ** 2 / np.maximum(sol.var_U[:, idx[k]], 1e-30)
mis_sta = np.sqrt(res2.reshape(-1, 3).sum(1) / dat2.reshape(-1, 3).sum(1))
chi_sta = np.sqrt(chi2.reshape(-1, 3).sum(1) / (3 * K))

stations = pd.DataFrame({
    "station": names, "lon": st_llh[:, 0], "lat": st_llh[:, 1], "height_km": st_llh[:, 2],
    "frac_misfit": mis_sta, "chi_rms": chi_sta,
    "dist_to_deep_band_km": DH[deep, :].min(0),
    "within_150km_of_deep_band": near_deep,
    "seg": seg_of(st_llh[:, 1]),
})
print("    station frac_misfit by segment (median): "
      + str(stations.groupby("seg")["frac_misfit"].median().round(3).to_dict()))

lat_rows = []
for lo in range(39, 51):
    m = (stations["lat"] >= lo) & (stations["lat"] < lo + 1)
    mt = (lat >= lo) & (lat < lo + 1)
    lat_rows.append({
        "lat_band": lo,
        "n_stations": int(m.sum()),
        "n_stations_within_150km_of_deep_band": int((m & stations["within_150km_of_deep_band"]).sum()),
        "median_frac_misfit": float(stations.loc[m, "frac_misfit"].median()) if m.sum() else np.nan,
        "median_chi_rms": float(stations.loc[m, "chi_rms"].median()) if m.sum() else np.nan,
        "n_tri": int(mt.sum()),
        "n_tri_deep": int((mt & deep).sum()),
        "median_nearest_station_km_deep_tri": (float(np.median(DH[mt & deep].min(1)))
                                               if (mt & deep).sum() else np.nan),
        "median_r30_deep_tri": (float(np.median(out.loc[mt & deep, "r30"]))
                                if (mt & deep).sum() else np.nan),
    })
latdf = pd.DataFrame(lat_rows)
print(latdf.round(3).to_string(index=False))

# ============================================================= A5 synthetic recovery
print("[6] A5 synthetic recovery through the real inversion (no 13-pattern restriction)")
opt = sol.options["inversion"]
lam = float(np.ravel(opt["lambda_dip"])[0])
lam0 = float(np.ravel(opt["lambda0"])[0])
which = str(np.ravel(opt["which_smoothing"])[0])
sigmas = np.ravel(np.asarray(opt["sigma"], dtype=float))
mis_grid = np.asarray(sol.misfit_comps, dtype=float)
ind_sig = np.nanargmin(mis_grid, axis=0)
sig0_used = sigmas[ind_sig]
pe.loc[: K - 1, "sigma0_used"] = sig0_used
print(f"    lambda {lam:.2f} km ({which}); per-component sigma0 {np.unique(sig0_used)}")

Cspatial = (np.exp(-0.5 * DTRI ** 2 / lam) if which.lower().startswith("gauss")
            else np.exp(-DTRI / lam))

SIG_SYN = 30.0
U30 = PATCH[SIG_SYN]
Mtrue = np.zeros((2 * NTRI, NTRI))
Mtrue[NTRI:, :] = -U30                      # dip-slip = -thrust
Dsyn = G @ Mtrue                            # 639 x 742 displacement per unit peak slip

POT = np.empty((K, NTRI))
PEAK = np.empty((K, NTRI))
POSTSD = np.empty((K, NTRI))
for k in range(K):
    e = sig0_used[k] * lam0 / lam
    Cm0 = np.zeros((2 * NTRI, 2 * NTRI))
    Cm0[:NTRI, :NTRI] = e ** 2 * Cspatial
    Cm0[NTRI:, NTRI:] = e ** 2 * Cspatial
    Cd = np.diag(sol.var_U[:, idx[k]])
    GC = G @ Cm0
    A_ = GC @ G.T + Cd
    X_ = np.linalg.solve(A_, GC)                       # 639 x 1484
    Mrec = GC.T @ np.linalg.solve(A_, Dsyn)
    thr = -Mrec[NTRI:, :]
    POT[k] = (area @ thr) / (area @ U30)
    PEAK[k] = np.diag(thr)
    post_var = np.diag(Cm0)[NTRI:] - np.einsum("ij,ij->j", GC[:, NTRI:], X_[:, NTRI:])
    POSTSD[k] = np.sqrt(np.maximum(post_var, 0.0))

det = 2.0 * POSTSD / np.maximum(PEAK, 1e-12)
det_rel = det / np.median(det[:, ST_N], axis=1, keepdims=True)   # per component, north = 1
out["res_pot_frac"] = np.median(POT, 0)
out["res_peak_frac"] = np.median(PEAK, 0)
out["det_threshold_rel"] = np.median(det_rel, 0)
out["max_stn_disp_mm_per_cm"] = 10.0 * np.linalg.norm(
    Dsyn.reshape(-1, 3, NTRI), axis=1).max(0)
out["n_stn_gt_0p5mm_per_cm"] = (10.0 * np.linalg.norm(
    Dsyn.reshape(-1, 3, NTRI), axis=1) > 0.5).sum(0)

for nm_, m in (("central TREMOR-ONLY", TO_C), ("north SLIP+TREMOR", ST_N)):
    print(f"    {nm_:22s} res_pot {np.median(out.loc[m,'res_pot_frac']):.3f}  "
          f"res_peak {np.median(out.loc[m,'res_peak_frac']):.3f}  "
          f"det_rel {np.median(out.loc[m,'det_threshold_rel']):.3f}  "
          f"max_disp {np.median(out.loc[m,'max_stn_disp_mm_per_cm']):.2f} mm/cm  "
          f"n_stn>0.5mm {np.median(out.loc[m,'n_stn_gt_0p5mm_per_cm']):.0f}")

# ============================================================== A6 class robustness
print("[7] A6 does correcting for representability dissolve the S7 TREMOR-ONLY class?")
out["rank_S_raw"] = np.nan
out["rank_S_corr"] = np.nan
for s, m in masks.items():
    out.loc[m, "rank_S_raw"] = rank_pct(out.loc[m, "S_tot_mm"].values)
    out.loc[m, "rank_S_corr"] = rank_pct(
        (out.loc[m, "S_tot_mm"] / out.loc[m, "pot_frac30"]).values)
print(f"    central TREMOR-ONLY: within-segment slip rank raw "
      f"{np.median(out.loc[TO_C,'rank_S_raw']):.1f} -> representability-corrected "
      f"{np.median(out.loc[TO_C,'rank_S_corr']):.1f}  "
      f"(range {out.loc[TO_C,'rank_S_corr'].min():.1f}-{out.loc[TO_C,'rank_S_corr'].max():.1f})")
up = (out.loc[TO_C, "rank_S_corr"] > out.loc[TO_C, "rank_S_raw"]).sum()
still = (out.loc[TO_C, "rank_S_corr"] < 25).sum()
print(f"    {up} of {TO_C.sum()} rise, {TO_C.sum()-up} fall (the four with pot_frac30 > 0.85 "
      f"lose ground to neighbours that gain more); the WORST case is rank "
      f"{out.loc[TO_C,'rank_S_corr'].max():.1f}, so {still} of {TO_C.sum()} stay below the "
      f"25th percentile even after the basis handicap is fully undone")
print(out.loc[TO_C, ["tri", "lat", "depth_km", "S_tot_mm", "pot_frac30", "r30",
                     "rank_S_raw", "rank_S_corr", "res_pot_frac", "det_threshold_rel",
                     "nearest_station_km", "n_stn_100km"]].round(3).to_string(index=False))

# ===================================================================== summary table
print("[8] writing tables")
summ = []


def block(name, m):
    if m.sum() == 0:
        return
    d = out[m]
    summ.append({
        "group": name, "n_tri": int(m.sum()),
        "median_lat": d["lat"].median(), "median_depth_km": d["depth_km"].median(),
        "median_r30": d["r30"].median(), "median_r50": d["r50"].median(),
        "median_pot_frac30": d["pot_frac30"].median(),
        "median_pot_frac50": d["pot_frac50"].median(),
        "median_leak30": d["leak30"].median(),
        "median_capacity_density_rel": d["capacity_density_rel"].median(),
        "median_nearest_station_km": d["nearest_station_km"].median(),
        "median_n_stn_100km": d["n_stn_100km"].median(),
        "median_res_pot_frac": d["res_pot_frac"].median(),
        "median_res_peak_frac": d["res_peak_frac"].median(),
        "median_det_threshold_rel": d["det_threshold_rel"].median(),
        "median_max_stn_disp_mm_per_cm": d["max_stn_disp_mm_per_cm"].median(),
        "median_D3nb_tot_days": d["D3nb_tot_days"].median(),
    })


for s in ("north", "central", "south"):
    block(f"segment:{s}", masks[s])
    block(f"segment:{s}|ETS band {ETS_DEPTH[0]:.0f}-{ETS_DEPTH[1]:.0f} km",
          masks[s] & band)
for cl in ("SLIP+TREMOR", "TREMOR-ONLY", "SLIP-ONLY", "QUIET", "MIXED", "varies"):
    block(f"class:{cl}", (out["class_s7"] == cl).values)
block("S7 central-Oregon TREMOR-ONLY cluster (9)", TO_C)
block("north SLIP+TREMOR (145)", ST_N)
block("north TREMOR-ONLY (mesh-edge cluster)",
      ((out["class_s7"] == "TREMOR-ONLY") & (out["seg"] == "north")).values)
block("ALL", np.ones(NTRI, bool))
sm = pd.DataFrame(summ)

(HERE / "results").mkdir(exist_ok=True)
(HERE / "figures").mkdir(exist_ok=True)
out.round(6).to_csv(HERE / "results/s12_representability_triangles.csv", index=False)
sm.round(4).to_csv(HERE / "results/s12_summary.csv", index=False)
pe.round(6).to_csv(HERE / "results/s12_pattern_energy.csv", index=False)
stations.round(4).to_csv(HERE / "results/s12_stations.csv", index=False)
latdf.round(4).to_csv(HERE / "results/s12_misfit_by_latitude.csv", index=False)
print(sm[["group", "n_tri", "median_r30", "median_pot_frac30",
          "median_capacity_density_rel", "median_nearest_station_km",
          "median_res_pot_frac", "median_det_threshold_rel"]].round(3).to_string(index=False))

# ============================================================================ figure
print("[9] figure")
LONLIM = (lon.min() - 0.4, lon.max() + 1.4)
LATLIM = (lat.min() - 0.4, lat.max() + 0.4)


def base_map(ax, seglines=True):
    ax.set_xlim(*LONLIM)
    ax.set_ylim(*LATLIM)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(45.0)))
    if seglines:
        for yv in (44.0, 46.0):
            ax.axhline(yv, color="k", lw=0.6, ls="--", alpha=0.55)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("lon", fontsize=8)


def mark_clusters(ax):
    ax.scatter(lon[TO_C], lat[TO_C], s=70, facecolors="none", edgecolors="k", lw=1.1,
               marker="o", zorder=5)
    ax.scatter(lon[ST_N], lat[ST_N], s=6, facecolors="none", edgecolors="0.25", lw=0.35,
               marker="s", zorder=4)


def scmap(ax, vals, title, cmap="viridis", vmin=None, vmax=None, cbl=""):
    base_map(ax)
    sc = ax.scatter(lon, lat, c=vals, s=16, cmap=cmap, vmin=vmin, vmax=vmax,
                    edgecolors="none")
    mark_clusters(ax)
    cb = plt.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.ax.tick_params(labelsize=6.5)
    cb.set_label(cbl, fontsize=7)
    ax.set_title(title, fontsize=8.5)
    return sc


fig = plt.figure(figsize=(19.5, 15.5))
gs = fig.add_gridspec(3, 3, hspace=0.24, wspace=0.30)

ax = fig.add_subplot(gs[0, 0])
scmap(ax, out["r30"], f"A1  representability r, sigma=30 km\n"
      f"(fraction of a unit thrust patch the 13 fixed\npatterns can express; random vector = "
      f"{K/NTRI:.3f})", vmin=0, vmax=1, cbl="r30")
ax.set_ylabel("lat", fontsize=8)

ax = fig.add_subplot(gs[0, 1])
scmap(ax, out["pot_frac30"], "A1  potency fraction, sigma=30 km\n"
      "(share of the patch's POTENCY the basis reproduces)",
      cmap="RdYlBu", vmin=0.4, vmax=1.1, cbl="pot_frac30")

ax = fig.add_subplot(gs[0, 2])
scmap(ax, out["r50"], "A1  representability r, sigma=50 km",
      vmin=0, vmax=1, cbl="r50")

ax = fig.add_subplot(gs[1, 0])
scmap(ax, out["capacity_density_rel"],
      "A2  model slip CAPACITY density\n(amplitude-weighted sum_k |thrust_k|, "
      "relative to\nthe median over north SLIP+TREMOR)",
      cmap="magma", vmin=0, vmax=2.0, cbl="capacity / north median")
ax.set_ylabel("lat", fontsize=8)

ax = fig.add_subplot(gs[1, 1])
base_map(ax)
sc = ax.scatter(lon, lat, c=out["nearest_station_km"], s=16, cmap="cividis_r",
                vmin=0, vmax=60, edgecolors="none")
ax.scatter(stations["lon"], stations["lat"], s=11, marker="^", c="crimson",
           edgecolors="k", lw=0.25, zorder=6)
mark_clusters(ax)
cb = plt.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
cb.ax.tick_params(labelsize=6.5)
cb.set_label("nearest GNSS station (km)", fontsize=7)
ax.set_title(f"A3  GNSS network ({len(names)} stations, red)\nand nearest-station distance "
             "per triangle", fontsize=8.5)

ax = fig.add_subplot(gs[1, 2])
scmap(ax, out["res_pot_frac"],
      "A5  synthetic recovery through the REAL inversion\n"
      "(GNSS + TV prior, no 13-pattern restriction):\nrecovered potency / true potency, "
      "sigma=30 km patch",
      cmap="RdYlBu", vmin=0.4, vmax=1.1, cbl="res_pot_frac")

# --- bar chart by segment / class
ax = fig.add_subplot(gs[2, 0])
groups = [("north\nSLIP+TREMOR", ST_N), ("central\nTREMOR-ONLY", TO_C),
          ("north\nall", masks["north"]), ("central\nall", masks["central"]),
          ("south\nall", masks["south"])]
metrics = [("r30", "r30", "#2166ac"), ("pot_frac30", "potency frac (s=30)", "#67a9cf"),
           ("res_pot_frac", "inversion recovery", "#b2182b"),
           ("capacity_density_rel", "capacity (rel. north)", "#ef8a62")]
w = 0.2
xs = np.arange(len(groups))
for i, (col, labl, c) in enumerate(metrics):
    vals = [np.median(out.loc[m, col]) for _, m in groups]
    ax.bar(xs + (i - 1.5) * w, vals, w, color=c, label=labl, edgecolor="k", lw=0.3)
    for x, v in zip(xs + (i - 1.5) * w, vals):
        ax.text(x, v + 0.015, f"{v:.2f}", ha="center", fontsize=5.8, rotation=90)
ax.set_xticks(xs)
ax.set_xticklabels([g for g, _ in groups], fontsize=7)
ax.axhline(1.0, color="k", lw=0.6, ls=":")
ax.set_ylim(0, 1.45)
ax.legend(fontsize=6.5, ncol=2, loc="upper center", framealpha=0.9)
ax.set_title("A1/A2/A5 by group (medians)\nbasis-limited (blue) vs network-limited (red)",
             fontsize=8.5)
ax.tick_params(labelsize=7)

# --- per-latitude bars
ax = fig.add_subplot(gs[2, 1])
lb = latdf[latdf["n_tri_deep"] > 0]
ax.bar(lb["lat_band"] + 0.5, lb["n_stations_within_150km_of_deep_band"], 0.8,
       color="#4393c3", edgecolor="k", lw=0.3, label="stations <150 km of the 30-45 km band")
ax.bar(lb["lat_band"] + 0.5, lb["n_stations"] - lb["n_stations_within_150km_of_deep_band"],
       0.8, bottom=lb["n_stations_within_150km_of_deep_band"], color="0.8",
       edgecolor="k", lw=0.3, label="other stations in the band")
ax.axvspan(44, 46, color="orange", alpha=0.15, zorder=0)
ax.set_xlabel("latitude band", fontsize=8)
ax.set_ylabel("GNSS stations", fontsize=8)
ax2 = ax.twinx()
ax2.plot(lb["lat_band"] + 0.5, lb["median_nearest_station_km_deep_tri"], "o-",
         color="crimson", ms=3.5, lw=1.1, label="median nearest-station distance, deep tri")
ax2.set_ylabel("km", fontsize=8, color="crimson")
ax2.tick_params(labelsize=7, colors="crimson")
ax2.set_ylim(0, 40)
ax.tick_params(labelsize=7)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=6.2, loc="upper left", framealpha=0.9)
ax.set_title("A3  GNSS coverage by latitude\n(orange = central Oregon 44-46N)", fontsize=8.5)

# --- misfit vs latitude
ax = fig.add_subplot(gs[2, 2])
ax.scatter(stations["frac_misfit"], stations["lat"], s=13, c="0.4", edgecolors="none",
           label="station")
mm = latdf[latdf["n_stations"] > 0]
ax.plot(mm["median_frac_misfit"], mm["lat_band"] + 0.5, "o-", color="crimson", ms=4,
        lw=1.2, label="1-deg median")
ax.axhspan(44, 46, color="orange", alpha=0.15, zorder=0)
ax.set_ylim(*LATLIM)
ax.set_xlabel("fractional misfit of the ICA spatial response\n"
              "(all 13 retained components, per station)", fontsize=8)
ax.set_ylabel("lat", fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=6.5, loc="lower left", framealpha=0.9)
ax.set_title("A4  data misfit vs latitude\n(no central-Oregon anomaly)", fontsize=8.5)

leg = [Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
              markeredgecolor="k", markersize=8,
              label="S7 central-Oregon TREMOR-ONLY cluster (9 tri)"),
       Line2D([0], [0], marker="s", color="none", markerfacecolor="none",
              markeredgecolor="0.25", markersize=5,
              label="north SLIP+TREMOR (145 tri)")]
fig.legend(handles=leg, fontsize=8, loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 0.945))
fig.suptitle("S12 — can the 13-pattern vbICA slip model place slip under the central-Oregon "
             "tremor cluster?", fontsize=12, y=0.965)
fig.savefig(HERE / "figures/s12_representability.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("    wrote figures/s12_representability.png")
print("DONE")
