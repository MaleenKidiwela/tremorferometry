#!/usr/bin/env python
"""Merlin steps 2-3 prerequisites (notes/2026-09-03_merlin_sse_plan_ruling.md).
STEP 2 reconciliation (before any stress test): for the Sep-2024 northern episode, the ON-PRE contrast of shear stress
  from (i) the stored trailing-frame field, (ii) E @ delta A_trail (must equal (i) to float32), (iii) E @ delta A_raw
  (un-averaged); the ratio (iii)/(ii) must match the same ratio on the slip series within 20 %, else build error.
STEP 3 predictors: strain Green's functions at EVERY cell at 0.5, 10 and 20 km above the mesh interface, stored as the
  13-vector per (cell, offset, 6 strain components) so any window/operator can be applied later; local receiver
  orientation = nearest mesh patch strike/dip (rake 90) -> shear/normal/Coulomb 13-vectors, plus the fixed 0/12/90.
Output: res_catalog_sse/strain_gf_offsets.npz  (E[off, cell, 6, 13], mu, lam per cell/offset, strike/dip per cell)
"""
import sys, numpy as np, pandas as pd
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
from sse_strain import ElasticModel, load_fault_ascii, load_solution
from sse_strain import strain as st
from sse_strain.green import strain_from_components
from sse_strain.megathrust import slab_depth_interpolator
from sse_strain.geodesy import llh2localxy
SD = "fault_tomography/inversion/res_catalog_sse"
OFFS = (0.5, 10.0, 20.0)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)
days = pd.to_datetime(zc["days"]); cells = zc["cell"]; clat, clon = zc["cell_lat"], zc["cell_lon"]; M = len(cells)
A_trail = np.nan_to_num(zc["A_trail"].astype(float)); A_raw = np.nan_to_num(zc["A_raw"].astype(float)); ip = zc["patch"]
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns(); slip_k = fault.slip_enu(Ls, Ld)
em = ElasticModel.from_csv(f"{REPO}/data/cascadia_casc16_profile.csv", nu_greens=sol.nu, source="casc1.6")
d_iface = np.asarray(zc["iface_depth_km"], float); valid = np.isfinite(d_iface)
x_km, y_km = llh2localxy(clat, clon, fault.origin)
strike = np.asarray(fault.strike)[ip]; dip = np.asarray(fault.dip)[ip]
print(f"local orientation: strike median {np.median(strike):.0f}, dip median {np.median(dip):.1f} (north: strike {np.median(strike[clat>47.3]):.0f}, dip {np.median(dip[clat>47.3]):.1f})")
E = np.zeros((len(OFFS), M, 6, A_raw.shape[1]), np.float32); MU = np.zeros((len(OFFS), M)); LAM = np.zeros((len(OFFS), M))
for k, off in enumerate(OFFS):
    d_obs = np.where(valid, np.maximum(d_iface - off, 0.5), 20.0)
    pts = np.column_stack([x_km * 1e3, y_km * 1e3, -d_obs * 1e3])
    E[k] = strain_from_components(pts, fault.vertices_m(), slip_k, sol.nu, progress=False)
    MU[k] = em.mu(d_obs * 1e3); LAM[k] = em.lam(d_obs * 1e3)
    print(f"  offset {off} km: GF done", flush=True)
np.savez_compressed(f"{SD}/strain_gf_offsets.npz", E=E, offsets=np.array(OFFS), mu=MU, lam=LAM, strike=strike, dip=dip, cell=cells,
                    note="E[off,cell,6,K]: strain per unit vbICA amplitude at (d_iface - off) km depth; strain = E @ A; stress = hooke(strain, mu, lam)")

# ---- STEP 2 reconciliation on the Sep-2024 episode (Merlin windows) ----
def win(a, b): return (days >= pd.Timestamp(a)) & (days <= pd.Timestamp(b))
ON, PRE = win("2024-08-20", "2024-10-05"), win("2024-06-01", "2024-07-15")
k10 = OFFS.index(10.0)
def resolved(Ecmp, A_on_minus_pre, k, local):
    eps = np.einsum("obk,k->ob", Ecmp[k], A_on_minus_pre)[None]                        # (1, M, 6)
    sig = st.hooke(eps, MU[k][None], LAM[k][None])
    if local:
        out = np.zeros((3, M))
        for j in range(M):
            c, s_, n_ = st.coulomb_stress(sig[:, j:j + 1], float(strike[j]), float(dip[j]), 90.0, friction=0.4)
            out[:, j] = float(c[0, 0]), float(s_[0, 0]), float(n_[0, 0])
        return out
    c, s_, n_ = st.coulomb_stress(sig, 0.0, 12.0, 90.0, friction=0.4); return np.vstack([c[0], s_[0], n_[0]])
dAt = A_trail[ON].mean(0) - A_trail[PRE].mean(0); dAr = A_raw[ON].mean(0) - A_raw[PRE].mean(0)
LDc = -Ld[ip]; slip_raw = (A_raw @ LDc.T); pk = slip_raw[ON].max(0) - slip_raw[PRE].mean(0)
sel = (clat > 47.3) & (clat < 49.5) & valid & (pk >= 0.004)
print(f"\nSTEP 2 reconciliation, Sep-2024, northern cells with peak slip >= 4 mm: {sel.sum()} cells")
stored = zc["shear_stress"].astype(float); c_stored = stored[ON].mean(0) - stored[PRE].mean(0)
c_trail = resolved(E, dAt, k10, False)[1]; c_raw = resolved(E, dAr, k10, False)[1]
s_trail = (LDc @ dAt); s_raw = (LDc @ dAr)
print(f"  shear (fixed receiver, 10 km): stored {np.median(c_stored[sel])/1e3:+.3f} kPa ({100*(c_stored[sel]<0).mean():.0f}% neg) | E@dA_trail {np.median(c_trail[sel])/1e3:+.3f} | E@dA_raw {np.median(c_raw[sel])/1e3:+.3f} ({100*(c_raw[sel]<0).mean():.0f}% neg)")
print(f"  max |stored - E@dA_trail| over sel = {np.abs(c_stored[sel]-c_trail[sel]).max():.2e} Pa (float32 check)")
r_stress = np.median(c_raw[sel]) / np.median(c_trail[sel]); r_slip = np.median(s_raw[sel]) / np.median(s_trail[sel])
print(f"  ratio raw/trail: stress {r_stress:.2f} vs slip {r_slip:.2f} -> {'CONSISTENT' if abs(r_stress/r_slip-1) < 0.2 else 'INCONSISTENT (build error?)'}")
print("  local-orientation shear/normal/Coulomb medians (E@dA_raw) by offset:")
for k, off in enumerate(OFFS):
    loc = resolved(E, dAr, k, True); fix = resolved(E, dAr, k, False)
    print(f"    {off:4.1f} km: LOCAL shear {np.median(loc[1][sel])/1e3:+.2f} kPa ({100*(loc[1][sel]<0).mean():.0f}% neg), normal {np.median(loc[2][sel])/1e3:+.2f} ({100*(loc[2][sel]>0).mean():.0f}% pos), Coulomb {np.median(loc[0][sel])/1e3:+.2f} | FIXED shear {np.median(fix[1][sel])/1e3:+.2f} ({100*(fix[1][sel]<0).mean():.0f}% neg)"
          f" | corr(shear, local slip) {np.corrcoef(loc[1][sel], pk[sel])[0,1]:+.2f}")
print(f"wrote {SD}/strain_gf_offsets.npz")
