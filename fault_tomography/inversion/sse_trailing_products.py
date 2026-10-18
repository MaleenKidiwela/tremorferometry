#!/usr/bin/env python
"""Slow-slip model products in OUR measurement frame (user 2026-09-02): the repo's 31-day CENTRED average on the vbICA
time functions is replaced by our dv/v operator, the 30-day TRAILING mean over (t-30, t], equal weight per day, and
the reference is our inversion gauge (zero mean over 2019-01-01..2024-12-31). Because every field is
sum_k pattern_k(x) * a_k(t), smoothing the K amplitude functions is exactly smoothing every field.

Uses cascadia-sse-strain (github MaleenKidiwela) unmodified: Gualandi 2026-08-18 solution, Slab2 mesh, casc1.6
elastic profile, receiver (strike 0, dip 12, rake 90, friction 0.4) as in its committed run.

Outputs (fault_tomography/inversion/res_catalog_sse/, new files only):
  components_at_cells_trailing30.npz  daily 2007-2026 at our 872 interface cells (10 km above the interface):
      strain_e{xx,yy,zz,xy,xz,yz} (1), stress_s{..} (Pa), dilatation, eps_eq, max_shear_strain, mean_stress, pressure,
      von_mises_stress, max_shear_stress, coulomb, shear_stress, normal_stress, slip_thrust_cum (m, thrust-positive,
      trailing frame), slip_thrust_inc30 (m per 30 d), slip_strike_cum, plus amplitudes A_trail (T x K)
  megathrust_strain_trailing30.nc     the repo's full interface-following grid (0.1 deg, weekly), same fields, our
      operator + gauge, masked + 21-km low-pass exactly as the repo's script does
  sse_trailing_check.png              latitude-time of thrust slip increment and Coulomb on the deep band (sanity)
"""
import sys, os, numpy as np, pandas as pd
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
from sse_strain import ElasticModel, dataset, load_fault_ascii, load_solution
from sse_strain import strain as st
from sse_strain.green import strain_from_components
from sse_strain.megathrust import make_slab_grid, resolution_filter, slab_depth_interpolator
from sse_strain.geodesy import llh2localxy
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
R0, R1 = np.datetime64("2019-01-01"), np.datetime64("2024-12-31")
WIN, OFFSET_KM, RECEIVER, FRICTION = 30, 10.0, (0.0, 12.0, 90.0), 0.4

fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt")
sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns(); slip_k = fault.slip_enu(Ls, Ld)          # (n_patch, 3, K)
em = ElasticModel.from_csv(f"{REPO}/data/cascadia_casc16_profile.csv", nu_greens=sol.nu,
                           source="casc1.6 forearc median profile (Stephenson et al., 2017)")

# ---- our operator on the K amplitude functions: daily grid, trailing 30-d mean, gauge 2019-24 ----
A = sol.component_amplitudes(); t = pd.to_datetime(np.asarray(sol.dates))
Ad = pd.DataFrame(A, index=t).groupby(level=0).mean()
days = pd.date_range(Ad.index.min(), Ad.index.max(), freq="D")
Ad = Ad.reindex(days).interpolate(limit_area="inside")
A_raw = Ad.values.copy()                                                             # daily, un-averaged
A_trail = Ad.rolling(f"{WIN}D", min_periods=WIN - 2).mean().values                # (t-30, t] equal weight per day
gauge = (days >= R0) & (days <= R1)
A_trail = A_trail - np.nanmean(A_trail[gauge], axis=0)
K = A_trail.shape[1]; T = len(days)
print(f"amplitudes: {T} days x {K} components | trailing {WIN}-d mean, gauge {str(R0)[:10]}..{str(R1)[:10]}")

# ---- our 872 cells as observation points (10 km above the mesh interface, mesh's own depth) ----
gz = np.load(f"{RD}/G.npz", allow_pickle=True)
cells, clat, clon = gz["cell"], gz["cell_lat"].astype(float), gz["cell_lon"].astype(float)
d_iface = slab_depth_interpolator(fault)(clon, clat)                              # km, positive down
d_iface = np.where(np.isfinite(d_iface), d_iface, np.nan)
valid = np.isfinite(d_iface)
print(f"cells: {len(cells)} | interface depth from mesh finite at {valid.sum()} | median {np.nanmedian(d_iface):.1f} km")
x_km, y_km = llh2localxy(clat, clon, fault.origin)
d_obs = np.where(valid, d_iface - OFFSET_KM, 20.0)
pts = np.column_stack([x_km * 1e3, y_km * 1e3, -np.maximum(d_obs, 0.5) * 1e3])
print("strain Green's functions at cells ...", flush=True)
E = strain_from_components(pts, fault.vertices_m(), slip_k, sol.nu, progress=False)  # (872, 6, K)
eps = np.einsum("obk,tk->tob", E, np.nan_to_num(A_trail))                          # (T, 872, 6)
mu = em.mu(d_obs * 1e3); lam = em.lam(d_obs * 1e3)
sig = st.hooke(eps, np.broadcast_to(mu, (T, len(cells))), np.broadcast_to(lam, (T, len(cells))))
inv = st.stress_invariants(sig)
dcfs, shear, normal = st.coulomb_stress(sig, *RECEIVER, friction=FRICTION)
out = {}
for i, c in enumerate(("exx", "eyy", "ezz", "exy", "exz", "eyz")):
    out[f"strain_{c}"] = eps[..., i]; out[f"stress_s{c[1:]}"] = sig[..., i]
out.update(dilatation=st.dilatation(eps), eps_eq=st.equivalent_shear_strain(eps), max_shear_strain=st.max_shear_strain(eps),
           mean_stress=inv["mean_stress"], pressure=inv["pressure"], von_mises_stress=inv["von_mises_stress"],
           max_shear_stress=inv["max_shear_stress"], coulomb=dcfs, shear_stress=shear, normal_stress=normal)
# slip at cells (nearest patch), trailing frame, thrust-positive
from scipy.spatial import cKDTree
cll = fault.centroid_llh; px, py = llh2localxy(cll[:, 1], cll[:, 0], fault.origin)
dist, ip = cKDTree(np.column_stack([px, py])).query(np.column_stack([x_km, y_km]))
ss_cum = (Ls @ A_trail.T)[ip]; ds_cum = (Ld @ A_trail.T)[ip]                       # (872, T) metres, trailing frame
# sign check: the model's ETS months (e.g. 2024-09) must be thrust-POSITIVE after the flip
sm = sol.slip_magnitude(apply_rake_flip=True)                                        # stored: ETS negative (checked)
sign = -1.0
thrust_cum = sign * ds_cum
raw_cum = sign * (Ld @ np.nan_to_num(A_raw).T)[ip]                                   # un-averaged cumulative thrust slip
thrust_inc = np.full_like(thrust_cum, np.nan); thrust_inc[:, WIN:] = raw_cum[:, WIN:] - raw_cum[:, :-WIN]   # slip within (t-30, t] = 30-d trailing mean of the rate
na = np.array([(days >= np.datetime64("2024-08-15")) & (days <= np.datetime64("2024-10-15"))]).ravel()
north = (clat > 47.3) & (clat < 48.8)
print(f"sign check: mean thrust 30-d increment, northern cells, Sep 2024 = {np.nanmean(thrust_inc[north][:, na]) * 1000:+.1f} mm "
      f"(must be positive)")
for k in list(out): out[k] = out[k].astype(np.float32)
np.savez_compressed(f"{SD}/components_at_cells_trailing30.npz", days=days.astype(str).values, cell=cells, cell_lat=clat, cell_lon=clon,
                    iface_depth_km=d_iface, obs_depth_km=d_obs, valid=valid, patch=ip, patch_dist_km=dist,
                    slip_thrust_cum=thrust_cum.T.astype(np.float32), slip_thrust_inc30=thrust_inc.T.astype(np.float32),
                    slip_strike_cum=ss_cum.T.astype(np.float32), A_trail=A_trail.astype(np.float32), A_raw=A_raw.astype(np.float32),
                    window_days=WIN, gauge=f"{str(R0)[:10]}..{str(R1)[:10]}", offset_km=OFFSET_KM, receiver=np.array(RECEIVER),
                    friction=FRICTION, note="trailing 30-d mean on vbICA amplitudes (all fields + slip_thrust_cum); slip_thrust_inc30 = raw slip accumulated within (t-30,t] = 30-d trailing mean of the slip rate; thrust = -dip_slip (ETS positive); units strain 1, stress Pa, slip m", **out)
print(f"wrote {SD}/components_at_cells_trailing30.npz  ({len(out)} fields + slip)")

# ---- full-grid netCDF (repo grid, weekly) in the trailing frame ----
grid = make_slab_grid(fault, spacing_deg=0.1, offset_km=OFFSET_KM, depth_limits_km=(25.0, 50.0))
print(f"grid {grid.lon.size} x {grid.lat.size}, valid {int(grid.valid.sum())} — Green's functions ...", flush=True)
Eg = strain_from_components(grid.points_m, fault.vertices_m(), slip_k, sol.nu, progress=False)
grid.depth_m = np.array([np.nanmean(grid.observation_depth_km) * 1e3])
keep = np.isfinite(A_trail).all(1)
ds = dataset.evaluate(grid, sol, em, Eg, amplitudes=A_trail[keep], times=days.values[keep], stride=7,
                      tensor_output=True, receiver=RECEIVER, friction=FRICTION)
mask = grid.valid[None, None, :, :]
for v in ds.data_vars:
    if {"lat", "lon"} <= set(ds[v].dims): ds[v] = ds[v].where(mask)
lam_km = 21.0
for v in list(ds.data_vars):
    if {"lat", "lon"} <= set(ds[v].dims):
        a = ds[v].values; sh = a.shape
        ds[v] = (ds[v].dims, resolution_filter(a.reshape((-1,) + sh[-2:]), grid, lam_km).reshape(sh).astype("float32"), ds[v].attrs)
ds = ds.assign_coords(interface_depth=(("lat", "lon"), grid.interface_depth_km), observation_depth=(("lat", "lon"), grid.observation_depth_km))
ds.attrs.update(measurement_operator=f"TRAILING {WIN}-day mean over (t-{WIN}, t], equal weight per day (dv/v stack operator)",
                reference_interval=f"{str(R0)[:10]}/{str(R1)[:10]} (dv/v inversion gauge)", offset_above_interface_km=OFFSET_KM,
                spatial_lowpass_km=lam_km, built_by="tremorferometry/fault_tomography/inversion/sse_trailing_products.py")
dataset.write(ds, f"{SD}/megathrust_strain_trailing30.nc")
print(f"wrote {SD}/megathrust_strain_trailing30.nc")

# ---- sanity figure: latitude-time on the deep band ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
deepc = valid & (d_iface >= 30) & (d_iface <= 45)
lb = np.arange(40, 50.01, 0.25); li = np.digitize(clat, lb)
def latbin(X):
    return np.array([np.nanmean(X[:, deepc & (li == k)], 1) if (deepc & (li == k)).sum() else np.full(T, np.nan) for k in range(1, len(lb))])
fig, ax = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
for a, (X, nm, cm, sym) in zip(ax, ((thrust_inc.T * 1000, "thrust slip increment per 30 d (mm), trailing frame", "magma", False),
                                    (out["coulomb"] / 1e3, "Coulomb stress (kPa), trailing frame, gauge 2019-24", "RdBu_r", True))):
    Z = latbin(X); v = np.nanpercentile(np.abs(Z), 99)
    im = a.pcolormesh(days, lb[:-1] + 0.125, Z, cmap=cm, vmin=-v if sym else 0, vmax=v, shading="auto"); a.set_ylabel("lat"); a.set_title(nm, fontsize=10); plt.colorbar(im, ax=a, pad=.01)
ax[-1].set_xlim(pd.Timestamp("2009-01-01"), pd.Timestamp("2026-09-01"))
fig.savefig(f"{SD}/sse_trailing_check.png", dpi=110, bbox_inches="tight"); print("wrote sse_trailing_check.png")
