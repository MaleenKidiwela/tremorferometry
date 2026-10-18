#!/usr/bin/env python
"""Build the matched (interface dbeta/beta  vs  modeled slow-slip) dataset on our 872 interface cells.

dv/v side: res_catalog/inversion_4d_daily.npz MODEL (days x cells, 2-4 s production joint inversion; value at date t
  is the 30-day TRAILING stack (t-30, t] of the days that had a daily stack).
slip side: res_catalog_sse/slip_at_cells.npz (Gualandi 2026-08-18 solution via cascadia-sse-strain; daily cumulative
  slip rel. long-term trend + 7-d slip rate, nearest mesh patch per cell).
TIMING CALIBRATION (user 2026-09-02): the slip predictors are passed through the SAME operator as the dv/v — for each
  cell and date t, the mean over the stack-days d in (t-30, t] (days on which at least one family of that cell has a
  row in family_days.parquet) of slip_rate(d) [transient predictor] and of slip_cum(d) [state predictor]. No fixed lag.
strain side: out/megathrust_strain.nc invariants (dilatation, max_shear_strain, coulomb; weekly, ALREADY 31-d centred
  averaged, referenced to 2025-08..2026-08) sampled at the nearest grid node, linearly interpolated to days and
  shifted by +15 d (centred -> trailing) — approximate, flagged.
Output: res_catalog_sse/matched.npz  (days, cells, MODEL, pred_rate, pred_cum, nstackdays, dil, mxs, coul, ...)
"""
import numpy as np, pandas as pd, os
RD = "fault_tomography/inversion/res_catalog"; SD = "fault_tomography/inversion/res_catalog_sse"
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
days = pd.to_datetime(zp["days"]); MODEL = np.asarray(zp["MODEL"], float); cells = zp["cell"]; T, M = MODEL.shape
zs = np.load(f"{SD}/slip_at_cells.npz", allow_pickle=True)
assert list(zs["cell"]) == list(cells)
t_cum = pd.to_datetime(zs["t_cum"]); t_rate = pd.to_datetime(zs["t_rate"])
grid = pd.date_range(min(t_cum.min(), days.min()) - pd.Timedelta(days=31), max(t_cum.max(), days.max()), freq="D")
G = len(grid); gi = {d: k for k, d in enumerate(grid)}
cum = pd.DataFrame(zs["slip_cum"].T, index=t_cum).groupby(level=0).mean().reindex(grid).interpolate(limit_area="inside").values      # (G, M)
rate = pd.DataFrame(zs["slip_rate"].T, index=t_rate).groupby(level=0).mean().reindex(grid).interpolate(limit_area="inside").values   # (G, M)
# stack-day mask per cell from the family-day table
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["cell", "date"]).drop_duplicates()
fd = fd[fd.cell.isin(set(cells))]
cpos = {c: i for i, c in enumerate(cells)}
mask = np.zeros((G, M), np.float32)
di = np.array([gi.get(pd.Timestamp(d), -1) for d in fd.date.values]); ci = np.array([cpos[c] for c in fd.cell.values])
ok = di >= 0; mask[di[ok], ci[ok]] = 1.0
# trailing-30 operator with the stack-day mask: pred(t) = sum_{d in (t-30, t]} m(d) x(d) / sum m(d)
def trailing(x):
    xm = np.nan_to_num(x) * mask; cs = np.vstack([np.zeros((1, M)), np.cumsum(xm, 0)]); cm = np.vstack([np.zeros((1, M)), np.cumsum(mask, 0)])
    out = np.full((G, M), np.nan)
    for k in range(30, G):
        n = cm[k + 1] - cm[k - 29]; s = cs[k + 1] - cs[k - 29]
        out[k] = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return out
pr, pc = trailing(rate), trailing(cum)
nsd = np.vstack([np.zeros((1, M)), np.cumsum(mask, 0)]); nsd = np.array([nsd[k + 1] - nsd[k - 29] if k >= 30 else np.zeros(M) for k in range(G)])
sel = np.array([gi[d] for d in days])
pred_rate, pred_cum, nstack = pr[sel], pc[sel], nsd[sel]
# also the plain (unmasked) trailing-30 mean and the centred 31-d mean for reference
def plain(x, w=30):
    cs = np.vstack([np.zeros((1, M)), np.cumsum(np.nan_to_num(x), 0)])
    return np.array([(cs[k + 1] - cs[k + 1 - w]) / w if k + 1 >= w else np.full(M, np.nan) for k in range(G)])
pred_rate_plain = plain(rate)[sel]
# strain invariants at cells (weekly centred -> daily, +15 d shift to trailing convention)
try:
    import xarray as xr
    ds = xr.open_dataset("/home/jovyan/cascadia-sse-strain/out/megathrust_strain.nc")
    lat, lon = ds.lat.values, ds.lon.values
    tt = pd.to_datetime(ds.time.values) + pd.Timedelta(days=15)
    out_inv = {}
    for v in ("dilatation", "max_shear_strain", "coulomb", "eps_eq"):
        A = ds[v].isel(depth=0).values                                 # (time, lat, lon)
        vals = np.full((len(tt), M), np.nan)
        for j in range(M):
            ia = int(np.argmin(np.abs(lat - zs["cell_lat"][j]))); io = int(np.argmin(np.abs(lon - zs["cell_lon"][j])))
            col = A[:, ia, io]
            if np.isnan(col).all():                                    # nearest non-nan node within 2 steps
                sub = A[:, max(ia - 2, 0):ia + 3, max(io - 2, 0):io + 3]
                col = np.nanmean(sub.reshape(len(tt), -1), axis=1)
            vals[:, j] = col
        out_inv[v] = pd.DataFrame(vals, index=tt).reindex(grid).interpolate(limit_area="inside").values[sel].astype(np.float32)
except Exception as e:
    print("strain sampling failed:", e); out_inv = {}
np.savez_compressed(f"{SD}/matched.npz", days=days.astype(str).values, cell=cells, cell_lat=zs["cell_lat"], cell_lon=zs["cell_lon"],
                    cell_depth=zs["cell_depth"], patch_dist_km=zs["patch_dist_km"], MODEL=MODEL.astype(np.float32),
                    pred_rate=pred_rate.astype(np.float32), pred_cum=pred_cum.astype(np.float32),
                    pred_rate_plain=pred_rate_plain.astype(np.float32), nstackdays=nstack.astype(np.int16),
                    **{f"inv_{k}": v for k, v in out_inv.items()})
fin = np.isfinite(MODEL) & np.isfinite(pred_rate)
print(f"matched: {T} days x {M} cells | finite pairs {fin.sum():,} | median stack-days per window {np.median(nstack[fin]):.0f} "
      f"| pred_rate (m/yr) p50/p99 {np.nanpercentile(pred_rate[fin], 50):.3f}/{np.nanpercentile(pred_rate[fin], 99):.3f} "
      f"| strain vars {list(out_inv)}")
print(f"wrote {SD}/matched.npz")
