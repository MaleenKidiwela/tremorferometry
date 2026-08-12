"""Per-year WHOLE-INTERFACE delta-beta/beta maps from the joint 4-D inversion (all cells, NO depth cut).
Actual grid cells (pcolormesh, no interpolation, no gridlines). Reads {RESDIR}/inversion_4d.npz."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20"); STEP = float(os.environ.get("STEP", "0.2"))
d = np.load(f"{RD}/inversion_4d.npz", allow_pickle=True)
wins = list(d["wins"]); MODEL = d["MODEL"]; clat = d["clat"]; clon = d["clon"]
ill = np.isfinite(MODEL).any(0)                                  # ALL illuminated interface cells (no depth cut)
di = np.where(ill)[0]; plon, plat = clon[di], clat[di]
ulon = np.round(np.arange(plon.min(), plon.max()+STEP/2, STEP), 3); ulat = np.round(np.arange(plat.min(), plat.max()+STEP/2, STEP), 3)
lox = {v: i for i, v in enumerate(ulon)}; lax = {v: i for i, v in enumerate(ulat)}
edx = np.append(ulon-STEP/2, ulon[-1]+STEP/2); edy = np.append(ulat-STEP/2, ulat[-1]+STEP/2); LON, LAT = np.meshgrid(edx, edy)
proj = ccrs.PlateCarree(); vlim = 0.4; yrs = list(range(2010, 2026))
fig, ax = plt.subplots(4, 4, figsize=(15, 17), subplot_kw={"projection": proj})
for k, y in enumerate(yrs):
    a = ax.flat[k]; wi = wins.index(y); Z = np.full((len(ulat), len(ulon)), np.nan)
    for ci in di:
        v = MODEL[wi][ci]
        if np.isfinite(v):
            iy = lax.get(round(clat[ci], 3)); ix = lox.get(round(clon[ci], 3))
            if iy is not None and ix is not None: Z[iy, ix] = v
    a.set_extent([plon.min()-0.6, plon.max()+0.6, plat.min()-0.6, plat.max()+0.6], crs=proj)
    a.add_feature(cfeature.LAND, facecolor="#f0f0f0"); a.add_feature(cfeature.OCEAN, facecolor="#dbeafe"); a.coastlines(resolution="50m", lw=.6)
    pc = a.pcolormesh(LON, LAT, np.ma.masked_invalid(Z), cmap="RdBu_r", vmin=-vlim, vmax=vlim, transform=proj, edgecolors="none")
    gl = a.gridlines(draw_labels=(k % 4 == 0 or k >= 12), lw=0); gl.top_labels = gl.right_labels = False; gl.xlines = gl.ylines = False
    a.set_title(f"{y}", fontsize=11)
fig.suptitle(f"Joint 4-D inversion — WHOLE-interface δβ/β anomaly by year (rel. 2019-24, {STEP}° cells, all 189 stations, no depth cut, no interp)", fontsize=13, y=.995)
cax = fig.add_axes([1.01, 0.3, 0.012, 0.4]); plt.colorbar(pc, cax=cax, label="δβ/β anomaly (%, model)")
fig.tight_layout(); fig.savefig(f"{RD}/inversion_4d_interface_maps.png", dpi=110, bbox_inches="tight")
print(f"wrote {RD}/inversion_4d_interface_maps.png ({len(di)} interface cells, {STEP}deg)")
