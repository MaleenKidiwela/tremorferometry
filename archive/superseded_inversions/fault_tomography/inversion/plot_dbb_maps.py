import warnings, os; warnings.filterwarnings("ignore")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g50")
STEP = float(os.environ.get("STEP", "0.2"))
d = np.load(f"{RD}/coarse_cm.npz", allow_pickle=True)
wins = d["wins"]; MODEL = d["MODEL"]; clat = d["clat"]; clon = d["clon"]; deep = d["deep"]
di = np.where(deep)[0]; plon, plat = clon[di], clat[di]
# regular cell mesh (NO interpolation; each cell shown at its true footprint, blank where no data)
ulon = np.round(np.arange(plon.min(), plon.max()+STEP/2, STEP), 3)
ulat = np.round(np.arange(plat.min(), plat.max()+STEP/2, STEP), 3)
lox = {v: i for i, v in enumerate(ulon)}; lax = {v: i for i, v in enumerate(ulat)}
edx = np.append(ulon-STEP/2, ulon[-1]+STEP/2); edy = np.append(ulat-STEP/2, ulat[-1]+STEP/2)
LON, LAT = np.meshgrid(edx, edy)
proj = ccrs.PlateCarree(); vlim = 0.4; yrs = list(range(2010, 2026))
fig, ax = plt.subplots(4, 4, figsize=(15, 16), subplot_kw={"projection": proj})
for k, y in enumerate(yrs):
    a = ax.flat[k]; wi = list(wins).index(y)
    Z = np.full((len(ulat), len(ulon)), np.nan)
    for ci in di:
        v = MODEL[wi][ci]
        if np.isfinite(v):
            iy = lax.get(round(clat[ci], 3)); ix = lox.get(round(clon[ci], 3))
            if iy is not None and ix is not None: Z[iy, ix] = v
    a.set_extent([plon.min()-0.6, plon.max()+0.6, plat.min()-0.6, plat.max()+0.6], crs=proj)
    a.add_feature(cfeature.LAND, facecolor="#f0f0f0"); a.add_feature(cfeature.OCEAN, facecolor="#dbeafe")
    a.coastlines(resolution="50m", lw=.6)
    pc = a.pcolormesh(LON, LAT, np.ma.masked_invalid(Z), cmap="RdBu_r", vmin=-vlim, vmax=vlim, transform=proj, edgecolors="none")
    a.set_title(f"{y}", fontsize=11)
    gl = a.gridlines(draw_labels=(k % 4 == 0 or k >= 12), lw=0); gl.top_labels = gl.right_labels = False; gl.xlines = gl.ylines = False
fig.suptitle(f"Deep plate-interface δβ/β by year ({STEP}° grid cells, common-mode inversion, mirror-free, >30 km) — no interpolation", fontsize=14, y=.995)
cax = fig.add_axes([1.01, 0.3, 0.012, 0.4]); plt.colorbar(pc, cax=cax, label="δβ/β (%, model)")
fig.tight_layout(); fig.savefig(f"{RD}/dbb_maps_by_year.png", dpi=110, bbox_inches="tight")
print(f"wrote {RD}/dbb_maps_by_year.png ({len(di)} deep cells, {STEP}deg)")
