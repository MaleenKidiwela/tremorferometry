"""Yearly SURFACE (near-station) dv/v maps = the joint 4-D inversion's per-station SITE TERMS, interpolated
across space by year. This is the near-RECEIVER shallow field the site terms absorbed (complement of the
near-source interface maps). Masked beyond RAD km of a station (no wild extrapolation). Reads {RESDIR}/inversion_4d.npz."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import cartopy.crs as ccrs, cartopy.feature as cfeature
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20"); RAD = float(os.environ.get("RAD", "45"))
d = np.load(f"{RD}/inversion_4d.npz", allow_pickle=True)
wins = list(d["wins"]); SITES = d["SITES"]; slat = d["sta_lat"]; slon = d["sta_lon"]
proj = ccrs.PlateCarree()
glon = np.arange(slon.min()-0.3, slon.max()+0.3+1e-9, 0.1); glat = np.arange(slat.min()-0.3, slat.max()+0.3+1e-9, 0.1)
GLO, GLA = np.meshgrid(glon, glat)
lat0, lon0 = slat.mean(), slon.mean()
def xy(la, lo): return np.column_stack([(np.asarray(lo)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(la)-lat0)*111])
gxy = xy(GLA.ravel(), GLO.ravel())
def surface(vals, la, lo):
    Z = griddata(np.column_stack([lo, la]), vals, (GLO, GLA), method="linear")
    dmin, _ = cKDTree(xy(la, lo)).query(gxy); Z[(dmin.reshape(GLO.shape) > RAD)] = np.nan
    return Z
allv = SITES[np.isfinite(SITES)]; VL = min(max(float(np.nanpercentile(np.abs(allv), 97)), 0.1), 0.6)
ext = [slon.min()-0.5, slon.max()+0.5, slat.min()-0.4, slat.max()+0.4]; yrs = list(range(2010, 2026))
fig, ax = plt.subplots(4, 4, figsize=(15, 17), subplot_kw={"projection": proj})
for k, y in enumerate(yrs):
    a = ax.flat[k]; wi = wins.index(y); v = SITES[wi]; ok = np.isfinite(v)
    a.set_extent(ext, crs=proj); a.add_feature(cfeature.LAND, facecolor="#f0f0f0"); a.add_feature(cfeature.OCEAN, facecolor="#dbeafe"); a.coastlines(resolution="50m", lw=.5)
    if ok.sum() >= 8:
        Z = surface(v[ok], slat[ok], slon[ok])
        pc = a.pcolormesh(glon, glat, np.ma.masked_invalid(Z), cmap="PuOr_r", vmin=-VL, vmax=VL, transform=proj, shading="auto")
    a.scatter(slon, slat, s=4, c="k", marker="^", transform=proj, zorder=5); a.set_title(f"{y} ({int(ok.sum())} stns)", fontsize=10)
fig.suptitle(f"Joint 4-D inversion — yearly SURFACE (near-station SITE-TERM) dv/v anomaly (rel. 2019-24, ±{RAD:.0f} km mask)", fontsize=13, y=.995)
cax = fig.add_axes([1.01, 0.3, 0.012, 0.4]); plt.colorbar(pc, cax=cax, label="site-term dv/v anomaly (%)")
fig.tight_layout(); fig.savefig(f"{RD}/inversion_4d_surface_maps.png", dpi=110, bbox_inches="tight")
print(f"wrote {RD}/inversion_4d_surface_maps.png")
