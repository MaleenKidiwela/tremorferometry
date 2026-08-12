"""inversion_ets_index.png for the COARSE common-mode grids (family-centroid fix, mirror-free).
Left: ETS-phase composite (ETS - inter-ETS) deep interface delta-beta/beta, actual grid cells (no interp).
Right: deep index vs year with the matched scrambled-year null band. Reads {RESDIR}/coarse_cm.npz."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20")
STEP = float(os.environ.get("STEP", "0.2"))
d = np.load(f"{RD}/coarse_cm.npz", allow_pickle=True)
clat = d["clat"]; clon = d["clon"]; deep = d["deep"]; wins = d["wins"]; idx = d["idx"]; NULL = float(d["null"])
ets = d["ets_map"] if "ets_map" in d.files else np.full(len(clat), np.nan)
mixed = d["mixed"] if "mixed" in d.files else np.zeros(len(clat), bool)
di = np.where(deep)[0]; plon, plat = clon[di], clat[di]

# regular cell mesh (NO interpolation; each cell at its true footprint, blank where no data)
ulon = np.round(np.arange(plon.min(), plon.max()+STEP/2, STEP), 3)
ulat = np.round(np.arange(plat.min(), plat.max()+STEP/2, STEP), 3)
lox = {v: i for i, v in enumerate(ulon)}; lax = {v: i for i, v in enumerate(ulat)}
edx = np.append(ulon-STEP/2, ulon[-1]+STEP/2); edy = np.append(ulat-STEP/2, ulat[-1]+STEP/2)
LON, LAT = np.meshgrid(edx, edy)
Z = np.full((len(ulat), len(ulon)), np.nan)
for ci in di:
    v = ets[ci]
    if np.isfinite(v):
        iy = lax.get(round(clat[ci], 3)); ix = lox.get(round(clon[ci], 3))
        if iy is not None and ix is not None: Z[iy, ix] = v
vlim = float(np.nanmax(np.abs(ets[deep]))) if np.isfinite(ets[deep]).any() else 0.3
vlim = min(max(vlim, 0.1), 0.6)

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(13, 6.6))
# --- left: ETS composite map ---
a0 = fig.add_axes([0.03, 0.08, 0.44, 0.84], projection=proj)
a0.set_extent([plon.min()-0.6, plon.max()+0.6, plat.min()-0.6, plat.max()+0.6], crs=proj)
a0.add_feature(cfeature.LAND, facecolor="#f0f0f0"); a0.add_feature(cfeature.OCEAN, facecolor="#dbeafe")
a0.coastlines(resolution="50m", lw=.6)
pc = a0.pcolormesh(LON, LAT, np.ma.masked_invalid(Z), cmap="RdBu_r", vmin=-vlim, vmax=vlim, transform=proj, edgecolors="none")
# outline depth-mixed cells (fuzzy near-30km boundary)
for ci in di:
    if mixed[ci] and np.isfinite(ets[ci]):
        a0.add_patch(plt.Rectangle((clon[ci]-STEP/2, clat[ci]-STEP/2), STEP, STEP, fill=False,
                                    edgecolor="k", lw=.6, ls=(0, (2, 1)), transform=proj, zorder=4))
a0.set_title(f"ETS-phase COMPOSITE  (ETS − inter-ETS)\ndeep interface δβ/β ({STEP}° cells, dashed = depth-mixed)", fontsize=11)
plt.colorbar(pc, ax=a0, fraction=.045, pad=.02, label="δβ/β (%, model)")
# --- right: deep index vs year + null ---
a1 = fig.add_axes([0.58, 0.13, 0.39, 0.74])
a1.axhspan(-NULL, NULL, color="grey", alpha=.2, label=f"matched null (±{NULL:.2f}%, 95th of scrambles)")
a1.plot(wins, idx, "o-", color="#1a1a2e", lw=1.8, label="deep network-mean anomaly")
a1.axhline(0, color="k", lw=.5)
a1.set_xlabel("year"); a1.set_ylabel("deep δβ/β anomaly (%, model)")
a1.set_title(f"Deep index vs year — max |idx| {np.nanmax(np.abs(idx)):.2f}%\n(within null → no resolvable secular deep change)", fontsize=11)
a1.legend(fontsize=8); a1.grid(alpha=.3)
fig.suptitle(f"Coarse common-mode DEEP inversion (family-centroid, mirror-free) — ETS composite + index [{os.path.basename(RD)}]", fontsize=12.5, y=.99)
fig.savefig(f"{RD}/inversion_ets_index.png", dpi=130, bbox_inches="tight")
print(f"wrote {RD}/inversion_ets_index.png (deep cells {int(deep.sum())}, {STEP}deg)")
