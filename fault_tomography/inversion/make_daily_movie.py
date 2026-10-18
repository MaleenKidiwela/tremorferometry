"""Movie of the DAILY joint 4-D inversion (inversion_4d_daily.npz -> inversion_4d_daily_movie.mp4).

Weekly-sampled frames of the interface delta-beta/beta field (true sign, rel. 2019-24 gauge) on the
~8 km production cell grid, with a deep-index (>=30 km illuminated-mean) timeline + cursor below.
Frames are the raw 30-d trailing product: per-frame posterior sigma_m ~0.7% median, so single-frame
texture is mostly noise — the movie is for coverage/QC and slow structure, not per-frame claims.
Sub-30-d structure = artifact by construction (I6)."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cartopy.crs as ccrs, cartopy.feature as cfeature

RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
STRIDE = int(os.environ.get("STRIDE", "7"))          # days between frames
FPS = int(os.environ.get("FPS", "12"))
VLIM = float(os.environ.get("VLIM", "1.0"))          # % color scale (frame sigma_m ~0.7%)

d = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
MODEL, SG = d["MODEL"], d["SIGM"]
clat, clon, dep = d["clat"], d["clon"], d["depth"]
days = pd.to_datetime(d["days"])
deep = dep >= 30.0

# CORRECTED=1 (2026-08-22, user decision): transmission-normalized movie, 2009-10 excluded.
# Cell-level T = per-band interpolation of the daily patch T between patch centroids —
# VISUALIZATION-GRADE per Merlin r2 item 6a (quantitative product = patch series in
# transmission_deconv.npz). Deep-index strip = deconvolved precision-weighted aggregate.
CORR = os.environ.get("CORRECTED", "0") == "1"
BAND = os.environ.get("BAND", "all")                 # all | deep | shallow
if CORR:
    zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
    zd = np.load(f"{RD}/transmission_deconv.npz", allow_pickle=True)
    xy = np.load(f"{RD}/G.npz", allow_pickle=True)["cxy"]
    ac = xy - xy.mean(0)
    u = np.linalg.eigh(ac.T @ ac)[1][:, -1]; s_along = ac @ u
    Tpat = np.asarray(zt["T"], float)
    psc = np.asarray(zt["patch_sc"], float); pband = zt["patch_band"]
    def cell_T(t):
        out = np.full(len(dep), np.nan)
        for bd, bmask in (("shallow", ~deep), ("deep", deep)):
            qi = np.where(pband == bd)[0]; o = np.argsort(psc[qi])
            out[bmask] = np.interp(s_along[bmask], psc[qi][o], Tpat[t, qi][o])
        return out
    keepf = days >= pd.Timestamp("2011-01-01")       # user decision: 2009-10 excluded
else:
    keepf = np.ones(len(days), bool)
ti = np.arange(0, len(days), STRIDE)
ti = ti[keepf[ti]]

# which cells are drawn
show = deep if BAND == "deep" else (~deep if BAND == "shallow" else np.ones(len(dep), bool))

# index strip: corrected -> precision-weighted deconvolved aggregate of the band; else naive mean
if CORR:
    A = np.asarray(zd["A"], float); pb = zd["patch_band"]
    if BAND == "shallow":
        qi = np.where(pb == "shallow")[0]
        S = np.asarray(zd["SDEC"], float)[1:].mean(0)[np.ix_(qi, qi)]   # mature epochs
        wi = np.linalg.solve(S, np.ones(len(qi))); w = wi / wi.sum()
        idx_raw = A[:, qi] @ w
        strip_lab = "shallow index (%)"
    else:
        idx_raw = np.asarray(zd["agg"], float)       # frozen deep weights
        strip_lab = "deep index (%)"
else:
    idx_raw = np.nanmean(np.where(np.isfinite(MODEL[:, deep]), MODEL[:, deep], np.nan), axis=1)
    strip_lab = "deep index (%)"
idx_s = pd.Series(idx_raw, index=days).rolling("90D", min_periods=30).mean()

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(7.2, 9.6))
axm = fig.add_axes([0.06, 0.24, 0.88, 0.72], projection=proj)
axt = fig.add_axes([0.10, 0.05, 0.84, 0.13])
axm.set_extent([clon.min() - 0.6, clon.max() + 0.6, clat.min() - 0.5, clat.max() + 0.5], crs=proj)
axm.add_feature(cfeature.LAND, facecolor="#f0f0f0")
axm.add_feature(cfeature.OCEAN, facecolor="#dbeafe")
axm.coastlines(resolution="50m", lw=0.6)

# marker size ~ one 7.8-km cell (0.07 deg) in points^2
bb = axm.get_window_extent(fig.canvas.get_renderer())
pts_per_deg = bb.height / (clat.max() - clat.min() + 1.0) / fig.dpi * 72
s_cell = max((pts_per_deg * 0.075) ** 2, 4)
sc = axm.scatter(clon, clat, c=np.zeros(len(clon)), s=s_cell, marker="s", cmap="RdBu_r",
                 vmin=-VLIM, vmax=VLIM, transform=proj, edgecolors="none", zorder=3)
axm.scatter(clon[deep], clat[deep], s=s_cell * 0.15, marker=".", c="none", zorder=4)  # placeholder layer
cb = fig.colorbar(sc, ax=axm, shrink=0.6, pad=0.02)
cb.set_label("interface δβ/β (%, rel. 2019–24)")
ttl = axm.set_title("", fontsize=12)

axt.plot(days, idx_s.values, lw=1.0, color="#334155")
axt.axhline(0, lw=0.5, color="#94a3b8")
axt.set_ylabel(strip_lab, fontsize=8)
axt.tick_params(labelsize=7)
lo = np.nanpercentile(idx_s, 1); hi = np.nanpercentile(idx_s, 99)
pad = 0.2 * (hi - lo + 1e-9); axt.set_ylim(lo - pad, hi + pad)
cursor = axt.axvline(days[0], color="#dc2626", lw=1.2)

def update(k):
    t = ti[k]
    v = MODEL[t].astype(float)
    if CORR:
        Tc = cell_T(t)
        v = v / np.clip(Tc, 0.1, None)
        v[~(Tc >= 0.1)] = np.nan                     # un-correctable cells masked
    good = np.isfinite(v) & show
    colors = np.where(good, v, np.nan)
    sc.set_array(colors)
    sc.set_alpha(None)
    fc = sc.cmap((np.clip(colors, -VLIM, VLIM) + VLIM) / (2 * VLIM))
    fc[~good & show] = (0.85, 0.85, 0.85, 0.5)       # lit-but-masked = gray
    fc[~show] = (0, 0, 0, 0)                         # other band = invisible
    sc.set_facecolors(fc)
    cursor.set_xdata([days[t], days[t]])
    tag = "transmission-corrected (viz-grade)" if CORR else "raw"
    ttl.set_text(f"Daily joint inversion [{BAND}, {tag}] — {days[t].date()}  "
                 f"(30-d trailing, {good.sum()} cells)")
    return sc, cursor, ttl

anim = FuncAnimation(fig, update, frames=len(ti), blit=False)
suffix = ("_corrected" if CORR else "") + ("" if BAND == "all" else f"_{BAND}")
out = f"{RD}/inversion_4d_daily_movie{suffix}.mp4"
anim.save(out, writer=FFMpegWriter(fps=FPS, bitrate=2500), dpi=110)
print(f"wrote {out} | {len(ti)} frames (stride {STRIDE} d) @ {FPS} fps = {len(ti)/FPS:.0f} s")
