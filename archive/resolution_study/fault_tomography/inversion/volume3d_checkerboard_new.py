#!/usr/bin/env python
"""3-D VOLUME checkerboard on the NEW full network (recreates the old smoke_volume3d_checkerboard with the
187-station / 6726-pair catalog + honest operator). Model = crustal VOLUME (lon,lat,depth), single-scatter
kernel (2-4 s), per-station site terms, REALISTIC per-datum noise (time-varying residual sigma). Sources =
interface cells at Slab2 depth; receivers = stations. Plants multi-scale 3-D checkers, inverts, scores recovery
by horizontal scale AND depth layer. CAVEAT (Merlin): one 2-4 s window => weak VERTICAL resolution; the figure
shows horizontal recovery per depth layer + the depth smearing, not clean depth separation.
Outputs res_catalog/volume3d_checkerboard_new.{csv,png}."""
import os, sys
import numpy as np, pandas as pd
from numpy.linalg import lstsq
sys.path.insert(0, "fault_tomography/kernels")
from kernel import kernel_singlescatter
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = "fault_tomography/inversion/res_catalog"
BETA, ELL, W1, W2 = 3.5, 40.0, 2.0, 4.0
HGRID, DEPTHS = 0.2, np.array([3., 12., 22., 32., 40., 48.])   # same layers as the old volume test
SCALES_KM = [70, 140, 250]
LAM = 0.5

pairs = pd.read_csv(f"{OUT}/pairs.csv")
cells = pd.read_csv(f"{OUT}/cells.csv").sort_values("cell").reset_index(drop=True)
cidx = {c: i for i, c in enumerate(cells.cell)}
lat0, lon0 = cells.cell_lat.mean(), cells.cell_lon.mean()
def proj(lat, lon):
    return np.column_stack([(np.asarray(lon)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(lat)-lat0)*111])
CXY = proj(cells.cell_lat.values, cells.cell_lon.values); CZ = cells.depth_km.values

# 3-D model grid
mlat = np.arange(cells.cell_lat.min()-0.3, cells.cell_lat.max()+0.3+1e-9, HGRID)
mlon = np.arange(cells.cell_lon.min()-0.3, cells.cell_lon.max()+0.3+1e-9, HGRID)
GLA, GLO = np.meshgrid(mlat, mlon, indexing="ij")
hlat = np.repeat(GLA.ravel(), len(DEPTHS)); hlon = np.repeat(GLO.ravel(), len(DEPTHS))
mz = np.tile(DEPTHS, GLA.size); mxy = proj(hlat, hlon)
Mv = len(mz); print(f"3-D volume cells: {Mv} ({GLA.size} cols x {len(DEPTHS)} layers)")

# realistic sigma per pair (time-varying residual) from tensor
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["site_mean"] = pm.groupby(["tag", "ym"]).dvv_month.transform("mean"); pm["resid"] = pm.dvv_month-pm.site_mean
pm["pair"] = pm.tag+"|"+pm.cell
tv = pm.groupby("pair").resid.apply(lambda r: (r-r.mean()).std() if len(r) >= 6 else np.nan)
TVM = float(np.nanmedian(tv)); sig_real = {p: (v if np.isfinite(v) and v > 0 else TVM) for p, v in tv.items()}

# forward operator over the volume (normalized single-scatter rows)
sxy = proj(pairs.sta_lat.values, pairs.sta_lon.values)
G = np.zeros((len(pairs), Mv), np.float32); sig = np.zeros(len(pairs))
for i, r in pairs.iterrows():
    ci = cidx[r.cell]
    K = kernel_singlescatter(mxy, mz, CXY[ci], float(CZ[ci]), sxy[i], BETA, W1, W2, ELL)
    G[i] = -K
    sig[i] = sig_real.get(f"{r.tag}|{r.cell}", TVM)
hit = (np.abs(G) > 0).sum(0); ill = hit >= 3
G = G[:, ill]; mxy_i, mz_i, hlat_i, hlon_i = mxy[ill], mz[ill], hlat[ill], hlon[ill]
M = G.shape[1]; print(f"illuminated volume cells (>=3 rays): {M}")

# site design + 3-D graph Laplacian (7-NN, depth stretched x3 so smoothing is weaker across layers)
tags = pairs.tag.values; ut = pd.unique(tags); tc = {t: k for k, t in enumerate(ut)}
S = np.zeros((len(pairs), len(ut)))
for r, t in enumerate(tags):
    S[r, tc[t]] = 1.0
P3 = np.column_stack([mxy_i, mz_i*3.0]); _, nb = cKDTree(P3).query(P3, k=8)
Lap = np.zeros((M, M))
for i in range(M):
    for j in nb[i, 1:]:
        Lap[i, i] += 1; Lap[i, j] -= 1

w = 1.0/sig
A = np.hstack([G, S]); Aw = A*w[:, None]
reg = np.zeros((A.shape[1], A.shape[1])); reg[:M, :M] = LAM*Lap
reg[np.arange(M), np.arange(M)] += 1e-3
si = np.arange(M, M+len(ut)); reg[si, si] += (LAM*0.125)**2
# stacked LS system solved once per checker via normal equations
AtA = Aw.T@Aw + reg; AtAinv = np.linalg.inv(AtA)


def checker3d(hkm, dkm=25.0):
    return np.sign(np.sin(np.pi*mxy_i[:, 0]/hkm)*np.sin(np.pi*mxy_i[:, 1]/hkm)*np.sin(np.pi*mz_i/dkm))*1.0


rows = []; grids = {}
for hk in SCALES_KM:
    mt = checker3d(hk); dclean = G@mt
    corr_all = []; bydep = {z: [] for z in DEPTHS}; mrec_acc = np.zeros(M)
    for rr in range(8):
        rng = np.random.RandomState(rr+hk)
        dn = dclean + sig*rng.randn(len(sig))
        x = AtAinv @ (Aw.T @ (w*dn)); mrec = x[:M]; mrec_acc += mrec
        corr_all.append(np.corrcoef(mt, mrec)[0, 1])
        for z in DEPTHS:
            m = np.abs(mz_i-z) < 1
            if m.sum() > 3:
                bydep[z].append(np.corrcoef(mt[m], mrec[m])[0, 1])
    grids[hk] = mrec_acc/8
    rec = dict(scale_km=hk, overall=np.nanmean(corr_all))
    for z in DEPTHS:
        rec[f"z{int(z)}"] = np.nanmean(bydep[z]) if bydep[z] else np.nan
    rows.append(rec)
    print(f"{hk}km checker: overall r={rec['overall']:.2f}  by depth " + " ".join(f"{int(z)}km:{rec[f'z{int(z)}']:.2f}" for z in DEPTHS))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/volume3d_checkerboard_new.csv", index=False)

# figure: (top) recovery corr vs scale per depth layer; (bottom) recovered 250km checker across depth layers
fig = plt.figure(figsize=(4.6*len(DEPTHS)/2.2, 8))
axc = fig.add_axes([0.08, 0.60, 0.36, 0.34])
for z in DEPTHS:
    axc.plot(R.scale_km, R[f"z{int(z)}"], "o-", label=f"{int(z)} km")
axc.axhline(0.7, ls=":", color="k", lw=1); axc.set_xlabel("horizontal checker size (km)"); axc.set_ylabel("recovery corr")
axc.set_title("recovery by depth layer & scale\n(new 187-stn network, honest noise)"); axc.legend(fontsize=7, title="depth"); axc.grid(alpha=.3); axc.set_ylim(-.1, 1.05)
axt = fig.add_axes([0.52, 0.60, 0.44, 0.34]); axt.axis("off")
axt.text(0, 0.5, "CAVEAT: single 2-4 s window\n=> weak VERTICAL resolution.\nDepth layers trade off (smearing),\nso 'by-depth' curves largely track\neach other. Horizontal multi-scale\nrecovery is the informative axis:\nlarge recovers, fine doesn't.", fontsize=10, va="center")
# bottom row: recovered 250km checker at each depth layer
mrec = grids[250]
for k, z in enumerate(DEPTHS):
    a = fig.add_axes([0.04+k*0.16, 0.06, 0.14, 0.40])
    m = np.abs(mz_i-z) < 1
    sc = a.scatter(hlon_i[m], hlat_i[m], c=mrec[m], cmap="RdBu_r", vmin=-1, vmax=1, s=10, marker="s")
    a.set_title(f"{int(z)} km", fontsize=9); a.set_xticks([]); a.set_yticks([])
fig.text(0.5, 0.49, "Recovered 250 km checker across depth layers (2-4 s; note depth smearing)", ha="center", fontsize=10)
cax = fig.add_axes([0.965, 0.08, 0.012, 0.36]); plt.colorbar(sc, cax=cax, label="dbeta/beta (%)")
fig.savefig(f"{OUT}/volume3d_checkerboard_new.png", dpi=130, bbox_inches="tight")
print(f"wrote {OUT}/volume3d_checkerboard_new.png, .csv")
print(R.round(2).to_string(index=False))
