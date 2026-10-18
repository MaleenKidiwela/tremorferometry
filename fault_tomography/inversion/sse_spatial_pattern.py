#!/usr/bin/env python
"""SPATIAL-PATTERN test (2026-09-04): per episode, the MAP of delta-beta/beta change (ON[+15,+60]-PRE[-90,-30] per cell,
all cells within 300 km of the footprint) is correlated across cells with the MAPS of modeled change of each stress/
strain variable (raw-amplitude ON-PRE through the Green's functions at 0.5/10/20 km, fixed + local receiver), the slip
map, and a pure GEOMETRY template (signed across-dip position relative to the slip centroid: up-dip +, down-dip -).
Pooled over episodes (model catalog N+C+S, and tremor catalog), null = rigid shift of the episode times (same maps,
same cells), max-statistic over variables.  Products: interface (deseasoned), family 2-4 s, family 3-5 s.
Outputs: res_catalog_sse/spatial_pattern.csv, figures/sse_spatial_pattern.png
"""
import numpy as np, pandas as pd, sys
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(123)
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src"); from sse_strain import strain as st
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; gi = pd.Series(np.arange(G), index=days); A_raw = np.nan_to_num(zc["A_raw"].astype(float)); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
zg = np.load(f"{SD}/strain_gf_offsets.npz", allow_pickle=True); EOFF = zg["E"].astype(float); MU = zg["mu"]; LAM = zg["lam"]; strike = zg["strike"]; dip = zg["dip"]; offs = zg["offsets"]
# products
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
F24 = np.load(f"{SD}/product_family24.npz", allow_pickle=True)["MODEL"].astype(float); F35 = np.load(f"{SD}/product_family35.npz", allow_pickle=True)["MODEL"].astype(float)
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C", "S"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
EPM = [dict(key=r.key, seg=r.seg, ion=int(gi[pd.Timestamp(r.fp_onset)]), ipeak=int(r.ipeak), fp=zf[f"fp_{r.key}"]) for _, r in okep.iterrows()]
tre = pd.read_csv(f"{SD}/tremor_episodes.csv"); EPT = []
for _, r in tre.iterrows():
    i = int(gi[pd.Timestamp(r.peak)]); on = int(gi[pd.Timestamp(r.onset)]); fp = np.where(np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (np.abs(clat - (r.band + 0.25)) <= 0.5))[0]
    if len(fp) >= 5: EPT.append(dict(key=r.key, seg=r.seg, ion=on, ipeak=i, fp=fp))
month = days.month.values; mask = np.zeros(G, bool)
for e in EPM: mask[max(e["ipeak"] - 45, 0):e["ipeak"] + 91] = True
def deseason(X):
    Xd = X.copy()
    for m in range(1, 13): sel = (month == m) & ~mask; Xd[month == m] -= np.nanmean(X[sel], 0)[None]
    return Xd
PROD = {"interface": deseason(Y), "family 2-4 s": deseason(F24), "family 3-5 s": deseason(F35)}
def model_maps(e, shift=0):
    i = wrap(e["ion"] + shift) if shift else e["ion"]; dA = A_raw[max(i + 15, 0):i + 61].mean(0) - A_raw[max(i - 90, 0):i - 29].mean(0); out = {}
    for k, off in enumerate(offs):
        eps_ = np.einsum("obk,k->ob", EOFF[k], dA)[None]; sig = st.hooke(eps_, MU[k][None], LAM[k][None])
        c, s_, n_ = st.coulomb_stress(sig, 0.0, 12.0, 90.0, friction=0.4)
        out[f"dilatation@{off:g}"] = eps_[0, :, :3].sum(1) * 1e8; out[f"Coulomb@{off:g}"] = c[0] / 1e3; out[f"shear@{off:g}"] = s_[0] / 1e3; out[f"normal@{off:g}"] = n_[0] / 1e3
        for b, cn in enumerate(("exx", "eyy", "ezz", "exy", "exz", "eyz")): out[f"{cn}@{off:g}"] = eps_[0, :, b] * 1e8
    out["slip map"] = inc[wrap(e["ipeak"] + shift) if shift else e["ipeak"]]
    return out
def geom_template(e):
    fp = e["fp"]; c = XY[fp].mean(0); v = XY - c; d = np.sqrt((v ** 2).sum(1)); return -v[:, 0] * 1.0   # west (up-dip) positive
def dvv_map(X, e, shift=0):
    i = wrap(e["ion"] + shift) if shift else e["ion"]; a = np.nanmean(X[max(i + 15, 0):i + 61], 0); b = np.nanmean(X[max(i - 90, 0):i - 29], 0); return a - b
def region(e): d = np.sqrt(((XY[:, None, :] - XY[e["fp"]][None, :, :]) ** 2).sum(-1)).min(1); return np.where(np.isfinite(dep) & (d <= 300))[0]
def scorr(a, b, cs):
    m = np.isfinite(a[cs]) & np.isfinite(b[cs]); return np.corrcoef(a[cs][m], b[cs][m])[0, 1] if m.sum() >= 30 and np.std(b[cs][m]) > 0 else np.nan
VARS = None; rows = []
for catname, EPS in (("model catalog", EPM), ("tremor catalog", EPT)):
    for pname, X in PROD.items():
        obs = {}; nul = {}
        for e in EPS:
            cs = region(e); mm = model_maps(e); mm["geometry (up-dip +)"] = geom_template(e); dm = dvv_map(X, e)
            if VARS is None: VARS = list(mm)
            for v in VARS: obs.setdefault(v, []).append(scorr(dm, mm[v], cs))
        obs = {v: np.array(x) for v, x in obs.items()}
        for d in range(200):
            sh = {s: int(rng.randint(90, span - 90)) for s in ("N", "C", "S", "SC")}; tmp = {v: [] for v in VARS}
            for e in EPS:
                cs = region(e); mm = model_maps(e, sh[e["seg"]]); mm["geometry (up-dip +)"] = geom_template(e); dm = dvv_map(X, e, sh[e["seg"]])   # shift BOTH -> tests whether dv/v-model map alignment exceeds what any random epoch gives
                for v in VARS: tmp[v].append(scorr(dm, mm[v], cs))
            for v in VARS: nul.setdefault(v, []).append(np.nanmean(tmp[v]))
        nulA = np.array([nul[v] for v in VARS]); zmax = np.nanmax(np.abs((nulA - nulA.mean(1, keepdims=True)) / nulA.std(1, keepdims=True)), 0)
        for v in VARS:
            o = np.nanmean(obs[v]); n_ = np.array(nul[v]); z = (o - n_.mean()) / n_.std(); p = float(np.mean(np.abs(n_ - n_.mean()) >= abs(o - n_.mean()))); pf = float(np.mean(zmax >= abs(z)))
            rows.append(dict(catalog=catname, product=pname, variable=v, n_ep=int(np.isfinite(obs[v]).sum()), mean_spatial_r=o, frac_pos=float(np.nanmean(obs[v] > 0)), null_mean=n_.mean(), null_sd=n_.std(), z=z, p=p, p_fwer=pf))
        top = sorted([r for r in rows if r["catalog"] == catname and r["product"] == pname], key=lambda r: r["p"])[:4]
        print(f"{catname} | {pname}: " + " ; ".join(f"{r['variable']} r={r['mean_spatial_r']:+.3f} z={r['z']:+.1f} p={r['p']:.3f} pF={r['p_fwer']:.2f}" for r in top), flush=True)
T = pd.DataFrame(rows); T.to_csv(f"{SD}/spatial_pattern.csv", index=False)
print("\nvariables with p_fwer<0.05:"); print(T[T.p_fwer < 0.05].round(3).to_string(index=False) if (T.p_fwer < 0.05).any() else "none")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(17, 7))
for k, catname in enumerate(("model catalog", "tremor catalog")):
    pv = T[T.catalog == catname].pivot(index="variable", columns="product", values="z").reindex(VARS)
    im = ax[k].imshow(pv.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto"); ax[k].set_xticks(range(pv.shape[1])); ax[k].set_xticklabels(pv.columns, fontsize=8); ax[k].set_yticks(range(len(VARS))); ax[k].set_yticklabels(VARS, fontsize=6); ax[k].set_title(f"{catname}: z of mean spatial correlation (dv/v change map vs model change map)", fontsize=9)
plt.colorbar(im, ax=ax, shrink=.6, label="z vs shifted-episode null"); fig.savefig("figures/sse_spatial_pattern.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_spatial_pattern.png, spatial_pattern.csv")
