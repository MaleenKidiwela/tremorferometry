#!/usr/bin/env python
"""Exploration round 2 (user: "explore more"), shared rigid-shift nulls + max-statistic family-wise control.
A. VARIABLES x PATCHES x LAGS: whole-record cross-correlation of patch-mean 30-d change of dv/v (5 products: interface
   deseasoned, family TOTAL 2-4 s, family site-removed, SITES, family TOTAL 3-5 s) with the patch-mean 30-d change of 12
   model variables (slip-in-window, cum slip, 6 strain, dilatation, Coulomb/shear/normal at 10 km, trailing frame),
   lags -120..+120 d step 10; z vs 400 shifts; FWER over (variables x lags) per patch-product and over everything.
B. STATION-LEVEL STRAIN RELATION: model strain at each station's own location (1 km depth) from the repo Green's
   functions; near-station dv/v (SITES) 30-d change vs station dilatation / max-shear / e_zz 30-d change; per-station
   correlation with shift null, pooled slope, sign map, and relation to distance from the deep patch.
C. STATION ETS-RESPONSE MAP: per station the mean episode contrast (ON[+15,+60]-PRE[-90,-30], admitted episodes whose
   footprint is within 150 km) of its certified families in 2-4 s, 3-5 s and MIRROR windows; sign/size vs instrument
   type, distance to the slip footprint, up-dip/down-dip position; consistency 2-4 vs 3-5 vs mirror.
Outputs: res_catalog_sse/explore2_*.csv, figures/sse_explore2.png
"""
import numpy as np, pandas as pd, os, sys, glob
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; REPO = "/home/jovyan/cascadia-sse-strain"
sys.path.insert(0, f"{REPO}/src"); from sse_strain import strain as st
rng = np.random.RandomState(31); ND = int(os.environ.get("ND", "400"))
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; gi = pd.Series(np.arange(G), index=days)
win = np.asarray((days >= "2011-01-01") & (days <= "2026-07-20"))
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True); PAT = np.asarray(zt["patterns"], float); plat = zt["patch_lat"]; band = zt["patch_band"]
deepq = [q for q in np.argsort(plat) if band[q] == "deep" and 40.5 < plat[q] < 49.5]
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); sta_tags = list(zp["sta_tags"])
Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=md, columns=sta_tags).reindex(days)
# family products at cells: 2-4 s TOTAL / site-removed (parquet) and 3-5 s TOTAL (CSV, certified families)
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell", "date", "dvv"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv; fd["dw"] = fd.v - fd.groupby(["tag", "date"]).v.transform("mean")
cd = fd.groupby(["ci", "date"]).agg(tot=("v", "mean"), dw=("dw", "mean")).reset_index(); okd = cd.date.isin(days)
def grid_from(cd, col):
    A = np.full((G, M), np.nan); A[gi[cd.date[okd]].values, cd.ci[okd].values] = cd[col][okd].values; return A
TOT, DW = grid_from(cd, "tot"), grid_from(cd, "dw")
famcell = fd[["tag", "fam", "ci"]].drop_duplicates(); sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); famset = {t: set(g.fam) for t, g in fd.groupby("tag")}
del fd, cd
rows = []
for f in glob.glob("data/daily_dvv_*_Z_3to5.csv"):
    tag = os.path.basename(f).split("daily_dvv_")[1].split("_Z_3to5")[0]
    if tag not in famset: continue
    x = pd.read_csv(f, usecols=["patch", "date", "dvv", "cc_max"]); x = x[(x.cc_max > 0.7) & (x.patch.isin(famset[tag]))]
    x["tag"] = tag; x["date"] = pd.to_datetime(x.date); x["v"] = -x.dvv * 100; rows.append(x[["tag", "patch", "date", "v"]].rename(columns={"patch": "fam"}))
f35 = pd.concat(rows).merge(famcell, on=["tag", "fam"]); c35 = f35.groupby(["ci", "date"]).v.mean().reset_index(); ok35 = c35.date.isin(days)
TOT35 = np.full((G, M), np.nan); TOT35[gi[c35.date[ok35]].values, c35.ci[ok35].values] = c35.v[ok35].values
print(f"3-5 s family product: {f35.tag.nunique()} stations, {len(f35):,} rows", flush=True)
def pm(X, q):
    w = PAT[q]; fin = np.isfinite(X); num = np.nansum(np.where(fin, X, 0) * w[None], 1); den = (fin * w[None]).sum(1); o = num / np.where(den > 0, den, np.nan)
    o[(fin * w[None]).sum(1) / w.sum() < 0.5] = np.nan; return o
def deseason(y):
    s = pd.Series(y, index=days); return (s - s.groupby(s.index.month).transform("mean")).values
def d30(y): o = np.full_like(y, np.nan); o[30:] = y[30:] - y[:-30]; return o
VARS = {"slip in window": zc["slip_thrust_inc30"].astype(float) * 1e3, "cum slip": zc["slip_thrust_cum"].astype(float) * 1e3,
        "dilatation": zc["dilatation"].astype(float) * 1e8, "Coulomb": zc["coulomb"].astype(float) / 1e3, "shear stress": zc["shear_stress"].astype(float) / 1e3,
        "normal stress": zc["normal_stress"].astype(float) / 1e3}
for c in ("exx", "eyy", "ezz", "exy", "exz", "eyz"): VARS[f"strain {c}"] = zc[f"strain_{c}"].astype(float) * 1e8
LAGL = np.arange(-120, 121, 10)
def xcorr_lags(a, b):
    out = []
    for L in LAGL:
        bb = np.roll(b, L); m = np.isfinite(a) & np.isfinite(bb) & win; out.append(np.corrcoef(a[m], bb[m])[0, 1] if m.sum() > 500 else np.nan)
    return np.array(out)
PRODS = {"interface": Y, "family 2-4s": TOT, "site-removed 2-4s": DW, "family 3-5s": TOT35}
resA = []; Zall = {}
for q in deepq:
    cs = np.where(PAT[q] > 0.3)[0]; sr = sta_rows[sta_rows.ci.isin(cs)].groupby("tag").n.sum(); sr = sr / sr.sum(); stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    resp = {k: d30(deseason(pm(X, q))) for k, X in PRODS.items()}; resp["SITES"] = d30(deseason(SITES[stas].mean(axis=1).values)) if stas else np.full(G, np.nan)
    pred = {k: pm(X, q) for k, X in VARS.items()}; pred = {k: (v if k == "slip in window" else d30(v)) for k, v in pred.items()}
    obs = {(p, v): xcorr_lags(resp[p], pred[v]) for p in resp for v in pred}
    null = np.zeros((ND, len(resp), len(pred), len(LAGL)))
    for d in range(ND):
        s = rng.randint(365, G - 365)
        for i, p in enumerate(resp):
            for j, v in enumerate(pred): null[d, i, j] = xcorr_lags(resp[p], np.roll(pred[v], s))
    sd = np.nanstd(null, 0); zmax_null = np.nanmax(np.abs(null / np.where(sd > 0, sd, np.nan)), axis=(2, 3))     # (ND, prods)
    for i, p in enumerate(resp):
        for j, v in enumerate(pred):
            z = obs[(p, v)] / np.where(sd[i, j] > 0, sd[i, j], np.nan); k = int(np.nanargmax(np.abs(z))) if np.isfinite(z).any() else 0
            resA.append(dict(patch=round(plat[q], 1), product=p, variable=v, best_lag=int(LAGL[k]), r=obs[(p, v)][k], z=z[k],
                             p_fwer_lags=float((zmax_null[:, i] >= abs(z[k])).mean()), r_lag0=obs[(p, v)][LAGL == 0][0]))
    Zall[q] = zmax_null
    print(f"patch {plat[q]:.1f}N done", flush=True)
A = pd.DataFrame(resA); A.to_csv(f"{SD}/explore2_A_varmap.csv", index=False)
zglob = np.concatenate([Zall[q].max(1) for q in deepq])                       # max over prods per draw per patch -> approx global
A["p_fwer_global"] = [float((np.max([Zall[qq].max(1) for qq in deepq], axis=0) >= abs(z)).mean()) for z in A.z]
print("\nA. TOP variable-patch-product cells (z at best lag; p_fwer over lags x variables within patch-product; global over all):")
print(A.sort_values("p_fwer_global").head(15).round(3).to_string(index=False))
print(f"   cells with p_fwer_global<0.05: {(A.p_fwer_global < 0.05).sum()} of {len(A)}")

# ---- B. station-level strain relation ----
from sse_strain import ElasticModel, load_fault_ascii, load_solution
from sse_strain.green import strain_from_components
from sse_strain.geodesy import llh2localxy
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns(); slip_k = fault.slip_enu(Ls, Ld)
pairs = pd.read_csv(f"{RD}/pairs.csv"); sll = pairs.groupby("tag")[["sta_lat", "sta_lon"]].first(); sll = sll[sll.index.isin(SITES.columns)]
sx, sy = llh2localxy(sll.sta_lat.values, sll.sta_lon.values, fault.origin)
Es = strain_from_components(np.column_stack([sx * 1e3, sy * 1e3, -1000.0 * np.ones(len(sll))]), fault.vertices_m(), slip_k, sol.nu, progress=False)
A_trail = np.nan_to_num(zc["A_trail"].astype(float)); eps_s = np.einsum("obk,tk->tob", Es, A_trail)                     # (G, nsta, 6)
dil_s = eps_s[..., :3].sum(-1) * 1e8; ezz_s = eps_s[..., 2] * 1e8
ev = np.linalg.eigvalsh(st.voigt_to_tensor(eps_s.reshape(-1, 6))).reshape(G, len(sll), 3); mxs_s = (ev[..., 2] - ev[..., 0]) / 2 * 1e8
resB = []; pooled = []
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111])
deep_xy = km(clat[np.isfinite(dep) & (dep >= 30) & (dep <= 45)], clon[np.isfinite(dep) & (dep >= 30) & (dep <= 45)])
for j, tag in enumerate(sll.index):
    y = d30(deseason(SITES[tag].values)); preds = {"dilatation": d30(dil_s[:, j]), "max shear": d30(mxs_s[:, j]), "e_zz": d30(ezz_s[:, j])}
    m = np.isfinite(y) & win
    if m.sum() < 1000: continue
    dist = np.sqrt(((km(np.array([sll.sta_lat.iloc[j]]), np.array([sll.sta_lon.iloc[j]])) - deep_xy) ** 2).sum(1)).min()
    rec = dict(tag=tag, n=int(m.sum()), dist_deep_km=float(dist), borehole=tag.startswith("B") and "p90f40" in tag, amp_dil=float(np.nanstd(preds["dilatation"][win])))
    for k, x in preds.items():
        mm = m & np.isfinite(x); r = np.corrcoef(y[mm], x[mm])[0, 1]; nl = [np.corrcoef(y[mm], np.roll(x, rng.randint(365, G - 365))[mm])[0, 1] for _ in range(200)]
        rec[f"r_{k}"] = r; rec[f"p_{k}"] = float(np.mean(np.abs(nl) >= abs(r))); rec[f"z_{k}"] = r / np.std(nl)
        if k == "dilatation": pooled.append(np.column_stack([y[mm], x[mm], np.full(mm.sum(), j)]))
    resB.append(rec)
B = pd.DataFrame(resB); B.to_csv(f"{SD}/explore2_B_station_strain.csv", index=False)
P = np.vstack(pooled); slope = np.polyfit(P[:, 1], P[:, 0], 1)[0]
nl = [np.polyfit(np.roll(P[:, 1], rng.randint(1000, len(P) - 1000)), P[:, 0], 1)[0] for _ in range(200)]
print(f"\nB. STATION dv/v vs modeled strain at the station ({len(B)} stations): per-station r(dilatation) median {B.r_dilatation.median():+.3f}, "
      f"stations with p<0.01: {(B.p_dilatation < 0.01).sum()} (expected {0.01*len(B):.1f}); positive r fraction {(B.r_dilatation > 0).mean():.2f} | "
      f"max shear p<0.01: {(B['p_max shear'] < 0.01).sum()} | e_zz p<0.01: {(B.p_e_zz < 0.01).sum()}")
print(f"   pooled slope d(dv/v %)/d(dilatation 1e-8) = {slope:+.4f} (null sd {np.std(nl):.4f}, z {slope/np.std(nl):+.1f}); literature CWI sensitivity ~ -1e3..-1e4 per strain => {-1e3*1e-8*100:+.4f}..{-1e4*1e-8*100:+.4f} % per 1e-8")
print(B.sort_values("p_dilatation").head(8)[["tag", "dist_deep_km", "borehole", "r_dilatation", "z_dilatation", "p_dilatation", "r_max shear", "r_e_zz"]].round(3).to_string(index=False))

# ---- C. station ETS-response map (coda 2-4, coda 3-5, mirror) ----
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C", "S"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
XY = km(clat, clon)
def fam_contrasts(df, ions):
    piv = df.pivot_table(index="date", columns="fam", values="v").reindex(days); piv = piv - piv.groupby(piv.index.month).transform("mean"); out = []
    for fam in piv.columns:
        y = piv[fam].values.astype(float); vals = []
        for i in ions:
            a = y[max(i + 15, 0):i + 61]; b = y[max(i - 90, 0):i - 29]
            if np.isfinite(a).mean() > .6 and np.isfinite(b).mean() > .6: vals.append(np.nanmean(a) - np.nanmean(b))
        if len(vals) >= 4: out.append(np.mean(vals))
    return np.array(out)
resC = []
for tag in sorted(famset):
    if tag not in sll.index: continue
    sxy = km(np.array([sll.sta_lat[tag]]), np.array([sll.sta_lon[tag]]))
    ions = [int(gi[pd.Timestamp(r.fp_onset)]) for _, r in okep.iterrows() if np.sqrt(((XY[zf[f"fp_{r.key}"]] - sxy) ** 2).sum(1)).min() <= 150]
    if len(ions) < 4: continue
    rec = dict(tag=tag, n_ep=len(ions), lat=float(sll.sta_lat[tag]), lon=float(sll.sta_lon[tag]), borehole=tag.startswith("B") and "p90f40" in tag)
    for kind, fn in (("c24", f"data/daily_dvv_{tag}_Z_2to4.csv"), ("c35", f"data/daily_dvv_{tag}_Z_3to5.csv"), ("mir", f"data/daily_dvv_{tag}_MIRROR_2to4.csv")):
        if not os.path.exists(fn): rec[kind] = np.nan; rec[f"n_{kind}"] = 0; continue
        x = pd.read_csv(fn, usecols=["patch", "date", "dvv", "cc_max"]); x = x[(x.cc_max > 0.7) & (x.patch.isin(famset[tag]))].rename(columns={"patch": "fam"})
        x["date"] = pd.to_datetime(x.date); x["v"] = -x.dvv * 100; c = fam_contrasts(x, ions); rec[kind] = float(np.mean(c)) if len(c) else np.nan; rec[f"n_{kind}"] = len(c); rec[f"fpos_{kind}"] = float((c > 0).mean()) if len(c) else np.nan
    resC.append(rec)
C = pd.DataFrame(resC); C.to_csv(f"{SD}/explore2_C_station_response.csv", index=False)
cc = C.dropna(subset=["c24"])
print(f"\nC. STATION ETS-RESPONSE (2-4 s coda contrast, {len(cc)} stations): mean {cc.c24.mean():+.3f}%, positive {100*(cc.c24>0).mean():.0f}%; "
      f"|c24|>0.1%: {(cc.c24.abs()>0.1).sum()} stations; boreholes mean {cc[cc.borehole].c24.mean():+.3f} vs broadband {cc[~cc.borehole].c24.mean():+.3f}")
both = C.dropna(subset=["c24", "c35"]); print(f"   2-4 vs 3-5 s across stations: Spearman {pd.Series(both.c24).corr(both.c35, method='spearman'):+.2f}, sign agreement {100*np.mean(np.sign(both.c24)==np.sign(both.c35)):.0f}% (n={len(both)})")
bm = C.dropna(subset=["c24", "mir"]); print(f"   coda vs MIRROR across stations (n={len(bm)}): Spearman {pd.Series(bm.c24).corr(bm.mir, method='spearman'):+.2f}; mirror mean {bm.mir.mean():+.3f}")
print(f"   latitude dependence of c24: Spearman(lat, c24) {pd.Series(cc.lat).corr(cc.c24, method='spearman'):+.2f}")
print(cc.sort_values("c24").iloc[np.r_[0:5, -5:0]][["tag", "lat", "n_ep", "c24", "fpos_c24", "c35", "mir"]].round(3).to_string(index=False))

# ---- figure ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig = plt.figure(figsize=(18, 11)); gs = fig.add_gridspec(2, 3, width_ratios=[1.6, 1, 1])
ax = fig.add_subplot(gs[0, :2]); pv = A[A["product"] == "interface"].pivot(index="variable", columns="patch", values="z")
im = ax.imshow(pv.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto"); ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels(pv.columns); ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels(pv.index, fontsize=8)
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]):
        r = A[(A["product"] == "interface") & (A.variable == pv.index[i]) & (A.patch == pv.columns[j])].iloc[0]
        ax.text(j, i, f"{r.best_lag:+d}", ha="center", va="center", fontsize=6, color="k" if abs(r.z) < 2.5 else "w")
plt.colorbar(im, ax=ax, shrink=.7, label="z at best lag (label = lag, d)"); ax.set_title("A. interface δβ/β 30-d change vs model variables (patch-mean, best lag −120..+120 d); FWER threshold ≈ |z| 3.3", fontsize=9)
ax = fig.add_subplot(gs[0, 2]); pv2 = A[A["product"] == "family 3-5s"].pivot(index="variable", columns="patch", values="z")
im = ax.imshow(pv2.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto"); ax.set_xticks(range(pv2.shape[1])); ax.set_xticklabels(pv2.columns, fontsize=7); ax.set_yticks(range(pv2.shape[0])); ax.set_yticklabels(pv2.index, fontsize=7); ax.set_title("same, family 3-5 s coda", fontsize=9)
ax = fig.add_subplot(gs[1, 0]); sc = ax.scatter(C.lon, C.lat, c=C.c24, cmap="RdBu_r", vmin=-0.3, vmax=0.3, s=28, edgecolor="k", linewidth=.3); ax.set_aspect(1 / np.cos(np.radians(47)))
plt.colorbar(sc, ax=ax, shrink=.7, label="station ETS-locked coda (2-4 s) contrast (%)"); ax.set_title("C. per-station response to admitted slip episodes within 150 km", fontsize=9)
ax = fig.add_subplot(gs[1, 1]); ax.scatter(both.c24, both.c35, s=14, c="#1d4ed8"); ax.axhline(0, color="k", lw=.5); ax.axvline(0, color="k", lw=.5); ax.plot([-.4, .4], [-.4, .4], "k--", lw=.5)
ax.set_xlabel("2-4 s contrast (%)"); ax.set_ylabel("3-5 s contrast (%)"); ax.set_title("C. stations: 2-4 vs 3-5 s (velocity change would fall on the diagonal)", fontsize=9)
ax = fig.add_subplot(gs[1, 2]); ax.scatter(B.dist_deep_km, B.z_dilatation, c=["#dc2626" if b else "#1d4ed8" for b in B.borehole], s=14); ax.axhline(0, color="k", lw=.5); ax.axhline(2.6, color="r", ls="--", lw=.6); ax.axhline(-2.6, color="r", ls="--", lw=.6)
ax.set_xlabel("station distance to deep interface band (km)"); ax.set_ylabel("z of r(near-station dv/v, modeled dilatation at station)"); ax.set_title("B. station-level strain relation (red = borehole)", fontsize=9)
fig.suptitle("Exploration round 2: model variables vs coda dv/v across patches, lags, windows, and stations", fontsize=11)
fig.savefig("figures/sse_explore2.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_explore2.png + explore2_*.csv")
