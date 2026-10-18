#!/usr/bin/env python
"""DELAYED responses + OVERLAP-aware lag profiles (user 2026-09-10: "beware of overlapping patches, how delayed responses").
For each post-2011 substantive patch: patch-mean corrected+deseasoned interface δβ/β from -150 to +300 d relative to the catalog
START, minus its PRE mean (-120..-31). Contamination mask: days on which ANOTHER substantive patch is active on this footprint
(or, stricter, anywhere within 200 km of the footprint) -> those days are dropped in the "clean" profile. Binned in 30-d lags.
Stacks by terrane class and by segment, for Mw>=6.2 and Mw>=6.0, all vs clean; null = rigid time shifts of the field (same shift
for every patch in a draw, 200 draws) -> z per lag bin and FWER over bins. Per patch: lag bin of the extreme |deviation| and its z.
Outputs: res_catalog_sse/patch_lag_profiles.csv (class x bin), patch_lag_perpatch.csv, figures/sse_patch_lag_profiles.png
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(5); NDRAW = 200
cat = pd.read_csv(f"{SD}/slip_patches_5mm.csv", parse_dates=["start", "end"]); zl = np.load(f"{SD}/slip_patches_5mm.npz", allow_pickle=True); L = zl["label"].T; days = pd.to_datetime(zl["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); cells = zc["cell"]; clat, clon, ip, valid = zc["cell_lat"], zc["cell_lon"], zc["patch"], zc["valid"]
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); Tc = pd.read_csv(f"{RD}/cell_T_patchmapped.csv").set_index("cell").T_cell.reindex([str(c) for c in cells]).values; Tc = np.where(np.isfinite(Tc) & (Tc > 0.2), Tc, 1.0)
Y = pd.DataFrame(np.asarray(zp["MODEL"], float) / Tc[None, :], index=pd.to_datetime(zp["days"])).reindex(days).values.copy()
active = (L[ip] > 0).T; month = days.month.values; Ym = np.where(active, np.nan, Y)
for m in range(1, 13): s = month == m; Y[s] -= np.nanmean(Ym[s], 0)[None]
sub = cat[~cat.minor & ~cat.broad_pulse]; nonminor = set(sub.id); terrane = lambda la: "Klamath/Franciscan (<43.2)" if la < 43.2 else ("Siletzia core (43.2-46.8)" if la < 46.8 else ("Olympic/Crescent (46.8-48.6)" if la < 48.6 else "Vancouver Is. (>48.6)"))
TERR = ["Klamath/Franciscan (<43.2)", "Siletzia core (43.2-46.8)", "Olympic/Crescent (46.8-48.6)", "Vancouver Is. (>48.6)"]; seg = lambda la: "north" if la >= 46 else ("central" if la >= 44 else "south")
REL = np.arange(-150, 301); BINS = list(range(-150, 301, 30)); nb = len(BINS) - 1; PRE = (REL >= -120) & (REL <= -31)
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
tri_ll = np.loadtxt("/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt", comments="#").reshape(-1, 3, 3)[:, :, :2].mean(1); TXY = km(tri_ll[:, 1], tri_ll[:, 0])
P = []
for _, r in sub[sub.start >= "2011-01-01"].iterrows():
    ts, te = int(gi[r.start]), int(gi[r.end]); fp = np.unique(np.where(L[:, ts:te + 1] == r.id)[0]); inside = np.where(valid & np.isin(ip, fp))[0]
    if len(inside) < 5: continue
    near = np.where(np.sqrt(((TXY[:, None, :] - TXY[fp][None, :, :]) ** 2).sum(-1)).min(1) <= 200)[0]
    idx = ts + REL; okt = (idx >= 0) & (idx < G); idxc = np.clip(idx, 0, G - 1)
    other_fp = np.array([(np.isin(L[fp, t], list(nonminor - {r.id}))).any() if okt[k] else True for k, t in enumerate(idxc)])
    other_near = np.array([(np.isin(L[near, t], list(nonminor - {r.id}))).any() if okt[k] else True for k, t in enumerate(idxc)])
    P.append(dict(id=r.id, ts=ts, te=te, dur=te - ts + 1, inside=inside, Mw=r.Mw, terrane=terrane(r.lat_centroid), seg=seg(r.lat_centroid), start=r.start, okt=okt, cont_fp=other_fp, cont_near=other_near))
print(f"{len(P)} patches; per-patch contamination (another substantive patch on footprint) covers a median {np.median([p['cont_fp'][p['okt']].mean() for p in P])*100:.0f} % of the -150..+300 d window; within 200 km: {np.median([p['cont_near'][p['okt']].mean() for p in P])*100:.0f} %")
def series(p, shift=0, mask=None):
    idx = (p["ts"] + REL + shift) % G; s = np.nanmean(Y[idx][:, p["inside"]], 1); s[~p["okt"]] = np.nan
    if mask is not None: s = np.where(mask, np.nan, s)
    b = np.nanmean(s[PRE]); return s - b
def binned(s): return np.array([np.nanmean(s[(REL >= BINS[k]) & (REL < BINS[k + 1])]) for k in range(nb)])
rows = []; prof = []; PROF = {}
for setname, mwmin in (("Mw>=6.2", 6.2), ("Mw>=6.0", 6.0)):
    Ps = [p for p in P if p["Mw"] >= mwmin]
    for clean in ("all", "clean_fp", "clean_near"):
        mk = lambda p: None if clean == "all" else (p["cont_fp"] if clean == "clean_fp" else p["cont_near"])
        S = np.array([binned(series(p, 0, mk(p))) for p in Ps]); groups = {"ALL": np.ones(len(Ps), bool)}
        groups.update({t: np.array([p["terrane"] == t for p in Ps]) for t in TERR}); groups.update({s: np.array([p["seg"] == s for p in Ps]) for s in ("north", "central", "south")})
        NUL = np.full((NDRAW, len(Ps), nb), np.nan)
        for d in range(NDRAW):
            sh = int(rng.randint(365, G - 365)); NUL[d] = np.array([binned(series(p, sh, mk(p))) for p in Ps])
        for gname, g in groups.items():
            if g.sum() < 3: continue
            obs = np.nanmean(S[g], 0); nul = np.nanmean(NUL[:, g, :], 1); z = (obs - np.nanmean(nul, 0)) / np.nanstd(nul, 0); zmax = np.nanmax(np.abs((nul - np.nanmean(nul, 0)) / np.nanstd(nul, 0)), 1)
            PROF[(setname, clean, gname)] = (obs, np.nanpercentile(nul, 2.5, 0), np.nanpercentile(nul, 97.5, 0), int(g.sum()))
            for k in range(nb): prof.append(dict(set=setname, mask=clean, group=gname, n=int(g.sum()), lag_bin_start=BINS[k], lag_bin_end=BINS[k + 1] - 1, mean=obs[k], null_sd=np.nanstd(nul[:, k]), z=z[k], p=float(np.nanmean(np.abs(nul[:, k] - np.nanmean(nul[:, k])) >= abs(obs[k] - np.nanmean(nul[:, k])))), p_fwer_bins=float(np.nanmean(zmax >= abs(z[k])))))
        if clean == "clean_fp":
            for i, p in enumerate(Ps):
                s = S[i]; nz = (s - np.nanmean(NUL[:, i, :], 0)) / np.nanstd(NUL[:, i, :], 0)
                if np.isfinite(nz).sum() < 3: continue
                k = int(np.nanargmax(np.abs(nz))); rows.append(dict(set=setname, id=p["id"], start=p["start"].date(), Mw=p["Mw"], terrane=p["terrane"], seg=p["seg"], dur_d=p["dur"], extreme_lag_bin_start=BINS[k], extreme_value=s[k], extreme_z=nz[k], bins_valid=int(np.isfinite(nz).sum()), frac_contaminated=float(p["cont_fp"][p["okt"]].mean())))
PR = pd.DataFrame(prof); PR.to_csv(f"{SD}/patch_lag_profiles.csv", index=False); PP = pd.DataFrame(rows); PP.to_csv(f"{SD}/patch_lag_perpatch.csv", index=False)
pd.set_option("display.width", 250)
print("\n=== lag bins with p_fwer_bins < 0.10 (any set/mask/group) ==="); h = PR[PR.p_fwer_bins < 0.10].sort_values("p_fwer_bins"); print(h[["set", "mask", "group", "n", "lag_bin_start", "mean", "null_sd", "z", "p", "p_fwer_bins"]].round(3).to_string(index=False) if len(h) else "none")
print("\n=== Mw>=6.2, clean_fp: mean profile by group (lag bin start -> mean %, z) ===")
for gname in ["ALL"] + TERR + ["north", "central", "south"]:
    s = PR[(PR.set == "Mw>=6.2") & (PR["mask"] == "clean_fp") & (PR.group == gname)]
    if len(s): print(f"{gname:30s} n={s.n.iloc[0]:3d} | " + " ".join(f"{int(b):+4d}:{m:+.2f}({z:+.1f})" for b, m, z in zip(s.lag_bin_start, s["mean"], s.z)))
print("\n=== per-patch extreme lag (Mw>=6.2, clean): distribution of the lag bin of the largest |z| ===")
pp = PP[PP.set == "Mw>=6.2"]; print(pp.groupby("extreme_lag_bin_start").size().to_string()); print(f"patches with |z|>=2.5 at the extreme: {(pp.extreme_z.abs() >= 2.5).sum()} of {len(pp)} (expected under null with {nb} bins ~ {len(pp) * nb * 0.012:.1f})")
print(pp[pp.extreme_z.abs() >= 2.5][["id", "start", "Mw", "terrane", "extreme_lag_bin_start", "extreme_value", "extreme_z", "frac_contaminated"]].round(2).to_string(index=False))
# ---------- figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = dict(zip(TERR, ["#ea580c", "#059669", "#1d4ed8", "#7c3aed"])); COL.update(north="#1d4ed8", central="#059669", south="#ea580c", ALL="k")
fig, ax = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True); xc = np.array(BINS[:-1]) + 15
for i, setname in enumerate(("Mw>=6.2", "Mw>=6.0")):
    for j, (clean, groups, ttl) in enumerate((("all", ["ALL", "north", "central", "south"], "all days, by segment"), ("clean_fp", ["ALL", "north", "central", "south"], "contaminated days masked (other patch on footprint), by segment"), ("clean_fp", TERR, "contaminated days masked, by terrane"))):
        a = ax[i, j]
        for g in groups:
            if (setname, clean, g) not in PROF: continue
            obs, lo, hi, n = PROF[(setname, clean, g)]; c = COL[g]; a.fill_between(xc, lo, hi, color=c, alpha=.07); a.plot(xc, obs, "o-", color=c, ms=4, lw=1.3, label=f"{g.split(' (')[0]} (n {n})")
        a.axvline(0, color="k", lw=.8); a.axhline(0, color="k", lw=.5); a.axvspan(-120, -31, color="#e5e7eb", alpha=.5); a.set_xlabel("days since patch start (30-d bins)"); a.set_ylabel("patch-mean dv/v rel. PRE (%)"); a.set_title(f"{setname}: {ttl}", fontsize=9); a.legend(fontsize=7)
fig.suptitle("Lag profiles of dv/v inside slip patches: delayed responses and overlap contamination (bands = shift-null 95 %; grey = PRE window)", fontsize=10)
fig.savefig("figures/sse_patch_lag_profiles.png", dpi=130); print("\nwrote figures/sse_patch_lag_profiles.png, patch_lag_profiles.csv, patch_lag_perpatch.csv")
