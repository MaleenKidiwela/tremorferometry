#!/usr/bin/env python
"""Per-patch MAPS of interface δβ/β inside each slip patch (user 2026-09-09): PRE = [start-60, start-1], DUR = [start, end],
POST = [end+31, end+60] (catalog start/end = first/last day the trailing-30-d slip window > 5 mm; dv/v is in the same trailing frame).
For every patch (post-2011, Mw >= MWMIN, not minor / broad pulse; unbridged catalog) and window: per-cell time-mean map of the
corrected interface δβ/β (deseasoned per cell unless DESEASON=0) for the 0.1° cells whose nearest mesh triangle is in the footprint;
patch mean, sd across cells, median, n cells; overlap flags: another component active on this footprint during PRE / POST
(footprint overlap with a non-minor component; minor blips counted separately) and the number of other non-minor components active anywhere in the margin during PRE / POST.
Outputs: res_catalog_sse/patch_dvv_maps_summary.csv; figures/sse_patch_dvv_maps_p{n}.png (6 patches x 3 windows per page).
Env: MWMIN (6.2), DESEASON (1), VLIM (1.5 %), CATSUF ("").
"""
import os, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; REPO = "/home/jovyan/cascadia-sse-strain"
MWMIN = float(os.environ.get("MWMIN", "6.2")); DES = os.environ.get("DESEASON", "1") == "1"; VLIM = float(os.environ.get("VLIM", "1.5")); SUF = os.environ.get("CATSUF", ""); TAG = "" if MWMIN == 6.2 else f"_Mw{MWMIN:g}"
cat = pd.read_csv(f"{SD}/slip_patches_5mm{SUF}.csv", parse_dates=["start", "end"]); zl = np.load(f"{SD}/slip_patches_5mm{SUF}.npz", allow_pickle=True); L = zl["label"].T; days = pd.to_datetime(zl["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); cells = zc["cell"]; clat, clon, ip, valid = zc["cell_lat"], zc["cell_lon"], zc["patch"], zc["valid"]
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); assert list(zp["cell"]) == list(cells)
Tc = pd.read_csv(f"{RD}/cell_T_patchmapped.csv").set_index("cell").T_cell.reindex([str(c) for c in cells]).values; Tc = np.where(np.isfinite(Tc) & (Tc > 0.2), Tc, 1.0)
Y = pd.DataFrame(np.asarray(zp["MODEL"], float) / Tc[None, :], index=pd.to_datetime(zp["days"])).reindex(days).values.copy()
if DES:
    active = (L[ip] > 0).T; month = days.month.values; Ym = np.where(active, np.nan, Y)
    for m in range(1, 13): s = month == m; Y[s] -= np.nanmean(Ym[s], 0)[None]
tri = np.loadtxt(f"{REPO}/data/cascadia_mesh.txt", comments="#").reshape(-1, 3, 3)[:, :, :2]
sel = cat[(cat.start >= "2011-01-01") & ~cat.minor & ~cat.broad_pulse & (cat.Mw >= MWMIN)].sort_values("start")
nonminor = set(cat[~cat.minor].id); seg = lambda la: "north" if la >= 46 else ("central" if la >= 44 else "south")
WIN = lambda ts, te: {"PRE (−60..−1 d)": np.arange(max(ts - 60, 0), ts), "DUR (start..end)": np.arange(ts, te + 1), "POST (+31..+60 d)": np.arange(te + 31, min(te + 61, G))}
rows = []; panels = []
for _, r in sel.iterrows():
    ts, te = int(gi[r.start]), int(gi[r.end]); fp = np.unique(np.where(L[:, ts:te + 1] == r.id)[0]); inside = np.where(valid & np.isin(ip, fp))[0]
    if len(inside) < 5: continue
    rec = dict(id=r.id, start=r.start.date(), end=r.end.date(), seg=seg(r.lat_centroid), Mw=r.Mw, n_cells=len(inside), n_triangles=len(fp)); maps = {}
    for w, idx in WIN(ts, te).items():
        blk = Y[idx][:, inside]; cm = np.nanmean(blk, 0); maps[w] = cm; k = w.split()[0]
        rec[f"{k}_mean"] = np.nanmean(cm); rec[f"{k}_sd_cells"] = np.nanstd(cm); rec[f"{k}_median"] = np.nanmedian(cm); rec[f"{k}_sd_celldays"] = np.nanstd(blk); rec[f"{k}_n_valid_cells"] = int(np.isfinite(cm).sum()); rec[f"{k}_ndays"] = len(idx)
        oth = L[fp][:, idx]; oth_ids = set(np.unique(oth[(oth > 0) & (oth != r.id)])) if len(idx) else set()
        anyw = set(np.unique(L[:, idx][L[:, idx] > 0])) - {r.id} if len(idx) else set()
        rec[f"{k}_overlap_footprint_ids"] = ";".join(str(int(x)) for x in sorted(oth_ids & nonminor)) if oth_ids & nonminor else ""; rec[f"{k}_overlap_footprint_nonminor"] = bool(oth_ids & nonminor); rec[f"{k}_overlap_footprint_minor_n"] = len(oth_ids - nonminor)
        rec[f"{k}_other_nonminor_active_anywhere"] = len(anyw & nonminor)
    rec["DUR_minus_PRE"] = rec["DUR_mean"] - rec["PRE_mean"]; rec["POST_minus_PRE"] = rec["POST_mean"] - rec["PRE_mean"]; rows.append(rec); panels.append((rec, inside, fp, maps))
S = pd.DataFrame(rows); S.to_csv(f"{SD}/patch_dvv_maps_summary{TAG}.csv", index=False)
pd.set_option("display.width", 250)
print(S[["id", "start", "end", "seg", "Mw", "n_cells", "PRE_mean", "PRE_sd_cells", "DUR_mean", "DUR_sd_cells", "POST_mean", "POST_sd_cells", "DUR_minus_PRE", "POST_minus_PRE", "PRE_overlap_footprint_nonminor", "POST_overlap_footprint_nonminor", "PRE_other_nonminor_active_anywhere", "POST_other_nonminor_active_anywhere"]].round(2).to_string(index=False))
print(f"\n{len(S)} patches | PRE footprint-overlap with another NON-MINOR patch: {S.PRE_overlap_footprint_nonminor.sum()} | POST: {S.POST_overlap_footprint_nonminor.sum()} | mean DUR−PRE north {S[S.seg=='north'].DUR_minus_PRE.mean():+.3f}, central {S[S.seg=='central'].DUR_minus_PRE.mean():+.3f}, south {S[S.seg=='south'].DUR_minus_PRE.mean():+.3f}")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import Rectangle; from matplotlib.collections import PolyCollection, PatchCollection
PER = 6; npg = int(np.ceil(len(panels) / PER)); outs = []
for pg in range(npg):
    sub = panels[pg * PER:(pg + 1) * PER]; fig, axs = plt.subplots(len(sub), 3, figsize=(15, 4.6 * len(sub)), constrained_layout=True); axs = np.atleast_2d(axs)
    for i, (rec, inside, fp, maps) in enumerate(sub):
        lo0, la0 = clon[inside].min() - 0.8, clat[inside].min() - 0.6; lo1, la1 = clon[inside].max() + 0.8, clat[inside].max() + 0.6
        box = np.where(valid & (clon >= lo0) & (clon <= lo1) & (clat >= la0) & (clat <= la1))[0]
        for j, (w, cm) in enumerate(maps.items()):
            a = axs[i, j]; k = w.split()[0]
            a.add_collection(PatchCollection([Rectangle((clon[c] - 0.05, clat[c] - 0.05), 0.1, 0.1) for c in box], facecolor="#eef2f7", edgecolor="none"))
            a.add_collection(PolyCollection(tri[fp], facecolor="none", edgecolor="#64748b", lw=0.7))
            ok = np.isfinite(cm); pc = PatchCollection([Rectangle((clon[c] - 0.05, clat[c] - 0.05), 0.1, 0.1) for c in inside[ok]], cmap="RdBu_r", edgecolor="k", lw=0.15); pc.set_array(cm[ok]); pc.set_clim(-VLIM, VLIM); a.add_collection(pc)
            a.set_xlim(lo0, lo1); a.set_ylim(la0, la1); a.set_aspect(1 / np.cos(np.radians(46))); a.tick_params(labelsize=7); a.grid(alpha=.2)
            flag = "no other patch on footprint" if k in ("PRE", "POST") else ""
            if k in ("PRE", "POST"):
                if rec[f"{k}_overlap_footprint_nonminor"]: ids_ = rec[f"{k}_overlap_footprint_ids"].split(";"); flag = f"⚠ overlap: patch #{', #'.join(ids_[:3])}{'…' if len(ids_) > 3 else ''}"
                elif rec[f"{k}_overlap_footprint_minor_n"]: flag = f"{rec[f'{k}_overlap_footprint_minor_n']} minor blip(s) on footprint"
                if rec[f"{k}_other_nonminor_active_anywhere"]: flag += f"; {rec[f'{k}_other_nonminor_active_anywhere']} other active elsewhere"
            a.set_title(f"{w}: mean {rec[f'{k}_mean']:+.2f} %, sd {rec[f'{k}_sd_cells']:.2f} (n {rec[f'{k}_n_valid_cells']}, {rec[f'{k}_ndays']} d)\n{flag}", fontsize=7.5)
            if j == 0: a.set_ylabel(f"#{rec['id']}  {rec['start']} → {rec['end']}\n{rec['seg']}, Mw {rec['Mw']:.2f}, {rec['n_cells']} cells", fontsize=8)
    fig.colorbar(pc, ax=axs.ravel().tolist(), location="bottom", shrink=0.4, pad=0.01, aspect=40, label=f"interface δβ/β, per-cell window mean (%, corrected{', deseasoned' if DES else ''})")
    fig.suptitle(f"δβ/β inside each slip patch — PRE [start−60, start−1] · DUR [start, end] · POST [end+31, end+60]  (Mw ≥ {MWMIN:g}, post-2011, unbridged; page {pg + 1}/{npg}; grey outline = model footprint triangles)", fontsize=10)
    out = f"figures/sse_patch_dvv_maps{TAG}_p{pg + 1}.png"; fig.savefig(out, dpi=130); plt.close(fig); outs.append(out)
print("wrote", ", ".join(outs), f"and {SD}/patch_dvv_maps_summary{TAG}.csv")
