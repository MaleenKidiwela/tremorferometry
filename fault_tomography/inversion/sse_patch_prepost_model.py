#!/usr/bin/env python
"""Companion to sse_patch_prepost.py (user 2026-09-09): how every SSE-model-derived variable changes INSIDE each slip patch
(same patches, same inside/ring cells, same PRE/DUR/POST/POST2 windows in the trailing-30-d frame), and whether the dv/v change
inside the patch scales with any variable's change across patches.
Fields (components_at_cells_trailing30.npz, 10 km above the interface, trailing frame): strain_e** (nanostrain), stress_s** (kPa),
dilatation / eps_eq / max_shear_strain (nanostrain), mean_stress / pressure / von_mises / max_shear / coulomb / shear / normal
(kPa), slip_thrust_cum (mm), slip_thrust_inc30 (mm per 30 d).
Per patch x field: window means over inside cell-days; changes DUR-PRE, POST-PRE, POST2-PRE; ring equivalents.
Cross-patch test: Spearman rho between the dv/v change (d_median_dur, d_median_post, d_median_dur_minus_ring; interface and
family 2-4 s) and each field change, per set; null = the dv/v shift-null draws saved by sse_patch_prepost.py correlated with the
same field changes (max-|rho| over fields -> FWER).
Outputs: res_catalog_sse/patch_prepost_model.csv (long), patch_prepost_model_summary.csv (typical change per field per set),
patch_prepost_model_corr.csv, figures/sse_patch_prepost_model.png
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.stats import spearmanr
SD = "fault_tomography/inversion/res_catalog_sse"
zn = np.load(f"{SD}/patch_prepost_null.npz", allow_pickle=True); pids = list(zn["pids"]); INS = list(zn["inside"]); RNG = list(zn["ring"]); TS = zn["ts"]; TE = zn["te"]; PRODS = list(zn["products"])
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days)
SCALE = {"strain": 1e9, "dilatation": 1e9, "eps_eq": 1e9, "max_shear_strain": 1e9, "stress": 1e-3, "mean_stress": 1e-3, "pressure": 1e-3, "von_mises_stress": 1e-3, "max_shear_stress": 1e-3, "coulomb": 1e-3, "shear_stress": 1e-3, "normal_stress": 1e-3, "slip_thrust_cum": 1e3, "slip_thrust_inc30": 1e3}
UNIT = {1e9: "nstrain", 1e-3: "kPa", 1e3: "mm"}
FIELDS = [k for k in zc.files if k.startswith(("strain_", "stress_")) or k in ("dilatation", "eps_eq", "max_shear_strain", "mean_stress", "pressure", "von_mises_stress", "max_shear_stress", "coulomb", "shear_stress", "normal_stress", "slip_thrust_cum", "slip_thrust_inc30")]
def scale(k): return SCALE.get(k, SCALE.get(k.split("_")[0], 1.0))
WIN = {"PRE": (-120, -31, "s"), "DUR": (0, 0, "se"), "POST": (31, 120, "e"), "POST2": (121, 240, "e")}
def widx(ts, te, w):
    a, b, ref = WIN[w]
    if ref == "se": return np.arange(ts, te + 1)
    r = ts if ref == "s" else te; return np.arange(max(r + a, 0), min(r + b, G - 1) + 1)
T = pd.read_csv(f"{SD}/patch_prepost_patches.csv"); cat = pd.read_csv(f"{SD}/slip_patches_5mm.csv")[["id", "Mw", "lat_centroid", "start"]].set_index("id")
rows = []
for k in FIELDS:
    X = zc[k].astype(float) * scale(k)
    for j, pid in enumerate(pids):
        rec = dict(id=pid, field=k, unit=UNIT[scale(k)])
        for w in WIN:
            idx = widx(TS[j], TE[j], w); rec[w] = np.nanmean(X[idx][:, INS[j]]); rec[w + "_ring"] = np.nanmean(X[idx][:, RNG[j]])
        rec["d_dur"] = rec["DUR"] - rec["PRE"]; rec["d_post"] = rec["POST"] - rec["PRE"]; rec["d_post2"] = rec["POST2"] - rec["PRE"]
        rec["ring_d_dur"] = rec["DUR_ring"] - rec["PRE_ring"]; rec["ring_d_post"] = rec["POST_ring"] - rec["PRE_ring"]; rows.append(rec)
MD = pd.DataFrame(rows); MD.to_csv(f"{SD}/patch_prepost_model.csv", index=False)
mw = cat.Mw.reindex(pids).values; lat = cat.lat_centroid.reindex(pids).values
SETS = {"Mw>=6.2 all": (mw >= 6.2), "Mw>=6.2 north": (mw >= 6.2) & (lat >= 46), "Mw>=6.2 south": (mw >= 6.2) & (lat < 46), "Mw>=6.0 all": (mw >= 6.0), "all post-2011": np.ones(len(pids), bool)}
# typical change of each field inside patches, per set
srows = []
for k in FIELDS:
    m = MD[MD.field == k].set_index("id").reindex(pids)
    for sname, sel in SETS.items():
        srows.append(dict(field=k, unit=m.unit.iloc[0], set=sname, n=int(sel.sum()), PRE_mean=np.nanmean(m.PRE[sel]), DUR_mean=np.nanmean(m.DUR[sel]), POST_mean=np.nanmean(m.POST[sel]), d_dur_mean=np.nanmean(m.d_dur[sel]), d_dur_median=np.nanmedian(m.d_dur[sel]), d_dur_frac_pos=np.nanmean(m.d_dur[sel] > 0),
                          d_post_mean=np.nanmean(m.d_post[sel]), d_post_frac_pos=np.nanmean(m.d_post[sel] > 0), ring_d_dur_mean=np.nanmean(m.ring_d_dur[sel])))
SM = pd.DataFrame(srows); SM.to_csv(f"{SD}/patch_prepost_model_summary.csv", index=False)
pd.set_option("display.width", 250)
print("typical change INSIDE the patch, set Mw>=6.2 north (DUR-PRE mean, frac positive | POST-PRE mean, frac pos | ring DUR-PRE):")
s = SM[SM.set == "Mw>=6.2 north"]; print(s[["field", "unit", "PRE_mean", "DUR_mean", "POST_mean", "d_dur_mean", "d_dur_frac_pos", "d_post_mean", "d_post_frac_pos", "ring_d_dur_mean"]].round(3).to_string(index=False))
# cross-patch correlation of dv/v change with field change, null from the dv/v shift draws
crows = []
for pn in PRODS:
    Tp = T[T["product"] == pn].set_index("id").reindex(pids)
    for metric in ("d_median_dur", "d_median_post", "d_median_dur_minus_ring"):
        obs_dv = Tp[metric].values; nul_dv = zn[f"{pn}|{metric}"]                      # (NDRAW, npatch)
        for sname, sel in SETS.items():
            if sel.sum() < 8: continue
            res = []
            for k in FIELDS:
                m = MD[MD.field == k].set_index("id").reindex(pids)
                for fch in ("d_dur", "d_post", "DUR"):
                    x = m[fch].values; ok = sel & np.isfinite(x) & np.isfinite(obs_dv)
                    if ok.sum() < 8 or np.nanstd(x[ok]) == 0: continue
                    rho = spearmanr(x[ok], obs_dv[ok]).correlation
                    nr = np.array([spearmanr(x[ok], nul_dv[d, ok], nan_policy="omit").correlation for d in range(nul_dv.shape[0])]); nr = nr[np.isfinite(nr)]
                    if len(nr) < 50: continue
                    res.append(dict(product=pn, metric=metric, set=sname, field=k, field_change=fch, n=int(ok.sum()), n_null=len(nr), rho=rho, null_sd=np.std(nr), p=float(np.mean(np.abs(nr) >= abs(rho))), _nr=nr))
            if not res: continue
            nmin = min(len(r["_nr"]) for r in res); NR = np.array([r["_nr"][:nmin] for r in res]); zmax = np.max(np.abs(NR) / NR.std(1, keepdims=True), 0)   # max |z| over fields per draw
            for r in res: r["p_fwer"] = float(np.mean(zmax >= abs(r["rho"]) / r["null_sd"])); r.pop("_nr"); crows.append(r)
CR = pd.DataFrame(crows); CR.to_csv(f"{SD}/patch_prepost_model_corr.csv", index=False)
print("\nbest cross-patch correlations per (product, metric, set) [Spearman rho, p vs dv/v shift null, FWER over fields]:")
best = CR.sort_values("p").groupby(["product", "metric", "set"]).head(1).sort_values(["product", "metric", "set"]); print(best[["product", "metric", "set", "field", "field_change", "n", "n_null", "rho", "null_sd", "p", "p_fwer"]].round(3).to_string(index=False))
print("\nany with p_fwer < 0.05:"); print(CR[CR.p_fwer < 0.05].round(3).to_string(index=False) if (CR.p_fwer < 0.05).any() else "none")
# figure: north Mw>=6.2, interface d_median_dur vs key field changes + bar of typical changes
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
KEY = [("slip_thrust_inc30", "DUR"), ("coulomb", "d_dur"), ("shear_stress", "d_dur"), ("normal_stress", "d_dur"), ("dilatation", "d_dur"), ("von_mises_stress", "d_dur"), ("slip_thrust_cum", "d_post"), ("coulomb", "d_post")]
pn = PRODS[0]; Tp = T[T["product"] == pn].set_index("id").reindex(pids); fig, ax = plt.subplots(2, 4, figsize=(18, 8.5), constrained_layout=True)
for a, (k, fch) in zip(ax.ravel(), KEY):
    m = MD[MD.field == k].set_index("id").reindex(pids)
    for sname, col in (("Mw>=6.2 north", "#1d4ed8"), ("Mw>=6.2 south", "#ea580c")):
        sel = SETS[sname]; a.plot(m[fch].values[sel], Tp.d_median_dur.values[sel], "o", color=col, ms=5, label=sname)
    r = CR[(CR["product"] == pn) & (CR.metric == "d_median_dur") & (CR.set == "Mw>=6.2 all") & (CR.field == k) & (CR.field_change == fch)]
    a.axhline(0, color="k", lw=.5); a.axvline(0, color="k", lw=.5); a.set_xlabel(f"{k} {fch} ({m.unit.iloc[0]})"); a.set_ylabel(f"dv/v DUR−PRE median (%), {pn}")
    if len(r): a.set_title(f"rho {r.rho.iloc[0]:+.2f}, p {r.p.iloc[0]:.2f}, pFWER {r.p_fwer.iloc[0]:.2f} (Mw≥6.2 all)", fontsize=9)
ax[0, 0].legend(fontsize=8); fig.suptitle("dv/v change inside the slip patch vs the model-variable change inside the same patch (per patch; unbridged catalog, post-2011)", fontsize=10)
fig.savefig("figures/sse_patch_prepost_model.png", dpi=130); print("wrote figures/sse_patch_prepost_model.png, patch_prepost_model*.csv")
