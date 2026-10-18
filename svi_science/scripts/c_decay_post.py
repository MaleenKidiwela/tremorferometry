"""C post-checks (added after the pre-registered run; labelled NOT pre-registered in RESULTS_C.md):
1. Is a patch's active duration just the time the front needs to cross the patch? transit = along-strike span (q10-q90 of
   the cell's events) / |episode front speed|.
2. Durations at single sources (SVI family, NW source) instead of patches: removes the patch-width transit.
3. Plateau fits: how many plateau winners have a real plateau (T > 2w and T > 0.5 d) vs an exponential-like tail.
"""
import json
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import c_decay as C

R = C.R; OUT = C.OUT
K = C.K.set_index(["region", "episode"])
res = {}

CT = pd.read_csv(f"{OUT}/cells_duration.csv")
E = pd.read_parquet(f"{OUT}/events_lag.parquet")
rows = []
for _, r in CT[CT.run.isin(["headline", "S13"])].dropna(subset=["D80"]).iterrows():
    g = E[(E.region == r.region) & (E.run == r.run) & (E.ep == r.ep) & (E.patch == r.patch)]
    span = np.quantile(g.as_km, 0.9) - np.quantile(g.as_km, 0.1)
    sp = abs(K.loc[(r.region, r.ep), "lfe_speed"])
    rows.append(dict(region=r.region, run=r.run, ep=r.ep, patch=r.patch, D80=r.D80, span_km=span, speed=sp,
                     transit_d=span / sp if sp > 0.5 else np.nan, dcl=r.dcl, dd=r.dd))
T = pd.DataFrame(rows); T["excess_d"] = T.D80 - T.transit_d
res["transit"] = {f"{k[0]}_{k[1]}": dict(n=int(len(g)), median_D80=float(g.D80.median()), median_transit=float(g.transit_d.median()),
                                         median_excess=float(g.excess_d.median()),
                                         frac_D80_gt_2x_transit=float((g.D80 > 2 * g.transit_d).mean()),
                                         rho_D80_transit=float(C.stats.spearmanr(g.D80, g.transit_d, nan_policy="omit")[0]))
                  for k, g in T.groupby(["region", "run"])}
T.to_csv(f"{OUT}/cells_transit.csv", index=False)

# single-source durations
Es = pd.read_parquet(f"{R}/s3/events_svi.parquet")
Es = Es[(Es.series == "data") & (Es.P_same >= 0.75) & Es.counted_once & Es.ep.isin(C.EPS["SVI"])].copy()
Es["patch"] = "F" + Es.unit.astype(str)
En = pd.read_parquet(f"{R}/s3/events_nw.parquet")
En = En[(En.series == "data") & (En.P_same >= 0.75) & En.selected_v2 & En.counted_once & En.ep.isin(C.EPS["NW"])].copy()
En["patch"] = En.source
res["single_source"] = {}; SS = []
for reg, X in [("SVI", Es), ("NW", En)]:
    X = X[["time", "ep", "as_km", "dd_km", "patch", "P_same"]].assign(region=reg, qualC=False)
    X = C.add_lag(C.prep(X, reg), reg)
    ct = C.cell_table(X, reg, "single_source").dropna(subset=["D80_n30"]); SS.append(ct)
    rho, p, ng = C.perm_within(ct, "D80_n30")
    res["single_source"][reg] = dict(n_cells=int(len(ct)), n_episodes=ng, rho=rho, p=p,
                                     median_D80n30_by_class=ct.groupby("dcl").D80_n30.median().to_dict(),
                                     median_D80_by_class=ct.groupby("dcl").D80.median().to_dict(),
                                     n_by_class=ct.groupby("dcl").size().to_dict(),
                                     median_D80n30=float(ct.D80_n30.median()))
pd.concat(SS).to_csv(f"{OUT}/cells_duration_single_source.csv", index=False)

CF = pd.read_csv(f"{OUT}/cell_fits.csv")
pl = CF[CF.winner == "plateau_stop"]
real = (pl.plateau_T > 2 * pl.plateau_w) & (pl.plateau_T > 0.5)
res["plateau_detail"] = {f"{k[0]}_{k[1]}": dict(n_plateau_winners=int(len(g)), n_real_plateau=int(real[g.index].sum()),
                                                median_T=float(g.plateau_T.median()), median_w=float(g.plateau_w.median()))
                         for k, g in pl.groupby(["region", "run"])}
om = CF[CF.winner == "omori"]
res["omori_detail"] = {f"{k[0]}_{k[1]}": dict(n=int(len(g)), p_values=g.omori_p.round(2).tolist()) for k, g in om.groupby(["region", "run"])}
json.dump(res, open(f"{OUT}/result_post.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
