#!/usr/bin/env python
"""Builder for res_catalog_sse/product_family24.npz and product_family35.npz (Merlin audit 2026-09-04 F6: provenance).
Family-level coda dv/v (= -eps, %) averaged per cell-day over certified families with cc_max > 0.7, on the model day grid,
packaged with the same keys as inversion_4d_daily.npz so sse_primary_test.py can consume it via DVV_PRODUCT.
2-4 s: from family_days.parquet (production certified rows, all gates); 3-5 s: from data/daily_dvv_<TAG>_Z_3to5.csv
restricted to the same certified families (184 stations have 3-5 s files; B011p90f40 does not)."""
import numpy as np, pandas as pd, glob, os
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells); gi = pd.Series(np.arange(G), index=days)
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"])
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=md, columns=list(zp["sta_tags"])).reindex(days).values
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell", "date", "dvv"]); cpos = {c: i for i, c in enumerate(cells)}; fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv
famset = {t: set(g.fam) for t, g in fd.groupby("tag")}; famcell = fd[["tag", "fam", "ci"]].drop_duplicates()
c24 = fd.groupby(["ci", "date"]).v.mean().reset_index(); ok = c24.date.isin(days); T24 = np.full((G, M), np.nan, np.float32); T24[gi[c24.date[ok]].values, c24.ci[ok].values] = c24.v[ok].values
rows = []
for f in glob.glob("data/daily_dvv_*_Z_3to5.csv"):
    tag = os.path.basename(f).split("daily_dvv_")[1].split("_Z_3to5")[0]
    if tag not in famset: continue
    x = pd.read_csv(f, usecols=["patch", "date", "dvv", "cc_max"]); x = x[(x.cc_max > 0.7) & (x.patch.isin(famset[tag]))]; x["tag"] = tag; x["date"] = pd.to_datetime(x.date); x["v"] = -x.dvv * 100; rows.append(x.rename(columns={"patch": "fam"})[["tag", "fam", "date", "v"]])
f35 = pd.concat(rows).merge(famcell, on=["tag", "fam"]); c35 = f35.groupby(["ci", "date"]).v.mean().reset_index(); ok = c35.date.isin(days); T35 = np.full((G, M), np.nan, np.float32); T35[gi[c35.date[ok]].values, c35.ci[ok].values] = c35.v[ok].values
for name, T in (("family24", T24), ("family35", T35)):
    np.savez_compressed(f"{SD}/product_{name}.npz", days=days.astype(str).values, MODEL=T, SIGM=np.full_like(T, np.nan), SITES=SITES.astype(np.float32), cell=cells, clat=zc["cell_lat"], clon=zc["cell_lon"], depth=zc["iface_depth_km"], sta_tags=zp["sta_tags"], lam_f=0.0, lam_s=0.0, inR=np.zeros(G, bool))
    print(name, "finite cell-days", int(np.isfinite(T).sum()))
