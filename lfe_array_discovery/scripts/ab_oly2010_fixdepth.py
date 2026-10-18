"""Olympic 2010 real families (site-count >=3, OLY2010_cvm/site_score_real.csv): Hyp2000 with depth FIXED to the McCrory 2012
plate interface under the family (USGS_NW model, Vp/Vs 1.823), same picks as ab_oly2010_cvm.py (2026-09-28).
Starting epicentre for the interface depth: the USGS free-depth solution if it exists, else the L9 location.
-> results/ab_royer/OLY2010_cvm/loc_USGS_NW_fixed.csv"""
import re
import os, sys, shutil
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import ab_hypo_cvm  # noqa: F401  registers USGS_NW
import ab_hypo as H
import ab_common as A
import ab_oly2010_cvm as O
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_oregon/scripts")
import xs_locate as XL
OUT = f"{O.RES}/ab_royer/OLY2010_cvm"
A.FS = 20.0
code = O.sites(); stas = sorted(set(code.values()))
P = pd.read_csv(f"{OUT}/picks.csv"); P = P[P.population == "real"].drop(columns="population")
good = set(pd.read_csv(f"{OUT}/site_score_real.csv").query("n_imp_sites>=3").grp)
C = pd.read_csv(f"{OUT}/compare.csv").set_index("grp")
rows = []
for g in sorted(good):
    r = C.loc[g]
    lat, lon = (r.USGS_NW_lat, r.USGS_NW_lon) if r.get("USGS_NW_located") is True and np.isfinite(r.USGS_NW_lat) else (r.lat, r.lon)
    if not (np.isfinite(lat) and np.isfinite(lon)):
        continue
    z = XL.iface_depth(lat, lon)
    wd = f"{OUT}/run_fixed/{g}"; shutil.rmtree(wd, ignore_errors=True)
    ids = H.write_hyp_inputs(wd, P[P.grp == g], stas, 1.823, O.RMIN, model="USGS_NW")
    if not ids:
        rows.append(dict(grp=g, located=False, z_iface=z, reason="too few picks")); continue
    s = re.sub(r"(?m)^ZTR .*$", f"ZTR {z:.1f} T", open(f"{wd}/run.hyp").read()); open(f"{wd}/run.hyp", "w").write(s)  # ab_hypo now writes "ZTR 35 F"
    Hy, _ = H.run_hyp(wd)
    if len(Hy):
        rows.append(dict(Hy.iloc[0].to_dict(), located=True, z_iface=z))
    else:
        rows.append(dict(grp=g, located=False, z_iface=z, reason="hyp2000 no solution"))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/loc_USGS_NW_fixed.csv", index=False)
L = D[D.located]
print(f"real {len(good)}; input {len(D)}; located fixed {len(L)}; erh median {L.erh.median():.1f}; erh<5 {(L.erh<5).sum()}")
print(D[~D.located].reason.value_counts().to_dict() if (~D.located).any() else "")
