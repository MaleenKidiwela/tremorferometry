"""ab_hypo.py for MENDO (2026-09-29): same picking, Hyp2000 + hypoDD (Table S1) and station bootstrap,
with the USGS CVM v1.7 Mendocino profile (scripts/mendo_cvm_profile.py ->
data/velocity_models/usgs_cascadia_v1.7/profiles/mendo_MTJ.csv) added as model USGS_MENDO and made the
HEADLINE (bootstrap + loc_headline.csv) at its 0-35 km path-average Vp/Vs (mendo_MTJ_summary.csv
pos_0_35km, 1.71).  Columns: Vp/Vs 1.78 (Shelly 2025's value) and the GSC model.  UPDATE 09-29: headline switched to the NCSN MEN model (see H.MODELS["MEN"] below).  The profile increases monotonically (no LVZ to merge, see
ab_hypo_cvm.py).
Usage as ab_hypo.py:  python scripts/ab_hypo_mendo.py --sub MENDO --det ... --out ...
                      (defaults here: --models USGS_MENDO,GSC --vpvs <pos>,1.78)"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_hypo as H  # noqa: E402

_V = "/home/jovyan/tremorferometry/data/velocity_models/usgs_cascadia_v1.7/profiles"
_L = pd.read_csv(f"{_V}/mendo_MTJ.csv")
zt, vp, vmax = [], [], 0.0
for z, v in zip(_L.ztop_km, _L.vp):
    if v < vmax:
        continue
    zt.append(float(z)); vp.append(float(v)); vmax = v
H.MODELS["USGS_MENDO"] = (zt, vp)
POS = round(float(pd.read_csv(f"{_V}/mendo_MTJ_summary.csv").pos_0_35km.iloc[0]), 2)
# MEN model = Yoon & Shelly (2024) TSR supplement Table S3 (1-D layered, adapted from NCSN men.crt), the model
# Shelly 2025 located with (added 2026-09-29; data/velocity_models/yoon_shelly_2024/YoonShelly2024_TSR_supplement.pdf).
# Their 23 layer tops are cut to 15 (hypoDD MAXLAY=15): their 0.6-km steps above 4.2 km taken every 1.2 km; all
# values are theirs unchanged; 4.2-24 km (incl. the 22 km 6.9 and 24 km 7.9 tops) exactly as tabulated.
H.MODELS["MEN"] = ([0.0, 1.2, 2.4, 3.6, 4.2, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0],
                   [3.5, 3.9429, 4.3857, 4.8286, 5.05, 5.2371, 5.4449, 5.6528, 5.8607, 6.0685, 6.2764,
                    6.4843, 6.6921, 6.9, 7.9])
H.HEAD = ("MEN", 1.78, 1.5)      # headline = Shelly 2025's model and Vp/Vs 1.78; Yoon & Shelly Vs column = 1.732; USGS_MENDO + GSC = columns
_boot = H.bootstrap


def _boot_safe(*a, **k):
    """ab_hypo.bootstrap returns a column-less frame when no family has >= 10 resamples, and the
    merge on 'grp' then fails; give it the columns (all NaN) instead."""
    B = _boot(*a, **k)
    if not len(B):
        B = pd.DataFrame(columns=["grp", "boot_km", "boot_dz", "boot_n"])
    return B


H.bootstrap = _boot_safe
if __name__ == "__main__":
    if "--models" not in sys.argv:
        sys.argv += ["--models", "MEN,USGS_MENDO,GSC"]
    if "--vpvs" not in sys.argv:
        sys.argv += ["--vpvs", f"1.78,1.73,{POS}"]
    H.main()
