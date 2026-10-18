"""ab_hypo.py with the USGS Cascadia CVM v1.7 regional 1-D models added to MODELS (2026-09-28).
Same picking, Hyp2000 + hypoDD and bootstrap as ab_hypo.py; only extra model names:
USGS_OR, USGS_NW, USGS_SW, USGS_SVI (layer tops + Vp from data/velocity_models/usgs_cascadia_v1.7/profiles/<region>.csv).
Usage as ab_hypo.py, e.g. --models GSC,USGS_SW --vpvs 1.73,1.785"""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import ab_hypo as H
_V = "/home/jovyan/tremorferometry/data/velocity_models/usgs_cascadia_v1.7/profiles"
for _k, _r in [("USGS_OR", "oregon_TB2015"), ("USGS_NW", "royer_NW_WA"), ("USGS_SW", "royer_SW_WA"), ("USGS_SVI", "south_VI")]:
    _L = pd.read_csv(f"{_V}/{_r}.csv")
    # Hypoinverse CRH rejects a layer slower than the one above ("BAD DATA FOR LAYER CRUST MODEL"): the WA/VI
    # low-velocity zone (~18.5-29.5 km) is merged into the faster layer above it (velocity held constant).
    zt, vp, vmax = [], [], 0.0
    for z, v in zip(_L.ztop_km, _L.vp):
        if v < vmax:
            continue
        zt.append(float(z)); vp.append(float(v)); vmax = v
    H.MODELS[_k] = (zt, vp)
if __name__ == "__main__":
    H.main()
