#!/usr/bin/env python
"""MENDO (2026-09-29): regional 1-D model from the USGS Cascadia CVM v1.7 for the Mendocino discovery.

The Yoon & Shelly (2024) 1-D layered model (the NCSN "MEN" model, used by Shelly et al. 2025) is NOT on
disk (layer values are not in the Shelly 2025 SI text), so -- per the project rule (memory
usgs-cvm17-velocity-model) -- the CVM v1.7 is used.  None of the existing profiles/ boxes is near
Mendocino (nearest oregon_TB2015, 44 N), so a new box is extracted with the SAME code
(data/velocity_models/usgs_cascadia_v1.7/scripts/extract_profiles.py: region_samples + layered +
write_crh, imported unchanged).  The CVM's southern edge is ~40.2 N, so the box is 40.2-40.7 N,
124.1-123.3 W: over the stations and just north of the 40.03-40.09 N LFE patch (PARTIAL coverage, like
royer_NVI: the LFE patch itself is outside the model).
-> data/velocity_models/usgs_cascadia_v1.7/profiles/mendo_MTJ{.csv,_raw.csv,.crh,_S.crh,_summary.csv}
"""
import os
import sys

V = "/home/jovyan/tremorferometry/data/velocity_models/usgs_cascadia_v1.7"
sys.path.insert(0, f"{V}/scripts")
import extract_profiles as E  # noqa: E402

NAME = "mendo_MTJ"
BOX = (40.2, 40.7, -124.1, -123.3, "Mendocino triple junction stations (PARTIAL: grid ends ~40.2 N)")


def main():
    if os.path.exists(f"{E.OUT}/{NAME}.csv"):
        print("exists", f"{E.OUT}/{NAME}.csv")
        return
    R = E.region_samples(BOX)
    R.to_csv(f"{E.OUT}/{NAME}_raw.csv", index=False, float_format="%.4f")
    if (R.n > 0).sum() == 0:
        raise SystemExit("NO COVERAGE")
    L, k, pos = E.layered(R)
    L.to_csv(f"{E.OUT}/{NAME}.csv", index=False)
    E.write_crh(f"{E.OUT}/{NAME}.crh", f"USGS CVM1.7 {NAME} P", L.ztop_km, L.vp)
    E.write_crh(f"{E.OUT}/{NAME}_S.crh", f"USGS CVM1.7 {NAME} S", L.ztop_km, L.vs)
    import pandas as pd
    pd.DataFrame([dict(region=NAME, n_layers=k, pos_0_35km=round(pos, 3),
                       ncells_deep=int(R[R.z_km > 20].n.median()),
                       lat_max_covered=round(R.lat_max_covered.max(), 3), box=str(BOX[:4]))]).to_csv(
        f"{E.OUT}/{NAME}_summary.csv", index=False)
    print(NAME, k, "layers, POS", round(pos, 3))
    print(L.to_string(index=False))


if __name__ == "__main__":
    main()
