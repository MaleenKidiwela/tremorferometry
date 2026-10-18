#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6): select earthquakes for the differential Wadati per 30-km along-strike zone.

  python s11_eq_select.py
writes tomo/eq/selection.csv (one row per selected event, priority order = download order) and
tomo/eq/selection_counts.csv (per region x zone x class).
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s11_tomo_wadati import along_strike  # noqa: E402

EPI_MAX = 50.0          # km, epicentre to nearest LFE family (crustal, intraslab)
EPI_LAYER = 30.0        # km, layer events (tier 2)
CAP_CRUST_WA = 120      # per (region, zone), largest magnitude first
CAP_LAYER = 50


def lfe_families():
    F = []
    for s, r in (("BV3", "WA_NW"), ("CV3", "WA_SW"), ("EallV3", "SVI")):
        P = pd.read_csv(f"{T}/tables/picks_{s}.csv", low_memory=False)
        f = P.groupby("fam")[["lat", "lon", "depth_km"]].first().dropna()
        f = f[f.depth_km >= 20]                  # mislocated shallow groups do not define the LFE area
        f["reg"] = r
        F.append(f)
    return pd.concat(F)


def main():
    F = lfe_families()
    d = pd.read_parquet(f"{ROOT}/svi_science/results/g2_local/tables/eq_local.parquet")
    d = d[(d.type == "earthquake") & (d.removed == "")].copy()           # duplicates already merged in the table
    # ComCat id carrying the phase-data: the row itself, or its merged ComCat twin for NRCan-primary rows
    d["cc_id"] = np.where(d.src == "ComCat", d.id,
                          np.where(d.dup_of.str.match(r"^(uw|us)"), d.dup_of, ""))
    lat0 = np.radians(47.5)
    D = np.hypot((d.lon.values[:, None] - F.lon.values[None]) * 111.19 * np.cos(lat0),
                 (d.lat.values[:, None] - F.lat.values[None]) * 111.19)
    d["epi_lfe"] = D.min(1)
    d["reg"] = F.reg.values[D.argmin(1)]
    d["along"] = along_strike(d.lat.values, d.lon.values)
    d["zone"] = (np.floor(d.along / 30) * 30).astype(int)
    crust = (d.depth <= 25) & (d.cls == "crustal")
    slab = d.cls == "intraslab"
    layer = (d.depth > 25) & (d.depth <= 45) & d.cls.isin(["overriding_deep", "near_interface"])
    d["kind"] = np.select([crust, slab, layer], ["crust", "slab", "layer"], "other")
    svi = d.reg == "SVI"
    sel_slab = slab & (d.epi_lfe <= EPI_MAX) & (d.mag >= 1.0)
    sel_crust_svi = crust & svi & (d.epi_lfe <= EPI_MAX) & (d.mag >= 1.0)
    c_wa = d[crust & ~svi & (d.epi_lfe <= EPI_MAX) & (d.mag >= 1.5)].copy()
    c_wa["rk"] = c_wa.groupby(["reg", "zone"]).mag.rank(ascending=False, method="first")
    lay = d[layer & (d.epi_lfe <= EPI_LAYER) & (d.mag >= 1.5)].copy()
    lay["rk"] = lay.groupby(["reg", "zone"]).mag.rank(ascending=False, method="first")
    parts = []
    a = d[sel_slab].copy(); a["tier"] = 1; a["rk"] = a.groupby(["reg", "zone"]).mag.rank(ascending=False, method="first")
    parts.append(a)
    b = d[sel_crust_svi].copy(); b["tier"] = 1; b["rk"] = b.groupby("zone").mag.rank(ascending=False, method="first")
    parts.append(b)
    c = c_wa[c_wa.rk <= CAP_CRUST_WA].copy(); c["tier"] = 1; parts.append(c)
    e = lay[lay.rk <= CAP_LAYER].copy(); e["tier"] = 2; parts.append(e)
    S = pd.concat(parts)
    S = S[~S.id.duplicated()]
    # download order: tier, then zone-balanced rank (round robin over zones, largest magnitude first)
    S = S.sort_values(["tier", "rk", "mag"], ascending=[True, True, False]).reset_index(drop=True)
    S["source"] = np.where(S.cc_id != "", "comcat", "none_public")
    cols = ["id", "cc_id", "src", "source", "tier", "kind", "reg", "zone", "along", "t", "t_PT", "lat", "lon", "depth",
            "mag", "magType", "cls", "slab2_km", "dz_slab2", "epi_lfe", "rk"]
    os.makedirs(f"{T}/eq", exist_ok=True)
    S[cols].to_csv(f"{T}/eq/selection.csv", index=False)
    cnt = S.groupby(["reg", "zone", "kind", "source"]).size().unstack(["kind", "source"], fill_value=0)
    cnt.to_csv(f"{T}/eq/selection_counts.csv")
    print(cnt.to_string())
    print("total", len(S), S.groupby(["tier", "source"]).size().to_dict())


if __name__ == "__main__":
    main()
