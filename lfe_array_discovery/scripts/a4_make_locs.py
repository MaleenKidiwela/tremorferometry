#!/usr/bin/env python
"""A4 location lists: real = Royer & Bostock (2014) Fig 4 Hyp2000 templates (NW, SW, NVI); nullR =
the same cloud rotated about its centroid by a random angle and shifted uniformly within a disc of
the cloud's RMS radius (the ab_truth_score null), one realisation per subarray, depths kept.
Missing depth -> subarray median (flag depth_filled); depth_spread_km > 1 -> depth_uncertain.
-> results/a4_seed/locs_all.csv, locs_pilot.csv (10 per subarray, seed 0, + their nullR images)"""
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(os.path.dirname(HERE), "results")
T = pd.read_csv(f"{R}/ab_royer/truth/royer2014_templates_digitized.csv")
T = T[(T.locator == "Hyp2000") & T.subarray.isin(["NW", "SW", "NVI"])].reset_index(drop=True)
rows = []
lid = 0
for sub in ["NW", "SW", "NVI"]:
    t = T[T.subarray == sub].copy()
    med = t.depth_km.median()
    t["depth_filled"] = t.depth_km.isna()
    t["depth"] = t.depth_km.fillna(med)
    t["depth_uncertain"] = (t.depth_spread_km > 1) | t.depth_filled
    lat0, lon0 = t.lat.mean(), t.lon.mean()
    kx = 111.19 * np.cos(np.radians(lat0))
    x, y = (t.lon - lon0) * kx, (t.lat - lat0) * 111.19
    rng = np.random.default_rng(abs(hash(sub)) % 2**32 if False else {"NW": 11, "SW": 12, "NVI": 13}[sub])
    th = rng.uniform(0, 2 * np.pi); Rr = np.sqrt(np.mean(x ** 2 + y ** 2)); rr = Rr * np.sqrt(rng.uniform()); ph = rng.uniform(0, 2 * np.pi)
    nx = np.cos(th) * x - np.sin(th) * y + rr * np.cos(ph)
    ny = np.sin(th) * x + np.cos(th) * y + rr * np.sin(ph)
    nlon, nlat = lon0 + nx / kx, lat0 + ny / 111.19
    D = np.hypot(nx.values[:, None] - x.values[None], ny.values[:, None] - y.values[None]).min(1)
    for i, r in enumerate(t.itertuples()):
        base = dict(subarray=sub, truth_row=int(r.Index), depth=r.depth, depth_uncertain=bool(r.depth_uncertain),
                    depth_filled=bool(r.depth_filled), thesis_id=r.template_id)
        rows.append(dict(lid=lid, cond="real", lon=r.lon, lat=r.lat, d_nearest_truth_km=0.0, **base)); lid += 1
        rows.append(dict(lid=lid, cond="nullR", lon=float(nlon.iloc[i]), lat=float(nlat.iloc[i]),
                         d_nearest_truth_km=float(D[i]), rot_deg=float(np.degrees(th)), shift_km=float(rr), **base)); lid += 1
L = pd.DataFrame(rows)
os.makedirs(f"{R}/a4_seed", exist_ok=True)
L.to_csv(f"{R}/a4_seed/locs_all.csv", index=False)
pil = []
for sub in ["NW", "SW", "NVI"]:
    tr = L[(L.subarray == sub) & (L.cond == "real")].truth_row.unique()
    pick = np.random.default_rng(0).choice(tr, 10, replace=False)
    pil.append(L[L.truth_row.isin(pick) & (L.subarray == sub)])
P = pd.concat(pil); P.to_csv(f"{R}/a4_seed/locs_pilot.csv", index=False)
print(L.groupby(["subarray", "cond"]).agg(n=("lid", "size"), dnt=("d_nearest_truth_km", "median")))
print(P.groupby(["subarray", "cond"]).size())
