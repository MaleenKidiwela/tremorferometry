#!/usr/bin/env python
"""s10_kelb_merge -- PLAN §9: the KELB(b)->KLNB gain-normalised check passed for family 87 only (tables/kelb_klnb_check.csv);
build the merged one-station series PO.KLNBm (KELB 2003-07-31..2004-06-04 + KLNB 2004-06-04..2006-06-22) for fam 87.
Windows are already RMS-normalised, so gain cancels. Writes win/PO.KLNBm.npz and tables/row_station_kelbm.parquet."""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s10_common as S
K = pd.read_csv(f"{S.TB}/kelb_klnb_check.csv")
fams = K[K["merge"].astype(bool)].fam.tolist()
RS = pd.read_parquet(f"{S.TB}/row_station.parquet"); X = pd.read_parquet(f"{S.TB}/rows.parquet", columns=["row_id", "fam"])
RS = RS.merge(X, on="row_id")
m = RS[(((RS.netsta == "PO.KELB") & (RS.epoch == "Eb")) | ((RS.netsta == "PO.KLNB") & (RS.epoch == "E1"))) & RS.fam.isin(fams)].copy()
Ws, rids = [], []
for ns in ("PO.KELB", "PO.KLNB"):
    z = np.load(f"{S.WIN}/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
    r = m[m.netsta == ns].row_id
    r = r[r.isin(pos)]
    Ws.append(z["W"][[pos[x] for x in r]]); rids.append(r.to_numpy())
np.savez(f"{S.WIN}/PO.KLNBm.npz", row_id=np.concatenate(rids), W=np.concatenate(Ws), tau=S.TAU)
m = m.assign(netsta="PO.KLNBm", epoch="Em").drop(columns=["fam"])
m.to_parquet(f"{S.TB}/row_station_kelbm.parquet", index=False)
print("merged families", fams, "rows", len(m))
