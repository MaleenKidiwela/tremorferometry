#!/usr/bin/env python
"""STAGE R1a -- build the two Ducellier reproduction BEDS and their download plans.

BED A (primary, fidelity): HIS EXACT permanent set, catalogs/ducellier_catalog/data/Ducellier/
       stations_permanent.txt, channels and epochs as he lists them.  ~70 % of it is
       SINGLE-COMPONENT VERTICAL (EHZ/SHZ); those stations are kept and each contributes its one
       component to the per-channel cc mean, exactly as find_all_LFEs_parallel.py does.
BED B (same-plus-more): bed A + every permanent station in the FDSN metadata within 80 km of one
       of his 66 family hypocentres that is alive in the window and was NOT in his set.
Window: the 2008 April southern ETS (2008-04-01 .. 2008-05-10) -- the richest bed in his
catalogue (64 of 66 families with >= 20 events in 30 days).
"""
import os,sys,numpy as np,pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
DUC=f"{ROOT}/catalogs/ducellier_catalog"
W0,W1="2008-04-01","2008-05-10"
def his_stations():
    rows=[]
    for ln in open(f"{DUC}/data/Ducellier/stations_permanent.txt"):
        p=ln.split()
        if len(p)<10: continue
        rows.append(dict(sta=p[0],net=p[1],chans=p[2],loc=("" if p[3]=="--" else p[3]),prov=p[4],
                         lat=float(p[5]),lon=float(p[6]),on=p[7],off=p[8],dt=float(p[9])))
    return pd.DataFrame(rows)
def main():
    H=his_stations()
    H["nch"]=H.chans.str.count(",")+1
    H["zonly"]=(H.nch==1).astype(int)
    H.to_csv(f"{RES}/duc_bedA_stations.csv",index=False)
    d0,d1=pd.Timestamp(W0),pd.Timestamp(W1)
    alive=H[(pd.to_datetime(H.on)<=d1)&(pd.to_datetime(H.off)>=d0)]
    print(f"BED A: {len(H)} rows / {H.sta.nunique()} stations; Z-only rows {H.zonly.sum()} "
          f"({100*H.zonly.mean():.0f} %); alive in {W0}..{W1}: {alive.sta.nunique()}")
    F=pd.read_csv(f"{RES}/families_south.csv"); F=F[F.source=="Ducellier2022"]
    C=pd.read_csv(f"{RES}/fdsn_channels_south.csv")
    C=C[(C.sr>=20)&C.cha.str[0].isin(list("BHES"))&C.net.isin({"NC","BK","PB","UW","CC","US","UO"})].copy()
    C["start"]=pd.to_datetime(C.start)
    C["end"]=pd.to_datetime(C.end.fillna("2100-01-01"),errors="coerce").fillna(pd.Timestamp("2100-01-01"))
    A=C[(C.start<=d1)&(C.end>=d0)]
    g=A.groupby(["net","sta"]).agg(lat=("lat","first"),lon=("lon","first"),
        chans=("cha",lambda x:",".join(sorted(set(x)))),sr=("sr","max")).reset_index()
    keep=[]
    for r in g.itertuples():
        dmin=min(np.sqrt(((r.lat-f.lat)*111.19)**2+((r.lon-f.lon)*111.19*np.cos(np.radians(f.lat)))**2
                         +f.depth**2) for f in F.itertuples())
        if dmin<=80: keep.append((r.net,r.sta,round(dmin,1),r.chans,r.sr))
    B=pd.DataFrame(keep,columns=["net","sta","min_hypo_km","chans","sr"])
    B["in_bedA"]=B.sta.isin(set(H.sta)).astype(int)
    B.to_csv(f"{RES}/duc_bedB_stations.csv",index=False)
    print(f"BED B: {len(B)} stations alive in the window within 80 km "
          f"({int(B.in_bedA.sum())} of them in bed A, {int((1-B.in_bedA).sum())} new)")
    days=[str(d.date()) for d in pd.date_range(W0,W1,freq="D")]
    rows=[]
    for r in alive.drop_duplicates("sta").itertuples():
        for d in days: rows.append((r.net,r.sta,d,"A"))
    for r in B[B.in_bedA==0].itertuples():
        for d in days: rows.append((r.net,r.sta,d,"B"))
    P=pd.DataFrame(rows,columns=["net","sta","day","bed"]).drop_duplicates(["net","sta","day"])
    P.to_csv(f"{RES}/duc_download_plan.csv",index=False)
    print(f"station-days: {len(P)} (bed A {int((P.bed=='A').sum())}, new in bed B "
          f"{int((P.bed=='B').sum())}), est {len(P)*0.025:.0f} GB")
if __name__=="__main__": main()
