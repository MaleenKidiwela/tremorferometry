#!/usr/bin/env python
"""Emit one MONTH's continuous download plan for the extension roster of a source."""
import os,sys,argparse
import pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True)
ap.add_argument("--month",required=True); ap.add_argument("--nsta",type=int,default=12)
ap.add_argument("--out",required=True)
a=ap.parse_args()
M=pd.read_csv(f"{RES}/templates_south.csv"); M=M[M.source==a.source]
# extension roster = the nsta stations nearest to the most families (a budget, not a filter:
# every family keeps every station it has, ranked by distance)
rank=M.groupby("netsta").hypo_km.median().sort_values()
stas=list(rank.index[:a.nsta])
m=pd.Timestamp(a.month)
days=pd.date_range(m,m+pd.offsets.MonthEnd(0),freq="D")
C=pd.read_csv(f"{RES}/fdsn_channels_south.csv")
C["start"]=pd.to_datetime(C.start); C["end"]=pd.to_datetime(C.end.fillna("2100-01-01"),errors="coerce").fillna(pd.Timestamp("2100-01-01"))
rows=[]
for ns in stas:
    net,sta=ns.split(".")
    c=C[(C.net==net)&(C.sta==sta)]
    if len(c)==0: continue
    if not ((c.start<=days[-1])&(c.end>=days[0])).any(): continue
    for d in days: rows.append((net,sta,str(d.date())))
pd.DataFrame(rows,columns=["net","sta","day"]).to_csv(a.out,index=False)
print(f"[plan] {a.source} {a.month}: {len(rows)} station-days over {len(stas)} stations")
