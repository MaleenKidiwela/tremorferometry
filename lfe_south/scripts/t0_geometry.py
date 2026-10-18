"""Stage 0: southern truth inventory + station geometry.
Writes results/families_south.csv and results/station_geometry.csv"""
import pandas as pd, numpy as np, os, re, io

ROOT="/home/jovyan/tremorferometry"
OUT=f"{ROOT}/lfe_south/results"

# ---------- families ----------
fams=[]
# Shelly 2025 (27), locations from the USGS metadata XML abstract
xml=open(f"{ROOT}/lfe_catalog_repair/results/external/other/Shelly2025_metadata.xml").read()
blk=re.search(r"depth \(km below mean surface\):(.*?)The LFE family locations may be updated", xml, re.S).group(1)
for line in blk.strip().splitlines():
    p=line.split()
    if len(p)==4:
        fams.append(dict(source="Shelly2025", fam=p[0], lat=float(p[1]), lon=float(p[2]), depth=float(p[3])))
# Ducellier/Plourde families (template -> location)
tl=pd.read_csv(f"{ROOT}/catalogs/ducellier_catalog/data/Plourde_2015/templates_list.txt",
               sep=r"\s+", skiprows=1, header=None,
               names=["tpl","famnum","lat","lon","depth","eH","eZ","ntpl"])
have=set(os.listdir(f"{ROOT}/catalogs/ducellier_catalog/data/Ducellier/catalogs"))
for _,r in tl.iterrows():
    if r.tpl in have:
        fams.append(dict(source="Ducellier2022", fam=r.tpl, lat=r.lat, lon=r.lon, depth=r.depth,
                         plourde_fam=r.famnum))
F=pd.DataFrame(fams)

# event counts / spans
def shelly_counts():
    df=pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
    df.columns=[c.strip('%') for c in df.columns]
    t=pd.to_datetime(dict(year=df.year,month=df.month,day=df.day,hour=df.hour,
                          minute=df.minute,second=df.second.astype(int))).astype("datetime64[ns]")
    df["t"]=t
    return df
S=shelly_counts()
sc=S.groupby("idnum").t.agg(["size","min","max"])
for i,r in sc.iterrows():
    F.loc[(F.source=="Shelly2025")&(F.fam==str(i)),["n_events","t0","t1"]]=[r["size"],str(r["min"]),str(r["max"])]

dd=f"{ROOT}/catalogs/ducellier_catalog/data/Ducellier/catalogs"
for f in F.loc[F.source=="Ducellier2022","fam"]:
    n=0; lo=None; hi=None
    for cf in ("catalog_2004_2011.pkl","catalog_2007_2009.pkl"):
        p=os.path.join(dd,f,cf)
        if not os.path.exists(p): continue
        d=pd.read_pickle(p)
        if len(d)==0: continue
        t=pd.to_datetime(dict(year=d.year,month=d.month,day=d.day,hour=d.hour,
                              minute=d.minute,second=d.second.astype(int))).astype("datetime64[ns]")
        n+=len(d)
        lo=t.min() if lo is None else min(lo,t.min()); hi=t.max() if hi is None else max(hi,t.max())
        F.loc[F.fam==f,"n_"+cf[8:12]]=len(d)
    F.loc[F.fam==f,["n_events","t0","t1"]]=[n,str(lo),str(hi)]
F.to_csv(f"{OUT}/families_south.csv",index=False)
print("families:", F.groupby("source").agg(n=("fam","size"), ev=("n_events","sum")).to_string())

# ---------- stations ----------
C=pd.read_csv(f"{OUT}/fdsn_channels_south.csv")
C=C[C.sr>=20]
# band codes usable for 1-8 Hz: B,H,E,S (broad/high/extremely short/short period)
C["band"]=C.cha.str[0]; C["comp"]=C.cha.str[2]
C=C[C.band.isin(list("BHES"))]
C["end"]=C["end"].fillna("3000-01-01")
# a station-epoch is 3-C if it has >=3 distinct components in one band+loc
g=C.groupby(["net","sta","loc","band"]).agg(ncomp=("comp","nunique"), sr=("sr","max"),
        lat=("lat","first"), lon=("lon","first"), elev=("elev","first"), dep=("dep","max"),
        start=("start","min"), end=("end","max")).reset_index()
g=g[g.ncomp>=3]
# prefer one row per station: keep the band with the highest sample rate, then longest span
g["span"]=(pd.to_datetime(g.end.clip(upper="2262-01-01"))-pd.to_datetime(g.start)).dt.days
g=g.sort_values(["net","sta","sr","span"],ascending=[1,1,0,0])
ST=g.drop_duplicates(["net","sta"]).reset_index(drop=True)
ST.to_csv(f"{OUT}/stations_3c_south.csv",index=False)
print("3-C station-epochs:",len(g),"distinct 3-C stations:",len(ST))
print(ST.net.value_counts().to_dict())

# ---------- geometry ----------
R=6371.0
def hypo(lat1,lon1,dep1,lat2,lon2,elev2):
    dlat=(lat2-lat1)*111.19
    dlon=(lon2-lon1)*111.19*np.cos(np.radians((lat1+lat2)/2))
    dz=dep1+elev2/1000.0
    return np.sqrt(dlat**2+dlon**2+dz**2), np.sqrt(dlat**2+dlon**2)

rows=[]
for _,f in F.iterrows():
    d,e=hypo(f.lat,f.lon,f.depth,ST.lat.values,ST.lon.values,ST.elev.values)
    for i in np.where(d<=80)[0]:
        s=ST.iloc[i]
        rows.append(dict(source=f.source, fam=f.fam, net=s.net, sta=s.sta, loc=s.loc, band=s.band,
                         sr=s.sr, hypo_km=round(float(d[i]),2), epi_km=round(float(e[i]),2),
                         sta_lat=s.lat, sta_lon=s.lon, start=s.start, end=s.end))
G=pd.DataFrame(rows)
G.to_csv(f"{OUT}/station_geometry.csv",index=False)
print("\nfamily-station pairs <=80 km:",len(G))
for src,sub in G.groupby("source"):
    print(f"  {src}: {sub.fam.nunique()} families, {sub.sta.nunique()} distinct stations, "
          f"median {sub.groupby('fam').size().median():.0f} stations/family")
    print("   stations:", sorted(set(sub.net+"."+sub.sta)))
