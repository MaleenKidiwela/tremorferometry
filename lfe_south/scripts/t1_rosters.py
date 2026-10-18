"""Stage-1 rosters + download plan.

Roster rule (geometry only, no quality gate):
  permanent networks (NC,BK,PB,UW,CC,US,UO), band B/H/E/S, sample rate >=20 Hz,
  hypocentral distance <= 80 km of the published family hypocentre,
  alive on the sampled detection days AND alive 2026-06-01 (so the same station can carry the
  template forward in stage 3),  nearest K_NEAR kept as the download budget.
Component count is RECORDED, not required: the southern permanent network is largely
vertical-only short-period in the 2004-2011 era and Ducellier's own catalogue used those stations.
"""
import pandas as pd, numpy as np, os, json
ROOT='/home/jovyan/tremorferometry'; OUT=f'{ROOT}/lfe_south/results'
K_NEAR=10; MAXD=80.0; TARGET=400; DAYCAP=80
PERM_NETS={'NC','BK','PB','UW','CC','US','UO'}
NOW=pd.Timestamp('2026-06-01')

C=pd.read_csv(f'{OUT}/fdsn_channels_south.csv')
C=C[(C.sr>=20)&C.cha.str[0].isin(list('BHES'))&C.net.isin(PERM_NETS)].copy()
C['band']=C.cha.str[0]; C['comp']=C.cha.str[2]
C['start']=pd.to_datetime(C.start)
C['end']=pd.to_datetime(C.end.fillna('2100-01-01'),errors='coerce').fillna(pd.Timestamp('2100-01-01'))
ALIVE_NOW=set(map(tuple,C[(C.start<=NOW)&(C.end>=NOW)][['net','sta']].drop_duplicates().values))

F=pd.read_csv(f'{OUT}/families_south.csv')

def alive(day):
    c=C[(C.start<=day)&(C.end>=day)]
    g=(c.groupby(['net','sta','band']).agg(ncomp=('comp','nunique'),sr=('sr','max'),
        lat=('lat','first'),lon=('lon','first'),elev=('elev','first'),
        locs=('loc',lambda x: ','.join(sorted(set(map(str,x))))),
        chas=('cha',lambda x: ','.join(sorted(set(x))))).reset_index())
    g=g.sort_values(['ncomp','sr'],ascending=False).drop_duplicates(['net','sta'])
    return g[[ (n,s) in ALIVE_NOW for n,s in zip(g.net,g.sta)]]

def hypo(f,lat,lon,elev):
    dlat=(lat-f.lat)*111.19; dlon=(lon-f.lon)*111.19*np.cos(np.radians(f.lat))
    return np.sqrt(dlat**2+dlon**2+(f.depth+elev/1000.0)**2)

det={}
S=pd.read_csv(f'{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv')
S.columns=[c.strip('%') for c in S.columns]
S['t']=(pd.to_datetime(dict(year=S.year,month=S.month,day=S.day,hour=S.hour,minute=S.minute,
        second=S.second.astype(int))).astype('datetime64[ns]')+pd.to_timedelta(S.second%1,unit='s'))
for i,g in S.groupby('idnum'): det[('Shelly2025',str(i))]=g[['t']].copy()
A=pd.read_parquet(f'{OUT}/ducellier_detections.parquet')
A=A[A.cat=='net']
for f,g in A.groupby('fam'): det[('Ducellier2022',f)]=g[['t']].copy()

rows=[];plan=set();famdays={};famtimes=[]
for _,f in F.iterrows():
    key=(f.source,str(f.fam))
    if key not in det: continue
    d=det[key].copy(); d['day']=d.t.dt.floor('D')
    per=d.groupby('day').size().sort_values(ascending=False)
    k=min(int((per.cumsum()<TARGET).sum())+1,len(per),DAYCAP)
    days=list(per.index[:k])
    famdays[f'{f.source}|{f.fam}']=[str(x.date()) for x in days]
    famtimes.append(d[d.day.isin(days)].assign(source=f.source,fam=f.fam))
    st=alive(days[len(days)//2]).copy()
    st['hypo']=hypo(f,st.lat.values,st.lon.values,st.elev.values)
    st=st[st.hypo<=MAXD].sort_values('hypo').head(K_NEAR)
    for _,s in st.iterrows():
        rows.append(dict(source=f.source,fam=f.fam,net=s.net,sta=s.sta,band=s.band,sr=s.sr,
                         ncomp=s.ncomp,chas=s.chas,locs=s.locs,hypo_km=round(s.hypo,1),
                         sta_lat=s.lat,sta_lon=s.lon,elev=s.elev,n_days=len(days),
                         n_det=int(per.iloc[:k].sum())))
        for dd in days: plan.add((s.net,s.sta,str(dd.date())))
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/stage1_rosters.csv',index=False)
pd.concat(famtimes,ignore_index=True).to_parquet(f'{OUT}/stage1_times.parquet')
json.dump(famdays,open(f'{OUT}/stage1_days.json','w'),indent=0)
P=pd.DataFrame(sorted(plan),columns=['net','sta','day'])
P.to_csv(f'{OUT}/stage1_download_plan.csv',index=False)
print(R.groupby('source').agg(fams=('fam','nunique'),stns=('sta','nunique'),med_hypo=('hypo_km','median'),
      med_det=('n_det','median'),med_days=('n_days','median'),med_ncomp=('ncomp','median')))
for src,sub in R.groupby('source'):
    ns=sub.groupby('fam').size()
    print(src,'stations/family: min',ns.min(),'median',int(ns.median()))
    print('  ',sorted(set(sub.net+'.'+sub.sta+'('+sub.ncomp.astype(str)+'c)')))
print('\nstation-days:',len(P),'unique days',P.day.nunique(),'unique stations',P.sta.nunique())
print('est GB @35MB/station-day:',round(len(P)*0.035,1)); print(P.groupby('net').size().to_dict())
