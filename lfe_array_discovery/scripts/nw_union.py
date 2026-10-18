"""NW union of real families across runs (2026-09-29, user: combine runs for more families).
Runs: NW_lag15 (final base), NW_mute, NW_fair. Real = site-count score >=3 (NW_mute: class_final lfe). Duplicate test between two families of different
runs on their 2-8 Hz stacks (hypo_28/stacks_loc.npz): per common station, N/E cross-correlation, lag of max within +-3 s; duplicate if >=2 stations with
cc>=0.5 and lag spread <=0.15 s (same rule as lfe_oregon/scripts/xs_pair_distinct.py). Candidates: pairs whose locations are <=10 km apart (or either unlocated).
Greedy: start from NW_lag15; add a family from another run only if it is not a duplicate of any family already in the union."""
import numpy as np, pandas as pd
from scipy.signal import correlate
R='results/ab_royer'; FS=40.0; ML=int(3*FS)
EXTRA={'NW_run1':dict(stk=f'{R}/NW/hypo_28/stacks_loc.npz',loc=f'{R}/NW/hypo_allsta/loc_GSC_1.73_R1.5.csv',lab=f'{R}/NW/hypo_28/review/eye_class.csv'),
       'NW_qc':dict(stk=f'{R}/NW/hypo_qc28/stacks_loc.npz',loc=f'{R}/NW/hypo_qc/loc_GSC_1.73_R1.5.csv',lab=f'{R}/NW/hypo_qc28/review/eye_class.csv')}
def real(run):
    if run in EXTRA:
        c=pd.read_csv(EXTRA[run]['lab']); return set(c[c['class']=='lfe'].grp)
    if run.endswith('_mute'):
        c=pd.read_csv(f'{R}/{run}/hypo_28/review/class_final.csv'); return set(c[c.cls=='lfe'].grp)
    s=pd.read_csv(f'{R}/{run}/hypo_28/review/auto_score.csv'); col='n_sites_3km' if 'n_sites_3km' in s else 'n_imp_sites'; return set(s[s[col]>=3].grp)
def locs(run):
    h=pd.read_csv(EXTRA[run]['loc'] if run in EXTRA else f'{R}/{run}/hypo/loc_GSC_1.73_R1.5.csv'); return h.set_index('grp')[['lat','lon','dep','erh']]
def stacks(run):
    z=np.load(EXTRA[run]['stk'] if run in EXTRA else f'{R}/{run}/hypo_28/stacks_loc.npz'); S={}
    for k in z.files:
        g,s=k.split('_',1); S.setdefault(int(g),{})[s]=z[k].astype(float)
    return S
def km(a,b,c,d):
    a,b,c,d=map(np.radians,(a,b,c,d)); q=np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2; return 6371*2*np.arcsin(np.sqrt(q))
def dup(A,B):
    lags,ccs=[],[]
    for s in set(A)&set(B):
        best=(0,0)
        for c in (0,1):
            x,y=A[s][c],B[s][c]; n=np.linalg.norm(x)*np.linalg.norm(y)
            if n==0: continue
            cc=correlate(x,y,'full')/n; mid=len(y)-1; w=cc[mid-ML:mid+ML+1]; k=int(np.argmax(w))
            if w[k]>best[0]: best=(w[k],(k-ML)/FS)
        if best[0]>=0.5: ccs.append(best[0]); lags.append(best[1])
    return len(ccs)>=2 and np.ptp(lags)<=0.15, len(ccs), (np.ptp(lags) if len(lags)>1 else np.nan)
import sys
runs=sys.argv[1].split(',') if len(sys.argv)>1 else ['NW_lag15','NW_mute','NW_fair','NW_run1','NW_qc']
TAG=sys.argv[2] if len(sys.argv)>2 else 'NW_union'; import os; os.makedirs(f'{R}/{TAG}',exist_ok=True)
RE={r:real(r) for r in runs}; LO={r:locs(r) for r in runs}; ST={r:stacks(r) for r in runs}
U=[(runs[0],g) for g in sorted(RE[runs[0]])]; rows=[]
for r in runs[1:]:
    for g in sorted(RE[r]):
        la=LO[r].lat.get(g,np.nan); lo=LO[r].lon.get(g,np.nan); isdup=None
        for (ur,ug) in U:
            ula=LO[ur].lat.get(ug,np.nan); ulo=LO[ur].lon.get(ug,np.nan)
            d=km(la,lo,ula,ulo) if np.isfinite([la,lo,ula,ulo]).all() else 0.0
            if d>10: continue
            ok,n,sp=dup(ST[r][g],ST[ur][ug])
            if ok: isdup=(ur,ug,round(d,1),n,sp); break
        rows.append(dict(run=r,grp=g,lat=la,lon=lo,dup_of=f'{isdup[0]}:{isdup[1]}' if isdup else '',dist_km=isdup[2] if isdup else np.nan))
        if not isdup: U.append((r,g))
T=pd.DataFrame(rows); T.to_csv(f'{R}/{TAG}/other_run_families.csv',index=False)
out=[]
for (r,g) in U:
    L=LO[r]; out.append(dict(run=r,grp=g,lat=L.lat.get(g,np.nan),lon=L.lon.get(g,np.nan),depth_km=L.dep.get(g,np.nan),erh=L.erh.get(g,np.nan)))
O=pd.DataFrame(out); O.to_csv(f'{R}/{TAG}/union_families.csv',index=False)
print('real per run:',{r:len(RE[r]) for r in runs})
print('added from other runs:',T[T.dup_of==''].groupby('run').size().to_dict(),' duplicates found:',T[T.dup_of!=''].groupby('run').size().to_dict())
print('UNION:',len(O),' located',O.lat.notna().sum(),' erh<5',(O.erh<5).sum())
