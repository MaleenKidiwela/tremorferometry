import sys, itertools, json, time
sys.path.insert(0,'/home/jovyan/tremorferometry/lfe_array_discovery/scripts')
import numpy as np, pandas as pd
import ab_common as A, ab_auto as AU
OUT='/home/jovyan/tremorferometry/lfe_array_discovery/results/ab_royer/NVI/diag_7of10'
import os; os.makedirs(OUT,exist_ok=True)
stas=A.SUBARRAY['NVI']['stas']
dates=sys.argv[1].split(',')
dc=A.DayCache(stas)
def nms(P, w=30):   # 2-D pair suppression on the union, same neighbourhood as ab_auto (+-15 s in i and j, 0.5 s grid)
    P=P.sort_values('cc',ascending=False).reset_index(drop=True)
    keep=[]; I=P.i.to_numpy(); J=P.j.to_numpy()
    taken=np.zeros(len(P),bool)
    for k in range(len(P)):
        if taken[k]: continue
        keep.append(k)
        taken |= (np.abs(I-I[k])<=w)&(np.abs(J-J[k])<=w)
    return P.iloc[keep]
rows=[]
for d in dates:
    for hh in range(24):
        H=dc.hour(d,hh); H={s:v for s,v in H.items() if v[1].mean()>0}
        ok3=[s for s in stas if s in H and (H[s][0][1:]**2).sum()>0]
        combos=list(itertools.combinations(ok3,7))
        allp=[]; per=[]
        for c in combos:
            R,st=AU.autocorr_hour({s:H[s] for s in c})
            if R is None: continue
            R['combo']=','.join(c); allp.append(R); per.append(len(R))
        U=pd.concat(allp) if allp else pd.DataFrame(columns=['i','j','cc'])
        first=per[0] if per else 0          # combos[0] = first 7 in the listed order = current reading
        un=len(nms(U)) if len(U) else 0
        R10,_=AU.autocorr_hour({s:H[s] for s in ok3})
        rows.append(dict(date=d,hour=hh,n3c=len(ok3),stas=','.join(ok3),n_combos=len(combos),
             first7=first,mean_combo=float(np.mean(per)) if per else 0,max_combo=max(per) if per else 0,
             union_nms=un,all_available=len(R10) if R10 is not None else 0))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(f'{OUT}/per_hour_{"_".join(dates)}.csv',index=False)
