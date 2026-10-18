import numpy as np, pandas as pd, sys
from scipy.signal import hilbert
sys.path.insert(0,'scripts'); import ab_common as A
SITE={'B013':'S1','W030':'S1','B943':'S1'}
def score(npz):
    z=np.load(npz); fam={}
    for k in z.files: g,s=k.split('_',1); fam.setdefault(int(g),[]).append(s)
    rows=[]
    for g,ss in fam.items():
        imp={}; tpk={}
        for s in ss:
            x=z[f'{g}_{s}'].astype(float); e=np.abs(hilbert(x,axis=1)).sum(0)
            k=int(np.argmax(e)); w=e[max(k-20,0):k+21].mean()
            imp[s]=e.max()/np.median(e); tpk[s]=k/40.0
            # duration: time above half max
            imp[s+'_dur']=(e>e.max()/2).sum()/40.0
        good=[s for s in ss if imp[s]>5 and imp[s+'_dur']<3]
        sites={SITE.get(s,s) for s in good}
        t=np.array([tpk[s] for s in good]) if good else np.array([])
        rows.append(dict(grp=g,n_imp_sites=len(sites),t_spread=(np.ptp(t) if len(t)>1 else np.nan),
                         med_imp=np.median([imp[s] for s in ss])))
    return pd.DataFrame(rows)
