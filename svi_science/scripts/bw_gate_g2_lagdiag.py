"""bw_gate_g2_lagdiag -- G2 diagnostic: per-event best lag (+-10 s) of the 2012 PGC template against PGC data around
Bostock event times (S-1..S+7 s window), per era; a peak away from 0 = a time-base shift of the catalogue (or data).
  python bw_gate_g2_lagdiag.py 2004-07,2008-05 [STA]
"""
import bw_gate_g2_epoch as E, numpy as np, pandas as pd
S=E.setup(); T=S['T']; V=S['V']; A4=S['A4']
Z=np.load(f"{E.G.REP}/results/e5r_templates_2012.npz", allow_pickle=True)
QC=S['QC']; QC=QC[(QC.sta=='PGC')&QC.usable&QC.template].set_index('bkey')
import sys
for win in sys.argv[1].split(','):
    year=int(E.WIN[win][0][:4]); T.CACHES=E.CACHES[year]
    days=T.load_days('CN.PGC',year); Ev=E.events(win)
    lags=[]
    fams=Ev.fam.value_counts().index[:12]
    for f in fams:
        pass
        if f not in QC.index: continue
        fr=S['F'].loc[f]; dx=T.epi_km(fr.lat,fr.lon,float(S['stab'].loc['PGC','lat']),float(S['stab'].loc['PGC','lon']))
        ttS=float(V.tt(dx,fr.dep,'S'))
        ref=np.asarray(Z[f"PGC|{int(f)}"],float); i_ref=400+int(round(QC.loc[f,'s_pick_s']*20))-20
        tpl=ref[:, i_ref:i_ref+160]; tpl=tpl/np.linalg.norm(tpl)
        g=Ev[Ev.fam==f].sort_values('mag',ascending=False).head(150)
        Sw,_=T.cut_windows(days,g.t_origin.values+ttS-20,800)
        if Sw is None: continue
        for s in Sw:
            # template expected at index 400-20=380 (S-1s); scan +-10 s
            best=(-9,0)
            for k in range(-200,201,1):
                j=380+k; seg=s[:, j:j+160]; n=np.linalg.norm(seg)
                if n>0:
                    c=(seg*tpl).sum()/n
                    if c>best[0]: best=(c,k)
            if best[0]>0.3: lags.append(best[1]/20)
    lags=np.array(lags)
    h=np.histogram(lags,bins=np.arange(-2,2.01,0.25))[0]; c=lags[np.abs(lags)<2]
    print(win,'core median',np.median(c),'core mean',c.mean().round(3))
    print(win, len(lags), 'median lag', np.median(lags) if len(lags) else None, 'frac |lag|<0.5', np.mean(np.abs(lags)<0.5) if len(lags) else None, h)
