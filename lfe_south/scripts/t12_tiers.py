#!/usr/bin/env python
"""STAGE 3b -- FA-ONLY tiers per network epoch, exactly the frozen rule of
lfe_extension/scripts/c_tiers_by_epoch.py, re-measured on the southern data.

  T_core = smallest T on a 0.25-MAD grid at or above 8.0 with quiet-day FA <= 0.1 / template-day
  T_ext  = smallest T with quiet-day FA <= 1.0 AND post-SSE FA <= 1.5
  quiet family-day  = observable, no PNSN tremor within BOX deg of the family that day, and >= 3
                      days from any such tremor day
  post-SSE day      = 1..10 days after the last tremor day of an episode (episode ends after >= 3
                      tremor-free days)
Thresholds are set by false alarms ONLY; recall is never used to pick T.  No family is dropped --
a family whose epoch has too few quiet days falls back to the epoch median T and is flagged.
"""
import os,sys,argparse,glob
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
BOX=0.3; GRID=np.arange(8.0,14.01,0.25)
def episodes(trem_days):
    td=sorted(pd.to_datetime(sorted(set(trem_days))))
    if not td: return set()
    ends=[];cur=td[0]
    for i in range(1,len(td)):
        if (td[i]-td[i-1]).days>=3: ends.append(td[i-1])
        cur=td[i]
    ends.append(td[-1])
    out=set()
    for e in ends:
        for k in range(1,11): out.add(e+pd.Timedelta(days=k))
    return out
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True)
    ap.add_argument("--det",required=True,help="glob of monthly detection parquets")
    ap.add_argument("--out",default=None)
    a=ap.parse_args()
    out=a.out or f"{RES}/tiers_{a.source}.csv"
    D=pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(a.det))],ignore_index=True)
    D["fam"]=D.fam.astype(str); D["day"]=D.t.dt.floor("D")
    COV=pd.concat([pd.read_csv(p) for p in sorted(glob.glob(a.det.replace(".parquet","_coverage.csv")))],
                  ignore_index=True)
    COV["day"]=pd.to_datetime(COV.day)
    COV=COV.sort_values("day").drop_duplicates("day")
    F=pd.read_csv(f"{RES}/families_south.csv"); F=F[F.source==a.source]; F["fam"]=F.fam.astype(str)
    TRE=pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv")
    TRE["t"]=pd.to_datetime(TRE.time).astype("datetime64[ns]"); TRE["day"]=TRE.t.dt.floor("D")
    M=pd.read_csv(f"{RES}/templates_south.csv"); M=M[M.source==a.source]
    ntpl=M.groupby(M.fam.astype(str)).size().to_dict()
    # network epochs = maximal runs with the same number of contributing stations
    COV["epoch"]=(COV.n_sta!=COV.n_sta.shift()).cumsum()
    ep=COV.groupby("epoch").agg(t0=("day","min"),t1=("day","max"),n_sta=("n_sta","first"),
                                n_days=("day","size")).reset_index()
    ep=ep[ep.n_days>=15].reset_index(drop=True)
    rows=[]
    for f in F.fam:
        fl=F[F.fam==f].iloc[0]
        tf=TRE[(TRE.lat.between(fl.lat-BOX,fl.lat+BOX))&(TRE.lon.between(fl.lon-BOX,fl.lon+BOX))]
        tdays=set(tf.day)
        post=episodes(tdays)
        far=set()
        for d in COV.day:
            if d in tdays: continue
            if min([abs((d-x).days) for x in tdays],default=999)>=3: far.add(d)
        Df=D[D.fam==f]
        nt=max(ntpl.get(f,1),1)
        for e in ep.itertuples():
            days=COV[(COV.day>=e.t0)&(COV.day<=e.t1)].day
            q=[d for d in days if d in far]; p=[d for d in days if d in post]
            de=Df[(Df.day>=e.t0)&(Df.day<=e.t1)]
            fa_q=lambda T: (de[(de.day.isin(q))&(de.net_mad>=T)].shape[0]/max(len(q),1))
            fa_p=lambda T: (de[(de.day.isin(p))&(de.net_mad>=T)].shape[0]/max(len(p),1))
            Tc=next((T for T in GRID if fa_q(T)<=0.1),np.nan)
            Te=next((T for T in GRID if fa_q(T)<=1.0 and fa_p(T)<=1.5),np.nan)
            rows.append(dict(source=a.source,fam=f,epoch=int(e.epoch),t0=str(e.t0.date()),
                t1=str(e.t1.date()),n_sta=int(e.n_sta),n_quiet_days=len(q),n_post_days=len(p),
                n_det=len(de),T_core=Tc,T_ext=Te,
                fa_quiet_at_8=fa_q(8.0),fa_post_at_8=fa_p(8.0),n_templates=nt,
                derived=int(len(q)>=10)))
    R=pd.DataFrame(rows)
    med=R.groupby("epoch")[["T_core","T_ext"]].median()
    R["T_core"]=R.apply(lambda r: r.T_core if np.isfinite(r.T_core) else med.T_core.get(r.epoch,np.nan),1)
    R["T_ext"]=R.apply(lambda r: r.T_ext if np.isfinite(r.T_ext) else med.T_ext.get(r.epoch,np.nan),1)
    R.to_csv(out,index=False)
    # per-family per-year catalogue
    key=R.set_index(["fam","epoch"])
    D["epoch"]=pd.cut(D.day,bins=[pd.Timestamp("1900-01-01")]+list(ep.t1)+[pd.Timestamp("2100-01-01")],
                      labels=False,include_lowest=True)
    def thr(r,col):
        try: return key.loc[(r.fam,int(ep.epoch.iloc[min(int(r.epoch),len(ep)-1)])),col]
        except Exception: return np.nan
    D["T_core"]=[thr(r,"T_core") for r in D.itertuples()]
    D["T_ext"]=[thr(r,"T_ext") for r in D.itertuples()]
    D["core"]=(D.net_mad>=D.T_core).astype(int); D["ext"]=(D.net_mad>=D.T_ext).astype(int)
    D["year"]=D.t.dt.year
    D.to_parquet(f"{RES}/catalog_{a.source}.parquet")
    py=D.groupby(["fam","year"]).agg(n_all=("t","size"),n_core=("core","sum"),n_ext=("ext","sum")).reset_index()
    # every family keeps a row for every year, even with zero detections
    full=pd.MultiIndex.from_product([sorted(F.fam),sorted(D.year.unique())],names=["fam","year"])
    py=py.set_index(["fam","year"]).reindex(full,fill_value=0).reset_index()
    py.to_csv(f"{RES}/catalog_{a.source}_per_family_year.csv",index=False)
    log(f"[tiers {a.source}] epochs {len(ep)}, families {R.fam.nunique()}, "
        f"median T_core {R.T_core.median():.2f} T_ext {R.T_ext.median():.2f}; "
        f"catalogue {int(D.core.sum())} core / {int(D.ext.sum())} ext of {len(D)} detections")
if __name__=="__main__": main()
