#!/usr/bin/env python
"""STAGE 1 -- the southern EXISTENCE TEST (Merlin's o4 spec, lfe_picker/results/o4_deferred.md),
run under the RANK-DON'T-DISCARD rule (user, 2026-09-23).

For every southern family, stack the modern permanent stations at the family's PUBLISHED detection
times, with a per-family/per-station constant offset tau searched over +-30 s, and compare the
S-window energy excess against the identical search over 200 random-time null stacks.

  statistic(tau) = 10*log10( E(S+0..S+5 s) / E(S-12..S-6 s) )   summed over components
  excess         = max over tau of statistic
  null_pctile    = % of the 200 null maxima below `excess`      <- the MEASUREMENT
  z_vs_null      = (excess - mean(null maxima)) / sd(null maxima)

NOTHING IS FILTERED.  Every family x station tested keeps its row.  Merlin's ">= 2 stations above
the null p99" and ">= 10 of 27 families" are recorded as REFERENCE LINES (columns `above_p99`,
`above_median`) on the distribution, never as gates.  Families are RANKED by evidence and every
family with any station above its own null median is carried to stage 2/3 -- the extension itself
is the test.

Alignment: t_detection + tt_S(station, family) + tau, tt_S from taup/iasp91 at the published
hypocentre (any constant error is absorbed by tau; only the inter-station moveout matters).
Stacks are accumulated twice -- per-window RMS-normalised (primary; every event weighted equally)
and raw linear (secondary) -- and the null gets the identical treatment.

  python scripts/t4_stack_test.py [--nnull 200] [--out results/stage1_existence.csv]
"""
import os,sys,time,argparse,warnings,json
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
PRE=45.0; POST=40.0; TAU=30.0; SWIN=(-1.0,5.0); NWIN=(-12.0,-6.0)
L=int(round((PRE+POST)*FS))
def taup_S(depth_km,epi_km):
    try:
        from obspy.taup import TauPyModel
        global _M
        try: _M
        except NameError: _M=TauPyModel(model="iasp91")
        a=_M.get_travel_times(source_depth_in_km=max(depth_km,0.1),
                              distance_in_degree=epi_km/111.19,
                              phase_list=["s","S","Sn","Sg","S+"])
        if a: return float(min(x.time for x in a))
    except Exception: pass
    return float(np.sqrt(epi_km**2+depth_km**2)/3.5)
def stat_curve(stk,taus_idx,ns,nn0,nn1):
    """stk (C,L) -> statistic at each tau index (vectorised over tau)."""
    e=np.cumsum(np.concatenate([np.zeros((stk.shape[0],1)),stk.astype(np.float64)**2],1),axis=1)
    def E(a,b): return (e[:,b]-e[:,a]).sum(0)
    s0=taus_idx+ns[0]; s1=taus_idx+ns[1]; n0=taus_idx+nn0; n1=taus_idx+nn1
    num=E(s0,s1); den=E(n0,n1)
    return 10*np.log10(np.maximum(num,1e-30)/np.maximum(den,1e-30))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--nnull",type=int,default=200)
    ap.add_argument("--out",default=f"{RES}/stage1_existence.csv")
    ap.add_argument("--rosters",default=f"{RES}/stage1_rosters.csv"); ap.add_argument("--seed",type=int,default=20260923)
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    R=pd.read_csv(a.rosters)
    F=pd.read_csv(f"{RES}/families_south.csv"); F["fam"]=F.fam.astype(str)
    T=pd.read_parquet(f"{RES}/stage1_times.parquet"); T["fam"]=T.fam.astype(str)
    T["day"]=T.t.dt.floor("D")
    R["fam"]=R.fam.astype(str)
    fmeta={(r.source,str(r.fam)):r for r in F.itertuples()}
    # travel times
    R["tt_s"]=[taup_S(fmeta[(r.source,r.fam)].depth,
                      np.sqrt(((r.sta_lat-fmeta[(r.source,r.fam)].lat)*111.19)**2 +
                              ((r.sta_lon-fmeta[(r.source,r.fam)].lon)*111.19*
                               np.cos(np.radians(fmeta[(r.source,r.fam)].lat)))**2))
             for r in R.itertuples()]
    taus=np.arange(-TAU,TAU+1e-9,1/FS)
    tidx=np.round((PRE+taus)*FS).astype(int)
    ns=(int(round(SWIN[0]*FS)),int(round(SWIN[1]*FS)))
    nn0,nn1=int(round(NWIN[0]*FS)),int(round(NWIN[1]*FS))
    rng=np.random.default_rng(a.seed)
    rows=[]; t00=time.time()
    for (net,sta),Rs in R.groupby(["net","sta"]):
        keys=[(r.source,r.fam) for r in Rs.itertuples()]
        tt={(r.source,r.fam):r.tt_s for r in Rs.itertuples()}
        sub=T[[ (s,f) in set(keys) for s,f in zip(T.source,T.fam)]]
        if len(sub)==0: continue
        # modal component count for this station over its sampled days (1-C stations are normal
        # in the southern permanent network before ~2012)
        ncs=[]
        for day in sub.day.unique():
            p0=npz_path(net,sta,day)
            if os.path.exists(p0):
                try: ncs.append(int(np.load(p0)["x"].shape[0]))
                except Exception: pass
        if not ncs: continue
        C=int(pd.Series(ncs).mode().iloc[0])
        acc={k:None for k in keys}          # (real_norm, real_raw, null_norm, null_raw, n)
        for day,Td in sub.groupby("day"):
            p=npz_path(net,sta,day)
            if not os.path.exists(p): continue
            try: z=np.load(p); X=z["x"].astype(np.float32); t0=float(z["t0"])
            except Exception: continue
            if X.shape[1]<0.5*86400*FS: continue
            if X.shape[0]!=C: continue
            if not np.any(X): continue
            for (s_,f_),Tf in Td.groupby(["source","fam"]):
                k=(s_,f_)
                if k not in acc: continue
                if acc[k] is None:
                    acc[k]=[np.zeros((C,L)),np.zeros((C,L)),
                            np.zeros((a.nnull,C,L)),np.zeros((a.nnull,C,L)),0]
                A=acc[k]
                st0=((Tf.t.values.astype("datetime64[ns]").astype(np.int64)/1e9)+tt[k]-t0-PRE)
                idx=np.round(st0*FS).astype(int)
                idx=idx[(idx>=0)&(idx+L<=X.shape[1])]
                if len(idx)==0: continue
                w=X[:,idx[:,None]+np.arange(L)]                       # (C,n,L)
                A[1]+=w.sum(1)
                r=np.sqrt((w.astype(np.float64)**2).mean(axis=(0,2)))+1e-20
                A[0]+=(w/r[None,:,None]).sum(1)
                A[4]+=len(idx)
                hi=X.shape[1]-L
                for j in range(a.nnull):
                    ridx=rng.integers(0,hi,size=len(idx))
                    wn=X[:,ridx[:,None]+np.arange(L)]
                    A[3][j]+=wn.sum(1)
                    rn=np.sqrt((wn.astype(np.float64)**2).mean(axis=(0,2)))+1e-20
                    A[2][j]+=(wn/rn[None,:,None]).sum(1)
            del X
        for k,A in acc.items():
            if A is None or A[4]<20: continue
            out={}
            for tag,ireal,inull in (("norm",0,2),("raw",1,3)):
                c=stat_curve(A[ireal],tidx,ns,nn0,nn1)
                ex=float(c.max()); bt=float(taus[int(np.argmax(c))])
                nm=np.array([stat_curve(A[inull][j],tidx,ns,nn0,nn1).max() for j in range(a.nnull)])
                out[tag]=dict(excess=ex,tau=bt,p50=float(np.median(nm)),p99=float(np.percentile(nm,99)),
                              mean=float(nm.mean()),sd=float(nm.std()),
                              pctile=100.0*float((nm<ex).mean()),
                              z=float((ex-nm.mean())/(nm.std()+1e-12)))
            o=out["norm"]; r=out["raw"]
            rows.append(dict(source=k[0],fam=k[1],net=net,sta=sta,ncomp=C,n_stack=A[4],
                tt_s=round(tt[k],2),hypo_km=float(Rs[(Rs.source==k[0])&(Rs.fam==k[1])].hypo_km.iloc[0]),
                excess_db=o["excess"],best_tau=o["tau"],null_pctile=o["pctile"],z_vs_null=o["z"],
                null_median=o["p50"],null_p99=o["p99"],null_mean=o["mean"],null_sd=o["sd"],
                above_p99=int(o["excess"]>o["p99"]),above_median=int(o["excess"]>o["p50"]),
                excess_db_raw=r["excess"],best_tau_raw=r["tau"],null_pctile_raw=r["pctile"],
                z_vs_null_raw=r["z"],above_p99_raw=int(r["excess"]>r["p99"]),
                above_median_raw=int(r["excess"]>r["p50"])))
        log(f"  {net}.{sta}: {len([1 for k,A in acc.items() if A is not None])} families, "
            f"{time.time()-t00:.0f}s")
    D=pd.DataFrame(rows); D.to_csv(a.out,index=False)
    if len(D):
        fam=D.groupby(["source","fam"]).agg(
              n_sta_tested=("sta","size"),n_stack=("n_stack","max"),
              max_excess_db=("excess_db","max"),max_z=("z_vs_null","max"),
              max_pctile=("null_pctile","max"),median_pctile=("null_pctile","median"),
              n_above_p99=("above_p99","sum"),n_above_median=("above_median","sum"),
              n_above_p99_raw=("above_p99_raw","sum"),
              min_hypo_km=("hypo_km","min"),med_hypo_km=("hypo_km","median")).reset_index()
        # EVIDENCE RANK -- an ordering, not a filter.  Every family keeps its row.
        fam["evidence"]=fam.n_above_p99*100+fam.max_z
        fam=fam.sort_values(["source","evidence"],ascending=[True,False])
        fam["rank_in_source"]=fam.groupby("source").cumcount()+1
        # carried forward = ANY station above its own null MEDIAN (user rule 4); flagged, never dropped
        fam["carry_stage2"]=(fam.n_above_median>=1).astype(int)
        fam.to_csv(a.out.replace(".csv","_family.csv"),index=False)
        for src,g in fam.groupby("source"):
            log(f"STAGE1 {src}: {len(g)} families tested, all retained | "
                f"median max_z {g.max_z.median():.2f}, median max_pctile {g.max_pctile.median():.1f} | "
                f"REFERENCE LINES: >=2 stations above null p99 -> {int((g.n_above_p99>=2).sum())}, "
                f">=1 station above null p99 -> {int((g.n_above_p99>=1).sum())}, "
                f"carried to stage 2 -> {int(g.carry_stage2.sum())}")
if __name__=="__main__": main()
