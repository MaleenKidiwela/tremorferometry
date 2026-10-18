#!/usr/bin/env python
"""STAGE 2 -- build a 12-s template for EVERY family, from EACH station's OWN stack.

A station's template is cut from that station's own stack of the family's published detection
times, at the offset tau that stage 1 measured for that (family, station).  A waveform is never
transferred between stations: instrument responses differ.
RANK-DON'T-DISCARD: a template is built for every (family, station) stage 1 tested; the stage-1
evidence columns travel with it so downstream steps can rank, never filter.

CHANNEL COUNT IS PER STATION, NEVER ASSUMED (fix 2026-09-23).  The southern permanent bed is
1-C (30 stations), 3-C (22) and -- where two co-located sensors share the preferred band, as at
BK.YBH (HH1/HH2/HH3 from one sensor plus HHE/HHN/HHZ from the other) -- 6-C (1 station).  The same
station can also change count between days (BK.YBH: 405 six-channel days, 13 three-channel days).
A template is a PER-CHANNEL stack, so every day that builds it must carry the same channel count:
we take the count stage 1 measured tau on (`ncomp` in stage1_existence.csv, its modal count), skip
the non-conforming days, and raise a message naming the station if one ever slips through.

  results/templates_south.npz : tpl {(source,fam,sta): (C,L) float32}, meta table
  python scripts/t9_templates.py [--pre 2.0 --len 12.0]
"""
import os,sys,argparse,warnings
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
def _chans(p):
    try: return [str(c) for c in np.load(p)["chans"]]
    except Exception: return []
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--pre",type=float,default=2.0)
    ap.add_argument("--len",type=float,default=12.0)
    ap.add_argument("--out",default=f"{RES}/aa_shelly/templates6_build.npz")
    ap.add_argument("--stations",default="",help="smoke test: restrict to NET.STA,NET.STA")
    ap.add_argument("--half",default="build",choices=["build","test","all"],
        help="split each family's sampled days by alternating rank so the reproduction score is "
             "measured on days the template never saw")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    E=pd.read_csv(f"{RES}/stage1_existence.csv"); E["fam"]=E.fam.astype(str)
    if a.stations:
        keep_ns=set(a.stations.split(","))
        E=E[(E.net+"."+E.sta).isin(keep_ns)].copy()
        log(f"[templates] restricted to {sorted(keep_ns)} -- {len(E)} (family,station) rows")
    T=pd.read_parquet(f"{RES}/stage1_times.parquet"); T["fam"]=T.fam.astype(str)
    T["day"]=T.t.dt.floor("D")
    if a.half!="all":
        keep=[]; DJ={"build":{},"test":{}}
        for (s_,f_),g in T.groupby(["source","fam"]):
            dd=sorted(g.day.unique())
            DJ["build"][f"{s_}|{f_}"]=[str(pd.Timestamp(x).date()) for x in dd[0::2]]
            DJ["test"][f"{s_}|{f_}"] =[str(pd.Timestamp(x).date()) for x in dd[1::2]]
            sel=set(dd[0::2]) if a.half=="build" else set(dd[1::2])
            keep.append(g[g.day.isin(sel)])
        T=pd.concat(keep,ignore_index=True)
        import json as _j
        # BOTH halves are written here, always: the held-out TEST list is what the Shelly
        # reproduction is scored on, and it must exist even though no template is built from it.
        for h in ("build","test"):
            _j.dump(DJ[h],open(f"{RES}/aa_shelly/stage2_days_{h}.json","w"),indent=0)
    R=pd.read_csv(f"{RES}/stage1_rosters.csv"); R["fam"]=R.fam.astype(str)
    # per-station channel count, exactly the one stage 1 stacked on (constant per station there)
    NCOMP={}
    if "ncomp" in E.columns:
        for (net,sta),g in E.groupby(["net","sta"]):
            v=g.ncomp.dropna()
            if len(v): NCOMP[(net,sta)]=int(pd.Series(v.astype(int)).mode().iloc[0])
    L=int(round(a.len*FS)); rows=[]; TPL={}; nmix=0
    for (net,sta),Es in E.groupby(["net","sta"]):
        sub=T[[(s,f) in set(zip(Es.source,Es.fam)) for s,f in zip(T.source,T.fam)]]
        tt={(r.source,r.fam):r.tt_s for r in Es.itertuples()}
        tau={(r.source,r.fam):r.best_tau for r in Es.itertuples()}
        C=NCOMP.get((net,sta))
        if C is None:                      # fall back to this station's own modal cached count
            ncs=[]
            for day in sub.day.unique():
                p0=npz_path(net,sta,day)
                if os.path.exists(p0):
                    try: ncs.append(int(np.load(p0)["x"].shape[0]))
                    except Exception: pass
            if not ncs:
                log(f"  {net}.{sta}: no cached day readable, skipped"); continue
            C=int(pd.Series(ncs).mode().iloc[0])
        acc={}; n={}; nday=0; nskip=0
        for day,Td in sub.groupby("day"):
            p=npz_path(net,sta,day)
            if not os.path.exists(p): continue
            try: z=np.load(p); X=z["x"].astype(np.float32); t0=float(z["t0"])
            except Exception: continue
            if X.shape[0]!=C:              # e.g. the day one of two co-located sensors is absent
                nskip+=1; continue
            nday+=1
            for (s_,f_),Tf in Td.groupby(["source","fam"]):
                k=(s_,f_)
                if k not in tt: continue
                st0=((Tf.t.values.astype("datetime64[ns]").astype(np.int64)/1e9)+tt[k]+tau[k]-a.pre-t0)
                idx=np.round(st0*FS).astype(int); idx=idx[(idx>=0)&(idx+L<=X.shape[1])]
                if len(idx)==0: continue
                w=X[:,idx[:,None]+np.arange(L)]
                r=np.sqrt((w.astype(np.float64)**2).mean(axis=(0,2)))+1e-20
                s=(w/r[None,:,None]).sum(1)
                if k in acc and acc[k].shape!=s.shape:
                    raise ValueError(
                        f"channel-count mismatch while stacking {net}.{sta} on {pd.Timestamp(day).date()} "
                        f"for {s_}|{f_}: the stack carries {acc[k].shape[0]} channels, this day has "
                        f"{s.shape[0]} (station channel count from stage 1 = {C}; cached channels "
                        f"{_chans(p)}).  A template is a per-channel stack -- a station's channel "
                        f"count must be constant across the days that build it.")
                acc[k]=s if k not in acc else acc[k]+s; n[k]=n.get(k,0)+len(idx)
        nmix+=nskip
        for k,s in acc.items():
            key=f"{k[0]}|{k[1]}|{net}.{sta}"
            TPL[key]=(s/max(n[k],1)).astype(np.float32)
            e=Es[(Es.source==k[0])&(Es.fam==k[1])].iloc[0]
            rows.append(dict(source=k[0],fam=k[1],net=net,sta=sta,netsta=f"{net}.{sta}",key=key,
                ncomp=int(s.shape[0]),ncomp_station=int(C),n_days=nday,n_days_offcount=nskip,
                n_stack=n[k],pre_s=a.pre,len_s=a.len,best_tau=e.best_tau,
                excess_db=e.excess_db,null_pctile=e.null_pctile,z_vs_null=e.z_vs_null,
                hypo_km=e.hypo_km,tt_s=e.tt_s))
        log(f"  {net}.{sta}: {len(acc)} templates, {C}-C, {nday} days"
            +(f", {nskip} days skipped (channel count != {C})" if nskip else ""))
    M=pd.DataFrame(rows); M.to_csv(a.out.replace(".npz",".csv"),index=False)
    np.savez_compressed(a.out,**TPL)
    log(f"[templates] {len(TPL)} templates over {M.groupby('source').fam.nunique().to_dict()} families"
        f"; channel counts {M.ncomp.value_counts().to_dict()}; {nmix} station-days skipped for an "
        f"off-modal channel count")
if __name__=="__main__": main()
