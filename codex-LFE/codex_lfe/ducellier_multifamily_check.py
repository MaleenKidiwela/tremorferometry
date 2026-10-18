"""Reuse one acquired day for seven published families; no new waveform requests."""
from pathlib import Path
import json,argparse
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
from .provenance import write_json,sha256
from .published_pickle import read_published_pickle
from .neural import orientation_matrix
from .waveforms import preprocess_stream
from .matched_filter import normalized_cc,extract_candidates,shifted_sites,robust_scale
from .evaluate import match_events

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',default='results/ducellier_continuous_november_day_v1')
    parser.add_argument('--output',default='results/ducellier_multifamily_november_v1')
    parser.add_argument('--holdout',default='')
    parser.add_argument('--threshold',type=float)
    args=parser.parse_args()
    raw=Path(args.input);out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output must stay inside codex-LFE')
    out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed scan is immutable')
    thresholds=[args.threshold] if args.threshold is not None else [6,8,10,12]
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    source=Path('results/southern_benchmark_inventory_v1')
    bank=json.loads((source/'comparison_bank_design.json').read_text())
    acquisition=json.loads((raw/'acquisition.json').read_text())
    design=json.loads((raw/'design.json').read_text());start=UTCDateTime(design['start']);duration=design['duration_s']
    base=Path('results/ducellier_roster_packets_v1/responses')
    inv=read_inventory(str(base/'BK.xml'))+read_inventory(str(base/'NC.xml'))+read_inventory(str(base/'PB.xml'))
    write_json(out/'design.json',dict(source_sha256=sha256(__file__),acquisition_sha256=sha256(raw/'acquisition.json'),bank_manifest_sha256=sha256(source/'comparison_bank_manifest.json'),threshold_mad=thresholds,withheld_station=args.holdout or None,family_ids=bank['families'],status='development_panel',time_basis='family-specific published template-start clocks; no cross-family origin-time merging',new_waveform_requests=0))
    family_sites={f:{} for f in bank['families']}
    for r in acquisition['records']:
        if r['status']!='received':continue
        station=r['station'];channels=r['channels'];saved=out/(station+'_correlations.npz');receipt=saved.with_suffix('.json')
        if saved.exists() and receipt.exists() and json.loads(receipt.read_text())['sha256']==sha256(saved):
            with np.load(saved) as arrays:
                for family in arrays.files:family_sites[family][station]=arrays[family]
            continue
        if sha256(r['path'])!=r['sha256']:raise ValueError('Raw checksum changed')
        stream=read(r['path'])
        for tr in stream:
            if str(inv.get_response(tr.id,tr.stats.starttime))!=str(inv.get_response(tr.id,tr.stats.endtime)):raise ValueError('Response epoch crossing')
            if inv.get_orientation(tr.id,tr.stats.starttime)!=inv.get_orientation(tr.id,tr.stats.endtime):raise ValueError('Orientation epoch crossing')
            tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80)
        w=preprocess_stream(stream,channels,station,float(start),duration+60,fs=20,band=(1.5,9))
        current=orientation_matrix([inv.get_orientation(c,start) for c in channels]);zne=np.linalg.solve(current,w.x)
        valid=np.broadcast_to(w.valid.all(axis=0),w.x.shape)
        station_arrays={}
        for family in bank['families']:
            tp=source/'source/src/templates_v1'/family
            paths=[tp/f'{station.split(".")[1]}_{c.split(".")[-1]}.pkl' for c in channels]
            if not all(p.exists() for p in paths):continue
            templates=[];angles=[]
            for p in paths:
                tr,angle=read_published_pickle(p)
                if tr.stats.sampling_rate!=20:raise ValueError('Template rate')
                templates.append(tr.data);angles.append(angle)
            x=orientation_matrix(angles)@zne
            cc,good=normalized_cc(x,np.asarray(templates),valid)
            cc=cc[:duration*20]
            if np.isfinite(cc).any():family_sites[family][station]=cc;station_arrays[family]=cc
        np.savez_compressed(saved,**station_arrays);write_json(receipt,dict(sha256=sha256(saved),raw_sha256=r['sha256'],families=list(station_arrays)))
        print(station,'processed',len(station_arrays),'families',flush=True)
    rng=np.random.default_rng(20260914)
    stations=sorted({s for panel in family_sites.values() for s in panel})
    shifts=[{s:int(rng.choice([-1,1])*rng.integers(60,601)*20) for s in stations} for _ in range(8)]
    allrows=[];summary=[]
    for family,all_sites in family_sites.items():
        held=all_sites.get(args.holdout) if args.holdout else None
        if args.holdout and held is None:raise ValueError('Withheld station unavailable for '+family)
        sites={s:x for s,x in all_sites.items() if s!=args.holdout}
        if len(sites)<3:continue
        controls=[shifted_sites(sites,shift) for shift in shifts]
        good=np.isfinite(np.asarray(list(sites.values()))).all(axis=0)
        for panel in controls:good &= np.isfinite(np.asarray(list(panel.values()))).all(axis=0)
        held_panels=[]
        if held is not None:
            held_panels=[held]+[shifted_sites({args.holdout:held},shift)[args.holdout] for shift in shifts]
            for h in held_panels:good &= np.isfinite(h)
        refs=read_published_pickle(source/'source/data/Ducellier/catalogs'/family/'catalog_2004_2011.pkl').copy()
        refs['time_s']=pd.to_datetime(refs[['year','month','day','hour','minute']],utc=True).map(lambda t:t.timestamp())+refs.second
        refs=refs[(refs.time_s>=float(start))&(refs.time_s<float(start)+duration)].copy()
        indices=np.rint((refs.time_s.to_numpy()-float(start))*20).astype(int)
        eligible=np.array([bool(good[i]) if 0<=i<len(good) else False for i in indices],dtype=bool)
        refs=refs.loc[eligible].copy();refs['family_id']=family;refs['time_basis']='published_template_start_'+family
        for threshold in thresholds:
            for i,panel in enumerate([sites]+controls):
                events,mask=extract_candidates(panel,20,float(start),threshold,min_sites=3,common_valid=good)
                rows=[dict(e,family_id=family,time_basis='published_template_start_'+family,threshold_mad=threshold,control=i-1) for e in events]
                if held_panels:
                    h=held_panels[i];scale=robust_scale(h[good])
                    if scale is None:raise ValueError('No withheld normalization exposure')
                    for row in rows:
                        index=int(round((row['time_s']-float(start))*20))
                        row['withheld_station']=args.holdout
                        row['withheld_cc']=float(h[index])
                        row['withheld_score_mad']=float((h[index]-scale[0])/scale[1])
                allrows.extend(rows)
                frame=pd.DataFrame(rows,columns=['time_s','family_id','time_basis'])
                matched=len(match_events(frame,refs,.5)) if i==0 else None
                summary.append(dict(family_id=family,threshold_mad=threshold,control=i-1,candidates=len(events),matched_reference=matched,eligible_reference=len(refs),exposure_s=float(mask.sum()/20),sites=len(sites)))
        print(family,'scored',len(refs),'eligible references',flush=True)
    pd.DataFrame(allrows).to_csv(out/'candidates.csv',index=False);pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_development_panel',production_validated=False,control_shifts_samples=shifts,output_sha256={p.name:sha256(p) for p in [out/'candidates.csv',out/'summary.csv']}))

if __name__=='__main__':main()
