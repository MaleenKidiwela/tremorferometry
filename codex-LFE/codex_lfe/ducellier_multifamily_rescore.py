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
from .response_validity import response_window_mask

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--correlations',required=True)
    parser.add_argument('--input',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    raw=Path(args.input);previous=Path(args.correlations);out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed report is immutable')
    original=json.loads((previous/'design.json').read_text())
    args.holdout=original.get('withheld_station') or ''
    thresholds=original['threshold_mad']
    source=Path('results/southern_benchmark_inventory_v1')
    acquisition=json.loads((raw/'acquisition.json').read_text())
    design=json.loads((raw/'design.json').read_text())
    start=UTCDateTime(design['start']);duration=design['duration_s']
    family_sites={f:{} for f in original['family_ids']}
    grouping=pd.read_csv(source/'comparison_bank_source_families.csv').set_index('template_id').published_family_id.to_dict()
    coverage=[]
    for r in acquisition['records']:
        if r['status']!='received':continue
        station=r['station'];saved=previous/(station+'_correlations.npz')
        receipt=json.loads(saved.with_suffix('.json').read_text())
        if sha256(saved)!=receipt['sha256']:raise ValueError('Correlation checksum changed')
        if sha256(r['path'])!=r['sha256']:raise ValueError('Raw checksum changed')
        stream=read(r['path'],headonly=True)
        with np.load(saved) as arrays:
            for family in arrays.files:
                path=source/'source/src/templates_v1'/family/f'{station.split(".")[1]}_{r["channels"][0].split(".")[-1]}.pkl'
                template,_=read_published_pickle(path)
                cc=arrays[family].copy()
                valid=response_window_mask(stream,r['channels'],float(start),len(cc),20,template.stats.npts)
                before=int(np.isfinite(cc).sum());cc[~valid]=np.nan
                family_sites[family][station]=cc
                coverage.append(dict(station=station,family_id=family,previous_valid_s=before/20,corrected_valid_s=float(np.isfinite(cc).sum()/20)))
    pd.DataFrame(coverage).to_csv(out/'response_edge_coverage.csv',index=False)
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    write_json(out/'design.json',dict(original_design_sha256=sha256(previous/'design.json'),source_sha256=sha256(__file__),response_mask_source_sha256=sha256(Path(__file__).with_name('response_validity.py')),threshold_mad=thresholds,withheld_station=args.holdout or None,edge_policy='Exclude default 5% total cosine taper plus 10 seconds per trace edge, requiring complete template support.',new_waveform_requests=0,new_correlation_evaluations=0))
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
                for row in rows:
                    row['source_template_id']=family
                    row['published_family_id']=grouping[family]
                allrows.extend(rows)
                frame=pd.DataFrame(rows,columns=['time_s','family_id','time_basis'])
                matched=len(match_events(frame,refs,.5)) if i==0 else None
                summary.append(dict(family_id=family,source_template_id=family,published_family_id=grouping[family],threshold_mad=threshold,control=i-1,candidates=len(events),matched_reference=matched,eligible_reference=len(refs),exposure_s=float(mask.sum()/20),sites=len(sites)))
        print(family,'scored',len(refs),'eligible references',flush=True)
    pd.DataFrame(allrows).to_csv(out/'candidates.csv',index=False);pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_development_panel',production_validated=False,control_shifts_samples=shifts,output_sha256={p.name:sha256(p) for p in [out/'candidates.csv',out/'summary.csv']}))

if __name__=='__main__':main()
