"""Conditional shared-reference comparison using training-derived model clocks.

Clocks transfer published stack labels, not evaluation detections. This is a
development reference-agreement diagnostic, never a precision estimate.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.interpolate import RegularGridInterpolator

from .evaluate import match_events
from .provenance import sha256,write_json
from .reference import read_bostock

BASIS='crosson_model_training_stack_label_transfer_origin'


def in_intervals(times, intervals):
    """Membership in the union; overlapping family intervals are supported."""
    if len(intervals)==0: return np.zeros(len(times),bool)
    intervals=np.asarray(intervals,float)
    order=np.argsort(intervals[:,0]);intervals=intervals[order]
    starts=intervals[:,0];ends=np.maximum.accumulate(intervals[:,1])
    j=np.searchsorted(starts,times,side='right')-1
    return (j>=0)&(times<ends[np.maximum(j,0)])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['config','matcher','association','clocks','output']: ap.add_argument('--'+key,required=True)
    ap.add_argument('--coverage',help='exact corrected neural prediction intervals; require all configured stations to cover predicted S +/-5 s')
    ap.add_argument('--travel-table',default='../lfe_catalog_repair/results/e2_ttable_crosson.npz')
    ap.add_argument('--station-metadata',default='../lfe_catalog_repair/results/e2_stations.csv')
    ap.add_argument('--daily-event-guard-s',type=float,default=0.)
    args=ap.parse_args();out=Path(args.output)
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load(Path(args.config).read_text());mf=Path(args.matcher);nn=Path(args.association)
    clocks={r['family_id']:r['origin_minus_published_s'] for r in
            json.loads((Path(args.clocks)/'family_clocks.json').read_text()) if r['accepted']}
    training=json.loads((mf/'training.json').read_text());anchors=training['anchors_from_published_s']
    clocks={f:v for f,v in clocks.items() if f in anchors}
    if not clocks: raise ValueError('no qualified model clocks')
    reference=read_bostock(cfg['reference_csv'])
    summaries=[];pairs=[];sources={}
    def record(path): sources[str(Path(path).resolve())]=sha256(path)
    for p in [args.config,cfg['reference_csv'],Path(args.clocks)/'family_clocks.json',
              Path(args.clocks)/'manifest.json',mf/'training.json',nn/'manifest.json',Path(__file__)]: record(p)
    if args.coverage:
        coverage=pd.read_parquet(args.coverage)
        station_meta=pd.read_csv(args.station_metadata)
        station_meta.index=station_meta.net+'.'+station_meta.sta
        with np.load(args.travel_table) as table:
            travel=RegularGridInterpolator((table['x'],table['z']),table['S'],bounds_error=True)
        for p in [args.coverage,args.travel_table,args.station_metadata]: record(p)
    for split in ['calibration','test']:
        intervals=[];mf_candidates=[]
        for date in cfg[split+'_dates']:
            ip=mf/f'{split}_{date}_eligible_intervals.parquet';cp=mf/f'{split}_{date}_candidates.parquet'
            record(ip);record(cp)
            intervals.append(pd.read_parquet(ip));mf_candidates.append(pd.read_parquet(cp))
        intervals=pd.concat(intervals,ignore_index=True)
        intervals=intervals[intervals.family_id.isin(clocks)].copy()
        delta=intervals.family_id.map(clocks)-intervals.family_id.map(anchors)
        # Additional 60 s erosion for whole-waveform context, before either arm's scoring.
        intervals['start_s']+=delta+60
        intervals['end_s']+=delta-60
        if args.daily_event_guard_s:
            day=pd.to_datetime(intervals.date,utc=True).dt.as_unit('ns').astype('int64')/1e9
            intervals['start_s']=np.maximum(intervals.start_s,day+args.daily_event_guard_s)
            intervals['end_s']=np.minimum(intervals.end_s,day+86400-args.daily_event_guard_s)
        intervals=intervals[intervals.end_s>intervals.start_s]
        intervals.to_parquet(out/f'{split}_waveform_intervals.parquet',index=False)
        union=intervals[['start_s','end_s']].to_numpy()
        ref=reference[reference.date.isin(cfg[split+'_dates'])&reference.family_id.isin(clocks)].copy()
        ref['time_s']+=ref.family_id.map(clocks);ref['time_basis']=BASIS
        keep=np.zeros(len(ref),bool)
        for family,ii in ref.groupby('family_id').indices.items():
            support=intervals[intervals.family_id==family][['start_s','end_s']].to_numpy()
            keep[ii]=in_intervals(ref.iloc[ii].time_s.to_numpy(),support)
        if args.coverage:
            for station in cfg['stations']:
                m=station_meta.loc[station]
                epi=np.hypot((ref.latitude-m.lat)*111.195,
                    (ref.longitude-m.lon)*111.195*np.cos(np.deg2rad((ref.latitude+m.lat)/2)))
                predicted=ref.time_s.to_numpy()+travel(np.column_stack([epi,ref.depth_km]))
                support=coverage[(coverage.station==station)&(coverage.phase=='S')][['start_s','end_s']].to_numpy().copy()
                support[:,0]+=5.;support[:,1]-=5.
                support=support[support[:,1]>support[:,0]]
                keep &= in_intervals(predicted,support)
        before=len(ref);ref=ref[keep].reset_index(drop=True)
        ref.to_parquet(out/f'{split}_reference.parquet',index=False)
        def evaluate(events,arm,setting,realization):
            events=events.copy()
            events['time_basis']=BASIS
            events=events[in_intervals(events.time_s.to_numpy(),union)].reset_index(drop=True)
            for tolerance in [1.,2.,5.]:
                match=match_events(events,ref,tolerance,same_family=False)
                spatial=[]
                for row in match.itertuples():
                    e,r=events.iloc[row.detection_index],ref.iloc[row.reference_index]
                    distance=float(np.hypot((e.latitude-r.latitude)*111.195,
                        (e.longitude-r.longitude)*111.195*np.cos(np.deg2rad(r.latitude))))
                    spatial.append(distance)
                    if realization=='observed':
                        pairs.append({'split':split,'arm':arm,'setting':setting,'tolerance_s':tolerance,
                            'event_id':str(e.event_id),'reference_id':r.reference_id,'family_id':r.family_id,
                            'time_residual_s':row.residual_s,'horizontal_separation_km':distance})
                summaries.append({'split':split,'arm':arm,'setting':setting,'realization':realization,
                    'tolerance_s':tolerance,'n_clock_qualified_reference_before_masks':before,
                    'n_common_reference':len(ref),'n_candidates_in_common_waveform_intervals':len(events),
                    'n_temporal_matches':len(match),'reference_recovery':len(match)/len(ref) if len(ref) else None,
                    'n_matches_within_10km':int(np.sum(np.asarray(spatial)<=10)),
                    'median_horizontal_separation_km':float(np.median(spatial)) if spatial else None})
        candidates=pd.concat(mf_candidates,ignore_index=True)
        candidates=candidates[candidates.family_id.isin(clocks)].copy()
        candidates['time_s']+=candidates.family_id.map(clocks)-candidates.family_id.map(anchors)
        locations=reference.groupby('family_id')[['latitude','longitude']].first()
        candidates['latitude']=candidates.family_id.map(locations.latitude)
        candidates['longitude']=candidates.family_id.map(locations.longitude)
        candidates['event_id']=candidates.index.astype(str)
        for threshold in [8.,9.,10.]:
            for realization in range(cfg['null_replicates']+1):
                subset=candidates[(candidates.threshold==threshold)&(candidates.realization==realization)]
                evaluate(subset,'matched_filter',f'MAD{threshold:g}',
                         'observed' if realization==0 else f'site_shift_{realization}')
        for threshold,minimum in [(.3,4),(.3,5),(.5,4),(.5,5)]:
            for mode,seed in [('observed',0)]+[(m,s) for m in ['station','station_hour'] for s in [1,2,3]]:
                path=nn/f't{threshold}_n{minimum}_{mode}_{seed}_events.parquet';record(path)
                events=pd.read_parquet(path)
                if len(events):
                    dates=pd.to_datetime(events.time,unit='s',utc=True).dt.strftime('%Y-%m-%d')
                    events=events[dates.isin(cfg[split+'_dates'])].copy()
                    events['time_s']=events.time;events['family_id']=None;events['event_id']=events.idx.astype(str)
                else:
                    events=pd.DataFrame(columns=['time_s','family_id','event_id','latitude','longitude'])
                evaluate(events,'neural_association',f'p{threshold}_n{minimum}',
                         'observed' if mode=='observed' else f'{mode}_{seed}')
    pd.DataFrame(summaries).to_csv(out/'reference_agreement.csv',index=False)
    pd.DataFrame(pairs).to_csv(out/'matched_pairs.csv',index=False)
    write_json(out/'manifest.json',{'sources':sources,'family_clocks':clocks,
        'status':'conditional_development_reference_agreement_not_method_ranking',
        'exact_neural_prediction_coverage_applied':bool(args.coverage),
        'daily_event_guard_s':args.daily_event_guard_s,
        'common_waveform_context_erosion_s':60.,'tolerances_s':[1.,2.,5.],
        'limitations':['Model-dependent origin clock transferred from published templates, not absolute truth',
          'Clock-qualified sources only; discovery outside them not scored',
          'Neural and matcher candidate inventories differ; no precision comparison',
          'Matcher locations inherit template coordinates; zero separation is not independent location validation',
          'Exact neural S coverage applied with +/-5 s guard' if args.coverage else 'Exact legacy inference coverage unavailable',
          'Control methods and exposure after shifts differ; report counts, not equal-FDR claims',
          'No human-adjudicated candidate labels; spatial distances diagnostic only']})
    d=pd.DataFrame(summaries)
    print(d[(d.realization=='observed')&(d.tolerance_s==2)].to_string(index=False))


if __name__=='__main__': main()
