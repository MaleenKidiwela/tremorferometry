"""Controlled discovery-seed association on the bounded central recent panel."""
from pathlib import Path
import json,argparse,shutil
import numpy as np
import pandas as pd
import pyocto
from .reuse_association import finish_events,site_groups
from .site_selection import representative_sites
from .provenance import write_json,sha256

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',default='results/central_recent_picker_panel_v1');parser.add_argument('--output',default='results/central_recent_association_v1');parser.add_argument('--selection')
    parser.add_argument('--control-shift-grid-s',type=float,default=0.)
    parser.add_argument('--reuse-observed',help='Prior identical observed associations; only controls are recomputed')
    args=parser.parse_args();source=Path(args.input);out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed association immutable')
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    input_design=json.loads((source/'design_before_acquisition.json').read_text())
    start=pd.Timestamp(input_design['start_utc']).timestamp();end=start+3600
    if args.selection:
        selection=json.loads(Path(args.selection).read_text());lat=selection['activity_center_latitude'];lon=selection['activity_center_longitude']
        lat_bounds=[lat-.6,lat+.6];lon_bounds=[lon-.7,lon+.7]
    else:lat_bounds=[43.7,45.];lon_bounds=[-124.3,-122.3]
    raw=pd.read_parquet(source/'picks.parquet');stations=pd.read_csv(source/'stations.csv')[['id','latitude','longitude','elevation']]
    stations,groups=representative_sites(stations,1.)
    raw=raw[raw.station.isin(stations.id)].copy()
    write_json(out/'physical_site_selection.json',dict(groups=groups,retained=stations.id.tolist(),rule='Lexicographically first station per 1-km proximity group, independent of waveform scores'))
    for station in stations.id:
        with np.load(source/(station+'_processed.npz')) as a:
            if not np.isfinite(a['probabilities'][:,1200:73200]).all():raise ValueError('Panel lacks full core inference coverage')
    rng=np.random.default_rng(20260914)
    if not np.isfinite(args.control_shift_grid_s) or args.control_shift_grid_s<0:raise ValueError('Invalid control shift grid')
    if args.control_shift_grid_s:
        allowed=np.arange(np.ceil(60/args.control_shift_grid_s),np.floor(180/args.control_shift_grid_s)+1)*args.control_shift_grid_s
        if not len(allowed):raise ValueError('No shifts on requested grid')
        shifts=[{s:float(rng.choice(allowed)*rng.choice([-1,1])) for s in sorted(stations.id)} for _ in range(8)]
    else:
        shifts=[{s:float(rng.uniform(60,180)*rng.choice([-1,1])) for s in sorted(stations.id)} for _ in range(8)]
    design=dict(source_sha256=sha256(__file__),picks_sha256=sha256(source/'picks.parquet'),station_sha256=sha256(source/'stations.csv'),thresholds=[.1,.3],minimum_S_sites=[3,4],lat=lat_bounds,lon=lon_bounds,depth_km=[20,50],velocity_model='Homogeneous Vp=6.0 km/s, Vs=3.5 km/s; diagnostic assumption, not a calibrated regional location model',pick_match_tolerance_s=1.2,edt_std_s=.5,core_event_guard_s=240,control_shifts_s=shifts,control_shift_grid_s=args.control_shift_grid_s,full_core_P_S_coverage_verified=True,status='unvalidated_discovery_seed_association')
    if args.reuse_observed:
        prior=Path(args.reuse_observed);old=json.loads((prior/'design.json').read_text())
        for key in ['picks_sha256','station_sha256','thresholds','minimum_S_sites','lat','lon','depth_km','pick_match_tolerance_s','edt_std_s','core_event_guard_s']:
            if old[key]!=design[key]:raise ValueError('Observed reuse mismatch: '+key)
        previous_sites=json.loads((prior/'physical_site_selection.json').read_text())
        if previous_sites['retained']!=stations.id.tolist():raise ValueError('Changed physical-site selection')
        design['reused_observed']={}
        for threshold in [.1,.3]:
            for minimum in [3,4]:
                for suffix in ['.json','_events.parquet','_assignments.parquet']:
                    name=f't{threshold}_n{minimum}_r0'+suffix;path=prior/name
                    shutil.copyfile(path,out/name);design['reused_observed'][str(path)]=sha256(path)
    write_json(out/'design.json',design)
    summaries=[]
    for threshold in [.1,.3]:
        base=raw[raw.model_score>=threshold].rename(columns={'time_s':'time','model_score':'probability'}).reset_index(drop=True)
        base['source_pick_id']=np.arange(len(base))
        for minimum in [3,4]:
            for realization in range(9):
                name=f't{threshold}_n{minimum}_r{realization}';done=out/(name+'.json')
                if done.exists():summaries.append(json.loads(done.read_text()));continue
                picks=base.copy()
                if realization:picks['time']+=picks.station.map(shifts[realization-1])
                picks=picks.sort_values('time').reset_index(drop=True)
                vm=pyocto.VelocityModel0D(6.,3.5,tolerance=2.,association_cutoff_distance=200.)
                assoc=pyocto.OctoAssociator.from_area(lat=lat_bounds,lon=lon_bounds,zlim=(20.,50.),velocity_model=vm,time_before=120.,n_picks=minimum,n_s_picks=minimum,n_p_picks=0,n_p_and_s_picks=0,pick_match_tolerance=1.2,edt_pick_std=.5,exponential_edt=False,min_interevent_time=3.,n_threads=4)
                st=stations.copy();assoc.transform_stations(st)
                events,assignments=assoc.associate(picks,st)
                if len(events):
                    assoc.transform_events(events)
                    assignments['source_pick_id']=picks.iloc[assignments.pick_idx.to_numpy()].source_pick_id.to_numpy()
                events,assignments=finish_events(events,assignments,minimum,start+240,end-240)
                events.to_parquet(out/(name+'_events.parquet'),index=False);assignments.to_parquet(out/(name+'_assignments.parquet'),index=False)
                row=dict(threshold=threshold,minimum_S_sites=minimum,realization=realization,candidates=len(events),exposure_s=3120,near_depth_boundary=int(events.near_depth_boundary.sum()) if len(events) else 0)
                write_json(done,row);summaries.append(row);print(name,len(events),'candidate seeds',flush=True)
    frame=pd.DataFrame(summaries);frame.to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_unvalidated_seed_association',production_validated=False,velocity_model_validated=False,phase_truth_available=False,exposure_s=3120))
    print(frame[frame.realization==0].to_string(index=False))

if __name__=='__main__':main()
