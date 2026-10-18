"""Transfer published stack S labels to training stacks by waveform similarity.

Exploratory, model-dependent timing evidence. No single-event phase truth is
created, and absolute origin-time accuracy is not established.
"""
import argparse
import json
from pathlib import Path
import re

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey
import pandas as pd
import yaml

from .matched_filter import normalized_cc
from .neural import orientation_matrix
from .provenance import sha256, write_json
from .reference import read_bostock


def response_correct(stack, channels, inventory, dates, fs):
    from obspy import UTCDateTime
    n = stack.shape[1]
    length = 1 << (2*n-1).bit_length()
    frequency = np.fft.rfftfreq(length,1/fs)
    band = (frequency>=.7)&(frequency<=9.)
    taper = np.zeros_like(frequency)
    taper[band] = 1.
    a=(frequency>=.7)&(frequency<1)
    taper[a]=.5-.5*np.cos(np.pi*(frequency[a]-.7)/.3)
    a=(frequency>8)&(frequency<=9)
    taper[a]=.5+.5*np.cos(np.pi*(frequency[a]-8))
    rows,angles=[],[]
    for channel,x in zip(channels,stack):
        epoch=UTCDateTime(dates[0])
        response=inventory.get_response(channel,epoch)
        angle=inventory.get_orientation(channel,epoch)
        for date in dates:
            for t in [UTCDateTime(date),UTCDateTime(date)+86399]:
                if str(inventory.get_response(channel,t))!=str(response) or inventory.get_orientation(channel,t)!=angle:
                    raise ValueError('response/orientation changes across training dates')
        transfer=response.get_evalresp_response_for_frequencies(frequency[band],output='VEL')
        if not np.isfinite(transfer).all() or (np.abs(transfer)==0).any():
            raise ValueError('invalid in-band response')
        spectrum=np.fft.rfft((x-x.mean())*tukey(n,.1),n=length)
        output=np.zeros_like(spectrum)
        output[band]=spectrum[band]/transfer*taper[band]
        rows.append(np.fft.irfft(output,n=length)[:n]);angles.append(angle)
    return np.linalg.solve(orientation_matrix(angles),np.asarray(rows))


def best_s_match(stack, template, fs=20., pre_phase_s=2.):
    cc,_=normalized_cc(stack,template)
    # Our stacks start at published time -10 s; template starts 2 s before S.
    s_times=-10+np.arange(len(cc))/fs+pre_phase_s
    allowed=(s_times>=0)&(s_times<=50)&np.isfinite(cc)
    if not allowed.any():
        raise ValueError('no finite correlation')
    index=int(np.argmax(np.where(allowed,cc,-np.inf)))
    alternatives=allowed&(np.abs(s_times-s_times[index])>=2)
    second=float(np.max(cc[alternatives])) if alternatives.any() else None
    return float(s_times[index]),float(cc[index]),second


def phase_site_clocks(group):
    """P and S from one site get one vote, and contradictory phases are rejected."""
    passing=group[group.accepted_waveform_transfer]
    by_site=passing.groupby('station').origin_minus_published_s.agg(['median','min','max'])
    return by_site[by_site['max']-by_site['min']<=1.]['median']


def interval_source_identity(source, families):
    """Unique compatibility with printed two-decimal source coordinates.

    Include either formatting outcome at decimal ties. This is inferred identity
    evidence, not an authoritative mapping; never choose a nearest nonmatch.
    """
    mask=((families.source_latitude_deg-source.latitude).abs()<=.00500001) & \
         ((families.source_longitude_deg-source.longitude).abs()<=.00500001) & \
         ((families.source_depth_km-source.depth_km).abs()<=.00500001)
    ids=families.index[mask]
    return ids[0] if len(ids)==1 else None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True);ap.add_argument('--stacks',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--phases',default='S',choices=['S','PS'])
    ap.add_argument('--source-identity',choices=['legacy_nearest','rounding_interval'],default='legacy_nearest')
    ap.add_argument('--require-date-stability',action='store_true')
    ap.add_argument('--published',default='../lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia')
    ap.add_argument('--inventory',default='results/map_network_metadata_2009_2026/CN.xml')
    ap.add_argument('--travel-table',default='../lfe_catalog_repair/results/e2_ttable_crosson.npz')
    args=ap.parse_args();out=Path(args.output);root=Path(args.stacks);published=Path(args.published)
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    from obspy import read_inventory
    cfg=yaml.safe_load(Path(args.config).read_text());fs=float(cfg['fs'])
    meta=pd.read_csv(published/'metadata.csv'); inv=read_inventory(args.inventory)
    ref=read_bostock(cfg['reference_csv']).groupby('family_id').first()
    families=meta.groupby('source_family_id')[['source_latitude_deg','source_longitude_deg','source_depth_km']].first()
    counts=pd.DataFrame(json.loads((root/'counts.json').read_text())).set_index(['family_id','station','half']).n_windows
    archive=np.load(root/'stacks.npz',allow_pickle=False)
    with np.load(args.travel_table) as v:
        travel={phase:RegularGridInterpolator((v['x'],v['z']),v[phase],bounds_error=True) for phase in args.phases}
    rows=[];mapping=[]
    with h5py.File(published/'waveforms.hdf5','r') as h5:
        for family in sorted({key.split('__')[0] for key in archive.files}):
            source=ref.loc[family]
            distance=np.sqrt(((families.source_latitude_deg-source.latitude)*111.195)**2+
                ((families.source_longitude_deg-source.longitude)*111.195*np.cos(np.deg2rad(source.latitude)))**2+
                (families.source_depth_km-source.depth_km)**2)
            rank=distance.sort_values();source_id=rank.index[0]
            pub_source=families.loc[source_id]
            reverse=np.sqrt(((ref.latitude-pub_source.source_latitude_deg)*111.195)**2+
                ((ref.longitude-pub_source.source_longitude_deg)*111.195*np.cos(np.deg2rad(source.latitude)))**2+
                (ref.depth_km-pub_source.source_depth_km)**2)
            reciprocal=reverse.idxmin()==family
            accepted_map=bool(rank.iloc[0]<1 and rank.iloc[1]-rank.iloc[0]>.25 and reciprocal)
            if args.source_identity=='rounding_interval':
                identity=interval_source_identity(source,families)
                accepted_map=identity is not None
                if accepted_map:
                    source_id=identity
                    pub_source=families.loc[source_id]
            mapping.append({'family_id':family,'published_source_id':float(source_id),
                'nearest_distance_km':float(rank.iloc[0]),'next_distance_km':float(rank.iloc[1]),
                'reciprocal_nearest':bool(reciprocal),'accepted_mapping':accepted_map,
                'identity_method':args.source_identity})
            if not accepted_map: continue
            for station in ['CN.PGC','CN.SNB','CN.VGZ']:
                candidates=meta[(meta.source_family_id==source_id)&(meta.station_code==station.split('.')[1])]
                if len(candidates)!=1: continue
                row=candidates.iloc[0]
                if row.trace_sampling_rate_hz!=fs: raise ValueError('unexpected published sampling rate')
                # PO is the archive's generic network label; resolve by code AND coordinates.
                epoch=__import__('obspy').UTCDateTime(cfg['train_dates'][0])
                coordinates=inv.get_coordinates(cfg['stations'][station]['channels'][0],epoch)
                station_error=np.hypot((coordinates['latitude']-row.station_latitude_deg)*111.195,
                    (coordinates['longitude']-row.station_longitude_deg)*111.195*np.cos(np.deg2rad(row.station_latitude_deg)))
                if station_error>.1: raise ValueError('published station identity/coordinates mismatch')
                match=re.fullmatch(r'(\w+)\$(\d+),:3,:(\d+)',row.trace_name)
                if not match: raise ValueError('unsupported published trace key')
                pub=h5['data'][match[1]][int(match[2]),:3,:int(match[3])]
                pub=sosfiltfilt(butter(4,[1,8],fs=fs,btype='bandpass',output='sos'),pub,axis=1)
                stacks={h:archive[f'{family}__{station}__{h}'] for h in ['early','late']
                        if f'{family}__{station}__{h}' in archive}
                if not stacks: continue
                total=sum(counts.loc[family,station,h] for h in stacks)
                if total<20: continue
                merged=sum(x*counts.loc[family,station,h] for h,x in stacks.items())/total
                aligned={h:response_correct(x,cfg['stations'][station]['channels'],inv,cfg['train_dates'],fs)
                         for h,x in {**stacks,'all':merged}.items()}
                epi=np.hypot((source.latitude-row.station_latitude_deg)*111.195,
                    (source.longitude-row.station_longitude_deg)*111.195*np.cos(np.deg2rad((source.latitude+row.station_latitude_deg)/2)))
                for phase in args.phases:
                    label=row[f'trace_{phase}_arrival_sample']
                    if pd.isna(label): continue  # Never substitute predicted phase labels.
                    before,duration,nchannels=(1.,4.,1) if phase=='P' else (2.,8.,3)
                    start=int(round(label-before*fs));template=pub[:nchannels,start:start+int(duration*fs)]
                    if start<0 or template.shape!=(nchannels,int(duration*fs)): continue
                    matches={h:best_s_match(x[:nchannels],template,fs,before) for h,x in aligned.items()}
                    time_s,cc,second=matches['all']
                    check_dates=all(h in stacks and counts.loc[family,station,h]>=10 for h in ['early','late'])
                    drift=abs(matches['early'][0]-matches['late'][0]) if check_dates else None
                    accepted=cc>=.6 and cc-second>=.1 and (not check_dates or drift<=.5)
                    if args.require_date_stability:accepted=accepted and check_dates
                    tt=float(travel[phase]([[epi,source.depth_km]])[0])
                    rows.append({'family_id':family,'published_source_id':float(source_id),'station':station,'phase':phase,
                        'n_training_windows':int(total),'correlation':cc,'alternative_peak_cc':second,
                        'phase_arrival_after_published_s':time_s,'model_travel_s':tt,'origin_minus_published_s':time_s-tt,
                        'date_stability_measurable':check_dates,'date_half_drift_s':drift,
                        'accepted_waveform_transfer':bool(accepted)})
    frame=pd.DataFrame(rows);frame.to_csv(out/'station_transfers.csv',index=False)
    write_json(out/'source_mapping.json',mapping)
    clocks=[]
    if len(frame):
        for family,g in frame.groupby('family_id'):
            v=phase_site_clocks(g)
            spread=float(v.max()-v.min()) if len(v) else None
            accepted=len(v)>=2 and spread<=1.
            clocks.append({'family_id':family,'n_independent_sites':len(v),'clock_spread_s':spread,
                'accepted':bool(accepted),'origin_minus_published_s':float(v.median()) if accepted else None})
    write_json(out/'family_clocks.json',clocks)
    files=[Path(args.config),root/'manifest.json',root/'counts.json',root/'stacks.npz',published/'metadata.csv',
        published/'waveforms.hdf5',Path(args.inventory),Path(args.travel_table),Path(cfg['reference_csv']),Path(__file__)]
    write_json(out/'manifest.json',{'sources':{str(p.resolve()):sha256(p) for p in files},
        'phases':args.phases,'P_window':'vertical only, -1 to +3 s from published observed P label',
        'source_identity_method':args.source_identity,'require_date_stability':args.require_date_stability,
        'S_window':'three components, -2 to +6 s from published observed S label',
        'source_component_convention':'published Z12 treated as ZNE following prior E4 polarization audit; raw sensors use StationXML rotation',
        'criteria':{'max_source_join_km':1.,'min_next_source_margin_km':.25,'min_correlation':.6,
                    'min_peak_advantage':.1,'max_measurable_date_drift_s':.5,'min_sites':2,'max_clock_spread_s':1.},
        'status':'exploratory_training_stack_label_transfer','absolute_origin_time_validated':False,
        'limitations':['Source IDs joined by coordinates and waveform consistency, not authoritative identity table',
                       'Published template construction may include these years; this is clock calibration, not independent picker evaluation',
                       'Model-dependent clock and search-selected correlation; no single-event phase labels',
                       'Temporal stability not measurable for every station/source pair']})
    print(frame.to_string(index=False),flush=True)
    print('Family clocks:',json.dumps(clocks),flush=True)


if __name__=='__main__': main()
