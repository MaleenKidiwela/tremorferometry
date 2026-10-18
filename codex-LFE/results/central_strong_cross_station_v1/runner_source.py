"""Bounded horizontal cross-station coherence diagnostic, not published-method replication."""
from pathlib import Path
from itertools import combinations
import json,argparse
import numpy as np
import pandas as pd
from scipy.fft import rfft,irfft,next_fast_len
from scipy.signal import find_peaks
from obspy.geodetics import gps2dist_azimuth
from .provenance import sha256,write_json


def pair_scores(a,b,fs=20,lag_s=2):
    """a: N x 2 x L; b: N x 2 x (L+2*lag_samples), normalized per overlap.

    Positive lag means the matching segment in b starts later than its central
    segment. Each row searches the four horizontal component pairs equally.
    """
    length=a.shape[-1];pad=round(lag_s*fs)
    if b.shape[-1]!=length+2*pad:raise ValueError('Incorrect lag padding')
    a=a-a.mean(axis=-1,keepdims=True)
    nfft=next_fast_len(a.shape[-1]+b.shape[-1]-1)
    af=rfft(a[...,::-1],n=nfft);bf=rfft(b,n=nfft)
    sums=np.pad(np.cumsum(b,axis=-1),((0,0),(0,0),(1,0)))
    sums2=np.pad(np.cumsum(b*b,axis=-1),((0,0),(0,0),(1,0)))
    benergy=np.maximum(sums2[...,length:]-sums2[...,:-length]-(sums[...,length:]-sums[...,:-length])**2/length,0)
    aenergy=np.sum(a*a,axis=-1)
    best=np.zeros(len(a));lags=np.zeros(len(a));components=np.zeros(len(a),dtype=int)
    for i in range(2):
        for j in range(2):
            cc=irfft(af[:,i]*bf[:,j],n=nfft)[:,length-1:length+2*pad]
            denom=np.sqrt(aenergy[:,i,None]*benergy[:,j])
            cc=np.divide(cc,denom,out=np.zeros_like(cc),where=denom>0)
            cc=np.clip(np.abs(cc),0,1);index=cc.argmax(axis=1);values=cc[np.arange(len(cc)),index];use=values>best
            best[use]=values[use];lags[use]=(index[use]-pad)/fs;components[use]=2*i+j
    return best,lags,components


def station_windows(array,start,fs,times,travel,shift,length_s,pad_s):
    indices=np.rint((times+travel+shift-pad_s-start)*fs).astype(int)
    length=round((length_s+2*pad_s)*fs)
    if indices.min()<0 or indices.max()+length>array.shape[-1]:raise ValueError('Requested windows exceed stored support')
    return np.lib.stride_tricks.sliding_window_view(array[1:],length,axis=-1)[:,indices,:].transpose(1,0,2).astype(np.float64)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',nargs='+',default=['central_activity_picker_panel_v1','central_activity_matched_quiet_v1'])
    parser.add_argument('--selection',default='results/central_tremor_context_v1/active_panel_selection_before_waveforms.json')
    parser.add_argument('--output',default='results/central_cross_station_v1')
    args=parser.parse_args();out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed diagnostic immutable')
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    selection=json.loads(Path(args.selection).read_text())
    design=dict(window_s=8,step_s=1,residual_lag_search_s=2,triangle_closure_tolerance_s=.2,minimum_pair_cc=[.4,.5,.6],refractory_s=8,depth_assumption_km=35,Vs_km_s=3.5,guard_s=240,thresholds_not_calibrated=True,method='Maximum absolute horizontal component-pair CC for all station pairs; triangles require lag closure; same search for shifted controls',purpose='Test cross-station discovery input on existing central panels without neural seeds',limitations=['Not a replication of Savard/Bostock','Source-centered homogeneous moveout is diagnostic','Triangle matches are not LFEs or validated origin times'],source_sha256=sha256(__file__),selection_sha256=sha256('results/central_tremor_context_v1/active_panel_selection_before_waveforms.json'))
    design['selection_sha256']=sha256(args.selection)
    write_json(out/'design_before_computation.json',design)
    summary=[]
    for panel in args.inputs:
        source=Path('results')/panel;dest=out/panel;dest.mkdir(exist_ok=True)
        stations=pd.read_csv(source/'stations.csv').sort_values('id');ids=stations.id.tolist();waves={}
        travel={}
        for row in stations.itertuples():
            with np.load(source/(row.id+'_processed.npz')) as a:
                if not a['valid'].all():raise ValueError('This diagnostic requires continuous valid input')
                waves[row.id]=(a['waveform'].copy(),float(a['start_s']),float(a['fs']))
            distance=gps2dist_azimuth(selection['activity_center_latitude'],selection['activity_center_longitude'],row.latitude,row.longitude)[0]/1000
            travel[row.id]=float(np.hypot(distance,35+row.elevation/1000)/3.5)
        design_source=json.loads((source/'design_before_acquisition.json').read_text());start=pd.Timestamp(design_source['start_utc']).timestamp();times=np.arange(start+240,start+3600-240,1.)
        rng=np.random.default_rng(20260914)
        shifts=[{station:float(rng.uniform(60,180)*rng.choice([-1,1])) for station in ids} for _ in range(8)]
        write_json(dest/'inputs.json',dict(travel_s=travel,control_shifts_s=shifts,source_hashes={p.name:sha256(p) for p in source.glob('*_processed.npz')}))
        for realization in range(9):
            core={};padded={}
            for station,(x,begin,fs) in waves.items():
                shift=0 if realization==0 else shifts[realization-1][station]
                core[station]=station_windows(x,begin,fs,times,travel[station],shift,8,0)
                padded[station]=station_windows(x,begin,fs,times,travel[station],shift,8,2)
            pairs={}
            for a,b in combinations(ids,2):pairs[a,b]=pair_scores(core[a],padded[b])
            score=np.zeros(len(times));triangle=np.full(len(times),'',dtype=object);closure_at_best=np.full(len(times),np.nan)
            for a,b,c in combinations(ids,3):
                ab,ac,bc=pairs[a,b],pairs[a,c],pairs[b,c]
                closure=np.abs(ab[1]+bc[1]-ac[1]);value=np.minimum.reduce([ab[0],ac[0],bc[0]])
                value[closure>.200001]=0;use=value>score;score[use]=value[use];triangle[use]='|'.join([a,b,c]);closure_at_best[use]=closure[use]
            np.savez_compressed(dest/f'r{realization}_scores.npz',time_s=times+4,score=score)
            for threshold in [.4,.5,.6]:
                peaks,_=find_peaks(score,height=threshold,distance=8)
                candidates=pd.DataFrame(dict(window_center_diagnostic_time_s=times[peaks]+4,score=score[peaks],stations=triangle[peaks],closure_s=closure_at_best[peaks],status='unvalidated_cross_station_seed'))
                candidates.to_csv(dest/f'r{realization}_cc{threshold}_seeds.csv',index=False)
                summary.append(dict(panel=panel,realization=realization,threshold=threshold,candidates=len(candidates),nominal_window_start_exposure_s=len(times)))
            print(panel,realization,'complete',flush=True)
        del core,padded,waves
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_cross_station_diagnostic',production_validated=False,phase_accuracy_measured=False))
if __name__=='__main__':main()
