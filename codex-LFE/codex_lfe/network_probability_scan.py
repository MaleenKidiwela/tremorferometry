"""Exploratory continuous-score network scan with grid-preserving site controls.

No station score is treated as a calibrated event probability. Source coordinates
use a homogeneous diagnostic model and are not validated earthquake locations.
"""
from pathlib import Path
import json,os
import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import find_peaks
from .matched_filter import robust_scale
from .site_selection import representative_sites
from .provenance import write_json,sha256


def score_grid(scores,offsets,indices,shifts,device):
    import torch
    x=torch.as_tensor(scores,dtype=torch.float32,device=device)
    base=torch.as_tensor(indices,dtype=torch.int64,device=device)
    offset=torch.as_tensor(offsets-shifts[None,:],dtype=torch.int64,device=device)
    best_all=torch.full((len(indices),),-torch.inf,device=device)
    best_gate=best_all.clone();best_node=torch.full((len(indices),),-1,dtype=torch.int64,device=device)
    station=torch.arange(scores.shape[0],device=device)[None,:,None]
    for first in range(0,len(offset),32):
        ix=base[None,None,:]+offset[first:first+32,:,None]
        values=x[station,ix]
        mean=values.clamp(-6,6).mean(dim=1)
        best_all=torch.maximum(best_all,mean.max(dim=0).values)
        gated=mean.masked_fill((values>=2).sum(dim=1)<4,-torch.inf)
        value,node=gated.max(dim=0);use=value>best_gate
        best_gate=torch.maximum(best_gate,value);best_node[use]=node[use]+first
    return best_all.cpu().numpy(),best_gate.cpu().numpy(),best_node.cpu().numpy()


def main():
    import torch
    root=Path(__file__).resolve().parents[1];out=root/'results/north_network_probability_scan_v1';out.mkdir(exist_ok=False)
    selection_path=root/'results/north_recent_activity_context_v1/selection_before_waveforms.json'
    selection=json.loads(selection_path.read_text());lat=selection['activity_center_latitude'];lon=selection['activity_center_longitude']
    # A deliberately bounded, fixed search centered on independent tremor metadata.
    nodes=np.array([(la,lo,depth) for la in lat+np.arange(-6,7)*.05 for lo in lon+np.arange(-8,9)*.05 for depth in [20.,35.,50.]])
    pd.DataFrame(nodes,columns=['latitude','longitude','depth_km']).to_csv(out/'source_grid.csv',index=False)
    write_json(out/'design_before_scan.json',dict(status='running',input_model='original jms, unchanged',panels=['quiet','active'],
        score='Per-station S log-score median/MAD normalization over the common core; clip standardized scores at +/-6; equal mean across eight physical sites; require >=4 sites at >=2 MAD.',
        detection='Maximum over all source nodes at each sample; normalize that maximum by its core median/MAD; thresholds 4/6/8 MAD, 3 s refractory; identical operation on observed and controls.',
        source_grid_nodes=len(nodes),grid_step_degrees=.05,depths_km=[20,35,50],Vs_km_s=3.5,core_guard_s=240,
        shifts='Eight independent physical-site shifts on the 30 s inference grid, signed 60/90/120/150/180 s; nonwrapping.',
        model_probabilities_calibrated=False,production_validated=False,selection_sha256=sha256(selection_path),source_sha256=sha256(__file__)))
    torch.set_num_threads(2);device='cuda' if torch.cuda.is_available() else 'cpu';summaries=[];inputs={}
    for panel in ['quiet','active']:
        source=root/f'results/north_recent_{panel}_picker_v1';stations=pd.read_csv(source/'stations.csv');stations,groups=representative_sites(stations,1.)
        if len(stations)!=8:raise ValueError('Eight independent sites required')
        design=json.loads((source/'design_before_acquisition.json').read_text());origin=pd.Timestamp(design['start_utc']).timestamp()
        scores=[];normalization=[];fs=20.;start=origin-60
        indices=np.arange(int((origin+240-start)*fs),int((origin+3600-240-start)*fs))
        for station in stations.id:
            path=source/f'{station}_processed.npz';inputs[str(path)]=sha256(path)
            with np.load(path) as a:
                if float(a['fs'])!=fs or float(a['start_s'])!=start:raise ValueError('Unequal station clocks')
                raw=a['probabilities'][1].astype(float)
            if not np.isfinite(raw[1200:73200]).all():raise ValueError('Incomplete core predictions')
            log=np.log(np.maximum(raw,1e-8));scale=robust_scale(log[indices])
            if scale is None:raise ValueError('Degenerate station score distribution')
            scores.append((log-scale[0])/scale[1]);normalization.append(dict(station=station,log_score_median=scale[0],log_score_mad=scale[1]))
        scores=np.asarray(scores,np.float32)
        offsets=np.empty((len(nodes),len(stations)),int)
        for j,station in enumerate(stations.itertuples()):
            for i,(la,lo,depth) in enumerate(nodes):
                horizontal=gps2dist_azimuth(la,lo,station.latitude,station.longitude)[0]/1000
                travel=np.hypot(horizontal,depth+station.elevation/1000)/3.5;offsets[i,j]=int(round(travel*fs))
        rng=np.random.default_rng(20260914)
        shifts=[np.zeros(len(stations),int)]+[rng.choice([60,90,120,150,180],len(stations))*rng.choice([-1,1],len(stations))*int(fs) for _ in range(8)]
        if min(indices)+offsets.min()-max(np.max(x) for x in shifts)<1200 or max(indices)+offsets.max()-min(np.min(x) for x in shifts)>=73200:
            raise ValueError('Guard insufficient for travel times and controls')
        write_json(out/f'{panel}_normalization.json',dict(stations=normalization,shifts_s=[(s/fs).tolist() for s in shifts],station_order=stations.id.tolist(),physical_groups=groups))
        for realization,shift in enumerate(shifts):
            all_scores,gated,node=score_grid(scores,offsets,indices,shift,device)
            scale=robust_scale(all_scores)
            if scale is None:raise ValueError('Degenerate network scores')
            z=(gated-scale[0])/scale[1]
            for threshold in [4.,6.,8.]:
                peaks=find_peaks(z,height=threshold,distance=int(3*fs))[0];rows=[]
                for k in peaks:
                    source_id=int(node[k]);la,lo,depth=nodes[source_id]
                    rows.append(dict(time_s=float(start+indices[k]/fs),score_mad=float(z[k]),network_score=float(gated[k]),latitude=float(la),longitude=float(lo),depth_km=float(depth),grid_node=source_id,quality_status='unvalidated_network_score_proposal'))
                pd.DataFrame(rows,columns=['time_s','score_mad','network_score','latitude','longitude','depth_km','grid_node','quality_status']).to_parquet(out/f'{panel}_r{realization}_t{threshold:g}_proposals.parquet',index=False)
                summaries.append(dict(panel=panel,realization=realization,threshold=threshold,count=len(rows),exposure_s=len(indices)/fs))
            np.savez_compressed(out/f'{panel}_r{realization}_scores.npz',network_all=all_scores,network_gated=gated,node=node,start_s=start+indices[0]/fs,fs=fs,median=scale[0],mad=scale[1])
            print(panel,'realization',realization,'counts',[r['count'] for r in summaries[-3:]],flush=True)
    pd.DataFrame(summaries).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_exploratory_network_score_scan',production_validated=False,inputs=inputs,raw_bytes_downloaded=0,
        limitations=['Scores are standardized model outputs, not calibrated likelihood ratios.','Hypocenters and origin times use an uncalibrated homogeneous model.','Discovery proposals require independent waveform and phase validation.']))

if __name__=='__main__':main()
