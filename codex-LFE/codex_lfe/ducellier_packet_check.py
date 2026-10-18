"""Template registration and neural input comparison on a small roster-matched panel."""
from pathlib import Path
import json,os
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
from .published_pickle import read_published_pickle
from .neural import orientation_matrix,verified_zne,infer_probabilities,picks_from_probabilities
from .templates import cut
from .waveforms import Waveform,preprocess_stream
from .provenance import write_json,sha256


def main():
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch
    import seisbench.models as sbm
    torch.set_num_threads(2)
    root=Path('results/ducellier_roster_packets_v1')
    out=Path('results/ducellier_roster_packet_check_v1');out.mkdir(exist_ok=True)
    manifest=json.loads((root/'manifest.json').read_text())
    inv=read_inventory(str(root/'responses/BK.xml'))+read_inventory(str(root/'responses/NC.xml'))+read_inventory(str(root/'responses/PB.xml'))
    tp=Path('results/southern_benchmark_inventory_v1/source/src/templates_v1/080424.05.060')
    weight=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    model=sbm.LFEDetect.load(str(weight),map_location='cpu',weights_only=True)
    device='cuda' if torch.cuda.is_available() else 'cpu'
    write_json(out/'design_before_check.json',dict(input_sha256=sha256(root/'manifest.json'),source_sha256=sha256(__file__),weights_sha256=sha256(weight.with_suffix('.pt')),template_band_hz=[1.5,9],model_band_hz=[1,8],model_score_threshold=.3,template_time_rule='Catalog time is tested as first sample of the stored truncated template. Report fixed-time per-channel mean CC and ±2s lag diagnostic; no clock adjustment fitted.',model_input_arms=['native_counts_gain_equalized','response_corrected_velocity'],purpose='Small reference-enriched transfer/input diagnostic, not event precision or unbiased recall. Template channels were independently RMS-normalized by source authors; no absolute amplitude calibration implied.'))
    rows=[];allp=[]
    for idx,record in enumerate(manifest['records']):
        if record['status']!='received':continue
        station=record['station'];sta=station.split('.')[1]
        done=out/f'{idx:03d}.json'
        if done.exists():rows.append(json.loads(done.read_text()));continue
        t=pd.Timestamp(record['reference_utc']).timestamp()+record['shift_s']
        channels=record['channels'];stream=read(record['path'])
        corrected=stream.copy()
        for tr in corrected:
            if tr.id not in channels:continue
            if str(inv.get_response(tr.id,tr.stats.starttime))!=str(inv.get_response(tr.id,tr.stats.endtime)):
                raise ValueError('Response epoch crosses packet')
            tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80)
        current_angles=[inv.get_orientation(c,UTCDateTime(t)) for c in channels]
        matrix=orientation_matrix(current_angles)
        template=[];angles=[]
        for channel in channels:
            tr,angle=read_published_pickle(tp/f'{sta}_{channel.split(".")[-1]}.pkl')
            if tr.stats.sampling_rate!=20:raise ValueError('Template rate')
            template.append(tr.data);angles.append(angle)
        template=np.asarray(template)
        w=preprocess_stream(corrected,channels,station,t-30,100,fs=20,band=(1.5,9))
        target=orientation_matrix(angles)@np.linalg.solve(matrix,w.x)
        wave=Waveform(target,w.valid,w.start_s,w.fs,tuple(channels),station)
        values=[]
        for lag in np.arange(-2,2.001,.05):
            packet=cut(wave,t+lag,template.shape[1],alignment='sinc')
            if packet is None:values.append(np.nan);continue
            a=packet-packet.mean(axis=1,keepdims=True);b=template-template.mean(axis=1,keepdims=True)
            den=np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)
            values.append(float(np.mean(np.sum(a*b,axis=1)/den)))
        values=np.asarray(values);best=int(np.nanargmax(values)) if np.isfinite(values).any() else None
        row=dict(station=station,reference_utc=record['reference_utc'],shift_s=record['shift_s'],template_samples=int(template.shape[1]),fixed_time_mean_cc=float(values[40]) if np.isfinite(values[40]) else None,max_mean_cc=float(values[best]) if best is not None else None,best_lag_s=float(-2+best*.05) if best is not None else None,raw_sha256=sha256(record['path']),model_arms={})
        archives={}
        for arm,data in [('native_counts_gain_equalized',stream),('response_corrected_velocity',corrected)]:
            w=preprocess_stream(data,channels,station,t-30,100,fs=20,band=(1,8))
            if arm=='native_counts_gain_equalized':w,orientation=verified_zne(w,inv)
            else:w=Waveform(np.linalg.solve(matrix,w.x).astype(np.float32),w.valid,w.start_s,w.fs,('Z','N','E'),station)
            prob,n=infer_probabilities(w,model,device)
            picks=picks_from_probabilities(prob,w.start_s,20,station,.05)
            picks['input_arm']=arm;picks['packet_index']=idx;allp.append(picks)
            inside=(np.arange(prob.shape[1])/20+w.start_s>=t)&(np.arange(prob.shape[1])/20+w.start_s<t+template.shape[1]/20)
            selected=picks[(picks.model_score>=.3)&(picks.time_s>=t)&(picks.time_s<t+template.shape[1]/20)]
            row['model_arms'][arm]=dict(n_windows=n,P_max_probability=float(np.nanmax(prob[0,inside])) if np.isfinite(prob[0,inside]).any() else None,S_max_probability=float(np.nanmax(prob[1,inside])) if np.isfinite(prob[1,inside]).any() else None,P_picks=int(selected.phase.eq('P').sum()),S_picks=int(selected.phase.eq('S').sum()))
            archives[arm+'_waveform']=w.x;archives[arm+'_probabilities']=prob
        np.savez_compressed(out/f'{idx:03d}_waveforms_probabilities.npz',**archives)
        picks_for_packet=pd.concat(allp[-2:],ignore_index=True);picks_for_packet.to_parquet(out/f'{idx:03d}_picks.parquet',index=False)
        write_json(done,row);rows.append(row)
        print(idx,station,record['shift_s'],'template CC',round(row['fixed_time_mean_cc'],3) if row['fixed_time_mean_cc'] is not None else None,flush=True)
    flat=[]
    for r in rows:
        for arm,a in r['model_arms'].items():flat.append({k:v for k,v in r.items() if k!='model_arms'}|dict(input_arm=arm,**a))
    pd.DataFrame(flat).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_small_packet_diagnostic',station_packets=len(rows),production_picker_validated=False))

if __name__=='__main__':main()
