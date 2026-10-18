"""Compare adapted and original pickers on source-held-out natural phase proxies.

Phase centers come from training-stack transfers, not reviewed single-event picks.
Only cached calibration dates are used; no threshold or weights are fitted here.
"""
from pathlib import Path
import json,os
import numpy as np
import pandas as pd
import yaml
from obspy import read_inventory
from .pilot import load_waves,date_start
from .published_stack_clock import response_correct
from .reference import read_bostock
from .neural import picks_from_probabilities
from .miss_diagnosis import phase_window
from .evaluate import match_events
from .provenance import sha256,write_json


def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/north_adapted_natural_reference_v1';out.mkdir(exist_ok=False)
    adaptation=root/'results/north_picker_adaptation_v1';transfers=root/'results/north_date_stability_completion_v1/station_transfers.csv'
    cfg=yaml.safe_load((root/'configs/north_published_phase_transfer_2012.yaml').read_text())
    splits=pd.read_csv(adaptation/'source_splits.csv');transfer=pd.read_csv(transfers)
    transfer=transfer[transfer.accepted_waveform_transfer].merge(splits,left_on='published_source_id',right_on='source_family_id')
    transfer=transfer[transfer.split=='test'].copy()
    if transfer.empty:raise ValueError('No transfer-supported held-out source groups')
    dates=cfg['calibration_dates'];reference_path=root/cfg['reference_csv'];reference=read_bostock(reference_path)
    reference=reference[reference.date.isin(dates)&reference.family_id.isin(transfer.family_id)]
    references=reference.merge(transfer[['family_id','station','phase','phase_arrival_after_published_s']],on='family_id')
    references['expected_phase_s']=references.time_s+references.phase_arrival_after_published_s
    references.to_parquet(out/'reference_phase_proxies.parquet',index=False)
    shifts=[0.,-300.,-120.,120.,300.];radii=[.5,1.]
    write_json(out/'design_before_inference.json',dict(status='running',dates=dates,n_reference_events=len(reference),n_station_phase_proxies=len(references),
        phase_timing_basis='Published event time plus fixed training-only station/phase transfer; unreviewed individual-arrival proxies',
        source_scope='Fine-tuning test groups only; source groups and all synthetic test examples were excluded from weight fitting and epoch selection',
        natural_split_scope='Existing 2012 calibration dates are disjoint from transfer training and from 2026 adaptation noise; not claimed untouched by every earlier project.',
        thresholds=[.3],radii_s=radii,control_shifts_s=shifts,
        preprocessing='Cached masked native 1-8 Hz/20 Hz; response correction and ZNE rotation on 90 s contexts, retain central 60 s; central 30 s prediction output on the original UTC half-minute grid.',
        selection='All phase proxies from the declared dates and supported test-source groups, not selected by model output.',
        production_validated=False,reference_sha256=sha256(reference_path),transfer_sha256=sha256(transfers),adapted_weights_sha256=sha256(adaptation/'adapted_state_dict.pt')))
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(root/'cache/seisbench'))
    import torch,seisbench.models as sbm
    torch.set_num_threads(2);device='cuda' if torch.cuda.is_available() else 'cpu'
    original=root.parent/'lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms'
    models={name:sbm.LFEDetect.load(str(original),map_location='cpu',weights_only=True).to(device).eval() for name in ['baseline','adapted']}
    models['adapted'].load_state_dict(torch.load(adaptation/'adapted_state_dict.pt',map_location='cpu',weights_only=True))
    inv_path=root/'cache/metadata/cn_pilot_2012.xml';inv=read_inventory(inv_path);all_inputs={};rows=[];summaries=[];matched_rows=[]
    for date in dates:
        day=references[references.date==date];config=cfg.copy();config['stations']={s:v for s,v in cfg['stations'].items() if s in set(day.station)}
        waves=load_waves(config,date,all_inputs)
        for station,wave in waves.items():
            target=day[day.station==station].reset_index(drop=True);n=wave.x.shape[1];fs=wave.fs;starts=set()
            for t in target.expected_phase_s:
                for shift in shifts:
                    for delta in [-2.,2.]:
                        position=(t+shift+delta-wave.start_s)*fs
                        start=int(np.floor((position-300)/600))*600
                        if start>=300 and start+1500<=n:starts.add(start)
            contexts=[];kept=[]
            for start in sorted(starts):
                lo,hi=start-300,start+1500
                if not wave.valid[:,lo:hi].all():continue
                zne=response_correct(wave.x[:,lo:hi],config['stations'][station]['channels'],inv,[date],fs)
                contexts.append(zne[:,300:1500].astype(np.float32));kept.append(start)
            if not contexts:continue
            x=np.asarray(contexts)
            for name,model in models.items():
                prob=np.full((2,n),np.nan,np.float32)
                with torch.inference_mode():
                    for first in range(0,len(x),128):
                        tensor=torch.from_numpy(x[first:first+128]).to(device);tensor=model.annotate_batch_pre(tensor,{})
                        pred=model(tensor).cpu().numpy()
                        for start,value in zip(kept[first:first+128],pred):prob[:,start+300:start+900]=value[:,300:900]
                picks=picks_from_probabilities(prob,wave.start_s,fs,station,.3)
                np.savez_compressed(out/f'{date}_{station}_{name}_probabilities.npz',probabilities=prob,start_s=wave.start_s,fs=fs)
                picks.to_parquet(out/f'{date}_{station}_{name}_picks.parquet',index=False)
                for phase_index,phase in enumerate(['P','S']):
                    tt=target[target.phase==phase].reset_index(drop=True);pp=picks[picks.phase==phase].sort_values('time_s').reset_index(drop=True)
                    for shift in shifts:
                        for radius in radii:
                            covered=[]
                            for ref in tt.itertuples():
                                hit,maximum,peak,residual=phase_window(prob[phase_index],wave.start_s,fs,ref.expected_phase_s+shift,radius,pp.time_s.to_numpy(),pp.model_score.to_numpy())
                                covered.append(hit);rows.append(dict(model=name,date=date,station=station,phase=phase,reference_id=ref.reference_id,family_id=ref.family_id,shift_s=shift,radius_s=radius,covered=bool(hit),max_probability=maximum,max_peak_score=peak,peak_residual_s=residual))
                            refs=tt[np.asarray(covered,dtype=bool)].copy();refs['time_s']=refs.expected_phase_s+shift;refs['time_basis']='station_phase_arrival_utc'
                            detections=pp.copy();detections['family_id']=None
                            matches=match_events(detections,refs,radius,same_family=False)
                            summaries.append(dict(model=name,date=date,station=station,phase=phase,shift_s=shift,radius_s=radius,n_reference=len(tt),n_covered=len(refs),n_matched=len(matches)))
                            for match in matches.itertuples():
                                matched_rows.append(dict(model=name,date=date,station=station,phase=phase,shift_s=shift,radius_s=radius,reference_id=refs.iloc[int(match.reference_index)].reference_id,pick_time_s=float(detections.iloc[int(match.detection_index)].time_s),residual_s=float(match.residual_s)))
            print(date,station,'contexts',len(kept),'proxies',len(target),flush=True)
    pd.DataFrame(rows).to_parquet(out/'proxy_windows.parquet',index=False);summary=pd.DataFrame(summaries);summary.to_csv(out/'summary_by_station_day.csv',index=False)
    totals=summary.groupby(['model','phase','shift_s','radius_s'])[['n_reference','n_covered','n_matched']].sum().reset_index();totals['proxy_match_fraction']=totals.n_matched/totals.n_covered;totals.to_csv(out/'summary.csv',index=False)
    pd.DataFrame(matched_rows).to_csv(out/'one_to_one_matches.csv',index=False);write_json(out/'waveform_inputs.json',all_inputs)
    write_json(out/'report.json',dict(status='complete_source_held_out_natural_proxy_check',production_validated=False,individual_phase_truth_available=False,
        interpretation='One-to-one agreement with transferred phase proxies; shifted windows measure accidental agreement, not a labeled false-positive rate.',
        sources={str(p):sha256(p) for p in [Path(__file__),root/'codex_lfe/published_stack_clock.py',root/'codex_lfe/neural.py',root/'codex_lfe/miss_diagnosis.py',root/'codex_lfe/evaluate.py',inv_path,transfers,adaptation/'source_splits.csv',adaptation/'adapted_state_dict.pt',original.with_suffix('.pt')]}))
    print(totals.to_string(index=False),flush=True)

if __name__=='__main__':main()
