"""Complete previously unmeasurable phase checks using count-balanced whole dates.

Reuse all full-stack correlation results and every already-measurable date check.
Only date membership counts choose new halves; evaluation waveforms are not read.
"""
from pathlib import Path
import json,re,itertools
import numpy as np
import pandas as pd
import yaml,h5py
from scipy.signal import butter,sosfiltfilt
from obspy import read_inventory
from .pilot import load_waves
from .templates import cut,stack_windows
from .published_stack_clock import response_correct,best_s_match,phase_site_clocks
from .reference import read_bostock
from .provenance import sha256,write_json


def balanced_dates(counts):
    dates=sorted(d for d,n in counts.items() if n>0)
    if len(dates)<2:return None
    choices=[]
    for bits in itertools.product([False,True],repeat=len(dates)-1):
        a=[dates[0]]+[d for d,b in zip(dates[1:],bits) if b]
        b=[d for d in dates if d not in a]
        if not b:continue
        na,nb=sum(counts[d] for d in a),sum(counts[d] for d in b)
        choices.append((abs(na-nb),tuple(a),tuple(b),na,nb))
    _,a,b,na,nb=min(choices)
    return dict(a=list(a),b=list(b),n_a=int(na),n_b=int(nb)) if min(na,nb)>=10 else None


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'results/north_date_stability_completion_v1';out.mkdir(exist_ok=False)
    config=root/'configs/north_published_phase_transfer_2012.yaml';cfg=yaml.safe_load(config.read_text())
    previous=root/'results/north_published_phase_transfers_v1';stacks=root/'results/north_published_phase_training_stacks_v1'
    published=root.parent/'lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia'
    transfer=pd.read_csv(previous/'station_transfers.csv');manifest=json.loads((stacks/'manifest.json').read_text())
    old=np.load(stacks/'stacks.npz',allow_pickle=False)
    reference=read_bostock(root/cfg['reference_csv']).set_index('reference_id')
    used=pd.DataFrame(manifest['used_reference_station_records']).join(reference[['family_id','date','time_s']],on='reference_id',validate='many_to_one')
    used=used[[f'{r.family_id}__{r.station}__{r.half}' in old for r in used.itertuples()]].copy()
    plans={};eligible=transfer[~transfer.date_stability_measurable & (transfer.correlation>=.6) & ((transfer.correlation-transfer.alternative_peak_cc)>=.1)]
    for (family,station),_ in eligible.groupby(['family_id','station']):
        rows=used[(used.family_id==family)&(used.station==station)]
        partition=balanced_dates(rows.groupby('date').size().to_dict())
        if partition:plans[(family,station)]=partition
    write_json(out/'design_before_waveform_checks.json',dict(status='running',criterion='minimize whole-date count imbalance; chronological lexicographic tie break; >=10 windows per half; same .6 correlation/.1 peak margin/.5 s drift criteria',
        plans=[dict(family_id=f,station=s,**v) for (f,s),v in plans.items()],
        reuse='Full-stack correlations and already-measurable date checks remain unchanged; original used cohort preserved, including omitted singleton-half exclusions.',
        previous_transfers_sha256=sha256(previous/'station_transfers.csv')))
    windows={(f,s,h):[] for f,s in plans for h in ['a','b']};inputs={}
    for date in cfg['train_dates']:
        wanted=used[(used.date==date)&pd.Series([(f,s) in plans for f,s in zip(used.family_id,used.station)],index=used.index)]
        if wanted.empty:continue
        local=cfg.copy();local['stations']={s:v for s,v in cfg['stations'].items() if s in set(wanted.station)}
        waves=load_waves(local,date,inputs)
        for row in wanted.itertuples():
            w=cut(waves[row.station],row.time_s-10,int(70*cfg['fs']),alignment='sinc')
            if w is None:raise ValueError('Previously used training window is no longer valid')
            half='a' if date in plans[row.family_id,row.station]['a'] else 'b'
            windows[row.family_id,row.station,half].append(w)
        print('Completed cached training date',date,flush=True)
    inv=read_inventory(root/'cache/metadata/cn_pilot_2012.xml');meta=pd.read_csv(published/'metadata.csv')
    transfer['date_check_source']='original_chronological_halves'
    records=[];arrays={}
    with h5py.File(published/'waveforms.hdf5') as h5:
        for (family,station),plan in plans.items():
            native={h:stack_windows(windows[family,station,h]) for h in ['a','b']}
            merged=sum(native[h]*len(windows[family,station,h]) for h in native)/(plan['n_a']+plan['n_b'])
            original_halves={h:old[f'{family}__{station}__{h}'] for h in ['early','late'] if f'{family}__{station}__{h}' in old}
            counts=used[(used.family_id==family)&(used.station==station)].groupby('half').size()
            original=sum(x*counts[h] for h,x in original_halves.items())/counts.sum()
            error=float(np.max(np.abs(merged-original)))
            if error>1e-10:raise ValueError('Full training stack changed')
            aligned={h:response_correct(x,cfg['stations'][station]['channels'],inv,cfg['train_dates'],cfg['fs']) for h,x in native.items()}
            for h,x in native.items():arrays[f'{family}__{station}__{h}']=x
            for index,row in eligible[(eligible.family_id==family)&(eligible.station==station)].iterrows():
                found=meta[(meta.source_family_id==row.published_source_id)&(meta.station_code==station.split('.')[1])]
                if len(found)!=1:raise ValueError('Published source/station identity ambiguous')
                m=found.iloc[0];key=re.fullmatch(r'(\w+)\$(\d+),:3,:(\d+)',m.trace_name)
                pub=h5['data'][key[1]][int(key[2]),:3,:int(key[3])]
                pub=sosfiltfilt(butter(4,[1,8],fs=cfg['fs'],btype='bandpass',output='sos'),pub,axis=1)
                before,duration,channels=(1.,4.,1) if row.phase=='P' else (2.,8.,3)
                start=int(round(m[f'trace_{row.phase}_arrival_sample']-before*cfg['fs']))
                template=pub[:channels,start:start+int(duration*cfg['fs'])]
                matches={h:best_s_match(x[:channels],template,cfg['fs'],before) for h,x in aligned.items()}
                drift=abs(matches['a'][0]-matches['b'][0])
                transfer.loc[index,['date_stability_measurable','date_half_drift_s','accepted_waveform_transfer','date_check_source']]=[True,drift,drift<=.5,'count_balanced_whole_dates']
                records.append(dict(family_id=family,station=station,phase=row.phase,drift_s=drift,n_a=plan['n_a'],n_b=plan['n_b'],full_stack_reproduction_max_error=error))
    transfer.to_csv(out/'station_transfers.csv',index=False);pd.DataFrame(records).to_csv(out/'completed_date_checks.csv',index=False)
    np.savez_compressed(out/'balanced_native_stacks.npz',**arrays)
    clocks=[]
    for family,g in transfer.groupby('family_id'):
        v=phase_site_clocks(g);spread=float(v.max()-v.min()) if len(v) else None;accepted=len(v)>=2 and spread<=1.
        clocks.append(dict(family_id=family,n_independent_sites=len(v),clock_spread_s=spread,accepted=bool(accepted),origin_minus_published_s=float(v.median()) if accepted else None))
    write_json(out/'family_clocks.json',clocks);write_json(out/'waveform_inputs.json',inputs)
    write_json(out/'report.json',dict(status='complete_training_date_check_completion',new_checks=len(records),
        accepted_transfers=int(transfer.accepted_waveform_transfer.sum()),accepted_model_dependent_family_clocks=sum(x['accepted'] for x in clocks),
        production_validated=False,absolute_origin_time_validated=False,full_stack_scores_recomputed=False,
        inputs={str(p):sha256(p) for p in [Path(__file__),previous/'station_transfers.csv',stacks/'manifest.json',stacks/'stacks.npz',config,published/'metadata.csv',published/'waveforms.hdf5',root/'cache/metadata/cn_pilot_2012.xml']}))
    print((out/'report.json').read_text(),flush=True)

if __name__=='__main__':main()
