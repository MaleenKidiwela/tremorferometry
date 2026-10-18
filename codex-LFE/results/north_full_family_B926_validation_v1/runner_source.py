"""Unused-site check with date-separated training QC and reusable checkpoints."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from .pilot import ROOT,date_start
from .reference import read_bostock
from .templates import cut,build_station_template
from .waveforms import load_day
from .matched_filter import normalized_cc,align_to_anchor
from .heldout import sample_series
from .provenance import write_json,sha256

OUT=Path('results/north_full_family_B926_validation_v1')
RUN=Path('results/north_full_family_frozen_evaluation_v1')
STATION='PB.B926';CHANNELS=['PB.B926..EH1','PB.B926..EH2','PB.B926..EHZ']

def load(cfg,date):
    d=pd.Timestamp(date);p=ROOT/cfg['raw_root']/STATION/str(d.year)/f'{d.dayofyear:03d}.mseed'
    return load_day(p,ROOT/'cache/waveforms',CHANNELS,STATION,date_start(date),fs=cfg['fs'],alignment='sinc')

def train(cfg):
    OUT.mkdir(exist_ok=False)
    training=json.loads((RUN/'training.json').read_text());anchors=training['anchors_from_published_s']
    ref=read_bostock(ROOT/cfg['reference_csv']);windows={f:[[],[]] for f in anchors};inputs={}
    design=dict(station=STATION,channels=CHANNELS,training_dates=cfg['train_dates'],test_dates=cfg['test_dates'],minimum_training_events=20,minimum_training_split_cc=.25,split_rule='Alternating whole UTC training dates; evenly sample each half to equal size, then preserve those halves in the existing split-stack checker',evaluation='Signed correlation at exact training-predicted lag; 100 date-level nonwrapping shifts between 300 and 1800 seconds; no lag optimization',source_sha256=sha256(__file__),detector_manifest_sha256=sha256(RUN/'manifest.json'))
    write_json(OUT/'design_before_training.json',design);(OUT/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    for k,date in enumerate(cfg['train_dates']):
        wave,manifest=load(cfg,date);inputs[date]=manifest
        for row in ref[(ref.date==date)&ref.family_id.isin(anchors)].itertuples():
            seg=cut(wave,row.time_s-10,int(70*cfg['fs']),alignment='sinc')
            if seg is not None:windows[row.family_id][k%2].append(seg)
        print('B926 training',date,flush=True)
    archive={};rows=[]
    for family,(a,b) in windows.items():
        n=min(len(a),len(b))
        balanced=[]
        if n:
            ai=np.linspace(0,len(a)-1,n,dtype=int);bi=np.linspace(0,len(b)-1,n,dtype=int)
            for i,j in zip(ai,bi):balanced.extend([a[i],b[j]])
        template,qc=build_station_template(balanced,cfg['fs'],template_s=cfg['template_s'],min_events=20,min_split_cc=.25)
        rows.append(dict(qc,family_id=family,half_a_windows=len(a),half_b_windows=len(b)))
        if template is not None:
            archive[family+'__template']=template;archive[family+'__offset_s']=np.asarray(qc['template_start_from_published_s']-anchors[family])
    np.savez_compressed(OUT/'templates.npz',**archive);pd.DataFrame(rows).to_csv(OUT/'training_qc.csv',index=False)
    write_json(OUT/'training_inputs.json',inputs)
    write_json(OUT/'training_complete.json',dict(status='training_complete_no_test_read',usable_groups=len(archive)//2,files_sha256={p.name:sha256(p) for p in [OUT/'templates.npz',OUT/'training_qc.csv',OUT/'training_inputs.json',OUT/'design_before_training.json']}))
    print('B926 coherent date-separated groups:',len(archive)//2,flush=True)

def evaluate(cfg):
    record=json.loads((OUT/'training_complete.json').read_text())
    for name,h in record['files_sha256'].items():
        if sha256(OUT/name)!=h:raise ValueError('Training artifact changed')
    if (OUT/'report.json').exists():raise ValueError('Completed validation immutable')
    with np.load(OUT/'templates.npz') as bank:
        templates={k.removesuffix('__template'):(bank[k],float(bank[k.removesuffix('__template')+'__offset_s'])) for k in bank.files if k.endswith('__template')}
    candidates=pd.read_parquet(RUN/'test_experimental_candidates.parquet');reference=pd.read_parquet(RUN/'test_eligible_reference.parquet')
    rows=[];inputs={};shifts_record={}
    for date in cfg['test_dates']:
        wave,manifest=load(cfg,date);inputs[date]=manifest
        rng=np.random.default_rng(np.random.SeedSequence([20260914,int(date.replace('-',''))]))
        shifts=rng.uniform(300,1800,100)*rng.choice([-1,1],100);shifts_record[date]=shifts.tolist()
        for family,(template,offset) in templates.items():
            cc,valid=normalized_cc(wave.x,template,wave.valid);series,_=align_to_anchor(cc,valid,offset,wave.fs,wave.x.shape[1])
            for kind,frame in [('candidate',candidates),('published_reference',reference)]:
                times=frame[(frame.date==date)&frame.family_id.eq(family)].time_s.to_numpy()
                observed=sample_series(series,times,wave.start_s,wave.fs)
                null=np.array([sample_series(series,times+shift,wave.start_s,wave.fs) for shift in shifts])
                for i,t in enumerate(times):
                    if not np.isfinite(observed[i]):continue
                    finite=np.isfinite(null[:,i]);values=null[finite,i]
                    rows.append(dict(date=date,family_id=family,kind=kind,time_s=float(t),signed_cc=float(observed[i]),valid_controls=int(finite.sum()),mean_control_cc=float(values.mean()) if len(values) else None,paired_cc_excess=float(observed[i]-values.mean()) if len(values) else None))
        print('B926 evaluation',date,flush=True)
    data=pd.DataFrame(rows);data.to_csv(OUT/'observations.csv',index=False)
    summary=[]
    if len(data):
        for (kind,family),table in data.groupby(['kind','family_id']):
            summary.append(dict(kind=kind,family_id=family,n=len(table),median_signed_cc=float(table.signed_cc.median()),mean_paired_cc_excess=float(table.paired_cc_excess.mean())))
    write_json(OUT/'report.json',dict(status='unused_site_diagnostic_not_event_precision',groups_with_usable_training_templates=len(templates),summary=summary,shifts_s=shifts_record,limitations=['No event-level reviewed truth','A weak unused station may be inconclusive','Positive aggregate waveform support does not certify every candidate']))
    write_json(OUT/'evaluation_inputs.json',inputs);print(pd.DataFrame(summary).to_string(index=False),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['train','evaluate'],required=True);args=parser.parse_args()
    cfg=json.loads((RUN/'manifest.json').read_text())['config']
    if STATION in cfg['stations'] or STATION in {s['site_id'] for s in cfg['stations'].values()}:raise ValueError('Validation site used in detection')
    (train if args.stage=='train' else evaluate)(cfg)

if __name__=='__main__':main()
