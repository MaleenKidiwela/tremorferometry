"""Consume completed validation days and check the unchanged unused-site bank."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from .pilot import ROOT
from .north_b926_validation import load,STATION
from .matched_filter import normalized_cc,align_to_anchor
from .heldout import sample_series
from .reference import read_bostock
from .provenance import sha256,write_json


def check_day(source,date,cfg,templates,anchors,reference,out):
    completion=json.loads((source/f'{date}_complete.json').read_text())
    for name,digest in completion['products'].items():
        if sha256(source/f'{date}_{name}.parquet')!=digest:raise ValueError('Completed source product changed')
    candidates=pd.read_parquet(source/f'{date}_candidates.parquet');candidates=candidates[candidates.realization==0]
    spans=pd.read_parquet(source/f'{date}_intervals.parquet')
    refs=reference[(reference.date==date)&reference.family_id.isin(templates)].copy()
    refs['time_s']+=refs.family_id.map(anchors)
    supported=candidates[candidates.family_id.isin(templates)]
    rng=np.random.default_rng(np.random.SeedSequence([20260915,int(date.replace('-',''))]))
    shifts=rng.uniform(300,1800,100)*rng.choice([-1,1],100)
    rows=[];null_rows=[];wave=None;manifest=None
    for family,(template,offset) in templates.items():
        rf=refs[refs.family_id==family].copy();eligible=np.zeros(len(rf),bool)
        for span in spans[spans.family_id==family].itertuples():eligible|=(rf.time_s>=span.start_s)&(rf.time_s<span.end_s)
        rf=rf[eligible];det=supported[supported.family_id==family]
        if det.empty and rf.empty:continue
        if wave is None:wave,manifest=load(cfg,date)
        cc,valid=normalized_cc(wave.x,template,wave.valid)
        series,_=align_to_anchor(cc,valid,offset,wave.fs,wave.x.shape[1])
        for kind,frame,key in [('candidate',det,'candidate_id'),('reference',rf,'reference_id')]:
            times=frame.time_s.to_numpy();ids=frame[key].to_numpy()
            observed=sample_series(series,times,wave.start_s,wave.fs)
            null=np.array([sample_series(series,times+s,wave.start_s,wave.fs) for s in shifts])
            for i,t in enumerate(times):
                values=null[:,i];finite=np.isfinite(values);obs=float(observed[i])
                mean=float(values[finite].mean()) if finite.any() else None
                rows.append(dict(date=date,family_id=family,kind=kind,record_id=ids[i],time_s=float(t),
                    observed_covered=bool(np.isfinite(obs)),signed_cc=obs,valid_controls=int(finite.sum()),
                    mean_control_cc=mean,paired_cc_excess=obs-mean if mean is not None else None))
                for j,value in enumerate(values):null_rows.append(dict(date=date,family_id=family,kind=kind,record_id=ids[i],realization=j+1,shift_s=float(shifts[j]),signed_cc=float(value)))
    columns=['date','family_id','kind','record_id','time_s','observed_covered','signed_cc','valid_controls','mean_control_cc','paired_cc_excess']
    frame=pd.DataFrame(rows,columns=columns);frame.to_parquet(out/f'{date}_observations.parquet',index=False)
    pd.DataFrame(null_rows,columns=['date','family_id','kind','record_id','realization','shift_s','signed_cc']).to_parquet(out/f'{date}_controls.parquet',index=False)
    if manifest is not None:write_json(out/f'{date}_waveform_input.json',manifest)
    report=dict(date=date,status='complete_unused_site_day',candidate_rows=len(candidates),
        candidates_with_training_supported_family=len(supported),candidates_without_supported_family=len(candidates)-len(supported),
        covered_candidate_checks=int(((frame.kind=='candidate')&frame.observed_covered).sum()),
        source_completion_sha256=sha256(source/f'{date}_complete.json'),shifts_s=shifts.tolist())
    write_json(out/f'{date}_complete.json',report)
    print(date,'unused-site check complete:',len(supported),'supported-family candidates of',len(candidates),flush=True)
    return report


def main():
    out=ROOT/'results/north_extended_B926_validation_v1';out.mkdir(exist_ok=False)
    training=ROOT/'results/north_full_family_B926_validation_v1'
    cfg=json.loads((ROOT/'results/north_full_family_bank_2012_v1/manifest.json').read_text())['config']
    if STATION in cfg['stations'] or STATION in {s['site_id'] for s in cfg['stations'].values()}:raise ValueError('Held-out station used by detector')
    record=json.loads((training/'training_complete.json').read_text())
    for name,digest in record['files_sha256'].items():
        if sha256(training/name)!=digest:raise ValueError('Frozen held-out template artifact changed')
    with np.load(training/'templates.npz') as bank:
        templates={k.removesuffix('__template'):(bank[k].copy(),float(bank[k.removesuffix('__template')+'__offset_s'])) for k in bank.files if k.endswith('__template')}
    anchors=json.loads((ROOT/'results/north_full_family_bank_2012_v1/training.json').read_text())['anchors_from_published_s']
    sources=[ROOT/'results/north_extended_family_validation_v1',ROOT/'results/north_interleaved_family_validation_v1']
    work={d:source for source in sources for d in json.loads((source/'design_before_scan.json').read_text())['dates']}
    reference=read_bostock(cfg['reference_csv'])
    write_json(out/'design_before_checks.json',dict(status='awaiting_completed_validation_days',production_validated=False,
        dates=sorted(work),station=STATION,training_groups=len(templates),training_repeated=False,
        procedure='Frozen date-separated B926 templates and exact predicted lag, without maximizing over lag or polarity. 100 common date-level nonwrapping shifts, signed 300-1800 s.',
        reporting='Retain covered and uncovered candidate checks, unsupported-family counts, all control scores, and reference checks separately.',
        limitations=['Station-level support is not individual-event truth.','Missing B926 support is not proof of a false detection.','Family hypotheses may refer to the same physical event; summary rows are not unique-event totals.'],
        raw_bytes_to_download=0,sources={str(p):sha256(p) for p in [Path(__file__),training/'training_complete.json',training/'templates.npz',ROOT/'codex_lfe/north_b926_validation.py',ROOT/'codex_lfe/heldout.py']}))
    done={};deadline=time.monotonic()+7200
    while len(done)<len(work):
        progress=False
        for date,source in sorted(work.items(),reverse=True):
            if date in done or not (source/f'{date}_complete.json').exists():continue
            done[date]=check_day(source,date,cfg,templates,anchors,reference,out);progress=True
        write_json(out/'progress.json',dict(completed_dates=sorted(done),total_dates=len(work),state='checking_or_waiting_for_detector_days'))
        if len(done)==len(work):break
        if time.monotonic()>deadline:raise TimeoutError('Detector days not complete within two hours; existing results preserved')
        if not progress:time.sleep(20)
    observations=pd.concat([pd.read_parquet(out/f'{d}_observations.parquet') for d in sorted(done)],ignore_index=True)
    observations.to_parquet(out/'observations.parquet',index=False)
    summary=[]
    for (kind,family),frame in observations[observations.observed_covered].groupby(['kind','family_id']):
        summary.append(dict(kind=kind,family_id=family,n=len(frame),median_signed_cc=float(frame.signed_cc.median()),mean_paired_cc_excess=float(frame.paired_cc_excess.mean())))
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_extended_unused_site_validation',production_validated=False,
        individual_events_certified=0,raw_bytes_downloaded=0,days=list(done.values()),summary=summary))


if __name__=='__main__':main()
