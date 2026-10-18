"""Fixed-threshold extended validation with explicit family-specific clocks.

Retain every family hypothesis. No cross-family event merging or physical-origin
claim is made. All controls share exactly the observed family's valid exposure.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, multiprocessing, time
import numpy as np
import pandas as pd
from .pilot import ROOT, load_waves, date_start
from .training_artifact import load_training_artifact
from .matched_filter import normalized_cc, align_to_anchor, combine_sites, shifted_sites, extract_candidates
from .reference import read_bostock
from .evaluate import match_events
from .waveforms import runs
from .provenance import sha256, write_json


def family_basis(family):
    return 'arrival_anchor:CN.PGC:training_stack_peak:source_group='+family


def score_family_rows(candidates, reference, intervals, anchors):
    """Score each source clock separately; never match a different family's peak."""
    rows=[];matches=[]
    for family in anchors:
        refs=reference[reference.family_id==family].copy()
        refs['time_s']+=anchors[family];refs['time_basis']=family_basis(family)
        spans=intervals[intervals.family_id==family]
        eligible=np.zeros(len(refs),bool)
        for span in spans.itertuples(): eligible|=(refs.time_s>=span.start_s)&(refs.time_s<span.end_s)
        refs=refs.loc[eligible].reset_index(drop=True)
        observed=candidates[(candidates.family_id==family)&(candidates.realization==0)].reset_index(drop=True)
        pairs=match_events(observed,refs,1.,True)
        rows.append(dict(family_id=family,n_candidates=len(observed),n_eligible_references=len(refs),n_reference_matches=len(pairs)))
        for pair in pairs.itertuples():
            matches.append(dict(family_id=family,candidate_id=observed.iloc[pair.detection_index].candidate_id,
                reference_id=refs.iloc[pair.reference_index].reference_id,residual_s=pair.residual_s))
    return pd.DataFrame(rows),pd.DataFrame(matches,columns=['family_id','candidate_id','reference_id','residual_s'])


def scan(task):
    cfg,templates,anchors,date,out_path=task;out=Path(out_path);started=time.monotonic();inputs={}
    waves=load_waves(cfg,date,inputs);fs=cfg['fs'];n=int(86400*fs)
    block_n=int(cfg['normalization_block_s']*fs);guard=int(cfg['block_guard_s']*fs)
    all_sites=sorted({cfg['stations'][s]['site_id'] for members in templates.values() for s in members})
    block_shifts={}
    for block in range(n//block_n):
        rng=np.random.default_rng(np.random.SeedSequence([cfg['seed'],int(date.replace('-','')),block]))
        block_shifts[block]=[{s:0 for s in all_sites}]+[
            {s:int(rng.choice([-1,1])*rng.integers(int(cfg['null_shift_min_s']*fs),int(cfg['null_shift_max_s']*fs)+1)) for s in all_sites}
            for _ in range(cfg['null_replicates'])]
    rows=[];exposure=[];intervals=[]
    for family_index,(family,members) in enumerate(templates.items()):
        aligned={}
        for station,(template,offset) in members.items():
            cc,valid=normalized_cc(waves[station].x,template,waves[station].valid)
            aligned[station],_=align_to_anchor(cc,valid,offset,fs,n)
        full=combine_sites(aligned,{s:cfg['stations'][s]['site_id'] for s in aligned})
        for block,first in enumerate(range(0,n,block_n)):
            last=min(n,first+block_n);local={s:v[first:last] for s,v in full.items()}
            panels=[local]+[shifted_sites(local,shift) for shift in block_shifts[block][1:]]
            common=np.ones(last-first,bool);common[:guard]=False;common[-guard:]=False
            for panel in panels:common&=np.isfinite(np.asarray(list(panel.values()))).all(axis=0)
            block_exposure=[]
            for realization,panel in enumerate(panels):
                found,good=extract_candidates(panel,fs,date_start(date)+first/fs,9.,
                    support_mad=cfg['site_support_mad'],min_sites=cfg['min_sites'],
                    refractory_s=cfg['refractory_s'],common_valid=common)
                block_exposure.append(float(good.sum()/fs))
                for event_index,event in enumerate(found):
                    rows.append(event|dict(candidate_id=f'{date}_{family}_b{block}_r{realization}_e{event_index}',
                        family_id=family,time_basis=family_basis(family),realization=realization,date=date,block=block,
                        quality_status='experimental_family_hypothesis',threshold_mad=9.))
                if realization==0:
                    for lo,hi in runs(good):
                        intervals.append(dict(date=date,family_id=family,start_s=date_start(date)+(first+lo)/fs,end_s=date_start(date)+(first+hi)/fs))
            if len(set(block_exposure))!=1:raise ValueError('Unequal observed/control exposure')
            exposure.append(dict(date=date,family_id=family,block=block,shared_exposure_s=block_exposure[0]))
        if (family_index+1)%5==0:
            write_json(out/f'{date}_progress.json',dict(date=date,families_completed=family_index+1,total_families=len(templates),elapsed_s=time.monotonic()-started))
            print(date,'families',family_index+1,'/',len(templates),flush=True)
    columns=['candidate_id','time_s','score_mad','network_cc','n_sites','n_supporting_sites','supporting_sites','family_id','time_basis','realization','date','block','quality_status','threshold_mad']
    candidates=pd.DataFrame(rows,columns=columns);spans=pd.DataFrame(intervals,columns=['date','family_id','start_s','end_s'])
    reference=read_bostock(cfg['reference_csv']);reference=reference[reference.date==date]
    family_scores,matches=score_family_rows(candidates,reference,spans,anchors)
    frames={'candidates':candidates,'exposure':pd.DataFrame(exposure),'intervals':spans,'family_scores':family_scores,'matches':matches}
    for name,frame in frames.items():frame.to_parquet(out/f'{date}_{name}.parquet',index=False)
    write_json(out/f'{date}_inputs.json',inputs)
    report=dict(date=date,status='complete_family_validation_day',elapsed_s=time.monotonic()-started,
        n_observed_family_rows=int((candidates.realization==0).sum()),n_eligible_references=int(family_scores.n_eligible_references.sum()),
        n_family_matches=int(family_scores.n_reference_matches.sum()),
        control_counts=[int((candidates.realization==r).sum()) for r in range(1,17)],
        shared_family_exposure_s=float(frames['exposure'].shared_exposure_s.sum()),
        products={name:sha256(out/f'{date}_{name}.parquet') for name in frames})
    write_json(out/f'{date}_complete.json',report);print(json.dumps(report),flush=True)
    return report


def main():
    out=ROOT/'results/north_extended_family_validation_v1';out.mkdir(exist_ok=False)
    base=ROOT/'results/north_full_family_bank_2012_v1';training=ROOT/'results/north_full_family_training_recovered_v1'
    original_manifest=json.loads((base/'manifest.json').read_text());cfg=original_manifest['config']
    templates,anchors,_=load_training_artifact(training,cfg,original_manifest['reference_sha256'])
    dates=[str(d.date()) for d in pd.date_range('2012-09-21','2012-09-30')]
    if set(dates)&set(cfg['train_dates']+cfg['calibration_dates']+cfg['test_dates']):raise ValueError('Date split overlap')
    # Every declared input must already exist. No downloader is called here.
    for date in dates:
        dt=pd.Timestamp(date)
        for station in cfg['stations']:
            if not (Path(cfg['raw_root'])/station/str(dt.year)/f'{dt.dayofyear:03d}.mseed').exists():raise ValueError('Missing cached input')
    write_json(out/'design_before_scan.json',dict(status='running_extended_validation',production_validated=False,
        dates=dates,selection='All ten consecutive dates after the previous evaluation, through the end of the same cached month; no activity-based date selection.',
        split_note='New evaluation block for this fixed template bank; not claimed untouched by every previous project.',
        templates='Original 47 frozen groups; no training or threshold retuning',threshold_mad=9.,controls=16,
        exposure='Per-family intersection of valid observed and all shifted-control samples; identical support/exposure for every realization.',
        clock_policy='Family-specific PGC training-peak clocks, possibly different physical phases. Preserve every family hypothesis; no cross-family merging.',
        original_evaluation='Retained unchanged; its failed criterion is not erased or reclassified by this extension.',
        reporting='One-to-one within-family reference agreement; sum of family hypotheses is not a unique event count or calibrated false-discovery rate.',
        descriptive_control_check='At least 20 observed family rows and maximum shifted/observed row ratio <=0.05; this alone cannot establish production validity.',
        raw_bytes_to_download=0,workers=3,config=cfg,
        inputs={str(p):sha256(p) for p in [Path(__file__),training/'templates.npz',training/'training_artifact.json',base/'selection_before_test.json',ROOT/'codex_lfe/matched_filter.py',ROOT/'codex_lfe/evaluate.py',ROOT/'codex_lfe/waveforms.py']}))
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes());completed=[]
    with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(scan,(cfg,templates,anchors,date,str(out))) for date in reversed(dates)]
        for future in as_completed(futures):
            completed.append(future.result());write_json(out/'progress.json',dict(completed_dates=sorted(r['date'] for r in completed),total_dates=len(dates)))
    completed.sort(key=lambda x:x['date']);counts=np.sum([r['control_counts'] for r in completed],axis=0)
    observed=sum(r['n_observed_family_rows'] for r in completed);eligible=sum(r['n_eligible_references'] for r in completed);matched=sum(r['n_family_matches'] for r in completed)
    pd.DataFrame(completed).to_csv(out/'daily_summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_extended_family_validation',production_validated=False,
        n_observed_family_hypotheses=observed,n_eligible_references=eligible,n_family_matches=matched,
        within_family_reference_recall=matched/eligible if eligible else None,control_counts=counts.tolist(),
        maximum_control_ratio=float(counts.max()/observed) if observed else None,
        descriptive_control_check_passed=bool(observed>=20 and counts.max()/observed<=.05),
        unique_event_count=None,raw_bytes_downloaded=0,days=completed))
    print((out/'report.json').read_text(),flush=True)


if __name__=='__main__':main()
