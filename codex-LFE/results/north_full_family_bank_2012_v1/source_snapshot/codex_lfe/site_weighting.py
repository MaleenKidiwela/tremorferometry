"""Calibration-only test of training-stack contrast weights at physical-site level.

Reuses the E5 amplitude-ratio/capped-weight idea, not its fitted weights or
claimed performance. Existing templates, site grouping, masks and controls stay
fixed. The chosen packet is not independently labeled S; contrast is a proxy.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .pilot import ROOT, load_waves, scan_date_worker, summarize_split
from .provenance import sha256, source_manifest, write_json
from .reference import read_bostock
from .roster import validate_membership
from .templates import cut, stack_windows


def training_weights(run, stacks, cfg, repair_inputs=None):
    training=json.loads((run/'training.json').read_text())
    sm=json.loads((stacks/'manifest.json').read_text())
    if set(sm['train_dates']) != set(cfg['train_dates']):
        raise ValueError('stack dates differ from training dates')
    if sm['fs'] != cfg['fs'] or sm['start_relative_s'] != -10:
        raise ValueError('unexpected stack sample clock')
    counts=pd.DataFrame(json.loads((stacks/'counts.json').read_text())).set_index(['family_id','station','half']).n_windows
    qc=pd.read_csv(run/'template_qc.csv').set_index(['family_id','station'])
    # Review plots omit halves with a single event. Restore exactly those windows
    # from cached raw data instead of quietly changing the template population.
    repair_inputs = {} if repair_inputs is None else repair_inputs
    reference=read_bostock(cfg['reference_csv']).set_index('reference_id')
    used=pd.DataFrame(sm['used_reference_station_records']).join(
        reference[['family_id','date','time_s']],on='reference_id',validate='many_to_one')
    restored_waves={}
    templates={};rows=[];fs=cfg['fs']
    with np.load(stacks/'stacks.npz') as archive, np.load(run/'templates.npz') as original:
        for family in training['families_usable']:
            templates[family]={}
            for station,info in cfg['stations'].items():
                key=f'{family}__{station}'
                if key+'__template' not in original: continue
                halves={h:archive[f'{key}__{h}'] for h in ['early','late'] if f'{key}__{h}' in archive}
                for half in ['early','late']:
                    if half in halves or counts.loc[family,station,half]==0: continue
                    missing=used[(used.family_id==family)&(used.station==station)&(used.half==half)]
                    if len(missing)!=1 or counts.loc[family,station,half]!=1:
                        raise ValueError('unexpected missing training-stack half')
                    event=missing.iloc[0]
                    if event.date not in cfg['train_dates']: raise ValueError('nontraining singleton')
                    if event.date not in restored_waves:
                        restored_waves[event.date]=load_waves(cfg,event.date,repair_inputs)
                    snippet=cut(restored_waves[event.date][station],event.time_s-10,int(70*fs),
                                alignment=cfg.get('template_alignment','linear'))
                    if snippet is None: raise ValueError('missing singleton waveform')
                    halves[half]=stack_windows([snippet])
                total=sum(counts.loc[family,station,h] for h in halves)
                row=qc.loc[family,station]
                if total != row.n_windows: raise ValueError('training stack event count mismatch')
                merged=sum(x*counts.loc[family,station,h] for h,x in halves.items())/total
                start=int(round((row.template_start_from_published_s+10)*fs))
                template=original[key+'__template'];offset=float(original[key+'__offset_s'])
                error=float(np.max(np.abs(merged[:,start:start+template.shape[1]]-template)))
                if error>1e-10: raise ValueError('training stack does not reproduce frozen template')
                templates[family][station]=(template,offset)
                first=int(round((row.peak_offset_from_published_s+10-.5)*fs)); n=int(4*fs)
                noise_first=first-int(15*fs)
                if noise_first<0 or first+n>merged.shape[1]: raise ValueError('contrast window outside training stack')
                signal=float(np.mean(merged[:,first:first+n]**2))
                background=float(np.mean(merged[:,noise_first:noise_first+n]**2))
                if not signal>0 or not background>0: raise ValueError('invalid training contrast')
                ratio=float(np.sqrt(signal/background))
                rows.append({'family_id':family,'station':station,'site':info['site_id'],
                    'n_training_windows':int(total),'stack_amplitude_ratio':ratio,
                    'stack_contrast_db':float(20*np.log10(ratio)),
                    'template_reproduction_max_abs_error':error})
    station=pd.DataFrame(rows)
    # A colocated pair cannot double a site's network weight by instrument count.
    site=station.groupby(['family_id','site']).stack_amplitude_ratio.mean().rename('uncapped_weight').reset_index()
    median=float(site.uncapped_weight.median())
    site['weight']=np.minimum(site.uncapped_weight,3*median)/median
    weights={f:g.set_index('site').weight.to_dict() for f,g in site.groupby('family_id')}
    return templates,weights,station,site


def operating_point(scores,cfg):
    for key in sorted(scores,key=float):
        s=scores[key]
        if (s['reference']['n_detections']>=cfg['calibration_min_detections'] and s['null'] and
            s['null']['control_to_observed_range'][1]<=cfg['calibration_max_control_ratio']):
            return float(key)
    return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',required=True);ap.add_argument('--stacks',required=True)
    ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args(); run,stacks,out=map(Path,[args.run,args.stacks,args.output])
    if not 1<=args.workers<=8: raise ValueError('workers must be 1..8')
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((run/'manifest.json').read_text())['config']
    roster=validate_membership(cfg['stations'],ROOT/cfg['station_roster'])
    repair_inputs={}
    templates,weights,station,site=training_weights(run,stacks,cfg,repair_inputs)
    cfg['site_weights']=weights
    station.to_csv(out/'training_station_contrast.csv',index=False)
    site.to_csv(out/'training_site_weights.csv',index=False)
    files=[run/'manifest.json',run/'training.json',run/'template_qc.csv',run/'templates.npz',
           run/'calibration_scores.json',stacks/'manifest.json',stacks/'counts.json',stacks/'stacks.npz',
           Path(cfg['reference_csv']),Path(__file__),
           ROOT.parent/'lfe_catalog_repair/scripts/e5r_templates.py',
           ROOT.parent/'lfe_catalog_repair/scripts/e4_template_test.py']
    design={'status':'frozen_before_calibration_scan','scope':'calibration_only_no_evaluation_scan',
        'config':cfg,'roster':roster,'sources':{str(p.resolve()):sha256(p) for p in files},
        'source_code':source_manifest(ROOT/'codex_lfe'),
        'restored_singleton_waveform_inputs':repair_inputs,
        'rule':'Station weight proxy = RMS(training stack peak -0.5..+3.5s) / RMS(same-length window 15s earlier). Average station proxies within each physical site. Cap at 3x global median site proxy and median-normalize. Fixed weights for observed and controls; four-site support unchanged.',
        'limitations':['Training peak is not independently labeled S; contrast is not calibrated single-event SNR',
          'Earlier stack window may contain signal; this proxy is tested, not asserted optimal',
          'Colocated station mean and site-level cap adapt E5; not an exact replication',
          'Previously examined calibration dates; no new generalization or true precision claim']}
    write_json(out/'design_before_scan.json',design)
    tick=time.monotonic()
    tasks=[(cfg,d,'calibration',templates,out,'numpy','cpu',2) for d in cfg['calibration_dates']]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        completed=list(pool.map(scan_date_worker,tasks))
    pieces=[r for r,_ in completed]
    candidates,exposure,intervals=[pd.concat([p[k] for p in pieces],ignore_index=True) for k in range(3)]
    summary,eligible=summarize_split(cfg,'calibration',candidates,exposure,intervals,
        read_bostock(cfg['reference_csv']),list(templates),json.loads((run/'training.json').read_text())['anchors_from_published_s'])
    eligible.to_parquet(out/'calibration_eligible_reference.parquet',index=False)
    write_json(out/'calibration_scores.json',summary)
    write_json(out/'input_manifest.json',[m for _,m in completed])
    old=json.loads((run/'calibration_scores.json').read_text())
    old_refs=pd.read_parquet(run/'calibration_eligible_reference.parquet')
    if set(eligible.reference_id)!=set(old_refs.reference_id): raise ValueError('reference eligibility changed')
    baseline_intervals=pd.concat([pd.read_parquet(run/f'calibration_{d}_eligible_intervals.parquet') for d in cfg['calibration_dates']],ignore_index=True)
    keys=['date','family_id','start_s','end_s']
    pd.testing.assert_frame_equal(intervals[keys].sort_values(keys).reset_index(drop=True),baseline_intervals[keys].sort_values(keys).reset_index(drop=True))
    comparison=[]
    for threshold,result in summary.items():
        for label,scores in [('baseline',old[threshold]),('training_site_weighted',result)]:
            comparison.append({'threshold':float(threshold),'variant':label,
                'n_detected':scores['reference']['n_detections'],
                'n_correct_group_matches':scores['reference']['n_matched_same_family'],
                'eligible_reference':scores['reference_eligible'],
                'control_counts':json.dumps(scores['null']['control_counts']),
                'max_control_ratio':scores['null']['control_to_observed_range'][1]})
    pd.DataFrame(comparison).to_csv(out/'comparison.csv',index=False)
    selected=operating_point(summary,cfg);baseline=operating_point(old,cfg)
    decision={'status':'complete_calibration_experiment_no_evaluation_scan','elapsed_s':time.monotonic()-tick,
        'baseline_operating_threshold':baseline,'weighted_operating_threshold':selected,
        'shared_reference_count':len(eligible),'identical_waveform_intervals':True,
        'baseline_matches_at_selected':old[f'{baseline:g}']['reference']['n_matched_same_family'] if baseline else None,
        'weighted_matches_at_selected':summary[f'{selected:g}']['reference']['n_matched_same_family'] if selected else None,
        'interpretation':'Compare calibration recovery under the same engineering control gate. Empirical control ratios are not true false-discovery rates. No evaluation scan or automatic production adoption.'}
    write_json(out/'decision.json',decision)
    print(pd.DataFrame(comparison).to_string(index=False),flush=True)
    print(json.dumps(decision),flush=True)


if __name__=='__main__': main()
