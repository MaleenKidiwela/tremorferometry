"""Freeze a calibration-selected method, scan a held-out block and export candidates."""
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
import numpy as np
import pandas as pd
from .pilot import ROOT, scan_date_worker, summarize_split
from .provenance import sha256, source_manifest, write_json
from .reference import read_bostock
from .roster import validate_membership


def main():
    baseline=Path('results/method_interleaved_mf_v1')
    experiment=Path('results/method_interleaved_short_packet_v1')
    decision=json.loads((experiment/'decision.json').read_text())
    use_short=decision['advance_to_fresh_validation']
    selected=experiment if use_short else baseline
    threshold=decision['short_packet_operating_threshold' if use_short else 'baseline_operating_threshold']
    if threshold is None: raise ValueError('No calibrated operating point')
    cfg=json.loads((baseline/'manifest.json').read_text())['config']
    prior_dates=set(cfg['train_dates']+cfg['calibration_dates']+cfg['test_dates'])
    dates=['2010-09-06','2010-09-07','2010-09-08','2010-09-09','2010-09-10']
    if prior_dates & set(dates): raise ValueError('Overlapping dates')
    cfg['test_dates']=dates
    cfg['thresholds']=[threshold]
    cfg['template_s']=4. if use_short else 8.
    cfg.pop('site_weights',None)
    training=json.loads((baseline/'training.json').read_text())
    templates={f:{} for f in training['families_usable']}
    with np.load(selected/'templates.npz') as bank:
        for f in templates:
            for station in cfg['stations']:
                key=f'{f}__{station}'
                if key+'__template' in bank:
                    templates[f][station]=(bank[key+'__template'].copy(),float(bank[key+'__offset_s']))
    out=Path('results/north_catalog_block_20100906_10_v1')
    out.mkdir(exist_ok=False)
    write_json(out/'selection_before_test.json',dict(config=cfg,selected_method=str(selected),threshold=threshold,
        roster=validate_membership(cfg['stations'],ROOT/cfg['station_roster']),
        date_selection='Contiguous cached five-day block absent from recorded codex-LFE training/calibration/test date lists; not claimed untouched by other projects.',
        decision_source_sha256=sha256(experiment/'decision.json'),template_sha256=sha256(selected/'templates.npz'),
        baseline_manifest_sha256=sha256(baseline/'manifest.json'),source_code=source_manifest(ROOT/'codex_lfe'),
        validation_rule='Report same-family reference recovery, per-family coverage and unchanged maximum control ratio <=0.05; no retuning on this block.',
        release_status='experimental_unverified; arrival-anchor times, not validated origin times'))
    tasks=[(cfg,d,'test',templates,out,'numpy','cpu',2) for d in dates]
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
        completed=list(pool.map(scan_date_worker,tasks))
    pieces=[x for x,_ in completed]
    candidates,exposure,intervals=[pd.concat([x[k] for x in pieces],ignore_index=True) for k in range(3)]
    reference=read_bostock(cfg['reference_csv'])
    scores,eligible=summarize_split(cfg,'test',candidates,exposure,intervals,reference,list(templates),training['anchors_from_published_s'])
    write_json(out/'test_scores.json',scores)
    write_json(out/'input_manifest.json',[m for _,m in completed])
    eligible.to_parquet(out/'eligible_reference.parquet',index=False)
    intervals.to_parquet(out/'coverage_intervals.parquet',index=False)
    exposure.to_parquet(out/'control_exposure.parquet',index=False)
    catalog=candidates[candidates.realization==0].sort_values('time_s').copy()
    catalog.insert(0,'candidate_id',[f'north201009_{i:06d}' for i in range(1,len(catalog)+1)])
    catalog['arrival_anchor_utc']=pd.to_datetime(catalog.time_s,unit='s',utc=True).astype(str)
    catalog['validation_status']='experimental_unverified'
    catalog.to_parquet(out/'candidate_catalog.parquet',index=False)
    catalog.to_csv(out/'candidate_catalog.csv',index=False)
    reference[['family_id','latitude','longitude','depth_km']].drop_duplicates().query('family_id in @templates').to_csv(out/'source_groups.csv',index=False)
    score=scores[f'{threshold:g}']
    write_json(out/'report.json',dict(selected_method=str(selected),frozen_threshold=threshold,
        candidate_count=len(catalog),test_scores=score,
        empirical_control_gate_passed=bool(score['null'] and score['null']['control_to_observed_range'][1] is not None and score['null']['control_to_observed_range'][1]<=0.05),
        release_status='experimental_unverified',
        limits=['Known northern coordinate groups only','No independent event review or validated origin times','No unknown-family completeness estimate','Control ratios are not true false-discovery rates']))
    print(json.dumps(score,indent=2),flush=True)

if __name__=='__main__': main()
