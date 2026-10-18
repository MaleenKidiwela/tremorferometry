"""Calibration-only four-second template experiment with unchanged site controls."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .pilot import ROOT, scan_date_worker, summarize_split
from .provenance import sha256, source_manifest, write_json
from .reference import read_bostock
from .roster import validate_membership
from .site_weighting import operating_point

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',required=True)
    ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args(); run,out=map(Path,[args.run,args.output])
    if not 1<=args.workers<=8: raise ValueError('workers must be 1..8')
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((run/'manifest.json').read_text())['config']
    roster=validate_membership(cfg['stations'],ROOT/cfg['station_roster'])
    templates = {}
    archive = {}
    with np.load(run/'templates.npz') as bank:
        for family in json.loads((run/'training.json').read_text())['families_usable']:
            templates[family] = {}
            for station in cfg['stations']:
                key = f'{family}__{station}'
                if key+'__template' not in bank:
                    continue
                original = bank[key+'__template']
                # Frozen eight-second template begins two seconds before training peak.
                first, last = int(1.5*cfg['fs']), int(5.5*cfg['fs'])
                if original.shape[1] != int(8*cfg['fs']):
                    raise ValueError('Expected eight-second baseline template')
                template = original[:, first:last].copy()
                offset = float(bank[key+'__offset_s']) + first/cfg['fs']
                templates[family][station] = (template, offset)
                archive[key+'__template'] = template
                archive[key+'__offset_s'] = np.asarray(offset)
    cfg['template_s'] = 4.0
    cfg.pop('site_weights', None)
    np.savez_compressed(out/'templates.npz', **archive)
    write_json(out/'design_before_scan.json', {
        'status':'frozen_before_calibration_scan', 'scope':'calibration_only',
        'config':cfg, 'roster':roster,
        'rule':'Crop frozen training templates to peak -0.5 through +3.5 seconds; preserve anchor clock, source/station membership, equal site weights, support gate and controls.',
        'sources':{str(p.resolve()):sha256(p) for p in [run/'templates.npz',run/'manifest.json',run/'calibration_scores.json',Path(__file__)]},
        'source_code':source_manifest(ROOT/'codex_lfe'),
        'decision_rule':'Advance only if selected-threshold same-family matches exceed baseline under the unchanged engineering gate and identical eligibility; fresh-date validation required.'})
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
    identical_intervals = intervals[keys].sort_values(keys).reset_index(drop=True).equals(baseline_intervals[keys].sort_values(keys).reset_index(drop=True))
    comparison=[]
    for threshold,result in summary.items():
        for label,scores in [('baseline',old[threshold]),('short_packet',result)]:
            comparison.append({'threshold':float(threshold),'variant':label,
                'n_detected':scores['reference']['n_detections'],
                'n_correct_group_matches':scores['reference']['n_matched_same_family'],
                'eligible_reference':scores['reference_eligible'],
                'control_counts':json.dumps(scores['null']['control_counts']),
                'max_control_ratio':scores['null']['control_to_observed_range'][1]})
    pd.DataFrame(comparison).to_csv(out/'comparison.csv',index=False)
    selected=operating_point(summary,cfg);baseline=operating_point(old,cfg)
    decision={'status':'complete_calibration_experiment_no_evaluation_scan','elapsed_s':time.monotonic()-tick,
        'baseline_operating_threshold':baseline,'short_packet_operating_threshold':selected,
        'shared_reference_count':len(eligible),'identical_waveform_intervals':identical_intervals,
        'baseline_matches_at_selected':old[f'{baseline:g}']['reference']['n_matched_same_family'] if baseline else None,
        'short_packet_matches_at_selected':summary[f'{selected:g}']['reference']['n_matched_same_family'] if selected else None,
        'interpretation':'Compare calibration recovery under the same engineering control gate. Empirical control ratios are not true false-discovery rates. No evaluation scan or automatic production adoption.'}
    decision['advance_to_fresh_validation'] = bool(selected is not None and decision['short_packet_matches_at_selected'] > decision['baseline_matches_at_selected'])
    write_json(out/'decision.json',decision)
    print(pd.DataFrame(comparison).to_string(index=False),flush=True)
    print(json.dumps(decision),flush=True)


if __name__=='__main__': main()
