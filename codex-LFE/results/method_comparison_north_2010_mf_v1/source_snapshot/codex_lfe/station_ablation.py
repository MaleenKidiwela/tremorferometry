"""Calibration-only station ablation from separated-training-date template QC."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .pilot import ROOT, scan_date, summarize_split
from .reference import read_bostock
from .provenance import write_json, sha256, source_manifest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',required=True)
    ap.add_argument('--qc',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    run,out=Path(args.run),Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((run/'manifest.json').read_text())['config']
    training=json.loads((run/'training.json').read_text())
    qc=pd.read_csv(args.qc)
    remove=qc[(qc.status=='measured') & (qc.separated_dates_cc < cfg['min_training_split_cc'])]
    excluded=set(zip(remove.family_id,remove.station))
    archive=np.load(run/'templates.npz')
    templates={}
    for family in training['families_usable']:
        members={}
        for station in cfg['stations']:
            key=f'{family}__{station}'
            if key+'__template' in archive and (family,station) not in excluded:
                members[station]=(archive[key+'__template'],float(archive[key+'__offset_s']))
        if cfg['anchor_station'] not in members or len({cfg['stations'][s]['site_id'] for s in members}) < cfg['min_sites']:
            raise ValueError('ablation would invalidate a source group; redesign explicitly')
        templates[family]=members
    design={'status':'running','scope':'calibration_only_no_test_scan','config':cfg,
        'source_code':source_manifest(ROOT/'codex_lfe'),'parent_manifest_sha256':sha256(run/'manifest.json'),
        'templates_sha256':sha256(run/'templates.npz'),'training_qc_sha256':sha256(args.qc),
        'excluded':remove[['family_id','station','separated_dates_cc']].to_dict('records'),
        'rule':'Remove measured separated-date QC failures only; retain groups with insufficient temporal coverage. Preserve original templates, clocks, thresholds, and controls.'}
    write_json(out/'design_before_scan.json',design)
    started=time.monotonic(); manifests={}
    def scan(date):
        return scan_date(cfg,date,'calibration',templates,out,manifests,'numpy','cpu')
    with ThreadPoolExecutor(max_workers=3) as pool:
        pieces=list(pool.map(scan,cfg['calibration_dates']))
    c,e,i=[pd.concat([p[k] for p in pieces],ignore_index=True) for k in range(3)]
    reference=read_bostock(Path(cfg['reference_csv']))
    summary,eligible=summarize_split(cfg,'calibration',c,e,i,reference,list(templates),training['anchors_from_published_s'])
    eligible.to_parquet(out/'calibration_eligible_reference.parquet',index=False)
    write_json(out/'calibration_scores.json',summary)
    write_json(out/'input_manifest.json',manifests)
    old=json.loads((run/'calibration_scores.json').read_text()); comparison=[]
    for threshold,result in summary.items():
        for label,scores in [('baseline',old[threshold]),('training_qc_ablation',result)]:
            comparison.append({'threshold':float(threshold),'variant':label,
                'n_detected':scores['reference']['n_detections'],
                'n_correct_group_matches':scores['reference']['n_matched_same_family'],
                'eligible_reference':scores['reference_eligible'],
                'max_control_ratio':scores['null']['control_to_observed_range'][1] if scores['null'] else None})
    pd.DataFrame(comparison).to_csv(out/'comparison.csv',index=False)
    write_json(out/'completion.json',{'status':'complete_calibration_experiment_not_test_evidence',
        'elapsed_s':time.monotonic()-started,'test_started':False})

if __name__=='__main__':
    main()
