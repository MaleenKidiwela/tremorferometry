"""Evaluate the completed northern bank without repeating training or calibration."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing,json,shutil
import pandas as pd
from .pilot import ROOT,scan_date_worker,summarize_split,validate_config
from .training_artifact import load_training_artifact
from .reference import read_bostock
from .provenance import sha256,write_json,source_manifest
from .roster import validate_membership

def main():
    baseline=Path('results/north_full_family_bank_2012_v1');recovered=Path('results/north_full_family_training_recovered_v1')
    out=Path('results/north_full_family_frozen_evaluation_v1');out.mkdir(exist_ok=False)
    manifest=json.loads((baseline/'manifest.json').read_text());cfg=manifest['config']
    selection=json.loads((baseline/'selection_before_test.json').read_text());threshold=selection['threshold']
    if threshold is None or not selection['calibration_gate_passed']:raise ValueError('No calibration-selected threshold')
    templates,anchors,_=load_training_artifact(recovered,cfg,manifest['reference_sha256'])
    cfg['thresholds']=[threshold];validate_config(cfg)
    write_json(out/'manifest.json',dict(config=cfg,status='running_frozen_evaluation',source_code=source_manifest(ROOT/'codex_lfe'),template_sha256=sha256(recovered/'templates.npz'),roster=validate_membership(cfg['stations'],ROOT/cfg['station_roster'])))
    write_json(out/'selection_before_test.json',dict(threshold=threshold,calibration_selection_sha256=sha256(baseline/'selection_before_test.json'),training_artifact_sha256=sha256(recovered/'training_artifact.json'),dates=cfg['test_dates'],training_rerun=False,calibration_rerun=False,date_scope='Previously declared evaluation dates for this bank; not claimed untouched by all prior projects'))
    for name in ['training.json','templates.npz','template_qc.csv']:shutil.copyfile(recovered/name,out/name)
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    tasks=[(cfg,date,'test',templates,out,'numpy','cpu',2) for date in cfg['test_dates']]
    with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:completed=list(pool.map(scan_date_worker,tasks))
    pieces=[x for x,_ in completed];c,e,i=[pd.concat([p[k] for p in pieces],ignore_index=True) for k in range(3)]
    reference=read_bostock(ROOT/cfg['reference_csv'])
    scores,eligible=summarize_split(cfg,'test',c,e,i,reference,list(templates),anchors)
    write_json(out/'test_scores.json',scores);write_json(out/'input_manifest.json',[m for _,m in completed])
    eligible.to_parquet(out/'test_eligible_reference.parquet',index=False)
    observed=c[c.realization==0].copy();observed['quality_status']='experimental_unverified'
    observed.to_parquet(out/'test_experimental_candidates.parquet',index=False)
    write_json(out/'report.json',dict(status='complete_frozen_evaluation_not_final_catalog',scores=scores,threshold=threshold,training_rerun=False,calibration_rerun=False))
    print(json.dumps(scores,indent=2),flush=True)

if __name__=='__main__':main()
