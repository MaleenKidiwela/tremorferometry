"""Complement the consecutive late-month test with unused interleaved dates.

All four gaps within the original September 7-20 design are included, with no
selection by detector results and no change to the frozen threshold or templates.
"""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import json,multiprocessing
import numpy as np
import pandas as pd
from .pilot import ROOT
from .training_artifact import load_training_artifact
from .north_extended_validation import scan
from .provenance import sha256,write_json


def main():
    out=ROOT/'results/north_interleaved_family_validation_v1';out.mkdir(exist_ok=False)
    base=ROOT/'results/north_full_family_bank_2012_v1';training=ROOT/'results/north_full_family_training_recovered_v1'
    manifest=json.loads((base/'manifest.json').read_text());cfg=manifest['config']
    templates,anchors,_=load_training_artifact(training,cfg,manifest['reference_sha256'])
    used=set(cfg['train_dates']+cfg['calibration_dates']+cfg['test_dates'])
    dates=[str(d.date()) for d in pd.date_range(min(cfg['train_dates']),max(cfg['test_dates'])) if str(d.date()) not in used]
    if dates!=['2012-09-09','2012-09-12','2012-09-15','2012-09-17']:raise ValueError('Unexpected unused-date design')
    for date in dates:
        dt=pd.Timestamp(date)
        for station in cfg['stations']:
            if not (Path(cfg['raw_root'])/station/str(dt.year)/f'{dt.dayofyear:03d}.mseed').exists():raise ValueError('Missing cached input')
    write_json(out/'design_before_scan.json',dict(status='running_interleaved_date_validation',production_validated=False,
        dates=dates,selection='Every unused calendar date within the original September 7-20 train/calibration/evaluation span; no detector-result or activity-count selection.',
        split_note='Excluded from fitting and threshold selection of this bank; not claimed untouched by other projects.',
        purpose='Complement late-month low-reference exposure without replacing it or tuning on its results.',
        threshold_mad=9.,controls=16,workers=2,raw_bytes_to_download=0,
        clock_policy='Explicit family-specific PGC peak clocks, no cross-family event merging.',
        exposure='Identical per-family observed/control support intersection.',
        original_evaluation_changed=False,config=cfg,
        sources={str(p):sha256(p) for p in [Path(__file__),ROOT/'codex_lfe/north_extended_validation.py',training/'templates.npz',training/'training_artifact.json',base/'selection_before_test.json']}))
    (out/'runner_source.py').write_bytes(Path(__file__).read_bytes());done=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        jobs=[pool.submit(scan,(cfg,templates,anchors,date,str(out))) for date in reversed(dates)]
        for future in as_completed(jobs):
            done.append(future.result());write_json(out/'progress.json',dict(completed_dates=sorted(r['date'] for r in done),total_dates=len(dates)))
    done.sort(key=lambda r:r['date']);observed=sum(r['n_observed_family_rows'] for r in done)
    eligible=sum(r['n_eligible_references'] for r in done);matched=sum(r['n_family_matches'] for r in done);controls=np.sum([r['control_counts'] for r in done],axis=0)
    pd.DataFrame(done).to_csv(out/'daily_summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_interleaved_family_validation',production_validated=False,
        n_observed_family_hypotheses=observed,n_eligible_references=eligible,n_family_matches=matched,
        within_family_reference_recall=matched/eligible if eligible else None,control_counts=controls.tolist(),
        maximum_control_ratio=float(controls.max()/observed) if observed else None,
        descriptive_control_check_passed=bool(observed>=20 and controls.max()/observed<=.05),
        unique_event_count=None,raw_bytes_downloaded=0,days=done))
    print((out/'report.json').read_text(),flush=True)


if __name__=='__main__':main()
