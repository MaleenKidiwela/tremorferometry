"""Reuse the established loss diagnostic on the expanded bank's calibration split."""
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd
from .miss_diagnosis import matcher_date
from .provenance import write_json,sha256


def main():
    root=Path(__file__).resolve().parents[1]
    run=root/'results/north_full_family_bank_2012_v1'
    out=root/'results/north_full_bank_calibration_loss_v1'
    out.mkdir(exist_ok=False)
    cfg=json.loads((run/'manifest.json').read_text())['config']
    threshold=json.loads((run/'selection_before_test.json').read_text())['threshold']
    write_json(out/'design.json',dict(status='running',split='calibration',dates=cfg['calibration_dates'],threshold=threshold,
        purpose='Locate loss before selecting a new development method; no threshold changes or holdout rescoring'))
    tasks=[(cfg,str(run),date,threshold) for date in cfg['calibration_dates']]
    with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
        result=list(pool.map(matcher_date,tasks))
    frame=pd.DataFrame([row for rows,_ in result for row in rows])
    frame.to_parquet(out/'reference_scores.parquet',index=False)
    frame.groupby(['date','status']).size().rename('n_references').reset_index().to_csv(out/'loss_by_date.csv',index=False)
    frame.groupby(['family_id','status']).size().rename('n_references').reset_index().to_csv(out/'loss_by_family.csv',index=False)
    write_json(out/'inputs.json',[m for _,m in result])
    report=dict(status='complete_calibration_loss_diagnostic',n_references=len(frame),
        counts=frame.status.value_counts().to_dict(),threshold=threshold,production_validated=False,
        inputs={str(p):sha256(p) for p in [Path(__file__),root/'codex_lfe/miss_diagnosis.py',run/'templates.npz',run/'calibration_eligible_reference.parquet',run/'manifest.json']})
    write_json(out/'report.json',report);print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
