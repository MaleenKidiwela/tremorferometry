"""Audit template repeatability across separated training dates, without lag fitting."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .pilot import ROOT, load_waves
from .reference import read_bostock
from .templates import cut, stack_windows, vector_cc
from .provenance import write_json, sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((run/'manifest.json').read_text())['config']
    dates = sorted(cfg['train_dates'])
    if len(dates) < 4:
        raise ValueError('at least four training dates required')
    early, late = dates[:len(dates)//2], dates[len(dates)//2:]
    write_json(out/'design.json', {'early_dates':early,'late_dates':late,
        'minimum_windows_per_half':10,'time_lag_search':False,
        'purpose':'Post-fit temporal stability audit of training templates; not independent generalization evidence.'})
    ref = read_bostock((ROOT/cfg['reference_csv']).resolve())
    qc = pd.read_csv(run/'template_qc.csv')
    bank, manifests = {}, {}
    for date in dates:
        print('[temporal QC]',date,flush=True)
        waves = load_waves(cfg,date,manifests)
        records = ref[ref.date == date]
        for row in qc.itertuples():
            if not row.usable:
                continue
            windows=[]
            for t in records.loc[records.family_id == row.family_id,'time_s']:
                w=cut(waves[row.station],t-10,int(70*cfg['fs']))
                if w is not None:
                    windows.append(w)
            if windows:
                bank[row.family_id,row.station,date]=(len(windows),stack_windows(windows))
    rows=[]; arrays={}
    for row in qc.itertuples():
        if not row.usable:
            continue
        halves=[]; counts=[]
        for group in (early,late):
            items=[bank[row.family_id,row.station,d] for d in group if (row.family_id,row.station,d) in bank]
            count=sum(n for n,_ in items); counts.append(count)
            halves.append(sum(n*s for n,s in items)/count if count else None)
        first=int(round((row.template_start_from_published_s+10)*cfg['fs']))
        last=first+int(round(cfg['template_s']*cfg['fs']))
        enough=min(counts)>=10
        cc=vector_cc(halves[0][:,first:last],halves[1][:,first:last]) if enough else None
        rows.append({'family_id':row.family_id,'station':row.station,'early_n':counts[0],'late_n':counts[1],
            'alternating_event_cc':row.split_stack_cc,'separated_dates_cc':cc,
            'passes_same_cc_bar':bool(cc>=cfg['min_training_split_cc']) if enough else False,
            'status':'measured' if enough else 'insufficient_windows'})
        if enough:
            for label,s in zip(('early','late'),halves):
                arrays[f'{row.family_id}__{row.station}__{label}']=s[:,first:last]
    frame=pd.DataFrame(rows); frame.to_csv(out/'temporal_qc.csv',index=False)
    np.savez(out/'separated_stacks.npz',**arrays)
    write_json(out/'manifest.json',{'source_sha256':sha256(__file__),
        'template_qc_sha256':sha256(run/'template_qc.csv'),'reference_sha256':sha256(ROOT/cfg['reference_csv']),
        'inputs':manifests,'n_measured':int((frame.status=='measured').sum()),
        'n_pass_same_bar':int(frame.passes_same_cc_bar.sum()),
        'limitation':'Peak windows were previously selected using all training dates. The two stacks use disjoint dates, but window selection is shared. No test or calibration data used.'})

if __name__=='__main__':
    main()
