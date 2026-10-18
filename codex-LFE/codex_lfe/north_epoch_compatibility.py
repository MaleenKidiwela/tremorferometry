"""Metadata-only compatibility of the frozen northern bank across station epochs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .provenance import write_json,sha256


def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/north_bank_epoch_compatibility_v1';out.mkdir(exist_ok=False)
    baseline=root/'results/north_full_family_bank_2012_v1/manifest.json'
    artifact=root/'results/north_full_family_training_recovered_v1/templates.npz'
    metadata=root/'results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv'
    cfg=json.loads(baseline.read_text())['config'];archive=np.load(artifact,allow_pickle=False)
    families=sorted({k.split('__')[0] for k in archive.files})
    members={f:[s for s in cfg['stations'] if f'{f}__{s}__template' in archive] for f in families}
    e=pd.read_csv(metadata);e['id']=e.network+'.'+e.station;e=e[e.id.isin(cfg['stations'])].copy()
    e['lo']=pd.to_datetime(e.start_utc,utc=True);e['hi']=pd.to_datetime(e.end_utc_exclusive,utc=True)
    boundaries=sorted(set(e.lo)|set(e.hi));rows=[];details=[]
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        active=e[(e.lo<=lo)&(e.hi>=hi)]
        if active.id.duplicated().any():raise ValueError('Overlapping selected epochs')
        channels=active.set_index('id').channels.to_dict();counts=[]
        for f,ss in members.items():
            current=[s for s in ss if s in channels]
            exact=[s for s in current if set(channels[s].split(','))==set(cfg['stations'][s]['channels'])]
            sites=lambda a:len({cfg['stations'][s]['site_id'] for s in a})
            row=dict(family_id=f,start_utc=lo.isoformat(),end_utc_exclusive=hi.isoformat(),
                     original_sites=sites(ss),current_station_identity_sites=sites(current),exact_triplet_sites=sites(exact),
                     every_original_triplet_present=len(exact)==len(ss),
                     anchor_triplet_present=cfg['anchor_station'] in exact,
                     three_sites_after_unvalidated_station_subset=sites(current)>=cfg['min_sites'])
            counts.append(row);details.append(row)
        rows.append(dict(start_utc=lo.isoformat(),end_utc_exclusive=hi.isoformat(),
                         families=len(families),complete_original_triplet_sets=sum(r['every_original_triplet_present'] for r in counts),
                         families_with_at_least_three_identity_sites=sum(r['three_sites_after_unvalidated_station_subset'] for r in counts),
                         families_with_original_anchor_triplet=sum(r['anchor_triplet_present'] for r in counts)))
    pd.DataFrame(rows).to_csv(out/'epoch_summary.csv',index=False)
    pd.DataFrame(details).to_parquet(out/'family_epoch_compatibility.parquet',index=False)
    latest=rows[-1]
    write_json(out/'report.json',dict(status='complete_metadata_compatibility_audit',latest=latest,
        production_validated=False,limitations=['Station metadata is not waveform availability.',
        'Matching channel IDs do not prove unchanged instrument response or orientation.',
        'Removing stations changes the detector and requires separate calibration.',
        'The original runner requests all configured stations and cannot run unchanged when one is absent.',
        'Cross-channel transfer must be validated before applying a historical bank to current channels.'],
        inputs={str(p):sha256(p) for p in [baseline,artifact,metadata,Path(__file__)]}))
    print(json.dumps(latest,indent=2))

if __name__=='__main__':main()
