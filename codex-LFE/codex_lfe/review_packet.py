"""Create a small development-only arrival-review packet from target-roster data."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from .pilot import ROOT, load_waves, date_start
from .reference import read_bostock
from .roster import validate_membership
from .templates import cut
from .provenance import write_json, sha256


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',required=True)
    args=ap.parse_args();out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load((ROOT/'configs/map_network_north_calibration.yaml').read_text())
    cfg['stations']={s:m for s,m in cfg['stations'].items() if s in ['CN.PGC','CN.SNB']}
    membership=validate_membership(cfg['stations'],ROOT/cfg['station_roster'])
    date='2012-09-08';manifests={};waves=load_waves(cfg,date,manifests)
    reference=read_bostock(ROOT/cfg['reference_csv']);reference=reference[reference.date==date]
    rng=np.random.default_rng(20260912);starts=[];seeds=[];windows=[];traces={};keys=[]
    duration=60.;context=30.;n=int((duration+2*context)*cfg['fs'])
    for sampling in ['reference_enriched','random']:
        accepted=0
        if sampling=='reference_enriched':
            # Select from a public source table, never from competing picker output.
            pool=reference.iloc[rng.permutation(len(reference))]
            candidates=[(float(x.time_s+rng.uniform(-5,5)),x.reference_id) for x in pool.itertuples()]
        else:
            candidates=[(float(date_start(date)+rng.uniform(180,86400-180)),None) for _ in range(1000)]
        for start,source_id in candidates:
            if any(abs(start-s)<180 for s in starts):continue
            cuts={s:cut(w,start-context,n) for s,w in waves.items()}
            if any(v is None for v in cuts.values()):continue
            block=f'block_{len(starts):03d}';starts.append(start)
            for station,w in waves.items():
                window_id=f'window_{len(windows):03d}'
                row={'window_id':window_id,'station':station,'start_s':start,'end_s':start+duration,
                    'split':'development','sampling':sampling,'review_status':'pending','reviewer':'',
                    'block_id':block,'notes':''}
                windows.append(row);traces[window_id]=cuts[station]
                keys.append({'window_id':window_id,'sampling':sampling,'source_reference_id':source_id,'block_id':block})
            accepted+=1
            if accepted==6:break
        if accepted!=6:raise RuntimeError('insufficient valid review windows')
    order=rng.permutation(len(windows));windows=[windows[i] for i in order]
    write_json(out/'review_template.json',{'windows':windows,'labels':[]})
    write_json(out/'sampling_key.json',keys)
    np.savez_compressed(out/'waveforms.npz',**traces)
    packet=[]
    for row in windows:
        # Hide source selection, dates and station identity during initial picking.
        x=traces[row['window_id']];scale=float(np.quantile(abs(x),.99)) or 1.
        packet.append({'window_id':row['window_id'],'start_s':row['start_s'],'duration_s':duration,
            'channels':[c.split('.')[-1] for c in waves[row['station']].channels],
            'x':np.round(x/scale,5).tolist()})
    manifest={'status':'unreviewed_development_packet_not_ground_truth','date':date,'seed':20260912,
        'n_windows':len(windows),'n_paired_time_blocks':len(starts),'scored_window_s':duration,
        'context_s':context,'fs':cfg['fs'],'station_roster':membership,'waveform_inputs':manifests,
        'source_code_sha256':sha256(__file__),'reference_sha256':sha256(ROOT/cfg['reference_csv']),
        'sampling':'6 reference-enriched and 6 random times, each at PGC and SNB; randomized display order.',
        'limitations':['Development data previously used in matcher calibration; not held-out evidence.',
                       'Northern broadband examples only; central/south and boreholes remain required.',
                       'Original source timestamps guide sampling only; no P/S labels inferred.',
                       'Random times are not assumed negative; all windows require review.',
                       'Whole paired time blocks must stay together in later splits.']}
    write_json(out/'manifest.json',manifest)
    template=(ROOT/'review/picker_review.html').read_text()
    payload={'packet':packet,'windows':windows,'fs':cfg['fs'],'context':context,'key':sha256(out/'waveforms.npz')}
    (out/'review.html').write_text(template.replace('__PAYLOAD__',json.dumps(payload,separators=(',',':'))))
    print(f'{len(windows)} unlabeled windows: {out}/review.html')

if __name__=='__main__':main()
