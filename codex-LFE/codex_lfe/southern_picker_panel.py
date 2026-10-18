"""Checkpointed southern picker transfer panel using existing development data."""
from pathlib import Path
import json
import os
import time
import numpy as np
import pandas as pd
from .neural import verified_zne,infer_probabilities,picks_from_probabilities,prediction_coverage
from .waveforms import load_day
from .provenance import sha256,write_json


def main():
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch
    import seisbench.models as sbm
    from obspy import read_inventory
    torch.set_num_threads(2)
    out=Path('results/southern_picker_transfer_v1');out.mkdir(exist_ok=True)
    cfg=json.loads(Path('results/southern_catalog_multifamily_2024_v1/manifest.json').read_text())['config']
    dates=['2024-02-01','2024-02-07','2024-02-14']
    weights=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    inventories=[Path(f'results/southern_picker_response_2024_v1/{n}.xml') for n in ['BK','PB']]
    design=dict(dates=dates,stations=cfg['stations'],date_rule='First, middle and last dates of existing development calibration interval; February 15–29 evaluation remains sealed.',
        model='jms',sources={str(p):sha256(p) for p in inventories+[weights.with_suffix('.pt'),weights.with_suffix('.json'),Path(__file__),Path('codex_lfe/neural.py')]},
        raw_root=cfg['raw_root'],minimum_score=.05,threshold_grid=[.05,.1,.2,.3,.5,.7,.9],
        purpose='New regional/instrument transfer diagnostic on cached southern data. Does not establish phase precision from published event timestamps.',
        inference_rule='Existing masked 1–8 Hz 20 Hz preprocessing, verified component gains and geographic rotation, fixed-grid central stitching, one-second peak guards; no new waveform downloads.')
    dp=out/'design_before_inference.json'
    if dp.exists() and json.loads(dp.read_text())!=design:raise ValueError('Changed design')
    write_json(dp,design)
    inventory=read_inventory(str(inventories[0]))+read_inventory(str(inventories[1]))
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True)
    device='cuda' if torch.cuda.is_available() else 'cpu'
    records=[]
    for date in dates:
        day=pd.Timestamp(date,tz='UTC')
        for station,meta in cfg['stations'].items():
            dst=out/f'{date}_{station}';dst.mkdir(exist_ok=True)
            done=dst/'manifest.json'
            if done.exists():
                receipt=json.loads(done.read_text())
                if all((dst/n).exists() and sha256(dst/n)==h for n,h in receipt['outputs'].items()):
                    records.append(receipt);continue
            raw=Path(cfg['raw_root'])/station/str(day.year)/f'{day.dayofyear:03d}.mseed'
            if not raw.exists():raise FileNotFoundError(f'No implicit acquisition allowed: {raw}')
            tick=time.monotonic()
            wave,wm=load_day(raw,Path('cache/southern_picker_waveforms'),meta['channels'],meta['site_id'],day.timestamp(),duration_s=86400)
            wave,orientation=verified_zne(wave,inventory)
            scores,n=infer_probabilities(wave,model,device,batch_size=256)
            picks=picks_from_probabilities(scores,wave.start_s,wave.fs,station,.05)
            coverage=prediction_coverage(scores,wave.start_s,wave.fs,station)
            picks.to_parquet(dst/'picks.parquet',index=False)
            coverage.to_parquet(dst/'coverage.parquet',index=False)
            np.savez_compressed(dst/'probabilities.npz',probabilities=scores,start_s=wave.start_s,fs=wave.fs)
            receipt=dict(date=date,station=station,status='inference_only_no_accuracy_claim',waveform=wm,orientation=orientation,n_windows=n,n_picks=len(picks),elapsed_s=time.monotonic()-tick,outputs={n:sha256(dst/n) for n in ['picks.parquet','coverage.parquet','probabilities.npz']})
            write_json(done,receipt);records.append(receipt)
            print(date,station,n,'windows',len(picks),'picks',round(receipt['elapsed_s'],1),'s',flush=True)
    write_json(out/'report.json',dict(status='complete_inference_only',device=device,station_days=len(records),records=records,production_picker_validated=False))

if __name__=='__main__':main()
