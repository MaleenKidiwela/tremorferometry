"""CPU date-parallel inference with verified orientations and explicit coverage."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
import yaml

from .neural import infer_probabilities,picks_from_probabilities,prediction_coverage,verified_zne
from .pilot import ROOT,load_waves,validate_config
from .provenance import sha256,write_json,environment
from .roster import validate_membership

MODEL=None
INVENTORY=None


def initialize(model_path,inventories,threads):
    global MODEL,INVENTORY
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(ROOT/'cache/seisbench'))
    import torch
    import seisbench.models as sbm
    from obspy import read_inventory
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    MODEL=sbm.LFEDetect.load(model_path,map_location='cpu',weights_only=True).eval()
    INVENTORY=read_inventory(inventories[0])
    for path in inventories[1:]: INVENTORY+=read_inventory(path)


def run_date(task):
    cfg,date,output=task;out=Path(output);manifest={};all_picks=[];all_coverage=[];rows=[]
    waves=load_waves(cfg,date,manifest)
    for station,wave in waves.items():
        start=time.monotonic()
        wave,orientation=verified_zne(wave,INVENTORY)
        scores,n=infer_probabilities(wave,MODEL,device='cpu',batch_size=64)
        picks=picks_from_probabilities(scores,wave.start_s,wave.fs,station,.05)
        coverage=prediction_coverage(scores,wave.start_s,wave.fs,station)
        picks['date']=date;coverage['date']=date
        all_picks.append(picks);all_coverage.append(coverage)
        name=f'{date}_{station}_probabilities.npz'
        np.savez_compressed(out/name,probabilities=scores)
        rows.append({'date':date,'station':station,'n_windows':n,'n_picks_min_score_0_05':len(picks),
                     'seconds':time.monotonic()-start,'orientation':orientation,
                     'probabilities_file':name,'probabilities_sha256':sha256(out/name)})
        print(f'[inference] {date} {station}: {n} windows, {len(picks)} peaks, {rows[-1]["seconds"]:.1f}s',flush=True)
    picks=pd.concat(all_picks,ignore_index=True);coverage=pd.concat(all_coverage,ignore_index=True)
    picks.to_parquet(out/f'{date}_picks.parquet',index=False)
    coverage.to_parquet(out/f'{date}_coverage.parquet',index=False)
    write_json(out/f'{date}_manifest.json',{'inputs':manifest,'stations':rows})
    return date,rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--model',default='../lfe_catalog_repair/results/external/other/seisbench_lfedetect/cascadia')
    ap.add_argument('--inventory',nargs='+',required=True)
    ap.add_argument('--workers',type=int,default=6);ap.add_argument('--torch-threads',type=int,default=2)
    args=ap.parse_args();out=Path(args.output)
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load(Path(args.config).read_text());validate_config(cfg)
    roster=validate_membership(cfg['stations'],cfg['station_roster'])
    dates=sorted(set(cfg['calibration_dates']+cfg['test_dates']))
    files=[args.config,args.model+'.pt',args.model+'.json',__file__,str(Path(__file__).with_name('neural.py'))]+args.inventory
    manifest={'status':'running','sources':{str(Path(p).resolve()):sha256(p) for p in files},
        'config':cfg,'dates':dates,'roster':roster,'environment':environment(),
        'workers':args.workers,'torch_threads_per_worker':args.torch_threads,
        'model_name':'cascadia','minimum_export_score':.05,
        'preprocessing':'masked 1-8 Hz at 20 Hz, verified component gains and ZNE rotation, central-window stitching',
        'peak_rule':'local maxima separated by at least 1 s, +/-1 s coverage guard',
        'limitation':'This is a corrected processing arm, not an isolated rotation ablation or exact legacy reproduction.'}
    write_json(out/'manifest.json',manifest);tick=time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn'),
                             initializer=initialize,initargs=(args.model,args.inventory,args.torch_threads)) as pool:
        completed=list(pool.map(run_date,[(cfg,d,str(out)) for d in dates]))
    picks=pd.concat([pd.read_parquet(out/f'{d}_picks.parquet') for d,_ in completed],ignore_index=True)
    coverage=pd.concat([pd.read_parquet(out/f'{d}_coverage.parquet') for d,_ in completed],ignore_index=True)
    picks.to_parquet(out/'picks.parquet',index=False);coverage.to_parquet(out/'coverage.parquet',index=False)
    # Explicit adapter into the audited association runner's saved-pick schema.
    adapted=pd.DataFrame({'set':'interleaved2010','model':'cascadia_corrected',
        'net':picks.station.str.split('.').str[0],'sta':picks.station.str.split('.').str[1],
        'day':picks.date,'phase':picks.phase,'time':pd.to_datetime(picks.time_s,unit='s',utc=True),
        'prob':picks.model_score,'thr_first':.05})
    adapted.to_parquet(out/'association_input.parquet',index=False)
    manifest.update(status='complete_inference_only',elapsed_s=time.monotonic()-tick,n_picks=len(picks),
                    outputs_sha256={n:sha256(out/n) for n in ['picks.parquet','coverage.parquet','association_input.parquet']})
    write_json(out/'manifest.json',manifest)
    print(f'Finished {len(dates)*len(cfg["stations"])} station-days in {manifest["elapsed_s"]:.1f}s',flush=True)


if __name__=='__main__': main()
