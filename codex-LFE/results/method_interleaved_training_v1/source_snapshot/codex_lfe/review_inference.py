"""Run a frozen picker on padded review waveforms, without reading human labels."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from .neural import infer_probabilities, picks_from_probabilities, prediction_coverage, verified_zne
from .waveforms import Waveform
from .provenance import write_json, sha256, environment


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',required=True)
    ap.add_argument('--model',required=True)
    ap.add_argument('--inventory',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--device',default='cpu')
    ap.add_argument('--minimum-score',type=float,default=.05)
    args=ap.parse_args();packet,out=Path(args.packet),Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    from .pilot import ROOT
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(ROOT/'cache/seisbench'))
    import torch
    import seisbench.models as sbm
    from obspy import read_inventory
    torch.set_num_threads(4)
    if args.device!='cpu' and not torch.cuda.is_available():raise RuntimeError('requested CUDA unavailable')
    meta=json.loads((packet/'manifest.json').read_text())
    # Template contains only initial pending status; reviewed exports are never read.
    windows=json.loads((packet/'review_template.json').read_text())['windows']
    archive=np.load(packet/'waveforms.npz');inv=read_inventory(args.inventory)
    model=sbm.LFEDetect.load(args.model,map_location='cpu',weights_only=True)
    all_picks=[];all_coverage=[];probabilities={};orientation={}
    for w in windows:
        x=archive[w['window_id']]
        wm=meta['waveform_inputs'][f"{meta['date']}/{w['station']}"]
        wave=Waveform(x,np.ones(x.shape,bool),w['start_s']-meta['context_s'],meta['fs'],tuple(wm['channels']),wm['site_id'])
        wave,orient=verified_zne(wave,inv);orientation[w['window_id']]=orient
        scores,nwindows=infer_probabilities(wave,model,args.device)
        picks=picks_from_probabilities(scores,wave.start_s,wave.fs,w['station'],args.minimum_score)
        picks['window_id']=w['window_id'];all_picks.append(picks)
        coverage=prediction_coverage(scores,wave.start_s,wave.fs,w['station']);coverage['window_id']=w['window_id'];all_coverage.append(coverage)
        probabilities[w['window_id']]=scores
    pd.concat(all_picks,ignore_index=True).to_parquet(out/'picks.parquet',index=False)
    pd.concat(all_coverage,ignore_index=True).to_parquet(out/'prediction_coverage.parquet',index=False)
    np.savez_compressed(out/'probabilities.npz',**probabilities)
    write_json(out/'manifest.json',{'status':'inference_only_no_accuracy_estimate','packet_manifest_sha256':sha256(packet/'manifest.json'),
        'packet_waveforms_sha256':sha256(packet/'waveforms.npz'),'packet_windows_sha256':sha256(packet/'review_template.json'),
        'model_config_sha256':sha256(args.model+'.json'),'weights_sha256':sha256(args.model+'.pt'),
        'inventory_sha256':sha256(args.inventory),'source_sha256':sha256(__file__),'minimum_score':args.minimum_score,
        'device':args.device,'environment':environment(),'orientation':orientation,'model_training':model.weights_docstring,
        'limitations':['Northern development packet only. No human labels scored.','Exact authors preprocessing and stitching are not reproduced; this runner uses the same preprocessing for both model variants.','Candidates below the minimum score are not retained in picks.parquet; probabilities are retained.']})
    print('Finished inference; no performance metrics until review is supplied.')

if __name__=='__main__':main()
