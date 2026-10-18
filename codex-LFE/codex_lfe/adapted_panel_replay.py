"""Replay cached natural waveforms with the experimental adapted P/S model."""
from pathlib import Path
import json,os,shutil
import numpy as np
import pandas as pd
from .neural import Waveform,infer_probabilities,picks_from_probabilities,prediction_coverage
from .provenance import sha256,write_json


def main():
    root=Path(__file__).resolve().parents[1];os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(root/'cache/seisbench'))
    import torch,seisbench.models as sbm
    torch.set_num_threads(2);device='cuda' if torch.cuda.is_available() else 'cpu'
    original=root.parent/'lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms'
    state=root/'results/north_picker_adaptation_v1/adapted_state_dict.pt'
    model=sbm.LFEDetect.load(str(original),map_location='cpu',weights_only=True)
    model.load_state_dict(torch.load(state,map_location='cpu',weights_only=True));model.eval()
    for panel in ['quiet','active']:
        source=root/f'results/north_recent_{panel}_picker_v1';out=root/f'results/north_recent_adapted_{panel}_v1';out.mkdir(exist_ok=False)
        design=json.loads((source/'design_before_acquisition.json').read_text());start=pd.Timestamp(design['start_utc']).timestamp()
        write_json(out/'design_before_acquisition.json',design|dict(purpose='Experimental adapted-model replay of already acquired natural data',acquisition_performed=False,waveform_source=str(source),adapted_weights_sha256=sha256(state)))
        shutil.copyfile(source/'stations.csv',out/'stations.csv');picks=[];coverage=[];inputs={}
        for path in sorted(source.glob('*_processed.npz')):
            station=path.name.removesuffix('_processed.npz')
            with np.load(path) as a:wave=Waveform(a['waveform'],a['valid'],float(a['start_s']),float(a['fs']),('Z','N','E'),station)
            probabilities,n=infer_probabilities(wave,model,device,batch_size=128)
            p=picks_from_probabilities(probabilities,wave.start_s,wave.fs,station,.05);p=p[(p.time_s>=start)&(p.time_s<start+3600)]
            picks.append(p);coverage.append(prediction_coverage(probabilities,wave.start_s,wave.fs,station))
            np.savez_compressed(out/path.name,probabilities=probabilities,start_s=wave.start_s,fs=wave.fs)
            inputs[str(path)]=sha256(path);print(panel,station,n,'windows',len(p),'peaks >=.05',flush=True)
        pd.concat(picks,ignore_index=True).to_parquet(out/'picks.parquet',index=False)
        pd.concat(coverage,ignore_index=True).to_parquet(out/'coverage.parquet',index=False)
        write_json(out/'report.json',dict(status='complete_adapted_natural_waveform_replay',production_validated=False,
            phase_accuracy_measured=False,raw_bytes_downloaded=0,source_sha256=sha256(__file__),weights_sha256=sha256(state),waveform_inputs=inputs,
            limitations=['The comparison-hour spectra supplied synthetic training noise; this is not an independent field validation of false-positive rate.','The activity hour was not used for adaptation, but has no event-level phase truth.','Processed NPZ files contain predictions only; original waveforms remain in waveform_source.']))

if __name__=='__main__':main()
