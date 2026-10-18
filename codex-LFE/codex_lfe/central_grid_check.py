"""Fixed 15-second window-grid perturbation; diagnostic only, no threshold tuning."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
from .neural import infer_probabilities,picks_from_probabilities
from .waveforms import Waveform
from .provenance import sha256,write_json

def main():
    out=Path('results/central_grid_stability_v1');out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed diagnostic immutable')
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch,seisbench.models as sbm
    torch.set_num_threads(2)
    weights=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    write_json(out/'design_before_inference.json',dict(offset_s=15,thresholds=[.1,.3],matching_tolerance_s=.5,one_to_one=True,panels=['central_activity_picker_panel_v1','central_activity_matched_quiet_v1'],purpose='Window-grid stability, not phase accuracy or a production method selection',source_sha256=sha256(__file__),weights_sha256=sha256(weights.with_suffix('.pt'))))
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True);device='cuda' if torch.cuda.is_available() else 'cpu'
    rows=[]
    for panel in ['central_activity_picker_panel_v1','central_activity_matched_quiet_v1']:
        source=Path('results')/panel;dest=out/panel;dest.mkdir(exist_ok=True)
        baseline=pd.read_parquet(source/'picks.parquet')
        for file in sorted(source.glob('*_processed.npz')):
            station=file.name.removesuffix('_processed.npz');a=np.load(file);start=float(a['start_s']);fs=float(a['fs']);offset=int(15*fs)
            wave=Waveform(a['waveform'][:,offset:],a['valid'][:,offset:],start+15,fs,('Z','N','E'),station)
            prob,n=infer_probabilities(wave,model,device,batch_size=128)
            core_start=start+60;core_end=core_start+3600
            assert np.isfinite(prob[:,int((core_start-wave.start_s)*fs):int((core_end-wave.start_s)*fs)]).all()
            picks=picks_from_probabilities(prob,wave.start_s,fs,station,.05)
            picks=picks[(picks.time_s>=core_start)&(picks.time_s<core_end)]
            picks.to_parquet(dest/(station+'_picks.parquet'),index=False)
            np.savez_compressed(dest/(station+'_probabilities.npz'),probabilities=prob,start_s=wave.start_s,fs=fs)
            for phase in ['P','S']:
                for threshold in [.1,.3]:
                    x=baseline[(baseline.station==station)&(baseline.phase==phase)&(baseline.model_score>=threshold)].time_s.to_numpy()
                    y=picks[(picks.phase==phase)&(picks.model_score>=threshold)].time_s.to_numpy()
                    # Ascending-time greedy matching yields maximum cardinality for equal tolerance intervals.
                    x=np.sort(x);y=np.sort(y);i=j=0;res=[]
                    while i<len(x) and j<len(y):
                        if y[j]<x[i]-.5:j+=1
                        elif y[j]>x[i]+.5:i+=1
                        else:res.append(y[j]-x[i]);i+=1;j+=1
                    rows.append(dict(panel=panel,station=station,phase=phase,threshold=threshold,baseline_count=len(x),shifted_count=len(y),matched_count=len(res),median_absolute_delta_s=float(np.median(np.abs(res))) if res else None))
            print(panel,station,'completed',flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(out/'comparison.csv',index=False)
    grouped=frame.groupby(['panel','phase','threshold'])[['baseline_count','shifted_count','matched_count']].sum();grouped.to_csv(out/'summary.csv')
    write_json(out/'report.json',dict(status='complete_grid_stability_diagnostic',production_validated=False,phase_accuracy_measured=False,summary_sha256=sha256(out/'summary.csv')))
    print(grouped.to_string(),flush=True)
if __name__=='__main__':main()
