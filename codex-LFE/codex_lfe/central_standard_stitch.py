"""Isolate library overlap averaging on existing processed panels; no new filtering."""
from pathlib import Path
import json,os,inspect
import numpy as np
import pandas as pd
from obspy import Stream,Trace,UTCDateTime
from .provenance import write_json,sha256
from .neural import picks_from_probabilities

def main():
    out=Path('results/central_standard_stitch_v1');out.mkdir(exist_ok=True)
    if (out/'report.json').exists():raise ValueError('Completed experiment immutable')
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch,seisbench,seisbench.models as sbm
    torch.set_num_threads(2)
    weights=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True).eval().to('cuda' if torch.cuda.is_available() else 'cpu')
    model.filter_args=None;model.filter_kwargs=None
    write_json(out/'design_before_inference.json',dict(method='Installed SeisBench annotate overlap=600 stacking=avg blinding=(0,0)',preprocessing='Reuse saved response-corrected ZNE 1–8 Hz 20 Hz; disable additional model filtering',grid_offsets_s=[0,15],thresholds=[.1,.3],purpose='Isolate stitching sensitivity; no accuracy claim; prior E2 used library stitching on different recordings',seisbench_version=seisbench.__version__,source_sha256=sha256(__file__),weights_sha256=sha256(weights.with_suffix('.pt'))))
    for panel in ['central_activity_picker_panel_v1','central_activity_matched_quiet_v1']:
        source=Path('results')/panel
        for offset in [0,15]:
            dest=out/(panel+f'_offset{offset}');dest.mkdir(exist_ok=True);all_picks=[]
            (dest/'stations.csv').write_bytes((source/'stations.csv').read_bytes());(dest/'design_before_acquisition.json').write_bytes((source/'design_before_acquisition.json').read_bytes())
            for file in sorted(source.glob('*_processed.npz')):
                station=file.name.removesuffix('_processed.npz');a=np.load(file);start=float(a['start_s']);fs=float(a['fs']);skip=int(offset*fs)
                if not a['valid'].all():raise ValueError('This bounded comparison requires complete processed inputs')
                net,sta=station.split('.')
                stream=Stream([Trace(x[skip:].copy(),header=dict(network=net,station=sta,channel='HH'+c,sampling_rate=fs,starttime=UTCDateTime(start+offset))) for c,x in zip('ZNE',a['waveform'])])
                ann=model.annotate(stream,overlap=600,stacking='avg',blinding=(0,0),strict=True,batch_size=128)
                prob=np.full_like(a['probabilities'],np.nan)
                for k,phase in enumerate('PS'):
                    selected=ann.select(channel='LFEDetect_'+phase)
                    if len(selected)!=1:raise ValueError('Unexpected annotation segmentation')
                    tr=selected[0];lo=round((float(tr.stats.starttime)-start)*fs);prob[k,lo:lo+len(tr.data)]=tr.data
                assert np.isfinite(prob[:,1200:73200]).all()
                picks=picks_from_probabilities(prob,start,fs,station,.05);picks=picks[(picks.time_s>=start+60)&(picks.time_s<start+3660)];all_picks.append(picks)
                np.savez_compressed(dest/(station+'_processed.npz'),probabilities=prob,start_s=start,fs=fs)
                print(panel,offset,station,'complete',flush=True)
            pd.concat(all_picks,ignore_index=True).to_parquet(dest/'picks.parquet',index=False)
    write_json(out/'report.json',dict(status='complete_stitching_comparison_inputs',production_validated=False))
if __name__=='__main__':main()
