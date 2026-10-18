"""Published labeled-stack timing diagnostic, not continuous-event accuracy."""
from pathlib import Path
import json
import os
import h5py
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from .neural import infer_probabilities, picks_from_probabilities, prediction_coverage
from .waveforms import Waveform
from .provenance import sha256, write_json


def main():
    os.environ.setdefault('SEISBENCH_CACHE_ROOT', str(Path('cache/seisbench').resolve()))
    import torch
    import seisbench.models as sbm
    torch.set_num_threads(2)
    source = Path('../lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia')
    model_path = Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    out = Path('results/picker_labeled_stacks_jms_v2')
    out.mkdir(exist_ok=True)
    (out/'checkpoints').mkdir(exist_ok=True)
    meta = pd.read_csv(source/'metadata.csv')
    design = dict(status='development_diagnostic', model='jms',
        model_training=json.loads(model_path.with_suffix('.json').read_text())['docstring'],
        thresholds=[.05,.1,.2,.3,.5,.7,.9], tolerances_s=[.25,.5,1.0],
        label_columns=['trace_P_arrival_sample','trace_S_arrival_sample'],
        labels_not_used=['trace_P_predicted_arrival_sample','trace_S_predicted_arrival_sample'],
        component_rule='Original Z12 interpreted as ZNE following existing e4_polarization_check; horizontal permutation sensitivity also reported, never selected by score.',
        station_rule='Retain original PO network labels; no unsupported remapping to modern roster stations.',
        limitations=['Stack labels are not natural continuous-event truth.', 'Missing phase labels are not negatives.',
                    'Original dataset split retained, but this is an already accessed development resource.',
                    'Unmatched peaks are not a measured false-pick rate.', 'Label uncertainties absent; residuals are relative to nominal published picks.'],
        sources={str(p):sha256(p) for p in [source/'metadata.csv',source/'waveforms.hdf5',model_path.with_suffix('.pt'),model_path.with_suffix('.json'),Path(__file__)]})
    dp=out/'design_before_inference.json'
    if dp.exists() and json.loads(dp.read_text()) != design:
        raise ValueError('Cannot resume changed design or source')
    write_json(dp,design)
    model = sbm.LFEDetect.load(str(model_path),map_location='cpu',weights_only=True).eval()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    sos = butter(4,[1,8],btype='bandpass',fs=20,output='sos')
    rows=[]
    parity=json.loads((out/'cpu_gpu_parity.json').read_text()) if (out/'cpu_gpu_parity.json').exists() else None
    with h5py.File(source/'waveforms.hdf5') as h:
        assert h['data_format/component_order'][()] == b'Z12'
        for count,r in enumerate(meta.itertuples()):
            checkpoint=out/'checkpoints'/f'{count:05d}.json'
            if checkpoint.exists():
                rows.extend(json.loads(checkpoint.read_text()))
                continue
            row_start=len(rows)
            if '$' in r.trace_name:
                bucket,selection = r.trace_name.split('$')
                ix,c,n = selection.split(',')
                x=np.asarray(h['data'][bucket][int(ix),:int(c[1:]),:int(n[1:])],dtype=np.float32)
            else:
                x=np.asarray(h['data'][r.trace_name],dtype=np.float32)
            if x.ndim != 2 or x.shape[0] != 3:
                raise ValueError('Expected three component waveform')
            if float(r.trace_sampling_rate_hz)!=20:raise ValueError('Unexpected sampling rate')
            x=sosfiltfilt(sos,x,axis=-1).astype(np.float32)
            for variant,order in [('native_Z12',[0,1,2]),('horizontal_swap_sensitivity',[0,2,1])]:
                xx=x[order]
                wave=Waveform(xx,np.isfinite(xx),0,20,('Z','N','E'),str(r.station_code))
                prob,nw=infer_probabilities(wave,model,device)
                if parity is None and nw:
                    cpu,_=infer_probabilities(wave,model,'cpu')
                    mask=np.isfinite(prob)
                    parity=dict(max_absolute_probability_difference=float(np.max(abs(cpu[mask]-prob[mask]))),identical_valid_mask=bool(np.array_equal(np.isfinite(cpu),mask)))
                    write_json(out/'cpu_gpu_parity.json',parity)
                picks=picks_from_probabilities(prob,0,20,str(r.station_code),.05)
                coverage=prediction_coverage(prob,0,20,str(r.station_code))
                for phase in ['P','S']:
                    sample=getattr(r,f'trace_{phase}_arrival_sample')
                    if not np.isfinite(sample):continue
                    label=float(sample)/20
                    cov=coverage[coverage.phase.eq(phase)]
                    covered=bool(((cov.start_s<=label)&(cov.end_s>label)).any())
                    for threshold in design['thresholds']:
                        p=picks[picks.phase.eq(phase)&(picks.model_score>=threshold)]
                        residuals=p.time_s.to_numpy()-label
                        residual=float(residuals[np.argmin(abs(residuals))]) if len(residuals) and covered else None
                        rows.append(dict(trace_name=r.trace_name,source_family_id=r.source_family_id,station=f'{r.station_network_code}.{r.station_code}',split=r.split,variant=variant,phase=phase,threshold=threshold,label_s=label,covered=covered,residual_s=residual,n_peaks=len(p)))
            write_json(checkpoint,rows[row_start:])
            if count%100==0:print('stack',count,'/',len(meta),flush=True)
    frame=pd.DataFrame(rows);frame.to_parquet(out/'arrival_results.parquet',index=False)
    summaries=[]
    for keys,g in frame.groupby(['split','variant','phase','threshold']):
        v=g[g.covered];res=v.residual_s.dropna()
        for tolerance in design['tolerances_s']:
            matched=res[abs(res)<=tolerance]
            summaries.append(dict(zip(['split','variant','phase','threshold'],keys),tolerance_s=tolerance,n_labels=len(g),n_covered=len(v),n_matched=len(matched),recall_within_tolerance=len(matched)/len(v) if len(v) else None,median_signed_residual_s=float(matched.median()) if len(matched) else None,p90_absolute_residual_s=float(abs(matched).quantile(.9)) if len(matched) else None))
    pd.DataFrame(summaries).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_stack_diagnostic_only',n_stacks=len(meta),device=device,parity=parity,design_sha256=sha256(out/'design_before_inference.json'),production_picker_validated=False))
    print('Complete: published-stack diagnostic; continuous validation still required.',flush=True)

if __name__=='__main__':main()
