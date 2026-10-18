"""Evaluate only stack S labels excluded by the original fixed-grid edge mask."""
from pathlib import Path
import os
import json
import h5py
import numpy as np
import pandas as pd
from scipy.signal import butter,sosfiltfilt
from .neural import infer_probabilities,picks_from_probabilities,prediction_coverage
from .waveforms import Waveform
from .provenance import sha256,write_json


def main():
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch
    import seisbench.models as sbm
    torch.set_num_threads(2)
    src=Path('../lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia')
    old=Path('results/picker_labeled_stacks_jms_v2')
    out=Path('results/picker_labeled_stacks_tail_v1');out.mkdir(exist_ok=True)
    (out/'checkpoints').mkdir(exist_ok=True)
    source=pd.read_parquet(old/'arrival_results.parquet')
    names=set(source.loc[~source.covered,'trace_name'])
    meta=pd.read_csv(src/'metadata.csv')
    weights=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    design=dict(rule='Only previously uncovered S labels; evaluate a final unpadded 60-second window ending at the actual trace end. Window placement uses record length, never arrival label position. Existing covered P/S rows are reused unchanged.',
        input_sha256=sha256(old/'arrival_results.parquet'),metadata_sha256=sha256(src/'metadata.csv'),waveforms_sha256=sha256(src/'waveforms.hdf5'),weights_sha256=sha256(weights.with_suffix('.pt')),source_sha256=sha256(__file__),variants=['native_Z12','horizontal_swap_sensitivity'],production_code_changed=False)
    dp=out/'design_before_inference.json'
    if dp.exists() and json.loads(dp.read_text())!=design:raise ValueError('Changed design')
    write_json(dp,design)
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True)
    device='cuda' if torch.cuda.is_available() else 'cpu'
    sos=butter(4,[1,8],fs=20,btype='bandpass',output='sos')
    additions=[]
    with h5py.File(src/'waveforms.hdf5') as h:
        for idx,r in enumerate(meta[meta.trace_name.isin(names)].itertuples()):
            cp=out/'checkpoints'/f'{idx:05d}.json'
            if cp.exists():additions.extend(json.loads(cp.read_text()));continue
            begin=len(additions)
            if '$' in r.trace_name:
                b,k=r.trace_name.split('$');i,c,n=k.split(',')
                x=np.asarray(h['data'][b][int(i),:int(c[1:]),:int(n[1:])],dtype=np.float32)
            else:x=np.asarray(h['data'][r.trace_name],dtype=np.float32)
            x=sosfiltfilt(sos,x,axis=-1).astype(np.float32)
            offset=max(0,x.shape[1]-1200)
            for variant,order in [('native_Z12',[0,1,2]),('horizontal_swap_sensitivity',[0,2,1])]:
                xx=x[order,offset:]
                w=Waveform(xx,np.isfinite(xx),offset/20,20,('Z','N','E'),r.station_code)
                prob,_=infer_probabilities(w,model,device)
                cov=prediction_coverage(prob,w.start_s,20,r.station_code)
                picks=picks_from_probabilities(prob,w.start_s,20,r.station_code,.05)
                subset=source[(source.trace_name==r.trace_name)&(source.variant==variant)&~source.covered]
                for a in subset.itertuples():
                    cv=cov[cov.phase==a.phase]
                    covered=bool(((cv.start_s<=a.label_s)&(cv.end_s>a.label_s)).any())
                    pp=picks[(picks.phase==a.phase)&(picks.model_score>=a.threshold)]
                    residuals=pp.time_s.to_numpy()-a.label_s
                    row={k:getattr(a,k) for k in source.columns}
                    row.update(covered=covered,residual_s=float(residuals[np.argmin(abs(residuals))]) if len(residuals) and covered else None,n_peaks=len(pp))
                    additions.append(row)
            write_json(cp,additions[begin:])
            if idx%200==0:print('tail check',idx,'/',len(names),flush=True)
    fixed=pd.concat([source[source.covered],pd.DataFrame(additions)],ignore_index=True)
    assert len(fixed)==len(source)
    assert not fixed.duplicated(['trace_name','variant','phase','threshold']).any()
    fixed.to_parquet(out/'arrival_results.parquet',index=False)
    rows=[]
    for key,g in fixed.groupby(['variant','phase','threshold']):
        v=g[g.covered]
        for tol in [.25,.5,1.0]:
            m=v[abs(v.residual_s)<=tol]
            rows.append(dict(zip(['variant','phase','threshold'],key),tolerance_s=tol,n_labels=len(g),n_covered=len(v),n_matched=len(m),recall_within_tolerance=len(m)/len(v) if len(v) else None,median_signed_residual_s=float(m.residual_s.median()) if len(m) else None))
    pd.DataFrame(rows).to_csv(out/'summary.csv',index=False)
    write_json(out/'report.json',dict(status='complete_stack_coverage_diagnostic',original_rows=len(source),reused_covered_rows=int(source.covered.sum()),tail_evaluated_rows=len(additions),remaining_uncovered_rows=int((~fixed.covered).sum()),production_picker_validated=False))
    print('Completed tail coverage check.',flush=True)

if __name__=='__main__':main()
