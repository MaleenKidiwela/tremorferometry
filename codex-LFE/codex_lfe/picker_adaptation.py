"""Bounded northern P/S adaptation with spatially grouped source holdouts.

Synthetic validation is developmental evidence only. Real waveform validation is
required separately; no production picker is replaced by this experiment.
"""
from pathlib import Path
import copy,hashlib,json,os,re,time
import h5py
import numpy as np
import pandas as pd
from scipy.signal import butter,sosfiltfilt,find_peaks
from scipy.sparse.csgraph import connected_components
from .provenance import write_json,sha256

ROOT=Path(__file__).resolve().parents[1]
FS=20.;N=1200


def source_splits(metadata):
    sources=metadata.groupby('source_family_id')[['source_latitude_deg','source_longitude_deg','source_depth_km']].first()
    x=sources.to_numpy();scale=np.array([111.195,111.195*np.cos(np.deg2rad(x[:,0].mean())),1.])
    distance=np.linalg.norm((x[:,None]-x[None,:])*scale,axis=-1)
    _,groups=connected_components(distance<=2.,directed=False)
    split={};rows=[]
    for group in sorted(set(groups)):
        ids=sources.index[groups==group].tolist()
        value=int(hashlib.sha256(('20260914:'+','.join(str(int(i)) for i in ids)).encode()).hexdigest()[:8],16)/2**32
        fold='train' if value<.6 else ('validation' if value<.8 else 'test')
        for source in ids:split[source]=fold;rows.append(dict(source_family_id=source,spatial_group=group,split=fold))
    return split,pd.DataFrame(rows)


def prepare_sources(metadata,path,splits):
    sources={k:[] for k in ['train','validation','test']};sos=butter(4,[1,8],fs=FS,btype='bandpass',output='sos')
    with h5py.File(path) as h5:
        for row in metadata.itertuples():
            if row.trace_sampling_rate_hz!=FS:raise ValueError('Unexpected stack sample rate')
            labels=np.array([row.trace_P_arrival_sample,row.trace_S_arrival_sample],float)
            known=np.isfinite(labels)
            if not known.any():continue
            key=re.fullmatch(r'(\w+)\$(\d+),:3,:(\d+)',row.trace_name)
            if key:
                x=h5['data'][key[1]][int(key[2]),:3,:int(key[3])]
            elif row.trace_name in h5['data']:
                x=h5['data'][row.trace_name][:]
                if x.ndim!=2 or x.shape[0]!=3:raise ValueError('Unsupported unbucketed trace shape')
            else:raise ValueError('Unsupported stack address')
            x=sosfiltfilt(sos,x,axis=1)
            lo=max(0,int(labels[known].min()-4*FS));hi=min(x.shape[1],int(labels[known].max()+10*FS))
            x=x[:,lo:hi].copy();taper=np.ones(x.shape[1]);ramp=np.sin(np.linspace(0,np.pi/2,int(FS)))**2;taper[:len(ramp)]=ramp;taper[-len(ramp):]=ramp[::-1];x*=taper
            rms=np.sqrt(np.mean(x*x))
            if not np.isfinite(rms) or rms<=0:continue
            sources[splits[row.source_family_id]].append(dict(waveform=(x/rms).astype(np.float32),labels=labels-lo,known=known,family=int(row.source_family_id)))
    if any(not v for v in sources.values()):raise ValueError('Empty source split')
    return sources


def prepare_noise(panel):
    result={k:[] for k in ['train','validation','test']}
    files=sorted(panel.glob('*_processed.npz'))
    for path in files:
        with np.load(path) as a:
            # Core starts after the 60 s processing pad. Nonoverlapping 60 s
            # waveform intervals: first/second/third twenty minutes are disjoint.
            for i in range(60):
                first=1200+i*N;x=a['waveform'][:,first:first+N]
                if x.shape!=(3,N) or not a['valid'][:,first:first+N].all():continue
                fold=['train','validation','test'][i//20]
                result[fold].append(np.fft.rfft(x-x.mean(axis=1,keepdims=True),axis=1))
    if any(not v for v in result.values()):raise ValueError('Empty noise split')
    return {k:np.asarray(v) for k,v in result.items()},files


def batch(sources,noise,rng,size):
    xx=[];yy=[];masks=[];truth=[]
    grid=np.arange(N)
    for _ in range(size):
        spec=noise[rng.integers(len(noise))]
        # A shared random phase preserves three-component cross spectra while
        # removing the original temporal earthquake signatures. Synthetic noise.
        phase=np.exp(1j*rng.uniform(-np.pi,np.pi,spec.shape[-1]));phase[0]=1;phase[-1]=1
        x=np.fft.irfft(spec*phase,n=N,axis=1);rms=np.sqrt(np.mean(x*x))
        if rms<=0:raise ValueError('Zero noise')
        x/=rms;y=np.zeros((2,N));mask=np.ones((2,N),bool);times=[[],[]];centers=[]
        n_events=rng.choice([0,1,2],p=[.3,.5,.2])
        for _event in range(n_events):
            s=sources[rng.integers(len(sources))];chosen=rng.choice(np.flatnonzero(s['known']))
            center=int(rng.integers(200,1001))
            if any(abs(center-c)<200 for c in centers):continue
            start=int(round(center-s['labels'][chosen]));labels=s['labels']+start
            if any(any(abs(t-v)<FS for v in times[p]) for p,t in enumerate(labels) if np.isfinite(t)):continue
            signal=s['waveform'].copy();angle=rng.uniform(-np.pi,np.pi)
            signal[1:]=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])@signal[1:]
            signal*=rng.choice([-1,1])*np.exp(rng.uniform(np.log(.25),np.log(4.)))
            lo,hi=max(0,start),min(N,start+signal.shape[1])
            if hi<=lo:continue
            x[:,lo:hi]+=signal[:,lo-start:hi-start];centers.append(center)
            for p,t in enumerate(labels):
                if not s['known'][p]:mask[p,:]=False
                elif 0<=t<N:
                    y[p]=np.maximum(y[p],np.exp(-.5*((grid-t)/(FS*.2))**2));times[p].append(float(t/FS))
        x-=x.mean(axis=1,keepdims=True);x/=max(np.max(np.abs(x)),1e-10)
        xx.append(x);yy.append(y);masks.append(mask);truth.append(times)
    return np.asarray(xx,np.float32),np.asarray(yy,np.float32),np.asarray(masks),truth


def cardinality(predicted,reference,tolerance=.5):
    a=sorted(predicted);b=sorted(reference);i=j=count=0
    while i<len(a) and j<len(b):
        if a[i]<b[j]-tolerance:i+=1
        elif a[i]>b[j]+tolerance:j+=1
        else:count+=1;i+=1;j+=1
    return count


def evaluate(model,data,device):
    import torch
    x,y,mask,truth=data;model.eval();prob=[]
    with torch.no_grad():
        for i in range(0,len(x),64):prob.append(model(torch.from_numpy(x[i:i+64]).to(device)).cpu().numpy())
    prob=np.concatenate(prob);rows=[]
    for p,phase in enumerate(['P','S']):
        tp=fp=fn=0;windows=0;noise_fp=0;noise_windows=0
        for k in range(len(x)):
            if not mask[k,p].all():continue
            windows+=1;peaks=find_peaks(prob[k,p],height=.3,distance=int(FS))[0]/FS
            matched=cardinality(peaks,truth[k][p]);tp+=matched;fp+=len(peaks)-matched;fn+=len(truth[k][p])-matched
            if not truth[k][p]:noise_windows+=1;noise_fp+=len(peaks)
        rows.append(dict(phase=phase,tp=tp,fp=fp,fn=fn,precision=tp/max(1,tp+fp),recall=tp/max(1,tp+fn),f1=2*tp/max(1,2*tp+fp+fn),scored_windows=windows,synthetic_noise_false_peaks=noise_fp,synthetic_noise_windows=noise_windows))
    return rows,prob


def main():
    import torch
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(ROOT/'cache/seisbench'))
    import seisbench.models as sbm
    out=ROOT/'results/north_picker_adaptation_v1';out.mkdir(exist_ok=True)
    if (out/'training_progress.json').exists() or (out/'report.json').exists():raise ValueError('Existing training run is immutable')
    external=ROOT.parent/'lfe_catalog_repair/results/external/other'
    metadata_path=external/'seisbench_lfestacks_cascadia/metadata.csv';waveforms_path=metadata_path.with_name('waveforms.hdf5')
    weights=external/'seisbench_lfedetect/jms';metadata=pd.read_csv(metadata_path)
    splits,split_rows=source_splits(metadata);split_rows.to_csv(out/'source_splits.csv',index=False)
    sources=prepare_sources(metadata,waveforms_path,splits)
    noise,noise_paths=prepare_noise(ROOT/'results/north_recent_quiet_picker_v1')
    design=dict(status='running_development_experiment',epochs=8,steps_per_epoch=32,batch_size=64,learning_rate=.0001,weight_decay=.0001,seed=20260914,
        source_split='2 km 3D connected source groups assigned by stable hash: 60/20/20 train/validation/test; no source group crosses splits',
        source_family_counts=split_rows.groupby('split').size().to_dict(),trace_counts={k:len(v) for k,v in sources.items()},
        noise='Component-coherent phase-randomized spectra from disjoint twenty-minute portions of the eight-station August 20 comparison hour; synthetic noise, not labeled natural negatives.',
        signal='1-8 Hz published P/S packets; random horizontal rotation, polarity, placement and packet/noise RMS ratio 0.25-4; 0/1/2 injections',
        unknown_phase_policy='Mask a whole output phase in any window containing a source without that phase label',
        selection='Best validation macro P/S F1 at fixed .3 threshold and .5 s tolerance; score test split once after selection',
        production_validated=False,limitations=['Synthetic accuracy cannot certify natural-waveform precision or recall.','Natural waveform evaluation is required before deployment.','The original jms model remains unchanged.'],
        sources={str(p):sha256(p) for p in [Path(__file__),metadata_path,waveforms_path,weights.with_suffix('.pt'),weights.with_suffix('.json'),*noise_paths]})
    write_json(out/'design_before_training.json',design)
    torch.manual_seed(20260914);torch.set_num_threads(2);device='cuda' if torch.cuda.is_available() else 'cpu'
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True).to(device)
    validation=batch(sources['validation'],noise['validation'],np.random.default_rng(100),512)
    baseline,_=evaluate(model,validation,device)
    opt=torch.optim.AdamW(model.parameters(),lr=.0001,weight_decay=.0001);rng=np.random.default_rng(20260914)
    history=[];best=-1;best_state=None;best_epoch=None
    for epoch in range(1,9):
        model.train();total=0
        for _ in range(32):
            x,y,mask,_=batch(sources['train'],noise['train'],rng,64)
            x=torch.from_numpy(x).to(device);y=torch.from_numpy(y).to(device);mask=torch.from_numpy(mask).to(device)
            opt.zero_grad(set_to_none=True);logits=model(x,logits=True)
            values=torch.nn.functional.binary_cross_entropy_with_logits(logits,y,reduction='none')
            loss=values[mask].mean();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step();total+=float(loss.detach())
        scores,_=evaluate(model,validation,device);metric=float(np.mean([r['f1'] for r in scores]));history.append(dict(epoch=epoch,loss=total/32,validation=scores))
        if metric>best:best=metric;best_state={k:v.detach().cpu().clone() for k,v in model.state_dict().items()};best_epoch=epoch
        write_json(out/'training_progress.json',dict(epoch=epoch,total_epochs=8,history=history));print('epoch',epoch,'loss',total/32,'validation macro F1',metric,flush=True)
    model.load_state_dict(best_state);torch.save(best_state,out/'adapted_state_dict.pt')
    test=batch(sources['test'],noise['test'],np.random.default_rng(101),1024)
    adapted,prob=evaluate(model,test,device);np.savez_compressed(out/'test_predictions.npz',probabilities=prob,targets=test[1],mask=test[2]);write_json(out/'test_truth_times.json',test[3])
    original=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True).to(device);baseline_test,bprob=evaluate(original,test,device);np.savez_compressed(out/'baseline_test_predictions.npz',probabilities=bprob)
    rows=[dict(model=name,**r) for name,data in [('baseline',baseline_test),('adapted',adapted)] for r in data];pd.DataFrame(rows).to_csv(out/'test_metrics.csv',index=False)
    report=dict(status='complete_synthetic_adaptation_test',best_epoch=best_epoch,baseline_validation=baseline,history=history,test=rows,
        production_validated=False,natural_waveform_validation_completed=False,weights_sha256=sha256(out/'adapted_state_dict.pt'))
    write_json(out/'report.json',report);print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':main()
