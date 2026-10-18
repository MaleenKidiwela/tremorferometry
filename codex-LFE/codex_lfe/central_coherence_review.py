"""Render strongest coherence seeds with unused sites and spectral diagnostics."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import welch
from .central_cross_station import station_windows,pair_scores
from .provenance import sha256,write_json

def main():
    out=Path('results/central_strong_waveform_review_v1');out.mkdir(exist_ok=True)
    rows=[];pages=[]
    for panel in ['central_strong_picker_panel_v1','central_strong_matched_quiet_v1']:
        source=Path('results')/panel;analysis=Path('results/central_strong_cross_station_v2')/panel
        inputs=json.loads((analysis/'inputs.json').read_text());sites=json.loads((analysis/'physical_site_selection.json').read_text())['retained'];waves={}
        for station in sites:
            a=np.load(source/(station+'_processed.npz'));waves[station]=(a['waveform'].copy(),float(a['start_s']),float(a['fs']))
            f,power=welch(a['waveform'][1:,1200:73200],fs=20,nperseg=2048,axis=-1);power=power.sum(axis=0);band=(f>=1)&(f<=8);total=power[band].sum()
            rows.append(dict(panel=panel,station=station,dominant_frequency_hz=float(f[band][np.argmax(power[band])]),fraction_power_1_2Hz=float(power[(f>=1)&(f<2)].sum()/total),fraction_power_2_4Hz=float(power[(f>=2)&(f<4)].sum()/total),fraction_power_4_8Hz=float(power[(f>=4)&(f<=8)].sum()/total)))
        seeds=pd.read_csv(analysis/'r0_cc0.6_seeds.csv').sort_values('score',ascending=False).head(3)
        for number,seed in enumerate(seeds.itertuples()):
            selected=seed.stations.split('|');t=seed.window_center_diagnostic_time_s;lag={s:0. for s in sites}
            reference=selected[0]
            x,begin,fs=waves[reference];core=station_windows(x,begin,fs,np.array([t-4]),inputs['travel_s'][reference],0,8,0)
            for station in selected[1:]:
                x,begin,fs=waves[station];padded=station_windows(x,begin,fs,np.array([t-4]),inputs['travel_s'][station],0,8,2)
                score,offset,components=pair_scores(core,padded);lag[station]=float(offset[0])
            fig,axes=plt.subplots(len(sites),1,figsize=(14,12),sharex=True);snippets={}
            for row,station in enumerate(sites):
                x,begin,fs=waves[station];center=t+inputs['travel_s'][station]+lag[station];lo=round((center-begin-20)*fs);hi=lo+round(40*fs);segment=x[1:,lo:hi]
                scale=np.percentile(np.abs(segment),95) or 1.;axis=np.arange(segment.shape[-1])/fs-20
                for channel in range(2):axes[row].plot(axis,segment[channel]/scale+3*channel,lw=.6)
                axes[row].axvspan(-4,4,color='orange',alpha=.15);axes[row].set_ylabel(station,color='black' if station in selected else 'gray');snippets[station]=segment
            axes[-1].set_xlabel('Seconds relative to predicted moveout plus selected-site correlation lag; unused sites have no fitted lag')
            fig.suptitle(f'{panel}: {pd.Timestamp(t,unit="s",tz="UTC")} — CC {seed.score:.3f}\nBlack station labels: selected triangle; gray: unused physical sites. Not phase labels or a confirmed LFE.')
            fig.tight_layout(rect=(0,0,1,.95));name=f'{panel}_{number:02d}';fig.savefig(out/(name+'.png'),dpi=140);plt.close(fig)
            np.savez_compressed(out/(name+'.npz'),**snippets);write_json(out/(name+'.json'),dict(seed=seed._asdict(),source=str(source),fitted_lag_s=lag,travel_s=inputs['travel_s'],review_status='pending'))
            pages.append(f'<h2>{name}</h2><img src="{name}.png" alt="Station waveforms around coherence seed">')
    pd.DataFrame(rows).to_csv(out/'spectral_diagnostics.csv',index=False)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Central waveform review</title><style>body{font:16px system-ui;max-width:1400px;margin:2em auto}img{width:100%}</style><h1>Central coherence seed review</h1><p>Three highest-score observed seeds per panel, selected for diagnostic inspection. These are not independent validation cases. Coherence and tremor context do not establish LFE identity.</p>'+''.join(pages))
    write_json(out/'manifest.json',dict(status='review_ready_labels_pending',source_sha256=sha256(__file__),phase_labels_created=0,production_validated=False))
if __name__=='__main__':main()
