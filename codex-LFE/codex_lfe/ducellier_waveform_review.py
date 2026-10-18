"""Inspectable waveform/template overlays for every next-day candidate row.

This is a template audit, not a blinded phase-label set. All review labels remain
pending. Existing raw recordings are reused; no network requests or neural runs.
"""
from pathlib import Path
import json,html
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .published_pickle import read_published_pickle
from .waveforms import preprocess_stream
from .neural import orientation_matrix
from .provenance import write_json,sha256

def rms_scaled(a):
    a=a-a.mean(axis=1,keepdims=True)
    return a/np.maximum(np.sqrt(np.mean(a*a,axis=1,keepdims=True)),1e-30)

def main():
    root=Path('results/ducellier_multifamily_november23_holdout_edges_v2')
    raw=Path('results/ducellier_november23_acquisition_v1');out=Path('results/ducellier_candidate_waveform_review_v1');out.mkdir(exist_ok=True)
    cases=pd.read_csv(root/'within_family_review_candidates.csv');cases=cases[cases.control==-1].reset_index(drop=True)
    cases['review_id']=[f'case_{i:03d}' for i in range(len(cases))];cases['review_status']='pending';cases.to_csv(out/'case_manifest.csv',index=False)
    source=Path('results/southern_benchmark_inventory_v1/source/src/templates_v1')
    base=Path('results/ducellier_roster_packets_v1/responses')
    inv=read_inventory(str(base/'BK.xml'))+read_inventory(str(base/'NC.xml'))+read_inventory(str(base/'PB.xml'))
    acquisition=json.loads((raw/'acquisition.json').read_text());packets={r.review_id:{} for r in cases.itertuples()}
    for record in acquisition['records']:
        if record['status']!='received':continue
        station=record['station'];channels=record['channels'];stream=read(record['path'])
        if sha256(record['path'])!=record['sha256']:raise ValueError('Raw checksum')
        for case in cases.itertuples():
            paths=[source/case.source_template_id/f'{station.split(".")[1]}_{c.split(".")[-1]}.pkl' for c in channels]
            if not all(p.exists() for p in paths):continue
            traces=[read_published_pickle(p) for p in paths]
            template=np.asarray([tr.data for tr,angle in traces]);n=template.shape[1];t=UTCDateTime(case.time_s)
            segment=stream.slice(t-60,t+n/20+60).copy()
            for tr in segment:
                tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80)
            wave=preprocess_stream(segment,channels,station,float(t)-5,n/20+10,fs=20,band=(1.5,9))
            current=orientation_matrix([inv.get_orientation(c,t) for c in channels]);target=orientation_matrix([angle for tr,angle in traces])
            x=target@np.linalg.solve(current,wave.x)
            valid=wave.valid.all(axis=0)
            x[:,~valid]=np.nan
            packets[case.review_id][station]=(x,template,channels)
        print(station,'review snippets prepared',flush=True)
    links=[];outputs={}
    for case in cases.itertuples():
        panel=packets[case.review_id];stations=sorted(panel)
        fig,axes=plt.subplots(len(stations),3,figsize=(14,2*len(stations)),squeeze=False,layout='constrained')
        archive={}
        for i,station in enumerate(stations):
            data,template,channels=panel[station];archive[station+'_data']=data;archive[station+'_template']=template
            data=rms_scaled(data);template=rms_scaled(template)
            for c in range(3):
                ax=axes[i,c];ax.plot(np.arange(data.shape[1])/20-5,data[c],color='.35',lw=.6,label='Recorded waveform')
                ax.plot(np.arange(template.shape[1])/20,template[c],color='#c33732',lw=.7,alpha=.85,label='Published template')
                label=channels[c]+(' (WITHHELD)' if station=='PB.B039' else '')
                ax.set_title(label,fontsize=9);ax.set_xlim(-5,template.shape[1]/20+5);ax.tick_params(labelsize=8)
                if i==len(stations)-1:ax.set_xlabel('Seconds relative to source-template start')
        axes[0,0].legend(fontsize=7)
        fig.suptitle(f'{case.review_id} | {UTCDateTime(case.time_s)} | template {case.source_template_id}\nComponent RMS scaled independently; pending template audit, not phase truth',fontsize=12)
        png=out/(case.review_id+'.png');fig.savefig(png,dpi=120);plt.close(fig)
        np.savez_compressed(out/(case.review_id+'.npz'),**archive)
        outputs[png.name]=sha256(png);links.append(f'<li><a href="{case.review_id}.png">{case.review_id}</a> — {html.escape(str(UTCDateTime(case.time_s)))} — pending</li>')
        print(case.review_id,'saved',flush=True)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Candidate waveform audit</title><h1>Next-day candidate waveform audit</h1><p>All reviews pending. Overlays use published templates and are not blinded phase labels. PB.B039 was withheld from detection. Candidate counts are not finalized unique events.</p><ul>'+''.join(links)+'</ul>')
    write_json(out/'verification.json',dict(status='packet_complete_labels_pending',candidate_rows=len(cases),labels_submitted=0,new_waveform_downloads=0,source_sha256=sha256(__file__),candidate_source_sha256=sha256(root/'within_family_review_candidates.csv'),image_sha256=outputs))

if __name__=='__main__':main()
