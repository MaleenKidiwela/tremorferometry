"""Acquire bounded southern training packets and test native-channel template viability."""
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
import urllib.request
import numpy as np
import pandas as pd
from obspy import read
from .provenance import sha256, write_json
from .templates import build_station_template, cut
from .waveforms import preprocess_stream


def main():
    plan=Path('results/catalog_regional_preflight_v1/southern_station_volume_plan.csv')
    refpath=Path('results/southern_reference_audit_v1/published_reference.parquet')
    stations=pd.read_csv(plan,dtype=str,keep_default_na=False)
    reference=pd.read_parquet(refpath)
    train=reference[reference.date.between('2024-01-01','2024-01-31')&~reference.date.between('2024-01-06','2024-01-10')].copy()
    counts=train.groupby('family_id').size()
    families=sorted(counts[counts>=20].index)
    train=train[train.family_id.isin(families)].sort_values('time_s')
    family='multiple_named_families'
    epochs=pd.concat([pd.read_csv('results/map_network_metadata_2009_2026/triplet_epochs.csv'),pd.read_csv('results/map_network_metadata_ncedc_retry/triplet_epochs.csv')])
    for row in stations.itertuples():
        e=epochs[(epochs.network==row.network)&(epochs.station==row.station)&(epochs.channels==row.channels)]
        good=(pd.to_datetime(e.start_utc,utc=True)<=pd.Timestamp('2024-01-01',tz='UTC'))&(pd.to_datetime(e.end_utc_exclusive,utc=True)>=pd.Timestamp('2024-02-01',tz='UTC'))
        if not good.any(): raise ValueError('Training expansion crosses an unverified channel epoch')
    out=Path('results/southern_training_multifamily_v1');out.mkdir(exist_ok=False)
    minsep=min(111.195*np.hypot(float(a.latitude)-float(b.latitude),(float(a.longitude)-float(b.longitude))*np.cos(np.deg2rad(40.))) for i,a in stations.iterrows() for j,b in stations.iterrows() if i<j)
    write_json(out/'design_before_acquisition.json',{'family_ids':families,'selection':'All named families with >=20 references on January training dates, excluding reserved January 6–10; station/template choices use training only','reference_ids':train.reference_id.tolist(),'training_dates':sorted(train.date.unique()),'station_plan_sha256':sha256(plan),'reference_sha256':sha256(refpath),'source_sha256':sha256(__file__),'packet_bounds_s':[-35,85],'packet_count_per_station':len(train),'stations':(stations.network+'.'+stations.station).tolist(),'minimum_station_separation_km':minsep,'response_cap_per_station_bytes':32*1024*1024,'test_data_requested':False,'estimated_uncompressed_float32_bytes':int(sum(stations.min_sampling_rate_hz.astype(float))*3*4*120*len(train)),'qc':'Eight-second training template; >=20 valid windows; split-stack CC >=0.25; min four physical sites for a prospective regional detector.'})
    def fetch(row):
        station=f'{row.network}.{row.station}'
        endpoint=('https://service.ncedc.org' if row.network in ('BK','NC') else 'https://service.earthscope.org')+'/fdsnws/dataselect/1/query'
        channels=row.channels.split(',');loc=channels[0].split('.')[2] or '--'
        lines=[]
        intervals=[]
        for event in train.itertuples():
            lo,hi=event.time_s-35,event.time_s+85
            if intervals and lo<=intervals[-1][1]:
                intervals[-1][1]=max(intervals[-1][1],hi)
            else:
                intervals.append([lo,hi])
        for lo,hi in intervals:
            for channel in channels:
                begin=pd.Timestamp(lo,unit='s',tz='UTC').strftime('%Y-%m-%dT%H:%M:%S.%f')
                end=pd.Timestamp(hi,unit='s',tz='UTC').strftime('%Y-%m-%dT%H:%M:%S.%f')
                lines.append(f'{row.network} {row.station} {loc} {channel.split(".")[-1]} {begin} {end}')
        body=('\n'.join(lines)+'\n').encode()
        (out/f'{station}_request.txt').write_bytes(body)
        record={'station':station,'endpoint':endpoint,'requested_windows':len(train),'status':'pending'}
        try:
            req=urllib.request.Request(endpoint,data=body,headers={'Content-Type':'text/plain'})
            with urllib.request.urlopen(req,timeout=45) as response: data=response.read(32*1024*1024+1)
            if len(data)>32*1024*1024: raise ValueError('response size cap exceeded')
            if not data: raise ValueError('empty waveform response')
            path=out/f'{station}.mseed';path.write_bytes(data)
            stream=read(io.BytesIO(data));windows={f:[] for f in families};valid_ids={f:[] for f in families}
            for event in train.itertuples():
                wave=preprocess_stream(stream,channels,station,event.time_s-15,80,fs=20,alignment='sinc')
                snippet=cut(wave,event.time_s-10,1400,alignment='sinc')
                if snippet is not None: windows[event.family_id].append(snippet);valid_ids[event.family_id].append(event.reference_id)
            qcs={}
            for f in families:
                template,qc=build_station_template(windows[f],20,min_events=20,min_split_cc=.25)
                qcs[f]=qc
                if template is not None: np.savez_compressed(out/f'{f}__{station}_template.npz',template=template,offset_from_published_s=qc['template_start_from_published_s'])
            record.update(status='downloaded_and_checked',bytes=len(data),sha256=sha256(path),valid_reference_ids=valid_ids,qc_by_family=qcs)
        except Exception as exc:
            record.update(status='failed_not_evidence_of_absence',error=f'{type(exc).__name__}: {exc}')
        write_json(out/f'{station}_result.json',record)
        print(station,record['status'],'usable families',sum(q.get('usable',False) for q in record.get('qc_by_family',{}).values()),record.get('error',''),flush=True)
        return record
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(fetch,stations.itertuples(index=False)))
    family_sites={f:[r['station'] for r in results if r.get('qc_by_family',{}).get(f,{}).get('usable')] for f in families}
    write_json(out/'report.json',{'family_sites':family_sites,'four_site_families':[f for f,s in family_sites.items() if len(s)>=4],'training_reference_count':len(train),'results':results,'interpretation':'Training-only source-bank extension; independent calibration required for this new bank.'})
    print('Four-site families',sum(len(s)>=4 for s in family_sites.values()),flush=True)

if __name__=='__main__':main()
