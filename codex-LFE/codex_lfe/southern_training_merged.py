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
    stations=stations[stations.network.isin(['BK','NC'])]
    reference=pd.read_parquet(refpath)
    train=reference[reference.date.between('2024-01-01','2024-01-05')]
    family=train.groupby('family_id').size().sort_values(ascending=False).index[0]
    train=train[train.family_id==family].sort_values('time_s')
    out=Path('results/southern_training_ncedc_merged_v1');out.mkdir(exist_ok=False)
    minsep=min(111.195*np.hypot(float(a.latitude)-float(b.latitude),(float(a.longitude)-float(b.longitude))*np.cos(np.deg2rad(40.))) for i,a in stations.iterrows() for j,b in stations.iterrows() if i<j)
    write_json(out/'design_before_acquisition.json',{'family_id':family,'selection':'Most frequent named family on the fixed five training dates only','reference_ids':train.reference_id.tolist(),'training_dates':sorted(train.date.unique()),'station_plan_sha256':sha256(plan),'reference_sha256':sha256(refpath),'source_sha256':sha256(__file__),'packet_bounds_s':[-35,85],'packet_count_per_station':len(train),'stations':(stations.network+'.'+stations.station).tolist(),'minimum_station_separation_km':minsep,'response_cap_per_station_bytes':32*1024*1024,'test_data_requested':False,'estimated_uncompressed_float32_bytes':int(sum(stations.min_sampling_rate_hz.astype(float))*3*4*120*len(train)),'qc':'Eight-second training template; >=20 valid windows; split-stack CC >=0.25; min four physical sites for a prospective regional detector.'})
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
            stream=read(io.BytesIO(data));windows=[];valid_ids=[]
            for event in train.itertuples():
                wave=preprocess_stream(stream,channels,station,event.time_s-15,80,fs=20,alignment='sinc')
                snippet=cut(wave,event.time_s-10,1400,alignment='sinc')
                if snippet is not None: windows.append(snippet);valid_ids.append(event.reference_id)
            template,qc=build_station_template(windows,20,min_events=20,min_split_cc=.25)
            record.update(status='downloaded_and_checked',bytes=len(data),sha256=sha256(path),valid_reference_ids=valid_ids,qc=qc)
            if template is not None: np.savez_compressed(out/f'{station}_template.npz',template=template,offset_from_published_s=qc['template_start_from_published_s'])
        except Exception as exc:
            record.update(status='failed_not_evidence_of_absence',error=f'{type(exc).__name__}: {exc}')
        write_json(out/f'{station}_result.json',record)
        print(station,record['status'],record.get('qc',record.get('error')),flush=True)
        return record
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(fetch,stations.itertuples(index=False)))
    usable=[r['station'] for r in results if r.get('qc',{}).get('usable')]
    write_json(out/'report.json',{'family_id':family,'training_reference_count':len(train),'usable_template_stations':usable,'four_site_template_viability':len(usable)>=4,'results':results,'interpretation':'Training consistency only; not a southern candidate catalog or independent validation.'})

if __name__=='__main__':main()
