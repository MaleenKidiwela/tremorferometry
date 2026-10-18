"""Complete multi-family training acquisition with individual disjoint-window requests."""
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
import urllib.request
import urllib.parse
import numpy as np
import pandas as pd
from obspy import read, Stream
from .provenance import sha256, write_json
from .waveforms import preprocess_stream
from .templates import build_station_template, cut


def main():
    original=Path('results/southern_training_multifamily_v1')
    design=json.loads((original/'design_before_acquisition.json').read_text())
    reference=pd.read_parquet('results/southern_reference_audit_v1/published_reference.parquet')
    reference=reference[reference.reference_id.isin(design['reference_ids'])].sort_values('time_s')
    families=design['family_ids']
    stations=pd.read_csv('results/catalog_regional_preflight_v1/southern_station_volume_plan.csv',dtype=str,keep_default_na=False)
    intervals=[]
    for event in reference.itertuples():
        lo,hi=event.time_s-35,event.time_s+85
        if intervals and lo<=intervals[-1][1]:intervals[-1][1]=max(intervals[-1][1],hi)
        else:intervals.append([lo,hi])
    out=Path('results/southern_training_complete_v1');out.mkdir(exist_ok=True)
    write_json(out/'design_before_acquisition.json',{'source_design':str(original/'design_before_acquisition.json'),'source_design_sha256':sha256(original/'design_before_acquisition.json'),'family_ids':families,'reference_ids':reference.reference_id.tolist(),'union_intervals_s':intervals,'reason':'Multi-window NCEDC POST response does not cover all requested packets. Request each disjoint union separately and verify masks. Reuse EarthScope packet files where available.','source_sha256':sha256(__file__)})
    def process(row):
        station=f'{row.network}.{row.station}';channels=row.channels.split(',');records=[]
        if row.network=='PB':
            path=original/f'{station}.mseed'
            while not (original/f'{station}_result.json').exists():
                import time;time.sleep(1)
            stream=read(str(path));records.append({'reused':str(path),'sha256':sha256(path)})
        else:
            stream=Stream();request_dir=out/station;request_dir.mkdir(exist_ok=True)
            for i,(lo,hi) in enumerate(intervals):
                path=request_dir/f'{i:04d}.mseed';receipt=path.with_suffix('.json')
                if path.exists() and receipt.exists():
                    record=json.loads(receipt.read_text())
                    if record['sha256']!=sha256(path):raise ValueError('Changed saved packet')
                    stream+=read(str(path));records.append(record);continue
                params={'net':row.network,'sta':row.station,'loc':channels[0].split('.')[2] or '--','cha':','.join(c.split('.')[-1] for c in channels),
                        'starttime':pd.Timestamp(lo,unit='s',tz='UTC').strftime('%Y-%m-%dT%H:%M:%S.%f'),
                        'endtime':pd.Timestamp(hi,unit='s',tz='UTC').strftime('%Y-%m-%dT%H:%M:%S.%f')}
                url='https://service.ncedc.org/fdsnws/dataselect/1/query?'+urllib.parse.urlencode(params)
                record={'url':url,'interval_s':[lo,hi]}
                try:
                    with urllib.request.urlopen(url,timeout=40) as response:data=response.read(16*1024*1024+1)
                    if len(data)>16*1024*1024:raise ValueError('packet cap exceeded')
                    piece=read(io.BytesIO(data));path.write_bytes(data);stream+=piece
                    record.update(status='received',bytes=len(data),sha256=sha256(path));write_json(receipt,record)
                except Exception as exc:
                    record.update(status='failed',error=f'{type(exc).__name__}: {exc}');write_json(request_dir/f'{i:04d}_failure.json',record)
                records.append(record)
                if (i+1)%20==0:print(station,'downloaded',i+1,'/',len(intervals),flush=True)
        windows={f:[] for f in families};ids={f:[] for f in families}
        for event in reference.itertuples():
            wave=preprocess_stream(stream,channels,station,event.time_s-15,80,fs=20,alignment='sinc')
            packet=cut(wave,event.time_s-10,1400,alignment='sinc')
            if packet is not None:windows[event.family_id].append(packet);ids[event.family_id].append(event.reference_id)
        qcs={}
        for family in families:
            template,qc=build_station_template(windows[family],20,min_events=20,min_split_cc=.25);qcs[family]=qc
            if template is not None:np.savez_compressed(out/f'{family}__{station}_template.npz',template=template,offset_from_published_s=qc['template_start_from_published_s'])
        record={'station':station,'qc_by_family':qcs,'valid_reference_ids':ids,'waveform_requests':records}
        write_json(out/f'{station}_result.json',record);print(station,'usable families',sum(q.get('usable',False) for q in qcs.values()),flush=True)
        return record
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(process,stations.itertuples(index=False)))
    sites={f:[r['station'] for r in results if r['qc_by_family'][f].get('usable')] for f in families}
    write_json(out/'report.json',{'family_sites':sites,'four_site_families':[f for f,s in sites.items() if len(s)>=4],'training_reference_count':len(reference),'results':results,'interpretation':'Training-only templates. Calibration required before using this expanded source bank.'})
    print('Four-site families',sum(len(s)>=4 for s in sites.values()),flush=True)

if __name__=='__main__':main()
