"""Acquire a fixed small positive/control panel for a roster-matched published family."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import urllib.request
import urllib.parse
import pandas as pd
from obspy import read
import io
from .provenance import sha256,write_json


def fetch(task):
    index,row,shift=task
    root=Path('results/ducellier_roster_packets_v1');p=root/f'{index:02d}_{row.station}_{shift:+d}.mseed';receipt=p.with_suffix('.json')
    if p.exists() and receipt.exists():
        d=json.loads(receipt.read_text())
        if d.get('sha256')==sha256(p):return d
    t=pd.Timestamp(row.reference_utc)+pd.Timedelta(seconds=shift)
    net,sta=row.station.split('.');channels=row.channels.split(',')
    host='https://service.earthscope.org' if net=='PB' else 'https://service.ncedc.org'
    params=dict(net=net,sta=sta,loc=channels[0].split('.')[2] or '--',cha=','.join(c.split('.')[-1] for c in channels),starttime=(t-pd.Timedelta(seconds=45)).strftime('%Y-%m-%dT%H:%M:%S.%f'),endtime=(t+pd.Timedelta(seconds=95)).strftime('%Y-%m-%dT%H:%M:%S.%f'))
    url=host+'/fdsnws/dataselect/1/query?'+urllib.parse.urlencode(params)
    record=dict(station=row.station,reference_utc=row.reference_utc,shift_s=shift,channels=channels,url=url)
    try:
        with urllib.request.urlopen(url,timeout=30) as r:b=r.read(2_000_001)
        if len(b)>2_000_000 or not b:raise ValueError('Empty or oversized response')
        stream=read(io.BytesIO(b));p.write_bytes(b)
        record.update(status='received',bytes=len(b),sha256=sha256(p),path=str(p),received_channels=sorted({tr.id for tr in stream}))
    except Exception as exc:record.update(status='failed_not_evidence_of_absence',error=f'{type(exc).__name__}: {exc}')
    write_json(receipt,record);print(row.station,row.reference_utc,shift,record['status'],flush=True)
    return record


def main():
    root=Path('results/ducellier_roster_packets_v1');root.mkdir(exist_ok=True)
    src=Path('results/southern_benchmark_inventory_v1/targeted_reference_packet_plan.csv')
    plan=pd.read_csv(src)
    tasks=[(i,r,shift) for i,r in enumerate(plan.itertuples()) for shift in [0,600]]
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(fetch,tasks))
    write_json(root/'manifest.json',dict(status='bounded_packet_acquisition_complete',plan_sha256=sha256(src),source_sha256=sha256(__file__),records=records,total_bytes=sum(r.get('bytes',0) for r in records),event_identity_validated=False))

if __name__=='__main__':main()
