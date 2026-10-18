"""Resumable southern calibration, frozen evaluation and candidate export."""
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import io
import json
import multiprocessing
from pathlib import Path
import urllib.parse
import urllib.request
import time
import numpy as np
import pandas as pd
from obspy import read
from .pilot import ROOT, scan_date_worker, summarize_split
from .provenance import sha256, source_manifest, write_json
from .roster import validate_membership
from .site_weighting import operating_point


def fetch_day(task):
    station,meta,date,root=task
    n,s=station.split('.')
    day=pd.Timestamp(date);path=root/station/str(day.year)/f'{day.dayofyear:03d}.mseed'
    receipt=path.with_suffix('.json')
    if path.exists() and receipt.exists():
        record=json.loads(receipt.read_text())
        if record.get('sha256')==sha256(path):return record
    path.parent.mkdir(parents=True,exist_ok=True)
    channels=meta['channels'];location=channels[0].split('.')[2] or '--'
    endpoint=('https://service.ncedc.org' if n in ('BK','NC') else 'https://service.earthscope.org')+'/fdsnws/dataselect/1/query'
    params={'net':n,'sta':s,'loc':location,'cha':','.join(c.split('.')[-1] for c in channels),
            'starttime':day.strftime('%Y-%m-%dT%H:%M:%S'),'endtime':(day+pd.Timedelta(days=1)).strftime('%Y-%m-%dT%H:%M:%S')}
    url=endpoint+'?'+urllib.parse.urlencode(params)
    record={'station':station,'date':date,'url':url}
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=45) as response:
                data=response.read(128*1024*1024+1)
                record['headers']=dict(response.headers)
            if not data or len(data)>128*1024*1024:raise ValueError('empty or oversized day response')
            stream=read(io.BytesIO(data));coverage={}
            for channel in channels:
                segments=sorted((max(float(tr.stats.starttime),day.timestamp()),min(float(tr.stats.endtime)+tr.stats.delta,(day+pd.Timedelta(days=1)).timestamp())) for tr in stream.select(id=channel))
                total=0;last=-np.inf
                for lo,hi in segments:total+=max(0,hi-max(lo,last));last=max(last,hi)
                coverage[channel]=total/86400
            tmp=path.with_suffix('.part');tmp.write_bytes(data);tmp.replace(path)
            record.update(status='received',bytes=len(data),sha256=sha256(path),raw_time_coverage_fraction=coverage,attempt=attempt+1)
            write_json(receipt,record)
            print('download',station,date,len(data),'bytes; minimum time coverage',round(min(coverage.values()),4),flush=True)
            return record
        except Exception as exc:
            record.update(status='failed',error=f'{type(exc).__name__}: {exc}',attempt=attempt+1)
            print('download retry',station,date,record['error'],flush=True)
    write_json(receipt,record)
    raise RuntimeError(f'Failed acquisition: {station} {date}: {record["error"]}')


def main():
    training_root=Path('results/southern_training_expanded_v1')
    training=json.loads((training_root/'report.json').read_text())
    family=training['family_id'];anchor='BK.RVIT'
    plan=pd.read_csv('results/catalog_regional_preflight_v1/southern_station_volume_plan.csv',dtype=str,keep_default_na=False)
    cfg=json.loads(Path('results/method_interleaved_mf_v1/manifest.json').read_text())['config']
    cfg.update(name='southern_one_family_2024_v1',anchor_station=anchor,family_count=1,
        train_dates=json.loads((training_root/'design_before_acquisition.json').read_text())['training_dates'],
        calibration_dates=pd.date_range('2024-02-01','2024-02-14').strftime('%Y-%m-%d').tolist(),
        test_dates=pd.date_range('2024-02-15','2024-02-29').strftime('%Y-%m-%d').tolist(),
        reference_csv=str((ROOT/'results/southern_reference_audit_v1/published_reference.parquet').resolve()),
        raw_root=str((ROOT/'traces/southern_2024').resolve()),stations={})
    for row in plan.itertuples():
        station=f'{row.network}.{row.station}'
        if station in training['usable_template_stations']:
            cfg['stations'][station]={'site_id':station,'channels':row.channels.split(',')}
    if len(cfg['stations'])<4:raise ValueError('Four training-qualified sites required')
    epochs=pd.concat([pd.read_csv('results/map_network_metadata_2009_2026/triplet_epochs.csv'),pd.read_csv('results/map_network_metadata_ncedc_retry/triplet_epochs.csv')])
    for station,meta in cfg['stations'].items():
        n,s=station.split('.');e=epochs[(epochs.network==n)&(epochs.station==s)&(epochs.channels==','.join(meta['channels']))]
        good=(pd.to_datetime(e.start_utc,utc=True)<=pd.Timestamp('2024-02-01',tz='UTC'))&(pd.to_datetime(e.end_utc_exclusive,utc=True)>=pd.Timestamp('2024-03-01',tz='UTC'))
        if not good.any():raise ValueError('Unverified February channel epoch')
    anchor_peak=next(r['qc']['peak_offset_from_published_s'] for r in training['results'] if r['station']==anchor)
    templates={family:{}};archive={}
    for station in cfg['stations']:
        with np.load(training_root/f'{station}_template.npz') as z:
            t=z['template'].copy();offset=float(z['offset_from_published_s'])-anchor_peak
        templates[family][station]=(t,offset);archive[f'{family}__{station}__template']=t;archive[f'{family}__{station}__offset_s']=offset
    out=Path('results/southern_catalog_2024_v1');out.mkdir(exist_ok=True)
    reference=pd.read_parquet(cfg['reference_csv'])
    design={'config':cfg,'family_id':family,'anchor_from_published_s':anchor_peak,
        'roster':validate_membership(cfg['stations'],ROOT/cfg['station_roster']),
        'date_rule':'First half of February calibration; second half sealed evaluation. January-only training. Calendar choice before candidate processing.',
        'calibration_published_reference_count':int(((reference.family_id==family)&reference.date.isin(cfg['calibration_dates'])).sum()),
        'calibration_float32_byte_estimate':int(14*86400*3*4*sum(float(plan[(plan.network+'.'+plan.station)==s].iloc[0].min_sampling_rate_hz) for s in cfg['stations'])),
        'gate':'Existing >=20 calibration candidates and maximum exposure-adjusted control ratio <=0.05; no evaluation scan on failure.',
        'sources':{str(p.resolve()):sha256(p) for p in [training_root/'report.json',Path(cfg['reference_csv'])]+[training_root/f'{s}_template.npz' for s in cfg['stations']]},
        'source_code':source_manifest(ROOT/'codex_lfe')}
    dp=out/'design_before_acquisition.json'
    if dp.exists():
        old=json.loads(dp.read_text())
        if old['config']!=cfg or old['sources']!=design['sources']:raise ValueError('Cannot resume changed inputs or config')
    else:
        write_json(dp,design);np.savez_compressed(out/'templates.npz',**archive)
    write_json(out/'manifest.json',{'config':cfg,'status':'running','design_sha256':sha256(dp)})
    for split in ['calibration','test']:
        if split=='test':
            selected=json.loads((out/'selection_before_test.json').read_text())['threshold']
            if selected is None:
                print('Calibration gate failed; evaluation remains sealed.',flush=True);break
            cfg['thresholds']=[selected]
        dates=cfg[split+'_dates']
        print('START',split,len(dates),'days',len(cfg['stations']),'stations',flush=True)
        tasks=[(s,m,d,Path(cfg['raw_root'])) for d in dates for s,m in cfg['stations'].items()]
        with ThreadPoolExecutor(max_workers=2) as pool:receipts=list(pool.map(fetch_day,tasks))
        write_json(out/f'{split}_acquisition.json',receipts)
        tasks=[(cfg,d,split,templates,out,'numpy','cpu',2) for d in dates if not (out/f'{split}_{d}_eligible_intervals.parquet').exists()]
        if tasks:
            with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
                for _,manifest in pool.map(scan_date_worker,tasks):
                    write_json(out/f'{split}_waveform_manifest_{len(list(out.glob(split+"_waveform_manifest_*.json"))):03d}.json',manifest)
        candidates,exposure,intervals=[pd.concat([pd.read_parquet(out/f'{split}_{d}_{suffix}.parquet') for d in dates],ignore_index=True) for suffix in ['candidates','exposure','eligible_intervals']]
        scores,eligible=summarize_split(cfg,split,candidates,exposure,intervals,reference,[family],{family:anchor_peak})
        write_json(out/f'{split}_scores.json',scores);eligible.to_parquet(out/f'{split}_eligible_reference.parquet',index=False)
        if split=='calibration':
            selected=operating_point(scores,cfg)
            write_json(out/'selection_before_test.json',{'threshold':selected,'scores_sha256':sha256(out/'calibration_scores.json'),'test_accessed':False})
        chosen=selected if selected is not None else cfg['diagnostic_threshold']
        c=candidates[(candidates.realization==0)&(candidates.threshold==chosen)].copy()
        c['validation_status']='experimental_unverified' if selected is not None else 'failed_calibration_diagnostic'
        c['arrival_anchor_utc']=pd.to_datetime(c.time_s,unit='s',utc=True).astype(str)
        c.to_csv(out/f'{split}_candidate_catalog.csv',index=False);c.to_parquet(out/f'{split}_candidate_catalog.parquet',index=False)
        print('COMPLETE',split,'selected threshold',selected,'candidate rows',len(c),'eligible refs',len(eligible),flush=True)
    write_json(out/'run_status.json',{'status':'complete','calibration_gate_passed':selected is not None,'test_completed':(out/'test_scores.json').exists(),'scientific_release_validated':False})

if __name__=='__main__':main()
