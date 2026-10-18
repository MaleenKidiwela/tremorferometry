"""Exercise verified monthly retention with existing, unvalidated August panels.

Only newly staged raw COPIES are deleted. Original evaluation data are preserved.
This is an eight-station diagnostic replay, never a full-month production scan.
"""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
from .monthly_storage import MonthlyStorage
from .provenance import sha256, write_json


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'results/monthly_panel_replay_v1'
    if out.exists(): raise ValueError('Replay output already exists')
    store=MonthlyStorage(out,25_000_000);store.begin_month('2026-08')
    units=[];originals={};summaries=[]
    for panel in ['quiet','active']:
        source=root/f'results/north_recent_{panel}_picker_v1'
        association=root/f'results/north_recent_{panel}_grid_controls_v1'
        design=json.loads((source/'design_before_acquisition.json').read_text())
        start=pd.Timestamp(design['start_utc']).timestamp();end=start+3600
        unit=pd.Timestamp(design['start_utc']).strftime('%Y%m%dT%H%M')
        units.append(unit);store.begin_unit(unit)
        retained=out/'retained'/unit;retained.mkdir(parents=True)
        artifacts={};evidence={};coverage=[]
        stations=pd.read_csv(source/'stations.csv').id.tolist()
        # Every retained file is a ledger artifact, not an unchecked manifest link.
        def retain(path, folder):
            destination=retained/folder/path.name;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,destination)
            artifacts[folder+'/'+path.name]=destination
            return destination
        staged_bytes=0
        for station in stations:
            raw=source/f'{station}.mseed';receipt=json.loads(raw.with_suffix('.mseed.json').read_text())
            digest=sha256(raw)
            if digest!=receipt['sha256'] or raw.stat().st_size!=receipt['bytes']:
                raise ValueError('Original acquisition receipt mismatch')
            originals[str(raw)]=digest
            store.stage(unit,raw.name,raw.read_bytes());staged_bytes+=raw.stat().st_size
            path=source/f'{station}_processed.npz'
            with np.load(path) as a:
                fs=float(a['fs']);i=int(round((start-float(a['start_s']))*fs));j=i+int(3600*fs)
                if not a['valid'][:,i:j].all() or not np.isfinite(a['probabilities'][:,i:j]).all():
                    raise ValueError('Incomplete waveform/prediction core')
                if a['waveform'][:,i:j].shape!=(3,int(3600*fs)):
                    raise ValueError('Wrong core extent')
            coverage.append(dict(station=station,start_s=start,end_s=end,status='processed',evidence_id=station))
            evidence[station]=dict(description='Retained processed waveform and finite P/S predictions inspected over complete one-hour core; original raw receipt hash verified.',source_sha256=sha256(path))
        for path in sorted(source.iterdir()):
            if path.is_file() and path.suffix!='.mseed':retain(path,'picker')
        for path in sorted(association.iterdir()):
            if path.is_file():retain(path,'association')
        for path in [root/'codex_lfe/monthly_storage.py',root/'codex_lfe/catalog_products.py',Path(__file__)]:retain(path,'code')
        for path in sorted((root.parent/'lfe_catalog_repair/results/external/other/seisbench_lfedetect').glob('jms.*')):
            if path.is_file():retain(path,'model')
        event_source=pd.read_parquet(association/'t0.3_n3_r0_events.parquet')
        assignment=pd.read_parquet(association/'t0.3_n3_r0_assignments.parquet');events=[]
        clock='homogeneous_model_origin_utc_unvalidated'
        for event in event_source.itertuples():
            picks=assignment[assignment.event_idx==event.idx]
            windows=[dict(station=s,start_s=float(g.time.min())-10,end_s=float(g.time.max())+10) for s,g in picks.groupby('station')]
            events.append(dict(event_id=f'diagnostic_{unit}_{event.idx}',time_s=float(event.time),time_basis=clock,
                family_id=None,quality_status='experimental_unverified',supporting_stations=','.join(w['station'] for w in windows),
                support_windows=json.dumps(windows),latitude=event.latitude,longitude=event.longitude,depth_km=event.depth))
        controls=pd.read_csv(association/'summary.csv');controls=controls[(controls.threshold==.3)&(controls.minimum_S_sites==3)&(controls.realization>0)]
        products=dict(events=pd.DataFrame(events,columns=['event_id','time_s','time_basis','family_id','quality_status','supporting_stations','support_windows','latitude','longitude','depth_km']),
            families=pd.DataFrame(columns=['family_id']),coverage=pd.DataFrame(coverage),controls=controls[['realization','candidates','exposure_s']].rename(columns={'candidates':'count'}))
        for name,frame in products.items():
            path=retained/(name+'.parquet');frame.to_parquet(path,index=False);artifacts[name]=path
        artifacts['provenance']=retained/'provenance.json';artifacts['verification']=retained/'verification.json'
        write_json(artifacts['provenance'],dict(time_basis=clock,coverage_evidence=evidence,
            scope='Eight-station one-hour diagnostic panel; remaining 175 roster stations were not part of this replay.',
            source_panel=str(source),source_association=str(association),production_validated=False,
            observed_threshold=.3,minimum_S_sites=3,control_exposure_s=3120,
            control_guard_note='Full waveform coverage is one hour; association and controls share a 3120 s interior.',
            raw_acquisition='Original cached receipt and bytes verified; no new network requests.',
            retention='Full processed waveform/prediction arrays, response XML, acquisition receipts, source picks, associations and all shifted controls retained; all files committed individually.'))
        report=store.commit_verified_unit(unit,artifacts,start,end,stations)
        # Exercise process restart from disk before cleanup.
        store=MonthlyStorage(out,25_000_000);store.cleanup(unit)
        for path,digest in originals.items():
            if sha256(path)!=digest:raise ValueError('Original evaluation waveform changed')
        summaries.append(dict(unit=unit,station_count=len(stations),experimental_events=len(events),
            staged_copy_bytes_deleted=staged_bytes,retained_file_count=len(artifacts),verification_passed=report['passed']))
        print(unit,summaries[-1],flush=True)
    store.finish_month(units)
    if any((out/'staged').rglob('*.mseed')):raise ValueError('Owned staging not empty')
    write_json(out/'report.json',dict(status='complete_real_data_retention_replay',production_validated=False,
        full_month_processed=False,full_roster_processed=False,raw_bytes_downloaded=0,
        original_waveforms_preserved=originals,units=summaries,
        note='Ledger completed_months applies only to this isolated two-unit diagnostic plan, not the production history queue.'))


if __name__=='__main__':main()
