"""Verify retained unit structure before deleting owned waveform staging.

This checks data integrity and traceability. Scientific validity is a separate gate;
passing these checks never promotes experimental detections into a final catalog.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .provenance import sha256, write_json

COVERAGE_STATUSES={'processed','no_metadata','provider_no_data','waveform_gap','insufficient_context'}


def verify_unit(artifacts, start_s, end_s, station_ids):
    """Require explicit, nonoverlapping coverage of the entire roster and unit."""
    if not np.isfinite([start_s,end_s]).all() or end_s<=start_s:
        raise ValueError('Invalid unit bounds')
    stations=set(station_ids)
    if not stations or any(not isinstance(x,str) or not x for x in stations):
        raise ValueError('Explicit station roster required')
    required={'events','families','coverage','controls','provenance'}
    if not required<=set(artifacts):raise ValueError('Missing products')
    paths={k:Path(artifacts[k]).resolve() for k in required}
    events=pd.read_parquet(paths['events']);families=pd.read_parquet(paths['families'])
    coverage=pd.read_parquet(paths['coverage']);controls=pd.read_parquet(paths['controls'])
    provenance=json.loads(paths['provenance'].read_text())
    for frame,columns,name in [
        (events,{'event_id','time_s','time_basis','family_id','quality_status','supporting_stations','support_windows'},'events'),
        (families,{'family_id'},'families'),
        (coverage,{'station','start_s','end_s','status','evidence_id'},'coverage'),
        (controls,{'realization','count','exposure_s'},'controls')]:
        if not columns<=set(frame):raise ValueError(f'Missing {name} columns')
    if families.family_id.isna().any() or families.family_id.duplicated().any():
        raise ValueError('Invalid family identities')
    if events.event_id.isna().any() or events.event_id.duplicated().any():
        raise ValueError('Invalid event identities')
    if not np.isfinite(events.time_s.to_numpy(float)).all() or not events.time_s.between(start_s,end_s,inclusive='left').all():
        raise ValueError('Event outside half-open core')
    if not set(events.family_id.dropna())<=set(families.family_id):
        raise ValueError('Unknown event family')
    if not provenance.get('time_basis') or not set(events.time_basis)<=set([provenance['time_basis']]):
        raise ValueError('Incompatible event clocks')
    if not set(events.quality_status)<=set(['experimental_unverified','reviewed_candidate']):
        raise ValueError('This integration path cannot certify final events')
    if set(coverage.station)!=stations or not set(coverage.status)<=COVERAGE_STATUSES:
        raise ValueError('Incomplete roster or unresolved coverage status')
    if not np.isfinite(coverage[['start_s','end_s']].to_numpy(float)).all():
        raise ValueError('Nonfinite coverage')
    evidence=provenance.get('coverage_evidence',{})
    if coverage.evidence_id.isna().any() or not set(coverage.evidence_id)<=set(evidence):
        raise ValueError('Coverage evidence missing')
    for key in set(coverage.evidence_id):
        item=evidence[key]
        if not isinstance(item,dict) or not item.get('description') or not item.get('source_sha256'):
            raise ValueError('Incomplete coverage evidence')
        digest=item['source_sha256']
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Invalid evidence checksum')
    for station,group in coverage.groupby('station'):
        g=group.sort_values('start_s')
        lo,hi=g.start_s.to_numpy(float),g.end_s.to_numpy(float)
        if (hi<=lo).any() or abs(lo[0]-start_s)>1e-6 or abs(hi[-1]-end_s)>1e-6:
            raise ValueError('Coverage does not span unit')
        if len(g)>1 and np.any(np.abs(lo[1:]-hi[:-1])>1e-6):
            raise ValueError('Coverage gap or overlap')
    for event in events.itertuples():
        support=event.supporting_stations
        if not isinstance(support,str):raise ValueError('Supporting station list required')
        selected=support.split(',')
        if not selected or len(set(selected))!=len(selected) or not set(selected)<=stations:
            raise ValueError('Invalid supporting stations')
        windows=json.loads(event.support_windows)
        if not isinstance(windows,list) or len(windows)!=len(selected) or {w['station'] for w in windows}!=set(selected):
            raise ValueError('Explicit waveform support windows required')
        for window in windows:
            lo,hi=window['start_s'],window['end_s']
            if not np.isfinite([lo,hi]).all() or hi<=lo:raise ValueError('Invalid support window')
            valid=coverage[(coverage.station==window['station'])&(coverage.status=='processed')].sort_values('start_s')
            cursor=lo
            for span in valid.itertuples():
                if span.start_s<=cursor<span.end_s:cursor=max(cursor,span.end_s)
                if cursor>=hi:break
            if cursor<hi:raise ValueError('Event support outside processed coverage')
    if controls.realization.duplicated().any() or not np.isfinite(controls[['realization','count','exposure_s']].to_numpy(float)).all():
        raise ValueError('Invalid control rows')
    if ((controls[['realization','count','exposure_s']]<0).any().any()
        or (controls.realization<1).any()
        or (controls.realization%1!=0).any() or (controls['count']%1!=0).any()):
        raise ValueError('Invalid control counts or exposure')
    if (coverage.status=='processed').any() and controls.empty:
        raise ValueError('Processed coverage requires retained controls')
    report=dict(passed=True,verification_kind='retained_product_integrity_only',
                scientific_validation_performed=False,production_validated=False,
                n_events=len(events),n_families=len(families),n_stations=len(stations),
                n_controls=len(controls),core_start_s=float(start_s),core_end_s=float(end_s),
                products={k:dict(path=str(p),sha256=sha256(p)) for k,p in paths.items()})
    return report


def write_verification(path, artifacts, start_s, end_s, station_ids):
    report=verify_unit(artifacts,start_s,end_s,station_ids)
    write_json(path,report)
    return report
