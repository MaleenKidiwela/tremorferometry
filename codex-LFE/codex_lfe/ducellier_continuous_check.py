"""Bounded one-hour development scan; published template clock, no phase claims."""
from pathlib import Path
import io,json,urllib.request,urllib.parse,argparse,math
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
from .provenance import write_json,sha256
from .published_pickle import read_published_pickle
from .neural import orientation_matrix
from .waveforms import preprocess_stream
from .matched_filter import normalized_cc,extract_candidates,shifted_sites

ROOT=Path('results/ducellier_continuous_development_v1')
START=UTCDateTime('2010-02-17T18:00:00'); DURATION=3600

def acquire(row):
    station=row.station; channels=row.channels.split(',')
    path=ROOT/'raw'/f'{station}.mseed'; receipt=path.with_suffix('.json')
    net,sta=station.split('.')
    host='https://service.earthscope.org' if net=='PB' else 'https://service.ncedc.org'
    params=dict(net=net,sta=sta,loc=channels[0].split('.')[2] or '--',cha=','.join(c.split('.')[-1] for c in channels),starttime=str(START-60),endtime=str(START+DURATION+95))
    url=host+'/fdsnws/dataselect/1/query?'+urllib.parse.urlencode(params)
    if receipt.exists() and path.exists():
        old=json.loads(receipt.read_text())
        if old.get('url')!=url:raise ValueError('Cached request differs; use a new output directory')
        if old.get('sha256')==sha256(path):return old
        raise ValueError('Cached waveform changed; preserve and inspect')
    rec=dict(station=station,channels=channels,url=url,path=str(path))
    try:
        cap=30_000_000*math.ceil(DURATION/3600)
        with urllib.request.urlopen(url,timeout=90) as r: data=r.read(cap+1)
        if not data or len(data)>cap:raise ValueError('empty or exceeds bounded request cap')
        read(io.BytesIO(data));path.write_bytes(data)
        rec.update(status='received',bytes=len(data),sha256=sha256(path))
    except Exception as e:rec.update(status='failed_not_evidence_of_absence',error=str(e))
    write_json(receipt,rec);print(station,rec['status'],flush=True)
    return rec

def main():
    global ROOT, START, DURATION
    parser=argparse.ArgumentParser()
    parser.add_argument('--start',default=str(START))
    parser.add_argument('--output',default=str(ROOT))
    parser.add_argument('--duration',type=int,default=3600)
    parser.add_argument('--acquire-only',action='store_true')
    args=parser.parse_args()
    ROOT=Path(args.output).resolve()
    if not ROOT.is_relative_to(Path.cwd().resolve()):raise ValueError('Output must stay inside codex-LFE')
    START=UTCDateTime(args.start)
    DURATION=args.duration
    if not 3600<=DURATION<=86400:raise ValueError('Duration must be one hour to one day')
    (ROOT/'raw').mkdir(parents=True,exist_ok=True)
    if (ROOT/'design.json').exists():
        previous=json.loads((ROOT/'design.json').read_text())
        if previous['start']!=str(START) or previous['duration_s']!=DURATION:
            raise ValueError('Existing design differs; use a new output directory')
        if (ROOT/'report.json').exists():raise ValueError('Completed scan is immutable; reuse its outputs')

    write_json(ROOT/'design.json',dict(status='development_only',source_sha256=sha256(__file__),start=str(START),duration_s=DURATION,family='080424.05.060',band_hz=[1.5,9],input='response-corrected velocity',threshold_mad=[6,8,10,12],min_sites=3,support_mad=2,controls='Eight nonwrapping independent station shifts; matched common exposure',selection='Bounded development/control interval; see start and duration. Threshold grid unchanged; no production validation claim',maximum_total_request_bytes=210000000*math.ceil(DURATION/3600)))
    plan=pd.read_csv('results/southern_benchmark_inventory_v1/targeted_reference_packet_plan.csv').drop_duplicates('station')
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(acquire,plan.itertuples()))
    write_json(ROOT/'acquisition.json',dict(records=records,total_bytes=sum(r.get('bytes',0) for r in records)))
    if args.acquire_only:return
    base=Path('results/ducellier_roster_packets_v1/responses')
    inv=read_inventory(str(base/'BK.xml'))+read_inventory(str(base/'NC.xml'))+read_inventory(str(base/'PB.xml'))
    tp=Path('results/southern_benchmark_inventory_v1/source/src/templates_v1/080424.05.060')
    sites={}; coverage=[]
    for r in records:
        if r['status']!='received':continue
        channels=r['channels']; st=read(r['path'])
        for tr in st:
            if str(inv.get_response(tr.id,tr.stats.starttime))!=str(inv.get_response(tr.id,tr.stats.endtime)):raise ValueError('response epoch crossing')
            tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80)
        templates=[];angles=[]
        for c in channels:
            tr,angle=read_published_pickle(tp/f'{r["station"].split(".")[1]}_{c.split(".")[-1]}.pkl')
            assert tr.stats.sampling_rate==20
            templates.append(tr.data);angles.append(angle)
        w=preprocess_stream(st,channels,r['station'],float(START),DURATION+36,fs=20,band=(1.5,9))
        current=orientation_matrix([inv.get_orientation(c,START) for c in channels])
        x=orientation_matrix(angles)@np.linalg.solve(current,w.x)
        # Every rotated component depends on the entire original triplet.
        valid=np.broadcast_to(w.valid.all(axis=0),x.shape)
        cc,good=normalized_cc(x,np.asarray(templates),valid)
        cc=cc[:DURATION*20];good=good[:DURATION*20]
        coverage.append(dict(station=r['station'],valid_seconds=float(good.sum()/20)))
        if good.any():sites[r['station']]=cc
        print(r['station'],'valid seconds',good.sum()/20,flush=True)
    np.savez_compressed(ROOT/'station_correlations.npz',**sites)
    write_json(ROOT/'coverage.json',coverage)
    if len(sites)<3:
        write_json(ROOT/'report.json',dict(status='insufficient_usable_sites',sites=list(sites),production_validated=False))
        raise RuntimeError('Fewer than three usable physical sites; no network scan')
    rng=np.random.default_rng(20260914)
    controls=[];shift_records=[]
    for i in range(8):
        shifts={s:int(rng.choice([-1,1])*rng.integers(60,601)*20) for s in sorted(sites)}
        controls.append(shifted_sites(sites,shifts));shift_records.append(shifts)
    common=np.isfinite(np.asarray(list(sites.values()))).all(axis=0)
    for c in controls:common &= np.isfinite(np.asarray(list(c.values()))).all(axis=0)
    rows=[];counts=[]
    for threshold in [6,8,10,12]:
        for i,panel in enumerate([sites]+controls):
            events,valid=extract_candidates(panel,20,float(START),threshold,min_sites=3,common_valid=common)
            rows.extend(dict(e,threshold_mad=threshold,control=i-1) for e in events)
            counts.append(dict(threshold_mad=threshold,control=i-1,candidates=len(events),exposure_s=float(valid.sum()/20)))
    pd.DataFrame(rows).to_csv(ROOT/'candidates.csv',index=False)
    pd.DataFrame(counts).to_csv(ROOT/'counts.csv',index=False)
    write_json(ROOT/'report.json',dict(status='complete_development_scan',sites=list(sites),control_shifts_samples=shift_records,common_exposure_s=float(common.sum()/20),production_validated=False,outputs_sha256={p.name:sha256(p) for p in [ROOT/'station_correlations.npz',ROOT/'candidates.csv',ROOT/'counts.csv',ROOT/'coverage.json']}))
    print(pd.DataFrame(counts).to_string(index=False),flush=True)

if __name__=='__main__':main()
