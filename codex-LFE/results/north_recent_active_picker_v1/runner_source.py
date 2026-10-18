"""Bounded current-network central Cascadia discovery-input panel; no accuracy claim."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import io,json,os,urllib.request,urllib.parse,argparse
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
from .provenance import write_json,sha256
from .neural import orientation_matrix,infer_probabilities,picks_from_probabilities,prediction_coverage
from .waveforms import preprocess_stream,Waveform

OUT=Path('results/central_recent_picker_panel_v1')
STATIONS=['UO.ALSE','UO.VITI','UO.BUCK','UO.BEER','UO.WLOO','UO.UCHR','UW.BABR','PB.B027']
START=UTCDateTime('2026-09-13T12:00:00')
REUSE=None

def fetch(station):
    probe=Path('results/recent_waveform_probes_v1')/(station+'.json')
    if not probe.exists():return dict(station=station,status='no_current_triplet_in_saved_metadata')
    prior=json.loads(probe.read_text())
    channels=prior['requested_channels'];net,sta=station.split('.')
    params=dict(net=net,sta=sta,loc=channels[0].split('.')[2] or '--',cha=','.join(c.split('.')[-1] for c in channels),starttime=str(START-120),endtime=str(START+3720))
    records={}
    for kind,service,extra,cap in [('raw','dataselect',{},30_000_000),('response','station',{'level':'response'},5_000_000)]:
        path=OUT/(station+('.mseed' if kind=='raw' else '.xml'))
        url='https://service.earthscope.org/fdsnws/'+service+'/1/query?'+urllib.parse.urlencode(params|extra)
        receipt=path.with_suffix(path.suffix+'.json')
        if receipt.exists() and path.exists():
            old=json.loads(receipt.read_text())
            if old['url']==url and old['sha256']==sha256(path):records[kind]=old;continue
            raise ValueError('Changed cached acquisition')
        if REUSE:
            old_path=REUSE/path.name;old_receipt=old_path.with_suffix(old_path.suffix+'.json')
            if old_path.exists() and old_receipt.exists():
                old=json.loads(old_receipt.read_text())
                if old.get('url')==url and old.get('sha256')==sha256(old_path):
                    old=old|dict(path=str(old_path),reused=True);records[kind]=old;continue
        with urllib.request.urlopen(url,timeout=60) as r:b=r.read(cap+1)
        if not b or len(b)>cap:raise ValueError('Empty or oversized response')
        (read if kind=='raw' else read_inventory)(io.BytesIO(b));path.write_bytes(b)
        record=dict(path=str(path),url=url,sha256=sha256(path),bytes=len(b));write_json(receipt,record);records[kind]=record
    print(station,'acquired',records['raw']['bytes'],'raw bytes',flush=True)
    return dict(station=station,status='received',channels=channels,**records)

def fetch_record(station):
    try:return fetch(station)
    except Exception as exc:
        record=dict(station=station,status='acquisition_failed_not_evidence_of_absence',error=f'{type(exc).__name__}: {exc}')
        write_json(OUT/(station+'_acquisition_failure.json'),record);print(station,record['status'],flush=True)
        return record

def main():
    global OUT, START, STATIONS, REUSE
    parser=argparse.ArgumentParser()
    parser.add_argument('--start',default=str(START));parser.add_argument('--output',default=str(OUT))
    parser.add_argument('--selection');parser.add_argument('--reuse')
    args=parser.parse_args();START=UTCDateTime(args.start);OUT=Path(args.output).resolve()
    if not OUT.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    REUSE=Path(args.reuse) if args.reuse else None
    selection=json.loads(Path(args.selection).read_text()) if args.selection else None
    if selection:STATIONS=selection['stations']
    OUT.mkdir(exist_ok=True)
    if (OUT/'report.json').exists():raise ValueError('Completed panel is immutable')
    (OUT/'runner_source.py').write_bytes(Path(__file__).read_bytes())
    roster=pd.read_csv('results/target_network_from_map/stations.csv')
    if not set(STATIONS).issubset(set(roster.network+'.'+roster.station)):raise ValueError('Station outside roster')
    write_json(OUT/'design_before_acquisition.json',dict(start_utc=str(START),core_duration_s=3600,requested_padding_s=120,stations=STATIONS,selection=selection or 'One fixed hour around the completed recent availability probe; no accuracy claim',maximum_raw_bytes=len(STATIONS)*30_000_000,model='jms',response_taper='Linear detrend and fixed 5-second Hann taper; remove_response taper=False; 60-second additional preprocessing padding',purpose='Bounded modern-roster discovery input panel; not production scanning or phase validation',source_sha256=sha256(__file__)))
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(fetch_record,STATIONS))
    write_json(OUT/'acquisition.json',dict(records=records,total_raw_bytes=sum(r.get('raw',{}).get('bytes',0) for r in records)))
    os.environ.setdefault('SEISBENCH_CACHE_ROOT',str(Path('cache/seisbench').resolve()))
    import torch,seisbench.models as sbm
    torch.set_num_threads(2)
    weights=Path('../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    model=sbm.LFEDetect.load(str(weights),map_location='cpu',weights_only=True);device='cuda' if torch.cuda.is_available() else 'cpu'
    stations=[];all_picks=[]
    for record in records:
        if record['status']!='received':continue
        station=record['station'];channels=record['channels'];stream=read(record['raw']['path']);inv=read_inventory(record['response']['path'])
        for tr in stream:
            if inv.get_orientation(tr.id,tr.stats.starttime)!=inv.get_orientation(tr.id,tr.stats.endtime):raise ValueError('Orientation epoch crossing')
            if str(inv.get_response(tr.id,tr.stats.starttime))!=str(inv.get_response(tr.id,tr.stats.endtime)):raise ValueError('Response epoch crossing')
            tr.detrend('linear');tr.taper(max_percentage=None,max_length=5,type='hann')
            tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80,taper=False)
        wave=preprocess_stream(stream,channels,station,float(START)-60,3720,fs=20,band=(1,8))
        matrix=orientation_matrix([inv.get_orientation(c,START) for c in channels]);zne=np.linalg.solve(matrix,wave.x).astype(np.float32)
        valid=np.broadcast_to(wave.valid.all(axis=0),zne.shape).copy();zne[~valid]=0
        wave=Waveform(zne,valid,wave.start_s,wave.fs,('Z','N','E'),station)
        prob,n=infer_probabilities(wave,model,device,batch_size=128)
        picks=picks_from_probabilities(prob,wave.start_s,wave.fs,station,.05)
        picks=picks[(picks.time_s>=float(START))&(picks.time_s<float(START)+3600)];all_picks.append(picks)
        picks.to_parquet(OUT/(station+'_picks.parquet'),index=False)
        prediction_coverage(prob,wave.start_s,wave.fs,station).to_parquet(OUT/(station+'_coverage.parquet'),index=False)
        np.savez_compressed(OUT/(station+'_processed.npz'),waveform=wave.x,valid=wave.valid,probabilities=prob,start_s=wave.start_s,fs=wave.fs)
        coordinates=inv.get_coordinates(channels[0],START)
        stations.append(dict(id=station,latitude=coordinates['latitude'],longitude=coordinates['longitude'],elevation=coordinates['elevation'],n_windows=n,n_picks_at_005=len(picks),n_picks_at_03=int((picks.model_score>=.3).sum())))
        print(station,'inferred',n,'windows;',len(picks),'peaks >= .05',flush=True)
    pd.DataFrame(stations).to_csv(OUT/'stations.csv',index=False);pd.concat(all_picks,ignore_index=True).to_parquet(OUT/'picks.parquet',index=False)
    write_json(OUT/'report.json',dict(status='complete_discovery_input_panel',waveform_bytes=sum(r.get('raw',{}).get('bytes',0) for r in records),model_weights_sha256=sha256(weights.with_suffix('.pt')),station_count=len(stations),phase_accuracy_measured=False,production_validated=False,outputs_sha256={p.name:sha256(p) for p in OUT.glob('*_processed.npz')}))

if __name__=='__main__':main()
