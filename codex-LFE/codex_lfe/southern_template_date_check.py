"""Check existing southern template families across nonoverlapping UTC-date halves."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json
import numpy as np
import pandas as pd
from obspy import read,Stream
from scipy.ndimage import gaussian_filter1d
from .templates import cut,stack_windows,vector_cc
from .waveforms import preprocess_stream
from .provenance import write_json,sha256


def worker(task):
    station,channels,output=task
    out=Path(output);done=out/f'{station}.json'
    if done.exists():return json.loads(done.read_text())
    original=Path('results/southern_training_complete_v1')
    previous=json.loads((original/f'{station}_result.json').read_text())
    ref=pd.read_parquet('results/southern_reference_audit_v1/published_reference.parquet').set_index('reference_id')
    if station.startswith('PB.'):
        paths=[Path('results/southern_training_multifamily_v1')/f'{station}.mseed']
    else:paths=sorted((original/station).glob('*.mseed'))
    stream=Stream()
    for p in paths:stream+=read(str(p))
    results=[]
    for family,ids in previous['valid_reference_ids'].items():
        halves=[[],[]];valid=[[],[]]
        for rid in ids:
            event=ref.loc[rid];t=float(event.time_s)
            # Keep raw/filter context away from UTC boundaries so halves cannot share samples.
            seconds=t%86400
            if seconds<120 or seconds>86280:continue
            half=int(np.floor(t/86400))%2
            wave=preprocess_stream(stream,channels,station,t-15,80,fs=20,alignment='sinc')
            packet=cut(wave,t-10,1400,alignment='sinc')
            if packet is not None:halves[half].append(packet);valid[half].append(rid)
        record=dict(station=station,family_id=family,half_event_counts=[len(x) for x in halves],previous_qc=previous['qc_by_family'][family],date_half_ids=valid)
        if min(len(x) for x in halves)<10:
            record.update(status='insufficient_independent_date_halves',usable=False)
        else:
            a,b=[stack_windows(x) for x in halves]
            energy=gaussian_filter1d(np.sum(a*a,axis=0),sigma=5)
            peak=200+int(np.argmax(energy[200:1200]))
            first=peak-40;last=first+160
            cc=vector_cc(a[:,first:last],b[:,first:last])
            record.update(status='checked',split_stack_cc=cc,usable=cc>=.25,peak_offset_from_published_s=-10+peak/20)
        results.append(record)
    result=dict(station=station,results=results,raw_sources={str(p):sha256(p) for p in paths},previous_result_sha256=sha256(original/f'{station}_result.json'))
    write_json(done,result);print(station,'date-split coherent families',sum(r['usable'] for r in results),flush=True)
    return result


def main():
    out=Path('results/southern_template_date_check_v1');out.mkdir(exist_ok=True)
    stations=pd.read_csv('results/catalog_regional_preflight_v1/southern_station_volume_plan.csv',dtype=str,keep_default_na=False)
    write_json(out/'design_before_check.json',dict(rule='UTC epoch-day parity, no event-lag realignment; exclude events within120s of UTC midnight; >=10 windows per half; choose 8s window on half A and test half B at same position, CC>=.25 as original check.',source_sha256=sha256(__file__),purpose='Assess alternate-event QC dependence on overlapping waveform intervals; no new downloads or production template replacement.'))
    with ProcessPoolExecutor(max_workers=2) as pool:
        records=list(pool.map(worker,[(f'{r.network}.{r.station}',r.channels.split(','),str(out)) for r in stations.itertuples()]))
    rows=[r for station in records for r in station['results']]
    pd.DataFrame([{k:v for k,v in r.items() if k not in ['previous_qc','date_half_ids']} for r in rows]).to_csv(out/'summary.csv',index=False)
    sites={f:[r['station'] for r in rows if r['family_id']==f and r['usable']] for f in sorted({r['family_id'] for r in rows})}
    write_json(out/'report.json',dict(status='complete_training_consistency_diagnostic',family_sites=sites,four_site_families=[f for f,s in sites.items() if len(s)>=4],production_templates_changed=False))
    print('Date-separated four-site families:',[f for f,s in sites.items() if len(s)>=4],flush=True)

if __name__=='__main__':main()
