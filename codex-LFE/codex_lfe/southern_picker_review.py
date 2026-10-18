"""Build blinded paired southern review windows from already scanned waveforms."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .waveforms import load_day
from .templates import cut
from .provenance import write_json,sha256


def main():
    root=Path('results/southern_picker_transfer_v1')
    design=json.loads((root/'design_before_inference.json').read_text())
    out=Path('results/southern_picker_review_v2');out.mkdir(exist_ok=False)
    refpath=Path('results/southern_reference_audit_v1/published_reference.parquet')
    ref=pd.read_parquet(refpath)
    families=json.loads(Path('results/southern_training_complete_v1/report.json').read_text())['four_site_families']
    rng=np.random.default_rng(20260914)
    windows=[];traces={};keys=[];manifests={};packets=[];used=[];exclusions=[]
    write_json(out/'design_before_sampling.json',dict(seed=20260914,dates=design['dates'],minimum_available_sites=3,sampling='One reference-enriched and one random block per date, without using model predictions; retain available station packets and explicitly record excluded stations.'))
    for date in design['dates']:
        day=pd.Timestamp(date,tz='UTC');waves={}
        for station,meta in design['stations'].items():
            raw=Path(design['raw_root'])/station/str(day.year)/f'{day.dayofyear:03d}.mseed'
            waves[station],manifests[f'{date}/{station}']=load_day(raw,Path('cache/southern_picker_waveforms'),meta['channels'],meta['site_id'],day.timestamp(),duration_s=86400)
        rr=ref[ref.date.eq(date)&ref.family_id.isin(families)]
        for sampling in ['reference_enriched','random']:
            if sampling=='reference_enriched':
                candidates=[(r.time_s-15+float(rng.uniform(-5,5)),r.reference_id) for r in rr.iloc[rng.permutation(len(rr))].itertuples()]
            else:candidates=[(day.timestamp()+float(rng.uniform(180,86220)),None) for _ in range(1000)]
            accepted=False
            for start,reference_id in candidates:
                if any(abs(start-x)<180 for x in used):continue
                data={s:cut(w,start-30,2400) for s,w in waves.items()}
                if sum(v is not None for v in data.values())<3:continue
                block=f'block_{len(used):03d}';used.append(start)
                for station,wave in waves.items():
                    if data[station] is None:
                        exclusions.append(dict(date=date,block_id=block,station=station,reason="missing_valid_waveform_packet"))
                        continue
                    wid=f'window_{len(windows):03d}'
                    row=dict(window_id=wid,station=station,start_s=start,end_s=start+60,split='development',sampling=sampling,review_status='pending',reviewer='',block_id=block,notes='',date=date)
                    windows.append(row);traces[wid]=data[station]
                    keys.append(dict(window_id=wid,sampling=sampling,source_reference_id=reference_id,block_id=block))
                    scale=float(np.quantile(abs(data[station]),.99)) or 1
                    packets.append(dict(window_id=wid,start_s=start,duration_s=60,channels=[c.split('.')[-1] for c in wave.channels],x=np.round(data[station]/scale,5).tolist()))
                accepted=True;break
            if not accepted:raise RuntimeError(f'No common valid paired block for {date}/{sampling}')
    order=rng.permutation(len(windows));windows=[windows[i] for i in order];packets=[packets[i] for i in order]
    np.savez_compressed(out/'waveforms.npz',**traces)
    write_json(out/'review_template.json',dict(windows=windows,labels=[]))
    write_json(out/'sampling_key.json',keys)
    write_json(out/'excluded_station_packets.json',exclusions)
    payload=dict(packet=packets,windows=windows,fs=20,context=30,key=sha256(out/'waveforms.npz'))
    template=Path('review/picker_review.html').read_text()
    (out/'review.html').write_text(template.replace('__PAYLOAD__',json.dumps(payload,separators=(',',':'))))
    write_json(out/'manifest.json',dict(status='unreviewed_development_packet',dates=design['dates'],fs=20,context_s=30,n_windows=len(windows),paired_time_blocks=len(used),waveform_inputs=manifests,source_sha256=sha256(__file__),reference_sha256=sha256(refpath),labels=0,scope='Three reference-enriched and three random paired time blocks at five southern roster stations, including surface and borehole instruments. Random windows are not labeled negative. Model outputs do not guide window selection. Independent review required.'))
    # Reuse existing prediction tables, without new inference or exposing scores in review UI.
    picks=[];coverage=[]
    for date in design['dates']:
        for station in design['stations']:
            w=[x for x in windows if x['date']==date and x['station']==station]
            folder=root/f'{date}_{station}'
            p=pd.read_parquet(folder/'picks.parquet');c=pd.read_parquet(folder/'coverage.parquet')
            for row in w:
                pp=p[(p.time_s>=row['start_s'])&(p.time_s<row['end_s'])].copy();pp['window_id']=row['window_id'];picks.append(pp)
            coverage.append(c)
    pd.concat(picks,ignore_index=True).to_parquet(out/'saved_model_picks.parquet',index=False)
    pd.concat(coverage,ignore_index=True).to_parquet(out/'saved_model_coverage.parquet',index=False)
    print(len(windows),'unlabeled paired station windows, existing predictions reused.',flush=True)

if __name__=='__main__':main()
