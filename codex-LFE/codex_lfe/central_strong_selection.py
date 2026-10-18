"""Select the strongest available 2026 central tremor day before waveform inspection."""
from pathlib import Path
import json,urllib.request,urllib.parse
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
from .provenance import write_json,sha256

def main():
    out=Path('results/central_strong_context_v1');out.mkdir(exist_ok=True)
    if (out/'selection_before_waveforms.json').exists():raise ValueError('Selection immutable')
    counts_path=Path('results/central_tremor_context_v1/central_day_counts.json');counts=json.loads(counts_path.read_text())
    day=max((d for d in counts if d.startswith('2026') and d<'2026-09-14'),key=lambda d:(counts[d],d))
    params=dict(starttime=day+'T00:00:00',endtime=(pd.Timestamp(day)+pd.Timedelta(days=1)).strftime('%Y-%m-%dT00:00:00'),lat_min=43.7,lat_max=45.,lon_min=-124.3,lon_max=-122.3)
    url='https://tremorapi.pnsn.org/api/v3.0/events?'+urllib.parse.urlencode(params)
    path=out/'events.geojson'
    if path.exists():
        receipt=json.loads((out/'receipt.json').read_text())
        if receipt['url']!=url or receipt['sha256']!=sha256(path):raise ValueError('Changed source metadata')
        data=path.read_bytes()
    else:
        with urllib.request.urlopen(url,timeout=45) as r:data=r.read(2_000_001)
        if len(data)>2_000_000:raise ValueError('Metadata response exceeds cap')
        path.write_bytes(data);write_json(out/'receipt.json',dict(url=url,bytes=len(data),sha256=sha256(path)))
    obj=json.loads(data);rows=[]
    for feature in obj['features']:
        props=feature['properties'];lon,lat=feature['geometry']['coordinates'][:2]
        rows.append(dict(time=pd.Timestamp(props['time']),latitude=lat,longitude=lon,energy=props['energy'],duration_s=props['duration']))
    frame=pd.DataFrame(rows);frame['hour']=frame.time.dt.floor('h');hour=frame.groupby('hour').energy.sum().idxmax();active=frame[frame.hour==hour]
    lat=float(active.latitude.median());lon=float(active.longitude.median())
    roster=pd.read_csv('results/target_network_from_map/stations.csv');candidates=[]
    for station in roster.itertuples():
        key=station.network+'.'+station.station;probe=Path('results/recent_waveform_probes_v1')/(key+'.json')
        if not probe.exists():continue
        record=json.loads(probe.read_text())
        # Accept only the previously successful, complete triplet probe records.
        if record.get('status')!='waveform_received' or not record.get('all_requested_channels_present'):continue
        candidates.append(dict(station=key,distance_km=gps2dist_azimuth(lat,lon,station.latitude,station.longitude)[0]/1000))
    candidates=pd.DataFrame(candidates).sort_values(['distance_km','station']);selected=candidates.head(8)
    if len(selected)!=8:raise ValueError('Fewer than eight eligible stations')
    candidates.to_csv(out/'station_distances.csv',index=False);frame.to_csv(out/'day_activity.csv',index=False)
    selection=dict(selection='Highest reported central tremor count day in cached 2026 metadata; ties choose latest date; highest total reported energy hour on that day; station choice before waveform inspection',activity_date=day,day_count=counts[day],selected_hour_features=len(active),start_utc=hour.isoformat(),core_duration_s=3600,activity_center_latitude=lat,activity_center_longitude=lon,stations=selected.station.tolist(),quiet_comparison_start_utc='2026-09-13T12:00:00Z',individual_LFE_truth=False,source_url=url,source_sha256=sha256(path),day_counts_sha256=sha256(counts_path),runner_sha256=sha256(__file__))
    write_json(out/'selection_before_waveforms.json',selection);print(json.dumps(selection,indent=2),flush=True)
if __name__=='__main__':main()
