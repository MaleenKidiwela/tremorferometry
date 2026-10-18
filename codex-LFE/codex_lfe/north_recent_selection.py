"""Select a bounded recent northern activity/control panel before waveforms."""
from pathlib import Path
import json,urllib.request,urllib.parse
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
from .reuse_association import site_groups
from .provenance import sha256,write_json


def main():
    out=Path('results/north_recent_activity_context_v1')
    if (out/'selection_before_waveforms.json').exists():raise ValueError('Selection already frozen')
    counts_path=out/'day_counts.json';counts=json.loads(counts_path.read_text())
    day=max(d for d,n in counts.items() if d.startswith('2026') and d<'2026-09-14' and n>=20)
    params=dict(starttime=day+'T00:00:00',endtime=(pd.Timestamp(day)+pd.Timedelta(days=1)).strftime('%Y-%m-%dT00:00:00'),lat_min=48.1,lat_max=48.7,lon_min=-124.3,lon_max=-123.0)
    url='https://tremorapi.pnsn.org/api/v3.0/events?'+urllib.parse.urlencode(params)
    with urllib.request.urlopen(url,timeout=45) as response:data=response.read(2_000_001)
    if len(data)>2_000_000:raise ValueError('Oversized metadata')
    path=out/'active_day_events.geojson';path.write_bytes(data)
    write_json(out/'active_day_receipt.json',dict(url=url,sha256=sha256(path),bytes=len(data)))
    rows=[]
    for f in json.loads(data)['features']:
        p=f['properties'];lon,lat=f['geometry']['coordinates'][:2]
        rows.append(dict(time=pd.Timestamp(p['time']),latitude=lat,longitude=lon,energy=p['energy']))
    events=pd.DataFrame(rows);events['hour']=events.time.dt.floor('h')
    energy=events.groupby('hour').energy.sum().sort_index();hour=energy.idxmax();active=events[events.hour==hour]
    lat,lon=float(active.latitude.median()),float(active.longitude.median())
    # Nearest later day with no reported windows in the same calendar month.
    quiet_day=next(d.strftime('%Y-%m-%d') for d in pd.date_range(pd.Timestamp(day)+pd.Timedelta(days=1),pd.Timestamp(day)+pd.offsets.MonthEnd()) if counts.get(d.strftime('%Y-%m-%d'),0)==0)
    quiet=pd.Timestamp(quiet_day,tz='UTC')+pd.Timedelta(hours=hour.hour)
    roster=pd.read_csv('results/target_network_from_map/stations.csv');roster['id']=roster.network+'.'+roster.station
    groups=site_groups(roster,1.)
    candidates=[]
    for row in roster.itertuples():
        probe=Path('results/recent_waveform_probes_v1')/(row.id+'.json')
        if not probe.exists():continue
        p=json.loads(probe.read_text())
        if p.get('status')!='waveform_received' or not p.get('all_requested_channels_present'):continue
        candidates.append(dict(station=row.id,site=groups[row.id],latitude=row.latitude,longitude=row.longitude,distance_km=gps2dist_azimuth(lat,lon,row.latitude,row.longitude)[0]/1000))
    candidates=pd.DataFrame(candidates).sort_values(['distance_km','station'])
    selected=['CN.PGC','CN.VGZ']
    if not set(selected)<=set(candidates.station):raise ValueError('Required published overlap sites unavailable')
    selected_sites={groups[s] for s in selected}
    for row in candidates.itertuples():
        if row.site not in selected_sites:
            selected.append(row.station);selected_sites.add(row.site)
        if len(selected)==8:break
    if len(selected)!=8 or len(selected_sites)!=8:raise ValueError('Eight physical sites required')
    candidates.to_csv(out/'station_selection_candidates.csv',index=False);events.to_csv(out/'active_day_events.csv',index=False)
    selection=dict(selection='Latest cached 2026 northern day with at least 20 reported tremor windows; highest total-energy hour, earliest tie. Two current published-overlap sites plus six nearest distinct physical sites with successful recent triplet probes.',
        activity_date=day,day_count=counts[day],selected_hour_features=len(active),start_utc=hour.isoformat(),core_duration_s=3600,
        activity_center_latitude=lat,activity_center_longitude=lon,stations=selected,physical_site_groups={s:groups[s] for s in selected},
        quiet_comparison_start_utc=quiet.isoformat(),quiet_selection='First later date in the same month with no reported windows in the queried bounds; same UTC hour; this is not negative LFE truth.',
        execution_order='Quiet panel first, then activity panel; both in August 2026.',
        individual_LFE_truth=False,production_validated=False,source_sha256=sha256(path),day_counts_sha256=sha256(counts_path),source_sha256_code=sha256(__file__),
        limitations=['CN.SNB has no selected metadata after March 2021; never substitute a third published-overlap site.',
                    'CN.PGC changed from BH to HH channels; response and orientation correction are required.',
                    'This is a two-hour bounded transfer/discovery diagnostic, not full-archive production.'])
    write_json(out/'selection_before_waveforms.json',selection);print(json.dumps(selection,indent=2))

if __name__=='__main__':main()
