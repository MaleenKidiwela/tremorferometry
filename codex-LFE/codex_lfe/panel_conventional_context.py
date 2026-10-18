"""Retain conventional-earthquake context for preselected diagnostic hours."""
from pathlib import Path
import argparse,json,urllib.request,urllib.parse
import pandas as pd
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from .provenance import write_json,sha256


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--selection',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    out.mkdir(exist_ok=False);selection=json.loads(Path(args.selection).read_text())
    model=TauPyModel(model='iasp91');predictions=[];receipts=[]
    for panel,start in [('active',selection['start_utc']),('comparison',selection['quiet_comparison_start_utc'])]:
        start=pd.Timestamp(start);end=start+pd.Timedelta(hours=1)
        for region,extra in [('global',{'minmagnitude':5.5}),('regional',{'minlatitude':40,'maxlatitude':51,'minlongitude':-129,'maxlongitude':-119})]:
            params=dict(format='geojson',starttime=(start-pd.Timedelta(hours=1)).isoformat(),endtime=end.isoformat(),**extra)
            url='https://earthquake.usgs.gov/fdsnws/event/1/query?'+urllib.parse.urlencode(params)
            path=out/f'{panel}_{region}.geojson'
            try:
                with urllib.request.urlopen(url,timeout=45) as response:data=response.read(3_000_001)
                if len(data)>3_000_000:raise ValueError('Context response exceeds cap')
                obj=json.loads(data);path.write_bytes(data)
                receipts.append(dict(panel=panel,scope=region,url=url,sha256=sha256(path),events=len(obj['features']),status='received'))
                for event in obj['features']:
                    lon,lat,depth=event['geometry']['coordinates'];p=event['properties'];origin=pd.Timestamp(p['time'],unit='ms',tz='UTC')
                    distance=locations2degrees(lat,lon,selection['activity_center_latitude'],selection['activity_center_longitude'])
                    arrivals=model.get_travel_times(source_depth_in_km=max(0.,depth),distance_in_degree=distance,phase_list=['P','p','S','s'])
                    for a in arrivals:
                        t=origin+pd.Timedelta(seconds=a.time)
                        if start-pd.Timedelta(minutes=2)<=t<=end+pd.Timedelta(minutes=2):
                            predictions.append(dict(panel=panel,event_id=event['id'],magnitude=p['mag'],origin_utc=origin.isoformat(),phase=a.name,predicted_arrival_utc=t.isoformat(),distance_deg=distance,depth_km=depth,source_url=p['url']))
            except Exception as exc:
                receipts.append(dict(panel=panel,scope=region,url=url,status='context_query_failed',error=str(exc)))
    pd.DataFrame(predictions,columns=['panel','event_id','magnitude','origin_utc','phase','predicted_arrival_utc','distance_deg','depth_km','source_url']).to_csv(out/'arrival_context.csv',index=False)
    write_json(out/'report.json',dict(status='complete_context_audit' if all(r['status']=='received' for r in receipts) else 'incomplete_context_audit',queries=receipts,arrival_rows=len(predictions),
        limitations=['Iasp91 at the panel center is approximate; arrival flags are not automatic event exclusions.',
                    'Global query excludes magnitudes below 5.5 and events earlier than one hour before the panel.',
                    'No catalog match does not certify an LFE.'],selection_sha256=sha256(args.selection),runner_sha256=sha256(__file__)))
    print(json.dumps(receipts,indent=2));print('Arrival context rows:',len(predictions))

if __name__=='__main__':main()
