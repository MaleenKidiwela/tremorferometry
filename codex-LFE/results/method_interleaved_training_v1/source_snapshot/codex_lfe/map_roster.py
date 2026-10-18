"""Extract station membership from the user-designated dv/v HTML map.

Read its JSON DATA object without executing JavaScript or importing map builders.
The large embedded time-series payload is not needed for station membership.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .provenance import sha256, write_json


def read_map_stations(path):
    text=''
    with Path(path).open() as stream:
        while ',DAILY=' not in text:
            chunk=stream.read(1024*1024)
            if not chunk or len(text)>32*1024*1024:
                raise ValueError('map DATA payload missing or exceeds 32 MB')
            text+=chunk
    start=text.index('const DATA=')+len('const DATA=')
    data=json.loads(text[start:text.index(',DAILY=',start)])
    if not isinstance(data,dict) or not data:
        raise ValueError('empty or invalid station metadata')
    return {k:{'latitude':float(v['lat']),'longitude':float(v['lon'])} for k,v in data.items()}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--map',required=True)
    ap.add_argument('--lookup',required=True,nargs='+',help='network lookup CSV files with net,sta columns')
    ap.add_argument('--anchor-inventory',required=True,help='StationXML to resolve stations absent from lookup')
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    stations=read_map_stations(args.map)
    lookup=pd.concat([pd.read_csv(p)[['net','sta']] for p in args.lookup]).drop_duplicates()
    from obspy import read_inventory
    inv=read_inventory(args.anchor_inventory)
    rows=[]
    for name,coordinates in sorted(stations.items()):
        hits=sorted(lookup.loc[lookup.sta==name,'net'].unique())
        origin='lookup_exact_station_code'
        if not hits:
            hits=sorted({n.code for n in inv for s in n if s.code==name and
                np.hypot(s.latitude-coordinates['latitude'],s.longitude-coordinates['longitude'])<0.01})
            origin='stationxml_name_and_coordinates'
        if len(hits)!=1:raise ValueError(f'ambiguous or missing network for map station {name}: {hits}')
        rows.append({'network':hits[0],'station':name,**coordinates,'network_mapping':origin})
    table=pd.DataFrame(rows);table.to_csv(out/'stations.csv',index=False)
    write_json(out/'manifest.json',{'map_path':str(Path(args.map).resolve()),'map_sha256':sha256(args.map),
        'network_lookup_sha256':{str(Path(p).resolve()):sha256(p) for p in args.lookup},
        'anchor_inventory_sha256':sha256(args.anchor_inventory),'source_sha256':sha256(__file__),
        'n_station_markers':len(stations),'n_network_station_pairs':len(table),
        'authorization':'User identified borehole_dvv_map.html as containing the target stations.',
        'interpretation':'Membership extracted from current map. Map coordinates and lookup identities are not a substitute for independently verified instrument epochs and waveform availability. No LFE certification labels imported.'})
    print(f'{len(table)} stations extracted from the designated map')

if __name__=='__main__':main()
