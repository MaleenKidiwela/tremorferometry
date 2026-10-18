"""Acquire channel metadata for exactly the target roster; no waveform downloads."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
from pathlib import Path
import urllib.parse
import urllib.request
import pandas as pd
from obspy import read_inventory
from .inventory import triplet_epochs
from .roster import read_roster
from .provenance import sha256, write_json


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--roster',required=True)
    ap.add_argument('--start',required=True)
    ap.add_argument('--end',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args(); allowed=read_roster(args.roster);out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    groups={n:sorted(s.split('.')[1] for s in allowed if s.split('.')[0]==n) for n in sorted({s.split('.')[0] for s in allowed})}
    def fetch(item):
        net,stations=item
        endpoint='https://service.ncedc.org/fdsnws/station/1/query' if net in ('BK','NC') else 'https://service.earthscope.org/fdsnws/station/1/query'
        params={'network':net,'station':','.join(stations),'starttime':pd.Timestamp(args.start).isoformat(),'endtime':pd.Timestamp(args.end).isoformat(),'channel':'?H?','level':'channel','format':'xml'}
        url=endpoint+'?'+urllib.parse.urlencode(params);record={'network':net,'requested_stations':stations,'url':url}
        try:
            with urllib.request.urlopen(url,timeout=45) as response:data=response.read(20*1024*1024+1)
            if len(data)>20*1024*1024:raise ValueError('metadata response exceeds 20 MB cap')
            path=out/f'{net}.xml';path.write_bytes(data)
            inv=read_inventory(io.BytesIO(data))
            # Enforce roster even if a provider returns extra inventory entries.
            for network in inv:network.stations=[s for s in network if f'{network.code}.{s.code}' in allowed]
            table=triplet_epochs(inv,pd.Timestamp(args.start,tz='UTC').timestamp(),pd.Timestamp(args.end,tz='UTC').timestamp())
            record.update(status='metadata_received',sha256=sha256(path),bytes=len(data),n_triplet_epochs=len(table),returned_stations=sorted({f'{n.code}.{s.code}' for n in inv for s in n}))
            print(net,record['status'],len(table),'triplet epochs',flush=True)
        except Exception as exc:
            table=pd.DataFrame();record.update(status='request_failed_not_evidence_of_absence',error=f'{type(exc).__name__}: {exc}')
            print(net,record['status'],flush=True)
        write_json(out/f'{net}_request.json',record)
        return table,record
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(fetch,groups.items()))
    tables=[t for t,_ in results if not t.empty]
    table=pd.concat(tables,ignore_index=True) if tables else pd.DataFrame()
    table.to_csv(out/'triplet_epochs.csv',index=False)
    write_json(out/'manifest.json',{'roster_sha256':sha256(args.roster),'start':args.start,'end':args.end,
        'source_sha256':sha256(__file__),'n_requested_stations':len(allowed),'requests':[r for _,r in results],
        'interpretation':'Channel metadata only. No waveforms downloaded. Epoch alternatives overlap and are not additive. Metadata does not establish usable signal or catalog completeness.'})

if __name__=='__main__':main()
