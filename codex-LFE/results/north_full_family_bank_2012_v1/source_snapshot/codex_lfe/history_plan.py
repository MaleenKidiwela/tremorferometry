"""Create a nonoverlapping all-history station-epoch work ledger from metadata."""
from pathlib import Path
import json
import sqlite3
import pandas as pd
from .provenance import sha256, write_json


def main():
    source=Path('results/full_history_metadata_v1/triplet_epochs.csv')
    roster=Path('results/target_network_from_map/stations.csv')
    raw=pd.read_csv(source,dtype={'location':str},keep_default_na=False)
    raw['lo']=pd.to_datetime(raw.start_utc,utc=True);raw['hi']=pd.to_datetime(raw.end_utc_exclusive,utc=True)
    rank={'BH':0,'HH':1,'EH':2,'SH':3,'VH':4}
    rows=[]
    for (network,station),g in raw.groupby(['network','station'],sort=True):
        edges=sorted(set(g.lo)|set(g.hi))
        for lo,hi in zip(edges[:-1],edges[1:]):
            active=g[(g.lo<=lo)&(g.hi>=hi)].copy()
            if active.empty:continue
            active['priority']=active.band_instrument.map(rank).fillna(10)
            selected=active.sort_values(['priority','orientation_condition','min_sampling_rate_hz','channels']).iloc[0]
            # Instrument alternatives must not multiply advertised station exposure.
            rate_bytes=selected.uncompressed_float32_bytes/selected.duration_days_metadata/86400
            row={'network':network,'station':station,'channels':selected.channels,'start_utc':lo.isoformat(),'end_utc_exclusive':hi.isoformat(),
                 'latitude':selected.latitude,'longitude':selected.longitude,'metadata_seconds':(hi-lo).total_seconds(),
                 'uncompressed_float32_bytes':int((hi-lo).total_seconds()*rate_bytes),'status':'awaiting_waveform_availability_and_regional_method'}
            rows.append(row)
    out=Path('results/full_history_work_plan_v1');out.mkdir(exist_ok=True)
    frame=pd.DataFrame(rows);frame.to_csv(out/'nonoverlapping_triplet_epochs.csv',index=False)
    conn=sqlite3.connect(out/'work.sqlite')
    conn.execute('CREATE TABLE IF NOT EXISTS waveform_blocks (id INTEGER PRIMARY KEY, network TEXT, station TEXT, channels TEXT, start_utc TEXT, end_utc_exclusive TEXT, status TEXT, error TEXT, UNIQUE(network,station,channels,start_utc,end_utc_exclusive))')
    jobs=[]
    for row in frame.itertuples():
        lo=pd.Timestamp(row.start_utc);end=pd.Timestamp(row.end_utc_exclusive)
        while lo<end:
            hi=min(end,lo.normalize()+pd.offsets.MonthBegin(1))
            jobs.append((row.network,row.station,row.channels,lo.isoformat(),hi.isoformat(),'awaiting_validation'))
            lo=hi
    conn.executemany('INSERT OR IGNORE INTO waveform_blocks(network,station,channels,start_utc,end_utc_exclusive,status) VALUES(?,?,?,?,?,?)',jobs)
    conn.commit();conn.close()
    table=frame.groupby(['network','station']).agg(first_metadata=('start_utc','min'),last_metadata=('end_utc_exclusive','max'),metadata_seconds=('metadata_seconds','sum'),uncompressed_float32_bytes=('uncompressed_float32_bytes','sum')).reset_index()
    table.to_csv(out/'station_history_summary.csv',index=False)
    write_json(out/'manifest.json',{'status':'planning_not_completed_processing','roster_sha256':sha256(roster),'metadata_sha256':sha256(source),'source_sha256':sha256(__file__),
        'stations':len(table),'nonoverlapping_instrument_epochs':len(frame),'monthly_blocks':len(jobs),'metadata_station_years':float(frame.metadata_seconds.sum()/86400/365.25),
        'uncompressed_float32_bytes':int(frame.uncompressed_float32_bytes.sum()),'instrument_preference':'BH then HH then EH then SH then VH; provisional planning choice, must validate sensitivity before release',
        'limits':['Metadata is not waveform availability','Single-component earlier epochs remain a separate validation requirement','Monthly blocks must be subdivided/padded for actual processing','No block is marked processed merely because it is in the ledger']})
    print((out/'manifest.json').read_text(),flush=True)

if __name__=='__main__':main()
