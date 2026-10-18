"""Audit and normalize the local Shelly 2025 release without treating it as our detections."""
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from .provenance import write_json, sha256


def read_release(catalog, metadata):
    abstract=ET.parse(metadata).findtext('.//abstract')
    locations=re.findall(r'^\s*(\d+)\s+(40\.\d+)\s+(-123\.\d+)\s+(\d+\.\d+)\s*$',abstract,re.M)
    locations=pd.DataFrame(locations,columns=['idnum','latitude','longitude','depth_km']).astype({'idnum':int,'latitude':float,'longitude':float,'depth_km':float})
    if len(locations)!=27 or locations.idnum.duplicated().any():
        raise ValueError('Expected the 27 unique release family coordinates')
    raw=pd.read_csv(catalog).rename(columns={'%year':'year'})
    if raw.duplicated().any(): raise ValueError('Duplicate release rows')
    mapping=raw[['ID','idnum']].drop_duplicates()
    if mapping.ID.duplicated().any() or mapping.idnum.duplicated().any():
        raise ValueError('Family ID mapping is not one-to-one')
    seconds=raw.hour*3600+raw.minute*60+raw.second
    if not np.allclose(seconds,raw.s_of_day,atol=.011,rtol=0):
        raise ValueError('Calendar and seconds-of-day clocks disagree')
    time=pd.to_datetime(raw[['year','month','day']],utc=True)+pd.to_timedelta(seconds,unit='s')
    events=pd.DataFrame({'reference_id':[f'shelly2025_row_{i}' for i in range(len(raw))],
        'family_id':'shelly2025_'+raw.ID,'published_family_id':raw.ID,'idnum':raw.idnum,
        'time_s':time.astype('datetime64[ns, UTC]').astype('int64')/1e9,
        'time_basis':'shelly2025_template_start_first_station_estimated_S_minus_1s',
        'date':time.dt.strftime('%Y-%m-%d')})
    events=events.merge(locations,on='idnum',how='left',validate='many_to_one')
    if events.latitude.isna().any(): raise ValueError('Unmapped family location')
    return events,mapping.merge(locations,on='idnum',validate='one_to_one')


def main():
    catalog=Path('../lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv')
    metadata=Path('../lfe_catalog_repair/results/external/other/Shelly2025_metadata.xml')
    out=Path('results/southern_reference_audit_v1');out.mkdir(exist_ok=False)
    events,families=read_release(catalog,metadata)
    events.to_parquet(out/'published_reference.parquet',index=False)
    families.to_csv(out/'published_families.csv',index=False)
    events.groupby('date').size().rename('reference_count').to_csv(out/'daily_reference_counts.csv')
    write_json(out/'manifest.json',{'status':'external_reference_not_codex_detections','event_count':len(events),'family_count':len(families),
        'first_date':events.date.min(),'last_date':events.date.max(),
        'source_doi':'10.5066/P1TCKK7G','source_files':{str(p.resolve()):sha256(p) for p in [catalog,metadata,Path(__file__)]},
        'time_basis':events.time_basis.iloc[0],'origin_times_validated':False,'magnitudes_available':False,
        'checks':['No duplicate rows','One-to-one named-family mapping','27 metadata coordinates','Calendar versus seconds-of-day within 0.011 seconds']})
    print(len(events),'external references;',len(families),'named families;',events.date.min(),events.date.max())

if __name__=='__main__':main()
