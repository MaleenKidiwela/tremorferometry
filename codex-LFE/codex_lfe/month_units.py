"""Plan bounded network-coherent daily units; metadata does not imply availability."""
from pathlib import Path
import argparse,json
import pandas as pd
from .provenance import write_json,sha256


def plan_units(epochs,month,cutoff,padding_s=120):
    if padding_s<0:raise ValueError('Negative padding')
    cutoff=pd.Timestamp(cutoff)
    if cutoff.tzinfo is None:raise ValueError('UTC-aware cutoff required')
    begin=pd.Timestamp(month+'-01',tz='UTC');end=min(begin+pd.offsets.MonthBegin(1),cutoff)
    if end<=begin:raise ValueError('Month is after cutoff')
    epochs=epochs.copy();epochs['lo']=pd.to_datetime(epochs.start_utc,utc=True);epochs['hi']=pd.to_datetime(epochs.end_utc_exclusive,utc=True)
    for _,group in epochs.groupby(['network','station']):
        group=group.sort_values('lo')
        if (group.hi<=group.lo).any():raise ValueError('Nonpositive epoch')
        if len(group)>1 and (group.lo.iloc[1:].to_numpy()<group.hi.iloc[:-1].to_numpy()).any():raise ValueError('Overlapping station epochs')
    rows=[];units=[];hi=end
    while hi>begin:
        lo=max(begin,(hi-pd.Timedelta(nanoseconds=1)).normalize())
        unit_id=lo.strftime('%Y%m%d');units.append(dict(unit_id=unit_id,core_start_utc=lo.isoformat(),core_end_utc_exclusive=hi.isoformat()))
        for row in epochs[(epochs.lo<hi)&(epochs.hi>lo)].itertuples():
            core_lo=max(lo,row.lo);core_hi=min(hi,row.hi)
            request_lo=max(row.lo,core_lo-pd.Timedelta(seconds=padding_s));request_hi=min(row.hi,core_hi+pd.Timedelta(seconds=padding_s),pd.Timestamp(cutoff))
            rows.append(dict(unit_id=unit_id,network=row.network,station=row.station,channels=row.channels,core_start_utc=core_lo.isoformat(),core_end_utc_exclusive=core_hi.isoformat(),request_start_utc=request_lo.isoformat(),request_end_utc_exclusive=request_hi.isoformat(),left_padding_s=(core_lo-request_lo).total_seconds(),right_padding_s=(request_hi-core_hi).total_seconds(),status='planned_not_requested'))
        hi=lo
    return pd.DataFrame(units),pd.DataFrame(rows,columns=[
        'unit_id','network','station','channels','core_start_utc',
        'core_end_utc_exclusive','request_start_utc','request_end_utc_exclusive',
        'left_padding_s','right_padding_s','status'])


def metadata_gaps(units,requests,station_ids):
    """Explicit complementary intervals, including partial station days."""
    missing=[]
    for unit in units.itertuples():
        r=requests[requests.unit_id==unit.unit_id] if len(requests) else requests
        for station in sorted(station_ids):
            cursor=pd.Timestamp(unit.core_start_utc);end=pd.Timestamp(unit.core_end_utc_exclusive)
            selected=r[(r.network+'.'+r.station)==station].sort_values('core_start_utc') if len(r) else r
            spans=[(pd.Timestamp(x.core_start_utc),pd.Timestamp(x.core_end_utc_exclusive)) for x in selected.itertuples()]
            for lo,hi in spans+[(end,end)]:
                if cursor<lo:missing.append(dict(unit_id=unit.unit_id,station=station,core_start_utc=cursor.isoformat(),core_end_utc_exclusive=lo.isoformat(),status='no_selected_three_component_metadata_not_an_LFE_absence'))
                cursor=max(cursor,hi)
    return pd.DataFrame(missing,columns=['unit_id','station','core_start_utc','core_end_utc_exclusive','status'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--month',required=True);parser.add_argument('--cutoff');parser.add_argument('--output',required=True)
    args=parser.parse_args();source=Path('results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv');epochs=pd.read_csv(source)
    cutoff=args.cutoff or pd.to_datetime(epochs.end_utc_exclusive,utc=True).max().isoformat()
    out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    out.mkdir(exist_ok=False)
    units,requests=plan_units(epochs,args.month,cutoff);units.to_csv(out/'units.csv',index=False);requests.to_csv(out/'station_epoch_requests.csv',index=False)
    roster=pd.read_csv('results/target_network_from_map/stations.csv');all_ids=set(roster.network+'.'+roster.station)
    assert set(requests.network+'.'+requests.station)<=all_ids
    # Retain an explicit entry for stations without a selected epoch in every unit.
    missing=metadata_gaps(units,requests,all_ids)
    missing.to_csv(out/'stations_without_selected_epoch.csv',index=False)
    write_json(out/'manifest.json',dict(status='planned_not_downloaded',month=args.month,cutoff_utc=cutoff,units=len(units),station_epoch_requests=len(requests),missing_metadata_station_units=len(missing),padding_s=120,ordering='Descending UTC days; half-open event cores; padding may cross month boundaries but remains owned by the active month',source_sha256=sha256(__file__),epochs_sha256=sha256(source),production_validated=False,limits=['Selected triplet metadata is not waveform availability or sensitivity','No requests have been submitted','Cutoff cannot exceed refreshed metadata coverage','Epoch boundaries may truncate padding; record effective processing coverage after acquisition']))
    print((out/'manifest.json').read_text())
if __name__=='__main__':main()
