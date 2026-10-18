"""Import the identified Chestler MAT layout while preserving unresolved clock/selection status."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.io import loadmat
from .provenance import sha256,write_json


def matlab_nominal_seconds(values):
    """MATLAB serial day to nominal Unix seconds; source timezone is not inferred."""
    values=np.asarray(values,dtype=np.float64).reshape(-1)
    if not np.isfinite(values).all():raise ValueError('Nonfinite source serial date')
    days=np.floor(values)
    return (days-719529)*86400+(values-days)*86400


def parse_reference(path):
    obj=loadmat(path,simplify_cells=True)
    if 'LFEs' not in obj:raise ValueError('Expected LFEs structure')
    raw=obj['LFEs'];items=[raw] if isinstance(raw,dict) else list(np.asarray(raw,dtype=object).reshape(-1))
    families=[];events=[];seen=set()
    for number,family in enumerate(items):
        if not isinstance(family,dict) or not {'name','peakTimesMo'}<=set(family):raise ValueError('Unknown MAT family schema')
        name=str(family['name']).strip()
        if not name or name in seen:raise ValueError('Missing or duplicate source family name')
        seen.add(name)
        dates=np.asarray(family['peakTimesMo'],dtype=np.float64).reshape(-1);seconds=matlab_nominal_seconds(dates)
        coords={k:float(np.asarray(family[k]).item()) if k in family else np.nan for k in ['lat','lon']}
        if np.isfinite(coords['lat']) and not -90<=coords['lat']<=90:raise ValueError('Invalid source latitude')
        if np.isfinite(coords['lon']) and not -180<=coords['lon']<=180:raise ValueError('Invalid source longitude')
        families.append(dict(source_family_id=name,source_row=number,latitude=coords['lat'],longitude=coords['lon'],source_detection_count=len(dates),final_published_selection_verified=False))
        for index,(serial,time_s) in enumerate(zip(dates,seconds)):
            events.append(dict(reference_id=f'chestler_source_row_{number:03d}_{index:07d}',source_family_id=name,source_matlab_datenum=serial,nominal_time_s_unverified=time_s,time_basis='peakTimesMo_MATLAB_datenum_nominal_UTC_unverified',quality_status='external_reference_clock_and_final_selection_unverified'))
    frame=pd.DataFrame(events,columns=['reference_id','source_family_id','source_matlab_datenum','nominal_time_s_unverified','time_basis','quality_status'])
    frame['nominal_time_utc_unverified']=pd.to_datetime(frame.nominal_time_s_unverified,unit='s',utc=True).astype(str)
    return pd.DataFrame(families),frame


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    source=Path(args.input);out=Path(args.output).resolve()
    if not out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output outside codex-LFE')
    families,events=parse_reference(source);out.mkdir(exist_ok=False)
    families.to_csv(out/'source_families.csv',index=False);events.to_parquet(out/'source_detections.parquet',index=False)
    write_json(out/'manifest.json',dict(status='imported_external_reference_unverified',input_path=str(source.resolve()),input_sha256=sha256(source),source_sha256=sha256(__file__),family_count=len(families),source_detection_count=len(events),same_family_exact_time_duplicates=int(events.duplicated(['source_family_id','source_matlab_datenum']).sum()),phase_labels_imported=0,time_convention_verified=False,final_published_family_selection_verified=False,production_validated=False))
    print((out/'manifest.json').read_text())
if __name__=='__main__':main()
